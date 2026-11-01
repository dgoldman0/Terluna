// Continuous-mass, equal-particle snow column. Local implementation; see
// event_column_methods.md for the generator, thinning bound and split operators.
// No mass grid, statistical weights, collision batches or particle-size cap.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdio>

namespace snow {
using Vec = std::array<double,3>;
struct Particle { double m, z, d, v, growth; };
struct Physics {
    const double *p;
    Particle particle(double m,double z,const double *a) const {
        if (!(m>0) || !std::isfinite(m)) throw std::runtime_error("Invalid particle mass");
        double d=std::pow(m/p[0],1/p[1]);
        double x=a[0]*m, delta2=p[2]*p[2];
        double root=4*std::sqrt(x)/(delta2*std::sqrt(p[3]));
        // Rationalization preserves the same drag formula at small Best number.
        double diff=root/(std::sqrt(1+root)+1);
        double re=delta2/4*diff*diff;
        double v=a[1]*re/d, sc=a[3]*std::sqrt(re);
        double vent=sc<1 ? 1+.14*sc*sc : .86+.28*sc;
        double growth=a[2]*d*vent;
        if (!std::isfinite(v) || !std::isfinite(growth) || !(v>0) || growth<0)
            throw std::runtime_error("Invalid particle properties");
        return {m,z,d,v,growth};
    }
    double lost(const Particle&a,const Particle&b) const {
        double small=std::min(a.m,b.m);
        double d=std::clamp(std::pow(small/p[0],1/p[1]),p[5],p[6]);
        double scale=p[7]*(1+p[8]/std::pow(d,1.5))*std::acos(-1.)*d*d;
        double energy=(a.m/(a.m+b.m))*b.m*std::pow(a.v-b.v,2);
        double count=std::min(scale*(-std::expm1(-std::pow(p[9]*energy/scale,p[10]))),p[11]);
        return std::min(count*p[4],.5*small);
    }
};
struct RNG {
    std::mt19937_64 engine;
    explicit RNG(uint64_t s):engine(s){}
    double unit(){return std::generate_canonical<double,53>(engine);}
    double exponential(double rate){return -std::log1p(-unit())/rate;}
    size_t integer(size_t n){return std::uniform_int_distribution<size_t>(0,n-1)(engine);}
};
// A positive sum tree. Parent sums are recomputed from nonnegative children;
// no rate-subtraction drift or tiny-negative clipping enters selection.
struct Population {
    std::vector<Particle> x;
    std::vector<Vec> tree;
    size_t capacity=1;
    long double fragments=0; // unallocated count, in [0,1); retained in ledgers
    Population():tree(2,Vec{0,0,0}){}
    static Vec weights(const Particle&v){double d2=v.d*v.d;return {d2,v.v,d2*v.v};}
    void rebuild(){
        while(capacity<x.size()) capacity*=2;
        tree.assign(2*capacity,Vec{0,0,0});
        for(size_t i=0;i<x.size();i++) tree[capacity+i]=weights(x[i]);
        for(size_t i=capacity;i-->1;) for(int k=0;k<3;k++) tree[i][k]=tree[2*i][k]+tree[2*i+1][k];
    }
    void leaf(size_t i,Vec value){
        i+=capacity; tree[i]=value;
        while(i>1){i/=2;for(int k=0;k<3;k++) tree[i][k]=tree[2*i][k]+tree[2*i+1][k];}
    }
    void add(Particle value){
        x.push_back(value);
        if(x.size()>capacity) rebuild(); else leaf(x.size()-1,weights(value));
    }
    void set(size_t i,Particle value){x[i]=value;leaf(i,weights(value));}
    void erase(size_t i){
        size_t last=x.size()-1;
        if(i!=last){x[i]=x[last];leaf(i,weights(x[i]));}
        leaf(last,Vec{0,0,0});x.pop_back();
    }
    size_t select(int k,RNG&rng) const {
        double target=rng.unit()*tree[1][k];
        if(target>=tree[1][k]) target=std::nextafter(tree[1][k],0.);
        size_t at=1;
        while(at<capacity){double left=tree[2*at][k];if(target<left) at*=2;else {target-=left;at=2*at+1;}}
        size_t index=at-capacity;
        if(index>=x.size() || !(tree[at][k]>0)) throw std::runtime_error("Invalid sum-tree draw");
        return index;
    }
};
struct Stats {
    uint64_t proposals=0, sticking=0,chipping=0, emitted=0;
    long double fragment_number=0, deposited=0, injected=0,exited=0;
    uint64_t injected_n=0, exited_n=0,peak_particles=0,steps=0;
    double max_acceptance=0,max_diameter=0, fragment_growth_bound=0;
    std::chrono::steady_clock::time_point start=std::chrono::steady_clock::now();
    double wall_limit=1800;
    uint64_t proposal_limit=UINT64_C(10000000000),particle_limit=UINT64_C(10000000);
    void budget(){
        double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if(elapsed>wall_limit) throw std::runtime_error("Declared wall-time budget exceeded");
        if(proposals>proposal_limit) throw std::runtime_error("Declared proposal budget exceeded");
        if(peak_particles>particle_limit) throw std::runtime_error("Declared particle budget exceeded");
    }
};
void collide(Population &pop,const Physics&physics,const double *air,double volume,double duration,
             double C,bool breakup,RNG&rng,Stats&stats,int mode=0,double constant=0){
    double t=0;
    while(pop.x.size()>1){
        size_t n=pop.x.size();
        double first=n*pop.tree[1][2],second=pop.tree[1][0]*pop.tree[1][1];
        double multiplier=breakup ? 1 : air[4];
        double rate=mode==1 ? constant*double(n)*double(n-1)/(2*volume) : 2*C*(first+second)*multiplier/volume;
        if(!(rate>=0) || !std::isfinite(rate)) throw std::runtime_error("Invalid total event rate");
        if(rate==0) break;
        double wait=rng.exponential(rate);
        if(wait>duration-t) break;
        if(t+wait==t) throw std::runtime_error("Event time fell below floating-point resolution");
        t+=wait;stats.proposals++;
        if(stats.proposals%65536==0) stats.budget();
        size_t i,j;
        if(mode==1){i=rng.integer(n);j=rng.integer(n-1);if(j>=i)j++;}
        else if(rng.unit()*(first+second)<first){i=rng.integer(n);j=pop.select(2,rng);}
        else {i=pop.select(0,rng);j=pop.select(1,rng);}
        if(i==j) continue;
        Particle a=pop.x[i],b=pop.x[j];
        double bound=2*(a.d*a.d+b.d*b.d)*(a.v+b.v);
        double acceptance=mode==1 ? 1 : (a.d+b.d)*(a.d+b.d)*std::abs(a.v-b.v)/bound;
        if(!(acceptance>=0 && acceptance<=1+1e-14)) throw std::runtime_error("Invalid thinning probability");
        stats.max_acceptance=std::max(stats.max_acceptance,acceptance);
        if(mode!=1 && rng.unit()>=acceptance) continue;
        if(!breakup || rng.unit()<air[4]){
            // Retain the larger particle's position, the collecting trajectory.
            double z=a.m>=b.m ? a.z : b.z;
            Particle product=physics.particle(a.m+b.m,z,air);
            pop.set(i,product);pop.erase(j);stats.sticking++;
            stats.max_diameter=std::max(stats.max_diameter,product.d);
        }else{
            size_t small=a.m<=b.m ? i : j;
            Particle old=pop.x[small]; double loss=physics.lost(a,b);
            pop.set(small,physics.particle(old.m-loss,old.z,air));
            pop.fragments+=(long double)loss/physics.p[4];
            stats.fragment_number+=(long double)loss/physics.p[4];
            stats.chipping++;
            uint64_t count=static_cast<uint64_t>(std::floor(pop.fragments));
            pop.fragments-=count;
            for(uint64_t q=0;q<count;q++) pop.add(physics.particle(physics.p[4],old.z,air));
            stats.emitted+=count;
        }
        stats.peak_particles=std::max<uint64_t>(stats.peak_particles,pop.x.size());
    }
}
void statistics(const Population&pop,const double *query,int nq,double *moments,double*cdf){
    long double m1=0,m2=0;
    for(const Particle &v:pop.x){m1+=v.m;m2+=(long double)v.m*v.m;}
    moments[0]=pop.x.size();moments[1]=m1;moments[2]=m2;
    moments[3]=0;for(const Particle&v:pop.x)moments[3]=std::max(moments[3],v.m);
    std::fill(cdf,cdf+2*nq,0.);
    for(const Particle&v:pop.x){
        int index=std::lower_bound(query,query+nq,v.m)-query;
        if(index<nq){cdf[index]+=1;cdf[nq+index]+=v.m;}
    }
    for(int j=1;j<nq;j++){cdf[j]+=cdf[j-1];cdf[nq+j]+=cdf[nq+j-1];}
}
void grow(std::vector<Population>&cells,const Physics&physics,const double*air,double volume,
          double production,double duration,Stats&stats){
    for(size_t k=0;k<cells.size();k++){
        auto &pop=cells[k]; const double *s=air+5*k;
        // Bound includes empty cells whose allocation remainder is still present.
        stats.fragment_growth_bound+=duration*physics.particle(physics.p[4],0,s).growth;
        if(pop.x.empty()) continue;
        long double demand=0; for(const auto&x:pop.x)demand+=x.growth;
        double scale=demand>0 ? std::min((long double)1,production*volume/demand) : 0;
        std::vector<double> first(pop.x.size()), second(pop.x.size());
        long double d2=0;
        for(size_t i=0;i<pop.x.size();i++){
            first[i]=pop.x[i].growth*scale;
            second[i]=physics.particle(pop.x[i].m+duration*first[i],pop.x[i].z,s).growth;
            d2+=second[i];
        }
        double scale2=d2>0 ? std::min((long double)1,production*volume/d2) : 0;
        for(size_t i=0;i<pop.x.size();i++){
            auto old=pop.x[i];
            double add=duration*.5*(first[i]+scale2*second[i]);
            pop.x[i]=physics.particle(old.m+add,old.z,s);
            stats.deposited+=(long double)pop.x[i].m-old.m;
            stats.max_diameter=std::max(stats.max_diameter,pop.x[i].d);
        }
        pop.rebuild();
    }
}
void move(std::vector<Population>&cells,const Physics&physics,const double*air,double dz,double duration,
          Stats&stats,std::vector<double>&arrivals,bool gather,RNG&rng,bool upwind){
    std::vector<std::vector<Particle>> moved(cells.size());
    for(size_t k=0;k<cells.size();k++) for(auto x:cells[k].x){
        double remain=duration;int layer=k;
        while(layer>=0){
            double crossing=upwind ? rng.exponential(x.v/dz) : (x.z-layer*dz)/x.v;
            if(crossing>remain){if(!upwind)x.z-=x.v*remain;moved[layer].push_back(x);break;}
            remain-=crossing;x.z=layer*dz;layer--;
            if(layer<0){stats.exited+=x.m;stats.exited_n++;if(gather)arrivals.push_back(x.m);break;}
            x=physics.particle(x.m,x.z,air+5*layer);
        }
    }
    for(size_t k=0;k<cells.size();k++){cells[k].x.swap(moved[k]);cells[k].rebuild();}
}
void inject(std::vector<Population>&cells,const Physics&physics,const double*air,double depth,
            double rate,double duration,double lo,double hi,RNG&rng,Stats&stats){
    uint64_t n=std::poisson_distribution<uint64_t>(rate*duration)(rng.engine);
    for(uint64_t i=0;i<n;i++){
        double d=std::exp(std::log(lo)+rng.unit()*std::log(hi/lo));
        double m=physics.p[0]*std::pow(d,physics.p[1]);
        cells.back().add(physics.particle(m,depth,air+5*(cells.size()-1)));
        stats.injected+=m;stats.max_diameter=std::max(stats.max_diameter,d);
    }
    stats.injected_n+=n;
}
void error(char*out,const char*message){std::snprintf(out,1024,"%s",message);}
} // namespace snow

extern "C" int snow_properties(const double*p,const double*a,const double*m,int n,double*out,char*err){
    try{snow::Physics physics{p};for(int i=0;i<n;i++){auto x=physics.particle(m[i],0,a);out[3*i]=x.d;out[3*i+1]=x.v;out[3*i+2]=x.growth;}return 0;}
    catch(const std::exception&e){snow::error(err,e.what());return 1;}
}
extern "C" int snow_fragment_loss(const double*p,const double*a,const double*left,const double*right,int n,double*out,char*err){
    try{snow::Physics physics{p};for(int i=0;i<n;i++)out[i]=physics.lost(physics.particle(left[i],0,a),physics.particle(right[i],0,a));return 0;}
    catch(const std::exception&e){snow::error(err,e.what());return 1;}
}
extern "C" int snow_box(const double*p,const double*a,const double*m,int n,double volume,
                        const double*times,int nt,const double*query,int nq,uint64_t seed,int mode,
                        double constant,double C,double *moments,double*cdfs,double*info,char*err){
    try{
        snow::Physics physics{p};snow::RNG rng(seed);snow::Stats stats;snow::Population pop;
        for(int i=0;i<n;i++)pop.add(physics.particle(m[i],0,a));
        double before=0;
        for(int t=0;t<nt;t++){
            snow::collide(pop,physics,a,volume,times[t]-before,C,false,rng,stats,mode,constant);
            snow::statistics(pop,query,nq,moments+4*t,cdfs+2*nq*t);before=times[t];
        }
        info[0]=stats.proposals;info[1]=stats.sticking;info[2]=stats.max_acceptance;return 0;
    }catch(const std::exception&e){snow::error(err,e.what());return 1;}
}
// cfg: layers,area,dz,max_dt,fill_end,end,production,number_flux,seed_lo,seed_hi,
//      encounter_C,breakup,wall_limit,collisions,growth,settling,upwind_settling.
// output: 24 scalars, 4 cumulative arrival moments per query, 4 height moments
// per layer, and 4 time blocks of [number,mass,second moment,area] on arrival.
extern "C" int snow_column(const double*p,const double*air,const double*cfg,const double*query,int nq,
                           uint64_t seed,double*out,char*err){
    try{
        snow::Physics physics{p};snow::RNG rng(seed);snow::Stats stats;stats.wall_limit=cfg[12];
        int layers=int(cfg[0]);double area=cfg[1],dz=cfg[2],step=cfg[3],fill=cfg[4],end=cfg[5];
        double volume=area*dz,depth=dz*layers;
        std::vector<snow::Population> cells(layers);std::vector<double> arrivals;
        std::vector<long double> profiles(4*layers,0),blocks(16,0);
        double t=0;int phase=fill==0 ? 1 : 0;
        std::array<double,5> boundaries={fill,fill+(end-fill)/4,fill+(end-fill)/2,fill+3*(end-fill)/4,end};
        while(t<end){
            while(phase<4 && t>=boundaries[phase])phase++;
            stats.budget();bool gathering=phase>0;
            int block=phase-1;
            double boundary=boundaries[phase];
            double h=std::min(step,boundary-t);
            if(!(h>0)) throw std::runtime_error("Invalid split time boundary");
            size_t old_arrivals=arrivals.size();
            snow::inject(cells,physics,air,depth,cfg[7]*area,h/2,cfg[8],cfg[9],rng,stats);
            if(cfg[14])snow::grow(cells,physics,air,volume,cfg[6],h/2,stats);
            if(cfg[15])snow::move(cells,physics,air,dz,h/2,stats,arrivals,gathering,rng,bool(cfg[16]));
            if(cfg[13])for(int k=0;k<layers;k++)
                snow::collide(cells[k],physics,air+5*k,volume,h,cfg[10],bool(cfg[11]),rng,stats);
            if(cfg[15])snow::move(cells,physics,air,dz,h/2,stats,arrivals,gathering,rng,bool(cfg[16]));
            if(cfg[14])snow::grow(cells,physics,air,volume,cfg[6],h/2,stats);
            snow::inject(cells,physics,air,depth,cfg[7]*area,h/2,cfg[8],cfg[9],rng,stats);
            uint64_t live=0;
            for(int k=0;k<layers;k++){
                live+=cells[k].x.size();
                if(gathering){
                    for(const auto&x:cells[k].x){
                        profiles[4*k]+=h;profiles[4*k+1]+=(long double)h*x.m;
                        profiles[4*k+2]+=(long double)h*x.m*x.m;
                        profiles[4*k+3]+=(long double)h*cfg[10]*x.d*x.d;
                    }
                }
            }
            stats.peak_particles=std::max(stats.peak_particles,live);
            if(gathering)for(size_t i=old_arrivals;i<arrivals.size();i++){
                double m=arrivals[i],d=std::pow(m/p[0],1/p[1]);
                blocks[4*block]++;blocks[4*block+1]+=m;blocks[4*block+2]+=(long double)m*m;
                blocks[4*block+3]+=cfg[10]*d*d;
            }
            t+=h;stats.steps++;
        }
        long double remaining=0,remainder=0,count=0;
        for(const auto&cell:cells){
            for(const auto&x:cell.x)remaining+=x.m;
            count+=cell.x.size();remainder+=cell.fragments;
        }
        long double residual=stats.injected+stats.deposited-stats.exited-remaining-remainder*p[4];
        out[0]=stats.injected;out[1]=stats.deposited;out[2]=stats.exited;out[3]=remaining;
        out[4]=remainder*p[4];out[5]=residual;out[6]=stats.injected_n;out[7]=stats.exited_n;
        out[8]=count;out[9]=remainder;out[10]=stats.sticking;out[11]=stats.fragment_number;
        out[12]=stats.proposals;out[13]=stats.chipping;out[14]=stats.emitted;
        out[15]=stats.peak_particles;out[16]=stats.steps;out[17]=t;
        out[18]=stats.max_diameter;out[19]=stats.max_acceptance;out[20]=stats.fragment_growth_bound;
        std::sort(arrivals.begin(),arrivals.end());
        long double total=0;for(auto m:arrivals)total+=m;
        long double sum=0;int q=0;double quantiles[3]={.5,.9,.99};
        for(auto m:arrivals){sum+=m;while(q<3 && sum>=quantiles[q]*total){out[21+q]=std::pow(m/p[0],1/p[1]);q++;}}
        for(auto m:arrivals){
            int i=std::lower_bound(query,query+nq,m)-query;
            if(i<nq){double d=std::pow(m/p[0],1/p[1]);out[24+4*i]++;out[25+4*i]+=m;
                     out[26+4*i]+=m*m;out[27+4*i]+=cfg[10]*d*d;}
        }
        for(int j=1;j<nq;j++)for(int k=0;k<4;k++)out[24+4*j+k]+=out[24+4*(j-1)+k];
        for(int k=0;k<4*layers;k++)out[24+4*nq+k]=profiles[k];
        for(int k=0;k<16;k++)out[24+4*nq+4*layers+k]=blocks[k];
        return 0;
    }catch(const std::exception&e){snow::error(err,e.what());return 1;}
}
