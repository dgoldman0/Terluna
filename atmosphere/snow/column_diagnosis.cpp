// Diagnostics around the frozen event-column physics. The included file owns
// particle properties, collisions, growth, injection and transport. This driver
// retains its operator order and output arithmetic; tests require identical
// scientific output at fixed seed. Observation consumes no random numbers.
#include "event_column.cpp"

namespace diagnosis {
constexpr int WIDTH=15;
struct Demand {long double sum=0,square=0;double largest=0;};
void add(Demand&d,double g){d.sum+=g;d.square+=(long double)g*g;d.largest=std::max(d.largest,g);}
void sample(size_t n,const Demand&d,double capacity,double dt,long double*out){
    out[0]+=dt;
    if(n==0)out[1]+=dt;
    if(n==1)out[2]+=dt;
    if(d.sum>capacity)out[3]+=dt;
    out[4]+=dt*d.sum;out[5]+=dt*std::min(d.sum,(long double)capacity);
    if(n==0)out[6]+=(long double)dt*capacity;
    else out[7]+=dt*std::max((long double)0,capacity-d.sum);
    out[8]+=(long double)dt*n;out[9]+=(long double)dt*n*n;
    if(d.square>0)out[10]+=dt*d.sum*d.sum/d.square;
    if(d.sum>0)out[11]+=dt*d.largest/d.sum;
}
void observe(const std::vector<snow::Population>&cells,const snow::Physics&physics,
             const double*air,double volume,double production,double duration,
             bool enabled,long double*out){
    double capacity=enabled ? production*volume : 0;
    for(size_t k=0;k<cells.size();k++){
        const auto&pop=cells[k];auto dest=out+k*WIDTH;
        Demand first,second;
        if(enabled){
            for(const auto&x:pop.x)add(first,x.growth);
            double scale=first.sum>0 ? std::min((long double)1,production*volume/first.sum) : 0;
            for(const auto&x:pop.x){
                double increment=x.growth*scale;
                add(second,physics.particle(x.m+duration*increment,x.z,air+5*k).growth);
            }
        }
        sample(pop.x.size(),first,capacity,duration/2,dest);
        sample(pop.x.size(),second,capacity,duration/2,dest);
    }
}
} // namespace diagnosis

// The original output array is unchanged. trace[phase,layer,field] records
// filling then four gathering blocks. Growth statistics use both Heun stages
// with their exact quadrature weights, not exact continuous residence times.
extern "C" int snow_column_diagnosis(const double*p,const double*air,const double*cfg,const double*query,int nq,
                           uint64_t seed,double*out,double*trace,char*err){
    try{
        snow::Physics physics{p};snow::RNG rng(seed);snow::Stats stats;stats.wall_limit=cfg[12];
        int layers=int(cfg[0]);double area=cfg[1],dz=cfg[2],step=cfg[3],fill=cfg[4],end=cfg[5];
        double volume=area*dz,depth=dz*layers;
        std::vector<snow::Population> cells(layers);std::vector<double> arrivals;
        std::vector<long double> profiles(4*layers,0),blocks(16,0);
        std::vector<long double> audit(5*layers*diagnosis::WIDTH,0);
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
            diagnosis::observe(cells,physics,air,volume,cfg[6],h/2,bool(cfg[14]),audit.data()+phase*layers*diagnosis::WIDTH);
            if(cfg[14])snow::grow(cells,physics,air,volume,cfg[6],h/2,stats);
            if(cfg[15])snow::move(cells,physics,air,dz,h/2,stats,arrivals,gathering,rng,bool(cfg[16]));
            if(cfg[13])for(int k=0;k<layers;k++){
                uint64_t sticking=stats.sticking,chipping=stats.chipping,emitted=stats.emitted;
                snow::collide(cells[k],physics,air+5*k,volume,h,cfg[10],bool(cfg[11]),rng,stats);
                auto dest=audit.data()+(phase*layers+k)*diagnosis::WIDTH;
                dest[12]+=stats.sticking-sticking;dest[13]+=stats.chipping-chipping;dest[14]+=stats.emitted-emitted;
            }
            if(cfg[15])snow::move(cells,physics,air,dz,h/2,stats,arrivals,gathering,rng,bool(cfg[16]));
            diagnosis::observe(cells,physics,air,volume,cfg[6],h/2,bool(cfg[14]),audit.data()+phase*layers*diagnosis::WIDTH);
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
        for(size_t k=0;k<audit.size();k++)trace[k]=audit[k];
        return 0;
    }catch(const std::exception&e){snow::error(err,e.what());return 1;}
}
