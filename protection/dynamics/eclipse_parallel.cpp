// Joint body/array occlusion at receiving-area quadrature points. Full and
// empty body masks retain the analytic parallel-square result.
#include "fast_parallel.cpp"
static bool sphere_hit(const double* d,const double* b,double r) {
  double dd=d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
  double db=d[0]*b[0]+d[1]*b[1]+d[2]*b[2];
  if(db<=0 || db>=dd) return false;
  double c0=b[1]*d[2]-b[2]*d[1],c1=b[2]*d[0]-b[0]*d[2],c2=b[0]*d[1]-b[1]*d[0];
  return c0*c0+c1*c1+c2*c2<=r*r*dd;
}
extern "C" void masked_moments(int n,const double* xyz,const double* sun,
 const double* body,const double* radii,double side,int grid,double* output) {
  parallel_moments(n,xyz,sun,side,output);
  double h=side/2,cell=side/grid,weight=cell*cell;
  std::vector<Rect> rect;rect.reserve(n);
  for(int i=0;i<n;i++) {
    double d[3]={sun[0]-xyz[3*i],sun[1]-xyz[3*i+1],sun[2]-xyz[3*i+2]};
    double ds=std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
    bool partial=false,full=false;
    for(int b=0;b<2;b++) {
      const double* v=body+6*i+3*b;
      double along=(v[0]*d[0]+v[1]*d[1]+v[2]*d[2])/ds;
      if(along<radii[b]+side || along>ds-radii[b]-side) continue;
      double c0=v[1]*d[2]-v[2]*d[1],c1=v[2]*d[0]-v[0]*d[2],c2=v[0]*d[1]-v[1]*d[0];
      double impact=std::sqrt(c0*c0+c1*c1+c2*c2)/ds;
      // Lipschitz enclosure for any point in the receiving square.
      double allowance=side/std::sqrt(2.)*(1+2*std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/(ds-side));
      if(impact<radii[b]-allowance) full=true;
      else if(impact<radii[b]+allowance) partial=true;
    }
    if(full) {output[3*i]=output[3*i+1]=output[3*i+2]=0;continue;}
    if(!partial) continue;
    rect.clear();double dzsun=d[2];
    for(int j=0;j<n;j++) {
      if(i==j)continue;
      double dz=xyz[3*j+2]-xyz[3*i+2],fraction=dz/dzsun;
      if(!(fraction>0&&fraction<1))continue;
      double stretch=dz/(sun[2]-xyz[3*j+2]),mag=1+stretch;
      double x=xyz[3*j]-xyz[3*i]+(xyz[3*j]-sun[0])*stretch;
      double y=xyz[3*j+1]-xyz[3*i+1]+(xyz[3*j+1]-sun[1])*stretch;
      double a=std::max(-h,x-h*mag),b=std::min(h,x+h*mag),c=std::max(-h,y-h*mag),e=std::min(h,y+h*mag);
      if(b>a&&e>c)rect.push_back({a,b,c,e});
    }
    double area=0,mx=0,my=0;
    for(int ix=0;ix<grid;ix++)for(int iy=0;iy<grid;iy++) {
      double x=-h+(ix+.5)*cell,y=-h+(iy+.5)*cell;
      double ray[3]={d[0]-x,d[1]-y,d[2]};bool blocked=false;
      for(int b=0;b<2;b++) {
        const double* v=body+6*i+3*b;double shifted[3]={v[0]-x,v[1]-y,v[2]};
        if(sphere_hit(ray,shifted,radii[b])){blocked=true;break;}
      }
      if(!blocked)for(auto r:rect)if(x>=r.a&&x<=r.b&&y>=r.c&&y<=r.d){blocked=true;break;}
      if(!blocked){area+=weight;mx+=weight*x;my+=weight*y;}
    }
    output[3*i]=area;output[3*i+1]=mx;output[3*i+2]=my;
  }
}
