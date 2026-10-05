// Exact parallel-square shadow union and first moments; no geometric sampling.
#include <algorithm>
#include <cmath>
#include <vector>
struct Rect { double a,b,c,d; };
extern "C" void parallel_moments(int n, const double* xyz, const double* sun,
                                  double side, double* output) {
  const double h=side/2;
  std::vector<Rect> rect; std::vector<double> xs;
  std::vector<std::pair<double,double>> ys;
  rect.reserve(n); xs.reserve(2*n); ys.reserve(n);
  for(int i=0;i<n;i++) {
    rect.clear(); xs.clear(); bool full=false;
    double dzsun=sun[2]-xyz[3*i+2];
    for(int j=0;j<n;j++) {
      if(i==j) continue;
      double dz=xyz[3*j+2]-xyz[3*i+2];
      double fraction=dz/dzsun;
      if(!(fraction>0 && fraction<1)) continue;
      double stretch=dz/(sun[2]-xyz[3*j+2]), mag=1+stretch;
      double x=xyz[3*j]-xyz[3*i]+(xyz[3*j]-sun[0])*stretch;
      double y=xyz[3*j+1]-xyz[3*i+1]+(xyz[3*j+1]-sun[1])*stretch;
      double a=std::max(-h,x-h*mag), b=std::min(h,x+h*mag);
      double c=std::max(-h,y-h*mag), d=std::min(h,y+h*mag);
      if(b<=a || d<=c) continue;
      if(a==-h && b==h && c==-h && d==h) {full=true;break;}
      rect.push_back({a,b,c,d}); xs.push_back(a);xs.push_back(b);
    }
    if(full) {output[3*i]=output[3*i+1]=output[3*i+2]=0;continue;}
    std::sort(xs.begin(),xs.end());
    double area=0,mx=0,my=0;
    for(size_t k=1;k<xs.size();k++) {
      double w=xs[k]-xs[k-1]; if(w<=0) continue;
      double x=(xs[k]+xs[k-1])/2; ys.clear();
      for(auto r:rect) if(x>r.a && x<r.b) ys.push_back({r.c,r.d});
      std::sort(ys.begin(),ys.end()); double top=-h,height=0,moment=0;
      for(auto y:ys) {
        double low=std::max(top,y.first);
        if(y.second>low) {height+=y.second-low;moment+=(y.second*y.second-low*low)/2;}
        top=std::max(top,y.second);
      }
      area+=w*height;mx+=w*x*height;my+=w*moment;
    }
    output[3*i]=std::max(0.,std::min(side*side,side*side-area));
    output[3*i+1]=-mx; output[3*i+2]=-my;
  }
}
