#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
using namespace Rcpp;

struct Entry {
  int taxon;
  double share;
  bool dominant;
};

using Sample = std::vector<Entry>;

struct Rectangle {
  std::vector<std::pair<int,double>> shared;
  std::vector<int> dominant;
  double richness, minimum, occupied, numerical, zeta;
  int dominant_count;
  double mean_dominance;
  bool defined;
};

double ratio(double x, double y) { return y > 0 ? x / y : NA_REAL; }

std::vector<Sample> communities(const NumericMatrix& a, const NumericMatrix& b,
                                const std::vector<std::vector<int>>& taxa,
                                const IntegerMatrix& mask, int coalition, double threshold) {
  std::vector<Sample> samples(a.nrow());
  for(int i=0;i<a.nrow();i++) {
    double total=0;
    auto& s=samples[i];
    for(int j:taxa[i]) {
      double q=(mask(i,j) & coalition) ? b(i,j) : a(i,j);
      if(q>0) { total+=q; s.push_back({j,q,false}); }
    }
    if(total==0) continue;
    std::vector<int> order(s.size());
    for(std::size_t j=0;j<s.size();j++) { s[j].share/=total; order[j]=j; }
    std::stable_sort(order.begin(),order.end(),[&](int x,int y){return s[x].share>s[y].share;});
    double cum=0;
    for(int j:order) { s[j].dominant=true; cum+=s[j].share; if(cum>=threshold-1e-12) break; }
  }
  return samples;
}

Rectangle rectangle(const std::vector<Sample>& samples,const IntegerVector& idx,int species) {
  Rectangle r{};
  std::vector<int> count(species,0),dom(species,0),touched;
  std::vector<double> sum(species,0.0),mn(species,R_PosInf);
  r.minimum=R_PosInf;
  r.defined=true;
  for(int id:idx) {
    const auto& s=samples[id-1];
    r.minimum=std::min(r.minimum,double(s.size()));
    r.richness+=s.size();
    r.defined &= !s.empty();
    for(const auto& v:s) {
      const int j=v.taxon;
      if(count[j]++==0) touched.push_back(j);
      sum[j]+=v.share;
      mn[j]=std::min(mn[j],v.share);
      if(v.dominant) {dom[j]++; r.mean_dominance++;}
    }
  }
  const double n=idx.size();
  r.richness/=n;
  r.mean_dominance/=n;
  // Sorted shared vectors facilitate a deterministic merge.
  std::sort(touched.begin(),touched.end());
  for(int j:touched) {
    const bool h=count[j]==n;
    if(h) {r.shared.push_back({j,mn[j]}); r.zeta+=mn[j];}
    if(dom[j]==n) r.dominant.push_back(j);
    r.occupied+=(count[j]/n-h)*mn[j];
    r.numerical+=std::max(0.0,sum[j]-count[j]*mn[j])/n;
  }
  r.dominant_count=r.dominant.size();
  return r;
}

template <typename A,typename B,typename Id,typename Weight>
double distance(const A& a,const B& b,Id id, Weight w) {
  std::size_t i=0,j=0;
  double sum=0,change=0;
  while(i<a.size() || j<b.size()) {
    if(j==b.size() || (i<a.size() && id(a[i])<id(b[j]))) {
      const double v=w(a[i++]); sum+=v; change+=v;
    } else if(i==a.size() || id(b[j])<id(a[i])) {
      const double v=w(b[j++]); sum+=v; change+=v;
    } else {
      const double x=w(a[i++]),y=w(b[j++]); sum+=x+y; change+=std::abs(x-y);
    }
  }
  return ratio(change,sum);
}

NumericMatrix response(const std::vector<Sample>& before,const std::vector<Sample>& after,
                       const std::vector<Rectangle>& baseline,const List& sets,int species) {
  NumericMatrix F(sets.size(),16);
  NumericMatrix local(before.size(),3);
  for(std::size_t i=0;i<before.size();i++) {
    const auto& a=before[i]; const auto& b=after[i];
    auto taxon=[](const Entry& z){return z.taxon;};
    local(i,0)=distance(a,b,taxon,[](const Entry&){return 1.0;});
    local(i,1)=(a.empty() || b.empty()) ? NA_REAL :
      distance(a,b,taxon,[](const Entry& z){return z.share;});
    local(i,2)=distance(a,b,taxon,[](const Entry& z){return double(z.dominant);});
  }
  for(int i=0;i<sets.size();i++) {
    IntegerVector idx=sets[i];
    const auto& a=baseline[i];
    const auto b=rectangle(after,idx,species);
    F(i,0)=b.shared.size();
    F(i,1)=b.defined ? b.zeta : NA_REAL;
    F(i,2)=b.dominant_count;
    F(i,3)=ratio(b.richness-b.shared.size(),b.richness);
    F(i,4)=ratio(b.minimum-b.shared.size(),b.minimum);
    F(i,5)=F(i,3)-F(i,4);
    F(i,6)=b.defined ? 1-b.zeta : NA_REAL;
    F(i,7)=b.defined ? b.occupied : NA_REAL;
    F(i,8)=b.defined ? b.numerical : NA_REAL;
    auto key=[](const std::pair<int,double>& x){return x.first;};
    F(i,9)=distance(a.shared,b.shared,key,[](const std::pair<int,double>&){return 1.0;});
    F(i,10)=(a.defined && b.defined) ?
      distance(a.shared,b.shared,key,[](const std::pair<int,double>& x){return x.second;}) : NA_REAL;
    for(int j=0;j<3;j++) {
      double d=0;
      for(int id:idx) d+=local(id-1,j);
      F(i,11+j)=d/idx.size();
    }
    F(i,14)=distance(a.dominant,b.dominant,[](int x){return x;},[](int){return 1.0;});
    F(i,15)=ratio(b.mean_dominance-b.dominant_count,b.mean_dominance);
  }
  return F;
}

// [[Rcpp::export]]
List st4_evaluate_cpp(NumericMatrix a,NumericMatrix b,List sets,List masks,
                      NumericMatrix weights,double threshold=0.8) {
  std::vector<std::vector<int>> taxa(a.nrow());
  for(int j=0;j<a.ncol();j++) for(int i=0;i<a.nrow();i++) if(a(i,j)>0 || b(i,j)>0) taxa[i].push_back(j);
  const int n=sets.size(), d=masks.size(), q=16;
  IntegerMatrix mask=masks[0];
  auto s0=communities(a,b,taxa,mask,0,threshold);
  auto s1=communities(a,b,taxa,mask,15,threshold);
  std::vector<Rectangle> baseline;
  for(int i=0;i<n;i++) baseline.push_back(rectangle(s0,as<IntegerVector>(sets[i]),a.ncol()));
  NumericMatrix f0=response(s0,s0,baseline,sets,a.ncol());
  NumericMatrix f1=response(s0,s1,baseline,sets,a.ncol());
  if(d==1 && weights.ncol()==0) return List::create(_["before"]=f0,_["after"]=f1);
  NumericVector phi(R_xlen_t(n)*q*4*d,NA_REAL);
  phi.attr("dim")=IntegerVector::create(n,q,4,d);
  NumericMatrix supported(n,q);
  double max_eff=0;
  for(int i=0;i<d;i++) {
    mask=as<IntegerMatrix>(masks[i]);
    std::vector<NumericMatrix> values(16);
    values[0]=f0; values[15]=f1;
    for(int s=1;s<15;s++) values[s]=response(s0,communities(a,b,taxa,mask,s,threshold),baseline,sets,a.ncol());
    for(int h=0;h<n;h++) for(int j=0;j<q;j++) {
      bool good=true;
      for(int s=0;s<16;s++) if(!std::isfinite(values[s](h,j))) {good=false;break;}
      if(!good) continue;
      supported(h,j)++;
      double total=0;
      for(int r=0;r<4;r++) {
        double x=0;
        for(int s=1;s<16;s++) if(weights(s,r)>0) {
          x+=weights(s,r)*(values[s](h,j)-values[s-(1<<r)](h,j));
        }
        phi[h+R_xlen_t(n)*(j+q*(r+4*i))]=x;
        total+=x;
      }
      max_eff=std::max(max_eff,std::abs(total-(f1(h,j)-f0(h,j))));
    }
  }
  return List::create(_["before"]=f0,_["after"]=f1,_["phi"]=phi,_["support"]=supported/d,
                      _["efficiency_error"]=max_eff);
}
