#include "FvTestMesh2D.hpp"
#include "cartmesh2d/fv/HybridHeat2D.hpp"
#include <iostream>
#include <numeric>
#include <map>
using namespace cartmesh2d;using namespace cartmesh2d::fv;
namespace {
void require(bool ok,const char* s){if(!ok)throw std::runtime_error(s);}
// These are dimensional/scale-normalized algebra regressions, not physical gates.
constexpr double budget=32768*std::numeric_limits<double>::epsilon();
std::vector<HeatBoundary2D> boundaries(const FvMesh2D& m,bool fixed,double scale){std::vector<HeatBoundary2D>b;for(std::size_t f=0;f<m.faces.size();++f)if(!m.faces[f].neighbour)b.push_back({f,fixed?HeatBoundaryKind2D::Temperature:HeatBoundaryKind2D::Insulated,fixed?3+.1*m.faces[f].centre.x/scale-.2*m.faces[f].centre.y/scale:0,{}});return b;}
void affine(double scale) {
 auto m=fv_test::rectangle(7,5,2.3,true);for(auto& c:m.cells){c.centre.x*=scale;c.centre.y*=scale;c.area*=scale*scale;}for(auto& f:m.faces){f.centre.x*=scale;f.centre.y*=scale;f.areaVector.x*=scale;f.areaVector.y*=scale;f.correction.x*=scale;f.correction.y*=scale;}
 auto b=boundaries(m,true,scale);HybridHeatOperator2D op(m,b,.37);auto v=op.solveSteady();
 for(std::size_t i=0;i<m.cells.size();++i)require(std::abs(v.temperature[i]-(3+.1*m.cells[i].centre.x/scale-.2*m.cells[i].centre.y/scale))<budget*3,"HMM affine cell value");
 for(std::size_t f=0;f<m.faces.size();++f){double q=-.37*(.1*m.faces[f].areaVector.x-.2*m.faces[f].areaVector.y)/scale;require(std::abs(q-v.faceHeatFlux[f])<budget*(1+std::abs(q)),"HMM affine face flux");}
 std::map<std::pair<std::size_t,std::size_t>,double> a;op.visitSteadyTraceMatrix([&](auto i,auto j,double x){a[{i,j}]=x;});for(const auto& [ij,x]:a)require(x==a[{ij.second,ij.first}],"HMM trace symmetry");
 // Cholesky each ACTUAL local matrix: failure is not repaired by shifting.
 std::map<std::size_t,std::map<std::pair<std::size_t,std::size_t>,double>> locals;op.visitLocalMatrices([&](auto c,auto i,auto j,double x){locals[c][{i,j}]=x;});
 for(auto& [c,mat]:locals){auto n=m.cells[c].faces.size();std::vector<double> l(n*n);for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<=i;++j){double q=mat[{i,j}];for(std::size_t h=0;h<j;++h)q-=l[i*n+h]*l[j*n+h];if(i==j){require(q>0,"HMM local energy not positive");l[i*n+j]=std::sqrt(q);}else l[i*n+j]=q/l[j*n+j];}}
 auto fixed=op.evaluateAtCells(v.temperature);for(std::size_t i=0;i<v.temperature.size();++i)require(std::abs(fixed.cellResidual[i]-v.cellResidual[i])<budget,"static trace/cell condensation disagree");
 auto original=v;for(auto& f:m.faces)if(f.neighbour){auto owner=f.owner;f.owner=*f.neighbour;f.neighbour=owner;f.areaVector={-f.areaVector.x,-f.areaVector.y};f.correction={-f.correction.x,-f.correction.y};f.neighbourWeight=1-f.neighbourWeight;}
 v=HybridHeatOperator2D(m,b,.37).solveSteady();for(std::size_t i=0;i<v.temperature.size();++i)require(std::abs(v.temperature[i]-original.temperature[i])<budget*3,"HMM owner flip changed solution");
}
void closed() {
 auto m=fv_test::rectangle(8,6,2,true);HybridHeatOperator2D op(m,boundaries(m,false,1),.37);std::vector<double> t(m.cells.size()),cv(t.size(),1.3);
 for(std::size_t i=0;i<t.size();++i)t[i]=3+.1*std::sin(double(i));
 auto actual=op.evaluateAtCells(t);double power=0;for(std::size_t i=0;i<t.size();++i)power+=t[i]*actual.cellResidual[i];
 for(std::size_t f=0;f<m.faces.size();++f)if(!m.faces[f].neighbour)power-=actual.trace[f]*actual.faceHeatFlux[f];
 require(power>0&&std::abs(power-actual.dissipation)<budget,"HMM global energy/flux identity");
 auto initial=t;double mean=0,volume=0;for(std::size_t i=0;i<t.size();++i){mean+=m.cells[i].area*t[i];volume+=m.cells[i].area;}mean/=volume;
 double energy=0;for(std::size_t i=0;i<t.size();++i)energy+=cv[i]*m.cells[i].area*(t[i]-mean)*(t[i]-mean)/2;
 for(int step=0;step<8;++step){auto v=op.backwardEuler(t,cv,.01);require(v.maximumStageToFinalBudgetRatio<=1,"HMM stage/final algebraic budget");for(std::size_t i=0;i<t.size();++i)require(v.temperature[i]==t[i]-.01*v.cellResidual[i]/(cv[i]*m.cells[i].area)||std::abs(v.temperature[i]-(t[i]-.01*v.cellResidual[i]/(cv[i]*m.cells[i].area)))<budget,"HMM final not shared-flux update");double next=0,heat=0;for(std::size_t i=0;i<t.size();++i){next+=cv[i]*m.cells[i].area*(v.temperature[i]-mean)*(v.temperature[i]-mean)/2;heat+=cv[i]*m.cells[i].area*(v.temperature[i]-initial[i]);}require(next<=energy+budget*energy,"HMM closed perturbation energy increased");require(std::abs(heat)<budget,"HMM closed heat not conserved");require(v.maximumFaceFluxJump<budget,"HMM action/reaction mismatch");energy=next;t=v.temperature;}
 bool failed=false;try{(void)op.solveSteady();}catch(const std::exception&){failed=true;}require(failed,"pure Neumann steady silently fixed");
 failed=false;auto old=t;try{(void)op.backwardEuler(t,cv,0);}catch(const std::exception&){failed=true;}require(failed&&old==t,"failed step mutated accepted state");
 m.cells[0].centre.x-=5;failed=false;try{HybridHeatOperator2D bad(m,boundaries(m,false,1),.37);}catch(const std::exception&){failed=true;}require(failed,"invalid geometric kernel accepted");
}
}
int main(){try{for(double scale:{1.,1e-4,1e4})affine(scale);closed();std::cout<<"HMM affine/scale/local SPD/shared flux/energy/transaction regressions passed\n";return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
