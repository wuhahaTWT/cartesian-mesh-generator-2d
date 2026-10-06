// Complete native Euler/Corrected annular Couette verification, not product moving-wall qualification.
#include "cartmesh2d/fv/Euler2D.hpp"
#include "cartmesh2d/fv/EulerCheckpoint2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <chrono>
#include <fstream>
#include <filesystem>
#include <csignal>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
  volatile std::sig_atomic_t stopRequested=0;
  void requestStop(int){stopRequested=1;}
  bool thermalReference=false;
  constexpr double pi=3.14159265358979323846;
  const IdealGas2D gas{
    1.4,1
  }
  ;
  const EulerTransport2D tr{
    .02*3.5/.72,.02
  }
  ;
  const double G=.02*.2*.2/tr.thermalConductivity,C=(.2+G*(.25-1))/std::log(2.),D=1+G;
  struct Q{
    std::vector<double>x,w;
    Q(int n){
      for(int i=0;
      i<n;
      ++i){
        double z=std::cos(pi*(i+.75)/(n+.5)),dp=0;
        for(int it=0;
        it<30;
        ++it){
          double a=1,b=0;
          for(int k=1;
          k<=n;
          ++k){
            double t=a;
            a=((2*k-1)*z*a-(k-1)*b)/k;
            b=t;
          }
          dp=n*(z*a-b)/(z*z-1);
          double next=z-a/dp;
          if(std::abs(next-z)<2e-16){
            z=next;
            break;
          }
          z=next;
        }
        x.push_back((z+1)/2);
        w.push_back(1/((1-z*z)*dp*dp));
      }
    }
  }
  ;
  const Q q4(4),q8(8),q16(16);
  double T(double r){
    if(thermalReference)return 1+.2*std::log(r)/std::log(2.);
    return C*std::log(r)+D-G/(r*r);
  }
  double v(double r){
    if(thermalReference)return 0;
    return .3*r+.2/r;
  }
  double pressure(double r,const Q&q=q16){
    if(thermalReference)return 1;
    double a=0;
    for(size_t i=0;
    i<q.x.size();
    ++i){
      double s=1+(r-1)*q.x[i];
      a+=q.w[i]*v(s)*v(s)/(T(s)*s);
    }
    return std::exp((r-1)*a);
  }
  EulerConservative2D exact(Point2D p){
    double r=std::hypot(p.x,p.y),pr=pressure(r),rho=pr/T(r);
    return eulerConservative2D({
      rho,-v(r)*p.y/r,v(r)*p.x/r,pr
    }
    ,gas);
  }
  EulerConservative2D triangle(Point2D a,Point2D b,Point2D c,const Q&q){
    EulerConservative2D out{
    }
    ;
    double cross=(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
    for(size_t i=0;
    i<q.x.size();
    ++i)for(size_t j=0;
    j<q.x.size();
    ++j){
      double u=q.x[i],z=q.x[j];
      auto f=exact({
        a.x+u*(b.x-a.x)+u*z*(c.x-b.x),a.y+u*(b.y-a.y)+u*z*(c.y-b.y)
      }
      );
      double w=q.w[i]*q.w[j]*u*cross;
      for(int k=0;
      k<4;
      ++k)out[k]+=w*f[k];
    }
    return out;
  }
}
int main(int argc,char**argv){
  try{
    if(argc==2&&std::string(argv[1])=="--reference-test"){
      if(std::abs(T(1)-1)>1e-14||std::abs(T(2)-1.2)>1e-14)throw std::runtime_error("temperature endpoints");
      for(double r:{
        1.,1.2,1.6,2.
      }
      ){
        double h=1e-4;
        auto logp=[](double x){
          return std::log(pressure(x));
        }
        ;
        double derivative=(-logp(r+2*h)+8*logp(r+h)-8*logp(r-h)+logp(r-2*h))/(12*h);
        if(std::abs(derivative-v(r)*v(r)/(T(r)*r))>1e-9)throw std::runtime_error("radial pressure balance");
      }
      double heat=-4*pi*tr.thermalConductivity*G*(.25-1),work=4*pi*.02*.2*.2*(.25-1);
      if(std::abs(heat+work)>1e-15)throw std::runtime_error("thermal mechanical balance");
      for(auto*q:{
        &q4,&q8,&q16
      }
      ){
        double sum=0;
        for(double w:q->w){
          if(w<=0)throw std::runtime_error("quadrature weight");
          sum+=w;
        }
        if(std::abs(sum-1)>1e-14)throw std::runtime_error("quadrature normalization");
      }
      std::cout<<"native reference identities passed\n";
      return 0;
    }
    if(argc<6 || argc>8)throw std::runtime_error("mesh prefix exact|rest end dt [budget] [checkpoint]");
    std::signal(SIGINT,requestStop);
    std::signal(SIGTERM,requestStop);
    auto begin=std::chrono::steady_clock::now();
    auto secs=[&]{
      return std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    }
    ;
    auto rd=readCm2dTopology(argv[1]);
    if(!rd.valid())throw std::runtime_error(rd.error);
    auto m=makeFvMesh2D(rd.topology);
    std::string prefix=argv[2],init=argv[3];
    const bool readOnly=init=="snapshot"||init=="thermal-snapshot";
    thermalReference=init=="thermal-exact"||init=="thermal-rest"||init=="thermal-snapshot";
    if(readOnly)init=thermalReference?"thermal-exact":"exact";
    if(thermalReference)init=init=="thermal-rest"?"rest":"exact";
    if(init!="exact"&&init!="rest")throw std::runtime_error("initialization");
    double end=std::stod(argv[4]),dt=std::stod(argv[5]),budget=argc>6?std::stod(argv[6]):180;
    if(!(end>0&&dt>0&&budget>0))throw std::runtime_error("positive horizon/step/budget required");
    std::vector<EulerConservative2D> avg(m.cells.size()),point(m.cells.size());
    double mass=0,quadError=0,pressureQuad=0,normalDefect=0,wallApprox=0;
    for(size_t c=0;
    c<m.cells.size();
    ++c){
      EulerConservative2D low{
      }
      ;
      auto centre=m.cells[c].centre;
      auto&poly=rd.topology.cells[c].vertices;
      for(size_t j=0;
      j<poly.size();
      ++j){
        auto a=rd.topology.vertices[poly[j]].point,b=rd.topology.vertices[poly[(j+1)%poly.size()]].point;
        auto h=triangle(centre,a,b,q8),l=triangle(centre,a,b,q4);
        for(int k=0;
        k<4;
        ++k){
          avg[c][k]+=h[k];
          low[k]+=l[k];
        }
      }
      for(int k=0;
      k<4;
      ++k){
        quadError=std::max(quadError,std::abs(avg[c][k]-low[k])/m.cells[c].area);
        avg[c][k]/=m.cells[c].area;
      }
      mass+=avg[c][0]*m.cells[c].area;
      point[c]=exact(centre);
      double r=std::hypot(centre.x,centre.y);
      pressureQuad=std::max(pressureQuad,std::abs(pressure(r)-pressure(r,q8)));
    }
    const double targetMass=3*pi,norm=targetMass/mass;
    double circleMass=0;
    for(size_t i=0;
    i<q16.x.size();
    ++i){
      double r=1+q16.x[i];
      circleMass+=2*pi*r*pressure(r)/T(r)*q16.w[i];
    }
    double circleNorm=targetMass/circleMass;
    for(auto*f:{
      &avg,&point
    }
    )for(auto&a:*f)for(double&x:a)x*=norm;
    std::vector<EulerBoundary2D> bc;
    std::ofstream boundaries(prefix+".wall-definition.csv");
    boundaries<<std::setprecision(17)<<"face,x,y,sx,sy,inner,u,v,temperature,normalDefect,projectionDifference\n";
    for(size_t f=0;
    f<m.faces.size();
    ++f)if(!m.faces[f].neighbour){
      auto&face=m.faces[f];
      double r=std::hypot(face.centre.x,face.centre.y),l=std::hypot(face.areaVector.x,face.areaVector.y);
      bool inner=r<1.5;
      double radius=inner?1:2;
      Vector2D n{
        face.areaVector.x/l,face.areaVector.y/l
      }
      ;
      Vector2D circular{
        -v(radius)*face.centre.y/r,v(radius)*face.centre.x/r
      }
      ;
      double vn=circular.x*n.x+circular.y*n.y;
      Vector2D tang{
        circular.x-vn*n.x,circular.y-vn*n.y
      }
      ;
      double defect=std::abs(tang.x*n.x+tang.y*n.y);
      normalDefect=std::max(normalDefect,defect);
      wallApprox=std::max(wallApprox,std::abs(vn));
      EulerBoundary2D b;
      b.face=f;
      b.kind=EulerBoundaryKind2D::NoSlipWall;
      b.thermalKind=HeatBoundaryKind2D::Temperature;
      b.thermalValue=thermalReference?T(r):(inner?1:1.2);
      b.wallVelocity=tang;
      b.name=inner?"inner":"outer";
      bc.push_back(b);
      boundaries<<f<<','<<face.centre.x<<','<<face.centre.y<<','<<face.areaVector.x<<','<<face.areaVector.y<<','<<inner<<','<<tang.x<<','<<tang.y<<','<<b.thermalValue<<','<<defect<<','<<std::abs(vn)<<'\n';
    }
    EulerState2D state;
    state.cells=avg;
    if(init=="rest")for(auto&a:state.cells){
      a[1]=a[2]=0;
      a[3]=a[0]*1.1/(gas.gamma-1);
    }
    if(argc>7){
      std::ifstream f(argv[7]);
      state=readEulerCheckpoint2D(f,m,bc,gas,"couette-reference",tr);
    }
    auto save=[&]{
      const auto temporary=prefix+".checkpoint.tmp";
      std::ofstream f(temporary);
      writeEulerCheckpoint2D(f,m,bc,gas,state,"couette-reference",tr);
      f.close();
      if(!f)throw std::runtime_error("checkpoint write failed");
      std::filesystem::rename(temporary,prefix+".checkpoint");
    }
    ;
    save();
    EulerStepper2D solver(m,bc,gas,tr);
    EulerStepControls2D ctl;
    ctl.integrator=EulerTimeIntegrator2D::Sdirk2;
    ctl.fluxScheme=EulerFluxScheme2D::Hllc;
    ctl.order=2;
    ctl.maximumStep=dt;
    ctl.endTime=end;
    ctl.interrupted=[&]{
      return stopRequested||secs()>budget;
    }
    ;
    std::ofstream referenceFaces(prefix+".reference-faces.csv");
    referenceFaces<<std::setprecision(17)<<"face,advX,advY,pressureX,pressureY,viscX,viscY,mass\n";
    for(std::size_t id=0;id<m.faces.size();++id){
      const auto& f=m.faces[id];const double len=std::hypot(f.areaVector.x,f.areaVector.y);
      const Vector2D tangent{-f.areaVector.y/len,f.areaVector.x/len};
      std::array<double,7> integral{};
      for(std::size_t j=0;j<q16.x.size();++j){
        const Point2D p{f.centre.x+(q16.x[j]-.5)*len*tangent.x,f.centre.y+(q16.x[j]-.5)*len*tangent.y};
        const auto primitive=eulerPrimitive2D(exact(p),gas);const double un=primitive.u*f.areaVector.x+primitive.v*f.areaVector.y;
        const double r2=p.x*p.x+p.y*p.y;
        const double ux=thermalReference?0:.4*p.x*p.y/(r2*r2),vy=-ux;
        const double shear=thermalReference?0:.4*(p.y*p.y-p.x*p.x)/(r2*r2);
        const std::array<double,7> value{norm*primitive.density*primitive.u*un,norm*primitive.density*primitive.v*un,norm*primitive.pressure*f.areaVector.x,norm*primitive.pressure*f.areaVector.y,-tr.dynamicViscosity*(2*ux*f.areaVector.x+shear*f.areaVector.y),-tr.dynamicViscosity*(shear*f.areaVector.x+2*vy*f.areaVector.y),norm*primitive.density*un};
        for(std::size_t k=0;k<7;++k)integral[k]+=q16.w[j]*value[k];
      }
      referenceFaces<<id;for(double value:integral)referenceFaces<<','<<value;referenceFaces<<'\n';
    }
    auto snapshot=[&](const std::string& label,const std::vector<EulerConservative2D>& values){
      EulerState2D probe;probe.cells=values;
      const auto op=solver.spatialSnapshot(probe,ctl);
      std::ofstream out(prefix+".snapshot-"+label+".csv");
      out<<std::setprecision(17)<<"face,owner,neighbour,x,y,sx,sy,mass,mx,my,energy,viscX,viscY,work,heat,pressureX,pressureY\n";
      for(std::size_t f=0;f<m.faces.size();++f){
        const auto& face=m.faces[f];out<<f<<','<<face.owner<<','<<(face.neighbour?std::to_string(*face.neighbour):"-1")<<','<<face.centre.x<<','<<face.centre.y<<','<<face.areaVector.x<<','<<face.areaVector.y;
        for(auto value:op.faceFlux[f])out<<','<<value;
        for(auto value:op.faceViscousFlux[f])out<<','<<value;
        out<<','<<op.faceHeatFlux[f]<<','<<op.facePressureFlux[f].x<<','<<op.facePressureFlux[f].y<<'\n';
      }
      std::ofstream residual(prefix+".residual-"+label+".csv");
      residual<<std::setprecision(17)<<"cell,area,mass,mx,my,energy\n";
      for(std::size_t c=0;c<m.cells.size();++c){residual<<c<<','<<m.cells[c].area;for(auto value:op.cellResidual[c])residual<<','<<value;residual<<'\n';}
    };
    snapshot("point",point);snapshot("mean",avg);snapshot("input",state.cells);
    std::ofstream hist(prefix+".history.csv");
    hist<<std::setprecision(17)<<"step,time,dt,residualRate,changeRate,wallOutputRate,innerHeat,outerHeat,innerWork,outerWork,innerTorque,outerTorque,balanceMass,balanceEnergy,maxCellBalance,rejected,spatial,newton,krylov\n";
    std::string status=readOnly?"read-only-space-snapshot":"target-time-not-steady",error;
    EulerStepResult2D last;
    double rate=solver.residualRate(state,ctl),initialRate=rate,change=1e300,outputRate=1e300;
    std::array<double,6> prev{
    }
    ;
    bool have=false;
    size_t stable=0,calls=0,reject=0,newton=0,krylov=0;
    const double tref=1/std::sqrt(1.4),gate=1e-5/tref;
    try{
      while(!readOnly&&state.time<end){
        auto old=state;
        auto r=solver.advance(state,ctl);
        calls+=r.spatialEvaluations;
        reject+=r.rejectedCandidates;
        newton+=r.nonlinearIterations;
        krylov+=r.linearIterations;
        state=r.state;
        last=r;
        change=0;
        for(size_t c=0;
        c<state.cells.size();
        ++c)for(int k=0;
        k<4;
        ++k)change=std::max(change,std::abs(state.cells[c][k]-old.cells[c][k])/(r.step*std::max(1.,std::abs(avg[c][k]))));
        rate=solver.residualRate(state,ctl);
        std::array<double,6> out{
        }
        ;
        for(auto&b:bc){
          size_t side=b.name=="inner"?0:1;
          auto f=b.face;
          out[side]+=r.faceHeatFlux[f];
          out[2+side]+=r.faceViscousFlux[f][2];
          out[4+side]+=m.faces[f].centre.x*r.faceViscousFlux[f][1]-m.faces[f].centre.y*r.faceViscousFlux[f][0];
        }
        outputRate=0;
        if(!have)outputRate=1e300;
        else for(size_t k=0;
        k<6;
        ++k)outputRate=std::max(outputRate,std::abs(out[k]-prev[k])/(r.step*.1));
        prev=out;
        have=true;
        hist<<state.steps<<','<<state.time<<','<<r.step<<','<<rate<<','<<change<<','<<outputRate;
        for(double x:out)hist<<','<<x;
        hist<<','<<r.balanceError[0]<<','<<r.balanceError[3]<<','<<r.maximumCellBalanceError<<','<<r.rejectedCandidates<<','<<r.spatialEvaluations<<','<<r.nonlinearIterations<<','<<r.linearIterations<<'\n';
        hist.flush();
        save();
        stable=rate<gate&&change<gate&&outputRate<gate?stable+1:0;
        if(stable>=10){
          status="discrete-steady-three-gates";
          break;
        }
      }
    }
    catch(const std::exception&e){
      status="failed-last-accepted";
      error=e.what();
      save();
    }
    snapshot("final",state.cells);
    std::ofstream cells(prefix+".cells.csv");
    cells<<std::setprecision(17)<<"cell,x,y,area,rho,u,v,p,T,pointRho,pointU,pointV,pointP,pointT,meanRho,meanU,meanV,meanP,meanT\n";
    double epoint[5]{
    }
    ,emean[5]{
    }
    ,vol=0,minT=1e300,minRho=1e300,minP=1e300;
    for(size_t c=0;
    c<m.cells.size();
    ++c){
      auto a=eulerPrimitive2D(state.cells[c],gas),b=eulerPrimitive2D(point[c],gas),h=eulerPrimitive2D(avg[c],gas);
      double x[]={
        a.density,a.u,a.v,a.pressure,a.pressure/a.density
      }
      ,y[]={
        b.density,b.u,b.v,b.pressure,b.pressure/b.density
      }
      ,z[]={
        h.density,h.u,h.v,h.pressure,h.pressure/h.density
      }
      ;
      cells<<c<<','<<m.cells[c].centre.x<<','<<m.cells[c].centre.y<<','<<m.cells[c].area;
      for(auto*f:{
        x,y,z
      }
      )for(int k=0;
      k<5;
      ++k)cells<<','<<f[k];
      cells<<'\n';
      for(int k=0;
      k<5;
      ++k){
        epoint[k]+=m.cells[c].area*std::abs(x[k]-y[k]);
        emean[k]+=m.cells[c].area*std::abs(x[k]-z[k]);
      }
      vol+=m.cells[c].area;
      minT=std::min(minT,x[4]);
      minRho=std::min(minRho,x[0]);
      minP=std::min(minP,x[3]);
    }
    std::ofstream ref(prefix+".analytic-wall.csv");
    ref<<std::setprecision(17)<<"wall,heat,work,torque\n";
    for(int side=0;side<2;++side){
      double heat=0;
      if(thermalReference){
        for(const auto& boundary:bc)if((boundary.name=="inner")==!side){
          const auto& f=m.faces[boundary.face];const double len=std::hypot(f.areaVector.x,f.areaVector.y);
          const Vector2D tangent{-f.areaVector.y/len,f.areaVector.x/len};
          for(std::size_t j=0;j<q16.x.size();++j){
            const double x=f.centre.x+(q16.x[j]-.5)*len*tangent.x,y=f.centre.y+(q16.x[j]-.5)*len*tangent.y;
            heat-=q16.w[j]*tr.thermalConductivity*.2/std::log(2.)*(x*f.areaVector.x+y*f.areaVector.y)/(x*x+y*y);
          }
        }
        ref<<(side?"outer":"inner")<<','<<heat<<",0,0\n";
      }else{
        const double r=side?2:1,n=side?1:-1;
        ref<<(side?"outer":"inner")<<','<<-2*pi*tr.thermalConductivity*n*(C+2*G/(r*r))<<','<<4*pi*tr.dynamicViscosity*.2*n*(.3+.2/(r*r))<<','<<4*pi*tr.dynamicViscosity*.2*n<<'\n';
      }
    }
    std::ofstream faces(prefix+".faces.csv");
    faces<<std::setprecision(17)<<"face,x,y,inner,mass,heat,work,mx,my,viscX,viscY\n";
    if(!last.faceFlux.empty())for(auto&b:bc){
      size_t f=b.face;
      faces<<f<<','<<m.faces[f].centre.x<<','<<m.faces[f].centre.y<<','<<(b.name=="inner")<<','<<last.faceFlux[f][0]<<','<<last.faceHeatFlux[f]<<','<<last.faceViscousFlux[f][2]<<','<<last.faceFlux[f][1]<<','<<last.faceFlux[f][2]<<','<<last.faceViscousFlux[f][0]<<','<<last.faceViscousFlux[f][1]<<'\n';
    }
    double actualMass=0,maxMach=0;
    for(size_t c=0;
    c<state.cells.size();
    ++c)actualMass+=state.cells[c][0]*m.cells[c].area;
    for(int i=0;
    i<=1000;
    ++i){
      double rr=1+i/1000.;
      maxMach=std::max(maxMach,v(rr)/std::sqrt(1.4*T(rr)));
    }
    double boundaryTorque=0,pressureTorque=0,internalConvectiveTorque=0,internalViscousTorque=0,boundaryLeverCorrection=0,cellAngularResidual=0;
    std::ofstream allfaces(prefix+".all-faces.csv");
    allfaces<<std::setprecision(17)<<"face,owner,neighbour,x,y,mx,my,viscX,viscY\n";
    if(!last.faceFlux.empty())for(size_t f=0;
    f<m.faces.size();
    ++f){
      auto&face=m.faces[f];
      auto a=m.cells[face.owner].centre;
      double mx=last.faceFlux[f][1],my=last.faceFlux[f][2],vx=last.faceViscousFlux[f][0],vy=last.faceViscousFlux[f][1];
      if(face.neighbour){
        auto b=m.cells[*face.neighbour].centre;
        internalConvectiveTorque+=(a.x-b.x)*(my-vy)-(a.y-b.y)*(mx-vx);
        internalViscousTorque+=(a.x-b.x)*vy-(a.y-b.y)*vx;
        cellAngularResidual+=(a.x-b.x)*my-(a.y-b.y)*mx;
      }
      else{
        boundaryTorque+=face.centre.x*my-face.centre.y*mx;
        pressureTorque+=face.centre.x*(my-vy)-face.centre.y*(mx-vx);
        boundaryLeverCorrection+=(face.centre.x-a.x)*my-(face.centre.y-a.y)*mx;
        cellAngularResidual+=a.x*my-a.y*mx;
      }
      allfaces<<f<<','<<face.owner<<','<<(face.neighbour?std::to_string(*face.neighbour):"-1")<<','<<face.centre.x<<','<<face.centre.y<<','<<mx<<','<<my<<','<<vx<<','<<vy<<'\n';
    }
    std::ofstream report(prefix+".json");
    report<<std::setprecision(17)<<"{\"referenceKind\":"<<std::quoted(thermalReference?"polygon-static-harmonic-heat":"circle-rotating-Couette")<<",\"pointPressureMassMapping\":\"polygon-normalized\",\"status\":"<<std::quoted(status)<<",\"error\":"<<std::quoted(error)<<",\"cells\":"<<m.cells.size()<<",\"time\":"<<state.time<<",\"steps\":"<<state.steps<<",\"residualRate\":"<<rate<<",\"changeRate\":"<<change<<",\"wallOutputRate\":"<<outputRate<<",\"steadyGateRatePerSecond\":"<<gate<<",\"quadrature8vs4ConservedMax\":"<<quadError*norm<<",\"pressureQuadrature16vs8\":"<<pressureQuad*norm<<",\"boundaryTotalTorque\":"<<boundaryTorque<<",\"boundaryPressureTorque\":"<<pressureTorque<<",\"internalCentroidConvectiveAngularFlux\":"<<internalConvectiveTorque<<",\"internalCentroidViscousAngularFlux\":"<<internalViscousTorque<<",\"boundaryCentroidLeverCorrection\":"<<boundaryLeverCorrection<<",\"cellCentroidAngularResidual\":"<<cellAngularResidual<<",\"initialResidualRate\":"<<initialRate<<",\"maximumReferenceMach\":"<<maxMach<<",\"actualMass\":"<<actualMass<<",\"targetMass\":"<<targetMass<<",\"polygonPressureScale\":"<<norm<<",\"circlePressureScale\":"<<circleNorm<<",\"maximumWallNormalDefect\":"<<normalDefect<<",\"maximumCircleVelocityRemovedNormal\":"<<wallApprox<<",\"minimumTemperature\":"<<minT<<",\"minimumDensity\":"<<minRho<<",\"minimumPressure\":"<<minP<<",\"spatialCalls\":"<<calls<<",\"rejected\":"<<reject<<",\"newton\":"<<newton<<",\"krylov\":"<<krylov<<",\"pointL1\":[";
    for(int k=0;
    k<5;
    ++k)report<<(k?",":"")<<epoint[k]/vol;
    report<<"],\"conservativeMeanPrimitiveL1\":[";
    for(int k=0;
    k<5;
    ++k)report<<(k?",":"")<<emean[k]/vol;
    report<<"],\"fullSecondsIncludingExport\":"<<secs()<<"}\n";
    std::cout<<status<<" time="<<state.time<<" residual="<<rate<<" seconds="<<secs()<<'\n';
    return status=="failed-last-accepted"?2:0;
  }
  catch(const std::exception&e){
    std::cerr<<e.what()<<'\n';
    return 1;
  }
}
