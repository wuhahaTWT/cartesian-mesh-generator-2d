// Pseudo-transient globalization of the unchanged steady native residual.
// Only cell P1 velocity receives an L2 mass. Face velocity and pressure are
// algebraic DAE variables; no artificial pressure time derivative is added.
// This is a steady solver path, not a qualified physical transient method.
#define CARTMESH_NEWTON_NO_MAIN
#include "native-laminar-newton.cpp"

double inversePseudoTime=0;
template<class Base> struct PseudoMass : Base {
    static constexpr bool iterationUsesState=true;
    PseudoMass(const Fixture& f,int t,const std::string& problem,double nu,bool nonlinear,int model,int order,const Vec& previous):
        Base(f,t,problem,nu,nonlinear,model,order,previous){
        if(inversePseudoTime==0)return;
        auto& e=this->e;const int m=e.a.m;const auto old=localState(f,t,m,previous);
        // [J + M/dt] delta = -F(u), written in absolute-state coordinates.
        // Use the native exact P1 mass in both velocity components. This
        // changes the iteration, while a fixed point still satisfies F=0.
        for(int c=0;c<2;++c)for(int i=0;i<3;++i)for(int j=0;j<3;++j){
            const double mass=inversePseudoTime*e.a.mass(i,j);
            e.matrix(c*m+i,c*m+j)+=mass;e.rhs[c*m+i]+=mass*old[c*m+j];
        }
        e.condense();
    }
};
int main(int argc,char** argv)try{
    if(argc<2)throw std::runtime_error("pseudo-time requires a mode");
    const std::string mode=argv[1];
    // All nonlinear residuals and physical diagnostics use the original
    // operator, never the regularized matrix used to construct a candidate.
    if(mode!="assemble"&&mode!="recover")return runNewton(argc,argv);
    if(argc!=12)throw std::runtime_error("usage: pseudo-time assemble|recover mesh n problem nu stokes|ns boundary order previous.state|zero prefix inverse_dt");
    std::size_t used=0;const std::string value=argv[11];inversePseudoTime=std::stod(value,&used);
    if(used!=value.size()||!std::isfinite(inversePseudoTime)||inversePseudoTime<0)throw std::runtime_error("finite nonnegative inverse pseudo step required");
    const std::string boundary=argv[7];
    if(boundary=="traction")return runOseen<NewtonTransport<PseudoMass<OpenTransport<OutletForm::ExactTraction>>>,ExplicitOpenBoundary<OutletForm::ExactTraction>>(11,argv);
    if(boundary=="pseudo-traction")return runOseen<NewtonTransport<PseudoMass<OpenTransport<OutletForm::PseudoTraction>>>,ExplicitOpenBoundary<OutletForm::PseudoTraction>>(11,argv);
    if(boundary=="normal-stress")return runOseen<NewtonTransport<PseudoMass<OpenTransport<OutletForm::NormalStress>>>,ExplicitOpenBoundary<OutletForm::NormalStress>>(11,argv);
    return runOseen<NewtonTransport<PseudoMass<Transport>>,StandardOseenBoundary>(11,argv);
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
