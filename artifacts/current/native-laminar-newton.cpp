// Research Newton solve of the SAME compatible conservative RT1 residual.
// Differentiate the advector as well as the transported velocity. On a face
// the sign branch is frozen at the current flux, and integration is split at
// its actual zero. The moving split has no boundary term because flux=0 there.
#define CARTMESH_OPEN_BOUNDARY_NO_MAIN
#include "native-laminar-open-boundary.cpp"

template<class Base> struct NewtonTransport : Base {
    static constexpr bool linearizedNewton=true;
    Mat advectorDerivative,picardCondensed;
    NewtonTransport(const Fixture& f,int t,const std::string& problem,double nu,bool nonlinear,int model,int order,const Vec& previous):
        Base(f,t,problem,nu,nonlinear,model,order,previous),advectorDerivative(2*this->e.a.m,2*this->e.a.m),picardCondensed(this->e.condensed().first){
        if(!nonlinear)return;
        auto& e=this->e;const int m=e.a.m;
        auto& lift=this->lift;auto& beta=this->beta;
        for(std::size_t k=0;k<lift.tri.size();++k){const auto& tr=lift.tri[k];
            for(auto q:tr.quadrature(order)){
                const auto z=tr.ref(q.p);const auto u=this->betaAt(int(k),q.p);
                std::array<Vector2D,8> dx{{{0,0},{1,0},{0,0},{0,0},{0,1},{0,0},{2*z.x,z.y},{z.y,0}}},dy{{{0,0},{0,0},{1,0},{0,0},{0,0},{0,1},{0,z.x},{z.x,2*z.y}}};
                std::vector<Vector2D> gx(2*m),gy(2*m),value(2*m);
                for(int j=0;j<2*m;++j){value[j]=lift.value(int(k),j,q.p);
                    for(int l=0;l<8;++l){
                        const auto physical=[&](double x,double y){Vector2D d{x*dx[l].x+y*dy[l].x,x*dx[l].y+y*dy[l].y};return Vector2D{(tr.j0.x*d.x+tr.j1.x*d.y)/tr.det,(tr.j0.y*d.x+tr.j1.y*d.y)/tr.det};};
                        const auto x=physical(tr.j1.y/(tr.h*tr.det),-tr.j0.y/(tr.h*tr.det));
                        const auto y=physical(-tr.j1.x/(tr.h*tr.det),tr.j0.x/(tr.h*tr.det));
                        const double c=lift.coefficients(8*k+l,j);gx[j].x+=c*x.x;gx[j].y+=c*x.y;gy[j].x+=c*y.x;gy[j].y+=c*y.y;
                    }}
                for(int i=0;i<2*m;++i)for(int j=0;j<2*m;++j)advectorDerivative(i,j)-=q.w*(value[j].x*dot(gx[i],u)+value[j].y*dot(gy[i],u));
            }
            const int prev=(int(k)+int(lift.tri.size())-1)%int(lift.tri.size());const auto d=tr.a-tr.c;Vector2D S{d.y,-d.x};
            for(auto [z,w]:upwindRule(order,dot(this->betaAt(int(k),tr.c),S),dot(this->betaAt(int(k),tr.a),S))){
                Point2D p{tr.c.x+z*d.x,tr.c.y+z*d.y};const double flux=dot(this->betaAt(int(k),p),S);
                const auto u=flux>=0?this->betaAt(int(k),p):this->betaAt(prev,p);std::vector<Vector2D> left(2*m),jump(2*m);
                for(int j=0;j<2*m;++j){left[j]=lift.value(int(k),j,p);jump[j]=difference(left[j],lift.value(prev,j,p));}
                for(int i=0;i<2*m;++i)for(int j=0;j<2*m;++j)advectorDerivative(i,j)+=w*dot(left[j],S)*dot(jump[i],u);
            }
            const int l=lift.faceLocal[k];const auto& face=f.mesh.faces[f.mesh.cells[t].faces[l]];const double sign=face.owner==std::size_t(t)?1:-1;const auto sf=face.areaVector;
            auto fluxEnd=[&](double s){return sign*((beta[3+2*l]+s*beta[4+2*l])*sf.x+(beta[m+3+2*l]+s*beta[m+4+2*l])*sf.y);};
            for(auto [z,w]:upwindRule(order,fluxEnd(-.5),fluxEnd(.5))){const double s=z-.5;Point2D p{face.centre.x-s*sf.y,face.centre.y+s*sf.x};const double flux=fluxEnd(s);
                Vector2D uf{beta[3+2*l]+s*beta[4+2*l],beta[m+3+2*l]+s*beta[m+4+2*l]};const auto upwind=flux>=0?this->betaAt(int(k),p):uf;std::vector<Vector2D> facev(2*m);
                facev[3+2*l].x=1;facev[4+2*l].x=s;facev[m+3+2*l].y=1;facev[m+4+2*l].y=s;
                for(int i=0;i<2*m;++i){const auto jump=difference(lift.value(int(k),i,p),facev[i]);
                    for(int j=0;j<2*m;++j){const double df=sign*dot(facev[j],sf);advectorDerivative(i,j)+=w*df*dot(jump,upwind);
                        if(boundaryKind(face,model,problem)==1)advectorDerivative(i,j)+=w*df*dot(facev[i],uf);}}
            }
        }
        // Solve J(u_k) u_candidate = f + D_beta C(u_k)[u_k] u_k.
        // This is exactly the increment equation J d = -F, with inhomogeneous
        // velocity traces kept in the ordinary retained-unknown elimination.
        for(int i=0;i<2*m;++i)for(int j=0;j<2*m;++j){e.matrix(i,j)+=advectorDerivative(i,j);e.rhs[i]+=advectorDerivative(i,j)*beta[j];}
        e.condense();
    }
};

template<class Boundary> std::vector<bool> freeStateDofs(const Fixture& f,int model,const std::string& problem){
    const int nc=int(f.mesh.cells.size());std::vector<bool> free(9*nc+4*f.mesh.faces.size(),true);
    if(model==0)free[9*(nc-1)+6]=false;
    for(std::size_t id=0;id<f.mesh.faces.size();++id){const auto& face=f.mesh.faces[id];const int kind=Boundary::kind(face,model,problem);
        for(int c=0;c<2;++c)for(int k=0;k<2;++k)free[9*nc+4*id+2*c+k]=face.neighbour||kind==1||(kind==2&&c==0);}
    return free;
}
template<class Boundary> void imposeTrace(const Fixture& f,int model,const std::string& problem,int order,Vec& state){
    const int nc=int(f.mesh.cells.size());const auto free=freeStateDofs<Boundary>(f,model,problem);
    for(std::size_t id=0;id<f.mesh.faces.size();++id){const auto& face=f.mesh.faces[id];if(face.neighbour)continue;std::array<double,4> known{};
        if(problem=="cylinder"){if(face.patch==BoundaryPatch2D::DomainBoundary&&Boundary::kind(face,model,problem)==0)known[0]=1;}
        else if(problem!="noslip"&&problem!="noslip-sheared")for(auto [z,w]:gauss(order)){const double s=z-.5;const auto sf=face.areaVector;const auto u=exactAt({face.centre.x-s*sf.y,face.centre.y+s*sf.x},problem).u;known[0]+=w*u.x;known[1]+=12*w*s*u.x;known[2]+=w*u.y;known[3]+=12*w*s*u.y;}
        for(int j=0;j<4;++j)if(!free[9*nc+4*id+j])state[9*nc+4*id+j]=known[j];}
    if(model==0){const double p=state[9*(nc-1)+6];for(int t=0;t<nc;++t)state[9*t+6]-=p;}
}
int rawIndex(const Fixture& f,int t,int m,int local){
    if(local>=2*m)return 9*t+6+local-2*m;
    const int c=local/m,k=local%m;if(k<3)return 9*t+3*c+k;
    return 9*int(f.mesh.cells.size())+4*int(f.mesh.cells[t].faces[(k-3)/2])+2*c+(k-3)%2;
}
template<class Base,class Boundary> Vec residualOrDerivative(const Fixture& f,const std::string& problem,double nu,int model,int order,const Vec& state,const Vec* direction){
    Vec result(state.size());
    for(std::size_t t=0;t<f.mesh.cells.size();++t){
        auto accumulate=[&](const auto& o,bool derivative){const auto& e=o.e;const int m=e.a.m;const auto local=localState(f,int(t),m,derivative?*direction:state);
            for(int i=0;i<2*m+3;++i){double r=derivative?0:-e.rhs[i];for(int j=0;j<2*m+3;++j)r+=e.matrix(i,j)*local[j];result[rawIndex(f,int(t),m,i)]+=r;}};
        if(direction){NewtonTransport<Base> o(f,int(t),problem,nu,true,model,order,state);accumulate(o,true);}
        else{Base o(f,int(t),problem,nu,true,model,order,state);accumulate(o,false);}
    }
    const auto free=freeStateDofs<Boundary>(f,model,problem);for(std::size_t i=0;i<result.size();++i)if(!free[i])result[i]=0;return result;
}
template<class Base,class Boundary> int dispatchNewton(int argc,char** argv){
    const std::string mode=argv[1];
    if(mode=="seed"){
        if(argc!=9)throw std::runtime_error("usage: newton seed mesh n problem boundary order initial.state|zero output.state");
        auto f=readFixture(argv[2],std::stoi(argv[3]));const std::string problem=argv[4];const int model=Boundary::model(argv[5]),order=std::stoi(argv[6]);
        if(!Boundary::allows(problem,model)||order<4||order>12||std::filesystem::exists(argv[8]))throw std::runtime_error("invalid seed options or existing output");
        auto state=readState(argv[7],f);imposeTrace<Boundary>(f,model,problem,order,state);writeState(argv[8],f,state);
        std::cout<<"{\"seedOnly\":true,\"physicalSolution\":false,\"stateSize\":"<<state.size()<<"}\n";return 0;
    }
    if(mode=="tangent"){
        if(argc!=9)throw std::runtime_error("usage: newton tangent mesh n problem nu boundary order epsilon");
        auto f=readFixture(argv[2],std::stoi(argv[3]));const std::string problem=argv[4];const double nu=std::stod(argv[5]),eps=std::stod(argv[8]);const int model=Boundary::model(argv[6]),order=std::stoi(argv[7]);
        if(!std::isfinite(nu)||!(nu>0)||!std::isfinite(eps)||!(eps>0)||order<4||order>12||!Boundary::allows(problem,model))throw std::runtime_error("invalid tangent options");
        Vec state(9*f.mesh.cells.size()+4*f.mesh.faces.size()),d(state.size());const auto free=freeStateDofs<Boundary>(f,model,problem);
        for(std::size_t i=0;i<state.size();++i){state[i]=.3*std::sin(.71*(i+1));d[i]=free[i]?std::cos(.37*(i+1)):0;}imposeTrace<Boundary>(f,model,problem,order,state);
        auto plus=state,minus=state;for(std::size_t i=0;i<state.size();++i){plus[i]+=eps*d[i];minus[i]-=eps*d[i];}
        const auto rp=residualOrDerivative<Base,Boundary>(f,problem,nu,model,order,plus,nullptr),rm=residualOrDerivative<Base,Boundary>(f,problem,nu,model,order,minus,nullptr),jd=residualOrDerivative<Base,Boundary>(f,problem,nu,model,order,state,&d);
        double err=0,norm=0,maximum=0;for(std::size_t i=0;i<state.size();++i){const double difference=(rp[i]-rm[i])/(2*eps)-jd[i];err+=difference*difference;norm+=jd[i]*jd[i];maximum=std::max(maximum,std::abs(difference));}
        std::cout<<std::setprecision(17)<<"{\"epsilon\":"<<eps<<",\"relativeDerivativeError\":"<<std::sqrt(err/norm)<<",\"absoluteDerivativeMax\":"<<maximum<<",\"directionDerivativeL2\":"<<std::sqrt(norm)<<"}\n";return 0;
    }
    // Never interpret a Newton recovery's linearized reaction as physical.
    // Full residual, energy, force and error checks use the original operator.
    if(mode=="check")return runOseen<Base,Boundary>(argc,argv);
    return runOseen<NewtonTransport<Base>,Boundary>(argc,argv);
}
int runNewton(int argc,char** argv)try{
    if(argc<2)throw std::runtime_error("newton requires a mode");
    const std::string mode=argv[1];
    if(mode=="blend"){
        if(argc!=8)throw std::runtime_error("usage: newton blend mesh n current.state candidate.state alpha output.state");
        auto f=readFixture(argv[2],std::stoi(argv[3]));auto v=readState(argv[4],f),candidate=readState(argv[5],f);const double alpha=std::stod(argv[6]);
        if(!std::isfinite(alpha)||!(alpha>0)||alpha>1||std::filesystem::exists(argv[7]))throw std::runtime_error("invalid blend or existing output");
        double maximum=0;for(std::size_t i=0;i<v.size();++i){const double change=alpha*(candidate[i]-v[i]);v[i]+=change;if(!std::isfinite(v[i]))throw std::runtime_error("nonfinite trial");maximum=std::max(maximum,std::abs(change));}
        writeState(argv[7],f,v);std::cout<<std::setprecision(17)<<"{\"stateCoefficientChange\":"<<maximum<<",\"alpha\":"<<alpha<<"}\n";return 0;
    }
    const int at=mode=="seed"?5:mode=="tangent"?6:7;if(argc<=at)throw std::runtime_error("missing boundary mode");const std::string boundary=argv[at];
    if(boundary=="traction")return dispatchNewton<OpenTransport<OutletForm::ExactTraction>,ExplicitOpenBoundary<OutletForm::ExactTraction>>(argc,argv);
    if(boundary=="pseudo-traction")return dispatchNewton<OpenTransport<OutletForm::PseudoTraction>,ExplicitOpenBoundary<OutletForm::PseudoTraction>>(argc,argv);
    if(boundary=="normal-stress")return dispatchNewton<OpenTransport<OutletForm::NormalStress>,ExplicitOpenBoundary<OutletForm::NormalStress>>(argc,argv);
    return dispatchNewton<Transport,StandardOseenBoundary>(argc,argv);
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}

#ifndef CARTMESH_NEWTON_NO_MAIN
int main(int argc,char** argv){return runNewton(argc,argv);}
#endif
