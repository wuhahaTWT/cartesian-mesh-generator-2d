// Research: conservative RT1-advected P1 hybrid Oseen/NS, static pressure.
// Native geometry/operators/load/condensation/recovery; no product changes.
#define CARTMESH_P1_STRESS_NO_MAIN
#include "native-laminar-p1-stress.cpp"
#include <cstring>

std::uint64_t meshKey(const Fixture& f){
    std::uint64_t h=1469598103934665603ULL;
    auto add=[&](auto x){const auto* p=reinterpret_cast<const unsigned char*>(&x);for(std::size_t i=0;i<sizeof x;++i){h^=p[i];h*=1099511628211ULL;}};
    add(f.mesh.cells.size());add(f.mesh.faces.size());
    for(const auto& c:f.mesh.cells){add(c.centre.x);add(c.centre.y);add(c.area);for(auto i:c.faces)add(i);}
    for(const auto& g:f.mesh.faces){add(g.centre.x);add(g.centre.y);add(g.areaVector.x);add(g.areaVector.y);add(g.owner);add(g.neighbour.value_or(std::size_t(-1)));add(g.patch);}return h;
}
Vec readState(const std::string& path,const Fixture& f){
    Vec v(9*f.mesh.cells.size()+4*f.mesh.faces.size());if(path=="zero")return v;
    std::ifstream in(path,std::ios::binary);std::uint64_t header[3]{};in.read(reinterpret_cast<char*>(header),sizeof header);
    if(header[0]!=0x504f5345454e3031ULL||header[1]!=meshKey(f)||header[2]!=v.size())throw std::runtime_error("state version/geometry/size mismatch");
    in.read(reinterpret_cast<char*>(v.data()),v.size()*sizeof(double));if(!in||in.peek()!=EOF)throw std::runtime_error("truncated/trailing state");
    for(double x:v)if(!std::isfinite(x))throw std::runtime_error("nonfinite state");
    return v;
}
void writeState(const std::string& path,const Fixture& f,const Vec& v){
    std::ofstream out(path+".tmp",std::ios::binary);std::uint64_t header[]{0x504f5345454e3031ULL,meshKey(f),v.size()};out.write(reinterpret_cast<const char*>(header),sizeof header);out.write(reinterpret_cast<const char*>(v.data()),v.size()*sizeof(double));out.close();if(!out)throw std::runtime_error("state write failed");std::filesystem::rename(path+".tmp",path);
}
Vec localState(const Fixture& f,int t,int m,const Vec& state){
    Vec v(2*m+3);for(int c=0;c<2;++c)for(int j=0;j<3;++j)v[c*m+j]=state[9*t+3*c+j];
    for(int p=0;p<3;++p)v[2*m+p]=state[9*t+6+p];
    for(std::size_t l=0;l<f.mesh.cells[t].faces.size();++l)for(int c=0;c<2;++c)for(int j=0;j<2;++j)v[c*m+3+2*l+j]=state[9*f.mesh.cells.size()+4*f.mesh.cells[t].faces[l]+2*c+j];
    return v;
}
// 0 prescribed velocity, 1 natural outlet, 2 symmetry (normal fixed).
int boundaryKind(const Face& face,bool open){
    if(face.neighbour||!open)return 0;
    if(face.patch==BoundaryPatch2D::EmbeddedBoundary)return 0;
    if(face.patch!=BoundaryPatch2D::DomainBoundary)throw std::runtime_error("explicit external patches required");
    auto S=face.areaVector;double round=128*std::numeric_limits<double>::epsilon()*std::hypot(S.x,S.y);
    if(std::abs(S.y)<=round)return S.x>0?1:0;
    if(std::abs(S.x)<=round)return 2;
    throw std::runtime_error("external control requires axis-aligned far boundary");
}
Vector2D forceAt(Point2D p,const std::string& problem,double nu,bool nonlinear){
    if(problem=="cylinder")return {0,0};
    auto ex=exactAt(p,problem);
    if(problem=="couette"||problem=="hydrostatic"||problem=="poiseuille"||problem=="rotation"){
        Vector2D gp=problem=="hydrostatic"?Vector2D{1,2}:problem=="poiseuille"?Vector2D{-8,0}:Vector2D{0,0};
        Vector2D result{nu*(ex.f.x-gp.x)+gp.x,nu*(ex.f.y-gp.y)+gp.y};
        if(nonlinear&&problem=="rotation"){result.x-=p.x;result.y-=p.y;}return result;
    }
    double s=problem=="noslip-sheared"?.7:0,pi=std::acos(-1.),xi=p.x-s*p.y;
    Vector2D gp{pi*std::cos(pi*xi)*std::sin(pi*p.y),pi*(std::sin(pi*xi)*std::cos(pi*p.y)-s*std::cos(pi*xi)*std::sin(pi*p.y))};
    Vector2D force{nu*(ex.f.x-gp.x)+gp.x,nu*(ex.f.y-gp.y)+gp.y};
    if(nonlinear){force.x+=ex.u.x*ex.gradient[0].x+ex.u.y*ex.gradient[0].y;force.y+=ex.u.x*ex.gradient[1].x+ex.u.y*ex.gradient[1].y;}return force;
}
struct Oseen {
    Element e;Lift lift;Mat convection;Vec beta;std::vector<Vec> bc;bool open;double nu;int order;const Fixture& f;int t;
    Oseen(const Fixture& F,int T,const std::string& problem,double viscosity,bool nonlinear,bool Open,int Order,const Vec& previous):
        e(F,T,problem,0,false,true,Order),lift(F,T,e.a,Order),convection(2*e.a.m,2*e.a.m),beta(localState(F,T,e.a.m,previous)),open(Open),nu(viscosity),order(Order),f(F),t(T){
        auto& a=e.a;int m=a.m;const auto& cell=f.mesh.cells[t];e.rhs.assign(e.rhs.size(),0);
        for(int i=0;i<2*m;++i)for(int j=0;j<2*m;++j)e.matrix(i,j)*=nu;
        for(std::size_t k=0;k<lift.tri.size();++k){Vec coeff(8);for(int l=0;l<8;++l)for(int j=0;j<2*m;++j)coeff[l]+=lift.coefficients(8*k+l,j)*beta[j];bc.push_back(coeff);
            for(auto q:lift.tri[k].quadrature(order)){auto phi=a.basis.phi(q.p);auto b=betaAt(int(k),q.p);auto force=forceAt(q.p,problem,nu,nonlinear);
                for(int j=0;j<2*m;++j)e.rhs[j]+=q.w*dot(force,lift.value(int(k),j,q.p));
                for(int c=0;c<2;++c)for(int i=1;i<3;++i)for(int j=0;j<3;++j)convection(c*m+i,c*m+j)-=q.w*phi[j]*(i==1?b.x:b.y)/a.basis.h;
            }}
        for(std::size_t l=0;l<cell.faces.size();++l){const auto& face=f.mesh.faces[cell.faces[l]];double sign=face.owner==std::size_t(t)?1:-1;
            for(auto [z,w]:gauss(order)){double s=z-.5;auto S=face.areaVector;Point2D p{face.centre.x-s*S.y,face.centre.y+s*S.x};auto phi=a.basis.phi(p);double flux=sign*((beta[3+2*l]+s*beta[4+2*l])*S.x+(beta[m+3+2*l]+s*beta[m+4+2*l])*S.y),positive=std::max(flux,0.);
                for(int c=0;c<2;++c){Vec test(m),trial(m);for(int j=0;j<3;++j){test[j]=phi[j];trial[j]=positive*phi[j];}for(int j=0;j<2;++j){test[3+2*l+j]=-(j?s:1);trial[3+2*l+j]=(flux-positive)*(j?s:1);}
                    for(int i=0;i<m;++i)for(int j=0;j<m;++j)convection(c*m+i,c*m+j)+=w*test[i]*trial[j];
                    if(boundaryKind(face,open)==1)for(int i=0;i<2;++i)for(int j=0;j<2;++j)convection(c*m+3+2*l+i,c*m+3+2*l+j)+=w*flux*(i?s:1)*(j?s:1);
                }}}
        for(int i=0;i<2*m;++i)for(int j=0;j<2*m;++j)e.matrix(i,j)+=convection(i,j);
        Mat ii(8,8);for(int i=0;i<8;++i)for(int j=0;j<8;++j)ii(i,j)=e.matrix(e.inside[i],e.inside[j]);DenseLU lu(ii.v,8);
        for(std::size_t j=0;j<e.outside.size();++j){Vec r(8);for(int i=0;i<8;++i)r[i]=e.matrix(e.inside[i],e.outside[j]);r=lu.solve(r);for(int i=0;i<8;++i)e.eliminated(i,int(j))=r[i];}
        Vec r(8);for(int i=0;i<8;++i)r[i]=e.rhs[e.inside[i]];e.loadInternal=lu.solve(r);
    }
    Vector2D betaAt(int k,Point2D p)const{auto shape=lift.tri[k].shape(p);Vector2D b{};for(int l=0;l<8;++l){b.x+=shape[l].x*bc[k][l];b.y+=shape[l].y*bc[k][l];}return b;}
    Vector2D liftedVelocity(const Vec& v,int k,Point2D p)const{Vector2D u{};for(int j=0;j<2*e.a.m;++j){auto r=lift.value(k,j,p);u.x+=v[j]*r.x;u.y+=v[j]*r.y;}return u;}
    double energyIdentity(const Vec& v)const{
        int m=e.a.m;double result=0;const auto& cell=f.mesh.cells[t];
        for(std::size_t k=0;k<lift.tri.size();++k)for(auto q:lift.tri[k].quadrature(order)){auto phi=e.a.basis.phi(q.p);double u=0,w=0,d=0;for(int j=0;j<3;++j){u+=phi[j]*v[j];w+=phi[j]*v[m+j];for(int l=0;l<m;++l)d+=phi[j]*(e.a.gx(j,l)*beta[l]+e.a.gy(j,l)*beta[m+l]);}result+=.5*q.w*d*(u*u+w*w);}
        for(std::size_t l=0;l<cell.faces.size();++l){const auto& face=f.mesh.faces[cell.faces[l]];double sign=face.owner==std::size_t(t)?1:-1;
            for(auto [z,w]:gauss(order)){double s=z-.5;auto S=face.areaVector;auto phi=e.a.basis.phi({face.centre.x-s*S.y,face.centre.y+s*S.x});double flux=sign*((beta[3+2*l]+s*beta[4+2*l])*S.x+(beta[m+3+2*l]+s*beta[m+4+2*l])*S.y);
                for(int c=0;c<2;++c){double u=0,uf=v[c*m+3+2*l]+s*v[c*m+4+2*l];for(int j=0;j<3;++j)u+=phi[j]*v[c*m+j];result+=.5*w*(std::abs(flux)*(u-uf)*(u-uf)+(boundaryKind(face,open)==1?1:-1)*flux*uf*uf);}}}return result;
    }
};
template<class Operator> int runOseen(int argc,char** argv)try{
    if(argc!=11)throw std::runtime_error("usage: p1-oseen assemble|recover|check mesh n noslip|noslip-sheared|cylinder nu stokes|ns closed|open order previous.state|zero prefix");
    std::string mode=argv[1],name=argv[2],problem=argv[4],prefix=argv[10];int n=std::stoi(argv[3]),order=std::stoi(argv[8]);double nu=std::stod(argv[5]);bool nonlinear=std::string(argv[6])=="ns",open=std::string(argv[7])=="open";
    if((mode!="assemble"&&mode!="recover"&&mode!="check")||!(nu>0)||order<4||order>12||(std::string(argv[6])!="ns"&&std::string(argv[6])!="stokes")||(std::string(argv[7])!="open"&&std::string(argv[7])!="closed"))throw std::runtime_error("invalid research options");
    if((problem=="noslip"&&name!="square")||(problem=="noslip-sheared"&&name!="sheared")||(problem!="noslip"&&problem!="noslip-sheared"&&problem!="cylinder"&&problem!="couette"&&problem!="hydrostatic"&&problem!="poiseuille"&&problem!="rotation")||(open&&problem!="cylinder"))throw std::runtime_error("unsupported manufactured/boundary pairing");
    auto start=std::chrono::steady_clock::now();auto f=readFixture(name,n);const auto& mesh=f.mesh;int nc=int(mesh.cells.size()),nf=int(mesh.faces.size()),raw=4*nf+nc,count=0;Vec previous=readState(argv[9],f),outState(previous.size()),zero(previous.size());
    std::vector<int> map(raw,-1);Vec known(raw);int outlets=0;
    for(int i=0;i<nf;++i){const auto& face=mesh.faces[i];int kind=boundaryKind(face,open);outlets+=kind==1;
        for(int c=0;c<2;++c)for(int j=0;j<2;++j)if(face.neighbour||kind==1||(kind==2&&c==0))map[4*i+2*c+j]=count++;
        if(!face.neighbour&&problem=="cylinder"&&face.patch==BoundaryPatch2D::DomainBoundary&&kind==0)known[4*i]=1.;
        if(!face.neighbour&&problem!="cylinder"&&problem!="noslip"&&problem!="noslip-sheared")for(auto [z,w]:gauss(order)){double ss=z-.5;auto S=face.areaVector;auto u=exactAt({face.centre.x-ss*S.y,face.centre.y+ss*S.x},problem).u;
            known[4*i]+=w*u.x;known[4*i+1]+=12*w*ss*u.x;known[4*i+2]+=w*u.y;known[4*i+3]+=12*w*ss*u.y;}
    }
    if(open&&!outlets)throw std::runtime_error("no pressure-reference natural outlet");
    for(int t=0;t<nc-(open?0:1);++t)map[4*nf+t]=count++;
    Vec rhs(count),solution;std::ofstream entries,cells,faces;std::uint64_t nnz=0;
    if(mode=="assemble")entries.open(prefix+".entries.tmp",std::ios::binary);
    else{
        if(mode=="recover"){std::ifstream in(prefix+".solution",std::ios::binary);solution.resize(count);in.read(reinterpret_cast<char*>(solution.data()),count*sizeof(double));if(!in||in.peek()!=EOF)throw std::runtime_error("bad sparse solution");for(int i=0;i<raw;++i)if(map[i]>=0)known[i]=solution[map[i]];}
        else{for(int i=0;i<4*nf;++i)known[i]=previous[9*nc+i];for(int t=0;t<nc;++t)known[4*nf+t]=previous[9*t+6];}
        cells.open(prefix+"."+mode+".cells.csv.tmp");cells<<std::setprecision(17)<<"cell,x,y,area,u0,uX,uY,v0,vX,vY,p0,pX,pY\n";
    }
    Vec reaction(4*nf),viscousReaction(4*nf),convectionReaction(4*nf),stressReaction(4*nf),exactReaction(4*nf);double area=0,urms=0,prms=0,gauge=0,pressureWork=0,visc=0,conv=0,identity=0,force=0,maxdiv=0,maxspeed=0,pmin=1e300,pmax=-1e300,internal=0,trace=0,delta=0,rtRms=0,rtMax=0,pVertexMin=1e300,pVertexMax=-1e300;
    for(int t=0;t<nc;++t){Operator o(f,t,problem,nu,nonlinear,open,order,nonlinear?previous:zero);auto& e=o.e;const auto& a=e.a;int m=a.m;area+=mesh.cells[t].area;trace=std::max(trace,o.lift.traceResidual);
        std::vector<int> ids;for(int j:e.outside)if(j==2*m)ids.push_back(4*nf+t);else{int c=j/m,k=j%m;ids.push_back(4*int(mesh.cells[t].faces[(k-3)/2])+2*c+(k-3)%2);}
        if(mode=="assemble"){auto [k,b]=e.condensed();for(std::size_t i=0;i<ids.size();++i){int row=map[ids[i]];if(row<0)continue;rhs[row]+=b[i];for(std::size_t j=0;j<ids.size();++j){int col=map[ids[j]];double v=k(int(i),int(j));if(col<0)rhs[row]-=v*known[ids[j]];else if(v!=0){Entry z{row,col,v};entries.write(reinterpret_cast<char*>(&z),sizeof z);++nnz;}}}}
        else{Vec ext;for(int id:ids)ext.push_back(known[id]);auto v=mode=="check"?localState(f,t,m,previous):e.recover(ext);P6 ru{},rv{};for(int l=0;l<6;++l)for(int j=0;j<m;++j){ru[l]+=a.potential(l,j)*v[j];rv[l]+=a.potential(l,j)*v[m+j];}
            for(int c=0;c<2;++c)for(int j=0;j<3;++j)outState[9*t+3*c+j]=v[c*m+j];
            for(int j=0;j<3;++j)outState[9*t+6+j]=v[2*m+j];
            for(auto q:a.q){auto phi=a.basis.phi(q.p);auto theta=a.basis.theta(q.p);auto ex=problem=="cylinder"?Exact{}:exactAt(q.p,problem);double u=0,w=0,p=0,d=0;for(int j=0;j<6;++j){u+=theta[j]*ru[j];w+=theta[j]*rv[j];}for(int j=0;j<3;++j){p+=phi[j]*v[2*m+j];for(int l=0;l<m;++l)d+=phi[j]*(a.gx(j,l)*v[l]+a.gy(j,l)*v[m+l]);}
                urms+=q.w*((u-ex.u.x)*(u-ex.u.x)+(w-ex.u.y)*(w-ex.u.y));gauge+=q.w*(p-ex.p);prms+=q.w*(p-ex.p)*(p-ex.p);maxspeed=std::max(maxspeed,std::hypot(u,w));maxdiv=std::max(maxdiv,std::abs(d));pmin=std::min(pmin,p);pmax=std::max(pmax,p);}
            for(std::size_t k=0;k<o.lift.tri.size();++k){
                for(auto q:o.lift.tri[k].quadrature(order)){auto u=o.liftedVelocity(v,int(k),q.p);auto ex=problem=="cylinder"?Exact{}:exactAt(q.p,problem);rtRms+=q.w*((u.x-ex.u.x)*(u.x-ex.u.x)+(u.y-ex.u.y)*(u.y-ex.u.y));rtMax=std::max(rtMax,std::hypot(u.x,u.y));}
                for(auto point:{o.lift.tri[k].c,o.lift.tri[k].a,o.lift.tri[k].b}){auto phi=a.basis.phi(point);double pressure=0;for(int j=0;j<3;++j)pressure+=phi[j]*v[2*m+j];pVertexMin=std::min(pVertexMin,pressure);pVertexMax=std::max(pVertexMax,pressure);auto u=o.liftedVelocity(v,int(k),point);rtMax=std::max(rtMax,std::hypot(u.x,u.y));}
            }
            identity+=o.energyIdentity(v);
            for(int i=0;i<2*m;++i){force+=v[i]*e.rhs[i];for(int j=0;j<2*m;++j){conv+=v[i]*o.convection(i,j)*v[j];visc+=v[i]*(e.matrix(i,j)-o.convection(i,j))*v[j];}for(int p=0;p<3;++p)pressureWork+=v[i]*e.matrix(i,2*m+p)*v[2*m+p];}
            for(int i=0;i<2*m+3;++i){double r=-e.rhs[i];for(int j=0;j<2*m+3;++j)r+=e.matrix(i,j)*v[j];if(std::find(e.inside.begin(),e.inside.end(),i)!=e.inside.end())internal=std::max(internal,std::abs(r));else if(i<2*m){int c=i/m,j=i%m;int id=4*int(mesh.cells[t].faces[(j-3)/2])+2*c+(j-3)%2;reaction[id]+=r;
                    for(int l=0;l<2*m;++l){viscousReaction[id]+=(e.matrix(i,l)-o.convection(i,l))*v[l];convectionReaction[id]+=o.convection(i,l)*v[l];}}}
            for(std::size_t l=0;l<mesh.cells[t].faces.size();++l){int id=int(mesh.cells[t].faces[l]);const auto& face=mesh.faces[id];if(face.neighbour)continue;
                for(auto [z,w]:gauss(order)){double ss=z-.5;auto S=face.areaVector;auto phi=a.basis.phi({face.centre.x-ss*S.y,face.centre.y+ss*S.x});double ux=0,uy=0,vx=0,vy=0,p=0;
                    for(int k=0;k<3;++k){p+=phi[k]*v[2*m+k];for(int j=0;j<m;++j){ux+=phi[k]*a.gx(k,j)*v[j];uy+=phi[k]*a.gy(k,j)*v[j];vx+=phi[k]*a.gx(k,j)*v[m+j];vy+=phi[k]*a.gy(k,j)*v[m+j];}}
                    double tx=(2*nu*ux-p)*S.x+nu*(uy+vx)*S.y,ty=nu*(uy+vx)*S.x+(2*nu*vy-p)*S.y;
                    for(int j=0;j<2;++j){stressReaction[4*id+j]+=w*(j?ss:1)*tx;stressReaction[4*id+2+j]+=w*(j?ss:1)*ty;}
                    if(problem=="noslip"||problem=="noslip-sheared"){auto ex=exactAt({face.centre.x-ss*S.y,face.centre.y+ss*S.x},problem);auto gu=ex.gradient[0],gv=ex.gradient[1];double exx=(2*nu*gu.x-ex.p)*S.x+nu*(gu.y+gv.x)*S.y,exy=nu*(gu.y+gv.x)*S.x+(2*nu*gv.y-ex.p)*S.y;
                        for(int j=0;j<2;++j){exactReaction[4*id+j]+=w*(j?ss:1)*exx;exactReaction[4*id+2+j]+=w*(j?ss:1)*exy;}}}
            }
            cells<<t<<','<<mesh.cells[t].centre.x<<','<<mesh.cells[t].centre.y<<','<<mesh.cells[t].area;for(int j=0;j<9;++j)cells<<','<<outState[9*t+j];cells<<'\n';
        }
        if((t+1)%2000==0)std::cerr<<mode<<" cells "<<t+1<<'/'<<nc<<'\n';
    }
    double boundary=0,faceResidual=0,wallFx=0,wallFy=0,wallNormal=0,wallPressureFx=0,wallPressureFy=0,wallViscousFx=0,wallConvectionFx=0,wallStressFx=0,wallStressFy=0,wallTractionDifference=0,massFlux=0,momentumX=0,momentumY=0,tractionX=0,tractionY=0,fluxWork=0,allFluxWork=0,tractionError=0,boundaryLength=0;
    if(mode=="assemble"){entries.close();std::ofstream out(prefix+".rhs.tmp",std::ios::binary);out.write(reinterpret_cast<char*>(rhs.data()),rhs.size()*sizeof(double));out.close();if(!entries||!out)throw std::runtime_error("matrix write failed");std::filesystem::rename(prefix+".entries.tmp",prefix+".entries");std::filesystem::rename(prefix+".rhs.tmp",prefix+".rhs");}
    else{faces.open(prefix+"."+mode+".faces.csv.tmp");faces<<std::setprecision(17)<<"face,x,y,Sx,Sy,owner,neighbour,u0,us,v0,vs,reactionX,reactionY\n";
        for(int i=0;i<nf;++i){const auto& face=mesh.faces[i];faces<<i<<','<<face.centre.x<<','<<face.centre.y<<','<<face.areaVector.x<<','<<face.areaVector.y<<','<<face.owner<<','<<(face.neighbour?int(*face.neighbour):-1);
            for(int j=0;j<4;++j){outState[9*nc+4*i+j]=known[4*i+j];faces<<','<<known[4*i+j];if(map[4*i+j]>=0)faceResidual=std::max(faceResidual,std::abs(reaction[4*i+j]));else boundary+=known[4*i+j]*reaction[4*i+j];}faces<<','<<reaction[4*i]<<','<<reaction[4*i+2]<<'\n';
            if(!face.neighbour){
                if(problem=="noslip"||problem=="noslip-sheared"){double length=std::hypot(face.areaVector.x,face.areaVector.y);boundaryLength+=length;
                    for(int c=0;c<2;++c)for(int j=0;j<2;++j){double exact=exactReaction[4*i+2*c+j]-(j?0:gauge/area*(c?face.areaVector.y:face.areaVector.x));double error=reaction[4*i+2*c+j]-exact;tractionError+=(j?12:1)*error*error/length;}}
                tractionX+=reaction[4*i];tractionY+=reaction[4*i+2];
                for(auto [z,w]:gauss(order)){double ss=z-.5;Vector2D u{known[4*i]+ss*known[4*i+1],known[4*i+2]+ss*known[4*i+3]};double mass=dot(u,face.areaVector),flux=nonlinear?((previous[9*nc+4*i]+ss*previous[9*nc+4*i+1])*face.areaVector.x+(previous[9*nc+4*i+2]+ss*previous[9*nc+4*i+3])*face.areaVector.y):0;
                    massFlux+=w*mass;momentumX+=w*flux*u.x;momentumY+=w*flux*u.y;allFluxWork+=w*flux*dot(u,u);
                    if(nonlinear&&boundaryKind(face,open)!=1){tractionX+=w*flux*u.x;tractionY+=w*flux*u.y;fluxWork+=w*flux*dot(u,u);}
                    if(face.patch==BoundaryPatch2D::EmbeddedBoundary){int t=int(face.owner);Basis b{mesh.cells[t].centre,f.diameter[t],{}};auto phi=b.phi({face.centre.x-ss*face.areaVector.y,face.centre.y+ss*face.areaVector.x});double p=0;for(int j=0;j<3;++j)p+=phi[j]*outState[9*t+6+j];wallPressureFx+=w*p*face.areaVector.x;wallPressureFy+=w*p*face.areaVector.y;}
                }
            }
            if(!face.neighbour&&face.patch==BoundaryPatch2D::EmbeddedBoundary){wallFx-=reaction[4*i];wallFy-=reaction[4*i+2];wallViscousFx-=viscousReaction[4*i];wallConvectionFx-=convectionReaction[4*i];wallStressFx-=stressReaction[4*i];wallStressFy-=stressReaction[4*i+2];
                for(int j=0;j<4;++j)wallTractionDifference=std::max(wallTractionDifference,std::abs(reaction[4*i+j]-stressReaction[4*i+j])/std::hypot(face.areaVector.x,face.areaVector.y));
                for(double s:{-.5,.5})wallNormal=std::max(wallNormal,std::abs((known[4*i]+s*known[4*i+1])*face.areaVector.x+(known[4*i+2]+s*known[4*i+3])*face.areaVector.y)/std::hypot(face.areaVector.x,face.areaVector.y));}}
        for(std::size_t i=0;i<outState.size();++i)delta=std::max(delta,std::abs(outState[i]-previous[i]));
        writeState(prefix+"."+mode+".state",f,outState);cells.close();faces.close();if(!cells||!faces)throw std::runtime_error("field write failed");std::filesystem::rename(prefix+"."+mode+".cells.csv.tmp",prefix+"."+mode+".cells.csv");std::filesystem::rename(prefix+"."+mode+".faces.csv.tmp",prefix+"."+mode+".faces.csv");
    }
    if(mode!="assemble"&&problem!="cylinder"){
        gauge/=area;prms=0;
        for(int t=0;t<nc;++t){Basis b{mesh.cells[t].centre,f.diameter[t],{}};
            for(auto q:cellQuadrature(mesh,t,order)){auto phi=b.phi(q.p);double pressure=0;for(int j=0;j<3;++j)pressure+=phi[j]*outState[9*t+6+j];double error=pressure-exactAt(q.p,problem).p-gauge;prms+=q.w*error*error;}}
    }
    std::cout<<std::setprecision(17)<<"{\"mode\":\""<<mode<<"\",\"cells\":"<<nc<<",\"faces\":"<<nf<<",\"unknowns\":"<<count<<",\"entries\":"<<nnz<<",\"area\":"<<area<<",\"nu\":"<<nu<<",\"velocityRms\":"<<optionalMetric(std::sqrt(urms/area),problem!="cylinder")<<",\"pressureRms\":"<<optionalMetric(std::sqrt(prms/area),problem!="cylinder")<<",\"maxSpeed\":"<<maxspeed<<",\"pressureRange\":"<<(mode=="assemble"?0:pVertexMax-pVertexMin)<<",\"quadraturePressureRange\":"<<(mode=="assemble"?0:pmax-pmin)<<",\"rt1VelocityRms\":"<<optionalMetric(std::sqrt(rtRms/area),problem!="cylinder")<<",\"rt1SampledMaxSpeed\":"<<rtMax<<",\"maxDivergence\":"<<maxdiv<<",\"pressureWork\":"<<pressureWork<<",\"viscousEnergy\":"<<visc<<",\"convectionEnergy\":"<<conv<<",\"convectionIdentity\":"<<identity<<",\"convectionIdentityError\":"<<conv-identity<<",\"forceWork\":"<<force<<",\"boundaryWork\":"<<boundary<<",\"energyBalance\":"<<visc+conv+pressureWork-force-boundary<<",\"internalResidual\":"<<internal<<",\"freeFaceResidual\":"<<faceResidual<<",\"wallFx\":"<<wallFx<<",\"wallFy\":"<<wallFy<<",\"wallPressureFx\":"<<wallPressureFx<<",\"wallViscousAndStabilizationFx\":"<<wallViscousFx<<",\"wallConvectionReconstructionFx\":"<<wallConvectionFx<<",\"stressBodyFx\":"<<wallStressFx<<",\"stressBodyFy\":"<<wallStressFy<<",\"wallTractionMomentDifferencePerLength\":"<<wallTractionDifference<<",\"wallPressureFy\":"<<wallPressureFy<<",\"netBoundaryMassFlux\":"<<massFlux<<",\"boundaryMomentumFluxX\":"<<momentumX<<",\"boundaryMomentumFluxY\":"<<momentumY<<",\"boundaryTractionX\":"<<tractionX<<",\"boundaryTractionY\":"<<tractionY<<",\"globalMomentumBalanceX\":"<<optionalMetric(tractionX-(nonlinear?momentumX:0),problem=="cylinder")<<",\"globalMomentumBalanceY\":"<<optionalMetric(tractionY-(nonlinear?momentumY:0),problem=="cylinder")<<",\"physicalConvectionEnergy\":"<<conv+fluxWork<<",\"physicalBoundaryTractionWork\":"<<boundary+fluxWork<<",\"upwindDissipationPlusDivergenceWork\":"<<conv+fluxWork-.5*allFluxWork<<",\"manufacturedTractionMomentRms\":"<<optionalMetric(std::sqrt(tractionError/std::max(boundaryLength,1e-300)),problem=="noslip"||problem=="noslip-sheared")<<",\"wallNormal\":"<<wallNormal<<",\"stateCoefficientChange\":"<<delta<<",\"rtTraceResidual\":"<<trace<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}

#ifndef CARTMESH_P1_OSEEN_NO_MAIN
int main(int argc,char** argv){return runOseen<Oseen>(argc,argv);}
#endif
