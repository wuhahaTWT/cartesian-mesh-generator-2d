#include "cartmesh2d/fv/detail/CompatibleFlowTransport2D.hpp"
namespace cartmesh2d::fv::detail::compatible {
namespace {
std::vector<std::pair<double,double>> upwind(int order,double left,double right) {
    std::vector<double> ends{0,1};
    if((left<0&&right>0)||(left>0&&right<0))ends.insert(ends.begin()+1,left/(left-right));
    std::vector<std::pair<double,double>> out;const auto rule=gauss(order);
    for(std::size_t i=1;i<ends.size();++i)for(auto [z,w]:rule)
        out.push_back({ends[i-1]+(ends[i]-ends[i-1])*z,(ends[i]-ends[i-1])*w});
    return out;
}
Vector2D difference(Vector2D a,Vector2D b){return {a.x-b.x,a.y-b.y};}
struct Advector {
    const Lift& lift;std::vector<Vec> coefficients;
    Advector(const Lift& l,const Vec& beta,std::size_t m):lift(l) {
        if(beta.size()!=2*m+3)throw std::invalid_argument("Compatible advector dimensions");
        for(double x:beta)if(!std::isfinite(x))throw std::runtime_error("Nonfinite compatible advector");
        for(std::size_t k=0;k<lift.tri.size();++k){Vec c(8);
            for(std::size_t i=0;i<8;++i)for(std::size_t j=0;j<2*m;++j)c[i]+=lift.coefficients(8*k+i,j)*beta[j];
            coefficients.push_back(std::move(c));}
    }
    Vector2D at(std::size_t k,Point2D p)const {
        const auto shape=lift.tri[k].shape(p);Vector2D b{};
        for(std::size_t l=0;l<8;++l){b.x+=shape[l].x*coefficients[k][l];b.y+=shape[l].y*coefficients[k][l];}return b;
    }
};
std::pair<std::array<Vector2D,8>,std::array<Vector2D,8>> derivatives(Vector2D z) {
    return {{{{0,0},{1,0},{0,0},{0,0},{0,1},{0,0},{2*z.x,z.y},{z.y,0}}},
            {{{0,0},{0,0},{1,0},{0,0},{0,0},{0,1},{0,z.x},{z.x,2*z.y}}}};
}
void check(const FvMesh2D& mesh,std::size_t cell,const P1Local& a,const Lift& lift,const std::vector<bool>& open) {
    if(cell>=mesh.cells.size()||open.size()!=mesh.cells[cell].faces.size()||a.m!=3+2*mesh.cells[cell].faces.size()||lift.tri.size()!=mesh.cells[cell].faces.size())
        throw std::invalid_argument("Compatible transport topology dimensions");
}
}
Mat transportMatrix(const FvMesh2D& mesh,std::size_t cell,const P1Local& a,const Lift& lift,const Vec& beta,const std::vector<bool>& open,int order) {
    check(mesh,cell,a,lift,open);const auto m=a.m;const Advector advector(lift,beta,m);Mat convection(2*m,2*m);
    for(std::size_t k=0;k<lift.tri.size();++k){const auto& tr=lift.tri[k];
        for(auto q:tr.quadrature(order)){const auto z=tr.ref(q.p),b=advector.at(k,q.p);Vec r0(2*m),r1(2*m),d0(2*m),d1(2*m);const auto shape=tr.shape(q.p);const auto [dx,dy]=derivatives(z);
            const double bx=(tr.j1.y*b.x-tr.j1.x*b.y)/(tr.h*tr.det),by=(-tr.j0.y*b.x+tr.j0.x*b.y)/(tr.h*tr.det);
            for(std::size_t l=0;l<8;++l){Vector2D d{bx*dx[l].x+by*dy[l].x,bx*dx[l].y+by*dy[l].y};d={(tr.j0.x*d.x+tr.j1.x*d.y)/tr.det,(tr.j0.y*d.x+tr.j1.y*d.y)/tr.det};
                for(std::size_t j=0;j<2*m;++j){const double c=lift.coefficients(8*k+l,j);r0[j]+=c*shape[l].x;r1[j]+=c*shape[l].y;d0[j]+=c*d.x;d1[j]+=c*d.y;}}
            for(std::size_t i=0;i<2*m;++i)for(std::size_t j=0;j<2*m;++j)convection(i,j)-=q.w*(d0[i]*r0[j]+d1[i]*r1[j]);
        }
        const auto prev=(k+lift.tri.size()-1)%lift.tri.size();const auto d=tr.a-tr.c;const Vector2D S{d.y,-d.x};
        for(auto [z,w]:upwind(order,dot(advector.at(k,tr.c),S),dot(advector.at(k,tr.a),S))){const Point2D p{tr.c.x+z*d.x,tr.c.y+z*d.y};const double flux=dot(advector.at(k,p),S);std::vector<Vector2D> left(2*m),right(2*m);
            for(std::size_t j=0;j<2*m;++j){left[j]=lift.value(k,j,p);right[j]=lift.value(prev,j,p);}
            for(std::size_t i=0;i<2*m;++i)for(std::size_t j=0;j<2*m;++j)convection(i,j)+=w*flux*dot(difference(left[i],right[i]),flux>=0?left[j]:right[j]);
        }
        const auto l=lift.faceLocal[k],id=mesh.cells[cell].faces[l];const auto& face=mesh.faces[id];const double sign=face.owner==cell?1:-1;
        auto fluxEnd=[&](double s){return sign*((beta[3+2*l]+s*beta[4+2*l])*face.areaVector.x+(beta[m+3+2*l]+s*beta[m+4+2*l])*face.areaVector.y);};
        for(auto [z,w]:upwind(order,fluxEnd(-.5),fluxEnd(.5))){const double s=z-.5;const auto sf=face.areaVector;const Point2D p{face.centre.x-s*sf.y,face.centre.y+s*sf.x};const double flux=fluxEnd(s);std::vector<Vector2D> inside(2*m),facev(2*m);
            for(std::size_t j=0;j<2*m;++j)inside[j]=lift.value(k,j,p);
            facev[3+2*l].x=1;facev[4+2*l].x=s;facev[m+3+2*l].y=1;facev[m+4+2*l].y=s;
            for(std::size_t i=0;i<2*m;++i)for(std::size_t j=0;j<2*m;++j){convection(i,j)+=w*flux*dot(difference(inside[i],facev[i]),flux>=0?inside[j]:facev[j]);if(open[l])convection(i,j)+=w*flux*dot(facev[i],facev[j]);}
        }
    }
    return convection;
}
Mat transportAdvectorDerivative(const FvMesh2D& mesh,std::size_t cell,const P1Local& a,const Lift& lift,const Vec& beta,const std::vector<bool>& open,int order) {
    check(mesh,cell,a,lift,open);const auto m=a.m;const Advector advector(lift,beta,m);Mat derivative(2*m,2*m);
    for(std::size_t k=0;k<lift.tri.size();++k){const auto& tr=lift.tri[k];
        for(auto q:tr.quadrature(order)){const auto z=tr.ref(q.p),u=advector.at(k,q.p);const auto [dx,dy]=derivatives(z);std::vector<Vector2D> gx(2*m),gy(2*m),value(2*m);
            for(std::size_t j=0;j<2*m;++j){value[j]=lift.value(k,j,q.p);
                for(std::size_t l=0;l<8;++l){
                    const auto physical=[&](double x,double y){Vector2D d{x*dx[l].x+y*dy[l].x,x*dx[l].y+y*dy[l].y};return Vector2D{(tr.j0.x*d.x+tr.j1.x*d.y)/tr.det,(tr.j0.y*d.x+tr.j1.y*d.y)/tr.det};};
                    const auto x=physical(tr.j1.y/(tr.h*tr.det),-tr.j0.y/(tr.h*tr.det));const auto y=physical(-tr.j1.x/(tr.h*tr.det),tr.j0.x/(tr.h*tr.det));
                    const double c=lift.coefficients(8*k+l,j);gx[j].x+=c*x.x;gx[j].y+=c*x.y;gy[j].x+=c*y.x;gy[j].y+=c*y.y;
                }}
            for(std::size_t i=0;i<2*m;++i)for(std::size_t j=0;j<2*m;++j)derivative(i,j)-=q.w*(value[j].x*dot(gx[i],u)+value[j].y*dot(gy[i],u));
        }
        const auto prev=(k+lift.tri.size()-1)%lift.tri.size();const auto d=tr.a-tr.c;const Vector2D S{d.y,-d.x};
        for(auto [z,w]:upwind(order,dot(advector.at(k,tr.c),S),dot(advector.at(k,tr.a),S))){const Point2D p{tr.c.x+z*d.x,tr.c.y+z*d.y};const double flux=dot(advector.at(k,p),S);const auto u=flux>=0?advector.at(k,p):advector.at(prev,p);std::vector<Vector2D> left(2*m),jump(2*m);
            for(std::size_t j=0;j<2*m;++j){left[j]=lift.value(k,j,p);jump[j]=difference(left[j],lift.value(prev,j,p));}
            for(std::size_t i=0;i<2*m;++i)for(std::size_t j=0;j<2*m;++j)derivative(i,j)+=w*dot(left[j],S)*dot(jump[i],u);
        }
        const auto l=lift.faceLocal[k],id=mesh.cells[cell].faces[l];const auto& face=mesh.faces[id];const double sign=face.owner==cell?1:-1;const auto sf=face.areaVector;
        auto fluxEnd=[&](double s){return sign*((beta[3+2*l]+s*beta[4+2*l])*sf.x+(beta[m+3+2*l]+s*beta[m+4+2*l])*sf.y);};
        for(auto [z,w]:upwind(order,fluxEnd(-.5),fluxEnd(.5))){const double s=z-.5;const Point2D p{face.centre.x-s*sf.y,face.centre.y+s*sf.x};const double flux=fluxEnd(s);const Vector2D uf{beta[3+2*l]+s*beta[4+2*l],beta[m+3+2*l]+s*beta[m+4+2*l]};const auto u=flux>=0?advector.at(k,p):uf;std::vector<Vector2D> facev(2*m);
            facev[3+2*l].x=1;facev[4+2*l].x=s;facev[m+3+2*l].y=1;facev[m+4+2*l].y=s;
            for(std::size_t i=0;i<2*m;++i){const auto jump=difference(lift.value(k,i,p),facev[i]);
                for(std::size_t j=0;j<2*m;++j){const double df=sign*dot(facev[j],sf);derivative(i,j)+=w*df*dot(jump,u);if(open[l])derivative(i,j)+=w*df*dot(facev[i],uf);}}
        }
    }
    return derivative;
}
}
