// Explicit static-stress and pseudo-traction outlets for the coupled solver.
// G_ij = d u_i / d x_j. Prescribed physical traction is (nu(G+G^T)-pI)n.
// A pseudo-traction (nu G-pI)n=-p_D n becomes physical traction
// nu G^T n-p_D n. The G^T term is implicit, and NOT silently zeroed.
// This weak vector condition is not separate strong p=p_D and d_n u=0.
#define CARTMESH_P1_TRANSPORT_NO_MAIN
#include "native-laminar-p1-transport.cpp"

enum class OutletForm { ExactTraction, PseudoTraction, NormalStress };
template<OutletForm form> struct ExplicitOpenBoundary : StandardOseenBoundary {
    static int model(const std::string& name){
        const std::string expected=form==OutletForm::ExactTraction?"traction":form==OutletForm::PseudoTraction?"pseudo-traction":"normal-stress";
        if(name!=expected)throw std::runtime_error("outlet operator/policy mismatch");
        return form==OutletForm::PseudoTraction?2:1;
    }
    static bool allows(const std::string&,int){return true;}
    static int kind(const Face& face,int outletModel,const std::string& problem){
        const int result=boundaryKind(face,outletModel,problem);
        // Channel controls prescribe full analytic velocity on the other
        // boundaries; the cylinder keeps its original horizontal symmetry.
        return result==2&&problem!="cylinder"?0:result;
    }
};
template<OutletForm form> struct OpenTransport : Transport {
    OpenTransport(const Fixture& f,int t,const std::string& problem,double nu,bool nonlinear,int outletModel,int order,const Vec& previous):
        Transport(f,t,problem,nu,nonlinear,outletModel,order,previous){
        const auto& a=e.a;const int m=a.m;
        for(std::size_t l=0;l<f.mesh.cells[t].faces.size();++l){const auto& face=f.mesh.faces[f.mesh.cells[t].faces[l]];
            if(face.neighbour||boundaryKind(face,outletModel,problem)!=1)continue;
            const auto S=face.areaVector;
            for(auto [z,w]:gauss(order)){const double s=z-.5;Point2D p{face.centre.x-s*S.y,face.centre.y+s*S.x};
                const auto exact=problem=="cylinder"?Exact{}:exactAt(p,problem);
                Vector2D traction{-exact.p*S.x,-exact.p*S.y};
                if constexpr(form==OutletForm::ExactTraction){const auto gu=exact.gradient[0],gv=exact.gradient[1];
                    traction.x+=nu*(2*gu.x*S.x+(gu.y+gv.x)*S.y);traction.y+=nu*((gu.y+gv.x)*S.x+2*gv.y*S.y);}
                for(int c=0;c<2;++c)for(int k=0;k<2;++k){const int row=c*m+3+2*l+k;const double test=w*(k?s:1);
                    outletLoad[row]+=test*(c?traction.y:traction.x);
                }
            }}
        // The shared Oseen operator already assembled the implicit G^T term
        // for model 2. Add only the prescribed load here, exactly once.
        for(int i=0;i<2*m;++i)e.rhs[i]+=outletLoad[i];
        // Only retained face-test rows changed. Interior block, interior-to-
        // trace block and interior RHS are unchanged, so their elimination
        // already computed by Transport remains exact for this operator.
    }
};
int main(int argc,char** argv)try{
    if(argc==4&&std::string(argv[1])=="write-warped"){
        const int n=std::stoi(argv[2]);const std::string path=argv[3];
        if(n<2||n>32||std::filesystem::exists(path))throw std::runtime_error("invalid small fixture or existing output");
        const double pi=std::acos(-1.);std::vector<std::vector<Point2D>> nodes(n+1,std::vector<Point2D>(n+1));
        for(int j=0;j<=n;++j)for(int i=0;i<=n;++i){const double x=double(i)/n,y=double(j)/n;nodes[j][i]={x,y};
            if(i>0&&i<n&&j>0&&j<n){nodes[j][i].x+=.09*std::sin(2*pi*x)*std::sin(2*pi*y);nodes[j][i].y+=.06*std::sin(pi*x)*std::sin(2*pi*y);}}
        std::vector<Polygon2D> polygons;
        for(int j=0;j<n;++j)for(int i=0;i<n;++i)polygons.push_back({{nodes[j][i],nodes[j][i+1],nodes[j+1][i+1],nodes[j+1][i]}});
        const auto topology=fromPolygons(polygons);const auto checked=makeFvMesh2D(topology);std::string error;
        if(!writeCm2dTopology(topology,path,&error))throw std::runtime_error(error);
        double area=0;for(const auto& c:checked.cells)area+=c.area;
        std::cout<<std::setprecision(17)<<"{\"cells\":"<<checked.cells.size()<<",\"faces\":"<<checked.faces.size()<<",\"area\":"<<area<<"}\n";return 0;
    }

    if(argc>7){const std::string name=argv[7];
        if(name=="traction")return runOseen<OpenTransport<OutletForm::ExactTraction>,ExplicitOpenBoundary<OutletForm::ExactTraction>>(argc,argv);
        if(name=="pseudo-traction")return runOseen<OpenTransport<OutletForm::PseudoTraction>,ExplicitOpenBoundary<OutletForm::PseudoTraction>>(argc,argv);
        if(name=="normal-stress")return runOseen<OpenTransport<OutletForm::NormalStress>,ExplicitOpenBoundary<OutletForm::NormalStress>>(argc,argv);
    }
    return runOseen<Transport>(argc,argv);
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
