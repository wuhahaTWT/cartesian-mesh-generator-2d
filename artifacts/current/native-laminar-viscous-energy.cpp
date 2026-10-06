#include "cartmesh2d/fv/FvMesh2D.hpp"
#include "cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

// Research restriction of the actual product diffusion flux to wall cells and
// whole graph rings. All omitted perturbation values are zero, not deleted
// physical cells. Homogeneous Dirichlet perturbations on every outer face;
// a compact wall-supported restriction does not reach the far-field boundary.
// Matrix B=M^(-1/2) K M^(-1/2), K=-integrated outward viscous momentum flux,
// at nu=1 m^2/s. Positive eigenvalues of (B+B^T)/2 certify instantaneous
// kinetic-energy production, NOT an incompressible eigenmode or branch verdict.
int main(int argc,char** argv) {
    using namespace cartmesh2d;
    using namespace cartmesh2d::fv;
    if(argc!=5) throw std::runtime_error("usage: viscous-energy mesh prefix rings symmetric|laplacian|orthogonal|pressure|affine");
    const auto read=readCm2dTopology(argv[1]);
    if(!read.valid()) throw std::runtime_error(read.error);
    const auto m=makeFvMesh2D(read.topology);
    const std::size_t n=m.cells.size(), nf=m.faces.size();
    const int rings=std::stoi(argv[3]); const std::string mode=argv[4];
    if(rings<0 || (mode!="symmetric" && mode!="laplacian" && mode!="orthogonal" && mode!="pressure" && mode!="affine"))
        throw std::runtime_error("invalid experiment control");
    std::vector<bool> chosen(n,false),fixed(nf,true);
    for(const auto& f:m.faces) if(!f.neighbour && f.patch==BoundaryPatch2D::EmbeddedBoundary) chosen[f.owner]=true;
    for(int ring=0;ring<rings;++ring) {
        auto next=chosen;
        for(const auto& f:m.faces) if(f.neighbour && (chosen[f.owner] || chosen[*f.neighbour]))
            next[f.owner]=next[*f.neighbour]=true;
        chosen=std::move(next);
    }
    std::vector<std::size_t> cells;
    for(std::size_t i=0;i<n;++i) if(chosen[i]) cells.push_back(i);
    if(cells.empty()) throw std::runtime_error("no embedded-wall cells");
    if(mode=="pressure" || mode=="affine") {
        double xmax=-1e300;
        for(const auto& f:m.faces)if(!f.neighbour)xmax=std::max(xmax,f.centre.x);
        for(std::size_t id=0;id<nf;++id)fixed[id]=!m.faces[id].neighbour && m.faces[id].centre.x==xmax;
    }
    const auto stencil=detail::buildFlowGradientStencil2D(m,fixed,mode=="pressure" || mode=="affine");
    if(mode=="affine") {
        std::ofstream out(std::string(argv[2])+".affine.csv");
        out<<std::setprecision(17)<<"axis,cell,area,productGradientError,centralGradientError\n";
        double area=0;for(const auto& c:m.cells)area+=c.area;
        std::cout<<std::setprecision(17)<<"[";
        for(int axis=0;axis<2;++axis) {
            std::vector<double> p(n),bc(nf);
            for(std::size_t i=0;i<n;++i)p[i]=axis?m.cells[i].centre.y:m.cells[i].centre.x;
            for(std::size_t id=0;id<nf;++id)bc[id]=axis?m.faces[id].centre.y:m.faces[id].centre.x;
            const auto gp=stencil.apply(p,bc);
            auto pf=detail::pressureFaceValues(m,p,gp,bc,fixed);
            const auto original=detail::conservativePressureGradient(m,pf);
            for(std::size_t id=0;id<nf;++id){const auto& f=m.faces[id];pf[id]=f.neighbour?.5*(p[f.owner]+p[*f.neighbour]):(fixed[id]?bc[id]:p[f.owner]);}
            const auto central=detail::conservativePressureGradient(m,pf);
            double rms[2]{},maximum[2]{};
            for(std::size_t i=0;i<n;++i) {
                const double e0=std::hypot(original[i].x-(axis==0),original[i].y-(axis==1));
                const double e1=std::hypot(central[i].x-(axis==0),central[i].y-(axis==1));
                rms[0]+=m.cells[i].area*e0*e0;rms[1]+=m.cells[i].area*e1*e1;
                maximum[0]=std::max(maximum[0],e0);maximum[1]=std::max(maximum[1],e1);
                out<<axis<<','<<i<<','<<m.cells[i].area<<','<<e0<<','<<e1<<'\n';
            }
            if(axis)std::cout<<",";
            std::cout<<"{\"axis\":"<<axis<<",\"productRms\":"<<std::sqrt(rms[0]/area)<<",\"productMax\":"<<maximum[0]
                <<",\"centralRms\":"<<std::sqrt(rms[1]/area)<<",\"centralMax\":"<<maximum[1]<<"}";
        }
        out.close();if(!out)throw std::runtime_error("output failed");std::cout<<"]\n";return 0;
    }
    std::ofstream map(std::string(argv[2])+".cells.csv"),out(std::string(argv[2])+".matrix.csv");
    if(!map || !out)throw std::runtime_error("output open failed");
    map<<std::setprecision(17)<<"local,cell,x,y,area\n";
    for(std::size_t k=0;k<cells.size();++k){const auto& c=m.cells[cells[k]];map<<k<<','<<cells[k]<<','<<c.centre.x<<','<<c.centre.y<<','<<c.area<<'\n';}
    out<<std::setprecision(17)<<"row,col,value\n";
    const std::vector<double> zero(nf);
    std::vector<double> u(n),v(n);
    if(mode=="pressure") {
        // Geometric pressure-flux stabilization in the product predicted flux,
        // with a uniform positive momentum response rAU=1 s. This is not the
        // solution-dependent full saddle-point Jacobian. S=div(interpolate(Gp)
        // dot Sface - compact normal pressure flux); positive S is dissipative
        // only when velocity divergence and pressure force are adjoints.
        for(std::size_t col=0;col<cells.size();++col) {
            const auto source=cells[col];u[source]=1/std::sqrt(m.cells[source].area);
            const auto gp=stencil.apply(u,zero);
            const auto gc=detail::conservativePressureGradient(m,detail::pressureFaceValues(m,u,gp,zero,fixed));
            std::vector<double> divergence(n);
            for(std::size_t id=0;id<nf;++id) {
                const auto& f=m.faces[id];const auto i=f.owner;
                if(!f.neighbour && !fixed[id])continue;
                auto g=gp[i],cg=gc[i];double other=0;
                if(f.neighbour) {const auto j=*f.neighbour;const double w=f.neighbourWeight;other=u[j];
                    g={g.x*(1-w)+gp[j].x*w,g.y*(1-w)+gp[j].y*w};
                    cg={cg.x*(1-w)+gc[j].x*w,cg.y*(1-w)+gc[j].y*w};}
                const double q=dot(cg,f.areaVector)-f.transmissibility*(other-u[i])-dot(g,f.correction);
                divergence[i]+=q;if(f.neighbour)divergence[*f.neighbour]-=q;
            }
            for(std::size_t row=0;row<cells.size();++row) {
                const auto i=cells[row];const double a=divergence[i]/std::sqrt(m.cells[i].area);
                if(a!=0)out<<row<<','<<col<<','<<a<<'\n';
            }
            u[source]=0;
        }
        map.close();out.close();if(!map || !out)throw std::runtime_error("output write failed");
        std::cout<<"cells="<<n<<" support="<<cells.size()<<" rings="<<rings<<" mode="<<mode<<'\n';
        return 0;
    }
    for(std::size_t col=0;col<2*cells.size();++col) {
        const auto source=cells[col/2];
        auto& field=(col%2)?v:u; field[source]=1/std::sqrt(m.cells[source].area);
        const auto gu=stencil.apply(u,zero),gv=stencil.apply(v,zero);
        const auto stress=mode=="symmetric"
            ?detail::symmetricViscousCorrection(m,u,v,gu,gv,zero,zero,fixed,fixed,fixed,fixed,1.)
            :std::vector<Vector2D>(nf);
        std::vector<Vector2D> integrated(n);
        for(std::size_t id=0;id<nf;++id) {
            const auto& f=m.faces[id]; const auto i=f.owner;
            const auto component=[&](const auto& value,const auto& gradient) {
                auto g=gradient[i];double other=0;
                if(f.neighbour){const auto j=*f.neighbour;other=value[j];
                    g={g.x*(1-f.neighbourWeight)+gradient[j].x*f.neighbourWeight,
                       g.y*(1-f.neighbourWeight)+gradient[j].y*f.neighbourWeight};}
                return f.transmissibility*(other-value[i])+(mode=="orthogonal"?0.:dot(g,f.correction));
            };
            const Vector2D flux{component(u,gu)-stress[id].x,component(v,gv)-stress[id].y};
            integrated[i].x+=flux.x; integrated[i].y+=flux.y;
            if(f.neighbour){integrated[*f.neighbour].x-=flux.x;integrated[*f.neighbour].y-=flux.y;}
        }
        for(std::size_t k=0;k<cells.size();++k) {
            const auto i=cells[k];const auto a=integrated[i]*(1/std::sqrt(m.cells[i].area));
            if(a.x!=0)out<<2*k<<','<<col<<','<<a.x<<'\n';
            if(a.y!=0)out<<2*k+1<<','<<col<<','<<a.y<<'\n';
        }
        field[source]=0;
    }
    map.close();out.close();if(!map || !out)throw std::runtime_error("output write failed");
        std::cout<<"cells="<<n<<" support="<<cells.size()<<" rings="<<rings<<" mode="<<mode<<'\n';
}
