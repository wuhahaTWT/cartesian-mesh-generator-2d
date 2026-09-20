#include "cartmesh2d/grid/RectilinearMesh2D.hpp"
#include "cartmesh2d/fv/FvMesh2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace cartmesh2d;
int main(int argc,char** argv){try{
    if(argc!=5)throw std::runtime_error("usage: rectilinear_probe nx ny stretch prefix");
    const auto parse=[](const char* text){std::size_t end=0;const double v=std::stod(text,&end);
        if(end!=std::string(text).size()||!std::isfinite(v))throw std::runtime_error("invalid number");return v;};
    const double xd=parse(argv[1]),yd=parse(argv[2]),a=parse(argv[3]);
    if(xd<2||yd<2||xd>256||yd>256||std::floor(xd)!=xd||std::floor(yd)!=yd||a<0||a>20)
        throw std::runtime_error("probe requires integer dimensions 2..256 and stretch 0..20");
    const auto nx=static_cast<std::size_t>(xd),ny=static_cast<std::size_t>(yd);
    std::vector<double>x(nx+1),y(ny+1);
    for(std::size_t i=0;i<=nx;++i)x[i]=double(i)/double(nx);
    for(std::size_t j=0;j<=ny;++j)y[j]=a==0?double(j)/double(ny):std::expm1(a*double(j)/double(ny))/std::expm1(a);
    const auto t=makeRectilinearMesh2D(x,y);const auto mesh=fv::makeFvMesh2D(t);(void)mesh;
    std::string error;const std::string p=argv[4];
    if(!writeCm2dTopology(t,p+".cm2d",&error)||!writeLegacyVtk2D(t,p+".vtk",&error))throw std::runtime_error(error);
    std::ofstream out(p+".quality.json");out<<solverQualityReportToJson(evaluateSolverQuality2D(t));
    if(!out)throw std::runtime_error("quality export failed");
    std::cout<<"cells="<<t.cells.size()<<" firstHeight="<<y[1]<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
