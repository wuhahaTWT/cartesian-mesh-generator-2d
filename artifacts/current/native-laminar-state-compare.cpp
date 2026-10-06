// Native geometry and polynomial field differences for the coupled research solver.
// No independent PDE assembly, acceptance threshold or pressure branch selection.
#define CARTMESH_P1_OSEEN_NO_MAIN
#include "native-laminar-p1-oseen.cpp"
int main(int argc,char** argv)try {
    if(argc<5)throw std::runtime_error("usage: state-compare geometry mesh n prefix | compare mesh n first_state second_state");
    const std::string mode=argv[1];const auto f=readFixture(argv[2],std::stoi(argv[3]));const auto& mesh=f.mesh;
    if(mode=="geometry"&&argc==5){
        const std::string prefix=argv[4];if(std::filesystem::exists(prefix+".cells.csv"))throw std::runtime_error("geometry output exists");
        std::ofstream out(prefix+".cells.csv");out<<std::setprecision(17)<<"cell,area\n";double area=0;
        for(std::size_t t=0;t<mesh.cells.size();++t){out<<t<<','<<mesh.cells[t].area<<'\n';area+=mesh.cells[t].area;}
        out.close();if(!out)throw std::runtime_error("geometry output failed");
        std::cout<<std::setprecision(17)<<"{\"meshKey\":"<<meshKey(f)<<",\"cells\":"<<mesh.cells.size()<<",\"faces\":"<<mesh.faces.size()<<",\"area\":"<<area<<"}\n";return 0;
    }
    if(mode!="compare"||argc!=6)throw std::runtime_error("invalid state comparison request");
    const auto a=readState(argv[4],f),b=readState(argv[5],f);double area=0,gauge=0,velocity=0,pressure=0,pressureGauged=0,faceError=0,length=0,faceMax=0,coefficientMax=0;
    for(std::size_t j=0;j<a.size();++j)coefficientMax=std::max(coefficientMax,std::abs(a[j]-b[j]));
    for(std::size_t t=0;t<mesh.cells.size();++t){Basis basis{mesh.cells[t].centre,f.diameter[t],{}};area+=mesh.cells[t].area;
        for(auto q:cellQuadrature(mesh,int(t),4)){auto phi=basis.phi(q.p);double du=0,dv=0,dp=0;
            for(int j=0;j<3;++j){du+=phi[j]*(a[9*t+j]-b[9*t+j]);dv+=phi[j]*(a[9*t+3+j]-b[9*t+3+j]);dp+=phi[j]*(a[9*t+6+j]-b[9*t+6+j]);}
            velocity+=q.w*(du*du+dv*dv);pressure+=q.w*dp*dp;gauge+=q.w*dp;}}
    gauge/=area;
    for(std::size_t t=0;t<mesh.cells.size();++t){Basis basis{mesh.cells[t].centre,f.diameter[t],{}};
        for(auto q:cellQuadrature(mesh,int(t),4)){auto phi=basis.phi(q.p);double dp=-gauge;for(int j=0;j<3;++j)dp+=phi[j]*(a[9*t+6+j]-b[9*t+6+j]);pressureGauged+=q.w*dp*dp;}}
    for(std::size_t i=0;i<mesh.faces.size();++i){const auto& face=mesh.faces[i];double l=std::hypot(face.areaVector.x,face.areaVector.y);length+=l;std::array<double,4> d{};
        for(int j=0;j<4;++j)d[j]=a[9*mesh.cells.size()+4*i+j]-b[9*mesh.cells.size()+4*i+j];
        faceError+=l*(d[0]*d[0]+d[2]*d[2]+(d[1]*d[1]+d[3]*d[3])/12);
        for(double s:{-.5,.5})faceMax=std::max(faceMax,std::hypot(d[0]+s*d[1],d[2]+s*d[3]));}
    std::cout<<std::setprecision(17)<<"{\"cells\":"<<mesh.cells.size()<<",\"velocityP1Rms\":"<<std::sqrt(velocity/area)<<",\"pressureP1RmsAbsolute\":"<<std::sqrt(pressure/area)
        <<",\"pressureP1RmsAfterGlobalGauge\":"<<std::sqrt(pressureGauged/area)<<",\"pressureMeanDifference\":"<<gauge<<",\"faceVelocityRms\":"<<std::sqrt(faceError/length)
        <<",\"faceVelocityMaximum\":"<<faceMax<<",\"stateCoefficientMaximum\":"<<coefficientMax<<"}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
