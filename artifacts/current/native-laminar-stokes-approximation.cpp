// Research-only actual-mesh P2 approximation comparison. Reuse the native
// reference, quadrature and potential; do not assemble a second PDE. The
// cellwise L2 optimum need not satisfy global velocity/flux constraints.
#define CARTMESH_POLYGON_STOKES_NO_MAIN
#include "native-laminar-polygon-stokes.cpp"
int main(int argc,char** argv)try{
    using namespace polygon_stokes;
    if(argc!=5)throw std::runtime_error("usage: approximation mesh state xy fresh_output");
    const auto fixture=readFixture(argv[1],0);const auto& mesh=fixture.mesh;const auto state=readState(argv[2],fixture);Reference reference(argv[3],.25,1);
    if(std::filesystem::exists(argv[4]))throw std::runtime_error("output exists");
    std::ofstream out(argv[4]);out<<std::setprecision(17)<<"cell,area,best_p2_error_squared,interpolated_potential_error_squared,actual_potential_error_squared,actual_minus_interpolated_squared\n";
    double area=0,bestSum=0,interpolatedSum=0,actualSum=0,diffSum=0,referenceSum=0,divergenceMax=0;
    for(std::size_t t=0;t<mesh.cells.size();++t){const auto& cell=mesh.cells[t];P1Local a(mesh,t,fixture.diameter[t],16);const auto m=a.m;Mat mass6(6,6);Vec bx(6),by(6),ix(3),iy(3),ux(m),uy(m);
        for(auto q:a.q){const auto e=reference.at(q.p);const auto theta=a.basis.theta(q.p);const auto phi=a.basis.phi(q.p);
            for(std::size_t j=0;j<6;++j){bx[j]+=q.w*theta[j]*e.u.x;by[j]+=q.w*theta[j]*e.u.y;for(std::size_t k=0;k<6;++k)mass6(j,k)+=q.w*theta[j]*theta[k];}
            for(std::size_t j=0;j<3;++j){ix[j]+=q.w*phi[j]*e.u.x;iy[j]+=q.w*phi[j]*e.u.y;}}
        DenseLU lu3(a.mass.v,3),lu6(mass6.v,6);ix=lu3.solve(ix);iy=lu3.solve(iy);bx=lu6.solve(bx);by=lu6.solve(by);
        for(std::size_t j=0;j<3;++j){ux[j]=ix[j];uy[j]=iy[j];}
        for(std::size_t l=0;l<cell.faces.size();++l){const auto& f=mesh.faces[cell.faces[l]];for(auto [z,w]:gauss(16)){const double s=z-.5;const auto e=reference.at({f.centre.x-s*f.areaVector.y,f.centre.y+s*f.areaVector.x});
            ux[3+2*l]+=w*e.u.x;ux[4+2*l]+=12*w*s*e.u.x;uy[3+2*l]+=w*e.u.y;uy[4+2*l]+=12*w*s*e.u.y;}}
        const auto solved=localState(fixture,int(t),int(m),state);P6 rx{},ry{},sx{},sy{};
        for(std::size_t j=0;j<6;++j)for(std::size_t k=0;k<m;++k){rx[j]+=a.potential(j,k)*ux[k];ry[j]+=a.potential(j,k)*uy[k];sx[j]+=a.potential(j,k)*solved[k];sy[j]+=a.potential(j,k)*solved[m+k];}
        double best=0,interp=0,actual=0,difference=0;
        for(auto q:a.q){const auto exact=reference.at(q.p);const auto theta=a.basis.theta(q.p);Vector2D b{},i{},s{};for(std::size_t j=0;j<6;++j){b.x+=theta[j]*bx[j];b.y+=theta[j]*by[j];i.x+=theta[j]*rx[j];i.y+=theta[j]*ry[j];s.x+=theta[j]*sx[j];s.y+=theta[j]*sy[j];}
            const auto be=b-exact.u,ie=i-exact.u,se=s-exact.u,de=s-i;best+=q.w*dot(be,be);interp+=q.w*dot(ie,ie);actual+=q.w*dot(se,se);difference+=q.w*dot(de,de);referenceSum+=q.w*dot(exact.u,exact.u);}
        for(auto id:cell.faces){const auto& f=mesh.faces[id];for(double s:{-.5,.5}){auto phi=a.basis.phi({f.centre.x-s*f.areaVector.y,f.centre.y+s*f.areaVector.x});double div=0;
            for(std::size_t j=0;j<3;++j)for(std::size_t k=0;k<m;++k)div+=phi[j]*(a.gx(j,k)*ux[k]+a.gy(j,k)*uy[k]);divergenceMax=std::max(divergenceMax,std::abs(div));}}
        area+=cell.area;bestSum+=best;interpolatedSum+=interp;actualSum+=actual;diffSum+=difference;out<<t<<','<<cell.area<<','<<best<<','<<interp<<','<<actual<<','<<difference<<'\n';
    }checkedClose(out);
    std::cout<<std::setprecision(17)<<"{\"cells\":"<<mesh.cells.size()<<",\"bestP2Rms\":"<<std::sqrt(bestSum/area)<<",\"interpolatedP2Rms\":"<<std::sqrt(interpolatedSum/area)<<",\"actualP2Rms\":"<<std::sqrt(actualSum/area)<<",\"actualMinusInterpolatedP2Rms\":"<<std::sqrt(diffSum/area)<<",\"exactVelocityRms\":"<<std::sqrt(referenceSum/area)<<",\"interpolatedWeakDivergenceMaximum\":"<<divergenceMax<<"}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
