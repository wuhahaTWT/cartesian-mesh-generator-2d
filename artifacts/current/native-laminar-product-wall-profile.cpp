// Research postprocessor for two saved compatible P1 flow states on the same
// fixed polygon wall. It compares every wall pressure trace after removing one
// volume-weighted pressure gauge from each complete fluid field. No cell or
// face is removed, and no pressure equation is assembled here.
#define CARTMESH_P1_OSEEN_NO_MAIN
#include "native-laminar-p1-oseen.cpp"

struct WallSample {
    const Face* face{};
    int owner{};
    double angle{},length{},constant{},slope{};
};

std::vector<WallSample> wallSamples(const Fixture& fixture,const Vec& state,double& gauge) {
    const auto& mesh=fixture.mesh;double area=0;
    gauge=0;
    for(std::size_t i=0;i<mesh.cells.size();++i){area+=mesh.cells[i].area;gauge+=mesh.cells[i].area*state[9*i+6];}
    gauge/=area;
    Point2D centre{};double lengthSum=0;
    for(const auto& face:mesh.faces)if(!face.neighbour&&face.patch==BoundaryPatch2D::EmbeddedBoundary){
        const double length=std::hypot(face.areaVector.x,face.areaVector.y);centre.x+=length*face.centre.x;centre.y+=length*face.centre.y;lengthSum+=length;
    }
    if(!(lengthSum>0))throw std::runtime_error("no embedded wall faces");
    centre.x/=lengthSum;centre.y/=lengthSum;std::vector<WallSample> result;
    for(const auto& face:mesh.faces)if(!face.neighbour&&face.patch==BoundaryPatch2D::EmbeddedBoundary){
        const int owner=int(face.owner);Basis basis{mesh.cells[owner].centre,fixture.diameter[owner],{}};const auto phi=basis.phi(face.centre);
        double constant=-gauge;for(int j=0;j<3;++j)constant+=phi[j]*state[9*owner+6+j];
        const double slope=(-face.areaVector.y*state[9*owner+7]+face.areaVector.x*state[9*owner+8])/basis.h;
        result.push_back({&face,owner,std::atan2(face.centre.y-centre.y,face.centre.x-centre.x),std::hypot(face.areaVector.x,face.areaVector.y),constant,slope});
    }
    std::sort(result.begin(),result.end(),[](const auto& a,const auto& b){return a.angle<b.angle;});return result;
}

int main(int argc,char** argv)try {
    if(argc!=8)throw std::runtime_error("usage: wall-profile meshA nA stateA meshB nB stateB output.csv");
    const std::string output=argv[7];if(std::filesystem::exists(output))throw std::runtime_error("wall profile output already exists");
    auto a=readFixture(argv[1],std::stoi(argv[2])),b=readFixture(argv[4],std::stoi(argv[5]));auto sa=readState(argv[3],a),sb=readState(argv[6],b);double ga=0,gb=0;auto wa=wallSamples(a,sa,ga),wb=wallSamples(b,sb,gb);
    const double geometryScale=std::max(1.,std::sqrt(std::accumulate(a.mesh.cells.begin(),a.mesh.cells.end(),0.,[](double x,const auto& c){return x+c.area;})));
    const double tolerance=512*std::numeric_limits<double>::epsilon()*geometryScale;
    double norm=0,maximum=0,overlapSum=0,lengthA=0,lengthB=0,forceA=0,forceB=0;std::size_t intervals=0;
    for(const auto& x:wa){lengthA+=x.length;forceA+=x.constant*x.face->areaVector.x;}
    for(const auto& x:wb){lengthB+=x.length;forceB+=x.constant*x.face->areaVector.x;}
    std::ofstream out(output+".tmp");out<<std::setprecision(17)<<"interval,length,x,y,pA0,pAs,pB0,pBs\n";
    auto rule=gauss(10);
    for(const auto& x:wa){const double lx=x.length,nx=x.face->areaVector.x/lx,ny=x.face->areaVector.y/lx,tx=-ny,ty=nx,ca=x.face->centre.x*tx+x.face->centre.y*ty;
        for(const auto& y:wb){const double ly=y.length,mx=y.face->areaVector.x/ly,my=y.face->areaVector.y/ly;
            if(1-(nx*mx+ny*my)>1024*std::numeric_limits<double>::epsilon())continue;
            const double dx=y.face->centre.x-x.face->centre.x,dy=y.face->centre.y-x.face->centre.y;if(std::abs(dx*nx+dy*ny)>tolerance)continue;
            const double cb=y.face->centre.x*tx+y.face->centre.y*ty,lo=std::max(ca-lx/2,cb-ly/2),hi=std::min(ca+lx/2,cb+ly/2);if(hi-lo<=tolerance)continue;
            const double length=hi-lo,mid=(lo+hi)/2,sa0=(mid-ca)/lx,sb0=(mid-cb)/ly,pa=x.constant+sa0*x.slope,pb=y.constant+sb0*y.slope,pas=x.slope*length/lx,pbs=y.slope*length/ly;
            for(auto [z,w]:rule){double s=z-.5,d=(pb+s*pbs)-(pa+s*pas);norm+=w*length*d*d;maximum=std::max(maximum,std::abs(d));}
            const double px=x.face->centre.x+(mid-ca)*tx,py=x.face->centre.y+(mid-ca)*ty;out<<intervals++<<','<<length<<','<<px<<','<<py<<','<<pa<<','<<pas<<','<<pb<<','<<pbs<<'\n';overlapSum+=length;
        }
    }
    if(std::abs(overlapSum-lengthA)>tolerance*wa.size()||std::abs(overlapSum-lengthB)>tolerance*wb.size())throw std::runtime_error("fixed wall coverage differs");
    out.close();if(!out)throw std::runtime_error("wall profile write failure");std::filesystem::rename(output+".tmp",output);
    std::cout<<std::setprecision(17)<<"{\"wallFacesA\":"<<wa.size()<<",\"wallFacesB\":"<<wb.size()<<",\"commonIntervals\":"<<intervals<<",\"wallLengthA\":"<<lengthA<<",\"wallLengthB\":"<<lengthB<<",\"commonLength\":"<<overlapSum<<",\"pressureGaugeA\":"<<ga<<",\"pressureGaugeB\":"<<gb<<",\"wallPressureDifferenceRms\":"<<std::sqrt(norm/overlapSum)<<",\"wallPressureDifferenceMaxAtQuadrature\":"<<maximum<<",\"pressureForceXA\":"<<forceA<<",\"pressureForceXB\":"<<forceB<<",\"pressureForceXDifference\":"<<forceB-forceA<<",\"geometryTolerance\":"<<tolerance<<"}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
