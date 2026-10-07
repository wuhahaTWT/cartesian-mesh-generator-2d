// Research: conservative P1 wall traction and pressure recovery for stationary
// no-slip Stokes walls. Reuses the native coupled operator and provided solution.
#define CARTMESH_P1_STRESS_NO_MAIN
#include "native-laminar-p1-stress.cpp"

struct WallTraction {
    // Coefficients of t(s)=constant+s*slope, s in [-1/2,1/2].
    std::array<Vector2D,2> value{},withoutFaceLoad{};
};
WallTraction recoverWallTraction(const Element& e,const Vec& state,int localFace,double length) {
    WallTraction result;
    for(int component=0;component<2;++component)for(int moment=0;moment<2;++moment){
        const int row=component*e.a.m+3+2*localFace+moment;
        double residual=-e.rhs[row];
        for(std::size_t column=0;column<e.matrix.nc;++column)residual+=e.matrix(row,column)*state[column];
        const double scale=(moment?12.:1.)/length;
        double& value=component?result.value[moment].y:result.value[moment].x;
        double& omitted=component?result.withoutFaceLoad[moment].y:result.withoutFaceLoad[moment].x;
        value=scale*residual;omitted=scale*(residual+e.rhs[row]);
    }
    return result;
}
struct WallRecord {
    int id,owner;Point2D centre;double length;Vector2D normal;
    WallTraction traction;std::array<double,2> cellPressure{},exactPressure{},normalViscous{};
    std::array<Vector2D,2> exactTraction{};
};
int main(int argc,char** argv)try {
    if(argc==5&&std::string(argv[1])=="write-ring"){
        const int radial=std::stoi(argv[2]),angular=std::stoi(argv[3]);const std::string path=argv[4];
        if(radial<1||angular<8||radial>16||angular>128||std::filesystem::exists(path))
            throw std::runtime_error("invalid small ring fixture or existing output");
        std::vector<std::vector<Point2D>> nodes(radial+1,std::vector<Point2D>(angular));
        for(int j=0;j<=radial;++j)for(int i=0;i<angular;++i){double radius=.5+.5*j/radial,angle=2*std::acos(-1.)*i/angular;
            nodes[j][i]={radius*std::cos(angle),radius*std::sin(angle)};}
        std::vector<Polygon2D> polygons;
        for(int j=0;j<radial;++j)for(int i=0;i<angular;++i){int next=(i+1)%angular;
            polygons.push_back({{nodes[j][i],nodes[j+1][i],nodes[j+1][next],nodes[j][next]}});}
        auto topology=fromPolygons(polygons);const auto checked=makeFvMesh2D(topology);std::string error;
        if(!writeCm2dTopology(topology,path,&error))throw std::runtime_error(error);
        std::cout<<"{\"cells\":"<<checked.cells.size()<<",\"faces\":"<<checked.faces.size()<<"}\n";return 0;
    }
    if(argc!=9)throw std::runtime_error("usage: wall-recovery mesh n problem lambda lift|cell symmetric|laplace order existing_prefix");
    const std::string name=argv[1],problem=argv[3],prefix=argv[8];
    const int n=std::stoi(argv[2]),order=std::stoi(argv[7]);const double lambda=std::stod(argv[4]);
    const bool lifted=std::string(argv[5])=="lift",symmetric=std::string(argv[6])=="symmetric";
    if((std::string(argv[5])!="lift"&&std::string(argv[5])!="cell")||
       (std::string(argv[6])!="symmetric"&&std::string(argv[6])!="laplace")||order<4||order>12||!std::isfinite(lambda))
        throw std::runtime_error("invalid recovery controls");
    if(problem!="hydrostatic"&&!(problem=="noslip"&&name=="square")&&!(problem=="noslip-sheared"&&name=="sheared"))
        throw std::runtime_error("normal traction is pressure only for the validated stationary no-slip fixtures");
    if(std::filesystem::exists(prefix+".wall-recovery.csv"))throw std::runtime_error("wall output already exists; preserve evidence");
    const auto start=std::chrono::steady_clock::now();const auto fixture=readFixture(name,n);const auto& mesh=fixture.mesh;
    const int nc=int(mesh.cells.size()),nf=int(mesh.faces.size());std::vector<int> map(4*nf+nc,-1);Vec known(map.size());int count=0;
    for(int i=0;i<nf;++i)if(mesh.faces[i].neighbour)for(int j=0;j<4;++j)map[4*i+j]=count++;
    for(int i=0;i<nc-1;++i)map[4*nf+i]=count++;
    for(int i=0;i<nf;++i)if(!mesh.faces[i].neighbour)for(auto [z,w]:gauss(order)){
        const double s=z-.5;
        // Match the coupled driver: these validated stationary traces are exactly zero.
        const Vector2D u{0,0};
        known[4*i]+=w*u.x;known[4*i+1]+=12*w*s*u.x;known[4*i+2]+=w*u.y;known[4*i+3]+=12*w*s*u.y;
    }
    Vec solution(count);std::ifstream input(prefix+".solution",std::ios::binary);
    input.read(reinterpret_cast<char*>(solution.data()),count*sizeof(double));
    if(!input||input.peek()!=EOF)throw std::runtime_error("solution absent, truncated or wrong size");
    for(int i=0;i<count;++i)if(!std::isfinite(solution[i]))throw std::runtime_error("nonfinite solution");
    for(std::size_t i=0;i<map.size();++i)if(map[i]>=0)known[i]=solution[map[i]];
    double area=0,gauge=0,wallVelocity=0;Vector2D force{},reaction{};std::vector<WallRecord> walls;
    for(int t=0;t<nc;++t){
        Element e(fixture,t,problem,lambda,lifted,symmetric,order);const auto& a=e.a;const int m=a.m;Vec ext;
        for(int j:e.outside){int id;if(j==2*m)id=4*nf+t;else{int component=j/m,k=j%m;id=4*int(mesh.cells[t].faces[(k-3)/2])+2*component+(k-3)%2;}ext.push_back(known[id]);}
        const auto state=e.recover(ext);area+=mesh.cells[t].area;
        for(auto q:a.q){auto phi=a.basis.phi(q.p);auto exact=manufactured(q.p,problem,lambda);double pressure=0;
            for(int j=0;j<3;++j)pressure+=phi[j]*state[2*m+j];
            gauge+=q.w*(pressure-exact.p);force.x+=q.w*exact.f.x;force.y+=q.w*exact.f.y;}
        for(std::size_t fi=0;fi<mesh.cells[t].faces.size();++fi){
            const int id=int(mesh.cells[t].faces[fi]);const auto& face=mesh.faces[id];if(face.neighbour)continue;
            const double length=std::hypot(face.areaVector.x,face.areaVector.y);const Vector2D normal=face.areaVector*(1/length);
            WallRecord wall{id,t,face.centre,length,normal,recoverWallTraction(e,state,int(fi),length),{}, {}, {}};
            auto centrePhi=a.basis.phi(face.centre);for(int j=0;j<3;++j)wall.cellPressure[0]+=centrePhi[j]*state[2*m+j];
            wall.cellPressure[1]=(-face.areaVector.y*state[2*m+1]+face.areaVector.x*state[2*m+2])/a.basis.h;
            for(int mode=0;mode<2;++mode){
                const P3 phi=mode?P3{0,-face.areaVector.y/a.basis.h,face.areaVector.x/a.basis.h}:centrePhi;
                double ux=0,uy=0,vx=0,vy=0;
                for(int l=0;l<3;++l)for(int j=0;j<m;++j){ux+=phi[l]*a.gx(l,j)*state[j];uy+=phi[l]*a.gy(l,j)*state[j];
                    vx+=phi[l]*a.gx(l,j)*state[m+j];vy+=phi[l]*a.gy(l,j)*state[m+j];}
                wall.normalViscous[mode]=(symmetric?2.:1.)*(normal.x*normal.x*ux+normal.x*normal.y*(uy+vx)+normal.y*normal.y*vy);
            }
            reaction.x+=length*wall.traction.value[0].x;reaction.y+=length*wall.traction.value[0].y;
            for(auto [z,w]:gauss(order)){const double s=z-.5;const Point2D point{face.centre.x-face.areaVector.y*s,face.centre.y+face.areaVector.x*s};
                auto exact=manufactured(point,problem,lambda);const auto gu=exact.gradient[0],gv=exact.gradient[1];
                const Vector2D traction{(2*gu.x-exact.p)*normal.x+(gu.y+gv.x)*normal.y,(gu.y+gv.x)*normal.x+(2*gv.y-exact.p)*normal.y};
                for(int j=0;j<2;++j){const double factor=w*(j?12*s:1);wall.exactPressure[j]+=factor*exact.p;wall.exactTraction[j].x+=factor*traction.x;wall.exactTraction[j].y+=factor*traction.y;}
                wallVelocity=std::max(wallVelocity,std::hypot(known[4*id]+s*known[4*id+1],known[4*id+2]+s*known[4*id+3]));
            }
            walls.push_back(wall);
        }
    }
    gauge/=area;double lengthSum=0,pCell=0,pTraction=0,pNoLoad=0,pProjection=0,pMoment=0,pMax=0,cellMax=0,pStress=0,pStressMoment=0,pStressMax=0,tractionError=0,tractionMoment=0,shearError=0;
    std::ofstream output(prefix+".wall-recovery.csv");if(!output)throw std::runtime_error("wall output open failure");
    output<<std::setprecision(17)<<"face,owner,x,y,length,nx,ny,p_cell0,p_cells,p_reaction0,p_reactions,p_no_load0,p_no_loads,p_exact0,p_exacts,p_stress0,p_stresss,normal_viscous0,normal_viscouss,tx0,txs,ty0,tys,exact_tx0,exact_txs,exact_ty0,exact_tys\n";
    for(auto wall:walls){
        const auto nvec=wall.normal;wall.cellPressure[0]-=gauge;
        wall.traction.value[0].x+=gauge*nvec.x;wall.traction.value[0].y+=gauge*nvec.y;
        wall.traction.withoutFaceLoad[0].x+=gauge*nvec.x;wall.traction.withoutFaceLoad[0].y+=gauge*nvec.y;
        const std::array<double,2> pressure{-dot(wall.traction.value[0],nvec),-dot(wall.traction.value[1],nvec)};
        const std::array<double,2> noLoad{-dot(wall.traction.withoutFaceLoad[0],nvec),-dot(wall.traction.withoutFaceLoad[1],nvec)};
        const std::array<double,2> stressPressure{pressure[0]+wall.normalViscous[0],pressure[1]+wall.normalViscous[1]};
        lengthSum+=wall.length;
        for(auto [z,w]:gauss(order)){const double s=z-.5,weight=w*wall.length;const Point2D point{wall.centre.x-wall.length*nvec.y*s,wall.centre.y+wall.length*nvec.x*s};
            auto exact=manufactured(point,problem,lambda);const auto gu=exact.gradient[0],gv=exact.gradient[1];
            const Vector2D traction{(2*gu.x-exact.p)*nvec.x+(gu.y+gv.x)*nvec.y,(gu.y+gv.x)*nvec.x+(2*gv.y-exact.p)*nvec.y};
            const double pc=wall.cellPressure[0]+s*wall.cellPressure[1],pr=pressure[0]+s*pressure[1],pn=noLoad[0]+s*noLoad[1],pe=wall.exactPressure[0]+s*wall.exactPressure[1];
            pCell+=weight*(pc-exact.p)*(pc-exact.p);pTraction+=weight*(pr-exact.p)*(pr-exact.p);pNoLoad+=weight*(pn-exact.p)*(pn-exact.p);pProjection+=weight*(pe-exact.p)*(pe-exact.p);pMoment+=weight*(pr-pe)*(pr-pe);
            const double ps=stressPressure[0]+s*stressPressure[1];pStress+=weight*(ps-exact.p)*(ps-exact.p);pStressMoment+=weight*(ps-pe)*(ps-pe);pStressMax=std::max(pStressMax,std::abs(ps-exact.p));
            pMax=std::max(pMax,std::abs(pr-exact.p));cellMax=std::max(cellMax,std::abs(pc-exact.p));
            const Vector2D recovered{wall.traction.value[0].x+s*wall.traction.value[1].x,wall.traction.value[0].y+s*wall.traction.value[1].y};
            const Vector2D projected{wall.exactTraction[0].x+s*wall.exactTraction[1].x,wall.exactTraction[0].y+s*wall.exactTraction[1].y};
            const Vector2D error{recovered.x-traction.x,recovered.y-traction.y},projectionError{recovered.x-projected.x,recovered.y-projected.y};tractionError+=weight*dot(error,error);tractionMoment+=weight*dot(projectionError,projectionError);
            const double shear=-nvec.y*error.x+nvec.x*error.y;shearError+=weight*shear*shear;
        }
        output<<wall.id<<','<<wall.owner<<','<<wall.centre.x<<','<<wall.centre.y<<','<<wall.length<<','<<nvec.x<<','<<nvec.y;
        for(auto pair:{wall.cellPressure,pressure,noLoad,wall.exactPressure,stressPressure,wall.normalViscous})for(double value:pair)output<<','<<value;
        output<<','<<wall.traction.value[0].x<<','<<wall.traction.value[1].x<<','<<wall.traction.value[0].y<<','<<wall.traction.value[1].y<<','<<wall.exactTraction[0].x<<','<<wall.exactTraction[1].x<<','<<wall.exactTraction[0].y<<','<<wall.exactTraction[1].y<<'\n';
    }
    output.close();if(!output)throw std::runtime_error("wall output write failure");
    const auto rms=[&](double value){return std::sqrt(value/lengthSum);};
    std::cout<<std::setprecision(17)<<"{\"cells\":"<<nc<<",\"unknowns\":"<<count<<",\"wallFaces\":"<<walls.size()<<",\"wallLength\":"<<lengthSum<<",\"area\":"<<area<<",\"pressureGauge\":"<<gauge
        <<",\"cellPressureTraceRms\":"<<rms(pCell)<<",\"reactionPressureRms\":"<<rms(pTraction)<<",\"reactionPressureProjectionRms\":"<<rms(pMoment)<<",\"bestP1PressureRms\":"<<rms(pProjection)
        <<",\"omittedFaceLoadPressureRms\":"<<rms(pNoLoad)<<",\"reactionPressureMaxAtQuadrature\":"<<pMax<<",\"cellPressureMaxAtQuadrature\":"<<cellMax
        <<",\"stressPressureRms\":"<<rms(pStress)<<",\"stressPressureProjectionRms\":"<<rms(pStressMoment)<<",\"stressPressureMaxAtQuadrature\":"<<pStressMax
        <<",\"tractionRms\":"<<rms(tractionError)<<",\"tractionProjectionRms\":"<<rms(tractionMoment)<<",\"shearRms\":"<<rms(shearError)<<",\"wallVelocityMax\":"<<wallVelocity
        <<",\"forceBalanceNorm\":"<<std::hypot(force.x+reaction.x,force.y+reaction.y)<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
