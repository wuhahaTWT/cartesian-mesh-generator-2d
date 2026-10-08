// Native linear response attribution of one production Stokes SIMPLE step.
// This diagnostic retains every operator; it never deletes a discretization term.
#define main cellStepTraceMain
#include "cell_step.cpp"
#undef main
#include "FlowPressure2D.hpp"
#include <sstream>

static std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> a;std::stringstream s(line);std::string x;
    while(std::getline(s,x,','))a.push_back(x);return a;
}
static FlowInitialGuess2D readMode(const std::filesystem::path& dir,const FvMesh2D& m) {
    FlowInitialGuess2D x;x.u.resize(m.cells.size());x.v.resize(m.cells.size());x.p.resize(m.cells.size());x.flux.resize(m.faces.size());
    std::ifstream a(dir/"input-modal-cells.csv"),f(dir/"input-modal-faces.csv");
    ensure(bool(a)&&bool(f),"Cannot read captured modal fields");std::string line;std::getline(a,line);std::size_t count=0;
    while(std::getline(a,line)){const auto z=split(line);const auto i=std::stoul(z[0]);ensure(i==count++&&z.size()==8,"Invalid captured cell order");x.u[i]=std::stod(z[5]);x.v[i]=std::stod(z[6]);x.p[i]=std::stod(z[7]);}
    ensure(count==m.cells.size(),"Incomplete captured cells");std::getline(f,line);count=0;
    while(std::getline(f,line)){const auto z=split(line);const auto i=std::stoul(z[0]);ensure(i==count++&&z.size()==4,"Invalid captured face order");x.flux[i]=std::stod(z[3]);}
    ensure(count==m.faces.size(),"Incomplete captured faces");return x;
}
struct Response {std::string name;FlowInitialGuess2D predictor,increment,final;Vec rhsU,rhsV,q;};
int main(int argc,char** argv) {
    try {
        ensure(argc==5,"Usage: source_response MESH BOUNDARIES TRACE_DIR OUTPUT_DIR");
        const auto input=readCm2dTopology(argv[1]);ensure(input.valid(),"Invalid mesh");const auto m=makeFvMesh2D(input.topology);const auto n=m.cells.size(),nf=m.faces.size();
        const std::filesystem::path dir=argv[4];std::filesystem::create_directories(dir);
        FlowControls2D c;c.scenario="custom";c.nu=.1;c.speed=.5;c.tolerance=1e-8;c.maxIterations=1;c.momentumInertia=0;c.pressureCorrectionPasses=4;
        c.pressurePreconditioner=PressurePreconditioner2D::IncompleteCholesky0;c.convection=ConvectionScheme2D::FaceLimitedLinearUpwind;c.viscousStress=ViscousStress2D::Symmetric;
        std::ifstream bcFile(argv[2]);ensure(bool(bcFile),"Cannot read boundaries");c.boundaryConditions=readFlowBoundaryConditions2D(bcFile,m,c);const auto b=boundaries(m,c);
        ensure(b.closed,"Source attribution requires closed fixed-velocity annulus");for(std::size_t id=0;id<nf;++id)ensure(m.faces[id].neighbour || (b.fixedU[id]&&b.fixedV[id]),"Source attribution requires all velocity wall traces fixed");
        auto homogeneous=b;std::fill(homogeneous.u.begin(),homogeneous.u.end(),0);std::fill(homogeneous.v.begin(),homogeneous.v.end(),0);std::fill(homogeneous.p.begin(),homogeneous.p.end(),0);
        FlowInitialGuess2D zero;zero.u.resize(n);zero.v.resize(n);zero.p.resize(n);zero.flux.resize(nf);const auto x=readMode(argv[3],m);
        std::ifstream summary(std::filesystem::path(argv[3])/"summary.csv");ensure(bool(summary),"Cannot read native trace summary");std::string summaryLine;std::getline(summary,summaryLine);std::getline(summary,summaryLine);const auto targetCell=std::stoul(split(summaryLine).at(2));ensure(targetCell<n,"Captured target outside mesh");
        Trace t,z;active=&t;const auto actual=solveIncompressibleFromGuess2D(m,c,x);active=&z;const auto base=solveIncompressibleFromGuess2D(m,c,zero);active=nullptr;
        ensure(t.momentum.size()==4&&z.momentum.size()==4&&t.rc.size()==1&&z.rc.size()==1,"Unexpected native production call sequence");
        for(std::size_t call=0;call<4;++call){ensure(t.momentum[call].diag==z.momentum[call].diag&&t.momentum[call].off==z.momentum[call].off,"Source attribution requires identical modal/base matrices");for(const auto* trace:{&t.momentum[call],&z.momentum[call]})for(const auto& source:trace->source)ensure(source.x==0&&source.y==0,"Unexpected nonzero volume source");}
        ensure(t.rc[0].ra==z.rc[0].ra&&t.rc[0].df==z.rc[0].df,"Source attribution requires identical pressure responses");
        std::vector<std::pair<std::size_t,std::size_t>> edges;for(const auto& f:m.faces)if(f.neighbour)edges.emplace_back(f.owner,*f.neighbour);
        const detail::SparsePattern2D pattern(n,edges);ensure(pattern.rows==t.momentum[0].rows&&pattern.columns==t.momentum[0].columns,"Captured matrix graph differs");
        System au(pattern),av(pattern),ap(pattern);detail::LinearWorkspace2D work(n);
        const auto pressureStencil=detail::buildFlowGradientStencil2D(m,b.fixedP,true),uStencil=detail::buildFlowGradientStencil2D(m,b.fixedU),vStencil=detail::buildFlowGradientStencil2D(m,b.fixedV);
        Vec zeros(nf),cellZeros(n);Grad gradientZeros(n);const auto pf=detail::pressureFaceValues(m,x.p,pressureStencil.apply(x.p,b.p),b.p,b.fixedP),pf0=detail::pressureFaceValues(m,cellZeros,gradientZeros,b.p,b.fixedP);
        const std::vector<std::string> names{"wall_nonorth_viscosity","wall_symmetric_stress","internal_nonorth_viscosity","internal_symmetric_stress","momentum_gauss_pressure","momentum_underrelaxation","rc_direct_gauss_pressure","rc_compact_pressure","rc_nonorth_pressure","rc_flux_underrelaxation","old_pressure_carry"};
        std::vector<Response> responses;for(const auto& name:names)responses.push_back({name,zero,zero,zero,Vec(n),Vec(n),Vec(nf)});
        for(int component=0;component<2;++component) {
            const auto& a=t.momentum[component];const auto& q=z.momentum[component];
            for(std::size_t i=0;i<n;++i) {
                auto& pressure=component?responses[4].rhsV:responses[4].rhsU;
                auto& relaxation=component?responses[5].rhsV:responses[5].rhsU;
                pressure[i]=-m.cells[i].area*((component?a.force[i].y:a.force[i].x)-(component?q.force[i].y:q.force[i].x));
                const double common=std::max(t.momentum[0].diag[i],t.momentum[1].diag[i])/c.velocityRelaxation;
                relaxation[i]=(common-a.diag[i])*a.field[i]-(common-q.diag[i])*q.field[i];
            }
            for(std::size_t id=0;id<nf;++id) {
                const auto& f=m.faces[id];
                for(const auto cell:{f.owner,f.neighbour?*f.neighbour:f.owner}) {
                    const auto terms=faceTerms(m,c,b,a,pf,id,cell),before=faceTerms(m,c,b,q,pf0,id,cell);
                    const auto nonorth=f.neighbour?2UL:0UL,stress=f.neighbour?3UL:1UL;
                    auto& nonorthSource=component?responses[nonorth].rhsV:responses[nonorth].rhsU;
                    auto& stressSource=component?responses[stress].rhsV:responses[stress].rhsU;
                    nonorthSource[cell]+=(terms.wall-before.wall)+(terms.nonorth-before.nonorth);
                    stressSource[cell]+=terms.stress-before.stress;
                    if(!f.neighbour)break;
                }
            }
        }
        for(std::size_t id=0;id<nf;++id) {
            const auto a=rcTerms(m,c,t.rc[0],id),q=rcTerms(m,c,z.rc[0],id);
            responses[6].q[id]=a.gauss-q.gauss;responses[7].q[id]=a.compact-q.compact;responses[8].q[id]=a.nonorth-q.nonorth;responses[9].q[id]=a.relaxation-q.relaxation;
        }
        auto passes=csv(dir/"pressure-pass-response.csv","source,pass,pressure_projection_dimensionless,pressure_norm_ratio_dimensionless,iterations");
        const auto projection=[&](const FlowInitialGuess2D& a,int kind) {
            long double dotProduct=0,size=0;
            if(kind==2) {for(std::size_t i=0;i<nf;++i){const double w=std::hypot(m.faces[i].areaVector.x,m.faces[i].areaVector.y);dotProduct+=w*static_cast<long double>(x.flux[i])*a.flux[i];size+=w*static_cast<long double>(x.flux[i])*x.flux[i];}}
            else for(std::size_t i=0;i<n;++i){const double w=m.cells[i].area;if(kind==0){dotProduct+=w*(static_cast<long double>(x.u[i])*a.u[i]+static_cast<long double>(x.v[i])*a.v[i]);size+=w*(static_cast<long double>(x.u[i])*x.u[i]+static_cast<long double>(x.v[i])*x.v[i]);}else{dotProduct+=w*static_cast<long double>(x.p[i])*a.p[i];size+=w*static_cast<long double>(x.p[i])*x.p[i];}}
            return static_cast<double>(dotProduct/size);
        };
        const auto normRatio=[&](const FlowInitialGuess2D& a,int kind) {
            long double sum=0,size=0;if(kind==2){for(std::size_t i=0;i<nf;++i){const double w=std::hypot(m.faces[i].areaVector.x,m.faces[i].areaVector.y);sum+=w*static_cast<long double>(a.flux[i])*a.flux[i];size+=w*static_cast<long double>(x.flux[i])*x.flux[i];}}
            else for(std::size_t i=0;i<n;++i){const double w=m.cells[i].area;if(kind==0){sum+=w*(static_cast<long double>(a.u[i])*a.u[i]+static_cast<long double>(a.v[i])*a.v[i]);size+=w*(static_cast<long double>(x.u[i])*x.u[i]+static_cast<long double>(x.v[i])*x.v[i]);}else{sum+=w*static_cast<long double>(a.p[i])*a.p[i];size+=w*static_cast<long double>(x.p[i])*x.p[i];}}
            return static_cast<double>(std::sqrt(sum/size));
        };
        for(std::size_t k=0;k<responses.size();++k) {
            auto& r=responses[k];
            if(k==10){r.predictor.p=x.p;r.final.p=x.p;continue;}
            if(k<6) {
                for(int component=0;component<2;++component) {
                    const auto& a=t.momentum[component];auto& matrix=component?av:au;matrix.diag=a.diag;matrix.off=a.off;matrix.rhs=component?r.rhsV:r.rhsU;
                    for(std::size_t i=0;i<n;++i)matrix.diag[i]=std::max(t.momentum[0].diag[i],t.momentum[1].diag[i])/c.velocityRelaxation;
                    matrix.solve(component?r.predictor.v:r.predictor.u,work,.01*c.tolerance*c.speed*c.velocityRelaxation,std::numeric_limits<double>::infinity(),detail::LinearSolveMethod2D::Jacobi,1e-11);
                }
                FlowResult2D state;state.u=r.predictor.u;state.v=r.predictor.v;state.p=cellZeros;state.flux=zeros;Vec df(nf);
                const auto gu=uStencil.apply(state.u,zeros),gv=vStencil.apply(state.v,zeros);
                r.q=rhieChowFlux(m,c,homogeneous,state,t.rc[0].ra,cellZeros,cellZeros,gradientZeros,gradientZeros,gu,gv,gradientZeros,gradientZeros,zeros,false,0,df);
            }
            r.predictor.flux=r.q;Vec pc(n);FlowResult2D scratch;std::size_t pass=0;
            const auto solve=[&](const System& matrix,Vec& field,bool pressure) {
                ensure(pressure,"Expected only pressure correction solve");const auto it=matrix.solvePressure(field,work,detail::LinearPressureMethod2D::IC0,1e-11);
                auto pstate=zero;pstate.p=field;passes<<r.name<<','<<++pass<<','<<projection(pstate,1)<<','<<normRatio(pstate,1)<<','<<it<<'\n';return it;
            };
            const auto lagged=solvePressureCorrection(m,c,b,scratch,ap,pc,t.rc[0].ra,t.rc[0].df,r.q,pressureStencil,zeros,solve);
            const auto gc=detail::conservativePressureGradient(m,detail::pressureFaceValues(m,pc,pressureStencil.apply(pc,zeros),zeros,b.fixedP));
            r.final=r.predictor;
            for(std::size_t i=0;i<n;++i){r.increment.u[i]=-t.rc[0].ra[i]*gc[i].x;r.increment.v[i]=-t.rc[0].ra[i]*gc[i].y;r.increment.p[i]=c.pressureRelaxation*pc[i];r.final.u[i]+=r.increment.u[i];r.final.v[i]+=r.increment.v[i];r.final.p[i]+=r.increment.p[i];}
            for(std::size_t id=0;id<nf;++id){const auto& f=m.faces[id];r.increment.flux[id]=(f.neighbour||b.fixedP[id])?t.rc[0].df[id]*(pc[f.owner]-(f.neighbour?pc[*f.neighbour]:0))+lagged[id]:0;r.final.flux[id]+=r.increment.flux[id];}
        }
        auto projected=csv(dir/"source-projection.csv","source,stage,velocity_projection_dimensionless,velocity_norm_ratio_dimensionless,pressure_projection_dimensionless,pressure_norm_ratio_dimensionless,flux_projection_dimensionless,flux_norm_ratio_dimensionless");
        auto target=csv(dir/"source-target.csv","source,stage,cell0,u_m_s,v_m_s,p_m2_s2");
        auto sum=zero,predSum=zero;Vec rhsUSum(n),rhsVSum(n);
        for(const auto& r:responses) {
            for(const auto& entry:{std::pair{"predictor",&r.predictor},std::pair{"increment",&r.increment},std::pair{"final",&r.final}}) {
                projected<<r.name<<','<<entry.first;for(int kind=0;kind<3;++kind)projected<<','<<projection(*entry.second,kind)<<','<<normRatio(*entry.second,kind);projected<<'\n';
                for(auto cell:{targetCell,targetCell+1})if(cell<n)target<<r.name<<','<<entry.first<<','<<cell<<','<<entry.second->u[cell]<<','<<entry.second->v[cell]<<','<<entry.second->p[cell]<<'\n';
            }
            for(std::size_t i=0;i<n;++i){sum.u[i]+=r.final.u[i];sum.v[i]+=r.final.v[i];sum.p[i]+=r.final.p[i];predSum.u[i]+=r.predictor.u[i];predSum.v[i]+=r.predictor.v[i];rhsUSum[i]+=r.rhsU[i];rhsVSum[i]+=r.rhsV[i];}
            for(std::size_t id=0;id<nf;++id){sum.flux[id]+=r.final.flux[id];predSum.flux[id]+=r.predictor.flux[id];}
        }
        fields(dir,"source-sum",m,sum);fields(dir,"replayed-modal",m,state(actual));fields(dir,"replayed-base",m,state(base));
        auto check=csv(dir/"response-closure.csv","quantity,max_abs_difference,units");double eu=0,ev=0,ep=0,ef=0,epu=0,epv=0,eq=0,eru=0,erv=0;
        for(std::size_t i=0;i<n;++i){eu=std::max(eu,std::abs(sum.u[i]-(actual.u[i]-base.u[i])));ev=std::max(ev,std::abs(sum.v[i]-(actual.v[i]-base.v[i])));ep=std::max(ep,std::abs(sum.p[i]-(actual.p[i]-base.p[i])));epu=std::max(epu,std::abs(predSum.u[i]-(t.rc[0].u[i]-z.rc[0].u[i])));epv=std::max(epv,std::abs(predSum.v[i]-(t.rc[0].v[i]-z.rc[0].v[i])));for(int component=0;component<2;++component){const auto& a=t.momentum[component];const auto& q=z.momentum[component];const double d=std::max(t.momentum[0].diag[i],t.momentum[1].diag[i])/c.velocityRelaxation;const double actualRhs=a.rhs[i]-q.rhs[i]+(d-a.diag[i])*a.field[i]-(d-q.diag[i])*q.field[i];auto& e=component?erv:eru;e=std::max(e,std::abs((component?rhsVSum[i]:rhsUSum[i])-actualRhs));}}
        for(std::size_t id=0;id<nf;++id){ef=std::max(ef,std::abs(sum.flux[id]-(actual.flux[id]-base.flux[id])));eq=std::max(eq,std::abs(predSum.flux[id]-(t.rc[0].predicted[id]-z.rc[0].predicted[id])));}
        check<<"relaxed_rhs_u,"<<eru<<",m3/s2\nrelaxed_rhs_v,"<<erv<<",m3/s2\npredictor_u,"<<epu<<",m/s\npredictor_v,"<<epv<<",m/s\npredicted_phi,"<<eq<<",m2/s\nfinal_u,"<<eu<<",m/s\nfinal_v,"<<ev<<",m/s\nfinal_p,"<<ep<<",m2/s2\nfinal_phi,"<<ef<<",m2/s\n";
        std::cout<<std::setprecision(17)<<"source_sum_velocity_projection="<<projection(sum,0)<<" pressure="<<projection(sum,1)<<" flux="<<projection(sum,2)<<" max_final_errors="<<eu<<','<<ev<<','<<ep<<','<<ef<<'\n';
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
