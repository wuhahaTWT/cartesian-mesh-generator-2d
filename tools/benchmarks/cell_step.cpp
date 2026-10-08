// Read-only GNU ld interceptors. Every actual solve is the production solve.
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include "FlowSolverDetail2D.hpp"
#include "symbols.hpp"
#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::fv::solver_detail;
using Grad=std::vector<Vector2D>;
struct MomentumTrace {
    Vec field,flux,diag,rhs,off;
    Grad gradient,force,source,stress;
    std::vector<std::size_t> rows,columns;
    bool y;
};
struct FluxTrace {
    Vec u,v,p,phi,ra,oldU,oldV,predicted,df;
    Grad gu,gv,gup,gvp,gp,force;
};
struct Trace {std::vector<MomentumTrace> momentum; std::vector<FluxTrace> rc;};
static Trace* active=nullptr;
#ifdef DIAG_WRAP
#define MOM_ARGS System& a,const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,const Vec& field,const Vec& flux,const Grad& gradField,const Grad& gp,const Grad& source,const Grad& stress,const Grad& faceVelocity,bool y,const Vec* previous,double timeStep
void realMomentum(MOM_ARGS) asm("__real_" MOMENTUM_SYMBOL);
void wrapMomentum(MOM_ARGS) asm("__wrap_" MOMENTUM_SYMBOL);
void wrapMomentum(MOM_ARGS) {
    realMomentum(a,m,c,b,field,flux,gradField,gp,source,stress,faceVelocity,y,previous,timeStep);
    if(active) {
        ensure(c.momentumInertia==0 && !previous,"Trace only supports steady production Stokes step");
        active->momentum.push_back({field,flux,a.diag,a.rhs,a.off,gradField,gp,source,stress,a.pattern.rows,a.pattern.columns,y});
    }
}
#define RC_ARGS const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,const FlowResult2D& r,const Vec& ra,const Vec& oldU,const Vec& oldV,const Grad& gu,const Grad& gv,const Grad& gup,const Grad& gvp,const Grad& gp,const Grad& force,const Vec& oldDefect,bool previous,double timeStep,Vec& df
Vec realRC(RC_ARGS) asm("__real_" RHIECHOWFLUX_SYMBOL);
Vec wrapRC(RC_ARGS) asm("__wrap_" RHIECHOWFLUX_SYMBOL);
Vec wrapRC(RC_ARGS) {
    auto predicted=realRC(m,c,b,r,ra,oldU,oldV,gu,gv,gup,gvp,gp,force,oldDefect,previous,timeStep,df);
    if(active) {
        ensure(!previous,"Trace only supports steady RC");
        active->rc.push_back({r.u,r.v,r.p,r.flux,ra,oldU,oldV,predicted,df,gu,gv,gup,gvp,gp,force});
    }
    return predicted;
}
#endif
static std::ofstream csv(const std::filesystem::path& p,const char* header) {
    std::ofstream o(p);ensure(bool(o),"Cannot write diagnostic file");o<<std::setprecision(17)<<header<<'\n';return o;
}
static std::uint64_t hash(const Vec& v) {
    std::uint64_t h=14695981039346656037ULL;
    for(double x:v) {auto u=std::bit_cast<std::uint64_t>(x);for(int k=0;k<8;++k){h^=(u>>(8*k))&255;h*=1099511628211ULL;}}
    return h;
}
static void hashes(std::ostream& o,const FlowResult2D& r) {
    o<<hash(r.u)<<','<<hash(r.v)<<','<<hash(r.p)<<','<<hash(r.flux);
}
static void fields(const std::filesystem::path& dir,const char* name,const FvMesh2D& m,const FlowInitialGuess2D& x) {
    auto o=csv(dir/(std::string(name)+"-cells.csv"),"cell0,cell1,x,y,area,u_m_s,v_m_s,p_m2_s2");
    for(std::size_t i=0;i<m.cells.size();++i)o<<i<<','<<i+1<<','<<m.cells[i].centre.x<<','<<m.cells[i].centre.y<<','<<m.cells[i].area<<','<<x.u[i]<<','<<x.v[i]<<','<<x.p[i]<<'\n';
    auto f=csv(dir/(std::string(name)+"-faces.csv"),"face0,owner0,neighbour0,phi_m2_s");
    for(std::size_t i=0;i<m.faces.size();++i)f<<i<<','<<m.faces[i].owner<<','<<(m.faces[i].neighbour?std::to_string(*m.faces[i].neighbour):"-1")<<','<<x.flux[i]<<'\n';
}
static FlowInitialGuess2D state(const FlowResult2D& r) {return {r.u,r.v,r.p,r.flux};}
struct FaceTerms {double pressure=0,wall=0,nonorth=0,stress=0,diagonal=0,offdiag=0;};
static FaceTerms faceTerms(const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,const MomentumTrace& t,const Vec& pf,std::size_t face,std::size_t cell) {
    const auto& f=m.faces[face];const double sign=f.owner==cell?1:-1;const auto& bc=t.y?b.v:b.u;
    const auto g=f.neighbour?interpolateGradient(f,t.gradient):t.gradient[f.owner];
    const double nu=faceNu(c,face),d=nu*f.transmissibility;
    FaceTerms a;
    a.pressure=-sign*pf[face]*(t.y?f.areaVector.y:f.areaVector.x);
    a.nonorth=sign*nu*dot(g,f.correction);
    a.stress=t.stress.empty()?0:-sign*(t.y?t.stress[face].y:t.stress[face].x);
    a.diagonal=d;
    if(f.neighbour)a.offdiag=-d;
    else a.wall=d*bc[face];
    return a;
}
struct RcTerms {double interpolation=0,skew=0,gauss=0,compact=0,nonorth=0,relaxation=0,actual=0,df=0;};
static RcTerms rcTerms(const FvMesh2D& m,const FlowControls2D& c,const FluxTrace& t,std::size_t id) {
    const auto& f=m.faces[id];RcTerms a;a.actual=t.predicted[id];a.df=t.df[id];if(!f.neighbour)return a;
    const auto i=f.owner,j=*f.neighbour;const double w=f.neighbourWeight,rf=interpolate(f,t.ra);
    const Point2D point{m.cells[i].centre.x*(1-w)+m.cells[j].centre.x*w,m.cells[i].centre.y*(1-w)+m.cells[j].centre.y*w};
    const auto skew=f.centre-point;
    a.interpolation=interpolate(f,t.u)*f.areaVector.x+interpolate(f,t.v)*f.areaVector.y;
    a.skew=dot(interpolateGradient(f,t.gup),skew)*f.areaVector.x+dot(interpolateGradient(f,t.gvp),skew)*f.areaVector.y;
    const Vector2D rag{(1-w)*t.ra[i]*t.force[i].x+w*t.ra[j]*t.force[j].x,(1-w)*t.ra[i]*t.force[i].y+w*t.ra[j]*t.force[j].y};
    a.gauss=dot(rag,f.areaVector);
    a.compact=-rf*f.transmissibility*(t.p[j]-t.p[i]);
    a.nonorth=-rf*dot(interpolateGradient(f,t.gp),f.correction);
    const double oldUf=interpolate(f,t.oldU)+dot(interpolateGradient(f,t.gu),skew),oldVf=interpolate(f,t.oldV)+dot(interpolateGradient(f,t.gv),skew);
    a.relaxation=(1-c.velocityRelaxation)*(t.phi[id]-oldUf*f.areaVector.x-oldVf*f.areaVector.y);
    return a;
}
static void exportTrace(const std::filesystem::path& dir,const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,const Trace& t,const Trace& base,const FlowInitialGuess2D& input,const FlowResult2D& output,const FlowResult2D& baseOutput,std::size_t target) {
    ensure(t.momentum.size()==4 && base.momentum.size()==4 && t.rc.size()==1 && base.rc.size()==1,"Expected exactly pre/post U/V plus one RC");
    const auto stencil=detail::buildFlowGradientStencil2D(m,b.fixedP,true);
    std::set<std::size_t> selected{target};if(target+1<m.cells.size())selected.insert(target+1);
    for(const auto& f:m.faces)if(f.neighbour && (f.owner==target || *f.neighbour==target)){selected.insert(f.owner);selected.insert(*f.neighbour);}
    std::vector<std::vector<std::size_t>> incident(m.cells.size());
    auto geometry=csv(dir/"face-geometry.csv","face0,owner0,neighbour0,x,y,Sx_m,Sy_m,correction_x_m,correction_y_m,transmissibility,neighbourWeight");
    for(std::size_t id=0;id<m.faces.size();++id){const auto& f=m.faces[id];incident[f.owner].push_back(id);if(f.neighbour)incident[*f.neighbour].push_back(id);geometry<<id<<','<<f.owner<<','<<(f.neighbour?std::to_string(*f.neighbour):"-1")<<','<<f.centre.x<<','<<f.centre.y<<','<<f.areaVector.x<<','<<f.areaVector.y<<','<<f.correction.x<<','<<f.correction.y<<','<<f.transmissibility<<','<<f.neighbourWeight<<'\n';}
    auto fo=csv(dir/"momentum-face-terms.csv","run,phase,component,cell0,cell1,face0,owner0,neighbour0,kind,orientation,pressure_rhs_m3_s2,wall_dirichlet_rhs_m3_s2,viscous_nonorth_rhs_m3_s2,symmetric_stress_rhs_m3_s2,compact_diagonal_m2_s,compact_offdiagonal_m2_s,compact_oldfield_action_m3_s2,complete_viscous_oldfield_action_m3_s2");
    auto oldAction=csv(dir/"momentum-old-field-action.csv","run,phase,component,cell0,wall_compact_m3_s2,wall_nonorth_m3_s2,wall_stress_m3_s2,wall_complete_m3_s2,internal_compact_m3_s2,internal_nonorth_m3_s2,internal_stress_m3_s2,internal_complete_m3_s2,pressure_m3_s2,source_m3_s2,sum_unrelaxed_defect_m3_s2,actual_rhs_minus_matrix_oldfield_m3_s2,alignment_m3_s2");
    auto ro=csv(dir/"momentum-rows.csv","run,phase,component,cell0,cell1,pressure_face_sum_m3_s2,pressure_actual_m3_s2,wall_rhs_m3_s2,nonorth_rhs_m3_s2,stress_rhs_m3_s2,source_rhs_m3_s2,reconstructed_rhs_m3_s2,actual_rhs_m3_s2,rhs_difference_m3_s2,compact_diagonal_m2_s,actual_diagonal_m2_s,compact_diagonal_difference_m2_s,actual_offdiagonal_sum_m2_s,reconstructed_offdiagonal_sum_m2_s,offdiag_difference_m2_s,relaxation_rhs_m3_s2,common_diagonal_rhs_m3_s2,relaxed_rhs_m3_s2,relaxed_diagonal_m2_s,diagonal_times_evaluated_u_m3_s2,offdiag_times_evaluated_u_m3_s2,row_residual_m3_s2");
    auto raw=csv(dir/"momentum-actual-cells.csv","run,call,phase,component,cell0,field_m_s,gradient_x_s_1,gradient_y_s_1,force_x_m_s2,force_y_m_s2,diagonal_m2_s,rhs_m3_s2");
    auto off=csv(dir/"momentum-actual-offdiagonal.csv","run,call,row0,column0,coefficient_m2_s");
    double maxRhsError=0,maxForceError=0,maxRcError=0,maxDiagonalError=0,maxOffError=0;
    for(int run=0;run<3;++run)for(std::size_t call=0;call<4;++call) {
        const auto& q=t.momentum[call];const auto& z=base.momentum[call];const bool delta=run==2;const auto& a=run==1?z:q;
        const bool post=call>=2;const auto& p=post?output.p:input.p;Vec pbase=post?baseOutput.p:Vec(m.cells.size());
        const auto pf=detail::pressureFaceValues(m,p,stencil.apply(p,b.p),b.p,b.fixedP);
        const auto pz=detail::pressureFaceValues(m,pbase,stencil.apply(pbase,b.p),b.p,b.fixedP);
        const auto& pa=run==1?pz:pf;const char* name=run==0?"modal":run==1?"base":"delta";
        const auto difference=[&](double x,double y){return delta?x-y:x;};
        for(std::size_t i=0;i<m.cells.size();++i) {
            raw<<name<<','<<call<<','<<(post?"post":"pre")<<','<<(a.y?'v':'u')<<','<<i<<','<<difference(a.field[i],z.field[i])<<','<<difference(a.gradient[i].x,z.gradient[i].x)<<','<<difference(a.gradient[i].y,z.gradient[i].y)<<','<<difference(a.force[i].x,z.force[i].x)<<','<<difference(a.force[i].y,z.force[i].y)<<','<<difference(a.diag[i],z.diag[i])<<','<<difference(a.rhs[i],z.rhs[i])<<'\n';
            for(auto k=a.rows[i];k<a.rows[i+1];++k)off<<name<<','<<call<<','<<i<<','<<a.columns[k]<<','<<difference(a.off[k],z.off[k])<<'\n';
        }
        for(std::size_t cell=0;cell<m.cells.size();++cell) {
            FaceTerms sums;double offsum=0,offEval=0,wallCompact=0,wallNonorth=0,wallStress=0,internalCompact=0,internalNonorth=0,internalStress=0,oldMatrixAction=0;
            const Vec& evaluate=post?(a.y?output.v:output.u):(a.y?t.rc[0].v:t.rc[0].u);
            const Vec& evalbase=post?(a.y?baseOutput.v:baseOutput.u):(a.y?base.rc[0].v:base.rc[0].u);
            const Vec& eval=run==1?evalbase:evaluate;
            for(const auto id:incident[cell]) {
                const auto& f=m.faces[id];
                auto terms=faceTerms(m,c,b,a,pa,id,cell);const auto zero=faceTerms(m,c,b,z,pz,id,cell);
                const auto& bc=a.y?b.v:b.u;const double d=faceNu(c,id)*f.transmissibility;
                const auto other=f.neighbour?(f.owner==cell?*f.neighbour:f.owner):cell;
                const double compact=d*((f.neighbour?a.field[other]:bc[id])-a.field[cell]);
                const double compactBase=d*((f.neighbour?z.field[other]:bc[id])-z.field[cell]);
                const double compactAction=difference(compact,compactBase);
                if(delta){terms.pressure-=zero.pressure;terms.wall-=zero.wall;terms.nonorth-=zero.nonorth;terms.stress-=zero.stress;terms.diagonal-=zero.diagonal;terms.offdiag-=zero.offdiag;}
                sums.pressure+=terms.pressure;sums.wall+=terms.wall;sums.nonorth+=terms.nonorth;sums.stress+=terms.stress;sums.diagonal+=terms.diagonal;sums.offdiag+=terms.offdiag;
                if(f.neighbour){internalCompact+=compactAction;internalNonorth+=terms.nonorth;internalStress+=terms.stress;}else{wallCompact+=compactAction;wallNonorth+=terms.nonorth;wallStress+=terms.stress;}
                if(selected.contains(cell))fo<<name<<','<<(post?"post":"pre")<<','<<(a.y?'v':'u')<<','<<cell<<','<<cell+1<<','<<id<<','<<f.owner<<','<<(f.neighbour?std::to_string(*f.neighbour):"-1")<<','<<(f.neighbour?"internal":"wall")<<','<<(f.owner==cell?1:-1)<<','<<terms.pressure<<','<<terms.wall<<','<<terms.nonorth<<','<<terms.stress<<','<<terms.diagonal<<','<<terms.offdiag<<','<<compactAction<<','<<compactAction+terms.nonorth+terms.stress<<'\n';
            }
            for(auto k=a.rows[cell];k<a.rows[cell+1];++k){offsum+=difference(a.off[k],z.off[k]);offEval+=difference(a.off[k]*eval[a.columns[k]],z.off[k]*evalbase[a.columns[k]]);}
            const double source=a.source.empty()?0:(a.y?a.source[cell].y:a.source[cell].x),sourcebase=z.source.empty()?0:(z.y?z.source[cell].y:z.source[cell].x);
            const double force=-m.cells[cell].area*difference(a.y?a.force[cell].y:a.force[cell].x,z.y?z.force[cell].y:z.force[cell].x);
            const double reconstructed=force+sums.wall+sums.nonorth+sums.stress+difference(source,sourcebase),actual=difference(a.rhs[cell],z.rhs[cell]);
            oldMatrixAction=difference(a.diag[cell]*a.field[cell],z.diag[cell]*z.field[cell]);for(auto k=a.rows[cell];k<a.rows[cell+1];++k)oldMatrixAction+=difference(a.off[k]*a.field[a.columns[k]],z.off[k]*z.field[a.columns[k]]);
            const double completeWall=wallCompact+wallNonorth+wallStress,completeInternal=internalCompact+internalNonorth+internalStress,defect=completeWall+completeInternal+force+difference(source,sourcebase);
            oldAction<<name<<','<<(post?"post":"pre")<<','<<(a.y?'v':'u')<<','<<cell<<','<<wallCompact<<','<<wallNonorth<<','<<wallStress<<','<<completeWall<<','<<internalCompact<<','<<internalNonorth<<','<<internalStress<<','<<completeInternal<<','<<force<<','<<difference(source,sourcebase)<<','<<defect<<','<<actual-oldMatrixAction<<','<<defect-(actual-oldMatrixAction)<<'\n';
            const double common=std::max(t.momentum[0].diag[cell],t.momentum[1].diag[cell])/c.velocityRelaxation,commonbase=std::max(base.momentum[0].diag[cell],base.momentum[1].diag[cell])/c.velocityRelaxation;
            const double relax=post?0:difference((a.diag[cell]/c.velocityRelaxation-a.diag[cell])*a.field[cell],(z.diag[cell]/c.velocityRelaxation-z.diag[cell])*z.field[cell]);
            const double commonRhs=post?0:difference(((run==1?commonbase:common)-a.diag[cell]/c.velocityRelaxation)*a.field[cell],(commonbase-z.diag[cell]/c.velocityRelaxation)*z.field[cell]);
            const double diagonal=post?difference(a.diag[cell],z.diag[cell]):difference(run==1?commonbase:common,commonbase);
            const double diagEval=post?difference(a.diag[cell]*eval[cell],z.diag[cell]*evalbase[cell]):difference((run==1?commonbase:common)*eval[cell],commonbase*evalbase[cell]);
            ro<<name<<','<<(post?"post":"pre")<<','<<(a.y?'v':'u')<<','<<cell<<','<<cell+1<<','<<sums.pressure<<','<<force<<','<<sums.wall<<','<<sums.nonorth<<','<<sums.stress<<','<<difference(source,sourcebase)<<','<<reconstructed<<','<<actual<<','<<reconstructed-actual<<','<<sums.diagonal<<','<<difference(a.diag[cell],z.diag[cell])<<','<<sums.diagonal-difference(a.diag[cell],z.diag[cell])<<','<<offsum<<','<<sums.offdiag<<','<<offsum-sums.offdiag<<','<<relax<<','<<commonRhs<<','<<actual+relax+commonRhs<<','<<diagonal<<','<<diagEval<<','<<offEval<<','<<actual+relax+commonRhs-diagEval-offEval<<'\n';
            maxRhsError=std::max(maxRhsError,std::abs(reconstructed-actual));maxForceError=std::max(maxForceError,std::abs(sums.pressure-force));maxDiagonalError=std::max(maxDiagonalError,std::abs(sums.diagonal-difference(a.diag[cell],z.diag[cell])));maxOffError=std::max(maxOffError,std::abs(offsum-sums.offdiag));
        }
    }
    auto rc=csv(dir/"rc-face-terms.csv","run,face0,owner0,neighbour0,kind,interpolation_m2_s,predictor_skew_m2_s,rag_gauss_m2_s,compact_pressure_m2_s,nonorth_pressure_m2_s,relaxation_defect_m2_s,sum_m2_s,actual_predicted_m2_s,sum_difference_m2_s,df_s");
    auto rcraw=csv(dir/"rc-actual-cells.csv","run,cell0,predictor_u_m_s,predictor_v_m_s,p_m2_s2,ra_s,old_u_m_s,old_v_m_s,gu_x_s_1,gu_y_s_1,gv_x_s_1,gv_y_s_1,gup_x_s_1,gup_y_s_1,gvp_x_s_1,gvp_y_s_1,gp_x_m_s2,gp_y_m_s2,force_x_m_s2,force_y_m_s2");
    auto pc=csv(dir/"pressure-correction.csv","cell0,pc_modal_m2_s2,pc_base_m2_s2,pc_delta_m2_s2,pressure_relaxation,predictor_delta_u_m_s,predictor_delta_v_m_s,corrected_delta_u_m_s,corrected_delta_v_m_s,correction_delta_u_m_s,correction_delta_v_m_s");
    for(std::size_t i=0;i<m.cells.size();++i) {
        const auto& q=t.rc[0];const auto& z=base.rc[0];
        pc<<i<<','<<(output.p[i]-q.p[i])/c.pressureRelaxation<<','<<(baseOutput.p[i]-z.p[i])/c.pressureRelaxation<<','<<((output.p[i]-q.p[i])-(baseOutput.p[i]-z.p[i]))/c.pressureRelaxation<<','<<c.pressureRelaxation<<','<<q.u[i]-z.u[i]<<','<<q.v[i]-z.v[i]<<','<<output.u[i]-baseOutput.u[i]<<','<<output.v[i]-baseOutput.v[i]<<','<<(output.u[i]-q.u[i])-(baseOutput.u[i]-z.u[i])<<','<<(output.v[i]-q.v[i])-(baseOutput.v[i]-z.v[i])<<'\n';
    }
    for(int run=0;run<3;++run) {
        const auto& z=base.rc[0];const auto& q=run==1?z:t.rc[0];const bool delta=run==2;
        const char* name=run==0?"modal":run==1?"base":"delta";const auto difference=[&](double x,double y){return delta?x-y:x;};
        for(std::size_t i=0;i<m.cells.size();++i) {
            rcraw<<name<<','<<i<<','<<difference(q.u[i],z.u[i])<<','<<difference(q.v[i],z.v[i])<<','<<difference(q.p[i],z.p[i])<<','<<difference(q.ra[i],z.ra[i])<<','<<difference(q.oldU[i],z.oldU[i])<<','<<difference(q.oldV[i],z.oldV[i]);
            for(const auto pair:{std::pair{&q.gu,&z.gu},std::pair{&q.gv,&z.gv},std::pair{&q.gup,&z.gup},std::pair{&q.gvp,&z.gvp},std::pair{&q.gp,&z.gp},std::pair{&q.force,&z.force}})rcraw<<','<<difference((*pair.first)[i].x,(*pair.second)[i].x)<<','<<difference((*pair.first)[i].y,(*pair.second)[i].y);
            rcraw<<'\n';
        }
        for(std::size_t id=0;id<m.faces.size();++id) {
            const auto& f=m.faces[id];auto a=rcTerms(m,c,q,id);const auto zero=rcTerms(m,c,z,id);
            if(delta){a.interpolation-=zero.interpolation;a.skew-=zero.skew;a.gauss-=zero.gauss;a.compact-=zero.compact;a.nonorth-=zero.nonorth;a.relaxation-=zero.relaxation;a.actual-=zero.actual;a.df-=zero.df;}
            const double sum=a.interpolation+a.skew+a.gauss+a.compact+a.nonorth+a.relaxation;
            maxRcError=std::max(maxRcError,std::abs(sum-a.actual));
            rc<<name<<','<<id<<','<<f.owner<<','<<(f.neighbour?std::to_string(*f.neighbour):"-1")<<','<<(f.neighbour?"internal":"wall")<<','<<a.interpolation<<','<<a.skew<<','<<a.gauss<<','<<a.compact<<','<<a.nonorth<<','<<a.relaxation<<','<<sum<<','<<a.actual<<','<<sum-a.actual<<','<<a.df<<'\n';
        }
    }
    auto audit=csv(dir/"alignment.csv","quantity,max_abs_difference,units");
    audit<<"momentum_rhs,"<<maxRhsError<<",m3/s2\npressure_face_vs_actual_force,"<<maxForceError<<",m3/s2\ncompact_diagonal,"<<maxDiagonalError<<",m2/s\ncompact_offdiag_sum,"<<maxOffError<<",m2/s\nrc_predicted,"<<maxRcError<<",m2/s\n";
    auto stages=csv(dir/"stage-projection.csv","stage,velocity_projection_dimensionless,velocity_norm_ratio_dimensionless,pressure_projection_dimensionless,pressure_norm_ratio_dimensionless,flux_projection_dimensionless,flux_norm_ratio_dimensionless");
    auto points=csv(dir/"stage-target-response.csv","stage,cell0,cell1,u_delta_m_s,v_delta_m_s,p_delta_m2_s2");
    FlowInitialGuess2D predictor=input,final=state(output),correction=input;
    for(std::size_t i=0;i<m.cells.size();++i){predictor.u[i]=t.rc[0].u[i]-base.rc[0].u[i];predictor.v[i]=t.rc[0].v[i]-base.rc[0].v[i];predictor.p[i]=t.rc[0].p[i]-base.rc[0].p[i];final.u[i]-=baseOutput.u[i];final.v[i]-=baseOutput.v[i];final.p[i]-=baseOutput.p[i];correction.u[i]=final.u[i]-predictor.u[i];correction.v[i]=final.v[i]-predictor.v[i];correction.p[i]=final.p[i]-predictor.p[i];}
    for(std::size_t id=0;id<m.faces.size();++id){predictor.flux[id]=t.rc[0].predicted[id]-base.rc[0].predicted[id];final.flux[id]-=baseOutput.flux[id];correction.flux[id]=final.flux[id]-predictor.flux[id];}
    const FlowInitialGuess2D* values[]={&input,&predictor,&correction,&final};const char* labels[]={"input","momentum_predictor_and_rc","pressure_correction_increment","final"};
    for(int stage=0;stage<4;++stage){const auto& a=*values[stage];long double vv=0,vi=0,vn=0,pp=0,pi=0,pn=0,ff=0,fi=0,fn=0;for(std::size_t i=0;i<m.cells.size();++i){const auto w=m.cells[i].area;vv+=w*(input.u[i]*input.u[i]+input.v[i]*input.v[i]);vi+=w*(input.u[i]*a.u[i]+input.v[i]*a.v[i]);vn+=w*(a.u[i]*a.u[i]+a.v[i]*a.v[i]);pp+=w*input.p[i]*input.p[i];pi+=w*input.p[i]*a.p[i];pn+=w*a.p[i]*a.p[i];}for(std::size_t id=0;id<m.faces.size();++id){const auto w=std::hypot(m.faces[id].areaVector.x,m.faces[id].areaVector.y);ff+=w*input.flux[id]*input.flux[id];fi+=w*input.flux[id]*a.flux[id];fn+=w*a.flux[id]*a.flux[id];}stages<<labels[stage]<<','<<static_cast<double>(vi/vv)<<','<<static_cast<double>(std::sqrt(vn/vv))<<','<<static_cast<double>(pi/pp)<<','<<static_cast<double>(std::sqrt(pn/pp))<<','<<static_cast<double>(fi/ff)<<','<<static_cast<double>(std::sqrt(fn/ff))<<'\n';for(const auto i:selected)points<<labels[stage]<<','<<i<<','<<i+1<<','<<a.u[i]<<','<<a.v[i]<<','<<a.p[i]<<'\n';}
}
int main(int argc,char** argv) {
    try {
        ensure(argc>=4 && argc<=6,"Usage: cell_step MESH BOUNDARIES OUTPUT_DIR [target0=3710] [power_steps=180]");
        const auto target=argc>4?std::stoul(argv[4]):3710UL,steps=argc>5?std::stoul(argv[5]):180UL;
        const std::filesystem::path dir=argv[3];std::filesystem::create_directories(dir);
        const auto input=readCm2dTopology(argv[1]);ensure(input.valid(),"Invalid mesh");const auto m=makeFvMesh2D(input.topology);ensure(target<m.cells.size(),"Target0 outside mesh");
        FlowControls2D c;c.scenario="custom";c.nu=.1;c.speed=.5;c.tolerance=1e-8;c.maxIterations=1;c.momentumInertia=0;c.pressureCorrectionPasses=4;
        c.pressurePreconditioner=detail::systemCholeskyAvailable2D()?PressurePreconditioner2D::SystemCholesky:PressurePreconditioner2D::IncompleteCholesky0;
        c.convection=ConvectionScheme2D::FaceLimitedLinearUpwind;c.viscousStress=ViscousStress2D::Symmetric;
        std::ifstream boundariesFile(argv[2]);ensure(bool(boundariesFile),"Cannot read boundaries");c.boundaryConditions=readFlowBoundaryConditions2D(boundariesFile,m,c);const auto b=boundaries(m,c);
        ensure(b.closed,"Requires closed annulus");for(std::size_t id=0;id<m.faces.size();++id)ensure(m.faces[id].neighbour || (b.fixedU[id]&&b.fixedV[id]),"Requires fixed velocity walls");
        FlowInitialGuess2D zero;zero.u.resize(m.cells.size());zero.v.resize(m.cells.size());zero.p.resize(m.cells.size());zero.flux.resize(m.faces.size());
        const auto base=solveIncompressibleFromGuess2D(m,c,zero);auto x=zero;
        for(std::size_t i=0;i<m.cells.size();++i){x.u[i]=std::sin(1.71*static_cast<double>(i));x.v[i]=std::cos(1.53*static_cast<double>(i));x.p[i]=std::sin(1.31*static_cast<double>(i));}x.p[0]=0;
        double area=0;for(const auto& cell:m.cells)area+=cell.area;
        const auto product=[&](const auto& a,const auto& other){long double sum=0;for(std::size_t i=0;i<m.cells.size();++i)sum+=m.cells[i].area/area*(a.u[i]*other.u[i]+a.v[i]*other.v[i]+a.p[i]*other.p[i]);return static_cast<double>(sum);};
        auto power=csv(dir/"power.csv","iteration,normRatio,rayleigh,cellEigenResidual,worstCell0,x,y");double lastRayleigh=0,lastResidual=0;
        for(std::size_t it=1;it<=steps;++it){const auto r=solveIncompressibleFromGuess2D(m,c,x);auto next=zero;for(std::size_t i=0;i<m.cells.size();++i){next.u[i]=r.u[i]-base.u[i];next.v[i]=r.v[i]-base.v[i];next.p[i]=r.p[i]-base.p[i];}for(std::size_t i=0;i<m.faces.size();++i)next.flux[i]=r.flux[i]-base.flux[i];const double size=std::sqrt(product(next,next));ensure(size>0&&std::isfinite(size),"Degenerate power iterate");lastRayleigh=product(x,next)/product(x,x);auto residual=next;std::size_t worst=0;for(std::size_t i=0;i<m.cells.size();++i){residual.u[i]-=lastRayleigh*x.u[i];residual.v[i]-=lastRayleigh*x.v[i];residual.p[i]-=lastRayleigh*x.p[i];if(std::hypot(next.u[i],next.v[i])>std::hypot(next.u[worst],next.v[worst]))worst=i;}lastResidual=std::sqrt(product(residual,residual))/size;power<<it<<','<<size/std::sqrt(product(x,x))<<','<<lastRayleigh<<','<<lastResidual<<','<<worst<<','<<m.cells[worst].centre.x<<','<<m.cells[worst].centre.y<<'\n';for(auto* f:{&next.u,&next.v,&next.p,&next.flux})for(auto& v:*f)v/=size;x=std::move(next);}
        // Trace is disabled throughout modal reconstruction; only F(x), F(0) below capture.
        Trace modalTrace,baseTrace;active=&modalTrace;const auto result=solveIncompressibleFromGuess2D(m,c,x);active=&baseTrace;const auto tracedBase=solveIncompressibleFromGuess2D(m,c,zero);active=nullptr;
        // Same input with recording disabled: wrappers call precisely the same real functions.
        const auto noTrace=solveIncompressibleFromGuess2D(m,c,x);
        ensure(result.u==noTrace.u&&result.v==noTrace.v&&result.p==noTrace.p&&result.flux==noTrace.flux,"Trace changed production output");
        auto delta=state(result);for(std::size_t i=0;i<m.cells.size();++i){delta.u[i]-=tracedBase.u[i];delta.v[i]-=tracedBase.v[i];delta.p[i]-=tracedBase.p[i];}for(std::size_t i=0;i<m.faces.size();++i)delta.flux[i]-=tracedBase.flux[i];
        fields(dir,"input-modal",m,x);fields(dir,"output-modal",m,state(result));fields(dir,"output-base",m,state(tracedBase));fields(dir,"output-delta",m,delta);
#ifdef DIAG_WRAP
        exportTrace(dir,m,c,b,modalTrace,baseTrace,x,result,tracedBase,target);
#endif
        auto summary=csv(dir/"summary.csv","cells,faces,target0,target1,powerSteps,alphaU,alphaP,rayleigh180,residual180,momentumCalls,rcCalls,output_u_hash,output_v_hash,output_p_hash,output_phi_hash");
        summary<<m.cells.size()<<','<<m.faces.size()<<','<<target<<','<<target+1<<','<<steps<<','<<c.velocityRelaxation<<','<<c.pressureRelaxation<<','<<lastRayleigh<<','<<lastResidual<<','<<modalTrace.momentum.size()<<','<<modalTrace.rc.size()<<',';hashes(summary,result);summary<<'\n';
        auto verify=csv(dir/"trace-invariance.csv","field,trace_enabled_hash,trace_disabled_hash,max_abs_difference,bitwise_equal");
        const Vec* aa[]={&result.u,&result.v,&result.p,&result.flux};const Vec* bb[]={&noTrace.u,&noTrace.v,&noTrace.p,&noTrace.flux};const char* names[]={"u","v","p","phi"};
        for(int k=0;k<4;++k)verify<<names[k]<<','<<hash(*aa[k])<<','<<hash(*bb[k])<<",0,"<<(hash(*aa[k])==hash(*bb[k])?"true":"false")<<'\n';
        std::cout<<std::setprecision(17)<<"rayleigh="<<lastRayleigh<<" residual="<<lastResidual<<" target0="<<target<<" target1="<<target+1<<" calls="<<modalTrace.momentum.size()<<'/'<<modalTrace.rc.size()<<" outputHashes=";hashes(std::cout,result);std::cout<<'\n';
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
