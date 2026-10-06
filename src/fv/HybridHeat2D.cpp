#include "cartmesh2d/fv/HybridHeat2D.hpp"
#include "cartmesh2d/fv/HybridCellGeometry2D.hpp"
#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"
#include <numeric>
#include <limits>
#include <stdexcept>
namespace cartmesh2d::fv {
namespace {
void check(bool v,const char* s){if(!v)throw std::runtime_error(std::string("hybrid heat: ")+s);}
constexpr std::size_t absent=std::numeric_limits<std::size_t>::max();
}
HybridHeatOperator2D::HybridHeatOperator2D(const FvMesh2D& mesh,const std::vector<HeatBoundary2D>& bc,double k) {
    validateFvMesh2D(mesh);check(std::isfinite(k)&&k>0,"invalid conductivity");
    boundary_.resize(mesh.faces.size());unknown_.assign(mesh.faces.size(),absent);
    std::vector<bool> assigned(mesh.faces.size());
    for(const auto& x:bc){check(x.face<assigned.size()&&!assigned[x.face]&&!mesh.faces[x.face].neighbour,"invalid/duplicate boundary");
        check(!x.partner&&x.kind!=HeatBoundaryKind2D::Periodic,"periodic traces not implemented");
        check(std::isfinite(x.value),"invalid boundary value");assigned[x.face]=true;
        auto& b=boundary_[x.face];b.fixed=x.kind==HeatBoundaryKind2D::Temperature;b.temperature=x.value;
        if(b.fixed)check(x.value>0,"nonpositive prescribed Kelvin");
        if(x.kind==HeatBoundaryKind2D::OutwardFlux)b.flux=x.value*std::hypot(mesh.faces[x.face].areaVector.x,mesh.faces[x.face].areaVector.y);
        else check(b.fixed||x.kind==HeatBoundaryKind2D::Insulated,"unsupported boundary");
        if(x.kind==HeatBoundaryKind2D::Insulated)check(x.value==0,"inactive insulated value");
    }
    for(std::size_t f=0;f<mesh.faces.size();++f){check(mesh.faces[f].neighbour||assigned[f],"missing boundary");
        owner_.push_back(mesh.faces[f].owner);neighbour_.push_back(mesh.faces[f].neighbour);
        if(!boundary_[f].fixed)unknown_[f]=count_++;
    }
    diagnostics_.minimumNormalDistance=std::numeric_limits<double>::infinity();
    for(std::size_t cell=0;cell<mesh.cells.size();++cell){const auto& c=mesh.cells[cell];Local l;l.faces=c.faces;l.area=c.area;
        const auto geometry=hybridCellGeometry2D(mesh,cell);const auto n=l.faces.size();
        const auto& d=geometry.offsets;const auto& g=geometry.gradientColumns;const auto& r=geometry.remainder;
        auto w=geometry.stabilizationWeights;for(auto& x:w)x*=k;
        diagnostics_.minimumNormalDistance=std::min(diagnostics_.minimumNormalDistance,geometry.minimumNormalDistance);
        diagnostics_.closureRoundoffRatio=std::max(diagnostics_.closureRoundoffRatio,geometry.closureRoundoffRatio);
        diagnostics_.momentRoundoffRatio=std::max(diagnostics_.momentRoundoffRatio,geometry.momentRoundoffRatio);
        l.a.resize(n*n);l.b.resize(n);
        for(std::size_t i=0;i<n;++i)for(std::size_t j=i;j<n;++j){double a=k*c.area*dot(g[i],g[j]);for(std::size_t q=0;q<n;++q)a+=w[q]*r[q*n+i]*r[q*n+j];check(std::isfinite(a),"nonfinite local matrix");l.a[i*n+j]=l.a[j*n+i]=a;}
        for(std::size_t i=0;i<n;++i){for(std::size_t j=0;j<n;++j)l.b[i]+=l.a[i*n+j];l.sum+=l.b[i];
            for(unsigned a=0;a<2;++a){double affine=0;for(std::size_t j=0;j<n;++j)affine+=r[i*n+j]*(a?d[j].y:d[j].x);
                double h=0;for(const auto& x:d)h=std::max(h,std::hypot(x.x,x.y));diagnostics_.maximumLocalAffineResidual=std::max(diagnostics_.maximumLocalAffineResidual,std::abs(affine)/h);}}
        check(std::isfinite(l.sum)&&l.sum>0,"nonpositive local cell pivot");local_.push_back(std::move(l));
    }
}
HybridHeatResult2D HybridHeatOperator2D::solveSteady(double tol) const {
    std::vector<bool> seen(local_.size());
    for(std::size_t seed=0;seed<local_.size();++seed)if(!seen[seed]) {
        bool referenced=false;std::vector<std::size_t> pending{seed};seen[seed]=true;
        while(!pending.empty()) {
            const auto cell=pending.back();pending.pop_back();
            for(auto f:local_[cell].faces) {
                referenced=referenced||boundary_[f].fixed;
                if(neighbour_[f]) {const auto other=owner_[f]==cell?*neighbour_[f]:owner_[f];if(!seen[other]){seen[other]=true;pending.push_back(other);}}
            }
        }
        check(referenced,"steady component without Dirichlet reference; prescribed mean not implemented");
    }
    return solve(std::vector<double>(local_.size()),std::vector<double>(local_.size()),tol);
}
HybridHeatResult2D HybridHeatOperator2D::backwardEuler(const std::vector<double>& t,const std::vector<double>& cv,double dt,double tol) const {
    check(t.size()==local_.size()&&cv.size()==t.size()&&std::isfinite(dt)&&dt>0,"invalid step dimensions/time");
    std::vector<double> mass(t.size());for(std::size_t i=0;i<t.size();++i){check(std::isfinite(t[i])&&t[i]>0&&std::isfinite(cv[i])&&cv[i]>0,"invalid previous temperature/capacity");mass[i]=cv[i]*local_[i].area/dt;check(std::isfinite(mass[i])&&mass[i]>0,"invalid cell mass/time");}return solve(t,mass,tol);
}
HybridHeatResult2D HybridHeatOperator2D::evaluateAtCells(const std::vector<double>& t,double tol) const {
    check(t.size()==local_.size(),"invalid fixed cell dimensions");for(double x:t)check(std::isfinite(x),"nonfinite fixed cell value");
    return solve(t,std::vector<double>(t.size()),tol,nullptr,true);
}
void HybridHeatOperator2D::visitSteadyTraceMatrix(const MatrixVisitor& v) const {solve(std::vector<double>(local_.size()),std::vector<double>(local_.size()),1e-12,&v);}
void HybridHeatOperator2D::visitLocalMatrices(const std::function<void(std::size_t,std::size_t,std::size_t,double)>& v) const {
    for(std::size_t c=0;c<local_.size();++c){const auto& l=local_[c];for(std::size_t i=0;i<l.faces.size();++i)for(std::size_t j=0;j<l.faces.size();++j)v(c,i,j,l.a[i*l.faces.size()+j]);}}
HybridHeatResult2D HybridHeatOperator2D::solve(const std::vector<double>& previous,const std::vector<double>& mass,double tol,const MatrixVisitor* visitor,bool fixedCells) const {
    check(std::isfinite(tol)&&tol>0&&tol<=1e-2,"invalid solve tolerance");
    std::vector<std::pair<std::size_t,std::size_t>> edges;
    for(const auto& l:local_)for(auto f:l.faces)for(auto h:l.faces)if(unknown_[f]!=absent&&unknown_[h]!=absent&&unknown_[f]!=unknown_[h])edges.emplace_back(unknown_[f],unknown_[h]);
    HybridHeatResult2D out;out.trace.resize(boundary_.size());out.traceUnknowns=count_;
    std::vector<double> x(count_);
    if(count_){detail::SparsePattern2D pattern(count_,edges);detail::SparseSystem2D sys(pattern);
        for(std::size_t c=0;c<local_.size();++c){const auto& l=local_[c];const auto n=l.faces.size();const double pivot=l.sum+mass[c];
            for(std::size_t i=0;i<n;++i){const auto row=unknown_[l.faces[i]];if(row==absent)continue;
                sys.rhs[row]+=l.b[i]*(fixedCells?previous[c]:mass[c]*previous[c]/pivot);
                for(std::size_t j=0;j<n;++j){const double a=l.a[i*n+j]-(fixedCells?0:l.b[i]*l.b[j]/pivot);const auto f=l.faces[j],col=unknown_[f];
                    if(col==absent)sys.rhs[row]-=a*boundary_[f].temperature;
                    else if(row==col)sys.diag[row]+=a;else sys.add(row,col,a);
                }
            }
        }
        for(std::size_t f=0;f<boundary_.size();++f)if(unknown_[f]!=absent)sys.rhs[unknown_[f]]-=boundary_[f].flux;
        if(visitor){for(std::size_t i=0;i<count_;++i){(*visitor)(i,i,sys.diag[i]);for(auto q=pattern.rows[i];q<pattern.rows[i+1];++q)(*visitor)(i,pattern.columns[q],sys.off[q]);}return out;}
        detail::LinearWorkspace2D work(count_);out.iterations=sys.solvePressure(x,work,detail::LinearPressureMethod2D::Jacobi,tol);
        std::vector<double> ax;sys.apply(x,ax);for(std::size_t i=0;i<count_;++i)ax[i]-=sys.rhs[i];out.traceResidualNorm=detail::linearNorm(ax);out.traceRhsNorm=detail::linearNorm(sys.rhs);
    }
    for(std::size_t f=0;f<out.trace.size();++f)out.trace[f]=unknown_[f]==absent?boundary_[f].temperature:x[unknown_[f]];
    out.temperature.resize(local_.size());out.faceHeatFlux.resize(boundary_.size());out.cellResidual.resize(local_.size());
    std::vector<double> otherFlux(boundary_.size());
    for(std::size_t c=0;c<local_.size();++c){const auto& l=local_[c];const auto n=l.faces.size();double rhs=0;for(std::size_t i=0;i<n;++i)rhs+=l.b[i]*(out.trace[l.faces[i]]-previous[c]);const double t=fixedCells?previous[c]:previous[c]+rhs/(l.sum+mass[c]);check(std::isfinite(t),"nonfinite cell solution");out.temperature[c]=t;
        double residual=mass[c]*(t-previous[c]);
        for(std::size_t i=0;i<n;++i){double q=0;for(std::size_t j=0;j<n;++j)q-=l.a[i*n+j]*(out.trace[l.faces[j]]-t);
            residual+=q;out.dissipation-=(out.trace[l.faces[i]]-t)*q;
            if(owner_[l.faces[i]]==c)out.faceHeatFlux[l.faces[i]]=q;else otherFlux[l.faces[i]]=q;
        }
        out.maximumLocalCellBalance=std::max(out.maximumLocalCellBalance,std::abs(residual));
    }
    for(std::size_t f=0;f<boundary_.size();++f){const double q=out.faceHeatFlux[f];out.cellResidual[owner_[f]]+=q;if(neighbour_[f]){out.cellResidual[*neighbour_[f]]-=q;out.maximumFaceFluxJump=std::max(out.maximumFaceFluxJump,std::abs(q+otherFlux[f]));}}
    out.stageTemperature=out.temperature;
    const double solverBudget=1e-13+tol*out.traceRhsNorm;
    for(std::size_t c=0;c<local_.size();++c) {
        const double stageBalance=mass[c]*(out.stageTemperature[c]-previous[c])+out.cellResidual[c];
        out.maximumSharedStageCellBalance=std::max(out.maximumSharedStageCellBalance,std::abs(stageBalance));
        if(mass[c]>0) {
            // Flux-form accepted value, never a clipped temperature. All cells
            // use the SAME actual-face flux and the same uniform physical dt.
            const double final=previous[c]-out.cellResidual[c]/mass[c];
            check(std::isfinite(final),"nonfinite conservative final value");
            double fluxScale=0;for(auto f:local_[c].faces)fluxScale+=std::abs(out.faceHeatFlux[f]);
            const double budget=std::sqrt(static_cast<double>(local_[c].faces.size()))*solverBudget/mass[c]+
                4096*std::numeric_limits<double>::epsilon()*(std::abs(previous[c])+std::abs(out.stageTemperature[c])+std::abs(final)+fluxScale/mass[c]);
            const double difference=std::abs(final-out.stageTemperature[c]);
            out.maximumStageToFinalDifferenceK=std::max(out.maximumStageToFinalDifferenceK,difference);
            out.maximumStageToFinalBudgetRatio=std::max(out.maximumStageToFinalBudgetRatio,difference/budget);
            check(difference<=budget,"stage/conservative-final difference exceeds algebraic budget");out.temperature[c]=final;
        }
        const double balance=mass[c]*(out.temperature[c]-previous[c])+out.cellResidual[c];
        out.maximumCellBalance=std::max(out.maximumCellBalance,std::abs(balance));
    }
    return out;
}
}
