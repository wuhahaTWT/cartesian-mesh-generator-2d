"""Illustrate existing native reports; no new solver or mesh acceptance gate."""
from pathlib import Path
import hashlib,json,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
from matplotlib.font_manager import FontProperties
root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root/'tools/flow'));import native_mesh as nm
font=FontProperties(fname='/Library/Fonts/Arial Unicode.ttf');plt.rcParams.update({'font.family':font.get_name(),'axes.unicode_minus':False,'font.size':12})
base=root/'outputs/annulus-cross';out=root/'artifacts/current'
fig,axes=plt.subplots(1,2,figsize=(11.4,6.3),layout='constrained');records=[]
for ax,stage,title in zip(axes,['raw','solver'],['GEOS 原始切割：6,008格','共用后处理结束：6,208格']):
 path=base/('geos-canonical128.'+stage+'.cm2d');q=base/('geos-canonical128.'+('raw-quality' if stage=='raw' else 'quality')+'.json')
 m=nm.read_cm2d(path);d=json.loads(q.read_text());metrics=d['metrics']
 polys=[[m.vertices[v] for v in c.vertices] for c in m.cells]
 ax.add_collection(PolyCollection(polys,facecolors='#e0edf5',edgecolors='#4c697b',linewidths=.85))
 ax.add_collection(LineCollection([(m.vertices[e.v0],m.vertices[e.v1]) for e in m.edges if e.neighbour<0],colors='#273e48',linewidths=1.4))
 ax.set_aspect('equal');ax.set_xlim(-.55,-.39);ax.set_ylim(.15,.31);ax.set_xlabel('x / m');ax.set_ylabel('y / m');ax.set_title(title,fontsize=15)
 ax.text(.02,-.22,f"全网格：最小内角 {metrics['min_interior_angle_deg']:.2f}°；最大非正交 {metrics['max_non_orthogonality_deg']:.2f}°\n最大长宽比 {metrics['max_cell_aspect']:.2f}；当前质量门失败项 {d['issue_count']}",transform=ax.transAxes,fontsize=11,va='top')
 records.append(dict(stage=stage,mesh=str(path.relative_to(root)),meshSha256=hashlib.sha256(path.read_bytes()).hexdigest(),quality=str(q.relative_to(root)),qualitySha256=hashlib.sha256(q.read_bytes()).hexdigest(),nativeQualityReport=d))
fig.suptitle('同一原多边形圆环：后处理改善部分指标，同时恶化形状指标',fontsize=17)
fig.savefig(out/'cutcell-postprocess-comparison.png',dpi=170)
record=dict(format='cutcell-treatment-review-evidence-v1',baseCommit='6223f2a',stages=records,processing=json.loads((base/'geos-canonical128.metrics.json').read_text()),notes=['Existing native metric reports only; this is not a new mesh/solver acceptance gate.','Raw input fails existing quality and is not presented as solver-ready.','Postprocessing includes small-cell agglomeration, convex partition and local repartition; this two-stage snapshot does not isolate the contribution of each operation.','Both stages retain original physical walls.','The original-versus-GEOS solver comparison shares native postprocessing and cannot exclude defects in that stage.'],policyHistory=dict(minInteriorAngleDeg=.5,minNeighbourAreaRatio=.01,originalCutFractionMergeTrigger=.1,solverQualityPolicyIntroduced='10fa0f7 (2026-08-24), export gate',areaRatioGateIntroduced='f604042 (2026-08-24)',smallCellPolicyIntroduced='1300ec2 (2026-08-19)'),codeFindings={'strictConvexFilter':'src/quality/SolverTopology2D.cpp:34-52; original 5c46daa, strengthened 4b0af3f', 'partition':'src/quality/SolverTopology2D.cpp:606-643; non-preserved non-rectangle non-strictly-convex sources are partitioned', 'geometryImplication':'The fluid-side boundary around a convex polygonal solid has reflex corners. Preserving multiple physical facets in one integration polygon conflicts with the current strict-convex default filter.', 'smallCellSelection':'src/stabilization/SmallCell2D.cpp:94-166; original cut area fraction selects cells; neighbour ranking uses non-small/shared length/alpha/area, not full post-merge solver conditioning.', 'scope':'Current-source mechanism plus existing before/after snapshots; not a per-diagonal lineage proof or a new flow run.'},status='Mesh treatment is now the first repair priority; no production algorithm or threshold changed by this literature review.')
(out/'cutcell-treatment-review.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n')
print(out/'cutcell-postprocess-comparison.png')
