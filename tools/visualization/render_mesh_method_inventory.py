from pathlib import Path
import sys, json, hashlib, re, html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Polygon
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'outputs/mesh-method-inventory'
PRE=ROOT/'outputs/airfoil-mesh-preview'
ART=ROOT/'artifacts/current'
STEM='mesh-method-inventory'
plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
plt.rcParams['axes.unicode_minus']=False
wall=np.loadtxt(PRE/'naca2412.xy')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
records=[]; meshes={}

def load_mesh(path):
    lines=path.read_text().splitlines(); nv=int(lines[1].split()[1]); k=2
    vertices=np.array([[float(v) for v in s.split()[1:3]] for s in lines[k:k+nv]]); k+=nv
    ne=int(lines[k].split()[1]); k+=1
    edges=np.array([[int(v) for v in s.split()] for s in lines[k:k+ne]]); k+=ne
    nc=int(lines[k].split()[1]); k+=1; cells=[]; areas=[]
    for line in lines[k:k+nc]:
        tok=line.split(); n=int(tok[4]); cells.append([int(x) for x in tok[5:5+n]]); areas.append(float(tok[3]))
    return vertices,edges,cells,areas,lines[k+nc]



def read_case(id,label,directory,hybrid=False,reused=False):
    candidates=['mesh.hybrid.solver.cm2d','mesh.fallback.solver.cm2d'] if hybrid else ['mesh.solver.cm2d','mesh.failed.solver.cm2d']
    run=json.loads((directory/'run.json').read_text()); log=(directory/'run.log').read_text()
    mesh=next((directory/f for f in candidates if (directory/f).exists()),None)
    rec=dict(id=id,label=label,directory=str(directory.relative_to(ROOT)).replace('\\','/'),command=[run['exe'],*run['args']],exit_code=run['code'],status=run['status'],reused_from_previous_turn=reused)
    if mesh is None:
        rec.update(result='预算结束，无完整网格' if run['status']=='timeout' else '失败，无完整网格',native_pass=False)
        records.append(rec);return
    v,e,c,a,audit=load_mesh(mesh)
    res=json.loads((directory/'mesh.resolution.json').read_text())
    layers=res.get('boundary_layers',{}); cols=layers.get('columns',[])
    layerids=sorted(set(i for x in cols for i in x['solver_cell_ids']))
    reduced=[x for x in cols if x['retained_layers']<x['requested_layers']]
    qpath=directory/('mesh.hybrid.solver-quality.json' if hybrid else 'mesh.failed.solver-quality.json')
    q=json.loads(qpath.read_text()) if qpath.exists() else {}
    passed=q.get('valid','solver_quality=PASS' in log)
    if hybrid and 'mesh_mode=pure_cutcell_fallback' in log:result='退回纯 Cut-cell，未保留混合层'
    else:result='本项目网格检查通过' if passed else '拒绝求解输出'
    rec.update(result=result,native_pass=passed,mesh=str(mesh.relative_to(ROOT)).replace('\\','/'),mesh_sha256=sha(mesh),cells=len(c),stored_topology_audit=audit,layer_cells=layers.get('retained_cells',0),first_layer_wall_fraction=layers.get('first_layer_wall_length_fraction'),full_layer_wall_fraction=layers.get('full_requested_layers_wall_length_fraction'),reduced_columns=len(reduced),issues=q.get('issues',[]),quality_metrics=q.get('metrics',{}))
    if hybrid:
        hr=directory/'mesh.hybrid.json'
        h=json.loads(hr.read_text()) if hr.exists() else {}
        keys=['zero_layer_column_count','q3_accepted_transaction_count','q32_accepted_transaction_count','q33_accepted_transaction_count','q41_accepted_construction_count','q51_outer_transition_radial_subdivision','q51_outer_transition_radial_declined','q52_termination_buffer_radial_committed','q52_termination_buffer_rows']
        rec['strategy_metrics']={k:h.get(k) for k in keys}
    meshes[id]=dict(v=v,e=e,c=c,layers=layerids,reduced=reduced,issues=rec['issues'],rec=rec)
    records.append(rec)

read_case('old-cut','旧默认：小格合并＋凸划分',PRE/'A',reused=True)
read_case('shape','新处理：形状聚合',BASE/'shape-fine')
read_case('shape-coarse','新处理：较粗背景',BASE/'shape-coarse')
read_case('hybrid-default','当前混合：共享构造',PRE/'C',True,True)
variants=[('hybrid-legacy','旧混合：legacy 构造'),('hybrid-q3','终止区合并'),('hybrid-q3-pair','终止区两格重划'),('hybrid-q3-group','终止区分组重划'),('hybrid-q4','构造阶段连接选优'),('hybrid-q5-outer','外过渡宽度匹配'),('hybrid-q5-buffer','终止缓冲宽度匹配'),('hybrid-combined','上述可选策略组合')]
for id,label in variants:read_case(id,label,BASE/id,True)

views={'body':(-.08,1.08,-.15,.19),'nose':(-.025,.12,-.075,.09),'tail':(.94,1.04,-.045,.055)}

def draw(ax,m,view='body'):
    bounds=views[view];v,e,c=m['v'],m['e'],m['c'];ax.set_facecolor('white')
    if m['layers']:ax.add_collection(PolyCollection([v[c[i]] for i in m['layers']],facecolors='#c8eee4',edgecolors='none'))
    ax.add_collection(LineCollection(v[e[:,1:3]],colors='#346b8b',linewidths=.42,alpha=.85))
    for cid in sorted(set(x['cell_id'] for x in m['issues'] if 'cell_id' in x and x['cell_id']<len(c))):
        ax.add_patch(Polygon(v[c[cid]],facecolor='#f9b3a8',edgecolor='#bd3322',linewidth=.9,zorder=3))
    ax.add_patch(Polygon(wall,closed=True,facecolor='#e4e8ed',edgecolor='#172f45',linewidth=1.05,zorder=4))
    if m['reduced']:ax.add_collection(LineCollection([x['wall_segment_endpoints'] for x in m['reduced']],colors='#e36b20',linewidths=2,zorder=5))
    ax.set_xlim(bounds[:2]);ax.set_ylim(bounds[2:]);ax.set_aspect('equal');ax.tick_params(labelsize=8);ax.set_xlabel('x/c',fontsize=9);ax.set_ylabel('y/c',fontsize=9)
    ax.spines[['top','right']].set_visible(False)

fig,axes=plt.subplots(2,2,figsize=(15,7),layout='constrained')
for ax,id in zip(axes.flat,['old-cut','shape','hybrid-legacy','hybrid-default']):
    m=meshes[id];draw(ax,m);r=m['rec'];ax.set_title(f"{r['label']} | {r['cells']:,} 格 | {r['result']}",loc='left',fontsize=12,color='#a43425' if not r['native_pass'] else '#172f45')
fig.suptitle('同一生成器的新旧处理策略 · 同一翼型与外域',fontsize=18)
fig.supxlabel('上排仅切换壁面处理，背景格完全相同；下排仅切换构造身份策略。红色为失败单元，橙色为降层物面。',fontsize=10)
fig.savefig(ART/f'{STEM}.png',dpi=145,facecolor='#f8fafc');plt.close(fig)
fig,axes=plt.subplots(4,2,figsize=(11,13),layout='constrained')
for row,id in enumerate(['old-cut','shape','hybrid-legacy','hybrid-default']):
    for col,view in enumerate(['nose','tail']):
        draw(axes[row,col],meshes[id],view);axes[row,col].set_title(meshes[id]['rec']['label']+' · '+('前缘' if col==0 else '尾缘'),loc='left',fontsize=11)
fig.suptitle('旧划分与新聚合的真实局部 · 失败网格明确标红',fontsize=16)
fig.savefig(ART/f'{STEM}-details.png',dpi=135,facecolor='#f8fafc');plt.close(fig)
fig,axes=plt.subplots(3,3,figsize=(16,10),layout='constrained')
for ax,id in zip(axes.flat,['hybrid-default',*[x[0] for x in variants]]):
    if id not in meshes:ax.text(.1,.5,next(r['result'] for r in records if r['id']==id),transform=ax.transAxes);ax.set_title(id);continue
    m=meshes[id];draw(ax,m,'tail');r=m['rec'];ax.set_title(f"{r['label']}\n{r['cells']:,} 格 · {r['result']}",loc='left',fontsize=10)
fig.suptitle('混合网格全部显式策略 · 同一 4 层请求、同一尾缘窗口',fontsize=17)
fig.savefig(ART/f'{STEM}-hybrid.png',dpi=125,facecolor='#f8fafc');plt.close(fig)

# Compare real segment endpoint tuples exactly, without nearest-neighbour matching.
def edge_geometry(m):
    result=set()
    for edge in m['e']:
        p=tuple(m['v'][edge[1]]);q=tuple(m['v'][edge[2]])
        result.add(tuple(sorted((p,q))))
    return result
reference_edges=edge_geometry(meshes['hybrid-default'])
fig,axes=plt.subplots(3,1,figsize=(13,10),layout='constrained')
for ax,id in zip(axes,['hybrid-q3-group','hybrid-q4','hybrid-combined']):
    m=meshes[id];edges=edge_geometry(m)
    added=sorted(edges-reference_edges);removed=sorted(reference_edges-edges)
    m['rec']['geometry_edge_changes']={'added':len(added),'removed':len(removed),'comparison':'exact endpoint-coordinate pairs, no proximity matching'}
    draw(ax,m)
    if removed:ax.add_collection(LineCollection(removed,colors='#c33c30',linewidths=2.1,linestyles='dashed',zorder=7))
    if added:ax.add_collection(LineCollection(added,colors='#007b65',linewidths=2.1,zorder=8))
    ax.set_title(f"{m['rec']['label']} | 新增{len(added)} / 移除{len(removed)}条几何连接",loc='left',fontsize=12)
fig.suptitle('真正改变连接的位置 · 与当前共享构造默认结果比较',fontsize=17)
fig.supxlabel('绿色：新连接；红色虚线：原连接。按精确端点比较，未做近邻猜配；几何变化不代表流场改善。',fontsize=10)
fig.savefig(ART/f'{STEM}-changes.png',dpi=125,facecolor='#f8fafc');plt.close(fig)

backgrounds=[]
fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
for ax,id in zip(axes,['background-uniform','background-adaptive']):
    data=json.loads((BASE/id/'mesh.background.json').read_text());cells=data['cells'];segments=[];polys=[];colors=[]
    for cell in cells:
        x0,y0,x1,y1=cell['bounds']
        if x1<-.1 or x0>1.1 or y1<-.16 or y0>.2:continue
        p=[[x0,y0],[x1,y0],[x1,y1],[x0,y1]];polys.append(p);colors.append(['#fff','#f6c27b','#bbdef2'][cell['classification']]);segments.extend([[p[i],p[(i+1)%4]] for i in range(4)])
    ax.add_collection(PolyCollection(polys,facecolors=colors,edgecolors='none'));ax.add_collection(LineCollection(segments,colors='#346b8b',linewidths=.45));ax.plot(*np.vstack([wall,wall[0]]).T,color='#172f45',lw=1.3)
    ax.set_xlim(-.1,1.1);ax.set_ylim(-.16,.2);ax.set_aspect('equal');ax.set_title(('均匀' if 'uniform' in id else '自适应')+f"完整背景 | {len(cells):,} 格");ax.set_xlabel('x/c');ax.set_ylabel('y/c')
    backgrounds.append(dict(id=id,cells=len(cells),classification_counts=data['counts'],solver_ready=data['solver_ready'],result='几何产品；不能进入既有 FvMesh 求解器'))
fig.suptitle('完整背景网格保留物体内单元 · 不属于流体求解网格',fontsize=16)
fig.savefig(ART/f'{STEM}-background.png',dpi=130,facecolor='#f8fafc');plt.close(fig)
default_hash=meshes['hybrid-default']['rec']['mesh_sha256']
for r in records:
    if r['id'].startswith('hybrid') and 'mesh_sha256' in r:r['same_file_as_default']=r['mesh_sha256']==default_hash
provenance=json.loads((BASE/'binary-provenance.json').read_text())
meta=dict(source_commit=provenance['head'],geometry_source='examples/complex/naca2412_dense.xy',geometry_sha256=sha(PRE/'naca2412.xy'),families=3,fluid_meshing_families=2,processing_configurations=4,background_modes=2,hybrid_option_forms=6,binaries=provenance,records=records,backgrounds=backgrounds,scope='No CFD, external checkMesh, default changes, removals, or new accuracy thresholds. Existing Release binaries were run; fresh rebuild blocked by local VS instance/file-tracker setup. Previous-turn A and C are explicitly reused with hashes.',research=[dict(method='近壁一级细化＋形状聚合',entry='tools/benchmarks/annulus_mesh_cross.cpp --wall-refine 1 --distance-index 1',scope='固定圆环/方域适配器；尚非通用翼型入口',keep='保留并通用化'),dict(method='GEOS 原始切割',entry='tools/benchmarks/geos_annulus.py',scope='圆环切割对照；后处理仍共用自研旧路径',keep='独立几何参照'),dict(method='Gmsh 映射四边形 / 三角形',entry='tools/benchmarks/pc_meshing/gmsh_reference.py',scope='断言两个256边圆环；不能直接用于单翼型',keep='第三方参照'),dict(method='cfMesh cartesian2DMesh',entry='tools/benchmarks/annulus_cross.py',scope='独立参考，物面有近似，非严格同壁面Cut-cell',keep='第三方参照'),dict(method='AMReX EB',entry='tools/benchmarks/amrex_cutcell',scope='直斜壁切割交叉诊断；不是产品网格核',keep='参照/历史证据')])
update_path=ROOT/'outputs/quality-advisory/result.json'
policy_update_html=''
if update_path.exists():
    update=json.loads(update_path.read_text(encoding='utf8'))
    meta['quality_policy_update']=update
    findings=update['actualFindings']
    policy_update_html=(
        '<section class="alert"><h2>当前行为：经验指标仅作建议</h2>'
        '<p>当前程序已移除经验质量数值对导出和原生读入的硬拦截；桌面也不再以其判定合格。'
        '下面图表保留旧程序的新旧处理对照，其中“拒绝输出”是旧版本行为，不能当作当前状态。</p>'
        f'<p>同一张翼型现在正常导出{update["cells"]:,}格，保留{update["advisoryCount"]}条诊断。'
        f'实际OpenFOAM2606标准检查报告{findings["standardFailedChecks"]}项失败，'
        f'扩展检查报告{findings["extendedFailedChecks"]}项失败；'
        f'涉及{findings["facePyramidIssueFaces"]}个面锥体方向问题和{findings["highSkewFaces"]}个高偏斜面。'
        '没有替换为项目自定meshQualityDict，也没有把工具退出码0当成网格合格。'
        '网格仍然保留，可供目标CFD继续处理。</p>'
        '<p>原始日志和继续计算输入：outputs/quality-advisory/airfoil/foam/；'
        '具体当前行为以docs/CURRENT_STATE_CN.md为准。</p></section>')

(ART/f'{STEM}.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf8')
esc=html.escape
rows=''
for r in records:
    sm=r.get('strategy_metrics',{}); effect=''
    if sm:
        if r.get('same_file_as_default'):effect='最终网格与基础配置逐字节相同。'
        else:
            delta=r.get('geometry_edge_changes',{});effect=f"与基础配置相比：新增{delta.get('added',0)}条连接，移除{delta.get('removed',0)}条；不代表流场更优。"
    detail=esc(json.dumps(r,ensure_ascii=False,indent=2))
    rows+=f"<tr><td>{esc(r['label'])}</td><td>{r.get('cells','—'):,}</td><td>{esc(r['result'])}</td><td>{effect}<details><summary>命令、哈希、原生问题及局部降层</summary><pre>{detail}</pre></details></td></tr>" if isinstance(r.get('cells'),int) else f"<tr><td>{esc(r['label'])}</td><td>—</td><td>{esc(r['result'])}</td><td><pre>{detail}</pre></td></tr>"
research_rows=''.join(f"<tr><td>{esc(x['method'])}</td><td><code>{esc(x['entry'])}</code></td><td>{esc(x['scope'])}</td><td>{esc(x['keep'])}</td></tr>" for x in meta['research'])
report=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>当前网格方法与版本清单</title><style>body{{font:16px/1.7 'Microsoft YaHei',sans-serif;background:#f3f6f9;color:#183146;margin:0}}main{{max-width:1250px;margin:auto;padding:32px}}section{{background:white;padding:24px;margin:20px 0;border-radius:12px}}h1{{font-size:30px}}h2{{font-size:22px}}img{{width:100%;height:auto}}table{{border-collapse:collapse;width:100%}}td,th{{border-bottom:1px solid #d6e1e8;padding:12px;text-align:left;vertical-align:top}}pre{{white-space:pre-wrap;word-break:break-all;font-size:12px}}code{{font-size:13px;word-break:break-all}}.alert{{border-left:5px solid #c45339;padding:12px;background:#fff1ed}}a{{color:#16708e}}</style><main><h1>当前网格方法与版本清单</h1><p>开发线 {esc(meta['source_commit'][:7])}。只生成、检查及绘制真实网格；未启动 CFD。网格方法与内部处理策略分开说明。</p>{policy_update_html}<section><h2>历史选型说明</h2><p class="alert">上一轮 A/B 使用当前程序里的旧默认 Cut-cell 分支，漏加 --shape-agglomeration；C/D 只是两组层数参数。四张图没有覆盖新旧算法。桌面 Cut-cell runtime 也落后于 build/Release：缺少新聚合选项；混合程序两处哈希相同。</p><p>源码和功能可调用不等于已经替换桌面默认。此次不删除方法、不改默认，也不把失败网格当成可求输入。</p></section><section><h2>到底有多少种</h2><p><b>流体网格按两类维护：纯 Cut-cell、贴壁层混合。</b>完整背景是第三个几何产品家族。旧凸划分/新形状聚合属于纯 Cut-cell 内部处理策略；旧交点构造/共享交点身份属于混合网格内部策略，不能分别算作独立版本。背景有均匀与自适应两种模式。</p><p>Q3/Q4/Q5只是开发编号：Q3指终止区合并/重划，Q4指构造阶段选择连接方式，Q5指过渡区宽度匹配。对应6个显式选项，属于同一混合生成器的局部策略，不能当成三个版本。4/8层、粗/细、均匀/自适应尺寸场、尾流带、背景相位和面积预算也不能各算一套独立算法。</p></section><section><h2>新旧处理策略对照（不是四个网格版本）</h2><img src="{STEM}.png"><p>新聚合在同一背景下去掉了大部分原壁面凸分解碎片，但细/粗两档仍各有5项质量问题，均涉及前缘和尾缘：有强反角、边界偏斜、粗档非正交以及非正面心法向距离。明确拒绝正常Solver/OpenFOAM输出。失败网格保留为下一步修复输入。</p><a href="{STEM}-details.png">打开前缘和尾缘等比例大图</a></section><section><h2>历史拒绝原因（当前已改为建议并实际运行checkMesh）</h2><p>“本项目网格检查”只是源码里的一组检查条件，不等于外部checkMesh通过、成熟求解器收敛或物理精度认证。必须分开看以下三层。</p><table><tr><th>检查层次</th><th>依据与局限</th></tr><tr><td>几何/离散适用条件</td><td>面积、连接、面方向须合法；当前有限体积离散要求单元中心至面沿外法向投影为正，否则扩散系数会变号。翼型尾缘两面为约-0.00317m，违反当前算子要求。不能据此外推所有CFD离散都无法处理该多边形。</td></tr><tr><td>项目质量限制</td><td>一些数值参考OpenFOAM默认检查和生成配置，但用途不完全相同。默认checkMesh的skew告警线4，与snappyHexMesh生成配置示例的边界skew20是不同设置。本项目把边界skew4作为阻止输出的条件，是项目政策，未证明它是普适CFD失败界线。原生形状比定义也不同，1000不能直接照搬为同一物理含义。</td></tr><tr><td>项目生成目标</td><td>最小角20°、原生形状比10、内部非正交45°是工程假设，并非行业标准或已完成通用CFD认证的界线。此翼型分别为24.13°/8.83/38.36°，全部满足，失败并非由这三项触发。</td></tr></table><p>细档只剩2个问题单元、5条记录：前缘1项边界偏斜5.26；尾缘凹性111.49°、偏斜7.40及两条负法向距离。把边界偏斜限制改为20也不能消除尾缘的几何问题，本次没有放宽限制。聚合以不移动物面为前提，只做有限范围并集；本例停止于没有可接受并集，尚无通用的局部细化/受控划分补救。</p><p>圆环采用许多小转角直边，本例凹性仅约1.406°，与翼型尾缘强凹口不同。9,076格形状聚合圆环已有自研SIMPLE3037轮、后续OpenFOAM SIMPLE3473轮、Fluent2026R1残差停止2027轮的既有记录；停止量定义不同，不据此比较快慢。该翼型失败网格尚未运行外部checkMesh或成熟CFD，不能假称成熟软件已经判失败。</p><p><a href="https://api.openfoam.com/2312/primitiveMeshCheck_8C_source.html">OpenFOAM默认检查源码</a> · <a href="https://doc.openfoam.com/2306/tools/pre-processing/mesh/generation/snappyhexmesh/meshquality/">OpenFOAM网格生成质量配置示例</a> · 圆环实算记录见原<a href="pc-mesh-study.json">数据JSON</a>。只纠正已有报告表述，不新增报告或修改网格/求解器算法。</p></section><section><h2>所有混合显式策略的实际结果</h2><img src="{STEM}-hybrid.png"><p>同一4层请求、首层0.0025c。局部降层橙色显示；有些开关接收数为0，有些改变中间构造或最终连接。最终文件哈希不同也可能包含身份/编号变化，不能直接推断流场差异。并行执行时间只记录，不作性能排名。</p><img src="{STEM}-changes.png"><p>终止区分组重划只改变3条连接；构造选优及组合的变化涉及中段上下侧。这些是同一生成器的小范围修复，不是独立版本。</p><table><tr><th>路径 / 策略</th><th>流体格数</th><th>本例结果</th><th>实际动作与完整证据</th></tr>{rows}</table></section><section><h2>完整背景：两个独立几何产品</h2><img src="{STEM}-background.png"><p>两者都保留Inside/Outside/Intersected完整单元，solver_ready=false。均匀粗格可能没有完全位于薄翼型内部的单元，这不表示固体被删除。独立MAC浸入边界求解器有自己的计算格和固体辅助未知量，不是把这份背景JSON直接送进原生FvMesh。</p></section><section><h2>已有研究 / 第三方参照，不冒充通用翼型入口</h2><table><tr><th>方法</th><th>入口</th><th>实际范围</th><th>建议</th></tr>{research_rows}</table><p><a href="shape-agglomeration.png">既有圆环新旧聚合图</a> · <a href="annulus-cross-meshes.png">既有第三方/混合参照图</a> · <a href="pc-mesh-study.html">近壁细化与Gmsh实验原报告</a></p></section><section><h2>建议如何保留和更迭</h2><table><tr><th>处理</th><th>对象与理由</th></tr><tr><td>退出常用入口候选，保留回归基线</td><td>旧纯Cut-cell凸分解流程；旧legacy交点构造。前者会制造尖碎片，后者已有共享身份构造替代。新聚合尚不能覆盖尾缘失败，所以现在不能直接删库或盲目切换默认。</td></tr><tr><td>优先继续维护</td><td>形状聚合＋有界近壁细化；当前共享构造混合网格。新聚合需解决强反角/尾缘受控划分；近壁细化研究需通用化。混合层需要保持真实降层与质量告警。</td></tr><tr><td>按用例保留研究开关</td><td>Q3/Q4/Q5和patch事务。单翼型上零接收不证明无用，原窄缝、终止区、喷管等触发例仍有价值。不能按阶段数字删除。</td></tr><tr><td>独立保留</td><td>两种完整背景产品、独立浸入边界研究、第三方几何参照。</td></tr><tr><td>只保留为内部工具 / 历史</td><td>边界层条带构造CLI是混合网格组件，非完整外流网格器；RectilinearMesh是算例/测试建格工具；H2顺序穷举参考用于H3等价测试。历史里程碑不移动。Q1已移除，不恢复。</td></tr></table></section><section><h2>验证与后续消费者</h2><p>现有Release的新聚合与Solver导出两项原生检查通过；每个生成命令的本项目检查结果、退出码和原始日志另列。没有全量测试、外部checkMesh、CFD、物理验证或安装包。尝试重新构建时，本机CMake的VS实例识别和直接MSBuild的路径追踪分别失败；本轮使用既有Release并记录哈希，没有把它描述成重新构建成功。</p><p>通过网格供用户选型后的工况配置；新聚合失败网格及输入供尾缘修复；日志和命令作复现材料；PNG/HTML与成功中间导出可重建。原始材料保留在 outputs/mesh-method-inventory 与 outputs/airfoil-mesh-preview，未删除或覆盖旧结果。</p><p><a href="{STEM}.json">完整JSON、源码入口、原生问题和二进制哈希</a></p></section></main></html>'''
(ART/f'{STEM}.html').write_text(report,encoding='utf8')
print(json.dumps([dict(id=r['id'],cells=r.get('cells'),result=r['result'],same_file_as_default=r.get('same_file_as_default')) for r in records],ensure_ascii=False,indent=2))
