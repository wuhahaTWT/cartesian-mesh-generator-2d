"""Build a local HTML/PNG report from actual retained runs only."""
from pathlib import Path
import argparse,csv,datetime,html,json,math,sys
from study import ROOT,BASE,read,dump,sha
sys.path.insert(0,str(ROOT/'tools/flow'));sys.path.insert(0,str(ROOT/'tools/benchmarks'))
import native_mesh as nm
from annulus_cross import errors

OUT=ROOT/'artifacts/current';STEM='pc-mesh-study'
def load(path):return read(path) if path.exists() else None
def collect():
    data=dict(updatedUtc=datetime.datetime.now(datetime.timezone.utc).isoformat(),environment=load(BASE/'environment.json'),
        meshes=[],solves=[],limitations=[
        'Exploratory single-machine measurements, not isolated or repeated performance certification. Mesh generation overlapped some setup/early solver work.',
        'Native strict combines momentum, continuity, iterate changes and final linear certification. Fluent/OpenFOAM/Firedrake use different stopping definitions.',
        'The original polygon with facet-wall motion differs from true circular Couette. Circular diagnostics include model and discretization differences.',
        'Firedrake uses the separate analytic-trace boundary family, P2/P1 triangles or Q2/Q1 quadrilaterals, integration errors and more DOFs per cell. Do not rank it as a same-discretization solver speedup.',
        'GPU and electrical energy are not measured by these CPU experiments. CPU-seconds and RSS are measured costs, not FLOP or kWh estimates.',
        'No product defaults, physical walls or solver kernels were changed. Research wall refinement is opt-in and not a general mesher qualification.'])
    for receipt in sorted((BASE/'runs').glob('*.json')):
        r=read(receipt);prefix=ROOT/r.get('prefix','');kind=r.get('kind');item=dict(run=r)
        if kind=='mesh':
            item.update(name=receipt.stem.removeprefix('mesh-'),metrics=load(Path(str(prefix)+'.metrics.json')),
                quality=load(Path(str(prefix)+'.shape-quality.json')) or load(Path(str(prefix)+'.quality.json')),
                gmsh=load(Path(str(prefix)+'.gmsh.json')))
            if item['gmsh']:
                gr=load(BASE/'runs'/('gmsh-'+item['name']+'.json'))
                item['completeGenerationSeconds']=r.get('wallSeconds',0)+(gr or {}).get('wallSeconds',0)
            else:item['completeGenerationSeconds']=r.get('wallSeconds')
            q=item.get('quality')
            if q and 'issues' in q:
                from collections import Counter
                events=q.pop('issues');q['issueCodeCounts']=dict(Counter(str(x['code']) for x in events))
                q['issueExamples']=events[:8]
                file=Path(str(prefix)+'.shape-quality.json')
                if not file.exists():file=Path(str(prefix)+'.quality.json')
                q['fullLocalReport']=str(file.relative_to(ROOT));q['fullReportSha256']=sha(file)
            data['meshes'].append(item)
        elif kind in ('native','openfoam','firedrake','fluent'):
            item.update(kind=kind,name=r['name'])
            if kind=='native':
                item.update(result=load(Path(str(prefix)+'.json')),performance=load(Path(str(prefix)+'.performance.json')))
                if r.get('exitCode')==0 and item['result'] and item['result'].get('converged'):
                    with Path(str(prefix)+'.cells.csv').open(encoding='utf-8-sig') as f:
                        rows=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
                    item['circleDiagnostics']=errors(rows,r.get('boundaryModel','facet-wall'))
            elif kind=='firedrake':item.update(result=load(prefix/'resource.json'))
            else:
                result=load(prefix/'result.json')
                if kind=='fluent' and result:
                    keys=['status','version','threads','gpuUsed','cells','iterations','residualConverged','finalResiduals',
                        'meshCheckPassed','launchSeconds','solveSeconds','totalSeconds','sessionClosed','caseDataWritten','error','circleErrors','cellMapping','finiteWholeVolumeField','coupling']
                    result={k:result.get(k) for k in keys}
                    result['sourceRecord']=str((prefix/'result.json').relative_to(ROOT))
                item.update(result=result,inner=load(prefix/'inner-resource.json'))
            if r.get('exitCode') not in (None,0):
                log=receipt.with_suffix('.log').read_text(encoding='utf-8',errors='replace')
                item['failureStage']='input-boundary' if 'cannot open explicit boundary' in log else 'runtime-launch' if r.get('exitCode')==127 or 'No such file or directory' in log else 'budget' if r['status']=='timeout' else 'iteration-limit' if kind=='native' and (item.get('result') or {}).get('iterations')==r.get('maximumIterations') else 'solver'
            data['solves'].append(item)
    files=[p for p in BASE.rglob('*') if p.is_file()]
    data['indexedEquivalence']=load(BASE/'indexed-equivalence.json')
    data['sameMeshFieldDifferences']=load(BASE/'same-mesh-fields.json')
    data['hardware']=load(BASE/'hardware.json')
    data['validation']=load(BASE/'validation.json')
    data['sourceSha256']={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'tools/benchmarks/annulus_mesh_cross.cpp',ROOT/'src/quality/SolverTopology2D.cpp',ROOT/'src/fv/Incompressible2D.cpp',*sorted(Path(__file__).parent.glob('*.py'))]}
    data['storage']=dict(totalBytes=sum(p.stat().st_size for p in files),files=len(files),
        largest=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size) for p in sorted(files,key=lambda p:p.stat().st_size,reverse=True)[:12]],
        roles=dict(continueCalculation='Chosen .solver.cm2d/.boundaries and accepted native fields/face flux or mature HDF5 checkpoint/case-data.',
            representativeEvidence='Report/summary, source hashes, commands, residuals and minimal failed inputs/reports.',
            rebuildable='Render PNG/HTML, VTK previews, successful raw construction stages, local experiment image and build binaries; review after a major version.',
            cleanup='No automatic deletion of raw fields, failed inputs or gas-turbine/base images. Keep a compact current report and representative inputs; inspect storage list before cleanup.'))
    dump(OUT/(STEM+'.json'),data);return data

def plots(data):
    import numpy as np
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection,LineCollection
    from matplotlib.colors import LogNorm
    names=['shape160','graded128','gmsh16k','gmsh65k'];fig,axes=plt.subplots(2,4,figsize=(16,8),layout='constrained')
    for col,name in enumerate(names):
        prefix=BASE/'meshes'/name;file=Path(str(prefix)+'.solver.cm2d')
        if not file.exists():
            for row in range(2):axes[row,col].set_title(name+' pending')
            continue
        mesh=nm.read_cm2d(file);areas=[c.stored_area for c in mesh.cells];metric=load(Path(str(prefix)+'.metrics.json'))
        polygons=[[mesh.vertices[v] for v in c.vertices] for c in mesh.cells]
        ratio=np.array(areas)/(metric['interiorMedianArea'] or np.median(areas))
        walls=[(mesh.vertices[e.v0],mesh.vertices[e.v1]) for e in mesh.edges if e.neighbour<0]
        for row in range(2):
            ax=axes[row,col]
            collection=PolyCollection(polygons,array=ratio,cmap='viridis',norm=LogNorm(vmin=max(float(ratio.min()),1e-12),vmax=float(ratio.max())),
                edgecolors='#37495a',linewidths=.16 if row==0 else .6)
            ax.add_collection(collection);ax.add_collection(LineCollection(walls,colors='#e14433',linewidths=1.2))
            ax.set_aspect('equal');ax.set(xlabel='x / m',ylabel='y / m')
            ax.set_xlim((-1.04,1.04) if row==0 else (-.6,-.32));ax.set_ylim((-1.04,1.04) if row==0 else (-.44,-.12))
            if row==0:ax.set_title(f'{name}: {len(mesh.cells):,} cells')
            else:ax.set_title(f'Wall max / bulk median\n{metric["boundaryMaxOverInteriorMedian"]:.3g}')
        fig.colorbar(collection,ax=axes[1,col],label='Cell area / interior median (log)')
    path=OUT/(STEM+'-meshes.png');fig.savefig(path,dpi=125);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
    for item in data['solves']:
        r=item['run']
        if item['kind']!='native':continue
        file=Path(str(ROOT/r['prefix'])+'.residuals.csv')
        if not file.exists():continue
        with file.open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
        if not rows:continue
        it=[float(x['iteration']) for x in rows]
        label=r['mesh']+' '+r['coupling']+' '+str(r['tolerance'])+' '+r['preconditioner']
        axes[0].semilogy(it,[max(float(x['momentumResidual']),1e-30) for x in rows],label=label)
        axes[1].semilogy(it,[max(float(x['velocityChange']),1e-30) for x in rows],label=label)
    for ax,title in zip(axes,['Normalized momentum residual','Velocity change / Uref']):
        ax.set(xlabel='Outer iteration',ylabel=title);ax.grid(alpha=.2)
        if ax.lines:ax.legend(fontsize=6)
    fig.savefig(OUT/(STEM+'-convergence.png'),dpi=130);plt.close(fig)

def number(v):
    if v is None:return '—'
    return f'{v:,.4g}' if isinstance(v,(int,float)) else html.escape(str(v))

def build_html(data):
    meshrows=[];solverrows=[]
    for item in data['meshes']:
        r=item['run'];m=item['metrics'] or {};q=item['quality'] or {};status=r['status']
        if status=='finished':status='通过' if r['exitCode']==0 else '质量拒绝' if r['exitCode']==3 else '失败'
        values=[item['name'],m.get('cells'),status,item['completeGenerationSeconds'],r.get('processTreeCpuSeconds'),
            r.get('peakProcessTreeRssMiB'),m.get('boundaryMaxOverInteriorMedian'),m.get('boundaryP95Area'),q.get('issue_count')]
        meshrows.append('<tr>'+''.join('<td>'+number(x)+'</td>' for x in values)+'</tr>')
    for item in data['solves']:
        r=item['run'];result=item.get('result') or {};perf=item.get('performance') or {};inner=item.get('inner') or {}
        converged=result.get('converged',False)
        if item['kind']=='fluent':converged=result.get('residualConverged',False)
        status='收敛' if converged and r.get('exitCode')==0 else item.get('failureStage',r['status'] if r['status']!='finished' else result.get('status','未确认收敛'))
        cpu=result.get('cpuSeconds') if item['kind']=='firedrake' else inner.get('cpuSeconds',r.get('processTreeCpuSeconds'))
        rss=result.get('peakRssMiB') if item['kind']=='firedrake' else inner.get('peakChildRssMiB',r.get('peakProcessTreeRssMiB'))
        if item['kind']=='openfoam':cpu=inner.get('cpuSeconds');rss=inner.get('peakChildRssMiB')
        diag=item.get('circleDiagnostics') or result.get('circleErrors') or {}
        values=[item['name'],result.get('cells',perf.get('cells')),result.get('totalDofs'),status,
            result.get('iterations',result.get('nonlinearIterations')),r.get('wallSeconds'),
            perf.get('solveSeconds',result.get('solveSeconds')),cpu,rss,
            diag.get('velocityL2OverUref',result.get('velocityRmsOverUref')),
            diag.get('pressureL2OverAnalyticPressureSpan',result.get('pressureRmsOverAnalyticSpan'))]
        solverrows.append('<tr>'+''.join('<td>'+number(x)+'</td>' for x in values)+'</tr>')
    def table(head,rows):return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+x+'</th>' for x in head)+'</tr></thead><tbody>'+''.join(rows)+'</tbody></table></div>'
    env=data.get('environment') or {};storage=data['storage']
    doc='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>PC 网格与成熟求解器实测</title><style>body{font:16px/1.65 system-ui,"Microsoft YaHei",sans-serif;background:#f2f5f8;color:#172737;margin:0}main{max-width:1500px;margin:auto;padding:28px}section{background:white;padding:24px;margin:20px 0;border-radius:12px}h1{font-size:30px}h2{font-size:22px}p{max-width:1100px}.scroll{overflow:auto}table{border-collapse:collapse;font-size:13px;width:100%;white-space:nowrap}td,th{padding:9px;border-bottom:1px solid #dce3ea;text-align:right}td:first-child,th:first-child{text-align:left}th{background:#edf3f9;position:sticky;top:0}img{width:100%;height:auto}code{background:#edf2f7;padding:2px 5px}.notice{border-left:4px solid #dca02e;padding-left:16px}a{color:#1767a3}small{color:#526577}</style><main>
<h1>PC 网格与成熟求解器：实际实验报告</h1>'''
    doc+=f'<p>更新 UTC：{html.escape(data["updatedUtc"])} · 当前源码基点 <code>{html.escape(env.get("sourceCommit","")[:10])}</code>，实验源码哈希及命令见 JSON。表格只读取真实产物，未完成项不计作收敛。</p>'
    doc+='<section><h2>已得到的结论</h2><p>近壁先加密可以控制大单元，但面积变小并不保证压力诊断更好。9076 格纯聚合的原生压力 RMS/真圆解析跨度约 1.87%，约 9000 格近壁加密方案约 3.08%；这两组背景分辨率不同，是接近单元预算的策略对照。面积上限从 4 收紧到 2 会拒绝原先通过的 160/256 档，不能代替改进生成方法。</p><p>Gmsh 映射圆环仅作几何专用成熟参照。104,968格Fluent虽然62轮停止，真圆压力诊断仍约15.91%、局部最大约17.72倍解析跨度，速度RMS/Uref约0.0393%。这些诊断含边界模型差异。原生高密度尝试的预算结束保留为未完成；成熟求解器已经完成的规模和实际成本见下表。不同停止量、边界族、阶数和后台负载不支持直接宣布通用提速或精度资格。</p></section>'
    doc+='<section><h2>实验问题和边界</h2><p>解决纯聚合后近壁单元面积增大的问题，比较旧流程、按形状聚合、近壁先加密再聚合、面积预算、网格相位及成熟 Gmsh 网格。固定原每圆 256 条直边、ν=0.1 m²/s、ρ=1 kg/m³、内壁逐边切向速度 0.5 m/s、外壁静止、零初场。</p><p class="notice">网格合法性、原生质量、外部 checkMesh、离散数值收敛和物理误差分别判断。真圆解析诊断包含多边形壁面模型差异。Firedrake 的解析 trace 边界族单列，不能与原逐边壁面族混算精度或提速。</p></section>'
    doc+='<section><h2>真实网格与近壁单元面积</h2><img src="'+STEM+'-meshes.png" alt="真实网格及近壁面积对照"><p>颜色为面积比，使用每张图自己的完整正值对数范围；红线是原物面。最大面积比直接列出，未裁掉离群单元。</p></section>'
    doc+='<section><h2>网格生成成本和失败</h2><p>完整生成时间含进程、切割、处理、拓扑、质量与导出。Gmsh 行合计成熟生成和原生无修复导入；CPU/RSS 列为原生导入进程树，Gmsh 自身 CPU/RSS 另见 JSON。超时前缀不算失败根因已定位。</p>'
    doc+=table(['组','流体格','状态','完整生成 s','CPU s','峰值 MiB','最大近壁面积比','近壁 P95 面积 m²','质量项'],meshrows)+'</section>'
    doc+='<section><h2>完整求解成本和收敛</h2><p>进程完整时间包括启动/读取/求解/导出，离线字段映射和本报告绘图另计；内核时间由求解器自己报告。CPU 秒是计算成本，不等于墙钟秒。Windows RSS 为同时进程树采样，Linux RSS 为最大子进程 rusage；两种口径不冒充完全相同。Fluent 使用 CPU，Firedrake 是成熟组装/PETSc上的自定义稳态变分参照，单 MPI rank（三角形 P2/P1，四边形 Q2/Q1），OpenFOAM 为单进程。</p>'
    doc+=table(['组','单元数','自由度','状态','外轮/Newton','完整 s','求解内核 s','CPU s','峰值 MiB','速度 RMS/Uref','压力 RMS/解析跨度'],solverrows)+'</section>'
    doc+='<section><h2>实际接受场的压力诊断</h2><img src="pc-mesh-study-fields.png" alt="实际压力偏差"><p>全场保留所有单元和离群压力，先移除最佳面积加权压力常数。对称对数色标仅用于显示，不作为验收门；真圆参照与逐直边运动模型不同。Fluent 读取同次保存的 case/data，并以完整节点和 cell-face 连接映射，ASCII 导出行号不作为单元身份。</p></section>'
    doc+='<section><h2>原生求解的实际迭代历史</h2><img src="'+STEM+'-convergence.png" alt="实际收敛曲线"><p>默认同对流、黏性离散及 0.6/0.25 松弛；参数变化在组名和命令中明确记录。残差达到上限、预算结束和失败候选均不算收敛。</p></section>'
    doc+='<section><h2>方法来源与当前限制</h2><ul>'+''.join('<li>'+html.escape(s)+'</li>' for s in data['limitations'])+'</ul><p>成熟参考：<a href="https://gmsh.info/doc/texinfo/gmsh.html">Gmsh 尺寸场/Transfinite</a> · <a href="https://www.firedrakeproject.org/documentation.html">Firedrake/PETSc</a> · <a href="https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_ug/flu_ug_GridCheck.html">Fluent 网格检查</a>。映射圆环是几何专用参考，不是一般 Cut-cell 生成器。</p></section>'
    doc+=f'<section><h2>产物去向与版本后清理</h2><p>当前实验目录约 {storage["totalBytes"]/1024**2:,.1f} MiB，{storage["files"]:,} 个文件。继续计算保留选定网格、边界、接受场/面通量或 HDF5。代表失败与日志保留；可重建图、VTK、中间成功网格、实验镜像和编译缓存可在大版本收尾时核对清理。没有删除原始场、失败材料或 gas-turbine 环境。</p><p><a href="{STEM}.json">完整数字、命令、资源与存储清单 JSON</a></p></section></main></html>'
    (OUT/(STEM+'.html')).write_text(doc,encoding='utf-8')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--no-plots',action='store_true');a=p.parse_args()
    data=collect()
    if not a.no_plots:plots(data)
    build_html(data);print(OUT/(STEM+'.html'))
