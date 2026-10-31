"""Plain-language front page; technical identifiers stay in expandable details."""
from pathlib import Path
import html,json
from study import BASE,ROOT,FLOW,read,dump,sha
import subprocess

def collect():
 out={'native':[],'fluent':[],'probe':None,'provenance':{'sourceBaseCommit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'workingTreeSnapshot':True,'flowBinarySha256':sha(FLOW),'gpuPluginSha256':sha(ROOT/'build/Release/cartmesh2d_torch_pressure.dll'),'scope':'single-machine exploratory measurements; full GPU only offloads pressure; device telemetry includes desktop'}}
 for p in sorted((BASE/'runs').glob('gpu-study-v2-*.json')):
  r=read(p);prefix=ROOT/r.get('prefix','')
  resultFile=Path(str(prefix)+'.json');result=read(resultFile) if resultFile.exists() else None
  if result is not None:
   result={k:v for k,v in result.items() if not isinstance(v,(dict,list))}
   result.update(fullLocalRecord=str(resultFile.relative_to(ROOT)),fullRecordSha256=sha(resultFile))
  out['native'].append({'run':r,'performance':read(Path(str(prefix)+'.performance.json')) if Path(str(prefix)+'.performance.json').exists() else None,'result':result})
 for p in sorted((BASE/'gpu-fluent').glob('*/result.json')):out['fluent'].append(read(p))
 p=BASE/'gpu-fluent/probe3d-result.json'
 if p.exists():out['fluent3dProbe']=read(p)
 p=BASE/'gpu-v2/validation.json'
 if p.exists():out['validation']=read(p)
 p=BASE/'gpu-v2/field-differences.json'
 if p.exists():out['fieldDifferences']=read(p)
 p=BASE/'gpu/probe2d/result.json'
 if p.exists():out['probe']=read(p)
 p=BASE/'gpu-v2/replay65k.json'
 if p.exists():
  try:out['pressureReplay']=read(p)
  except json.JSONDecodeError:out['pressureReplayStatus']='incomplete; see retained run receipt'
 return out

def front(data):
 g=data.get('gpuFollowup') or collect()
 def esc(x):return html.escape(str(x))
 def num(x):return '未取得' if x is None else f'{x:,.2f}' if isinstance(x,(float,int)) else esc(x)
 def table(head,rows):return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+esc(x)+'</th>' for x in head)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'
 text='<section><h2>先看实测结论</h2><ul><li><b>GPU已实际运行，但现有压力模块更慢：</b>同65,536格、同SIMPLE/IC(0)、同5轮，求解CPU 12.81秒，GPU压力路径74.28秒。两边场相同到浮点舍入量级。</li><li><b>600秒慢在方程内核：</b>原生coupled的第一轮就约125秒；它不是花10分钟读网格。</li><li><b>最少10轮有冗余，但不是主要瓶颈：</b>9,076格可在6轮通过全部原门，完整时间10.78→9.81秒。</li><li><b>上一轮不是“6万格求解5分钟”：</b>Fluent方程迭代约39秒，完整流程约326秒；配置与导出流程需要改进。</li></ul></section><section><h2>这次究竟在比较什么</h2><p>每行是一项真实试验。<b>背景分辨率是切割前的方格数；流体单元数才是实际求解规模。</b>例如160×160背景切掉固体后，最终只有9,076个流体单元。失败、超时和短程测试不算完整收敛。</p>'
 text+=table(['要回答的问题','应该固定什么','只改变什么','先看哪几列'],[
 ['GPU是否更快','同一网格、初值、停止条件、数值模型','CPU / GPU执行路径','是否完整收敛 → 求解秒 → 场差异 → 启动和保存秒'],
 ['网格分布是否更好','同一物面、求解器和控制；接近单元预算','近壁细分、聚合、贴壁参考','质量通过 → 壁面大格 → 速度/压力诊断'],
 ['求解算法哪里慢','同一网格、边界和误差要求','SIMPLE / coupled / 线性加速方法','每轮成本、线性迭代数、总轮数'],
 ['为什么超过600秒','固定已失败输入','1轮/2轮前缀与阶段计时','最后完成轮、线性方程耗时；不能仅凭超时判发散']])
 text+='<p>下面的CPU/GPU补测是串行执行。最少轮数试验只改变等待轮数，不放宽任何残差、场变化、守恒或最终线性检查；产品默认仍为10轮。</p></section>'
 text+='<section><h2>名字翻译：不用先记代号</h2>'
 text+=table(['报告代号','中文含义','它改变了什么'],[
 ['shape160','原物面切割后，按形状合并；160×160背景，最终9,076格','网格方案'],
 ['graded128','靠墙先细分一级再合并；128×128背景，最终9,000格','网格方案；与上一行并非同背景单因素'],
 ['graded-indexed512','靠墙先细分再合并；512×512背景，最终104,968格','高密度方案；indexed仅加速几何距离查询'],
 ['gmsh65k','Gmsh贴壁圆环参考，65,536格','成熟、圆环专用参考网格'],
 ['simple','速度和压力分开求解，再反复修正','求解算法'],
 ['coupled','速度和压力联立求解','求解算法；每轮工作更多'],
 ['ic0 / aggregation','不完全分解 / 聚合多层预条件','线性方程加速方法，不是网格'],
 ['strict / adaptive','每步严格求解 / 前期较松、最终严格认证','线性方程精度安排'],
 ['1e-6 / 1e-8','数值停止阈值；后者更严格','不代表物理误差一定达到此数值']])+'</section>'
 text+='<section><h2>完整时间、求解时间：各在等什么</h2><p><b>完整时间</b>是从启动实验进程到退出。<b>求解时间</b>是实际迭代方程阶段。其余包括网格转换、软件启动、读取、设置边界、初始化、保存检查点和退出。这些步骤有用途，但重复转换和逐壁面远程设置可以避免。上一轮“launchSeconds”还含网格转换，因此不能直接叫纯启动时间。</p><p>上一轮65,536格Fluent：完整约326秒，其中迭代约39秒。它说明整个实验流程慢，不能据此说6万格方程本身要5分钟。新对照复用配置好的case，下面逐段列出。</p>'
 rows=[]
 for r in g['fluent']:
  if r['name'].startswith('prepare'):continue
  stage=r.get('stages',{});runs=r.get('runs',[])
  for i,s in enumerate(runs or [{}]):
   label=('单层3D参考' if r['dimension']==3 else '原二维')+' / '+('GPU' if r['gpuRequested'] else str(r['cpuProcesses'])+'核CPU')+(' / 会话内重算' if i else '')
   rows.append([label,num(stage.get('launch') if i==0 else 0),num(stage.get('readAndSetup') if i==0 else 0),num(stage.get('initialize'+str(i))),num(s.get('solveSeconds')),num(stage.get('saveCaseData') if i==len(runs)-1 else None),'收敛' if s.get('converged') else esc(r.get('error','未完整收敛'))])
 text+=table(['参考与硬件','启动秒','读入/设置秒','初始化秒','求解秒','保存秒','状态'],rows) if rows else '<p>新阶段实测待生成。</p>'
 text+='<p>本轮新的3D GPU命令行启动停在许可服务器连接阶段，未进入可观察的GPU求解；新的2D case复用会话也未连通。均保留为启动失败，未计算CPU/GPU提速；既有成熟CPU结果仍有效。</p>'
 text+='<p>单层3D使用原平面物面、原壁面速度，前后为对称面；它是Fluent GPU限制下的独立对照，不冒充二维原生求解器。二维GPU启动探测结果：'+esc((g.get('probe') or {}).get('error','尚未记录'))+'</p></section>'
 text+='<section><h2>自研 CPU/GPU 与最少轮数补测</h2><p>自研GPU目前只接管压力线性方程，动量、组装和其他逻辑仍在CPU。下面的短程试验用于定位成本，<b>不能把“跑5轮”当成完成整个CFD</b>。内存、命令和全部数值见JSON。</p>'
 groups={k:[] for k in ['CPU/GPU：每对只改执行位置','最少轮数：只改等待下限','CPU加速方法：同精度对照','超时定位：大网格前缀','更多短程记录']}
 for item in g['native']:
  r=item['run'];p=item['performance'] or {};result=item['result'] or {}
  if r.get('kind')!='gpu-study':continue
  labels={'shape160':'9,076格形状合并','gmsh65k':'65,536格Gmsh贴壁参考','graded128':'9,000格近壁先细分'}
  label=labels.get(r['mesh'],r['mesh'])+' / '+('GPU压力' if r['backend']=='torch-cuda' else 'CPU')+' / '+{'simple':'分开修正','coupled':'联立求解'}[r['coupling']]+' / '+{'ic0':'不完全分解','aggregation':'多层加速'}[r['preconditioner']]
  status='完整收敛' if result.get('converged') else '预算结束' if r.get('status')=='timeout' else '短程结束，未收敛' if r.get('exitCode')==2 else '失败'
  category='CPU/GPU：每对只改执行位置' if r['preconditioner']=='ic0' and (r['mesh']=='gmsh65k' and r['maximumIterations']==5 or r['mesh']=='shape160' and r['coupling']=='coupled') else '最少轮数：只改等待下限' if r['mesh']=='shape160' and r['coupling']=='coupled' else '超时定位：大网格前缀' if r['mesh']=='gmsh65k' and r['coupling']=='coupled' else 'CPU加速方法：同精度对照' if r['preconditioner']=='aggregation' else '更多短程记录'
  groups[category].append([esc(label),str(r['minimumIterations']),str(r['maximumIterations']),num(result.get('iterations')),num(r.get('wallSeconds')),num(p.get('solveSeconds')),num(p.get('coupledLinearSolveSeconds') if r['coupling']=='coupled' else p.get('pressureLinearSolveSeconds')),status])
 for title,rows in groups.items():
  text+='<h3>'+title+'</h3>'+table(['网格 / 执行位置 / 算法','最少轮数','最多轮数','实际轮数','完整秒','求解秒','联立/压力线性秒','结果'],rows)
 text+='<p>9,000格SIMPLE的多层CPU完整对照仍1333轮：完整166.42秒、求解165.06秒；上一轮同网格IC(0)为429.88/427.39秒、同为1333轮，速度差约1.4e-12 Uref。后台负载不同，不当成隔离重复提速认证。9,076格完整coupled/IC(0)配对都6轮：CPU求解19.00秒，GPU压力路径163.85秒，速度差约8.7e-15 Uref。</p>'
 text+='<p><b>已定位的慢处：</b>65,536格coupled第一轮求解约125.46秒，其中联立线性124.56秒，内部压力预处理61.00秒；读取网格约0.65秒。197次Krylov/5,610次内部压力迭代说明它主要在忙着算方程，不是读取卡死。600秒预算没有证明必然发散，也没有给出完整收敛时长。</p>'
 text+='<p><b>同5轮的直接CPU/GPU对照：</b>65,536格SIMPLE+IC(0)，CPU求解12.81秒，GPU压力路径74.28秒。GPU实际执行15次压力求解、5次矩阵上传；启动14.30秒、迭代57.27秒、上传0.185秒。两边压力迭代均3,440次，最大速度差约1.3e-15 Uref。当前实现的瓶颈主要是细碎GPU内核和同步/三角求解成本，不能归咎于PCIe上传。这是短程成本诊断，不是完整CFD收敛对比。CPU改用聚合多层预条件后同5轮求解1.80秒，但属于另一加速方法。</p>'
 text+='<p>压力线性时间在coupled联立时间之内，不能重复相加。达到600秒只表示预算耗尽；大网格的前缀、线性迭代与失败原因在展开记录中保留。</p></section>'
 text+='<section><h2>怎样判断数据有没有进步</h2><ol><li>先看是否网格合格、是否真正完整收敛。超时或短程数据仅用于诊断。</li><li>同一问题比较时间时，先固定网格、边界、初值、精度和保存范围；CPU/GPU达到各自停止门还要核对接受场。</li><li>网格更多、残差更小都不自动等于物理更准。本报告压力百分比是对真圆解析参照的诊断，包含多边形壁面模型差异。</li><li>没有统一“10万格应当多少秒”：在本机该圆环上，成熟Fluent的CPU方程阶段已实测约54秒；不能外推到瞬态、湍流或其他网格。</li></ol><p>成熟做法：CPU处理几何/输入输出，GPU长期保留大数组并执行适合并行的数值内核，避免每轮来回拷贝和同步。Fluent官方建议单GPU先配一个CPU进程，更多CPU主要帮助设置和后处理；不是把CPU/GPU百分比手工分成两半。<a href="https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_ug/flu_ug_sec_gpu_solver_starting.html">Fluent官方说明</a>。下一步GPU候选为成熟<a href="https://github.com/NVIDIA/AMGX">NVIDIA AMGX</a>：支持GPU多层预条件、标量/块系统；尚未集成或获得本项目速度资格</p></section>'
 return text
