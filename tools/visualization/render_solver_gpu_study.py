#!/usr/bin/env python3
"""Render measured native SIMPLE and optional GPU pressure results."""
import argparse
import hashlib
import html
import json
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools/benchmarks'))
from simple_convergence_study import fields, differences


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--simple', type=Path, default=ROOT/'outputs/simple-pc-study/summary.json')
    parser.add_argument('--gpu', type=Path, default=ROOT/'outputs/gpu-cfd-study/summary.json')
    parser.add_argument('--pressure', type=Path, default=ROOT/'outputs/gpu-pressure-replays')
    parser.add_argument('--output', type=Path, default=ROOT/'artifacts/current/simple-gpu-pc')
    parser.add_argument('--validation', type=Path, default=ROOT/'outputs/solver-gpu-regression/validation.json')
    args = parser.parse_args(); prefix = args.output.resolve(); prefix.parent.mkdir(parents=True, exist_ok=True)
    simple, gpu = load(args.simple), load(args.gpu)
    sr = {x['name']: x for x in simple['runs']}; gr = {x['name']: x for x in gpu['runs']}
    block_gpu, block_cpu = gr['annulus-small-coupled-torch-cuda'], gr['annulus-small-coupled-cpu']
    if block_gpu['accepted'] and block_cpu['accepted']:
        block_gpu['differenceFromMatchedCPU'] = differences(fields(block_gpu['prefix']), fields(block_cpu['prefix']),
            .5, .25+.05/block_gpu['summary']['domainHeight'])
    validation = load(args.validation) if args.validation.exists() else {'status': 'pending final checks'}
    pressure = []
    for path in args.pressure.glob('*.json'):
        item = load(path)
        if item.get('format') == 'cartmesh2d-pressure-replay-v1':
            item['evidenceFile'] = str(path.relative_to(ROOT)); pressure.append(item)
    pressure.sort(key=lambda x: x['cells'])
    tol_names = ['dfg-research-tol-0.05', 'dfg-research-tol-0.01', 'dfg-research-tol-0.001',
                 'dfg-research-tol-1e-06', 'dfg-research-tol-1e-08']
    tolerances = []
    for name in tol_names:
        run = sr[name]; body = next(x for x in run['summary']['namedWallLoads'] if x['name'] == 'cylinder')
        tolerances.append({'tolerance': run['tolerance'], 'iterations': run['summary']['iterations'],
            'wallSeconds': run['wallSeconds'], 'Cd': body['forceX']/.002, 'Cl': body['forceY']/.002,
            'velocityMaxDifferencePercentOfUref': 100*run['differenceFromSameResearchStrictSolution']['maxVelocityComponentOverUref'],
            'pressureMaxDifferencePercentOfScale': 100*run['differenceFromSameResearchStrictSolution']['maxPressureOverReferenceScale']})
    def compact(run):
        return {k: v for k, v in run.items() if k not in ('summary', 'performance')} | {
            'summary': {k: v for k, v in run.get('summary', {}).items() if k in (
                'converged', 'status', 'iterations', 'tolerance', 'momentumResidual', 'continuity',
                'velocityChange', 'pressureChange', 'pressureBackend', 'strictLinearFinal', 'namedWallLoads')},
            'performance': run.get('performance', {})}
    sources = [ROOT/'src/fv/TorchPressureModule2D.cpp', ROOT/'src/fv/FlowGpuPressure2D.cpp',
               ROOT/'include/cartmesh2d/fv/detail/FlowLinearSystem2D.hpp',
               ROOT/'include/cartmesh2d/fv/detail/FlowConvergence2D.hpp', Path(__file__)]
    data = {'format': 'cartmesh2d-simple-gpu-pc-evidence-v1',
            'sourceBaseCommit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'sourceSha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            'hardware': 'RTX 4060 Laptop 8188 MiB, driver 572.70, existing PyTorch 2.12.1+cu126',
            'scope': 'Windows serial native runs; stopping, same-discretization fields and frozen pressure solves; not new physical or large full-CFD qualification',
            'toleranceComparison': tolerances, 'simpleRuns': [compact(x) for x in simple['runs']],
            'gpuFullCFDRuns': [compact(x) for x in gpu['runs']], 'pressureReplay': pressure,
            'pressureStudyManifest': load(args.pressure/'manifest.json'),
            'researchVariantBuild': load(ROOT/'outputs/simple-pc-variants/comparison.json'),
            'validation': validation}
    prefix.with_suffix('.json').write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
    labels = ['5%', '1%', '0.1%', '1e-6', '1e-8']
    ax[0].bar(labels, [x['iterations'] for x in tolerances], color='#28758c')
    ax[0].set(title='Same DFG research discretization', ylabel='SIMPLE outer iterations', xlabel='Normalized stopping tolerance')
    for i, x in enumerate(tolerances): ax[0].text(i, x['iterations']+7, str(x['iterations']), ha='center')
    ax[1].plot(labels[:-1], [x['velocityMaxDifferencePercentOfUref'] for x in tolerances[:-1]], 'o-', color='#ba5d2f')
    ax[1].set(yscale='log', ylabel='Max velocity-component difference / Uref (%)', xlabel='Normalized stopping tolerance',
              title='Difference from same-mesh 1e-8 field')
    ax[1].grid(axis='y', alpha=.2)
    fig.savefig(prefix.with_name(prefix.name+'-stopping.png'), dpi=170); plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    names = ['annulus-small-cpu-jacobi', 'annulus-small-gpu-jacobi', 'annulus-small-cpu-ic0',
             'annulus-small-gpu-ic0', 'annulus-small-cpu-aggregation']
    labels = ['CPU\nJacobi', 'GPU\nJacobi', 'CPU\nIC0', 'GPU\nIC0', 'CPU\nmultigrid']
    values = [gr[n]['wallSeconds'] for n in names]
    ax[0].bar(labels, values, color=['#28758c', '#d27a40', '#28758c', '#d27a40', '#267448'])
    ax[0].set(yscale='log', ylabel='Complete fresh-process CFD seconds', title='484-cell annulus: identical stopping gates')
    for i, value in enumerate(values): ax[0].text(i, value*1.12, f'{value:.3g}', ha='center')
    large = [x for x in pressure if x['cells'] >= 100000]
    for backend, pc, label, color in [('cpu', 'jacobi', 'CPU Jacobi', '#5b82b7'),
                                    ('torch-cuda', 'jacobi', 'GPU Jacobi', '#d27a40'),
                                    ('cpu', 'ic0', 'CPU IC0', '#28758c'),
                                    ('torch-cuda', 'ic0', 'GPU IC0', '#ad5353'),
                                    ('cpu', 'aggregation', 'CPU multigrid', '#267448')]:
        rows = [next(r for r in x['runs'] if r['backend'] == backend and r['preconditioner'] == pc) for x in large]
        ax[1].plot([x['cells'] for x in large], [r['warmMedianSeconds'] for r in rows], 'o-', label=label, color=color)
    ax[1].set(yscale='log', ylabel='Frozen pressure solve seconds (median of 3 warm runs)',
              xlabel='Actual channel CFD cells (one outer step supplies matrix)', title='Linear subproblem, not complete large CFD')
    ax[1].legend(fontsize=8); ax[1].grid(axis='y', alpha=.2)
    fig.savefig(prefix.with_name(prefix.name+'-gpu.png'), dpi=170); plt.close(fig)

    tolerance_rows = ''.join(f"<tr><td>{x['tolerance']:g}</td><td>{x['iterations']}</td><td>{x['velocityMaxDifferencePercentOfUref']:.5g}%</td><td>{x['Cd']:.6f}</td><td>{x['Cl']:.6f}</td></tr>" for x in tolerances)
    work_names = ['annulus-small-product-simple', 'annulus-small-product-coupled', 'annulus-fine-product-coupled',
                  'dfg-research-tol-1e-08', 'dfg-research-anderson-adaptive', 'dfg-research-coupled']
    work_rows = ''.join(f"<tr><td>{html.escape(n)}</td><td>{sr[n]['summary']['iterations']}</td><td>{sr[n]['performance'].get('coupledKrylovIterations',0):,}</td><td>{sr[n]['performance']['pressureIterations']:,}</td><td>{sr[n]['wallSeconds']:.3f}</td></tr>" for n in work_names)
    pressure_rows = ''.join(f"<tr><td>{x['cells']:,}</td><td>{r['backend']}</td><td>{r['preconditioner']}</td><td>{r['iterations']}</td><td>{r['warmMedianSeconds']:.6g}</td><td>{r['trueResidualNorm']:.3g} / {x['nativeAbsoluteTarget']:.3g}</td></tr>" for x in pressure for r in x['runs'] if r['accepted'])
    small, fine = sr['annulus-small-product-coupled'], sr['annulus-fine-product-coupled']
    plain, accelerated, block = sr['dfg-research-tol-1e-08'], sr['dfg-research-anderson-adaptive'], sr['dfg-research-coupled']
    largest = max(pressure, key=lambda x: x['cells'])
    pr = {(x['backend'],x['preconditioner']): x for x in largest['runs']}
    jacobi_speedup = pr['cpu','jacobi']['warmMedianSeconds']/pr['torch-cuda','jacobi']['warmMedianSeconds']
    validation_text = html.escape(json.dumps(validation, indent=2, ensure_ascii=False))
    report = f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>SIMPLE 与 GPU 压力求解实测</title>
<style>body{{font:16px/1.7 system-ui,"Microsoft YaHei",sans-serif;margin:0;background:#f3f5f4;color:#172e36}}main{{max-width:1150px;margin:auto;padding:32px}}h1{{font-size:30px}}h2{{margin-top:38px}}section,.lead{{background:white;padding:24px;border-radius:12px;margin:20px 0}}.lead{{border-left:5px solid #28758c}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{padding:9px;text-align:left;border-bottom:1px solid #d9e1e2}}th{{background:#edf3f3}}img{{width:100%;height:auto}}code{{background:#edf3f3;padding:2px 5px}}.note{{color:#526771}}.warn{{color:#94481e}}.scroll{{overflow-x:auto}}a{{color:#11657a}}</style>
<main><h1>SIMPLE 为什么慢，GPU 到底加速了什么</h1><p class="note">Windows PC · RTX 4060 Laptop · 复用现有 CUDA 12.6 / PyTorch · 真实原生方程与输出</p>
<div class="lead"><b>结论：</b>1%/5% 单步停止量不能保证解已经接近收敛场；coupled 每轮内部做大量线性工作；GPU Jacobi 在大压力矩阵上有效，但当前 CPU 多重网格更快，小网格完整 CFD 用 GPU 明显更慢。CPU 默认、原生质量门和停止门均保留。</div>
<section><h2>1. 当前软件如何停止</h2><p>界面默认 <code>1e-6</code>，研究常用 <code>1e-8</code>。strict 要求至少 10 轮，最大归一化动量残差、速度变化和压力变化都小于该容差，局部连续性和全局流量不平衡小于 <code>1e-8</code>，最后还有严格线性求解检查。coupled 另检查面通量一致性。</p>
<p>速度变化是最大单元速度增量除以 Uref；压力变化除以 <code>Uref² + ν·Uref/H</code>；动量残差按动量对角系数和 Uref 缩放。这不是“最终物理误差百分比”，也不等于其他软件同名 residual。</p>
<p class="warn">下表将本项目归一化门槛直接设为 0.05/0.01；如果老师指“残差相对初值”或“监测量相对变化”，定义另有不同。连续性门未放宽。严格解仍有已知曲壁物理误差，只用作同网格迭代参照。</p>
<div class="scroll"><table><tr><th>停止容差</th><th>外迭代</th><th>速度分量最大差/Uref</th><th>圆柱 Cd</th><th>圆柱 Cl</th></tr>{tolerance_rows}</table></div>
<img src="{prefix.name}-stopping.png"><p>原默认 SIMPLE、只改速度修正、只改压力组合均失败；两项研究组合后，upwind 为 {sr['dfg-ablation-both']['summary']['iterations']} 轮，limited-linear 为 {plain['summary']['iterations']} 轮。研究压力组合改变空间离散，尚未替换默认。1e-6 相对 1e-8 迭代场已很接近，但这不修复原有压力与升力误差。</p></section>
<section><h2>2. coupled 的“10 轮”里包含什么</h2><p>SIMPLE 先解速度，再用近似的局部速度响应修正压力，通常还以 0.6/0.25 松弛逐步更新。这里的 coupled 联立速度两分量和压力，在一轮内用块 FGMRES 反复求解完整冻结响应，并嵌套压力预条件求解；它不使用 SIMPLE 的这两个松弛因子，而用线搜索控制更新。</p>
<p>至少 10 轮是共用停止函数中的下限。小圆环的 {small['summary']['iterations']} 轮实际有 {small['performance']['coupledKrylovIterations']:,} 次块 Krylov 和 {small['performance']['pressureIterations']:,} 次压力迭代；细圆环为 {fine['performance']['coupledKrylovIterations']:,} / {fine['performance']['pressureIterations']:,}。外轮数不能直接代表总工作量或耗时。研究 DFG 的 coupled {block['summary']['iterations']} 轮耗时约 {block['wallSeconds']:.2f} 秒，SIMPLE {plain['summary']['iterations']} 轮约 {plain['wallSeconds']:.2f} 秒；后者加 Anderson/adaptive 为 {accelerated['summary']['iterations']} 轮、约 {accelerated['wallSeconds']:.2f} 秒。它们达到接近的同离散场，但每轮工作和内部线性策略不同。</p>
<div class="scroll"><table><tr><th>实际运行（product / research 分开）</th><th>外轮数</th><th>块 Krylov</th><th>压力迭代</th><th>完整秒数</th></tr>{work_rows}</table></div><p class="note">压力耗时包含在 coupled 线性耗时内，不重复相加。此处是本机串行单次完整成本，不能外推通用性能倍数。</p></section>
<section><h2>3. GPU 正确性与性能</h2><p>可选原生 C++ 模块直接复用已安装 Torch 的 CUDA 运算，Jacobi/IC(0) 用双精度，返回后仍用原生矩阵检查原来的真实残差；不静默转回 CPU。IC(0) 因子在 CPU 构造，三角求解在 GPU。GPU 多重网格尚未实现，表中的 aggregation 是现有 CPU 对照。</p>
<img src="{prefix.name}-gpu.png"><p>完整 484 格圆环：CPU IC(0) 约 {gr['annulus-small-cpu-ic0']['wallSeconds']:.2f} 秒，GPU IC(0) 约 {gr['annulus-small-gpu-ic0']['wallSeconds']:.2f} 秒，均为 {gr['annulus-small-gpu-ic0']['summary']['iterations']} 轮；速度分量差约 {gr['annulus-small-gpu-ic0']['differenceFromMatchedCPU']['maxVelocityComponentOverUref']:.2g} Uref。profile 表明主要额外成本在细碎迭代运算/同步，不只是上传。{largest['cells']:,} 未知量压力子问题，GPU Jacobi 相对 CPU Jacobi 约 {jacobi_speedup:.2f} 倍，但 CPU 多重网格仍略快。不能据此声称大规模完整 CFD 已提速。</p>
<div class="scroll"><table><tr><th>未知量</th><th>后端</th><th>预条件</th><th>迭代</th><th>缓存后中位秒数</th><th>原生真实残差 / 原门槛</th></tr>{pressure_rows}</table></div>
<p class="note">每个冻结矩阵/方法含一次冷解与三次缓存后重复，表列后者中位数；完整 CFD 单独含进程和 GPU 初始化。大规则通道只运行一个外迭代生成真实压力矩阵，并未收敛，也未取得高密度物理精度资格。</p></section>
<section><h2>4. 后续采用范围</h2><p>保留 CPU 默认与 GPU 实验入口。优先继续修复默认 SIMPLE 的压力离散/曲壁一致性，并减少外迭代；大规模线性求解继续研究多重网格及 GPU 上的预条件实现。当前数据不支持把现有 GPU 原型设为通用默认。</p>
<p>完整数字、输入/二进制哈希、命令与实际范围：<a href="{prefix.name}.json">机器可读证据</a>。接受场、原始失败与日志保存在本机 outputs。原有历史研究文件未删除。</p><details><summary>本次软件检查范围</summary><pre>{validation_text}</pre></details></section></main></html>'''
    prefix.with_suffix('.html').write_text(report, encoding='utf-8')
    print(prefix.with_suffix('.html'))


if __name__ == '__main__':
    main()
