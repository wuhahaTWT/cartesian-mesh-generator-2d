'use strict';

function exportGuide({ result, rasterImport, flow, thermal, euler, background }) {
  if(background) return `# 完整笛卡尔背景网格

本次 ${result.counts.cells} 个完整单元，保留物体内部。分类为外部、内部、相交；相交单元没有被裁切。

- mesh-preview.png：真实网格预览与几何轮廓。
- *.background.json：完整计算域、输入轮廓、矩形单元坐标/层级/分类。
- *.background.vtk：ParaView四边形网格，classification为0外部、1内部、2相交，level为细化层级。
- *.xy、*.sizing.json、result.json：输入、加密分布和实际参数。

这不是流体求解拓扑。自适应粗细交界可有悬挂节点，物体轮廓仅作为几何叠加；没有OpenFOAM算例或流场结果。需要浸入边界等数值处理后才能用于物体绕流，不能直接输入当前层流求解器。
`;

  const duct = (flow?.summary?.case || thermal?.summary?.flowCase) === 'duct';
  const transient = flow?.summary?.temporalDiscretization === 'backward-euler';
  const cells = Number(result.counts.cells).toLocaleString('en-US');
  const gate = value => (value?.pass ?? value?.valid) === true ? '通过' : (value?.pass ?? value?.valid) === false ? '未通过' : '未检查';
  const convection = flow?.summary?.convection || 'upwind';
  const convectionLabel = convection === 'face-limited-linear' ? '线性迎风（面方向限制）'
    : convection === 'limited-linear' ? '线性迎风（限制重构）' : '一阶迎风';
  const convectionNote = flow && (flow.summary.convectionInferred || flow.summary.convection === undefined)
    ? '（旧摘要缺少格式字段，按一阶迎风推断）' : '';
  const pressureDiscretization = flow?.summary?.pressureDiscretization === 'shared-face-gauss'
    ? '共享面压力' : '旧结果未记录';
  const pressurePreconditioner = flow?.summary?.pressurePreconditionerInferred
    || flow?.summary?.pressurePreconditioner === 'legacy-unspecified'
    || flow?.summary?.pressurePreconditioner === undefined ? '旧结果未记录'
    : flow?.summary?.pressurePreconditioner === 'aggregation'
    ? '多重网格（试验）' : flow?.summary?.pressurePreconditioner === 'ic0'
      ? '标准（IC0）' : '旧结果未记录';
  const viscousStress = flow?.summary?.viscousStress === 'symmetric'
    ? '共享面牛顿应力' : '旧结果未记录';
  const convectionQualification = convection === 'face-limited-linear'
    ? '速度沿面法向、切向进行局部有界重构；这不保证整个耦合解无振荡或所有工况都更准确。'
    : convection === 'limited-linear'
    ? '有界限制重构可减少数值扩散，但不保证所有工况都更准确，也不代表通用二阶或新的工程合格结论。'
    : '一阶迎风为默认稳健格式。';
  return `# 先看这里

本次网格：**${cells} 个单元**。先打开 **mesh-preview.png**，看全景和壁面局部。

| 文件 | 用来做什么 |
| --- | --- |
| mesh-preview.png | 本次真实网格的图片，可查看、汇报；不是仿真输入文件。 |
| *-openfoam 文件夹 | OpenFOAM 算例。网格在 constant/polyMesh；计算前设置边界条件、物性和求解器，再运行 checkMesh。 |
| *.cm2d | CartMesh2D 原生网格。混合模式以 *.hybrid.solver.cm2d 为最终求解网格。 |
| *.vtk | 用 ParaView 查看。混合模式以 *.hybrid.solver.vtk 为最终网格，其余可能是中间过程。 |
| result.json / selection.json | 本次参数、检查结果和自动选参记录。 |
| *.resolution.json / *.solver-quality.json | 实际尺寸、边界层和内部质量的详细数据。 |
${rasterImport ? '| source-image.* / image-outline.png / image-import.json | 原图、确认的轮廓叠加图及尺寸标定记录。 |\n' : ''}| *.xy 及其他 JSON | 输入轮廓或过程诊断，排查问题时保留。 |
${flow ? '| *.flow.json / *.flow.fields.json | 自研二维层流摘要与按最终 CM2D cell id 对齐的速度、运动学压力场。 |\n| *.flow.vtk / *.flow.cells.csv / *.flow.faces.csv | ParaView 流场、逐单元数值及面通量；faces.csv 还含压力与动量通量列。 |\n| *.flow.residuals.csv | SIMPLE 内迭代历史；非定常时仅含最后一个时间步。 |\n' : ''}${transient ? '| *.flow.time-history.csv | 本次全部物理时间步、CFL、动能及受力监测。 |\n| *.flow.checkpoint | 最后接受状态，可配合同一最终网格继续计算；需要保持物性与离散格式。 |\n' : ''}${flow?.summary?.timeStepControl === 'adaptive-cfl-retry' ? '| *.flow.attempt-history.csv | 自动步长全部试算，包括拒绝原因、重试步长、实际 CFL 与收敛指标；time-history 只含已接受步。 |\n' : ''}| flow-incomplete-* | 取消/失败诊断；候选CSV不能当作完成的流场。目录内完整 .checkpoint 可用于续算，忽略 .tmp。 |

内部拓扑：${gate(result.gates?.topology)}；内部 Solver：${gate(result.gates?.solver)}。
${duct ? '本次曲壁通道／喷管工况：最左侧竖直端面为均匀速度入口，最右侧竖直端面为运动学 p=0 出口，其余曲壁及障碍物无滑移。CM2D 的几何物面标签并不全是流动壁面；这些入口/出口条件由原生求解器的 duct 工况指定，使用 OpenFOAM 时仍需分别设置。\n' : ''}
${flow?.summary?.case === 'custom' ? '本次使用命名边界。*.flow.boundaries 保存实际输入，绑定最终网格面的位置、外法向与 owner；摘要含同一组逐面条件。参考速度仅用于归一化，实际入口速度在各边界记录中；出口压力为运动学压力。续算必须保持原名称、类型和数值。当前不支持自定义滑移或出口回流；OpenFOAM 边界仍须另外设置。\n' : ''}
${flow?.summary?.namedWallLoads ? '命名壁面载荷在摘要 namedWallLoads 中：压力与黏性力分项及总力，单位 m³/s²；力矩单位 m⁴/s²，基准原点 (0,0)，逆时针为正。它们表示流体作用于壁面的力/力矩除以密度和单位厚度，乘实际密度与厚度可换成 N 和 N·m；开放壁面分组的压力载荷依赖压力基准。\n' : ''}
${flow?.summary?.timeStepControl === 'adaptive-cfl-retry' ? '自动步长仅控制 CFL 与失败重试，不是时间误差估计；需独立做步长细化。续算目标为绝对物理时间。\n' : ''}
${flow?.summary?.initialVortex ? `初始局部涡 compact-cubic-v1：中心 (${flow.summary.initialVortex.centre.join(', ')}) m，支撑半径 ${flow.summary.initialVortex.radius} m，带符号峰值速度 ${flow.summary.initialVortex.peakSpeed} m/s（正值逆时针）。仅在零时刻施加，整个支撑圆盘须位于流体内；不是持续源项。*.flow.initial.checkpoint 保存真实初始场，可用于独立复核或从零重算；*.flow.checkpoint 继续最后接受的状态，续算不再施加扰动。\n` : ''}
${flow ? `自研${transient ? '非定常' : '稳态'}层流：${transient ? `本次时间推进完成，已接受到 t=${flow.summary.acceptedTime} s；本次 ${flow.summary.completedSteps} 步，最后一步 dt=${flow.summary.dt} s` : flow.summary.converged ? '已收敛' : '到达迭代上限，未收敛'}；流动停止容差 ${flow.summary.tolerance ?? "未记录"}（动量残差与速度/压力变化），连续性门 1e-8；工况 ${flow.summary.case}，${transient ? '最后一步内' : ''}迭代 ${flow.summary.iterations} 次。对流格式：${convectionLabel}${convectionNote}；压力求解：${pressurePreconditioner}；压力离散：${pressureDiscretization}；黏性应力：${viscousStress}；压力 p 的单位是 m²/s²。${convectionQualification}${flow.summary.viscousStress === 'symmetric' ? '压力力和黏性力在同一组共享面上积分。' : ''}\n` : ''}
${flow ? `线性迭代精度：${flow.summary.linearPolicy==='adaptive' || flow.summary.adaptiveLinear===true ? '自适应（最终仍须严格复核）' : '固定精度'}；速度松弛系数：${flow.summary.velocityRelaxation ?? '旧结果未记录'}。这些是数值迭代设置，不改变物理时间步或最终停止容差。\n` : ''}
${thermal ? `## 温度结果

最后完整温度场：t=${thermal.summary.acceptedTime} s，${thermal.summary.minValue.toPrecision(6)}–${thermal.summary.maxValue.toPrecision(6)} K。先看 **temperature-preview.png**。

当前结果目录：${thermal.files['.json'].replace(/thermal\.json$/, '')}

- thermal.vtk / thermal.cells.csv：ParaView 场与逐格温度（K）。
- thermal.thermal-history.csv：每个已接受时间步的残差、温度积分和守恒平衡。
- thermal.thermal.checkpoint：温度与流动的联合续算状态；carrier.checkpoint 仅供诊断。
- boundary.csv / desktop-state.json：真实逐面热边界、请求参数与完成/失败状态。
- thermal-run-*：各次计算的独立目录；以 desktop-state.json 的 complete 为完整结果，failed/cancelled 目录中的候选场不能冒充完成结果。

D=k/(rho cp)，源项为 Q/(rho cp)，通量为向外 q/(rho cp)。这是单向恒物性热输运，没有浮力、温度反馈、辐射或共轭传热。
` : ''}
${euler ? `## 可压 Euler 结果

最后完整可压场：t=${euler.summary.time} s，本次${euler.summary.acceptedSteps}个接受步。先看 **euler-preview.png**。

当前完整结果清单：${euler.manifest}

- euler.cells.csv / euler.fields.json / euler.vtk：密度kg/m³、绝对压力Pa、速度m/s、温度K、Mach数；rhoE为总能量密度J/m³。
- euler.faces.csv：实际面上质量、两分量动量和总能量通量，按owner外向、每单位厚度；检查点保存守恒场。
- euler.history.csv：每个已接受时间步的物理时间、声学CFL、密度/压力下限、质量与能量。
- desktop-state.json + euler.checkpoint：恢复相同网格后可在App中载入的续算清单与状态；两文件需要保留在同一目录。
- euler-run-* 各自独立；只有 desktop-state.json 中 complete 是到达目标时间的完整结果。failed/cancelled 中保存的接受状态可续算，不能冒充完成结果。

当前为一阶Rusanov、前向欧拉无黏理想气体模型，物面自由滑移。没有黏性、热传导或湍流，光滑涡精度研究尚未通过。可压绝对压力不能当成不可压的运动学压力，目标时间不表示稳态。
` : ''}
**导出不等于通过外部检查**：本次打包没有运行 OpenFOAM / Fluent 检查。${flow ? '上面的收敛状态只适用于本次自研层流离散工况，不代表其他流动工况。' : thermal ? '同步流动保存在温度联合状态中；本次没有独立流场显示包。' : euler ? '上述可压结果只适用于本次明确的理想气体无黏工况。' : '本次没有自研流场结果。'}使用哪个 CFD 软件，就在该软件中检查导入后的网格。
`;
}

module.exports = { exportGuide };
