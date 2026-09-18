'use strict';

function exportGuide({ result, rasterImport, flow }) {
  const cells = Number(result.counts.cells).toLocaleString('en-US');
  const gate = value => (value?.pass ?? value?.valid) === true ? '通过' : (value?.pass ?? value?.valid) === false ? '未通过' : '未检查';
  const convection = flow?.summary?.convection || 'upwind';
  const convectionLabel = convection === 'limited-linear' ? '线性迎风（限制重构）' : '一阶迎风';
  const convectionNote = flow && (flow.summary.convectionInferred || flow.summary.convection === undefined)
    ? '（旧摘要缺少格式字段，按一阶迎风推断）' : '';
  const pressureDiscretization = flow?.summary?.pressureDiscretization === 'shared-face-gauss'
    ? '共享面压力' : '旧结果未记录';
  const convectionQualification = convection === 'limited-linear'
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
${flow ? '| *.flow.json / *.flow.fields.json | 自研二维稳态层流摘要与按最终 CM2D cell id 对齐的速度、运动学压力场。 |\n| *.flow.vtk / *.flow.cells.csv / *.flow.faces.csv | ParaView 流场、逐单元数值及面通量；faces.csv 还含压力与动量通量列。 |\n| *.flow.residuals.csv | SIMPLE 迭代连续性、速度变化与动量残差历史。 |\n' : ''}| flow-incomplete-* | 被取消或失败求解保留的诊断目录；不是有效求解结果。 |

内部拓扑：${gate(result.gates?.topology)}；内部 Solver：${gate(result.gates?.solver)}。
${flow ? `自研稳态层流：${flow.summary.converged ? '已收敛' : '到达迭代上限，未收敛'}；工况 ${flow.summary.case}，迭代 ${flow.summary.iterations} 次。对流格式：${convectionLabel}${convectionNote}；压力离散：${pressureDiscretization}；压力 p 的单位是 m²/s²。${convectionQualification}\n` : ''}
**导出不等于通过外部检查**：本次打包没有运行 OpenFOAM / Fluent 检查。${flow ? '上面的收敛状态只适用于本次自研层流离散工况，不代表其他流动工况。' : '本次没有自研流场结果。'}使用哪个 CFD 软件，就在该软件中检查导入后的网格。
`;
}

module.exports = { exportGuide };
