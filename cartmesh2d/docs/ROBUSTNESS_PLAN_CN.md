# 2D-R 复杂几何鲁棒性修复

目标：解决正式 2D-0~V 之后压力测试暴露出的复杂单环边界问题，而不降低既有门禁。

## 已知失败基线

- gear_star：Cut-cell Unsupported（大量）
- serpentine_body：Cut-cell Unsupported（少量）
- naca2412_dense：source global topology audit fail
- superellipse_24：source global topology audit fail
- nozzle_profile：PASS

## 修复原则

1. Intersected leaf 不再依赖“全局凹 polygon 直接 Sutherland-Hodgman 裁到 AABB 后仍是单 polygon”的假设。
2. 从局部 embedded-boundary fragments 与 Cartesian cell perimeter fluid intervals 构造有向边界图。
3. 只有闭合、非分叉、正面积的单一 fluid loop 才形成普通 CutCell2D。
4. 真正多个 fluid components 仍显式 Unsupported；不得制造假桥。
5. AABB 裁剪交点在 tolerance 内吸附到精确 cell side 坐标，降低相邻 cell topology 漂移。
6. 原 15 项 Stage 0~6 acceptance 必须保持全绿。
7. 复杂压力目标：gear_star / serpentine_body / nozzle_profile / naca2412_dense / superellipse_24 全部完整 CLI + topology + SVG PASS。

当前分支：`agent/native-2d-robustness`。
