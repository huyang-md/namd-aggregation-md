# Changelog

本文件记录所有对外可见的变更。格式遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)，
版本号遵循语义化版本（破坏行为 → 主版本，新功能 → 次版本，修复 → 补丁）。

## [Unreleased]

### Added
- `aggregate.py --vmd-out`：生成按簇着色的 VMD Tcl 脚本——一簇一色、单体灰色、
  水/离子暗色小点，在 VMD 里 source 即可直接"看见"聚集体；demo 已默认生成
- `aggregate.py --plots`：可选 matplotlib 出图（簇尺寸分布、尺寸与 Rg 时间线、
  径向密度），未安装 matplotlib 时自动跳过；CSV 仍是主输出
- `examples/demo_stacking.py`：堆叠几何层（`analyze.py`）的零依赖 demo——
  合成 J 型堆积数据，2 分钟跑通 J/H 判别、层间距与 RDF
- Issue 模板（bug / feature request）

### Changed
- `.gitignore` 忽略规则推广为 `examples/*-case/`：任何 demo 的产物目录自动覆盖，
  不再逐个手写

## [1.0.1] - 2026-11-22

### Added
- `CITATION.cff`：仓库页出现 "Cite this repository" 按钮

### Changed
- `.gitignore` 显式忽略 `examples/demo-case/`，不再依赖 `*.dcd`/`*.psf` 通配符覆盖

### Fixed
- `examples/demo.py` 结尾提示改为绝对路径，从任何工作目录复制粘贴都能运行

## [1.0.0] - 2026-11-22

### Added
- 首次公开发布。体系搭建：`pipeline.py` 一条命令走完定盒子 → Packmol → psfgen →
  加水加离子 → 生成 conf → 体检；`scaffold.py` / `inspect_molecule.py` / `preflight.py` /
  `restart.py` / `inspect_traj.py` / `aggregate.py` / `analyze.py` / `compare_systems.py` /
  `study_report.py`
- 通用自组装分析（任何溶质）与可选堆叠几何分析（J/H 判别、层间距、滑移角）
- 13 步可复现教程（`docs/tutorial.md`）与可打勾清单（`docs/checklist.md`）
- 62 个自合成数据的单元测试与 GitHub Actions CI（Python 3.10 / 3.12）
- 零依赖 demo（`examples/demo.py`）
