# namd-aggregation-md

NAMD 小分子自组装模拟的**准备 · 体检 · 续跑 · 分析**工具集。

适用场景：几十个分子用 Packmol 随机分散到水里，跑 NAMD 自发组装。
**分析分两层**：

- **通用层**（`aggregate.py`）——任何自组装都能用：聚集体多大、什么形状、什么组成、
  内部是均一还是核壳。胶束、囊泡、纤维、凝胶、多组分共组装、无机杂化都适用。
- **堆叠几何层**（`analyze.py`，可选）——只对有芳香平面的体系：J/H 判别、层间距、倾角。

不适用：蛋白-配体、膜体系、GROMACS 输入（但 DCD/XTC 之外的读取器可以复用）。

许可证：**GPL-3.0-or-later**（见 `LICENSE`）。你可以自由使用、修改、再分发，
基于本工具的衍生作品也须以同一许可证开源。

> **第一次做？先看 [`docs/tutorial.md`](docs/tutorial.md)。**
> 那份 13 步教程从"定电荷态"讲到"能不能下结论"，每步都有验收标准和常见错误。
> 本 README 是命令速查，教程是完整流程。
>
> 想直接照着打勾？用 [`docs/checklist.md`](docs/checklist.md)，标好了哪几步必须你亲手做。

---

## 1. 安装

把整个 `namd-aggregation-md` 文件夹放到下面任一处：

| 位置 | 生效范围 |
|---|---|
| `<项目>\.dsh\skills\` | 只对该项目 |
| `~\.dsh\skills\` | 该用户所有项目 |
| `<项目>\.agents\skills\` 或 `~\.agents\skills\` | agents 约定 |

**依赖分两档：**

| 用途 | 需要什么 | 怎么装 |
|---|---|---|
| 只读轨迹、做分析 | Python 3.10+ 和 numpy | `pip install numpy` |
| 从零搭体系（`pipeline`） | 再加 **Packmol** 和 **VMD** | `pip install packmol`；VMD 从官网装 |

Packmol 从 v100 起官方提供了 pip 包，自带 Windows/macOS/Linux 二进制，**不需要 Julia，也不需要 conda**：

```powershell
pip install packmol          # 提供 packmol 命令
packmol --version
```

如果 `packmol` 不在 PATH 上，把它所在的目录加到环境变量 `MD_TOOLS_DIR`：

```powershell
setx MD_TOOLS_DIR "C:\Users\<you>\AppData\Roaming\Python\Python313\Scripts"
```

分析部分不需要 NAMD、VMD、GROMACS——只读文件、只写文本和 CSV。

验证安装：

```powershell
python scripts/study_report.py --help
python -m unittest discover -s tests        # 62 个测试，应全过
```

---

## 2. 五分钟上手

假设你的体系目录长这样（`psf`/`pdb`/`toppar` 在一起）：

```
my-case/
├─ ionized.psf
├─ ionized.pdb
├─ toppar/
└─ 04_prod.dcd
```

```powershell
$S = "<skill 目录>\scripts"

# 想一条命令走完（需要 Packmol 和 VMD）
python "$S\pipeline.py" --dir my-case --name mycase --ligand-pdb mol.pdb `
       --molecules 50 --target-conc-mM 50 --prod-ns 100

# 或者分步做：
# 1) 建立 case.json 和四个阶段的 conf
python "$S\scaffold.py" --dir my-case

# 2) 自动识别分子的芳香环系，写入分析配置
python "$S\inspect_molecule.py" --dir my-case --write

# 3) 跑前体检
python "$S\preflight.py" --dir my-case

# 4a) 通用分析：聚集尺寸、形状、组成、径向密度（任何自组装体系）
python "$S\aggregate.py" --dir my-case --traj 04_prod.dcd --last-ns 50

# 4b) 可选：堆叠几何（只有芳香平面体系才需要）
python "$S\analyze.py" --dir my-case --traj 04_prod.dcd --last-ns 50 --planes ring1 ring1
```

第 4a 步会在 `my-case\analysis\` 下生成：

- `summary.txt` —— 人先看这个
- `clusters.csv` / `clusters_timeseries.csv` —— 聚集尺寸分布与演化
- `shape.csv` —— 最大簇的 Rg、各向异性、形状分类
- `density_profile.csv` —— 径向密度（判断核壳结构）
- `composition.csv` —— 多组分体系的组成（仅多组分时生成）

第 4b 步另外生成：
- `planes_*_nearest.csv` —— 相邻分子中心距
- `clusters.csv` / `clusters_timeseries.csv` —— 聚集尺寸与演化

---

## 3. 换一个分子要改什么

**通用分析（`aggregate.py`）不需要任何配置。** 它自动找出所有非水非离子的残基当溶质，
多组分体系也一样处理。

**堆叠几何分析（`analyze.py`）才需要指定平面组。** `inspect_molecule.py` 会读拓扑、找环、
按共享边合并成环系、算平面性，把可用的环系写进 `case.json`：

```json
"analysis": {
  "plane_groups": {"ring1": ["C1","C2",...], "ring3": [...]},
  "plane_pair": ["ring1","ring1"],
  "rdf_pairs": [["ring:ring1","ring:ring1"]]
}
```

需要人工确认的只有三件事：

1. **平面性 rms** 应 < 0.2 Å。大于这个值说明那个环系不是平面的（含 sp3 碳），
   拿它当平面会算出没意义的倾角。脚本会自动剔除，但你要确认剔得对。
2. **`plane_pair`** 默认选最大的平面环系与自己做同基团堆叠。想比较不同基团
   （比如"给体 vs 受体"），改成 `["ring1","ring3"]`。
3. **`rdf_pairs`** 可以用原子名（`"C5"`）或**环质心**（`"ring:ring1"`）。
   环质心的好处是不依赖具体原子命名，换分子直接能用。

写完想手改 `case.json`，改完直接重跑 `analyze.py` 即可，不用重新生成 conf。

---

## 4. 命令一览

| 命令 | 干什么 |
|---|---|
| `pipeline.py` | **一条命令走完搭建**：定盒子 → Packmol → psfgen → 加水加离子 → 生成 conf → 体检 |
| `scaffold.py` | 从 `case.json` 生成 min/nvt/npt/prod 四个 conf；`--packmol` 还能按目标浓度算盒子、生成 Packmol 和 VMD solvate 输入 |
| `inspect_molecule.py` | 自动识别溶质与芳香环系，写入分析配置 |
| `preflight.py` | 跑前体检：参数文件、盒子、跨 conf 一致性、`firsttimestep` 与 `.xsc`、磁盘够不够、配置有没有被手改 |
| `restart.py` | 断点续跑：从 `.xsc` 取真实步数、改输出名、算剩余步数和磁盘 |
| `inspect_traj.py` | 轨迹审计：帧数、覆盖区间、缺口/重叠、日志是否正常结束 |
| `aggregate.py` | **通用自组装分析**：聚集尺寸、形状（Rg/各向异性）、组成、径向密度 |
| `analyze.py` | 堆叠几何（可选）：RDF / 滑移角 / 层间距 |
| `compare_systems.py` | 多体系可比性：浓度、帧间隔、时长、力场旋钮 |
| `study_report.py` | 把上面几项汇总成"现在能下哪些结论" |

---

## 5. 三条自检线（都很便宜，能省几小时）

1. **平面性 rms < 0.2 Å** —— 大了说明原子组选错。
2. **相邻分子中心距 > 2 Å** —— 小于这个物理上不可能，说明坐标系或原子选择错了
   （`analyze.py` 会自动报 `SANITY FAIL`）。
3. **堆叠层间距中位数 ≥ 3 Å** —— 小于 3 Å 通常是原子组互相穿插，不是真堆积。

---

## 6. 输入约定

脚本假设 NAMD/CHARMM 的常规布局：

- PSF + PDB + `toppar/` 在同一目录
- 阶段命名 `01_min` / `02_nvt` / `03_equil_npt` / `04_prod`（不匹配就用 `scaffold.py --force` 重建）
- 轨迹 `.dcd` + 日志 `.log` 同名
- 重启文件 `X.restart.coor/.vel/.xsc` 或 `X.coor/.vel/.xsc`

DCD 读取器兼容：NAMD 的 6 值单元胞记录、9 值记录、以及 VMD 导出的无步长信息轨迹
（后者需要 `--frame-ps` 指定帧间隔）。

---

## 7. 测试

```powershell
python -m unittest discover -s tests        # 62 个测试
```

测试不需要任何 MD 数据——DCD/PSF/PDB/XSC 都由测试自己合成。覆盖：

| 测试文件 | 覆盖什么 |
|---|---|
| `test_core.py` | DCD/PSF/PDB/XSC 解析、周期边界、分子解缠绕 |
| `test_generators.py` | 盒子反算、Packmol/VMD 脚本、四个 conf 的生成与回读 |
| `test_checks.py` | conf 解析、漂移检测、`firsttimestep` 校验、可比性判定 |
| `test_analysis.py` | 回转半径与形状分类、跨边界聚簇、簇判据 |
| `test_docs.py` | 编码/BOM、语法、无个人数据、每个工具的 `--help` |

其中几条是**回归测试**，专门盯住真出现过的 bug：

- `test_unwrap_molecule_stays_in_global_frame` —— 分子解缠绕曾把每个分子重新以自身为原点，
  导致所有分子间距离作废
- `test_gather_cluster_rejoins_across_boundary` —— 聚集体跨周期边界时 Rg 曾算出 55 Å（应为 24 Å）
- `test_import` 与 `test_hand_written_conf_is_not_drift` —— 老项目接入后不应产生误报

GitHub Actions 会在每次 push 时于 Linux 上跑 Python 3.10 和 3.12（见 `.github/workflows/tests.yml`）。

---

## 8. 边界

- 脚本**不启动** NAMD，只生成输入和读结果。运行请自行提交。
- 聚集尺寸依赖判据（接触/质心、cutoff），换判据数字会变；报告时写明判据。
- 跨体系比较尺寸前，必须确认浓度、帧间隔、时长、力场旋钮一致——
  `compare_systems.py` 和 `study_report.py` 就是干这个的。

---

## 许可证

Copyright (C) 2026 huyang-md

本程序是自由软件：你可以依据自由软件基金会发布的 GNU 通用公共许可证第 3 版
（或你选择的任何更新版本）的条款重新分发和/或修改它。

本程序的分发是希望它有用，但不提供任何担保，甚至不提供适销性或特定用途适用性的
默示担保。详见 LICENSE（GNU GPL v3 全文）。

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. See LICENSE for the full text.
