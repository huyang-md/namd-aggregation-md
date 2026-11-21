# 小分子自组装 MD：完整可复现教程

从"一个分子"到"一句能写进论文的结论"，13 个步骤，每步都有**验收标准**和**常见错误**。

适用于：几十个染料/探针分子随机分散到水中，自发组装成聚集体，统计堆积模式与尺寸。
不适用于：蛋白-配体、膜体系。

> 教程里的每个坑都来自真实项目：配体电荷无化学依据（CGenFF 罚分 103）、
> Packmol 导致的浓度漂移、生产阶段力场旋钮误写、分析坐标系错误
> 导致层间距从 3.6 Å 变成 9.9 Å。这些坑现在都有对应的自动检查。

---

## 总览

```
① 定电荷态 ──► ② CGenFF 参数 ──► ③ 定浓度/盒子 ──► ④ Packmol
                                                        │
   ⑬ 结论闸门 ◄── ⑫ 分析 ◄── ⑪ 认分子 ◄── ⑩ 轨迹审计 ◄── ⑨ 运行
        ▲                                                      ▲
        └── ⑬' 可比性检查 ◄────────────────────────────────────┘

                     ⑤ psfgen ► ⑥ 加水/离子 ► ⑦ 生成 conf ► ⑧ 跑前体检
```

四个**闸门**，不过就停下修，不要往下走：

| 闸门 | 在哪 | 不过的后果 |
|---|---|---|
| 罚分 | ② CGenFF | 参数不可信，所有结论都悬空 |
| 预检 | ⑧ preflight | 跑到一半炸，浪费几十小时 |
| 审计 | ⑩ inspect_traj | 在残缺数据上做统计 |
| 可比 | ⑬' compare_systems | 比的是浓度/参数差异，不是化学 |

---

## 命令约定

下面所有命令都用这两个变量，先设好再复制：

```powershell
$py = '<你的 python 路径>'
$S  = '<skill 目录>\scripts'
```

- `$py` 换成你自己的 Python 3.10+（带 numpy 即可）
- `$S` 换成 skill 实际所在目录
- `& $py "$S\xxx.py"` 是 PowerShell 写法；bash 下写成 `python "$S/xxx.py"`

---

## 准备：记录版本

可复现的第一步是**知道你用了什么**。开一个 `RUNLOG.md`，先抄下：

```
NAMD                  : namd2 --version  → 例如 Git-2025-12-04 Linux-x86_64-multicore-CUDA
CGenFF program        : 见 lig.rtf 第二行的 "version 4.0"
CHARMM36 toppar       : 记录来源与下载日期
Packmol               : packmol < /dev/null  → 版本号
VMD                   : vmd -h
Python / numpy        : python -c "import numpy; print(numpy.__version__)"
分析阈值              : slip 32°、lateral 6 Å、stack 8 Å、cluster 4 Å（见第 ⑫ 步）
```

把这些写进论文的方法部分或 SI，是审稿人会查的东西。

---

## ① 定电荷态（最关键的一步）

**这一步错了，后面全错，而且不会报错。**

### 要决定什么

1. 分子是中性还是离子？
2. 如果带电，带在哪个原子上？

依据是**实验用的盐形式**：氯化物盐 → 阳离子；中性化合物 → 净电荷 0。

### 做法

1. 在 ChemDraw 里画出**带电形式**，形式电荷显式标出（`N+`、`O-`）
2. 导出 3D mol2 时**必须保留形式电荷**
3. 验收：mol2 里所有原子电荷之和 = 你设定的净电荷

```powershell
# 检查 mol2 的电荷和
python -c "a=[l.split() for l in open('mol.mol2') if len(l.split())>8 and l.split()[0].isdigit()]; print('net charge =', round(sum(float(x[8]) for x in a),3))"
```

### 常见错误

| 现象 | 原因 |
|---|---|
| 所有原子电荷都是 `0.000` | 导出工具丢了形式电荷（ChemDraw → OpenBabel 常见） |
| 几十个碳全叫 `C` | 2D 结构直接导出，原子名没唯一化 |
| 原子类型是 `C`/`N`/`O` 而非 `C.ar`/`C.3` | 同上，缺少 SYBYL 类型，CGenFF 只能靠距离猜键级 |

**这三种情况都会让 CGenFF 的类比质量大幅下降。** 出现任一，回 ChemDraw 重做。

---

## ② 生成 CGenFF 参数并验收罚分

### 做法

- 上传 mol2 到 **CHARMM-GUI Ligand Reader & Modeler**（`charmm-gui.org/input/ligand`），
  或直接 **cgenff.paramchem.org**
- 填**净电荷**（必须与 ① 一致）、自旋多重度
- 只需要参数时选 "Just generate parameters"

产出：`lig.rtf`（拓扑）、`lig.prm`（参数）。

### 验收标准

打开 `lig.rtf`，第一行 RESI：

```
RESI lig  1.000 ! param penalty= ... ; charge penalty= ...
```

| 罚分 | 含义 | 处理 |
|---|---|---|
| < 10 | 类比合理 | 直接用 |
| 10–50 | 需要基本验证 | 至少做一次二面角 QM 对比 |
| **> 50** | **类比很差** | 必须优化，或明确记录并论证 |

**罚分 > 50 时，先查是不是 ① 的电荷问题**：如果分子有明确的阳离子中心，
但 `lig.rtf` 里所有 N 都是中性类型（`NG2S0`/`NG301`/`NG2R52`），
而净电荷又是 +1，那就是电荷被摊到了环上——回 ① 重做。

### 记录

把 `param penalty` 和 `charge penalty` 写进 RUNLOG，这是方法学的一部分。

---

## ③ 定浓度与盒子尺寸

**浓度决定了聚集尺寸，是跨体系比较的前提。**

### 公式

```
V(Å³) = n_molecules / (c_mol_per_L × 6.022×10⁻⁴)
a(Å)  = V^(1/3)                      # 立方盒
```

例：50 个分子、50 mM → V = 1.66×10⁶ Å³ → a = 118.4 Å

### 关键决定

- 所有体系**必须同浓度**（相对偏差 < 5%）
- 记录最终盒子边长，后面每一步都要核对

### 常见错误

Packmol 的 `inside box 0 0 0 120 120 120` 只规定**摆放范围**，
输出的 CRYST1 是**已摆放坐标的外接框**。小分子占 109 Å、大分子占 121 Å，
再用 VMD `-t 8` 外扩，浓度就差出 1.5 倍。

**正确做法**：把最终盒子钉死，Packmol 的摆放盒按水层厚度收缩。

---

## ④ Packmol 摆放

### 输入文件 `wa.inp`

```
tolerance 2.0
add_box_sides 1.2
seed 12345                                # 固定随机种子！否则无法复现
output mycase_com.pdb

structure mol.pdb
  number 50
  inside box 0. 0. 0. 102.4 102.4 102.4   # = 最终盒子 − 2×8 Å
end structure
```

### 运行

```bash
packmol < wa.inp
```

### 验收

- 输出文件存在，残差（objective function）收敛
- 分子数正确：`grep -c '^ATOM'` ÷ 每分子原子数 = 50

### 用 skill 生成

```powershell
& $py "$S\scaffold.py" --dir <新目录> --name mycase --packmol `
      --molecules 50 --target-conc-mM 50 --ligand-pdb mol.pdb
```

会按目标浓度反算盒子，写出 `wa.inp` 和 `solvate.tcl`。

---

## ⑤ psfgen 建 PSF

需要拓扑：`top_all36_cgenff.rtf` + `lig.rtf`。

```tcl
package require psfgen
topology top_all36_cgenff.rtf
topology lig.rtf
pdbalias residue LIG lig
segment S { pdb mycase_com.pdb }
coordpdb mycase_com.pdb S
guesscoord
writepsf mycase_proc.psf
writepdb mycase_proc.pdb
```

### 验收

```
NATOM = 50 × 每分子原子数
```

查 PSF 第 6 行的 `!NATOM` 数字。

---

## ⑥ 加水与加离子

### 加水（盒子钉死）

```tcl
package require solvate
solvate mycase_proc.psf mycase_proc.pdb \
    -minmax {{0.0 0.0 0.0} {118.4 118.4 118.4}} \
    -o mycase_solv
```

`-minmax` 而不是 `-t 8`：前者把盒子钉死到指定尺寸，后者按溶质外接框外扩，
体积会随分子大小漂移。

### 加离子

只加到**中和**为止：`n_ion = |净电荷 × 分子数|`。

### 验收（两条都要）

1. 体系总电荷 = 0
2. 记下离子的浓度 → 离子强度

**离子强度是个被忽略的量。** 50 个 +1 染料配 50 个 Cl⁻，
在 2.4×10⁶ Å³ 的盒子里就是 35 mM 的氯离子——会显著屏蔽静电作用。

如果实验是在 10 mM 缓冲液里、染料只有 µM 级，模拟的离子强度高了上千倍。
这时候"静电差异导致聚集不同"这类结论会被屏蔽效应淹没。

**在同浓度的前提下，离子强度与染料浓度绑定，无法单独降低**（除非让分子保持中性）。
报告时写明这一点。

---

## ⑦ 生成 NAMD 输入

```powershell
& $py "$S\scaffold.py" --dir <体系目录> --name mycase --prod-ns 100 --dcdfreq-ps 10
```

生成 `case.json` + 四个 conf。`case.json` 是**唯一真源**，改参数改它再 `--force` 重建。

### 必须一致的参数（跨体系）

| 参数 | 值 | 说明 |
|---|---|---|
| `oneFourScaling` | **1.0**（CHARMM36/CGenFF） | 写成 `1-4scaling` 会触发弃用警告 |
| `timestep` | 2.0 | |
| `cutoff` / `switchdist` | 12.0 / 10.0 | |
| `dcdfreq` | 统一 | 100 ps 和 10 ps 混用会让统计权重失衡 |

### 常见错误

| 现象 | 后果 |
|---|---|
| 平衡用 1.0、生产用 1.07 | 势能面在交接处跳变，速度对应旧哈密顿量 |
| 中途改 `dcdfreq` | 粗采样段权重缩水 |
| 生产阶段 `outputName` 与已有轨迹同名 | 追加写入，两段数据混在一个文件里 |

---

## ⑧ 跑前体检

```powershell
& $py "$S\preflight.py" --dir <体系目录>
```

### 验收

**FAIL 数量 = 0**。WARN 要逐条看。

检查内容：参数文件是否存在、PSF/PDB 与原子数、盒子 vs CRYST1、
跨 conf 的物理参数一致性、`firsttimestep` 与 `.xsc`、`outputName` 会不会追加、
磁盘够不够、`case.json` 有没有被手改。

### 如果 FAIL

按提示改 `case.json`，`--force` 重建，再跑一次预检。**不要手改 conf。**

---

## ⑨ 运行与断点续跑

```bash
namd2 mycase.conf > mycase.log 2>&1 &
```

### 启动后第一件事

看日志确认实际生效的物理参数：

```bash
grep "1-4 ELECTROSTATICS" mycase.log
grep "Configuration file is" mycase.log
```

### 中断后续跑

```powershell
& $py "$S\restart.py" --dir <目录> --from <中断的outputName> --target-ns 100 --write
```

它从 `.xsc` 读真实步数、改输出名、去掉 `cellBasisVector`、算剩余步数和磁盘。

**不要从一条旧轨迹的末尾直接换参数接着跑**——要在同一步、同参数下续。

---

## ⑩ 轨迹审计

```powershell
& $py "$S\inspect_traj.py" --dir <体系目录>
```

### 验收

- 帧数与 DCD 头的 NSET 一致（不一致说明被截断）
- 覆盖时长达到目标，**无 GAP**
- 日志正常 `End of program`（或知道为什么没有）
- 帧间隔记录在案

### 关键认知

**日志不能作为"跑完了"的证据。** 缓冲写入会让日志落后于文件；
`End of program` 也不说明模拟了多少时间。DCD 头才是准的。

---

## ⑪ 定义分析口径

```powershell
& $py "$S\inspect_molecule.py" --dir <体系目录> --write
```

读拓扑、找环、并成环系、算平面性，把可用的写进 `case.json`。

### 验收

**平面性 rms < 0.2 Å**。大了说明那个环系含 sp3 碳，不能当平面用。

输出示例：

```
group     atoms  aromatic  planarity rms  usable
ring1        16     12/16          0.097     yes
ring2        14     10/14          0.325      no      ← 含 sp3，剔除
ring3         6       6/6          0.001     yes
```

### 需要人工决定的

- `plane_pair`：默认最大环系与自身（同基团堆叠）。要比较给体 vs 受体就改。
- `rdf_pairs`：可以用原子名（`C5`）或**环质心**（`ring:ring1`）。
  环质心的好处是不依赖原子命名，换分子直接能用。

---

## ⑫ 分析

```powershell
& $py "$S\analyze.py" --dir <体系目录> --traj 04_prod.dcd --last-ns 50
```

### 输出的量

| 文件 | 内容 | 验收 |
|---|---|---|
| `summary.txt` | 全部摘要 | 先看这个 |
| `rdf_*.csv` | 径向分布 | 峰位与层间距自洽 |
| `planes_*_slip.csv` | 滑移角分布 | J 型占比 > 50% 才算"滑移堆积主导" |
| `planes_*_stacked_d.csv` | 堆叠层间距 | 中位数应在 3.2–4.0 Å |
| `planes_*_nearest.csv` | 最近邻中心距 | **必须 > 2 Å**，否则坐标系或原子选择错了 |
| `clusters_timeseries.csv` | 聚集演化 | 趋势应稳定，不是"一直长大" |

### 三条自检线

1. **最近邻中心距 > 2 Å** —— 小于这个物理上不可能
2. **层间距中位数 ≥ 3 Å** —— 小于 3 Å 通常是原子组穿插
3. **平面性 rms < 0.2 Å** —— 大了说明原子组选错

前两条脚本会自动报 `SANITY FAIL`。

### 判据必须写明

聚集尺寸依赖判据（接触/质心、cutoff）。同一份数据换个 cutoff，
"最大簇"可以从 20 个分子变成 50 个。论文里必须写清楚用了什么。

---

## ⑬ 可比性与结论闸门

```powershell
& $py "$S\compare_systems.py" --root <父目录>
& $py "$S\study_report.py"  --root <父目录> --analysis-root <分析输出目录>
```

### 可比性的五关

1. 浓度一致（< 5% 偏差）
2. 帧间隔一致
3. 生产时长一致
4. 力场旋钮一致（`oneFourScaling`/`cutoff`/`timestep`）
5. 分析口径一致（cutoff、阈值）

### 结论闸门

`study_report.py` 把每条结论判成 PASS / BLOCKED：

```
[PASS]    J-type slip stacking exists     all systems > 50% J-type
[PASS]    Interplanar distance <= 4 A     median d <= 4.0 A in every system
[BLOCKED] Assembly has equilibrated       still growing in ...
[BLOCKED] Substituent comparison          systems differ in concentration, spacing
```

**有 BLOCKED 就不要写跨体系排序。** 那是在比浓度和参数，不是比化学。

---

## 完全自动化：能到哪一步

```powershell
& $py "$S\pipeline.py" --dir <case> --name <name> --ligand-pdb mol.pdb `
      --molecules 50 --target-conc-mM 50 --prod-ns 100
```

一条命令跑完 ③ 定盒子 → ④ Packmol → ⑤ psfgen → ⑥ 加水加离子 → ⑦ 生成 conf → ⑧ 体检，
每步打印 `[OK]` / `[SKIP]` / `[FAIL]`，缺工具时说清楚缺什么而不是半途而废。
用 `--skip-pack` / `--skip-build` / `--skip-solvate` 可以从中间续。

实测一次完整运行：

```
[OK  ] ligand: mol.pdb (86 atoms), 50 copies at 50 mM
[OK  ] box: final 118.42 A (1,661 nm3), packing 102.42 A
[OK  ] psfgen: wrote demo_proc.psf
[OK  ] solvate + neutralize: 159,069 atoms, net charge -0.0000, ions {'CLA': 50}
[OK  ] NAMD .conf written
[OK  ] preflight passed
```

### 不能自动化的三步，和为什么

| 步骤 | 为什么必须人做 |
|---|---|
| ① 定电荷态 | **是化学事实，不是优化问题。** 需要知道实验用的是氯化物盐还是中性化合物、正电荷在哪个原子上。脚本能检查自洽性，但决定不了 |
| ② CGenFF 参数 | 需要 CHARMM-GUI / ParamChem **网页**，或本地的 CGenFF 程序。没有公开 API，无法脚本调用 |
| ⑨ 跑 MD | NAMD 在服务器上，本地通常没有 |

如果以后拿到了**本地 CGenFF 程序**（`cgenff` 命令行版），② 也能自动化；① 仍然需要你回答一次问题。

### 依赖

| 工具 | 用途 | 缺失时 |
|---|---|---|
| Python 3.10+ / numpy | 全部脚本 | 必装 |
| **Packmol** | 随机摆放 | `pip install packmol`（v100+ 自带二进制）；缺失时 `pipeline` 停在 ④，仍会写好 `wa.inp` |
| **VMD** | psfgen、solvate、autoionize | `pipeline` 停在 ⑤，会写好 tcl |
| NAMD | 跑 MD | 不在本地也行，把目录拷到集群 |

VMD 必须支持 `-dispdev text`（无界面）。验证：

```powershell
vmd -dispdev text -e check.tcl     # check.tcl 里 package require psfgen / solvate / autoionize
```

---

## 已有项目如何接入（不必重做）

手上已经跑完一批模拟、conf 是手写的，也能用这套工具，**不用重建任何文件**：

```powershell
# 只写 case.json，现有 .conf 一个都不动
& $py "$S\scaffold.py" --dir <已有体系目录> --case-only
& $py "$S\preflight.py" --dir <已有体系目录>
```

`--case-only` 只写 `case.json`（记录原子数、盒子、分析配置），
然后 `preflight` 就会检查跨 conf 一致性、`firsttimestep`、磁盘，
但**不会**把手写的 conf 判成"漂移"——只有本工具生成的 conf 才做漂移检查。

这样一批旧数据立刻能享受：

- `compare_systems.py` 的浓度/帧间隔/参数对照
- `study_report.py` 的结论闸门
- `analyze.py` 的自动环系识别

唯一做不了的是"改 case.json 重建 conf"——那需要 `--force`，会覆盖手写文件。
老项目要重建 conf 时，先备份。

---

## 其他自组装体系怎么用

这份教程以"芳香分子堆叠成 J-聚集体"为例，但**大部分内容与体系无关**。
换体系时，只有第 ⑪⑫ 步的分析口径要换。

### 哪些步骤完全不变

①②③④⑤⑥⑦⑧⑨⑩ —— 定电荷、生成参数、搭盒子、跑 MD、轨迹审计，对所有体系都一样。

### 第 ⑪⑫ 步怎么换

**用 `aggregate.py` 代替 `analyze.py`。** 它不需要任何分子知识，自动把
"不是水、不是离子"的残基都当溶质：

```powershell
& $py "$S\aggregate.py" --dir <体系目录> --traj 04_prod.dcd --last-ns 50
```

产出聚集尺寸分布、最大簇的 Rg、形状分类（球/棒/片）、多组分组成、径向密度。

### 不同体系看什么

| 体系 | 主要看 | 判读 |
|---|---|---|
| **胶束** | 聚集数分布、Rg、径向密度 | 密度在中心最高 = 实心核；聚集数分布应该窄 |
| **囊泡** | 径向密度 | 密度剖面在 r=0 附近**低**、中间高 = 空腔 |
| **纤维** | Rg、形状、各向异性 | rel_anisotropy 0.5–0.8 = 棒状；看它随时间是变长还是变粗 |
| **凝胶网络** | 簇数、最大簇 | 最大簇接近全部分子且持续增长 = 渗流 |
| **多组分共组装** | `composition.csv` | 组成随时间漂移 = 分相；稳定 = 共组装 |
| **π-堆叠发色团** | 另加 `analyze.py` | 见第 ⑫ 步的滑移角与层间距 |

### 形状描述符怎么读

`shape.csv` 里有三个量：

- **Rg**（回转半径）：聚集体有多大。判断紧密程度时看 Rg 与 N^(1/3) 是否同步增长。
- **rel_anisotropy**：接近 0 = 球形；0.5–0.8 = 棒状；接近 1 = 极度拉伸。
- **asphericity / acylindricity**：分别衡量"偏离球"和"偏离圆柱"。

### 判据仍然要写明

`--cutoff` 和 `--mode`（contact / centroid）会直接改变"最大簇"的数字。
同一份数据换个 cutoff，最大簇可以从 20 个分子变成 50 个。
**论文里必须写清楚用了哪个判据、哪个 cutoff。**

---

## 可复现性记录表

填完这张表，别人（或三个月后的你）才能重跑：

| 项目 | 记录什么 |
|---|---|
| 净电荷与带电位置 | +1 在哪个原子，依据是什么 |
| CGenFF 罚分 | param / charge 两个数，CGenFF 版本 |
| 分子数 / 目标浓度 | 如 50 / 50 mM |
| 最终盒子 | 实测 CRYST1 |
| 离子 | 种类、数量、离子强度 |
| Packmol seed / tolerance | 如 12345 / 2.0 |
| NAMD 版本与参数 | 见 RUNLOG |
| 生产时长 / 帧间隔 | 如 100 ns / 10 ps |
| 分析阈值 | slip、lateral、stack、cluster |
| 分析用的原子组 | 由 `inspect_molecule` 生成的那份 |

---

## 常见错误速查

| 现象 | 根因 | 查什么 |
|---|---|---|
| CGenFF 罚分 > 50 | 电荷态错 / 结构残缺 | mol2 电荷和、原子类型是否 SYBYL |
| 各体系浓度差 > 5% | Packmol 外接框 | 最终盒子 vs 目标 |
| 层间距 ~10 Å | 坐标系错误 | 最近邻距离是否 < 2 Å |
| 层间距 ~2 Å | 原子组穿插 | plane_groups 是否选错环 |
| 最大簇恒为全部分子 | 判据渗透 | 改用质心判据或更小 cutoff |
| 团聚趋势一直"growing" | 没平衡 | 加长时间，别比较尺寸 |
| 生产日志有 `No space left` | 磁盘 | `preflight` 的磁盘预估 |
| 轨迹帧数少于预期 | 崩溃/截断 | DCD 头 NSET vs 实际 |

---

## 附：阈值从哪来

本教程的默认阈值不是权威标准，而是**合理起点**：

| 参数 | 默认 | 依据 |
|---|---|---|
| 滑移角 J/H 分界 | 32° | 常用经验分界，30–35° 都出现过，设后要说明 |
| 横向偏移上限 | 6 Å | 分子横向尺寸的一半左右，按分子大小调 |
| 相邻中心距上限 | 8 Å | 层间距 3.5 Å + 横向 6 Å 的合成 |
| 接触判据 cutoff | 4 Å | 重原子间范德华接触 |
| 平面性上限 | 0.2 Å | 芳香环典型 rms < 0.1 Å，0.2 留余量 |

**换分子、换体系时要重新核对这些值是否合适**，并把最终采用的值写进方法。