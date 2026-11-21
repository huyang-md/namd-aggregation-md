# 操作清单

照着打勾。`👤` = 必须你自己做，`🤖` = 一条命令，`⛔` = 不过就停。

```powershell
$py = '<你的 python>'
$S  = '<skill 目录>\scripts'
```

---

## 阶段 1　分子参数　👤

- [ ] **确认盐形式**：实验用的是氯化物盐（阳离子）还是中性化合物？
- [ ] **确认电荷位置**：若带电，+ 在哪个原子上？（查论文/SI/合成路线）
- [ ] **ChemDraw 画出带电形式**，形式电荷显式标注
- [ ] **导出 3D mol2**，保留形式电荷
- [ ] ⛔ **mol2 电荷和 = 设定净电荷**（不是 0.000）
- [ ] **上传 CHARMM-GUI Ligand Reader**，填净电荷，下载 `lig.rtf` + `lig.prm`
- [ ] ⛔ **看罚分**：`RESI lig 1.000 ! param penalty= ... charge penalty= ...`
      - `< 10` 通过 ｜ `10–50` 需验证 ｜ `> 50` 回上一步查电荷
- [ ] 记录到 RUNLOG：净电荷、带电原子、CGenFF 版本、两个罚分

---

## 阶段 2　搭体系　👤 + 🤖

- [ ] 👤 决定：分子数（如 50）与目标浓度（如 50 mM），**所有体系必须相同**
- [ ] 👤 装 Packmol（免安装 exe，放 PATH 或记住路径）
- [ ] 🤖 一条命令：

```powershell
& $py "$S\pipeline.py" --dir <体系目录> --name <体系名> `
      --ligand-pdb mol.pdb --molecules 50 --target-conc-mM 50 `
      --prod-ns 100 --vmd "<vmd.exe 路径>"
```

- [ ] ⛔ 输出最后一行是 `preflight passed`（0 FAIL）
- [ ] ⛔ 体系净电荷 `0.0000`
- [ ] 记录：最终盒子边长、离子数与种类、离子强度

---

## 阶段 3　跑 MD　👤

- [ ] 目录拷到服务器
- [ ] 提交 `namd2 04_prod.conf > 04_prod.log`
- [ ] ⛔ 启动后确认生效参数：

```bash
grep "1-4 ELECTROSTATICS" 04_prod.log    # 必须是 1.0
grep "Configuration file is" 04_prod.log
```

- [ ] 中断时用 `restart.py --from <outputName> --target-ns <目标>` 续跑，不要手改 conf

---

## 阶段 4　分析　🤖

- [ ] 🤖 审计：`inspect_traj.py --dir <目录>` —— 帧数、覆盖、无 GAP
- [ ] 🤖 认分子：`inspect_molecule.py --dir <目录> --write`
- [ ] 🤖 分析：`analyze.py --dir <目录> --traj 04_prod.dcd --last-ns 50`
- [ ] ⛔ 人工看三个数（在 `summary.txt`）：

| 数 | 合格线 |
|---|---|
| 平面性 rms | < 0.2 Å |
| 最近邻中心距 | > 2 Å |
| 层间距中位数 | ≥ 3 Å |

---

## 阶段 5　下结论　👤

- [ ] 🤖 `compare_systems.py --root <父目录>`
- [ ] 🤖 `study_report.py --root <父目录> --analysis-root <分析目录>`
- [ ] ⛔ **没有 BLOCKED 才能写跨体系排序**
- [ ] ⛔ **聚集趋势必须 stable**，还在 growing 就不比尺寸
- [ ] 把最终阈值（滑移角、cutoff、接触判据）写进方法部分

---

## 交付前最后确认

- [ ] 所有体系的浓度偏差 < 5%
- [ ] 所有体系的帧间隔、时长、力场旋钮一致
- [ ] 每个体系的 `param penalty` 已记录
- [ ] 离子种类与浓度已记录
- [ ] 分析阈值已记录
- [ ] 可复现性记录表填完（见 `docs/tutorial.md` 同名章节）