# synth2 R Workflow Requirements

**Version:** 1.0
**Status:** Draft
**Date:** 2026-08-16

## Problem Frame

应用研究者用经典合成控制（Abadie–Diamond–Hainmueller）写一篇比较案例论文时，真正耗时的不是求解 \(W\)，而是把嵌套优化、空间安慰剂、时间安慰剂、leave-one-out、平衡表和论文图收成一份可辩护的附录。Stata 的 `synth2`（Yan & Chen, *Stata Journal* 2023）把这件事收成一条命令；R 里零件散在 `Synth`、`tidysynth`、`SCtools` 里，而且 R / Stata 的权重经常对不齐。

本项目的问题不是「R 没有合成控制」，而是：从 Stata / 陈强教材 / 连享会迁过来的人，无法用一个 R 入口得到与 `synth2` 同构、且在标准数据上数值可对齐的完整报告。gsynth、SDID、scpi 是另一条方法线，不在这个楔子里。

## Requirements

| ID | Requirement | Priority | Notes |
|----|-------------|----------|-------|
| R1 | 估计经典单处理单位、单处理时点的 Abadie SCM：\(W \ge 0,\ \sum W=1\)，对角 \(V \ge 0\)。 | Must Have | 与 `synth` / `synth2` 同一估计量，不是 gsynth / augsynth / SDID。 |
| R2 | 支持默认回归型 \(V\)、`nested` 嵌套优化、`allopt` 多起点；可传入 `customV`。 | Must Have | 对应 `synth2` 的 `nested` / `allopt` / `customV`。 |
| R3 | 预测变量时间窗口语法与 `synth2` 对齐：`xperiod`、`mspeperiod`、`var(t)`、`var(t1&t2&t3)`、整段预处理平均。 | Must Have | 复现 Prop 99 规格 `cigsale(1988) cigsale(1980) cigsale(1975)` 必须不变形。 |
| R4 | 长面板输入：单位、时间、结果、协变量；可选 donor 子集（`counit`）、`preperiod` / `postperiod`。 | Must Have | 等价于 `xtset` 后的 long data，不要求用户先 `dataprep()`。 |
| R5 | 空间安慰剂：对 donor 轮流当作处理单位；输出 pre/post MSPE、post/pre 比、Fisher \(p\)；支持 `cut` 过滤预处理 MSPE 过大的单位；输出逐年（左/双侧）\(p\)。 | Must Have | 对应 `placebo(unit cut(#))`。 |
| R6 | 时间安慰剂：指定一个假处理时点，重估并报告假处理窗与真处理窗的差距。 | Must Have | 对应 `placebo(period(t))`。 |
| R7 | Leave-one-out：每次去掉一个非零权重 donor，报告逐期效应的 min/max 带。 | Must Have | 对应 `loo`。 |
| R8 | 数值对齐：在锁定的金标准数据上，\(W\)、\(V\)、预处理 RMSE、逐期效应、ATT、MSPE 比、Fisher \(p\) 与 Stata `synth2` 对齐（见 Success Criteria）。 | Must Have | 本机无 Stata；金标准来自已发表表 + 外部 Stata 日志夹具。 |
| R9 | 报告对象包含 synth2 核心表：预处理拟合（RMSE、R²、N donor、N covariate）、平衡表（Treated / Synthetic / Average control 及 % bias）、\(V\) 权重、\(W\) 权重（含零权重单位可查询）。 | Must Have | 对应 `e(rmse)`、`e(r2)`、`e(bal)`、`e(V_wt)`、`e(U_wt)`。 |
| R10 | 默认产出与 `synth2` 同构的图：路径拟合、效应、\(V\) 条形图、\(W\) 条形图；若跑了安慰剂/LOO，再出 spaghetti、MSPE 比排序、逐年 \(p\)、LOO 带。 | Must Have | 对应 `savegraph` / `e(graph)`；图形设备用 ggplot2。 |
| R11 | 单一入口函数跑完「估计 →（可选）安慰剂 →（可选）LOO → 表/图/结果对象」，不要求用户手写循环。 | Must Have | 这是相对 `Synth` + `SCtools` 的产品差异。 |
| R12 | 结果对象可打印、可 `summary()`、可取出矩阵/表；安慰剂与 LOO 结果挂在同一对象上，而不是散落的全局副作用。 | Must Have | 对标 `e()`，但用 R 的 S3/S7，不用改全局环境。 |
| R13 | 混合安慰剂（Yan & Chen, *Economics Letters* 2023）：假单位 + 假时点可在同一次调用里组合。 | Should Have | `synth2` 已支持；v1 若时间紧可次于 R5/R6，但接口要预留。 |
| R14 | 优化调参与 Stata 插件对齐暴露：`margin`、`maxiter`、`sigf`、`bound`。 | Should Have | 没有这些，R8 很难调到金标准。 |
| R15 | 中英文 vignette：完整复现 smoking / Prop 99，对照已发表数字；中文 vignette 按连享会/`synth2` 教学顺序写。 | Should Have | 第一批用户是教材/课程迁移，不是 CRAN 猎奇。 |
| R16 | 将金标准夹具（Stata 日志或手工抄录的表）放进 `tests/testthat/`，CI 对 R8 做回归。 | Must Have | 没有 Stata 的机器也必须能测「还对齐不对齐」。 |
| R17 | 包名、函数名不与 CRAN 上 `seewave::synth2` 冲突；LICENSE 与依赖兼容（`synth2` SSC 为 GPL-3）。 | Must Have | 仓库可叫 synthAIO；CRAN 包名另定。 |
| R18 | 多处理单位、交错采纳、交互固定效应、矩阵补全、合成 DID、scpi 预测区间。 | Nice to Have | 明确推迟。见 Scope Boundaries。 |

## Success Criteria

在 **Abadie–Diamond–Hainmueller (2010) smoking** 数据、规格

`cigsale ~ lnincome + age15to24 + retprice + beer + cigsale(1988) + cigsale(1980) + cigsale(1975)`，
`trunit = California`，`trperiod = 1989`，`xperiod = 1980:1988`，`nested + allopt`

上，相对 Stata `synth2`（或 Yan & Chen 2023 / 已核对的教学复现）满足：

1. 非零 \(W\) 的单位集合一致；每个非零权重绝对误差 \(\le 10^{-3}\)（或相对误差 \(\le 1\%\)，取较松者，并在夹具里写死实际容差）。
2. 对角 \(V\) 的主导两项（通常是 `age15to24` 与 `cigsale(1975)`）排序一致，数值绝对误差 \(\le 10^{-3}\)。
3. 预处理 RMSE 相对误差 \(\le 1\%\)；平均 ATT 绝对误差 \(\le 0.05\) 包/人（smoking 量纲）。
4. 空间安慰剂：California 的 post/pre MSPE 比排名第一；未过滤 Fisher \(p\) 与 \(1/N_{\text{units}}\) 一致（smoking 上约为 \(1/39 = 0.0256\)）。
5. 一次调用产出 R9 的全部表和 R10 的核心图，无需用户再写循环。
6. `R CMD check` 在关闭建议依赖的情况下无 ERROR / WARNING。

未取得作者/用户提供的 Stata 日志前，以 Carlos Mendez 教学复现与 SJ 论文已发表数字为临时金标准，并在夹具头注释来源与版本（`synth2` 2.1.0，2023-10-05）。拿到 Stata 日志后替换为更高精度夹具，容差只许收紧不许放宽。

## Scope Boundaries

**In scope:**

- 经典单处理 SCM 的估计、嵌套优化、空间/时间安慰剂、LOO、平衡表、论文图、结果对象、金标准测试。
- 与 `synth2` 对等的调参与预测变量窗口语法。
- 教学级 Prop 99 复现（中英文）。

**Out of scope:**

- gsynth（交互固定效应 / 多处理 / 交错）。
- `synthdid`、`augsynth`、`scpi`、`microsynth`、`MSCMT` 多元结果作为用户可见估计量。
- 多个处理单位的联合估计（用户若先加总再当一个单位，那是数据预处理，不是本包功能）。
- 网页应用、Shiny、Stata 互操作桥（读 `.dta` 可以，调 Stata 不可以）。
- 把现有 R 包再导出一层无对齐保证的 wrapper。

## Key Decisions

| Decision | Chosen | Rationale | Alternatives Considered |
|----------|--------|-----------|-------------------------|
| 产品形态 | 经典 SCM 的 synth2 级工作流包，仓库可仍叫 synthAIO | 用户选定楔子 A；全家桶是海 | B 统一 API；C 只写手册 |
| 验收 | 数值对齐 Stata `synth2`，不是「看起来像」 | 否则复现教材时无法辩护 | 仅工作流对齐 |
| 估计量 | 只做 Abadie 2010 SCM | gsynth 是另一篇论文、另一组假设 | 顺手接入 gsynth |
| 求解器 | 为实现 R8，允许自研或深度封装二次规划；禁止「调用 `tidysynth` 就算完成」 | R/Stata `synth` 权重已知不一致 | 薄包装 `Synth` / `tidysynth` |
| 金标准 | 夹具入库；本机无 Stata 不阻塞开发 | 已发表表可启动；Stata 日志后替换 | 开发机必须装 Stata |
| 用户接口 | 一个主函数 + 长数据 + 单位/时间/处理标识；预测变量窗口用公式或专用 helper | Stata 迁移成本最低 | 只提供 tidy pipe；只提供 `dataprep` 列表 |
| 包名 | 不叫 `synth2` | CRAN 上 `seewave::synth2` 已占用 | 强行 `synth2` |
| 许可 | GPL-3 或与依赖兼容的更严 copyleft | 上游 SSC `synth2` 为 GPL-3；若只重实现方法、不抄 ado，也可 GPL-2+/MIT，但需法律复核后再定 | 未复核就 MIT |

## Outstanding Questions

| # | Question | Impact if Wrong | Owner |
|---|----------|-----------------|-------|
| Q1 | 能否拿到 smoking 规格下 `synth2, nested allopt` 的完整 Stata 日志（\(W,V\)、逐期效应、安慰剂矩阵）？ | 夹具精度不够，R8 会和真正的 Stata 用户对不上 | 用户 |
| Q2 | CRAN 包名用 `synthaio`、`synth2r` 还是其他？ | 改名成本在提交前低、提交后高 | 用户 |
| Q3 | 求解器自研（osqp / quadprog / 自写内点）还是先试 `MSCMT`/`Synth` 再对残差？ | 选错会浪费一周仍对不齐 | 实现阶段；默认：先对照，对不齐再自研内层 QP |
| Q4 | v1 是否必须交出 R13 混合安慰剂，还是接口预留即可？ | 多一周实现与测试 | 用户；文档默认 Should Have |
| Q5 | 目标 R 版本下限（建议 ≥ 4.2）与是否进 CRAN（相对 GitHub-only）？ | 影响 check 严格度和依赖选择 | 用户；文档默认按 CRAN 质量写，发布渠道后定 |
