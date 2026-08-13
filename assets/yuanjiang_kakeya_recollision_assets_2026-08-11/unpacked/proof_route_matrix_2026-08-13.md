# 圆酱四维挂谷—重碰撞母猜想：证明路线穷举矩阵

更新时间：2026-08-13。这里的“成立”只表示该局部命题已有证明；不表示四维挂谷猜想已解决。

| 路线 | 想要的桥 | 当前最硬结果 | 守门反例／缺口 | 状态 |
|---|---|---|---|---|
| 1. worldline 直译 | 三维自由飞行变成四维直线 | `(t,x(t))=(t,b+tv)` 完全严格 | 真实硬球 world-tubes 在同一时刻不重叠；挂谷需要高重叠 | 桥成立，朴素等同失败 |
| 2. 把管当硬球扫掠 | 管交就是碰撞 | 两轨迹近遇条件可精确写成横向投影小 | `ε≈δ` 时排斥消掉重叠；`ε≪δ` 时多数管交只是 ghost encounter | 单尺度方案否定 |
| 3. 偶环／C12 | 高 incidence 强迫闭环 | Naor--Verstraëte 对所有 `s<4` 强迫 `(1-o(1))Q^4` 个不同 C12；可抽 `Ω(Q^3)` vertex-disjoint C12 | 单 C12 有显式一般位置反例，不强迫 plane/bush/chord | 组合轮子完成，几何兑换开放 |
| 4. C12 拼刚性块 | 多环共享路径杀掉滑动核 | 四维 `C6` event polygon 恰多 1 flex；`Theta(3,3,3)` generic parallel rigid | 很多 C12 不强迫共享半圈；广义六边形型图可躲 | 新局部证书，全球计数不够 |
| 5. quadric net | 每个六边形落在二次曲面族 | P4 中任意六边形线并至少有 3 维 quadric net | 每环的 net 可完全不同；共享 ≤4 边通常无共同 quadric | 局部真，同步开放 |
| 6. exact hubs + polynomial | 把每格压成真共点，再用插值 | exact point--line 版本由 Zhang/多项式法排除 `s<4` | `δ`-近点不是真零点；Chebyshev 与 cube-center 例否定任何无条件固定幂稳定化 | 轮子清楚，连续缺口正是 Kakeya |
| 7. bearing／parallel rigidity | 近遇图谱隙迫使全局 bush | 若 normalized bearing gap 大，则 roots 近 `a+αv`，tube union 有常数体积 | gap 小时一般 matrix-Cheeger 无反向定理 | 严格条件定理，soft-mode 分类开放 |
| 8. incidence sheaf | 同一 coboundary 两种 Schur 消元 | 先消 cubes 得 Wang bearing；先消 tubes 得 DHM event paths；`r_cl+H1=2m` | 现成 sheaf 谱论不处理这张非正则 worldline sheaf | 母猜想最精确共同框架 |
| 9. exact central swaps | 碰撞交换速度，让 lanes 继续直走 | 若 swap graph parallel rigid，全部碰撞必同一时空 bush | 非 bush 的严格时序 swap 网络必须有 bearing flex；exact swap 是零测 | 严格 no-go／二分 |
| 10. nontrivial collision gates | 出射接到另外两根 tubes | 四速度必须是 velocity rectangle，即抛物面加法四元组 | 单 cube 可拆成互不相干 rectangles；需 tube-consistent weighted inverse theorem | 新接口，未闭合 |
| 11. BG 随机细化 | 每个粗近遇以 `(ε/δ)^2=Q^-1` 成为微碰撞 | 角度均匀时 ghost collisions 期望 `ND`; 一般角能量重现 `s<13/4` | 真碰撞散射会删未来 ghost encounters；局部窄方向帽正是 planebrush 难点 | 解释已知阈值，未到 4 |
| 12. cycle small-ball | 每个 fresh closure 付概率幂 | fixed graph：增益指数等于 quantitative closure-matrix rank；DHM Prop.7.11 给合法 molecule `ε^{cρ}` | 环数不等于 fresh rank；tube law 没有独立 quotient entropy | 严格 toy lemma + DHM 接口 |
| 13. relative entropy | lane law 编码便宜，重碰撞在 smooth gas 中极贵 | lane boxes 相对熵 `O(N log Q)`；若 global tail为 `exp(-cR logQ)` 且 `R≈ND`，所有 `s<4` 矛盾 | 尚无 lane→真实 `ND` recollisions 编译；DHM 无现成 global N-particle tail | 最短条件证明之一 |
| 14. collision holonomy／Möbius | 环上反射积的非平凡扭转收费 | 固定 history 可在 reflected sheets 展直；六球证书有非平凡 quotient holonomy | tree trivialization、变动法向与 Jacobian 不能丢；holonomy 不自动给体积 | 严格局部结构，全球量化开放 |
| 15. velocity paraboloid／decoupling | 能量动量守恒连接 Fourier/Kakeya | 碰撞矩形就是抛物面加法；与 Bourgain--Demeter 框架同型 | 这是现有 restriction/decoupling 大轮子，不是绕过 Kakeya 的捷径 | 查重完成 |
| 16. soft potential／virial | 用排斥势的时空作用量统计所有重叠 | 若轨迹保持贴管，virial 预算会给 `ND≲N` | 保持贴管要求耦合 `≲1/D`，正好消掉增益；力会删未来重叠 | 条件端点，裸路线抵消 |
| 17. projective/Möbius 放大角度 | 把局部窄方向帽吹开再碰撞 | projective/parabolic rescaling可把帽变宽 | tube 厚度、时间窗和 Jacobian 同时缩放，碰撞概率不变；硬球律也不具一般 projective invariance | 无免费增益，回到 induction on scales |
| 18. finite-field/generalized hexagon | 构造或排除高围长坏图 | 抽象尺度匹配反例覆盖 `11/3≤s<4`；exact real realization被多项式法排除 | `δ`-近欧氏实现若存在就是 Kakeya 反例 | 精确说明必须用欧氏连续几何 |

## 当前只剩的两个核心接口

1. **Worldline-specific Hodge/grain theorem**：对 tube--cube incidence sheaf，证明 near-section 要么给足 fresh coexact rank，要么低频 sections 同步成 plany/sticky grains。
2. **Collision compilation/cost + global tail**：把坏管族变成 `R≈ND` 个可收费的真实／条件碰撞 gates；或证明失败必落入前一个 carrier 分支。随后把 DHM 的 molecule 小因子提升为所需的 global large-deviation bound。

这两个接口若都无附加假设地完成，基本就是完整四维 Kakeya 证明的主体；不能把它们包装成“技术细节”。
