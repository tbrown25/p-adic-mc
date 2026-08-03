# p-adic-mc — a p-adic analog of Brownian noise & Monte Carlo, in C

**A sketch / design doc.** A C project (novice → intermediate) that builds
random walks and diffusion on the **p-adic numbers** — the non-Archimedean analog
of Brownian motion and Monte Carlo — and visualizes them. Chosen because it has
the most *research* upside of the C ideas: it sits directly on the number-theory
program (P10) and connects to probability, physics of disordered systems, and the
QBism/foundations thread.

## The idea, in one breath

Over the reals, "close" is Archimedean and Brownian motion has Gaussian
increments. Over **Qₚ** the metric is **ultrametric** — the strong triangle
inequality `|x+y|ₚ ≤ max(|x|ₚ, |y|ₚ)` — and its geometry is a **tree** (the tree
of Zₚ-balls / the Bruhat–Tits tree). Diffusion there is governed not by the usual
Laplacian but by the **Vladimirov operator** `Dᵅ` (a pseudo-differential
operator), and the resulting random walk is **hierarchical**: a particle escapes
small balls fast and large balls slowly, giving **anomalous, log-periodic**
diffusion. That structure is exactly why p-adic models describe spin glasses and
protein energy landscapes (Parisi ultrametricity; Avetisov et al.).

So: **p-adic noise = jumps weighted by ultrametric (tree) distance**, and
**p-adic Monte Carlo = sampling the Haar measure on Zₚ** to estimate p-adic
integrals. Both are cheap, discrete, and tree-shaped — a natural fit for C.

## What gets built (C)

- **`padic.h/.c`** — a p-adic number to precision `k`: digit array base `p`
  (`x = Σ aᵢ pⁱ`, `0 ≤ aᵢ < p`). Operations: `add`, `mul` (with carry mod pᵏ),
  `valuation` (vₚ), `abs` (`|x|ₚ = p^{-vₚ(x)}`), `dist` (`|x−y|ₚ`). Unit tests
  (Hensel-style sanity: `1/(1−p)` = `1 + p + p² + …`, etc.).
- **`tree.c`** — the ball/tree view: a node is a residue mod pⁿ; the ultrametric
  distance is the depth of the last common ancestor. One **ultrametric random
  walk**: from the current ball, jump to another with probability decaying in tree
  distance (transition kernel `∝ p^{-(α+1)·dist}`), the discrete Vladimirov kernel.
- **`diffusion.c`** — an ensemble of walkers → **MSD(t)** in the p-adic metric;
  show the anomalous / **log-periodic oscillation** signature (the tell that the
  geometry is a tree, not a line). Optionally integrate the p-adic heat equation
  `∂ₜu + Dᵅu = 0` on the truncated tree and watch the heat kernel spread.
- **`mc.c`** — **Monte Carlo over Zₚ** (Haar measure = uniform digits): estimate
  a p-adic integral (e.g. `∫_{Zₚ} |x|ₚ^s dx = (1−p⁻¹)/(1−p^{−s−1})`) and validate
  against the closed form — the p-adic analog of "sample uniformly, average."
- **`viz`** — emit the tree + a highlighted walk, and MSD(t), as SVG/PPM from C
  (ties the "data-viz in C" thread), or export coordinates for the Backporch site
  (a `/pillars` or `/proofs` demo: the tree, the walk, the log-periodic MSD).

## Milestones (each a shippable unit)

1. Qₚ arithmetic in C + tests (`+`, `×`, vₚ, `|·|ₚ`, distance).
2. The tree + one ultrametric random walk; dump the trajectory.
3. Ensemble → MSD(t); reproduce the log-periodic anomalous-diffusion curve.
4. Monte Carlo integral over Zₚ; match the closed form to N⁻¹ᐟ² error.
5. Visualization (SVG from C, or exported to the site).

## Why it's *research*, not just a portfolio piece

- **P10 (number theory).** The p-adic world is where p-adic L-functions, p-adic
  modular forms, and Iwasawa theory live — this is a computational, visual entry
  to it, and the natural home for the artifact on the site.
- **Probability.** A concrete non-Archimedean Wiener/Monte-Carlo analog, with a
  falsifiable signature (log-periodicity) to check the code against theory.
- **Physics / QBism.** Ultrametric diffusion models disordered systems; p-adic
  QM (Vladimirov–Volovich) sits near the foundations thread.
- **Open wedge (honest):** is there a computable "p-adic noise" clean enough to
  drive a p-adic Monte Carlo for number-theoretic quantities (p-adic L-values,
  measures)? Exploratory — but the discrete tree makes it *testable*.

## References

- Kochubei, *Pseudo-Differential Equations and Stochastics over Non-Archimedean
  Fields*.
- Vladimirov, Volovich, Zelenov, *p-adic Analysis and Mathematical Physics*.
- Avetisov, Bikulov, Kozyrev, *Application of p-adic analysis to models of
  ultrametric diffusion* (protein dynamics, spin glasses).
- Parisi, Sourlas — ultrametricity in disordered systems.

## Status

**Milestone 1 — done & verified.** `include/padic/padic.h` + `src/padic.c` +
`tests/test_padic.c`, dependency-free C, `make test` (15/15 pass): Z/pᵏ arithmetic
(add/sub/mul with carry, `neg` as base-p complement), valuation, `|·|ₚ`, distance,
the ultrametric strong-triangle inequality, and the Hensel headline
`(1−p)·(1+p+p²+…) ≡ 1`. Numbers are base-p digit strings (little-endian).

**Next — milestone 2:** the ball/tree view + one ultrametric random walk
(`tree.c`): tree distance = depth of the last common ancestor; jump kernel
`∝ p^{−(α+1)·dist}`.
