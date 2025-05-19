# PhD Thesis in Multivariate Conditional Density Estimation with Likelihood Inference and Regression Analysis: Comparing Transformation Forest Models, Copulas and Normalizing Flows

This repository collects the R notebooks and scripts used in the associated thesis.  All experiments are run directly in R to study multivariate conditional densities with mathematical rigour.

**Important:** Read `AGENTS.md` before running any scripts or using the assistant. It explains how to keep this README in mind when analyzing outputs or generating code.

## Table of contents
1. [Quickstart](#quickstart)
2. [Scientific approach](#scientific-approach)
3. [Implementation](#implementation)
4. [Theory](Theory.md)
5. [Notation](Notation.md)

## Quickstart
Execute scripts with `Rscript` and run the sanity checks via

```sh
./run_checks.sh
```

Refer to `AGENTS.md` for further guidance.

## Scientific approach
This repository is a research notebook rather than a software project. All scripts are executed directly in R to study
problems in mathematical statistics and probability theory. The focus is on
precise and reproducible experiments with **Mathematical Rigour**, not on
standard software engineering workflows or packaging. The choice of R reflects
the statistical setting; Python is intentionally avoided.

The triangular transport methodology is described in the preprint
**arXiv:2503.21673v1 [stat.CO], 27 Mar 2025**.
<https://arxiv.org/abs/2503.21673>

Our \(S_k(x_1,\ldots,x_k)\) functions are simply transformation forests for
the regression \(x_k \sim x_1,\ldots,x_{k-1}\). Monotonicity in \(x_k\) is
therefore automatically ensured.



## Theory
The theoretical background has been moved to [Theory.md](Theory.md) for clarity.
## Implementation
* **Triangular transport map:** Implement by choosing distribution \$d\_k\$ and functional form for \$\theta\_k\$ for each dimension. The algorithm computes \$S^{-1}\$ and \$\ell(x)\$ by iterating \$k=1\$ to \$K\$. *Sampling:* draw \$z \sim \mathcal{N}(0,I\_K)\$ and set \$x\_k = q\_{d\_k}(\Phi(z\_k) \mid x\_{\<k}, x\_{\text{cov}})\$ sequentially to get a sample \$x\$. *Density evaluation:* given \$x\$, compute \$z=S(x)\$ and accumulate \$\log f\_{d\_k}(x\_k \mid x\_{\<k},x\_{\text{cov}}) - \log \varphi(z\_k)\$ over \$k\$. This yields exact likelihood values and allows fast simulation from \$\pi\$.
* **Transformation forest approach:** Instead of parametric \$\theta\_k\$, fit each conditional distribution with a nonparametric forest. For each \$k\$, treat \$Y\_k\$ as response and \$(Y\_{\<k}, X)\$ as covariates in a transformation forest. This produces an empirical CDF \$\hat F\_k(y\_k \mid y\_{\<k},x)\$ that is monotonic in \$y\_k\$. Use \$\hat F\_k\$ (and its inverse) in place of \$F\_{d\_k}\$ in the triangular map to obtain a fully data-driven \$S\$. The forest automatically captures interactions via tree splits, without explicit basis functions.
- Let \(\hat F_k(x_k\mid x_{<k})\) be the conditional CDF estimated by a transformation forest fitting \(X_k\) on \(X_{<k}\).  Then
  \[
    S_k(x_1,\dots,x_k)
      \;=\;\Phi^{-1}\bigl(\hat F_k(x_k\mid x_{<k})\bigr).
  \]

- For each fixed \(x_{<k}\), \(\hat F_k(\cdot\mid x_{<k})\) is non-decreasing in \(x_k\) by construction of the forest-based CDF estimator; hence
  \[
    \partial_{x_k}S_k(x_1,\dots,x_k)
      = \frac{\partial_{x_k}\hat F_k(x_k\mid x_{<k})}{\varphi\bigl(\Phi^{-1}(\hat F_k(x_k\mid x_{<k}))\bigr)}
      \;\ge\;0
  \]
  almost everywhere, so each \(S_k\) is monotone in its last argument and \(S\) is invertible :contentReference[oaicite:0]{index=0}:contentReference[oaicite:1]{index=1}.

* **Copula method:** Two-stage fitting: (1) Estimate each marginal CDF \$\hat F\_j(y\_j \mid x)\$ for \$j=1,\dots,d\$ (e.g. using parametric models or empirical ranks). (2) Transform outputs to \$u\_{ij} = \hat F\_j(y\_{ij}\mid x\_i)\$ and fit a copula \$C(u\_{1},\dots,u\_{d})\$ on $\[0,1]^d\$ (e.g. a vine copula). Then approximate the joint CDF by
  $F(y_1,\dots,y_d \mid x) \approx C\!\big(\hat F_1(y_1|x),\,\dots,\,\hat F_d(y_d|x)\big).$
  This approach separates marginal and dependence modeling. If the dependence between \$Y\$ components changes with \$x\$, a static copula may be inadequate without allowing copula parameters to vary with \$x\$.
  * **Conditional normalizing flow:** Use an invertible neural network to model \$p(y|x)\$. For example, implement a coupling-layer or autoregressive flow with parameters conditioned on \$X\$. The network \$T\$ is trained such that \$z = T(y; x) \sim \mathcal{N}(0,I)\$. Training maximizes \$\sum\_i \log p\_{T}(y\_i \mid x\_i)\$ by gradient descent. This black-box model provides great flexibility (universal function approximation) at the cost of many parameters and less interpretability. It serves as a baseline for comparison, highlighting the trade-off between flexibility and interpretability.


# Workflow
The interaction of the main scripts is summarised in
[workflow.md](workflow.md). It visualises the sequence executed by
`run_all.R` and highlights which objects are passed between scripts.

# Notation
The full notation tables are available in [Notation.md](Notation.md).
