# PhD Thesis in Multivariate Conditional Density Estimation with Likelihood Inference and Regression Analysis: Comparing Transformation Forest Models, Copulas and Normalizing Flows

**Important:** Read `AGENTS.md` before running any scripts or using the assistant. It explains how to keep this README in mind when analyzing outputs or generating code.
## Motivation
* **Goal:** Estimate a flexible conditional density \$p(y \mid x)\$ for multivariate continuous outcomes \$Y \in \mathbb{R}^d\$ given features \$X \in \mathbb{R}^p\$, with full likelihood inference (not just point predictions).
* **Context:** Distributional regression - interest in the entire outcome distribution (beyond mean or variance). Relevant for heteroscedastic or multi-modal responses.
* **Challenges:** High-dimensional \$Y\$ with complex interdependencies; \$p(y|x)\$ may be multi-modal, non-Gaussian, and vary strongly with \$x\$. Needs models that are both flexible (to fit complex distributions) and tractable (for likelihood evaluation and sampling).
* **Approaches:** Consider both statistical and ML methods:

  * *Structured invertible transformation* (triangular transport map) for exact likelihood modeling.
  * *Copulas* for modular dependency modeling (with separate marginal fits).
  * *Transformation forests* for nonparametric conditional distribution estimation.
  * *Normalizing flows* (neural invertible networks) for maximum flexibility.

## Scientific approach
This repository is a research notebook rather than a software project. We will
never build an R package. All scripts are executed directly in R to study
problems in mathematical statistics and probability theory. The focus is on
precise and reproducible experiments with **Mathematical Rigour**, not on
standard software engineering workflows or packaging. The choice of R reflects
the statistical setting; Python is intentionally avoided.

## Theory
* **Triangular factorization:** Model a joint or conditional density via sequential univariate factors. For target variables \$X\_1,\dots,X\_K\$ (components of \$Y\$) conditioned on external covariates \$x\_{\text{cov}}\$:
  $\pi(x_1,\dots,x_K \mid x_{\text{cov}}) = \prod_{k=1}^K f_{d_k}\!\big(x_k \mid x_{1:k-1},\,x_{\text{cov}}\big),$
  where each \$f\_{d\_k}\$ is a chosen univariate density family with parameters \$\theta\_k(x\_{\<k},x\_{\text{cov}})\$.
* **Transport map \$S\$:** Define each conditional CDF \$F\_{d\_k}(x\_k \mid x\_{\<k},x\_{\text{cov}})\$. The *transport map* \$S: \mathbb{R}^K \to \mathbb{R}^K\$ sends \$x=(x\_1,\dots,x\_K)\$ to latent \$z=(z\_1,\dots,z\_K)\$ by
  $z_k = S_k(x_k \mid x_{<k},x_{\text{cov}}) = \Phi^{-1}\!\Big( F_{d_k}\!\big(x_k \mid x_{<k},x_{\text{cov}}\big) \Big),$
  where \$\Phi^{-1}\$ is the standard normal quantile. Each \$S\_k\$ is monotonic in \$x\_k\$, making \$S\$ invertible (triangular Jacobian).
* **Inverse transform:** Invert \$S\$ sequentially for sampling or inference. Given \$z \sim \mathcal{N}(0,I\_K)\$, set \$u\_k = \Phi(z\_k)\$ and then
  $x_k = F_{d_k}^{-1}(u_k \mid x_{<k}, x_{\text{cov}}), \qquad k=1,\dots,K,$
  to obtain \$x \sim \pi(\cdot \mid x\_{\text{cov}})\$. Conversely, any \$x\$ maps to \$z=S(x)\$.
* **Likelihood formula:** By change of variables,
  $\pi(x \mid x_{\text{cov}}) = \eta\!\big(S(x)\big)\;\prod_{k=1}^K \frac{f_{d_k}(x_k \mid x_{<k},x_{\text{cov}})}{\varphi(z_k)},$
  with \$\eta(z)=\prod\_{k=1}^K \varphi(z\_k)\$ the standard normal density and \$\varphi\$ the \$\mathcal{N}(0,1)\$ PDF. The log-likelihood for \$x\$ is
  $\ell(x) = -\tfrac{1}{2}\|S(x)\|^2 + \sum_{k=1}^K \Big[\log f_{d_k}\!\big(x_k \mid x_{<k},x_{\text{cov}}\big) - \log \varphi(z_k)\Big].$
  Maximize \$\sum\_i \ell(x^{(i)})\$ to estimate the functions \$\theta\_k\$ (maximum likelihood estimation).
* **Sparsity:** Impose conditional independence by restricting each \$\theta\_k\$ to a subset of \$x\_{\<k}\$. If \$j \notin J\_k \subseteq {1,\dots,k-1}\$, then \$X\_j\$ is excluded from \$f\_{d\_k}(x\_k|\cdot)\$ (implying \$X\_k \perp X\_j \mid {X\_i: i\in J\_k},,x\_{\text{cov}}\$ in the model). Sparse dependencies reduce complexity.
* **Flexible parameterization:** Specify \$\theta\_k(x\_{\<k},x\_{\text{cov}})\$ in a rich function class. For example, use a basis expansion
  $\theta_k(x_{<k},x_{\text{cov}}) \approx \sum_{m} c_{mk}\,\psi_{mk}(x_{<k},x_{\text{cov}}),$
  with \${\psi\_{mk}}\$ basis functions (polynomials, splines, tree indicators, etc). This captures nonlinear effects. The monotonic CDF link ensures the conditional density is valid.
* **Learning algorithm:** Train by maximizing log-likelihood. Use gradient-based optimization or boosting. A *map adaptation* strategy can start with an identity map and add basis terms to \$\theta\_k\$ iteratively (greedy selection by largest log-lik improvement). Triangular structure permits sequential or block-wise training of components.
* **Unifying cases:** This framework recovers many known models as special cases:

  * All \$d\_k\$ Gaussian with linear \$\theta\_k\$ \$;\to;\$ Gaussian copula (multivariate normal with regression).
  * Tree-based \$\theta\_k\$ \$;\to;\$ transformation forest models for each conditional.
  * Neural network \$\theta\_k\$ \$;\to;\$ autoregressive normalizing flow.

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
  almost everywhere, so each \(S_k\) is monotone in its last argument and \(S\) is invertible.

* **Copula method:** Two-stage fitting: (1) Estimate each marginal CDF \$\hat F\_j(y\_j \mid x)\$ for \$j=1,\dots,d\$ (e.g. using parametric models or empirical ranks). (2) Transform outputs to \$u\_{ij} = \hat F\_j(y\_{ij}\mid x\_i)\$ and fit a copula \$C(u\_{1},\dots,u\_{d})\$ on $\[0,1]^d\$ (e.g. a vine copula). Then approximate the joint CDF by
  $F(y_1,\dots,y_d \mid x) \approx C\!\big(\hat F_1(y_1|x),\,\dots,\,\hat F_d(y_d|x)\big).$
  This approach separates marginal and dependence modeling. If the dependence between \$Y\$ components changes with \$x\$, a static copula may be inadequate without allowing copula parameters to vary with \$x\$.
  * **Conditional normalizing flow:** Use an invertible neural network to model \$p(y|x)\$. For example, implement a coupling-layer or autoregressive flow with parameters conditioned on \$X\$. The network \$T\$ is trained such that \$z = T(y; x) \sim \mathcal{N}(0,I)\$. Training maximizes \$\sum\_i \log p\_{T}(y\_i \mid x\_i)\$ by gradient descent. This black-box model provides great flexibility (universal function approximation) at the cost of many parameters and less interpretability. It serves as a baseline for comparison, highlighting the trade-off between flexibility and interpretability.
  
  
# Notation 

| Symbol | Meaning (paper context) | Recommended R naming |
|--------|-------------------------|----------------------|
| **S** | Target->Reference map $S:\mathbb{R}^K \rightarrow \mathbb{R}^K$ | `S` / `S_map()` |
| $S_k$ | $k$-th map component, monotone in $x_k$ | `S_k()` / `S_list[[k]]` |
| $S^{-1}$ | Inverse map (Reference->Target) | `S_inv()` |
| $S^{-1}_k$ | Sequential inversion of the $k$-th component | `S_inv_k()` |
| $S^{\\sharp}\\pi$ | Push-forward distribution (image of $\\pi$ under $S$) | `push_fwd` |
| $S_{\\sharp}\\eta$ | Pull-back distribution (pre-image of $\\eta$ under $S^{-1}$) | `pull_back` |
| **R** | Reference->Target map (Sec. 3.2.2) | `R` / `R_map()` |
| $\\nabla S$ | Jacobian matrix of the map | `J_S` |
| $\\det\\nabla S$ | Determinant (triangular -> product of diagonals) | `det_J` |
| $\\partial_{x_k}S_k$ | Diagonal Jacobian entry | `dSdx_k` |
| **$\\pi$** | Target distribution (orange variables) | `pi_dist` |
| $\\pi(x)$ | Target density function | `pi_density` |
| **$\\eta$** | Reference distribution (green vars), usually $\\mathcal N(0,I)$ | `eta_dist` |
| $\\eta(z)$ | Reference density function | `eta_density` |
| **x** | Random vector $x \\sim \\pi$ | `X` (matrix N x K) |
| **z** | Random vector $z \\sim \\eta$ | `Z` |
| $X_i$ | $i$-th realisation of $x$ | `X[i, ]` |
| $Z_i$ | $i$-th realisation of $z$ | `Z[i, ]` |
| orange variable | Any quantity in the $\\pi$ domain | suffix `_pi` |
| green variable | Any quantity in the $\\eta$ domain | suffix `_eta` |
| $K$ | Number of target dimensions | `K` |
| $N$ | Ensemble / sample size | `N` |
| $p(\\cdot)$ | Generic pdf | `p_fun` |
| $a, b$ | Placeholder random variables | `a`, `b` |
| $y^*$ | Conditioning variable | `y_star` |
| $x^*$ | Conditional target variable | `x_star` |
| $L$ | Likelihood $\\eta(S(x))\\,|\\det\\nabla S|$ | `likelihood` |
| $\\ell$ | Log-likelihood | `loglik` |
| $\\|S(x)\\|^2$ | Squared norm of map outputs | `S_norm2` |
| $\\mathrm{KL}(\\pi\\,\\|\\,S_{\\sharp}\\eta)$ | KL divergence (objective) | `kl_loss` |
| $\\mathrm{KL}_{\\text{batch}}$ | Mini-batch KL estimator | `kl_batch` |
| $f(\\cdot)$ | Monotone sub-function | `f_mono()` |
| $g(\\cdot)$ | Non-monotone sub-function | `g_nonmono()` |
| $g(x_{1:k-1})$ | Conditional shift term | `g_shift()` |
| $f(x_k)$ | Univariate monotone basis | `f_u()` |
| $f(x_{1:k-1},x_k)$ | Cross-term basis | `f_cross()` |
| $c$ | Basis-function coefficient | `c_coef` |
| $c_{jk}$ | Individual fit coefficient | `c[j,k]` |
| $r(\\cdot)$ | Rectifier $r:\\mathbb{R}\\rightarrow\\mathbb{R}_{+}$ | `rectifier()` |
| $\\mathrm{He}_j(x)$ | Hermite polynomial of order $j$ | `hermite_j()` |
| $\\mathcal{H}$ | Basis-function space | `basis_set` |
| $\\tau$ | $2\\pi = 6.283185...$ | `tau <- 2*pi` |
| $\\Phi$ | Standard-normal CDF | `pnorm` |
| $\\varphi$ | Standard-normal PDF | `dnorm` |
| $U$ | Rank variable $u_{ik}=\\text{rank}/(N+1)$ | `U` |
| $Z = \\Phi^{-1}(U)$ | Probit transform | `Z_probit` |
| $\\mathrm{condSample}(\\cdot)$ | Draw from $\\pi(x_{k+1:K}\\mid x_{1:k}^*)$ | `cond_sample()` |

# Packages

## Transformation Forests
* trtf
* tram::traforest

## Reference

The triangular transport methodology is described in the preprint  
**arXiv:2503.21673v1 [stat.CO], 27 Mar 2025**.  
<https://arxiv.org/abs/2503.21673>

Our \(S_k(x_1,\ldots,x_k)\) functions are simply transformation forests for
the regression \(x_k \sim x_1,\ldots,x_{k-1}\). Monotonicity in \(x_k\) is
therefore automatically ensured.

The repository no longer stores the large PDF artifacts. Download the preprint
from the arXiv link above if needed.

## License

This project is licensed under the [MIT License](LICENSE).
