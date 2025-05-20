# PhD Thesis in Multivariate Conditional Density Estimation with Likelihood Inference and Regression Analysis: Comparing Transformation Forest Models, Copulas and Normalizing Flows

## Motivation
* **Goal:** Estimate a flexible conditional density \$p(y \mid x)\$ for multivariate continuous outcomes \$Y \in \mathbb{R}^d\$ given features \$X \in \mathbb{R}^p\$, with full likelihood inference (not just point predictions).
* **Context:** Distributional regression - interest in the entire outcome distribution (beyond mean or variance). Relevant for heteroscedastic or multi-modal responses.
* **Challenges:** High-dimensional \$Y\$ with complex interdependencies; \$p(y|x)\$ may be multi-modal, non-Gaussian, and vary strongly with \$x\$. Needs models that are both flexible (to fit complex distributions) and tractable (for likelihood evaluation and sampling).
* **Approaches:** Consider both statistical and ML methods:

  * *Structured invertible transformation* (triangular transport map) for exact likelihood modeling.
  * *Copulas* for modular dependency modeling (with separate marginal fits).
  * *Transformation forests* for nonparametric conditional distribution estimation.
  * *Normalizing flows* (neural invertible networks) for maximum flexibility.

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
  almost everywhere, so each \(S_k\) is monotone in its last argument and \(S\) is invertible .
\]



