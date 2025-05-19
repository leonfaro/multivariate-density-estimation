## Theory
* **Triangular factorization:** Model a joint or conditional density via sequential univariate factors. For target variables \(X_1,\dots,X_K\) (components of \(Y\)) conditioned on external covariates \(x_{\text{cov}}\):
  $\pi(x_1,\dots,x_K \mid x_{\text{cov}}) = \prod_{k=1}^K f_{d_k}\!\big(x_k \mid x_{1:k-1},\,x_{\text{cov}}\big),$
  where each \(f_{d_k}\) is a chosen univariate density family with parameters \(\theta_k(x_{<k},x_{\text{cov}})\).
* **Transport map $S$:** Define each conditional CDF $F_{d_k}(x_k \mid x_{<k},x_{\text{cov}})$. The *transport map* $S: \mathbb{R}^K \to \mathbb{R}^K$ sends $x=(x_1,\dots,x_K)$ to latent $z=(z_1,\dots,z_K)$ by
  $z_k = S_k(x_k \mid x_{<k},x_{\text{cov}}) = \Phi^{-1}\!\Big( F_{d_k}\!\big(x_k \mid x_{<k},x_{\text{cov}}\big) \Big),$
  where $\Phi^{-1}$ is the standard normal quantile. Each $S_k$ is monotonic in $x_k$, making $S$ invertible (triangular Jacobian).
* **Inverse transform:** Invert $S$ sequentially for sampling or inference. Given $z \sim \mathcal{N}(0,I_K)$, set $u_k = \Phi(z_k)$ and then
  $x_k = F_{d_k}^{-1}(u_k \mid x_{<k}, x_{\text{cov}}), \qquad k=1,\dots,K,$
  to obtain $x \sim \pi(\cdot \mid x_{\text{cov}})$. Conversely, any $x$ maps to $z=S(x)$.
* **Likelihood formula:** By change of variables,
  $\pi(x \mid x_{\text{cov}}) = \eta\!\big(S(x)\big)\;\prod_{k=1}^K \frac{f_{d_k}(x_k \mid x_{<k},x_{\text{cov}})}{\varphi(z_k)},$
  with $\eta(z)=\prod_{k=1}^K \varphi(z_k)$ the standard normal density and $\varphi$ the $\mathcal{N}(0,1)$ PDF. The log-likelihood for $x$ is
  $\ell(x) = -\tfrac{1}{2}\|S(x)\|^2 + \sum_{k=1}^K \Big[\log f_{d_k}\!\big(x_k \mid x_{<k},x_{\text{cov}}\big) - \log \varphi(z_k)\Big].$
  Maximize $\sum_i \ell(x^{(i)})$ to estimate the functions $\theta_k$ (maximum likelihood estimation).
* **Sparsity:** Impose conditional independence by restricting each $\theta_k$ to a subset of $x_{<k}$. If $j \notin J_k \subseteq \{1,\dots,k-1\}$, then $X_j$ is excluded from $f_{d_k}(x_k|\cdot)$ (implying $X_k \perp X_j \mid \{X_i: i\in J_k\},x_{\text{cov}}$ in the model). Sparse dependencies reduce complexity.
* **Flexible parameterization:** Specify $\theta_k(x_{<k},x_{\text{cov}})$ in a rich function class. For example, use a basis expansion
  $\theta_k(x_{<k},x_{\text{cov}}) \approx \sum_{m} c_{mk}\,\psi_{mk}(x_{<k},x_{\text{cov}}),$
  with $\{\psi_{mk}\}$ basis functions (polynomials, splines, tree indicators, etc). This captures nonlinear effects. The monotonic CDF link ensures the conditional density is valid.
* **Learning algorithm:** Train by maximizing log-likelihood. Use gradient-based optimization or boosting. A *map adaptation* strategy can start with an identity map and add basis terms to $\theta_k$ iteratively (greedy selection by largest log-lik improvement). Triangular structure permits sequential or block-wise training of components.
* **Unifying cases:** This framework recovers many known models as special cases:

  * All $d_k$ Gaussian with linear $\theta_k$ $\;\to\;$ Gaussian copula (multivariate normal with regression).
  * Tree-based $\theta_k$ $\;\to\;$ transformation forest models for each conditional.
  * Neural network $\theta_k$ $\;\to\;$ autoregressive normalizing flow.
* **Model inputs:** The three estimators—transformation forest, D-vine copula, and autoregressive normalizing flow—are trained and tested solely on the observation vectors $x^{(i)}\in\mathbb{R}^3$. Let $(x^{(i)})_{i=1}^{N_{\text{train}}}$ denote the training sample and $(x^{(i)})_{i=1}^{N_{\text{test}}}$ the test sample. Each algorithm processes the coordinates sequentially in the order $k=1,2,3$ without any side information about $\pi$ or the conditional factors. Parameter estimation relies only on these data.
