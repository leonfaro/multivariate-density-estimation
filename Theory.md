## Theory

### Triangular factorisation

Consider target variables $X_1,\ldots,X_K$. For each index $k$ choose a univariate density family $f_{d_k}$ with parameter function $\theta_k(x_{1:k-1}, x_{\text{cov}})$. The conditional joint density is
\[
  \pi(x_1,\ldots,x_K \mid x_{\text{cov}})
    = \prod_{k=1}^K f_{d_k}(x_k \mid x_{1:k-1}, x_{\text{cov}}).
\]

### Transport map

Define
\[
  F_{d_k}(x_k \mid x_{1:k-1}, x_{\text{cov}})
    = \int_{-\infty}^{x_k} f_{d_k}(t \mid x_{1:k-1}, x_{\text{cov}})\,\mathrm dt.
\]
The triangular map $S\colon\mathbb{R}^K\rightarrow\mathbb{R}^K$ sends $x=(x_1,\ldots,x_K)$ to $z=(z_1,\ldots,z_K)$ via
\[
  z_k = S_k(x_1,\ldots,x_k)
      = \Phi^{-1}\bigl(F_{d_k}(x_k \mid x_{1:k-1}, x_{\text{cov}})\bigr).
\]
Each component is monotone in $x_k$, hence $S$ is invertible.

### Inverse transform

Let $z\sim\mathcal N(0,I_K)$. Set $u_k=\Phi(z_k)$ and compute
\[
  x_k = F_{d_k}^{-1}(u_k \mid x_{1:k-1}, x_{\text{cov}}),\quad k=1,\ldots,K.
\]
This yields $x\sim\pi(\cdot\mid x_{\text{cov}})$. Conversely, any $x$ maps to $z=S(x)$.

### Likelihood

The change-of-variables formula gives
\[
  \pi(x\mid x_{\text{cov}})
    = \eta(S(x)) \prod_{k=1}^K \frac{f_{d_k}(x_k \mid x_{1:k-1}, x_{\text{cov}})}{\varphi(z_k)},
\]
where $\eta$ is the $K$-variate standard normal density and $\varphi$ is its one-dimensional version. The log-likelihood reads
\[
  \ell(x) = -\tfrac{1}{2}\|S(x)\|^2
             + \sum_{k=1}^K \bigl[\log f_{d_k}(x_k \mid x_{1:k-1}, x_{\text{cov}}) - \log \varphi(z_k)\bigr].
\]
Maximising $\sum_i \ell(x^{(i)})$ estimates the functions $\theta_k$.

### Sparsity and functional form

Conditional independence enters by restricting $\theta_k$ to depend only on $J_k\subseteq\{1,\ldots,k-1\}$. Excluding $X_j$ with $j\notin J_k$ yields $X_k\perp X_j\mid\{X_i:i\in J_k\},x_{\text{cov}}$. A sparse choice of $J_k$ reduces complexity.

We may expand
\[
  \theta_k(x_{1:k-1}, x_{\text{cov}})
    \approx \sum_{m} c_{mk}\,\psi_{mk}(x_{1:k-1}, x_{\text{cov}}),
\]
where $\{\psi_{mk}\}$ is a set of basis functions. The monotonic CDF link ensures validity of the conditional densities.

### Learning algorithm

Gradient-based optimisation or boosting maximises the log-likelihood. One strategy starts with the identity map and adds basis terms that increase the objective most. The triangular structure allows sequential or block-wise training.

### Special cases

- Gaussian $f_{d_k}$ with linear $\theta_k$ yield a Gaussian copula.
- Tree-based parameter functions lead to transformation forests.
- Neural network parameter functions give autoregressive normalising flows.

### Model inputs

The algorithms use only the observation vectors $x^{(i)}\in\mathbb{R}^3$. Let $\{x^{(i)}\}_{i=1}^{N_{\text{train}}}$ denote the training set and $\{x^{(i)}\}_{i=1}^{N_{\text{test}}}$ the test set. Each method processes the coordinates in the order $k=1,2,3$ without additional information.
