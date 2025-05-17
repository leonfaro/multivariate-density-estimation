# Notation

| Symbol | Meaning (paper context) | Recommended R naming |
|--------|-------------------------|----------------------|
| **S** | Target->Reference map $S:\mathbb{R}^K \rightarrow \mathbb{R}^K$ | `S` / `S_map()` |
| $S_k$ | $k$-th map component, monotone in $x_k$ | `S_k()` / `S_list[[k]]` |
| $S^{-1}$ | Inverse map (Reference->Target) | `S_inv()` |
| $S^{-1}_k$ | Sequential inversion of the $k$-th component | `S_inv_k()` |
| $S^{\sharp}\pi$ | Push-forward distribution (image of $\pi$ under $S$) | `push_fwd` |
| $S_{\sharp}\eta$ | Pull-back distribution (pre-image of $\eta$ under $S^{-1}$) | `pull_back` |
| **R** | Reference->Target map (Sec. 3.2.2) | `R` / `R_map()` |
| $\nabla S$ | Jacobian matrix of the map | `J_S` |
| $\det\nabla S$ | Determinant (triangular -> product of diagonals) | `det_J` |
| $\partial_{x_k}S_k$ | Diagonal Jacobian entry | `dSdx_k` |
| **$\pi$** | Target distribution (orange variables) | `pi_dist` |
| $\pi(x)$ | Target density function | `pi_density` |
| **$\eta$** | Reference distribution (green vars), usually $\mathcal N(0,I)$ | `eta_dist` |
| $\eta(z)$ | Reference density function | `eta_density` |
| **x** | Random vector $x \sim \pi$ | `X` (matrix N x K) |
| **z** | Random vector $z \sim \eta$ | `Z` |
| $X_i$ | $i$-th realisation of $x$ | `X[i, ]` |
| $Z_i$ | $i$-th realisation of $z$ | `Z[i, ]` |
| orange variable | Any quantity in the $\pi$ domain | suffix `_pi` |
| green variable | Any quantity in the $\eta$ domain | suffix `_eta` |
| $K$ | Number of target dimensions | `K` |
| $N$ | Ensemble / sample size | `N` |
| $p(\cdot)$ | Generic pdf | `p_fun` |
| $a, b$ | Placeholder random variables | `a`, `b` |
| $y^*$ | Conditioning variable | `y_star` |
| $x^*$ | Conditional target variable | `x_star` |
| $L$ | Likelihood $\eta(S(x))\,|\det\nabla S|$ | `likelihood` |
| $\ell$ | Log-likelihood | `loglik` |
| $\|S(x)\|^2$ | Squared norm of map outputs | `S_norm2` |
| $\mathrm{KL}(\pi\,\|\,S_{\sharp}\eta)$ | KL divergence (objective) | `kl_loss` |
| $\mathrm{KL}_{\text{batch}}$ | Mini-batch KL estimator | `kl_batch` |
| $f(\cdot)$ | Monotone sub-function | `f_mono()` |
| $g(\cdot)$ | Non-monotone sub-function | `g_nonmono()` |
| $g(x_{1:k-1})$ | Conditional shift term | `g_shift()` |
| $f(x_k)$ | Univariate monotone basis | `f_u()` |
| $f(x_{1:k-1},x_k)$ | Cross-term basis | `f_cross()` |
| $c$ | Basis-function coefficient | `c_coef` |
| $c_{jk}$ | Individual fit coefficient | `c[j,k]` |
| $r(\cdot)$ | Rectifier $r:\mathbb{R}\rightarrow\mathbb{R}_{+}$ | `rectifier()` |
| $\mathrm{He}_j(x)$ | Hermite polynomial of order $j$ | `hermite_j()` |
| $\mathcal{H}$ | Basis-function space | `basis_set` |
| $\tau$ | $2\pi = 6.283185...$ | `tau <- 2*pi` |
| $\Phi$ | Standard-normal CDF | `pnorm` |
| $\varphi$ | Standard-normal PDF | `dnorm` |
| $U$ | Rank variable $u_{ik}=\text{rank}/(N+1)$ | `U` |
| $Z = \Phi^{-1}(U)$ | Probit transform | `Z_probit` |
| $\mathrm{condSample}(\cdot)$ | Draw from $\pi(x_{k+1:K}\mid x_{1:k}^*)$ | `cond_sample()` |
