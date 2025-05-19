# Multivariate Conditional Density Estimation: Transformation Forests, Copulas, and Normalizing Flows

## Triangular Factorisation

For target variables \(X_1,\dots,X_K\) we factorise the joint density
\[
\pi(x_1,\ldots,x_K)
  = \prod_{k=1}^{K} f_{d_k}(x_k \mid x_{1:k-1}),
\]
where each univariate density \(f_{d_k}\) has parameter function \(\theta_k(x_{1:k-1})\).


## Transport Map
Define conditional CDFs
\[

F_{d_k}(x_k \mid x_{1:k-1})
  = \int_{-\infty}^{x_k} f_{d_k}(t \mid x_{1:k-1})\, \mathrm dt.
\]
The monotone triangular map \(S:\mathbb R^K\to\mathbb R^K\) is
\[
z_k = S_k(x_1,\dots,x_k) = \Phi^{-1}\!\bigl(F_{d_k}(x_k \mid x_{1:k-1})\bigr),

\qquad k=1,\dots,K.
\]

## Inverse Transform (Sampling)
For \(z\sim\mathcal N(0,I_K)\) set \(u_k=\Phi(z_k)\) and compute sequentially
\[

x_k = F_{d_k}^{-1}(u_k \mid x_{1:k-1}),
\]
yielding \(x\sim\pi\).


## Likelihood
Change of variables gives
\[

\pi(x)=
  \eta\!\bigl(S(x)\bigr)\,
  \prod_{k=1}^{K}
    \frac{f_{d_k}(x_k \mid x_{1:k-1})}{\varphi(z_k)},

\]
and the log‑likelihood
\[
\ell(x)=
  -\tfrac12\|S(x)\|^2
  + \sum_{k=1}^{K}\!\bigl[\log f_{d_k}(x_k \mid x_{1:k-1})-\log\varphi(z_k)\bigr].

\]
Maximising \(\sum_i \ell(x^{(i)})\) estimates \(\theta_k\).

## Sparsity
Impose conditional independence by restricting
\(\theta_k\) to depend only on \(J_k\subset\{1,\dots,k-1\}\):

\(X_k\perp X_j\mid\{X_i:i\in J_k\}\) for \(j\notin J_k\).

Sparsity lowers complexity.

## Flexible Parameterisation
Approximate the parameter functions via basis expansion
\[

\theta_k(x_{1:k-1})
  \approx \sum_m c_{mk}\,\psi_{mk}(x_{1:k-1}),

\]
with a monotone link to keep valid conditional densities.

## Special Cases
* **Gaussian copula** – Gaussian \(f_{d_k}\) with linear \(\theta_k\).
* **Transformation forest** – tree‑based \(\theta_k\).
* **D-vine Copula with TF marginals** – each transformation forest estimates the conditional CDF \(\hat F_k(y_k\mid x)\). D-vine copula models only the dependence among uniform scores \(U_{ik}=\hat F_k(y_{ik}\mid x_i)\), or their probit transforms \(Z_{ik}=\Phi^{-1}(U_{ik})\). conditional CDF \(\rightarrow\) ranks \(U\) \(\rightarrow\) copula.

* **Autoregressive normalising flow** – neural‑network \(\theta_k\).

## Model Inputs
Algorithms use only observation vectors \(x^{(i)}\in\mathbb R^3\) split into training \((i=1,\dots,N_{\text{train}})\) and test samples \((i=1,\dots,N_{\text{test}})\). Processing order is \(k=1,2,3\).

## Implementation Details
### Triangular Transport Map

* **Sampling:** draw \(z\sim\mathcal N(0,I_K)\), set \(x_k=q_{d_k}(\Phi(z_k)\mid x_{<k})\) sequentially.
* **Density Evaluation:** given \(x\), compute \(z=S(x)\) and accumulate \(\log f_{d_k}(x_k\mid x_{<k})-\log\varphi(z_k)\).

### Transformation Forest Variant
For each \(k\) fit a transformation forest of \(X_k\) on \(X_{<k}\) to obtain a monotone empirical CDF \(\hat F_k\). Substitute \(\hat F_k\) and its inverse for \(F_{d_k}\) in the map:
\[
S_k(x_1,\dots,x_k)=\Phi^{-1}\!\bigl(\hat F_k(x_k\mid x_{<k})\bigr).

### D-vine Copula with TF marginals** – 
For each transformation forest estimates the conditional CDF \(\hat F_k(y_k\mid x)\). D-vine copula models only the dependence among uniform scores \(U_{ik}=\hat F_k(y_{ik}\mid x_i)\), or their probit transforms \(Z_{ik}=\Phi^{-1}(U_{ik})\). conditional CDF \(\rightarrow\) ranks \(U\) \(\rightarrow\) copula.
\]



