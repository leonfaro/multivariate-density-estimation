# Theory Workflow

The diagram below summarises the algorithmic structure implied by [Theory.md](Theory.md). Boxes represent mathematical objects or operations.

```mermaid
flowchart TD
    XData["Training and test data $x^{(i)}$"] -->|"split"| DataSplit["$X_{\text{train}}, X_{\text{test}}$"]
    DataSplit --> Param["Parameter functions $\theta_k(x_{1:k-1},x_{\text{cov}})$"]
    Param --> Density["Conditional densities $f_{d_k}(x_k\mid x_{1:k-1}, x_{\text{cov}})$"]
    Density --> CDF["Conditional CDFs $F_{d_k}(x_k\mid x_{1:k-1}, x_{\text{cov}})$"]
    CDF --> MapS["Monotone map $S$ with components $S_k=\Phi^{-1}\circ F_{d_k}$"]
    MapS --> Inverse["Inverse transform $S^{-1}$ for sampling"]
    Inverse --> Samples["Samples $x \sim \pi(\cdot\mid x_{\text{cov}})$"]
    MapS -->|"$z=S(x)$"| Likelihood["Log-likelihood $\ell(x)$"]
    Likelihood --> Estim["Maximise $\sum_i \ell(x^{(i)})$"]
    Estim --> Param
    Param --> Basis["Basis expansion $\theta_k\approx\sum_m c_{mk}\psi_{mk}$"]
    Param --> Sparse["Sparsity sets $J_k\subset\{1,\dots,k-1\}$"]
    CDF -->|"replace $F_{d_k}$"| Forest["Transformation forest variant"]
    Forest --> MapS
```

The workflow traces how data inform the parameter functions, which in turn define the conditional densities and the triangular transport map. Sampling uses the inverse map, while likelihood evaluation feeds back into parameter estimation. Basis expansions and sparsity restrict the parameter functions. A non-parametric alternative replaces the parametric CDFs by transformation forests.
