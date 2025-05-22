# Parametric Baseline Extension

This repository implements simple experiments for multivariate conditional density estimation.
All notation follows `Theory.md` and the instructions in `AGENTS.md`.

## Background

For a random vector $x = (x_1,\ldots,x_K)$ with density $\pi(x)$ we study a
lower--triangular transport $S$ such that $z = S(x)$ has independent standard
normal components.  Densities relate via

$$\pi(x) = \eta(S(x)) |\det \nabla_x S(x)|,$$
where $\eta$ denotes the standard Gaussian.

The parametric baseline fits one--dimensional conditionals
$$\pi(x_k\mid x_{k-1})$$
using maximum likelihood.  Each conditional belongs to one of eight canonical
families listed below.  Parameters are linear in the previous coordinate and
transformed by either the identity or the softplus link
\[\operatorname{softplus}(z) = \log(1+e^z).\]

## Available Families

1. **Normal** $(\mu,\sigma)$
2. **Exponential** $(\lambda)$
3. **Gamma** $(\text{shape}, \text{rate})$
4. **Weibull** $(\text{shape}, \text{scale})$
5. **Lognormal** $(\mu_{\log}, \sigma_{\log})$
6. **Poisson** $(\lambda)$
7. **Beta** $(\alpha,\beta)$
8. **Logistic** $(\text{loc}, \text{scale})$

For each parameter $\theta$ we fit
$$\theta(x_{k-1}) = \ell\big(\beta_0 + \beta_1 x_{k-1}\big),$$
where $\ell$ is `identity` for unconstrained parameters and `softplus` for
positive ones.
Optimization uses `optim` with the BFGS method on the negative log-likelihood.

## Configuration Syntax

The configuration `config` is a list of length $K`.  Each entry has a field
`distr` specifying the family.  The order of `config` determines the factorization
of $\pi(x)$.

## Tests

Unit tests are written with **testthat** and executed via `run_checks.sh`:

```
Rscript basic_tests.R
Rscript -e "testthat::test_dir('tests/testthat')"
bash lint.sh
```

The new tests verify positivity of `softplus`, successful recovery of synthetic
Gamma parameters and monotonic improvement of the likelihood.

## Example

The script `example_bivariate_fit.R` performs a small bivariate fit and prints a
summary table.
