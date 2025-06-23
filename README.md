# Parametric Baseline Extension
Read `Theory.md` and 'roadmap.md' before starting any task. Adjust ALGORITHM_SPEC.md after every task if necessary.
This repository implements simple experiments for multivariate conditional density estimation.
All notation follows `Theory.md` and the instructions in `AGENTS.md`.

## Background

For a random vector $x = (x_1,\ldots,x_K)$ with density $\pi(x)$ we study a
lower--triangular transport $S$ such that $z = S(x)$ has independent standard
normal components.  Densities relate via

$$\pi(x) = \eta(S(x)) |\det \nabla_x S(x)|,$$
where $\eta$ denotes the standard Gaussian.

Optimization uses `optim` with the BFGS method on the negative log-likelihood.
All densities and likelihoods are computed in log-space from `00_setup.R` on;
only the final presentation of summary metrics converts them back to standard
scale.

## Configuration Syntax

The configuration `config` is a list of length $K`.  Each entry has a field
`distr` specifying the family.  

## End_of_readme

