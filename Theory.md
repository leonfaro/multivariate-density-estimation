# Triangular Transport Methods (TTM)

We summarise the method in three concise blocks—mapping, sampling, and optimising—while keeping the notation from the README.

## Mapping

**Goal:** learn a lower-triangular map `S` so that `z = S(x)` with `z ~ eta` and `x ~ pi`.

| symbol | meaning |
|--------|------------------------------------------------|
| `pi`   | target distribution |
| `eta`  | reference distribution, e.g. `N(0,I)` |
| `S`    | map `R^K -> R^K` |
| `S_k`  | component `R^k -> R` |
| `S^-1` | inverse map |
| `K`    | dimension |

We enforce monotonicity in each `x_k`. The change of variables reads
`pi(x) = eta(S(x)) |det \nabla_x S(x)|`. Inversion is sequential:
`x_k = S_k^-1(z_k; x_{1:k-1})`. Sparsity comes from conditional independence, reducing the arguments of `S_k`.

## Sampling

1. Draw `Z ~ eta`.
2. Compute `X = S^{-1}(Z)` by solving each `S_k^-1` in turn.
3. Optional: fix the first `k` coordinates to sample from a conditional of `pi`.

## Optimising

Given data from `pi` or its unnormalised density, adjust the coefficients of `S` to make the pushforward close to `eta`—for instance by maximising the likelihood of `Z = S(X)` or minimising `\mathrm{KL}(\pi \| S_\# \eta)`.

