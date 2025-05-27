# Research Roadmap

This roadmap summarises the next steps for extending the multivariate conditional density experiments.  The notation follows `Theory.md`.

## 1. Additional Likelihood-Based Models
- Identify two further transformation models from the cited paper and write their joint log-densities
  $$\ell^{(j)}(x;\theta_j)=\sum_{n=1}^N \log p^{(j)}(x^{(n)};\theta_j)\,,\quad j\in\{A,B\}. $$
- Estimate parameters $\hat\theta_j$ by maximum likelihood.
- Evaluate the resulting log-likelihoods with the metric in §6 and append the new $\Delta^{(j)}$ to the summary table.

## 2. Triangular Order Permutation
- Select several non-trivial permutations $\sigma\in S_K$.
- Train the Transformation Forest and kernel estimator on each permuted order.
- Compute the degradation
  $$\Delta_\sigma=\mathbb E_{\hat\pi}[\ell^{(\mathrm{true})}(X)-\ell_\sigma(X)]\,.$$
- Compare these values with the baseline to interpret sparsity loss.

## 3. Correct Bivariate and Trivariate Baselines
- For $k\in\{2,3\}$ compute the exact conditional densities
  $$p_k(x_k\mid x_{1:k-1};\theta_k)$$
  implied by the chosen family.
- Form the joint log-likelihood
  $$\ell_k(\theta_k)=\sum_{n=1}^N\log p_k(x_k^{(n)}\mid x_{1:k-1}^{(n)};\theta_k).$$
- Optimise $\theta_k$ with a standard optimiser and verify on synthetic data that $\Delta^{(\mathrm{param})}\approx0$ when the model is correct.

## 4. Fixed Marginal Transformation
- Keep the default Box--Cox transform within the Transformation Forest and document this choice.

## 5. Remove "TF + Copula"
- The TF already outputs approximately independent components.
- A subsequent vine copula on those components is redundant and can reduce performance.
- If a copula experiment is desired, fit univariate marginals first and apply the vine directly.

## 6. Evaluation Protocol
- For each model $\mathcal M$ compute the test log-density matrix
  $$L^{(\mathcal M)}_{n,k}=\log p^{(\mathcal M)}(x^{(n)}_k\mid x^{(n)}_{1:k-1};\hat\Theta_{\mathcal M}).$$
- The average per-dimension log-likelihood is
  $$\bar\ell^{(\mathcal M)}_k=\frac{1}{N_{\mathrm{test}}}\sum_{n=1}^{N_{\mathrm{test}}}L^{(\mathcal M)}_{n,k}.$$
- The gap to truth is
  $$\Delta^{(\mathcal M)}_k=\bar\ell^{(\mathrm{true})}_k-\bar\ell^{(\mathcal M)}_k\,,\qquad\Delta^{(\mathcal M)}=\sum_{k=1}^{K}\Delta^{(\mathcal M)}_k.$$
- These $\Delta$ values remain the sole ranking criterion.

## Execution Order
1. Repair the parametric baseline.
2. Run permutation experiments on the corrected baseline and non-parametric models.
3. Add the new transformation-based models and recompute all $\Delta$.
4. Eliminate the invalid copula variant; freeze TF hyper-parameters.
5. Maintain the evaluation protocol and update the summary table.

