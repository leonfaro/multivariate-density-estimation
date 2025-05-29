\begin{algorithm}[H]
\caption{DataSetup: Conditional 4-D Sample Generation and Pre-processing}
\begin{algorithmic}[1]
\textbf{Input:} sample size $N$; split ratio $\rho\!\in\!(0,1)$; configuration $\mathcal{C}$ (four pairs $\bigl(\text{distr}_k,\text{parm}_k\bigr)$); random seed $s$\\
\textbf{Output:} training matrix $\mathbf{X}^{\text{train}}\!\in\!\mathbb{R}^{\lfloor\rho N\rfloor\times4}$, test matrix $\mathbf{X}^{\text{test}}\!\in\!\mathbb{R}^{(N-\lfloor\rho N\rfloor)\times4}$, standardisation stats $\{(\mu_k,\sigma_k)\}_{k=1}^{4}$

\State Set pseudorandom generator seed $\gets s$
\For{$i \gets 1$ \textbf{to} $N$}                                    \Comment{simulate one observation}
    \State Sample $x_{i1} \sim \mathcal{N}(0,1)$
    \State $\theta_2 \gets \textsc{EvalParams}\bigl(\text{parm}_2;\,x_{i1}\bigr)$
    \State Sample $x_{i2} \sim \operatorname{Exp}\!\bigl(\text{rate}=\theta_2\bigr)$
    \State $\theta_3 \gets \textsc{EvalParams}\bigl(\text{parm}_3;\,x_{i1},x_{i2}\bigr)$
    \State Sample $x_{i3} \sim \operatorname{Beta}\!\bigl(\theta_{3,1},\theta_{3,2}\bigr)$
    \State $\theta_4 \gets \textsc{EvalParams}\bigl(\text{parm}_4;\,x_{i1},x_{i2},x_{i3}\bigr)$
    \State Sample $x_{i4} \sim \Gamma\!\bigl(\text{shape}=\theta_{4,1},\text{scale}=\theta_{4,2}\bigr)$
\EndFor
\State Assemble $\mathbf{X}\in\mathbb{R}^{N\times4}$ from rows $(x_{i1},x_{i2},x_{i3},x_{i4})$
\State $n_{\text{train}}\gets\lfloor\rho N\rfloor$
\State $\mathbf{X}^{\text{train}}\gets\mathbf{X}_{1:n_{\text{train}},\,\cdot}$;\;
       $\mathbf{X}^{\text{test}}\gets\mathbf{X}_{n_{\text{train}}+1:N,\,\cdot}$
\For{$k \gets 1$ \textbf{to} 4}                                         \Comment{z-score standardisation}
    \State $\mu_k\gets\operatorname{mean}\!\bigl(\mathbf{X}^{\text{train}}_{\cdot,k}\bigr)$
    \State $\sigma_k\gets\operatorname{sd}\!\bigl(\mathbf{X}^{\text{train}}_{\cdot,k}\bigr)$
    \ForAll{rows $j$ in $\mathbf{X}^{\text{train}}$ and $\mathbf{X}^{\text{test}}$}
        \State $x_{jk}\gets\bigl(x_{jk}-\mu_k\bigr)/\sigma_k$
    \EndFor
\EndFor
\State \Return $\mathbf{X}^{\text{train}},\mathbf{X}^{\text{test}},\{(\mu_k,\sigma_k)\}_{k=1}^{4}$
// Time: $O(NK)$, Space: $O(NK)$
\\[4pt]
\Procedure{EvalParams}{parm, $\mathbf{z}$}                    \Comment{helper: evaluate conditional parameters}
    \If{parm $\equiv$ NULL} \Return NULL \EndIf
    \State \Return $parm(\mathbf{z})$
\EndProcedure
\end{algorithmic}
\end{algorithm}

