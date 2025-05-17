source('Code/05_joint_evaluation.R', chdir = TRUE)
result <- all.equal(sum(ell_forest), sum(ll_test), tol = 1e-1)
if (!isTRUE(result)) stop(result)
cat('Test passed: log-likelihood sums agree within tolerance\n')

