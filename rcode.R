# generates time-series "raw" data that fits to the given correlation matrix(=functional connectome). 

fcon <- as.matrix(read_csv("~/fcondata/fcon.csv", col_names = FALSE))
make.positive.definite(fcon) -> fcon_pd 
fcon_sym <- (fcon_pd+t(fcon_pd))/2    # make symmetric as covariance

nobs = 2000                           # T scale while simulating, include tau 
nvars = dim(chol(fcon_sym))[1]        # "nodes", 30

rv_fcon = t(chol(fcon_sym)) %*% matrix(rnorm(nvars*nobs), nrow=nvars, ncol=nobs)    # RVs that follows cnxn corr(cov) matrix
rv_fcon_data = as.data.frame(t(rv_fcon))                                            # use t() if needed 

VAR(rv_fcon_data, p=1) -> fcon_raw                     # vector AR for time-series linear regression 
Bcoef(fcon_raw) -> fcon_A                              # regression coefficient matrix 
cov(residuals(VAR(rv_fcon_data, p=1))) -> fcon_E       # covariance matrix of residuals

View(fcon_sym)
View(summary(fcon_raw)$cov)                            # should be almost identical with 'correlation(connectivity) matrix' we used

write.csv(fcon_A, file = "fcon_A.csv")
write.csv(fcon_E, file = "fcon_E.csv")
