library(AER)
library(hdm)
library(glmnet)
library(ranger)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

set.seed(2026)
sig <- function(x) 1 / (1 + exp(-x))

sim_data <- function(n, gamma, dgp) {
  X <- matrix(rnorm(n * 5), n, 5); colnames(X) <- paste0("X", 1:5)
  Z <- rbinom(n, 1, sig(X[,1] + 0.5 * X[,2]))
  D <- ifelse(Z == 1, rbinom(n, 1, sig(-1 + 0.5 * X[,1] + gamma)),
                      rbinom(n, 1, sig(-1 + 0.5 * X[,1])))
  g <- switch(dgp,
    linear    = X[,1] + 0.6*X[,2] - 0.4*X[,3] + 0.3*X[,4] - 0.2*X[,5],
    nonlinear = 0.5*X[,1]^2 + sin(pi*X[,2]) + X[,1]*X[,2] - 0.4*X[,3]^2,
    sparse    = X[,1])
  list(Y = D + g + rnorm(n), D = D, Z = Z, X = X)
}

get_nuisance <- function(Y, D, Z, X, K = 3) {
  n   <- length(Y)
  fld <- sample(rep(1:K, length.out = n))
  g <- e <- d <- mY1 <- mY0 <- mD1 <- mD0 <- numeric(n)
  lam <- exp(seq(log(0.5), log(0.001), length.out = 20))
  for (k in 1:K) {
    tr  <- fld != k;  te <- fld == k
    nms <- c(paste0("X", 1:5), "Z")
    XZ  <- setNames(as.data.frame(cbind(X[tr,], Z = Z[tr])), nms)
    X1  <- setNames(as.data.frame(cbind(X[te,], Z = 1)),     nms)
    X0  <- setNames(as.data.frame(cbind(X[te,], Z = 0)),     nms)
    g[te] <- predict(cv.glmnet(X[tr,], Y[tr], nfolds = 3, lambda = lam),
                     X[te,], s = "lambda.min")
    d[te] <- predict(cv.glmnet(X[tr,], D[tr], nfolds = 3, lambda = lam),
                     X[te,], s = "lambda.min")
    e[te] <- pmax(pmin(
               predict(cv.glmnet(X[tr,], Z[tr], nfolds = 3, lambda = lam,
                                 family = "binomial"),
                       X[te,], type = "response", s = "lambda.min"),
               0.95), 0.05)
    fy <- suppressWarnings(
      ranger(Y ~ ., data = cbind(Y = Y[tr], XZ),
             num.trees = 100, num.threads = 1))
    fd <- suppressWarnings(
      ranger(D ~ ., data = cbind(D = D[tr], XZ),
             num.trees = 100, num.threads = 1))
    mY1[te] <- predict(fy, X1)$predictions
    mY0[te] <- predict(fy, X0)$predictions
    mD1[te] <- predict(fd, X1)$predictions
    mD0[te] <- predict(fd, X0)$predictions
  }
  list(g = g, e = e, d = d, mY1 = mY1, mY0 = mY0, mD1 = mD1, mD0 = mD0)
}

e_2sls <- function(Y, D, Z, X) {
  df <- as.data.frame(cbind(Y = Y, D = D, Z = Z, X))
  v  <- paste(paste0("X", 1:5), collapse = "+")
  ft <- ivreg(as.formula(paste("Y ~ D +", v, "| Z +", v)), data = df)
  s  <- summary(ft)
  c(tau = unname(coef(ft)["D"]), se = coef(s)["D", "Std. Error"])
}

e_pl <- function(Y, D, Z, X) {
  sel_y <- which(coef(rlasso(X, Y))[-1] != 0)
  sel_d <- which(coef(rlasso(X, D))[-1] != 0)
  sel   <- union(sel_y, sel_d)
  if (length(sel) == 0) sel <- 1L
  X_sel <- X[, sel, drop = FALSE]
  df    <- as.data.frame(cbind(Y = Y, D = D, Z = Z, X_sel))
  v     <- paste(colnames(X_sel), collapse = "+")
  ft    <- tryCatch(
    ivreg(as.formula(paste("Y ~ D +", v, "| Z +", v)), data = df),
    error = function(e) NULL)
  if (is.null(ft)) return(c(tau = NA, se = NA))
  s <- summary(ft)
  c(tau = unname(coef(ft)["D"]), se = coef(s)["D", "Std. Error"])
}

e_dml <- function(Y, D, Z, nu) {
  Yt    <- Y - nu$g;  Dt <- D - nu$d;  Zt <- Z - nu$e
  denom <- mean(Zt * Dt)
  if (abs(denom) < 1e-4) return(c(tau = NA, se = NA))
  tau <- mean(Zt * Yt) / denom
  psi <- Zt * (Yt - tau * Dt)
  c(tau = tau, se = sqrt(mean(psi^2) / denom^2 / length(Y)))
}

e_aipw <- function(Y, D, Z, nu) {
  ch  <- pmax(nu$mD1 - nu$mD0, 0.01)
  w   <- (Z - nu$e) / (nu$e * (1 - nu$e))
  pY  <- w * (Y - ifelse(Z == 1, nu$mY1, nu$mY0)) + nu$mY1 - nu$mY0
  pD  <- w * (D - ifelse(Z == 1, nu$mD1, nu$mD0)) + ch
  denom <- mean(pD)
  if (abs(denom) < 0.01) return(c(tau = NA, se = NA))
  tau <- mean(pY) / denom
  c(tau = tau, se = sd((pY - tau * pD) / denom) / sqrt(length(Y)))
}

e_tmle <- function(Y, D, Z, nu) {
  e <- nu$e
  H <- (Z - e) / (e * (1 - e))

  mYz  <- ifelse(Z == 1, nu$mY1, nu$mY0)
  epsY <- sum(H * (Y - mYz)) / sum(H^2)
  mY1s <- nu$mY1 + epsY / e
  mY0s <- nu$mY0 - epsY / (1 - e)

  mDz  <- ifelse(Z == 1, nu$mD1, nu$mD0)
  epsD <- sum(H * (D - mDz)) / sum(H^2)
  mD1s <- nu$mD1 + epsD / e
  mD0s <- nu$mD0 - epsD / (1 - e)
  ch_s <- pmax(mD1s - mD0s, 0.01)

  denom <- mean(ch_s)
  if (abs(denom) < 0.01) return(c(tau = NA, se = NA))

  tau  <- mean(mY1s - mY0s) / denom
  mYzs <- ifelse(Z == 1, mY1s, mY0s)
  mDzs <- ifelse(Z == 1, mD1s, mD0s)
  pY   <- H * (Y - mYzs) + mY1s - mY0s
  pD   <- H * (D - mDzs) + ch_s
  c(tau = tau, se = sd((pY - tau * pD) / denom) / sqrt(length(Y)))
}

one_rep <- function(n, gamma, dgp) {
  d  <- sim_data(n, gamma, dgp)
  nu <- with(d, get_nuisance(Y, D, Z, X))
  r1 <- tryCatch(with(d, e_2sls(Y, D, Z, X)),  error = function(e) c(NA, NA))
  r2 <- tryCatch(with(d, e_pl(Y, D, Z, X)),    error = function(e) c(NA, NA))
  r3 <- tryCatch(with(d, e_dml(Y, D, Z, nu)),  error = function(e) c(NA, NA))
  r4 <- tryCatch(with(d, e_aipw(Y, D, Z, nu)), error = function(e) c(NA, NA))
  r5 <- tryCatch(with(d, e_tmle(Y, D, Z, nu)), error = function(e) c(NA, NA))
  rbind(r1, r2, r3, r4, r5)
}

run_sim <- function(n, gamma, dgp, R = 100) {
  reps <- replicate(R, one_rep(n, gamma, dgp), simplify = FALSE)
  taus <- sapply(reps, function(r) r[, "tau"])
  ses  <- sapply(reps, function(r) r[, "se"])
  data.frame(
    method   = c("2SLS", "Post-Lasso", "DML-IV", "AIPW-IV", "IV-TMLE"),
    n = n, gamma = gamma, dgp = dgp,
    bias     = round(rowMeans(taus - 1,                na.rm = TRUE), 3),
    rmse     = round(sqrt(rowMeans((taus-1)^2,         na.rm = TRUE)), 3),
    coverage = round(rowMeans(abs(taus-1) <= 1.96*ses, na.rm = TRUE) * 100, 1),
    ci_width = round(rowMeans(2*1.96*ses,              na.rm = TRUE), 3),
    na_frac  = round(rowMeans(is.na(taus)), 3),
    row.names = NULL)
}

grid <- expand.grid(
  n     = c(500, 1000, 2000),
  gamma = c(0.5, 1.0, 2.0),
  dgp   = c("linear", "nonlinear", "sparse"),
  stringsAsFactors = FALSE)
grid <- grid[order(-grid$n), ]

n_cores <- max(1L, detectCores() - 1L)
cl      <- makeCluster(n_cores)
clusterExport(cl, c("sig", "sim_data", "get_nuisance",
                    "e_2sls", "e_pl", "e_dml", "e_aipw", "e_tmle",
                    "one_rep", "run_sim", "grid"))
invisible(clusterEvalQ(cl, {
  library(AER); library(hdm); library(glmnet); library(ranger)
}))

results <- do.call(rbind, parLapply(cl, 1:nrow(grid), function(i) {
  run_sim(grid$n[i], grid$gamma[i], grid$dgp[i], R = 100)
}))
stopCluster(cl)

mcols <- c("2SLS"       = "#D62728",
           "Post-Lasso" = "#FF7F0E",
           "DML-IV"     = "#2CA02C",
           "AIPW-IV"    = "#1F77B4",
           "IV-TMLE"    = "#9467BD")

mshapes <- c("2SLS"       = 16,
             "Post-Lasso" = 17,
             "DML-IV"     = 15,
             "AIPW-IV"    = 18,
             "IV-TMLE"    = 8)

dlabs <- c(linear = "Linear DGP", nonlinear = "Nonlinear DGP", sparse = "Sparse DGP")
glabs <- c("0.5" = "Weak\n(\u03b3=0.5)", "1" = "Moderate\n(\u03b3=1.0)",
           "2"   = "Strong\n(\u03b3=2.0)")

results_trim <- results |>
  filter(!(n == 500 & gamma == 0.5)) |>
  mutate(dgp    = factor(dlabs[dgp], levels = unname(dlabs)),
         method = factor(method,      levels = names(mcols)))

p_rmse <- results_trim |>
  filter(gamma == 2.0) |>
  ggplot(aes(n, rmse, color = method, shape = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 3.5) +
  facet_wrap(~dgp) +
  scale_color_manual(values = mcols) +
  scale_shape_manual(values = mshapes) +
  scale_x_continuous(breaks = c(500, 1000, 2000)) +
  labs(x = "Sample size (n)", y = "RMSE", color = NULL, shape = NULL,
       title = "RMSE by sample size  (\u03b3 = 2.0,  R = 100)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "gray93"),
        panel.grid.minor = element_blank())

min_cov <- floor(min(results$coverage, na.rm = TRUE) / 5) * 5

p_cov <- results |>
  filter(n == 2000) |>
  mutate(g_lab  = factor(glabs[as.character(gamma)], levels = unname(glabs)),
         dgp    = factor(dlabs[dgp],                 levels = unname(dlabs)),
         method = factor(method,                      levels = names(mcols))) |>
  ggplot(aes(g_lab, coverage, color = method, group = method, shape = method)) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "gray40") +
  geom_line(linewidth = 0.9) + geom_point(size = 3.5) +
  facet_wrap(~dgp) +
  scale_color_manual(values = mcols) +
  scale_shape_manual(values = mshapes) +
  scale_y_continuous(limits = c(min_cov, 103), breaks = seq(min_cov, 100, 10)) +
  labs(x = "Instrument strength", y = "Coverage (%)", color = NULL, shape = NULL,
       title = "95% CI coverage by instrument strength  (n = 2000,  R = 100)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "gray93"),
        panel.grid.minor = element_blank())

p_bias <- results_trim |>
  filter(gamma == 2.0) |>
  ggplot(aes(n, abs(bias), color = method, shape = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 3.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray60") +
  facet_wrap(~dgp) +
  scale_color_manual(values = mcols) +
  scale_shape_manual(values = mshapes) +
  scale_x_continuous(breaks = c(500, 1000, 2000)) +
  labs(x = "Sample size (n)", y = "|Bias|", color = NULL, shape = NULL,
       title = "Absolute bias by sample size  (\u03b3 = 2.0,  R = 100)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "gray93"),
        panel.grid.minor = element_blank())

p_ciwidth <- results |>
  filter(gamma == 2.0) |>
  mutate(dgp    = factor(dlabs[dgp], levels = unname(dlabs)),
         method = factor(method,      levels = names(mcols))) |>
  ggplot(aes(n, ci_width, color = method, shape = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 3.5) +
  facet_wrap(~dgp) +
  scale_color_manual(values = mcols) +
  scale_shape_manual(values = mshapes) +
  scale_x_continuous(breaks = c(500, 1000, 2000)) +
  labs(x = "Sample size (n)", y = "Mean 95% CI width",
       color = NULL, shape = NULL,
       title = "Mean 95% CI width by sample size  (\u03b3 = 2.0,  R = 100)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "gray93"),
        panel.grid.minor = element_blank())

print(p_rmse)
print(p_cov)
print(p_bias)
print(p_ciwidth)

