nma_additive <- function(TE, W,
                         treat1, treat2,
                         studlab,
                         level,
                         X, C.matrix, B.matrix,
                         df.Q.additive, n, sep.trts) {
  
  
  m <- length(TE)
  #
  # Drop Matrix attributes
  #
  W <- as.matrix(W)
  class(W) <- "matrix"
  #
  # Laplacian matrix and pseudoinverse of L
  #
  L <- t(X) %*% W %*% X
  ##
  Lplus <- ginv(L) # = Cov matrix of beta (components)
  Lplus[is_zero(Lplus)] <- 0
  colnames(Lplus) <- colnames(L)
  rownames(Lplus) <- rownames(L)
  #
  # H matrix
  #
  H <- X %*% Lplus %*% t(X) %*% W
  H[is_zero(H)] <- 0
  
  
  ##
  ## beta = effects of components
  ##
  beta <- as.vector(Lplus %*% t(X) %*% W %*% TE)
  se.beta <- sqrt(diag(Lplus))
  names(beta) <- names(se.beta)
  ##
  ## theta = estimates for combination treatments
  ##
  theta <- as.vector(C.matrix %*% beta)
  se.theta <- sqrt(diag(C.matrix %*% Lplus %*% t(C.matrix)))
  names(theta) <- names(se.theta)
  ##
  ## delta = treatment estimates for observed comparisons
  ##
  delta <- as.vector(X %*% beta) # = B.matrix %*% theta = H %*% TE
  se.delta <- unname(sqrt(diag(X %*% Lplus %*% t(X))))
  ##
  ## delta.all(.matrix) = all direct and indirect treatment estimates
  ##
  B.full <- createB(ncol = n)
  X.all <- B.full %*% C.matrix
  colnames(X.all) <- colnames(C.matrix)
  ##
  labels <- colnames(B.matrix)
  ##
  k <- 0
  lab <- vector(mode = "numeric", length = choose(n, 2))
  ##
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      k <- k + 1
      lab[k] <- paste(labels[i], labels[j], sep = sep.trts)
    }
  }
  ##
  rownames(X.all) <- lab
  #
  # Variance-covariance matrix for all comparisons
  #
  Cov <- X.all %*% Lplus %*% t(X.all)
  Cov[is_zero(Cov)] <- 0
  #
  delta.all <- as.vector(X.all %*% beta)
  se.delta.all <- sqrt(diag(X.all %*% Lplus %*% t(X.all)))
  names(delta.all) <- names(se.delta.all)
  ##
  delta.all.matrix <- se.delta.all.matrix <- matrix(0, ncol = n, nrow = n)
  ##
  k <- 0
  ##
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      k <- k + 1
      delta.all.matrix[i, j] <-  delta.all[k]
      delta.all.matrix[j, i] <- -delta.all[k]
      se.delta.all.matrix[i, j] <-
        se.delta.all.matrix[j, i] <- se.delta.all[k]
    }
  }
  ##
  colnames(delta.all.matrix) <- rownames(delta.all.matrix) <-
    colnames(se.delta.all.matrix) <- rownames(se.delta.all.matrix) <- labels
  
  
  comparisons <- c(list(studlab = studlab, treat1 = treat1, treat2 = treat2),
                   ci(delta, se.delta, level = level))
  ##
  all.comparisons <- ci(delta.all.matrix, se.delta.all.matrix,
                        level = level)
  ##
  components <- ci(beta, se.beta, level = level)
  ##
  combinations <- ci(theta, se.theta, level = level)
  ##
  ## Test of total heterogeneity / inconsistency:
  ##
  Q.additive <- as.vector(t(delta - TE) %*% W %*% (delta - TE))
  ##
  if (is.na(df.Q.additive) | df.Q.additive == 0)
    pval.Q.additive <- NA
  else
    pval.Q.additive <- 1 - pchisq(Q.additive, df.Q.additive)
  ##
  ## Heterogeneity variance
  ##
  I <- diag(m)
  E <- matrix(0, nrow = m, ncol = m)
  for (i in 1:m)
    for (j in 1:m)
      E[i, j] <- as.numeric(studlab[i] ==  studlab[j])
  ##
  if (df.Q.additive == 0) {
    tau2 <- NA
    tau <- NA
    I2 <- NA
    lower.I2 <- NA
    upper.I2 <- NA
  }
  else {
    tau2 <- max(0, (Q.additive - df.Q.additive) /
                   sum(diag((I - H) %*%
                            (B.matrix %*% t(B.matrix) * E / 2) %*% W)))
    tau <- sqrt(tau2)
    ci.I2 <- isquared(Q.additive, df.Q.additive, level)
    I2 <- ci.I2$TE
    lower.I2 <- ci.I2$lower
    upper.I2 <- ci.I2$upper
  }
  
  H.matrix <- H
  #
  rownames(H.matrix) <- colnames(H.matrix) <- studlab
  #
  res <- list(comparisons = comparisons,
              all.comparisons = all.comparisons,
              components = components,
              combinations = combinations,
              ##
              Q.additive = Q.additive,
              df.Q.additive = df.Q.additive,
              pval.Q.additive = pval.Q.additive,
              ##
              tau = tau,
              I2 = I2, lower.I2 = lower.I2, upper.I2 = upper.I2,
              ##
              L.matrix = L,
              Lplus.matrix = Lplus,
              H.matrix = H.matrix,
              #
              Cov = Cov
              )
  
  
  res
}
