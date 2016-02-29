pls1gm.fit <- function (X, Y, ncomp, ...) {
  Y <- as.matrix(Y)
  dnX <- dimnames(X)[[2]]
  dnY <- dimnames(Y)
  nobj <- dim(X)[1]
  npred <- dim(X)[2]
  nresp <- dim(Y)[2]
  Xmeans <- colMeans(X)
  X <- X - rep(Xmeans, each = nobj)
  Ymean <- mean(Y)
  Y <- Y - Ymean
  A <- ncomp
  Tau <- matrix(0, nrow = nrow(X), ncol = A)
  TT <- matrix(0, nrow = nrow(X), ncol = A)
  V <- matrix(0, nrow = ncol(X), ncol = A)
  W <- matrix(0, ncol(X), A)
  P <- matrix(0, ncol(X), A)
  D2 <- matrix(0, nrow = ncol(X), ncol = A)
  iD2 <- matrix(0, nrow = ncol(X), ncol = A)
  nw <- matrix(nrow = 1, ncol = A)
  nt <- matrix(nrow = 1, ncol = A)
  iB <- matrix(0, ncol(X), A)
  iPreds <- matrix(0, nrow(X), A)
  iResids <- matrix(0, nrow(X), A)
  Q <- matrix(0, nrow = A, ncol = 1)
  tQ <- matrix(0, nrow = A, ncol = 1)
if(A == 1) {
  V[, 1] <- t(X) %*% Y
  nw[1] <- sqrt(c(crossprod(V[, 1]))) 
  W[, 1] <- V[, 1] / nw[1]
  Tau[, 1] <- X %*% W[, 1]
  nt[1] <- sqrt(c(crossprod(Tau[, 1]))) 
  TT[, 1] <- Tau[, 1] / nt[1]
  P[, 1] <- crossprod(X, TT[, 1])
  D2[1, 1] <- nt[1]
  iD2[1, 1] <- 1 / D2[1, 1]
  R <- W * iD2[1, 1]
  q1 <- crossprod(TT, Y)
  tQ[1, ] <- q1
  Y <- Y - TT[, 1] %*% t(q1)
  Betas <- R[, 1] %*% tQ
  row.names(Betas) <- dnX
  row.names(W) <- dnX
  fitted.pre <- X %*% Betas
  fitted <- fitted.pre + Ymean
  residuals <- Y - fitted.pre
  Yactual <- Y + Ymean
  Xdata <- as.data.frame(X)
  yloadings <- as.matrix((c(tQ) / diag(D2))[1])
  list(loadings = P, weights = W, D2 = D2[1:A, 1:A], iD2 = iD2[1:A, 1:A], 
       Ymean = Ymean,
       Xmeans = Xmeans, coefficients = Betas, y.loadings = yloadings, 
       scores = TT %*% D2[1, 1], R = R, Y =  Y, Yactual = Yactual, fitted = fitted, 
       residuals = residuals, Xdata = Xdata, iPreds = fitted, y.loadings2 = tQ)
} else {
  V[, 1] <- t(X) %*% Y
  nw[1] <- sqrt(c(crossprod(V[, 1])))
  W[, 1] <- V[, 1] / nw[1]
  Tau[, 1] <- X %*%W[, 1]
  nt[1] <- sqrt(c(crossprod(Tau[, 1]))) 
  TT[, 1] <- Tau[, 1] / nt[1]
  P[, 1] <- crossprod(X, TT[, 1])
  q1 <- crossprod(TT[, 1], Y)
  tQ[1, ] <- q1
  out <- NULL
  for(a in 2:A) {
    V[, a] <- nw[a - 1] * (W[, a - 1] - (P[, a - 1] / nt[a - 1]))
    V[, a] <- V[, a] - (W[, 1:(a - 1)] %*% (t(W[, 1:(a - 1)]) %*% V[, a]))
    nw[a] <- sqrt(c(crossprod(V[, a]))) 
    W[, a] <- V[, a] / nw[a]
    Tau[, a] <- X %*% W[, a]
    Tau[, a] <- Tau[, a] - (TT[, 1:(a - 1)] %*% (crossprod(TT[, 1:(a - 1)], as.matrix(Tau[, a]))))
    nt[a] <- sqrt(c(crossprod(Tau[, a]))) 
    TT[, a] <- Tau[, a] / nt[a]
    P[, a] <- crossprod(X, TT[, a])
    q1 <- crossprod(TT[, a], Y)
    tQ[a, ] <- q1
    out <- list(W, D2[1:A, 1:A], TT, iD2[1:A, 1:A], tQ, P)
  }
  W <- out[[1]]
  D2 <- as.matrix(out[[2]])
  TT <- out[[3]]
  iD2 <- as.matrix(out[[4]])
  tQ <- as.matrix(out[[5]])
  P <- crossprod(X, TT)
  d2 <- -(nt[1:a-1] * (nw[2:A] / nw[1:a-1]))
  D2 <- diag(as.vector(nt))
  D2[seq((A + 1), (A + 1) * (A - 1), (A + 1))] <- d2
  iD2 <- backsolve(D2, diag(A))
  R <- W %*% iD2
  q1 <- crossprod(TT, Y)
  D2.b <- D2
  for(a in 1:A) {
    P[, a] <- P[, a] / D2.b[a, a]
  }
  R <- W %*% iD2
  if(A == 1) {
    iB[, 1] <- W[, 1] * iD2[1, 1] * tQ[1, 1]
    out.iB <- iB
  } else {
    for(a in 2:A) {
      iB[, 1] <- W[, 1] * iD2[1, 1] * tQ[1, 1]
      iB[, a] <- R[, 1:a] %*% tQ[1:a, 1]
      out.iB <- iB
    }
  }
  Betas <- R %*% q1
  fitted <- X %*% Betas
  residuals <- Y - fitted
  Yactual <- Y + Ymean
  Xdata <- as.data.frame(X)
  yloadings <- as.matrix(c(tQ) / diag(D2))
  if(A == 1) {
    iPreds[, 1] <- X %*% W[, 1] * iD2[1, 1] * tQ[1, 1] + Ymean
    out.iPreds <- iPreds
  } else {
    for(a in 2:A) {
      iPreds[, 1] <- X %*% W[, 1] * iD2[1, 1] * tQ[1, 1] + Ymean
      iPreds[, a] <- X %*% R[, 1:a] %*% tQ[1:a, 1] + Ymean
      out.iPreds <- iPreds
    }
  }
  if(A == 1) {
    iResids[, 1] <- Y - (X %*% out.iB[, 1])
    out.iResids <- iResids
  } else {
    for(a in 2:A) {
      iResids[, 1] <- Y - (X %*% out.iB[, 1])
      iResids[, a] <- Y - (X %*% out.iB)[, a]
      out.iResids <- iResids
    }
  }
  colnames(out.iB) <- paste("ncomp", 1:A, sep = ".")
  colnames(out.iPreds) <- paste("ncomp", 1:A, sep = ".")
  colnames(out.iResids) <- paste("ncomp", 1:A, sep = ".")
  colnames(R) <- paste("ncomp", 1:A, sep = ".")
  class(TT) <- "scores"
  class(P) <- class(q1) <- class(yloadings) <- "loadings"
  class(W) <- "weights"
  class(out.iPreds) <- "predict"
  class(out.iB) <- "coefficients"
  class(Xmeans) <- "vector"
  row.names(P) <- dnX
  row.names(W) <- dnX
  row.names(out.iB) <- dnX
  row.names(Betas) <- dnX
  list(loadings = P, weights = W, D2 = D2, iD2 = iD2, Ymean = Ymean,
       Xmeans = Xmeans, coefficients = out.iB, y.loadings = yloadings, 
       scores = TT %*% diag(diag(D2)), R = R, Y = Y, Yactual = Yactual, fitted = fitted, 
       residuals = out.iResids, Xdata = Xdata, iPreds = out.iPreds, y.loadings2 = tQ)
  }
}
