## NOTE: this is a copy of ../R/pca_compress.R so that the Shiny app directory is
## self-contained when it is deployed to shinyapps.io. Keep the two in sync.

## SDS 365 - Homework 2
## PCA-based image compression

## Rank-k PCA approximation of an image matrix.
##
## X: an n x p numeric matrix (the image)
## k: rank of the approximation (1 <= k <= p)
##
## Returns a list with
##   approx: the rank-k approximation X^(k) = z_1 u_1' + ... + z_k u_k'
##   error : ||X - X^(k)||_F, the Frobenius norm of the approximation error
pca_compress <- function(X, k) {
  X <- as.matrix(X)
  p <- ncol(X)

  if (!is.numeric(k) || length(k) != 1 || k < 1 || k > p) {
    stop("k must be a single number between 1 and ncol(X) = ", p)
  }

  V <- crossprod(X)                  # V = X'X  (p x p)
  eig <- eigen(V, symmetric = TRUE)  # V = U Lambda U'
  U <- eig$vectors                   # p x p eigenvectors, decreasing eigenvalue
  Z <- X %*% U                       # n x p principal component scores

  idx <- seq_len(k)
  Xk <- Z[, idx, drop = FALSE] %*% t(U[, idx, drop = FALSE])

  list(approx = Xk, error = norm(X - Xk, type = "F"))
}

## Same decomposition, but computed once so that many values of k can be
## evaluated without repeating the eigendecomposition. Used by the analysis
## script and the Shiny app; pca_compress() above is the stand-alone answer.
pca_decompose <- function(X) {
  X <- as.matrix(X)
  eig <- eigen(crossprod(X), symmetric = TRUE)
  list(X = X, U = eig$vectors, Z = X %*% eig$vectors, values = eig$values)
}

pca_rank_k <- function(fit, k) {
  idx <- seq_len(k)
  Xk <- fit$Z[, idx, drop = FALSE] %*% t(fit$U[, idx, drop = FALSE])
  list(approx = Xk, error = norm(fit$X - Xk, type = "F"))
}
