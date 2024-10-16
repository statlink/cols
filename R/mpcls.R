mpcls <- function(y, x) {

  dm <- dim(x)
  n <- dm[1]   ;   p <- dm[2]
  d <- dim(y)[2]
  dvec <-  2 * crossprod(x, y)

  xx <- 2 * crossprod(x)
  A <- cbind(1, diag(p) )
  bvec <- c(1, rep(0, p) )
  mse <- rep( sum(y^2), d)
  be <- matrix(nrow = p, ncol = d)

  for ( j in 1:d ) {
    f <- quadprog::solve.QP(Dmat = xx, dvec = dvec[, j], Amat = A, meq = 1, bvec = bvec)
    be[, j] <- f$solution
    mse[j] <- f$value
  }
  mse <- ( mse + Rfast::colsums(y^2) ) / n
  rownames(be) <- colnames(x)
  list(be = be, mse = mse)
}
