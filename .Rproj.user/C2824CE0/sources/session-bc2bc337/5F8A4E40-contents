mcls <- function(y, x, R, ca) {

  dm <- dim(x)
  n <- dm[1]  ;  p <- dm[2]
  d <- dim(y)[2]
  mse <- rep( sum(y^2), d)

  xxs <- solve( crossprod(x) )
  bols <- xxs %*% crossprod(x, y)
  be <- bols - xxs %*% R %*% solve(R %*% xxs %*% R, R %*% bols - ca)
  mse <- Rfast::colmeans( (y - x %*% be)^2 )

  list(be = be, mse = mse)
}
