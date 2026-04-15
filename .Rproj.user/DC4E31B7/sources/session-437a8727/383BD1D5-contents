cls <- function(y, x, R, ca) {
  mod <- Rfast2::cls(y, x, R, ca)
  be <- mod$bcls
  mse <- sum( (y - x %*% be)^2 ) / dim(x)[1]

  list(be = be, mse = mse)
}
