construct_block <- function(A1, A2, A3, A4) {
    block <- rbind(cbind(A1, A2), cbind(A3, A4))
    return(block)
}

#####################################################
construct_matrices <- function(X, Z, y, GLAM = FALSE, w = NULL) {
  if(is.matrix(X)) X <- Matrix(X)
  if(is.matrix(Z)) Z <- Matrix(Z)
  if(is.matrix(y)) y <- Matrix(y)
  if(is.numeric(y)) y <- Matrix(as.matrix(y))

  # if(GLAM) {
  #   XtX <- XtX(X,w)
  #   XtZ <- XtZ(X,Z,w)
  #   ZtX <- t(XtZ)
  #   ZtZ <- ZtZ(Z,w)
  #   Zty = Zty(Z,z,w)
  #   Xty = Xty(X,z,w)
  #   yty <- sum((z^2)*w)
  #   ZtXtZ = rbind(XtZ, ZtZ)
  #   u <- c(Xty,Zty)
  # } else {
  #XtW <- t(X * w)
  XtX <-crossprod(X)
  Xty <-crossprod(X, y)
  XtZ <-crossprod(X, Z)
  ZtX <- t(XtZ)
  #ZtW <- t(Z * w)
  ZtZ <-crossprod(Z)
  Zty <-crossprod(Z, y)
  ZtXtZ <- rbind(XtZ, ZtZ)
  u <- rbind(Xty, Zty)
  #yty <- sum((y^2)*w)
  yty <- as.numeric(crossprod(y))
  #}
  res <- list(XtX = XtX, XtZ = XtZ, ZtX = ZtX, ZtZ = ZtZ,
              Xty = Xty, Zty = Zty, yty = yty, ZtXtZ = ZtXtZ, u = u)
}

#########################################################

H <- function(X,A) {
    d <- dim(A)
    M <- matrix(A, nrow = d[1])
    XM <- X%*%M
    array(XM, c(nrow(XM),d[-1]))
}
#########################################################

Rotate <- function(A) {
    d <- 1:length(dim(A))
    d1 <- c(d[-1],d[1])
    aperm(A, d1)
}
#########################################################

RH <- function(X,A) { Rotate(H(X,A)) }
#########################################################

Rten2 <- function(X1,X2) {
        one.1 <- matrix(1,1,ncol(X1))
        one.2 <- matrix(1,1,ncol(X2))
        kronecker(X1,one.2)*kronecker(one.1,X2)
}
#########################################################

#########################################################
# construct_matrices <- function(X, Z, z, w, GLAM) {
#   if(GLAM) {
#
#     XtX <- XtX(X,w)
#     XtZ <- XtZ(X,Z,w)
#     ZtX <- t(XtZ)
#     ZtZ <- ZtZ(Z,w)
#     Zty = Zty(Z,z,w)
#     Xty = Xty(X,z,w)
#     yty <- sum((z^2)*w)
#     ZtXtZ = rbind(XtZ, ZtZ)
#     u <- c(Xty,Zty)
#   } else {
#     XtW = t(X * w)
#     XtX = XtW %*% X
#     Xty = XtW%*%z
#     XtZ = XtW %*% Z
#     ZtX = t(XtZ)
#     ZtW = t(Z * w)
#     ZtZ = ZtW %*% Z
#     Zty = ZtW %*% z
#     ZtXtZ = rbind(XtZ, ZtZ)
#     u <- c(Xty,Zty)
#     yty <- sum((z^2)*w)
#   }
#   res <- list(XtX = XtX, XtZ = XtZ, ZtX = ZtX, ZtZ = ZtZ,
#               Xty = Xty, Zty = Zty, yty = yty, ZtXtZ = ZtXtZ, u = u)
# }
