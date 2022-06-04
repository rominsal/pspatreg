MM_basis <- function (x, xl, xr, ndx, bdeg, pord, decom = 2) {
	Bb <- bspline(x,xl,xr,ndx,bdeg)
	knots <- Bb$knots
	B <- Bb$B
	m <- ncol(B)
	n <- nrow(B)
	D <- diff(diag(m), differences=pord)
	P.svd <- svd(crossprod(D))
	U <- (P.svd$u)[,1:(m-pord)] # eigenvectors
	d <- (P.svd$d)[1:(m-pord)]  # eigenvalues
	Z <- B %*% U
	if(decom == 1) {
		X <- B %*% ((P.svd$u)[,-(1:(m-pord))])
	} else if (decom == 2){
		X <- NULL
		for(i in 0:(pord-1)){ X <- cbind(X,x^i) }
	} else if(decom == 3) {
		Xf <- NULL
		for(i in 0:(pord-1)){
	    	Xf <- cbind(Xf,knots[-c((1:pord),
	    	                  (length(knots)- pord + 1):length(knots))]^i)
		}
		X <- B %*% Xf
	}
	list(X = X, Z = Z, d = d, B = B, m = m, D = D, U = P.svd$u, knots = knots)
}
