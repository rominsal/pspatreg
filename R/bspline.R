bspline <- function(x, xl, xr, ndx, bdeg){
	dx <- (xr - xl)/ndx
	knots <- seq(xl - bdeg*dx, xr + bdeg*dx, by=dx)
	B <- spline.des(knots, x, bdeg + 1, 0*x)$design
	res <- list(B = B, knots = knots)
	res
}
