#' @name summary.impactspar.pspatreg
#' @rdname summary.impactspar.pspatreg
#'
#' @title Summary method for object of class impactspar.pspatreg.
#'
#' @description This method summarizes direct, indirect and total effects (or impacts)
#'   for continous parametric covariates in semiparametric spatial regression models.
#'
#' @param object \emph{impactspar} object fitted using \code{\link{pspatfit}} function.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class \emph{summary.impactspar.pspatreg}
#'
#' @author Roman Minguez \email{roman.minguez@@uclm.es}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{impactspar}} Compute direct, indirect and total
#'     effects (or impacts)
#'     for continous parametric covariates in PS-SAR, PS-SEM, PS-SARAR, PS-SLX or PS-DURBIN regression models.
#'   \item \code{\link{print.summary.impactspar.pspatreg}} print objects of
#'     class \emph{summary.pspatreg}
#' }
#'
#' @examples
#'  # See examples for \code{\link{impactspar}} function.
#' @export
summary.impactspar.pspatreg <- function(object, ...) {
 z <- object
if (z$type %in% c("sdm", "sar", "sarar")) {
  tot <- z$mimpactstot
  dir <- z$mimpactsdir
  ind <- z$mimpactsind
  varpar <- rownames(tot)
  nrep <- ncol(tot)
  mean_dir <- apply(dir, 1, mean)
  mean_tot <- apply(tot, 1, mean)
  mean_ind <- apply(ind, 1, mean)
  sd_dir <- apply(dir, 1, sd)
  sd_tot <- apply(tot, 1, sd)
  sd_ind <- apply(ind, 1, sd)
  t_dir <- mean_dir / sd_dir
  t_tot <- mean_tot / sd_tot
  t_ind <- mean_ind / sd_ind
  z$tot_table <- cbind(mean_tot, sd_tot, t_tot,
                       2 * pnorm(abs(t_tot), 
                                 mean = 0, sd = 1, 
                                 lower.tail = FALSE))
  colnames(z$tot_table) <- c("Estimate", "Std. Error", 
                             "t value", "Pr(>|t|)")
  rownames(z$tot_table) <- varpar
  z$dir_table <- cbind(mean_dir, sd_dir, t_dir,
                       2 * pnorm(abs(t_dir), 
                                 mean = 0, sd = 1,
                                 lower.tail = FALSE))
  colnames(z$dir_table) <- c("Estimate", "Std. Error", 
                             "t value", "Pr(>|t|)")
  rownames(z$dir_table) <- varpar
  z$ind_table <- cbind(mean_ind, sd_ind, t_ind,
                       2 * pnorm(abs(t_ind), 
                                 mean = 0, sd = 1,
                                 lower.tail = FALSE))
  colnames(z$ind_table) <- c("Estimate", "Std. Error", 
                             "t value", "Pr(>|t|)")
  rownames(z$ind_table) <- varpar
 }
 class(z) <- c("summary.impactspar.pspatreg", class(z))
 z
}
