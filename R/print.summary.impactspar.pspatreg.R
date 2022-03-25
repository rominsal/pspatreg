#' @name print.summary.impactspar.pspatreg
#' @rdname print.summary.impactspar.pspatreg
#'
#' @title Print method for objects of class summary.impactspar.pspatreg
#'
#' @param x object of class \emph{summary.impactspar.pspatreg}.
#' @param digits number of digits to show in printed tables.
#'   Default: max(3L, getOption("digits") - 3L).
#' @param ... further arguments passed to or from other methods.
#'
#' @author 
#' \tabular{ll}{ 
#'   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Roberto Basile \tab \email{roberto.basile@@univaq.it} \cr Maria Durban \tab
#'   \email{mdurban@@est-econ.uc3m.es} \cr Gonzalo Espana-Heredia \tab
#'   \email{gehllanza@@gmail.com} \cr 
#'  }
#'  
#' @seealso
#' \itemize{
#'   \item \code{\link{impactspar}} Compute direct, indirect and
#'     total effects (or impacts).
#'     for continous parametric covariates in PS-SAR, PS-SEM, PS-SARAR, PS-SLX or PS-DURBIN regression models.
#'   \item \code{\link{summary.impactspar.pspatreg}} Summary method
#'     for \emph{summary.pspatreg} objects.
#' }
#'
#' @examples
#' # See examples for \code{\link{impactspar}} function.
#' @export
print.summary.impactspar.pspatreg <- function(x, 
      digits = max(3L, getOption("digits") - 3L), ...) {
  if (x$type %in% c("sar", "sdm", "sarar")) {
    if(!is.null(x$tot_table)) {
      cat(paste("\n Total Parametric Impacts (",x$type,") \n", sep = ""))
      printCoefmat( x$tot_table, P.values = FALSE, has.Pvalue = FALSE)
    }
    if(!is.null(x$dir_table)) {
      cat(paste("\n Direct Parametric Impacts (",x$type,") \n", sep = ""))
      printCoefmat( x$dir_table, P.values = FALSE, has.Pvalue = FALSE)
    }
    if(!is.null(x$ind_table)) {
      cat(paste("\n Indirect Parametric Impacts (",x$type,") \n", sep = ""))
      printCoefmat( x$ind_table, P.values = FALSE, has.Pvalue = FALSE)
    }
  }
  if (x$type %in% c("slx", "sdem")) {
    cat(paste("\n Parametric Impacts (",x$type,") \n", sep = ""))
    print(x$mimpacts, digits = digits)
    cat("\n Standard errors: \n")
    print(x$semimpacts, digits = digits)
    cat("\n Z-values: \n")
    print(x$zvalmimpacts, digits = digits)
    cat("\n p-values: \n")
    print(x$pvalmimpacts, digits = digits)
  }
  invisible(x)
}
