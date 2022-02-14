#' @name methods_pspatreg
#' @title Methods for class pspatreg
#' @description The \code{anova()} function provides tables of fitted 
#'  spsur models including information criteria (AIC and BIC), 
#'  log-likelihood and degrees of freedom of each fitted model. The 
#'  argument \code{lrtest} allows to perform LR tests between nested models.
#'  The \code{print()} function is used to print short tables including the 
#'  values of beta and spatial coefficients as well as p-values of significance test for each 
#'  coefficient. This can be used as an alternative to 
#'  \code{\link{summary.pspatreg}} when a brief output is needed. 
#'  The rest of methods works in the usual way. 
#'       
#' @param object a \code{pspatreg} object created by 
#'   \code{\link{pspatfit}}.
#' @param x similar to \code{object} argument for \code{print()} 
#'  and \code{plot} functions.
#' @param digits number of digits to show in printed tables.
#'   Default: max(3L, getOption("digits") - 3L).
#' @param lrtest logical value to compute likelihood ratio
#'   test for nested models in `anova` method. Default = \code{TRUE}
#' @param REML logical value to get restricted log-likelihood 
#'   instead of the usual log-likelihood. Default = \code{FALSE}
#' @param bayesian logical value to get bayesian or frequentist  
#'   covariance matrix for parametric terms. Default = \code{FALSE}
#' @param ... further arguments passed to or from other methods.  
#' @examples
#' ## INCLUDE SOME EXAMPLES OF ANOVA AND PERHAPS ANOTHER 
#' ## ADDITIONAL METHODS
#' @author
#'   \tabular{ll}{
#'   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   }

NULL

#' @name anova
#' @rdname methods_pspatreg
#' @export
#' 
anova.pspatreg <- function(object, ..., lrtest = TRUE) {
  object <- list(object, ...)
  cl <- match.call()
  ns <- sapply(object, function(x) length(x$residuals))
  if (any(ns != ns[1L]))
    stop("models were not all fitted to the same size of dataset")
  nmodels <- length(object)
  vtypes <- character(length = nmodels)
  vLL <- vRLL <- vAIC <- vBIC <- vdf <- numeric(nmodels) 
  vresdf <- vresdev <- vsigmasq <- numeric(nmodels)
  if (lrtest) {
    cat("Warning: LR test should be used only for nested models.\n")
    vlrtest <- vpval <- vector(mode = "numeric", 
                               length = nmodels)
    vlrtest[1] <- vpval[1] <- NA
  }
  for (i in 1:nmodels) {
    if (!inherits(object[[i]], "pspatreg"))
      stop("object not a fitted pspatreg model")
    vresdf[i] <- object[[i]]$df.residual
    vsigmasq[i] <- sum(object[[i]]$residuals^2) / vresdf[i]
    llik_i <- logLik(object[[i]], REML = FALSE)
    rllik_i <- logLik(object[[i]], REML = TRUE)
    vLL[i] <- as.numeric(llik_i)
    vRLL[i] <- as.numeric(rllik_i)
    vresdev[i] <- -2*vLL[i]
    vAIC[i] <- as.numeric(AIC(llik_i))
    vBIC[i] <- as.numeric(BIC(rllik_i))
    vdf[i] <- attr(llik_i, "df")
    vtypes[i] <- as.character(cl[[i + 1]]) #      paste("model ", i, sep = "")
    if (lrtest) {
      if(i>1) {
       vlrtest[i] <- 2*(vRLL[i] - vRLL[i-1])
       vpval[i] <- pchisq(vlrtest[i],
                          df = vdf[i] - vdf[i-1],
                          lower.tail = FALSE) 
      }
    }
    rm(llik_i, rllik_i)
  } 
  if (lrtest) {
       res <- data.frame(logLik = vLL,
                         rlogLik = vRLL,
                         edf = vdf,
                         AIC = vAIC, 
                         BIC = vBIC,
                         LRtest = vlrtest,
                         `p-val` = vpval,
                         row.names = vtypes)
     } else {
       res <- data.frame(`logLik` = vLL,
                         rlogLik = vRLL,
                         edf = vdf,
                         AIC = vAIC, 
                         BIC = vBIC,
                         row.names = vtypes)
     }
    class(res) <- c("anova", "data.frame")
    return(res)
}

#' @name coef
#' @rdname methods_pspatreg
#' @export
coef.pspatreg <- function(object, ...) {
  ret <- NULL
  if (!is.null(object$rho))
    ret <- c(ret, object$rho)
  if (!is.null(object$delta))
    ret <- c(ret, object$delta)
  if (!is.null(object$phi))
    ret <- c(ret, object$phi)
  ret <- c(ret, object$bfixed)
  names(ret) <- gsub("fixed_", "", names(ret))
  ret
}

#' @name fitted
#' @rdname methods_pspatreg
#' @export
fitted.pspatreg <- function(object, ...)
{
  if (is.null(object$na.action))
    res <- object$fitted.values
  else res <- napredict(object$na.action, object$fitted.values)
  res
}

#' @name logLik
#' @rdname methods_pspatreg
#' @export
logLik.pspatreg <- function(object, ..., REML = FALSE) {
  edftot <- object$edftot
  N <- length(object$residuals)
  N0 <- N
  if (REML) {
    N <- N - edftot
    LL <- object$llik_reml
  } else LL <- object$llik
  attr(LL, "nall") <- N0
  attr(LL, "nobs") <- N
  attr(LL, "df") <- edftot
  class(LL) <- "logLik"
  LL
}

#' @name residuals
#' @rdname methods_pspatreg
#' @export
residuals.pspatreg <- function(object, ...) {
  if (is.null(object$na.action))
    res <- object$residuals
  else res <- napredict(object$na.action, object$residuals)
  res
}

#' @name vcov
#' @rdname methods_pspatreg
#' @export
vcov.pspatreg <- function(object, ..., bayesian = FALSE) {
  if (!bayesian) 
    res <- as.matrix(object$vcov_fr)
  else res <- as.matrix(object$vcov_by)
  idxsigma <- grepl("fixed", rownames(res))
  res <- res[idxsigma, idxsigma]
  rownames(res) <- gsub("fixed_", "", rownames(res))
  colnames(res) <- gsub("fixed_", "", colnames(res))
  res
}

#' @name print
#' @rdname methods_pspatreg
#' @export
print.pspatreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                        ...) {
  if (!inherits(x, "pspatreg")) 
    stop("Argument must be a pspatreg object")
  summx <- summary(x)
  mtablepar <- summx$coef_par_table
  print(round(mtablepar, digits = digits))
  invisible(x)
}


