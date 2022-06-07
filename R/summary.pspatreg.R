#' @name summary.pspatreg
#' @rdname summary.pspatreg
#'
#' @title Summary method for objects of class pspatreg.
#'
#' @description This method summarizes both spatial (2-dimension) and
#'   spatio-temporal (3-dimension) \emph{pspatreg} objects.
#'   The tables include information of:
#'   \itemize{
#'      \item The spatial (or spatio-temporal) trends. When the model is ANOVA
#'        the trend is decomposed in main and interaction effects.
#'      \item The parametric and non-parametric covariates.
#'      \item The \eqn{\rho} parameter when the model is SAR.
#'      \item The \eqn{\phi} parameter when the model is spatio-temporal
#'        with a first-order autorregressive in the noise.
#'  }
#'
#' @param object \emph{pspatreg} object fitted using \code{\link[pspatreg]{pspatfit}} function.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class \emph{summary.pspatreg}
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
#'   \item \code{\link{pspatfit}} estimate spatial or spatio-temporal semiparametric 
#'     regression models.
#'   \item \code{\link{print.summary.pspatreg}} print objects of class \emph{summary.pspatreg}
#' }
#'
#' @examples
#'  # See examples for \code{\link{pspatfit}} function.
#'
#' @export
summary.pspatreg <- function(object,...) {
 z <- object 
 names_var <- labels(z$terms)
 names_varspt <- names_var[grepl("pspt", names_var)]
 nvarspt <- length(names_varspt)
 #names_varnopar <- names_var[grepl("pspl", names_var)]
 names(z$bfixed) <- gsub("fixed_", "", names(z$bfixed))
 names(z$se_bfixed) <- gsub("fixed_", "", names(z$se_bfixed))
 names_varpar <- names(z$bfixed) 
 names_varpar <- names_varpar[!grepl("pspl", names_varpar) & 
                              !grepl("pspt", names_varpar)]
 nvarpar <- length(names_varpar)
 names_varnopar <- names(z$bfixed)
 names_varnopar <- names_varnopar[grepl("pspl", names_varnopar)]
 names_varnopar <- gsub("fixed_", "", names_varnopar)
 names_varnopar <- gsub(").1", ")", names_varnopar)
 rdf <- z$df.residual
 r <- z$residuals
 f <- z$fitted.values
 rss <- sum(r^2)
 resvar <- rss/rdf
 #names_varpar <- gsub("fixed_", "", names_varpar)
 # match_names_varpar <- unique(grep(paste(c("spt", "_main", 
 #                                           "_int"),
 #                                   collapse = "|"),
 #                              names_varpar, value = TRUE))
 # names_varpar <- names_varpar[ !(names_varpar %in% 
 #                                    match_names_varpar)]
 # if ("Intercept" %in% names_varpar)
 #   names_varpar <- c("Intercept", names_varpar[!(names_varpar 
 #                                                %in% "Intercept")])
 # est_par <- z$bfixed[names_varpar]
 # se_par <- z$se_bfixed[names_varpar]
 est_par <- z$bfixed
 se_par <- z$se_bfixed
 tval_par <- est_par/se_par
 z$coef_par_table <- cbind(est_par, se_par, tval_par,
                      2 * pt(abs(tval_par),rdf, lower.tail = FALSE))
 colnames(z$coef_par_table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
 #rownames(z$coef_par_table) <- names_varpar
 if(!(is.null(z$rho))) {
   t_rho <- z$rho / z$se_rho
   pval_rho <- 2 * pt(abs(t_rho), rdf, lower.tail = FALSE)
   z$coef_par_table <- rbind(z$coef_par_table,
                         c(z$rho, z$se_rho, t_rho, pval_rho))
   n <- length(rownames(z$coef_par_table))
   rownames(z$coef_par_table)[n] <- c("rho")
 }
 if(!(is.null(z$delta))) {
    t_delta <- z$delta / z$se_delta
    pval_delta <- 2 * pt(abs(t_delta), rdf, lower.tail = FALSE)
    z$coef_par_table <- rbind(z$coef_par_table,
                              c(z$delta, z$se_delta, 
                                t_delta, pval_delta))
    n <- length(rownames(z$coef_par_table))
    rownames(z$coef_par_table)[n] <- c("delta")
 }
 if (!is.null(z$time)){
   if(z$cor == "ar1")
   {
     t_phi <- z$phi / z$se_phi
     pval_phi <- 2 * pt(abs(t_phi),rdf, lower.tail = FALSE)
     z$coef_par_table <- rbind(z$coef_par_table,
                               c(z$phi,z$se_phi,t_phi,pval_phi))
     n <- length(rownames(z$coef_par_table))
     rownames(z$coef_par_table)[n] <- c("phi")
   }
 }
 if (!is.null(z$edfspt)) {
   if (z$psanova) { # is.null(time)
     if (is.null(z$time)) {
       z$coef_spttrend_table <- matrix(0, nrow = 3 , ncol = 1)
       z$coef_spttrend_table[1, 1] <- z$edfspt[c("f1_main")]
       z$coef_spttrend_table[2, 1] <- z$edfspt[c("f2_main")]
       z$coef_spttrend_table[3, 1] <- sum(z$edfspt[c("f12.1",
                                                     "f12.2")])
       colnames(z$coef_spttrend_table) <- c("EDF")
       rownames(z$coef_spttrend_table) <- c("f1", "f2", "f12")
     } else { # !is.null(time)
       z$coef_spttrend_table <- matrix(0,nrow = 7, ncol = 1)
       z$coef_spttrend_table[1, 1] <- z$edfspt[c("f1_main")]
       z$coef_spttrend_table[2, 1] <- z$edfspt[c("f2_main")]
       z$coef_spttrend_table[3, 1] <- z$edfspt[c("ft_main")]
       z$coef_spttrend_table[4, 1] <- sum(z$edfspt[c("f12.1", 
                                                     "f12.2")])
       z$coef_spttrend_table[5, 1] <- sum(z$edfspt[c("f1t.1", 
                                                     "f1t.2")])
       z$coef_spttrend_table[6, 1] <- sum(z$edfspt[c("f2t.1", 
                                                     "f2t.2")])
       z$coef_spttrend_table[7, 1] <- sum(z$edfspt[c("f12t.1", 
                                                       "f12t.2", 
                                                       "f12t.3")])
       colnames(z$coef_spttrend_table) <- c("EDF")
       rownames(z$coef_spttrend_table) <- c("f1", "f2", "ft", 
                                            "f12", "f1t", "f2t", 
                                            "f12t")
     }
   } else { # PSANOVA = FALSE
     if (is.null(z$time)){
       z$coef_spttrend_table <- matrix(0, nrow = 1, ncol = 1)
       z$coef_spttrend_table[1, 1] <- sum(z$edfspt[c("sp1", "sp2")])
       colnames(z$coef_spttrend_table) <- c("EDF")
       rownames(z$coef_spttrend_table) <- c("f(sp1, sp2)")
     } else {
       z$coef_spttrend_table <- matrix(0, nrow = 1, ncol = 1)
       z$coef_spttrend_table[1, 1] <- sum(z$edfspt[c("sp1", "sp2",
                                                     "time")])

       colnames(z$coef_spttrend_table) <- c("EDF")
       rownames(z$coef_spttrend_table) <- c("f(sp1,sp2,time)")
     }
   }
 } else { z$coef_spttrend_table <- NULL}

 if( !is.null(z$coef_spttrend_table)){
   if (any(is.na(z$coef_spttrend_table))) {
     index_na_table <- which(is.na(z$coef_spttrend_table))
     z$coef_spttrend_table[c(index_na_table)] <- 0
   }
 }
 if(!is.null(z$edfnopar)) {
   z$coef_nopar_table <- matrix(z$edfnopar, ncol = 1)
   colnames(z$coef_nopar_table) <- c("EDF")
   rownames(z$coef_nopar_table) <- names_varnopar
 }
 z$sigma <- sqrt(resvar)
 z$edftot <- z$edftot
 class(z) <- c("summary.pspatreg",class(z))
 z
}
