#' @name plot_terms
#' @rdname plot_terms
#'
#' @title Plot terms of non-parametric covariates in PS-SAR, PS-SEM, PS-SARAR, PS-SLX or PS-DURBIN 
#'   regression models.
#'        
#' @description For each non-parametric covariate the plot of the term
#'   includes confidence intervals and the decomposition in fixed and 
#'   random part when the term is reparameterized as a mixed model.  
#'
#' @param fitterms object returned from \code{\link{fit_terms}} function.
#' @param data dataframe or sf with the data. 
#' @param conflevel numerical value for the confidence interval of the
#'   term. Default 0.95.
#' @param listw used to compute spatial lags for Durbin specifications. 
#'   Default =  \code{NULL}
#' @param dynamic Logical value to set a dynamic model.
#'   Dynamic models include a temporal lag of the dependent
#'   variable in the right-hand side of the equation.
#'   Default = \code{FALSE}.
#' @param nt  Number of temporal periods. It is needed
#'   for dynamic models.  
#'
#' @return list with the plots of the terms for each non-parametric 
#'   covariate included in the object returned from \code{\link{fit_terms}}.
#'                                 
#' @author Roman Minguez \email{roman.minguez@@uclm.es}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{fit_terms}} compute smooth functions for non-parametric
#'                                 continuous covariates.
#'   \item \code{\link{impactsnopar}} plot the effects functions 
#'     of non-parametric covariates.
#'   \item \code{\link[mgcv]{vis.gam}} plot the terms fitted by 
#'     \code{\link[mgcv]{gam}} function in \pkg{mgcv} package.   
#' }
#' 
#' @references \itemize{ 
#'   \item Wood, S.N. (2017). \emph{Generalized Additive Models. 
#'   An Introduction with \code{R}} (second edition). CRC Press, Boca Raton.
#'  }
#'         
#' @examples
#'  ###################### unemployment for 103 Italian provinces in period 1996-2014.
#' library(pspatreg)
#' data(unemp_it)
#' 
#' ######################  No Spatial Trend: PSAR including a spatial 
#' ######################  lag of the dependent variable
#' form1 <- unrate ~ partrate + agri + cons +
#'                  pspl(serv,nknots=15) +
#'                  pspl(empgrowth,nknots=20) 
#' gamsar <- pspatfit(form1, data = unemp_it, type = "sar", listw = Wsp_it)
#' summary(gamsar)
#' ######################  Fit non-parametric terms (spatial trend must be name "spttrend")
#' list_varnopar <- c("serv", "empgrowth")
#' terms_nopar <- fit_terms(gamsar, list_varnopar)
#' ######################  Plot non-parametric terms
#' #adjust plot margins
#' par(mar = c(1, 1, 1, 1))
#' plot_terms(terms_nopar, unemp_it)
#'  
#' @export
plot_terms <- function(fitterms, data, conflevel = 0.95, 
                       listw = NULL,  
                       dynamic = FALSE, 
                       nt = NULL) {
  if (inherits(data, "sf"))
    data <- st_drop_geometry(data)
  nfull <- nrow(data)
  if (dynamic) {
    if (is.null(nt)) 
      stop("plot_impactsnopar function needs nt as argument for dynamic models")
    idxyear1 <- seq(from = 1, to = nfull, by = nt)
    data <- data[-idxyear1, ]
  }
  fit <- fitterms$fitted_terms
  se_fit <- fitterms$se_fitted_terms
  fit_fixed <- fitterms$fitted_terms_fixed
  se_fit_fixed <- fitterms$se_fitted_terms_fixed
  fit_random <- fitterms$fitted_terms_random
  se_fit_random <- fitterms$se_fitted_terms_random
  variables <- colnames(fit)
  crval <- qnorm((1 - conflevel)/2, mean = 0, 
                 sd = 1, lower.tail = FALSE)
  for (i in 1:length(variables)) {
    name_var <- variables[i]
    fit_var <- matrix(fit[,c(name_var)], ncol = 1)
    colnames(fit_var) <- name_var
    se_fit_var <- matrix(se_fit[,c(name_var)], ncol = 1)
    colnames(se_fit_var) <- name_var
    up_fit_var <- fit_var + crval*se_fit_var
    colnames(up_fit_var) <- name_var
    low_fit_var <- fit_var - crval*se_fit_var
    colnames(low_fit_var) <- name_var
    fit_var_fixed <- matrix(fit_fixed[,c(name_var)], ncol = 1)
    colnames(fit_var_fixed) <- name_var
    se_fit_var_fixed <- matrix(se_fit_fixed[,c(name_var)], ncol = 1)
    colnames(se_fit_var_fixed) <- name_var
    up_fit_var_fixed <- fit_var_fixed + crval*se_fit_var_fixed
    colnames(up_fit_var_fixed) <- name_var
    low_fit_var_fixed <- fit_var_fixed - crval*se_fit_var_fixed
    colnames(low_fit_var_fixed) <- name_var
    fit_var_random <- matrix(fit_random[, c(name_var)], ncol = 1)
    colnames(fit_var_random) <- name_var
    se_fit_var_random <- matrix(se_fit_random[,c(name_var)], ncol = 1)
    colnames(se_fit_var_random) <- name_var
    up_fit_var_random <- fit_var_random + crval*se_fit_var_random
    colnames(up_fit_var_random) <- name_var
    low_fit_var_random <- fit_var_random - crval*se_fit_var_random
    colnames(low_fit_var_random) <- name_var
    # Check for Wlag.var
    if (grepl("Wlag", name_var)) {
      idx_name_var <- str_detect(name_var, colnames(data))
      var <- as.matrix(data[, idx_name_var])
      Wsp <- listw2mat(listw)
      # spatio-temporal data 
      if (nrow(var) > nrow(Wsp)) {
        nt <- nrow(var) %/% nrow(Wsp)
        It <- Diagonal(nt)
        Wspt <- kronecker(Wsp, It)
        listw <- mat2listw(as.matrix(Wspt))
      }
      colnames(var) <- name_var
      var <- create_WX(var, listw = listw, 
                       zero.policy = TRUE, 
                       prefix = "")
    } else var <- as.matrix(data[, c(name_var)])
    colnames(var) <- name_var
    ord <- order(var)
    par(mfrow = c(2, 1))
    plot(var[ord], fit_var[ord], 
         type = "l",
         ylab = paste("f(", name_var, ")"), 
         xlab = name_var,
         ylim = c(min(low_fit_var), max(up_fit_var)), 
         cex.lab = 1.0, 
         col = 2, 
         lty = 1, 
         lwd = 2, 
         cex.main = 1.0,
         main = paste("Term: ", name_var),
         sub = "Confidence Intervals in dashed lines")
    lines(var[ord], up_fit_var[ord], 
          xlab = "", 
          ylab = "", 
          type = "l", 
          col = 2, 
          lty = 2, 
          lwd = 1.5)
    lines(var[ord], low_fit_var[ord], 
          xlab = "",  
          ylab = "", 
          type = "l", 
          col = 2, 
          lty = 2, 
          lwd = 1.5)
    abline(a = 0, b = 0)
    plot(var[ord], fit_var[ord], 
         type = "l",
         ylab = paste("f(",name_var,")"), 
         xlab = name_var,
         ylim = c(min(low_fit_var), max(up_fit_var)), 
         cex.lab = 1.0, 
         col = 2, 
         lty = 1, 
         lwd = 2,
         cex.main = 1.0,
         main = paste("Global (red), Fixed (green) and Random (blue) terms"))
    lines(var[ord], fit_var_fixed[ord], 
          xlab = "", 
          ylab = "", 
          type = "l", 
          col = 3, 
          lty = 2, 
          lwd = 2)
    lines(var[ord], fit_var_random[ord], 
          xlab = "", 
          ylab = "", 
          type = "l", 
          col = 4, 
          lty = 3, 
          lwd = 2)
    abline(a = 0, b = 0)
    readline(prompt = "Press [enter] to continue")
  }
  par(mfrow = c(1 ,1))
}




