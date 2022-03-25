#' @name fit_terms
#' @rdname fit_terms
#'
#' @title Compute terms for non-parametric spatial (2d) or 
#'   spatio-temporal (3d) trends and for smooth functions of continous 
#'   non-parametric covariates in PS-SAR, PS-SEM, PS-SARAR, PS-SLX or PS-DURBIN regression models.
#'        
#' @description The \code{fit_terms} function compute both:
#' \itemize{
#'   \item Non-parametric spatial (2d) or spatio-temporal (3d) trends 
#'     including the decomposition in main and interaction trends 
#'     when the model is ANOVA.
#'     \item Smooth functions \eqn{f(x_i)} for non-parametric covariates 
#'       in semiparametric models. It also includes standard errors and the 
#'       decomposition of each non-parametric
#'       term in fixed and random parts.
#' }
#'        
#' @param object object fitted using \code{\link{pspatfit}} function. 
#' @param variables vector including names of non-parametric covariates. 
#'   To fit the terms of non-parametric spatial (2d) or spatio-temporal 
#'   (3d) trend this argument must be set equal to \emph{spttrend}. See \code{examples}.
#'   
#' @return A list including:
#'   \tabular{ll}{
#'     \emph{fitted_terms} \tab Matrix including terms in columns. \cr
#'     \emph{se_fitted_terms} \tab Matrix including standard errors of terms 
#'       in columns. \cr
#'     \emph{fitted_terms_fixed} \tab Matrix including fixed part of terms 
#'       in columns. \cr
#'     \emph{se_fitted_terms_fixed} \tab Matrix including standard errors of 
#'       fixed part of terms in columns. \cr
#'     \emph{fitted_terms_random} \tab Matrix including random part of terms 
#'       in columns. \cr
#'     \emph{se_fitted_terms_random} \tab Matrix including standard errors of 
#'       random part of terms in columns.\cr
#'  }
#'  This object can be used as an argument of \code{\link{plot_terms}} function
#'  to make plots of both non-parametric trends and smooth functions of 
#'  covariates. See \emph{examples} below. 
#'  
#' @author Roman Minguez \email{roman.minguez@@uclm.es} 
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{pspatfit}} estimate spatial or spatio-temporal 
#'     semiparametric PS-SAR, PS-SEM, PS-SARAR, PS-SLX or PS-DURBIN regression models.
#'   \item \code{\link{plot_terms}} plot smooth functions of non-parametric
#'     covariates.
#' }
#' 
#' @references \itemize{ 
#'  \item Lee, D. and Durbán, M. (2011). P-Spline ANOVA Type Interaction 
#'         Models for Spatio-Temporal Smoothing. \emph{Statistical Modelling}, 
#'         (11), 49-69.
#'         
#'   \item Wood, S.N. (2017). \emph{Generalized Additive Models. 
#'   An Introduction with \code{R}} (second edition). CRC Press, Boca Raton.
#'  } 
#'        
#' @examples
#' ###################### Examples using a panel data of rate of unemployment 
#' ###################### in 103 Italian provinces during the period 1996-2014.
#' library(pspatreg)
#' data(unemp_it)
#' unemp_it_short <- unemp_it[unemp_it$year == 2019, ]
#' ######################  No Spatial Trend: PSAR including a spatial 
#' ######################  lag of the dependent variable
#' form1 <- unrate ~ partrate + agri + cons +
#'                  pspl(serv, nknots = 15) +
#'                  pspl(empgrowth, nknots = 20) 
#' gamsar <- pspatfit(form1, data = unemp_it_short, 
#'                    type = "sar", listw = Wsp_it)
#' summary(gamsar)
#' 
#' ######################  Fit non-parametric terms 
#' ######################  (spatial trend must be name "spttrend")
#' list_varnopar <- c("serv", "empgrowth")
#' terms_nopar <- fit_terms(gamsar, list_varnopar)
#' 
#' ######################  Plot non-parametric terms
#' #adjust plot margins
#' par(mar = c(1, 1, 1, 1))
#' plot_terms(terms_nopar, unemp_it)
#'   
#' ###############################################
#' ############### Examples of terms corresponding to spatial (2d) or 
#' ###############       spatio-temporal (3d) trends
#' #################################################  
#' # Spatial (2d) semiparametric ANOVA model without spatial lag
#' # Interaction term f12 with nested basis
#' form3 <- unrate ~ partrate + agri + cons +
#'                   pspl(serv,nknots = 15) + 
#'                   pspl(empgrowth, nknots = 20) +
#'                   pspt(long, lat, nknots = c(20, 20), psanova = TRUE,
#'                   nest_sp1 = c(1, 2), 
#'                   nest_sp2 = c(1, 2))
#' # Spatial trend fixed for period 1996-2019
#' geospanova <- pspatfit(form3, data = unemp_it_short)
#' summary(geospanova)
#' ### Plot spatial trend and interaction (ANOVA)
#' ### Create sf object to make the plot 
#' library(sf)
#' unemp_it_sf_short <- st_as_sf(dplyr::left_join(
#'                               unemp_it_short, 
#'                               map_it,  
#'                         by = c("prov" = "COD_PRO")))
#' plot_sp2d(geospanova, 
#'           data = unemp_it_sf_short, 
#'           addmain = TRUE, 
#'           addint = TRUE)
#' ###############################################
#' # Spatio-temporal (3d) semiparametric ANOVA model without spatial lag
#' # Interaction terms f12,f1t,f2t and f12t with nested basis
#' # Remark: It is necessary to include ntime as argument
#' # Remark: nest_sp1, nest_sp2 and nest_time must be divisors of nknots
#' form4 <- unrate ~ partrate + agri + cons +
#'                   pspl(serv, nknots = 15) + 
#'                   pspl(empgrowth, nknots = 20) +
#'                   pspt(long, lat, year, nknots = c(18, 18, 8), 
#'                        psanova = TRUE,
#'                   nest_sp1 = c(1, 2, 3), 
#'                   nest_sp2 = c(1, 2, 3),
#'                   nest_time = c(1, 2, 2))
#' sptanova <- pspatfit(form4, data = unemp_it, 
#'                      control = list(tol = 1e-1))
#' summary(sptanova)
#' 
#' ####### Plot spatio-temporal trends
#' library(sf)
#' unemp_it_sf <- st_as_sf(dplyr::left_join(
#'                                 unemp_it, map_it,  
#'                         by = c("prov" = "COD_PRO")))
#' 
#' plot_sp3d(sptanova, data = unemp_it_sf, 
#'           time_var = "year", 
#'           time_index = c(1996, 2019),
#'           addmain = TRUE, addint = FALSE)
#'  
#' @export
fit_terms <- function(object, variables){
  X <- object$X
  Z <- object$Z
  bfixed <- object$bfixed
  brandom <- object$brandom
  psanova <- object$psanova
  cov_b <-object$vcov_b
  row_cov_fixed <- c(grepl("fixed", rownames(cov_b)))
  col_cov_fixed <- c(grepl("fixed", colnames(cov_b)))
  cov_b_fixed <- cov_b[row_cov_fixed,col_cov_fixed]
  row_cov_random <- c(grepl("random", rownames(cov_b)))
  col_cov_random <- c(grepl("random", colnames(cov_b)))
  cov_b_random <- cov_b[row_cov_random,col_cov_random]
  fitted_terms_fixed <- fitted_terms_random <- fitted_terms <- NULL
  se_fitted_terms <- se_fitted_terms_fixed <- se_fitted_terms_random <- NULL
  for (i in 1:length(variables)) {
    var_name <- variables[i]
    if (grepl("spttrend",var_name)) {
       time <- object$time
       if (psanova) {# psanova=TRUE
         if (!is.null(time)) { # 3d
           eff_spttrend_psanova <- c("Intercept","f1_main","f2_main","ft_main",
                                     "f12_int", "f1t_int", "f2t_int", 
                                     "f12t_int")
         } else { # 2d
           eff_spttrend_psanova <- c("Intercept","f1_main","f2_main","f12_int")
         }
         for (j in 1:length(eff_spttrend_psanova)) {
           eff_spttrend_psanova_j <- eff_spttrend_psanova[j]
           Xj <- as.matrix(X[,grepl(eff_spttrend_psanova_j, 
                                    colnames(X))])
           bfixed_j <- bfixed[grepl(eff_spttrend_psanova_j, 
                                    names(bfixed))]
           Zj <- as.matrix(Z[,grepl(eff_spttrend_psanova_j, 
                                    colnames(Z))])
           brandom_j <- brandom[grepl(eff_spttrend_psanova_j, 
                                      names(brandom))]
           term_fixed_j <- Xj %*% bfixed_j
           term_random_j <- Zj %*% brandom_j
           term_j <- term_fixed_j + term_random_j
           names(term_fixed_j) <- names(term_random_j) <- eff_spttrend_psanova_j
           names(term_j) <- eff_spttrend_psanova_j
           fitted_terms_fixed <- cbind(term_fixed_j, 
                                       fitted_terms_fixed)
           fitted_terms_random <- cbind(term_random_j, 
                                        fitted_terms_random)
           fitted_terms <- cbind(term_j, fitted_terms)
           colnames(fitted_terms_fixed)[1] <- eff_spttrend_psanova_j
           colnames(fitted_terms_random)[1] <- eff_spttrend_psanova_j
           colnames(fitted_terms)[1] <- eff_spttrend_psanova_j
           row_cov_j <- c(grepl(eff_spttrend_psanova_j, 
                                rownames(cov_b)))
           col_cov_j <- c(grepl(eff_spttrend_psanova_j, 
                                colnames(cov_b)))
           cov_b_j <- cov_b[row_cov_j,col_cov_j]
           se_term_j <- rowSums( (cbind(Xj,Zj) %*% cov_b_j)
                                        * cbind(Xj,Zj) )^0.5
           se_fitted_terms <- cbind(se_term_j,se_fitted_terms)
           colnames(se_fitted_terms)[1] <- eff_spttrend_psanova_j
           row_cov_j_fixed <- c(grepl(eff_spttrend_psanova_j,
                                      rownames(cov_b_fixed)))
           col_cov_j_fixed <- c(grepl(eff_spttrend_psanova_j,
                                      colnames(cov_b_fixed)))
           cov_b_j_fixed <- cov_b_fixed[row_cov_j_fixed,col_cov_j_fixed]
           se_term_j_fixed <- rowSums((Xj %*% cov_b_j_fixed)
                                               * Xj)^0.5
           se_fitted_terms_fixed <- cbind(se_term_j_fixed,
                                          se_fitted_terms_fixed)
           colnames(se_fitted_terms_fixed)[1] <- eff_spttrend_psanova_j
           row_cov_j_random <- c(grepl(eff_spttrend_psanova_j,
                                       rownames(cov_b_random)))
           col_cov_j_random <- c(grepl(eff_spttrend_psanova_j,
                                       colnames(cov_b_random)))
           cov_b_j_random <- cov_b_random[row_cov_j_random, 
                                          col_cov_j_random]
           se_term_j_random <- rowSums((Zj %*% cov_b_j_random)
                                               * Zj)^0.5
           se_fitted_terms_random <- cbind(se_term_j_random,
                                          se_fitted_terms_random)
           colnames(se_fitted_terms_random)[1] <- eff_spttrend_psanova_j
         }
         match_fixed <- unique(grep(paste(eff_spttrend_psanova,
                                               collapse="|"),
                                         names(bfixed), value=TRUE))
         bfixed_spt <- bfixed[match_fixed]
         match_random <- unique(grep(paste(eff_spttrend_psanova,
                                          collapse="|"),
                                    names(brandom), value=TRUE))
         brandom_spt <- brandom[match_random]
         match_X <- unique(grep(paste(eff_spttrend_psanova,
                                          collapse="|"),
                                    colnames(X), value=TRUE))
         X_spt <- X[,match_X]
         match_Z <- unique(grep(paste(eff_spttrend_psanova,
                                      collapse="|"),
                                colnames(Z), value=TRUE))
         Z_spt <- Z[,match_Z]
         fitted_term_fixed_spt <- X_spt %*% bfixed_spt
         fitted_term_random_spt <- Z_spt %*% brandom_spt
         fitted_term_spt <- fitted_term_fixed_spt + fitted_term_random_spt

         fitted_terms_fixed <- cbind(fitted_term_fixed_spt,
                                     fitted_terms_fixed)
         fitted_terms_random <- cbind(fitted_term_random_spt,
                                      fitted_terms_random)
         fitted_terms <- cbind(fitted_term_spt, fitted_terms)
         colnames(fitted_terms_fixed)[1] <- "spttrend"
         colnames(fitted_terms_random)[1] <- "spttrend"
         colnames(fitted_terms)[1] <- "spttrend"
         match_cov_bspt <- unique(grep(paste(eff_spttrend_psanova,
                                     collapse="|"),
                               colnames(cov_b), value=TRUE))
         cov_b_spt <- cov_b[match_cov_bspt, match_cov_bspt]
         se_term_spt <- rowSums((cbind(X_spt,Z_spt) %*% cov_b_spt)
                                        * cbind(X_spt,Z_spt))^0.5
         se_fitted_terms <- cbind(se_term_spt,
                                  se_fitted_terms)
         colnames(se_fitted_terms)[1] <- "spttrend"
         match_cov_bspt_fixed <- unique(grep(paste(eff_spttrend_psanova,
                                             collapse="|"),
                                       colnames(cov_b_fixed), 
                                       value=TRUE))
         cov_b_spt_fixed <- cov_b[match_cov_bspt_fixed, 
                                  match_cov_bspt_fixed]
         se_term_spt_fixed <- rowSums((X_spt %*% cov_b_spt_fixed)
                                              * X_spt)^0.5
         se_fitted_terms_fixed <- cbind(se_term_spt_fixed,
                                        se_fitted_terms_fixed)
         colnames(se_fitted_terms_fixed)[1] <- "spttrend"
         match_cov_bspt_random <- unique(grep(paste(eff_spttrend_psanova,
                                                   collapse="|"),
                                             colnames(cov_b_random), 
                                             value=TRUE))
         cov_b_spt_random <- cov_b[match_cov_bspt_random, 
                                   match_cov_bspt_random]
         se_term_spt_random <- rowSums( (Z_spt %*% cov_b_spt_random)
                                               * Z_spt )^0.5
         se_fitted_terms_random <- cbind(se_term_spt_random,
                                         se_fitted_terms_random)
         colnames(se_fitted_terms_random)[1] <- "spttrend"
      } else { # psanova=FALSE
        Xi <- as.matrix(X[,grepl("spt",colnames(X))])
        bfixed_i <- bfixed[grepl("spt",names(bfixed))]
        Zi <- as.matrix(Z[,grepl("spt",colnames(Z))])
        brandom_i <- brandom[grepl("spt",names(brandom))]
        term_fixed_i <- Xi %*% bfixed_i
        term_random_i <- Zi %*% brandom_i
        term_i <- term_fixed_i + term_random_i
        fitted_terms_fixed <- cbind(term_fixed_i,fitted_terms_fixed)
        fitted_terms_random <- cbind(term_random_i, fitted_terms_random)
        fitted_terms <- cbind(term_i,fitted_terms)
        colnames(fitted_terms_fixed)[1] <- "spttrend"
        colnames(fitted_terms_random)[1] <- "spttrend"
        colnames(fitted_terms)[1] <- "spttrend"
        row_cov_i <- c(grepl("spt",rownames(cov_b)))
        col_cov_i <- c(grepl("spt",colnames(cov_b)))
        cov_b_i <- cov_b[row_cov_i,col_cov_i]
        se_term_i <- rowSums((cbind(Xi,Zi) %*% cov_b_i)
                                     * cbind(Xi,Zi))^0.5
        se_fitted_terms <- cbind(se_term_i,
                                 se_fitted_terms)
        colnames(se_fitted_terms)[1] <- "spttrend"
        row_cov_i_fixed <- c(grepl("spt", 
                                   rownames(cov_b_fixed)))
        col_cov_i_fixed <- c(grepl("spt", 
                                   colnames(cov_b_fixed)))
        cov_b_i_fixed <- cov_b_fixed[row_cov_i_fixed, 
                                     col_cov_i_fixed]
        se_term_i_fixed <- rowSums((Xi %*% cov_b_i_fixed)
                                           * Xi)^0.5
        se_fitted_terms_fixed <- cbind(se_term_i_fixed,
                                       se_fitted_terms_fixed)
        colnames(se_fitted_terms_fixed)[1] <- "spttrend"
        row_cov_i_random <- c(grepl("spt", 
                                    rownames(cov_b_random)))
        col_cov_i_random <- c(grepl("spt", 
                                    colnames(cov_b_random)))
        cov_b_i_random <- cov_b_random[row_cov_i_random, 
                                       col_cov_i_random]
        se_term_i_random <- rowSums((Zi %*% cov_b_i_random)
                                            * Zi)^0.5
        se_fitted_terms_random <- cbind(se_term_i_random,
                                       se_fitted_terms_random)
        colnames(se_fitted_terms_random)[1] <- "spttrend"
      }
    } else { # No spttrend
      Xi <- X[, grepl(var_name, colnames(X)), drop = FALSE]
      bfixed_i <- bfixed[grepl(var_name, 
                               names(bfixed))]
      Zi <- Z[,grepl(var_name, colnames(Z)), drop = FALSE]
      brandom_i <- brandom[grepl(var_name, 
                                 names(brandom))]
      # Divide var and Wlag_var in Durbin case
      k <- 1
      if (any(grepl("Wlag", colnames(Xi)))) {
        k <- 2
        Wlag_Xi <- Xi[, grepl("Wlag", colnames(Xi)), drop = FALSE]
        Xi <- Xi[, !grepl("Wlag", colnames(Xi)), drop = FALSE]
        Wlag_Zi <- Zi[, grepl("Wlag", colnames(Zi)), drop = FALSE]
        Zi <- Zi[, !grepl("Wlag", colnames(Zi)), drop = FALSE]
        Wlag_bfixed_i <-  bfixed_i[grepl("Wlag", 
                                       names(bfixed_i))]
        bfixed_i <- bfixed_i[!grepl("Wlag", 
                                 names(bfixed_i))]
        Wlag_brandom_i <-  brandom_i[grepl("Wlag", 
                                       names(brandom_i))]
        brandom_i <-  brandom_i[!grepl("Wlag", 
                                        names(brandom_i))]
      }
      for (j in 1:k) { # If necessary repeat the process for Wlagvar (k = 2)
        if (j == 2) {
          var_name <- paste0("Wlag.", var_name)
          Xi <- Wlag_Xi
          Zi <- Wlag_Zi
          bfixed_i <- Wlag_bfixed_i
          brandom_i <- Wlag_brandom_i
        }
        term_fixed_i <- as.matrix(Xi %*% bfixed_i)
        term_random_i <- as.matrix(Zi %*% brandom_i)
        term_i <- as.matrix(term_fixed_i + term_random_i)
        colnames(term_i) <- var_name
        colnames(term_fixed_i) <- var_name
        colnames(term_random_i) <- var_name
        fitted_terms_fixed <- cbind(fitted_terms_fixed, 
                                    term_fixed_i)
        fitted_terms_random <- cbind(fitted_terms_random, 
                                     term_random_i)
        fitted_terms <- cbind(fitted_terms, term_i)
        cov_b_i <- cov_b[c(names(bfixed_i), names(brandom_i)),
                         c(names(bfixed_i), names(brandom_i))]
        se_term_i <- as.matrix( rowSums((cbind(Xi, Zi) %*% cov_b_i)
                             * cbind(Xi, Zi))^0.5 )
        colnames(se_term_i) <- var_name
        se_fitted_terms <- cbind(se_fitted_terms, se_term_i)
        cov_b_i_fixed <- cov_b_fixed[names(bfixed_i), 
                                     names(bfixed_i)]
        se_term_i_fixed <- as.matrix( rowSums( (Xi %*% cov_b_i_fixed)
                                    * Xi )^0.5 )
        colnames(se_term_i_fixed) <- var_name
        se_fitted_terms_fixed <- cbind(se_fitted_terms_fixed, 
                                       se_term_i_fixed)
        cov_b_i_random <- cov_b_random[names(brandom_i), 
                                       names(brandom_i)]
        se_term_i_random <- as.matrix( rowSums((Zi %*% cov_b_i_random)
                                    * Zi)^0.5 )
        colnames(se_term_i_random) <- var_name
        se_fitted_terms_random <- cbind(se_fitted_terms_random, 
                                        se_term_i_random)
      }
    }
  } # end for (i in 1:length(variables))
  res <- list(fitted_terms = as.matrix(fitted_terms),
              se_fitted_terms = as.matrix(se_fitted_terms),
              fitted_terms_fixed = as.matrix(fitted_terms_fixed),
              se_fitted_terms_fixed = as.matrix(se_fitted_terms_fixed),
              fitted_terms_random = as.matrix(fitted_terms_random),
              se_fitted_terms_random = as.matrix(se_fitted_terms_random) )
  res
}
