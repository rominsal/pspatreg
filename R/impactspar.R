#' @name impactspar
#' @rdname impactspar
#'
#' @title Compute direct, indirect and total impacts
#'   for continous parametric covariates.
#'        
#' @description Compute direct, indirect and total 
#'   impacts for parametric covariates included in a semiparametric spatial
#'   or spatio-temporal model. The models can be of type \emph{ps-sar}, 
#'   \emph{ps-sarar}, \emph{ps-sdm}, \emph{ps-sdem} or
#'   \emph{ps-slx}.       
#'
#' @param obj A `pspatreg` object created by \code{\link{pspatfit}}.
#' @inheritParams spatialreg::impacts 
#'                      
#' @details 
#'   This function is similar to the 
#'   \code{\link[spatialreg]{impacts}} method used in \pkg{spatialreg}.
#'    package. The function
#'   \code{\link{impactspar}} obtains the three type of impacts (total, direct
#'   and indirect) together with a measure of statistical 
#'   significance, according to the randomization approach described in 
#'   LeSage and Pace (2009). Briefly, they suggest to obtain a 
#'   sequence of \emph{nsim} random matrices using a multivariate normal 
#'   distribution N(0; \emph{Sigma}), being \emph{Sigma} the estimated 
#'   covariance matrix of the fitted \emph{beta} for parametric covariates
#'   and spatial parameters of the model. 
#'   These random matrices, combined with the values of the fitted 
#'   \emph{beta} for parametric covariates and the estimated 
#'   values of the spatial parameters, are used to obtain simulated values. 
#'   The function \code{\link{impactspar}} obtains the standard 
#'   deviations using the \emph{nsim} simulated impacts in the randomization 
#'   procedure, which are used to test the significance of the estimated 
#'   impacts for the original data.
#'   Finally, if the spatial model is type = "slx" or "sdem", then there is no
#'   need to simulate to make inference of the impacts. The standard errors of the 
#'   impacts are computed directly using the \emph{Sigma} matrix of the estimated
#'   covariances of \emph{beta} and spatial parameters. 
#'   
#' @return 
#'   An object of class \emph{impactspar.pspatreg}. Can be printed
#'   with \code{summary}.
#'         
#'  If \code{type = "sar", "sdm", "sarar"}, the object returned is a list 
#'  with 4 objects including the type of model and three matrices including the simulated
#'  total, direct and indirect impacts:         
#'  \tabular{ll}{
#'    \emph{type} \tab Type of spatial econometric model. \cr
#'    \emph{mimpactstot} \tab Matrix including simulated total impacts 
#'                            for each variable in rows. \cr
#'    \emph{mimpactsdir} \tab Matrix including simulated direct impacts 
#'                            for each variable in rows. \cr
#'    \emph{mimpactsind} \tab Matrix including simulated indirect impacts 
#'                            for each variable in rows. \cr
#'    }
#'    If \code{type = "slx", "sdem"} the object returned is a list 
#'    with 5 objects including the type of model and four matrices including 
#'    the computed total, direct and indirect impacts, the standard errors,
#'    the z-values and p-values of each type of impact:
#'  \tabular{ll}{
#'    \emph{type} \tab Type of spatial econometric model. \cr
#'    \emph{mimpact} \tab Matrix including computed total, direct and indirect
#'                        impacts for each variable in rows. \cr
#'    \emph{semimpact} \tab Matrix including standard errors of total, direct 
#'                          and indirect impacts for each variable 
#'                          in rows. \cr
#'    \emph{zvalmimpact} \tab Matrix including z-values of total, direct 
#'                          and indirect impacts for each variable 
#'                          in rows. \cr
#'    \emph{pvalmimpact} \tab Matrix including p-values of total, direct 
#'                          and indirect impacts for each variable 
#'                          in rows. \cr
#'    }
#'     
#' @seealso
#' \itemize{
#'   \item \code{\link{pspatfit}} estimate spatial or spatio-temporal 
#'           semiparametric ps-sar, ps-sem, ps-sarar, ps-slx or 
#'           ps-durbin regression models.
#'   \item \code{\link{impactsnopar}} compute total, direct and indirect impact
#'           functions for non-parametric continuous covariates.
#'   \item \code{\link{fit_terms}} compute smooth term functions for 
#'     non-parametric continuous covariates.
#'   \item \code{\link[spdep]{impacts}} similar function in \pkg{spdep} 
#'           package to compute impacts in spatial parametric econometric 
#'           models.
#' }
#' 
#' @references \itemize{ 
#'   \item LeSage, J. and Pace, K. (2009). \emph{Introduction to 
#'         Spatial Econometrics}. CRC Press, Boca Raton.
#'  }       
#'
#' @examples
#' ################################################
#'  ###################### Examples using a panel data of rate of
#'  ###################### unemployment for 103 Italian provinces in period 1996-2014.
#' library(pspatreg)
#' data(unemp_it, package = "pspatreg")
#' ## Wsp_it is a matrix. Create a neighboord list 
#' lwsp_it <- spdep::mat2listw(Wsp_it)
#' ## short sample for spatial pure case (2d)
#' unemp_it_short <- unemp_it[unemp_it$year == 2019, ]  
#' ########  No Spatial Trend: PSAR including a spatial 
#' ########  lag of the dependent variable
#' form1 <- unrate ~ partrate + agri + cons +
#'                  pspl(serv, nknots = 15) +
#'                  pspl(empgrowth, nknots = 20)
#' ### example with type = "sar"                   
#' gamsar <- pspatfit(form1, 
#'                    data = unemp_it_short, 
#'                    type = "sar", 
#'                    listw = lwsp_it)
#' summary(gamsar)
#' ###### Parametric Total, Direct and Indirect Effects
#' imp_parvar <- impactspar(gamsar, listw = lwsp_it)
#' summary(imp_parvar)
#' 
#' ### example with type = "slx"                   
#' 
#' gamslx <- pspatfit(form1, 
#'                    data = unemp_it_short, 
#'                    type = "slx", 
#'                    listw = lwsp_it)
#'                    
#' summary(gamslx)
#' ###### Parametric Total, Direct and Indirect Effects
#' imp_parvarslx <- impactspar(gamslx, listw = lwsp_it)
#' summary(imp_parvarslx) 
#' \donttest{
#' ######################   PSAR-ANOVA with spatio-temporal trend and 
#' ######################   temporal autorregresive noise
#' form2 <- unrate ~ partrate + agri + cons +
#'                   pspl(serv, nknots = 15) + 
#'                   pspl(empgrowth, nknots = 20) +
#'                   pspt(long, lat, year, 
#'                        nknots = c(18, 18, 8), 
#'                        psanova = TRUE,
#'                        nest_sp1 = c(1, 2, 3),
#'                        nest_sp2 = c(1, 2, 3),
#'                        nest_time = c(1, 2, 2))
#' sptanova_sar_ar1 <- pspatfit(form2, 
#'                              data = unemp_it, 
#'                              listw = lwsp_it, 
#'                              type = "sar", 
#'                              cor = "ar1", 
#'                              control = list(tol = 1e-1))
#' summary(sptanova_sar_ar1)
#' ##### Parametric Total, Direct and Indirect Effects
#' imp_parvar2 <- impactspar(sptanova_sar_ar1, 
#'                           listw = lwsp_it)
#' summary(imp_parvar2)
#' }
#'
#' @export
impactspar <- function(obj, ..., tr = NULL, 
                       R = 1000, listw = NULL, 
                       tol = 1e-06, 
                       Q = NULL) {
  if (!inherits(obj, "pspatreg"))
    stop("object must be of pspatreg class")
  if (obj$type %in% c("sim", "sem")) 
    stop("impact measures not for this model")
  # if (is.null(listw) && !is.null(obj$listw_style) && 
  #     obj$listw_style != "W") 
  #   stop("Only row-standardised weights supported")
  type <- obj$type
  if (is.null(listw)) 
    listw <- obj$listw
  if (inherits(listw, "listw")) {
    Wsp <- Matrix(listw2mat(listw))
  } else if (inherits(listw, "matrix")) {
    Wsp <- Matrix(listw)
  } else if (inherits(listw, "Matrix")) 
    Wsp <- listw   
  nfull <- obj$nfull
  nsp <- obj$nsp
  nt <- obj$nt
  bfixed <- obj$bfixed
  names_var <- names(bfixed)
  names_var <- names_var[!grepl("spt", names_var) & 
                         !grepl("pspl", names_var) &
                         !grepl("Intercept", names_var) &
                         !grepl("fixed_f1_main", names_var) &
                         !grepl("fixed_f2_main", names_var) &
                         !grepl("fixed_ft_main", names_var) &
                         !grepl("fixed_f12", names_var) & 
                         !grepl("fixed_f1t", names_var) & 
                         !grepl("fixed_f2t", names_var) &
                         !grepl("fixed_f12t", names_var) & 
                         !grepl("fixed_tlag", names_var)   ]
  bfixed <- bfixed[names_var]
  names_var <- gsub("fixed_", "", names_var)
  names(bfixed) <- names_var
  vcovb <- vcov(obj)[names_var, names_var]
  se_bfixed <- sqrt(diag(vcovb))
  if (type %in% c("slx", "sdem")) { # No simulation...
    variables <- names(bfixed)[!grepl("Wlag", names_var)]
    mimpacts <- matrix(NA, nrow = length(variables), ncol = 3)
    rownames(mimpacts) <- variables
    colnames(mimpacts) <- c("Direct", "Indirect", "Total")
    se_mimpacts <- mimpacts
    bfixedlag <- bfixed[grepl("Wlag", names(bfixed))]
    se_bfixedlag <- se_bfixed[names(bfixedlag)]
    for (i in 1:length(variables)) {
      vari <- variables[i]
      mimpacts[i, c("Direct")] <- bfixed[vari]
      se_mimpacts[i, c("Direct")] <- se_bfixed[vari]
      idxvarlagi <- grepl(vari, names(bfixedlag))
      if (any(idxvarlagi)) {
        varilag <- names(bfixedlag)[idxvarlagi]
        mimpacts[i, c("Indirect")] <- bfixedlag[varilag]
        se_mimpacts[i, c("Indirect")] <- se_bfixed[varilag]
        mimpacts[i, c("Total")] <- mimpacts[i, c("Direct")] +
                                   mimpacts[i, c("Indirect")]
        varbi <- vcovb[vari, vari]
        varblagi <- vcovb[varilag, varilag]
        covbiblagi <- vcovb[vari, varilag]
        
        se_mimpacts[i, c("Total")] <- sqrt(varbi + varblagi + 
                                             2*covbiblagi)
      } else {
        mimpacts[i, c("Total")] <- mimpacts[i, c("Direct")]
        se_mimpacts[i, c("Total")] <- se_mimpacts[i, c("Direct")]
      }
    }
    zval_mimpacts <- mimpacts / se_mimpacts
    pval_mimpacts <- pnorm(abs(zval_mimpacts), mean = 0, 
                           sd = 1, lower.tail = FALSE)
    res <- list(  type = type,
                  mimpacts = mimpacts,
                  semimpacts = se_mimpacts,
                  zvalmimpacts = zval_mimpacts,
                  pvalmimpacts = pval_mimpacts )
   
  } 
  if (type %in% c("sar", "sdm", "sarar")) {
    names_beta <- names(bfixed)[!grepl("Wlag", names(bfixed))]
    beta <- bfixed[names_beta]
    if (type == "sdm") {
      names_theta <- names(bfixed)[grepl("Wlag", names(bfixed))]
      theta <- bfixed[names_theta]
      # Checking if theta < beta
      for (i in 1:length(beta)) {
        namei <- names_beta[i]
        if (!any(grepl(namei, names_theta))) {
          theta <- c(theta, 0)
          name_thetai <- paste("Wlag.", namei, sep = "")
          names_theta <- c(names_theta, name_thetai)
          names(theta) <- names_theta
          rowvcovb <- rownames(vcovb)
          colvcovb <- colnames(vcovb)
          vcovb <- rbind(vcovb, 0)
          vcovb <- cbind(vcovb, 0)
          rownames(vcovb) <- c(rowvcovb, name_thetai)
          colnames(vcovb) <- c(colvcovb, name_thetai)
          }
      }
    }
    rho <- obj$rho
    se_rho <- obj$se_rho
    if(!is.null(tr)){
      trWsp <- tr
    } else {
      if (is.null(Q)) Q <- 30
      trWsp <- trW(as(Wsp,"CsparseMatrix"), 
                   m = Q, p = 16, type = "MC")
    }
    rho_sim <- rnorm(R, rho, se_rho)
    rho_sim[rho_sim > 1] <- 0.99
    if (type == "sdm") {
      beta_theta_par_sim <- mvrnorm(R, c(beta, theta), 
                              vcovb, tol = tol)
      beta_par_sim <- beta_theta_par_sim[, names_beta]
      theta_par_sim <- beta_theta_par_sim[, names_theta]
      rm(beta_theta_par_sim)
    }
    else beta_par_sim <- mvrnorm(R, beta, 
                                 vcovb, tol = tol)
    beta_par_sim <- t(beta_par_sim)
    if (type == "sdm")
      theta_par_sim <- t(theta_par_sim)
    q <- length(trWsp)
    nsp <- nrow(Wsp)
    if (type == "sdm") {
      q <- q - 1
      mT <- matrix(c(1, trWsp[1:q]/nsp, trWsp[1:(q+1)]/nsp),
                   nrow = 2, byrow = TRUE)
    }
    else mT <- matrix(c(1, trWsp[1:q]/nsp), nrow = 1)
    a <- matrix(1, nrow = q + 1, ncol = 1)
    #   G <- diag(g)
    mP <- list()
    k <- length(beta)
    mimpactsdir <- matrix(NA, nrow = k, ncol = R)
    mimpactstot <- mimpactsdir
    for (i in 1:R) {
      if (type == "sdm")
        mP[[i]] <- matrix(c(beta_par_sim[, i], 
                            theta_par_sim[, i]),
                            ncol = 2, byrow = FALSE)
      else
        mP[[i]] <- matrix(beta_par_sim[, i], ncol = 1)
      g <- 1
      for (j in 1:q)  
        g <- c(g, rho_sim[i]^j) 
      G <- diag(g)
      mimpactsdir[, i] <- mP[[i]] %*% mT %*% G %*% a
      mimpactstot[, i] <- matrix(rowSums(mP[[i]]), ncol = 1) %*%
        (matrix(g,nrow=1) %*% a)
    }
    mimpactsind <- mimpactstot - mimpactsdir
    rownames(mimpactsind) <- names(beta) 
    rownames(mimpactstot) <- names(beta) 
    rownames(mimpactsdir) <- names(beta)
    res <- list( type = type,
                 mimpactstot = mimpactstot, 
                 mimpactsdir = mimpactsdir, 
                 mimpactsind = mimpactsind )
  }
  class(res) <- ("impactspar.pspatreg")
  res
}

