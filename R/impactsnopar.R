#' @name impactsnopar
#' @rdname impactsnopar
#'
#' @title Compute direct, indirect and total impacts functions 
#'   for continous non-parametric covariates in semiparametric spatial 
#'   regression models.
#'        
#' @description Compute and plot direct, indirect and total impect functions 
#'   for non-parametric covariates included in a semiparametric spatial
#'   or spatio-temporal econometric model. This model must include a 
#'   spatial lag of the dependent variable and/or non-parametric covariates, 
#'   to have indirect impacts different from 0, otherwise, total and direct 
#'   function impacts are the same.         
#'
#' @param obj \emph{pspatfit} object fitted using \code{\link{pspatfit}} function. 
#' @param listw should be a spatial neighbours list object created for example by \code{nb2listw} from \code{spdep} package. 
#' It can also be a spatial weighting matrix of order (NxN) instead of a listw neighbours list object.
#' @param conflevel numerical value for the confidence interval of the
#'    impact functions. Default 0.95.
#' @param viewplot Default `TRUE` to plot impacts. If FALSE use \code{\link{plot_impactsnopar}} to plot impacts
#' @param smooth Default `TRUE`. Whether to smooth fitted impacts or not.
#' @param span span for the kernel of the smoothing (see \code{\link{loess}} 
#'             for details). Default c(0.1, 0.1, 0.2) 
#'             
#' @details
#'   To compute the impact functions of the non-parametric covariates, first 
#'   it is used the function 
#'   \code{\link{fit_terms}} to get fitted values of the terms and 
#'   standard errors of the fitted values for each non-parametric covariate. 
#'   Then, the intervals for the fitted term are computed as \cr
#'   
#'   \emph{ fitted_values plus/minus quantile*standard errors } \cr
#'   
#'   where \emph{quantile} is the corresponding quantile of the N(0,1) 
#'   distribution. The total impact function is computed as: \cr
#'   
#'   \code{solve(kronecker((I_N - rho*W_N), It), fitted_values) } \cr
#'   
#'   where \emph{(I_N - rho*W_N)} matrix is the spatial lag matrix and 
#'   \emph{It} is an identity matrix of order equals to the temporal periods 
#'   (\emph{t}). Obviously, \emph{t = 1} for pure spatial econometric models. 
#'   The upper and lower bounds of the total impact functions are computed 
#'   using the previous formula but using 
#'   \emph{fitted_values plus/minus quantile*standard errors} instead
#'   of \emph{fitted_values}. \cr
#'   
#'   The direct impacts function is computed using the formula: \cr
#'   
#'   \code{diag(solve(kronecker((I_N - rho*W_N), It),
#'                               diag(fitted_values))} \cr
#'                               
#'   that is, the fitted values are put in the main diagonal of 
#'   a diagonal matrix and, afterwards, the spatial lag is applied over 
#'   this diagonal matrix. Finally, the main diagonal of the resulting
#'   matrix is considered the direct impact function. 
#'   The upper and lower bounds of the direct impact functions are computed 
#'   using the previous formula but using 
#'   \emph{fitted_values plus/minus quantile*standard errors} instead
#'   of \emph{fitted_values}. \cr
#'   
#'   Eventually, the indirect impacts function are computed as the 
#'   difference between both total and direct impact functions, that is: \cr
#'   
#'    \emph{indirect impact function = total impacts function -
#'                                     direct impacts function} \cr
#'                                     
#'  In this way we can get both, the indirect impact functions and upper and
#'  lower bounds of the indirect impact functions. \cr
#'  
#'  It is important to remark that, usually, the indirect impact functions 
#'  are very wiggly. To get ride of this problem, the argument \code{smooth} 
#'  (default = `TRUE`) allows to smooth the impacts function using the 
#'  \code{\link[stats]{loess}} 
#'  function available in \pkg{stats}. This is very convenient when the 
#'  indirect impacts function is plotted.                                              
#' 
#'
#' @return A list including
#'   \tabular{ll}{
#'     \emph{impnopar_tot} \tab Matrix including total impacts in columns. \cr
#'     \emph{impnopar_dir} \tab Matrix including direct impacts in columns. \cr
#'     \emph{impnopar_ind} \tab Matrix including indirect 
#'                              impacts in columns. \cr
#'     \emph{impnopar_tot_up} \tab Matrix including upper bounds of 
#'                                 total impacts in columns. \cr
#'     \emph{impnopar_dir_up} \tab Matrix including upper bounds of 
#'                                 direct impacts in columns. \cr
#'     \emph{impnopar_ind_up} \tab Matrix including upper bounds of 
#'                                 indirect impacts in columns. \cr
#'     \emph{impnopar_tot_low} \tab Matrix including lower bounds of 
#'                                  total impacts in columns. \cr
#'     \emph{impnopar_dir_low} \tab Matrix including lower bounds of 
#'                                  direct impacts in columns. \cr
#'     \emph{impnopar_ind_low} \tab Matrix including lower bounds of 
#'                                  indirect impacts in columns. \cr
#'  }
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
#'   \item \code{\link{pspatfit}} estimate spatial or spatio-temporal semiparametric PS-SAR
#'     regression models.
#'   \item \code{\link{impactspar}} compute and simulate total, direct and indirect impect
#'     (or impacts) for parametric continuous covariates.
#'   \item \code{\link{fit_terms}} compute terms for smooth functions for non-parametric
#'     continuous covariates and for non-parametric trends.
#'   \item \code{\link{plot_impactsnopar}} plot the non-parametric impacts functions
#'     allowing for previous smoothing.
#' }
#' 
#' @references 
#' \itemize{
#'     \item Basile, R.; Durban, M.; Minguez, R.; Montero, J. M.; and 
#'     Mur, J. (2014). Modeling regional economic dynamics: Spatial
#'     dependence, spatial heterogeneity and nonlinearities. 
#'     \emph{Journal of Economic Dynamics and Control}, (48), 229-245.
#'     <doi: 10.1016/j.jedc.2014.06.011>
#'         
#'  \item Eilers, P. and Marx, B. (2021). \emph{Practical Smoothing. 
#'        The Joys of P-Splines}. Cambridge University Press.
#'     
#'  \item Fahrmeir, L.; Kneib, T.;  Lang, S.; and Marx, B. (2013). 
#'        \emph{Regression. Models, Methods and Applications}.
#'        Springer.         
#'         
#'   \item LeSage, J. and Pace, K. (2009). \emph{Introduction to 
#'         Spatial Econometrics}. CRC Press, Boca Raton.
#'         
#'   \item Minguez, R.; Basile, R. and Durban, M. (2020). An Alternative 
#'     Semiparametric Model for Spatial Panel Data. \emph{Statistical Methods and Applications},
#'     (29), 669-708. <doi:	10.1007/s10260-019-00492-8>
#'
#'   \item Montero, J., Minguez, R., and Durban, M. (2012). SAR models 
#'     with nonparametric spatial trends: A P-Spline approach. 
#'     \emph{Estadistica Espanola}, (54:177), 89-111.
#
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
#' ######  No Spatial Trend: PSAR including a spatial 
#' ######  lag of the dependent variable
#' form1 <- unrate ~ partrate + agri + cons +
#'                  pspl(serv, nknots = 15) +
#'                  pspl(empgrowth, nknots = 20) 
#' gamsar <- pspatfit(form1, 
#'                     data = unemp_it_short, 
#'                     type = "sar", 
#'                     listw = lwsp_it)
#'  summary(gamsar)
#'  ###### Non-Parametric Total, Direct and Indirect impacts
#'  ## adjust plot margins
#'  par(mar = c(1, 1, 1, 1))
#'  imp_nparvar <- impactsnopar(gamsar, 
#'                              listw = lwsp_it, 
#'                              viewplot = TRUE)
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
#' ###### Non-Parametric Total, Direct and Indirect impacts
#' ## adjust plot margins
#' par(mar = c(1, 1, 1, 1))
#' imp_nparvar2 <- impactsnopar(sptanova_sar_ar1, 
#'                              listw = lwsp_it, 
#'                              viewplot = TRUE)
#'  }                            
#'
#' @export
impactsnopar <- function(obj, listw = NULL, conflevel = 0.95, 
                         viewplot = TRUE, smooth = TRUE, 
                         span = c(0.1, 0.1, 0.2)) {
  type <- obj$type
  if (is.null(listw)) 
    listw <- obj$listw
  if (inherits(listw, "listw")) {
    Wsp <- Matrix(listw2mat(listw))
  } else if (inherits(listw, "matrix")) {
    Wsp <- Matrix(listw)
  } else if (inherits(listw, "Matrix")) 
    Wsp <- listw
  nsp <- obj$nsp
  nt <- obj$nt
  bfixed <- obj$bfixed
  dynamic <- obj$dynamic
  if (dynamic) nt <- nt - 1 #Remove first observ.
  varnopar <- names(bfixed)[grepl("pspl", names(bfixed))]
  # Prevent repetition of names for Wlag variables...
  varnopar <- varnopar[!grepl("Wlag", varnopar)]
  if (!(length(varnopar) > 0))
    stop("there is no any nonparametric variable in this model")
  # Build list of nonparametric variables
  variables <- varnopar
  for (i in 1:length(varnopar)) {
    idx_varnopari <- str_detect(varnopar[i], colnames(obj$data))
    variables[i] <- colnames(obj$data)[idx_varnopari]
  }
  fitsall <- fit_terms(obj, variables)
  fits <- fitsall$fitted_terms
  sefits <- fitsall$se_fitted_terms
  crval <- qnorm((1 - conflevel)/2, mean = 0, 
                 sd = 1, lower.tail = FALSE)
  fitsup <- fits + crval*sefits
  fitslow <- fits - crval*sefits
  if (type %in% c("slx", "sdm", "sdem")) {
    idxvarWlag <- str_detect(colnames(fits), "Wlag")
    namesWlag <- colnames(fits)[idxvarWlag]
    namesvar <- colnames(fits)[!idxvarWlag]
    
    impactsind <- fits[, namesWlag, drop = FALSE]
    impactsindup <- fitsup[, namesWlag, drop = FALSE]
    impactsindlow <- fitslow[, namesWlag, drop = FALSE]
    
    impactsdir <- fits[, namesvar, drop = FALSE]
    impactsdirup <- fitsup[, namesvar, drop = FALSE]
    impactsdirlow <- fitslow[, namesvar, drop = FALSE]
    
    # Add zeros to impactsind when there is no Wlag variable
    for (i in 1:length(namesvar)) {
      idxWlagi <- grepl(namesvar[i], namesWlag)
      if (!any(idxWlagi)) {
         varWlagi <- paste("Wlag.", namesvar[i], sep = "")
         namesimpactsind <- c(colnames(impactsind),
                              varWlagi)
         impactsind <- cbind(impactsind, 0)
         impactsindup <- cbind(impactsindup, 0)
         impactsindlow <- cbind(impactsindlow, 0)
         colnames(impactsind) <- namesimpactsind
         colnames(impactsindup) <- namesimpactsind
         colnames(impactsindlow) <- namesimpactsind
      }
    }
    colnames(impactsind) <- gsub("Wlag.", "", 
                                 colnames(impactsind))
    colnames(impactsindup) <- gsub("Wlag.", "", 
                                   colnames(impactsindup))
    colnames(impactsindlow) <- gsub("Wlag.", "", 
                                   colnames(impactsindlow))
    impactstot <- impactsdir[, namesvar] +
                  impactsind[, namesvar]
    impactstotup <- impactsdirup[, namesvar] +
                    impactsindup[, namesvar]
    impactstotlow <- impactsdirlow[, namesvar] +
                     impactsindlow[, namesvar]
  } else {
    impactstot <- fits
    impactstotup <- fitsup
    impactstotlow <- fitslow
    impactsdir <- fits
    impactsdirup <- fitsup
    impactsdirlow <- fitslow
  }
  if (type %in% c("sar", "sdm", "sarar"))
    rho <- obj$rho
  else rho <- 0
  A <- Diagonal(nsp) - rho*Wsp
  It <- Diagonal(nt)
  impactstot <- solve(kronecker(A, It), impactstot)
  impactstotup <- solve(kronecker(A, It), impactstotup)
  impactstotlow <- solve(kronecker(A, It), impactstotlow)
  for (i in 1:ncol(impactsdir)) {
    fitsi <- impactsdir[, i]
    impactsdir[, i] <- diag(solve(kronecker(A, It),
                                  diag(fitsi)))
    fitsupi <- impactsdirup[, i]
    impactsdirup[, i] <- diag(solve(kronecker(A, It),
                                    diag(fitsupi)))
    fitslowi <- impactsdirlow[, i]
    impactsdirlow[, i] <- diag(solve(kronecker(A, It),
                                     diag(fitslowi)))
  }
  impactsind <- impactstot - impactsdir
  impactsindup <- impactstotup - impactsdirup
  impactsindlow <- impactstotlow - impactsdirlow
  res <- list( impnopar_tot = impactstot,
               impnopar_dir = impactsdir,
               impnopar_ind = impactsind,
               impnopar_tot_up = impactstotup,
               impnopar_dir_up = impactsdirup,
               impnopar_ind_up = impactsindup,
               impnopar_tot_low = impactstotlow,
               impnopar_dir_low = impactsdirlow,
               impnopar_ind_low = impactsindlow)
  if (viewplot) {
   plot_impactsnopar(res, data = obj$data, 
                     smooth = smooth,
                     span = span,
                     dynamic = obj$dynamic,
                     nt = obj$nt)
  }
 invisible(res)
}



