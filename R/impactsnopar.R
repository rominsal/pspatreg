#' @name impactsnopar
#' @rdname impactsnopar
#'
#' @title Compute direct, indirect and total impacts functions 
#'   for continous non-parametric covariates in semiparametric spatial regression models.
#'        
#' @description Compute and plot direct, indirect and total impect functions for 
#'   non-parametric covariates included in a semiparametric spatial
#'   or spatio-temporal SAR model. This model must include a spatial
#'   lag of the dependent variable (SAR) to have indirect impacts 
#'   different from 0, otherwise, total and direct function impacts 
#'   are the same.         
#'
#' @param obj \emph{pspatfit} object fitted using \code{\link{pspatfit}} function. 
#' @param listw should be a spatial neighbours list object created for example by \code{nb2listw} from \code{spdep} package. 
#' It can also be a spatial weighting matrix of order (NxN) instead of a listw neighbours list object.
#' @param conflevel numerical value for the confidence interval of the
#'    impact functions. Default 0.95.
#' @param viewplot Default \code{TRUE} to plot impacts. If FALSE use \code{\link{plot_impactsnopar}} to plot impacts
#' @param smooth Default \code{TRUE}. Whether to smooth fitted impacts or not.
#' @param span span for the kernel of the smoothing (see \code{\link{loess}} 
#'             for details). Default c(0.1, 0.1, 0.2). 
#'             
#' @details DESCRIBE ALGORITHM TO COMPUTE impECT FUNCTIONS AND THE 
#'          SMOOTHING TO PLOT        
#'
#' @return A list including
#'   \tabular{ll}{
#'     \emph{impnopar_tot} \tab Matrix including total impects in columns. \cr
#'     \emph{impnopar_dir} \tab Matrix including direct impects in columns. \cr
#'     \emph{impnopar_ind} \tab Matrix including indirect impects in columns. \cr
#'     \emph{impnopar_tot_up} \tab Matrix including upper bounds of total impects in columns. \cr
#'     \emph{impnopar_dir_up} \tab Matrix including upper bounds of direct impects in columns. \cr
#'     \emph{impnopar_ind_up} \tab Matrix including upper bounds of indirect impects in columns. \cr
#'     \emph{impnopar_tot_low} \tab Matrix including lower bounds of total impects in columns. \cr
#'     \emph{impnopar_dir_low} \tab Matrix including lower bounds of direct impects in columns. \cr
#'     \emph{impnopar_ind_low} \tab Matrix including lower bounds of indirect impects in columns. \cr
#'  }
#'         
#' @author Roman Minguez \email{roman.minguez@@uclm.es}
#'
#' @family Direct, Indirect and Total impects.
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{pspatfit}} estimate spatial or spatio-temporal semiparametric PS-SAR
#'     regression models.
#'   \item \code{\link{impactspar}} compute and simulate total, direct and indirect impect
#'     (or impacts) for parametric continuous covariates.
#'   \item \code{\link{fit_terms}} compute terms for smooth functions for non-parametric
#'     continuous covariates and for non-parametric trends.
#'   \item \code{\link{plot_impactsnopar}} plot the non-parametric impects functions
#'     allowing for previous smoothing.
#' }
#' 
#' @references 
#' \itemize{ 
#'   \item Basile, R., Durbán, M., Mínguez, R., Montero, J.
#'         M., and Mur, J. (2014). Modeling regional economic 
#'         dynamics: Spatial dependence, spatial heterogeneity and 
#'         nonlinearities. \emph{Journal of Economic Dynamics and 
#'         Control}, (48), 229-245.
#'         
#'   \item LeSage, J. and Pace, K. (2009). \emph{Introduction to 
#'         Spatial Econometrics}. CRC Press, Boca Raton.
#'         
#'   \item Mínguez, R.; Basile, R. and Durbán, M. (2018). An Alternative Semiparametric Model
#'         for Spatial Panel Data. Under evaluation in \emph{Statistical
#'         Methods and Applications}.
#'  }       
#'
#' @examples
#' ################################################
#'  ###################### Examples using a panel data of rate of
#'  ###################### unemployment for 103 Italian provinces in period 1996-2014.
#' library(pspatreg)
#' data(unemp_it); Wsp <- Wsp_it
#' 
#' ######################  No Spatial Trend: PSAR including a spatial 
#' ######################  lag of the dependent variable
#' form1 <- unrate ~ partrate + agri + cons +
#'                  pspl(serv,nknots=15) +
#'                  pspl(empgrowth,nknots=20) 
#'  gamsar <- pspatfit(form1, data = unemp_it, type = "sar", listw = Wsp_it)
#'  summary(gamsar)
#'  ###### Non-Parametric Total, Direct and Indirect impacts
#'  imp_nparvar <- impactsnopar(gamsar, listw = Wsp_it, viewplot = TRUE)
#'  
#' ######################   PSAR-ANOVA with spatial trend
#' form2 <- unrate ~ partrate + agri + cons +
#'                   pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
#'                   pspt(long,lat,nknots=c(20,20),psanova=TRUE,
#'                   nest_sp1=c(1,2),nest_sp2=c(1,2))
#' ##### Spatial trend fixed for period 1996-2014
#' geospanova_sar <- pspatfit(form2, data=unemp_it, listw = Wsp_it, type = "sar", 
#'                            control=list(tol = 1e-1))
#' summary(geospanova_sar)
#'  ###### Non-Parametric Total, Direct and Indirect impacts
#'  imp_nparvar2 <- impactsnopar(geospanova_sar, listw = Wsp_it, viewplot = TRUE)
#'  
#' ######################   PSAR-ANOVA with spatio-temporal trend and 
#' ######################   temporal autorregresive noise
#'  form3 <- unrate ~ partrate + agri + cons +
#'                    pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
#'                    pspt(long,lat,year,nknots=c(18,18,8),psanova=TRUE,
#'                    nest_sp1=c(1,2,3),nest_sp2=c(1,2,3),
#'                    nest_time=c(1,2,2),ntime=19)
#' sptanova_sar_ar1 <- pspatfit(form3, data = unemp_it, listw = Wsp_it, 
#'                               type = "sar", cor = "ar1",
#'                               control=list(tol=1e-1))
#' summary(sptanova_sar_ar1)
#'  ###### Non-Parametric Total, Direct and Indirect impacts
#'  imp_nparvar3 <- impactsnopar(geospanova_sar, listw = Wsp_it, viewplot = TRUE)
#'
#' @keywords Indirect impects, Direct impects, SAR, non-parametric covariates.
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



