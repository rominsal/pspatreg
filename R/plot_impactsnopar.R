#' @name plot_impactsnopar
#' @rdname plot_impactsnopar
#'
#' @title Plot direct, indirect and total impacts  functions 
#'   for continous non-parametric covariates in semiparametric spatial regression models.
#'        
#' @description Plot direct, indirect and total impacts functions for 
#'   non-parametric covariates included in a semiparametric spatial
#'   or spatio-temporal SAR model. This model must include a spatial
#'   lag of the dependent variable (SAR) to have indirect effects 
#'   different from 0, otherwise, total and direct function effects 
#'   are the same. The effect functions can be smoothed to overcome 
#'   the instabilities created by the premultiplication of matrix
#'   \eqn{(I - \rho W)^{-1}} 
#'
#' @param impactsnopar object returned from \code{\link{impactsnopar}} function.
#' @param data dataframe with the data. 
#' @param smooth logical value to choose smoothing of the effects function
#'               prior to plot. Default TRUE.
#' @param span span for the kernel of the smoothing (see \code{\link{loess}} 
#'             for details). Default c(0.1, 0.1, 0.2). 
#' @param dynamic Logical value to set a dynamic model.
#'   Dynamic models include a temporal lag of the dependent
#'   variable in the right-hand side of the equation.
#'   Default = `FALSE`.
#' @param nt  Number of temporal periods. It is needed
#'   for dynamic models.  
#'  
#' @return plot of the direct, indirect and total impacts  function for each non-parametric
#'   covariate included in the object returned from \code{\link{impactsnopar}}.
#'                                 
#' @author 
#' \tabular{ll}{ 
#'   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Roberto Basile \tab \email{roberto.basile@@univaq.it} \cr 
#'   Maria Durban \tab \email{mdurban@@est-econ.uc3m.es} \cr 
#'   Gonzalo Espana-Heredia \tab \email{gehllanza@@gmail.com} \cr 
#'  }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{impactsnopar}} compute total, direct and indirect effect
#'           functions for non-parametric continuous covariates.
#'   \item \code{\link{fit_terms}} compute smooth functions for non-parametric
#'                                 continuous covariates.
#'   \item \code{\link{plot_terms}} plot the terms of non-parametric covariates.
#' }
#' 
#' @references \itemize{ 
#'     \item Basile, R.; Durban, M.; Minguez, R.; Montero, J. M.; and 
#'     Mur, J. (2014). Modeling regional economic dynamics: Spatial
#'     dependence, spatial heterogeneity and nonlinearities. 
#'     \emph{Journal of Economic Dynamics and Control}, (48), 229-245.
#'     <doi:10.1016/j.jedc.2014.06.011>               
#'  }
#'         
#' @examples
#' ################################################
#' # Examples using spatial data of Ames Houses.
#' ###############################################
#' # Getting and preparing the data
#' library(pspatreg)
#' library(spdep)
#' library(sf)
#' ames <- AmesHousing::make_ames() # Raw Ames Housing Data
#' ames_sf <- st_as_sf(ames, coords = c("Longitude", "Latitude"))
#' ames_sf$Longitude <- ames$Longitude
#' ames_sf$Latitude <- ames$Latitude
#' ames_sf$lnSale_Price <- log(ames_sf$Sale_Price)
#' ames_sf$lnLot_Area <- log(ames_sf$Lot_Area)
#' ames_sf$lnTotal_Bsmt_SF <- log(ames_sf$Total_Bsmt_SF+1)
#' ames_sf$lnGr_Liv_Area <- log(ames_sf$Gr_Liv_Area)
#' ########### Constructing the spatial weights matrix
#' ames_sf1 <- ames_sf[(duplicated(ames_sf$Longitude) == FALSE), ]
#' coord_sf1 <- cbind(ames_sf1$Longitude, ames_sf1$Latitude)
#' ID <- row.names(as(ames_sf1, "sf"))
#' col_tri_nb <- tri2nb(coord_sf1)
#' soi_nb <- graph2nb(soi.graph(col_tri_nb, 
#'                             coord_sf1), 
#'                    row.names = ID)
#' lw_ames <- nb2listw(soi_nb, style = "W", 
#'                     zero.policy = FALSE)
#'                     
#' form1 <- lnSale_Price ~ Fireplaces + Garage_Cars +
#'           pspl(lnLot_Area, nknots = 20) + 
#'           pspl(lnTotal_Bsmt_SF, nknots = 20) +
#'           pspl(lnGr_Liv_Area, nknots = 20)    
#' 
#' gamsar <- pspatfit(form1, data = ames_sf1, 
#'                    type = "sar", listw = lw_ames,
#'                    method = "Chebyshev")
#' summary(gamsar)
#' \donttest{ 
#' nparimpacts <- impactsnopar(gamsar, listw = lw_ames, viewplot = FALSE)
#' plot_impactsnopar(nparimpacts, data = ames_sf1, smooth = TRUE)
#' ###### Examples using a panel data of rate of
#' ###### unemployment for 103 Italian provinces in period 1996-2014.
#' library(pspatreg)
#' data(unemp_it)
#' ## Wsp_it is a matrix. Create a neighboord list 
#' lwsp_it <- spdep::mat2listw(Wsp_it)
#' ## short sample for spatial pure case (2d)
#' ########  No Spatial Trend: PSAR including a spatial 
#' ########  lag of the dependent variable
#' form1 <- unrate ~ partrate + agri + cons + empgrowth +
#'                   pspl(serv, nknots = 15) 
#' gamsar <- pspatfit(form1, data = unemp_it, 
#'                    type = "sar", 
#'                    listw = lwsp_it)
#' summary(gamsar)
#' ###### Non-Parametric Total, Direct and Indirect impacts
#' imp_nparvar <- impactsnopar(gamsar, alpha = 0.05,
#'                             listw = lwsp_it, 
#'                             viewplot = TRUE)  
#' ##### This returns the same result but using plot_impactsnopar()
#' imp_nparvar <- impactsnopar(gamsar, listw = lwsp_it, alpha = 0.05,
#'                             viewplot = FALSE)
#' plot_impactsnopar(imp_nparvar, data = unemp_it, 
#'                    smooth = TRUE)
#' }                    
#'
#' @export
plot_impactsnopar <- function(impactsnopar, data, smooth = TRUE, 
                           span = c(0.1, 0.1, 0.2),
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
 tot <- impactsnopar$impnopar_tot
 if (nrow(tot) != nrow(data)) 
   stop("Dimensions of impacts and data disagree. 
        Likely the model is dynamic and arguments dynamic = TRUE and nt value
        are needed in plot_impactsnopar")
 uptot <- impactsnopar$impnopar_tot_up
 lowtot <- impactsnopar$impnopar_tot_low
 dir <- impactsnopar$impnopar_dir
 updir <- impactsnopar$impnopar_dir_up
 lowdir <- impactsnopar$impnopar_dir_low
 ind <- impactsnopar$impnopar_ind
 upind <- impactsnopar$impnopar_ind_up
 lowind <- impactsnopar$impnopar_ind_low
 for (i in 1:ncol(tot)) {
    name_var <- colnames(tot)[i]
    var <- as.matrix(data[, c(name_var)])
    colnames(var) <- name_var
    ord <- order(var)
    tot_i <- matrix(tot[, c(name_var)], ncol = 1)
    colnames(tot_i) <- name_var
    uptot_i <- matrix(uptot[, c(name_var)], ncol = 1)
    colnames(uptot_i) <- name_var
    lowtot_i <- matrix(lowtot[, c(name_var)], ncol = 1)
    colnames(lowtot_i) <- name_var
    dir_i <- matrix(dir[, c(name_var)], ncol = 1)
    colnames(dir_i) <- name_var
    updir_i <- matrix(updir[,c(name_var)], ncol = 1)
    colnames(updir_i) <- name_var
    lowdir_i <- matrix(lowdir[,c(name_var)], ncol = 1)
    colnames(lowdir_i) <- name_var
    ind_i <- matrix(ind[,c(name_var)], ncol = 1)
    colnames(ind_i) <- name_var
    upind_i <- matrix(upind[,c(name_var)], ncol = 1)
    colnames(upind_i) <- name_var
    lowind_i <- matrix(lowind[,c(name_var)], ncol = 1)
    colnames(lowind_i) <- name_var
    if (smooth) {
      span_tot <- span[1]
      span_dir <- span[2]
      span_ind <- span[3]
      tot_i_smooth <- predict(loess(tot_i ~ var, span = span_tot),
                            method = "loess()")
      uptot_i_smooth <- predict(loess(uptot_i ~ var, span = span_tot),
                              method = "loess()")
      lowtot_i_smooth <- predict(loess(lowtot_i ~ var, span = span_tot),
                                method = "loess()")
      if (sum(is.nan(tot_i_smooth)) > 1 ||
          sum(is.nan(uptot_i_smooth)) > 1 ||
          sum(is.nan(lowtot_i_smooth)) > 1 ) {
        warning(paste("Smoothing of total impacts with variable ", name_var,
                  " produces NaN. ",
                   "This variable is not suited for smoothing \n", sep = ""))
      } else {
        tot_i <- tot_i_smooth
        uptot_i <- uptot_i_smooth
        lowtot_i <- lowtot_i_smooth
      }
      dir_i_smooth <- predict(loess(dir_i ~ var, span = span_dir),
                              method = "loess()")
      updir_i_smooth <- predict(loess(updir_i ~ var, span = span_dir),
                                method = "loess()")
      lowdir_i_smooth <- predict(loess(lowdir_i ~ var, span = span_dir),
                                 method = "loess()")
      if (sum(is.nan(dir_i_smooth)) > 1 ||
          sum(is.nan(updir_i_smooth)) > 1 ||
          sum(is.nan(lowdir_i_smooth)) > 1 ) {
        warning(paste("Smoothing of direct impacts with variable ", name_var,
                  " produces NaN. ",
                  "This variable is not suited for smoothing \n", sep = ""))
      } else {
        dir_i <- dir_i_smooth
        updir_i <- updir_i_smooth
        lowdir_i <- lowdir_i_smooth
      }      
      ind_i_smooth <- predict(loess(ind_i ~ var, span = span_ind),
                              method = "loess()")
      upind_i_smooth <- predict(loess(upind_i~var, span = span_ind),
                                method = "loess()")
      lowind_i_smooth <- predict(loess(lowind_i ~ var, span = span_ind),
                                 method = "loess()")
      if (sum(is.nan(ind_i_smooth)) > 1 ||
          sum(is.nan(upind_i_smooth)) > 1 ||
          sum(is.nan(lowind_i_smooth)) > 1 ) {
        warning(paste("Smoothing of indirect impacts with variable ", name_var,
                  " produces NaN. ",
                  "This variable is not suited for smoothing \n", sep = ""))
      } else {
        ind_i <- ind_i_smooth
        upind_i <- upind_i_smooth
        lowind_i <- lowind_i_smooth
      }       
    }
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(3, 1))
    plot(var[ord], tot_i[ord], 
         type = "l",
         ylab = paste("f(", name_var, ")"), 
         xlab = name_var,
         ylim = c(min(lowtot_i), max(uptot_i)), 
         cex.lab = 1.0, 
         col = 2, 
         lty = 1, 
         lwd = 2, 
         cex.main = 1.0, 
         main = paste("Total Impacts: ", name_var))
    lines(var[ord], uptot_i[ord], 
          xlab = "", ylab = "", 
          type = "l", col = 2, lty = 2, lwd = 1.5)
    lines(var[ord],lowtot_i[ord], 
          xlab = "", 
          ylab = "",
          type = "l", 
          col = 2, 
          lty = 2, 
          lwd = 1.5)
    abline(a = 0, b = 0)
    plot(var[ord], dir_i[ord], 
         type = "l",
         ylab = paste("f(", name_var, ")"), 
         xlab = name_var,
         ylim = c(min(lowdir_i), max(updir_i)), 
         cex.lab = 1.0, 
         col = 3,
         lty = 1, 
         lwd = 2,
         cex.main = 1.0, 
         main = paste("Direct Impacts:", name_var))
    lines(var[ord], updir_i[ord], 
          xlab = "", 
          ylab = "", 
          type = "l", 
          col = 3, 
          lty = 2, 
          lwd = 1.5)
    lines(var[ord],lowdir_i[ord], 
          xlab = "", 
          ylab = "", 
          type = "l", 
          col = 3, 
          lty = 2, 
          lwd = 1.5)
    abline(a = 0, b = 0)
    plot(var[ord], ind_i[ord], 
         type = "l",
         ylab = paste("f(", name_var, ")"), 
         xlab = name_var,
         ylim = c(min(lowind_i), max(upind_i)),
         cex.lab = 1.0, 
         col = 4,
         lty = 1,
         lwd = 2,
         cex.main = 1.0, 
         main = paste("Indirect Impacts: ", name_var))
    lines(var[ord], upind_i[ord], 
          xlab = "",
          ylab = "",
          type = "l",
          col = 4, 
          lty = 2, 
          lwd = 1.5)
    lines(var[ord], lowind_i[ord], 
          xlab = "", 
          ylab = "", 
          type = "l",
          col = 4, 
          lty = 2,
          lwd = 1.5)
    abline(a = 0, b = 0)
    readline(prompt="Press [enter] to continue")
 }
}
