#' @name plot_sp3d 
#' @rdname plot_sp3d
#'
#' @title Plot and mapping spatio-temporal trends        
#' @description DESCRIBE THE FUNCTION...
#' @param object object returned from \code{\link{pspatfit}} 
#' @param data sf object. 
#' @param time_var name of the temporal variable in data.
#' @param time_index vector of time points to plot. 
#' @param coordinates coordinates matrix if *data* is not an sf object.
#' @param npoints number of points to use in the interpolation.
#' @param addcontour Logical value to add contour lines.
#' @param addpoints Logical value to add spatial points to the graphics.
#' @param addmain Add f1_main and f2_main plots in psanova case.
#' @param addint Add f12_int in psanova case.
#' @return plots and maps of the spatial trends                                  
#' @author Roman Minguez \email{roman.minguez@@uclm.es}
#' @examples
#' 
#' library(pspatreg)
#' data("unemp_it_new_polyg_sf")
#' data("unemp_it")
#' unemp_it <- unemp_it_new_polyg_sf
#' rm(unemp_it_new_polyg_sf) # Only to be consistent with the names...
#' colnames(Wsp_it) <- rownames(Wsp_it)
#' Wsp_it <- as.matrix(Wsp_it)
#' Using spdep we transform the spatial weights matrix in a list of neighbours object:
#'  
#'  lwsp_it <- spdep::mat2listw(Wsp_it, 
#'                              row.names = rownames(Wsp_it))
#' 
#' summary(lwsp_it)
#' ######## FORMULA of the model
#' form3d_psanova_restr <- unrate ~ partrate + agri + cons +
#' pspl(serv, nknots = 15) + 
#'  pspl(empgrowth, nknots = 20) +
#'  pspt(long, lat, year, 
#'       nknots = c(18,18,8),
#'       psanova = TRUE, 
#'       nest_sp1 = c(1, 2, 3), 
#'       nest_sp2 = c(1, 2, 3),
#'       nest_time = c(1, 2, 2),
#'       f1t = FALSE, f2t = FALSE)
#' 
#' ####### FIT the model
#' sp3danovasarar1 <- pspatfit(form3d_psanova_restr, 
#' data = unemp_it, 
#' listw = lwsp_it, 
#' method = "Chebyshev", 
#' type = "sar",
#' cor = "ar1")
#' 
#' ###### PLOT spatio-temporal trends
#' plot_sp3d(sp3danovasarar1, data = unemp_it, 
#' time_var = "year", 
#' time_index = c(1996, 2005, 2019),
#' addmain = FALSE, addint = FALSE)
#' 
#' @export
plot_sp3d <- function(object, data,
                      time_var,
                      time_index,
                      addmain = TRUE,
                      addint = TRUE) {
  ## VIP: Assumption sf object...
## The function needs this arguments:
  # 1. A fitted object of class "pspatreg"
  # 2. A database that could be an sf object or a dataframe.
  # 3. In the last case (dataframe) the names of the
  ##   spatial coordinates need to be supplied.
  ## THE FUNCTION NEED TO BE DOCUMENTED IN THE USUAL WAY INCLUDING EXAMPLES.
  ## THE EXAMPLES CAN BE EXTRACTED FROM "Examples_plots_maps_spatialtrends.Rmd"
  ######################################################################
  if (!(inherits(object, "pspatreg"))) 
    stop("object must be of class pspatreg")
  if (!(inherits(data, "sf"))) 
    stop("data must be an sf object")
  if (!(time_var %in% names(data)))
    stop("time_var must be between the names of data")
  dynamic <- object$dynamic
  nfull <- nrow(data)
  nt <- object$nt
  if (dynamic) {
    idxyear1 <- seq(from = 1, to = nfull, by = nt)
    data <- data[-idxyear1, ]
  }  
  pos_time_var <- which(time_var == names(data))
  time_var <- data[, pos_time_var]
  if (inherits(time_var, "sf"))
    time_var <- st_drop_geometry(time_var)
  sp3dfitl <- fit_terms(object, "spttrend") # list object
  data$sp3dtrend <- sp3dfitl$fitted_terms[, "spttrend"]
  if (object$psanova) {
    data$f1_main <- sp3dfitl$fitted_terms[, "f1_main"]
    data$f2_main <- sp3dfitl$fitted_terms[, "f2_main"]
    data$f12_int <- sp3dfitl$fitted_terms[, "f12_int"] +
                    sp3dfitl$fitted_terms[, "Intercept"]
  }
  for (i in seq_along(time_index)) {
    year_i <- time_index[i]
    data_i <- data[time_var == year_i, ]
    df_i <- data_i[, c("sp3dtrend")]
    plot(df_i, main = paste("Spatial Trend for : ", year_i))
    if (object$psanova && addmain) {
      readline(prompt="Press [enter] to continue")
      df_i <- data_i[, c("f1_main")]
      plot(df_i, 
           main = paste("Spat. Trend: f1_main for : ", year_i))
      readline(prompt="Press [enter] to continue")
      df_i <- data_i[, c("f2_main")]
      plot(df_i, 
           main =  paste("Spat. Trend: f2_main for : ", year_i))
    }
    if (object$psanova && addint) {
      readline(prompt="Press [enter] to continue")
      df_i <- data_i[, c("f12_int")]
      plot(df_i, 
           main = paste("Spat. Trend: f12_int for : ", year_i))
    }
    readline(prompt = "Press [enter] to continue")
  }
}