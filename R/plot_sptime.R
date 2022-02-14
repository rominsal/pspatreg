#' @name plot_sptime 
#' @rdname plot_sptime 
#'
#' @title Plot of time trend for PS-ANOVA models in 3d. 
#' @description DESCRIBE THE FUNCTION...
#' @param object object returned from \code{\link{pspatfit}} 
#' @param data either sf or dataframe with the data. 
#' @param time_var name of the temporal variable in data. 
#' @param reg_var name of the regional variable in data.
#' @return time series plots of the temporal trend for each region                                  
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
#' 
#' ###### FORMULA OF THE MODEL
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
#' ######## PLOT OF TEMPORAL TREND FOR EACH PROVINCE 
#' plot_sptime(sp3danovasarar1, data = unemp_it, 
#' time_var = "year", 
#' reg_var = "prov")
#' 
#' @export
plot_sptime <- function(object, data, time_var, reg_var) {
  if (!(inherits(object, "pspatreg"))) 
    stop("object must be of class pspatreg")
  if (!(object$psanova))
    stop("object must include a psanova spatio-temporal trend")
  # if (!(inherits(data, "sf"))) 
  #   stop("data must be an sf object")
  if (!(time_var %in% names(data)))
    stop("time_var must be between the names of data")
  if (!(reg_var %in% names(data)))
    stop("reg_var must be between the names of data")
  pos_time_var <- which(time_var == names(data))
  time_var <- data[, pos_time_var]
  if (inherits(time_var, "sf"))
    time_var <- st_drop_geometry(time_var)
  if (inherits(time_var, "data.frame"))
    time_var <- time_var[, 1, drop = TRUE]
  pos_reg_var <- which(reg_var == names(data))
  reg_var <- data[, pos_reg_var]
  if (inherits(reg_var, "sf"))
    reg_var <- st_drop_geometry(reg_var)
  if (inherits(reg_var, "data.frame"))
    reg_var <- reg_var[, 1, drop = TRUE]
  time_var <- as.factor(time_var)
  nt <- nlevels(time_var)  
  dynamic <- object$dynamic
  reg_var <- as.factor(reg_var)
  nsp <- nlevels(reg_var)
  if (dynamic) {
    # Eliminate first year...
    year1 <- levels(time_var)[1]
    idxyear1 <- time_var == year1
    time_var <- time_var[!idxyear1]
    reg_var <- reg_var[!idxyear1]
    nt <- nt - 1
  }
  sp3dfitl <- fit_terms(object, "spttrend") # list object
  ft <- sp3dfitl$fitted_terms[, "ft_main"] +
        sp3dfitl$fitted_terms[, "f1t_int"] +
        sp3dfitl$fitted_terms[, "f2t_int"] +
        sp3dfitl$fitted_terms[, "f12t_int"]
  df <- data.frame(region = reg_var,
                   time = time_var,
                   ft = ft)
  regions <- levels(reg_var)
  matdata <- matrix(NA, nrow = nt, ncol = nsp)
  for (i in seq_along(levels(reg_var))) {
    reg_i <- levels(reg_var)[i]
    matdata[, i] <- df[df$region == reg_i, "ft"] 
  }
  colnames(matdata) <- levels(reg_var)
  if (dynamic) 
    time <- as.integer(levels(time_var)[2:nlevels(time_var)])
  else
    time <- as.integer(levels(time_var))
  mattime <- matrix(rep(time, nsp), nrow = nt, 
                     ncol = nsp, byrow = FALSE)
  colnames(mattime) <- rep("Time", nsp)
  matplot(mattime, matdata, type = "l",
          xlab = "Time", ylab = "ft", 
          main = "Temporal Trend for each region")
}
