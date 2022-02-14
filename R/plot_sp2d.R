#' @name plot_sp2d 
#' @rdname plot_sp2d 
#'
#' @title Plot and mapping spatial trends        
#' @description  DESCRIBE THE FUNCTION...
#' @param object object returned from \code{\link{pspatfit}} 
#' @param data either sf or dataframe with the data. 
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
#' library(AmesHousing)# The Ames housing data
#' library(tidymodels)
#' library(tune)
#' library(spdep)
#' library(ggplot2)
#' library(dplyr)
#' library(pspatreg)
#' library(tripack)
#' library(rgeos)
#' library(dbscan)
#' 
#' ######## getting and preparing the data
#' ames <- make_ames()# Raw Ames Housing Data
#' ames_sf <- sf::st_as_sf(ames, coords = c("Longitude", "Latitude"))
#' ames_sf$Longitude <- ames$Longitude
#' ames_sf$Latitude <- ames$Latitude
#' ames_sf$lnSale_Price <- log(ames_sf$Sale_Price)
#' ames_sf$lnLot_Area <- log(ames_sf$Lot_Area)
#' ames_sf$lnTotal_Bsmt_SF <- log(ames_sf$Total_Bsmt_SF+1)
#' ames_sf$lnGr_Liv_Area <- log(ames_sf$Gr_Liv_Area)
#' 
#' ########### Constructing the spatial weights matrix
#' ames_sf1 <- ames_sf[duplicated(ames_sf$Longitude)==FALSE,]
#' coord_sf1 <- cbind(ames_sf1$Longitude,ames_sf1$Latitude)
#' ID <- row.names(as(ames_sf1, "sf"))
#' col.tri.nb <- tri2nb(coord_sf1)
#' soi_nb <- graph2nb(soi.graph(col.tri.nb,coord_sf1), row.names=ID)
#' is.symmetric.nb(soi_nb,verbose=T,force=F)
#' listW <- nb2listw(soi_nb, style="W",zero.policy=F)
#' 
#' 
#' ######## formula of the model
#' form2d <- lnSale_Price ~ Fireplaces + Garage_Cars +
#' pspl(lnLot_Area, nknots = 20) + 
#'  pspl(lnTotal_Bsmt_SF, nknots = 20) +
#'  pspl(lnGr_Liv_Area, nknots = 20) +
#'  pspt(Longitude,Latitude, nknots = c(10, 10), 
#'       psanova = FALSE)
#' 
#' ######## fit the model
#' sp2dsar <- pspatfit(form2d, data = ames_sf1, 
#' listw = listW, method = "Chebyshev", 
#' type = "sar", control = list(trace = FALSE))
#'
#' summary(sp2dsar)
#' anova(sp2dsar)
#' 
#' ####### plot spatial trend for spatial point coordinate
#' plot_sp2d(sp2dsar, data = ames_sf1)
#' 
#' ###### MODEL WITH ANOVA DESCOMPOSITION
#'  form2d_psanova <- lnSale_Price ~ Fireplaces + Garage_Cars +
#'  pspl(lnLot_Area, nknots = 20) + 
#'  pspl(lnTotal_Bsmt_SF, nknots = 20) +
#'  pspl(lnGr_Liv_Area, nknots = 20) +
#'  pspt(Longitude,Latitude, nknots = c(10, 10), 
#'       psanova = TRUE)
#'       
#' sp2danovasar <- pspatfit(form2d_psanova, data = ames_sf1, 
#'                         listw = listW,method = "Chebyshev", 
#'                         type = "sar", control = list(trace = FALSE))
#'                         
#' ###### PLOT ANOVA DESCOMPOSITION MODEL
#' plot_sp2d(sp2danovasar, data = ames_sf1, addmain = TRUE, 
#' addint = TRUE)
#' 
#' @export
plot_sp2d <- function(object, data, 
                      coordinates = NULL,
                      npoints = 300, 
                      cexpoints = 0.25,
                      addcontour = TRUE,
                      addpoints = TRUE,
                      addmain = TRUE,
                      addint = TRUE) {
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
  if (!(inherits(data, "sf")) && is.null(coordinates)) 
    stop("data must be an sf object or 
          spatial coordinates must be provided as a matrix")
  if (!(inherits(data, "sf")) && (!is.null(coordinates))) {
    # Build an sf object with points geometry
    data$X <- coordinates[, 1]
    data$Y <- coordinates[, 2]
    data <- st_as_sf(data, coords = c("X", "Y"))
  }
  nfull <- nrow(data)
  nt <- object$nt
  dynamic <- object$dynamic
  if (dynamic) {
    idxyear1 <- seq(from = 1, to = nfull, by = nt)
    data <- data[-idxyear1, ]
  }  
  sp2dfitl <- fit_terms(object, "spttrend") # list object
  geom_data <- st_geometry_type(data)
  if(sum(geom_data == "POINT") == length(geom_data))
    type_geom <- "POINT" else type_geom <- "NO_POINT"
  if (type_geom == "POINT") { ## Case 1: Spatial points
    spco <- st_coordinates(data)
    sp1 <- spco[, 1]
    sp2 <- spco[, 2]
    if (!(object$psanova)) {
      sp2dtrend <- sp2dfitl$fitted_terms
      X <- seq(min(sp1), max(sp1), length = npoints)
      Y <- seq(min(sp2), max(sp2), length = npoints)
      dfintp <- interp(sp1, sp2, sp2dtrend,
                              xo = X, yo = Y,
                              linear = FALSE,
                              duplicate = "median",
                              extrap = FALSE)
      image(X, Y, dfintp$z)
      if (addcontour)
        contour(X, Y, dfintp$z, add = TRUE)
      if (addpoints)
        points(sp1, sp2, cex = cexpoints)
      title(main = "Spatial Trend")
    } else { # psanova case
      sp2dtrend <- sp2dfitl$fitted_terms[, "spttrend"]
      f1_main <- sp2dfitl$fitted_terms[, "f1_main"]
      f2_main <- sp2dfitl$fitted_terms[, "f2_main"]
      f12_int <- sp2dfitl$fitted_terms[, "f12_int"] +
                 sp2dfitl$fitted_terms[, "Intercept"]
      X <- seq(min(sp1), max(sp1), length = npoints)
      Y <- seq(min(sp2), max(sp2), length = npoints)
      df <- data.frame(sp1 = sp1, sp2 = sp2,
                       sp2dtrend = sp2dtrend,
                       f1_main = f1_main,
                       f2_main = f2_main,
                       f12_int = f12_int)
      dfintp <- interp(sp1, sp2, sp2dtrend,
                              xo = X, yo = Y,
                              linear = FALSE,
                              duplicate = "median",
                              extrap = FALSE)
      image(X, Y, dfintp$z)
      if (addcontour)
        contour(X, Y, dfintp$z, add = TRUE)
      if (addpoints)
        points(sp1, sp2, cex = cexpoints)
      title(main = "Spatial Trend")
      if (addmain) {
        readline(prompt="Press [enter] to continue")
        f1plot <- ggplot(data = df, 
                         mapping = aes(x = sp1,
                                       y = f1_main)) +
                  geom_line() +
                  ggtitle("Spat. Trend: f1_main")
        print(f1plot)
        readline(prompt="Press [enter] to continue")
        f2plot <- ggplot(data = df,
                         mapping = aes(x = sp2,
                                       y = f2_main)) +
                  geom_line() +
                  ggtitle("Spat. Trend: f2_main")
        print(f2plot)
      }
      if (addint) {
        readline(prompt="Press [enter] to continue")
        dfintp <- interp(sp1, sp2, f12_int,
                                xo = X, yo = Y,
                                linear = FALSE,
                                duplicate = "median",
                                extrap = FALSE)
        image(X, Y, dfintp$z)
        if (addcontour)
          contour(X, Y, dfintp$z, add = TRUE)
        if (addpoints)
          points(sp1, sp2, cex = cexpoints)
        title(main = "Spat. Trend: f12_int")
      }
     }
  } else { # NO_POINT
    if (!(object$psanova)) {
      data$sp2dtrend <- sp2dfitl$fitted_terms
      df <- data[, c("sp2dtrend")]
      plot(df, main = "Spatial Trend")
    } else {
      data$sp2dtrend <- sp2dfitl$fitted_terms[, "spttrend"]
      df <- data[, c("sp2dtrend")]
      plot(df, main = "Spatial Trend")
      if (addmain) {
        readline(prompt="Press [enter] to continue")
        data$f1_main <- sp2dfitl$fitted_terms[, "f1_main"]
        df <- data[, c("f1_main")]
        plot(df, main = "Spat. Trend: f1_main")
        readline(prompt="Press [enter] to continue")
        data$f2_main <- sp2dfitl$fitted_terms[, "f2_main"]
        df <- data[, c("f2_main")]
        plot(df, main = "Spat. Trend: f2_main")

      }
      if (addint) {
        readline(prompt="Press [enter] to continue")
        data$f12_int <- sp2dfitl$fitted_terms[, "f12_int"] +
          sp2dfitl$fitted_terms[, "Intercept"]
        df <- data[, c("f12_int")]
        plot(df, main = "Spat. Trend: f12_int")
      }
    }
  }
}
