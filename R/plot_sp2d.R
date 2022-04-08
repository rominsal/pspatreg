#' @name plot_sp2d 
#' @rdname plot_sp2d 
#'
#' @title Plot and mapping spatial trends.   
#'      
#' @description  Make plots and maps of the spatial trends 
#'   in 2d of the objects fitted with \code{\link{pspatfit}} function.
#'   
#' @param object object returned from \code{\link{pspatfit}} 
#' @param data either sf or dataframe with the data. 
#' @param coordinates coordinates matrix if \code{data} is not an sf object.
#' @param npoints number of points to use in the interpolation.
#' @param cexpoints size of the points. Default = 0.25 
#' @param addcontour Logical value to add contour lines.
#' @param addpoints Logical value to add spatial points to the graphics.
#' @param addmain Add f1_main and f2_main plots in psanova case.
#' @param addint Add f12_int in psanova case.
#' 
#' @return plots and maps of the spatial trends                                  
#' @author 
#' \tabular{ll}{ 
#'   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Roberto Basile \tab \email{roberto.basile@@univaq.it} \cr Maria Durban \tab
#'   \email{mdurban@@est-econ.uc3m.es} \cr Gonzalo Espana-Heredia \tab
#'   \email{gehllanza@@gmail.com} \cr 
#'  }
#'  
#' @references \itemize{ 
#'   \item Lee, D. and Durban, M. (2011). P-Spline ANOVA Type Interaction 
#'     Models for Spatio-Temporal Smoothing. \emph{Statistical Modelling}, 
#'     (11), 49-69. <doi:10.1177/1471082X1001100104>
#'     
#'   \item Eilers, P. and Marx, B. (2021). \emph{Practical Smoothing. 
#'   The Joys of P-Splines}. Cambridge University Press.
#'     
#'   \item Fahrmeir, L.; Kneib, T.;  Lang, S.; and Marx, B. (2021). 
#'     \emph{Regression. Models, Methods and Applications (2nd Ed.)}.
#'      Springer.
#'         
#'   \item Wood, S.N. (2017). \emph{Generalized Additive Models. 
#'   An Introduction with \code{R}} (second edition). CRC Press, Boca Raton.
#'  }   
#'  
#'@examples
#' library(pspatreg, package = "pspatreg")
#' ##### EXAMPLE 2D WITH UNEMPLOYMENT IN ITALIAN PROVINCES
#' ## 103 Italian provinces. Period 1996-2019
#' data(unemp_it, package = "pspatreg")
#' ## Wsp_it is a matrix. Create a neighboord list 
#' lwsp_it <- spdep::mat2listw(Wsp_it)
#' ## short sample for spatial pure case (2d)
#' unemp_it_short <- unemp_it[unemp_it$year == 2019, ]
#' ## Spatial trend with anova decomposition
#' formit2d_anova <- unrate ~ partrate + agri + cons +
#'                            pspl(serv, nknots = 15) + 
#'                            pspl(empgrowth, nknots = 20) +
#'                            pspt(long, lat, nknots = c(20, 20), 
#'                            psanova = TRUE)
#' \donttest{
#' geospit2d_anova <- pspatfit(formit2d_anova, 
#'                             data = unemp_it_short)
#' summary(geospit2d_anova)                             
#' ### Create sf object to make the plot 
#' library(sf)
#' unemp_it_sf_short <- st_as_sf(dplyr::left_join(
#'                               unemp_it_short, 
#'                               map_it,  
#'                         by = c("prov" = "COD_PRO")))
#' plot_sp2d(geospit2d_anova, data = unemp_it_sf_short, 
#'           addmain = TRUE, addint = TRUE)
#' ## Spatial trend without anova decomposition
#' formit2d <- unrate ~ partrate + agri + cons +
#'                      pspl(serv, nknots = 15) + 
#'                      pspl(empgrowth, nknots = 20) +
#'                      pspt(long, lat, nknots = c(20, 20))
#' geospit2d <- pspatfit(formit2d, 
#'                       data = unemp_it_short)
#' summary(geospit2d_anova) 
#' plot_sp2d(geospit2d, data = unemp_it_sf_short)
#' ######## EXAMPLE 2D WITH AMES DATA 
#' ######## getting and preparing the data
#' library(spdep)
#' library(dbscan)
#' ames <- AmesHousing::make_ames() # Raw Ames Housing Data
#' ames_sf <- st_as_sf(ames, coords = c("Longitude", "Latitude"))
#' ames_sf$Longitude <- ames$Longitude
#' ames_sf$Latitude <- ames$Latitude
#' ames_sf$lnSale_Price <- log(ames_sf$Sale_Price)
#' ames_sf$lnLot_Area <- log(ames_sf$Lot_Area)
#' ames_sf$lnTotal_Bsmt_SF <- log(ames_sf$Total_Bsmt_SF+1)
#' ames_sf$lnGr_Liv_Area <- log(ames_sf$Gr_Liv_Area)
#' 
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
#' ######## formula of the model IN AMES
#' form2d <- lnSale_Price ~ Fireplaces + Garage_Cars +
#'           pspl(lnLot_Area, nknots = 20) + 
#'           pspl(lnTotal_Bsmt_SF, nknots = 20) +
#'           pspl(lnGr_Liv_Area, nknots = 20) +
#'           pspt(Longitude, Latitude, 
#'                nknots = c(10, 10), 
#'                psanova = FALSE)
#' 
#' ######## fit the model
#' sp2dsar <- pspatfit(form2d, data = ames_sf1, 
#'                     listw = lw_ames, 
#'                     method = "Chebyshev", 
#'                     type = "sar")
#' summary(sp2dsar)
#' 
#' ####### plot spatial trend for spatial point coordinate
#' plot_sp2d(sp2dsar, data = ames_sf1)
#' 
#' ###### MODEL WITH ANOVA DESCOMPOSITION
#'  form2d_psanova <- lnSale_Price ~ Fireplaces + Garage_Cars +
#'                    pspl(lnLot_Area, nknots = 20) + 
#'                    pspl(lnTotal_Bsmt_SF, nknots = 20) +
#'                    pspl(lnGr_Liv_Area, nknots = 20) +
#'                    pspt(Longitude, Latitude, 
#'                         nknots = c(10, 10), 
#'                         psanova = TRUE)
#'       
#' sp2danovasar <- pspatfit(form2d_psanova, 
#'                         data = ames_sf1, 
#'                         listw = lw_ames,
#'                         method = "Chebyshev", 
#'                         type = "sar")
#' summary(sp2danovasar)                         
#'                         
#' ###### PLOT ANOVA DESCOMPOSITION MODEL
#' plot_sp2d(sp2danovasar, data = ames_sf1, 
#'           addmain = TRUE, addint = TRUE)
#' }           
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
    sp2dtrend <- sp2dfitl$fitted_terms[, "spttrend"]
    sp2dtrend <- sp2dtrend - mean(sp2dtrend)
    min_i <- min(sp2dtrend) 
    max_i <- max(sp2dtrend)
    if (object$psanova) {
      f1_main <- sp2dfitl$fitted_terms[, "f1_main"]
      f2_main <- sp2dfitl$fitted_terms[, "f2_main"]
      f12_int <- sp2dfitl$fitted_terms[, "f12_int"]
      intercept <- sp2dfitl$fitted_terms[, "Intercept"]
      if (addmain) {
        min_i <- min(c(min_i, f1_main, f2_main))
        max_i <- max(c(max_i, f1_main, f2_main))
      }
      if (addint) {
        min_i <- min(c(min_i, f12_int))
        max_i <- max(c(max_i, f12_int))
      }
    }
    range_i <- c(min_i - 0.01, max_i + 0.01)
    breaks_i <- seq(min_i - 0.01, max_i + 0.01, 
                    by = diff(range(range_i))/15)
    if (!(object$psanova)) {
      dfsptrend2d <- data.frame(sp1 = sp1, sp2 = sp2, 
                              sp2dtrend = sp2dtrend)
      interp_trend2d <- mba.surf(dfsptrend2d, 
                                 no.X = npoints,
                                 no.Y = npoints,
                                 extend = TRUE)$xyz
      fields::image.plot(interp_trend2d,
                        breaks = breaks_i,
                        col = heat.colors(
                        n = (length(breaks_i)-1)))
      if (addcontour)
        graphics::contour(interp_trend2d, add = TRUE)
      if (addpoints)
        points(sp1, sp2, cex = cexpoints)
      title(main = "Spatial Trend (centered)")
    } else { # psanova case
      df <- data.frame(sp1 = sp1, sp2 = sp2,
                       sp2dtrend = sp2dtrend,
                       f1_main = f1_main,
                       f2_main = f2_main,
                       f12_int = f12_int,
                       intercept = intercept)
      dfsptrend2d <- df[, c("sp1", "sp2", "sp2dtrend")]
      interp_trend2d <- mba.surf(dfsptrend2d, 
                                 no.X = npoints,
                                 no.Y = npoints,
                                 extend = TRUE)$xyz      
      fields::image.plot(interp_trend2d,
                        breaks = breaks_i,
                        col = heat.colors(
                        n = (length(breaks_i)-1)))
      if (addcontour)
        contour(interp_trend2d, add = TRUE)
      if (addpoints)
        points(sp1, sp2, cex = cexpoints)
      title(main = "Spatial Trend (centered)")
      if (addmain) {
        readline(prompt="Press [enter] to continue")
        f1plot <- ggplot(data = df) +
                  geom_line(mapping = aes(x = sp1,
                                          y = f1_main)) +
                  ggplot2::ylim(min_i, max_i) +
                  ggtitle("Spat. Trend: f1_main")
        print(f1plot)
        readline(prompt="Press [enter] to continue")
        f2plot <- ggplot(data = df,
                         mapping = aes(x = sp2,
                                       y = f2_main)) +
                  geom_line() +
                  ggplot2::ylim(min_i, max_i) +
                  ggtitle("Spat. Trend: f2_main")
        print(f2plot)
      }
      if (addint) {
        readline(prompt="Press [enter] to continue")
        df_f12_int <- df[, c("sp1", "sp2", "f12_int")]
        interp_f12_int <- mba.surf(df_f12_int, 
                                   no.X = npoints,
                                   no.Y = npoints,
                                   extend = TRUE)$xyz
        fields::image.plot(interp_f12_int,
                           breaks = breaks_i,
                           col = heat.colors(
                             n = (length(breaks_i)-1)))
        if (addcontour)
          contour(interp_f12_int, add = TRUE)
        if (addpoints)
          points(sp1, sp2, cex = cexpoints)
        title(main = "Spat. Trend: f12_int")
      }
     }
  } else { # NO_POINT
    data$sp2dtrend <- sp2dfitl$fitted_terms[, "spttrend"]
    data$sp2dtrend <- data$sp2dtrend - 
                        mean(data$sp2dtrend)
    df <- data[ ,"sp2dtrend"]
    min_i <- min(df$sp2dtrend)
    max_i <- max(df$sp2dtrend)
    if (!(object$psanova)) {
      plot(df, main = "Spatial Trend (centered)")
    } else { # object$psanova == TRUE
      data$f1_main <- sp2dfitl$fitted_terms[, "f1_main"]
      data$f2_main <- sp2dfitl$fitted_terms[, "f2_main"]
      data$f12_int <- sp2dfitl$fitted_terms[, "f12_int"]
      data$intercept <- sp2dfitl$fitted_terms[, "Intercept"]
      # Recompute min_i and max_i
      if (addmain) {
        min_i <- min(c(min_i, data$f1_main, 
                       data$f2_main))
        max_i <- max(c(max_i, data$f1_main, 
                       data$f2_main))
      } 
      if (addint) {
        min_i <- min(c(min_i, data$f12_int))
        max_i <- max(c(max_i, data$f12_int))
      }
      range_i <- c(min_i - 0.01, max_i + 0.01)
      breaks_i <- seq(min_i - 0.01, max_i + 0.01, 
                      by = diff(range(range_i))/5)
      plot(df, breaks = breaks_i,
           main = "Spatial Trend (centered)")
      if (addmain) {
        readline(prompt="Press [enter] to continue")
        df <- data[, c("f1_main")]
        plot(df, breaks = breaks_i ,
             main = "Spat. Trend: f1_main")
        readline(prompt="Press [enter] to continue")
        df <- data[, c("f2_main")]
        plot(df, breaks = breaks_i ,
             main = "Spat. Trend: f2_main")
      }
      if (addint) {
        readline(prompt="Press [enter] to continue")
        df <- data[, c("f12_int")]
        plot(df, breaks = breaks_i ,
             main = "Spat. Trend: f12_int")
      }
    }
  }
}
