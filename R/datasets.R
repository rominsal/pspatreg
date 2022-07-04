#' Regional unemployment rates Italian provinces
#'
#' A panel dataset containing unemployment rates and other economic
#' variables for Italian NUTS-3 provinces during the years 1996-2019.
#'
#' @format A data frame with 2472 rows and 17 variables:
#' \describe{
#'   \item{prov}{province (NUTS-3) coded as a number.}
#'   \item{name}{province (NUTS-3) coded as a name.}
#'   \item{reg}{region (NUTS-2) coded as a name.}
#'   \item{year}{year.}
#'   \item{area}{area of the province (km~2~).}
#'   \item{unrate}{unemployment rate (percentage).}
#'   \item{agri}{share of employment in agriculture (percentage).}
#'   \item{ind}{share of employment in industry (percentage).}
#'   \item{cons}{share of employment in construction (percentage).}
#'   \item{serv}{share of employment in services (percentage).}
#'   \item{popdens}{population density.}
#'   \item{partrate}{labor force participation rate, i.e. the
#'                   ratio between the total labor force and the
#'                   working population.}
#'   \item{empgrowth}{employment growth rate (percentage).}
#'   \item{long}{longitude of the centroid of the province.}
#'   \item{lat}{latitude of the centroid of the province.}
#'   \item{South}{dummy variable with unit value for southern provinces.}
#'   \item{ln_popdens}{logarithm of population density.}
#' }
#'
#' @source Italian National Institute of Statistics (ISTAT)
#'         \emph{https://www.istat.it}
"unemp_it"

#' Spatial weight matrix for Italian provinces
#' 
#' A spatial weight matrix row-standardized for Italian NUTS-3 provinces
#' 
#' @format A row-standardized squared matrix with 103 rows and columns. 
#' The rows and columns follow the same order than provinces included in 
#' \emph{unemp_it} data frame.
#' 
#' @source Italian National Institute of Statistics (ISTAT)
#'         \emph{https://www.istat.it}
"Wsp_it"

#' map of Italian provinces
#' 
#' An sf object including a map of Italian NUTS-3 provinces
#' 
#' @format An sf object with 103 rows and 2 columns:
#' \describe{
#'   \item{COD_PRO}{province (NUTS-3) coded as a number.}
#'   \item{geometry}{geometry (polygons) of the sf object.}
#'  }
#' @source Italian National Institute of Statistics (ISTAT)
#'         \emph{https://www.istat.it}
"map_it" 

#' Productivity growth and internal net migration - Italian provinces
#'
#' A spatial dataframe including a map of Italian NUTS-3 provinces and 
#' cross-sectional dataset on provincial labor productivity growth rates,
#' internal net migration rates, and other economic variables.
#'
#' @format An sf object with 107 rows and 9 columns:
#' \describe{
#'   \item{COD_PROV}{province (NUTS-3) coded as a number.}
#'   \item{DEN_PROV}{province (NUTS-3) coded as a name.}
#'   \item{longitude}{longitude of the centroid of the province.}
#'   \item{latitude}{latitude of the centroid of the province.}
#'   \item{lnPROD_0}{log of labor productivity in 2002
#'                   (measured as gross value added per worker).}
#'   \item{growth_PROD}{Average annual growth rate of labor productivity 
#'                      over the period 2002-2018.}
#'   \item{lnoccgr}{Average annual growth rate of employment over the period 2002-2018.}
#'   \item{net}{Average annual provincial internal net migration rate 
#'              (computed as the difference between internal immigration and 
#'              emigration flows of the working-age population, i.e. 
#'              people aged 15-65, divided by the total population).}
#'   \item{geometry}{geometry (polygons) of the sf object.}
#' }
#'
#' @source Italian National Institute of Statistics (ISTAT)
#'         \emph{https://www.istat.it}
"prod_it"

#' Spatial weight matrix for Italian provinces
#' 
#' A spatial weight matrix row-standardized for Italian NUTS-3 provinces
#' 
#' @format A row-standardized squared matrix with 107 rows and columns. 
#' The rows and columns follow the same order than provinces included in 
#' \emph{unemp_it} data frame.
#' 
#' @source Italian National Institute of Statistics (ISTAT)
#'         \emph{https://www.istat.it}
"lwsp_it"
