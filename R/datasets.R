#' Regional unemployment rates Italian provinces
#'
#' A panel dataset containing unemployment rates and other economic
#' variables for Italian NUTS-3 provinces during the years 1996-2014.
#'
#' @format A data frame with 1957 rows and 17 variables:
#' \describe{
#'   \item{prov}{province (NUTS-3) coded as a name.}
#'   \item{reg}{region (NUTS-2) coded as a name.}
#'   \item{area}{area of the province (km~2~).}
#'   \item{year}{year.}
#'   \item{unrate}{unemployment rate (percentage).}
#'   \item{agri}{share of employment in agriculture (percentage).}
#'   \item{ind}{share of employment in industry (percentage).}
#'   \item{cons}{share of employment in construction (percentage).}
#'   \item{serv}{share of employment in services (percentage).}
#'   \item{popdens}{population density.}
#'   \item{empgrowth}{employment growth rate (percentage).}
#'   \item{partrate}{labor force participation rate, i.e. the
#'                   ratio between the total labor force and the
#'                   working population.}
#'   \item{long}{longitude of the centroid of the province.}
#'   \item{lat}{latitude of the centroid of the province.}
#'   \item{South}{dummy variable with unit value for southern provinces.}
#' }
#'
#' @source Italian National Institute of Statistics (ISTAT)
#'         \url{https://www.istat.it/}
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
#'         \url{https://www.istat.it/}
"Wsp_it" 
