#' @docType package
#' @name pspatreg-package
#' @rdname pspatreg-package
#'
#' @title Spatio-Temporal Semiparametric Regression Models with Spatial Lags.
#'
#' @description
#'  \pkg{pspatreg} offers the user a collection of functions to estimate 
#'  and make inference of geoadditive spatial or spatio-temporal 
#'  semiparametric regression models of type \emph{ps-sim}, \emph{ps-sar}, 
#'  \emph{ps-sem}, \emph{ps-sarar}, \emph{ps-sdm}, \emph{ps-sdem} or 
#'  \emph{ps-slx}. These type of specifications are very
#'   general and they can include parametric and non-parametric 
#'   covariates, spatial or spatio-temporal non-parametric
#'   trends and spatial lags of the dependent and independent variables and/or 
#'   the noise of the model. The non-parametric terms (either trends or 
#'   covariates) are modeled using P-Splines. The non-parametric 
#'   trend can be decomposed in an ANOVA way including main and 
#'   interactions effects of 2nd and 3rd order. The estimation 
#'   method can be restricted maximum likelihood (REML)
#'   or maximum likelihood (ML).
#'
#' @details
#'  Some functionalities that have been included in \pkg{pspatreg} package are:
#'
#' @section 1. Estimation of the geoadditive spatial or spatio-temporal semiparametric regression model:
#'  \pkg{pspatreg} allows the estimation of geoadditive spatial or spatio-temporal 
#'  semiparametric regression models which could include:
#'  \itemize{
#'    \item An spatial or spatio-temporal trend, that is,
#'      a geoadditive model either for cross-section data or
#'      for panel data. This trend can be decomposed in main 
#'      and interaction functions in an ANOVA way. The spatial 
#'      (or spatio-temporal) trend gather the potential spatial 
#'      heterogeneity of the data. 
#'    \item Parametric covariates as usual in regression models.
#'    \item Non-parametric covariates in which the functional relationship 
#'      is estimated from the data. Both the trends and non-parametric
#'      covariates are modelled using P-splines.
#'    \item Spatial dependence adding spatial lags of the dependent
#'      and independent variables as usual in spatial econometric models.
#'      These models gather the potential spatial spillovers.
#'   }
#'   Once specified, the whole model can be estimated using either 
#'   restricted maximum-likelihood (REML) or maximum likelihood (ML).
#'   The spatial econometric specifications allowed in \pkg{pspatreg}
#'   are the following ones:  
#'  \itemize{
#'    \item \emph{ps-sim}: geoadditive semiparametric model without 
#'      spatial effects (in addition to the spatial or spatio-temporal trend, 
#'      if it is included).
#'      \deqn{ y = f(s_1,s_2,\tau_{t}) y 
#'                 + X \beta +
#'                 + \sum_{i=1}^k g(z_i) 
#'                 + \epsilon }
#'      where:
#'      \itemize{
#'        \item \eqn{f(s_1,s_2,\tau_t)} is a smooth spatio-temporal trend
#'          of the spatial coordinates \eqn{s1,s_2} and of the temporal
#'          coordinates \eqn{\tau_t}.
#'        \item \eqn{X} is a matrix including values of parametric covariates.
#'        \item \eqn{g(z_i)} are non-parametric smooth functions of the
#'          covariates \eqn{z_i}.
#'        \item \eqn{\epsilon ~ N(0,R)} where \eqn{ R = \sigma^2 I_T} if errors
#'          are uncorrelated or it follows an AR(1) temporal autoregressive 
#'          structure for serially correlated errors.
#'      }   
#'    \item \emph{ps-slx}: geoadditive semiparametric model with 
#'      spatial lags of the regresors (either parametric or non-parametric):
#'      \deqn{ y =  f(s_1,s_2,\tau_{t}) 
#'                + X \beta 
#'                + (W_{N} \otimes I_T) X \theta  
#'                + \sum_{i =1}^k g(z_i) 
#'                + \sum_{i = 1}^k g((\gamma_i*W_{N} \otimes I_T) z_i)  
#'                + \epsilon }
#'      where:
#'      \itemize{
#'        \item \eqn{W_N} is the spatial weights matrix.
#'        \item \eqn{I_T} is an identity matrix of order 
#'          \eqn{T} (\emph{T = 1} for pure spatial data).
#'      }
#'    \item \emph{ps-sar}: geoadditive semiparametric model with spatial 
#'      lag of the dependent variable
#'      \deqn{ y = (\rho*W_{N} \otimes I_T) y 
#'                + f(s_1,s_2,\tau_{t}) 
#'                + X \beta 
#'                + \sum_{i =1}^k g(z_i) 
#'                + \epsilon }
#'              
#'    \item \emph{ps-sem}: geoadditive semiparametric model 
#'      with a spatial lag of the noise of the model
#'      \deqn{ y = f(s_1,s_2,\tau_{t}) 
#'                 + X \beta 
#'                 + \sum_{i =1}^k g(z_i) 
#'                 + u }
#'      \deqn{ u = (\delta*W_{N} \otimes I_T) u 
#'                  + \epsilon }
#'    \item \emph{ps-sdm}: geoadditive semiparametric model 
#'      with spatial lags of the endogenous variable and 
#'      of the regressors (spatial durbin model)
#'      \deqn{ y = (\rho*W_{N} \otimes I_T) y 
#'                + f(s_1,s_2,\tau_{t}) 
#'                + X \beta 
#'                + (W_{N} \otimes I_T) X \theta  
#'                + \sum_{i = 1}^k g(z_i) 
#'                + \sum_{i = 1}^k g((\gamma_i*W_{N} \otimes I_T) z_i)  
#'                + \epsilon }
#'    \item \emph{ps-sdem}: geoadditive semiparametric model
#'      with spatial errors and spatial lags of 
#'      the endogenous variable and of the regressors
#'      \deqn{ y = f(s_1,s_2,\tau_{t}) 
#'                 + X \beta 
#'                 + (W_{N} \otimes I_T) X \theta  
#'                 + \sum_{i = 1}^k g(z_i) 
#'                 + \sum_{i = 1}^k g((\gamma_i*W_{N} \otimes I_T) z_i)  
#'                 + u }
#'      \deqn{ u = (\delta*W_{N} \otimes I_T) u 
#'                  + \epsilon }
#'    \item \emph{ps-sarar}: geoadditive semiparametric model 
#'      with a spatial lag for: both dependent variable 
#'      and errors
#'      \deqn{ y = (\rho*W_{N} \otimes I_T) y 
#'                + f(s_1,s_2,\tau_{t}) 
#'                + X \beta 
#'                + (W_{N} \otimes I_T) X \theta  
#'                + \sum_{i = 1}^k g(z_i) 
#'                + \sum_{i = 1}^k g((\gamma_i*W_{N} \otimes I_T) z_i)  
#'                + u }
#'      \deqn{ u = (\delta*W_{N} \otimes I_T) u 
#'                  + \epsilon }
#'    } 
#' @section 2. Plot of the spatial and spatio-temporal trends:
#' 
#' @section 3. Impacts spillovers:
#' 
#' @section 4. Additional methods:
#'
#' @section Datasets:
#'   \pkg{pspatreg} includes a spatio-temporal panel database including 
#'   observations of unemployment,  economic variables and spatial coordinates
#'   of 103 Italian provinces between 1996-2019 (yearly data). 
#'   This database is provided as R data (Rmd format) and can be 
#'   load using the command \code{data(unemp_it, package = "pspatreg")}.
#'   The database also includes a \emph{W} spatial neighborhood matrix
#'   of the Italian provinces (based on queen criterium). Furthermore, 
#'   the shapefile is also included and can be used to plot spatial and 
#'   spatio-temporal trends estimated for each province. A lot of examples 
#'   are included in the help of the funcions (see, for example, 
#'   \code{?pspatfit}).
#'   For the spatial pure case (2d) the examples use the database Ames 
#'   included in \pkg{AmesHousing}. 
#'   See the help of \code{?AmesHousing::make_ames} for the variables included 
#'   in this database. Examples of hedonic models including geoadditive 
#'   spatial econometric regressions are included in the examples 
#'   of \pkg{pspatreg} package. 
#'   CONTINUAR AQUÍ... 
#'   
#'   illustrate the capabilities of different functions. Briefly, their 
#'   main characteristics are the following \cr
#'   \itemize{
#'     \item The \emph{spc} dataset (Spatial Phillips-Curve) is a 
#'      classical dataset taken from Anselin (1988, p. 203), of small 
#'      dimensions.
#'      \item The \emph{NCOVR} dataset comprises Homicides and a list of 
#'      selected socio-economic variables for continental U.S. counties 
#'      in four decennial census years: 1960, 1970, 1980 and 1990. 
#'      It is freely available from
#'        \url{https://geodacenter.github.io/data-and-lab/ncovr/}. 
#'        \emph{NCOVR} is a typical spatial panel dataset \emph{(G=1)}.
#'      \item The \emph{spain.covid} dataset comprises Within and Exit mobility index 
#'      together with the weeklly incidence COVID-19 at Spain provinces from 
#'      February 21 to May 21 2020. 
#'      \url{https://www.mitma.gob.es/ministerio/covid-19/evolucion-movilidad-big-data}
#'    }
#'
#' @references 
#'   \itemize{ 
#'     \item Basile, R.; Durban, M.; Minguez, R.; Montero, J. M.; and 
#'     Mur, J. (2014). Modeling regional economic dynamics: Spatial
#'     dependence, spatial heterogeneity and nonlinearities. 
#'     \emph{Journal of Economic Dynamics and Control}, (48), 229-245.
#'     <doi: 10.1016/j.jedc.2014.06.011>
#'
#'   \item Eilers, P. and Marx, B. (1996). Flexible Smoothing with 
#'     B-Splines and Penalties. \emph{Statistical Science}, (11), 89-121.
#'     
#'   \item Fahrmeir, L.; Kneib, T.;  Lang, S.; and Marx, B. (2013). 
#'     \emph{Regression. Models, Methods and Applications}.
#'      Springer.
#'     
#'   \item Lee, D. and Durban, M. (2011). P-Spline ANOVA Type Interaction 
#'     Models for Spatio-Temporal Smoothing. \emph{Statistical Modelling}, 
#'     (11), 49-69. <doi: 10.1177/1471082X1001100104>
#'
#'   \item Lee, D. J., Durban, M., and Eilers, P. (2013). Efficient
#'     two-dimensional smoothing with P-spline ANOVA mixed models 
#'     and nested bases. \emph{Computational Statistics & Data Analysis}, 
#'     (61), 22-37. <doi: 10.1016/j.csda.2012.11.013>
#'
#'   \item LeSage, J. and Pace, K. (2009). \emph{Introduction to 
#'     Spatial Econometrics}. CRC Press, Boca Raton.
#'
#'   \item Minguez, R.; Basile, R. and Durban, M. (2020). An Alternative 
#'     Semiparametric Model for Spatial Panel Data. \emph{Statistical Methods and Applications},
#'     (29), 669-708. <doi:	10.1007/s10260-019-00492-8>
#'
#'   \item Montero, J., Minguez, R., and Durban, M. (2012). SAR models 
#'     with nonparametric spatial trends: A P-Spline approach. 
#'     \emph{Estadistica Espanola}, (54:177), 89-111.
#'
#'   \item Rodriguez-Alvarez, M. X.; Kneib, T.; Durban, M.; Lee, D.J.
#'     and Eilers, P. (2015). Fast smoothing parameter separation 
#'     in multidimensional generalized P-splines: the SAP algorithm.
#'     \emph{Statistics and Computing} 25 (5), 941-957. 
#'     <doi: 10.1007/s11222-014-9464-2>
#' }
#'
#' @author
#'   \tabular{ll}{
#'   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#'     Roberto Basile \tab \email{roberto.basile@@univaq.it} \cr
#'     Maria Durban \tab \email{mdurban@@est-econ.uc3m.es} \cr
#'     Gonzalo Espana-Heredia \tab \email{gehllanza@@gmail.com} \cr
#'   }
#'
#' @importFrom AmesHousing make_ames
#' @importFrom akima interp
#' @importFrom dplyr left_join
#' @importFrom fields image.plot
#' @importFrom ggplot2 ggplot geom_line ggtitle labs aes xlim ylim
#' @importFrom graphics image contour matplot title points
#' @importFrom graphics par abline lines
#' @importFrom grDevices heat.colors 
#' @importFrom MASS ginv mvrnorm
#' @importFrom Matrix bandSparse bdiag crossprod determinant  
#' @importFrom Matrix diag Diagonal kronecker Matrix 
#' @importFrom Matrix rowSums solve t tcrossprod 
#' @importFrom methods as 
#' @importFrom minqa bobyqa
#' @importFrom numDeriv hessian
#' @importFrom plm plm pdata.frame Within
#' @importFrom sf st_as_sf st_drop_geometry st_coordinates 
#' @importFrom sf st_geometry_type 
#' @importFrom spatialreg get.ZeroPolicyOption create_WX   
#' @importFrom spatialreg can.be.simmed jacobianSetup do_ldet 
#' @importFrom spatialreg intImpacts lmSLX invIrW trW
#' @importFrom spdep listw2mat mat2listw nb2listw    
#' @importFrom spdep tri2nb graph2nb soi.graph is.symmetric.nb  
#' @importFrom splines spline.des  
#' @importFrom stats var sd loess predict vcov 
#' @importFrom stats model.response as.formula .getXlevels
#' @importFrom stats pchisq pnorm pt rnorm qnorm
#' @importFrom stats coefficients fitted residuals printCoefmat
#' @importFrom stats model.frame model.matrix terms
#' @importFrom stats anova coef formula logLik AIC BIC
#' @importFrom stats na.action napredict update
#' @importFrom stringr str_detect str_replace str_split
NULL
