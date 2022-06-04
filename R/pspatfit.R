#' @name pspatfit
#' @rdname pspatfit
#'
#' @title Estimate spatial or spatio-temporal semiparametric 
#'   regression models from a spatial econometric perspective.
#'
#' @description Estimate geoadditive spatial or spatio-temporal 
#'   semiparametric regression models of type \emph{ps-sar}, \emph{ps-sem}, \emph{ps-sarar}, 
#'   \emph{ps-sdm}, \emph{ps-sdem} or
#'   \emph{ps-slx}. These type of specifications are very
#'   general and they can include parametric and non-parametric 
#'   covariates, spatial or spatio-temporal non-parametric
#'   trends and spatial lags of the dependent and independent variables and/or 
#'   the noise of the model. 
#'   The non-parametric terms (either trends or covariates) 
#'   are modeled using P-Splines. 
#'   The non-parametric trend can be decomposed in an ANOVA way 
#'   including main and interactions effects of 2nd and 3rd order. 
#'   The estimation method can be restricted maximum likelihood (REML)
#'   or maximum likelihood (ML).
#'
#' @param formula A formula similar to GAM specification including
#'   parametric and non-parametric terms. Parametric covariates
#'   are included in the usual way and non-parametric P-spline smooth terms are
#'   specified using \code{pspl(.)} and \code{pspt(.)} for the non-parametric 
#'   covariates and spatial   or spatio-temporal trend, respectively.
#'   More details in \emph{Details} and \emph{Examples}.
#' @param data  A data frame containing the parametric and non-parametric
#'   covariates included in the model. Also, if a \code{pspt(.)} term is 
#'   included in formula, the data frame must include the spatial and temporal
#'   coordinates specified in \code{pspt(.)}. In this case coordinates
#'   must be ordered choosing time as fast index and spatial coordinates
#'   as low indexes. See \code{head(unemp_it)} as an example.
#' @param na.action	A function (default \code{options("na.action")}),
#'  can also be `na.omit` or `na.exclude` with consequences 
#'   for residuals and fitted values. It may be necessary to set 
#'   `zero.policy` to `TRUE` because this subsetting may 
#'   create no-neighbour observations.    
#' @param listw Default = `NULL` This will create a model with no spatial dependency.
#'   To include spatial dependency, \code{listw} should be a spatial neighbours list object 
#'   created for example by \code{\link[spdep]{nb2listw}} from \pkg{spdep}
#'   package; if \code{\link[spdep]{nb2listw}} not given, set to 
#'   the same spatial weights as the \code{listw} argument. It can
#'   also be a spatial weighting matrix of order \emph{(NxN)} instead of
#'   a \code{listw} neighbours list object. 
#' @param method Similar to the corresponding parameter of 
#'   \code{\link[spatialreg]{lagsarlm}} function in \pkg{spatialreg} package. 
#'   "eigen" (default) - the Jacobian is computed as the product of 
#'   (1 - rho*eigenvalue) using \code{\link[spatialreg]{eigenw}} from package
#'   spatialreg. For big samples (> 500) method = "eigen" is not recommended.
#'   Use "spam" or "Matrix_J" for strictly symmetric weights lists of 
#'   styles "B" and "C", or made symmetric by similarity 
#'   (Ord, 1975, Appendix C) if possible for styles "W" and "S", 
#'   using code from the spam or Matrix packages to calculate the 
#'   determinant; "Matrix" and "spam_update" provide updating Cholesky 
#'   decomposition methods; "LU" provides an alternative sparse matrix 
#'   decomposition approach. In addition, there are "Chebyshev" and 
#'   Monte Carlo "MC" approximate log-determinant methods; 
#'   the Smirnov/Anselin (2009) trace approximation is available 
#'   as "moments". Three methods: "SE_classic", "SE_whichMin", 
#'   and "SE_interp" are provided experimentally, the first to 
#'   attempt to emulate the behaviour of Spatial Econometrics 
#'   toolbox ML fitting functions. All use grids of log determinant 
#'   values, and the latter two attempt to ameliorate some features 
#'   of "SE_classic".
#' @param interval Search interval for autoregressive parameter.
#'   Default = `NULL`.
#' @param trs Similar to the corresponding parameter of 
#'   \code{\link[spatialreg]{lagsarlm}} function in \pkg{spatialreg} package.
#'   Default `NULL`, if given, a vector of powered spatial weights 
#'   matrix traces output by \code{\link[spdep]{trW}}.
#' @param zero.policy Similar to the corresponding parameter of 
#'   \code{\link[spatialreg]{lagsarlm}} function in \pkg{spatialreg} package. 
#'   If `TRUE` assign zero to the lagged value of zones without 
#'   neighbours, if `FALSE` assign `NA` - causing 
#'   \code{pspatfit()} to terminate with an error. Default = `NULL`. 
#' @param type Type of spatial model specification following 
#'   the usual spatial econometric terminology. 
#'   Default = \code{"sim"} this creates
#'   a model with no type of spatial dependency. 
#'   Types of spatial models available (similar to 
#'   \pkg{spsur} package):
#'   \code{"sar"}, \code{"sem"}, \code{"sdm"}, 
#'   \code{"sdem"}, \code{"sarar"}, or \code{"slx"}. 
#'   When creating a \code{"slx"}, \code{"sdem"}
#'   or \code{"sdm"} model, it is necessary to include the formula of the Durbin part in 
#'   the \code{Durbin} argument in the same way than 
#'   \pkg{spsur} or \pkg{spatialreg} packages.
#'   There are examples on how to create these models 
#'   in \emph{Examples} section.
#' @param Durbin Default = `NULL`. 
#'   If model is of \code{type = "sdm"}, \code{"sdem"} or 
#'   \code{"slx"} then this argument should be a formula
#'    of the subset of explanatory variables to be 
#'    spatially lagged in the right hand side part of 
#'    the model. See \code{\link[spsur]{spsurml}} for 
#'    a similar argument.
#' @param cor Type of temporal correlation for temporal data. Possible values 
#'   are \code{"none"} (default) or \code{"ar1"}.
#' @
#' param demean Logical value to include a demeaning 
#'   for panel data. Default = `FALSE`. 
#'   The demeaning is done previously to the estimation for
#'   both parametric and nonparametric terms. It is not possible
#'   to set \code{demean = TRUE} when spatio-temporal trends
#'   are included. 
#' @param eff_demean Type of demeaning for panel data.
#'   Possible values are \code{"individual"} (default),
#'   \code{"time"} and \code{"twoways"}.
#' @param index Vector of variables indexing panel data. 
#'   First variable corresponds to individuals and second 
#'   variable corresponds to temporal coordinate (fast index). 
#'   It follows the same rules than \code{\link[plm]{plm}} function
#'   in package \pkg{plm}.        
#' @param dynamic Logical value to set a dynamic model.
#'   Dynamic models include a temporal lag of the dependent
#'   variable in the right-hand side of the equation.
#'   Default = `FALSE`. 
#' @param control List of extra control arguments. See 
#'   \emph{Control Arguments} section below.
#'
#' @details
#'   Function to estimate the model:
#'   \deqn{ y = (\rho*W_{N} \otimes I_T) y 
#'     + f(s_1,s_2,\tau_{t}) 
#'     + X \beta 
#'     + (W_{N} \otimes I_T) X \theta  
#'     + \sum_{i = 1}^k g(z_i) 
#'     + \sum_{i = 1}^k g((\gamma_i*W_{N} \otimes I_T) z_i)  
#'     + \epsilon }
#'  where:
#'   \itemize{
#'     \item \eqn{f(s_1,s_2,\tau_t)} is a smooth spatio-temporal trend
#'     of the spatial coordinates \eqn{s1,s_2} and of the temporal
#'     coordinates \eqn{\tau_t}.
#'     \item \eqn{X} is a matrix including values of parametric covariates.
#'     \item \eqn{g(z_i)} are non-parametric smooth functions of the
#'     covariates \eqn{z_i}.
#'     \item \eqn{W_N} is the spatial weights matrix.
#'     \item \eqn{\rho} is the spatial spillover parameter.
#'     \item \eqn{I_T} is an identity matrix of order \eqn{T} (\emph{T=1}
#'       for pure spatial data).
#'     \item \eqn{\epsilon ~ N(0,R)} where \eqn{R = \sigma^2 I_T} if errors
#'       are uncorrelated or it follows an AR(1) temporal autoregressive 
#'       structure for serially correlated errors.
#'   }
#' \describe{
#'   \item{Including non-parametric terms}{
#'   The non-parametric terms are included in \code{formula} using
#'   \code{pspt(.)} for spatial or spatio-temporal trends and 
#'   \code{pspl(.)} for other non-parametric smooth additive terms.
#'   For example, if a model includes:
#'    \itemize{
#'      \item An spatio-temporal trend with 
#'         variables \emph{long} and \emph{lat} 
#'         as spatial coordinates,and \emph{year} 
#'         as temporal coordinate.
#'      \item Two non-parametric covariates named 
#'        \emph{empgrowth} and \emph{serv}.
#'      \item  Three parametric covariates named 
#'        \emph{partrate}, \emph{agri} and \emph{cons}.
#'    }
#'    Then, the formula should be written as (choosing default values
#'    for each term): \cr
#'    
#'     \code{ unrate ~ partrate + agri + cons +
#'                    pspl(serv) + pspl(empgrowth) +
#'                    pspt(long,lat,year) } \cr
#'
#'   For a spatial trend case, the term \code{pspt(.)} does not include a 
#'   temporal coordinate, that is, in the previous example would be 
#'   specified as \code{pspt(long,lat)}.
#'   }
#'   \item{How to use \code{pspl()} and \code{pspt()}}{   
#'   Note that both in \code{pspl(.)} and \code{pspt(.)}, we have 
#'   to include the number of knots, named \code{nknots}, 
#'   which is the dimension of the basis used to represent the smooth term.
#'   The value of \code{nknots} should not be less than the dimension of the null space of the penalty
#'   for the term, see \code{\link[mgcv]{null.space.dimension}} 
#'   and \code{\link[mgcv]{choose.k}} from \pkg{mgcv} 
#'   package to know how to choose \code{nknots}.
#'   
#'   In \code{pspl(.)} the default is \code{nknots = 10}, see the help of \code{\link{pspl}} function.
#'   In this term we can only include single variables, so if we want more than one 
#'   non-parametric variable we will use a \code{pspl(.)} term for each nonparametric variable.
#'   
#'   On the other hand, \code{pspt(.)} is used for spatial smoothing 
#'   (when temporal coordinate is `NULL`)  or 
#'   spatio-temporal smoothing (when a variable is provided for
#'   the temporal coordinate). 
#'   The default for the temporal coordinate is \code{time = NULL}, 
#'   see the help of \code{\link{pspt}}, and the default number of knots
#'   are \code{nknots = c(10, 10, 5)}. If only 
#'   include spatial smoothing, \code{nknots} will be a length 2 vector 
#'   indicating the basis for each spatial coordinate. 
#'   For spatio-temporal smoothing, it will be a length 3 vector.
#'   }
#'   \item{ANOVA descomposition}{
#'   In many situations the  spatio-temporal trend, given by 
#'   \eqn{f(s_1,s_2,\tau_t)}, can be very complex and the use of a 
#'   multidimensional smooth function may not be flexible enough to 
#'   capture the structure in the data. Furthermore, the estimation of 
#'   this trend can become computationally intensive especially for 
#'   large databases.\cr
#'   To solve this problem, Lee and Durban (2011) proposed an ANOVA-type
#'   decomposition of this spatio-temporal trend where spatial and 
#'   temporal main effects, and second- and third-order interaction 
#'   effects can be identified as:
#'   
#'   \deqn{ f(s_1, s_2, \tau_t) = f_1(s_1) + f_2(s_2) + f_t(\tau_t) +
#'          f_{1,2}(s_1, s_2) +  f_{1,t}(s_1, \tau_t) +
#'          f_{2,t}(s_2, \tau_t) + f_{1,2,t}(s_1, s_2, \tau_t) }
#'   
#'   In this equation the decomposition of the spatio-temporal trend 
#'   is as follows:
#'   \itemize{
#'     \item Main effects given by the functions 
#'       \eqn{f_1(s_1), f_2(s_2)} and \eqn{f_t(\tau_t)}.
#'     \item Second-order interaction effects given by the functions 
#'       \eqn{f_{1,2}(s_1,s_2), f_{1,t}(s_1,\tau_t)} 
#'       and \eqn{f_{2,t}(s_2,\tau_t)}.
#'     \item Third-order interaction effect given by the 
#'       function \eqn{f_{1,2,t}(s_1,s_2,\tau_t)}.    
#'   }
#'   
#'   In this case, each effect can have its own 
#'   degree of smoothing allowing a greater flexibility for the 
#'   spatio-temporal trend. The ANOVA decomposition of the trend
#'   can be set as an argument in \code{pspt(.)} terms choosing
#'   \code{psanova = TRUE}.
#'   
#'   For example to choose an ANOVA decomposition in the 
#'   previous case we can set: \cr
#'   
#'    \code{pspt(long, lat, year, nknots = c(18,18,8),
#'   psanova = TRUE)}
#'   \cr
#'   
#'   In most empirical cases main effects functions are more flexible than 
#'   interaction effects functions and therefore, the number of knots in B-Spline 
#'   bases for interaction effects do not need to be as big as the 
#'   number of knots for main effects. 
#'   \emph{Lee et al.}, (2013) proposed a nested basis procedure
#'   in which the number of knots for the interaction effects functions are 
#'   reduced using \emph{divisors} such that the space spanned by 
#'   B-spline bases used for interaction effects are a subset of the 
#'   space spanned by B-spline bases used for main effects. 
#'   The \emph{divisors} can be specified as an argument in 
#'   \code{pspt(.)} terms. \cr 
#'   To do this, there are three
#'   arguments available inside \code{pspt()} to define the divisors.
#'   These arguments are named \code{nest_sp1}, \code{nest_sp2} and 
#'   \code{nest_time}, respectively. 
#'   The value for these arguments
#'   are vector parameters including divisors of
#'   the \code{nknots} values. \cr
#'   
#'   For example, if we set \code{nest_sp1 =
#'   c(1,2,2)} between the arguments of \code{pspl(.)}, 
#'   we will have all knots for main effect of \emph{s_1}, 
#'   \emph{18/2=9} knots for each second-order effect including \emph{s_1}, 
#'   and \emph{8/2=4} knots for the third order effect including \emph{s_1}. It
#'   is important that the vector of numbers will be integer divisors 
#'   of the values in \code{nknots}.
#'   See section \emph{Examples} for more details.
#'   
#'   Eventually, any effect function can be excluded of the ps-anova
#'   spatio-temporal trend. To exclude main effects, the arguments  
#'   \code{f1_main}, \code{f2_main} or \code{ft_main} have to be set to
#'   `FALSE` (default=`TRUE`).
#'   We can also exclude the second- and third-order
#'   effects functions setting to `FALSE` the arguments \code{f12_int}, \code{f1t_int}, 
#'   \code{f2t_int} or \code{f12t_int} in \code{pspl(.)}.
#'   }
#'  }
#'  All the terms included in the model are jointly fitted using Separation of Anisotropic 
#'   Penalties (SAP) algorithm (see \emph{Rodriguez-Alvarez et al., (2015)}) 
#'   which allows to the mixed model reparameterization of the model. 
#'   For type of models \code{"sar", "sem", "sdm", "sdem", "sarar"}  or 
#'   \code{cor = "ar1"}, the parameters \eqn{\rho}, \eqn{\lambda} and \eqn{\phi} 
#'   are numerically estimated using 
#'   \code{\link[minqa]{bobyqa}} function implemented in package \pkg{minqa}.
#'   In these cases, an iterative process between SAP and numerical 
#'   optimization of \eqn{\rho}, \eqn{\lambda} and \eqn{\phi} is applied until
#'   convergence. See details in \emph{Minguez et al.}, (2018).
#'  \describe{   
#'  \item{Plotting non-parametric terms}{
#'    To plot the non-linear functions corresponding to 
#'    non-parametric terms we need to compute the fitted values,
#'   and standard erros, using \code{fit_terms()} function 
#'   and, afterwards, use \code{plot_terms()} function to 
#'   plot the non-linear functions. \cr 
#'   An example of how plot the functions of non-parametric 
#'   terms given by \code{"var1"} and \code{"var2"} variables is given by 
#'   the next lines of code (it is assumed that a previous 
#'   model has been fitted using \code{pspatfit(.)} and 
#'   saved as an object named \code{model}): \cr
#'   
#'   \code{list_varnopar <- c("var1", "var2")} \cr
#'   \code{terms_nopar <- fit_terms(model, list_varnopar)} \cr
#'   \code{plot_terms(terms_nopar, data)} 
#'   \cr
#'   
#'   The \code{data} argument of \code{plot_terms()} usually 
#'   corresponds to the dataframe used to fitted the model 
#'   although a different database can be used to plot the 
#'   non-parametric terms.}
#'  \item{Spatial impacts}{
#'    For the spatial models given by  \code{type = "sar"}, 
#'    \code{"sdm"}, \code{"sdem"}, \code{"sarar"} 
#'    or \code{"slx"} it is possible to compute spatial 
#'    spillovers as usual in spatial econometric specifications. 
#'    Nevertheless, in this case we need to distinguish between 
#'    parametric and non-parametric covariates when computing spatial 
#'    impacts.}
#'    \itemize{
#'      \item spatial impacts for parametric covariates \cr
#'       In this case, the spatial impacts are computed in the 
#'       usual way using simulation. See LeSage and Page (2009) 
#'       for computational details. The function \code{impactspar()}
#'       computes the direct, indirect and total impacts for 
#'       parametric covariates and return and object similar to 
#'       the case of \pkg{spatialreg} and \pkg{spsur} packages.
#'       The inference for \code{"sar"}, \code{"sdm"}, 
#'       and \code{"sarar"} types is based on simulations 
#'       and for \code{"slx"} and \code{"sdem"} types the 
#'       standard errors or total impacts are computed using 
#'       the variance-covariance matrix of the fitted model. 
#'       The \code{summary()} method can be used to present the 
#'       the complete table of spatial impacts in this parametric case.
#'       See the help of \code{\link{impactspar}} to know the 
#'       additional arguments of the function. A little example 
#'       is given in the next lines of code:\cr
#'       
#'       \code{imp_parvar <- impactspar(MODEL, listw = W)} \cr
#'       \code{summary(imp_parvar)}
#'   \item spatial impacts for non-parametric covariates \cr
#'     In this case direct, indirect and total 
#'     \emph{spatial impacts functions} are 
#'     obtained using \code{impactsnopar}. The details of 
#'     computation and inference can be obtained from the help 
#'     of \code{\link{impactsnopar}}. 
#'     The argument \code{viewplot} of \code{impactsnopar} 
#'     have to be set as `TRUE` to plot the spatial impacts 
#'     functions. Another way to get the same plots is using 
#'     \code{plot_impactsnopar} function with the output 
#'     of \code{impactsnopar}. 
#'     Next lines give an example of both cases: \cr
#'       
#'     \code{imp_nparvar <- impactsnopar(MODEL, listw = W, viewplot = TRUE)} \cr
#'     \code{imp_nparvar <- impactsnopar(MODEL, listw = W, viewplot = FALSE)} \cr
#'     \code{plot_impactsnopar(imp_nparvar, data = DATA)} \cr
#'  }
#' }  
#'
#' @return A list object of class \emph{pspatreg}
#' \tabular{ll}{
#'  \code{call} \tab Matched call. \cr
#'  \code{terms} \tab The terms object used. \cr
#'  \code{contrasts} \tab (only where relevant) the contrasts used
#'              for parametric covariates. \cr
#'  \code{xlevels} \tab (only where relevant) a record of the levels
#'              of the parametric factors used in fitting. \cr
#'  \code{data} \tab dataframe used as database. \cr
#'  \code{nsp} \tab number of spatial observations. \cr
#'  \code{nt} \tab number of temporal observations. It is set
#'    to \code{nt=1} for spatial data. \cr
#'  \code{nfull} \tab total number of observations. \cr
#'  \code{edftot} \tab Equivalent degrees of freedom for the whole model. \cr
#'  \code{edfspt} \tab Equivalent degrees of freedom for smooth
#'              spatio-temporal or spatial trend. \cr
#'  \code{edfnopar} \tab Equivalent degrees of freedom for
#'              non-parametric covariates. \cr
#'  \code{psanova} \tab \emph{TRUE} if spatio-temporal or spatial trend is 
#'    PS-ANOVA. \cr
#'  \code{type} \tab Value of \code{type} argument in the call to \code{pspatfit}. \cr
#'  \code{listw} \tab Value of \code{listw} argument in the call to \code{pspatfit}. \cr
#'  \code{Durbin} \tab Value of \code{Durbin} argument in the call to \code{pspatfit}. \cr
#'  \code{cor} \tab Value of \code{cor} argument in the call to \code{pspatfit}. \cr
#'  \code{dynamic} \tab Value of \code{dynamic} argument in the call to \code{pspatfit}. \cr
#'  \code{demean} \tab Value of \code{demean} argument in the call to \code{pspatfit}. \cr
#'  \code{eff_demean} \tab Value of \code{eff_demean} argument in the call to \code{pspatfit}. \cr
#'  \code{index} \tab Value of \code{index} argument in the call to \code{pspatfit}. \cr
#'  \code{bfixed} \tab Estimated betas corresponding to fixed effects in
#'              mixed model. \cr
#'  \code{se_bfixed} \tab Standard errors of fixed betas. \cr
#'  \code{brandom} \tab Estimated betas corresponding to random effects
#'              in mixed model. \cr
#'  \code{se_brandom}\tab Standard errors of random betas. \cr
#'  \code{vcov_fr} \tab Covariance matrix of fixed and random
#'              effects using frequentist or sandwich method. \cr
#'  \code{vcov_by} \tab Covariance matrix of fixed and random
#'              effects using bayesian method. \cr
#'  \code{rho} \tab Estimated rho for spatial lag of the
#'    dependent variable. \cr
#'  \code{se_rho} \tab Standard error of \eqn{rho}. \cr
#'  \code{delta} \tab Estimated delta for spatial error models. \cr
#'  \code{se_delta} \tab Standard error of \eqn{delta}. \cr
#'  \code{phi} \tab Estimated phi. If \code{cor="none"} always \eqn{phi=0}. \cr
#'  \code{se_phi} \tab Standard error of \eqn{phi}. \cr
#'  \code{fitted.values} \tab Vector of fitted values of the dependent
#'              variable. \cr
#'  \code{se_fitted.values} \tab Vector of standard errors of
#'       \code{fitted.values}. \cr
#'  \code{fitted.values_Ay} \tab Vector of fitted values of the spatial lag of
#'      dependent variable: \eqn{(\rho*W_N \otimes I_T) y}. \cr
#'  \code{se_fitted.values_Ay} \tab Vector of standard errors of
#'       \code{fitted.values_Ay}. \cr
#'  \code{residuals} \tab Vector of residuals. \cr
#'  \code{df.residual} \tab Equivalent degrees of freedom for \code{residuals}. 
#'    \cr
#'  \code{sig2}  \tab Residual variance computed as SSR/df.residual. \cr
#'  \code{llik} \tab Log-likelihood value. \cr
#'  \code{llik_reml} \tab Restricted log-likelihood value. \cr
#'  \code{aic} \tab Akaike information criterion. \cr
#'  \code{bic} \tab Bayesian information criterion. \cr
#'  \code{sp1} \tab First spatial coordinate. \cr
#'  \code{sp2} \tab Second spatial coordinate. \cr
#'  \code{time} \tab Time coordinate. \cr
#'  \code{y} \tab Dependent variable. \cr
#'  \code{X} \tab Model matrix for fixed effects. \cr
#'  \code{Z} \tab Model matrix for random effects. \cr
#' }
#'
#' @section Control Arguments:
#' \tabular{ll}{
#'   \code{optim} \tab method of estimation: \code{"llik_reml"} (default) or
#'     \code{"llik"}. \cr
#'   \code{typese} \tab method to compute 
#'     standard errors. \code{"sandwich"}  or \code{"bayesian"} (default).
#'     See Fahrmeir et al, pp. 375 for details of computations. \cr   
#'   \code{vary_init} \tab Initial value of the noise variance in the model.
#'     Default = `NULL`. \cr
#'   \code{trace} \tab A logical value set to \emph{TRUE} to show 
#'     intermediate results during the estimation process. 
#'     Default = \emph{FALSE}. \cr
#'   \code{tol1} \tab Numerical value for the tolerance of convergence
#'     of penalization parameters during the estimation process. 
#'     Default 1e-6. This tolerance is only used 
#'     for small samples (<= 500 observations). \cr
#'   \code{tol2} \tab Numerical value for the tolerance of convergence
#'     of total estimated degrees of freedom ("edftot") during the 
#'     estimation process. Default 5e-3. This tolerance is used for 
#'     medium or big samples (> 500 observations). \cr
#'   \code{tol3} \tab Numerical value for the tolerance of convergence
#'     of spatial and correlation parameters during the 
#'     estimation process. Default 1e-2.  \cr
#'   \code{maxit} \tab An integer value for the maximum number of 
#'     iterations until convergence. Default = 200. \cr
#'   \code{rho_init} \tab An initial value for \eqn{rho} parameter. 
#'     Default 0. \cr
#'   \code{delta_init} \tab An initial value for \eqn{delta} parameter. 
#'     Default 0. \cr
#'   \code{phi_init} \tab An initial value for \eqn{phi} parameter. 
#'     Default 0. \cr
#'   \code{Imult} \tab default 2; used for preparing the Cholesky 
#'       decompositions for updating in the Jacobian function \cr
#'   \code{super} \tab  if `NULL` (default), set to `FALSE` to use 
#'       a simplicial decomposition for the sparse Cholesky decomposition 
#'       and method "Matrix_J", set to `as.logical(NA)` for method "Matrix", if 
#'       `TRUE`, use a supernodal decomposition \cr
#'     \code{cheb_q} \tab default 5; highest power of the approximating 
#'       polynomial for the Chebyshev approximation \cr
#'     \code{MC_p} \tab default 16; number of random variates \cr
#'     \code{MC_m} \tab default 30; number of products of random variates 
#'       matrix and spatial weights matrix \cr
#'     \code{spamPivot} \tab  default "MMD", alternative "RCM" \cr
#'     \code{in_coef} \tab default 0.1, coefficient value for initial Cholesky 
#'       decomposition in "spam_update" \cr
#'     \code{type} \tab default "MC", used with method "moments"; alternatives 
#'       "mult" and "moments", for use if trs is missing \cr 
#'     \code{correct} \tab default `TRUE`, used with method "moments" to 
#'       compute the Smirnov/Anselin correction term \cr
#'     \code{trunc} \tab default `TRUE`, used with method "moments" to 
#'       truncate the Smirnov/Anselin correction term \cr
#'     \code{SE_method} \tab default "LU", may be "MC" \cr
#'     \code{nrho} \tab default 200, as in SE toolbox; the size of the first 
#'       stage lndet grid; it may be reduced to for example 40 \cr
#'     \code{interpn} \tab default 2000, as in SE toolbox; the size of the 
#'       second stage lndet grid \cr
#'     \code{SElndet} \tab default `NULL`, may be used to pass a 
#'       pre-computed SE toolbox style matrix of coefficients and their lndet 
#'       values to the "SE_classic" and "SE_whichMin" methods \cr
#'     \code{LU_order} \tab default `FALSE`; used in "LU_prepermutate", 
#'       note warnings given for lu method \cr
#'     \code{pre_eig} \tab default `NULL`; may be used to pass a 
#'       pre-computed vector of eigenvalues \cr
#'  } 
#'  
#' @author 
#'   \tabular{ll}{
#'     Roman Minguez \tab \email{roman.minguez@@uclm.es} \cr
#'     Roberto Basile \tab \email{roberto.basile@@univaq.it} \cr
#'     Maria Durban \tab \email{mdurban@@est-econ.uc3m.es} \cr
#'     Gonzalo Espana-Heredia \tab \email{gehllanza@@gmail.com} \cr
#'   }   
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{impactspar}} compute total, direct and indirect effect
#'                                functions for parametric continuous 
#'                                covariates.
#'   \item \code{\link{impactsnopar}} compute total, direct and indirect effect
#'                                 functions for non-parametric continuous 
#'                                 covariates.
#'   \item \code{\link{fit_terms}} compute smooth functions for non-parametric
#'                                 continuous covariates.
#'    \item\code{\link[mgcv]{gam}} well-known alternative of estimation
#'                                 of semiparametric models in \pkg{mgcv} 
#'                                 package.
#' }
#'
#' @references 
#'   \itemize{ 
#'     \item Basile, R.; Durban, M.; Minguez, R.; Montero, J. M.; and 
#'     Mur, J. (2014). Modeling regional economic dynamics: Spatial
#'     dependence, spatial heterogeneity and nonlinearities. 
#'     \emph{Journal of Economic Dynamics and Control}, (48), 229-245.
#'     <doi:10.1016/j.jedc.2014.06.011>
#'
#'   \item Eilers, P. and Marx, B. (1996). Flexible Smoothing with 
#'     B-Splines and Penalties. \emph{Statistical Science}, (11), 89-121.
#'     
#'   \item Eilers, P. and Marx, B. (2021). \emph{Practical Smoothing. 
#'   The Joys of P-Splines}. Cambridge University Press.
#'     
#'   \item Fahrmeir, L.; Kneib, T.;  Lang, S.; and Marx, B. (2021). 
#'     \emph{Regression. Models, Methods and Applications (2nd Ed.)}.
#'      Springer.
#'     
#'   \item Lee, D. and Durban, M. (2011). P-Spline ANOVA Type Interaction 
#'     Models for Spatio-Temporal Smoothing. \emph{Statistical Modelling}, 
#'     (11), 49-69. <doi:10.1177/1471082X1001100104>
#'
#'   \item Lee, D. J., Durban, M., and Eilers, P. (2013). Efficient
#'     two-dimensional smoothing with P-spline ANOVA mixed models 
#'     and nested bases. \emph{Computational Statistics & Data Analysis}, 
#'     (61), 22-37. <doi:10.1016/j.csda.2012.11.013>
#'
#'   \item LeSage, J. and Pace, K. (2009). \emph{Introduction to 
#'     Spatial Econometrics}. CRC Press, Boca Raton.
#'
#'   \item Minguez, R.; Basile, R. and Durban, M. (2020). An Alternative 
#'     Semiparametric Model for Spatial Panel Data. \emph{Statistical Methods and Applications},
#'     (29), 669-708. <doi:10.1007/s10260-019-00492-8>
#'
#'   \item Montero, J., Minguez, R., and Durban, M. (2012). SAR models 
#'     with nonparametric spatial trends: A P-Spline approach. 
#'     \emph{Estadistica Espanola}, (54:177), 89-111.
#'
#'   \item Rodriguez-Alvarez, M. X.; Kneib, T.; Durban, M.; Lee, D.J.
#'     and Eilers, P. (2015). Fast smoothing parameter separation 
#'     in multidimensional generalized P-splines: the SAP algorithm.
#'     \emph{Statistics and Computing} 25 (5), 941-957. 
#'     <doi:10.1007/s11222-014-9464-2>
#'     
#'   \item Wood, S.N. (2017). \emph{Generalized Additive Models. 
#'   An Introduction with \code{R}} (second edition). CRC Press, Boca Raton.
#' }
#'
#' @examples 
#' ################################################
#' ###################### Examples using a panel data of rate of
#' ###################### unemployment for 103 Italian provinces in 1996-2019.
#' 
#' ##########################
#' library(pspatreg)
#' ## load spatial panel and Wsp_it
#' ## 103 Italian provinces. Period 1996-2019
#' data(unemp_it, package = "pspatreg")
#' ## Wsp_it is a matrix. Create a neighboord list 
#' lwsp_it <- spdep::mat2listw(Wsp_it)
#' ## short sample for spatial pure case (2d)
#' unemp_it_short <- unemp_it[unemp_it$year == 2019, ]
#' ####  GAM pure with pspatreg
#' form1 <- unrate ~ partrate + agri + cons +
#'                  pspl(serv, nknots = 15) +
#'                  pspl(empgrowth, nknots = 20)
#' gampure <- pspatfit(form1, data = unemp_it_short)
#' summary(gampure)
#' 
#' ######################  Get Non-parametric terms of GAM with pspatreg
#' list_varnopar <- c("serv", "empgrowth")
#' terms_nopar <- fit_terms(gampure, list_varnopar)
#' ######################  Plot non-parametric terms
#' plot_terms(terms_nopar, unemp_it_short)
#' 
#' 
#' #####################  GAM + SAR Model
#' gamsar <- pspatfit(form1, data = unemp_it_short, 
#'                    type = "sar", listw = lwsp_it)
#' summary(gamsar)
#' ### Test nested models
#' anova(gampure, gamsar)
#' \donttest{
#' ######### Non-Parametric Total, Direct and Indirect impacts
#' ### with impactsnopar(viewplot = TRUE)
#' imp_nparvar <- impactsnopar(gamsar, 
#'                             listw = lwsp_it, 
#'                             viewplot = TRUE)
#' 
#' ############ Non-Parametric Total, Direct and Indirect impacts
#' ### with impactsnopar(viewplot = FALSE) and using plot_impactsnopar()
#' imp_nparvar <- impactsnopar(gamsar, listw = lwsp_it, viewplot = FALSE)
#' plot_impactsnopar(imp_nparvar, data = unemp_it_short, smooth = TRUE)
#' 
#' ###################### Parametric Total, Direct and Indirect impacts
#' imp_parvar <- impactspar(gamsar, listw = lwsp_it)
#' summary(imp_parvar)
#'
#' ###############################################
#' ### Spatial semiparametric model without spatial lag
#' ### Add pspt(long, lat) to include non_parametric spatial variables
#' ### psanova = FALSE to not decompose the effects
#' form2 <- unrate ~ partrate + agri + cons +
#'                 pspl(serv, nknots = 15) + 
#'                 pspl(empgrowth, nknots = 20) +
#'                 pspt(long, lat, nknots = c(20, 20), 
#'                      psanova = FALSE)
#' geosp <- pspatfit(form2, data = unemp_it_short)                        
#' summary(geosp)
#' ## Plot of spatial trend fitted
#' ### Create sf object to make the plot 
#' library(sf)
#' unemp_it_sf_short <- st_as_sf(dplyr::left_join(
#'                               unemp_it_short, 
#'                               map_it,  
#'                         by = c("prov" = "COD_PRO")))
#' plot_sp2d(geosp, data = unemp_it_sf_short)
#' 
#' 
#' ###############################################
#' ### Spatial semiparametric model with spatial lag
#' ### Type = "sar" for spatial lag of the dependent variable
#' geospsar <- pspatfit(form2, 
#'                      data = unemp_it_short, 
#'                      listw = lwsp_it, 
#'                      type = "sar")
#' summary(geospsar)
#' ## Plot of spatial trend fitted
#' plot_sp2d(geospsar, data = unemp_it_sf_short)
#' 
#' ### Type = "sem" for spatial lag of the noise
#' geospsem <- pspatfit(form2, 
#'                      data = unemp_it_short, 
#'                      listw = lwsp_it, 
#'                      type = "sem")
#' summary(geospsem)                  
#' 
#' #### Non-Parametric Total, Direct and Indirect impacts 
#' #### for spatial sar. First with smooth, second without smoothing
#' imp_nparvar_smooth <- impactsnopar(geospsar, 
#'                                    listw = lwsp_it, 
#'                                    viewplot = TRUE, 
#'                                    smooth = TRUE)
#' imp_nparvar_no_smooth <- impactsnopar(geospsar, 
#'                                       listw = lwsp_it, 
#'                                       viewplot = TRUE, 
#'                                       smooth = FALSE)
#' ###### Parametric Total, Direct and Indirect impacts
#' list_varpar <- c("partrate","agri","cons")
#' imp_parvar <- impactspar(geospsar, 
#'                          listw = lwsp_it)
#' summary(imp_parvar)
#' 
#' 
#' ############################################
#' ### Spatial Durbin Model with GAM
#' ### Apart from defining type = "sdm" we have to include a
#' ### formula in parameter Durbin for the Durbin part
#' 
#' durbinform <- ~ partrate + agri + pspl(serv, nknots = 13) 
#' gamsdm <- pspatfit(form1, 
#'                    data = unemp_it_short,
#'                    listw = lwsp_it,
#'                    type = "sdm",
#'                    Durbin = durbinform)
#' summary(gamsdm)
#' 
#' ############################################
#' ### Spatial Durbin Error Model with GAM
#' 
#' gamsdem <- pspatfit(form1, 
#'                     data = unemp_it_short,
#'                     listw = lwsp_it,                     type = "sdem",
#'                     Durbin = durbinform)
#' 
#' summary(gamsdem)
#' 
#' ########################################
#' ### SLX model with GAM
#' gamslx <- pspatfit(form1, 
#'                    data = unemp_it_short,
#'                    listw = lwsp_it,
#'                    type = "slx",
#'                    Durbin = durbinform)
#' summary(gamslx)
#' 
#'
#' ###############################################
#' ### Spatial semiparametric ANOVA model without spatial lag
#' ### Parameters nest_sp1 and nest_sp2 indicate the following:
#' ### the first value of the vector is the divisor of the first value of
#' ### nknots un pspt(). Gives the nknots for the first-order interaction.
#' ### The second value is the divisor for the second value in nknots giving
#' ### the the nknots per each second-order interaction. They have to be divisible
#' ### Interaction term f12 with nested basis
#' form3 <- unrate ~ partrate + agri + cons +
#'                   pspl(serv, nknots = 15) + 
#'                   pspl(empgrowth, nknots = 20) +
#'                   pspt(long, lat, nknots = c(20, 20), 
#'                        psanova = TRUE,
#'                        nest_sp1 = c(1, 2), 
#'                        nest_sp2 = c(1, 2))
#'                        
#' geospanova <- pspatfit(form3, data = unemp_it_short)
#' summary(geospanova)
#' ### Plot spatial trend and interaction (ANOVA)
#' plot_sp2d(geospanova, 
#'          data = unemp_it_sf_short, 
#'          addmain = TRUE, 
#'          addint = TRUE)
#'          
#' ## In plot_sp2d() function, if data is not an sf object, 
#' ## you must indicate the spatial coordinates variables
#' ## first plot is the whole spatial trend,
#' ## f1_main and f2_main are plots of main functions (in anova)
#' ## f12_int is the interaction function (in anova)
#' plot_sp2d(geospanova, 
#'          data = unemp_it_short, 
#'          addmain = TRUE, 
#'          addint = TRUE, 
#'          coordinates = unemp_it_short[, c("long", "lat")])
#'
#' ################################
#' ### Spatial semiparametric ANOVA model with spatial lag
#' ### Same as the previous but include type = "sar" and parameter litsw
#' ### Interaction term f12 with nested basis
#' geospanova_sar <- pspatfit(form3, 
#'                        data = unemp_it_short, 
#'                        listw = lwsp_it, 
#'                        type = "sar")
#' summary(geospanova_sar)
#'
#'
#'  ###############################################
#'  ### Spatio-temporal semiparametric ANOVA model without spatial lag
#'  ### Interaction terms f12,f1t,f2t and f12t with nested basis
#'  ### Remark: nest_sp1, nest_sp2 and nest_time must be divisors of nknots
#'  form4 <- unrate ~ partrate + agri + cons +
#'                    pspl(serv, nknots = 15) + 
#'                    pspl(empgrowth, nknots = 20) +
#'                    pspt(long, lat, year, 
#'                         nknots = c(18, 18, 12),
#'                         psanova = TRUE, 
#'                         nest_sp1 = c(1, 2, 3), 
#'                         nest_sp2 = c(1, 2, 3),
#'                         nest_time = c(1, 2, 2))
#'  sptanova <- pspatfit(form4, data = unemp_it)
#'  summary(sptanova)
#'  
#' ### Create sf object to make the plot 
#' ### of spatio-temporal trends
#' library(sf)
#' unemp_it_sf <- st_as_sf(dplyr::left_join(
#'                               unemp_it, 
#'                               map_it,  
#'                         by = c("prov" = "COD_PRO")))
#' ###### Plot spatio-temporal trends for different years
#' plot_sp3d(sptanova, data = unemp_it_sf, 
#'           time_var = "year", 
#'           time_index = c(1996, 2005, 2019),
#'           addmain = FALSE, addint = FALSE)
#' ###### Plot of spatio-temporal trend, main effects 
#' ######      and interaction effect for a year
#' plot_sp3d(sptanova, data = unemp_it_sf, 
#'           time_var = "year", 
#'           time_index = c(2019),
#'           addmain = TRUE, addint = TRUE)
#'           
#' ###### Plot of temporal trends for each province
#' plot_sptime(sptanova, 
#'             data = unemp_it, 
#'             time_var = "year", 
#'             reg_var = "prov")
#' 
#'  ### Spatio-temporal semiparametric ANOVA model with spatial lag
#'  sptanova_sar <- pspatfit(form4, data = unemp_it,
#'                           listw = lwsp_it, 
#'                           type = "sar")
#'  summary(sptanova_sar)
#'  
#'  
#'  #########################################
#'  ### Spatio-temporal semiparametric ANOVA model with spatial lag
#'  ### and temporal autorregresive noise
#'  sptanova_sar_ar1 <- pspatfit(form4, data = unemp_it, 
#'                               listw = lwsp_it, 
#'                               type = "sar",
#'                               cor = "ar1")
#'  summary(sptanova_sar_ar1)
#'  ###### Non-Parametric Total, Direct and Indirect Impacts
#'  list_varnopar <- c("serv", "empgrowth")
#'  imp_nparvar <- impactsnopar(sptanova_sar_ar1, 
#'                              listw = lwsp_it, 
#'                              viewplot = TRUE)
#'  ###### Parametric Total, Direct and Indirect Effects
#'  list_varpar <- c("partrate","agri","cons")
#'  imp_parvar <- impactspar(sptanova_sar_ar1, listw = lwsp_it)
#'  summary(imp_parvar)
#'  
#'  
#'  ###############################################
#'  ### Spatio-temporal semiparametric ANOVA model without spatial lag
#'  ### Now we repeat previous spatio-temporal model but restricting some interactions
#'  ### Interaction terms f12,f1t and f12t with nested basis
#'  ### Interaction term f2t restricted to 0
#'   form5 <- unrate ~ partrate + agri + cons +
#'                   pspl(serv, nknots=15) + pspl(empgrowth, nknots=20) +
#'                   pspt(long, lat, year, 
#'                        nknots = c(18, 18, 6), 
#'                        psanova = TRUE,
#'                        nest_sp1 = c(1, 2, 3), 
#'                        nest_sp2 = c(1, 2, 3),
#'                        nest_time = c(1, 2, 2), 
#'                        f2t_int = FALSE)
#'  sptanova2 <- pspatfit(form5, data = unemp_it)
#'  summary(sptanova2)
#'  
#' ######################  Demeaning (Within Estimators)
#' formpar <- unrate ~ partrate + agri + cons
#' param <- pspatfit(formpar, data = unemp_it, listw = lwsp_it)
#' param_dem <- pspatfit(formpar, data = unemp_it, 
#'                       listw = lwsp_it,
#'                       demean = TRUE,
#'                       index = c("prov", "year") )
#' summary(param_dem)
#' # Compare results with plm package
#' param_plm <- plm::plm(formula = formpar,
#'                       data = unemp_it,
#'                       index = c("prov", "year"),
#'                       model = "within")
#' summary(param_plm)                                              
#' param_dem_twoways <- pspatfit(formpar, 
#'                       data = unemp_it, 
#'                       listw = lwsp_it,
#'                       demean = TRUE,
#'                       eff_demean = "twoways",
#'                       index = c("prov", "year"))
#' summary(param_dem_twoways)
#' param_plm_twoways <- plm::plm(formula = formpar,
#'                       data = unemp_it,
#'                       index = c("prov", "year"),
#'                       effect = "twoways",
#'                       model = "within")
#' summary(param_plm_twoways) 
#' param_dem_time <- pspatfit(formpar, 
#'                       data = unemp_it, 
#'                       listw = lwsp_it,
#'                       demean = TRUE,
#'                       eff_demean = "time",
#'                       index = c("prov", "year"))
#' summary(param_dem_time)
#' param_plm_time <- plm::plm(formula = formpar,
#'                       data = unemp_it,
#'                       index = c("prov", "year"),
#'                       effect = "time",
#'                       model = "within")
#' summary(param_plm_time)
#' ##### Demeaning with nonparametric covariates
#' formgam <- unrate ~ partrate + agri + cons +
#'                     pspl(serv, nknots = 15) + 
#'                     pspl(empgrowth, nknots = 20)
#' gam_dem <- pspatfit(formula = formgam,
#'                       data = unemp_it,
#'                       demean = TRUE,
#'                       index = c("prov", "year"))
#' summary(gam_dem)   
#' # Compare with GAM pure without demeaning                    
#' gam <- pspatfit(formula = formgam,
#'                  data = unemp_it)
#' summary(gam)
#' # Plot of terms for GAM (without and with demeaning)
#' fit_gam <- fit_terms(gam, c("serv", "empgrowth"))
#' plot_terms(fit_gam, unemp_it)
#' fit_gam_dem <- fit_terms(gam_dem, c("serv", "empgrowth"))
#' plot_terms(fit_gam_dem, unemp_it)
#' 
#' ## Demeaning with type = "sar" model
#' gamsar_dem <- pspatfit(formula = formgam,
#'                       data = unemp_it,
#'                       type = "sar", 
#'                       listw = lwsp_it,
#'                       demean = TRUE,
#'                       index = c("prov", "year"))
#' summary(gamsar_dem)
#' }
#'                
#' @export
pspatfit <- function(formula, data, na.action,
                     listw = NULL, 
                     type = "sim", 
                     method = "eigen", 
                     Durbin = NULL,
                     zero.policy = NULL, 
                     interval = NULL, 
                     trs = NULL,
                     cor = "none",
                     dynamic = FALSE,
                     demean = FALSE,
                     eff_demean = "individual",
                     index = NULL,
                     control = list()) {
  con <- list(tol1 = 1e-6, tol2 = 5e-3, tol3 = 1e-2, 
              maxit = 200, trace = FALSE,
              optim = "llik_reml", 
              typese = "sandwich", 
              ## Values:llik_reml, score_llik_reml, llik, score_llik
              vary_init = NULL, bold = NULL,
              rho_init = 0, delta_init = 0, phi_init = 0,
              Imult = 2, cheb_q = 5, MC_p = 16L, MC_m = 30L, super = NULL,
              spamPivot = "MMD", in_coef = 0.1, type = "MC", correct = TRUE,
              trunc = TRUE, SE_method = "LU", nrho = 200, interpn = 2000,
              SElndet = NULL, LU_order = FALSE,
              pre_eig = NULL) 
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))  
  if (!(type == "sim")) { ## Change cheking when adding additional spatial models
    if (is.null(listw) || !inherits(listw, 
                                    c("listw", "Matrix", 
                                      "matrix")))
      stop("listw format unknown or NULL")
    if (inherits(listw, "listw")) {
      if (is.null(formula) || is.null(data)) {
        Wsp <- Matrix(listw2mat(listw))
      }  
    }
    if (inherits(listw, "matrix")) {
      Wsp <- Matrix(listw)
      listw <- mat2listw(Wsp)
    }  
    if (inherits(listw, "Matrix")) {
      Wsp <- listw
      listw <- mat2listw(as.matrix(Wsp))
    } 
  } else Wsp <- NULL
  if (is.null(zero.policy))
    zero.policy <- get.ZeroPolicyOption()
  can.sim <- FALSE
  if (!(is.null(listw)) && listw$style %in% c("W", "S")) {
    can.sim <- can.be.simmed(listw)
  }
  cl <- match.call()
  if (!is.null(index)) {
    if (demean) {
      if (any(grepl("pspt", formula))) 
        stop("pspt terms are not allowed with demeaning")
      formula <- update(formula, . ~ . - 1)
    }
    var_form <- all.vars(formula)
    lhs <- var_form[1]
    rhs <- paste(var_form[2:length(var_form)],
                 collapse = " + ")
    pformula <- as.formula(paste(lhs, rhs, 
                                 sep = " ~ "))
    pmodel <- plm::plm(formula = pformula,
                       data = data, 
                       index = index,
                       effect = eff_demean,
                       model = "within")
    pdata <- pmodel$model
    if (demean) {
      for (i in 1:ncol(pdata)) {
        if (!dynamic) {
          vector_meansi <- plm::Between(pdata[[i]], 
                                        effect = eff_demean)
        } else {
          if (eff_demean == "individual") {
            vector_meansi <- plm::Between(lag(pdata[[i]]), 
                                          na.rm = TRUE, 
                                          effect = eff_demean)
            # Impute NA
            vector_meansi[which(is.na(vector_meansi))] <- 
              vector_meansi[which(is.na(vector_meansi)) + 1]
          }
          else if (eff_demean == "time") {
            vector_meansi <- plm::Between(pdata[[i]], 
                                          effect = eff_demean)
          } else {
            stop("In dynamic models with demeaning, the value of 
                 argument eff_demean must be individual or time")
          } 
        }
        pdata[[i]] <- pdata[[i]] - vector_meansi
      }
    }
  }
  mf <- match.call(expand.dots = TRUE)
  m <- match(c("formula", "data", "offset"), 
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  if (is.null(index))  
    mf <- eval(mf, envir = parent.frame())
  else {
    mf <- pdata
  }
  na.act <- attr(mf, "na.action")
  y <- as.numeric(model.response(mf, "numeric"))
  nfull <- length(y)
  Xfull <- Zfull <- cfull <- NULL
  pordfull <- bdegfull <- nknotsfull <- decomfull <- NULL
  dfull_list <- list()
  if (!(type == "sim")) {
    if (!is.null(na.act)) {
      subset <- !(1:length(listw$neighbours) %in% na.act)
      listw <- subset(listw, subset, 
                      zero.policy = zero.policy)
    }
    Wsp <- Matrix(listw2mat(listw))
  }
  #### ASIGNACIONES DE DATOS A NEW ENVIRONMENT ##############
  similar <- FALSE
  env <- new.env()
  assign("nfull", nfull, envir = env)
  assign("type", type, envir = env) # type can be "sim"
  assign("cor", cor, envir = env) # temporal correlation
  assign("listw", listw, envir = env)
  if (!(is.null(listw))) {
    assign("n", length(listw$neighbours), envir = env)
    assign("similar", FALSE, envir = env)
    assign("can.sim", can.sim, envir = env)
    assign("verbose", con$trace, envir = env)
    assign("family", "SAR", envir = env)
    if (!(type == "sim")) {
      interval <- jacobianSetup(method, env, con, 
                                pre_eig = con$pre_eig,
                                trs = trs, 
                                interval = interval)
      assign("interval", interval, envir = env)
      
    }
  }
  assign("Wsp", Wsp, envir = env)
  mt <- terms(formula, specials = c("pspl", "pspt"))
  names_var <- labels(mt)
  # Careful: The dataset could include factors...
  Xmodel <- model.matrix(mt, mf)
  if (attr(mt, "intercept") == 1) {
    Xpar <- Xmodel[, c("(Intercept)"), drop = FALSE]
    names_varpar <- c("(Intercept)")
  } else {
    Xpar <- NULL
    names_varpar <- NULL
  }
  
  names_varspt <- names_var[grepl("pspt", names_var)]
  nvarspt <- length(names_varspt)
  names_varnopar <- names_var[grepl("pspl", names_var)]
  nvarnopar <- length(names_varnopar)
  names_varpar_nointercept <- names_var[!grepl("pspl", names_var) & 
                              !grepl("pspt", names_var)]
  if (length(names_varpar_nointercept) > 0) {
    names_varpar <- c(names_varpar, names_varpar_nointercept)
    nvarpar <- length(names_varpar)
    Xpar <- Xmodel[, c(names_varpar), drop = FALSE]
    colnames(Xpar) <- names_varpar
  } 
  nvarpar <- length(names_varpar)
  Xfull <- Xpar
  if (nvarspt > 0) {
    varspt <- names_varspt
    Bi <- mf[, c(varspt)]
    sp1 <- attr(Bi, "sp1")
    sp2 <- attr(Bi, "sp2")
    nsp <- length(sp1)
    time <- attr(Bi, "time")
    ntime <- attr(Bi, "ntime")
    nknotsspt <- attr(Bi, "nknots")
    if (length(nknotsspt) == 2) 
      names(nknotsspt) <- c("sp1", "sp2")
    if (length(nknotsspt) == 3) 
      names(nknotsspt) <- c("sp1", "sp2", "time")
    nknotsfull <- c(nknotsfull, nknotsspt)
    bdegspt <- attr(Bi, "bdeg")
    if (length(bdegspt) == 2) 
      names(bdegspt) <- c("sp1", "sp2")
    if (length(bdegspt) == 3) 
      names(bdegspt) <- c("sp1", "sp2", "time")
    bdegfull <- c(bdegfull, bdegspt)
    pordspt <- attr(Bi, "pord")
    if (length(pordspt) == 2) 
      names(pordspt) <- c("sp1","sp2")
    if (length(pordspt) == 3) 
      names(pordspt) <- c("sp1", "sp2", "time")
    pordfull <- c(pordfull, pordspt)
    decomspt <- attr(Bi, "decom")
    names(decomspt) <- c("spt")
    decomfull <- c(decomfull, decomspt)
    psanova <- attr(Bi, "psanova")
    nest_sp1 <- attr(Bi, "nest_sp1")
    nest_sp2 <- attr(Bi, "nest_sp2")
    nest_time <- attr(Bi, "nest_time")
    f1_main <- attr(Bi, "f1_main")
    f2_main <- attr(Bi, "f2_main")
    ft_main <- attr(Bi, "ft_main")
    f12_int <- attr(Bi, "f12_int")
    f1t_int <- attr(Bi, "f1t_int")
    f2t_int <- attr(Bi, "f2t_int")
    f12t_int <- attr(Bi, "f12t_int")
    Bsptfull <- Bspt(sp1 = sp1, sp2 = sp2, time = time, 
                     nfull = nfull, ntime = ntime, 
                     psanova = psanova, Bi = Bi,
                     bdegspt = bdegspt)
    sp1_short <- Bsptfull$sp1 
    sp2_short <- Bsptfull$sp2
    nsp_short <- length(sp1_short)
    time_short <- Bsptfull$time
    Bsptlist <- Bsptfull$Bsptlist
    #rm(Bsptfull, Bi)
    XZsptlist <- B_XZ_spt(sp1 = sp1_short, 
                          sp2 = sp2_short, 
                          time = time_short, 
                          pordspt = pordspt, 
                          psanova = psanova,
                          decomspt = decomspt, 
                          f1_main = f1_main,
                          f2_main = f2_main, 
                          ft_main = ft_main,
                          f12_int = f12_int, 
                          f1t_int = f1t_int,
                          f2t_int = f2t_int, 
                          f12t_int = f12t_int,
                          Bsptlist = Bsptlist)
    Xspt <- XZsptlist$Xspt
    if (attr(mt, "intercept") == 1) 
      Xspt <- Xspt[, -c(1), drop = FALSE] # Remove intercept
    Zspt <- XZsptlist$Zspt
    dsptlist <- XZsptlist$dsptlist
    cspt <- unlist(XZsptlist$csptlist)
    Xfull <- cbind(Xfull, Xspt)
    Zfull <- cbind(Zfull, Zspt)
    dfull_list <- c(dfull_list, dsptlist)
    cfull <- c(cfull, cspt)
  } else {
    sp1 <- sp2 <- time <- NULL
    Xspt <- Zspt <- dsptlist <- cspt <- NULL
    nknotsspt <- pordspt <- bdegspt <- decomspt <- NULL
  } # if (!is.null(names_varspt))
  if (inherits(mf, "pdata.frame")) {
    # It is a panel...
    nsp <- length(unique(plm::index(mf)[, 1]))
    nt <- length(unique(plm::index(mf)[, 2]))
  } else if (!is.null(Wsp)) {
    nsp <- nrow(Wsp)
  } else if (!is.null(sp1)) {
    if (!is.null(time)) #3d case
      nsp <- length(unique(sp1)) 
    else nsp <- length(sp1) #2d case
  } else nsp <- nfull
  if (!is.null(time) && (!inherits(mf, "pdata.frame"))) 
    nt <- length(unique(time))  
  else if (nfull > nsp) {
    ##  spatio-temporal data but working in 2d...
    nt <- nfull %/% nsp
    if ((nfull %% nsp) > 0) 
      stop("Dimensions of Wsp matrix and data do not agree")
  } else nt <- 1
  assign("nsp", nsp, envir = env)
  assign("nt", nt, envir = env)  
  if (type %in% c("slx", "sdm", "sdem")) {
    if (is.null(Durbin)) Durbin <- update(formula, NULL ~ . )
    # Add intercept in Durbin formula to use it in factor case...
    # After it is removed in Xpar
    Durbin <- update(Durbin, NULL ~ . + 1) 
    mtDurbin <- terms(Durbin, specials = c("pspl"))
    names_var_Wlag <- labels(mtDurbin)
    names_varnopar_Wlag <- names_var_Wlag[grepl("pspl", 
                                                names_var_Wlag)]
    nvarnoparWlag <- length(names_varnopar_Wlag)
    # if (any(!(names_var_Wlag %in% names_var))) 
    #   warning("Some of the variables in Durbin formula are not 
    #           in the general formula")
    # Build a dataframe for Durbin formula
    vindex <- rep(0, length(colnames(data)))
    for (i in 1:length(names_var_Wlag)) {
      name_i <- names_var_Wlag[i]
      idx_i <- str_detect(name_i, colnames(data))
      vindex <- vindex + as.numeric(idx_i)
    }
    vindex <- as.logical(vindex)
    dfDurbin <- data[, vindex]
    if (inherits(dfDurbin, "sf"))
      dfDurbin <- st_drop_geometry(dfDurbin)
    # Take care of spatio-temporal data with Wsp
    It <- Diagonal(nt)
    Wspt <- kronecker(Wsp, It)
    listwspt <- mat2listw(as.matrix(Wspt))
    if (length(names_var_Wlag) > length(names_varnopar_Wlag)) {
      # There could be factors between parametric variables in Durbin formula
      XparWlag <- model.matrix(mtDurbin, dfDurbin)
      names_varpar_Wlag <- colnames(XparWlag)[!grepl("pspl", 
                                                     colnames(XparWlag)) &
                                                !grepl("Intercept",
                                                       colnames(XparWlag))]
      nvarparWlag <- length(names_varpar_Wlag)
      XparWlag <- XparWlag[, names_varpar_Wlag, drop = FALSE]
      XparWlag <- create_WX(XparWlag, listw = listwspt , 
                            zero.policy = zero.policy,
                            prefix = "Wlag")
      Xpar <- cbind(Xpar, XparWlag)
      nvarpar <- ncol(Xpar)
      names_varpar <- colnames(Xpar)
    }
    if (nvarnoparWlag > 0) {
      for (i in 1:nvarnoparWlag) {
        name_i <- names_varnopar_Wlag[i]
        idx_i <- str_detect(name_i, colnames(dfDurbin))
        var_i <- dfDurbin[, idx_i]
        if (inherits(var_i, "sf")) var_i <- st_drop_geometry(var_i)
        var_i <- as.matrix(var_i)
        colnames(var_i) <- colnames(dfDurbin)[idx_i]
        Wvar_i <- create_WX(var_i, listw = listwspt , 
                            zero.policy = zero.policy,
                            prefix = "Wlag")
        dfDurbin[, idx_i]  <- Wvar_i
      }
      colnames(dfDurbin) <- gsub("Wlag.", "", colnames(dfDurbin))
      mfDurbin <- model.frame(Durbin, 
                              data = dfDurbin, 
                              drop.unused.levels = TRUE) 
      mfDurbin_nopar <- mfDurbin[, grepl("pspl", colnames(mfDurbin)), 
                                 drop = FALSE]
      expr1 <- paste0("Wlag.pspl","\\(", collapse = "")
      expr2 <- paste0("pspl","\\(","Wlag.", collapse = "")
      
      colnames(mfDurbin_nopar) <- paste0("Wlag.", colnames(mfDurbin_nopar))
      colnames(mfDurbin_nopar) <- gsub(expr1, expr2, 
                                              colnames(mfDurbin_nopar))
      names_varnopar_Wlag <- paste0("Wlag.", names_varnopar_Wlag)
      names_varnopar_Wlag <- gsub(expr1, expr2,
                                  names_varnopar_Wlag)
      names_varnopar <- c(names_varnopar, names_varnopar_Wlag)
      nvarnopar <- nvarnopar + nvarnoparWlag
      mf <- cbind(mf, mfDurbin_nopar)
      rm(expr1, expr2)
    }
  }
  if ((dynamic == TRUE) && nt > 1) {
    my <- matrix(y, nrow = nt, ncol = nsp, byrow = FALSE)
    tlag_my <- rbind(my[2:(nt), ], NA)
    tlagy <- matrix(tlag_my, nrow = nsp*nt, ncol = 1, byrow = TRUE)
    colnames(tlagy) <- paste("tlag.", all.vars(formula)[1], sep = "")
    idxnatlagy <- !is.na(tlagy[, 1])
    Xpar <- cbind(Xpar, tlagy)
  }
  #if (!is.null(Xpar)) Xfull <- cbind(Xfull, Xpar)
  if (nvarnopar > 0) {
    Xnopar <- Znopar <-  cnopar <- NULL
    nknotsnopar <- bdegnopar <- pordnopar <- decomnopar <- NULL
    dnoparlist <- list()
    dnoparlistnames <- vector()
    for (i in 1:length(names_varnopar)) {
      varnopar <- names_varnopar[i]
      Bi <- mf[, c(varnopar)]
      colnames(Bi) <- paste(names_varnopar[i], 1:ncol(Bi), sep = ".")
      nknots_i <- attr(Bi, "nknots")
      pord_i <- attr(Bi, "pord")
      bdeg_i <- attr(Bi, "bdeg")
      decom_i <- attr(Bi, "decom")
      x_i <- attr(Bi, "x")
      names(nknots_i) <- names(pord_i) <- varnopar
      names(bdeg_i) <- names(decom_i) <- varnopar
      nknotsnopar <- c(nknotsnopar, nknots_i)
      pordnopar <- c(pordnopar, pord_i)
      bdegnopar <- c(bdegnopar, bdeg_i)
      decomnopar <- c(decomnopar, decom_i)
      BtoXZ <- B_XZ(Bi, x = x_i, pord = pord_i, 
                    decom = decom_i)
      Xi <- as.matrix(BtoXZ$X) 
      if (attr(mt, "intercept") == 1) 
        Xi <- BtoXZ$X[, -c(1), drop = FALSE] # Remove intercept
      colnames(Xi) <- paste(names_varnopar[i], 1:ncol(Xi), sep = ".")
      Zi <- BtoXZ$Z
      colnames(Zi) <- paste(names_varnopar[i], 1:ncol(Zi), sep = ".")
      dnoparlist[[i]] <- BtoXZ$d
      dnoparlistnames <- c(dnoparlistnames,
                           paste("d_", names_varnopar[i], sep = ""))
      ci <- ncol(Bi)
      names(ci) <- paste(names_varnopar[i], sep = "")
      Xnopar <- cbind(Xnopar, Xi)
      Znopar <- cbind(Znopar, Zi)
      cnopar <- c(cnopar, ci)
      rm(BtoXZ, Bi, Xi, Zi, ci, nknots_i, pord_i, bdeg_i,
         decom_i, varnopar)
    }
    names(dnoparlist) <- dnoparlistnames
    rm(dnoparlistnames)
    Xfull <- cbind(Xfull, Xnopar)
    Zfull <- cbind(Zfull, Znopar)
    dfull_list <- c(dfull_list, dnoparlist)
    cfull <- c(cfull, cnopar)
    nknotsfull <- c(nknotsfull, nknotsnopar)
    pordfull <- c(pordfull, pordnopar)
    bdegfull <- c(bdegfull, bdegnopar)
    decomfull <- c(decomfull, decomnopar)
  } else {
    Xnopar <- Znopar <- dnoparlist <- cnopar <- NULL
    nknotsnopar <- pordnopar <- bdegnopar <- decomnopar <- NULL
  }    # end if(!is.null(names_varnopar))
  if (dynamic == TRUE) {
    # Remove NA's
    y <- y[idxnatlagy]
    Xfull <- Xfull[idxnatlagy, ]
    Zfull <- Zfull[idxnatlagy, ]
    # Change nt to nt-1
    nt <- nt - 1
    nfull <- nsp*nt
    assign("nt", nt, envir = env)
    assign("nfull", nfull, envir = env)
    assign("idxnatlagy", idxnatlagy, envir = env)
  } else assign("idxnatlagy", NULL, envir = env)
  assign("y", y, envir = env)
  assign("Xfull", Xfull, envir = env)
  assign("Zfull", Zfull, envir = env)
  assign("nvarpar", nvarpar, envir = env)
  assign("names_varpar", names_varpar, envir = env)
  assign("nvarnopar", nvarnopar, envir = env)
  assign("names_nvarnopar", names_varnopar, 
         envir = env)
  assign("nvarspt", nvarspt, envir = env)
  if (nvarspt > 0) {
    assign("sp1", sp1, envir = env)
    assign("sp2", sp1, envir = env)
    assign("time", time, envir = env)
    assign("Bsptlist", Bsptlist, envir = env)
    assign("cspt", cspt, envir = env)
    assign("dsptlist", dsptlist, envir = env)
    assign("bdegspt", bdegspt, envir = env)
    assign("pordspt", pordspt, envir = env)
    assign("nknotsspt", nknotsspt, envir = env)
    assign("psanova", psanova, envir = env)
    assign("f1_main", f1_main, envir = env)
    assign("f2_main", f2_main, envir = env)
    assign("ft_main", ft_main, envir = env)
    assign("f12_int", f12_int, envir = env)
    assign("f1t_int", f1t_int, envir = env)
    assign("f2t_int", f2t_int, envir = env)
    assign("f12t_int", f12t_int, envir = env)
  }
  if (nvarnopar > 0) {
    assign("cnopar", cnopar, envir = env)
    assign("dnoparlist", dnoparlist, envir = env)
    assign("bdegnopar", bdegnopar, envir = env)
    assign("pordnopar", pordnopar, envir = env)
    assign("nknotsnopar", nknotsnopar, envir = env)
    assign("names_varnopar", names_varnopar, envir = env)
  }
  if (is.null(con$vary_init)) con$vary_init <- var(y)
  message("\nFitting Model...\n")
  model_fit <- fit_pspat(env, con)
  mt_terms <- attr(mt, "term.labels")
  model_fit$contrasts <- attr(Xpar, "contrasts")
  model_fit$xlevels <- .getXlevels(mt, mf)
  model_fit$call <- cl
  model_fit$terms <- mt
  model_fit$data <- data
  model_fit$type <- type
  model_fit$cor <- cor
  model_fit$demean <- demean
  model_fit$eff_demean <- eff_demean
  model_fit$index <- index
  model_fit$dynamic <- dynamic
  if (!is.null(listw)) model_fit$listw <- listw
  else model_fit$listw <- NULL
  model_fit$Durbin <- Durbin
  model_fit$X <- Xfull
  model_fit$Z <- Zfull
  model_fit$y <- y
  if (dynamic) {
    # Add the last observation
    nt <- nt + 1
    nfull <- nsp*nt
  }
  model_fit$nfull <- nfull
  model_fit$nsp <- nsp
  model_fit$nt <- nt
  model_fit$df.residual <- length(y) - model_fit$edftot
  model_fit$fitted.values_Ay <- model_fit$fit_A1y
  model_fit$fit_A1y <- NULL
  if (con$typese == "bayesian") {
    model_fit$se_bfixed <- model_fit$sefr_bfixed
    model_fit$se_brandom <- model_fit$sefr_brandom
    model_fit$se_fitted.values <- model_fit$sefr_fitted.values
    model_fit$se_fitted.values_Ay <- model_fit$sefr_fit_A1y
  } else {
    model_fit$se_bfixed <- model_fit$seby_bfixed
    model_fit$se_brandom <- model_fit$seby_brandom
    model_fit$se_fitted.values <- model_fit$seby_fitted.values
    model_fit$se_fitted.values_Ay <- model_fit$seby_fit_A1y
  }
  model_fit$sefr_bfixed <- NULL
  model_fit$seby_bfixed <- NULL
  model_fit$sefr_brandom <- NULL
  model_fit$seby_brandom <- NULL
  model_fit$seby_fitted.values <- NULL
  model_fit$sefr_fitted.values <- NULL
  model_fit$sefr_fit_A1y <- NULL
  model_fit$seby_fit_A1y <- NULL
  model_fit$tauspt <- NULL
  model_fit$taunopar <- NULL
  class(model_fit) <- c("pspatreg", "lm")
  model_fit
}
