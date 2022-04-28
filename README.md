# pspatreg

<!-- badges: start -->

\[![CRAN status](https://www.r-pkg.org/badges/version-ago/pspatreg)
\[![CRAN
downloads-last-month](https://cranlogs.r-pkg.org/badges/last-month/pspatreg)
\[![CRAN
downloads-grand-total](https://cranlogs.r-pkg.org/badges/grand-total/pspatreg)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

The goal of **pspatreg** is to provide tools for the estimation and inference of semiparametric
spatial and spatio-temporal regression models with spatial lags.

## Installation

You can install the released version of **pspatreg** from
[CRAN](https://CRAN.R-project.org) with:

``` r
#install.packages("pspatreg")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
#devtools::install_github("rominsal/pspatreg")
```

## Abstract

In recent times there is a growing interest en the development of flexible semiparametric models
espacially for the modelling of spatial and spatio-temporal data. The package **pspatreg** is intended 
to provide some tools, based on P-Splines methodology, for the estimation and inference with these type of 
models from a spatial econometrics perspective. 

In this sense, the type of spatial econometric models included are ps-sar, ps-sem, 
ps.sarar, ps-slx, or ps-Durbin. All the specifications can include parametric and non-parametric covariates, 
spatial or spatio-temporal non-parametric trends and spatial lags of the dependent variable and/or 
the noise of the model. 

The non-parametric terms (either trends or covariates) are modeled using P-Splines. Furthermore, 
the non-parametric trend can be decomposed in an ANOVA way including main and interactions 
effects of 2nd and 3rd order (for spatio-temporal trends). 

For every model, the estimation method can be Restricted Maximum Likelihood (REML) or Maximum Likelihood (ML).

In addition to estimation and inference procedures, he package also provides functions to compute:

- Non-parametric and parametric impacts (including inference about them).
- Plots of non-parametric impacts and non-parametric terms. 
- Plots of the spatial or spatio-temporal trends (including main and interaction terms for ANOVA specifications).
- Usual methods for summary, anova, print, cov, residuals, ...

## Vignettes

The package includes three vignettes:

- Vignette A: Introduction to the use of the package and the statistical theory behind the models.
- Vignette B: Includes several examples using the well-known AMES database for pure spatial data. The examples include both spatial trends and spatial lags (as usual in spatial econometric specifications). The vignette B also compare results of **pspatreg** with the well-kown **spatialreg** package for parametric specifications. 
- Vignette C: Include examples with spatial panel data (spatio-temporal data). The examples include spatio-temporal trends, non-parametric covariates and spatial lags. The vignette C also compare results of **pspatreg** with **splm** package for parametric specifications.  
