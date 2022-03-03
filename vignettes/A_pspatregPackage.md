An introduction to pspatreg
================
Román Mínguez (UCLM), Roberto Basile (UNIVAQ), María Durbán (UC3M)
2022-03-03

# Introduction

**pspatreg** is a package that fits penalized spline (PS) semiparametric
static and dynamic spatial autoregressive models via Restricted (or
Residual) Maximum Likelihood (REML) and Maximum Likelihood (ML). This
approach combines penalized regression spline methods (Eilers, Marx, and
Durbán 2015) with standard spatial autoregressive models (such as SAR,
SEM and SDM). These types of models are thoroughly discussed in
Mı́nguez, Basile, and Durbán (2020); see also Montero, Mı́nguez, and
Durbán (2012), Basile et al. (2014), and Hoshino (2018).

These models are very flexible since they make it possible to include
within the same specification: *i*) spatial autoregressive terms
(i.e. spatial lags of dependent and independent variables as well as
spatial error terms) to capture spatial interaction or network effects;
*ii*) time lags of the dependent variable to capture persistence
effects; *iii*) parametric and nonparametric (smooth) terms to identify
nonlinear relationships between the response variable and the
covariates; *iv*) a spatio-temporal trend, i.e. a smooth interaction
between the spatial coordinates and the time trend, to capture
site-specific nonlinear time trends.

The proposed method also allows the user to apply an ANOVA decomposition
of the spatio-temporal trend into several components (spatial and
temporal main effects, and second- and third-order interactions between
them), which gives further insights into the dynamics of the data. Thus,
we use the acronym PS-ANOVA-SAR (SEM, SDM, SLX) for the new data
generating process (DGP) proposed. The use of nested B-spline bases for
the interaction components of the spatio-temporal trend contributes to
the efficiency of the fitting procedure without compromising the
goodness of fit of the model. Finally, we also consider an extension of
the PS-ANOVA-SAR (SEM, SDM, SLX) including the time lag of the dependent
variable (dynamic spatial model) and/or a first-order time series
autoregressive term process (AR1) in the noise to accommodate residual
serial correlation.

## The semiparametric spatial autoregressive model

In a very general form, the semiparametric spatial autoregressive model
reads as:

  
![y\_{it}=\\alpha y\_{it-1}+\\rho \\sum\_{j=1}^N w\_{ij,N} y\_{jt} +
\\pi \\sum\_{j=1}^N w\_{ij,N} y\_{jt-1}+ \\\\&#10;\\sum\_{k=1}^K
\\beta^\*\_k x^\*\_{k,it} + \\sum\_{k=1}^K \\gamma^\*\_k \\sum\_{j=1}^N
w\_{ij,N} x^\*\_{k,jt}+&#10;\\sum\_{\\delta=1}^\\Delta
g\_\\delta(x\_{\\delta, it}) + \\sum\_{\\delta=1}^\\Delta
h\_\\delta\\left(\\sum\_{j=1}^N w\_{ij,N} x\_{\\delta,jt}\\right) +
\\\\&#10;\\widetilde{
f}(s\_{1i},s\_{2i},\\tau\_t)+\\epsilon\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bit%7D%3D%5Calpha%20y_%7Bit-1%7D%2B%5Crho%20%5Csum_%7Bj%3D1%7D%5EN%20w_%7Bij%2CN%7D%20y_%7Bjt%7D%20%2B%20%5Cpi%20%5Csum_%7Bj%3D1%7D%5EN%20w_%7Bij%2CN%7D%20y_%7Bjt-1%7D%2B%20%5C%5C%0A%5Csum_%7Bk%3D1%7D%5EK%20%5Cbeta%5E%2A_k%20x%5E%2A_%7Bk%2Cit%7D%20%2B%20%5Csum_%7Bk%3D1%7D%5EK%20%5Cgamma%5E%2A_k%20%5Csum_%7Bj%3D1%7D%5EN%20w_%7Bij%2CN%7D%20x%5E%2A_%7Bk%2Cjt%7D%2B%0A%5Csum_%7B%5Cdelta%3D1%7D%5E%5CDelta%20g_%5Cdelta%28x_%7B%5Cdelta%2C%20it%7D%29%20%2B%20%5Csum_%7B%5Cdelta%3D1%7D%5E%5CDelta%20h_%5Cdelta%5Cleft%28%5Csum_%7Bj%3D1%7D%5EN%20w_%7Bij%2CN%7D%20x_%7B%5Cdelta%2Cjt%7D%5Cright%29%20%2B%20%5C%5C%0A%5Cwidetilde%7B%20f%7D%28s_%7B1i%7D%2Cs_%7B2i%7D%2C%5Ctau_t%29%2B%5Cepsilon_%7Bit%7D
"y_{it}=\\alpha y_{it-1}+\\rho \\sum_{j=1}^N w_{ij,N} y_{jt} + \\pi \\sum_{j=1}^N w_{ij,N} y_{jt-1}+ \\\\
\\sum_{k=1}^K \\beta^*_k x^*_{k,it} + \\sum_{k=1}^K \\gamma^*_k \\sum_{j=1}^N w_{ij,N} x^*_{k,jt}+
\\sum_{\\delta=1}^\\Delta g_\\delta(x_{\\delta, it}) + \\sum_{\\delta=1}^\\Delta h_\\delta\\left(\\sum_{j=1}^N w_{ij,N} x_{\\delta,jt}\\right) + \\\\
\\widetilde{ f}(s_{1i},s_{2i},\\tau_t)+\\epsilon_{it}")  

  
![\\epsilon\_{it}=\\theta \\sum\_{j=1}^N w\_{ij,N}\\epsilon\_{jt}+\\phi
\\epsilon\_{it-1}+u\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bit%7D%3D%5Ctheta%20%5Csum_%7Bj%3D1%7D%5EN%20w_%7Bij%2CN%7D%5Cepsilon_%7Bjt%7D%2B%5Cphi%20%5Cepsilon_%7Bit-1%7D%2Bu_%7Bit%7D
"\\epsilon_{it}=\\theta \\sum_{j=1}^N w_{ij,N}\\epsilon_{jt}+\\phi \\epsilon_{it-1}+u_{it}")  
  
![u\_{it} \\sim
i.i.d.(0,\\sigma^2\_u)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u_%7Bit%7D%20%5Csim%20i.i.d.%280%2C%5Csigma%5E2_u%29
"u_{it} \\sim i.i.d.(0,\\sigma^2_u)")  

where
![(y\_{it},x^\*\_{1,it},...,x^\*\_{K,it},x\_{1,it},...,,x\_{\\Delta,it})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28y_%7Bit%7D%2Cx%5E%2A_%7B1%2Cit%7D%2C...%2Cx%5E%2A_%7BK%2Cit%7D%2Cx_%7B1%2Cit%7D%2C...%2C%2Cx_%7B%5CDelta%2Cit%7D%29
"(y_{it},x^*_{1,it},...,x^*_{K,it},x_{1,it},...,,x_{\\Delta,it})") are
data observed for a sample of spatial units
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i
"i")
![(i=1,...,N)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28i%3D1%2C...%2CN%29
"(i=1,...,N)") in each time period
![t=1,...,T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t%3D1%2C...%2CT
"t=1,...,T"). The terms
![(s\_{1i},s\_{2i})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28s_%7B1i%7D%2Cs_%7B2i%7D%29
"(s_{1i},s_{2i})") indicate the spatial coordinates (latitude and
longitude) of the spatial unit
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i
"i") (either a spatial point or the centroid of a spatial polygon), and
![W\_{ij}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W_%7Bij%7D
"W_{ij}") is an element of a row-standardized spatial weights matrix.
The terms included into the model are:

  - ![y\_{it-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bit-1%7D
    "y_{it-1}") = the time lag of
    ![y\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bit%7D
    "y_{it}");

  - ![\\sum\_{j=1}^N w\_{ij,N}
    y\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csum_%7Bj%3D1%7D%5EN%20w_%7Bij%2CN%7D%20y_%7Bit%7D
    "\\sum_{j=1}^N w_{ij,N} y_{it}") = the contemporaneous spatial lag
    of
    ![y\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bit%7D
    "y_{it}");

  - ![\\sum\_{j=1}^N w\_{ij,N}
    y\_{it-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csum_%7Bj%3D1%7D%5EN%20w_%7Bij%2CN%7D%20y_%7Bit-1%7D
    "\\sum_{j=1}^N w_{ij,N} y_{it-1}") = the time lag of the spatial lag
    of
    ![y\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bit%7D
    "y_{it}");

  - ![\\alpha](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha
    "\\alpha"),
    ![\\rho](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Crho
    "\\rho"),
    ![\\pi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpi
    "\\pi"),
    ![\\beta^\*\_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta%5E%2A_k
    "\\beta^*_k"), and
    ![\\gamma^\*\_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cgamma%5E%2A_k
    "\\gamma^*_k") = fixed parameters;

  - ![\\sum\_{k=1}^K \\beta^\*\_k
    x^\*\_{k,it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csum_%7Bk%3D1%7D%5EK%20%5Cbeta%5E%2A_k%20x%5E%2A_%7Bk%2Cit%7D
    "\\sum_{k=1}^K \\beta^*_k x^*_{k,it}") and ![\\sum\_{k=1}^K
    \\gamma^\*\_k \\sum\_{j=1}^N w\_{ij,N}
    x^\*\_{k,jt}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csum_%7Bk%3D1%7D%5EK%20%5Cgamma%5E%2A_k%20%5Csum_%7Bj%3D1%7D%5EN%20w_%7Bij%2CN%7D%20x%5E%2A_%7Bk%2Cjt%7D
    "\\sum_{k=1}^K \\gamma^*_k \\sum_{j=1}^N w_{ij,N} x^*_{k,jt}") =
    parametric linear terms of some covariates
    ![x\_{k,it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_%7Bk%2Cit%7D
    "x_{k,it}") and of their spatial lags;

  - ![g\_\\delta(.)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g_%5Cdelta%28.%29
    "g_\\delta(.)") and
    ![h\_\\gamma(.)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;h_%5Cgamma%28.%29
    "h_\\gamma(.)") = nonparametric smooth functions of other covariates
    and of their spatial lags (they can also accommodate varying
    coefficient terms, smooth interaction between covariates,
    factor-by-curve intercations, and so on);

  - ![\\widetilde{
    f}(s\_{1i},s\_{2i},\\tau\_t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidetilde%7B%20f%7D%28s_%7B1i%7D%2Cs_%7B2i%7D%2C%5Ctau_t%29
    "\\widetilde{ f}(s_{1i},s_{2i},\\tau_t)") = an unknown nonparametric
    spatio-temporal trend;

  - ![\\epsilon\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bit%7D
    "\\epsilon_{it}") = an idiosyncratic error term that can follow a
    spatial autoregressive process ![\\epsilon\_{it}=\\theta
    \\sum\_{j=1}^N
    w\_{ij,N}\\epsilon\_{it}+u\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bit%7D%3D%5Ctheta%20%5Csum_%7Bj%3D1%7D%5EN%20w_%7Bij%2CN%7D%5Cepsilon_%7Bit%7D%2Bu_%7Bit%7D
    "\\epsilon_{it}=\\theta \\sum_{j=1}^N w_{ij,N}\\epsilon_{it}+u_{it}")
    with ![u\_{it}\\sim
    N(0,\\sigma^2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u_%7Bit%7D%5Csim%20N%280%2C%5Csigma%5E2%29
    "u_{it}\\sim N(0,\\sigma^2)") (SEM), a time series autoregressive
    AR(1) process, i.e., ![\\epsilon\_{it}=\\phi
    \\epsilon\_{it-1}+u\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bit%7D%3D%5Cphi%20%5Cepsilon_%7Bit-1%7D%2Bu_%7Bit%7D
    "\\epsilon_{it}=\\phi \\epsilon_{it-1}+u_{it}") with ![u\_{it}\\sim
    N(0,\\sigma^2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u_%7Bit%7D%5Csim%20N%280%2C%5Csigma%5E2%29
    "u_{it}\\sim N(0,\\sigma^2)"), or both.

Obviously, this very general specification is hardly suitable in real
data applications. However, it is worth noticing that it nests most of
the spatial additive models which can be used in practice. For example,
a more suitable semiparametric static model reads as:

  
![y\_{it}=\\rho \\sum\_{j=1}^N w\_{ij,N} y\_{jt} + &#10;\\sum\_{k=1}^K
\\beta^\*\_k x^\*\_{k,it} + &#10;\\sum\_{\\delta=1}^\\Delta
g\_\\delta(x\_{\\delta, it}) +&#10;\\widetilde{
f}(s\_{1i},s\_{2i},\\tau\_t)+\\epsilon\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bit%7D%3D%5Crho%20%5Csum_%7Bj%3D1%7D%5EN%20w_%7Bij%2CN%7D%20y_%7Bjt%7D%20%2B%20%0A%5Csum_%7Bk%3D1%7D%5EK%20%5Cbeta%5E%2A_k%20x%5E%2A_%7Bk%2Cit%7D%20%2B%20%0A%5Csum_%7B%5Cdelta%3D1%7D%5E%5CDelta%20g_%5Cdelta%28x_%7B%5Cdelta%2C%20it%7D%29%20%2B%0A%5Cwidetilde%7B%20f%7D%28s_%7B1i%7D%2Cs_%7B2i%7D%2C%5Ctau_t%29%2B%5Cepsilon_%7Bit%7D
"y_{it}=\\rho \\sum_{j=1}^N w_{ij,N} y_{jt} + 
\\sum_{k=1}^K \\beta^*_k x^*_{k,it} + 
\\sum_{\\delta=1}^\\Delta g_\\delta(x_{\\delta, it}) +
\\widetilde{ f}(s_{1i},s_{2i},\\tau_t)+\\epsilon_{it}")  

  
![\\epsilon\_{it}=\\phi
\\epsilon\_{it-1}+u\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bit%7D%3D%5Cphi%20%5Cepsilon_%7Bit-1%7D%2Bu_%7Bit%7D
"\\epsilon_{it}=\\phi \\epsilon_{it-1}+u_{it}")  
  
![u\_{it} \\sim
i.i.d.(0,\\sigma^2\_u)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u_%7Bit%7D%20%5Csim%20i.i.d.%280%2C%5Csigma%5E2_u%29
"u_{it} \\sim i.i.d.(0,\\sigma^2_u)")  
This semiparametric SAR model turns out to be extremely useful to
capture interactive spatial and temporal unobserved heterogeneity when
the last one is smoothly distributed over space and time (Mı́nguez,
Basile, and Durbán 2020). The dynamic extension (including
![y\_{it-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bit-1%7D
"y_{it-1}") and ![\\sum\_{j=1}^N w\_{ij,N}
y\_{it-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csum_%7Bj%3D1%7D%5EN%20w_%7Bij%2CN%7D%20y_%7Bit-1%7D
"\\sum_{j=1}^N w_{ij,N} y_{it-1}")) is also very promising and merits
further theoretical investigation. Finally, the following semiparametric
SAR model is very useful for modeling cross-setional spatial data taking
into account nonlinearities, spatial dependence and spatial
heterogeneity:

  
![y\_{i}=\\rho \\sum\_{j=1}^N w\_{ij,N} y\_{j} + \\sum\_{k=1}^K
\\beta^\*\_k x^\*\_{k,i} + \\sum\_{\\delta=1}^\\Delta
g\_\\delta(x\_{\\delta, i}) + &#10;\\widetilde{
f}(s\_{1i},s\_{2i})+\\epsilon\_{i}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%7D%3D%5Crho%20%5Csum_%7Bj%3D1%7D%5EN%20w_%7Bij%2CN%7D%20y_%7Bj%7D%20%2B%20%5Csum_%7Bk%3D1%7D%5EK%20%5Cbeta%5E%2A_k%20x%5E%2A_%7Bk%2Ci%7D%20%2B%20%5Csum_%7B%5Cdelta%3D1%7D%5E%5CDelta%20g_%5Cdelta%28x_%7B%5Cdelta%2C%20i%7D%29%20%2B%20%0A%5Cwidetilde%7B%20f%7D%28s_%7B1i%7D%2Cs_%7B2i%7D%29%2B%5Cepsilon_%7Bi%7D
"y_{i}=\\rho \\sum_{j=1}^N w_{ij,N} y_{j} + \\sum_{k=1}^K \\beta^*_k x^*_{k,i} + \\sum_{\\delta=1}^\\Delta g_\\delta(x_{\\delta, i}) + 
\\widetilde{ f}(s_{1i},s_{2i})+\\epsilon_{i}")  

  
![\\epsilon\_{i}\\sim
i.i.d.(0,\\sigma^2\_\\epsilon)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bi%7D%5Csim%20i.i.d.%280%2C%5Csigma%5E2_%5Cepsilon%29
"\\epsilon_{i}\\sim i.i.d.(0,\\sigma^2_\\epsilon)")  
In many situations the spatio-temporal trend to be estimated can be
complex, and the use of a single multidimensional smooth function may
not be flexible enough to capture the structure in the data. To solve
this problem, an ANOVA-type decomposition of ![\\widetilde{
f}(s\_{1i},s\_{2i},\\tau\_t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidetilde%7B%20f%7D%28s_%7B1i%7D%2Cs_%7B2i%7D%2C%5Ctau_t%29
"\\widetilde{ f}(s_{1i},s_{2i},\\tau_t)") can be used, where spatial and
temporal main effects, and second- and third-order interactions between
them can be identified:

  
![\\widetilde{
f}(s\_{1i},s\_{2i},\\tau\_t)=f\_1(s\_{1i})+f\_2(s\_{2i})+f\_{\\tau}(\\tau\_t)+
\\\\
f\_{1,2}(s\_{1i},s\_{2i})+f\_{1,\\tau}(s\_{1i},\\tau\_t)+f\_{2,\\tau}+(s\_{2i},\\tau\_t)+f\_{1,2,\\tau}(s\_{1i},s\_{2i},\\tau\_t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidetilde%7B%20f%7D%28s_%7B1i%7D%2Cs_%7B2i%7D%2C%5Ctau_t%29%3Df_1%28s_%7B1i%7D%29%2Bf_2%28s_%7B2i%7D%29%2Bf_%7B%5Ctau%7D%28%5Ctau_t%29%2B%20%5C%5C%20f_%7B1%2C2%7D%28s_%7B1i%7D%2Cs_%7B2i%7D%29%2Bf_%7B1%2C%5Ctau%7D%28s_%7B1i%7D%2C%5Ctau_t%29%2Bf_%7B2%2C%5Ctau%7D%2B%28s_%7B2i%7D%2C%5Ctau_t%29%2Bf_%7B1%2C2%2C%5Ctau%7D%28s_%7B1i%7D%2Cs_%7B2i%7D%2C%5Ctau_t%29
"\\widetilde{ f}(s_{1i},s_{2i},\\tau_t)=f_1(s_{1i})+f_2(s_{2i})+f_{\\tau}(\\tau_t)+ \\\\ f_{1,2}(s_{1i},s_{2i})+f_{1,\\tau}(s_{1i},\\tau_t)+f_{2,\\tau}+(s_{2i},\\tau_t)+f_{1,2,\\tau}(s_{1i},s_{2i},\\tau_t)")  

First, the geoadditive terms given by
![f\_1(s\_{1i}),f\_2(s\_{2i}),f\_{1,2}(s\_{1i},s\_{2i})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f_1%28s_%7B1i%7D%29%2Cf_2%28s_%7B2i%7D%29%2Cf_%7B1%2C2%7D%28s_%7B1i%7D%2Cs_%7B2i%7D%29
"f_1(s_{1i}),f_2(s_{2i}),f_{1,2}(s_{1i},s_{2i})") work as control
functions to filter the spatial trend out of the residuals, and transfer
it to the mean response in a model specification. Thus, they make it
possible to capture the shape of the spatial distribution of
![y\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bit%7D
"y_{it}"), conditional on the determinants included in the model. These
control functions also isolate stochastic spatial dependence in the
residuals, that is<span style="color: blue;">,</span> spatially
autocorrelated unobserved heterogeneity. Thus, they can be regarded as
an alternative to the use of individual regional dummies to capture
unobserved heterogeneity, as long as such heterogeneity is smoothly
distributed over space. Regional dummies peak at significantly higher
and lower levels of the mean response variable. If these peaks are
smoothly distributed over a two-dimensional surface (i.e., if unobserved
heterogeneity is spatially autocorrelated), the smooth spatial trend is
able to capture them.

Second, the smooth time trend,
![f\_{\\tau}(\\tau\_t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f_%7B%5Ctau%7D%28%5Ctau_t%29
"f_{\\tau}(\\tau_t)"), and the smooth interactions between space and
time -
![f\_{1,\\tau}(s\_{1i},\\tau\_t),f\_{2,\\tau},(s\_{2i},\\tau\_t),f\_{1,2,\\tau}(s\_{1i},s\_{2i},\\tau\_t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f_%7B1%2C%5Ctau%7D%28s_%7B1i%7D%2C%5Ctau_t%29%2Cf_%7B2%2C%5Ctau%7D%2C%28s_%7B2i%7D%2C%5Ctau_t%29%2Cf_%7B1%2C2%2C%5Ctau%7D%28s_%7B1i%7D%2Cs_%7B2i%7D%2C%5Ctau_t%29
"f_{1,\\tau}(s_{1i},\\tau_t),f_{2,\\tau},(s_{2i},\\tau_t),f_{1,2,\\tau}(s_{1i},s_{2i},\\tau_t)")
- work as control functions to capture the heterogeneous effect of
common shocks. Thus, conditional on a smooth distribution of the
spatio-temporal heterogeneity, the PS-ANOVA-SAR (SDM, SEM, SLX) model
works as an alternative to the models proposed by Bai and Li (2013), Shi
and Lee (2018), Pesaran and Tosetti (2011), Bailey, Holly, and Pesaran
(2016) and Vega and Elhorst (2016) based on extensions of common factor
models to accommodate both strong cross-sectional dependence (through
the estimation of the spatio-temporal trend) and weak cross-sectional
dependence (through the estimation of spatial autoregressive
parameters).

Furthermore, this framework is also flexible enough to control for the
linear and nonlinear functional relationships between the dependent
variable and the covariates as well as the heterogeneous effects of
these regressors across space. The model inherits all the good
properties of penalized regression splines, such as coping with missing
observations by appropriately weighting them, and straightforward
interpolation of the smooth functions.

## Computational approach to estimation

All the non-parametric terms are modeled using Penalized-splines
(P-Splines, Eilers and Marx (1996)). This methodology assumes that the
unknown function to be estimated is smooth, and can be represented as a
linear combination of basis functions,

  
![f(x)=\\sum\_j \\theta\_j
B\_j(x),](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29%3D%5Csum_j%20%5Ctheta_j%20B_j%28x%29%2C
"f(x)=\\sum_j \\theta_j B_j(x),")  
a popular election is the use of cubic B-spline basis (De Boor 1977). In
the case of multidimensional functions, the basis is calculated as
tensor products of the marginal basis functions in each dimension. The
smoohness of the curve/surface is controlled by adding to the likelihood
a penalty tunned by a smoothing parameter,
![\\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda
"\\lambda") (the number of smoothing parameters used is equal to the
dimension of the function to be estimated). The penalty can take many
forms, depending of the prior knowledge on the curve/surface. The most
common choice is to use second order differences on adjacent B-spline
coefficients (see Mı́nguez, Basile, and Durbán 2020). In the case of a
univariate function,
![f(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29
"f(x)"), the penalty is:

  
![ Penalty= \\lambda\\sum\_j \\left (
\\theta\_j-2\\theta\_{j-1}+\\theta\_{j-2}\\right
)^2.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%20Penalty%3D%20%5Clambda%5Csum_j%20%5Cleft%20%28%20%5Ctheta_j-2%5Ctheta_%7Bj-1%7D%2B%5Ctheta_%7Bj-2%7D%5Cright%20%29%5E2.
" Penalty= \\lambda\\sum_j \\left ( \\theta_j-2\\theta_{j-1}+\\theta_{j-2}\\right )^2.")  

As we mentioned above, the model is estimated via REML or ML. In order
to be able to estimate simultaneously all parameters in the model
(including the smothing parameters) we make use of the fact that a
P-spline can be reparameterized as a mixed model, and so, all the
methodology debeloped in this context can be use for estimation and
inference. Furthermore, the mixed model setting allows us to impose the
necessary constraints so that all terms in the model are identifiable,
as well as, to add spatial (or other) random effects if they are
necessary.

To speed up computations, we use a modification of the SOP (Separation
of Anisotropic Penalties) algorithm derived by Rodriguez-Alvarez et al.
(2015) (for variance components estimation). Also, the use of nested
B-spline bases (Lee, Durban, and Eilers 2013) for the interaction
components of the spatio-temporal trend contributes to the efficiency of
the fitting procedure without compromising the goodness of fit of the
model. See Mı́nguez, Basile, and Durbán (2020) for more details.

## Installation

This package is available in GitHub
(<https://github.com/rominsal/pspatreg>) and can be installed in the
usual way.\[1\]

# Basic information on the package

The package is accompanied by other two vignettes
([B\_Examples\_pspatreg\_CS\_data.html](file:///C:/Users/rober/Dropbox/spatialmodels/B_Examples_pspatreg_CS_data.html)
and
[C\_Examples\_pspatreg\_Panel\_data.html](file:///C:/Users/rober/Dropbox/spatialmodels/C_Examples_pspatreg_Panel_data.html)),
introducing the application of `pspatreg` to cross-sectional and panel
data. Here, we introduce some basic general information about the
package.

## The function `pspatfit()`

The main function in the `pspatreg` package is `pspatfit()`, which
estimates spatio-temporal pernalized spline spatial regression models
using either the Restricted Maximum Likelihood (REML) method or the
Maximum Likelihood (ML) method. In its generic form `pspatfit()` appears
as:

`pspatfit(formula, data, na.action, listw = NULL,type = "sim", method =
"eigen", Durbin = NULL, zero.policy = NULL,` `interval = NULL,trs =
NULL, cor = "none", dynamic = FALSE, control = list())`

The function `pspatfit()` returns a list of quantities of class `pspat`,
including coefficients of the parametric terms and their standard
errors, estimated coefficients corresponding to random effects in mixed
model and their standard errors, equivalent degrees of freedom,
residuals, fitted values, etc. A wide range of standard methods is also
available for the `pspat` objects, including `print()`, `summary()`,
`coef()`, `vcov()`, `anova()`, `fitted()`, `residuals()`, and `plot()`.

### The argument `formula`

The argument `formula` within the function `pspatfit()` is formula
similar to GAM specification including parametric and non-parametric
terms. Parametric covariates are included in the usual way and
non-parametric p-spline smooth terms are specified using `pspl(.)` and
`pspat(.)` for the non-parametric covariates and spatial or
spatio-temporal trends, respectively. For example

``` r
formula <- y ~ x1 + x2 + pspl(x3, nknots = 15) + pspl(x4, nknots = 20) +
                  pspt(long, lat, year, nknots = c(18,18,8),
                       psanova = TRUE, 
                       nest_sp1 = c(1, 2, 3), 
                       nest_sp2 = c(1, 2, 3),
                       nest_time = c(1, 2, 2))
```

In the example above, the model includes two parametric terms, two
nonparametric terms, and a spatio-temporal trend (with long and lat as
spatial coordinates, and year as temporal coordinate). The dimension of
the basis function both in `pspl(.)` and `pspt(.)` is defined by
`nknots`. This term should not be less than the dimension of the null
space of the penalty for the term (see `null.space.dimension` and
`choose.k` from package `mgcv` to know how to choose `nknots`). The
default number of `nknots` in `pspl(.)` is 10, but in this example we
have chosen 15 `nknots` for `g_1(x_3)` and 20 `nknots` for `g_2(x_4)`.
The default number of `nknots` in `pspt(.)` is `c(10,10,5)`, but we have
chosen `c(18,18,8)`.

In this example we also adopt an ANOVA decomposition of the
spatio-temporal trend (choosing `psanova = TRUE`). Each effect has its
own degree of smoothing allowing a greater flexibility for the
spatio-temporal trend. Calculating up to third-order interactions can be
computationally expensive. To address this problem, we can select
subgroups of interaction effects for the second- and third-order
effects. To define these subgroups, we use three parameters available in
`pspt()`: `nest_sp1`, `nest_sp2`, and `nest_time`. These parameters
indicate the divisors of the `nknots` parameters. For example, if we set
`nest_sp1 = c(1,2,3)`, we will have all knots for the `s_1` effect, 18/2
for each second-order effects with `s_1`, and 18/3 nots for the third
order effect with `s_1`.\[2\]

If we want to set any main effect to `0`, we must set the parameters
`f1_main`, `f2_main` or `ft_main` to `FALSE`, The default is `TRUE`. We
can also exclude second- or third-order effects setting `f12_int`,
`f1t_int`, `f2t_int`, `f12t_int` to `FALSE`.

### The argument `Type`

Using the argument `Type` we can choose different spatial model
specifications: `"sar"`, `"sem"`, `"sdm"`, `"sdem"`, `"sarar"`, or
`"slx"`. When creating a `"slx"`, `"sdem"` or `"sdm"` model, we need to
include the formula of the durbin part in the Durbin parameter.

### Data structure

The argument `data` must contain all the variables included in
parametric and non-parametric terms of the model. If a `pspat(.)` term
is included in `formula`, the data must contain the spatial and temporal
coordinates specified in `pspat(.)`. In this case, the coordinates must
be ordered choosing time as fast index and spatial coordinates as slow
indexes.

Both `data.frame` and `sf` class objects can be used as `data`
inputs.\[3\] `sf` objects are recommended since they allow the user to
map spatial trends. In our demos we use two datasets in `sf` version.

## Plotting smmoth terms

Plotting the estimated non-parametric smooth terms represents an
important step in semiparametric regression analyses. First, the
function `fit_terms()` computes estimated non-parametric smooth terms.
Then, the functions `plot_sp2d()` and `plot_sp3d()` are used to plot and
map spatial and spatio-temporal trends, respectively, while
`plot_sptime()` is used to plot the time trend for PS-ANOVA models in
3d; finally, and `plot_terms()` is used to plot smooth non-parametric
terms.

## Conpunting marginal effects

In the case of a semiparametric model without the spatial lag of the
dependent variable (PS model), if all regressors are manipulated
independently of the errors, ![\\widehat{g}\_\\delta(x\_{\\delta,
it})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidehat%7Bg%7D_%5Cdelta%28x_%7B%5Cdelta%2C%20it%7D%29
"\\widehat{g}_\\delta(x_{\\delta, it})") can be interpreted as the
conditional expectation of
![y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y
"y") given
![x\_{\\delta}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_%7B%5Cdelta%7D
"x_{\\delta}") (net of the effect of the other regressors). @blu.pow.03
use the term Average Structural Function (ASF) with reference to these
functions. Instead, in PS-SAR, PS-SDM or in PS-SARAR model, when
![\\rho](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Crho
"\\rho") is different from zero, the estimated smooth functions cannot
be interpreted as ASF. Taking advantage of the results obtained for
parametric SAR, we can compute the total smooth effect (total–ASF) of
![x\_{\\delta}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_%7B%5Cdelta%7D
"x_{\\delta}") as:  
  
![\\widehat{g}\_{\\delta}^{T}\\left(x\_{\\delta}\\right)=\\Sigma\_{q}&#10;\\left\[\\textbf{I}\_{n}-\\widehat{\\rho}\\textbf{W}\_{n}\\right\]^{-1}\_{ij}
b\_{\\delta q}(x\_{\\delta})\\widehat{\\beta}\_{\\delta
q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidehat%7Bg%7D_%7B%5Cdelta%7D%5E%7BT%7D%5Cleft%28x_%7B%5Cdelta%7D%5Cright%29%3D%5CSigma_%7Bq%7D%0A%5Cleft%5B%5Ctextbf%7BI%7D_%7Bn%7D-%5Cwidehat%7B%5Crho%7D%5Ctextbf%7BW%7D_%7Bn%7D%5Cright%5D%5E%7B-1%7D_%7Bij%7D%20b_%7B%5Cdelta%20q%7D%28x_%7B%5Cdelta%7D%29%5Cwidehat%7B%5Cbeta%7D_%7B%5Cdelta%20q%7D
"\\widehat{g}_{\\delta}^{T}\\left(x_{\\delta}\\right)=\\Sigma_{q}
\\left[\\textbf{I}_{n}-\\widehat{\\rho}\\textbf{W}_{n}\\right]^{-1}_{ij} b_{\\delta q}(x_{\\delta})\\widehat{\\beta}_{\\delta q}")  

where ![b\_{\\delta
q}(x\_{\\delta})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b_%7B%5Cdelta%20q%7D%28x_%7B%5Cdelta%7D%29
"b_{\\delta q}(x_{\\delta})") are P-spline basis functions, and
![\\widehat{\\beta}\_{\\delta
q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidehat%7B%5Cbeta%7D_%7B%5Cdelta%20q%7D
"\\widehat{\\beta}_{\\delta q}") the corresponding estimated parameters.

We can also compute direct and indirect (or spillover) effects of smooth
terms in the PS-SAR case as:

  
![\\widehat{g}\_{\\delta}^{D}\\left(x\_{\\delta}\\right)=\\Sigma\_{q}&#10;\\left\[\\textbf{I}\_{n}-\\widehat{\\rho}\\textbf{W}\_{n}\\right\]^{-1}\_{ii}
b\_{\\delta q}(x\_{k})\\widehat{\\beta}\_{\\delta
q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidehat%7Bg%7D_%7B%5Cdelta%7D%5E%7BD%7D%5Cleft%28x_%7B%5Cdelta%7D%5Cright%29%3D%5CSigma_%7Bq%7D%0A%5Cleft%5B%5Ctextbf%7BI%7D_%7Bn%7D-%5Cwidehat%7B%5Crho%7D%5Ctextbf%7BW%7D_%7Bn%7D%5Cright%5D%5E%7B-1%7D_%7Bii%7D%20b_%7B%5Cdelta%20q%7D%28x_%7Bk%7D%29%5Cwidehat%7B%5Cbeta%7D_%7B%5Cdelta%20q%7D
"\\widehat{g}_{\\delta}^{D}\\left(x_{\\delta}\\right)=\\Sigma_{q}
\\left[\\textbf{I}_{n}-\\widehat{\\rho}\\textbf{W}_{n}\\right]^{-1}_{ii} b_{\\delta q}(x_{k})\\widehat{\\beta}_{\\delta q}")  

  
![\\widehat{g}\_{\\delta}^{I}\\left(x\_{\\delta}\\right)=\\widehat{g}\_{\\delta}^{T}\\left(x\_{\\delta}\\right)-\\widehat{g}\_{\\delta}^{D}\\left(x\_{\\delta}\\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidehat%7Bg%7D_%7B%5Cdelta%7D%5E%7BI%7D%5Cleft%28x_%7B%5Cdelta%7D%5Cright%29%3D%5Cwidehat%7Bg%7D_%7B%5Cdelta%7D%5E%7BT%7D%5Cleft%28x_%7B%5Cdelta%7D%5Cright%29-%5Cwidehat%7Bg%7D_%7B%5Cdelta%7D%5E%7BD%7D%5Cleft%28x_%7B%5Cdelta%7D%5Cright%29
"\\widehat{g}_{\\delta}^{I}\\left(x_{\\delta}\\right)=\\widehat{g}_{\\delta}^{T}\\left(x_{\\delta}\\right)-\\widehat{g}_{\\delta}^{D}\\left(x_{\\delta}\\right)")  

Similar expressions can be provided for the direct, indirect and total
effects of the PS-SDM.

The function `impactspar()` computes direct, indirect and total impacts
for continuous parametric covariates using the standard procedure for
their computation (LeSage and Pace 2009).

The function `impactsnopar()` computes direct, indirect and total
impacts functions for continuous non-parametric covariates, while the
function `plot_impactsnopar()` is used to plot these impacts functions.
It is worth noticing that total, direct and indirect effects are never
smooth over the domain of the variable
![x\_\\delta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_%5Cdelta
"x_\\delta") due to the presence of the spatial multiplier matrix in the
algorithm for their computation. Indeed, a wiggly profile of direct,
indirect and total effects would appear even if the model were linear.
Therefore, in the spirit of the semiparametric approach, we included the
possibility of applying a spline smoother to obtain smooth curves (using
the argument `smooth=TRUE` in the function `plot_impactsnopar()`).

# References

<div id="refs" class="references">

<div id="ref-bai2015dynamic">

Bai, J., and K. Li. 2013. “Spatial Panel Data Models with Common
Shocks.” MPRA Paper 52786. University of Munich, Germany.
<https://ideas.repec.org/p/pra/mprapa/52786.html>.

</div>

<div id="ref-bailey2015two">

Bailey, Natalia, Sean Holly, and M. Hashem Pesaran. 2016. “A Two-Stage
Approach to Spatio-Temporal Analysis with Strong and Weak
Cross-Sectional Dependence.” *Journal of Applied Econometrics* 31 (1):
249–80. <https://doi.org/10.1002/jae.2468>.

</div>

<div id="ref-basdurminmonmur14">

Basile, Roberto, María Durbán, Román Mínguez, Jose María Montero, and
Jesús Mur. 2014. “Modeling Regional Economic Dynamics: Spatial
Dependence, Spatial Heterogeneity and Nonlinearities.” *Journal of
Economic Dynamics and Control* 48: 229–45.
<https://doi.org/http://dx.doi.org/10.1016/j.jedc.2014.06.011>.

</div>

<div id="ref-boor77">

De Boor, C. 1977. “Package for Calculating with B-Splines.” *Journal of
Numerical Analysis* 14: 441–72.

</div>

<div id="ref-eilers2015twenty">

Eilers, Paul HC, Brian D Marx, and Maria Durbán. 2015. “Twenty Years of
P-Splines.” *SORT-Statistics and Operations Research Transactions* 39
(2): 149–86.

</div>

<div id="ref-eil.mar.96">

Eilers, P. H., and B. D. Marx. 1996. “Flexible Smoothing with B-Splines
and Penalties.” *Statistical Science* 11: 89–121.

</div>

<div id="ref-Hoshino2018">

Hoshino, T. 2018. “Semiparametric Spatial Autoregressive Models with
Endogenous Regressors: With an Application to Crime Data.” *Journal of
Business & Economic Statistics* 36: 160–72.

</div>

<div id="ref-Lee2013">

Lee, D. J., M. Durban, and P. Eilers. 2013. “Efficient Two-Dimensional
Smoothing with P-Spline ANOVA Mixed Models and Nested Bases.”
*Computational Statistics & Data Analysis* 61: 22–37.

</div>

<div id="ref-les.pac.09">

LeSage, J., and K. Pace. 2009. *Introduction to Spatial Econometrics*.
Boca Raton: CRC Press.

</div>

<div id="ref-minguez2020alternative">

Mı́nguez, Román, Roberto Basile, and Marı́a Durbán. 2020. “An
Alternative Semiparametric Model for Spatial Panel Data.” *Statistical
Methods & Applications* 29 (4): 669–708.

</div>

<div id="ref-montero2012sar">

Montero, J, R Mı́nguez, and M Durbán. 2012. “SAR Models with
Nonparametric Spatial Trends. A P-Spline Approach.” *Estadı́stica
Española* 54 (177): 89–111.

</div>

<div id="ref-pesaran2011large">

Pesaran, M Hashem, and Elisa Tosetti. 2011. “Large Panels with Common
Factors and Spatial Correlation.” *Journal of Econometrics* 161 (2):
182–202.

</div>

<div id="ref-Rodriguez2015">

Rodriguez-Alvarez, M. X., T. Kneib, M. Durban, D. J. Lee, and P. Eilers.
2015. “Fast Smoothing Parameter Separation in Multidimensional
Generalized P-Splines: The SAP Algorithm.” *Statistics and Computing* 25
(5): 941–57.

</div>

<div id="ref-shi2018spatial">

Shi, Wei, and Lung-fei Lee. 2018. “A Spatial Panel Data Model with Time
Varying Endogenous Weights Matrices and Common Factors.” *Regional
Science and Urban Economics* 72: 6–34.

</div>

<div id="ref-vega2016regional">

Vega, Solmaria Halleck, and J Paul Elhorst. 2016. “A Regional
Unemployment Model Simultaneously Accounting for Serial Dynamics,
Spatial Dependence and Common Factors.” *Regional Science and Urban
Economics* 60: 85–95.

</div>

</div>

1.  To install any R package from GitHub you need to have previously
    installed `devtools` package from CRAN. Then execute the commands
    `library(devtools)`, to load `devtools`, and install
    `github("rominsal/pspatreg")` to install `pspatreg` package:
    `devtools::install_github("rominsal/pspatreg")`.

2.  In most empirical cases, the main effects are more flexible than
    interaction effects and therefore the number of knots in B-Spline
    bases for interaction effects do not need to be as large as the
    number of knots for the main effects (Lee, Durban, and Eilers 2013).

3.  `sf` means simple features of spatial vector objects. The geographic
    vector data model is based on points located within a coordinate
    reference system (CRS). Points can represent self-standing features
    (e.g., the location of a house) or they can be linked together to
    form more complex geometries such as lines and polygons. Most point
    geometries contain only two dimensions `x` and `y` (3-dimensional
    CRSs contain an additional `z` value, typically representing height
    above sea level). `sf` objects provide both a *geometry*
    information, describing where on Earth the feature is located, and
    *attributes* information, describing other properties (like the
    population of the region, the unemployment rate, etc.). `data.frame`
    objects store only attributes information.
