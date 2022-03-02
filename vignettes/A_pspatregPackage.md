An introduction to pspatreg
================
Román Mínguez (UCLM), Roberto Basile (UNIVAQ), María Durbán (UC3M)
2022-03-02

# 1 Introduction

**pspatreg** is a package that fits penalized spline (PS) semiparametric
static and dynamic spatial autoregressive models via Restricted (or
Residual) Maximum Likelihood (REML) and Maximum Likelihood (ML). This
approach combines penalized regression spline methods (Paul HC Eilers,
Marx, and Durbán 2015) with standard spatial autoregressive models (such
as SAR, SEM and SDM). These types of models are thoroughly discussed in
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

## 1.1 The semiparametric spatial autoregressive model

In a very general form, the semiparametric spatial autoregressive model
reads as:

<img src="https://render.githubusercontent.com/render/math?math= 
$$y_{it}=\alpha y_{it-1}+\rho \sum_{j=1}^N w_{ij,N} y_{jt} + \pi \sum_{j=1}^N w_{ij,N} y_{jt-1}+ \\
\sum_{k=1}^K \beta^*_k x^*_{k,it} + \sum_{k=1}^K \gamma^*_k \sum_{j=1}^N w_{ij,N} x^*_{k,jt}+
\sum_{\delta=1}^\Delta g_\delta(x_{\delta, it}) + \sum_{\delta=1}^\Delta h_\delta\left(\sum_{j=1}^N w_{ij,N} x_{\delta,jt}\right) + \\
\widetilde{ f}(s_{1i},s_{2i},\tau_t)+\epsilon_{it}$$
">
$$\\epsilon\_{it}=\\theta \\sum\_{j=1}^N w\_{ij,N}\\epsilon\_{jt}+\\phi \\epsilon\_{it-1}+u\_{it}$$
*u*<sub>*i**t*</sub> ∼ *i*.*i*.*d*.(0,*σ*<sub>*u*</sub><sup>2</sup>)

where
(*y*<sub>*i**t*</sub>,*x*<sub>1, *i**t*</sub><sup>\*</sup>,...,*x*<sub>*K*, *i**t*</sub><sup>\*</sup>,*x*<sub>1, *i**t*</sub>,...,,*x*<sub>*Δ*, *i**t*</sub>)
are data observed for a sample of spatial units *i* (*i* = 1, …, *N*) in
each time period *t* = 1, …, *T*. The terms
(*s*<sub>1*i*</sub>,*s*<sub>2*i*</sub>) indicate the spatial coordinates
(latitude and longitude) of the spatial unit *i* (either a spatial point
or the centroid of a spatial polygon), and *W*<sub>*i**j*</sub> is an
element of a row-standardized spatial weights matrix. The terms included
into the model are:

-   *y*<sub>*i**t* − 1</sub> = the time lag of *y*<sub>*i**t*</sub>;

-   $\\sum\_{j=1}^N w\_{ij,N} y\_{it}$ = the contemporaneous spatial lag
    of *y*<sub>*i**t*</sub>;

-   $\\sum\_{j=1}^N w\_{ij,N} y\_{it-1}$ = the time lag of the spatial
    lag of *y*<sub>*i**t*</sub>;

-   *α*, *ρ*, *π*, *β*<sub>*k*</sub><sup>\*</sup>, and
    *γ*<sub>*k*</sub><sup>\*</sup> = fixed parameters;

-   $\\sum\_{k=1}^K \\beta^\*\_k x^\*\_{k,it}$ and
    $\\sum\_{k=1}^K \\gamma^\*\_k \\sum\_{j=1}^N w\_{ij,N} x^\*\_{k,jt}$
    = parametric linear terms of some covariates
    *x*<sub>*k*, *i**t*</sub> and of their spatial lags;

-   *g*<sub>*δ*</sub>(.) and *h*<sub>*γ*</sub>(.) = nonparametric smooth
    functions of other covariates and of their spatial lags (they can
    also accommodate varying coefficient terms, smooth interaction
    between covariates, factor-by-curve intercations, and so on);

-   *f̃*(*s*<sub>1*i*</sub>,*s*<sub>2*i*</sub>,*τ*<sub>*t*</sub>) = an
    unknown nonparametric spatio-temporal trend;

-   *ϵ*<sub>*i**t*</sub> = an idiosyncratic error term that can follow a
    spatial autoregressive process
    $\\epsilon\_{it}=\\theta \\sum\_{j=1}^N w\_{ij,N}\\epsilon\_{it}+u\_{it}$
    with *u*<sub>*i**t*</sub> ∼ *N*(0,*σ*<sup>2</sup>) (SEM), a time
    series autoregressive AR(1) process, i.e.,
    *ϵ*<sub>*i**t*</sub> = *ϕ**ϵ*<sub>*i**t* − 1</sub> + *u*<sub>*i**t*</sub>
    with *u*<sub>*i**t*</sub> ∼ *N*(0,*σ*<sup>2</sup>), or both.

Obviously, this very general specification is hardly suitable in real
data applications. However, it is worth noticing that it nests most of
the spatial additive models which can be used in practice. For example,
a more suitable semiparametric static model reads as:

$$y\_{it}=\\rho \\sum\_{j=1}^N w\_{ij,N} y\_{jt} + 
\\sum\_{k=1}^K \\beta^\*\_k x^\*\_{k,it} + 
\\sum\_{\\delta=1}^\\Delta g\_\\delta(x\_{\\delta, it}) +
\\widetilde{ f}(s\_{1i},s\_{2i},\\tau_t)+\\epsilon\_{it}$$

*ϵ*<sub>*i**t*</sub> = *ϕ**ϵ*<sub>*i**t* − 1</sub> + *u*<sub>*i**t*</sub>
*u*<sub>*i**t*</sub> ∼ *i*.*i*.*d*.(0,*σ*<sub>*u*</sub><sup>2</sup>)
This semiparametric SAR model turns out to be extremely useful to
capture interactive spatial and temporal unobserved heterogeneity when
the last one is smoothly distributed over space and time (Mı́nguez,
Basile, and Durbán 2020). The dynamic extension (including
*y*<sub>*i**t* − 1</sub> and $\\sum\_{j=1}^N w\_{ij,N} y\_{it-1}$) is
also very promising and merits further theoretical investigation.
Finally, the following semiparametric SAR model is very useful for
modeling cross-setional spatial data taking into account nonlinearities,
spatial dependence and spatial heterogeneity:

$$y\_{i}=\\rho \\sum\_{j=1}^N w\_{ij,N} y\_{j} + \\sum\_{k=1}^K \\beta^\*\_k x^\*\_{k,i} + \\sum\_{\\delta=1}^\\Delta g\_\\delta(x\_{\\delta, i}) + 
\\widetilde{ f}(s\_{1i},s\_{2i})+\\epsilon\_{i}$$

*ϵ*<sub>*i*</sub> ∼ *i*.*i*.*d*.(0,*σ*<sub>*ϵ*</sub><sup>2</sup>)
In many situations the spatio-temporal trend to be estimated can be
complex, and the use of a single multidimensional smooth function may
not be flexible enough to capture the structure in the data. To solve
this problem, an ANOVA-type decomposition of
*f̃*(*s*<sub>1*i*</sub>,*s*<sub>2*i*</sub>,*τ*<sub>*t*</sub>) can be
used, where spatial and temporal main effects, and second- and
third-order interactions between them can be identified:

$$\\widetilde{ f}(s\_{1i},s\_{2i},\\tau_t)=f_1(s\_{1i})+f_2(s\_{2i})+f\_{\\tau}(\\tau_t)+ \\\\ f\_{1,2}(s\_{1i},s\_{2i})+f\_{1,\\tau}(s\_{1i},\\tau_t)+f\_{2,\\tau}+(s\_{2i},\\tau_t)+f\_{1,2,\\tau}(s\_{1i},s\_{2i},\\tau_t)$$

First, the geoadditive terms given by
*f*<sub>1</sub>(*s*<sub>1*i*</sub>), *f*<sub>2</sub>(*s*<sub>2*i*</sub>), *f*<sub>1, 2</sub>(*s*<sub>1*i*</sub>,*s*<sub>2*i*</sub>)
work as control functions to filter the spatial trend out of the
residuals, and transfer it to the mean response in a model
specification. Thus, they make it possible to capture the shape of the
spatial distribution of *y*<sub>*i**t*</sub>, conditional on the
determinants included in the model. These control functions also isolate
stochastic spatial dependence in the residuals, that is<span
style="color: blue;">,</span> spatially autocorrelated unobserved
heterogeneity. Thus, they can be regarded as an alternative to the use
of individual regional dummies to capture unobserved heterogeneity, as
long as such heterogeneity is smoothly distributed over space. Regional
dummies peak at significantly higher and lower levels of the mean
response variable. If these peaks are smoothly distributed over a
two-dimensional surface (i.e., if unobserved heterogeneity is spatially
autocorrelated), the smooth spatial trend is able to capture them.

Second, the smooth time trend, *f*<sub>*τ*</sub>(*τ*<sub>*t*</sub>), and
the smooth interactions between space and time -
*f*<sub>1, *τ*</sub>(*s*<sub>1*i*</sub>,*τ*<sub>*t*</sub>), *f*<sub>2, *τ*</sub>, (*s*<sub>2*i*</sub>,*τ*<sub>*t*</sub>), *f*<sub>1, 2, *τ*</sub>(*s*<sub>1*i*</sub>,*s*<sub>2*i*</sub>,*τ*<sub>*t*</sub>)
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

## 1.2 Computational approach to estimation

All the non-parametric terms are modeled using Penalized-splines
(P-Splines, P. H. Eilers and Marx (1996)). This methodology assumes that
the unknown function to be estimated is smooth, and can be represented
as a linear combination of basis functions,

*f*(*x*) = ∑<sub>*j*</sub>*θ*<sub>*j*</sub>*B*<sub>*j*</sub>(*x*),
a popular election is the use of cubic B-spline basis (De Boor 1977). In
the case of multidimensional functions, the basis is calculated as
tensor products of the marginal basis functions in each dimension. The
smoohness of the curve/surface is controlled by adding to the likelihood
a penalty tunned by a smoothing parameter, *λ* (the number of smoothing
parameters used is equal to the dimension of the function to be
estimated). The penalty can take many forms, depending of the prior
knowledge on the curve/surface. The most common choice is to use second
order differences on adjacent B-spline coefficients (see Mı́nguez,
Basile, and Durbán 2020). In the case of a univariate function,
*f*(*x*), the penalty is:

*P**e**n**a**l**t**y* = *λ*∑<sub>*j*</sub>(*θ*<sub>*j*</sub>−2*θ*<sub>*j* − 1</sub>+*θ*<sub>*j* − 2</sub>)<sup>2</sup>.

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

## 1.3 Installation

This package is available in GitHub
(<https://github.com/rominsal/pspatreg>) and can be installed in the
usual way.[1]

# 2 Basic information on the package

The package is accompanied by other two vignettes
([B_Examples_pspatreg_CS_data.html](file:///C:/Users/rober/Dropbox/spatialmodels/B_Examples_pspatreg_CS_data.html)
and
[C_Examples_pspatreg_Panel_data.html](file:///C:/Users/rober/Dropbox/spatialmodels/C_Examples_pspatreg_Panel_data.html)),
introducing the application of `pspatreg` to cross-sectional and panel
data. Here, we introduce some basic general information about the
package.

## 2.1 The function `pspatfit()`

The main function in the `pspatreg` package is `pspatfit()`, which
estimates spatio-temporal pernalized spline spatial regression models
using either the Restricted Maximum Likelihood (REML) method or the
Maximum Likelihood (ML) method. In its generic form `pspatfit()` appears
as:

`pspatfit(formula, data, na.action, listw = NULL,type = "sim", method = "eigen", Durbin = NULL, zero.policy = NULL,`
`interval = NULL,trs = NULL, cor = "none", dynamic = FALSE, control = list())`

The function `pspatfit()` returns a list of quantities of class `pspat`,
including coefficients of the parametric terms and their standard
errors, estimated coefficients corresponding to random effects in mixed
model and their standard errors, equivalent degrees of freedom,
residuals, fitted values, etc. A wide range of standard methods is also
available for the `pspat` objects, including `print()`, `summary()`,
`coef()`, `vcov()`, `anova()`, `fitted()`, `residuals()`, and `plot()`.

### 2.1.1 The argument `formula`

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
order effect with `s_1`.[2]

If we want to set any main effect to `0`, we must set the parameters
`f1_main`, `f2_main` or `ft_main` to `FALSE`, The default is `TRUE`. We
can also exclude second- or third-order effects setting `f12_int`,
`f1t_int`, `f2t_int`, `f12t_int` to `FALSE`.

### 2.1.2 The argument `Type`

Using the argument `Type` we can choose different spatial model
specifications: `"sar"`, `"sem"`, `"sdm"`, `"sdem"`, `"sarar"`, or
`"slx"`. When creating a `"slx"`, `"sdem"` or `"sdm"` model, we need to
include the formula of the durbin part in the Durbin parameter.

### 2.1.3 Data structure

The argument `data` must contain all the variables included in
parametric and non-parametric terms of the model. If a `pspat(.)` term
is included in `formula`, the data must contain the spatial and temporal
coordinates specified in `pspat(.)`. In this case, the coordinates must
be ordered choosing time as fast index and spatial coordinates as slow
indexes.

Both `data.frame` and `sf` class objects can be used as `data`
inputs.[3] `sf` objects are recommended since they allow the user to map
spatial trends. In our demos we use two datasets in `sf` version.

## 2.2 Plotting smmoth terms

Plotting the estimated non-parametric smooth terms represents an
important step in semiparametric regression analyses. First, the
function `fit_terms()` computes estimated non-parametric smooth terms.
Then, the functions `plot_sp2d()` and `plot_sp3d()` are used to plot and
map spatial and spatio-temporal trends, respectively, while
`plot_sptime()` is used to plot the time trend for PS-ANOVA models in
3d; finally, and `plot_terms()` is used to plot smooth non-parametric
terms.

## 2.3 Conpunting marginal effects

In the case of a semiparametric model without the spatial lag of the
dependent variable (PS model), if all regressors are manipulated
independently of the errors,
*ĝ*<sub>*δ*</sub>(*x*<sub>*δ*, *i**t*</sub>) can be interpreted as the
conditional expectation of *y* given *x*<sub>*δ*</sub> (net of the
effect of the other regressors). @blu.pow.03 use the term Average
Structural Function (ASF) with reference to these functions. Instead, in
PS-SAR, PS-SDM or in PS-SARAR model, when *ρ* is different from zero,
the estimated smooth functions cannot be interpreted as ASF. Taking
advantage of the results obtained for parametric SAR, we can compute the
total smooth effect (total–ASF) of *x*<sub>*δ*</sub> as:  
*ĝ*<sub>*δ*</sub><sup>*T*</sup>(*x*<sub>*δ*</sub>) = *Σ*<sub>*q*</sub>\[**I**<sub>*n*</sub>−*ρ̂***W**<sub>*n*</sub>\]<sub>*i**j*</sub><sup>−1</sup>*b*<sub>*δ**q*</sub>(*x*<sub>*δ*</sub>)*β̂*<sub>*δ**q*</sub>

where *b*<sub>*δ**q*</sub>(*x*<sub>*δ*</sub>) are P-spline basis
functions, and *β̂*<sub>*δ**q*</sub> the corresponding estimated
parameters.

We can also compute direct and indirect (or spillover) effects of smooth
terms in the PS-SAR case as:

*ĝ*<sub>*δ*</sub><sup>*D*</sup>(*x*<sub>*δ*</sub>) = *Σ*<sub>*q*</sub>\[**I**<sub>*n*</sub>−*ρ̂***W**<sub>*n*</sub>\]<sub>*i**i*</sub><sup>−1</sup>*b*<sub>*δ**q*</sub>(*x*<sub>*k*</sub>)*β̂*<sub>*δ**q*</sub>

*ĝ*<sub>*δ*</sub><sup>*I*</sup>(*x*<sub>*δ*</sub>) = *ĝ*<sub>*δ*</sub><sup>*T*</sup>(*x*<sub>*δ*</sub>) − *ĝ*<sub>*δ*</sub><sup>*D*</sup>(*x*<sub>*δ*</sub>)

Similar expressions can be provided for the direct, indirect and total
effects of the PS-SDM.

The function `impactspar()` computes direct, indirect and total impacts
for continuous parametric covariates using the standard procedure for
their computation (LeSage and Pace 2009).

The function `impactsnopar()` computes direct, indirect and total
impacts functions for continuous non-parametric covariates, while the
function `plot_impactsnopar()` is used to plot these impacts functions.
It is worth noticing that total, direct and indirect effects are never
smooth over the domain of the variable *x*<sub>*δ*</sub> due to the
presence of the spatial multiplier matrix in the algorithm for their
computation. Indeed, a wiggly profile of direct, indirect and total
effects would appear even if the model were linear. Therefore, in the
spirit of the semiparametric approach, we included the possibility of
applying a spline smoother to obtain smooth curves (using the argument
`smooth=TRUE` in the function `plot_impactsnopar()`).

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-bai2015dynamic" class="csl-entry">

Bai, J., and K. Li. 2013. “Spatial Panel Data Models with Common
Shocks.” MPRA Paper 52786. University of Munich, Germany.
<https://ideas.repec.org/p/pra/mprapa/52786.html>.

</div>

<div id="ref-bailey2015two" class="csl-entry">

Bailey, Natalia, Sean Holly, and M. Hashem Pesaran. 2016. “A Two-Stage
Approach to Spatio-Temporal Analysis with Strong and Weak
Cross-Sectional Dependence.” *Journal of Applied Econometrics* 31 (1):
249–80. <https://doi.org/10.1002/jae.2468>.

</div>

<div id="ref-basdurminmonmur14" class="csl-entry">

Basile, Roberto, María Durbán, Román Mínguez, Jose María Montero, and
Jesús Mur. 2014. “Modeling Regional Economic Dynamics: Spatial
Dependence, Spatial Heterogeneity and Nonlinearities.” *Journal of
Economic Dynamics and Control* 48: 229–45.
https://doi.org/<http://dx.doi.org/10.1016/j.jedc.2014.06.011>.

</div>

<div id="ref-boor77" class="csl-entry">

De Boor, C. 1977. “Package for Calculating with B-Splines.” *Journal of
Numerical Analysis* 14: 441–72.

</div>

<div id="ref-eil.mar.96" class="csl-entry">

Eilers, P. H., and B. D. Marx. 1996. “Flexible Smoothing with B-Splines
and Penalties.” *Statistical Science* 11: 89–121.

</div>

<div id="ref-eilers2015twenty" class="csl-entry">

Eilers, Paul HC, Brian D Marx, and Maria Durbán. 2015. “Twenty Years of
p-Splines.” *SORT-Statistics and Operations Research Transactions* 39
(2): 149–86.

</div>

<div id="ref-Hoshino2018" class="csl-entry">

Hoshino, T. 2018. “Semiparametric Spatial Autoregressive Models with
Endogenous Regressors: With an Application to Crime Data.” *Journal of
Business & Economic Statistics* 36: 160–72.

</div>

<div id="ref-Lee2013" class="csl-entry">

Lee, D. J., M. Durban, and P. Eilers. 2013. “Efficient Two-Dimensional
Smoothing with P-Spline ANOVA Mixed Models and Nested Bases.”
*Computational Statistics & Data Analysis* 61: 22–37.

</div>

<div id="ref-les.pac.09" class="csl-entry">

LeSage, J., and K. Pace. 2009. *Introduction to Spatial Econometrics*.
Boca Raton: CRC Press.

</div>

<div id="ref-minguez2020alternative" class="csl-entry">

Mı́nguez, Román, Roberto Basile, and Marı́a Durbán. 2020. “An Alternative
Semiparametric Model for Spatial Panel Data.” *Statistical Methods &
Applications* 29 (4): 669–708.

</div>

<div id="ref-montero2012sar" class="csl-entry">

Montero, J, R Mı́nguez, and M Durbán. 2012. “SAR Models with
Nonparametric Spatial Trends. A P-Spline Approach.” *Estadı́stica
Española* 54 (177): 89–111.

</div>

<div id="ref-pesaran2011large" class="csl-entry">

Pesaran, M Hashem, and Elisa Tosetti. 2011. “Large Panels with Common
Factors and Spatial Correlation.” *Journal of Econometrics* 161 (2):
182–202.

</div>

<div id="ref-Rodriguez2015" class="csl-entry">

Rodriguez-Alvarez, M. X., T. Kneib, M. Durban, D. J. Lee, and P. Eilers.
2015. “Fast Smoothing Parameter Separation in Multidimensional
Generalized P-Splines: The SAP Algorithm.” *Statistics and Computing* 25
(5): 941–57.

</div>

<div id="ref-shi2018spatial" class="csl-entry">

Shi, Wei, and Lung-fei Lee. 2018. “A Spatial Panel Data Model with Time
Varying Endogenous Weights Matrices and Common Factors.” *Regional
Science and Urban Economics* 72: 6–34.

</div>

<div id="ref-vega2016regional" class="csl-entry">

Vega, Solmaria Halleck, and J Paul Elhorst. 2016. “A Regional
Unemployment Model Simultaneously Accounting for Serial Dynamics,
Spatial Dependence and Common Factors.” *Regional Science and Urban
Economics* 60: 85–95.

</div>

</div>

[1] To install any R package from GitHub you need to have previously
installed `devtools` package from CRAN. Then execute the commands
`library(devtools)`, to load `devtools`, and install
`github("rominsal/pspatreg")` to install `pspatreg` package:
`devtools::install_github("rominsal/pspatreg")`.

[2] In most empirical cases, the main effects are more flexible than
interaction effects and therefore the number of knots in B-Spline bases
for interaction effects do not need to be as large as the number of
knots for the main effects (Lee, Durban, and Eilers 2013).

[3] `sf` means simple features of spatial vector objects. The geographic
vector data model is based on points located within a coordinate
reference system (CRS). Points can represent self-standing features
(e.g., the location of a house) or they can be linked together to form
more complex geometries such as lines and polygons. Most point
geometries contain only two dimensions `x` and `y` (3-dimensional CRSs
contain an additional `z` value, typically representing height above sea
level). `sf` objects provide both a *geometry* information, describing
where on Earth the feature is located, and *attributes* information,
describing other properties (like the population of the region, the
unemployment rate, etc.). `data.frame` objects store only attributes
information.
