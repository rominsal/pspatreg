# pspatreg
Estimate **geoadditive spatial or spatio-temporal semiparametric regression models**
from an econometric perspective.
The type of econometric models included are PS-SAR, PS-SEM, PS-SARAR, PS-SLX or PS-Durbin
specifications.

These specifications could include parametric and non-parametric covariates, spatial or 
spatio-temporal non-parametric trends and spatial lags of the dependent variable and/or 
the noise of the model. 

The non-parametric terms (either trends or covariates) are modeled using **P-Splines**. 
The non-parametric trend can be decomposed in an ANOVA way including main and interactions 
effects of 2nd and 3rd order (for spatio-temporal trends). 

The estimation method can be Restricted Maximum Likelihood (REML) or Maximum Likelihood (ML).

There are specific functions to compute:

- Non-parametric and parametric impacts (including inference about them).
- Plots of non-parametric impacts and non-parametric terms. 
- Plots of the spatial or spatio-temporal trends (including main and interaction terms for ANOVA specifications).
- Usual methods for summary, anova, print, cov, residuals, ...

