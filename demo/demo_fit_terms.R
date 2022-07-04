### demo with the whole set of examples of fit_terms() ########
###################### Examples using a panel data of rate of unemployment 
###################### in 103 Italian provinces during the period 1996-2014.
library(pspatreg)
data(unemp_it, package = "pspatreg")
lwsp_it <- spdep::mat2listw(Wsp_it)
#######  No Spatial Trend: PSAR including a spatial 
#######  lag of the dependent variable
form1 <- unrate ~ partrate + agri + cons + 
                  pspl(serv, nknots = 15) +
                  pspl(empgrowth, nknots = 20)  
gamsar <- pspatfit(form1, data = unemp_it, 
                   type = "sar", listw = lwsp_it)
summary(gamsar)

######  Fit non-parametric terms 
######  (spatial trend must be name "spttrend")
list_varnopar <- c("serv", "empgrowth")
terms_nopar <- fit_terms(gamsar, list_varnopar)

######################  Plot non-parametric terms
plot_terms(terms_nopar, unemp_it)
 
