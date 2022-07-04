### demo with the whole set of examples of impactspar() ########
################################################
#### Examples using a panel data of rate of
##### unemployment for 103 Italian provinces in period 1996-2014.
library(pspatreg)
data(unemp_it, package = "pspatreg")
## Wsp_it is a matrix. Create a neighboord list 
lwsp_it <- spdep::mat2listw(Wsp_it)
## short sample for spatial pure case (2d)
########  No Spatial Trend: PSAR including a spatial 
########  lag of the dependent variable
form1 <- unrate ~ partrate + agri + cons + empgrowth +
                 pspl(serv, nknots = 15) 
### example with type = "sar"                   
gamsar <- pspatfit(form1, 
                   data = unemp_it, 
                   type = "sar", 
                   listw = lwsp_it)
summary(gamsar)
###### Parametric Total, Direct and Indirect Effects
imp_parvar <- impactspar(gamsar, listw = lwsp_it)
summary(imp_parvar)

### example with type = "slx"                   

gamslx <- pspatfit(form1, 
                   data = unemp_it, 
                   type = "slx", 
                   listw = lwsp_it)
                   
summary(gamslx)
###### Parametric Total, Direct and Indirect Effects
imp_parvarslx <- impactspar(gamslx, listw = lwsp_it)
summary(imp_parvarslx)

