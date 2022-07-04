### demo with the whole set of examples of pspatfit() ########
##########################
library(pspatreg)
###############################################
# Examples using spatial data of Ames Houses.
###############################################
# Getting and preparing the data
library(spdep)
library(sf)
ames <- AmesHousing::make_ames() # Raw Ames Housing Data
ames_sf <- st_as_sf(ames, coords = c("Longitude", "Latitude"))
ames_sf$Longitude <- ames$Longitude
ames_sf$Latitude <- ames$Latitude
ames_sf$lnSale_Price <- log(ames_sf$Sale_Price)
ames_sf$lnLot_Area <- log(ames_sf$Lot_Area)
ames_sf$lnTotal_Bsmt_SF <- log(ames_sf$Total_Bsmt_SF+1)
ames_sf$lnGr_Liv_Area <- log(ames_sf$Gr_Liv_Area)
########### Constructing the spatial weights matrix
ames_sf1 <- ames_sf[(duplicated(ames_sf$Longitude) == FALSE), ]
coord_sf1 <- cbind(ames_sf1$Longitude, ames_sf1$Latitude)
ID <- row.names(as(ames_sf1, "sf"))
col_tri_nb <- tri2nb(coord_sf1)
soi_nb <- graph2nb(soi.graph(col_tri_nb, 
                            coord_sf1), 
                   row.names = ID)
lw_ames <- nb2listw(soi_nb, style = "W", 
                    zero.policy = FALSE)

####  GAM pure with pspatreg
form1 <- lnSale_Price ~ Fireplaces + Garage_Cars +
          pspl(lnLot_Area, nknots = 20) + 
          pspl(lnTotal_Bsmt_SF, nknots = 20) +
          pspl(lnGr_Liv_Area, nknots = 20)    
gampure <- pspatfit(form1, data = ames_sf1)
summary(gampure)

######################  Get Non-parametric terms of GAM with pspatreg
list_varnopar <- c("lnLot_Area", "lnTotal_Bsmt_SF", 
"lnGr_Liv_Area")
terms_nopar <- fit_terms(gampure, list_varnopar, intercept = TRUE)
######################  Plot non-parametric terms
plot_terms(terms_nopar, ames_sf1)

#####################  GAM + SAR Model
gamsar <- pspatfit(form1, data = ames_sf1, 
                   type = "sar", listw = lw_ames,
                   method = "Chebyshev")
summary(gamsar)
######### Non-Parametric Total, Direct and Indirect impacts
### with impactsnopar(viewplot = TRUE)
nparimpacts <- impactsnopar(gamsar, 
                            listw = lw_ames, 
                            viewplot = TRUE)
############ Non-Parametric Total, Direct and Indirect impacts
### with impactsnopar(viewplot = FALSE) and using plot_impactsnopar()
nparimpacts <- impactsnopar(gamsar, listw = lw_ames, viewplot = FALSE)
plot_impactsnopar(nparimpacts, data = ames_sf1, smooth = TRUE)

###################### Parametric Total, Direct and Indirect impacts
parimpacts <- impactspar(gamsar, listw = lw_ames)
summary(parimpacts)

###############################################
### Models with 2d spatial trend
form2 <- lnSale_Price ~ Fireplaces + Garage_Cars +
          pspl(lnLot_Area, nknots = 20) + 
          pspl(lnTotal_Bsmt_SF, nknots = 20) +
          pspl(lnGr_Liv_Area, nknots = 20) +
          pspt(Longitude, Latitude, 
               nknots = c(10, 10), 
               psanova = FALSE)
#####################  GAM + GEO Model
gamgeo2d <- pspatfit(form2, data = ames_sf1)
summary(gamgeo2d)

gamgeo2dsar <- pspatfit(form2, data = ames_sf1,
                        type = "sar", 
                        listw = lw_ames, 
                        method = "Chebyshev")
summary(gamgeo2dsar)
####### plot spatial trend for spatial point coordinate
plot_sp2d(gamgeo2dsar, data = ames_sf1)
### Models with psanova 2d spatial trend
form3 <- lnSale_Price ~ Fireplaces + Garage_Cars +
          pspl(lnLot_Area, nknots = 20) + 
          pspl(lnTotal_Bsmt_SF, nknots = 20) +
          pspl(lnGr_Liv_Area, nknots = 20) +
          pspt(Longitude, Latitude, 
               nknots = c(10, 10), 
               psanova = TRUE)
gamgeo2danovasar <- pspatfit(form3, data = ames_sf1,
                        type = "sar", 
                        listw = lw_ames, method = "Chebyshev")
summary(gamgeo2danovasar)
####### plot spatial trend for spatial point coordinate
plot_sp2d(gamgeo2danovasar, data = ames_sf1, 
addmain = TRUE, addint = TRUE)

## Comparison between models
anova(gampure, gamsar, gamgeo2d, gamgeo2dsar,
gamgeo2danovasar, lrtest = FALSE)

###############################################
###################### Examples using a panel data of rate of
###################### unemployment for 103 Italian provinces in 1996-2019.
###############################################
## load spatial panel and Wsp_it
## 103 Italian provinces. Period 1996-2019
data(unemp_it, package = "pspatreg")
## Wsp_it is a matrix. Create a neighboord list 
lwsp_it <- spdep::mat2listw(Wsp_it)
 ### Models with spatio-temporal trend
 ### Spatio-temporal semiparametric ANOVA model without spatial lag
 ### Interaction terms f12,f1t,f2t and f12t with nested basis
 ### Remark: nest_sp1, nest_sp2 and nest_time must be divisors of nknots
 form4 <- unrate ~ partrate + agri + cons +
                   pspl(serv, nknots = 15) + 
                   pspl(empgrowth, nknots = 20) +
                   pspt(long, lat, year, 
                        nknots = c(18, 18, 12),
                        psanova = TRUE, 
                        nest_sp1 = c(1, 2, 3), 
                        nest_sp2 = c(1, 2, 3),
                        nest_time = c(1, 2, 2))
 sptanova <- pspatfit(form4, data = unemp_it)
 summary(sptanova)
 
### Create sf object to make the plot 
### of spatio-temporal trends
library(sf)
unemp_it_sf <- st_as_sf(dplyr::left_join(
                              unemp_it, 
                              map_it,  
                        by = c("prov" = "COD_PRO")))
###### Plot spatio-temporal trends for different years
plot_sp3d(sptanova, data = unemp_it_sf, 
          time_var = "year", 
          time_index = c(1996, 2005, 2019),
          addmain = FALSE, addint = FALSE)
###### Plot of spatio-temporal trend, main effects 
######      and interaction effect for a year
plot_sp3d(sptanova, data = unemp_it_sf, 
          time_var = "year", 
          time_index = c(2019),
          addmain = TRUE, addint = TRUE)
          
###### Plot of temporal trends for each province
plot_sptime(sptanova, 
            data = unemp_it, 
            time_var = "year", 
            reg_var = "prov")

 
 ###############################################
 ### Spatio-temporal semiparametric ANOVA model without spatial lag
 ### Now we repeat previous spatio-temporal model but 
 ### restricting some interactions
 ### Interaction terms f12,f1t and f12t with nested basis
 ### Interaction term f2t restricted to 0
 
  form5 <- unrate ~ partrate + agri + cons + empgrowth +
                  pspl(serv, nknots = 15) + 
                  pspt(long, lat, year, 
                       nknots = c(18, 18, 6), 
                       psanova = TRUE,
                       nest_sp1 = c(1, 2, 3), 
                       nest_sp2 = c(1, 2, 3),
                       nest_time = c(1, 2, 2), 
                       f2t_int = FALSE)
 ## Add sar specification and ar1 temporal correlation 
 sptanova2_sar_ar1 <- pspatfit(form5, data = unemp_it, 
                              listw = lwsp_it, 
                              type = "sar",
                              cor = "ar1")
summary(sptanova2_sar_ar1)                     
################ Comparison with parametric panels            
######################  Demeaning (Within Estimators)
formpar <- unrate ~ partrate + agri + cons
# Not demeaning model
param <- pspatfit(formpar, data = unemp_it, listw = lwsp_it)
summary(param)
# Demeaning model
param_dem <- pspatfit(formpar, data = unemp_it,
                      demean = TRUE,
                      index = c("prov", "year"),
                      eff_demean = "individual" )
summary(param_dem)
# Compare results with plm package
param_plm <- plm::plm(formula = formpar,
                      data = unemp_it,
                      index = c("prov", "year"),
                      effect = "individual",
                      model = "within")
summary(param_plm)                                              
param_dem_time <- pspatfit(formpar, 
                      data = unemp_it, 
                      listw = lwsp_it,
                      demean = TRUE,
                      eff_demean = "time",
                      index = c("prov", "year"))
summary(param_dem_time)
param_plm_time <- plm::plm(formula = formpar,
                      data = unemp_it,
                      index = c("prov", "year"),
                      effect = "time",
                      model = "within")
summary(param_plm_time)
param_dem_twoways <- pspatfit(formpar, data = unemp_it,
                      demean = TRUE,
                      eff_demean = "twoways",
                      index = c("prov", "year") )
summary(param_dem_twoways)
param_plm_twoways <- plm::plm(formula = formpar,
                      data = unemp_it,
                      index = c("prov", "year"),
                      effect = "twoways",
                      model = "within")
summary(param_plm_twoways)
##### Demeaning with nonparametric covariates
formgam <- unrate ~ partrate + agri + cons +  
                    pspl(serv, nknots = 15) +
                    pspl(empgrowth, nknots = 20)
                    
gam_dem <- pspatfit(formula = formgam,
                      data = unemp_it,
                      demean = TRUE,
                      index = c("prov", "year"))
summary(gam_dem)   
# Compare with GAM pure without demeaning                    
gam <- pspatfit(formula = formgam,
                 data = unemp_it)
summary(gam)

## Demeaning with type = "sar" model
gamsar_dem <- pspatfit(formula = formgam,
                      data = unemp_it,
                      type = "sar", 
                      listw = lwsp_it,
                      demean = TRUE,
                      index = c("prov", "year"))
summary(gamsar_dem)
