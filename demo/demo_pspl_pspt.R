### demo with the whole set of examples of pspl() and pspt() ###
###############################################
# Examples using spatial data of Ames Houses.
###############################################
library(pspatreg)
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
#########  GAM pure with spatial trend (2d)
#####################  GAM + SAR Model
gamsar <- pspatfit(form1, data = ames_sf1, 
                   type = "sar", listw = lw_ames,
                   method = "Chebyshev")
summary(gamsar)
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
 ###############################################
###################### Examples using a panel data of rate of
###################### unemployment for 103 Italian provinces in 1996-2019.
 ###############################################
## load spatial panel and Wsp_it
## 103 Italian provinces. Period 1996-2019
data(unemp_it, package = "pspatreg")
## Wsp_it is a matrix. Create a neighboord list 
lwsp_it <- spdep::mat2listw(Wsp_it)
 ### Spatio-temporal semiparametric ANOVA model 
 ### Interaction terms f12,f1t,f2t and f12t with nested basis
 ### Remark: nest_sp1, nest_sp2 and nest_time must be divisors of nknots
 form4 <- unrate ~ partrate + agri + cons +
                   pspl(serv, nknots = 15) + 
                   pspl(empgrowth, nknots = 20) +
                   pspt(long, lat, year, 
                        nknots = c(18, 18, 8), 
                        psanova = TRUE, 
                        nest_sp1 = c(1, 2, 2), 
                        nest_sp2 = c(1, 2, 2),
                        nest_time = c(1, 2, 2))
 sptanova <- pspatfit(form4, data = unemp_it)
 summary(sptanova)
 

 ################################################  
 ### Interaction terms f1t not included in ANOVA decomposition
 form5 <- unrate ~ partrate + agri + cons +
                   pspl(serv, nknots = 15) + 
                   pspl(empgrowth, nknots=20) +
                   pspt(long, lat, year, 
                        nknots = c(18, 18, 8),
                        psanova = TRUE, 
                        nest_sp1 = c(1, 2, 3), 
                        nest_sp2 = c(1, 2, 3),
                        nest_time = c(1, 2, 2), 
                        f1t_int = FALSE)
 ## Add sar specification and ar1 temporal correlation 
 sptanova2_sar_ar1 <- pspatfit(form5, data = unemp_it, 
                              listw = lwsp_it, 
                              type = "sar",
                              cor = "ar1")
summary(sptanova2_sar_ar1)                                 

