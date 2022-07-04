### demo with the whole set of examples of plot_impactsnopar() ########
################################################
# Examples using spatial data of Ames Houses.
###############################################
# Getting and preparing the data
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
form1 <- lnSale_Price ~ Fireplaces + Garage_Cars +
          pspl(lnLot_Area, nknots = 20) + 
          pspl(lnTotal_Bsmt_SF, nknots = 20) +
          pspl(lnGr_Liv_Area, nknots = 20)    

gamsar <- pspatfit(form1, data = ames_sf1, 
                   type = "sar", listw = lw_ames,
                   method = "Chebyshev")
summary(gamsar)
nparimpacts <- impactsnopar(gamsar, listw = lw_ames, viewplot = FALSE)
plot_impactsnopar(nparimpacts, data = ames_sf1, smooth = TRUE)
###### Examples using a panel data of rate of
###### unemployment for 103 Italian provinces in period 1996-2014.
library(pspatreg)
data(unemp_it)
## Wsp_it is a matrix. Create a neighboord list 
lwsp_it <- spdep::mat2listw(Wsp_it)
## short sample for spatial pure case (2d)
########  No Spatial Trend: PSAR including a spatial 
########  lag of the dependent variable
form1 <- unrate ~ partrate + agri + cons + empgrowth +
                  pspl(serv, nknots = 15) 
gamsar <- pspatfit(form1, data = unemp_it, 
                   type = "sar", 
                   listw = lwsp_it)
summary(gamsar)
###### Non-Parametric Total, Direct and Indirect impacts
imp_nparvar <- impactsnopar(gamsar, alpha = 0.05,
                            listw = lwsp_it, 
                            viewplot = TRUE)  
##### This returns the same result but using plot_impactsnopar()
imp_nparvar <- impactsnopar(gamsar, listw = lwsp_it, alpha = 0.05,
                            viewplot = FALSE)
plot_impactsnopar(imp_nparvar, data = unemp_it, 
                   smooth = TRUE)

