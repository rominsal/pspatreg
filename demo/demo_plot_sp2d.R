### demo with the whole set of examples of plot_sp2d() ########
######## EXAMPLE 2D WITH AMES DATA 
######## getting and preparing the data
library(pspatreg, package = "pspatreg")
library(spdep)
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
######## formula of the model 
form2d <- lnSale_Price ~ Fireplaces + Garage_Cars +
          pspl(lnLot_Area, nknots = 20) + 
          pspl(lnTotal_Bsmt_SF, nknots = 20) +
          pspl(lnGr_Liv_Area, nknots = 20) +
          pspt(Longitude, Latitude, 
               nknots = c(10, 10), 
               psanova = FALSE)
######## fit the model
sp2dsar <- pspatfit(form2d, data = ames_sf1, 
                    listw = lw_ames, 
                    method = "Chebyshev", 
                    type = "sar")
summary(sp2dsar)
####### plot spatial trend for spatial point coordinate
plot_sp2d(sp2dsar, data = ames_sf1)
###### MODEL WITH ANOVA DESCOMPOSITION
form2d_psanova <- lnSale_Price ~ Fireplaces + Garage_Cars +
                   pspl(lnLot_Area, nknots = 20) + 
                   pspl(lnTotal_Bsmt_SF, nknots = 20) +
                   pspl(lnGr_Liv_Area, nknots = 20) +
                   pspt(Longitude, Latitude, 
                        nknots = c(10, 10), 
                        psanova = TRUE)
      
sp2danovasar <- pspatfit(form2d_psanova, 
                        data = ames_sf1, 
                        listw = lw_ames,
                        method = "Chebyshev", 
                        type = "sar")
summary(sp2danovasar)                         
            
###### PLOT ANOVA DESCOMPOSITION MODEL
plot_sp2d(sp2danovasar, data = ames_sf1, 
          addmain = TRUE, addint = TRUE)
