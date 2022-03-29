## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE,message=FALSE, cache = FALSE)

## -----------------------------------------------------------------------------
library(spdep)
library(sf)
library(ggplot2)
library(dplyr)
library(dbscan)
ames <- AmesHousing::make_ames()# Raw Ames Housing Data
ames_sf <- st_as_sf(ames, 
                    coords = c("Longitude", "Latitude"))
ames_sf$Longitude <- ames$Longitude
ames_sf$Latitude <- ames$Latitude

#devtools::load_all()


## ----fig.align='center'-------------------------------------------------------
ggplot(data=ames_sf) + geom_histogram(mapping = aes(x = Sale_Price))

ggplot(data=ames_sf) + geom_histogram(mapping = aes(x = log(Sale_Price)))

ames_sf$lnSale_Price <- log(ames_sf$Sale_Price)
ames_sf$lnLot_Area <- log(ames_sf$Lot_Area)
ames_sf$lnTotal_Bsmt_SF <- log(ames_sf$Total_Bsmt_SF+1)
ames_sf$lnGr_Liv_Area <- log(ames_sf$Gr_Liv_Area)


## -----------------------------------------------------------------------------
ames_sf1 <- ames_sf[duplicated(ames_sf$Longitude)==FALSE,]
coord_sf1 <- cbind(ames_sf1$Longitude,ames_sf1$Latitude)
# Delauney triangulation neighbours (symmetric)
ID <- row.names(as(ames_sf1, "sf"))
col.tri.nb <- tri2nb(coord_sf1)
soi_nb <- graph2nb(soi.graph(col.tri.nb,coord_sf1), row.names=ID)
is.symmetric.nb(soi_nb,verbose=T,force=F)
listW <- nb2listw(soi_nb, style="W",zero.policy=F)

## -----------------------------------------------------------------------------
# Linear Model

formlin <- lnSale_Price ~ lnLot_Area+lnTotal_Bsmt_SF+
  lnGr_Liv_Area+Garage_Cars+Fireplaces

durbinformlin <- ~ lnLot_Area+lnTotal_Bsmt_SF+
  lnGr_Liv_Area+Garage_Cars+Fireplaces

# Semiparametric model without spatial trend
formgam <- lnSale_Price ~ Fireplaces + Garage_Cars +
  pspl(lnLot_Area, nknots = 20) + 
  pspl(lnTotal_Bsmt_SF, nknots = 20) +
  pspl(lnGr_Liv_Area, nknots = 20) 

# Semiparametric model with spatial trend in 2d

form2d <- lnSale_Price ~ Fireplaces + Garage_Cars +
  pspl(lnLot_Area, nknots = 20) + 
  pspl(lnTotal_Bsmt_SF, nknots = 20) +
  pspl(lnGr_Liv_Area, nknots = 20) +
  pspt(Longitude,Latitude, nknots = c(10, 10), 
       psanova = FALSE)

# Semiparametric model with PS-ANOVA spatial trend in 2d

  form2d_psanova <- lnSale_Price ~ Fireplaces + Garage_Cars +
    pspl(lnLot_Area, nknots = 20) + 
    pspl(lnTotal_Bsmt_SF, nknots = 20) +
    pspl(lnGr_Liv_Area, nknots = 20) +
  pspt(Longitude,Latitude, nknots = c(10, 10), 
       psanova = TRUE)

  durbinformnonlin <- ~ Fireplaces + Garage_Cars +
  pspl(lnLot_Area, nknots = 20) 


## -----------------------------------------------------------------------------
# Install pspatreg if needed:
# devtools::install_github("rominsal/pspatreg")
library(pspatreg)
linsar <- pspatfit(formlin, data = ames_sf1,
                   listw = listW, 
                   method = "eigen", type = "sar")
summary(linsar)

## -----------------------------------------------------------------------------
coef(linsar)

## ----fig.align='center'-------------------------------------------------------
fits <- fitted(linsar)
plot(fits, ames_sf1$lnSale_Price)
resids <- residuals(linsar)
plot(fits, resids)

## -----------------------------------------------------------------------------
logLik(linsar)
logLik(linsar, REML = TRUE)

## -----------------------------------------------------------------------------
vcov(linsar)
vcov(linsar, bayesian = TRUE)

## -----------------------------------------------------------------------------
print(linsar)

## -----------------------------------------------------------------------------
eff_parvar_sar <- impactspar(linsar, list_varpar)
summary(eff_parvar_sar)


## -----------------------------------------------------------------------------
spatregsar <- spatialreg::lagsarlm(formlin, data = ames_sf1,
                                   listw = listW, 
                                   method = "eigen") 
summary(spatregsar)
W <- as(listW, "CsparseMatrix")
trMatc <- spatialreg::trW(W, type="mult")
set.seed(1)
spatialreg::impacts(spatregsar,listw=listW)
SAR.impact <- spatialreg::impacts(spatregsar, tr=trMatc, R=200)
list_varpar <- as.character(names(summary(linsar)$bfixed[-1]))
eff_parvar <- impactspar(linsar, list_varpar)
summary(eff_parvar)
# Let's compare direct impacts
round(data.frame(spatialreg_direct=summary(SAR.impact, zstats=TRUE, short=TRUE)$res$direct,sptpreg_direct=summary(eff_parvar_sar)$dir_table[,1]),3)
# Let's compare indirect impacts
round(data.frame(spatialreg_indirect=summary(SAR.impact, zstats=TRUE, short=TRUE)$res$indirect,sptpreg_indirect=summary(eff_parvar_sar)$ind_table[,1]),3)
# Let's compare total impacts
round(data.frame(spatialreg_total=summary(SAR.impact, zstats=TRUE, short=TRUE)$res$total,sptpreg_total=summary(eff_parvar_sar)$tot_table[,1]),3)


## -----------------------------------------------------------------------------
linslx <- pspatfit(formlin, data = ames_sf1,
                   listw = listW, 
                   method = "eigen", 
                   type = "slx", 
                   Durbin = durbinformlin)
summary(linslx)

## -----------------------------------------------------------------------------
eff_parvar_slx <- impactspar(linslx, listw = listW)
summary(eff_parvar_slx)

## -----------------------------------------------------------------------------
anova(linsar, linslx, lrtest = FALSE)

## -----------------------------------------------------------------------------
spatregslx <- spatialreg::lmSLX(formlin, data = ames_sf1,
                                listw = listW) 
SLX.impact <- spatialreg::impacts(spatregslx)
# Let's compare direct impacts
round(data.frame(spatialreg_direct=summary(SLX.impact)$impacts$direct,
                 sptpreg_direct=summary(eff_parvar_slx)$mimpacts[,1]),3)
# Let's compare indirect impacts
round(data.frame(spatialreg_indirect=summary(SLX.impact)$impacts$indirect,
                 sptpreg_indirect=summary(eff_parvar_slx)$mimpacts[,2]),3)
# Let's compare indirect impacts
round(data.frame(spatialreg_total=summary(SLX.impact)$impacts$total,
                 sptpreg_total=summary(eff_parvar_slx)$mimpacts[,3]),3)

## -----------------------------------------------------------------------------
linsdm <- pspatfit(formlin, data = ames_sf1,
                   listw = listW, 
                   method = "eigen", type = "sdm")
summary(linsdm)

## -----------------------------------------------------------------------------
anova(linsar, linsdm, lrtest = TRUE)
anova(linslx, linsdm, lrtest = TRUE)

## -----------------------------------------------------------------------------
eff_parvar_sdm <- impactspar(linsdm, list_varpar)
summary(eff_parvar_sdm)

