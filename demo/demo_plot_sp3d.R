### demo with the whole set of examples of plot_sp3d() ########
library(pspatreg)
library(sf)
data(unemp_it, package = "pspatreg")
lwsp_it <- spdep::mat2listw(Wsp_it)
unemp_it_sf <- st_as_sf(dplyr::left_join(
                               unemp_it, map_it,  
                        by = c("prov" = "COD_PRO")))
######## FORMULA of the model
form3d_psanova_restr <- unrate ~ partrate + agri + cons +
                        pspl(serv, nknots = 15) + 
                        pspl(empgrowth, nknots = 20) +
                        pspt(long, lat, year, 
                             nknots = c(18, 18, 8),
                             psanova = TRUE, 
                             nest_sp1 = c(1, 2, 3), 
                             nest_sp2 = c(1, 2, 3),
                             nest_time = c(1, 2, 2),
                             f1t = FALSE, f2t = FALSE)

####### FIT the model
sp3danova <- pspatfit(form3d_psanova_restr, 
                      data = unemp_it_sf)
summary(sp3danova)                          

###### Plot spatio-temporal trends for different years
plot_sp3d(sp3danova, data = unemp_it_sf, 
          time_var = "year", 
          time_index = c(1996, 2005, 2019),
          addmain = FALSE, addint = FALSE)
###### Plot of spatio-temporal trend, main effects 
######      and interaction effect for a year
plot_sp3d(sp3danova, data = unemp_it_sf, 
          time_var = "year", 
          time_index = c(2019),
          addmain = TRUE, addint = TRUE)
          
#### Plot of temporal trend for each province
plot_sptime(sp3danova, 
            data = unemp_it, 
            time_var = "year", 
            reg_var = "prov") 
