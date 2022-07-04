### demo with the whole set of examples of plot_sptime() ########
library(pspatreg)
data(unemp_it, package = "pspatreg")
lwsp_it <- spdep::mat2listw(Wsp_it)
###### FORMULA OF THE MODEL
form3d_psanova <- unrate ~ partrate + agri + cons +
                  pspl(serv, nknots = 15) + 
                  pspl(empgrowth, nknots = 20) +
                  pspt(long, lat, year, 
                       nknots = c(18, 18, 8),
                       psanova = TRUE, 
                       nest_sp1 = c(1, 2, 3), 
                       nest_sp2 = c(1, 2, 3),
                       nest_time = c(1, 2, 2))
####### FIT the model
sp3danova <- pspatfit(form3d_psanova, 
                      data = unemp_it, 
                      listw = lwsp_it, 
                      method = "Chebyshev") 

summary(sp3danova)

######## Plot of temporal trend for each province 
plot_sptime(sp3danova, 
            data = unemp_it, 
            time_var = "year", 
            reg_var = "prov")
