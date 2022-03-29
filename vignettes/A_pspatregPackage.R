## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, cache = FALSE)

## -----------------------------------------------------------------------------
formula <- y ~ x1 + x2 + pspl(x3, nknots = 15) + pspl(x4, nknots = 20) +
                  pspt(long, lat, year, nknots = c(18,18,8),
                       psanova = TRUE, 
                       nest_sp1 = c(1, 2, 3), 
                       nest_sp2 = c(1, 2, 3),
                       nest_time = c(1, 2, 2))

