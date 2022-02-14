# Functions to compute vectors of variance component matrix of random effects
# in models 2d and 3d
##############################################################################
par_var_comp2d <- function(la, env) {
  np_fixed <- ncol(env$Xfull) 
  pord <- env$pordspt
  if (env$nvarspt > 0) { 
    if(!env$psanova) { # spatial trend no psanova
      tau_init = rep(1, 2)
      c1 <- env$cspt[1] 
      c2 <- env$cspt[2]
      d1 <- env$dsptlist$sp1
      d2 <- env$dsptlist$sp2
      g1u <- as.numeric(rep(1, pord[2]) %x% d1)
      g2u <-  as.numeric(d2 %x% rep(1, pord[1]))
      g1b <- as.numeric(rep(1, c2 - pord[2]) %x% d1)
      g2b <- as.numeric(d2 %x% rep(1, c1 - pord[1]))
      if (!(length(g1b) == length(g2b))) 
        stop("lengths g1b and g2b different...")
      # Number of parameters in each part
      np <-  c(np_fixed, length(g1u), length(g2u), length(g1b))
      names(np) <- c("fixed", "(c1-pord[1])*pord[2]",
                     "(c2-pord[2])*pord[1]",
                     "(c1-pord[1])*(c2-pord[2])")
      np_eff <- np
      # VIP: SPATIAL TREND 2d. NO PS-ANOVA 
      # DEFINE Zspt <- cbind(Rten2(X2, Z1), Rten2(Z2, X1), Rten2(Z2, Z1))
      G1inv.n <- c(g1u, rep(0, np[3]), g1b)
      G2inv.n <- c(rep(0, np[2]), g2u, g2b)
      # Initialize variance components
      la <- c(la, tau_init)
      res <- list(np = np, np_eff = np_eff,
                  g1u = g1u, g2u = g2u, g1b = g1b, g2b = g2b,
                  G1inv.n = G1inv.n, G2inv.n = G2inv.n)
    } else { # spatial trend psanova case
      tau_init = rep(1, 4)
      f1_main <- env$f1_main
      f2_main <- env$f2_main
      f12_int <- env$f12_int
      c1 <- env$cspt[c(grepl("sp1", names(env$cspt)))]
      c2 <- env$cspt[c(grepl("sp2", names(env$cspt)))]
      d1 <- env$dsptlist$sp1_main
      d2 <- env$dsptlist$sp2_main
      d1.2 <- env$dsptlist$sp1_int2ord
      d2.2 <- env$dsptlist$sp2_int2ord
      # # Main Effects
      g1u <- d1
      g2u <- d2
      # # Second Order Interaction
      g1v <- rep(1, pord[2]- 1) %x% d1.2
      g2v <- d2.2 %x% rep(1, pord[1]- 1)
      g1b <- rep(1,c2[2]-pord[2]) %x% d1.2
      g2b <- d2.2 %x% rep(1,c1[2] - pord[1])
      np_f1 <- (c1[1] - pord[1])
      np_f2 <- (c2[1] - pord[2])
      np_f12 <- c((c1[2] - pord[1]) * (pord[2]-1),
                  (c2[2] - pord[2]) * (pord[1]-1),
                  (c1[2] - pord[1]) * (c2[2]-pord[2]))
      # Number of parameters in each part
      np <- c(np_fixed, np_f1, np_f2, np_f12)
      np_eff <- np_fixed
      if (f1_main) np_eff <- c(np_eff, np_f1)
      if (f2_main) np_eff <- c(np_eff, np_f2)
      if (f12_int) np_eff <- c(np_eff, np_f12)
      # # Tau1, Tau2 (Main Effects)
      G1inv.n <- c(g1u, rep(0, sum(np[3:6])))
      G2inv.n <- c(rep(0, np[2]), g2u, rep(0, sum(np[4:6])))
      # # Tau3, Tau4 (f12 interaction)
      G3inv.n <- c(rep(0, sum(np[c(2,3)])), g1v, 
                   rep(0, np[5]), g1b)
      G4inv.n <- c(rep(0, sum(np[c(2,3,4)])), 
                   g2v, g2b)
      # # Change Ginv.n and tau = 0 if any
      # # interaction term is not included in psanova
      if (env$f1_main) { 
        la <- c(la, tau_init[1]) } else { 
          G1inv.n <- rep(0, length(G1inv.n))
          la <- c(la, 0) }
      if (env$f2_main) { 
        la <- c(la,tau_init[2]) } else { 
          G2inv.n <- rep(0, length(G2inv.n)) 
          la <- c(la, 0) }
      if (env$f12_int) { 
        la <- c(la, tau_init[3:4]) } else { 
          G3inv.n <- rep(0, length(G3inv.n)) 
          G4inv.n <- rep(0, length(G4inv.n)) 
          la <- c(la, 0, 0) }
      res <- list(g1u = g1u, g2u = g2u, 
                  g1v = g1v, g2v = g2v,
                  g1b = g1b, g2b = g2b,
                  G1inv.n = G1inv.n, G2inv.n = G2inv.n,
                  G3inv.n = G3inv.n, G4inv.n = G4inv.n)
    }
  } else { # Without spatial trend
    res <- list()
    np <- np_eff <- np_fixed
  }
  if (env$nvarnopar > 0) { # Add parameters corresponding to Znopar
    np <- c(np, env$cnopar - env$pordnopar)
    np_eff <- c(np_eff, env$cnopar - env$pordnopar)
    # add smoothing parameters of non-parametic covariates
    #if (is.null(tau_nopar_init)){
    la <- c(la, rep(1, env$nvarnopar))
    #   } else {
    #    la<-c(res$la,tau_nopar_init) }
    if (env$nvarspt > 0) {
      np_nvarnopar <- sum(np[-c(1:(length(np) - env$nvarnopar))])
      res$G1inv.n <- c(res$G1inv.n, rep(0, np_nvarnopar))
      res$G2inv.n <- c(res$G2inv.n, rep(0, np_nvarnopar))
      if(env$psanova){
        res$G3inv.n <- c(res$G3inv.n, rep(0, np_nvarnopar))
        res$G4inv.n <- c(res$G4inv.n, rep(0, np_nvarnopar))
      }
    }
  }
  res$np <- np 
  res$np_eff <- np_eff
  res$la <- la
  return(res)
}

###############################################################################
par_var_comp3d <- function(la, env) {
  np_fixed <- ncol(env$Xfull)
  pord <- env$pordspt
  if (env$nvarspt > 0) {
    if(!env$psanova) {
      tau_init <- rep(1,3)
      c1 <- env$cspt[1]
      c2 <- env$cspt[2]
      c3 <- env$cspt[3]
      d1 <- env$dsptlist$sp1
      d2 <- env$dsptlist$sp2
      d3 <- env$dsptlist$time
      g1u <- d1 %x% rep(1, pord[2]) %x% rep(1, pord[3])
      g2u <- rep(1, pord[1]) %x% d2 %x% rep(1, pord[3])
      g3u <- rep(1, pord[1]) %x% rep(1, pord[2]) %x% d3
      g11b <- d1 %x% rep(1, c2 - pord[2]) %x% rep(1, pord[3])
      g12b <- d1 %x% rep(1, pord[2]) %x% rep(1, c3 - pord[3])
      g21b <- rep(1, c1 - pord[1]) %x% d2 %x% rep(1, pord[3])
      g22b <- rep(1, pord[1]) %x% d2 %x% rep(1, c3 - pord[3])
      g31b <- rep(1, c1 - pord[1]) %x% rep(1, pord[2]) %x% d3
      g32b <- rep(1, pord[1]) %x% rep(1, c2 - pord[2]) %x% d3
      g1t <- d1 %x% rep(1, c2 - pord[2]) %x% rep(1, c3 - pord[3])
      g2t <- rep(1, c1 - pord[1]) %x% d2 %x% rep(1, c3 - pord[3])
      g3t <- rep(1, c1 - pord[1]) %x% rep(1, c2 - pord[2]) %x% d3
      # Number of parameters in each part
      np <- c(np_fixed,
              (c1 - pord[1])*prod(pord[2:3]),
              (c2 - pord[2])*prod(pord[c(1,3)]),
              (c3 - pord[3])*prod(pord[1:2]),
              (c1 - pord[1])*(c2 - pord[2])*pord[3],
              (c1 - pord[1])*(c3 - pord[3])*pord[2],
              (c2 - pord[2])*(c3 - pord[3])*pord[1],
              (c1 - pord[1])*(c2 - pord[2])*(c3 - pord[3]))
      names(np) <- c("fixed","(c1-pord[1])*prod(pord[2:3])",
                     "(c2-pord[2])*prod(pord[c(1,3)])",
                     "(c3-pord[3])*prod(pord[1:2])",
                     "(c1-pord[1])*(c2-pord[2])*pord[3]",
                     "(c1-pord[1])*(c3-pord[3])*pord[2]",
                     "(c2-pord[2])*(c3-pord[3])*pord[1]",
                     "(c1-pord[1])*(c2-pord[2])*(c3-pord[3])")
      np_eff <- np
      G1inv.n <- c(g1u, rep(0, sum(np[3:4])), g11b, g12b, 
                   rep(0, np[7]), g1t)
      G2inv.n <- c(rep(0, np[2]), g2u, rep(0, np[4]), g21b, 
                   rep(0, np[6]),
                   g22b, g2t)
      G3inv.n <- c(rep(0, sum(np[2:3])), g3u, rep(0, np[5]), 
                   g31b, g32b, g3t)
      la <- c(la, tau_init)
      res <- list(np = np, np_eff = np_eff,
                  g1u = g1u, g2u = g2u, g3u = g3u,
                  g11b = g11b, g12b = g12b,
                  g21b = g21b, g22b = g22b,
                  g31b = g31b, g32b = g32b,
                  g1t = g1t, g2t = g2t, g3t = g3t,
                  G1inv.n = G1inv.n, G2inv.n = G2inv.n, G3inv.n = G3inv.n)
    } else { #psanova case
      tau_init <- rep(1, 12)
      c1 <- env$cspt[c(grepl("sp1", names(env$cspt)))]
      c2 <- env$cspt[c(grepl("sp2", names(env$cspt)))]
      c3 <- env$ cspt[c(grepl("time", names(env$cspt)))]
      d1 <- env$dsptlist$sp1_main
      d2 <- env$dsptlist$sp2_main
      d3 <- env$dsptlist$time_main
      d1.2 <- env$dsptlist$sp1_int2ord 
      d2.2 <- env$dsptlist$sp2_int2ord
      d3.2 <- env$dsptlist$time_int2ord
      d1.3 <- env$dsptlist$sp1_int3ord
      d2.3 <- env$dsptlist$sp2_int3ord
      d3.3 <- env$dsptlist$time_int3ord
      # Main Effects
      g1u <- d1
      g2u <- d2
      g3u <- d3
      # Second Order Interaction
      # Interactions Z_i,1_j,x_k
      g12u <- d1.2 %x% rep(1, pord[2] -1)
      g13u <- d1.2 %x% rep(1, pord[3] - 1)
      g21u <- rep(1, pord[1] - 1) %x% d2.2
      g23u <- d2.2 %x% rep(1, pord[3] - 1)
      g31u <- rep(1, pord[1] - 1) %x% d3.2
      g32u <- rep(1, pord[2] - 1) %x% d3.2
      # Interactions Z_i,Z_j,1_k
      g12b <- d1.2 %x% rep(1, c2[2] - pord[2])
      g21b <- rep(1, c1[2] - pord[1]) %x% d2.2
      g13b <- d1.2 %x% rep(1, c3[2] - pord[3])
      g31b <- rep(1, c1[2] - pord[1]) %x% d3.2
      g23b <- d2.2 %x% rep(1, c3[2] - pord[3])
      g32b <- rep(1, c2[2] - pord[2]) %x% d3.2
      # Interactions third order
      # Interactions Z_i,x_j,x_k
      g123u <- d1.3 %x% rep(1, pord[2] - 1) %x% rep(1, pord[3] - 1)
      g213u <- rep(1, pord[1] - 1) %x% d2.3 %x% rep(1, pord[3] - 1)
      g321u <- rep(1, pord[1] - 1) %x% rep(1, pord[2] - 1) %x% d3.3
      # Interactions Z_i,Z_j,x_k
      g123b <- d1.3 %x% rep(1, c2[3] - pord[2]) %x% 
        rep(1, pord[3] - 1)
      g213b <- rep(1, c1[3] - pord[1]) %x% d2.3 %x% 
        rep(1, pord[3] - 1)
      g132b <- d1.3 %x% rep(1, pord[2] - 1) %x% 
        rep(1, c3[3] - pord[3])
      g312b <- rep(1, c1[3] - pord[1]) %x% 
        rep(1, pord[2] - 1) %x% d3.3
      g231b <- rep(1, pord[1] - 1) %x% d2.3 %x% 
        rep(1, c3[3] - pord[3])
      g321b <- rep(1, pord[1] - 1) %x% 
        rep(1, c2[3] - pord[2]) %x% d3.3
      # Interactions Z_i,Z_j,Z_k
      g1t <- d1.3 %x% rep(1, c2[3] - pord[2]) %x% 
        rep(1, c3[3] - pord[3])
      g2t <- rep(1, c1[3] - pord[1]) %x% d2.3 %x% 
        rep(1, c3[3] - pord[3])
      g3t <- rep(1, c1[3] - pord[1]) %x% 
        rep(1, c2[3] - pord[2]) %x% d3.3
      np_f1 <- c1[1]-pord[1]
      np_f2 <- c2[1]-pord[2]
      np_ft <- c3[1]-pord[3]
      np_f12 <- c((c1[2]-pord[1])*(pord[2]-1),
                  (c2[2]-pord[2])*(pord[1]-1),
                  (c1[2]-pord[1])*(c2[2]-pord[2]))
      np_f1t <- c((c1[2]-pord[1])*(pord[3]-1),
                  (c3[2]-pord[3])*(pord[1]-1),
                  (c1[2]-pord[1])*(c3[2]-pord[3]))
      np_f2t <- c((c2[2]-pord[2])*(pord[3]-1),
                  (c3[2]-pord[3])*(pord[2]-1),
                  (c2[2]-pord[2])*(c3[2]-pord[3]))
      np_f12t <- c((c1[3]-pord[1])*(pord[2]-1)*(pord[3]-1),
                   (c2[3]-pord[2])*(pord[1]-1)*(pord[3]-1),
                   (c3[3]-pord[3])*(pord[2]-1)*(pord[1]-1),
                   (c1[3]-pord[1])*(c2[3]-pord[2])*(pord[3]-1),
                   (c1[3]-pord[1])*(c3[3]-pord[3])*(pord[2]-1),
                   (c3[3]-pord[3])*(c2[3]-pord[2])*(pord[1]-1),
                   (c1[3]-pord[1])*(c2[3]-pord[2])*(c3[3]-pord[3]))
      np <- c(np_fixed, np_f1, np_f2, np_ft, 
              np_f12, np_f1t, np_f2t, np_f12t)
      #names(np) <- c("fixed","f1","f2","f3","f12","f13","f23","f123")
      np_eff <- np_fixed
      if (env$f1_main) np_eff <- c(np_eff, np_f1)
      if (env$f2_main) np_eff <- c(np_eff, np_f2)
      if (env$ft_main) np_eff <- c(np_eff, np_ft)
      if (env$f12_int) np_eff <- c(np_eff, np_f12)
      if (env$f1t_int) np_eff <- c(np_eff, np_f1t)
      if (env$f2t_int) np_eff <- c(np_eff, np_f2t)
      if (env$f12t_int) np_eff <- c(np_eff, np_f12t)
      # Tau1, Tau2 y Tau3 (Main Effects)
      G1inv.n <- c(g1u, rep(0, sum(np[3:20])))
      G2inv.n <- c(rep(0, np[2]), g2u, rep(0, sum(np[4:20])))
      G3inv.n <- c(rep(0, sum(np[2:3])), g3u, 
                   rep(0,sum(np[5:20])))
      # Tau4, Tau5 (f12 interaction)
      G4inv.n <- c(rep(0, sum(np[2:4])), g12u, rep(0, np[6]),
                   g12b, rep(0, sum(np[8:20])))
      G5inv.n <- c(rep(0, sum(np[2:5])), g21u, g21b, 
                   rep(0, sum(np[8:20])))
      # Tau6, Tau7 (f13 interaction)
      G6inv.n <- c(rep(0,sum(np[2:7])), g13u, rep(0,np[9]),
                   g13b, rep(0, sum(np[11:20])))
      G7inv.n <- c(rep(0, sum(np[2:8])), g31u, g31b, 
                   rep(0, sum(np[11:20])))
      # Tau8, Tau9 (f23 interaction)
      G8inv.n <- c(rep(0, sum(np[2:10])), g23u, rep(0, np[12]),
                   g23b, rep(0, sum(np[14:20])))
      G9inv.n <- c(rep(0, sum(np[2:11])), g32u, g32b, 
                   rep(0, sum(np[14:20])))
      # Tau10, Tau11 y Tau 12 (f123 interaction)
      G10inv.n <- c(rep(0, sum(np[2:13])), g123u, 
                    rep(0, sum(np[15:16])),
                    g123b, g132b, rep(0,np[19]), g1t)
      G11inv.n <- c(rep(0,sum(np[2:14])), g213u, rep(0,np[16]),
                    g213b, rep(0,np[18]), g231b, g2t)
      G12inv.n <- c(rep(0,sum(np[2:15])), g321u, rep(0,np[17]), 
                    g312b, g321b, g3t)
      #  change Ginv.n and tau=0 if any
      # interaction term is not included in psanova
      if (env$f1_main) { 
        la <- c(la, tau_init[1]) 
      } else {
        G1inv.n <- rep(0,length(G1inv.n))
        la <- c(la,0) 
      }
      if (env$f2_main) { 
        la <- c(la, tau_init[2]) 
      } else {
        G2inv.n <- rep(0,length(G2inv.n))
        la <- c(la,0) 
      }
      if (env$ft_main) { 
        la <- c(la, tau_init[3]) 
      } else {
        G3inv.n <- rep(0,length(G3inv.n))
        la <- c(la,0) 
      }
      if (env$f12_int) { 
        la <- c(la, tau_init[4:5]) 
      } else {
        G4inv.n <- rep(0,length(G4inv.n))
        G5inv.n <- rep(0,length(G5inv.n))
        la <- c(la,0,0) 
      }
      if (env$f1t_int) { 
        la <- c(la, tau_init[6:7]) 
      } else {
        G6inv.n <- rep(0,length(G6inv.n))
        G7inv.n <- rep(0,length(G7inv.n))
        la <- c(la,0,0) }
      if (env$f2t_int) { 
        la <- c(la, tau_init[8:9]) 
      } else {
        G8inv.n <- rep(0,length(G8inv.n))
        G9inv.n <- rep(0,length(G9inv.n))
        la <- c(la,0,0) }
      if (env$f12t_int) { 
        la <- c(la, tau_init[10:12]) 
      } else {
        G10inv.n <- rep(0,length(G10inv.n))
        G11inv.n <- rep(0,length(G11inv.n))
        G12inv.n <- rep(0,length(G12inv.n))
        la <- c(la,0,0,0) }
      res <- list(g1u = g1u, g2u = g2u, g3u = g3u,
                  g12u = g12u, 
                  g13u = g13u,
                  g21u = g21u, 
                  g23u = g23u,
                  g31u = g31u, 
                  g32u = g32u,
                  g12b = g12b, 
                  g21b = g21b,
                  g13b = g13b, 
                  g31b = g31b,
                  g23b = g23b, 
                  g32b = g32b,
                  g123u = g123u, 
                  g213u = g213u,
                  g321u = g321u,
                  g123b = g123b, 
                  g213b = g213b,
                  g132b = g132b, 
                  g312b = g312b,
                  g231b = g231b, 
                  g321b = g321b,
                  g1t = g1t, 
                  g2t = g2t, 
                  g3t = g3t,
                  G1inv.n = G1inv.n, 
                  G2inv.n = G2inv.n, 
                  G3inv.n = G3inv.n,
                  G4inv.n = G4inv.n, 
                  G5inv.n = G5inv.n, 
                  G6inv.n = G6inv.n,
                  G7inv.n = G7inv.n, 
                  G8inv.n = G8inv.n, 
                  G9inv.n = G9inv.n,
                  G10inv.n = G10inv.n, 
                  G11inv.n = G11inv.n,
                  G12inv.n = G12inv.n)
    }
  } else { #env$nvarspt == 0
    np <- np_eff <- np_fixed
    res <- list()
  }
  if (env$nvarnopar > 0) { # Add parameters corresponding to Znopar
      np <- c(np, env$cnopar - env$pordnopar)
      np_eff <- c(np_eff, env$cnopar - env$pordnopar)
      # add smoothing parameters of non-parametic covariates
      #if (is.null(tau_nopar_init)){
          la <- c(la, rep(1, env$nvarnopar))
         #   } else {
         #    la<-c(res$la,tau_nopar_init) }
        np_nvarnopar <- sum(np[-c(1:(length(np) - env$nvarnopar))])
        res$G1inv.n <- c(res$G1inv.n, rep(0, np_nvarnopar))
        res$G2inv.n <- c(res$G2inv.n, rep(0, np_nvarnopar))
        res$G3inv.n <- c(res$G3inv.n, rep(0, np_nvarnopar))
        if(env$psanova){
            res$G4inv.n <- c(res$G4inv.n, rep(0, np_nvarnopar))
            res$G5inv.n <- c(res$G5inv.n, rep(0, np_nvarnopar))
            res$G6inv.n <- c(res$G6inv.n, rep(0, np_nvarnopar))
            res$G7inv.n <- c(res$G7inv.n, rep(0, np_nvarnopar))
            res$G8inv.n <- c(res$G8inv.n, rep(0, np_nvarnopar))
            res$G9inv.n <- c(res$G9inv.n, rep(0, np_nvarnopar))
            res$G10inv.n <- c(res$G10inv.n, rep(0, np_nvarnopar))
            res$G11inv.n <- c(res$G11inv.n, rep(0, np_nvarnopar))
            res$G12inv.n <- c(res$G12inv.n, rep(0, np_nvarnopar))
        }
    }
    res$np <- np
    res$np_eff <- np_eff
    res$la <- la
    return(res)
}

