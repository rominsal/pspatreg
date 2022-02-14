# Function to build G matrix
build_G2d <- function(la, lg, env) {
  if (env$nvarspt > 0) {
    # if (!is.null(g1u) && !is.null(g2u) && 
    #     !is.null(g1b) && !is.null(g2b))     
    g1u <- lg$g1u
    g2u <- lg$g2u
    g1b <- lg$g1b
    g2b <- lg$g2b
    if(!env$psanova) { # no psanova
      Ginv <- c((1/la[2])*g1u, (1/la[3])*g2u, 
                (1/la[2])*g1b + (1/la[3])*g2b)
      Ginv_eff <- Ginv
    } else { # psanova
      g1v <- lg$g1v
      g2v <- lg$g2v
      Ginv_eff <- NULL
      if (env$f1_main) { 
        Ginv_f1 <- (1/la[2])*g1u
        Ginv_eff <- c(Ginv_eff, Ginv_f1)
      } else { Ginv_f1 <- rep(0, length(g1u)) }
      if (env$f2_main) {
        Ginv_f2 <- (1/la[3])*g2u
        Ginv_eff <- c(Ginv_eff,Ginv_f2) } else { 
          Ginv_f2 <- rep(0,length(g2u)) }
      if (env$f12_int) {
        Ginv_f12 <- c((1/la[4])*g1v, (1/la[5])*g2v, 
                      (1/la[4])*g1b + (1/la[5])*g2b)
        Ginv_eff <- c(Ginv_eff,Ginv_f12) } else {
          g12u <- lg$g12u
          g21u <- lg$g21u
          g12b <- lg$g12b 
          Ginv_f12 <- rep(0, length(g12u) + 
                            length(g21u) + length(g12b)) }
      Ginv <- c(Ginv_f1, Ginv_f2, Ginv_f12)
    }
  } else {  Ginv <- Ginv_eff <- NULL } 
  if (env$nvarnopar>0) {
    for (k in 1:env$nvarnopar) {
      Ginv <- c(Ginv,
                (1/la[length(la) - env$nvarnopar + k]) * 
                  env$dnoparlist[[k]])
      Ginv_eff <- c(Ginv_eff,
                    (1/la[length(la) - env$nvarnopar + k]) * 
                      env$dnoparlist[[k]]) } 
    }
  G <- ifelse(Ginv!=0, 1/Ginv, 0) # Keep zeros in same positions
  G_eff <- 1/Ginv_eff
  res <- list(G = G, Ginv = Ginv, G_eff = G_eff, 
              Ginv_eff = Ginv_eff) }
##########################################################################

# Function to build G matrix
build_G3d <- function(la, lg, env) {
    g1u <- lg$g1u 
    g2u <- lg$g2u
    g3u <- lg$g3u
    g11b <- lg$g11b
    g21b <- lg$g21b
    g12b <- lg$g12b
    g31b <- lg$g31b
    g22b <- lg$g22b
    g32b <- lg$g32b
    g1t <- lg$g1t
    g2t <- lg$g2t
    g3t <- lg$g3t
    if(!env$psanova){
        Ginv <- c((1/la[2])*g1u,
                  (1/la[3])*g2u,
                  (1/la[4])*g3u,
                  (1/la[2])*g11b + (1/la[3])*g21b,
                  (1/la[2])*g12b + (1/la[4])*g31b,
                  (1/la[3])*g22b + (1/la[4])*g32b,
                  (1/la[2])*g1t + (1/la[3])*g2t + (1/la[4])*g3t)
        Ginv_eff <- Ginv
    } else {
        # NUEVA DEFINICIÓN DE Ginv Y Ginv_eff EN FUNCIÓN DE LAS INTERACCIONES INCLUIDAS
        g12u <- lg$g12u
        g21u <- lg$g21u
        g13u <- lg$g13u
        g31u <- lg$g31u
        g13b <- lg$g13b
        g23u <- lg$g23u
        g32u <- lg$g32u
        g23b <- lg$g23b
        g123u <- lg$g123u
        g213u <- lg$g213u
        g321u <- lg$g321u
        g123b <- lg$g123b
        g213b <- lg$g213b
        g132b <- lg$g132b
        g312b <- lg$g312b
        g231b <- lg$g231b
        g321b <- lg$g321b
        Ginv_eff <- NULL
        if(env$f1_main){
            Ginv_f1 <- (1/la[2])*g1u
            Ginv_eff <- c(Ginv_eff, Ginv_f1)
        } else { Ginv_f1 <- rep(0, length(g1u)) }
        if (env$f2_main) {
            Ginv_f2 <- (1/la[3])*g2u
            Ginv_eff <- c(Ginv_eff,Ginv_f2) } else {
                Ginv_f2 <- rep(0,length(g2u)) }
        if (env$ft_main) {
            Ginv_ft <- (1/la[4])*g3u
            Ginv_eff <- c(Ginv_eff,Ginv_ft) } else {
                Ginv_ft <- rep(0,length(g3u)) }
        if (env$f12_int) {
            Ginv_f12 <- c((1/la[5])*g12u, (1/la[6])*g21u,
                          (1/la[5])*g12b + (1/la[6])*g21b)
            Ginv_eff <- c(Ginv_eff,Ginv_f12) } else {
                Ginv_f12 <- rep(0,length(g12u) + length(g21u) + 
                                  length(g12b))
            }
        if (env$f1t_int) {
            Ginv_f1t <- c((1/la[7])*g13u, (1/la[8])*g31u,
                          (1/la[7])*g13b + (1/la[8])*g31b)
            Ginv_eff <- c(Ginv_eff, Ginv_f1t) } else {
                Ginv_f1t <- rep(0, length(g13u) + length(g31u) + 
                                  length(g13b))
            }
        if (env$f2t_int) {
            Ginv_f2t <- c((1/la[9])*g23u, (1/la[10])*g32u,
                          (1/la[9])*g23b + (1/la[10])*g32b)
            Ginv_eff <- c(Ginv_eff, Ginv_f2t) } else {
                Ginv_f2t <- rep(0, length(g23u) + length(g32u) + 
                                  length(g23b))
            }
        if (env$f12t_int) {
            Ginv_f12t <- c((1/la[11])*g123u, (1/la[12])*g213u,
                           (1/la[13])*g321u,
                           (1/la[11])*g123b + (1/la[12])*g213b,
                           (1/la[11])*g132b + (1/la[13])*g312b,
                           (1/la[12])*g231b + (1/la[13])*g321b,
                           (1/la[11])*g1t + (1/la[12])*g2t +
                               (1/la[13])*g3t)
            Ginv_eff <- c(Ginv_eff,Ginv_f12t) } else {
            Ginv_f12t <- rep(0,length(g123u) + length(g213u) +
                               length(g321u) + length(g123b) +
                               length(g132b) + length(g231b) +
                               length(g1t))
            }
        # CONSTRUCCIÓN Ginv COMPLETA (INCLUYENDO 0's SI LOS HAY)
        Ginv <- c(Ginv_f1, Ginv_f2, Ginv_ft, Ginv_f12, Ginv_f1t,
                  Ginv_f2t, Ginv_f12t)
    }
    if (env$nvarnopar > 0) {
        for (k in 1:env$nvarnopar) {
            Ginv <- c(Ginv,
                      (1/la[length(la) - (env$nvarnopar+2)+k]) * 
                        env$dnoparlist[[k]])
            Ginv_eff <- c(Ginv_eff,
                          (1/la[length(la) - (env$nvarnopar+2)+k]) * 
                            env$dnoparlist[[k]])
        }
    }
    G <- ifelse(Ginv != 0, 1/Ginv, 0) # Keep zeros in same positions
    G_eff <- 1/Ginv_eff
    res <- list(G = G, Ginv = Ginv, G_eff = G_eff, 
                Ginv_eff = Ginv_eff)
}
