XZ3d <-  function(x1,x2,x3,decom=2,psanova=FALSE,
           f1_main=TRUE,f2_main=TRUE,f3_main=TRUE,
           f12_int=TRUE,f13_int=TRUE,f23_int=TRUE,f123_int=TRUE,
           div_sp=c(1,1),div_tm=c(1,1),
           x1lim=NULL,x2lim=NULL,x3lim=NULL,
           knots=c(10,10,10),pord=c(2,2,2),bdeg=c(3,3,3),
           Xpar=NULL,Xnopar=NULL,
           knots_nopar=NULL,bdeg_nopar=NULL,pord_nopar=NULL) {

    # Function to get X (fixed effects) and Z (random effects) matrices
    # of mixed model representation corresponding to semiparametric functions
    # including an ANOVA spatio-temporal trend and, potentially, additional
    # parametric and non-parametric covariates.

    # Assumptions: Spatio-Temporal trend ANOVA allowing main effects;
    # 2nd and 3rd order interactions. It is also allowed nested basis for
    # spatio-temporal trend (in 2nd and 3rd interactions) and
    # the addition of parametric and non-parametric covariates.

    # Inputs:
      # x1, x2 spatial coordinates (same n) and x3 temporal coord.
      # Options for spatial and temporal coordinates included in knots (vector
      # including knots for each spatial and temporal coordinates),
      # pord (order of differences for penalty matrices, quadratic by default)
      # and bdeg (degree of b-splines, cubic by default).
      # Logical variables for anova main and interaction functions (f1_main,
      # f2_main,f3_main,f12_int,f13_int,f23_int,f123_int)
      # Divisors for nested basis of 2nd and 3rd order in div_sp (spatial)
      # and div_tm (temporal)
      # Parametric covariates included as matrix in Xpar
      # Non-parametric covariates included as matrix in Xnopar
      # Options for non-parametric covariates chosen in knots_nopar, bdeg_nopar and
      # pord_nopar

   # Output:
    # X = design matrix for the whole set of fixed effects
    # Z = design matrix for the whole set of random effects
    # X_spt = design matrix for the fixed effects of spatio-temporal trend
    # Z_spt = design matrix for the random effects of spatio-temporal trend
    # X1 = design matrix for the fixed effects of first spatial coordinate
    # Z1 = design matrix for the random effects of first spatial coordinate
    # c1 = number of columns of B-spline basis matrix of first spat. coord.
    # X2 = design matrix for the fixed effects of second spatial coordinate
    # Z2 = design matrix for the random effects of second spatial coordinate
    # c2 = number of columns of B-spline basis matrix of second spat. coord.
    # X3 = design matrix for the fixed effects of temporal coordinate
    # Z3 = design matrix for the random effects of temporal coordinate
    # c3 = number of columns of B-spline basis matrix of temporal coordinate
    # d1 = eigenvalues of SVD decomposition of Penalty Matrix: crossprod(D1)
    # d2 = eigenvalues of SVD decomposition of Penalty Matrix: crossprod(D2)
    # d3 = eigenvalues of SVD decomposition of Penalty Matrix: crossprod(D3)
    # nvar_par = number of parametric covariates
    # nvar_nopar = number of non-parametric covariates
    # Xpar = design matrix of parametric covariates (all included as fixed effects)
    # Xnopar = design matrix of fixed effects for the whole set of
    #         non-parametric covariates.
    # Znopar = design matrix of random effects for the whole set of
    #           non-parametric covariates.
    # c_nopar = number of columns of B-spline basis matrix for the whole set of
    #           non-parametric covariates
    # d_nopar = eigenvalues of SVD decomposition of Penalty Matrix for the whole
    #           set of non-parametric covariates: crossprod(D_nopar)

  # length of spatial and temporal coordinates
  nsp <- length(x1); ntime <- length(x3)

  # Build design matrix for B-spline basis of x1, x2 and x3.
  # Unpenalized part as X=B*U_n if decom=1 or X=[1|x|...|x^(d-1)] if decom=2
  if(is.null(x1lim))  x1lim <- c(min(x1)-0.01, max(x1)+0.01)
  MM1 <- MM_basis(x1, x1lim[1], x1lim[2], knots[1], bdeg[1], pord[1],
                  decom=decom)
  X1 <- MM1$X; Z1<-MM1$Z; B1<-MM1$B; c1 <- ncol(B1); d1 <- MM1$d

  if(is.null(x2lim))  x2lim <- c(min(x2)-0.01, max(x2)+0.01)
  MM2 <- MM_basis(x2, x2lim[1], x2lim[2], knots[2], bdeg[2], pord[2],
                  decom=decom)
  X2 <- MM2$X; Z2<-MM2$Z;	B2<-MM2$B; c2 <- ncol(B2); d2 <- MM2$d

  rm(MM1,MM2)

  if ( ntime>1 ) {
    if(is.null(x3lim))  x3lim <- c(min(x3)-0.01, max(x3)+0.01)
    MM3 <- MM_basis(x3, x3lim[1], x3lim[2], knots[3], bdeg[3], pord[3],
                    decom=decom)
    X3 <- MM3$X; Z3<-MM3$Z;	B3<-MM3$B
    c3 <- ncol(B3); d3 <- MM3$d
    rm(MM3)
  } else {
   X3 <- Z3 <- B3 <- Matrix(1)
   c3 <- d3 <- 1
  }

  if(!psanova){
      X_spt <- kronecker(Rten2(X1,X2),X3)
      Z_spt <- cbind(kronecker(Rten2(Z1,X2), X3),
                     kronecker(Rten2(X1,Z2), X3),
                     kronecker(Rten2(X1,X2), Z3),
                     kronecker(Rten2(Z1,Z2), X3),
                     kronecker(Rten2(Z1,X2), Z3),
                     kronecker(Rten2(X1,Z2), Z3),
                     kronecker(Rten2(Z1,Z2), Z3))
  } else {
      # Use nested basis for interactions of 2nd and 3rd order
      MM1.2 <- MM_basis(x1, x1lim[1], x1lim[2], knots[1]/div_sp[1], bdeg[1],
                        pord[1],decom=decom)
      X1.2 <- MM1.2$X; Z1.2<-MM1.2$Z; B1.2<-MM1.2$B
      c1.2 <- ncol(B1.2); d1.2 <- MM1.2$d
      MM1.3 <- MM_basis(x1, x1lim[1], x1lim[2], knots[1]/div_sp[2], bdeg[1],
                        pord[1],decom=decom)
      X1.3 <- MM1.3$X; Z1.3<-MM1.3$Z; B1.3<-MM1.3$B
      c1.3 <- ncol(B1.3); d1.3 <- MM1.3$d
      c1 <- c(c1,c1.2,c1.3)
      d1 <- list(d1=d1,d1.2=d1.2,d1.3=d1.3)
      rm(MM1.2,MM1.3)

      MM2.2 <- MM_basis(x2, x2lim[1], x2lim[2], knots[2]/div_sp[1], bdeg[1],
                        pord[1],decom=decom)
      X2.2 <- MM2.2$X; Z2.2<-MM2.2$Z; B2.2<-MM2.2$B
      c2.2 <- ncol(B2.2); d2.2 <- MM2.2$d
      MM2.3 <- MM_basis(x2, x2lim[1], x2lim[2], knots[2]/div_sp[2], bdeg[1],
                        pord[1],decom=decom)
      X2.3 <- MM2.3$X; Z2.3<-MM2.3$Z; B2.3<-MM2.3$B
      c2.3 <- ncol(B2.3); d2.3 <- MM2.3$d
      c2 <- c(c2,c2.2,c2.3)
      d2 <- list(d2=d2,d2.2=d2.2,d2.3=d2.3)
      rm(MM2.2,MM2.3)

      MM3.2 <- MM_basis(x3, x3lim[1], x3lim[2], knots[3]/div_tm[1], bdeg[1],
                        pord[1],decom=decom)
      X3.2 <- MM3.2$X; Z3.2<-MM3.2$Z; B3.2<-MM3.2$B
      c3.2 <- ncol(B3.2); d3.2 <- MM3.2$d
      MM3.3 <- MM_basis(x3, x3lim[1], x3lim[2], knots[3]/div_tm[2], bdeg[1],
                        pord[1],decom=decom)
      X3.3 <- MM3.3$X; Z3.3<-MM3.3$Z; B3.3<-MM3.3$B
      c3.3 <- ncol(B3.3); d3.3 <- MM3.3$d
      c3 <- c(c3,c3.2,c3.3)
      d3 <- list(d3=d3,d3.2=d3.2,d3.3=d3.3)
      rm(MM3.2,MM3.3)

      one1. <- X1[,1,drop=FALSE]
      one2. <- X2[,1,drop=FALSE]
      one3. <- X3[,1,drop=FALSE]
      x1. <- X1[,-1,drop=FALSE]
      x2. <- X2[,-1,drop=FALSE]
      x3. <- X3[,-1,drop=FALSE]

      X.int <- kronecker(Rten2(one1.,one2.),one3.)
      X.f1 <-  kronecker(Rten2(x1.,one2.),one3.)
      X.f2 <-  kronecker(Rten2(one1.,x2.),one3.)
      X.f3 <-  kronecker(Rten2(one1.,one2.),x3.)
      X.f12 <- kronecker(Rten2(x1.,x2.),one3.)
      X.f13 <- kronecker(Rten2(x1.,one2.),x3.)
      X.f23 <- kronecker(Rten2(one1.,x2.),x3.)
      X.f123 <- kronecker(Rten2(x1.,x2.),x3.)
      X_spt <- X.int
      if(f1_main) X_spt <- cbind(X_spt,X.f1)
      if(f2_main) X_spt <- cbind(X_spt,X.f2)
      if(f3_main) X_spt <- cbind(X_spt,X.f3)
      if(f12_int) X_spt <- cbind(X_spt,X.f12)
      if(f13_int) X_spt <- cbind(X_spt,X.f13)
      if(f23_int) X_spt <- cbind(X_spt,X.f23)
      if(f123_int) X_spt <- cbind(X_spt,X.f123)

      Z.f1 <- kronecker(Rten2(Z1, one2.), one3.) # g1u
      Z.f2 <- kronecker(Rten2(one1., Z2), one3.) # g2u
      Z.f3 <- kronecker(Rten2(one1., one2.), Z3) # g3u
      # g12u | g21u | g12b+g21b
      Z.f12 <- cbind(kronecker(Rten2(Z1.2, x2.), one3.),
                     kronecker(Rten2(x1., Z2.2), one3.),
                     kronecker(Rten2(Z1.2, Z2.2), one3.))
      # g13u | g31u | g13b+g31b
      Z.f13 <- cbind(kronecker(Rten2(Z1.2, one2.), x3.),
                     kronecker(Rten2(x1., one2.), Z3.2),
                     kronecker(Rten2(Z1.2, one2.), Z3.2))
      # g23u | g32u | g23b+g32b
      Z.f23 <- cbind(kronecker(Rten2(one1., Z2.2), x3.),
                     kronecker(Rten2(one1., x2.), Z3.2),
                     kronecker(Rten2(one1., Z2.2), Z3.2))
      # g123u | g213u | g321u | g123b+g213b | g132b+g312b |
      # g231b+g321b | g1t+g2t+g3t
      Z.f123 <- cbind(kronecker(Rten2(Z1.3, x2.), x3.),
                      kronecker(Rten2(x1., Z2.3), x3.),
                      kronecker(Rten2(x1., x2.), Z3.3),
                      kronecker(Rten2(Z1.3, Z2.3), x3.),
                      kronecker(Rten2(Z1.3, x2.), Z3.3),
                      kronecker(Rten2(x1., Z2.3), Z3.3),
                      kronecker(Rten2(Z1.3, Z2.3), Z3.3))
      Z_spt <- NULL
      if(f1_main) Z_spt <- cbind(Z_spt,Z.f1)
      if(f2_main) Z_spt <- cbind(Z_spt,Z.f2)
      if(f3_main) Z_spt <- cbind(Z_spt,Z.f3)
      if(f12_int) Z_spt <- cbind(Z_spt,Z.f12)
      if(f13_int) Z_spt <- cbind(Z_spt,Z.f13)
      if(f23_int) Z_spt <- cbind(Z_spt,Z.f23)
      if(f123_int) Z_spt <- cbind(Z_spt,Z.f123)
  }
  rm(X1,X2,X3,Z1,Z2,Z3)
  X <- X_spt
  Z <- Z_spt
  nvar_par <- 0
  if(!is.null(Xpar)){
      nvar_par <- ncol(Xpar)
      X <- cbind(Xpar,X)
  }
  nvar_nopar <- 0;  Znopar <- NULL
  d_nopar <- as.list(NULL); c_nopar <- NULL
  if(!is.null(Xnopar)){
    nvar_nopar <- ncol(Xnopar)
    if(is.null(knots_nopar)) knots_nopar <- rep(10,nvar_nopar)
    if(is.null(bdeg_nopar)) bdeg_nopar <- rep(3,nvar_nopar)
    if(is.null(pord_nopar)) pord_nopar <- rep(2,nvar_nopar)
    for(i in 1:nvar_nopar){
      var_np <- Xnopar[,i]
      MM.var <- MM_basis(var_np,min(var_np)-0.01,max(var_np)+0.01,
                         knots_nopar[i],bdeg_nopar[i],pord_nopar[i],decom=1)
      X.var <- MM.var$X; Z.var <- MM.var$Z
      B.var <- MM.var$B; d.var <- MM.var$d
      #Xnopar <- cbind(Xnopar,X.var[,c(-1)]) # Quita columna intercepto
      Xnopar <- cbind(Xnopar,X.var[,-c(1)]) # Quita columna intercepto
      Znopar <- cbind(Znopar,Z.var)
      d_nopar[[i]] <- d.var
      c_nopar <- c(c_nopar,ncol(B.var))
    }
    X <- cbind(X,Xnopar)
    Z <- cbind(Z,Znopar)
  }

  res <- list(X=X,Z=Z,X_spt=X_spt,Z_spt=Z_spt,
              c1=c1,c2=c2,c3=c3,d1=d1,d2=d2,d3=d3,
              knots=knots,pord=pord,bdeg=bdeg,
              nvar_par=nvar_par,nvar_nopar=nvar_nopar,
              Xpar=Xpar,Xnopar=Xnopar,Znopar=Znopar,
              c_nopar=c_nopar,d_nopar=d_nopar,knots_nopar=knots_nopar,
              bdeg_nopar=bdeg_nopar,pord_nopar=pord_nopar)
    res
}


buildXZ <-
  function(x1,x2,x3=NULL,
           ANOVA=TRUE,ANOVA_part=FALSE,
           f1_main=TRUE,f2_main=TRUE,f3_main=TRUE,
           f12_int=TRUE,f13_int=TRUE,f23_int=TRUE,f123_int=TRUE,
           div_sp=c(1,1),div_tm=c(1,1),
           x1lim=NULL, x2lim=NULL,x3lim=NULL,
           knots=c(10,10,10),pord=c(2,2,2),bdeg=c(3,3,3),
           Xpar=NULL,Xnopar=NULL,
           knots_nopar=NULL,bdeg_nopar=NULL,pord_nopar=NULL) {
    # VIP: Las coordenadas espaciales x1 y x2 deben tener distinta
    # longitud en 2 y 3 dimensiones.
    if (!is.null(x1)){
      n <- length(x1)
    } else {
      if(!is.null(x2))
      {
        n <- length(x2)
      } else {
        if(!is.null(Xnopar)){
          n <- nrow(Xnopar)
        } else {n <- nrow(Xpar)}
      }
    }
    if(ANOVA){decom=2} else {decom=1}
    #decom=1
      # Unpenalized part as X=B*U_n if decom=1 or X=[1|x|...|x^(d-1)] if decom=2
    if(!is.null(x1))
    {
      if(is.null(x1lim))  x1lim <- c(min(x1)-0.01, max(x1)+0.01)
      MM1 <- MM_basis(x1, x1lim[1], x1lim[2], knots[1], bdeg[1], pord[1],decom=decom)
      X1 <- MM1$X; Z1<-MM1$Z; B1<-MM1$B
      c1 <- ncol(B1); d1 <- MM1$d
    } else { X1<-NULL; Z1<-NULL; B1<-NULL; c1<-NULL; d1<-NULL}

    if(!is.null(x2))
    {
      if(is.null(x2lim))  x2lim <- c(min(x2)-0.01, max(x2)+0.01)
      MM2 <- MM_basis(x2, x2lim[1], x2lim[2], knots[2], bdeg[2], pord[2],decom=decom)
      X2 <- MM2$X; Z2<-MM2$Z;	B2<-MM2$B
      c2 <- ncol(B2); d2 <- MM2$d
    } else { X2<-NULL; Z2<-NULL; B2<-NULL; c2<-NULL; d2<-NULL}

    if(!is.null(x3)){
      if(is.null(x3lim))  x3lim <- c(min(x3)-0.01, max(x3)+0.01)
      t <- length(x3)
      if(length(knots)==2) knots <- c(knots,10)
      if(length(bdeg)==2) bdeg <- c(bdeg,3)
      if(length(pord)==2) pord <- c(pord,2)
      MM3 <- MM_basis(x3, x3lim[1], x3lim[2], knots[3], bdeg[3], pord[3],decom=decom)
      X3 <- MM3$X; Z3<-MM3$Z;	B3<-MM3$B
      c3 <- ncol(B3); d3 <- MM3$d
    } else {X3 <- NULL; Z3<-NULL;	B3<-NULL; c3 <- NULL; d3 <- NULL}

    if(ANOVA){
      MM1.2 <- MM_basis(x1, x1lim[1], x1lim[2], knots[1]/div_sp[1],
                        bdeg[1], pord[1],decom=decom)
      X1.2 <- MM1.2$X; Z1.2<-MM1.2$Z; B1.2<-MM1.2$B
      c1.2 <- ncol(B1.2); d1.2 <- MM1.2$d
      MM1.3 <- MM_basis(x1, x1lim[1], x1lim[2], knots[1]/div_sp[2],
                        bdeg[1], pord[1],decom=decom)
      X1.3 <- MM1.3$X; Z1.3<-MM1.3$Z; B1.3<-MM1.3$B
      c1.3 <- ncol(B1.3); d1.3 <- MM1.3$d
      c1 <- c(c1,c1.2,c1.3)
      d1 <- list(d1=d1,d1.2=d1.2,d1.3=d1.3)

      MM2.2 <- MM_basis(x2, x2lim[1], x2lim[2], knots[2]/div_sp[1],
                        bdeg[1], pord[1],decom=decom)
      X2.2 <- MM2.2$X; Z2.2<-MM2.2$Z; B2.2<-MM2.2$B
      c2.2 <- ncol(B2.2); d2.2 <- MM2.2$d
      MM2.3 <- MM_basis(x2, x2lim[1], x2lim[2], knots[2]/div_sp[2],
                        bdeg[1], pord[1],decom=decom)
      X2.3 <- MM2.3$X; Z2.3<-MM2.3$Z; B2.3<-MM2.3$B
      c2.3 <- ncol(B2.3); d2.3 <- MM2.3$d
      c2 <- c(c2,c2.2,c2.3)
      d2 <- list(d2=d2,d2.2=d2.2,d2.3=d2.3)

      MM3.2 <- MM_basis(x3, x3lim[1], x3lim[2], knots[3]/div_tm[1],
                        bdeg[1], pord[1],decom=decom)
      X3.2 <- MM3.2$X; Z3.2<-MM3.2$Z; B3.2<-MM3.2$B
      c3.2 <- ncol(B3.2); d3.2 <- MM3.2$d
      MM3.3 <- MM_basis(x3, x3lim[1], x3lim[2], knots[3]/div_tm[2],
                        bdeg[1], pord[1],decom=decom)
      X3.3 <- MM3.3$X; Z3.3<-MM3.3$Z; B3.3<-MM3.3$B
      c3.3 <- ncol(B3.3); d3.3 <- MM3.3$d
      c3 <- c(c3,c3.2,c3.3)
      d3 <- list(d3=d3,d3.2=d3.2,d3.3=d3.3)

      one1. <- X1[,1,drop=FALSE]
      one2. <- X2[,1,drop=FALSE]
      x1. <- X1[,-1,drop=FALSE]
      x2. <- X2[,-1,drop=FALSE]
      one3. <- X3[,1,drop=FALSE]
      x3. <- X3[,-1,drop=FALSE]
    }

    if(!is.null(x1) && !is.null(x2))
    {
      if(!is.null(x3))
      {
        X.int <- kronecker(Rten2(one1.,one2.),one3.)
        X.f1 <-  kronecker(Rten2(x1.,one2.),one3.)
        X.f2 <-  kronecker(Rten2(one1.,x2.),one3.)
        X.f3 <-  kronecker(Rten2(one1.,one2.),x3.)
        X.f12 <- kronecker(Rten2(x1.,x2.),one3.)
        X.f13 <- kronecker(Rten2(x1.,one2.),x3.)
        X.f23 <- kronecker(Rten2(one1.,x2.),x3.)
        X.f123 <- kronecker(Rten2(x1.,x2.),x3.)
        X_spt <- X.int
        if(f1_main)
        {
          X_spt <- cbind(X_spt,X.f1)
        }
        if(f2_main)
        {
          X_spt <- cbind(X_spt,X.f2)
        }
        if(f3_main)
        {
          X_spt <- cbind(X_spt,X.f3)
        }
        if(f12_int)
        {
          X_spt <- cbind(X_spt,X.f12)
        }
        if(f13_int)
        {
          X_spt <- cbind(X_spt,X.f13)
        }
        if(f23_int)
        {
          X_spt <- cbind(X_spt,X.f23)
        }
        if(f123_int)
        {
          X_spt <- cbind(X_spt,X.f123)
        }
        #X_spt similar a la función fit.sap3D.ANOVA del paquete SAP
        #X_spt <- kronecker(Rten2(X1,X2),X3)
        #X_spt <- kronecker(Rten2(X2,X1),X3)
        # if(ANOVA_part)
        # {
        #   n <- length(x1)
        #   t <- length(x3)
        #   X_spt <- cbind(rep(1,n*t),x1%x%rep(1,t),x2%x%rep(1,t),
        #                  rep(1,n)%x%x3,
        #                  (x1*x2)%x%rep(1,t),
        #                  x1%x%x3,x2%x%x3)
        # }
        if(ANOVA)
        {
          # VIP: DEFINIR LAS MATRICES Z EN FUNCIÓN DE EFECTOS PRINCIPALES
          # E INTERACCIONES
          # VIP**: OJO CON EL ORDEN PARA LUEGO ESTABLECER LA MATRIZ
          # DE PENALIZACIONES Y PARÁMETROS NP.
          Z.f1 <- kronecker(Rten2(Z1, one2.), one3.) # g1u
          Z.f2 <- kronecker(Rten2(one1., Z2), one3.) # g2u
          Z.f3 <- kronecker(Rten2(one1., one2.), Z3) # g3u

          Z.f12 <- cbind(kronecker(Rten2(Z1.2, x2.), one3.),
                         kronecker(Rten2(x1., Z2.2), one3.),
                         kronecker(Rten2(Z1.2, Z2.2), one3.))
          # g12u | g21u | g12b+g21b
          Z.f13 <- cbind(kronecker(Rten2(Z1.2, one2.), x3.),
                         kronecker(Rten2(x1., one2.), Z3.2),
                         kronecker(Rten2(Z1.2, one2.), Z3.2))
          # g13u | g31u | g13b+g31b
          Z.f23 <- cbind(kronecker(Rten2(one1., Z2.2), x3.),
                         kronecker(Rten2(one1., x2.), Z3.2),
                         kronecker(Rten2(one1., Z2.2), Z3.2))
          # g23u | g32u | g23b+g32b
          Z.f123 <- cbind(kronecker(Rten2(Z1.3, x2.), x3.),
                          kronecker(Rten2(x1., Z2.3), x3.),
                          kronecker(Rten2(x1., x2.), Z3.3),
                          kronecker(Rten2(Z1.3, Z2.3), x3.),
                          kronecker(Rten2(Z1.3, x2.), Z3.3),
                          kronecker(Rten2(x1., Z2.3), Z3.3),
                          kronecker(Rten2(Z1.3, Z2.3), Z3.3))
          # g123u | g213u | g321u | g123b+g213b | g132b+g312b |
          # g231b+g321b | g1t+g2t+g3t
# NUEVA DEFINICIÓN DE Z_spt EN FUNCIÓN DE LAS INTERACCIONES INCLUIDAS
          Z_spt <- NULL
          if(f1_main)
          {
            Z_spt <- cbind(Z_spt,Z.f1)
          }
          if(f2_main)
          {
            Z_spt <- cbind(Z_spt,Z.f2)
          }
          if(f3_main)
          {
            Z_spt <- cbind(Z_spt,Z.f3)
          }
          if(f12_int)
          {
            Z_spt <- cbind(Z_spt,Z.f12)
          }
          if(f13_int)
          {
            Z_spt <- cbind(Z_spt,Z.f13)
          }
          if(f23_int)
          {
            Z_spt <- cbind(Z_spt,Z.f23)
          }
          if(f123_int)
          {
            Z_spt <- cbind(Z_spt,Z.f123)
          }

        } else {
          Z_spt <- cbind(kronecker(Rten2(Z1,X2), X3),
                         kronecker(Rten2(X1,Z2), X3),
                         kronecker(Rten2(X1,X2), Z3),
                         kronecker(Rten2(Z1,Z2), X3),
                         kronecker(Rten2(Z1,X2), Z3),
                         kronecker(Rten2(X1,Z2), Z3),
                         kronecker(Rten2(Z1,Z2), Z3))
        }
      } else {
          # HABRÍA QUE MODIFICAR ESTE CASO PARA INCLUIR EL PS-ANOVA SÓLO ESPACIAL
        X_spt <- Rten2(X2, X1)
        Z_spt <- cbind(Rten2(X2, Z1), Rten2(Z2, X1),Rten2(Z2, Z1))
      }
    } else {
      X_spt <- matrix(1,nrow=n,ncol=1) # Sólo se incluye intercepto
      Z_spt<-NULL
    }
    X <- X_spt
    Z <- Z_spt

    nvar_par <- 0; X_par <- NULL
    if(!is.null(Xpar))
    {
      nvar_par <- ncol(Xpar)
      X_par <- as.matrix(Xpar)
      X <- cbind(X_par,X)
    }
    nvar_nopar<-0; X_nopar <- NULL; Z_nopar <- NULL
    d_nopar <- as.list(NULL); c_nopar <- NULL
    if(!is.null(Xnopar)){
      nvar_nopar <- ncol(Xnopar)
      if(is.null(knots_nopar)) knots_nopar <- rep(10,nvar_nopar)
      if(is.null(bdeg_nopar)) bdeg_nopar <- rep(3,nvar_nopar)
      if(is.null(pord_nopar)) pord_nopar <- rep(2,nvar_nopar)
      for(i in 1:nvar_nopar)
      {
        var_np<-Xnopar[,i]
        MM.var <- MM_basis(var_np,min(var_np)-0.01,max(var_np)+0.01,
                   knots_nopar[i],bdeg_nopar[i],pord_nopar[i],decom=1)
        X.var <- MM.var$X; Z.var <- MM.var$Z
        B.var <- MM.var$B; d.var <- MM.var$d
        X_nopar <- cbind(X_nopar,X.var[,c(-1)]) # Quita columna intercepto
        Z_nopar <- cbind(Z_nopar,Z.var)
        d_nopar[[i]] <- d.var
        c_nopar <- c(c_nopar,ncol(B.var))
      }
      X <- cbind(X,X_nopar)
      Z <- cbind(Z,Z_nopar)
    }

    res <- list(X=X,Z=Z,X_spt=X_spt,Z_spt=Z_spt,
                X1=X1,Z1=Z1,c1=c1,X2=X2,Z2=Z2,c2=c2,X3=X3,Z3=Z3,c3=c3,
                d1=d1,d2=d2,d3=d3,
                nvar_par=nvar_par,nvar_nopar=nvar_nopar,
                X_par=X_par,
                X_nopar=X_nopar,Z_nopar=Z_nopar,
                c_nopar=c_nopar,d_nopar=d_nopar)
    res
  }
