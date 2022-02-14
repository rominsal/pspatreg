build_Omega_ar1 <- function(phi, nt) {
    # The noise of the model is assumed to follow N(0,Var(U)) distribution.  It is supposed
  # that Var(U)=I_n %x% sigma_u^2*Omega where sigma_u^2= sigma_eps^2/(1-phi^2).  Inputs: phi
  # parameter for AR1 time process, time for sample temporal size Ouput: Omega=Matrix of
  # Var-Cov in time for the noise (without sigma_u^2), Omegainv=Omega^(-1),
  # der_sig2u_Omega_phi=\frac{\partial \sigma_u^2*Omega}{\partial phi},
  # der2_sig2u_Omega_phi=\frac{\partial^2 \sigma_u^2*Omega}{\partial phi^2} Dependencies:
  # library Matrix

  diags_Omega <- list(rep(1, nt))
  diags_der_sig2u_Omega_phi <- list(rep(2 * phi, nt),
                                    rep(1 + phi^2, nt),
                                    rep(2 * phi,nt))
  for (l in 2:nt) {
    diags_Omega[[l]] <- rep(phi^(l - 1), nt)
    if (l > 3)
      diags_der_sig2u_Omega_phi[[l]] <- rep((l - 1) * phi^(l - 2) - (l - 3) * phi^l,
                                            nt)
  }
  Omega <- as.matrix(bandSparse(n = nt, m = nt, k = -c(0:(nt - 1)),
                                 diagonals = diags_Omega,
                                 symmetric = TRUE))
  # Assumption nt>=3
  diags_Omegainv <- list(c(1, rep(1 + phi^2, nt - 2), 1), rep(-phi, nt - 1))
  Omegainv <- 1/(1 - phi^2) * bandSparse(n = nt, m = nt, k = -c(0:1),
                                         diagonals = diags_Omegainv,
                                         symmetric = TRUE, giveCsparse = TRUE)
  # First derivative of Omega with respect to phi
  der_sig2u_Omega_phi <- 1/((1 - phi^2)^2) *
                         as(bandSparse(n = nt, m = nt, k = -c(0:(nt - 1)),
                                               diagonals = diags_der_sig2u_Omega_phi,
                                               symmetric = TRUE), "dsyMatrix")
  # Second derivative of Omega with respect to phi
  part1_der2_sig2u_Omega_phi <- (4 * phi/((1 - phi^2)^3)) *
                                as(bandSparse(n = nt, m = nt,
                                                      k = -c(0:(nt - 1)),
                                                      diagonals = diags_der_sig2u_Omega_phi,
                                                      symmetric = TRUE), "dsyMatrix")
  diags_part2_der2_sig2u_Omega_phi <- list(rep(2, nt), rep(2 * phi, nt), rep(2, nt))
  for (l in 3:nt) {
    diags_part2_der2_sig2u_Omega_phi[[l]] <- rep((l - 1) * (l - 2) *
                                                   phi^(l - 3) - (l - 3) * l *
                                                   phi^(l - 1), nt)
  }
  part2_der2_sig2u_Omega_phi <- 1/((1 - phi^2)^2) *
                        as(bandSparse(n = nt, m = nt, k = -c(0:(nt - 1)),
                                              diagonals = diags_part2_der2_sig2u_Omega_phi,
                                              symmetric = TRUE), "dsyMatrix")

  der2_sig2u_Omega_phi <- part1_der2_sig2u_Omega_phi + part2_der2_sig2u_Omega_phi

  res <- list(Omega = as.matrix(Omega), Omegainv = as.matrix(Omegainv),
              der_sig2u_Omega_phi = as.matrix(der_sig2u_Omega_phi),
              der2_sig2u_Omega_phi = as.matrix(der2_sig2u_Omega_phi))
  return(res)
}
