odin_ode_2dose_bsa_v1 <- odin::odin({
  
  initial(S[])      <- y0[1,i]
  initial(E[])      <- y0[2,i]
  initial(Inbsa[])      <- y0[3,i]
  initial(Ibsa_pre[])   <- y0[4,i]
  initial(Ibsa_post[])  <- y0[5,i]
  initial(Rnbsa[])      <- y0[6,i]
  initial(Rbsa[])      <- y0[7,i]
  initial(R[])      <- y0[8,i]
  initial(Dnbsa[])      <- y0[9,i]
  initial(Dbsa[])      <- y0[10,i]
  initial(D[])      <- y0[11,i]
  initial(P1[])     <- y0[12,i]
  initial(N1[])     <- y0[13,i]
  initial(P2[])     <- y0[14,i]
  initial(N2[])     <- y0[15,i]
  initial(cumV1[])  <- y0[16,i]
  initial(cumV2[])  <- y0[17,i]
  initial(cumV[])   <- y0[18,i]
  initial(cumI[])   <- y0[19,i]
  initial(I[])  <- y0[20,i]
  
  beta_matrix[,] <- contacts[i,j]*q[as.integer(t),j]*(Ibsa_pre[j]+Ibsa_post[j]+Inbsa[j])
  beta[] <- sum(beta_matrix[i,])
  
  va1[] <- 1/((1 + exp(ta[i] - t))*(1 + exp(50*(cumV1[i] - vstop[i]))))
  
  v1 <- new_v1[as.integer(t)]
  v2 <- new_v2[as.integer(t)]
  
  bsa_infec <- delta1[as.integer(t)]
  bsa_death <- delta2[as.integer(t)]
  
  # ODE equations are here:
  deriv(S[])       <- -(beta[i] + constantrisk)*S[i] + phi[i]*R[i] - va1[i]*v1*S[i]/(S[i]+R[i]*vrf)
  deriv(E[])       <-  (beta[i] + constantrisk)*(S[i]+N1[i]+N2[i]) - gamma1[i]*E[i]
  deriv(Inbsa[])   <-  (1-pbsa)*gamma1[i]*E[i] - gamma2[i]*Inbsa[i]
  deriv(Ibsa_pre[])    <-  pbsa*gamma1[i]*E[i] - gammaT[i]*Ibsa_pre[i]
  deriv(Ibsa_post[])   <-  gammaT[i]*Ibsa_pre[i] - bsa_infec*1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i]
  deriv(Rnbsa[])       <-  (1-pdeath[i])*gamma2[i]*Inbsa[i] - phi[i]*Rnbsa[i] - va1[i]*v1*Rnbsa[i]*vrf/(S[i]+R[i])
  deriv(Rbsa[])       <-  (1-bsa_death*pdeath[i])*(bsa_infec*1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i]) - phi[i]*Rbsa[i] - va1[i]*v1*Rbsa[i]*vrf/(S[i]+R[i])
  deriv(R[])            <- (1-pdeath[i])*gamma2[i]*Inbsa[i] - phi[i]*Rnbsa[i] - va1[i]*v1*Rnbsa[i]*vrf/(S[i]+R[i])+
    (1-bsa_death*pdeath[i])*(bsa_infec*1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i]) - phi[i]*Rbsa[i] - va1[i]*v1*Rbsa[i]*vrf/(S[i]+R[i])
  deriv(Dnbsa[])       <-  pdeath[i]*(gamma2[i]*Inbsa[i])
  deriv(Dbsa[])       <-  bsa_death*pdeath[i]*(bsa_infec*1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i])
  deriv(D[])            <- pdeath[i]*(gamma2[i]*Inbsa[i])+bsa_death*pdeath[i]*(bsa_infec*1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i])
  
  deriv(P1[])      <-  va1[i]*v1*(e1[i]*S[i]+R[i]*vrf)/(S[i]+R[i]*vrf) - v2*P1[i] - kappa1[i]*P1[i]
  deriv(N1[])      <-  va1[i]*v1*(1-e1[i])*S[i]/(S[i]+R[i]*vrf) - v2*N1[i] - (beta[i] + constantrisk)*N1[i] + kappa1[i]*P1[i]
  deriv(P2[])      <-  e2[i]*v2*(P1[i]+N1[i]) - kappa2[i]*P2[i]
  deriv(N2[])      <-  (1-e2[i])*v2*(P1[i]+N1[i]) + kappa2[i]*P2[i] - (beta[i] + constantrisk)*N2[i] 
  
  deriv(cumV1[])   <-  va1[i]*v1
  deriv(cumV2[])   <-  va1[i]*v2
  deriv(cumV[])    <-  va1[i]*v1 + va1[i]*v2
  deriv(cumI[])    <-  gamma1[i]*E[i]
  deriv(I[]) <- gamma1[i]*E[i] - gamma2[i]*Inbsa[i] - bsa_infec*1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i]
  
  # PARAMETERS: general
  Ngroups          <- user()
  Ndays            <- user()
  Nc               <- user()
  pbsa             <- user()
  constantrisk     <- user()
  vrf              <- user()
  delta1[]           <- user()
  delta2[]           <- user()
  y0[,]            <- user()
  contacts[,]      <- user()
  vstop[]          <- user()
  ta[]             <- user()
  new_v1[]         <- user()
  new_v2[]         <- user()
  full_period[]    <- user()
  q[,]              <- user()
  
  # Parameters: ODE system
  gamma1[]       <- user()
  gamma2[]       <- user()
  gammaT[]       <- user()
  kappa1[]       <- user()
  kappa2[]       <- user()
  e1[]           <- user()
  e2[]           <- user()
  phi[]          <- user()
  pdeath[]       <- user()
  
  # Parameters: ODE system
  dim(e1)           <- Ngroups
  dim(e2)           <- Ngroups
  dim(gamma1)       <- Ngroups
  dim(gamma2)       <- Ngroups
  dim(gammaT)       <- Ngroups
  dim(kappa1)       <- Ngroups
  dim(kappa2)       <- Ngroups
  dim(phi)          <- Ngroups
  dim(pdeath)       <- Ngroups
  
  dim(vstop)       <- Ngroups
  dim(ta)          <- Ngroups
  dim(delta1)      <- Ndays
  dim(delta2)      <- Ndays
  dim(va1)         <- Ngroups
  dim(q)           <- c(Ndays, Ngroups)
  dim(beta)        <- Ngroups
  dim(beta_matrix) <- c(Ngroups, Ngroups)
  dim(y0)          <- c(Nc, Ngroups)
  dim(contacts)    <- c(Ngroups, Ngroups)
  
  dim(S)    <- Ngroups
  dim(E)    <- Ngroups
  dim(Inbsa)    <- Ngroups
  dim(Ibsa_pre) <- Ngroups
  dim(Ibsa_post) <- Ngroups
  dim(Rnbsa)    <- Ngroups
  dim(Rbsa)    <- Ngroups
  dim(Dnbsa)    <- Ngroups
  dim(Dbsa)    <- Ngroups
  dim(D)    <- Ngroups
  dim(R)    <- Ngroups
  dim(N1)   <- Ngroups
  dim(N2)   <- Ngroups
  dim(P1)   <- Ngroups
  dim(P2)   <- Ngroups
  dim(cumV) <- Ngroups
  dim(cumV1) <- Ngroups
  dim(cumV2) <- Ngroups
  dim(cumI) <- Ngroups
  dim(I)    <- Ngroups
  
  dim(new_v1)      <- Ndays
  dim(new_v2)      <- Ndays
  dim(full_period) <- Ndays
})
