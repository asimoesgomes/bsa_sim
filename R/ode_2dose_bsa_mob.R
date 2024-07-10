odin_ode_2dose_bsa_mob_v1 <- odin::odin({
  
  initial(S[,])      <- y0[1,i,j]
  initial(E[,])      <- y0[2,i,j]
  initial(Inbsa[,])      <- y0[3,i,j]
  initial(Ibsa_pre[,])   <- y0[4,i,j]
  initial(Ibsa_post[,])  <- y0[5,i,j]
  initial(Rnbsa[,])      <- y0[6,i,j]
  initial(Rbsa[,])      <- y0[7,i,j]
  initial(R[,])      <- y0[8,i,j]
  initial(Dnbsa[,])      <- y0[9,i,j]
  initial(Dbsa[,])      <- y0[10,i,j]
  initial(D[,])      <- y0[11,i,j]
  initial(P1[,])     <- y0[12,i,j]
  initial(N1[,])     <- y0[13,i,j]
  initial(P2[,])     <- y0[14,i,j]
  initial(N2[,])     <- y0[15,i,j]
  initial(cumV1[,])  <- y0[16,i,j]
  initial(cumV2[,])  <- y0[17,i,j]
  initial(cumV[,])   <- y0[18,i,j]
  initial(cumI[,])   <- y0[19,i,j]
  initial(I[,])  <- y0[20,i,j]
  
  bsa_infec <- delta1[as.integer(t)]
  bsa_death <- delta2[as.integer(t)]
  
  beta_matrix[,,] <- contacts[i,j]*q[as.integer(t),j]*(Ibsa_pre[j,k]+Ibsa_post[j,k]*(1-(bsa_infec-1))+Inbsa[j,k])#k is country
  beta[,] <- sum(beta_matrix[i,,j])#j is country
  
  imig_s[,,] <- mob[k,j]*S[i,k]#I need the matrix of initial values to be proportional to the pop dist across countries
  emig_s[,,] <- mob[j,k]*S[i,j]
  imig_e[,,] <- mob[k,j]*E[i,k]
  emig_e[,,] <- mob[j,k]*E[i,j]
  
  imig_s_agg[,] <- sum(imig_s[i,j,])
  emig_s_agg[,] <- sum(emig_s[i,j,])
  imig_e_agg[,] <- sum(imig_e[i,j,])
  emig_e_agg[,] <- sum(emig_e[i,j,])
  
  va1[] <- 1#1/((1 + exp(ta[i] - t))*(1 + exp(50*(cumV1[i] - vstop[i])))) NOT USING FOR BSA
  
  v1[] <- new_v1[as.integer(t),i]#I need to set number of vaccines per country and relative to the world pop.
  v2[] <- new_v2[as.integer(t),i]
  
  # ODE equations are here:
  deriv(S[,])       <- -(beta[i,j] + constantrisk)*S[i,j] + phi[i]*R[i,j] - va1[i]*v1[j]*S[i,j]/(S[i,j]+R[i,j]*vrf)+imig_s_agg[i,j]-emig_s_agg[i,j]
  deriv(E[,])       <-  (beta[i,j] + constantrisk)*(S[i,j]+N1[i,j]+N2[i,j]) - gamma1[i]*E[i,j]+imig_e_agg[i,j]-emig_e_agg[i,j]
  deriv(Inbsa[,])   <-  (1-pbsa)*gamma1[i]*E[i,j] - gamma2[i]*Inbsa[i,j]
  deriv(Ibsa_pre[,])    <-  pbsa*gamma1[i]*E[i,j] - gammaT[i]*Ibsa_pre[i,j]
  deriv(Ibsa_post[,])   <-  gammaT[i]*Ibsa_pre[i,j] - 1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i,j]
  deriv(Rnbsa[,])       <-  (1-pdeath[i])*gamma2[i]*Inbsa[i,j] - phi[i]*Rnbsa[i,j] - va1[i]*v1[j]*Rnbsa[i,j]*vrf/(S[i,j]+R[i,j])
  deriv(Rbsa[,])       <-  (1-bsa_death*pdeath[i])*(1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i,j]) - phi[i]*Rbsa[i,j] - va1[i]*v1[j]*Rbsa[i,j]*vrf/(S[i,j]+R[i,j])
  deriv(R[,])            <- (1-pdeath[i])*gamma2[i]*Inbsa[i,j] - phi[i]*Rnbsa[i,j] - va1[i]*v1[j]*Rnbsa[i,j]*vrf/(S[i,j]+R[i,j])+
    (1-bsa_death*pdeath[i])*(1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i,j]) - phi[i]*Rbsa[i,j] - va1[i]*v1[j]*Rbsa[i,j]*vrf/(S[i,j]+R[i,j])
  deriv(Dnbsa[,])       <-  pdeath[i]*(gamma2[i]*Inbsa[i,j])
  deriv(Dbsa[,])       <-  bsa_death*pdeath[i]*(1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i,j])
  deriv(D[,])            <- pdeath[i]*(gamma2[i]*Inbsa[i,j])+bsa_death*pdeath[i]*(1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i,j])
  
  #Previous version of protection against transmission
  # deriv(S[,])       <- -(beta[i,j] + constantrisk)*S[i,j] + phi[i]*R[i,j] - va1[i]*v1[j]*S[i,j]/(S[i,j]+R[i,j]*vrf)+imig_s_agg[i,j]-emig_s_agg[i,j]
  # deriv(E[,])       <-  (beta[i,j] + constantrisk)*(S[i,j]+N1[i,j]+N2[i,j]) - gamma1[i]*E[i,j]+imig_e_agg[i,j]-emig_e_agg[i,j]
  # deriv(Inbsa[,])   <-  (1-pbsa)*gamma1[i]*E[i,j] - gamma2[i]*Inbsa[i,j]
  # deriv(Ibsa_pre[,])    <-  pbsa*gamma1[i]*E[i,j] - gammaT[i]*Ibsa_pre[i,j]
  # deriv(Ibsa_post[,])   <-  gammaT[i]*Ibsa_pre[i,j] - bsa_infec*1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i,j]
  # deriv(Rnbsa[,])       <-  (1-pdeath[i])*gamma2[i]*Inbsa[i,j] - phi[i]*Rnbsa[i,j] - va1[i]*v1[j]*Rnbsa[i,j]*vrf/(S[i,j]+R[i,j])
  # deriv(Rbsa[,])       <-  (1-bsa_death*pdeath[i])*(bsa_infec*1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i,j]) - phi[i]*Rbsa[i,j] - va1[i]*v1[j]*Rbsa[i,j]*vrf/(S[i,j]+R[i,j])
  # deriv(R[,])            <- (1-pdeath[i])*gamma2[i]*Inbsa[i,j] - phi[i]*Rnbsa[i,j] - va1[i]*v1[j]*Rnbsa[i,j]*vrf/(S[i,j]+R[i,j])+
  #   (1-bsa_death*pdeath[i])*(bsa_infec*1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i,j]) - phi[i]*Rbsa[i,j] - va1[i]*v1[j]*Rbsa[i,j]*vrf/(S[i,j]+R[i,j])
  # deriv(Dnbsa[,])       <-  pdeath[i]*(gamma2[i]*Inbsa[i,j])
  # deriv(Dbsa[,])       <-  bsa_death*pdeath[i]*(bsa_infec*1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i,j])
  # deriv(D[,])            <- pdeath[i]*(gamma2[i]*Inbsa[i,j])+bsa_death*pdeath[i]*(bsa_infec*1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i,j])
  
  deriv(P1[,])      <-  va1[i]*v1[j]*(e1[i]*S[i,j]+R[i,j]*vrf)/(S[i,j]+R[i,j]*vrf) - v2[j]*P1[i,j] - kappa1[i]*P1[i,j]
  deriv(N1[,])      <-  va1[i]*v1[j]*(1-e1[i])*S[i,j]/(S[i,j]+R[i,j]*vrf) - v2[j]*N1[i,j] - (beta[i,j] + constantrisk)*N1[i,j] + kappa1[i]*P1[i,j]
  deriv(P2[,])      <-  e2[i]*v2[j]*(P1[i,j]+N1[i,j]) - kappa2[i]*P2[i,j]
  deriv(N2[,])      <-  (1-e2[i])*v2[j]*(P1[i,j]+N1[i,j]) + kappa2[i]*P2[i,j] - (beta[i,j] + constantrisk)*N2[i,j] 
  
  deriv(cumV1[,])   <-  va1[i]*v1[j]
  deriv(cumV2[,])   <-  va1[i]*v2[j]
  deriv(cumV[,])    <-  va1[i]*v1[j] + va1[i]*v2[j]
  deriv(cumI[,])    <-  gamma1[i]*E[i,j]
  deriv(I[,]) <- gamma1[i]*E[i,j] - gamma2[i]*Inbsa[i,j] - 1/(1/gamma2[i]-1/gammaT[i])*Ibsa_post[i,j]
  
  # PARAMETERS: general
  Ngroups          <- user()
  Ncountries       <- user()
  Ndays            <- user()
  Nc               <- user()
  pbsa             <- user()
  constantrisk     <- user()
  vrf              <- user()
  delta1[]           <- user()
  delta2[]           <- user()
  y0[,,]            <- user()
  contacts[,]      <- user()
  vstop[]          <- user()
  ta[]             <- user()
  new_v1[,]         <- user()
  new_v2[,]         <- user()
  full_period[]    <- user()
  q[,]              <- user()
  mob[,]            <- user()
  
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
  dim(v1) <- Ncountries
  dim(v2) <- Ncountries
  dim(imig_s) <- c(Ngroups, Ngroups, Ncountries)
  dim(imig_e) <- c(Ngroups, Ngroups, Ncountries)
  dim(emig_s) <- c(Ngroups, Ngroups, Ncountries)
  dim(emig_e) <- c(Ngroups, Ngroups, Ncountries)
  dim(imig_s_agg)  <- c(Ngroups, Ncountries)
  dim(emig_s_agg)  <- c(Ngroups, Ncountries)
  dim(imig_e_agg)  <- c(Ngroups, Ncountries)
  dim(emig_e_agg)  <- c(Ngroups, Ncountries)
  dim(q)           <- c(Ndays, Ngroups)
  dim(mob)          <- c(Ncountries, Ncountries)
  dim(beta)        <- c(Ngroups, Ncountries)
  dim(beta_matrix) <- c(Ngroups, Ngroups, Ncountries)
  dim(y0)          <- c(Nc, Ngroups, Ncountries)
  dim(contacts)    <- c(Ngroups, Ngroups)
  
  dim(S)    <- c(Ngroups,Ncountries)
  dim(E)    <- c(Ngroups,Ncountries)
  dim(Inbsa)    <- c(Ngroups,Ncountries)
  dim(Ibsa_pre) <- c(Ngroups,Ncountries)
  dim(Ibsa_post) <- c(Ngroups,Ncountries)
  dim(Rnbsa)    <- c(Ngroups,Ncountries)
  dim(Rbsa)    <- c(Ngroups,Ncountries)
  dim(Dnbsa)    <- c(Ngroups,Ncountries)
  dim(Dbsa)    <- c(Ngroups,Ncountries)
  dim(D)    <- c(Ngroups,Ncountries)
  dim(R)    <- c(Ngroups,Ncountries)
  dim(N1)   <- c(Ngroups,Ncountries)
  dim(N2)   <- c(Ngroups,Ncountries)
  dim(P1)   <- c(Ngroups,Ncountries)
  dim(P2)   <- c(Ngroups,Ncountries)
  dim(cumV) <- c(Ngroups,Ncountries)
  dim(cumV1) <- c(Ngroups,Ncountries)
  dim(cumV2) <- c(Ngroups,Ncountries)
  dim(cumI) <- c(Ngroups,Ncountries)
  dim(I)    <- c(Ngroups,Ncountries)
  
  dim(new_v1)      <- c(Ndays,Ncountries)
  dim(new_v2)      <- c(Ndays,Ncountries)
  dim(full_period) <- Ndays
})
