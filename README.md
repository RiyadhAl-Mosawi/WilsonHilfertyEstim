  To use this package you need to specify the following
  
  par=true parameter vector
  
  R=censoring scheme
  
  m=no. of failure times
  
  n=sample size
  
  k=size of group
  
  L=lower specification limit
  
  U=upper specification limit
  
  P0=minimum permitted proportion of conformance
  
  sim_no=number of replications
  
  boot=size of boot sample
  
  MCMC_size=size of MCMC sample
  
  Burn_in=size of burn-in sample
  
  q=parameter of general entropy loss function
  
  c=parameter of LINEX loss function
  
  randomseed= seed number
  
  dir=location for the results to be saved
  
  To compute the MLE and MPS along with their confidence inetrval, write
  
  calssical(par,m,n,k,R,L,U,P0,sim_no,boot,MCMC_size,Burn_in,q,c,randomseed,dir)
  
  and to compute the Mayes estimates, write
  
  Bayes(par,m,n,k,R,L,U,P0,sim_no,MCMC_size,Burn_in,q,c,randomseed,dir)

