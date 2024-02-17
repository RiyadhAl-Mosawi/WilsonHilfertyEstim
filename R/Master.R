# Compute the MSE with respect to the loss function
MSE.fun=function(eta,eta.hat, type="SEL",q){
  if(type=="SEL") return(mean(apply(as.matrix(1:length(eta.hat)),1,function(i) (eta.hat[i]-eta)^2)))
  if(type=="GEL") return(mean(apply(as.matrix(1:length(eta.hat)),1,function(i) (eta.hat[i]/eta)^q-q*log(eta.hat[i]/eta)-1)))
  if(type=="LINEX") return(mean(apply(as.matrix(1:length(eta.hat)),1,function(i) exp(q*(eta.hat[i]-eta))-q*(eta.hat[i]-eta)-1)))
}

Classical=function(para=par,m,n,k,R,L,U,P0,sim_no,boot,
             MCMC_size,Burn_in,q,c,randomseed,dir){
  cpc_true=cpc_fun(para,L,U,P0)
  para=c(para,cpc_true)
  it=1
  set.seed(randomseed)
  while(it<=sim_no){
    if(it==1){
      MLE_alp=MLE_lam=MLE_cpc=NULL
      MPS_alp=MPS_lam=MPS_cpc=NULL
      Boot_MLE1_alp=Boot_MLE1_lam=Boot_MLE1_cpc=NULL
      Boot_MLE2_alp=Boot_MLE2_lam=Boot_MLE2_cpc=NULL
      Boot_MPS1_alp=Boot_MPS1_lam=Boot_MPS1_cpc=NULL
      Boot_MPS2_alp=Boot_MPS2_lam=Boot_MPS2_cpc=NULL
      #TK_MLE_alp=TK_MLE_lam=TK_MLE_cpc=NULL
      #TK_MPS_alp=TK_MPS_lam=TK_MPS_cpc=NULL
      #MH_MLE_alp=MH_MLE_lam=MH_MLE_cpc=NULL
      #MH_MPS_alp=MH_MPS_lam=MH_MPS_cpc=NULL
      
      start_time=date()
    }  
    
    X=try(GenDataRcpp(c(0.1,100),para[1:2],R,k),silent=TRUE)
    if(is.character(X)) next
    if(X[m]==100) X[m]=X[m-1]+0.1
    
    #################################################################
    #------------- MLE using likelihood function
    ############s#####################################################
    
    mle_res=suppressWarnings(try(Estim(para[1:2],X,R,L,U,P0,k,"LF"),silent=TRUE))
    # Check convergence
    if(is.character(mle_res)) next
    if((mle_res$Estim_para[1]<=0)|(mle_res$Estim_para[1]>=3*para[1])) next
    if((mle_res$Estim_para[2]<=0)|(mle_res$Estim_para[2]>=3*para[2])) next
    if((mle_res$Estim_Cpc<=0)|(mle_res$Estim_Cpc>=3*para[3])) next
    
    mle_est=c(as.numeric(mle_res$Estim_para),as.numeric(mle_res$Estim_Cpc))
    mle_ese=c(as.numeric(mle_res$SE_para),as.numeric(mle_res$SE_Cpc))
    mle_aci=matrix(nr=3,nc=3,NA)
    mle_aci[1,1:2]=c(max(0,mle_res$CI_para[1]),mle_res$CI_para[2])
    mle_aci[1,3]=ifelse((mle_aci[1,2]>=para[1])&(mle_aci[1,1]<=para[1]),1,0)
    mle_aci[2,1:2]=c(max(0,mle_res$CI_para[3]),mle_res$CI_para[4])
    mle_aci[2,3]=ifelse((mle_aci[2,2]>=para[2])&(mle_aci[2,1]<=para[2]),1,0)
    mle_aci[3,1:2]=c(max(0,mle_res$CI_Cpc[1]),mle_res$CI_Cpc[2])
    mle_aci[3,3]=ifelse((mle_aci[3,2]>=para[3])&(mle_aci[3,1]<=para[3]),1,0)
    
    #################################################################
    #------------- MPS using product spacing function
    ############s#####################################################
    
    mps_res=suppressWarnings(try(Estim(para[1:2],X,R,L,U,P0,k,"SPF"),silent=TRUE))

    # Check convergence
    if(is.character(mps_res)) next
    if((mps_res$Estim_para[1]<=0)|(mps_res$Estim_para[1]>=3*para[1])) next
    if((mps_res$Estim_para[2]<=0)|(mps_res$Estim_para[2]>=3*para[2])) next
    if((mps_res$Estim_Cpc<=0)|(mps_res$Estim_Cpc>=3*para[3])) next
    
    mps_est=c(as.numeric(mps_res$Estim_para),as.numeric(mps_res$Estim_Cpc))
    mps_ese=c(as.numeric(mps_res$SE_para),as.numeric(mps_res$SE_Cpc))
    mps_aci=matrix(nr=3,nc=3,NA)
    mps_aci[1,1:2]=c(max(0,mps_res$CI_para[1]),mps_res$CI_para[2])
    mps_aci[1,3]=ifelse((mps_aci[1,2]>=para[1])&(mps_aci[1,1]<=para[1]),1,0)
    mps_aci[2,1:2]=c(max(0,mps_res$CI_para[3]),mps_res$CI_para[4])
    mps_aci[2,3]=ifelse((mps_aci[2,2]>=para[2])&(mps_aci[2,1]<=para[2]),1,0)
    mps_aci[3,1:2]=c(max(0,mps_res$CI_Cpc[1]),mps_res$CI_Cpc[2])
    mps_aci[3,3]=ifelse((mps_aci[3,2]>=para[3])&(mps_aci[3,1]<=para[3]),1,0)
    
    ## Save the results in Matrices
    MLE_alp=rbind(MLE_alp,c(mle_est[1],(mle_est[1]-para[1])^2,mle_ese[1],mle_aci[1,]))
    MLE_lam=rbind(MLE_lam,c(mle_est[2],(mle_est[2]-para[2])^2,mle_ese[2],mle_aci[2,]))
    MLE_cpc=rbind(MLE_cpc,c(mle_est[3],(mle_est[3]-para[3])^2,mle_ese[3],mle_aci[3,]))
  
    MPS_alp=rbind(MPS_alp,c(mps_est[1],(mps_est[1]-para[1])^2,mps_ese[1],mps_aci[1,]))
    MPS_lam=rbind(MPS_lam,c(mps_est[2],(mps_est[2]-para[2])^2,mps_ese[2],mps_aci[2,]))
    MPS_cpc=rbind(MPS_cpc,c(mps_est[3],(mps_est[3]-para[3])^2,mps_ese[3],mps_aci[3,]))
    
    cat("====== Scheme ",sch,"=== Iteration no.", it," ===========================\n")

    #---- Bootstrap method
    MLE_Boot=MPS_Boot=NULL
    j=1
    #pb = txtProgressBar(min = 0, max = boot, initial = 0, width=50,char="*", style =2)
    pb = progress_bar$new(format = "Bootstrap estimates :: (:spin) :percent [Elapsed time: :elapsedfull] ",total = boot,complete = "=",   incomplete = "-", current = ">",   clear = TRUE, width = 100)
    while(j<=boot){
      Y=try(GenDataRcpp(c(0.1,100),mle_est,R,k),silent=TRUE)
      if(is.character(Y)) next
      if(Y[m]==100) Y[m]=Y[m-1]+0.1
      
      mle_res1=suppressWarnings(try(Estim(para[1:2],Y,R,L,U,P0,k,"LF"),silent=TRUE))
      
      if(is.character(mle_res1)) next
      if((mle_res1$Estim_para[1]<=0)|(mle_res1$Estim_para[1]>=3*para[1])) next
      if((mle_res1$Estim_para[2]<=0)|(mle_res1$Estim_para[2]>=3*para[2])) next
      if((mle_res1$Estim_Cpc<=0)|(mle_res1$Estim_Cpc>=3*para[3])) next
      
      MLE_Boot=rbind(MLE_Boot,c(mle_res1$Estim_para,mle_res1$Estim_Cpc))
      
      mps_res1=suppressWarnings(try(Estim(para[1:2],Y,R,L,U,P0,k,"SPF"),silent=TRUE))

      if(is.character(mps_res1)) next
      if((mps_res1$Estim_para[1]<=0)|(mps_res1$Estim_para[1]>=3*para[1])) next
      if((mps_res1$Estim_para[2]<=0)|(mps_res1$Estim_para[2]>=3*para[2])) next
      if((mps_res1$Estim_Cpc<=0)|(mps_res1$Estim_Cpc>=3*para[3])) next
      
      MPS_Boot=rbind(MPS_Boot,c(mps_res1$Estim_para,mps_res1$Estim_Cpc))
      
      j=j+1
      #setTxtProgressBar(pb,j)
      pb$tick()
    }

    Boot_MLE1_alp=rbind(Boot_MLE1_alp,c(max(0,as.numeric(quantile(MLE_Boot[,1],0.025))),
                                        as.numeric(quantile(MLE_Boot[,1],0.975)),
                                        (as.numeric(quantile(MLE_Boot[,1],0.025))<para[1])*(para[1]<as.numeric(quantile(MLE_Boot[,1],0.975)))))
    Boot_MLE1_lam=rbind(Boot_MLE1_lam,c(max(0,as.numeric(quantile(MLE_Boot[,2],0.025))),
                                        as.numeric(quantile(MLE_Boot[,2],0.975)),
                                        (as.numeric(quantile(MLE_Boot[,2],0.025))<para[2])*(para[2]<as.numeric(quantile(MLE_Boot[,2],0.975)))))
    Boot_MLE1_cpc=rbind(Boot_MLE1_cpc,c(max(0,as.numeric(quantile(MLE_Boot[,3],0.025))),
                                        as.numeric(quantile(MLE_Boot[,3],0.975)),
                                        (as.numeric(quantile(MLE_Boot[,3],0.025))<para[3])*(para[3]<as.numeric(quantile(MLE_Boot[,3],0.975)))))

    Boot_MLE2_alp=rbind(Boot_MLE2_alp,c(max(0,mean(MLE_Boot[,1])-1.96*sd(MLE_Boot[,1])),
                                        mean(MLE_Boot[,1])+1.96*sd(MLE_Boot[,1]),
                                        ((mean(MLE_Boot[,1])-1.96*sd(MLE_Boot[,1]))<para[1])*(para[1]<(mean(MLE_Boot[,1])+1.96*sd(MLE_Boot[,1])))))
    Boot_MLE2_lam=rbind(Boot_MLE2_lam,c(max(0,mean(MLE_Boot[,2])-1.96*sd(MLE_Boot[,2])),
                                        mean(MLE_Boot[,2])+1.96*sd(MLE_Boot[,2]),
                                        ((mean(MLE_Boot[,2])-1.96*sd(MLE_Boot[,2]))<para[2])*(para[2]<(mean(MLE_Boot[,2])+1.96*sd(MLE_Boot[,2])))))
    Boot_MLE2_cpc=rbind(Boot_MLE2_cpc,c(max(0,mean(MLE_Boot[,3])-1.96*sd(MLE_Boot[,3])),
                                        mean(MLE_Boot[,3])+1.96*sd(MLE_Boot[,3]),
                                        ((mean(MLE_Boot[,3])-1.96*sd(MLE_Boot[,3]))<para[3])*(para[3]<(mean(MLE_Boot[,3])+1.96*sd(MLE_Boot[,3])))))
    
    Boot_MPS1_alp=rbind(Boot_MPS1_alp,c(as.numeric(quantile(MPS_Boot[,1],0.025)),
                                        as.numeric(quantile(MPS_Boot[,1],0.975)),
                                        (as.numeric(quantile(MPS_Boot[,1],0.025))<para[1])*(para[1]<as.numeric(quantile(MPS_Boot[,1],0.975)))))
    Boot_MPS1_lam=rbind(Boot_MPS1_lam,c(as.numeric(quantile(MPS_Boot[,2],0.025)),
                                        as.numeric(quantile(MPS_Boot[,2],0.975)),
                                        (as.numeric(quantile(MPS_Boot[,2],0.025))<para[2])*(para[2]<as.numeric(quantile(MPS_Boot[,2],0.975)))))
    Boot_MPS1_cpc=rbind(Boot_MPS1_cpc,c(as.numeric(quantile(MPS_Boot[,3],0.025)),
                                        as.numeric(quantile(MPS_Boot[,3],0.975)),
                                        (as.numeric(quantile(MPS_Boot[,3],0.025))<para[3])*(para[3]<as.numeric(quantile(MPS_Boot[,3],0.975)))))
    
    Boot_MPS2_alp=rbind(Boot_MPS2_alp,c(max(0,mean(MPS_Boot[,1])-1.96*sd(MPS_Boot[,1])),
                                        mean(MPS_Boot[,1])+1.96*sd(MPS_Boot[,1]),
                                        ((mean(MPS_Boot[,1])-1.96*sd(MPS_Boot[,1]))<para[1])*(para[1]<(mean(MPS_Boot[,1])+1.96*sd(MPS_Boot[,1])))))
    Boot_MPS2_lam=rbind(Boot_MPS2_lam,c(max(0,mean(MPS_Boot[,2])-1.96*sd(MPS_Boot[,2])),
                                        mean(MPS_Boot[,2])+1.96*sd(MPS_Boot[,2]),
                                        ((mean(MPS_Boot[,2])-1.96*sd(MPS_Boot[,2]))<para[2])*(para[2]<(mean(MPS_Boot[,2])+1.96*sd(MPS_Boot[,2])))))
    Boot_MPS2_cpc=rbind(Boot_MPS2_cpc,c(max(0,mean(MPS_Boot[,3])-1.96*sd(MPS_Boot[,3])),
                                        mean(MPS_Boot[,3])+1.96*sd(MPS_Boot[,3]),
                                        ((mean(MPS_Boot[,3])-1.96*sd(MPS_Boot[,3]))<para[3])*(para[3]<(mean(MPS_Boot[,3])+1.96*sd(MPS_Boot[,3])))))
    
    ss=NULL
    ss=rbind(ss,c(mean(MLE_alp[,1]),mean(MLE_alp[,2]),mean(MLE_alp[,3]),mean(MLE_alp[,5])-mean(MLE_alp[,4]),mean(MLE_alp[,6]),
                  mean(MLE_lam[,1]),mean(MLE_lam[,2]),mean(MLE_lam[,3]),mean(MLE_lam[,5])-mean(MLE_lam[,4]),mean(MLE_lam[,6]),
                  mean(MLE_cpc[,1]),mean(MLE_cpc[,2]),mean(MLE_cpc[,3]),mean(MLE_lam[,5])-mean(MLE_lam[,4]),mean(MLE_lam[,6])))
    ss=rbind(ss,c(mean(MPS_alp[,1]),mean(MPS_alp[,2]),mean(MPS_alp[,3]),mean(MPS_alp[,5])-mean(MPS_alp[,4]),mean(MPS_alp[,6]),
                  mean(MPS_lam[,1]),mean(MPS_lam[,2]),mean(MPS_lam[,3]),mean(MPS_lam[,5])-mean(MPS_lam[,4]),mean(MPS_lam[,6]),
                  mean(MPS_cpc[,1]),mean(MPS_cpc[,2]),mean(MPS_cpc[,3]),mean(MPS_lam[,5])-mean(MPS_lam[,4]),mean(MPS_lam[,6])))
    ss=rbind(ss,c(NA,NA,NA,mean(Boot_MLE1_alp[,2])-mean(Boot_MLE1_alp[,1]),mean(Boot_MLE1_alp[,3]),
                  NA,NA,NA,mean(Boot_MLE1_lam[,2])-mean(Boot_MLE1_lam[,1]),mean(Boot_MLE1_lam[,3]),
                  NA,NA,NA,mean(Boot_MLE1_cpc[,2])-mean(Boot_MLE1_cpc[,1]),mean(Boot_MLE1_cpc[,3])))
    ss=rbind(ss,c(NA,NA,NA,mean(Boot_MPS1_alp[,2])-mean(Boot_MPS1_alp[,1]),mean(Boot_MPS1_alp[,3]),
                  NA,NA,NA,mean(Boot_MPS1_lam[,2])-mean(Boot_MPS1_lam[,1]),mean(Boot_MPS1_lam[,3]),
                  NA,NA,NA,mean(Boot_MPS1_cpc[,2])-mean(Boot_MPS1_cpc[,1]),mean(Boot_MPS1_cpc[,3])))
    ss=rbind(ss,c(NA,NA,NA,mean(Boot_MLE2_alp[,2])-mean(Boot_MLE2_alp[,1]),mean(Boot_MLE2_alp[,3]),
                  NA,NA,NA,mean(Boot_MLE2_lam[,2])-mean(Boot_MLE2_lam[,1]),mean(Boot_MLE2_lam[,3]),
                  NA,NA,NA,mean(Boot_MLE2_cpc[,2])-mean(Boot_MLE2_cpc[,1]),mean(Boot_MLE2_cpc[,3])))
    ss=rbind(ss,c(NA,NA,NA,mean(Boot_MPS2_alp[,2])-mean(Boot_MPS2_alp[,1]),mean(Boot_MPS2_alp[,3]),
                  NA,NA,NA,mean(Boot_MPS2_lam[,2])-mean(Boot_MPS2_lam[,1]),mean(Boot_MPS2_lam[,3]),
                  NA,NA,NA,mean(Boot_MPS2_cpc[,2])-mean(Boot_MPS2_cpc[,1]),mean(Boot_MPS2_cpc[,3])))
    rownames(ss)=c("MLE","MPS","MLE.BT1","MPS.BT1","MLE.BT2","MPS.BT2")
    colnames(ss)=c("Est.a","MSE.a","ESE.a","Len.a","CP.a","Est.l","MSE.l","ESE.l","Len.l","CP.l","Est.c","MSE.c","ESE.c","Len.c","CP.c")
    print(ss)
    it=it+1
  } 
  
    return(list(MLE_alp=MLE_alp,MLE_lam=MLE_lam,MLE_cpc=MLE_cpc,
                MPS_alp=MPS_alp,MPS_lam=MPS_lam,MPS_cpc=MPS_cpc,
                Boot_MLE1_alp=Boot_MLE1_alp,Boot_MLE1_lam=Boot_MLE1_lam,Boot_MLE1_cpc=Boot_MLE1_cpc,
                Boot_MLE2_alp=Boot_MLE2_alp,Boot_MLE2_lam=Boot_MLE2_lam,Boot_MLE2_cpc=Boot_MLE2_cpc,
                Boot_MPS1_alp=Boot_MPS1_alp,Boot_MPS1_lam=Boot_MPS1_lam,Boot_MPS1_cpc=Boot_MPS1_cpc,
                Boot_MPS2_alp=Boot_MPS2_alp,Boot_MPS2_lam=Boot_MPS2_lam,Boot_MPS2_cpc=Boot_MPS2_cpc))
}

Bayes=function(para=par,m,n,k,R,L,U,P0,sim_no,
             MCMC_size,Burn_in,q,c,randomseed,dir){
  cpc_true=cpc_fun(para,L,U,P0)
  para=c(para,cpc_true)
  it=1
  set.seed(randomseed)
  while(it<=sim_no){
    if(it==1){
      TK_MLE_alp=TK_MLE_lam=TK_MLE_cpc=NULL
      TK_MPS_alp=TK_MPS_lam=TK_MPS_cpc=NULL
      MH_MLE_alp=MH_MLE_lam=MH_MLE_cpc=NULL
      MH_MPS_alp=MH_MPS_lam=MH_MPS_cpc=NULL
      start_time=date()
    }  
    
    X=try(GenDataRcpp(c(0.1,100),para[1:2],R,k),silent=TRUE)
    if(is.character(X)) next
    if(X[m]==100) X[m]=X[m-1]+0.1
    
    #################################################################
    #------------- Bayes using TK and LK
    ############s#####################################################
    
    TK_Res_mle=suppressWarnings(try(TK(para[1:2],X,R,q,c,L,U,P0,k,"LF"),silent=TRUE))
    if(is.character(TK_Res_mle)) next
    if(is.nan(TK_Res_mle$SEL[3])|is.na(TK_Res_mle$SEL[3]))  next
    if(is.nan(TK_Res_mle$GEL[3])|is.na(TK_Res_mle$GEL[3]))  next
    if(is.nan(TK_Res_mle$LINEX[3])|is.na(TK_Res_mle$LINEX[3]))  next
    if((TK_Res_mle$SEL[1]<=0)|(TK_Res_mle$SEL[1]>3*para[1])) next
    if((TK_Res_mle$SEL[2]<=0)|(TK_Res_mle$SEL[2]>3*para[2])) next
    if((TK_Res_mle$SEL[3]<=0)|(TK_Res_mle$SEL[3]>3*para[3])) next
    if((TK_Res_mle$GEL[1]<=0)|(TK_Res_mle$GEL[1]>3*para[1])) next
    if((TK_Res_mle$GEL[2]<=0)|(TK_Res_mle$GEL[2]>3*para[2])) next
    if((TK_Res_mle$GEL[3]<=0)|(TK_Res_mle$GEL[3]>3*para[3])) next
    if((TK_Res_mle$LINEX[1]<=0)|(TK_Res_mle$LINEX[1]>3*para[1])) next
    if((TK_Res_mle$LINEX[2]<=0)|(TK_Res_mle$LINEX[2]>3*para[2])) next
    if((TK_Res_mle$LINEX[3]<=0)|(TK_Res_mle$LINEX[3]>3*para[3])) next
    
    TK_MLE_alp=rbind(TK_MLE_alp,c(TK_Res_mle$SEL[1],(TK_Res_mle$SEL[1]-para[1])^2,
                                  TK_Res_mle$GEL[1],(TK_Res_mle$GEL[1]-para[1])^2,
                                  TK_Res_mle$LINEX[1],(TK_Res_mle$LINEX[1]-para[1])^2))
    TK_MLE_lam=rbind(TK_MLE_lam,c(TK_Res_mle$SEL[2],(TK_Res_mle$SEL[2]-para[2])^2,
                                  TK_Res_mle$GEL[2],(TK_Res_mle$GEL[2]-para[2])^2,
                                  TK_Res_mle$LINEX[2],(TK_Res_mle$LINEX[2]-para[2])^2))
    TK_MLE_cpc=rbind(TK_MLE_cpc,c(TK_Res_mle$SEL[3],(TK_Res_mle$SEL[3]-para[3])^2,
                                  TK_Res_mle$GEL[3],(TK_Res_mle$GEL[3]-para[3])^2,
                                  TK_Res_mle$LINEX[3],(TK_Res_mle$LINEX[3]-para[3])^2))
    
    #################################################################
    #------------- Bayes using TK and PSF
    ############s#####################################################
    
    TK_Res_mps=suppressWarnings(try(TK(para[1:2],X,R,q,c,L,U,P0,k,"SPF"),silent=TRUE))
    if(is.character(TK_Res_mps)) next
    if(is.nan(TK_Res_mps$SEL[3])|is.na(TK_Res_mps$SEL[3]))  next
    if(is.nan(TK_Res_mps$GEL[3])|is.na(TK_Res_mps$GEL[3]))  next
    if(is.nan(TK_Res_mps$LINEX[3])|is.na(TK_Res_mps$LINEX[3]))  next
    
    if((TK_Res_mps$SEL[1]<=0)|(TK_Res_mps$SEL[1]>3*para[1])) next
    if((TK_Res_mps$SEL[2]<=0)|(TK_Res_mps$SEL[2]>3*para[2])) next
    if((TK_Res_mps$SEL[3]<=0)|(TK_Res_mps$SEL[3]>3*para[3])) next
    if((TK_Res_mps$GEL[1]<=0)|(TK_Res_mps$GEL[1]>3*para[1])) next
    if((TK_Res_mps$GEL[2]<=0)|(TK_Res_mps$GEL[2]>3*para[2])) next
    if((TK_Res_mps$GEL[3]<=0)|(TK_Res_mps$GEL[3]>3*para[3])) next
    if((TK_Res_mps$LINEX[1]<=0)|(TK_Res_mps$LINEX[1]>3*para[1])) next
    if((TK_Res_mps$LINEX[2]<=0)|(TK_Res_mps$LINEX[2]>3*para[2])) next
    if((TK_Res_mps$LINEX[3]<=0)|(TK_Res_mps$LINEX[3]>3*para[3])) next
    
    TK_MPS_alp=rbind(TK_MPS_alp,c(TK_Res_mps$SEL[1],MSE.fun(para[1],TK_Res_mps$SEL[1],"SEL",0),
                                  TK_Res_mps$GEL[1],MSE.fun(para[1],TK_Res_mps$GEL[1],"GEL",q),
                                  TK_Res_mps$LINEX[1],MSE.fun(para[1],TK_Res_mps$LINEX[1],"LINEX",c)))
    TK_MPS_lam=rbind(TK_MPS_lam,c(TK_Res_mps$SEL[2],MSE.fun(para[2],TK_Res_mps$SEL[2],"SEL",0),
                                  TK_Res_mps$GEL[2],MSE.fun(para[2],TK_Res_mps$GEL[2],"GEL",q),
                                  TK_Res_mps$LINEX[2],MSE.fun(para[2],TK_Res_mps$LINEX[2],"LINEX",c)))
    TK_MPS_cpc=rbind(TK_MPS_cpc,c(TK_Res_mps$SEL[3],MSE.fun(para[3],TK_Res_mps$SEL[3],"SEL",0),
                                  TK_Res_mps$GEL[3],MSE.fun(para[3],TK_Res_mps$GEL[3],"GEL",q),
                                  TK_Res_mps$LINEX[3],MSE.fun(para[3],TK_Res_mps$LINEX[3],"LINEX",c)))
    
    
    TK_Res=NULL
    TK_Res=rbind(TK_Res,c(colMeans(TK_MLE_alp),colMeans(TK_MPS_alp)))
    TK_Res=rbind(TK_Res,c(colMeans(TK_MLE_lam),colMeans(TK_MPS_lam)))
    TK_Res=rbind(TK_Res,c(colMeans(TK_MLE_cpc),colMeans(TK_MPS_cpc)))
    
    colnames(TK_Res)=c("MLE:SEL","MSE", "MLE:GEL","MSE","MLE:LINEX","MSE","MPS:SEL","MSE","MPS:GEL","MSE","MPS:LINEX","MSE")
    rownames(TK_Res)=c("alp","lam","cpc")
  
    ## Save the results in Matrices

    cat("====== Scheme ",sch,"=== Iteration no.", it," ===========================\n")
    
    cat("\n Generating MH samples based on likelihood function...\n");
    mle_res=suppressWarnings(try(Estim(para[1:2],X,R,L,U,P0,k,"LF"),silent=TRUE))
    # Check convergence
    if(is.character(mle_res)) next
    if((mle_res$Estim_para[1]<=0)|(mle_res$Estim_para[1]>=3*para[1])) next
    if((mle_res$Estim_para[2]<=0)|(mle_res$Estim_para[2]>=3*para[2])) next
    if((mle_res$Estim_Cpc<=0)|(mle_res$Estim_Cpc>=3*para[3])) next
    
    mle_est=c(as.numeric(mle_res$Estim_para),as.numeric(mle_res$Estim_Cpc))
    mle_ese=c(as.numeric(mle_res$SE_para),as.numeric(mle_res$SE_Cpc))
    
    Bayes_MH_MLE=MH_sample("LF", para[1:2],mle_ese,R,X,k,L,U,P0,MCMC_size,Burn_in,q,c,0,T)
    
    if(is.nan(Bayes_MH_MLE$est_SEL[1])) next
    if(is.nan(Bayes_MH_MLE$est_SEL[2])) next
    if(is.nan(Bayes_MH_MLE$est_SEL[3])) next
    if(is.nan(Bayes_MH_MLE$est_GEL[1])) next
    if(is.nan(Bayes_MH_MLE$est_GEL[2])) next
    if(is.nan(Bayes_MH_MLE$est_GEL[3])) next
    if(is.nan(Bayes_MH_MLE$est_LINEX[1])) next
    if(is.nan(Bayes_MH_MLE$est_LINEX[2])) next
    if(is.nan(Bayes_MH_MLE$est_LINEX[3])) next
    
    if(is.na(Bayes_MH_MLE$est_SEL[1])) next
    if(is.na(Bayes_MH_MLE$est_SEL[2])) next
    if(is.na(Bayes_MH_MLE$est_SEL[3])) next
    if(is.na(Bayes_MH_MLE$est_GEL[1])) next
    if(is.na(Bayes_MH_MLE$est_GEL[2])) next
    if(is.na(Bayes_MH_MLE$est_GEL[3])) next
    if(is.na(Bayes_MH_MLE$est_LINEX[1])) next
    if(is.na(Bayes_MH_MLE$est_LINEX[2])) next
    if(is.na(Bayes_MH_MLE$est_LINEX[3])) next
    
    MH_MLE_alp =rbind(MH_MLE_alp ,c(Bayes_MH_MLE$est_SEL[1],MSE.fun(para[1],Bayes_MH_MLE$est_SEL[1],"SEL",0),
                                    Bayes_MH_MLE$HPD_alp ,(Bayes_MH_MLE$HPD_alp[1]<para[1])*(para[1]<Bayes_MH_MLE$HPD_alp[2]),
                                    Bayes_MH_MLE$est_GEL[1],MSE.fun(para[1],Bayes_MH_MLE$est_GEL[1],"GEL",q),
                                    Bayes_MH_MLE$est_LINEX[1],MSE.fun(para[1],Bayes_MH_MLE$est_LINEX[1],"LINEX",c)))
    
    MH_MLE_lam =rbind(MH_MLE_lam ,c(Bayes_MH_MLE$est_SEL[2],MSE.fun(para[2],Bayes_MH_MLE$est_SEL[2],"SEL",0),
                                    Bayes_MH_MLE$HPD_lam ,(Bayes_MH_MLE$HPD_lam[1]<para[2])*(para[2]<Bayes_MH_MLE$HPD_lam[2]),
                                    Bayes_MH_MLE$est_GEL[2],MSE.fun(para[2],Bayes_MH_MLE$est_GEL[2],"GEL",q),
                                    Bayes_MH_MLE$est_LINEX[2],MSE.fun(para[2],Bayes_MH_MLE$est_LINEX[2],"LINEX",c)))
    
    MH_MLE_cpc =rbind(MH_MLE_cpc ,c(Bayes_MH_MLE$est_SEL[3],MSE.fun(para[3],Bayes_MH_MLE$est_SEL[3],"SEL",0),
                                    Bayes_MH_MLE$HPD_cpc ,(Bayes_MH_MLE$HPD_cpc[1]<para[3])*(para[3]<Bayes_MH_MLE$HPD_cpc[2]),
                                    Bayes_MH_MLE$est_GEL[3],MSE.fun(para[3],Bayes_MH_MLE$est_GEL[3],"GEL",q),
                                    Bayes_MH_MLE$est_LINEX[3],MSE.fun(para[3],Bayes_MH_MLE$est_LINEX[3],"LINEX",c)))
    
    
    cat("\n Generating MH samples based on product of spacing function...\n");
    Bayes_MH_MPS=MH_sample("SPF", para[1:2],mle_ese,R,X,k,L,U,P0,MCMC_size,Burn_in,q,c,0,T)
    
    if(is.nan(Bayes_MH_MPS$est_SEL[1])) next
    if(is.nan(Bayes_MH_MPS$est_SEL[2])) next
    if(is.nan(Bayes_MH_MPS$est_SEL[3])) next
    if(is.nan(Bayes_MH_MPS$est_GEL[1])) next
    if(is.nan(Bayes_MH_MPS$est_GEL[2])) next
    if(is.nan(Bayes_MH_MPS$est_GEL[3])) next
    if(is.nan(Bayes_MH_MPS$est_LINEX[1])) next
    if(is.nan(Bayes_MH_MPS$est_LINEX[2])) next
    if(is.nan(Bayes_MH_MPS$est_LINEX[3])) next
    
    if(is.na(Bayes_MH_MPS$est_SEL[1])) next
    if(is.na(Bayes_MH_MPS$est_SEL[2])) next
    if(is.na(Bayes_MH_MPS$est_SEL[3])) next
    if(is.na(Bayes_MH_MPS$est_GEL[1])) next
    if(is.na(Bayes_MH_MPS$est_GEL[2])) next
    if(is.na(Bayes_MH_MPS$est_GEL[3])) next
    if(is.na(Bayes_MH_MPS$est_LINEX[1])) next
    if(is.na(Bayes_MH_MPS$est_LINEX[2])) next
    if(is.na(Bayes_MH_MPS$est_LINEX[3])) next
    
    MH_MPS_alp =rbind(MH_MPS_alp ,c(Bayes_MH_MPS$est_SEL[1],MSE.fun(para[1],Bayes_MH_MPS$est_SEL[1],"SEL",0),
                                    Bayes_MH_MPS$HPD_alp ,(Bayes_MH_MPS$HPD_alp[1]<para[1])*(para[1]<Bayes_MH_MPS$HPD_alp[2]),
                                    Bayes_MH_MPS$est_GEL[1],MSE.fun(para[1],Bayes_MH_MPS$est_GEL[1],"GEL",q)^2,
                                    Bayes_MH_MPS$est_LINEX[1],MSE.fun(para[1],Bayes_MH_MPS$est_LINEX[1],"LINEX",c)))
    
    MH_MPS_lam =rbind(MH_MPS_lam ,c(Bayes_MH_MPS$est_SEL[2],MSE.fun(para[2],Bayes_MH_MPS$est_SEL[2],"SEL",0),
                                    Bayes_MH_MPS$HPD_lam ,(Bayes_MH_MPS$HPD_lam[1]<para[2])*(para[2]<Bayes_MH_MPS$HPD_lam[2]),
                                    Bayes_MH_MPS$est_GEL[2],MSE.fun(para[2],Bayes_MH_MPS$est_GEL[2],"GEL",q),
                                    Bayes_MH_MPS$est_LINEX[2],MSE.fun(para[2],Bayes_MH_MPS$est_LINEX[2],"LINEX",c)))
    
    MH_MPS_cpc =rbind(MH_MPS_cpc ,c(Bayes_MH_MPS$est_SEL[3],MSE.fun(para[3],Bayes_MH_MPS$est_SEL[3],"SEL",0),
                                    Bayes_MH_MPS$HPD_cpc ,(Bayes_MH_MPS$HPD_cpc[1]<para[3])*(para[3]<Bayes_MH_MPS$HPD_cpc[2]),
                                    Bayes_MH_MPS$est_GEL[3],MSE.fun(para[3],Bayes_MH_MPS$est_GEL[3],"GEL",q),
                                    Bayes_MH_MPS$est_LINEX[3],MSE.fun(para[3],Bayes_MH_MPS$est_LINEX[3],"LINEX",c)))
    
    print(TK_Res)
    if(nrow(MH_MLE_alp)>1){
      ss=NULL
      ss=rbind(ss,c(colMeans(MH_MLE_alp[,1:4],na.rm=T),mean(MH_MLE_alp[,6],na.rm=T),colMeans(MH_MPS_alp[,1:4],na.rm=T),mean(MH_MPS_alp[,6],na.rm=T)))
      ss=rbind(ss,c(colMeans(MH_MLE_alp[,7:8],na.rm=T),NA,NA,NA,NA,colMeans(MH_MPS_alp[,7:8],na.rm=T),NA,NA,NA,NA))
      ss=rbind(ss,c(colMeans(MH_MLE_alp[,9:10],na.rm=T),NA,NA,NA,NA,colMeans(MH_MPS_alp[,9:10],na.rm=T),NA,NA,NA,NA))
      ss=rbind(ss,c(colMeans(MH_MLE_lam[,1:4],na.rm=T),mean(MH_MLE_lam[,6],na.rm=T),colMeans(MH_MPS_lam[,1:4],na.rm=T),mean(MH_MPS_lam[,6],na.rm=T)))
      ss=rbind(ss,c(colMeans(MH_MLE_lam[,7:8],na.rm=T),NA,NA,NA,NA,colMeans(MH_MPS_lam[,7:8],na.rm=T),NA,NA,NA,NA))
      ss=rbind(ss,c(colMeans(MH_MLE_lam[,9:10],na.rm=T),NA,NA,NA,NA,colMeans(MH_MPS_lam[,9:10],na.rm=T),NA,NA,NA,NA))
      ss=rbind(ss,c(colMeans(MH_MLE_cpc[,1:4],na.rm=T),mean(MH_MLE_cpc[,6],na.rm=T),colMeans(MH_MPS_cpc[,1:4],na.rm=T),mean(MH_MPS_cpc[,6],na.rm=T)))
      ss=rbind(ss,c(colMeans(MH_MLE_cpc[,7:8],na.rm=T),NA,NA,NA,NA,colMeans(MH_MPS_cpc[,7:8],na.rm=T),NA,NA,NA,NA))
      ss=rbind(ss,c(colMeans(MH_MLE_cpc[,9:10],na.rm=T),NA,NA,NA,NA,colMeans(MH_MPS_cpc[,9:10],na.rm=T),NA,NA,NA,NA))
      colnames(ss)=c("M.Est","M.MSE","M.LCI","M.UCI","M.CP", "P.Est","P.MSE","P.LCI","P.UCI","P.CP")
      rownames(ss)=c("alp.sel","alp.gel","alp.lines","lam.sel","lam.gel","lam.linex","cpc.sel","cpc.gel","cpc.linex")
      print(ss)    
      
    }
    it=it+1
  } 
  
  return(list(TK_MLE_alp=TK_MLE_alp,TK_MLE_lam=TK_MLE_lam,TK_MLE_cpc=TK_MLE_cpc,
              TK_MPS_alp=TK_MPS_alp,TK_MPS_lam=TK_MPS_lam,TK_MPS_cpc=TK_MPS_cpc,
              MH_MLE_alp=MH_MLE_alp,MH_MLE_lam=MH_MLE_lam,MH_MLE_cpc=MH_MLE_cpc,
              MH_MPS_alp=MH_MPS_alp,MH_MPS_lam=MH_MPS_lam,MH_MPS_cpc=MH_MPS_cpc))
}

