X_sampling <- function(X_pool, snp_region,N_sample){
  ref_id_num <- nrow(X_pool)
  id = sample.int(ref_id_num, replace=TRUE, N_sample)
  X_sim <- X_pool[id, snp_region]
  return(X_sim)
}


sim_center_binom <- function(sim_data){
  sim_data$y.train = sim_data$y.train - mean(sim_data$y.train)
  sim_data$y1.train = sim_data$y1.train - mean(sim_data$y1.train)
  sim_data$y2.train = sim_data$y2.train - mean(sim_data$y2.train)
  
  sim_data$y.test = sim_data$y.test - mean(sim_data$y.test)
  sim_data$y1.test = sim_data$y1.test - mean(sim_data$y1.test)
  sim_data$y2.test = sim_data$y2.test - mean(sim_data$y2.test)
  
  return(sim_data)
}


DataGen_newhap_binom <- function(mc=2, n.train, n.test, p, b0=0, 
                                 nonzero_beta1, nonzero_beta2, gammas, 
                                 var1, var2, seed=NULL,group = "heter xy", X_pool=X_pool, snp_region=snp_region){
  
  gammas=matrix(gammas, ncol=2, byrow=F)
  
  if (is.null(seed)==F) {set.seed(seed)}
  pi.train=rbeta(n.train,2,3) # composition 40% of y1
  pi.test=rbeta(n.test,2,3)
  
  data=foreach (i = 1:mc)  %dorng% {
    
    # parameters
    beta1=beta2=rep(0, 1+p)
    beta1[1]=b0
    idx1=sample(2:(p+1), 2)
    idx2=2*(rbinom(2, 1, 0.5)-0.5)
    #idx2=rep(1, 2)
    beta1[idx1[1]]= idx2[1]* nonzero_beta1
    if (group == "heter xy"){
      beta2[idx1[2]]= idx2[2]* nonzero_beta2
    }
    if (group == "homo xy"){
      beta2[idx1[1]]= idx2[1]* nonzero_beta2 # homo
    }
    
    
    x.train = data.matrix(X_sampling(X_pool, snp_region, n.train))
    
    design=cbind(1, x.train)
    y1.train=design%*%beta1+rnorm(n.train, mean=0, sd=sqrt(var1))
    y2.train=design%*%beta2+rnorm(n.train, mean=0, sd=sqrt(var2))
    y.train=pi.train*y1.train + (1-pi.train)*y2.train
    
    Disease.train=NULL
    for (j in 1:nrow(gammas)) {
      a=-mean(gammas[j,1]*y1.train+gammas[j,2]*y2.train)
      pD=1/(1+ 1/exp(a+ gammas[j,1]*y1.train +gammas[j,2]* y2.train))
      D=rbinom(n.train, 1, pD)
      Disease.train=cbind(Disease.train, D)
    }
    
    # test
    x.test = data.matrix(X_sampling(X_pool, snp_region, n.test))
    design=cbind(1, x.test)
    y1.test=design%*%beta1+rnorm(n.test, mean=0, sd=sqrt(var1))
    y2.test=design%*%beta2+rnorm(n.test, mean=0, sd=sqrt(var2))
    y.test=pi.test*y1.test + (1-pi.test)*y2.test
    
    Disease.test=NULL
    for (j in 1:nrow(gammas)) {
      a=-mean(gammas[j,1]*y1.test+gammas[j,2]*y2.test)
      pD=1/(1+ 1/exp(a+ gammas[j,1]*y1.test +gammas[j,2]* y2.test ))
      D=rbinom(n.test, 1, pD)
      Disease.test=cbind(Disease.test, D)
    }
    
    list(y.train=y.train, y1.train =y1.train, y2.train=y2.train,
         x.train=x.train, pi.train=pi.train, D.train=Disease.train,
         y.test=y.test, y1.test =y1.test, y2.test=y2.test,
         x.test=x.test, pi.test=pi.test, D.test=Disease.test)
  }
  return(data)
}

