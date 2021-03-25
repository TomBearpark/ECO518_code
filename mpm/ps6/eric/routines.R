# Non-parametric bootstrap
np_boot = function(form, data){
  n = nrow(data)
  dataBoot = data[sample(n, replace = TRUE),]
  lmBoot   = lm(form, dataBoot)
  SE      = coeftest(lmBoot, vcovHC(lmBoot, type = "HC1"))[,"Std. Error"]
  
  # Organize output
  out        = list(lmBoot$coefficients, SE)
  names(out) = c('coeff', 'SE')
  
  return(out)    
}

# Parametric bootstrap
para_boot = function(form, data, model){
  dataBoot = data
  depVar   = strsplit(form, '~')[[1]][1]
  
  dataBoot[,depVar] = predict(model) + sample(model$residuals, replace = TRUE)
  
  lmBoot = lm(form, dataBoot)
  SE      = coeftest(lmBoot, vcovHC(lmBoot, type = "HC1"))[,"Std. Error"]  
  
  out        = list(lmBoot$coefficients, SE)
  names(out) = c('coeff', 'SE')
  
  return(out)    
  
}

# Clustered bootstrap
cluster_boot = function(form, data, clusterVar){
  
  nClusters = length(unlist( unique(data[,clusterVar])))
    

  
  dataBoot = data %>% group_by(stateid) %>% group_nest() %>% 
    sample_n(nClusters, replace = TRUE) %>% unnest(data) %>% ungroup()
  
  
  lmBoot = lm(form, dataBoot)
  SE     = coeftest(lmBoot, vcovCL(lmBoot, cluster= dataBoot[,clusterVar] ))[, "Std. Error"]
  out    = list(lmBoot$coefficients, SE)
  names(out) = c('coeff', 'SE')
  
  return(out)    
  
}

# Simulation of question 3
runSim = function(j, nBoot, N, alpha0, beta0, df){
  
  if(j%%5 == 0){
  print(paste0('Iteration ', j, '...'))    
  }

  
  qStar = .95
  
  # Draw data ----
  set.seed(j)
  X   = rnorm(N, mean = 0, sd = sqrt(4))
  eps = rt(N, df) 
  
  
  # Part i, make dataset, run regression, run bootstrap ----
  u        = eps
  Y1       = alpha0 + beta0 + u
  data1    = as.data.frame(Y1)
  data1$X  = X
  form1    = 'Y1~X'    
  lm1      = lm(form1, data1)
  betaHat1 = lm1$coefficients['X']
    
  # Get boot
  Res_np1   = replicate(nBoot, np_boot(form1, data1))
  Res_para1 = replicate(nBoot, para_boot(form1, data1, lm1))
  
  
  # Part ii ----
  u        = sigma0(X) * eps
  Y2       = alpha0 + beta0 + u
  data2    = as.data.frame(Y2)
  data2$X  = X
  form2    = 'Y2~X' 
  lm2      = lm(form2, data2)
  betaHat2 = lm2$coefficients['X']
  
  # Get boot
  Res_np2   = replicate(nBoot, np_boot(form2, data2))
  Res_para2 = replicate(nBoot, para_boot(form2, data2, lm2))
   
  # Add code that stores (a)-(e) ----
  
  # Get coefficients
  beta_np1    = sapply(Res_np1['coeff',]  , '[[', 'X')
  beta_para1  = sapply(Res_para1['coeff',], '[[', 'X')  
  beta_np2    = sapply(Res_np2['coeff',]  , '[[', 'X')
  beta_para2  = sapply(Res_para2['coeff',], '[[', 'X')  
  
  # Get standard errors
  se_np1    = sapply(Res_np1['SE',]  , '[[', 'X')
  se_para1  = sapply(Res_para1['SE',], '[[', 'X')  
  se_np2    = sapply(Res_np2['SE',]  , '[[', 'X')
  se_para2  = sapply(Res_para2['SE',], '[[', 'X')  

  # t-statistic
  tStar1 = abs(coeftest(lm1, vcovHC(lm1, type = "HC1"))['X', "t value"])
  tStar2 = abs(coeftest(lm2, vcovHC(lm2, type = "HC1"))['X', "t value"])
        
  # Asymptotic normal critical value (a)
  t_asymp = abs(qnorm((1-qStar)/2))
    
  # Nonparametric bootstrap critical value (b)
  t_np1   = quantile(  abs((beta_np1    - betaHat1)/se_np1), qStar)
  t_np2   = quantile(  abs((beta_np2    - betaHat2)/se_np2), qStar)
  
  # Residual bootstrap critical values (c)
  t_para1 = quantile(abs((beta_para1  - betaHat1)/se_para1), qStar)
  t_para2 = quantile(abs((beta_para2  - betaHat2)/se_para2), qStar)

  # t-tilde
  tTilde1 = abs(sqrt(N) * betaHat1)
  tTilde2 = abs(sqrt(N) * betaHat2)
    
  # Non-parametric BS t-tilde (d)
  tTilde_np1   = quantile(abs(sqrt(N) * (beta_np1 - betaHat1)), qStar)
  tTilde_np2   = quantile(abs(sqrt(N) * (beta_np2 - betaHat2)), qStar)
  
  # Non-parametric BS t-tilde (e)
  tTilde_para1   = quantile(abs(sqrt(N) * (beta_para1 - betaHat1)), qStar)
  tTilde_para2   = quantile(abs(sqrt(N) * (beta_para2 - betaHat2)), qStar)
  
  
  # Output test results ---
  ex1 = as.numeric(c(tStar1 > t_asymp, 
                     tStar1 > t_np1, 
                     tStar1 > t_para1, 
                     tTilde1 > tTilde_np1, 
                     tTilde1 > tTilde_para1))

  names(ex1) = c('a', 'b', 'c', 'd', 'e')  
  
  ex2 = as.numeric(c(tStar2 > t_asymp, 
                     tStar2 > t_np2, 
                     tStar2 > t_para2, 
                     tTilde2 > tTilde_np2, 
                     tTilde2 > tTilde_para2))
  names(ex2) = c('a', 'b', 'c', 'd', 'e')    
  out = list(ex1, ex2)
  names(out) = c('ex1', 'ex2')  
  return(out)
  
}


sigma0 = function(X){
  0.5*X
}
