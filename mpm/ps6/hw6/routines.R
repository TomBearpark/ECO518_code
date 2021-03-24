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
  out    = lmBoot$coefficients

  
    return(out)    
}

# Clustered bootstrap
cluster_boot = function(form, data, clusterVar){
  
  nClusters = length(unlist( unique(data[,clusterVar])))
    

  
  dataBoot = data %>% group_by(stateid) %>% group_nest() %>% 
    sample_n(nClusters, replace = TRUE) %>% unnest(data) %>% ungroup()
  
  
  lmBoot = lm(form, dataBoot)
  SE     = coeftest(lmBoot, vcovCL(lmBoot, cluster= data[,clusterVar] ))[, "Std. Error"]
  out    = list(lmBoot$coefficients, SE)
  names(out) = c('coeff', 'SE')
  
  return(out)    
  
}