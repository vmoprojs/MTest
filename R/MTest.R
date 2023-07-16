####################################################
### File name: MTest.r
####################################################
MTest <- function(object, nboot = 100,
                  nsam = NULL,trace = FALSE,seed = NULL,
                  valor_vif = 0.9)
{
  datos <- object$model
  ff <- formula(object)
  if(is.null(nsam)){nsam = nrow(datos)}
  
  vals <- 1:nrow(datos)
  
  if(!is.null(seed)) {set.seed(seed)}
  
  sol.rsq <- NULL
  sol.vif <- NULL
  i = 1
  
  tt <- proc.time()
  while(i <=nboot)
  {
    sam <- sample(vals,nsam,replace = TRUE)
    aux <- datos[sam,]
    maux <- lm(ff,data = aux)
    sm  <- summary(maux)
    if(any(attr(terms(object),"order")>1))
    {
      vif.vals <- suppressMessages(car::vif(maux,type = "predictor"))
      vif.vals <- vif.vals[,3]
    }else{
      vif.vals <- car::vif(maux,type = "terms")
    }
    
    Raux <- (vif.vals-1)/vif.vals
    
    s1 <- c(sm$r.squared,Raux)
    sol.rsq <- rbind(sol.rsq,s1)
    
    sol.vif <- rbind(sol.vif,vif.vals)
    
    if(trace)
    {
      cat("Iteration",i,"out of ",nboot,"\n")
    }
    i = i+1
  }
  tt <- proc.time()-tt
  print(tt)
  
  
  
  
  
  pval_vif <- NULL
  
  for(j in 2:ncol(sol.rsq))
  {
    pval_vif <- c(pval_vif,sum(sol.rsq[,j]>valor_vif)/nboot)
  }
  names(pval_vif) <- colnames(sol.rsq)[2:ncol(sol.rsq)]
  pval_klein <- NULL
  for(z in 2:ncol(sol.rsq))
  {
    pval_klein <- c(pval_klein,sum(sol.rsq[,1]<sol.rsq[,z])/nboot)
  }
  names(pval_klein) <- colnames(sol.rsq)[2:ncol(sol.rsq)]
  
  colnames(sol.rsq) <- c("global",paste(names(datos)[-1],sep =""))
  rownames(sol.rsq) <- 1:nrow(sol.rsq)
  return(list(Bvals= sol.rsq,pval_vif = pval_vif,pval_klein=pval_klein))
}