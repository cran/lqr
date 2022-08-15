Log.best.lqr = function(formula,data = NULL,
                        subset = NULL,
                        p=0.5,a=0,b=1,epsilon = 0.001,precision = 10^-6,criterion = "AIC")
{
  cat('\n')
  call <- match.call()
  cat("Call:\n")
  print(call)
  cat('\n')
  
  if(is.null(subset)){subset = ""}
  xy = get_xy(formula0 = formula,
              data0 = data,
              subset0 = subset)
  x = xy$x
  y.original = xy$y
  
  y = log((y.original - a + epsilon)/(b - y.original + epsilon))
  
  if(length(p)>1) stop("The function Log.best.lqr is only available for one quantile.")
  ## Verify error at parameters specification
  envelope = TRUE
  
  #No data
  if( (length(x) == 0) | (length(y) == 0)) stop("All parameters must be provided.")
  
  #Validating if exists NA's
  if(sum(y[is.na(y)==TRUE]) > 0) stop("There are some NA's values in y")
  if(sum(y[is.na(x)==TRUE]) > 0) stop("There are some NA's values in x")
  
  #Validating dims data set
  if(ncol(as.matrix(y)) > 1) stop("y must have just one column")
  if( length(y) != nrow(as.matrix(x)) ) stop("x variable does not have the same number of lines than y")
  
  #Validating supports
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1)")
  if(precision <= 0) stop("precision must be a positive value (suggested to be small)")
  if(criterion != "" && criterion != "AIC" && criterion != "BIC" && criterion != "HQ" && criterion != "loglik") stop("The criterion must be AIC, BIC, HQ or loglik.")
  
  #Running the algorithm
  #out <- suppressWarnings(EM(y,x,p,dist,nu,gamma,precision,envelope))
  
  vdist = c("normal","t","laplace","slash","cont")
  vDIST = c("Normal","Student-t","Laplace","Slash","Cont. Normal")
  obj.out  = vector("list", 5)
  
  #######################
  #oma margenes totales
  
  par(mfrow=c(2,3),oma=c(0,1,1,0),mar=c(5, 4, 4, 2.5))
  for(k in 1:5)
  {
    obj.out[[k]] <- EM2(y = y,x = x,p = p,dist = vdist[k],precision = precision,envelope=TRUE)
  }
  plot.new()
  par(mfrow=c(1,1))
  #mtext("Histogram of residuals and fitted densities", side = 3, line = 1, outer = TRUE,cex=1.3)
  
  
  #HISTOGRAM AND FITTED DENSITIES
  
  resall = unlist(lapply(X = obj.out,FUN = residuals))
  par(mfrow=c(2,3),oma=c(0,1,0,1),mar=c(5, 4, 5, 0.5))
  seqq = seq(from=min(resall),to = max(resall),length.out = 1000)
  up = dSKD(y = 0,mu = 0,sigma = obj.out[[3]]$theta[dim(x)[2]+1],p = p,dist = vdist[3])
  
  folginha = 0.1
  for(k in 1:5)
  {
    hist(x = obj.out[[k]]$residuals,freq=FALSE,breaks=sqrt(length(y)),xlab = "residuals",
         main = vDIST[k],cex.main=1.5,ylim=c(0,(1+folginha)*up))
    dens = dSKD(y = seqq,mu = 0,sigma = obj.out[[k]]$theta[dim(x)[2]+1],p = p,dist = vdist[k],
                nu=obj.out[[k]]$nu,gamma = obj.out[[k]]$gamma)
    lines(seqq,dens,lwd=1.5,col="blue")
  }
  plot.new()
  par(mfrow=c(1,1))
  #mtext("Histogram of residuals and fitted densities", side = 3, line = 1, outer = TRUE,cex=1.3)
  
  RES = matrix(data = NA,nrow = 4,ncol = 5)
  for(k in 1:5)
  {
    RES[1,k] = (obj.out[[k]])$AIC
    RES[2,k] = (obj.out[[k]])$BIC
    RES[3,k] = (obj.out[[k]])$HQ
    RES[4,k] = (obj.out[[k]])$loglik
  }
  colnames(RES) = c(vDIST[1:4],"C. Normal")
  rownames(RES) = c("AIC","BIC","HQ","loglik")
  
  if(criterion == "AIC")
  {
    index = which(RES[1,] == min(RES[1,]))
  }
  if(criterion == "BIC")
  {
    index = which(RES[2,] == min(RES[2,])) 
  }
  if(criterion == "HQ")
  {
    index = which(RES[3,] == min(RES[3,])) 
  }
  if(criterion == "loglik")
  {
    index = which(RES[4,] == max(RES[4,])) 
  }
  
  out  = obj.out[[index]]
  dist = vdist[index] 
  cat('\n')
  cat("Criterion:",criterion,'\n')
  cat("Best fit:",vDIST[index],'\n')
  cat("Quantile:",p,'\n')
  cat('\n')
  cat('Model Likelihood-Based criteria:\n')
  cat('\n')
  print(RES)
  #cat("Iterations =",out$iter)
  cat('\n')
  cat('Estimates:\n')
  cat('\n')
  print(out$table)
  cat('---\n')
  cat('Signif. codes:  0 "***" 0.001 "**" 0.01 "*" 0.05 "." 0.1 " " 1\n')
  cat('\n')
  cat('sigma =',out$theta[ncol(as.matrix(x))+1],'\n')
  if(dist == "normal" || dist == "laplace"){
    cat('\n')
  }
  if(dist == "t" || dist == "slash"){
    cat('nu    =',out$nu,'\n')
    cat('\n')
  }
  if(dist == "cont"){
    cat('nu    =',out$nu,'\n')
    cat('gamma =',out$gamma,'\n')
    cat('\n')
  }
  
  if(dist == "normal" || dist == "laplace"){
    obj.out = list(iter = out$iter,criteria = out$criterio,
                   beta = out$theta[1:ncol(as.matrix(x))],
                   sigma= out$theta[ncol(as.matrix(x))+1],
                   SE=out$SE,table = out$table,loglik=out$loglik,
                   AIC=out$AIC,BIC=out$BIC,HQ=out$HQ,
                   fitted.values = out$fitted.values,residuals=out$residuals)
  }
  
  if(dist == "t" || dist == "slash"){
    obj.out = list(iter = out$iter,criteria = out$criterio,
                   beta = out$theta[1:ncol(as.matrix(x))],
                   sigma= out$theta[ncol(as.matrix(x))+1],nu = out$nu,
                   SE=out$SE,table = out$table,loglik=out$loglik,
                   AIC=out$AIC,BIC=out$BIC,HQ=out$HQ,
                   fitted.values = out$fitted.values,residuals=out$residuals)
  }
  
  if(dist == "cont"){
    obj.out = list(iter = out$iter,criteria = out$criterio,
                   beta = out$theta[1:ncol(as.matrix(x))],
                   sigma= out$theta[ncol(as.matrix(x))+1],nu = out$nu,
                   gamma = out$gamma,
                   SE=out$SE,table = out$table,loglik=out$loglik,
                   AIC=out$AIC,BIC=out$BIC,HQ=out$HQ,
                   fitted.values = out$fitted.values,residuals=out$residuals)
  }
  class(obj.out)  =  "qr"
  invisible(obj.out)
  
  
  
  
  
  
  
  
  
  
  
  pred = function(predlog,a,b)
  {
    return((b*exp(predlog)+a)/(1+exp(predlog)))
  }
  
  #message("The interpretation of the regression coefficients is analogous to the interpretation of the coefficients of a logistic regression for binary outcomes. For references, please check Galarza, C.M., Zhang P. and Lachos, V.H. (2020). Logistic Quantile Regression for Bounded Outcomes Using a Family of Heavy-Tailed Distributions. Sankhya B.")
  
  obj.out$fitted.values = pred(predlog = obj.out$fitted.values,a,b)
  obj.out$residuals = y.original - obj.out$fitted.values
  invisible(obj.out)
}
