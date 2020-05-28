### PCG cv & specify df

Est_qiu_theo_lambda=function(df, alpha=0.05, c0=0.25){ 
  require(glmnet)
  n = dim(df)[1]; p = dim(df)[2]
  t0=2
  IndMatrix = matrix(1, p, p) - diag(rep(1, p))
  Eresidual = matrix(0, n, p) # regression residuals matrix n*p
  CoefMatrix = matrix(0, p, p - 1) # regression coefficient matrix p*p-1
  meanX = colMeans(df)
  X = t(t(df) - meanX)
  XS = matrix(0, n, p)
  # XS: Standardized X
  shat=sqrt(n/(log(p)^3))
  lambda=sqrt(2*(2+0.01)*log(p/shat)/n)
    
  for (i in 1 : p){
    XS[, i] = X[, i] / sd(X[, i])
  }
  colnames(X)=colnames(df)
  colnames(XS)=colnames(df)
  
  for (i in 1 : p){
    out=glmnet(XS[, -i], X[, i], family = "gaussian", lambda = lambda)

    Coef=coef(out, s=lambda)[-1]
    CoefMatrix[i,] = Coef / apply(X[, -i], 2, sd)
    Eresidual[,i] = X[, i] - Predict 
  }
  
  # Get Est
  {
    CovRes = t(Eresidual) %*% Eresidual / n # residuals covariance
    Est = matrix(1, p, p) # estimated partial correlation (rho hat in the paper )
    
    for (i in 1 : (p - 1)){
      for (j in (i + 1) : p){
        temp = Eresidual[, i] * Eresidual[, j] + Eresidual[, i]^2 * CoefMatrix[j, i] + Eresidual[, j]^2 * CoefMatrix[i, j - 1]
        Est[i, j] = mean(temp) / sqrt(diag(CovRes)[i] * diag(CovRes)[j])
        Est[j, i] = Est[i, j]
      }
    }
  } # end Est
  
  # inference
  {
    #EstThresh = Est * ( abs(Est) >= (t0 * sqrt(log(p) / n) * IndMatrix) ) 
    EstThresh = Est * ( abs(Est) >= (sqrt(log(p) / n) * IndMatrix) )
    EstThresh[abs(EstThresh)>1]=1
    kappa = (n / 3) * mean( colSums(Eresidual^4) / (colSums(Eresidual^2))^2 )  # forth moment, a number 
    
    tau = seq(0, 3.5, 0.01); smax = n / 2; lentau = length(tau) 
    
    resprop = list() # selected edges with different tau's, a list of 351 elements
    rejectprop = c()
    rejectpropC0 = c() 
    for (i in 1 : lentau){ # tau vary from 0 to 3.50 by 0.01, length=351
      Threshold = tau[i] * sqrt(kappa * log(p) / n) * (1 - EstThresh^2)
      
      # c=0
      SRec = 1 * (abs(Est) > Threshold) # selected edge (matrix with 0 & 1) at tau[i]
      NoNSRec = 1 * (SRec == 0)
      resprop[[i]] = which(SRec == 1, arr.ind = TRUE) # select those significant edges at tau[i], off-diagonal elements, first columns, then second columns
      rejectprop = c(rejectprop, max(1, (sum(SRec) - p))) 
      
      # c-level
      SRecC0 = 1 * (abs(Est) - c0 > Threshold)
      rejectpropC0 = c(rejectpropC0, max(1, (sum(SRecC0) - p))) # number of selected edges? corresponding to tau
    } # end loop
    
    # c=0
    FDPprop = 2 * (p * (p - 1)) * ( 1 - pnorm( tau * sqrt(log(p)) ) ) / rejectprop # FDP corresponding to each tau (page 10)
    
    FDPresprop = c()
    
    # determine thresholding parameter tau by controling FDP 
    if (sum(FDPprop <= alpha) > 0) tauprop = min(c(2, tau[FDPprop <= alpha])) 
    if (sum(FDPprop <= alpha) == 0) tauprop = 2
    Threshold = tauprop * sqrt(kappa * log(p) / n) * (1 - EstThresh^2)
    SRec = 1 * (abs(Est) > Threshold); NoNSRec = 1 * (SRec == 0) # SRec is a matrix (0-1 matrix)
    FDPresprop = which(SRec == 1, arr.ind = TRUE) # selected edge location
    
    sigs=FDPresprop[which(FDPresprop[,1]!=FDPresprop[,2]),]%>% as.data.frame()
    colnames(sigs)=c("node1","node2")
    
    # c=0.15
    FDPpropC0 = 2 * ( sqrt(n) * p * ( 1 - pnorm( tau * sqrt(log(p)) ) ) + p * (p - 1 - sqrt(n)) * ( 1 - pnorm( sqrt(n) * c0 / sqrt(kappa) + tau * sqrt(log(p)) ) ) ) / rejectpropC0
    
    if (sum(FDPpropC0 <= alpha) > 0) taupropC0 = min(c(2, tau[FDPpropC0 <= alpha]))
    else taupropC0 = 2
    ThresholdC0 = taupropC0 * sqrt(kappa * log(p) / n) * (1 - EstThresh^2)
    SRecC0 = 1 * (abs(Est) - c0 > ThresholdC0); NoNSRecC0 = 1 * (SRecC0 == 0)
    FDPrespropC0 = which(SRecC0 == 1, arr.ind = TRUE) # selected edge location
    
    sigs0=FDPrespropC0[which(FDPrespropC0[,1]!=FDPrespropC0[,2]),]%>% as.data.frame()
    colnames(sigs0)=c("node1","node2")
  } # end Inference
  return(list(Est=Est, Sigs=sigs, Sigs0=sigs0))
  
}# end function 