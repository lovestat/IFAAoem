
#-------------------------------------------------------------------
## function for cross validation using OEM package
#-------------------------------------------------------------------

runOEM=function(
  x,
  y,
  nPredics,
  nfolds=10,
  nLam=100,
  standardize,
  intercept=TRUE,
  zeroSDCut=0,
  seed,
  seedi
){
  
  results=list()

  ## Shang: calculate ncol for big.matrix
  nBeta=ncol(x)
  nObsAll=length(y)

  # remove near constant x columns
  ##: Shang: calculate column standard deviations in big.matrix
  
  sdX=biganalytics::colsd(x)
  
  xWithNearZeroSd=which(sdX<=zeroSDCut) 
  if(length(xWithNearZeroSd)>0){
    x=x[,-xWithNearZeroSd]
  }
  rm(sdX)
  
  if(standardize)x=x/(apply(x,2,sd))

  x=x
  y=as.vector(y)
  
  # calculate lambda max
  
  ## Shang: x*y in big.matrix
  
  xtimesy <- biganalytics::apply(x, 2, function(x) x * as.vector(y))
  
  lamMax=max(abs(colSums(xtimesy)))/nObsAll
  lamVec=seq(lamMax,0,length=(nLam+1))[1:nLam]

  cvResul=cvOEM(x=x,y=y,lamList=lamVec,nfolds=nfolds,zeroSDCut=zeroSDCut,
                intercept=intercept,seed=seed,seedi=seedi)
    
  lamOpi=as.numeric(cvResul$optiLam)
  rm(cvResul)

  ## Shang: big.oem
  finalMCPrun=big.oem(x=x,y=y,penalty="mcp",gamma = 3,intercept=intercept,
      standardize=FALSE,lambda=lamOpi)
  
  rm(x,y)
  
  OverallIntercp=finalMCPrun$beta$mcp[1,]
  finalLassoRunBeta=as.vector(finalMCPrun$beta$mcp[-1,])
  rm(finalMCPrun)
  
  # convert back to the full beta if there near constant x columns
  if(length(xWithNearZeroSd)>0){
    betaTrans=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                  unSelectList=sort(xWithNearZeroSd),
                                  newBetaNoInt=finalLassoRunBeta)
    beta=betaTrans$finalBeta
    rm(betaTrans)
  } else {
    beta=finalLassoRunBeta
  }

  rm(xWithNearZeroSd)
  
  # return 
  results$betaNoInt=beta[-seq(1,length(beta),by=(nPredics+1))]
  results$betaInt=beta
  rm(beta)
  results$overalIntercp=OverallIntercp
  
  return(results)
}

# runOEM(x=x,y=y,nPredics)