#---------------------------------------------------------------------------------------
## bootstrap results function for Lasso OLS HDCI
#---------------------------------------------------------------------------------------

bootResuHDCI=function(
  data,
  refTaxa,
  maxDimension=434*5*10^4,
  #maxDimension=100*10^4,
  bootB,
  bootLassoAlpha,
  binPredInd,
  covsPrefix,
  Mprefix,
  standardize,
  seed
){
  
  results=list()
  # load data info
  basicInfo=dataInfo(data=data,binPredInd=binPredInd,
                     covsPrefix=covsPrefix,Mprefix=Mprefix)
  taxaNames=basicInfo$taxaNames
  ii=which(basicInfo$taxaNames%in%refTaxa)
  
  nTaxa=basicInfo$nTaxa
  nPredics=basicInfo$nPredics
  rm(basicInfo)
  gc()
  
  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa
  
  countOfSelec=rep(0,nAlphaSelec)
  resultsByRefTaxon=list()
  
  # inital Lasso OLS estimate
  dataForEst=dataRecovTrans(data=data,ref=refTaxa,
                            Mprefix=Mprefix,covsPrefix=covsPrefix)
  
  x=as(as.matrix(dataForEst$xTildalong),"sparseMatrix")
  y=as(dataForEst$UtildaLong,"sparseVector")
  rm(dataForEst)
  print("dim(x):")
  print(dim(x))
  print("length(y):")
  print(length(y))
  
  xCol=ncol(x)
  
  maxSubSamplSiz=floor(maxDimension/xCol)
  nToSamplFrom=length(y)
  subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
  if(subSamplK==1)maxSubSamplSiz=nToSamplFrom
  cat("HDCI subSamplK:",subSamplK,"\n")
  
  nRuns=(ceiling(subSamplK/2))
  cat("HDCI nRuns:",nRuns,"\n")
  
  for(k in 1:nRuns){
    rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
    xSub=as((x[rowToKeep,]),"sparseMatrix")
    ySub=as((y[rowToKeep]),"sparseVector")
    if(k==1){
      print("dim(xSub):")
      print(dim(xSub))
      print("length(ySub):")
      print(length(ySub))
      cat("class(xSub): ",class(xSub),"\n")
      cat("class(ySub): ",class(ySub),"\n")
    }
    penal=runBootLassoHDCI(x=xSub,y=ySub,nPredics=nPredics,nTaxa=nTaxa,
                           refTaxaPosition=ii,bootLassoAlpha=bootLassoAlpha,
                           bootB=bootB,standardize=standardize,
                           seed=seed)
    rm(xSub,ySub)
    gc()
    finalBetaEst.k=penal$beta
    CIvecLow.k=penal$betaCIlow
    CIvecUp.k=penal$betaCIhi
    
    finalBetaEst.LPR.k=penal$beta.LPR
    CIvecLow.LPR.k=penal$betaCIlow.LPR
    CIvecUp.LPR.k=penal$betaCIhi.LPR
    rm(penal)
    if(k==1){
      finalBetaEst=finalBetaEst.k
      CIvecLow=CIvecLow.k
      CIvecUp=CIvecUp.k
      
      finalBetaEst.LPR=finalBetaEst.LPR.k
      CIvecLow.LPR=CIvecLow.LPR.k
      CIvecUp.LPR=CIvecUp.LPR.k
    }
    if(k>1){
      finalBetaEst=finalBetaEst+finalBetaEst.k
      CIvecLow=CIvecLow+CIvecLow.k
      CIvecUp=CIvecUp+CIvecUp.k
      
      finalBetaEst.LPR=finalBetaEst.LPR+finalBetaEst.LPR.k
      CIvecLow.LPR=CIvecLow.LPR+CIvecLow.LPR.k
      CIvecUp.LPR=CIvecUp.LPR+CIvecUp.LPR.k
    }
    if(k>0 & (k%%(ceiling(nRuns/10))==0)){
      message(floor(100*k/nRuns)," percent of the estimation analsis for refrence taxon ",refTaxa,
              " have been done.")
    }
  }
  rm(x,y)
  gc()
  
  results$finalBetaEst=finalBetaEst/nRuns
  results$CIvecLow=CIvecLow/nRuns
  results$CIvecUp=CIvecUp/nRuns
  
  results$finalBetaEst.LPR=finalBetaEst.LPR/nRuns
  results$CIvecLow.LPR=CIvecLow.LPR/nRuns
  results$CIvecUp.LPR=CIvecUp.LPR/nRuns
  
  return(results)
}