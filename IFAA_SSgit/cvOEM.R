
cvOEM=function(
  x,
  y,
  lamList,
  nfolds,
  zeroSDCut,
  intercept,
  seed,
  seedi
){

  results=list()

  nLam=length(lamList)
    
  nObsAll=length(y)
  nBeta=ncol(x)

  # partition the data randomly into nfolds partitions
  parSize=floor(nObsAll/nfolds)

  if(length(seed)>0)set.seed(as.numeric(seed)+10^7+seedi)

  randomShuf=sample(nObsAll, nObsAll)
  sampleInd=list()
  for (i in 1:nfolds){
    if(i<nfolds){
      sampleInd[[i]]=randomShuf[(1+(i-1)*parSize):(i*parSize)]
    }else{
      sampleInd[[i]]=randomShuf[(1+(i-1)*parSize):nObsAll]
    }
  }
  rm(randomShuf)

  # cross validation
  cvPara=matrix(NA,nrow=nLam,ncol=nfolds)

  for(i in 1:nfolds){
    # check if there is zero-variance x in the partitions
    startT.i=proc.time()[3]
    xi=x[-sampleInd[[i]],]
    sdX.i=apply(xi,2,sd)
    xWithNearZeroSd.i=which(sdX.i<=zeroSDCut)
    rm(sdX.i)

    # remove near constant columns in x
    if(length(xWithNearZeroSd.i)>0) {
      xi=xi[,-xWithNearZeroSd.i]
    }
    nearZeroSd.i=length(xWithNearZeroSd.i)
    xi <- as.big.matrix(xi)
    yi=as.vector(y[-sampleInd[[i]]])

    ## Shang: use big.oem
    cv.i=big.oem(x=xi,y=yi,penalty=c("mcp"),gamma=3,intercept=intercept,
             standardize=FALSE,lambda=lamList,tol=1e-10)
    rm(xi,yi)

    if(length(xWithNearZeroSd.i)>0) {
      missPositions=xWithNearZeroSd.i
      for(iLam in 2:nLam){
        missPositions=c(missPositions,xWithNearZeroSd.i+(iLam-1)*nBeta)
      }
      rm(xWithNearZeroSd.i)
      gc()
      betaa=as.vector(cv.i$beta$mcp[-1,])
      betaTrans=groupBetaToFullBeta(nTaxa=(nBeta*nLam),nPredics=1,
                                    unSelectList=sort(missPositions),newBetaNoInt=betaa)
      rm(betaa)
      gc()
      betaMatrix=Matrix(betaTrans$finalBeta,nrow=nBeta,ncol=nLam)
      rm(betaTrans)
    } else {betaMatrix=cv.i$beta$mcp[-1,]}

    # training is done, start testing
    interc=cv.i$beta$mcp[1,]
    if(sum(is.na(interc))>0){
      interc[is.na(interc)]=mean(interc[!is.na(interc)])
    }
    if(sum(is.na(interc))>nLam*0.15){
      warnings("Many intercept estimates are missing in picasso estimates")
    }
    intcep=rep(1,length(sampleInd[[i]]))%*%t(as.matrix(interc))
    rm(cv.i,interc)
    
    ## Shang: submatrix of x in big.matrix format and the product
    
    yHatTest.i=x[sampleInd[[i]],]%*%betaMatrix+intcep
    rm(betaMatrix,intcep)
    gc()

    resiVecs.i=(yHatTest.i-y[sampleInd[[i]]])
    rm(yHatTest.i)
    cvPara[,i]=apply(resiVecs.i,2,function(x)sum(x^2))
    rm(resiVecs.i)
  }
  rm(sampleInd,x,y)

  SSE=Matrix::rowSums(cvPara)
  rm(cvPara)

  optiLamInd=tail(which(SSE==min(SSE)),n=1)
  optiLam=lamList[optiLamInd]

  # return results

  results$optiLam=optiLam
  rm(optiLam)

  return(results)
}

