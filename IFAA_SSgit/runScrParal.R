
#---------------------------------------------------------------------------------------
## parallel screen function to select independent reference taxa
#---------------------------------------------------------------------------------------

runScrParal=function(
  method=c("mcp"),
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  nRef,
  paraJobs,
  x1permut,
  nPermu,
  refTaxa,
  maxDimensionScr=0.8*434*10*10^4,
  #maxDimensionScr=1000000,
  standardize,
  sequentialRun,
  allFunc=allFunc,
  refReadsThresh,
  SDThresh,
  SDquantilThresh,
  balanceCut,
  Mprefix,
  covsPrefix,
  binPredInd,
  seed){
  
  results=list()
  
  # load data info
  basicInfo=dataInfo(data=data,qualifyRefTax=T,
                     refReadsThresh=refReadsThresh,
                     SDThresh=SDThresh,SDquantilThresh=SDquantilThresh,
                     balanceCut=balanceCut,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binPredInd)
  
  taxaNames=basicInfo$taxaNames
  nTaxa=basicInfo$nTaxa
  nPredics=basicInfo$nPredics
  nSub=basicInfo$nSub
  predNames=basicInfo$predNames
  
  results$goodRefTaxaCandi=basicInfo$goodRefTaxaCandi
  rm(basicInfo)
  gc()
  
  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa
  
  # make reference taxa list
  if(length(refTaxa)<nRef){
    if(length(seed)>0){
      set.seed(as.numeric(seed))
    }
    refTaxa_extra=sample((taxaNames[!(taxaNames%in%refTaxa)]),(nRef-length(refTaxa)))
    refTaxa=c(refTaxa,refTaxa_extra)
    results$refTaxa=refTaxa
  }
  
  if(length(refTaxa)>=nRef){
    if(length(seed)>0){
      set.seed(as.numeric(seed))
    }
    refTaxa=sample(refTaxa,nRef)
    results$refTaxa=refTaxa
  }
  
  #
  ## run original data screen
  #
  screen1=originDataScreen(data=data,testCovInd=testCovInd,
                           nRef=nRef,refTaxa=refTaxa,
                           paraJobs=paraJobs,
                           method=method,allFunc=allFunc,Mprefix=Mprefix,
                           covsPrefix=covsPrefix,
                           binPredInd=binPredInd,
                           standardize=standardize,
                           sequentialRun=sequentialRun,
                           seed=seed)
  
  results$countOfSelecForAPred=screen1$countOfSelecForAPred
  yTildLongList=screen1$yTildLongList
  results$testCovCountMat=screen1$testCovCountMat
  Alist=screen1$Alist
  Wlist=screen1$Wlist
  twoList=screen1$twoList
  lengthTwoList=screen1$lengthTwoList
  rm(screen1)
  gc()
  print("x1permut:")
  print(x1permut)
  maxMemUsedInMb=sum(gc()[,6])
  message("maximum memory used after origin screen: ",maxMemUsedInMb," Mb")
  
  if(!x1permut){
    # generate Xbeta and residues for reduced model permutation
    XbetResi=XbetaAndResidu(data=data,testCovInd=testCovInd,
                            nRef=nRef,refTaxa=refTaxa,paraJobs=paraJobs,
                            method=method,allFunc=allFunc,Mprefix=Mprefix,
                            covsPrefix=covsPrefix,standardize=standardize,
                            binPredInd=binPredInd,
                            sequentialRun=sequentialRun,
                            seed=seed)
    xBetaList=XbetResi$xBetaList
    residuList=XbetResi$residuList
    xTildLong=XbetResi$xTildLong
    rm(XbetResi)
  }
  #
  ## start to run permutation screen
  #
    startT=proc.time()[3]
    message("start to run permutation")
    
    # permut the exposure variable
    if(length(seed)>0)set.seed(as.numeric(seed)+10^6)
    
    permutOrder=lapply(rep(nSub,nPermu),sample)
    if(!x1permut){
      residuPermuOrder=lapply(rep(length(residuList[[1]]),nPermu),sample)
    }
    
    screenStartTime = proc.time()[3]
    
    EName=testCovInNewNam
    EVar=data[,EName,drop=F]
    
    totNumOfLoops=nRef*nPermu
    
    if(length(paraJobs)==0){
      availCores=availableCores()
      if(is.numeric(availCores))paraJobs=max(1,availableCores()-2)
      if(!is.numeric(availCores))paraJobs=1
    }
    
    if(!sequentialRun){
      message(paraJobs, " parallel jobs are registered for the permutation analysis in Phase 1b")
    }
    
    batch=paraJobs
    forLoopN=ceiling(totNumOfLoops/batch)
    print("permut forLoopN:")
    print(forLoopN)
    
    for (jj in 1:forLoopN){
      cl <- parallel::makeCluster(paraJobs,outfile="")

      parallel::clusterExport(cl=cl, varlist=allFunc,envir=parent.env(environment()))
      doParallel::registerDoParallel(cl)
      
      if(sequentialRun){foreach::registerDoSEQ()}
      if(forLoopN>1 & jj<forLoopN){
        refResu.j=foreach(i=((jj-1)*batch+1):(jj*batch),
                          .multicombine=T,
                          .packages=c("picasso","Matrix","oem"),
                          .errorhandling="pass") %dopar% {
            # for(i in ((jj-1)*paraJobs+1):(jj*paraJobs)){
                            if((i%%10)==0){
                              cat("permutation i:",i,"\n")
                            }
                            if(i>totNumOfLoops | i<1){
                              cat("i is out of range.","\n")
                            }
                            permut.i=1+(i-1)%%nPermu
                            ref.i=1+floor((i-1)/nPermu)
                            ii=which(taxaNames==refTaxa[ref.i])
                            if(i==1){
                              print("ref.i:")
                              print(ref.i)
                            }
                            if(x1permut){
                              permutX1=EVar[permutOrder[[permut.i]],,drop=F]
                              newData=data
                              newData[,EName]=permutX1
                              rm(permutX1)
                              
                              xLong.i=dataRecovTrans(data=newData,ref=refTaxa[ref.i],
                                                     Mprefix=Mprefix,covsPrefix=covsPrefix,xOnly=TRUE)
                              xLongTild.i=xLong.i$xTildalong
                              rm(newData,xLong.i)
                              
                              if(i==1){
                                cat("class(xLongTild.i): ",class(xLongTild.i),"\n")
                                cat("Post xTildL function in permutation","\n")
                              }
                              maxMemUsedInMb=sum(gc()[,6])
                              if(i%%100==1)cat("maximum memory used after xTildL function: ",maxMemUsedInMb," Mb","\n")
                              gc()
                              maxSubSamplSiz=floor(maxDimensionScr/ncol(xLongTild.i))
                              nToSamplFrom=nrow(xLongTild.i)
                              
                              if(method=="mcp") {
                                subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                                if(subSamplK==1)maxSubSamplSiz=nToSamplFrom
                                
                                nRuns=ceiling(subSamplK/2)
                                
                                if(i==1){
                                  cat("Permutation subSamplK:",subSamplK,"\n")
                                  cat("nToSamplFrom:",nToSamplFrom,"\n")
                                  cat("maxSubSamplSiz:",maxSubSamplSiz,"\n")
                                  cat("nRuns:",nRuns,"\n")
                                }
                                
                                for (k in 1:nRuns){
                                  rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                                  ## Shang: convert to big.matrix
                                  x=as.big.matrix(as.matrix(xLongTild.i[rowToKeep,]))
                                  # x=as((xLongTild.i[rowToKeep,]),"sparseMatrix")
                                  y=yTildLongList[[ref.i]][rowToKeep]
                                  # Penal.i=runPicasso(x=x,y=y,
                                  #                    nPredics=nPredics,
                                  #                    method="mcp",permutY=permutY,
                                  #                    standardize=standardize,
                                  #                    seed=seed,seedi=i)
                                  Penal.i=runOEM(x=x,y=y,
                                                 nPredics=nPredics,standardize=standardize,
                                                 seed=seed,seedi=i)
                                  
                                  BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                                  rm(Penal.i)
                                  gc()
                                  if(k==1)BetaNoInt.i=BetaNoInt.k
                                  if(k>1)BetaNoInt.i=BetaNoInt.i+BetaNoInt.k
                                  rm(BetaNoInt.k)
                                }
                                rm(k,x,y,xLongTild.i)
                                gc()
                                BetaNoInt.i=BetaNoInt.i/nRuns
                              }
                              gc()
                            }else{
                              resid=residuList[[ref.i]][residuPermuOrder[[permut.i]]]
                              yTildLong.i=xBetaList[[ref.i]]+resid
                              rm(resid)

                              maxSubSamplSiz=floor(maxDimensionScr/ncol(xTildLong))
                              nToSamplFrom=nrow(xTildLong)
                              if(method=="mcp") {
                                subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                                if(subSamplK==1)maxSubSamplSiz=nToSamplFrom
                                cat("subSamplK:",subSamplK,"\n")
                                
                                nRuns=ceiling(subSamplK/2)
                                cat("nRuns:",nRuns,"\n")
                                
                                for (k in 1:nRuns){
                                  rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                                  x=as((xTildLong[rowToKeep,]),"sparseMatrix")
                                  y=as((yTildLong.i[rowToKeep]),"sparseVector")
                                  # Penal.i=runPicasso(x=x,y=y,
                                  #                    nPredics=nPredics,
                                  #                    method="mcp",permutY=permutY,
                                  #                    standardize=standardize,
                                  #                    seed=seed,seedi=i)
                                  
                                  Penal.i=runOEM(x=x,y=y,
                                                 nPredics=nPredics,standardize=standardize,
                                                 seed=seed,seedi=i)
                                  BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                                  rm(Penal.i)
                                  gc()
                                  if(k==1)BetaNoInt.i=BetaNoInt.k
                                  if(k>1)BetaNoInt.i=BetaNoInt.i+BetaNoInt.k
                                  rm(BetaNoInt.k)
                                }
                                rm(k,x,y,xTildLong)
                                gc()
                                BetaNoInt.i=BetaNoInt.i/nRuns
                              }
                              rm(yTildLong.i)
                            }
                            selection.i=as(rep(0,nAlphaSelec),"sparseVector")
                            if (ii==1){
                              selection.i[-seq(1,nPredics)]=BetaNoInt.i
                            }
                            if (ii==nTaxa) {
                              selection.i[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=BetaNoInt.i
                            }
                            if ((ii>1) & (ii<nTaxa)) {
                              selection.i[1:(nPredics*(ii-1))]=BetaNoInt.i[1:(nPredics*(ii-1))]
                              selection.i[(nPredics*ii+1):nAlphaSelec]=BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]
                            }
                            rm(BetaNoInt.i)
                            # create return vector
                            recturnVec=as(rep(0,nAlphaSelec),"sparseVector")
                            recturnVec=selection.i
                            rm(selection.i)
                            return(recturnVec)
                          }
        # rm(yTildLongList)
        parallel::stopCluster(cl)
        gc()
        if(jj==1)refResu=refResu.j
        if(jj>1)refResu=do.call(c,list(refResu,refResu.j))
      }
      if(jj==forLoopN){
        refResu.j=foreach(i=((forLoopN-1)*batch+1):totNumOfLoops,
                          .multicombine=T,
                          .packages=c("picasso","Matrix","oem"),
                          .errorhandling="pass") %dopar% {
                            
                            # for(i in ((forLoopN-1)*batch+1):totNumOfLoops){
                            if((i%%10)==0){
                              cat("permutation i:",i,"\n")
                            }
                            if(i>totNumOfLoops | i<1){
                              cat("i is out of range.","\n")
                            }
                            permut.i=1+(i-1)%%nPermu
                            ref.i=1+floor((i-1)/nPermu)
                            ii=which(taxaNames==refTaxa[ref.i])
                            
                            if(x1permut){
                              permutX1=EVar[permutOrder[[permut.i]],,drop=F]
                              newData=data
                              newData[,EName]=permutX1
                              rm(permutX1)
                              
                              xLong.i=dataRecovTrans(data=newData,ref=refTaxa[ref.i],
                                                     Mprefix=Mprefix,covsPrefix=covsPrefix,xOnly=TRUE)
                              xLongTild.i=xLong.i$xTildalong
                              rm(newData,xLong.i)
                              
                              if(i==1){
                                cat("class(xLongTild.i): ",class(xLongTild.i),"\n")
                                cat("Post xTildL function in permutation","\n")
                              }
                              gc()
                              maxMemUsedInMb=sum(gc()[,6])
                              if(i%%100==1)cat("maximum memory used after xTildL function: ",maxMemUsedInMb," Mb","\n")
                              
                              maxSubSamplSiz=floor(maxDimensionScr/ncol(xLongTild.i))
                              nToSamplFrom=nrow(xLongTild.i)
                              if(method=="mcp") {
                                subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                                if(subSamplK==1)maxSubSamplSiz=nToSamplFrom
                                nRuns=ceiling(subSamplK/2)
                                
                                if(i==1){
                                  cat("Permutation subSamplK:",subSamplK,"\n")
                                  cat("Permutation nRuns:",nRuns,"\n")
                                  cat("nToSamplFrom:",nToSamplFrom,"\n")
                                  cat("maxSubSamplSiz:",maxSubSamplSiz,"\n")
                                }
                                
                                for (k in 1:nRuns){
                                  rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                                  x=as.big.matrix(as.matrix(xLongTild.i[rowToKeep,]))
                                  y=yTildLongList[[ref.i]][rowToKeep]
                                  # Penal.i=runPicasso(x=x,y=y,
                                  #                    nPredics=nPredics,
                                  #                    method="mcp",permutY=permutY,
                                  #                    standardize=standardize,
                                  #                    seed=seed,seedi=i)
                                  
                                  Penal.i=runOEM(x=x,y=y,
                                                 nPredics=nPredics,standardize=standardize,
                                                 seed=seed,seedi=i)
                                  
                                  BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                                  rm(Penal.i)
                                  gc()
                                  if(k==1)BetaNoInt.i=BetaNoInt.k
                                  if(k>1)BetaNoInt.i=BetaNoInt.i+BetaNoInt.k
                                  rm(BetaNoInt.k)
                                }
                                rm(k,x,y,xLongTild.i)
                                gc()
                                BetaNoInt.i=BetaNoInt.i/nRuns
                              }
                              gc()
                            }else{
                              resid=residuList[[ref.i]][residuPermuOrder[[permut.i]]]
                              yTildLong.i=xBetaList[[ref.i]]+resid
                              rm(resid)
                              
                              maxSubSamplSiz=floor(maxDimensionScr/ncol(xTildLong))
                              nToSamplFrom=nrow(xTildLong)
                              if(method=="mcp") {
                                subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                                if(subSamplK==1)maxSubSamplSiz=nToSamplFrom
                                cat("subSamplK:",subSamplK,"\n")
                                
                                nRuns=ceiling(subSamplK/2)
                                cat("nRuns:",nRuns,"\n")
                                
                                for (k in 1:nRuns){
                                  rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                                  x=xTildLong[rowToKeep,]
                                  y=yTildLong.i[rowToKeep]
                                  # Penal.i=runPicasso(x=x,y=y,
                                  #                    nPredics=nPredics,
                                  #                    method="mcp",permutY=permutY,
                                  #                    standardize=standardize,
                                  #                    seed=seed,seedi=i)
                                  
                                  Penal.i=runOEM(x=x,y=y,
                                                 nPredics=nPredics,standardize=standardize,
                                                 seed=seed,seedi=i)
                                  BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                                  rm(Penal.i)
                                  gc()
                                  if(k==1)BetaNoInt.i=BetaNoInt.k
                                  if(k>1)BetaNoInt.i=BetaNoInt.i+BetaNoInt.k
                                  rm(BetaNoInt.k)
                                }
                                rm(k,x,y)
                                gc()
                                BetaNoInt.i=BetaNoInt.i/nRuns
                              }
                              rm(xTildLong,yTildLong.i)
                            }
                            
                            selection.i=as(rep(0,nAlphaSelec),"sparseVector")
                            if (ii==1){
                              selection.i[-seq(1,nPredics)]=BetaNoInt.i
                            }
                            if (ii==nTaxa) {
                              selection.i[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=BetaNoInt.i
                            }
                            if ((ii>1) & (ii<nTaxa)) {
                              selection.i[1:(nPredics*(ii-1))]=BetaNoInt.i[1:(nPredics*(ii-1))]
                              selection.i[(nPredics*ii+1):nAlphaSelec]=BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]
                            }
                            rm(BetaNoInt.i)
                            # create return vector
                            recturnVec=as(rep(0,nAlphaSelec),"sparseVector")
                            recturnVec=selection.i
                            rm(selection.i)
                            return(recturnVec)
                          }
        # rm(yTildLongList)
        parallel::stopCluster(cl)
        gc()
        if(forLoopN==1)refResu=refResu.j
        if(forLoopN>1)refResu=do.call(c,list(refResu,refResu.j))
      }
      
      if(jj>0 & (jj%%(ceiling(forLoopN/10))==0)){
        message(round(100*jj/forLoopN,0), " percent of phase 1b analysis has been done.")
      }
    }
    rm(Alist,Wlist,yTildLongList,twoList,lengthTwoList)
    gc()
    maxMemUsedInMb=sum(gc()[,6])
    message("maximum memory used after permutation screen: ",maxMemUsedInMb," Mb")
    refResu<- lapply(refResu, as, "sparseMatrix")
    refResu=do.call(cbind, refResu)
    
    rm(data,permutOrder,EVar)
    
    endT=proc.time()[3]
    
    message("Permutation analysis done and took ", (endT-startT)/60," minutes")
    
    # obtain the maximum vector
    permuColInd=1+seq(0,totNumOfLoops-1)%%nPermu
    permutResuMat=matrix(NA,nrow=nTaxa,ncol=nPermu)
    permutTestCovList=list()
    nTestCov=length(testCovInd)
    
    for(i in 1:nPermu){
      matrix.i.permu=refResu[,(permuColInd==i),drop=F]
      vec.i=as(rep(0,nTaxa),"sparseVector")
      
      for(j in 1:nRef){
        vector.j.ref=matrix.i.permu[,j]
        matrix.i.j=as(matrix(vector.j.ref,nrow=nPredics),"sparseMatrix")
        #testCovVec.j.ref=as((Matrix::colSums(matrix.i.j[testCovInd,,drop=F])>0)+0,"sparseVector")
        testCovVec.j.ref=as((Matrix::colSums(matrix.i.j[testCovInd,,drop=F])),"sparseVector")
        
        vec.i=vec.i+testCovVec.j.ref
      }
      
      permutResuMat[,i]=as.vector(vec.i)
    }
    rm(vec.i,matrix.i.permu,vector.j.ref,matrix.i.j,testCovVec.j.ref)
    
    if(nTestCov>1){
      for(i in 1:nPermu){
        matrix.i.permu=refResu[,(permuColInd==i),drop=F]
        allCovCountMat.i=as(matrix(Matrix::rowSums(matrix.i.permu),nrow=nPredics),"sparseMatrix")
        permutTestCovList[[i]]=allCovCountMat.i[testCovInd,,drop=F]
      }
      rm(matrix.i.permu,allCovCountMat.i)
    }
    
    # print("permutResuMat:")
    # print(permutResuMat)
    
    maxVec=rep(NA,nPermu)
    for(i in 1:nPermu){
      maxVec[i]=max(permutResuMat[,i])
    }
    rm(permutResuMat)
    
    results$nTestCov=nTestCov
    
    if(nTestCov>1){
      MaxMatTestCovByPermu=matrix(NA,nrow=nTestCov,ncol=nPermu)
      for(k in 1:nTestCov){
        for(i in 1:nPermu){
          MaxMatTestCovByPermu[k,i]=max(permutTestCovList[[i]][k,])
        }
      }
      rm(permutTestCovList)
      
      results$MaxMatTestCovByPermu=MaxMatTestCovByPermu
      rm(MaxMatTestCovByPermu)
    }
    
    # obtain the null binomial distribution for each taxa
    totSeleCount=Matrix::rowSums(refResu)
    rm(refResu)
    binomPar=totSeleCount/totNumOfLoops
    rm(totSeleCount,totNumOfLoops)
    
    results$maxVec=maxVec
    rm(maxVec)
    
    results$binomPar=binomPar
    rm(binomPar)
    results$nTaxa=nTaxa
    results$nPredics=nPredics

  results$taxaNames=taxaNames
  rm(taxaNames)
  return(results)
}

