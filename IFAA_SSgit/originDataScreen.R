
#---------------------------------------------------------------------------------------
## Original screen function to select independent reference taxa
#---------------------------------------------------------------------------------------
originDataScreen=function(
  method,
  data,
  testCovInd,
  nRef,
  paraJobs,
  refTaxa,
  maxDimensionScr=0.8*434*10*10^4,
  #maxDimensionScr=1000000,
  standardize,
  sequentialRun,
  allFunc,
  Mprefix,
  covsPrefix,
  binPredInd,
  seed){
  
  results=list()
  
  # load data info
  basicInfo=dataInfo(data=data,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binPredInd)
  
  taxaNames=basicInfo$taxaNames
  nTaxa=basicInfo$nTaxa
  nPredics=basicInfo$nPredics
  rm(basicInfo)
  gc()
  
  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa
  
  countOfSelec=rep(0,nAlphaSelec)
  resultsByRefTaxon=list()
  
  # overwrite nRef if the reference taxon is specified
  nRef=length(refTaxa)
  
  startT1=proc.time()[3]
  message("start Original screen")
  if(length(paraJobs)==0){
    availCores=availableCores()
    if(is.numeric(availCores))paraJobs=max(1,availableCores()-2)
    if(!is.numeric(availCores))paraJobs=1
  }
  batch=paraJobs
  forLoopN=ceiling(nRef/batch)

  if(!sequentialRun){
    message(paraJobs, " parallel jobs are registered for analyzing ", nRef, " reference taxa in Phase 1a")
  }
  
  for (jj in 1:forLoopN){

    cl<-parallel::makeCluster(paraJobs,outfile="")
    
    parallel::clusterExport(cl=cl, varlist=allFunc,envir=parent.env(environment()))
    doParallel::registerDoParallel(cl)
    
    if(sequentialRun){foreach::registerDoSEQ()}
    
    #start parallel computing
    if(forLoopN>1 & jj<forLoopN){
      # scr1Resu.j=foreach(i=((jj-1)*batch+1):(jj*batch),
      #                    .multicombine=T,
      #                    .packages=c("picasso","Matrix","oem"),
      #                    .errorhandling="pass") %dopar% {
                           
       for(i in 1:((jj-1)*paraJobs+1):(jj*paraJobs)){
                           if(i==1){
                             cat("originScreen i:", i,"\n")
                             if(i>nRef | i<1){
                               cat("i is out of range.","\n")
                             }
                             cat("refTaxa[i]:",refTaxa[i],"\n")
                           }
                           
                           ii=which(taxaNames==refTaxa[i])
                           dataForEst=dataRecovTrans(data=data,ref=refTaxa[i],Mprefix=Mprefix,
                                                     covsPrefix=covsPrefix)
                           
                           ## Shang: xTildLongTild.i needs to be big.matrix format
                           
                           xTildLongTild.i=dataForEst$xTildalong
                           yTildLongTild.i=dataForEst$UtildaLong
                           rm(dataForEst)
                           gc()
                           if(i==1){
                             cat("class(xTildLongTild.i): ",class(xTildLongTild.i),"\n")
                             cat("class(yTildLongTild.i): ",class(yTildLongTild.i),"\n")
                             print("dim(xTildLongTild.i):")
                             print(dim(xTildLongTild.i))
                             print("length(yTildLongTild.i):")
                             print(length(yTildLongTild.i))
                           }
                           maxSubSamplSiz=floor(maxDimensionScr/ncol(xTildLongTild.i))
                           nToSamplFrom=nrow(xTildLongTild.i)
                           if(i==1){
                             print("maxSubSamplSiz:")
                             print(maxSubSamplSiz)
                             print("nToSamplFrom:")
                             print(nToSamplFrom)
                           }
                           
                           if(method=="mcp") {
                             subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                             if(subSamplK==1)maxSubSamplSiz=nToSamplFrom
                             if(i==1)cat("Origin Screen subSamplK:",subSamplK,"\n")
                             
                             nRuns=ceiling(subSamplK/2)
                             if(i==1)cat("Origin Screen nRuns:",nRuns,"\n")
                             
                             for (k in 1:nRuns){
                               rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                               
                               ## Shang: x needs to be in big.matrix format
                               
                               # x=as((xTildLongTild.i[rowToKeep,]),"sparseMatrix")
                               x=as.big.matrix(xTildLongTild.i[rowToKeep,])
                               y=as.vector(yTildLongTild.i[rowToKeep])
                               # Penal.i=runPicasso(x=x,y=y,
                               #                    nPredics=nPredics,
                               #                    method="mcp",permutY=F,
                               #                    standardize=standardize,
                               #                    seed=seed,seedi=i)
                               
                               Penal.i=runOEM(x=x,y=y,
                                              nPredics=nPredics,standardize=standardize,
                                              seed=seed,seedi=i)
                               
                               # print("Penal.i$betaNoInt:")
                               # print(Penal.i$betaNoInt)
                               
                               BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                               rm(Penal.i)
                               if(k==1)BetaNoInt.i=BetaNoInt.k
                               if(k>1)BetaNoInt.i=BetaNoInt.i+BetaNoInt.k
                               rm(BetaNoInt.k)
                             }
                             rm(k,x,y,xTildLongTild.i)
                             gc()
                             BetaNoInt.i=BetaNoInt.i/nRuns
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
                           recturnlist=list()
                           recturnlist[[1]]=selection.i
                           recturnlist[[2]]=yTildLongTild.i
                           rm(selection.i,yTildLongTild.i)
                           return(recturnlist)
                         }
      parallel::stopCluster(cl)
      gc()
      if(jj==1)scr1Resu=scr1Resu.j
      if(jj>1)scr1Resu=do.call(c,list(scr1Resu,scr1Resu.j))
    }
    
    if(jj==forLoopN){
      scr1Resu.j=foreach(i=((forLoopN-1)*batch+1):nRef,
                         .multicombine=T,
                         .packages=c("picasso","Matrix","oem"),
                         .errorhandling="pass") %dopar% {
                           
                           # for(i in ((forLoopN-1)*paraJobs+1):nRef){
                           
                           if(i==1)cat("originScreen i:", i,"\n")
                           
                           ii=which(taxaNames==refTaxa[i])
                           dataForEst=dataRecovTrans(data=data,ref=refTaxa[i],Mprefix=Mprefix,
                                                     covsPrefix=covsPrefix)
                           
                           xTildLongTild.i=dataForEst$xTildalong
                           yTildLongTild.i=dataForEst$UtildaLong
                           rm(dataForEst)
                           gc()
                           if(i==1){
                             cat("class(xTildLongTild.i): ",class(xTildLongTild.i),"\n")
                             cat("class(yTildLongTild.i): ",class(yTildLongTild.i),"\n")
                             print("dim(xTildLongTild.i):")
                             print(dim(xTildLongTild.i))
                             print("length(yTildLongTild.i):")
                             print(length(yTildLongTild.i))
                           }
                           maxSubSamplSiz=floor(maxDimensionScr/ncol(xTildLongTild.i))
                           nToSamplFrom=nrow(xTildLongTild.i)
                           if(i==1){
                             print("maxSubSamplSiz:")
                             print(maxSubSamplSiz)
                             print("nToSamplFrom:")
                             print(nToSamplFrom)
                           }
                           if(method=="mcp") {
                             
                             subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                             if(subSamplK==1)maxSubSamplSiz=nToSamplFrom
                             if(i==1)cat("Origin Screen subSamplK:",subSamplK,"\n")
                             
                             nRuns=ceiling(subSamplK/2)
                             if(i==1)cat("Origin Screen nRuns:",nRuns,"\n")
                             
                             for (k in 1:nRuns){
                               rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                               # Shang
                               x=as.big.matrix(xTildLongTild.i[rowToKeep,])
                               y=as.vector(yTildLongTild.i[rowToKeep])
                               # Penal.i=runPicasso(x=x,y=y,
                               #                   nPredics=nPredics,
                               #                   method="mcp",permutY=F,
                               #                   standardize=standardize,
                               #                   seed=seed,seedi=i)
                               
                               Penal.i=runOEM(x=x,y=y,
                                              nPredics=nPredics,standardize=standardize,
                                              seed=seed,seedi=i)
                               
                               # print("Penal.i$betaNoInt:")
                               # print(Penal.i$betaNoInt)
                               BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                               rm(Penal.i)
                               gc()
                               if(k==1)BetaNoInt.i=BetaNoInt.k
                               if(k>1)BetaNoInt.i=BetaNoInt.i+BetaNoInt.k
                               rm(BetaNoInt.k)
                             }
                             rm(k,x,y,xTildLongTild.i)
                             gc()
                             BetaNoInt.i=BetaNoInt.i/nRuns
                           }
                           
                           selection.i=as(rep(0,nAlphaSelec),"sparseVector")
                           if (ii==1){selection.i[-seq(1,nPredics)]=BetaNoInt.i
                           }
                           if (ii==nTaxa) {selection.i[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=BetaNoInt.i
                           }
                           if ((ii>1) & (ii<nTaxa)) {
                             selection.i[1:(nPredics*(ii-1))]=BetaNoInt.i[1:(nPredics*(ii-1))]
                             selection.i[(nPredics*ii+1):nAlphaSelec]=BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]
                           }
                           rm(BetaNoInt.i)
                           # create return vector
                           recturnlist=list()
                           recturnlist[[1]]=selection.i
                           recturnlist[[2]]=yTildLongTild.i
                           rm(selection.i,yTildLongTild.i)
                           return(recturnlist)
                         }
      parallel::stopCluster(cl)
      if(forLoopN==1)scr1Resu=scr1Resu.j
      if(forLoopN>1)scr1Resu=do.call(c,list(scr1Resu,scr1Resu.j))
      gc()
    }
    # print("scr1Resu:")
    # print(scr1Resu)
    if(jj>0 & (jj%%(ceiling(forLoopN/10))==0)){
      message(round(100*jj/forLoopN,0), "percent of phase 1a analysis has been done")
    }
  }
  rm(data)
  endT=proc.time()[3]
  
  message("Original screen done and took",(endT-startT1)/60,"minutes")
  
  selecList=list()
  for(i in 1:nRef){
    selecList[[i]]=scr1Resu[[i]][[1]]
  }
  
  results$yTildLongList=list()
  for(i in 1:nRef){
    results$yTildLongList[[i]]=scr1Resu[[i]][[2]]
  }
  print(selecList)
  rm(scr1Resu)
  
  selecList<- lapply(selecList, as, "sparseMatrix")
  scr1ResuSelec=do.call(cbind, selecList)
  rm(selecList)
  
  results$scr1ResuSelec=scr1ResuSelec
  
  # create count of selection for individual testCov
  countOfSelecForAllPred=as(matrix(Matrix::rowSums(scr1ResuSelec),nrow=nPredics),"sparseMatrix")
  testCovCountMat=countOfSelecForAllPred[testCovInd,,drop=F]
  rm(testCovInd,countOfSelecForAllPred)
  
  # create overall count of selection for all testCov as a whole
  countOfSelecForAPred=as(matrix(rep(0,nTaxa),nrow=1),"sparseMatrix")
  for (tax in 1:nTaxa){
    countMatForTaxni=scr1ResuSelec[(1+(tax-1)*nPredics):(tax*nPredics),,drop=F]
    totCountVecForTaxoni=Matrix::colSums(countMatForTaxni)
    #countOfSelecForAPred[1,tax]=sum(totCountVecForTaxoni>0)
    countOfSelecForAPred[1,tax]=sum(totCountVecForTaxoni)
    
  }
  rm(tax,scr1ResuSelec,countMatForTaxni,totCountVecForTaxoni)
  gc()
  
  # print("countOfSelecForAPred:")
  # print(countOfSelecForAPred)
  
  colnames(countOfSelecForAPred)=taxaNames
  rm(taxaNames)
  
  # return results
  results$testCovCountMat=testCovCountMat
  rm(testCovCountMat)
  results$countOfSelecForAPred=countOfSelecForAPred
  rm(countOfSelecForAPred)
  return(results)
}
