champ.DMR<-function(beta.norm=NULL, pd=NULL, lassoStyle="max", lassoRadius=2000, minSigProbesLasso=3, minDmrSep=1000, 
                    minDmrSize=0, adjPVal=0.05, adjust.method="BH", DMRpval=0.05, 
                    minPVal=0.001, # added this
                    noSigContinue=TRUE # added this
){

# Assuming the following about the input beta:
# 1. No XY chromosomes
# 2. "Sample_Group" in the phenotype/samplesheet file has to be "C" for controls and "T" for cases

  data(probe.features)
  data(probe.450K.VCs.af)
  
  message("Run Probe Lasso DMR Hunter")
  resultsFile=champ.MVP(beta.norm=beta.norm, pd=pd, adjPVal=adjPVal, 
                        bedFile=bedFile, adjust.method=adjust.method, 
                        resultsDir=getwd())
  colnames(resultsFile)
  
  if(sum(resultsFile$adj.P.Val<=adjPVal)>0){
    message(paste("There are ", sum(resultsFile$adj.P.Val<=adjPVal), " FDR significant hits. Ending Lasso analysis"))
    return(emptyDMRList)
  } else if(!noSigContinue){
    message("Not continuing with Lasso per input.")
  } else if(noSigContinue){
    message("Continuing with Lasso using unadjusted p-values")
    colnames(resultsFile)
    resultsFile$origAdjP<-resultsFile$adj.P.Val
    resultsFile$adj.P.Val<-resultsFile$P.Value
    adjPVal<-minPVal
  }
  
  myResults <- data.frame("P.Value" = resultsFile$P.Value, "adj.P.Val" = resultsFile$adj.P.Val, probe.features[match(rownames(resultsFile), rownames(probe.features)), ], row.names = rownames(resultsFile))
  myResults$CHR <- as.character(myResults$CHR)
  
  
  library(plyr)
  
  emptyDMRList=as.data.frame(setNames(replicate(21,numeric(0),simplify=F),c("probeID","adj.P.Val","CHR","MAPINFO","arm","gene.1","feature","cgi","feat.cgi","pol.af.f","pol.af.r","lasso.radius","dmr.no","dmr.start","dmr.end","dmr.size","dmr.core,start","dmr.core.end","dmr.core.size","dmr.p","deltaBeta")))
  
  # Checkign that there are still DMRs at the regular p-value level
  if(dim(myResults)[1]==0){
    message("Your dataset is empty and champ.lasso() cannot proceed")
    return(emptyDMRList)
  } else if(min(myResults$adj.P.Val)>adjPVal) {
    
    message("The adusted p-values in your dataset exceed the cutoff you have chosen champ.lasso() cannot proceed. You might like to rerun champ.lasso with the normalised beta values without batchCorrect.")
    return(emptyDMRList)
  }else if(count(myResults$adj.P.Val<adjPVal)[which(count(myResults$adj.P.Val<adjPVal)$x==TRUE),][1,2] < 3){
    message("There are not enough MVPs in your dataset for champ.lasso() to proceed.")
    return(emptyDMRList)
  }
  myResults <- myResults[order(myResults$CHR, myResults$MAPINFO),]
  
  ###Probe spacing and quantile derivation
  myResults.split <- split(myResults, paste(myResults$CHR, myResults$arm), drop=T) #split by chromosome & arm
  probe.spacing <- lapply(myResults.split, function(x) apply(cbind(c(diff(x$MAPINFO), tail(diff(x$MAPINFO), n = 1)), c(head(diff(x$MAPINFO), n = 1), diff(x$MAPINFO))), 1, min))	# nearest probe calculations
  myResults <- data.frame(do.call(rbind, myResults.split), "nrst.probe" = unlist(probe.spacing))
  rm(myResults.split, probe.spacing)
  
  lasso.quantiles <- do.call(rbind, lapply(split(myResults$nrst.probe, myResults$feat.cgi), function(x) ecdf(x)(lassoRadius)))
  if(lassoStyle == "max")
  {
    value.lasso.quantile <- min(lasso.quantiles)
  }else{
    value.lasso.quantile <- max(lasso.quantiles)
  }
  rm(lasso.quantiles)
  
  lasso.radius <- round(do.call(rbind, lapply(split(myResults$nrst.probe, myResults$feat.cgi), function(x) quantile(x, value.lasso.quantile, na.rm=T))))
  
  myResults$lasso.radius <- lasso.radius[match(myResults$feat.cgi, rownames(lasso.radius))]
  
  ### Filter out insignificant results
  myResults.sig <- myResults[myResults$adj.P.Val < adjPVal,]
  
  ### Determine probes that overlap based on the lasso radius
  lasso.gr <- GRanges(seqnames=paste("chr", myResults.sig$CHR, sep = ""), ranges=IRanges(start=myResults.sig$MAPINFO - myResults.sig$lasso.radius, end=myResults.sig$MAPINFO + myResults.sig$lasso.radius))
  probe.gr <- GRanges(seqnames=paste("chr", myResults.sig$CHR, sep = ""), ranges=IRanges(start=myResults.sig$MAPINFO, end=myResults.sig$MAPINFO))
  lasso.probe.countOverlap <- countOverlaps(lasso.gr, probe.gr)
  myResults.sig <- data.frame(myResults.sig, lasso.probe.countOverlap)
  probeKeepers <- which(lasso.probe.countOverlap >= minSigProbesLasso)
  
  if(length(probeKeepers)==0){
    message("No Overlapping probes")
    return(emptyDMRList)
  }
  
  #save.image("ChAMP_DMR_Test.Rdata")
  
  myResults.sig.cap <- myResults.sig[probeKeepers,] 
  myResults.sig.cap <- myResults.sig.cap[order(myResults.sig.cap$CHR, myResults.sig.cap$MAPINFO),] # orders object by chromosome then position
  rm(lasso.gr, probe.gr, lasso.probe.countOverlap, myResults.sig, probeKeepers)
  
  #big code replacement before this point...
  
  ### Start and End Point of the Lasso
  lasso.coord <- data.frame("CHR" = myResults.sig.cap$CHR, "arm" = myResults.sig.cap$arm, "lasso.start"=myResults.sig.cap$MAPINFO - myResults.sig.cap$lasso.radius, "lasso.end"=myResults.sig.cap$MAPINFO + myResults.sig.cap$lasso.radius)
  lasso.coord <- lasso.coord[order(lasso.coord$CHR, lasso.coord$lasso.start),]
  lasso.seq <- split(lasso.coord, paste(lasso.coord$CHR, lasso.coord$arm)) # Lasso coordinates split by Chromosome and Arm
  
  rm(lasso.coord)
  lasso.bp <- vector("list", length(lasso.seq)) # genomic coordinates of every base within all lassos, by chromosome
  names(lasso.bp) <- names(lasso.seq)
  for (k in 1:length(lasso.bp)) # k = CHR
  {
    dd <- vector("list", nrow(lasso.seq[[k]]))
    for (i in 1:length(dd)) # i = lasso bounds
    {
      dd[[i]] <- seq(lasso.seq[[k]][i, 3], lasso.seq[[k]][i, 4])
    }
    lasso.bp[[k]] <- sort(do.call(c, dd))
  }
  rm(lasso.seq)
  
  lasso.bp <- lapply(lasso.bp, unique) # leaves unique genome coordinates. NB names are factors-based
  diffs <- unlist(lapply(lasso.bp, function(x) c(FALSE, diff(x) <= minDmrSep)))
  chr <- sapply(strsplit(names(lasso.bp), " "), "[[", 1)
  chr.un <- unlist(chr)
  chr.un.nu <- as.numeric(chr.un)
  rm(chr, chr.un)
  len1 <- lapply(lasso.bp, length)
  len.un <- unlist(len1)
  len.un.vec <- as.vector(len.un)
  rm(len1, len.un)
  chr <- rep(chr.un.nu, len.un.vec)
  rm(chr.un.nu)
  bp <- unlist(lasso.bp)
  dmr.start.index <- which(diffs == FALSE)
  dmr.end.index <- c(dmr.start.index[-1] - 1, length(bp))
  dmrs <- data.frame("CHR" = chr[dmr.start.index],
                     "dmr.start" = bp[dmr.start.index],
                     "dmr.end" = bp[dmr.end.index])
  dmrs <- dmrs[order(dmrs$CHR, dmrs$dmr.start), ]
  rm(chr, diffs, bp, dmr.start.index, dmr.end.index)
  dmrs$dmr.size <- (dmrs$dmr.end - dmrs$dmr.start) + 1
  dmrs <- dmrs[dmrs$dmr.size >= minDmrSize,]
  rownames(dmrs) <- 1:nrow(dmrs) # renames rows after removing small DMRs
  dmrs$dmr.no <- 1:nrow(dmrs) # renames DMRs after removing small DMRs
  
  dmrs.gr <- GRanges(seqnames=paste("chr", dmrs$CHR, sep = ""), ranges=IRanges(start=dmrs$dmr.start, end=dmrs$dmr.end))
  probe.gr <- GRanges(seqnames=paste("chr", myResults$CHR, sep = ""), ranges=IRanges(start=myResults$MAPINFO, end=myResults$MAPINFO))
  dmr.probes.gr <- findOverlaps(dmrs.gr, probe.gr)
  dmr.probes <- as.data.frame(dmr.probes.gr)
  rm(dmrs.gr, probe.gr, dmr.probes.gr)
  
  toMatch <- c("dmr.no", "dmr.start", "dmr.end", "dmr.size")
  dmr.col.keep <- match(toMatch, colnames(dmrs))
  myDf <- data.frame(myResults[dmr.probes$subjectHits, ], dmrs[match(dmr.probes$queryHits, rownames(dmrs)), dmr.col.keep ])
  core.start <- lapply(split(myDf$MAPINFO, myDf$dmr.no), min)
  core.end <- lapply(split(myDf$MAPINFO, myDf$dmr.no), max)
  core.start.un <- unlist(core.start)
  core.end.un <- unlist(core.end)
  len <- as.vector(table(as.factor(myDf$dmr.no)))
  myDf$dmr.core.start <- rep(core.start.un, len)
  myDf$dmr.core.end <- rep(core.end.un, len)
  myDf$dmr.core.size <- (myDf$dmr.core.end - myDf$dmr.core.start) + 1
  rm(core.start, core.start.un, core.end, core.end.un, len)
  
  #don't think this works - this is from the lasso guys
  toMatch <- paste(c("p.", "q."), collapse = "|")
  rownames(myDf) <- sapply(strsplit(rownames(myDf), toMatch), "[[", 2)
  rm(dmrs, toMatch)
  
  myDf$CHR <- as.factor(as.character(myDf$CHR))
  dmr.beta.means=resultsFile[match(rownames(myDf),rownames(resultsFile)),]
  myDf <- data.frame(myDf, dmr.beta.means[,17:19,])
  rm(dmr.beta.means)
  
  dmr.beta <- split(as.data.frame(beta.norm[match(rownames(myDf), rownames(beta.norm)),]), myDf$dmr.no)
  corel <- lapply(dmr.beta, function(x) cor(t(x)))
  weights <- lapply(corel, function(x) 1/apply(x^2,1,sum))
  dmr.ind.p <- split(myDf$adj.P.Val, myDf$dmr.no)
  dmr.qp <- lapply(dmr.ind.p, qnorm)
  dmr.qp.w <- mapply("*", dmr.qp, weights)
  
  if(class(dmr.qp.w) == "matrix")
  {
    dmr.stat <- sum(dmr.qp.w)
  }else
  {
    dmr.stat <- lapply(dmr.qp.w, sum)
  }
  
  dmr.sd <- lapply(weights, function(x) sqrt(sum(x^2)))
  dmr.p <- mapply(function(x,y) pnorm(x,0, sd=y), dmr.stat, dmr.sd)
  dmr.rep <- as.numeric(summary(dmr.ind.p)[,1])
  
  #changed this
  dmr.p <- p.adjust(dmr.p, method = "BH")
  dmr.p.rank <- rank(dmr.p, ties.method="min")
  
  myDf$dmr.p <- rep(dmr.p, dmr.rep)
  myDf$dmr.p.rank <- rep(dmr.p.rank, dmr.rep)
  colKeep1 <- grep("AVG", colnames(myDf))
  dmr.means <- do.call(rbind, lapply(split(myDf[, colKeep1], myDf$dmr.no), colMeans))
  deltaBeta=dmr.means[,2]-dmr.means[,1]
  toMatch<-c("CHR","arm")
  if("gene" %in% colnames(myDf))
  {
    toMatch[length(toMatch)+1]="gene"
    
  }else{
    toMatch[length(toMatch)+1]="gene.1"
  }
  toMatch[4:12]<-c("dmr.no", "dmr.start", "dmr.end", "dmr.size", "dmr.core.start", "dmr.core.end", "dmr.core.size","dmr.p", "dmr.p.rank")
  colKeep2 <- c(match(toMatch, colnames(myDf)))
  myDfTrunc <- do.call(rbind, lapply(split(myDf[, colKeep2], myDf$dmr.no), function(x) x[1, ]))
  index <- match(c("dmr.p", "dmr.p.rank"), colnames(myDfTrunc))
  myDfTrunc <- data.frame(myDfTrunc[, 1:(index[1]-1)], dmr.means, deltaBeta, myDfTrunc[, index])
  rm(dmr.beta, corel, weights, dmr.ind.p, dmr.qp, dmr.qp.w, dmr.stat, dmr.sd, dmr.p, dmr.rep, dmr.p.rank, colKeep1, dmr.means, toMatch, colKeep2, index)
  
  message("You have found ",max(myDf$dmr.no)," DMRs.")
  
  dmrs <- myDfTrunc[order(myDfTrunc$dmr.no),]
  
  sig.dmrs=myDfTrunc[which(myDfTrunc$dmr.p < DMRpval),]
  message("You have found ",length(unique(myDfTrunc$dmr.no))," significant DMRs with a dmr.p < ",DMRpval,".")
  myDf<-myDf[order(myDf$dmr.p.rank), ]
  
  return("dmr.probes" = myDf)
}