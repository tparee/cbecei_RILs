library(data.table)
library(ggplot2)
library(limSolve)
library(readr)
setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/")






runlsei = function(predictors, Y){
  require(limSolve)
  
  
  predictors = as.matrix(predictors)
  Y = as.numeric(Y)
  
  #ncols = ncol(predictors)
  #groups = cutree(hclust(dist(t(predictors))), h = 0)
  #if(max(groups)==1){return(rep(NA,ncols))}
  #predictors = predictors[,!duplicated(groups)]
  
  d = ncol(predictors)
  A = predictors
  E = t(matrix(rep(1,d)))
  Ff = 1
  G = diag(rep(1,d))
  H = matrix(rep(0.004,d))
  out = lsei(A=A,B=Y,E=E,F=Ff,G=G,H=H,verbose=F)
  
  if(out$IsError){out = rep(NA,d)}
  else{
    out = out$X
    out = out/sum(out)
  }
  
  return(out)
}


winsize = 2000
stepsize = 500


for(CHR in c('II','III','IV',"V", "X")){
  
  snps = as.data.frame( fread( paste0("genotypes/",CHR,"_becei_variantInfo_founders&Rils.csv.gz")) )
  FGT = as.data.frame( fread(paste0("genotypes/",CHR,"_becei_genotypes_founders.csv.gz")) )
  
  AD <- fread(
    paste0("suppl/private/alignment&variantCalling/VCF/RILs/unfiltered/",CHR,"_becei_variants_ad.csv.gz"),
    select = c("CHROM", "POS", "REF" , "ALT" , "POPA-1_AD", "POPA-2_AD", "POPA-3_AD","POPB-1_AD", "POPB-2_AD", "POPB-3_AD")
  )
  
  AD$ID = paste0(AD$CHROM, ":", AD$POS, ":", AD$REF, ":", AD$ALT)
  AD = AD[match(snps$ID, AD$ID),]
  
  AD = as.matrix(AD)[,grepl("POP", colnames(AD))]
  REF = do.call(cbind, lapply(1:ncol(AD), function(j){as.numeric(tstrsplit(AD[,j],",")[[1]])}))
  ALT = do.call(cbind, lapply(1:ncol(AD), function(j){as.numeric(tstrsplit(AD[,j],",")[[2]])}))
  colnames(REF) = colnames(ALT) = colnames(AD)
  
  REF = cbind(alpha=apply(REF[,1:3],1,sum),beta=apply(REF[,4:6],1,sum))
  ALT = cbind(alpha=apply(ALT[,1:3],1,sum),beta=apply(ALT[,4:6],1,sum))
  FREQ = ALT/(REF+ALT)
  
  
  NEWFREQ = matrix(NA, ncol=ncol(FREQ), nrow=nrow(FREQ))
  colnames(NEWFREQ) = colnames(FREQ)
  
  for(THISPANEL in c("alpha","beta")){
    
    if(THISPANEL == "alpha"){
      founders = c("FA.g1","FA.g2","FM.g1","FM.g2")
    }else{
      founders = c("FB.g1","FB.g2","FM.g1","FM.g2")
    }
    
    
    ix = seq(winsize/2, nrow(FGT)-stepsize, stepsize)
    
    foundersfreq = do.call(rbind, lapply(ix, function(i){
      #print(i)
      win = max(c(1,(i-(winsize/2)))):min(c(nrow(FGT),(i+(winsize/2))))
      p = runlsei(predictors = FGT[win,founders], Y = FREQ[win,THISPANEL])
      p = c(pos1 = min(win), pos2=max(win),p)
      return(p)
    }))
    
    foundersfreq = as.data.frame(foundersfreq)
    save(foundersfreq,file = paste0("analysis/temp/",CHR,"_foundersfreq_G5_",THISPANEL, ".Rdata"))
    
    midpos = (foundersfreq$pos1[2:nrow(foundersfreq)] +  foundersfreq$pos2[1:(nrow(foundersfreq)-1)])/2
  
    foundersfreq$pos1[2:nrow(foundersfreq)] = midpos +1
    foundersfreq$pos2[1:(nrow(foundersfreq)-1)] = midpos
    
    newfreq = unlist(lapply(1:nrow(foundersfreq), function(i){
      win = foundersfreq[i,"pos1"]:foundersfreq[i,"pos2"]
      correctedfreq = apply(t(as.matrix(FGT[win,founders])) * as.numeric(foundersfreq[i,founders]), 2, sum)
      return(correctedfreq)
    }))
    
    #plot(newfreq,FREQ[,THISPANEL])
    
    NEWFREQ[,THISPANEL] = newfreq
  }
  
  
  
  FREQ = cbind(FREQ, NEWFREQ)
  colnames(FREQ) = c("POPA-AF-obs", "POPB-AF-obs","POPA-AF-corrected", "POPB-AF-corrected")
  OUTPUT = cbind(snps, AD, FREQ)
 
  write_csv( as.data.frame(OUTPUT), paste0("suppl/populations_AF/",CHR,"_alleleFrequencies_G5populations.csv.gz"))
  
}




