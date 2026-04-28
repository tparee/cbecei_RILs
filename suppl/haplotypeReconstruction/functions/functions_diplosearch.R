if(!require("limSolve")){install.packages("limSolve")}
if(!require("dplyr")){install.packages("dplyr")}
if(!require("data.table")){install.packages("data.table")}
if(!require("readr")){install.packages("readr")}

require(limSolve)
require(dplyr)
library(data.table)
library(readr)

# Function returning a list of window
# Input:
# size: the total size (all window cumul)
# winsize: the window size
# minsize: the minimum window size, i.e. the last window is smaller than winsize when size is not a multiple of winsize
#          => In this case, if the last window is smaller, fuse the two last windows
# Example:
# input: size = 45; winsize=10, minsize = 10
# output: list(1:10, 11:20, 21:30, 31:45)
get.win = function(size, winsize, minsize){
  
  if(winsize < size){
    
    nblocks = size/winsize
    fullblocks = floor(nblocks)
    
    ## split column numbers into 'nblocks' groups
    SPLIT <- split(1:(fullblocks*winsize), rep(1:fullblocks, each = winsize))
    if(nblocks>fullblocks){ 
      dblock = (nblocks-fullblocks)*winsize
      SPLIT[[length(SPLIT)+1]] = 1:dblock + max(fullblocks*winsize)
    }
    
    if(length(SPLIT[[length(SPLIT)]]) < minsize & length(SPLIT) > 1){
      SPLIT[[length(SPLIT)-1]] = c(SPLIT[[length(SPLIT)-1]],SPLIT[[length(SPLIT)]])
      SPLIT = SPLIT[1:(length(SPLIT)-1)]
    }
    
    
  }else{
    
    SPLIT = list(1:size)
  }
  
  return(SPLIT)
}


foundermatching = function(gx,predictors, hetmode = F){
  # gx = vector with single ril genotype with 0/0.5/1 for ref hom / het / alt hom
  # predictors in the founders genotype matrix for the same window
  # The function return a vector of length = ncol(predictors) with founder proportion
  # Example:
  # c(0,0,0,1) => 100% match with the 4th founder
  # c(0,0.2,0.8,0) => likely mean that there is a part of this genomic window corresponding to founder 1 and another to founder 2
  # c(0,0.5,0.5,0) can also mean that founder 2 and 3 are both an exact match (founders are identicals)
  
  # First, SNP with at least missing values
  predictNotMissing = apply(cbind(predictors, gx),1,function(x) sum(is.na(x))==0) 
  if(sum(predictNotMissing)<2){return(rep(NA, ncol(predictors)))}
  predictors = as.matrix(predictors[predictNotMissing,])
  gx = gx[predictNotMissing]
  
  #If a founder is heterozygous at some loci, consider it is the genotype of the ril at this loci
  # at SNP 10, ril is 0 and founder A is 0.5, change founder A SNP 10 in 0
  # It assumes that the founder is truely heterozygous (no genotyping error) 
  # And that whatever is the ril genotype at the herozygous snp should thus be a match
  # Or that we just don't know and still assume a match
  # You may want to question this choice because heterozygous allele is supposedly homozygous founder can be genotyping error
  # So you may just want change to founders[founders==0.5]=NA
  # If you do that a SNP which is heterozygous is one founder will be also considered as SNP with one missing value a be eliminated (above)
  # So you would lose information
  # Anyway, if the typical haplotype blocks contains a lot of different snp between founders it should not matters much
  # (If some founder are only different from a couple of SNP, it likely matters)
  
  if(!hetmode){
    heterozygous = which(predictors==0.5)
    if(length(heterozygous)>0){
      nsnp = nrow(predictors)
      rowhet=sapply(heterozygous, function(h){
        #wcol = floor(h/nsnp)+1
        wrow = (h/nsnp - floor(h/nsnp))*nsnp
        if(wrow == 0) wrow = nsnp#; wcol = wcol-1
        as.integer(round(wrow,digits = 1))
      })
      
      predictors[heterozygous]=gx[rowhet]
    }
    
    predictors = round(predictors)
    predictors[predictors>1]=1 
  }
  
  # First look if there is a perfect match with a founder
  # If yes, the function return c(0,1,0,0) (perfect match with founder 2)
  # Alternatively, it would return something like  c(0.333, 0.333, 0.333, 0) if there is perfect match with multiple founders
  matchgeno = sapply(1:ncol(predictors), function(founder){ sum(predictors[,founder] == gx)/nrow(predictors)})
  matchgeno =  round(matchgeno,digits = 7) == 1
  if(sum(matchgeno)>0){out=ifelse(matchgeno,1/sum(matchgeno), 0); return(out)}
  
  #If not, use lsei to solves a least squares problem under both equality and inequality constraints
  # Basically, we get a vector of founder frequency that would explain well SNP frequencies in the ril
  d = ncol(predictors)
  A = predictors
  E = t(matrix(rep(1,d)))
  Ff = 1
  G = diag(rep(1,d))
  H = matrix(rep(0.000,d))
  Y = gx
  
  out = lsei(A=A,B=Y,E=E,F=Ff,G=G,H=H,verbose=F)
  
  if(out$IsError){out = rep(NA,d)}else{out = out$X}
  return(out)
  
}




diplotype_foundermatching = function(win, ril, founders, phetdiff_threshold=0.05, phomdfiff_threshold = 0.01){
  gx = as.numeric(ril[win])
  predictors = founders[win,] # predictors = founders in the window
  foundergroups = cutree(hclust(dist(t(predictors))), h = 0)
  if(max(foundergroups)==1){return(NULL)}
  predictors = predictors[,!duplicated(foundergroups)]
  p = foundermatching(gx=gx,predictors = predictors, hetmode=T)
  p=as.numeric(p)
  
  bestfounder = names(foundergroups)[foundergroups == which.max(p)][1]
  maxmatch = 1-sum(calculateHamming(win, founderid1=bestfounder, founderid2=bestfounder, ril=ril, founders=founders))/(2*length(win))
  
  # perfect match with a single founder = homozygous diplotype without genotype error
  isperfectmatch = maxmatch == 1
  # If missmatch:
  # small missmatch: possibly genotyping error or recombination breakpont at the hedge of the window, or mix of two very similar haplotypes
  ishighmatch = maxmatch >= (1 - phomdfiff_threshold) & maxmatch < 1
  # Likely recombination breakpont or heterozygous diplotypes
  ismultiplefounder = sum(p>0.05) > 1
  
  if(!isperfectmatch & !ishighmatch & !ismultiplefounder){return(NULL)}
  if(ishighmatch & ismultiplefounder){ishighmatch=F}
  
  if(isperfectmatch){
    # There is a perfect match
    p= p[foundergroups]
    names(p)=colnames(founders)
    
    pdiplo = rbind(p,p)
    het = F
    new_win = win
    missmatches=NULL
    hamming = 0
  }
  
  
  if(ishighmatch){
    het = F
    missmatches = abs(predictors[,which.max(p)] - gx)
    hamming = sum(missmatches,na.rm=T)*2
    missmatches = win[which(missmatches>0)]
    
    p=as.numeric(p == max(p))
    p= p[foundergroups]
    names(p)=colnames(founders)
    pdiplo = rbind(p,p)
  }
  
  
  if(ismultiplefounder){
    # strong support for mix of founders, heterozygous or recombination?
    twobestmatches = order(p, decreasing = T)[1:2]
    diplogt = cbind(predictors[,twobestmatches], (predictors[,twobestmatches[1]]+predictors[,twobestmatches[2]])/2)
    hammings = as.numeric(apply(diplogt,2, function(x){ sum(abs(gx-x), na.rm=T)}))*2
    
    if(which.min(hammings) == 3 & hammings[3] < ceiling(phetdiff_threshold*length(win)*2)){
      het = T
      missmatches = abs(gx - (predictors[,twobestmatches[1]]+predictors[,twobestmatches[2]])/2)
      hamming =  sum(missmatches,na.rm=T)*2
      #print(hamming)
      missmatches = win[which(missmatches > 0)]
      pdiplo = matrix(0,ncol=ncol(predictors),nrow=2)
      pdiplo[cbind(1:2,twobestmatches)]=1
      pdiplo= pdiplo[,foundergroups]
      colnames(pdiplo)=colnames(founders)
    }else{return(NULL)}
  }
  
  out = list(win=win, ishet=het, pdiplo=pdiplo, hamming = hamming, missmatches=missmatches)
  
  return(out)
}



multiwin_foundersmatching = function(ril,founders){
  
  # 1) First we try to find regions with good match with some founders
  # We look within many different windows with varying size for correspondence with founders
  # This is done with runlsei out with return a vector of weight of lenght = ncol(founders)
  # 0,1,0,0 Perfect match with founder 2
  # several numbers > 0 & < 1, several founders possible or combination of founders
  # see runlsei for details
  winsizes = ceiling(length(ril)*c(0.5, 0.35, 0.26, 0.2, 0.17, 0.15, 0.12, 0.11, 0.09, 0.07, 0.049, 0.031, 0.018, 0.01, 0.001))
  winsizes = winsizes[winsizes >= 50]
  windows = do.call(c, lapply(winsizes, function(wsize){get.win(nrow(founders), winsize=wsize, minsize=50)}))
  names(windows)=NULL
  
  rangewindows = do.call(rbind, lapply(windows, range))
  windows = windows[order(rangewindows[,1], rangewindows[,2])]
  
  # For each window size:
  OUTPUT = lapply(1:length(windows), function(i){
    #print(i)
    win = windows[[i]]
    diplotype_foundermatching(win,ril, founders)
  })
  
  OUTPUT = OUTPUT[-which(unlist(lapply(OUTPUT, length)) == 0)]
  
  return(OUTPUT)
}


pdiplo_to_founderid = function(pdiplo){
  paste( sort(unique(apply(pdiplo, 1, function(y){ paste( sort(colnames(pdiplo)[y==1]),collapse=";")}))), collapse = "&")
}

getdiplowin_index = function(diplowin){
  do.call(rbind, lapply(diplowin, function(x){
    data.frame(pos1=min(x$win),
               pos2 = max(x$win),
               ishet = x$ishet,
               hamming = x$hamming,
               founder = pdiplo_to_founderid(x$pdiplo))
  }))
}



supressEncased = function(DIPLO){
  
  # Search for inferred haplotype encased in another
  # i.e. hap1 go from 1 to 150 & hap2 go from 10 to 140
  # and choose the larger haplo
  
  index = getdiplowin_index(DIPLO)
  order_vector = order(index$pos1, index$pos2)
  index = index[order_vector,]
  DIPLO = DIPLO[order_vector]
  
  h = 1
  while(h <= nrow(index)){
    
    pos1=index[h,1];pos2=index[h,2]
    
    isEncased = index[,1] <= pos1 & index[,2] >= pos2
    if(sum(isEncased)>1){DIPLO=DIPLO[-h]; index = index[-h,]}else{h=h+1}
  }
  
  return(DIPLO)
}


iscommonfounder = function(pdiplo1,pdiplo2){
 sum(apply(pdiplo1+pdiplo2, 1, max) == 2) == 2 | sum(apply(pdiplo1+pdiplo2[2:1,], 1, max) == 2) == 2
}

fuseAdjacent = function(DIPLO, ril, founders, mustoverlap = T, phetdiff_threshold=0.05, phomdfiff_threshold = 0.01){
  #HAPLO=HAPLO2
  index = getdiplowin_index(DIPLO)
  order_vector = order(index$pos1, index$pos2)
  index = index[order_vector,]
  DIPLO = DIPLO[order_vector]
  
  i = 1
  while(i < nrow(index)){
    #print(i)
    j=i+1
    
    overlapwithnext = index[j,1] <= index[i,2]+1
    
    d1 = DIPLO[[i]]$pdiplo
    d2 = DIPLO[[j]]$pdiplo
    
    commonfounder = sum(apply(d1+d2, 1, max) == 2) == 2 | sum(apply(d1+d2[2:1,], 1, max) == 2) == 2
    
    if( (overlapwithnext & commonfounder) | (!mustoverlap & commonfounder) ){
      b1 = index[i,"pos1"]
      b2 = index[j,"pos2"]
      x = diplotype_foundermatching(win=b1:b2, ril, founders,phetdiff_threshold=phetdiff_threshold, phomdfiff_threshold = phomdfiff_threshold)
      if(!is.null(x)){
        DIPLO[[i]] = x
        index[i,]=data.frame(pos1=x$win[1],pos2=max(x$win), ishet=x$ishet, hamming=x$hamming,  founder=pdiplo_to_founderid(x$pdiplo))
        DIPLO = DIPLO[-j]
        index = index[-j,]
      }else{i = i+1}
    }else{i = i+1}
  }
  
  return(DIPLO)
}



calculateHamming = function(win, founderid1, founderid2, ril, founders){
  hamming = abs(0.5*(founders[win,founderid1] + founders[win,founderid2]) - ril[win])*2
  hamming[is.na(hamming)]=0
  hamming
}

findbreak = function(pos1,pos2,push, direction, founderid1, founderid2, ril, founders, max_excess_hamming = 4){
  if( !(direction %in% c(-1,1))){stop("direction mush be -1 for left or 1 for right")}
  if(direction == -1){extendedwin = pos2:(pos1-push)}
  if(direction == 1){extendedwin = pos1:(pos2+push)}
  #range(extendedwin)
  currentwin = c(rep(T, pos2-pos1+1), rep(F, push))
  hamming = calculateHamming(win=extendedwin, founderid1=founderid1, founderid2=founderid2, ril=ril, founders=founders)
  
  cumhamming = cumsum(hamming)
  baseslope = abs((cumhamming[min(which(currentwin))] - cumhamming[max(which(currentwin))])/(pos1-pos2))
  baseslope = rep(baseslope, length(extendedwin))
  baseslope[baseslope < 2e-4 & !currentwin] = 2e-4 # One wrong (biallelic) genotype every 10^4 snps
  expectedhamming = cumsum(baseslope)
  excesshamming = cumhamming - expectedhamming
  
  #ggplot()+
  #  geom_point(data=data.frame(x=extendedwin, y=cumhamming, cw=currentwin), aes(x,y, color = cw ))+
  #  geom_line(data=data.frame(x=extendedwin, y= expectedhamming), aes(x,y))+
   # coord_cartesian(ylim=c(0,10))
  
  breakpoint = max(c(which(excesshamming <= max_excess_hamming & !currentwin), max(which(currentwin))))
  
  if(hamming[breakpoint] > 0){breakpoint = breakpoint - 1}
  breakpoint = extendedwin[breakpoint]
 
  return(breakpoint)
}


HammingChangebyTrimming = function(hammingvector, ptrim = 0.01, direction=1){
  untrim_hamming_tot = sum(hammingvector, na.rm=T)
  if(direction == -1){trimwin = ceiling(length(hammingvector)*ptrim):length(hammingvector)}
  if(direction == 1){trimwin = 1:floor(length(hammingvector)*(1-ptrim)) }
  trim_hamming_tot = sum(hammingvector[trimwin], na.rm=T)
  data.frame(trimpos1 = min(trimwin), trimpos2 = max(trimwin),
             untrim_hamming_tot = untrim_hamming_tot, trim_hamming_tot=trim_hamming_tot,
             hammingFoldChangeReduction = untrim_hamming_tot/trim_hamming_tot)
}


find_diplotype_limits = function(diplowin, ril, founders, max_excess_hamming = 4){
  
  push = round(0.075*length(ril))
  pos1 = min(diplowin$win)
  pos2 = max(diplowin$win)
  
  founderids = lapply(1:2, function(i){colnames(diplowin$pdiplo)[which(diplowin$pdiplo[i,] > 0)]})
  foundercomb = expand.grid(founderids[[1]], founderids[[2]])
  
  newlim = do.call(rbind, lapply(1:nrow(foundercomb), function(n){
    #print(n)
    p1 = pos1; p2 = pos2
    for(dir in c(-1,1)){
    hamming= calculateHamming(win=p1:p2, founderid1=as.character(foundercomb[n,1]), founderid2=as.character(foundercomb[n,2]), ril=ril, founders=founders)
    if(sum(hamming) > 0){
       hct = HammingChangebyTrimming(hammingvector = hamming, direction = dir)
       if(hct$hammingFoldChangeReduction >= 2 & dir == 1){p2 = (p1:p2)[hct$trimpos2]}
       if(hct$hammingFoldChangeReduction >= 2 & dir == -1){p1 = (p1:p2)[hct$trimpos1]}
     }
    }
    
    
    extendedlim = NULL
    for(dir in c(-1,1)){
      
      if(dir == 1){P = min(push, (length(ril)-pos2))}else{P = min(push, (pos1-1))}
      if(P==0 & dir == 1){l=p2}
      if(P==0 & dir == -1){l=p1}
      if(P>0){
        l = findbreak(pos1=p1,pos2=p2,
                      push=P, direction=dir,
                      founderid1=as.character(foundercomb[n,1]), founderid2=as.character(foundercomb[n,2]), ril=ril,
                      founders=founders,
                      max_excess_hamming=max_excess_hamming)
      }
      extendedlim = c(extendedlim, l)
    }
    
    return(extendedlim)
  }))
  
  newlim = newlim[which.max(newlim[,2]-newlim[,1]),]
  return(newlim)
}


refine_breakpoints = function(DIPLO, ril, founders){
  
  if(length(DIPLO) == 1){return(DIPLO)}
  
  index = getdiplowin_index(DIPLO)
  order_vector = order(index$pos1, index$pos2)
  index = index[order_vector,]
  DIPLO = DIPLO[order_vector]
  
  NEWDIPLOLIM = index[,1:2]
  
  if(NEWDIPLOLIM[1,1] > 1){ 
    x = diplotype_foundermatching(win=1:NEWDIPLOLIM[1,2], ril=ril, founders=founders)
    if(!is.null(x)){NEWDIPLOLIM[1,1]=1}
  }
  
  if(NEWDIPLOLIM[nrow(NEWDIPLOLIM),2] < length(ril)){ 
    x= diplotype_foundermatching(win=NEWDIPLOLIM[nrow(NEWDIPLOLIM),2]:length(ril), ril=ril, founders=founders)
    if(!is.null(x)){NEWDIPLOLIM[nrow(NEWDIPLOLIM),2] = length(ril)}
  }
  
  for(i in 1:(length(DIPLO)-1)){
    j = i+1
    
    pos1_i = NEWDIPLOLIM[i,1]
    pos2_j = NEWDIPLOLIM[j,2]
    WIN = pos1_i:pos2_j
    
    cum_hammings_distances = lapply(c(i,j), function(n){
      fnames = colnames(DIPLO[[n]]$pdiplo)
      founderids = lapply(1:2, function(g){fnames[which(DIPLO[[n]]$pdiplo[g,] > 0)]})
      foundercomb = as.data.frame(expand.grid(founderids[[1]], founderids[[2]]))
      
      out = lapply(1:nrow(foundercomb), function(m){
        hamming =  calculateHamming(win=WIN,
                                    founderid1=as.character(foundercomb[m,1]),
                                    founderid2=as.character(foundercomb[m,2]),
                                    ril=ril, founders=founders)
        
        if(n == i){cumham = cumsum(hamming)}
        if(n==j){cumham = rev(cumsum(rev(hamming)))}
        return(cumham)
      })
      return(out)
      
    })
    
    
    combihamm = expand.grid(1:length(cum_hammings_distances[[1]]), 1:length(cum_hammings_distances[[2]]))
    
    
    #ggplot()+theme_classic()+
    #geom_point(data=data.frame(x=WIN, y=cum_hammings_distances[[1]][[1]]), aes(x,y), color = 'salmon', shape=1, alpha=0.4)+
     #geom_point(data=data.frame(x=WIN, y=cum_hammings_distances[[2]][[1]]), aes(x,y), color='skyblue', shape=1, alpha=0.4)+
     #geom_line(data=data.frame(x=WIN, y=totham), aes(x,y), color='yellow',size=1)#+
     #coord_cartesian(ylim = c(0,5000), xlim = c(5000,15000))+
     #geom_vline(xintercept = c(7188,9262))
      
    
   #0.03988915
    breakpoint = do.call(rbind, lapply(1:nrow(combihamm), function(n){
      totham = cum_hammings_distances[[1]][[combihamm[n,1]]] + cum_hammings_distances[[2]][[combihamm[n,2]]]
      mintotham = min(totham)
      #mintotham/mean(totham)
      out = matrix(c(range(WIN[which(totham == mintotham)]),mintotham), nrow=1)
      colnames(out) = c("pos1","pos2","total_hammingdistance")
      as.data.frame(out)
    }))
    
    breakpoint = breakpoint[which.min(breakpoint$total_hammingdistance),]
    
    d1 = diplotype_foundermatching(win=NEWDIPLOLIM[i,1]:(breakpoint$pos2-1), ril=ril, founders=founders, phetdiff_threshold=0.08, phomdfiff_threshold = 0.02)
    d2 = diplotype_foundermatching(win=breakpoint$pos1:NEWDIPLOLIM[j,2], ril=ril, founders=founders, phetdiff_threshold=0.08, phomdfiff_threshold = 0.02)
    if(!is.null(d1) & !is.null(d2)){
      NEWDIPLOLIM[i,2] = (breakpoint$pos2-1)
      NEWDIPLOLIM[j,1] = breakpoint$pos1
    }
    
  }
  
  
  NEWDIPLO = lapply(1:nrow(NEWDIPLOLIM), function(n){
    diplotype_foundermatching(win=NEWDIPLOLIM[n,1]:NEWDIPLOLIM[n,2], ril=ril, founders=founders, phetdiff_threshold=0.08, phomdfiff_threshold = 0.02)
  })
  
  #failed = which(unlist(lapply(NEWDIPLO, is.null)))
  #if(length(failed) > 0){NEWDIPLO[failed] = DIPLO[failed]}
  
  return(NEWDIPLO)
}

runs_test_z <- function(x) {
  # Wald–Wolfowitz runs test for binary vector (0/1)
  # Z < 0: clustering; Z > 0: dispersion; Z ~ N(0,1) under H0
  x=x[!is.na(x)]
  n <- length(x)
  runs <- 1 + sum(x[-1] != x[-n])  # count runs
  n1 <- sum(x)
  n0 <- n - n1
  
  if (n0 == 0 || n1 == 0){NA}
  
  expected_runs <- 1 + (2 * n0 * n1) / n
  var_runs <- (2 * n0 * n1 * (2 * n0 * n1 - n)) / (n^2 * (n - 1))
  
  Z <- (runs - expected_runs) / sqrt(var_runs)
  
  return(Z)
}



#thisril = "QG3234"
#if(subset(meta, rilname == thisril)$panel == "α") foundernames = c("FA.g1","FA.g2", "FM.g1", "FM.g2")
#if(subset(meta, rilname == thisril)$panel == "β") foundernames = c("FB.g1","FB.g2", "FM.g1", "FM.g2")
#diplosearch(ril = rils[,thisril], founders=phasedfounderhaplotypes[,foundernames])


diplosearch = function(ril, founders, mask_suspicious_snps = F){
  
  #if(maskhet){ril[ril == 0.5] = NA}
  
  # First we generate many windows
  # and test wether there is a founder that match well ( > 0.99)
  # Or a combination of founder in case of heterozygous 
  # The algorithm return an output if a good match is found
  # DIPLO is a list, each element is a list itself containing:
  # $win (a "which" vector corresonding to the wondo), $ishet, $pdiplo, $hamming, $missmatches
  # $ishet (T or F; was a heterozygous diplotype block inferred)
  # $pdiplo: a matrix of dims: 2 x ncol(founders) containing 0 or 1 value
  # 1 indicate a match with the corresponding founder; rows are the two genome copy and are identical in case of homozygous
  # $hamming is the total hamming distance between the observed genotype and the haplotype genotype at the window. +1 is one diverging allele (i.e., alt-ref = 2; het-ref=1)
  # missmatches indicate the SNPs that have hamming distance > 0
  DIPLO = multiwin_foundersmatching(ril,founders)

  # Here, we extract window info (including hamming distance)
  index = cbind(getdiplowin_index(DIPLO),
                # and compute Z, whose negative value indicate that there are many missmatch concentrated in one place (suspicious)
                Z= unlist(lapply(DIPLO, function(x){
                  if(length(x$missmatches) == 0){return(0)}
                  return(runs_test_z(as.numeric(x$win %in% x$missmatches)))
                })),
                S = unlist(lapply(DIPLO, function(x){
                  if(length(x$missmatches) == 0){return(0)}
                  return(cluster_mismatch_score(as.numeric(x$win %in% x$missmatches)))
                })))
  
  
  
  
  # We delete suspicious windows (high hamming and low Z value)
  # Except for heterozygous window where is it expected that we have more error (het can be observed at 0 or 1 by chance, especially at lower depth)
  DIPLO = DIPLO[which(!(index$hamming > 10 & (index$Z < -90 | index$S > 3) & !index$ishet))]
  
  if(mask_suspicious_snps){
    sus = sort(unique(unlist(lapply(DIPLO, function(x){x$missmatches}))))
    index = getdiplowin_index(DIPLO)
    index$wsize = index$pos2 - index$pos1
    sus = unlist(lapply(sus, function(i){
      w = (index$pos1 + 0.01*index$wsize) <= i & (index$pos2 - 0.01*index$wsize) >= i
      nwin = sum(w)
      if(nwin==0){return(NULL)}
      pmissmatch = sum(unlist(lapply(DIPLO[w],function(x){i %in% x$missmatches})))/nwin
      if(pmissmatch==1){return(i)}else{return(NULL)}
    }))
    ril[sus] = NA
  }
  
  
  
  
  #Here, we delete window where that are already contained in a bigger window
  # to avoid redundant information
  DIPLO = supressEncased(DIPLO)
  #sus = sort(unique(unlist(lapply(DIPLO, function(x){x$missmatches})))) # suspicious snps (snps that don't match a founder despite strong support of all other snps)
  #rilclean = ril; if(length(sus) > 0){rilclean[sus] = NA} # we ignore them just for the next step
  
  # Then, we fuse together overlapping windows were the inferred founders match
  DIPLO = fuseAdjacent(DIPLO, ril, founders)
  getdiplowin_index(DIPLO)
  # Here we extend the limit of the window
  # We do that by finding when hamming value start to increases when pushing windows limit
  # This step is not perfect and we will refine the diplotype breakpoint below
  DIPLO = lapply(1:length(DIPLO), function(i){
    #print(i)
    pass = F
    maxexcesshamming = 4
    while(!pass & maxexcesshamming > -4){
      newlim = find_diplotype_limits(diplowin = DIPLO[[i]],ril=ril, founders=founders, max_excess_hamming = maxexcesshamming)
      newdiplo = diplotype_foundermatching(win=newlim[1]:newlim[2], ril=ril, founders=founders)
      if(!is.null(newdiplo)){pass = T}
      maxexcesshamming = maxexcesshamming -2
    }
    if(!pass){newdiplo = DIPLO[[i]]}
    return(newdiplo)
  })
  
  # We fuse windows that need to be fused (case where extending their limit makes them overlap and their founder is compatible)
  DIPLO = supressEncased(DIPLO)
  DIPLO = fuseAdjacent(DIPLO, ril, founders)
  
  # Try to fill gaps
  index = getdiplowin_index(DIPLO)
  gaps =  data.frame(pos1 = c(1,index[2,nrow(index)]), pos2 = c(index[1,1],length(ril)))
  gaps = data.frame(pos1=index[1:(nrow(index)-1),2],pos2=index[2:nrow(index),1])
  gaps = gaps[ which(gaps[,2] - gaps[,1] >= 50),]
  if(nrow(gaps)>0){
    DIPLOGAP = do.call(c,lapply(1:nrow(gaps), function(i){
      gap = unlist(gaps[i,])
      win = gap[1]:gap[2]
      nsplit = 1
      wsize = round(length(win)/nsplit)
      pass = F
      while(!pass & wsize >= 50){
        
        subwindows = get.win(size=length(win),winsize = wsize, minsize = 50)
        
        d = lapply(subwindows, function(subwin){
          diplotype_foundermatching(win=win[subwin], ril=ril, founders=founders)
        })
        
        names(d)=NULL
        d = d[lapply(d,length) > 0]
        if(length(d) > 0){pass = T}else{ nsplit = nsplit*3}
        wsize = round(length(win)/nsplit)
      }
      
      if(pass){
        return(d)
      }else{
        return(NULL)
      }
    }))
    
  }else{DIPLOGAP=NULL}
  
  
  if(length(DIPLOGAP) > 0){
    DIPLO = c(DIPLO, DIPLOGAP);
    index = getdiplowin_index(DIPLO)
    DIPLO = DIPLO[order(index$pos1, index$pos2)]
    index = index[order(index$pos1, index$pos2),]
    }
  
  DIPLO = supressEncased(DIPLO)
  DIPLO = fuseAdjacent(DIPLO, ril, founders, mustoverlap = F,phetdiff_threshold=0.08, phomdfiff_threshold = 0.02)
  
  # We refine the breakpoint by finding the position between two adjacent window that minimize the total hamming distance
  DIPLO = refine_breakpoints(DIPLO=DIPLO, ril=ril, founders=founders)
  DIPLO = supressEncased(DIPLO)
  
 
  ril[unique(unlist(lapply(DIPLO,function(x){x$missmatches})))] = NA
  DIPLO = fuseAdjacent(DIPLO, ril, founders, mustoverlap = F,phetdiff_threshold=0.08, phomdfiff_threshold = 0.02)
  #DIPLO2 = fuseAdjacent(DIPLO, ril, founders, mustoverlap = F,phetdiff_threshold=0.08, phomdfiff_threshold = 0.02)
  
  #getdiplowin_index(DIPLO)
  #getdiplowin_index(DIPLO2)
  #plot(cumsum(as.numeric(DIPLO[[3]]$win %in% DIPLO[[3]]$missmatches)))
  
  diplotable = getdiplowin_index(DIPLO)[,c("pos1","pos2","ishet","founder")]
  colnames(diplotable)[1:2] = c("whichsnp1","whichsnp2")
  return(diplotable)
}



# test = cumsum(as.numeric(DIPLO[[1]]$win %in% DIPLO[[1]]$missmatches))
# xx = 1:length(test)
# 
# plot(xx,test)
# summary(lm(test~xx))
# plot(cumsum(astestplot(cumsum(as.numeric(DIPLO[[3]]$win %in% DIPLO[[3]]$missmatches)))))
# x = DIPLO[[1]]
# x = as.numeric(x$win %in% x$missmatches)


cluster_mismatch_score <- function(x) {
  # variance-based
  r <- rle(x)
  var_run0 <- var(r$lengths[r$values == 0])
  var_run1 <- var(r$lengths[r$values == 1])
  
  
  n0 <- sum(x == 0)
  n1 <- sum(x == 1)
  n  <- n0 + n1
  
  p0 <- n0/ n
  p1 <- n1 / n
  
  # geometric variance
  var_geom_0 <- p0 / p1^2
  #var_geom_1 <- p1 / p0^2
  
  s <- log(var_run0 / var_geom_0)
  #s1 <- log(var_run1 / var_geom_1)
  
  s
}


# var(rle(sample(x))$lengths[r$values == 0])
# 
# 
# 
# score <- log(obs_var / exp_var)
# 
# 
# acf(x, lag.max = 50, plot = T)$acf[2]
# 
# # 
#x = diplotype_foundermatching(win=8401:28366, ril=ril2, founders=founders)$missmatches
# # diplotype_foundermatching(win=66179:78778, ril=ril, founders=founders)
# 
#(8401:28366)[is.na(ril[8401:28366])]
#sum(is.na(ril2[8401:28366]))
# 
# 
# View(cbind(founders[win,], ril[win]))
# 
# 
# 
#diplotype_foundermatching(win=50437:54454, ril=ril, founders=founders, phetdiff_threshold=0.08, phomdfiff_threshold = 0.03)

# 
# 
# 


