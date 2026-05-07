#source('/Users/tomparee/Documents/Documents - MacBook Pro de tom/basics_TP.R')
library(readr)
library(data.table)
relative = function(x){(x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))}

#rils = fread("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/genomic/haploResolvedGeno/chrI_beceiPanels_geno_RILs_haplotypeResolved_v2.csv.gz")
#colnames(rils)


setwd("/Users/Sol/Documents/simu_panel_derivation/")

slimfile = "scripts/becei_panels_deriv_neutral.slim"
nsim = 50
nloci = 500

param = expand.grid(panel = c('A','B'),N = c(400,4000), nG_outcrossing = c(4,5,6,7))
param = param[-15,]

for(p in 1:nrow(param)){
  panel = param$panel[p]
  N = param$N[p]
  nG_outcrossing = param$nG_outcrossing[p]
  
  #panel = 'A'
  #N=4000
  #nG_outcrossing = 7
  
  nlines = ifelse(panel=='A',158,139)
  outdir = paste0("output/panel",panel,'_N',N,'_nGoutcross',nG_outcrossing,'_nloci',nloci)
  
  
  ###############################################
  ########### Run the simulations ###############
  
  # Create the directories to output simulations results for each chromosome
  if(!file.exists(outdir)){
    system(command = paste0("mkdir ", outdir))
    for(CHR in c("I",'II','III','IV','V','X')){
      system(command = paste0("mkdir ", outdir, "/", CHR))
    }
  }
  
  
  gg = lapply(1:nsim, function(i){
    
    print(i)
    
    for(CHR in c("I",'II','III','IV','V','X')){
      
      slimcmd =paste("slim",
                     paste0("-d nlines=",nlines),
                     paste0("-d nloci=",nloci),
                     paste0("-d N=",N),
                     paste0("-d nG_outcross=",nG_outcrossing),
                     paste0("-d chromType=\"'", ifelse(CHR=='X','X','A'), "'\""),
                     paste0("-d outpath=\"'", paste0(outdir,"/",CHR), "'\""),
                     slimfile,
                     " ")
      
      
      slim_out = system( slimcmd, intern=T)
      
    }
    
    return(NULL)
  })
  
  
}








###################################################
########### Analyse the simulations ###############


extractOutputMS=function(filename, nloci, outputpath){
  
  nloci = nloci*2
  
  outputMS <- fread(paste0(outputpath, "/", filename))
  
  relpos = as.character(outputMS[2,1])
  relpos = as.numeric(strsplit(relpos," ")[[1]][-1])
  relpos = relpos*(nloci-1) + 1
  relpos = as.integer(round(relpos, digits = 0))
  if(sum(duplicated(relpos))){stop("duplicated position")}
  
  geno = outputMS[-1:-2,1]
  
  length(as.numeric(unlist(strsplit(as.character(geno[1]),""))))
  
  geno = sapply(as.vector(geno), function(x){as.numeric(unlist(strsplit(as.character(x),"")))})
  geno = t(matrix(geno, ncol = length(relpos), byrow = T))
  
  
  #mx = match(round(1:nloci, digits = 0), round(relpos, digits = 0))   
  mx = match(1:nloci, relpos)   
  geno = geno[mx,]
  geno[is.na(geno)]=0
  
  return(geno)
}

getStats = function(thisfile,nloci,outputpath){
  
  # binary code for haplotype founders' haploid genome 1,2,3,4
  fmarkers = c(f1 = '11', # female founder g1
               f2 = '01', # female founder g2
               f3 = '10', # male founder g1
               f4 = '00') # male founder g2
  
  gt = extractOutputMS(filename=thisfile, nloci=nloci, outputpath=outputpath)
  lines = paste0('L',rep(1:(ncol(gt)/4),each = 4))
  
  gt = do.call(rbind, lapply(1:nloci, function(l){
    l = l*2 + c(-1,0)
    x = paste0(gt[l,][1,], gt[l,][2,])
    match(x,fmarkers)
  }))
  
  
  missingbins = which(is.na(match(levels(snps$cMbin),unique(snps$cMbin)))) # no snps there
  gt = gt[-missingbins,]
  
  linestats = lapply(unique(lines), function(thisline){
    #print(thisline)
    x = gt[,lines == thisline]
    x = apply(x,1,function(y){paste(unique(y), collapse = ";")})
    x = data.frame(cMbin=unique(snps$cMbin), founder = x)
    x$nfounder = (nchar(x$founder)+1)/2
    
    if(sum(x$nfounder > 1) == 0){nbinhet = 0}else{
      xhet = subset(x, nfounder > 1)
      nbinhet = sum(unlist(lapply(1:nrow(xhet), function(i){
        thesefounders = as.numeric(strsplit(xhet[i,]$founder,";")[[1]])
        hapgt = fgt[snps$cMbin == xhet[i,]$cMbin, thesefounders]
        
        ndiff = sum(apply(hapgt,1,function(y){length(unique(y))}) > 1)
        ndiff > 0
      })))
    }
    
    wherebreaks =  c(1,cumsum(rle(x$founder)$lengths))
    wherebreaks =  sort(unique(wherebreaks))
    
    breakpoints_position = do.call(rbind,lapply(wherebreaks[c(-1,-length(wherebreaks))], function(b){
      
      #print(b)
      weight = 1/max(x$nfounder[c(b,b+1)])
      fhap1 = as.numeric(strsplit(x$founder[b], ";")[[1]])
      fhap2 = as.numeric(strsplit(x$founder[b+1], ";")[[1]])
      hap1_lim = c(max(wherebreaks[wherebreaks<b]),b)
      hap2_lim = c(b+1, min(wherebreaks[wherebreaks>b]))
      
      whichcMbins1 = x$cMbin[hap1_lim[1]:hap1_lim[2]]
      whichcMbins2 = x$cMbin[hap2_lim[1]:hap2_lim[2]]
      
      
      #############################################################################################
      ###### let's find the divergent snps between hap1 and hap2 the closest to the breakpoint #### 
      ##### i.e., where we can observe the breakpoint
      
      
      # the fwo founder haplotype at the limit of hap1
      if(weight == 1){
        hapgt1 = fgt[snps$cMbin %in% whichcMbins1, c(fhap1,fhap2)]
      }else{
        hapgt1 = cbind(apply(fgt[snps$cMbin %in% whichcMbins1, fhap1, drop = F],1,mean),
                       apply(fgt[snps$cMbin %in% whichcMbins1, fhap2, drop = F],1,mean))
      }
      
      hapdiff1 = hapgt1[,1] != hapgt1[,2] # wether the two founder share a snp at each position
      last_different_snp1 = max(which(hapdiff1)) # the last snp that differ
      obspos1 = snps$POS[snps$cMbin %in% whichcMbins1][last_different_snp1] # the position of this last snp
      
      
      if(weight == 1){
        hapgt2 = fgt[snps$cMbin %in% whichcMbins2, c(fhap1,fhap2)]
      }else{
        hapgt2 = cbind(apply(fgt[snps$cMbin %in% whichcMbins2, fhap1, drop = F],1,mean),
                       apply(fgt[snps$cMbin %in% whichcMbins2, fhap2, drop = F],1,mean))
      }
      
      
      hapdiff2 = hapgt2[,1] != hapgt2[,2]
      first_different_snp2 = min(which(hapdiff2))
      obspos2 = snps$POS[snps$cMbin %in% whichcMbins2][first_different_snp2]
      
      realpos = max(snps$POS[snps$cMbin %in% whichcMbins1])
      
      if(is.na(obspos1) | is.na(obspos2)){obspos1=obspos2=NA}
      
      data.frame(realpos, obspos1, obspos2, line = thisline, weight = weight)
      
    }))

    list(nbinhet, breakpoints_position)
    
  })
  
  
  
  hetbins = unlist(lapply(linestats, function(x){x[[1]]}))
  recombinationbreaks = do.call(rbind, lapply(linestats, function(x){ x[[2]] }))
  
  list(hetbins=hetbins,recombinationbreaks=recombinationbreaks)
}








for(outdir in paste0('output/',list.files('output/'))){
  
  panel=tstrsplit(tstrsplit(outdir,'_')[[1]],'panel')[[2]]
  nloci = as.numeric(tstrsplit(tstrsplit(outdir,'_')[[4]],'nloci')[[2]])
  
  sim.results = list(I=NULL,II=NULL,III=NULL,IV=NULL,V=NULL,X=NULL)
  for(CHR in c("I",'II','III','IV','V','X')){
    
    snps = read_csv(paste0("/Users/Sol/Documents/becei/rils_genotable2026/",CHR,"_becei_variantInfo_founders&Rils.csv.gz"))
    fgt = read_csv(paste0("/Users/Sol/Documents/becei/rils_genotable2026/",CHR,"_becei_genotypes_founders.csv.gz"))
    fgt = as.data.frame(fgt)
    
    if(panel == 'A') fgt = fgt[,c("FA.g1","FA.g2","FM.g1","FM.g2")]
    if(panel == 'B') fgt = fgt[,c("FB.g1","FB.g2","FM.g1","FM.g2")]
    
    snps$cM[is.na(snps$cM) & snps$POS < mean(snps$POS, na.rm = T)] = min(snps$cM,na.rm=T)
    snps$cM[is.na(snps$cM) & snps$POS > mean(snps$POS, na.rm = T)] = max(snps$cM,na.rm=T)
    snps$cM = relative(snps$cM)
    
    cuthere = seq(0 ,1, 1/nloci)
    cuthere[which.min(cuthere)]  = -1e-10
    cuthere[which.max(cuthere)]  = 1 + 1e-10
    snps$cMbin = cut(snps$cM, cuthere)
    
    files = list.files(paste0(outdir,"/",CHR))
    
    output = parallel::mclapply(files, mc.cores = 2, function(thisfile){
      #print(thisfile)
      #thisfile = files[1]
      
      out=getStats(thisfile=thisfile,nloci=nloci,outputpath=paste0(outdir,"/",CHR))
      out$CHR = CHR
      out$nloci=nloci
      out$file = thisfile
      out
    })
    
    sim.results[[CHR]] = output
    
    
    
  }
  
  
  save(sim.results, file = paste0(outdir,"/","sim.results.Rdata"))
  
}




simstats = NULL
rmaps=NULL
for(outdir in paste0('output/',list.files('output/'))){
  
  load(paste0(outdir,"/","sim.results.Rdata"))
  print(outdir)
  
  
  nsim = min(unlist(lapply(sim.results, length)))
  
  if(nsim>0){
    res=do.call(rbind,lapply(1:nsim, function(i){
      h=do.call(rbind,lapply(sim.results, function(x){x[[i]]$hetbins}))
      nhetlines=sum(apply(h,2,function(x){sum(x>0)>0}))
      hetsize= apply(h,2,function(x){sum(x*(c(50,50,50,50,50,25)/nloci))})
      hetsize_all = mean(hetsize)
      hetsize_het = mean(hetsize[hetsize>0])
      
      co = sum(unlist(lapply(sim.results, function(x){
        b=x[[i]]$recombinationbreaks
        b=b[!is.na(b$obspos1),]
        
        mean(aggregate(obspos1~line, b, length)$obspos1)
        
      })))
      
      
      data.frame(outdir = outdir, nCO=co, nHet=nhetlines, sizeHetall=hetsize_all, sizeHet_het = hetsize_het)
      
    }))
    
    
    
    rr=do.call(rbind,lapply(1:nsim, function(i){
      
      do.call(rbind,lapply(sim.results, function(x){
        b=x[[i]]$recombinationbreaks
        b=b[!is.na(b$obspos1),]
        b$chrom = x[[i]]$CHR
        b$nsim = i
        b$outdir = outdir 
        b
      }))
      
    }))
    
    
    rmaps = rbind(rmaps,rr)
    
    
    simstats = rbind(simstats,res)
  }
  
  
  
}


simstats$nGoutcross = as.numeric(tstrsplit(tstrsplit(simstats$outdir,'_')[[3]],'nGoutcross')[[2]])
simstats$cross=tstrsplit(tstrsplit(simstats$outdir,'_')[[1]],'panel')[[2]]
simstats$N=as.numeric(tstrsplit(tstrsplit(simstats$outdir,'_')[[2]],'N')[[2]])

rmaps$nGoutcross = as.numeric(tstrsplit(tstrsplit(rmaps$outdir,'_')[[3]],'nGoutcross')[[2]])
rmaps$cross=tstrsplit(tstrsplit(rmaps$outdir,'_')[[1]],'panel')[[2]]
rmaps$N=as.numeric(tstrsplit(tstrsplit(rmaps$outdir,'_')[[2]],'N')[[2]])

save(simstats, rmaps, file = "./simu_deriv_results_COinterference.Rdata")


ggplot(simstats, aes(nGoutcross, nCO, color=as.factor(N)))+geom_point(alpha=0.5,shape=1)+stat_summary()+facet_wrap(~cross)
ggplot(simstats, aes(nGoutcross, nHet, color=as.factor(N)))+geom_point(alpha=0.5,shape=1)+stat_summary()+facet_wrap(~cross)
ggplot(simstats, aes(nGoutcross, sizeHet_het, color=as.factor(N)))+geom_point(alpha=0.5,shape=1)+stat_summary()+facet_wrap(~cross)
ggplot(simstats, aes(nGoutcross, sizeHetall, color=as.factor(N)))+geom_point(alpha=0.5,shape=1)+stat_summary()+facet_wrap(~cross)

View(simstats)



mean(unlist(lapply(split(recombinationbreaks,recombinationbreaks$line), nrow)))

r1 = recombinationbreaks[order(recombinationbreaks$realpos),]
r1$cumbreaks = 1:nrow(r1)


r2 = subset(recombinationbreaks, !is.na(obspos1))
r2 = r2[order(r2$obspos1+r2$obspos2),]
r2$cumbreaks = 1:nrow(r2)



ggplot()+
  geom_point(data=r1, aes(realpos, cumbreaks), color = 'grey')+ # real
  geom_point(data=r2, aes( (obspos1+obspos2)/2, cumbreaks)) # observed

