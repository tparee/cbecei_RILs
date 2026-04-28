setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/")

#######################
##### Parameters ######

min_nsnp = 20
min_psize_bp = 5e4
min_gsize_cM = 0.5


foundinghaplotype_homology = NULL
for(CHR in c("I","II", "III", 'IV', "V", "X")){
  print(CHR)
  phasedfounderhaplotypes = as.data.frame(fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/genotypes/",CHR,"_becei_genotypes_founders.csv.gz")))
  snps =  as.data.frame(fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/genotypes/",CHR,"_becei_variantInfo_founders&Rils.csv.gz")))
  snps$pos = snps$POS

  foundernames = c("FA.g1","FA.g2", "FM.g1", "FM.g2", "FB.g1","FB.g2")
  nfounders = length(foundernames)
  founderspairs = expand.grid(1:nfounders,1:nfounders)
  founderspairs = founderspairs[founderspairs[,1]<founderspairs[,2],]
  

  IBD_regions = do.call(rbind,lapply(1:nrow(founderspairs), function(i){
    print(i)
    f1 = foundernames[founderspairs[i,1]]
    f2 = foundernames[founderspairs[i,2]]
    #founders = phasedfounderhaplotypes[,foundernames]
    gt1 = phasedfounderhaplotypes[,f1]
    gt2 = phasedfounderhaplotypes[,f2]
    
    runs = rle(as.numeric(gt1 == gt2))
    sort(runs$lengths, decreasing = T)[1:10]
    
    
    group = unlist(lapply(1:length(runs$lengths), function(r){rep(r,runs$lengths[r])}))
    homology = unlist(lapply(1:length(runs$lengths), function(r){rep(runs$values[r],runs$lengths[r])}))
    
    ibd = cbind(whichsnp=1:length(gt1), group)[homology==1,]
    ibd = as.data.frame(ibd)
    

    ibd = do.call(rbind,lapply( split(ibd,ibd$group), function(x){
      keep = nrow(x) > min_nsnp &  diff(snps$POS[range(x$whichsnp)]) > min_psize_bp & diff(snps$cM[range(x$whichsnp)]) > min_gsize_cM
      if(keep){
        w1 = min(x$whichsnp); w2 = max(x$whichsnp)
        return(data.frame(
          nsnp = length(x$whichsnp),
          ppos1 = snps$POS[w1],
          ppos2 = snps$POS[w2],
          gpos1 =  snps$cM[w1],
          gpos2 =  snps$cM[w2]
        ))
      }else{return(NULL)}
    }))
    
    if(is.null(ibd)){return(NULL)}
    if(nrow(ibd)==0){return(NULL)}
    
    ibd$chrom = CHR
    ibd$founder1 = f1
    ibd$founder2 = f2
    return(ibd)
  }))
  
  foundinghaplotype_homology = rbind(foundinghaplotype_homology, IBD_regions)
}



#readr::write_csv(foundinghaplotype_homology,
#          "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/haplotypes/IBD_between_founding_haploid_genomes.csv")

read_csv("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/haplotypes/IBD_between_founding_haploid_genomes.csv")



library(circlize)
min_gsize_cM = 2
CHR = "V"
snps =  as.data.frame(fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/genotypes/",CHR,"_becei_variantInfo_founders&Rils.csv.gz")))
foundinghaplotype_homology = as.data.frame(fread("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/haplotypes/IBD_between_founding_haploid_genomes.csv"))
foundinghaplotype_homology = subset(foundinghaplotype_homology, chrom == CHR)    
foundinghaplotype_homology = foundinghaplotype_homology[foundinghaplotype_homology$gpos2 - foundinghaplotype_homology$gpos1 > min_gsize_cM,]

if(CHR == 'X'){
  foundercolors = c('#2A6218', '#9CC891','#0E1A52',"#C1443E", '#DFB640')
  foundernames = c('FA.g2', 'FA.g1','FM.g1','FB.g2', 'FB.g1')
  
  foundinghaplotype_homology = subset(foundinghaplotype_homology, !grepl("FM.g2",founder1) & !grepl("FM.g2",founder2))
  
}else{
  foundercolors = c('#2A6218', '#9CC891','#0E1A52', '#6DB9D2',"#C1443E", '#DFB640')
  foundernames = c('FA.g2', 'FA.g1','FM.g2', 'FM.g1','FB.g2', 'FB.g1')
}


#foundercolors = c('#2A6218', '#9CC891','#415384', '#6DB9D2','#9F4D3A', '#DFB640')
library(circlize)
sectorLim = data.frame(founder = rep(foundernames, each = 4),
                       foundercolors = rep(foundercolors, each = 4),
                       pos = rep(range(snps$POS), length(foundernames)*2), 
                       y = rep(c(0,0,0.3,0.3), each = length(foundernames)))


par(mar=rep(0,4))
circos.clear()
### Basic circos graphic parameters
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.03), start.degree = 90, gap.degree =2, track.height = 0.1)

circos.initialize(sectorLim$founder, x = sectorLim$pos/1e6)

circos.trackPlotRegion(sectorLim$founder, y = sectorLim$y,
                       panel.fun = function(x, y) {
                         thisfounder = CELL_META$sector.index
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                     col = foundercolors[foundernames==thisfounder], border=foundercolors[foundernames==thisfounder])
                         
                         circos.text(CELL_META$xcenter, 
                                     CELL_META$cell.ylim[2] + mm_y(5), 
                                     CELL_META$sector.index)
                         #circos.axis(labels.cex = 0.6)
                         
                       })

# grey for within founder homology (i.e., ROH)
# beige for between founder homology (i.e., IBD)
foundinghaplotype_homology$col_link = ifelse(tstrsplit(foundinghaplotype_homology$founder1, '.', fixed=T)[[1]]==tstrsplit(foundinghaplotype_homology$founder2, '.', fixed=T)[[1]], "#B2B2B280", "#DDD28380" )


for(i in 1:nrow(foundinghaplotype_homology)){
  
  circos.link(sector.index1=foundinghaplotype_homology$founder1[i],
              point1=c(foundinghaplotype_homology$ppos1[i], foundinghaplotype_homology$ppos2[i])/1e6,
              sector.index2=foundinghaplotype_homology$founder2[i],
              point2=c(foundinghaplotype_homology$ppos1[i], foundinghaplotype_homology$ppos2[i])/1e6,
              col = foundinghaplotype_homology$col_link[i])
  
}

text(0, 0, CHR, cex = 2.5, font = 2, col = "black")

