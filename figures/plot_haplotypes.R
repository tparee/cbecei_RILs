library(readr)
library(ggplot2)
library(gridExtra)
library(data.table)

setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/")
meta = read_csv("suppl/RILs_sequencing_metadata.csv")
source("Utils.R")

# foundercolors = c(`FA.g1` = '#9CC891',
#                   `FA.g2` = '#2A6218',
#                   `FM.g1` = '#6DB9D2',
#                   `FM.g2` = '#415384',
#                   `FB.g1` = '#DFB640',
#                   `FB.g2` = "#C1443E")

######################################################################################
#########################  MAIN FIGURE WITH ALL RILs  ################################
######################################################################################

map = as.data.frame(read_csv("analysis/temp/cumbreaks_sampled.csv"))

coord = do.call(rbind, lapply(split(map, map$chrom), function(x){
  data.frame(chrom = x$chrom[1], pos=seq(0,max(x$pos), 4e6))
}))
########################################
#### Import haplotype blocks data  #####

toadd = 0
RIls_FounderHaplotypeBlocks=NULL
for(CHR in c("I","II", "III", 'IV', "V", "X")){
  
  x=read_csv(paste0("haplotypes/",CHR,"_rils_foundingHaplotypesBlocks.csv"))
  x = as.data.frame(x)
  x$chrom = CHR
  maxpos=max(x$pos2)
  x$pos1=x$pos1+toadd 
  x$pos2=x$pos2+toadd 
  map$pos[map$chrom == CHR] = map$pos[map$chrom == CHR] + toadd
  coord$cumpos[coord$chrom == CHR] = coord$pos[coord$chrom == CHR] + toadd
  toadd = toadd+maxpos + 5e5
  RIls_FounderHaplotypeBlocks = rbind(RIls_FounderHaplotypeBlocks,x)
}


RIls_FounderHaplotypeBlocks$founder[RIls_FounderHaplotypeBlocks$chrom == "X"] = gsub("FM.g1", "FM.g2", RIls_FounderHaplotypeBlocks$founder[RIls_FounderHaplotypeBlocks$chrom == "X"])
RIls_FounderHaplotypeBlocks$panel = meta$panel[match(RIls_FounderHaplotypeBlocks$rilname,meta$rilname)]



##############################################
#### Get the cumulative break positions  #####

pseudomarkers = seq(min(RIls_FounderHaplotypeBlocks$pos1), max(RIls_FounderHaplotypeBlocks$pos2), 20000)


foundersfreq = do.call(rbind, lapply(pseudomarkers, function(i){
  
  do.call(rbind, lapply(unique(RIls_FounderHaplotypeBlocks$panel), function(cross){
    x=subset(RIls_FounderHaplotypeBlocks, pos1 <= i & pos2 > i & panel == cross)
    if(nrow(x)==0){return(NULL)}
    chr = x$chr[1]
    p=table(x$founder)
    
    p = unlist(lapply(names(p), function(np){
      if(grepl(";|&",np)){
        np2 = unlist(strsplit(np, ';|&'))
        y = rep(p[np]/2, length(np2))
        names(y) = np2
        return(y)
      }else{
        return(p[np])
      }
    }))
    
    p=data.frame(founder=names(p), count = p)
    p = aggregate(count~founder, p, sum)
    p$freq = p$count/sum(p$count)
    p$pos = i
    p$cross = cross
    p$chrom = chr
    
    if(cross=='α'){
      forder = c("FM.g2", "FM.g1", "FA.g2", "FA.g1")
    }else{
      forder = c("FB.g1", "FB.g2", "FM.g1", "FM.g2")
    }
    p = p[order(match(p$founder,forder)),]
    p$cumfreq1 = c(0,cumsum(p$freq)[-nrow(p)])
    p$cumfreq2 = cumsum(p$freq)
    p
  }))
  
}))


ggplot(foundersfreq, aes(pos, ymin=cumfreq1, ymax=cumfreq2, fill=founder,color=founder, group=paste0(chrom, founder)))+
  theme_perso()+
  geom_ribbon()+facet_grid(cross~.)+
  scale_color_manual(breaks = names(foundercolors), values = foundercolors)+
  scale_fill_manual(breaks = names(foundercolors), values = foundercolors)+
  coord_cartesian(expand=0)






RIls_FounderHaplotypeBlocks=do.call(rbind,lapply(split(RIls_FounderHaplotypeBlocks, paste0(RIls_FounderHaplotypeBlocks$rilname, RIls_FounderHaplotypeBlocks$chrom)), function(HAPLO){
  
  print(HAPLO$rilname[1])
  
  if(nrow(HAPLO) == 1){return(HAPLO)}
  #HAPLO=split(Rilshaplotypes, Rilshaplotypes$rilname)[[1]]
  #HAPLO = subset(Rilshaplotypes,rilname=="A_QG3408")
  HAPLO = HAPLO[order(HAPLO$pos1, HAPLO$pos2),]
  
  midpos = (HAPLO$pos2[1:(nrow(HAPLO)-1)] + HAPLO$pos1[2:nrow(HAPLO)])/2
  HAPLO$pos2[1:(nrow(HAPLO)-1)] =  HAPLO$pos1[2:nrow(HAPLO)] = midpos
  #fm = do.call(rbind, getFoundersVector(founders=HAPLO$founders, foundernames=foundercolnames))
  #colnames(fm) = foundercolnames
  #HAPLO = cbind(HAPLO, fm)
  
  HAPLO
  
}))


haploPlotFromat = function(Rilshaplotypes, widthRILs = 0.95){
  
  haplop = do.call(rbind, lapply(split(Rilshaplotypes, Rilshaplotypes$rilname), function(haplotype){
   
    #haplotype = split(Rilshaplotypes, Rilshaplotypes$rilname)[[5]]
    #haplotype=divideOverlapHapBlocksByFounderID(haplotype, foundercolnames=foundercolnames,info=info)
    
   
    haplotype$nfounder = unlist(lapply(strsplit(haplotype$founder, ""), function(x){ sum(x %in% c(";","&")) + 1}))
    haplotype=do.call(rbind,lapply(1:nrow(haplotype), function(i){
      #print(i)
      x=haplotype[i,]
      if(x$nfounder>1){
        founders=strsplit(x$founder,";|&")[[1]]
        rely = seq(0,1,1/length(founders))
        rely1 = rely[1:(length(rely)-1)]
        rely2 = rely[2:length(rely)]
        x=x[rep(1,length(rely1),),]
        x$rely1=rely1
        x$rely2=rely2
        x$founder=founders
      }else{
        x$rely1=0;x$rely2=1
      }
      x
    }))
    
    haplotype
  }))
  
  
  haplop$y1 = as.numeric(factor(haplop$rilname))-1 + (haplop$rely1)*widthRILs
  haplop$y2 = as.numeric(factor(haplop$rilname))-1 + (haplop$rely2)*widthRILs
  haplop
  
}


meta = meta[order(meta$panel, meta$rilname),]
meta = subset(meta, rilname %in% RIls_FounderHaplotypeBlocks$rilname)
RIls_FounderHaplotypeBlocks$rilname = factor(RIls_FounderHaplotypeBlocks$rilname, levels = meta$rilname)
haplop = haploPlotFromat(Rilshaplotypes = RIls_FounderHaplotypeBlocks, widthRILs = 0.95)
haplop$haplocolor=foundercolors[match(haplop$founder, names(foundercolors))]
#haplop$haplocolor[is.na(haplop$haplocolor)]="black"



chrompos = aggregate(pos1~chrom, haplop, min)
chrompos$pos2 =  aggregate(pos2~chrom, haplop, max)$pos2
chrompos$midpos = (chrompos$pos1 + chrompos$pos2)/2

haplop$y1[haplop$panel == "β"] = haplop$y1[haplop$panel == "β"] + 3.5
haplop$y2[haplop$panel == "β"] = haplop$y2[haplop$panel == "β"] + 3.5


facettitle_pos = aggregate(y1 ~ panel, haplop, function(x){mean(range(x))})
facettitle_pos$x = max(haplop$pos2) + 1e5

pmosaic = ggplot(haplop)+
  theme_Publication3(base_size = 10)+
  geom_text(data = facettitle_pos, aes(x = 1.5 + x/1e6, y = y1, label = panel), fontface = "bold", size = 10/2.8)+
  geom_rect(aes(xmin=pos1/1e6,xmax=pos2/1e6, ymin=y1, ymax=y2,fill=founder), color = NA)+
  scale_fill_manual(values = foundercolors, breaks = names(foundercolors))+
  #theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0), breaks = chrompos$midpos/1e6, labels = chrompos$chrom,position = 'top')+
  scale_y_reverse(expand = c(0,0))+
  theme(axis.text.y =element_blank(),
        axis.line.y = element_line(color = "white"),
        axis.ticks.y = element_line(color = "white"),
        axis.ticks.x = element_line(color = "white"),
        axis.title = element_blank(),
        axis.line  = element_blank(),
        #axis.text.x = element_text(size = 15),
        panel.background = element_rect(fill = "white", colour=NA),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5,20,2,16), "pt"))+
  coord_cartesian(clip = "off", xlim = range( c(haplop$pos1, haplop$pos2 ))/1e6)+
  ylab("")+xlab("")


pfreq = ggplot(foundersfreq, aes(pos, ymin=cumfreq1, ymax=cumfreq2, fill=founder,color=founder, group=paste0(chrom, founder)))+
  theme_perso()+
  geom_ribbon()+facet_grid(cross~.)+
  scale_color_manual(breaks = names(foundercolors), values = foundercolors)+
  scale_fill_manual(breaks = names(foundercolors), values = foundercolors)+
  scale_y_continuous(breaks = c(0,1))+
  coord_cartesian(expand=0, xlim = c(min(haplop$pos1), max(haplop$pos2)))+
  #theme(legend.position = "none")+
  theme(axis.text.x =element_blank(),
        #axis.line.y = element_line(color = "white"),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 7, face = 'plain'),
        #axis.text.y =element_text(),
        panel.background = element_rect(fill = "white", colour=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, size = 10),
        plot.margin = unit(c(2,0,3,10), "pt"))+
  ylab('p haplo')



pmaps = ggplot()+
  theme_classic()+
  geom_point(data=map,aes(pos, cumbreaks),size=0.55/.pt)+
  coord_cartesian(expand = 0,ylim=c(0, max(map$cumbreaks)*1.11),xlim = c(min(haplop$pos1), max(haplop$pos2)))+
  geom_rect(data=chrompos, aes(xmin=pos1, xmax=pos2, ymin=0, ymax=max(map$cumbreaks)*1.1, group=chrom), fill=NA, size=0.5, color='black')+
  xlab('Physical position (Mb)')+
  ylab('cum. breaks')+
  scale_x_continuous(breaks = coord$cumpos, labels = coord$pos/1e6)+
  theme( axis.text = element_text(size=6),
         axis.title = element_text(size=7),
         plot.margin = unit(c(2,85,10,5), "pt"))

P=grid.arrange(pmosaic, pfreq, pmaps, heights = c(10,1.65, 1.65))

ggsave(P, file="figures/haplotypes_main.jpg", width=6.5, height=7, dpi=1200)


###############################################################################################
#########################  SUPPLEMENTARY FIGURE WITH HET RILs  ################################
###############################################################################################

diplorils = unique(subset(RIls_FounderHaplotypeBlocks, ishet == T)$rilname)
Rilshaplotypes_diplo = subset(RIls_FounderHaplotypeBlocks, rilname %in% diplorils)
Rilshaplotypes_diplo = droplevels(Rilshaplotypes_diplo)

Rilshaplotypes_diplo  = do.call(rbind, lapply(split(Rilshaplotypes_diplo, Rilshaplotypes_diplo$rilname), function(x){
  
  founders = strsplit(x$founder, "&")
  g1 = sapply(founders, function(x){x[1]})
  g2 = sapply(founders, function(x){x[length(x)]})
  
  xg1 = xg2 = x
  xg1$founder = g1
  xg2$founder = g2
  xg1$genome = "g1"
  xg2$genome = "g2"
  xg2 = subset(xg2, chrom %in% subset(xg2, ishet == T)$chrom)
  rbind(xg1, xg2)
  
}))


Rilshaplotypes_diplo$rilname = paste0(Rilshaplotypes_diplo$rilname,'_', Rilshaplotypes_diplo$genome)

haplop = haploPlotFromat(Rilshaplotypes = Rilshaplotypes_diplo, widthRILs = 0.95)


haplop_diplo=haploPlotFromat(Rilshaplotypes_diplo,
                             widthRILs = 0.5)


haplop_diplo$rely1[haplop_diplo$genome=='g1'] = haplop_diplo$rely1[haplop_diplo$genome=='g1']*0.4
haplop_diplo$rely2[haplop_diplo$genome=='g1'] = haplop_diplo$rely2[haplop_diplo$genome=='g1']*0.4

haplop_diplo$rely1[haplop_diplo$genome=='g12'] = haplop_diplo$rely1[haplop_diplo$genome=='g12']*0.4
haplop_diplo$rely2[haplop_diplo$genome=='g12'] = haplop_diplo$rely2[haplop_diplo$genome=='g12']*0.4

haplop_diplo$rely1[haplop_diplo$genome=='g2'] = haplop_diplo$rely1[haplop_diplo$genome=='g2']*0.4 + 0.5
haplop_diplo$rely2[haplop_diplo$genome=='g2'] = haplop_diplo$rely2[haplop_diplo$genome=='g2']*0.4 + 0.6

haplop_diplo$rilname = tstrsplit(haplop_diplo$rilname, '_g')[[1]]

haplop_diplo = lapply(split(haplop_diplo, haplop_diplo$panel), function(x){
  x$y1 = as.numeric(factor(x$rilname))-1 + (x$rely1)*0.7
  x$y2 = as.numeric(factor(x$rilname))-1 + (x$rely2)*0.7
  x
})




#"α" "β"
pd = ggplot(haplop_diplo[[2]])+theme_Publication3(base_size = 10)+
  geom_rect(aes(xmin=pos1/1e6,xmax=pos2/1e6, ymin=y1, ymax=y2,fill=founder), color = NA)+
  scale_fill_manual(values = foundercolors, breaks = names(foundercolors))+
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0), breaks = chrompos$midpos/1e6, labels = chrompos$chrom)+
  scale_y_continuous(expand = c(0,0), breaks = rily$y1+0.1, labels = rily$rilname)+
  theme(#axis.text.y =element_blank(),
    axis.line.y = element_line(color = "white"),
    #axis.ticks.y = element_line(color = "white"),
    axis.title = element_blank(),
    #axis.text.x = element_text(size = 15),
    panel.background = element_rect(fill = "white", colour=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(10,10,10,10), "pt"))+
  ylab("")+xlab("")


ggsave(pd, file="figures/haplotypes_diploid_HetRilsB.jpg", width=7, height=9, dpi=1200)


