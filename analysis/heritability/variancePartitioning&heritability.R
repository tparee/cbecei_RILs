library(sommer)
setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/")
meta = read_csv("suppl/RILs_sequencing_metadata.csv")
############################################################################
############################################################################
#################### GROWTH RATE  ##########################################
############################################################################
############################################################################

growthrates = read.csv( "phenotypes/growthrates.csv", sep = " ")
growthrates$log_hours_to_starve = log(growthrates$hours_to_starve)
growthrates$strain = tstrsplit(growthrates$strain,"_")[[2]]
growthrates = subset(growthrates, strain %in% meta$rilname[meta$panel == "alpha"])

snps = as.data.frame(read_csv("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/genotypes/beceiPanels_variantsInfo_pruned0.999.csv.gz"))
genotypes <- as.matrix(fread("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/genotypes/beceiPanels_geno_RILs_pruned0.999.csv.gz"))
genotypes = (genotypes-0.5)*2
#genotypes[genotypes == 0]=1

growthrates = subset(growthrates, strain %in% colnames(genotypes))
genotypes = genotypes[,colnames(genotypes) %in% growthrates$strain]

af = apply(genotypes,1,function(x){sum(x+1, na.rm=T)/(2*sum(!is.na(x)))})
genotypes = genotypes[af > 0 & af < 1,]
snps = snps[af > 0 & af < 1,]



############################################################################
#################### SNP-based GRM  ########################################
############################################################################

get.K_ASV = function(GT){
  ASV = scale(GT,center=T,scale=F) %*% t(scale(GT,center=T,scale=F))
  ASV = ASV / (psych::tr(ASV)/(nrow(GT)-1))
  ASV
}

######################################
########## Heritability ##############

growthrates$id = paste0(growthrates$strain,'_', 1:nrow(growthrates))
GT = t(genotypes)
GT = GT[match(growthrates$strain,row.names(GT)),]
row.names(GT) = growthrates$id
#GT[is.na(GT)]=0

#growthrates$id = factor(growthrates$strain, levels = rownames(GT))
K=get.K_ASV(GT)

modh2 <- mmer(
  log_hours_to_starve ~ 1,
  random = ~  block + vs(id, Gu = K),
  data = growthrates
)


var_comp <- summary(modh2)$varcomp

#                                                 VarComp   VarCompSE    Zratio   Constraint
#block.log_hours_to_starve-log_hours_to_starve 0.003785181 0.002205150  1.716519   Positive
#u:id.log_hours_to_starve-log_hours_to_starve  0.029750386 0.005628033  5.286107   Positive
#units.log_hours_to_starve-log_hours_to_starve 0.017463510 0.001289625 13.541546   Positive

Va <- var_comp[2,1]
Ve <- var_comp[3,1]
SE_Va <- var_comp[2,2]
SE_Ve <-var_comp[3,2]

# Derivatives
d_Va <- Ve / (Va + Ve)^2
d_Ve <- -Va / (Va + Ve)^2

# SE of h2 (delta method)
h2 = Va/(Va+Ve)
SE_h2 <- sqrt(d_Va^2 * SE_Va^2 + d_Ve^2 * SE_Ve^2)
print(c(h2,SE_h2)) #  0.63011927 0.04733113


######################################################
########## Genetic variance partitioning ##############

# *******************************
# *** Per recobination domain ***

Karm = get.K_ASV(GT[,snps$domain == 'arm'])
Kcenter = get.K_ASV(GT[,snps$domain == 'center'])

growthrates$id2 = growthrates$id

modh2 <- mmer(
  log_hours_to_starve ~ 1,
  random = ~  block + vs(id, Gu = Karm) +  vs(id2, Gu = Kcenter),
  data = growthrates
)



var_comp <- summary(modh2)$varcomp

var_comp = data.frame(trait='growthrate',domain = c('arm','center'), Va = var_comp[2:3,1],
                      SE_Va = var_comp[2:3,2])

var_comp$Va_rel = var_comp$Va/sum(var_comp$Va)
var_comp$SE_Va_rel = var_comp$SE_Va/sum(var_comp$Va)

write.table(var_comp, file = "analysis/heritability/growthrate_gvariance_partition_recDomains.csv")

# ***********************
# *** Per chromosomes ***

KI = get.K_ASV(GT[,snps$chrom == 'I'])
KII = get.K_ASV(GT[,snps$chrom == 'II'])
KIII = get.K_ASV(GT[,snps$chrom == 'III'])
KIV = get.K_ASV(GT[,snps$chrom == 'IV'])
KV = get.K_ASV(GT[,snps$chrom == 'V'])
KX = get.K_ASV(GT[,snps$chrom == 'X'])

growthrates$id2 = growthrates$id3 =  growthrates$id4 =  growthrates$id5 =  growthrates$id6 =  growthrates$id

modh2 <- mmer(
  log_hours_to_starve ~ 1,
  random = ~ vs(id, Gu = KI) 
  +  vs(id2, Gu = KII)
  +  vs(id3, Gu = KIII) 
  + vs(id4, Gu = KIV)
  + vs(id5, Gu = KV)
  + vs(id6, Gu = KX)
  + block,
  data = growthrates
)

var_comp <- summary(modh2)$varcomp

var_comp = data.frame(trait='growthrate',chrom = c("I","II","III",'IV',"V", "X"), Va = var_comp[1:6,1],
                      SE_Va = var_comp[1:6,2])

var_comp$Va_rel = var_comp$Va/sum(var_comp$Va)
var_comp$SE_Va_rel = var_comp$SE_Va/sum(var_comp$Va)

write.table(var_comp, file = "analysis/heritability/growthrate_gvariance_partition_chrom.csv")

ggplot(var_comp, aes(x=chrom, y = Va_rel, ymin = Va_rel - SE_Va_rel, ymax = Va_rel + SE_Va_rel))+
  geom_pointrange()+xlab("Chromosome")+ylab("Proportion of genomic variance explained")




############################################################################
#################### haplotype-based GRM  ##################################
############################################################################
Ghap = as.matrix(read_csv("analysis/heritability/haplotype_relatedness_matrix/Ghap_haplotypeRelatednessMatrix.csv"))
rownames(Ghap) = colnames(Ghap)
mm=match(growthrates$strain,row.names(Ghap))
Ghap = Ghap[mm,mm]
colnames(Ghap) = row.names(Ghap) = growthrates$id

modh2 <- mmer(
  log_hours_to_starve ~ 1,
  random = ~  block + vs(id, Gu = Ghap),
  data = growthrates
)


var_comp <- summary(modh2)$varcomp

Va <- var_comp[2,1]
Ve <- var_comp[3,1]
SE_Va <- var_comp[2,2]
SE_Ve <-var_comp[3,2]

D = mean(diag(Ghap)) - mean(Ghap)
Va = D*Va
SE_Va = D*SE_Va
# Heritability
h2 <- Va / (Va + Ve)
h2
# Delta method: approximate SE of h²
SE_h2 <- sqrt((SE_Va^2 * Ve^2 + SE_Ve^2 * Va^2) / (Va + Ve)^4)
SE_h2



#################################################################################
# ------------------------------------------------------------------------------#
# --------------------------     SIZE      -------------------------------------#
# ------------------------------------------------------------------------------#
#################################################################################

#"*************************** SNP - BASED ***************************************"

size = as.data.frame(fread("phenotypes/size_summarystat.csv"))
size = subset(size, strain %in% subset(meta, panel == "alpha")$rilname)

snps = as.data.frame(read_csv("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/genotypes/beceiPanels_variantsInfo_pruned0.999.csv.gz"))
genotypes <- as.matrix(fread("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/genotypes/beceiPanels_geno_RILs_pruned0.999.csv.gz"))
genotypes = (genotypes-0.5)*2

genotypes = genotypes[,colnames(genotypes) %in% size$strain]

size = droplevels(subset(size, strain %in% colnames(genotypes) & !is.na(weight_rotated)))
size$weight_rotated = size$weight_rotated/mean(size$weight_rotated)
size$weight_f = size$weight_f/mean(size$weight_f)


afs = apply(genotypes, 1, function(x){mean(x+1,na.rm=T)/2})
genotypes = genotypes[afs > 0 & afs < 1,]
snps = snps[afs > 0 & afs < 1,]

GT = t(genotypes)
#GT[is.na(GT)]=0
K=get.K_ASV(GT)


# *** -----------------------------------------***
# ***  Sexually convergent variance component. ***
# *** -----------------------------------------***

modh2 <- mmer(
  conv ~ 1,
  random = ~ block + vs(strain, Gu = K),
  data = size, weights = weight_rotated
)

vcomp = summary(modh2)$varcomp


#                   VarComp   VarCompSE   Zratio Constraint
#block.conv-conv    0.539196 0.3683031 1.464001   Positive
#u:strain.conv-conv 0.360911 0.1936264 1.863956   Positive
#units.conv-conv    1.058958 0.1671784 6.334298   Positive

# h2 = 0.2736115
Va = vcomp[2,1]
Ve = vcomp[3,1]
SE_Va <- vcomp[2,2]
SE_Ve <- vcomp[3,2]

# Derivatives
d_Va <- Ve / (Va + Ve)^2
d_Ve <- -Va / (Va + Ve)^2

# SE of h2 (delta method)
SE_h2 <- sqrt(d_Va^2 * SE_Va^2 + d_Ve^2 * SE_Ve^2)
h2 = Va/(Va+Ve)
print(c(h2,SE_h2)) #  0.2541862 0.1060181


# *** ----------------------------------------***
# ***  Sexually divergent variance component. ***
# *** ----------------------------------------***

modh2 <- mmer(
  div ~ 1,
  random = ~ block + vs(strain, Gu = K),
  data = size, weights = weight_rotated
)

vcomp = summary(modh2)$varcomp

#                    VarComp  VarCompSE   Zratio Constraint
#block.div-div    0.01933485 0.01626140 1.189003   Positive
#u:strain.div-div 0.08197348 0.02950911 2.777904   Positive
#units.div-div    0.08595788 0.01491149 5.764539   Positive

Va = vcomp[2,1]
Ve = vcomp[3,1]
SE_Va <- vcomp[2,2]
SE_Ve <- vcomp[3,2]

# Derivatives
d_Va <- Ve / (Va + Ve)^2
d_Ve <- -Va / (Va + Ve)^2

# SE of h2 (delta method)
SE_h2 <- sqrt(d_Va^2 * SE_Va^2 + d_Ve^2 * SE_Ve^2)
h2 = Va/(Va+Ve)
print(c(h2,SE_h2)) #  0.48813680 0.09984421


######################################################
########## Genetic variance partitioning ##############

# *******************************
# *** Per recobination domain ***

Karm = get.K_ASV(GT[,snps$domain == 'arm'])
Kcenter = get.K_ASV(GT[,snps$domain == 'center'])

size$id2 = size$id = size$strain
modh2 <- mmer(
  conv ~ 1,
  random = ~ block + vs(id, Gu = Karm) + vs(id2, Gu = Kcenter),
  data = size, weights = weight_rotated
)



var_comp <- summary(modh2)$varcomp

var_comp = data.frame(trait='size',domain = c('arm','center'), Va = var_comp[2:3,1],
                      SE_Va = var_comp[2:3,2])

var_comp$Va_rel = var_comp$Va/sum(var_comp$Va)
var_comp$SE_Va_rel = var_comp$SE_Va/sum(var_comp$Va)

write.table(var_comp, file = "analysis/heritability/sizeConv_gvariance_partition_recDomains.csv")

####
modh2 <- mmer(
  div ~ 1,
  random = ~ block + vs(id, Gu = Karm) + vs(id2, Gu = Kcenter),
  data = size, weights = weight_rotated
)



var_comp <- summary(modh2)$varcomp

var_comp = data.frame(trait='size',domain = c('arm','center'), Va = var_comp[2:3,1],
                      SE_Va = var_comp[2:3,2])

var_comp$Va_rel = var_comp$Va/sum(var_comp$Va)
var_comp$SE_Va_rel = var_comp$SE_Va/sum(var_comp$Va)

write.table(var_comp, file = "analysis/heritability/sizeDiv_gvariance_partition_recDomains.csv")

# ***********************
# *** Per chromosomes ***

KI = get.K_ASV(GT[,snps$chrom == 'I'])
KII = get.K_ASV(GT[,snps$chrom == 'II'])
KIII = get.K_ASV(GT[,snps$chrom == 'III'])
KIV = get.K_ASV(GT[,snps$chrom == 'IV'])
KV = get.K_ASV(GT[,snps$chrom == 'V'])
KX = get.K_ASV(GT[,snps$chrom == 'X'])

size$id2 = size$id3 =  size$id4 =  size$id5 =  size$id6 =  size$id

#### Sexually-convergent 

modh2 <- mmer(
  conv ~ 1,
  random = ~ block + vs(id, Gu = KI) 
  +  vs(id2, Gu = KII)
  +  vs(id3, Gu = KIII) 
  + vs(id4, Gu = KIV)
  + vs(id5, Gu = KV)
  + vs(id6, Gu = KX),
  data = size, weights = weight_rotated
)

var_comp <- summary(modh2)$varcomp

var_comp = data.frame(trait='size',chrom = c("I","II","III",'IV',"V", "X"), Va = var_comp[2:7,1],
                      SE_Va = var_comp[2:7,2])

var_comp$Va_rel = var_comp$Va/sum(var_comp$Va)
var_comp$SE_Va_rel = var_comp$SE_Va/sum(var_comp$Va)

write.table(var_comp, file = "analysis/heritability/sizeConv_gvariance_partition_chrom.csv")

ggplot(var_comp, aes(x=chrom, y = Va_rel, ymin = Va_rel - SE_Va_rel, ymax = Va_rel + SE_Va_rel))+
  geom_pointrange()+xlab("Chromosome")+ylab("Proportion of genomic variance explained")


#### Sexually-Divergent
modh2 <- mmer(
  div ~ 1,
  random = ~ block + vs(id, Gu = KI) 
  +  vs(id2, Gu = KII)
  +  vs(id3, Gu = KIII) 
  + vs(id4, Gu = KIV)
  + vs(id5, Gu = KV)
  + vs(id6, Gu = KX),
  data = size, weights = weight_rotated
)

var_comp <- summary(modh2)$varcomp

var_comp = data.frame(trait='size',chrom = c("I","II","III",'IV',"V", "X"), Va = var_comp[2:7,1],
                      SE_Va = var_comp[2:7,2])

var_comp$Va_rel = var_comp$Va/sum(var_comp$Va)
var_comp$SE_Va_rel = var_comp$SE_Va/sum(var_comp$Va)

write.table(var_comp, file = "analysis/heritability/sizeDiv_gvariance_partition_chrom.csv")

ggplot(var_comp, aes(x=chrom, y = Va_rel, ymin = Va_rel - SE_Va_rel, ymax = Va_rel + SE_Va_rel))+
  geom_pointrange()+xlab("Chromosome")+ylab("Proportion of genomic variance explained")


#"*************************** HAPLOTYPE - BASED ***************************************"

Ghap = as.matrix(read_csv("analysis/heritability/haplotype_relatedness_matrix/Ghap_haplotypeRelatednessMatrix.csv"))
rownames(Ghap) = colnames(Ghap)
mm=match(size$strain,row.names(Ghap))
Ghap = Ghap[mm,mm]

### Sexually convergent
modh2 <- mmer(
  conv ~ 1,
  random = ~ block + vs(strain, Gu = Ghap),
  data = size, weights = weight_rotated
)

var_comp <- summary(modh2)$varcomp
Va <- var_comp[2,1]
Ve <- var_comp[3,1]
SE_Va <- var_comp[2,2]
SE_Ve <-var_comp[3,2]

D = mean(diag(Ghap)) - mean(Ghap)
Va = D*Va
SE_Va = D*SE_Va
# Heritability
h2 <- Va / (Va + Ve)
h2
# Delta method: approximate SE of h²
sqrt((SE_Va^2 * Ve^2 + SE_Ve^2 * Va^2) / (Va + Ve)^4)



### Sexually Divergent
modh2 <- mmer(
  div ~ 1,
  random = ~ block + vs(strain, Gu = Ghap),
  data = size, weights = weight_rotated
)

var_comp <- summary(modh2)$varcomp
Va <- var_comp[2,1]
Ve <- var_comp[3,1]
SE_Va <- var_comp[2,2]
SE_Ve <-var_comp[3,2]


D = mean(diag(Ghap)) - mean(Ghap)
Va = D*Va
SE_Va = D*SE_Va
# Heritability
h2 <- Va / (Va + Ve)
h2

# Derivatives
d_Va <- Ve / (Va + Ve)^2
d_Ve <- -Va / (Va + Ve)^2

# SE of h2 (delta method)
SE_h2 <- sqrt(d_Va^2 * SE_Va^2 + d_Ve^2 * SE_Ve^2)

# h2 = 0.4530909, SE = 0.0995689



############################################################################
#################### PLOTS  ################################################
############################################################################
library(gridExtra)
source("utils.R")
varpartition_chrom = rbind(as.data.frame(read.csv("analysis/heritability/sizeConv_gvariance_partition_chrom.csv", sep=" ")),
                           as.data.frame(read.csv("analysis/heritability/sizeDiv_gvariance_partition_chrom.csv", sep=" ")),
                           as.data.frame(read.csv("analysis/heritability/growthrate_gvariance_partition_chrom.csv", sep=" ")))

varpartition_recDomains = rbind(as.data.frame(read.csv("analysis/heritability/sizeConv_gvariance_partition_recDomains.csv", sep=" ")),
                                as.data.frame(read.csv("analysis/heritability/sizeDiv_gvariance_partition_recDomains.csv", sep=" ")),
                                as.data.frame(read.csv("analysis/heritability/growthrate_gvariance_partition_recDomains.csv", sep=" ")))


varpartition_chrom$grp = c(rep('1', 6),rep('2', 12))
varpartition_recDomains$grp = c(rep('1', 2),rep('2', 4))


p1 = ggplot(varpartition_chrom, aes(x=chrom, y = Va_rel, ymin = Va_rel - SE_Va_rel, ymax = Va_rel + SE_Va_rel, group = grp, color = grp))+
  geom_pointrange(size=0.1, position = position_dodge(width = 0.3))+
  xlab("Chromosomes")+ylab("Proportion of \n genomic variance explained")+
  facet_grid(.~trait)+
  theme_Publication3()+
  geom_vline(data = data.frame(trait = "size"), aes(xintercept = -Inf), size=0.75)+
  coord_cartesian(ylim=c(0,0.75))+
  scale_color_manual(values = c('grey65', 'black'))+
  theme(legend.position = "none")

p2 = ggplot(varpartition_recDomains, aes(x=domain, y = Va_rel, ymin = Va_rel - SE_Va_rel, ymax = Va_rel + SE_Va_rel, group = grp, color = grp))+
  geom_pointrange(size=0.1, position = position_dodge(width = 0.3))+
  xlab("Domain")+ylab("")+
  facet_grid(.~trait)+
  theme_Publication3()+
  geom_vline(data = data.frame(trait = "size"), aes(xintercept = -Inf), size=0.75)+
  coord_cartesian(ylim=c(0,1.1))+
  scale_color_manual(values = c('grey65', 'black'))+
  theme(legend.position = "none")


p=grid.arrange(p1,p2, nrow=1, widths = c(1,0.7))


ggsave(p, file="figures/fig_variancePartition.png", width=3.7, height=1.35, dpi=1200)







