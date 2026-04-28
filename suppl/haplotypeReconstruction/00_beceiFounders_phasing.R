library(readr)
library(data.table)


CHR = 'III'

for(CHR in c("IV","V",'X')){
source("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/haplotypeReconstruction/functions/functions_beceiFounders_phasing.R")
meta = as.data.frame(readxl::read_xlsx(path = "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/private/meta_V2.xlsx"))
snps = read_csv(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/genotype_stringentFiltering/",CHR,"_snps_becei_RILs&Pools_stringentFiltering.csv.gz"))
rils = read_csv(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/genotype_stringentFiltering/",CHR,"_geno_becei_RILs_stringentFiltering.csv.gz"))
foundersgeno_ll = as.data.frame(read_csv(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/pools/",CHR,"_founder_genotype_loglikelihood.csv")))
rils = as.matrix(rils)
snps = as.data.frame(snps)
if(sum(foundersgeno_ll$POS == snps$POS) != nrow(snps)){stop("snps info in founders genotype likelihood and rils genotypes do no match")}
foundersgeno_ll = foundersgeno_ll[,-1:-4]
foundersgeno_ll = as.matrix(foundersgeno_ll)



# code_to_gt = colnames(foundersgeno_ll)
# code_to_gt = do.call(rbind,lapply(strsplit(code_to_gt,";"), function(x){as.numeric(tstrsplit(x,":")[[2]])}))
# rownames(code_to_gt) =  colnames(foundersgeno_ll)
# 
# ml_foundersgeno = do.call(rbind, lapply(1:nrow(foundersgeno_ll), function(i){
#   gt = code_to_gt[which.max(foundersgeno_ll[i,]),]
#   gq = diff(sort(foundersgeno_ll[i,], decreasing = T)[1:2])*-1
#   out = matrix(c(gt,gq), nrow=1)
#   colnames(out) = c("FA","FB","FM","GQ")
#   out
# }))
# 
# ml_foundersgeno  = as.data.frame(ml_foundersgeno)
# 
# exp(0)/(1+exp(0))
# 
# 
# foundersgeno_ll2 = foundersgeno_ll[match(phasedfounderhaplotypes$POS, snps$POS),]
# 
# test = do.call(rbind,lapply(1:nrow(foundersgeno_ll2), function(i){
#   if(i %% 1000 == 0){print(i)}
#   code = paste0("FA:", mean( unlist(phasedfounderhaplotypes[i,c("FA.g1","FA.g2")])),";",
#                 "FB:", mean( unlist(phasedfounderhaplotypes[i,c("FB.g1","FB.g2")])),";",
#                 "FM:", mean( unlist(phasedfounderhaplotypes[i,c("FM.g1","FM.g2")])))
#   
#   mll = foundersgeno_ll2[i,which.max(foundersgeno_ll2[i,])]
#   ll = foundersgeno_ll2[i,code]
#   
#   hamm = sum(abs(code_to_gt[names(mll),] - code_to_gt[code,]))
#   
#   data.frame(ml_geno = names(mll), geno = code, mll = mll, ll = ll, deltall = mll - ll, hammdist = hamm )
# }))
# 
# test$cumsumhammdist = cumsum(test$hammdist)
# ggplot(cbind(pos=phasedfounderhaplotypes$POS, cm = phasedfounderhaplotypes$cM ,index = 1:nrow(test),test),
#        aes(pos, cumsumhammdist))+geom_point()
# 
# ggplot(cbind(pos=phasedfounderhaplotypes$POS, cm = phasedfounderhaplotypes$cM ,index = 1:nrow(test),test),
#        aes(pos, deltall))+geom_point()


#####################################
#### format the input data ##########

# Split rils by cross A and cross B
rilsA = rils[,colnames(rils) %in% subset(meta, panel == "α")$rilname]
rilsB = rils[,colnames(rils) %in% subset(meta, panel == "β")$rilname]

info = snps

# First, InferFoundersHaploBlocks look at the different possible haplotypes in each genomic window
# the most frequent haplotypes are assumed to be acestral
# each window need to be big enough to diferenciate founders  (default is 500SNP)
# And not too large in genetid distance => if too much recombination, the ancestral haplotype are too broken to be recognize
# If it is the case, the window is subdivided
# Then the major haplotypes are attributed to each founders in the way that minimize the change from the initially inferred genotypes (given in founder1 and founder2)
# The function return the four founder hap order for each window, in a list
founderhaplotypes = InferFoundersHaploBlocks(foundersgeno_ll, rilsA, rilsB, info,nsnpWin = 500, maxSizeCM=3)
#This function phase the haplotype block inferred above in a way that minimize the number of breakpoints
phasedfounderhaplotypes = phaseHaplotypes(founderhaplotypes, info)
phasedfounderhaplotypes = cbind(info[phasedfounderhaplotypes$whichsnp,],phasedfounderhaplotypes)
#=> More details of function are provided in "functions_beceiFounders_phasing.R"


### Remove snps where allelic count of crossPool are unlikely given founders haplotype inferred genotypes
### At this stage, we want to have founding haplotype "backbones" that are robust and with no doubt
### Later, we will add all the snps filtered out so that are coherent with these backbones and RILs
bf_threshold = 0.05/(nrow(phasedfounderhaplotypes)*3) # Alpha / number of test
exp_p = list(CrossApool = apply(phasedfounderhaplotypes[,c("FA.g1","FA.g2","FM.g1", "FM.g2")],1,mean),
              CrossBpool = apply(phasedfounderhaplotypes[,c("FB.g1","FB.g2","FM.g1", "FM.g2")],1,mean),
              CrossCpool = apply(phasedfounderhaplotypes[,c("FM.g1", "FM.g2")],1,mean)/2)

AD = fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/pools/",CHR,"_AllelicDepth_becei_CrossPools_stringentFiltering.csv.gz"))
AD = as.data.frame(AD)
AD = AD[phasedfounderhaplotypes$whichsnp,]
AD = as.matrix(AD[,grepl("Cross", colnames(AD))])
REF = do.call(cbind, lapply(1:ncol(AD), function(j){as.numeric(tstrsplit(AD[,j],",")[[1]])}))
ALT = do.call(cbind, lapply(1:ncol(AD), function(j){as.numeric(tstrsplit(AD[,j],",")[[2]])}))
colnames(REF) = colnames(ALT) = tstrsplit(colnames(AD),"_AD")[[1]]
nsnp = nrow(REF)
pvals = do.call(cbind, lapply(names(exp_p), function(thiscross){
  ref = REF[,thiscross]
  alt = ALT[,thiscross]
  dp = ref+alt
  p = exp_p[[thiscross]]
  p[p==0] = 10e-6
  p[p==1] = 1-10e-6
  sapply(1:nsnp, function(i) binom.test(alt[i], dp[i], p[i])$p.value)
}))

keep = apply(pvals, 1, function(x){sum(x < bf_threshold) == 0})
#which(!keep)
phasedfounderhaplotypes = phasedfounderhaplotypes[keep,]


#check across 4cM windows that the founding haplotypes are the most common ones
foundermajor = lapply(seq(1,max(phasedfounderhaplotypes$cM)-1,1), function(gpos){
  print(gpos)
  win = which(phasedfounderhaplotypes$cM > (gpos - 2) & phasedfounderhaplotypes$cM < (gpos + 2))
  #gisze=diff(info$cM[phasedfounderhaplotypes[c(min(win),max(win)),1]])

  p = phasedfounderhaplotypes[win,c("whichsnp","FA.g1","FA.g2","FB.g1","FB.g2","FM.g1", "FM.g2" )]

  gtrils=cbind(rilsA, rilsB)[p$whichsnp,]
  gtrils[gtrils==0.5]=NA
  gx = cbind(p[,-1], gtrils)
  gx = gx[,apply(gx,2, function(y){mean(is.na(y))}) < 0.3]
  groups = cutree(hclust(dist(t(gx))), h = 0)
  npergroups = table(groups)
  foundergroups = unique(groups[1:(ncol(p)-1)])
  nfounderhaplo=length(foundergroups)
  

  #Here we verify that the founders are the major haplotypes
  foundersAreMajorHaplo=sum((names(sort(npergroups, decreasing = T)) %in% foundergroups)[1:nfounderhaplo])==nfounderhaplo
  notfounder = npergroups[!(names(npergroups) %in% foundergroups)]
  nbroken=sum(notfounder)
  frequenthaps = sort(notfounder[notfounder >= 5],decreasing = T)
  sort(notfounder,decreasing = T)
  # If an hap is frequent and it's not a founding hap, it's potentially suspicious 
  # (that imply that the exact same CO happens several independent times or that a specific recombinant haplotype spread over generation)
  # If these haplotype are only supported by a few snps (let's say <= 10 that differ from a inferred founder haplotype),
  # that may be some snps that are prone to genotyping error or alignment error

  suspicious_snps = lapply(as.numeric(names(frequenthaps)), function(g){
    hx = apply(gx[,groups==g],1,mean,na.rm=T)
    divsnp = sapply(2:ncol(p), function(j) which(p[,j] != hx))
    lowdiv = sapply(divsnp, length)  <= 10
    sort(unlist(divsnp[lowdiv]))
  })
  
  suspicious_snps = unique( unlist(suspicious_snps) )
  

  list(foundersAreMajorHaplo=foundersAreMajorHaplo,
             nbroken=nbroken,
       suspicious_snps = win[suspicious_snps])
})

sus = unique(unlist(lapply(foundermajor, function(x){x$suspicious_snps})))
print(sus)
if(length(sus) > 0){phasedfounderhaplotypes = phasedfounderhaplotypes[-sus,]}
#unlist(lapply(foundermajor, function(x){x$foundersAreMajorHaplo}))


phasedfounderhaplotypes=phasedfounderhaplotypes[,-which(colnames(phasedfounderhaplotypes) == "whichsnp")]
save(phasedfounderhaplotypes, file= paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/haplotypeReconstruction/temp/phased_founders_base/",CHR,"_beceiFounderPhasedHaplotypes_base.Rdata"))
}



fgt_old = as.data.frame(fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/becei_rils/genotypes/chr",CHR,"_beceiPanels_geno_founders_haplotypeResolved_V1.csv.gz")))
snps_old = as.data.frame(fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/becei_rils/genotypes/chr",CHR,"_beceiPanels_snps_founders&Rils_haplotypeResolved_V1.csv.gz")))
mm = match(phasedfounderhaplotypes$POS, snps_old$POS)
fgt_old = fgt_old[mm,]

which(!(fgt_old[,"FA.g2"] == phasedfounderhaplotypes[,"FA.g2"]))
i=5264
cbind(fgt_old[i + -5:5,], phasedfounderhaplotypes[i + -5:5,])
# wins = get.win(nrow(pvals),winsize = 100,minsize = 25)
# pvals_comb = do.call(rbind,lapply(wins, function(win){
#   meanpos = mean(info$POS[win])
#   P = pvals[win,]
#   p_comb = apply(P, 2, function(p){
#     X <- -2 * sum(log(p))
#     pchisq(X, df = 2*length(p), lower.tail = FALSE)
#   })
#   c(pos=meanpos, p_comb)
# }))

# pvals_comb = as.data.frame(pvals_comb)
# ggplot(pvals_comb, aes(pos, -log10(p_CrossBpool)))+geom_point()


















