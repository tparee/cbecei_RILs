############################################################
######## Extract recombination brekpoints ##################
#### i.e., region between hap_j_start to hap_i_end #########

recombination_breakpoints = NULL
for(CHR in c("I","II","III","IV","V","X")){
  rilhaplotypes = as.data.frame(readr::read_csv(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/haplotypes/",CHR,"_rils_foundingHaplotypesBlocks2.csv")))
  breakpoints = do.call(rbind, lapply(split(rilhaplotypes, rilhaplotypes$rilname), function(x){
    if(nrow(x) > 1){
      breakpos = t(apply(rbind(x$pos1[-1],x$pos2[-nrow(x)]),2,sort))
      nbreak_per_haploid_genome = ifelse(x$ishet[-1] | x$ishet[-nrow(x)], 1, 2)
      out = cbind(breakpos, nbreak_per_haploid_genome)
      colnames(out) = c("pos1", "pos2", "nbreak_per_haploid_genome")
      out = as.data.frame(out)
      out$chrom = CHR
      out$rilname = x$rilname[1]
      return(out)
    }else{
      return(NULL)
    }
  }))
  
  recombination_breakpoints = rbind(recombination_breakpoints, breakpoints)
}

#write_csv(recombination_breakpoints, file="/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/haplotypes/recombination_breakpoints.csv")

################################################
######### Plot the map #########################

# the breakpoint is between a range of position (pos1 to pos2)
# Because it founder can be homologous for a small interval they overlap so the brekpoint cannot be assigned to a unique position
# Although intervals are small and don't impact much the map, we can brekpoints between their possible position (pos1:pos2)
# Here, we sample with a uniform probability distribution
# In reality, the true breakpoint position follow a probability ~ recombination map, so we could use the reference map
# But here, we don't do that because:
# i) we don't want to bias our map to look like the reference ones, because we gonna compare them
# ii) The intervals are so small and recombination in Caenorhabditis mostly vary at broadscale level and not fine-scale, so it won't make a noticable difference anyway

set.seed(123)
nsample = 20
sampled_map_pos = t(sapply(1:nrow(recombination_breakpoints), function(i){
  round(runif(nsample, min = recombination_breakpoints$pos1[i], max=recombination_breakpoints$pos2[i]))
}))

colnames(sampled_map_pos) = paste0("sample",1:nsample)
rownames(sampled_map_pos) = 1:nrow(sampled_map_pos)
sampled_map_pos = reshape2::melt(sampled_map_pos)
map = cbind(sampled_map_pos, recombination_breakpoints[as.numeric(sampled_map_pos$Var1), c("chrom",'nbreak_per_haploid_genome')])[,-1]
colnames(map) = c("id",'pos',"chrom",'nbreak_per_haploid_genome')
map$id = paste0(map$id, "_", map$chrom)

map = do.call(rbind, lapply(split(map, map$id), function(breakpoints){
  breakpoints =  breakpoints[order(breakpoints$pos),]
  breakpoints$cumbreaks = cumsum(breakpoints$nbreak_per_haploid_genome)/2
  return(breakpoints)
}))

ggplot(map, aes(pos/1e6, cumbreaks, group = id))+geom_line(alpha=0.3)+facet_wrap(~chrom, nrow=1)

write_csv(map[,c(1:3,5)], file = "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/analysis/temp/cumbreaks_sampled.csv")
