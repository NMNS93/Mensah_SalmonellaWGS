source('scripts/functions.R')
# funcholder
merger = function(main,low,up){
  mainr = round(main,2)
  lowr = round(low,2)
  upr = round(up,2)
  output = as.character(c(mainr," (",lowr,"-",upr,")"))
  #return(paste(output,sep='',collapse=''))
  return(paste(mainr,' (',lowr, '-', upr,')',sep=''))
}

# tables.R - Create figures for the paper

# Figure 1 - Genotype concordance
# Layout: Antimicrobials, [Ecoff, nwt, wt, genoconcord], [Clinbp, nwt, wt, genoconcord]

ecoff_data = read.csv('results/stats/ecoff_stats.csv') %>% select(amr, ecopos=POS, econeg=NEG, TP, TN) %>% mutate(ecoconc=(TP+TN)*100/(ecopos+econeg))
clin_data = read.csv('results/stats/resc_stats.csv') %>% select(amr, clinpos=POS, clinneg=NEG, TP, TN) %>% mutate(clinconc=(TP+TN)*100/(clinpos+clinneg))
table1 = full_join(ecoff_data, clin_data, by = c("amr", "TP", "TN"))


# Table 2 - Sensspecdata
agar_data = read.csv('results/stats/agar_stats.csv') %>% mutate(agarsens=merger(sens,sens_l,sens_u)) %>% mutate(agarspec=merger(spec,spec_l,spec_u)) %>% select(amr, agarsens, agarspec)
disk_data = read.csv('results/stats/disk_stats.csv') %>% mutate(disksens=merger(sens,sens_l,sens_u)) %>% mutate(diskspec=merger(spec,spec_l,spec_u)) %>% select(amr, disksens, diskspec)
ecoff_data2 = read.csv('results/stats/resc_stats.csv') %>% mutate(ressens=merger(sens,sens_l,sens_u)) %>% mutate(resspec=merger(spec,spec_l,spec_u)) %>% select(amr, ressens, resspec)
table2 = full_join(disk_data, agar_data, by="amr") %>% full_join(ecoff_data2, by="amr")

# Table 3 - Resgenes tally
class = read.csv('files/amr_class.csv')
points = read.csv('data/genotype/pointfinder_pointmutations.csv')
refs = read.csv('data/genotype/ariba_resfinder.csv') %>% mutate(sample=clnsample(sample)) %>% select(sample, contains('ref_seq')) %>% left_join(select(points,sample,gene), by="sample")
intermed = refs %>% gather(gene,ref,-sample) %>% drop_na() %>% select(ref)
table3 = intermed %>% group_by(ref) %>% summarise(counts=length(ref))


write.csv(table1, 'results/tables/wgsconcordance.csv')
write.csv(table2, 'results/tables/ASTcomparison.csv')
write.csv(table3, 'results/tables/genotype_counts.csv')
