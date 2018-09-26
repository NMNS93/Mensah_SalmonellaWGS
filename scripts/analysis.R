# analysis.R
# Statistical analyses for manuscript

# Import functions
source('scripts/functions.R')

# Read cleaned data
brothe = readcsv_gather('data/clean/broth_ecoff.csv')
brothc = readcsv_gather('data/clean/broth_clinbp.csv')
diskc = readcsv_gather('data/clean/disk_clinbp.csv')
agarc = readcsv_gather('data/clean/agar_clinbp.csv')
rese = readcsv_gather('data/clean/genotype_ecoff.csv')
resc = readcsv_gather('data/clean/genotype_clinbp.csv')

# Create dataframes matching ecoff and clinical data against broth
ecoff_match = mframe(brothe, rese, "rese")
clin_match = mframe(brothc, diskc, "diskc") %>% merge(mframe(brothc, agarc, "agarc"), by=c("sample", "amr")) %>%
  merge(mframe(brothc, resc, "resc"), by=c("sample", "amr")) %>% gather(test, result, -sample, -amr)

# Tables of resistance gene counts for paperout
# Vector of all binary count input data
gathered_data = list(list(brothe, "brothe"),  list(rese,"rese"), list(brothc, "brothc"), list(diskc, "diskc"), list(agarc, "agarc"), list(resc, "resc"))
# Gather counts of resistance for each AST table
count_data = lapply(gathered_data, function(x){newname = x[[2]]; count(x[[1]],amr,result) %>% filter(result==1) %>% select(amr, n) %>% rename(!!newname:=n)})
# Join resistance count tables and write out
result_table1 = reduce(count_data, full_join, by="amr")
write.csv(result_table1, 'results/dataframes/resistance_counts.csv')

# ECOFF comparison statistics
ecoff_stats = select(ecoff_match, result=rese, everything())
write.csv(ecoff_match, 'results/dataframes/genotype_match.csv')
write.csv(get_stats(ecoff_stats), 'results/stats/ecoff_stats.csv')

# ECOFF comparison statistics, excluding sulfonamides (sulfisoxazole and sulfamethoxazole)
ecoff_stats_nosulfsmx = brothe %>% filter(!amr =="SUL") %>% filter(!amr=="SMX") %>% 
  mframe(rese, "rese") %>% select(result=rese, everything()) %>% get_stats()
write.csv(ecoff_stats_nosulfsmx, 'results/stats/ecoff_nosulfsmx_stats.csv')

# CBP comparison statistics
write.csv(clin_match, 'results/dataframes/clin_match.csv', row.names=FALSE)
disk_match = clin_match %>% filter(test=="diskc") %>% select(sample, amr, result)
write.csv(get_stats(disk_match), 'results/stats/disk_stats.csv')
agar_match = clin_match %>% filter(test=="agarc") %>% select(sample, amr, result)
write.csv(get_stats(agar_match), 'results/stats/agar_stats.csv')
resc_match = clin_match %>% filter(test=="resc") %>% select(sample, amr, result)
write.csv(get_stats(resc_match), 'results/stats/resc_stats.csv')

# CBP comparison statistics, excluding sulfamethoxazole and streptomycin
resc_stats_nosmxstr = resc_match %>% filter(!amr=="SMX") %>% filter(!amr=="STR") %>% get_stats
write.csv(resc_stats_nosmxstr, 'results/stats/clinbpres_nosmxstr_stats.csv')