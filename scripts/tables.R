source('scripts/functions.R')
# tables.R - Create figures for the paper

# Table 1 - Genotype concordance
# Layout: Antimicrobials, [Ecoff, nwt, wt, genoconcord], [Clinbp, nwt, wt, genoconcord]
ecoff_data = read.csv('results/stats/ecoff_stats.csv') %>% select(amr, ecopos=POS, econeg=NEG, TP, TN) %>% mutate(ecoconc=(TP+TN)*100/(ecopos+econeg))
clin_data = read.csv('results/stats/resc_stats.csv') %>% select(amr, clinpos=POS, clinneg=NEG, TP, TN) %>% mutate(clinconc=(TP+TN)*100/(clinpos+clinneg))
table1 = full_join(ecoff_data, clin_data, by = c("amr", "TP", "TN"))

# Table 2 - Sensitivity/specificity data
agar_data = read.csv('results/stats/agar_stats.csv') %>% mutate(agarsens=merger(sens,sens_l,sens_u)) %>% mutate(agarspec=merger(spec,spec_l,spec_u)) %>% select(amr, agarsens, agarspec)
disk_data = read.csv('results/stats/disk_stats.csv') %>% mutate(disksens=merger(sens,sens_l,sens_u)) %>% mutate(diskspec=merger(spec,spec_l,spec_u)) %>% select(amr, disksens, diskspec)
ecoff_data2 = read.csv('results/stats/resc_stats.csv') %>% mutate(ressens=merger(sens,sens_l,sens_u)) %>% mutate(resspec=merger(spec,spec_l,spec_u)) %>% select(amr, ressens, resspec)
table2 = full_join(disk_data, agar_data, by="amr") %>% full_join(ecoff_data2, by="amr")

# Table 3 - Resistance Genes Tally - Counts of AMR reference sequences (includes genes identified by ARIBA as partial, incomplete, fragmented)
class = read.csv('files/amr_class.csv')
points = read.csv('data/genotype/pointfinder_pointmutations.csv')
refs = read.csv('data/genotype/ariba_resfinder.csv') %>% mutate(sample=clnsample(sample)) %>% select(sample, contains('ref_seq')) %>% left_join(select(points,sample,gene), by="sample")
intermed = refs %>% gather(gene,ref,-sample) %>% drop_na() %>% select(ref)
table3 = intermed %>% group_by(ref) %>% summarise(counts=length(ref))

# Table 4 - Present resistance genes tally - Counts of AMR reference sequences only where ARIBA match column is "yes"
# Get columns from ARIBA tables ending in 'match'
refs = read.csv('data/genotype/ariba_resfinder.csv') %>% mutate(sample=clnsample(sample))
match_cols = refs %>% select(contains("match")) %>% colnames()
# Create data list.
datalist = list()
# Loop over indexes of reference columns
 for (i in seq(length(match_cols))){
  # Get string of current reference column
  current_col = match_cols[i]
  # Get string of matching 'match' column
  current_col_ref = str_replace(current_col, "match", "ref_seq")
  # Create dataframe of sample, reference sequence where match is "yes"
  df = refs %>% select(sample, first=current_col, second=current_col_ref) %>% filter(str_detect(first, "yes")) %>%
    select(sample,ref=second)
  # Append to data list
  datalist[[i]] = df
 }
# Merge all tables in datalist
big_table = do.call(rbind, datalist)
# Add gyrA results from pointfinder
points = read.csv('data/genotype/pointfinder_pointmutations.csv') %>% select(sample,ref=gene)
table4 = big_table %>% full_join(points, by=c("sample","ref")) %>% select(ref) %>% group_by(ref) %>% summarise(counts=length(ref))

# Produce a final table of reference sequence data
table5 = table4 %>% select(ref, yes_counts=counts) %>% full_join(table3, by="ref")

# Write out data
write.csv(table1, 'results/tables/wgsconcordance.csv')
write.csv(table2, 'results/tables/ASTcomparison.csv')
write.csv(table3, 'results/tables/genotype_counts.csv')
write.csv(table4, 'results/tables/genotype_counts_yesrefs.csv')
write.csv(table5, 'results/tables/genotype_counts_full.csv')
