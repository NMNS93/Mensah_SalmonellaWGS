# clean_genotypes.R
# Create binary table of gene presence/absence for each isolate
## TODO fix variable names

# Import functions
source("scripts/functions.R")

# Load resfinder genotype data
resfdf = read.csv('data/genotype/ariba_resfinder.csv', row.names=1) %>% clndata()
# Load pointfinder genotype data
pointmutdf = read.csv('data/genotype/pointfinder_pointmutations.csv')
# Merge resfinder and pointfinder data (only gyrA calls were found in this dataset)
indiv_gyra = pointmutdf %>% select(sample, gene) %>% mutate(value=1) %>% spread(gene, value)
genotemp1 = full_join(resfdf, indiv_gyra, by="sample") %>% mutate_all(funs(replace(.,is.na(.),0)))
# Merge columns for AMR gene families (except str)
mygenes = c("qnrb", "aadA", "aph_3", "aph_4", "blaCARB", "blaTEM", "cat", "cml", "dfrA", "floR", "mph_B", "sul", "tet", "gyrA")
genotemp2 = genotemp1
# Note: The function mergebin searches for columns matching strings in 'mygenes' and merges
for (gene in mygenes){
  genotemp2 = mergebin(genotemp2, gene, gene)
}
# Merge strA and strB with '1' only where both genes present 
genotype_df = mergestr(genotemp2, "str", "str_AB")
# Write full binary genotype table
write.csv(genotype_df, file="data/genotype/ariba_genotype_tally.csv", row.names=FALSE)

# Create genotype table for ECOFF comparisons
# Load genotype mapping data
resamr = read.csv('files/genotype_conversion.csv')
rese = map_amrs(genotype_df, select(resamr, amr, ecc), "ecc")
write.csv(rese, 'data/clean/genotype_ecoff.csv', row.names = FALSE)
# Create genotype clinbp
resclin = map_amrs(genotype_df, select(resamr, amr, clin), "clin")
write.csv(resclin, 'data/clean/genotype_clinbp.csv', row.names = FALSE)
# Create genotype counts table
count_output = genotype_df  %>% gather(gene, result, -sample) %>% count(gene,result) %>% filter(result == 1) %>% select(gene,count=n)
write.csv(count_output, 'results/dataframes/tallycount.csv')
# Write out original genotype dataframe
write.csv(genotemp1, 'results/dataframes/genotype_full_dataframe.csv', row.names=FALSE)

# Create genotype accession column
gen_access = read.csv('data/genotype/ariba_resfinder.csv', row.names=1) %>% 
  select(contains("ref")) %>% rownames_to_column %>% mutate(rowname = clnsample(rowname)) %>%
  select(sample=rowname, everything()) %>% unite(refseqs, contains("ref"), sep=",") %>% 
  mutate(refseqs = str_replace_all(refseqs, c(",NA"="","NA,"="","NA"=""),""))
write.csv(gen_access, 'results/dataframes/accessions.csv')