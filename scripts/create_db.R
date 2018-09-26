# create_db.R
# Produces a csv file that can be queried for isolate-level results. Contains following headers:
# Sample, AMR, broth, agar, disk, ECOFF, CBP, genotype, ariba_partial, ariba_interrupted, stats_results...

# Import functions
source('scripts/functions.R')

# Read phenotype data for database
input_samples = read.csv('files/sample_list.csv')
ecoff_settings = read.csv('files/epidemiological_cutoffs.csv') %>% `colnames<-`(c("antimicrobial", "amr", "ecoff_conc", "source", "comment"))
cbp_settings = read.csv('files/clinical_breakpoints.csv') %>% `colnames<-`(c("antimicrobial","amr","cbp_conc","cbp_disk_zone","reference","source","comment","comment2","comment3"))
disk_raw = read.csv('data/phenotype/disk.csv') %>% gather(amr, disk_zone, -sample)
broth_raw = read.csv('data/phenotype/broth.csv') %>% gather(amr, broth_conc, -sample)
agar_raw = read.csv('data/phenotype/agar.csv') %>% gather(amr, agar_conc, -sample)

# Read genotype data for database
genotype_ariba_res = read.csv('results/dataframes/genotype_full_dataframe.csv')
ariba_genotype = merge_geno_wrapper(genotype_ariba_res, 1, "genotype")
ariba_data = read.csv('data/genotype/ariba_resfinder.csv') %>% mutate(sample=clnsample(sample))
ariba_partial =merge_geno_wrapper(ariba_data, "partial", "ariba_partial")
ariba_interrupted = merge_geno_wrapper(ariba_data, "interrupted", "ariba_interrupted")
ariba_fragmented = merge_geno_wrapper(ariba_data, "fragmented", "ariba_fragmented")

# Get genotype refseqs for database
geno_refseqs = read.csv('data/genotype/ariba_resfinder.csv', row.names=1) %>% 
  select(contains("ref")) %>% rownames_to_column %>% mutate(rowname = clnsample(rowname)) %>%
  select(sample=rowname, everything()) %>% unite(refseqs, contains("ref"), sep=",") %>% 
  mutate(refseqs = str_replace_all(refseqs, c(",NA"="","NA,"="","NA"=""),""))

# Read stats results
clin_match = read.csv('results/dataframes/clin_match.csv') %>% spread(test,result) %>% select(sample, amr, agar_clin=agarc, disk_clin=diskc, geno_clin=resc)
ecoff_match = read.csv('results/dataframes/genotype_match.csv') %>% select(sample, amr, geno_ecoff=rese)

# Get agar and broth dilution result symbols
broth_for_sym = read.csv('data/phenotype/broth.csv') 
agar_for_sym = read.csv('data/phenotype/agar.csv')
broth_sym = phenotype_symbol(broth_for_sym) %>% select(broth_symbol=symbol, everything())
agar_sym = phenotype_symbol(agar_for_sym) %>% select(agar_symbol=symbol, everything())

# Build database
database = input_samples %>% left_join(broth_raw, by="sample") %>% left_join(broth_sym, by=c("sample","amr")) %>% 
  left_join(disk_raw, by=c("sample","amr")) %>% left_join(agar_raw, by=c("sample","amr")) %>% left_join(agar_sym, by=c("sample","amr")) %>%
  left_join(select(ecoff_settings, amr, ecoff_conc), by="amr") %>% 
  left_join(select(cbp_settings, amr, cbp_conc, cbp_disk_zone), by="amr") %>%
  left_join(clin_match, by=c("sample","amr")) %>% left_join(ecoff_match, by=c("sample","amr")) %>%
  left_join(ariba_partial, by="sample") %>% left_join(ariba_interrupted, by="sample") %>% left_join(ariba_fragmented, by="sample") %>% 
  left_join(ariba_genotype, by="sample") %>% left_join(geno_refseqs, by="sample")

# Filter database amr values to antimicrobials of interest
db_filtered = filter(database, amr %in% ecoff_settings$amr)

# Example queries:
## genotype-negative results where CBPs were within 1 dilution of AST MIC
# db_filtered %>% mutate(broth_conc=as.numeric(broth_conc)) %>% filter(geno_clin=="FN" | geno_clin=="FP") %>% mutate(cblo=cbp_conc/2, cbhi=cbp_conc*2) %>% filter(broth_conc<=cbhi & broth_conc >=cblo)
## genotype-negative results where ECOFFs were within 1 dilution of AST MIC
# db_filtered %>% mutate(broth_conc=as.numeric(broth_conc)) %>% filter(geno_ecoff=="FN" | geno_ecoff=="FP") %>% mutate(cblo=ecoff_conc/2, cbhi=ecoff_conc*2) %>% filter(broth_conc>=cblo & broth_conc<=cbhi) %>% mutate(amr=as.factor(amr))

# Write database
write.csv(db_filtered, 'results/database.csv', row.names=FALSE)
