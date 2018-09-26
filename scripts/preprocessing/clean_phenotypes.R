# clean_phenotypes.R 
# Apply ECOFFs and CBPs to raw phenotype data, creating binary categories: 
# Resistant/Non-wild-type = 1, Susceptible/Wild-type = 0.

# Import required functions
source("scripts/functions.R")

# Read the list of sample IDs
samples = read.table("files/sample_list.csv", header=TRUE, sep=",")

# Read the phenotype data
broth = read.table("data/phenotype/broth.csv", header=TRUE, sep=",")
agar = read.table("data/phenotype/agar.csv", header=TRUE, sep=",")
disk = read.table("data/phenotype/disk.csv", header=TRUE, sep=",")

# Read the epidemiological cutoffs (ecoff) and clinical breakpoints (CBPs)
ecoff = read.table("files/epidemiological_cutoffs.csv",header=TRUE, sep=",")
clinbp = read.table("files/clinical_breakpoints.csv", header=TRUE, sep=",")

# Limit broth data to the samples of interest (broth.csv contains control strains)
broth_temp = broth %>% filter(broth$sample %in% samples$sample)
# Increment result of any antimicrobial with max concentration as MIC (SUL, STR and SXT) for accurate assignment as NWT by the apply_mic_breakpoints function.
broth_temp_wsul = broth_temp %>% mutate(SUL = ifelse(X.15 == ">", 257, SUL)) %>% 
  mutate(SXT = ifelse(X.14 == ">", 5, SXT)) %>% mutate(STR = ifelse(X.12 == ">", 65, STR))
# Set broth.csv values to 0 (NWT) and 1 (WT) based on ecoff for each antimicrobial
broth_ecoff = apply_mic_breakpoints(broth_temp_wsul, ecoff, samples)
# Set broth.csv values to 0 (S) and 1 (R) based on CBPs for each antimicrobial
broth_clinbp = apply_mic_breakpoints(broth_temp_wsul, clinbp, samples)

# Set agar.csv values to 0 (S) and 1 (R) based on CBPs for each antimicrobial
agar_temp = agar
agar_clinbp = apply_mic_breakpoints(agar_temp, clinbp, samples)

# Set disk.csv values to 0 (S) and 1 (R) based on CBPs for each antimicrobial
disk_temp = disk
disk_clinbp = apply_zone_breakpoints(disk_temp, clinbp, samples)

# Write data to output files
write.csv(broth_ecoff, file="data/clean/broth_ecoff.csv", row.names=FALSE)
write.csv(broth_clinbp, file="data/clean/broth_clinbp.csv", row.names=FALSE)
write.csv(agar_clinbp, file="data/clean/agar_clinbp.csv", row.names=FALSE)
write.csv(disk_clinbp, file="data/clean/disk_clinbp.csv", row.names=FALSE)