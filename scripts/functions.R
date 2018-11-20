# Install required packages
#install.packages('pacman')
pacman::p_load(tidyr, tibble, stringr, bdpv, scales, ggplot2)

# Verify required directories structure for all script outputs
directories = c("data/clean", "results/", "results/stats", "results/dataframes", "results/tables")
lapply(directories, FUN=dir.create, showWarnings=FALSE)

# Turn off warnings - Many warnings are printed as R performs type conversions automatically
options(warn=-1)
sink(type="message")
# Create all directories required for processing
## preprocessing/phenotype ##

# Creates the broth_ecoff.csv, broth_clinbp.csv and agar_clinbp.csv, using ECOFF or CBP values. Calculates R/NWT > breakpoint.
apply_mic_breakpoints = function(phenotype, breakpoints, sample_list){
  # Create an output dataframe that is updated by the function
  outdf = data.frame(sample_list)
  # For every antimicrobial (3-letter code) in the ECOFF table
  for (antimic in breakpoints$code){
    # Find the index of the current antimicrobial in the phenotype table
    anti_index = match(antimic,names(phenotype)) 
    # Subset the sample and antimicrobial result from the phenotype table to a temporary dataframe
    temp = select(phenotype, sample, anti_index)
    # Set the ecoff value for this antimicrobial as a column in the temporary dataframe
    bpoint_val = filter(breakpoints, code == antimic)
    temp$ecoff = as.numeric(rep(bpoint_val$mic, length(temp)))
    # Determine the NWT/WT status of the samples for this antimicrobial
    temp$result = ifelse(temp[2] > temp$ecoff,1,0)
    # Set the new column name to the antimicrobial (header of temp[2])
    temp[2] = as.numeric(temp$result)
    # Merge result with output dataframe
    outdf = inner_join(outdf, select(temp, sample, names(phenotype)[anti_index]), by="sample")
  }
  return(outdf)
}

# Creates the disk_clinbp.csv using CBP values. Calculates R < zone breakpoint.
apply_zone_breakpoints = function(disk, clinbps, sample_list){
  # Create an output dataframe that is updated by the function
  outdf = data.frame(sample_list)
  # For every antimicrobial (3-letter code) in the CLINBP table
  for (antimic in clinbps$code){
    # Find the index of the current antimicrobial in the disk table
    anti_index = match(antimic,names(disk)) 
    # Subset the sample and antimicrobial result from the disk table to a temporary dataframe
    temp = select(disk, sample, anti_index)
    # Set the cbp value for this antimicrobial as a column in the temporary dataframe
    cbp_val = filter(clinbps, code == antimic)
    temp$cbp = as.numeric(rep(cbp_val$zone, length(temp)))
    # Determine the NWT/WT status of the samples for this antimicrobial
    temp$result = ifelse(temp[2] < temp$cbp,1,0)
    # Set the new column name to the antimicrobial (header of temp[2])
    temp[2] = as.numeric(temp$result)
    # Merge result with output dataframe
    outdf = inner_join(outdf, select(temp, sample, names(disk)[anti_index]), by="sample")
  }
  return(outdf)
}


## preprocessing/genotype ##

# Retrieve columns from ARIBA genotype dataframe containing the string 'match'.
matchcols <- function(df) {
  return(df[,grep('match',colnames(df))])
}

# Clean ARIBA sample name strings
clnsample <- function(x) {
  return (x %>% str_replace("out.run","") %>% str_replace("/.*",""))
}

# Get 'match' columns from ARIBA summary dataframe
clndata <- function(x) {
  # Create a dataframe with samples in the first column and ARIBA matches in the remainder
  cmatch = x %>% matchcols() %>% rownames_to_column %>% mutate(rowname = clnsample(rowname))
  colnames(cmatch)[1] = "sample"
  # Recode matches as 1, non-matches as 0 
  card_cln = cmatch %>%  lapply(function(y) recode(y, no=0, yes=1)) %>% 
    as.data.frame() %>% mutate(sample=cmatch$sample)
  return(card_cln)
}

# Merge columns in a binary dataframe that match a given string (contstr)
mergebin <- function(df, contstr, outnm) {
  # Create an output dataframe fro editing
  outdf = df
  # Create a regular expression from all strings in `constr`, seperated by the OR operator (|)
  cstr = paste(contstr,collapse="|")
  # Creating a new column merging columns matching constr by summing across the rows
  mid = select(df,matches(cstr)) %>% rowSums()
  # Normalise new column to binary values 
  mid[mid>0] = 1
  # Add the new column named by the string stored in `outnm`
  outdf[, outnm] = mid
  return(outdf)
}

# Merge streptomycin columns, producing '1' for strA and strB matches only
mergestr <- function(df, contstr, outnm) {
  # Create an output dataframe fro editing
  stroutdf = df
  # Create a regular expression from all strings in `constr`, seperated by the OR operator (|)
  cstr = paste(contstr,collapse="|")
  # Creating a new column merging columns matching constr by summing across the rows
  mid = select(df,matches(cstr)) %>% rowSums()
  # Normalise new column to binary values 
  mid[mid<2] = 0
  mid[mid>0] = 1
  # Add the new column named by the string stored in `outnm`
  stroutdf[, outnm] = mid
  return(stroutdf)
}

# Merge binary columns within dataframe (-sample) that match contstr
binarise = function(df) {select(df, -sample) %>% rowSums %>% sapply(function(x) if_else(x > 0, 1, 0))}

# Function: Return DF of all genotype columns for each gene merged by the schema in mapping dataset 
# Input : Genotype dataframe [sample, gene_name... ]; Mapping data set [ gene_name, amrs ]
map_amrs = function(df, amr_map, fname){
  # set colnames of amr_map for use in this function
  colnames(amr_map) = c("gene","amrs")
  # Add sample IDs to a dataframe for joining
  outdf = as.data.frame(df$sample) %>% `colnames<-`(c("sample"))
  # Extract the set of amr_map AMRs for lookup
  amr_set = amr_map$amr %>% stringi::stri_remove_empty() %>% unique()
  
  
  # Loop through the input gene set
  for (amr in amr_set) {
    ## Filter amr_map genotypes with am_
    genecols = filter(amr_map, amrs==amr)$gene %>% as.character
    ## Select this list of genes from the resdf
    tobin = df %>% select(sample, matches(paste(genecols, collapse="|"))) 
    ## Merge binary values
    tomerge = tobin %>% mutate(amerge= binarise(.))
    colnames(tomerge)[which(colnames(tomerge)=="amerge")] = amr
    ## Write out for future audit
    write.csv(tomerge, paste('results/dataframes/genotype_',fname,'_',amr,'.csv', sep=""), row.names =FALSE)
    ## Add to a dataframe
    outdf = full_join(outdf, select(tomerge, sample, amr), by="sample")
  }
return(outdf)
}

## ANALYSIS ##

# Generate table of TP,FP,FN,TN counts
mframe = function(truthset, testset, cname){
  # rename results column of testset before merge
  names(testset)[which(names(brothe) == "result")] = "result_test"
  # merge truthset (broth) and test data (phenotype/genotype test) by sample and amr
  comp_df = left_join(truthset, testset, by=c("sample","amr")) %>% replace_na(list(result_test = 0)) %>%
    # calculate matches for statistics
    mutate( !!cname := case_when(result == result_test & result == 0 ~ "TN", 
                                 result == result_test & result == 1  ~ "TP", 
                                 result != result_test & result == 0 ~ "FP", 
                                 result != result_test & result == 1  ~ "FN"))
  return(select(comp_df, -contains("result")))
}

# Read input binary CSV in long format (each row is a sample, AMR and its result)
readcsv_gather = function(infile){ return(read.csv(infile) %>% gather(amr,result, -sample)) }

# Get statistics
get_stats = function(mframe) {
  
  df = mframe %>% select(amr, result) %>% count(amr, result) %>% spread(result,n) %>% 
    mutate_all(funs(replace(., is.na(.), 0))) %>% mutate(POS = TP + FN) %>% mutate(NEG = TN + FP) %>%
    mutate(PVPOS = TP + FP) %>% mutate(PVNEG = TN + FN) %>%
    bind_rows(summarise_all(., funs(if(is.numeric(.)) sum(.) else "Total")))
  # Account for gentamicin 0 pos causing binom.test error in PVPOS comparisons (replace PVPOS 0 with 1)
  pvpos_zero = which(df$PVPOS==0)
  df$PVPOS = replace(df$PVPOS, pvpos_zero, 1)
  
  outdf = df %>% rowwise() %>%
    mutate(sens = binom.test(TP,POS,conf.level=0.95)[[5]]) %>%
    mutate(sens_l = binom.test(TP,POS,conf.level=0.95)[[4]][1]) %>%
    mutate(sens_u = binom.test(TP,POS,conf.level=0.95)[[4]][2]) %>%
    mutate(sens_pval = binom.test(TP,POS,conf.level=0.95)[[3]]) %>%
    mutate(spec = binom.test(TN,NEG,conf.level=0.95)[[5]]) %>%
    mutate(spec_l = binom.test(TN,NEG,conf.level=0.95)[[4]][1]) %>%
    mutate(spec_u = binom.test(TN,NEG,conf.level=0.95)[[4]][2]) %>%
    mutate(spec_pval = binom.test(TN,NEG,conf.level=0.95)[[3]]) %>%
    mutate(ppv = binom.test(TP,PVPOS,conf.level=0.95)[[5]]) %>%
    mutate(ppv_l = binom.test(TP,PVPOS,conf.level=0.95)[[4]][1]) %>% # PROBLEM
    mutate(ppv_u = binom.test(TP,PVPOS,conf.level=0.95)[[4]][2]) %>% # PROBLEM
    mutate(ppv_pval = binom.test(TP,PVPOS,conf.level=0.95)[[3]]) %>% # PROBLEM
    mutate(npv = binom.test(TN,PVNEG,conf.level=0.95)[[5]]) %>%
    mutate(npv_l = binom.test(TN,PVNEG,conf.level=0.95)[[4]][1]) %>%
    mutate(npv_u = binom.test(TN,PVNEG,conf.level=0.95)[[4]][2]) %>%
    mutate(npv_pval = binom.test(TN,PVNEG,conf.level=0.95)[[3]])
  
  return(outdf)
}

## Functions for create_db.R
# Get genotype string from 'binarised' ariba genotype data. for each row, filter to cols == 1 and return 'sample, pasted colnames'. 
get_genotype_cln_list = function(inlist){
  out_string = inlist[inlist=="1"] %>% names %>% paste(collapse=",") %>% stringi::stri_replace_all("", fixed="match") %>% stringi::stri_replace_all("", fixed=".") %>% stringi::stri_replace_all("", fixed="_")
  return(out_string)
}

# Merge genotype string matching value
merge_genotype_string = function(input, value){
  outstring = names(input)[input==value] %>% na.omit() %>% paste(collapse=",") %>%
    stringi::stri_replace_all("", fixed="match") %>% stringi::stri_replace_all("", fixed=".") %>% 
    stringi::stri_replace_all("", fixed="_") %>% stringi::stri_replace_all("", fixed="assembled")
  return(outstring)
}

# Wrapper function for merging genotype strings across rows for each sample
merge_geno_wrapper = function(indf, value, colname){
  cname = colname
  output = indf %>% rowwise %>% 
    do({ result = as_data_frame(.); result$s = merge_genotype_string(result, value); result}) %>% 
    select(sample, !!cname:=s)
  return(output)
}

# Function for returning AMR symbols from agar and broth dilution data
phenotype_symbol = function(df){
  # Length of df without initial 'sample' column
  df_len = length(df)
  # Create dataframe containing AMR symbol columns (stepwise in 2 from column 2)
  subset = select(df, sample, seq(2,df_len,2))
  # Rename with AMR genes (stepwise in 2 from column 3). 
  # names(x)[-1] excludes the intial 'sample' column.
  names(subset)[-1] = names(df)[seq(3,df_len,2)]
  # Return in longformat for merging with database
  output_symbol_df = gather(subset, amr, symbol, 2:length(subset))
  return(output_symbol_df)
}

# Merge values from bdpv outputs
merger = function(main,low,up){
  mainr = round(main,2)
  lowr = round(low,2)
  upr = round(up,2)
  output = as.character(c(mainr," (",lowr,"-",upr,")"))
  #return(paste(output,sep='',collapse=''))
  return(paste(mainr,' (',lowr, '-', upr,')',sep=''))
}