# .PHONY target. Runs to make the files described below
all : output
	
clean_phenotype_prereqs = files/sample_list.csv data/phenotype/broth.csv data/phenotype/agar.csv data/phenotype/disk.csv files/epidemiological_cutoffs.csv
clean_genotype_prereqs = data/genotype/ariba_resfinder.csv data/genotype/pointfinder_pointmutations.csv files/genotype_conversion.csv

# Clean phenotypes
clean_phenotype : $(clean_phenotype_prereqs)
	Rscript scripts/preprocessing/clean_phenotypes.R
	
# Clean genotypes
clean_genotype : $(clean_genotype_prereqs)
	Rscript scripts/preprocessing/clean_genotypes.R
	
# Analysis
output : clean_phenotype clean_genotype
	Rscript scripts/analysis.R
	Rscript scripts/create_db.R
	Rscript scripts/tables.R