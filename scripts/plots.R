# plots.R
# Generate figures for manuscript. 
#     Inputs: None
#     Outputs: 
#         plots/gene_freqs.png: A bar graph of gene frequencies

# Import functions and libraries
source('scripts/functions.R')

# Import plot data, getting only sample and reference sequences where ARIBA match column was "yes"
geno_counts = read.csv('results/tables/genotype_counts_full.csv', colClasses=c("NULL",NA,NA,"NULL")) %>% drop_na()

# Clean reference sequences to produce column containing gene groups. Used as plot legend to colour data.
df_groups = geno_counts %>% mutate(new = str_replace(ref,"_.+\\..*","")) %>% mutate(new = str_replace(new,"\\..*","")) %>% mutate(new= str_replace(new," \\(.*","")) %>% mutate(new=str_replace(new,"aadA.*","aadA")) %>% mutate(new=str_replace(new,"dfrA.*","dfrA")) %>% mutate(new=str_replace(new,"sul.*","sul"))
df_groups$new <- factor(df_groups$new)

# Generate plot
out_plot = ggplot(df_groups, aes(x=reorder(ref, desc(ref)), y=yes_counts, fill=new)) + 
  geom_bar(stat='identity') + 
  theme(axis.text.y=element_text(size=8), legend.key.size=unit(0.5,"line")) +  # Edit element sizes
  labs(x="Reference sequence assembled by ARIBA", y="Frequency of isolates (n=102)", fill="Gene") + 
  theme(panel.grid.minor = element_line(colour="white")) +
  ylim(0,40)  + coord_flip() # + geom_text(aes(label=yes_counts), size=3)

# Save plots 
ggsave('plots/figure1.png', plot=out_plot, dpi=500, limitsize=TRUE)
