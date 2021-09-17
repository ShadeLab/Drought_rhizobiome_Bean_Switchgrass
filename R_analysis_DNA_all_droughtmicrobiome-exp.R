#load libraries
library(ggplot2)
library(phyloseq)
library(dplyr)
library(vegan)
library(reshape2)
library(scales)
library(grid)

otu = read.csv("otu_table.csv", sep=",", row.names=1)
tax = read.csv("taxonomy.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("sample-metadata-bean.csv", sep=",", row.names=1)
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)
##rownames(meta) <- meta$sampleid (this is not needed to merge)
##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged
sample_names(phyloseq_merged)

##check column names of taxonomy table
colnames(tax_table(phyloseq_merged))

##remove mito, chloroplast and archaea, eukarya
phyloseq_merged_clean <- phyloseq_merged %>%
  subset_taxa(
    Domain == "d__Bacteria" &
    Family   != " f__Chloroplast" &
   #Order   != "o_Mitochondria" &
    #Class  != "c__Chloroplast"
   # OTU     != "0a3cf58d4ca062c13d42c9db4ebcbc53" &
   # OTU     != "22f1fa0bdcc19746dee080bcc12a1840"
    Family  != " f__Mitochondria"
)
phyloseq_merged_clean
##one of the chloroplast id'd to OTU 2383b3a3e7fc6bd37f44efe7104b3fe4
#was not removed after this step, others were removed

#phyloseq_merged_clean <- phyloseq_merged_clean1 %>%
  #subset_taxa(
     # Family  != "f__Chloroplast"
  #)


#phyloseq_merged_clean <- phyloseq_merged %>%
  #subset_taxa( Domain == "d__Bacteria" &
              # Family!= "f__Mitochondria" | is.na(Family) & 
              # Family!= "f__Chloroplast" | is.na(Family) )

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(phyloseq_merged_clean))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
View(sample_sum_df)
#write.csv(sample_sum_df, "samplesumdf_controls_mitochlo_removed.csv")
smin <- min(sample_sums(phyloseq_merged_clean))
smean <- mean(sample_sums(phyloseq_merged_clean))
smax <- max(sample_sums(phyloseq_merged_clean))
smin
smean
smax

sum <- sum(sample_sum_df$sum
           )
sum

##20526147 reads to start with (check from table.qzv file interactive plot (minus the 8 positive controls and negative controls)
##19940569 reads remain after filtering out mitochondria and chloroplast, input file is from qiime post dada2 denoise
##10205803 reads total for switchgrass only samples after filtering chloroplast and mitochondria
##9734766 reads total for bean only samples after filtering chloroplast and mitochondria

sample_data(phyloseq_merged_clean)

##rarefaction curve
rarecurve(t(otu_table(phyloseq_merged_clean)), step=50)

##rarefy to 25,000 reads, removes a few samples based on where I set the sample.size to be
phyloseq_rarefy<- rarefy_even_depth(phyloseq_merged_clean, sample.size = 25000, rngseed = TRUE, trimOTUs=FALSE)
##bean 
8 samples removedbecause they contained fewer reads than `sample.size`.
Up to first five removed samples are: 
  
  DNA_129DNA_146DNA_154DNA_155DNA_163


##switchgrass
6 samples removedbecause they contained fewer reads than `sample.size`.
Up to first five removed samples are: 
  
  DNA_116DNA_13DNA_37DNA_38DNA_61

##all experimental
14 samples removedbecause they contained fewer reads than `sample.size`.
Up to first five removed samples are: 
  
  DNA_116DNA_129DNA_13DNA_146DNA_154





##check histogram
sample_sum_df <- data.frame(sum = sample_sums(phyloseq_rarefy))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
View(sample_sum_df)

# melt to long format (for ggploting) 
# prune out phyla below 2% in each sample

phyloseq_class <- phyloseq_rarefy %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class) # Sort data frame alphabetically by phylum/family etc

phyloseq_class$Class[phyloseq_class$Abundance < 0.01] <- "< 1% abund."


# Set colors for plotting
family_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black", "blue", "red", "violet", "yellow", 
  "purple", "pink", "magenta", "darkblue", "darkorchid", "green","darkgreen", "brown","khaki", "grey",
  "cornsilk4","white","cornflowerblue"
)

#OR 
install.packages("Polychrome")
library(Polychrome)
P52<- createPalette(52,  c("#010101", "#ff0000"))
swatch(P52)
View(phyloseq_order)
View(phyloseq_family)
View(phyloseq_class)
# Plot 
plot<- ggplot(phyloseq_class, aes(x= treatment, y = Abundance, fill =Class)) + 
  facet_wrap(drought~planted, scales="free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(P52))+
  #scale_x_discrete(
    #breaks = c("7/8", "8/4", "9/2", "10/6"),
    #labels = c("Jul", "Aug", "Sep", "Oct"), 
    #drop = FALSE
  #) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (class > 1%) \n") +
  ggtitle("Class Composition")
plot
ggsave(filename = "barplot_class_exp_bean_drought-planted.tiff", plot = plot,
        width = 50,
       height = 30, units = c("cm"),
       dpi = 300)
##NMDS

set.seed(1)
# Ordinate
phyloseq_nmds <- ordinate(
  physeq = phyloseq_rarefy, 
  method = "NMDS", 
  distance = "bray"
)
plot1<- plot_ordination(
  physeq = phyloseq_rarefy ,
  ordination = phyloseq_nmds,
  color = "drought",
  shape = "planted"
  
)+scale_shape_manual(values=c(15,19,17,18,20,23,25))+scale_fill_manual(values=c("#a65628", "red", "#ffae19",
                                                                                "#4daf4a", "#1919ff", "darkorchid3", "magenta"))+
  theme_classic()+
  geom_point(aes(color = drought), size = 4) #+stat_ellipse(aes(color = drought, group=drought),type="norm")
plot1
ggsave(filename="NMDS-bean-drought-planted.TIFF", plot=plot1, width=6.8, height=6, units="in", dpi=720)


##adonis
set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_rarefy, method = "bray")
phyloseq_bray
##cant use the sample data from phyloseq merged clean because rarefaction removed 5-6 samples so there
##will be discrepancy in the naming. hence need to remove those samples from sampledf dataframe
sampledf <- data.frame(sample_data(phyloseq_merged_clean))


# Adonis test

##
sample_df <- data.frame(sample_data(phyloseq_rarefy))
View(sample_df)

adonis(phyloseq_bray ~ drought, data = sample_df)

##Bean by drought and planted treatments

Call:
  adonis(formula = phyloseq_bray ~ drought, data = sample_df) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
drought     2    0.6172 0.30860  2.5071 0.04821  0.001 ***
  Residuals  99   12.1858 0.12309         0.95179           
Total     101   12.8030                 1.00000           
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Call:
  adonis(formula = phyloseq_bray ~ planted, data = sample_df) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model    R2 Pr(>F)    
planted     1    2.5734  2.5734  25.157 0.201  0.001 ***
  Residuals 100   10.2296  0.1023         0.799           
Total     101   12.8030                 1.000           
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> 




##Swicthgrass by planted and drought treatments
  Call:
  adonis(formula = phyloseq_bray ~ drought, data = sample_df) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
drought     2    0.3072 0.15361  1.4031 0.02703  0.012 *
  Residuals 101   11.0573 0.10948         0.97297         
Total     103   11.3645                 1.00000         
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Call:
  adonis(formula = phyloseq_bray ~ planted, data = sample_df) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
planted     1    0.4397 0.43970  4.1053 0.03869  0.001 ***
  Residuals 102   10.9248 0.10711         0.96131           
Total     103   11.3645                 1.00000           
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




##switchgrass and bean by crop
Call:
  adonis(formula = phyloseq_bray ~ crop, data = sample_df) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
crop        1    19.424 19.4242  164.04 0.44571  0.001 ***
  Residuals 204    24.156  0.1184         0.55429           
Total     205    43.580                 1.00000           
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
pairwise.adonis2(phyloseq_bray ~ drought, data = sample_df)

##switchgrass (planted vs unplanted)
$parent_call
[1] "phyloseq_bray ~ planted , strata = Null , permutations 999"

$planted_vs_unplanted
Df SumOfSqs      R2      F Pr(>F)    
planted    1   0.4397 0.03869 4.1053  0.001 ***
  Residual 102  10.9248 0.96131                  
Total    103  11.3645 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

attr(,"class")
[1] "pwadstrata" "list"    

##switchgrass (dr vs ww vs pre dr)

$parent_call
[1] "phyloseq_bray ~ drought , strata = Null , permutations 999"

$`pre-drought_vs_well-watered`
Df SumOfSqs      R2      F Pr(>F)  
drought   1   0.1500 0.02315 1.3509  0.051 .
Residual 57   6.3292 0.97685                
Total    58   6.4792 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`pre-drought_vs_drought`
Df SumOfSqs     R2      F Pr(>F)
drought   1   0.1077 0.0184 0.9935  0.436
Residual 53   5.7466 0.9816              
Total    54   5.8544 1.0000              

$`well-watered_vs_drought`
Df SumOfSqs      R2      F Pr(>F)   
drought   1   0.1827 0.01788 1.6746  0.005 **
  Residual 92  10.0387 0.98212                 
Total    93  10.2215 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

attr(,"class")
[1] "pwadstrata" "list"

##bean planted versus unplanted
$parent_call
[1] "phyloseq_bray ~ planted , strata = Null , permutations 999"

$planted_vs_unplanted
Df SumOfSqs    R2      F Pr(>F)    
planted    1   2.5734 0.201 25.157  0.001 ***
  Residual 100  10.2296 0.799                  
Total    101  12.8030 1.000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

attr(,"class")
[1] "pwadstrata" "list"      
> 
##bean (dr vs ww vs predr)
  
  $parent_call
[1] "phyloseq_bray ~ drought , strata = Null , permutations 999"

$`pre-drought_vs_drought`
Df SumOfSqs      R2      F Pr(>F)   
drought   1   0.3132 0.04556 2.4346  0.008 **
  Residual 51   6.5605 0.95444                 
Total    52   6.8737 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`pre-drought_vs_well-watered`
Df SumOfSqs      R2      F Pr(>F)    
drought   1   0.3555 0.05192 3.0667  0.001 ***
  Residual 56   6.4910 0.94808                  
Total    57   6.8464 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`drought_vs_well-watered`
Df SumOfSqs      R2      F Pr(>F)  
drought   1   0.2768 0.02386 2.2247  0.014 *
  Residual 91  11.3202 0.97614                
Total    92  11.5969 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

attr(,"class")
[1] "pwadstrata" "list"      
> 

#pcoa
  
phyloseq_pcoa <- ordinate(
    physeq = phyloseq_rarefy, 
    method = "PCoA", 
    distance = "bray"
  )

# Plot 
plot_ordination(
  physeq = phyloseq_rarefy,
  ordination = phyloseq_pcoa,
  color = "drought",
  shape = "planted",
  title = "PCoA"
) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                "#4daf4a", "#1919ff", "darkorchid3", "magenta", "black")
  ) +
  geom_point(aes(color = drought), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 


