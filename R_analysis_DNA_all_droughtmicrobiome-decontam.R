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
sample_sum_df <- data.frame(sum = sample_sums(phyloseq_merged_clean_decontam_final))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
View(sample_sum_df)
#write.csv(sample_sum_df, "samplesumdf_controls_mitochlo_removed.csv")
smin <- min(sample_sums(phyloseq_merged_clean_decontam_final))
smean <- mean(sample_sums(phyloseq_merged_clean_decontam_final))
smax <- max(sample_sums(phyloseq_merged_clean_decontam_final))
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

##19938535 reads total for bean and SW samples after removing mito, chloro, and contaminant OTUs and neg control samples from phyloseq object
#versus 19940569 reads that remained before removing contaminant OTUs
#10204741 reads for SW ater removing contaminant reads, mito, chloro, neg con
##9733794 reads for bean after removing contaminant reads, mito, chloro, neg con
##4139117 reads from bean planted samples minus pre drought and controls (n=50)-- time and treatment interaction tests
##4794706 reads from bean unplanted samples minus pre drought and controls (n=50)-- time and treatment interaction tests

sample_data(phyloseq_merged_clean_decontam_final)

##rarefaction curve
rarecurve(t(otu_table(phyloseq_merged_clean_decontam_final)), step=50)

##rarefy to 25,000 reads, removes a few samples based on where I set the sample.size to be
phyloseq_rarefy<- rarefy_even_depth(phyloseq_merged_clean_decontam_final, sample.size = 25000, rngseed = TRUE, trimOTUs=FALSE)
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
ggsave(filename = "barplot_class_decontam-bean-planted-drought.tiff", plot = plot,
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
  geom_point(aes(color = drought), size = 4) #+stat_ellipse(aes(color = crop, group=crop),type="norm")
plot1
ggsave(filename="NMDS-decontam-bean-drought-planted.TIFF", plot=plot1, width=6.8, height=6, units="in", dpi=720)


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

adonis(phyloseq_bray ~ planted, data = sample_df)

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



-------
##decontam
library(phyloseq)
library(ggplot2)
library(tidyverse)
#install.packages("vctrs")
#update.packages("tidyverse")
#update.packages("scales")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

library(decontam)

otu = read.csv("otu_table.csv", sep=",", row.names=1)
tax = read.csv("taxonomy.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("sample-metadata-decontam.csv", sep=",", row.names=1)
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

head(sample_data(phyloseq_merged_clean))

##inspect library sizes
df <- as.data.frame(sample_data(phyloseq_merged_clean)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phyloseq_merged_clean)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

##identify contaminants by prevalence. in this method, the distribution 
#of the frequency of each sequence feature as a function of the 
#prevalence is used to identify contaminants. In this method, 
#the prevalence (presence/absence across samples) of each sequence feature 
#in true positive samples is compared to the prevalence in negative controls 
#to identify contaminants.
sample_data(phyloseq_merged_clean)$is.neg <- sample_data(phyloseq_merged_clean)$Sample_or_Control == "Control Sample"
contamdf.prev0.5 <- isContaminant(phyloseq_merged_clean, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev0.5$contaminant)
FALSE  TRUE 
78160    38 


# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(phyloseq_merged_clean, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev0.5$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
write.csv(df.pa, "contaminant-table.csv")

write.csv(contamdf.prev0.5, "contaminant-prev-0.5.csv")

##removing contaminants from phyloseq object
df.pa <- read.csv("contaminant-table.csv")
View(df.pa)
subset.df <- subset(df.pa, contaminant== "FALSE")
View(subset.df)
keep.taxa <- as.vector(subset.df$X)
keep.taxa

phyloseq_merged_clean_decontam <- prune_taxa(keep.taxa, phyloseq_merged_clean)
phyloseq_merged_clean_decontam

phyloseq_merged_clean_decontam
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 78160 taxa and 268 samples ]
sample_data() Sample Data:       [ 268 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 78160 taxa by 7 taxonomic ranks ]

##filter metadata to true samples only then merge with phyloseq decontam (with removed OTU contaminants)
metadata<- read.csv("sample-metadata-decontam.csv")
subset.metadata<- subset(metadata, Sample_or_Control=="True Sample")
keep.samples <- as.vector(subset.metadata$sample_id)
keep.samples

phyloseq_merged_clean_decontam_final <- prune_samples(keep.samples, phyloseq_merged_clean_decontam)
phyloseq_merged_clean_decontam_final

phyloseq_merged_clean_decontam_final
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 78160 taxa and 220 samples ]
sample_data() Sample Data:       [ 220 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 78160 taxa by 7 taxonomic ranks ]
##use this to check workflow and stats
sample_sum_df <- data.frame(sum = sample_sums(phyloseq_merged_clean_decontam_final))
sum <- sum(sample_sum_df$sum
)
sum

##rarefy to 25,000 reads, removes a few samples based on where I set the sample.size to be
phyloseq_rarefy<- rarefy_even_depth(phyloseq_merged_clean_decontam_final, sample.size = 25000, rngseed = TRUE, trimOTUs=FALSE)


##differences by crop after decontam
set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_rarefy, method = "bray")
phyloseq_bray

# Adonis test
##cant use the sample data from phyloseq merged clean because rarefaction removed 5-6 samples so there
##will be discrepancy in the naming. hence need to remove those samples from sampledf dataframe using phyloseq rarefy object as input

sample_df <- data.frame(sample_data(phyloseq_rarefy))
View(sample_df)

adonis(phyloseq_bray ~ crop, data = sample_df)

Call:
  adonis(formula = phyloseq_bray ~ crop, data = sample_df) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
crop        1    19.409 19.4088  163.98 0.44562  0.001 ***
  Residuals 204    24.146  0.1184         0.55438           
Total     205    43.555                 1.00000           
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##prune to SW samples only

metadata<- read.csv("sample-metadata-switchgrass.csv")
#subset.metadata<- subset(metadata, crop=="switchgrass")
keep.samples <- as.vector(metadata$sample_id)
keep.samples

phyloseq_merged_clean_decontam_final <- prune_samples(keep.samples, phyloseq_merged_clean_decontam)
phyloseq_merged_clean_decontam_final


##results SW -decontam

phyloseq_merged_clean_decontam_final
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 78160 taxa and 110 samples ]
sample_data() Sample Data:       [ 110 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 78160 taxa by 7 taxonomic ranks ]

##do the same steps as above to check reads and rarefy, then generate BC distance object
adonis(phyloseq_bray ~ drought, data = sample_df)

Call:
  adonis(formula = phyloseq_bray ~ drought, data = sample_df) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
drought     2    0.3061 0.15307  1.3976 0.02693  0.012 *
  Residuals 101   11.0617 0.10952         0.97307         
Total     103   11.3678                 1.00000         
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> adonis(phyloseq_bray ~ planted, data = sample_df)

Call:
  adonis(formula = phyloseq_bray ~ planted, data = sample_df) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
planted     1     0.444 0.44395  4.1453 0.03905  0.001 ***
  Residuals 102    10.924 0.10710         0.96095           
Total     103    11.368                 1.00000           
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> library(pairwiseAdonis)
Loading required package: cluster
> pairwise.adonis2(phyloseq_bray ~ drought, data = sample_df)
$parent_call
[1] "phyloseq_bray ~ drought , strata = Null , permutations 999"

$`pre-drought_vs_well-watered`
Df SumOfSqs      R2      F Pr(>F)  
drought   1   0.1488 0.02294 1.3385  0.043 *
  Residual 57   6.3364 0.97706                
Total    58   6.4852 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`pre-drought_vs_drought`
Df SumOfSqs      R2      F Pr(>F)
drought   1   0.1074 0.01834 0.9903  0.443
Residual 53   5.7460 0.98166              
Total    54   5.8534 1.00000              

$`well-watered_vs_drought`
Df SumOfSqs      R2      F Pr(>F)   
drought   1   0.1825 0.01785 1.6721  0.005 **
  Residual 92  10.0409 0.98215                 
Total    93  10.2234 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

attr(,"class")
[1] "pwadstrata" "list"  


##prune to bean samples only

metadata<- read.csv("sample-metadata-bean.csv")
#subset.metadata<- subset(metadata, crop=="bean")
keep.samples <- as.vector(metadata$sample_id)
keep.samples

phyloseq_merged_clean_decontam_final <- prune_samples(keep.samples, phyloseq_merged_clean_decontam)
phyloseq_merged_clean_decontam_final

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 78160 taxa and 110 samples ]
sample_data() Sample Data:       [ 110 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 78160 taxa by 7 taxonomic ranks ]

##do the same steps as above to check reads and rarefy, then generate BC distance object

Call:
  adonis(formula = phyloseq_bray ~ drought, data = sample_df) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
drought     2     0.619 0.30948  2.5141 0.04834  0.001 ***
  Residuals  99    12.186 0.12309         0.95166           
Total     101    12.805                 1.00000           
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> adonis(phyloseq_bray ~ planted, data = sample_df)

Call:
  adonis(formula = phyloseq_bray ~ planted, data = sample_df) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
planted     1     2.576 2.57603  25.183 0.20117  0.001 ***
  Residuals 100    10.229 0.10229         0.79883           
Total     101    12.805                 1.00000           
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> pairwise.adonis2(phyloseq_bray ~ drought, data = sample_df)
$parent_call
[1] "phyloseq_bray ~ drought , strata = Null , permutations 999"

$`pre-drought_vs_drought`
Df SumOfSqs     R2      F Pr(>F)  
drought   1   0.3138 0.0457 2.4425  0.011 *
  Residual 51   6.5513 0.9543                
Total    52   6.8650 1.0000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`pre-drought_vs_well-watered`
Df SumOfSqs      R2      F Pr(>F)    
drought   1   0.3554 0.05187 3.0635  0.001 ***
  Residual 56   6.4976 0.94813                  
Total    57   6.8530 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`drought_vs_well-watered`
Df SumOfSqs      R2     F Pr(>F)   
drought   1   0.2784 0.02399 2.237  0.009 **
  Residual 91  11.3239 0.97601                
Total    92  11.6023 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

attr(,"class")
[1] "pwadstrata" "list"      


-----------------
  
##bean samples mantel test correlation: weighted unifrac and BC
  
set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_rarefy
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 78160 taxa and 102 samples ]
sample_data() Sample Data:       [ 102 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 78160 taxa by 7 taxonomic ranks ]
phyloseq_bray <- phyloseq::distance(phyloseq_rarefy, method = "bray")
phyloseq_bray

library("ape")
random_tree = rtree(ntaxa(phyloseq_rarefy), rooted=TRUE, tip.label=taxa_names(phyloseq_rarefy))
plot(random_tree)
phyloseq_tree = merge_phyloseq(phyloseq_rarefy, random_tree)
phyloseq_tree

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 78160 taxa and 102 samples ]
sample_data() Sample Data:       [ 102 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 78160 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 78160 tips and 78159 internal nodes ]


#calculate weighted unifrac distance matrix
phyloseq_wunifrac <- phyloseq::distance(phyloseq_tree, method = "wunifrac" )

##mantel tests to test correlation between wunifrac and BC (https://jkzorz.github.io/2019/07/08/mantel-test.html)
library(scales)
library(reshape2)
library(grid)
set.seed(1)
mantel_cor = mantel(phyloseq_bray, phyloseq_wunifrac, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_cor
class(mantel_cor)
mantel <- as.matrix(mantel_cor)
mantel
mantel.tibble<-mantel %>% 
  as_tibble()
mantel.df <- as.data.frame(mantel)

Mantel statistic based on Spearmans rank correlation rho

Call:
mantel(xdis = phyloseq_bray, ydis = phyloseq_wunifrac, method = "spearman",      permutations = 9999, na.rm = TRUE) 

Mantel statistic r: 0.9549 
      Significance: 1e-04 

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0583 0.0756 0.0912 0.1086 
Permutation: free
Number of permutations: 9999

##if used set.seed(1)
##result is only slightly different
Mantel statistic based on Spearmans rank correlation rho

Call:
mantel(xdis = phyloseq_bray, ydis = phyloseq_wunifrac, method = "spearman",      permutations = 9999, na.rm = TRUE) 

Mantel statistic r: 0.9435 
      Significance: 1e-04 

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0614 0.0804 0.0964 0.1169 
Permutation: free
Number of permutations: 9999




plot(mantel_cor$xdis,mantel_cor$ydis)

##mantel.test(phyloseq_bray, phyloseq_wunifrac, plot = TRUE, method = c("spearman"))

install.packages("ecodist")
library(ecodist)
set.seed(1)
mantel(phyloseq_bray~phyloseq_wunifrac, nperm = 1000,
       mrank = TRUE, nboot = 500, pboot = 0.9, cboot = 0.95)

mantelr      pval1      pval2      pval3      llim.2.5%    ulim.97.5% 
  0.9548856  0.0010000  1.0000000  0.0010000  0.9486073     0.9636574 
> 
  
#if used set.seed(1)
mantel(phyloseq_bray~phyloseq_wunifrac, nperm = 1000,
         +        mrank = TRUE, nboot = 500, pboot = 0.9, cboot = 0.95)
mantelr      pval1      pval2      pval3  llim.2.5% ulim.97.5% 
  0.9435211  0.0010000  1.0000000  0.0010000  0.9366578  0.9540723 

#plot(mgram(phyloseq_bray,phyloseq_wunifrac))
ggplot_mantel<- function(mant){

  df <- data.frame(x = mant$xdis,
                   y = mant$ydis)
ggplot(df, aes(x,y)) + 
    geom_smooth(method = "lm", alpha = 0.2, colour = "black") + 
    geom_point(aes(colour = black), size = 4) +
    labs(y = "weighted unifrac", x = "bray curtis similarity") + 
    theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
        axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
        axis.title= element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        legend.title = element_text(size =12, face = "bold", colour = "black"),
        legend.text = element_text(size = 10, face = "bold", colour = "black")) +
    scale_colour_continuous(high = "navy", low = "salmon")

}

ggplot_mantel(mantel_cor)
plot(mantel_cor)

install.packages("mpmcorrelogram")
library(mpmcorrelogram)
mpmcorrelogram(phyloseq_bray, phyloseq_wunifrac, method = "spearman",
               simil = TRUE,
               plot = TRUE, print = TRUE)
##plots well
plot<- plot(phyloseq_bray,phyloseq_wunifrac,pch=16,cex=0.5,col="black",bty="l")
ggsave(filename = "correlation-bray-wunifrac.tiff", plot = plot,
       width = 50,
       height = 30, units = c("cm"),
       dpi = 300)


--------------------
  
##bean planted time versus treatment interaction (using permanova)


##prune to bean planted samples only
  
metadata<- read.csv("bean-time-treat-nopredr.csv")
subset.metadata<- subset(metadata, planted=="planted")
keep.samples <- as.vector(subset.metadata$sample_id)
keep.samples

phyloseq_merged_clean_decontam_final <- prune_samples(keep.samples, phyloseq_merged_clean_decontam)
phyloseq_merged_clean_decontam_final

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 78160 taxa and 50 samples ]
sample_data() Sample Data:       [ 50 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 78160 taxa by 7 taxonomic ranks ]


##do the same steps as above to check reads and rarefy, then generate BC distance object
## rarefaction removed 5 samples DNA_154DNA_155DNA_163DNA_164DNA_245
phyloseq_rarefy<- rarefy_even_depth(phyloseq_merged_clean_decontam_final, sample.size = 25000, rngseed = TRUE, trimOTUs=FALSE)
set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_rarefy, method = "bray")
phyloseq_bray
##cant use the sample data from phyloseq merged clean because rarefaction removed 5-6 samples so there
##will be discrepancy in the naming. hence need to remove those samples from sampledf dataframe by using phyloseq rarefy object as input

# Adonis test
sample_df <- data.frame(sample_data(phyloseq_rarefy))
View(sample_df)

adonis(phyloseq_bray ~ drought*harvest, by="margin", data = sample_df) ##time in days (as numeric) gave error # Error in model.frame.default(formula, data, drop.unused.levels = TRUE) : 
##invalid type (closure) for variable 'time', hence did with harvest as categorical variable
##planted tests results bean 
Call:
  adonis(formula = phyloseq_bray ~ drought * harvest, data = sample_df,      by = "margin") 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                 Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
drought          1    0.3348 0.33485  3.5828 0.07144  0.001 ***
harvest          4    0.5783 0.14457  1.5469 0.12338  0.001 ***
drought:harvest  4    0.5027 0.12568  1.3448 0.10726  0.010 ** 
Residuals       35    3.2711 0.09346         0.69791           
Total           44    4.6870                 1.00000           
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



## rarefaction removed 2 samples removedbecause they contained fewer reads than `sample.size`.
#removed samples are: DNA_146DNA_168


##unplanted test results bean
Call:
  adonis(formula = phyloseq_bray ~ drought * harvest, data = sample_df,      by = "margin") 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                 Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
drought          1    0.1329 0.132923  1.5264 0.02999  0.004 ** 
harvest          4    0.5716 0.142903  1.6410 0.12895  0.001 ***
drought:harvest  4    0.4190 0.104759  1.2030 0.09453  0.008 ** 
Residuals       38    3.3091 0.087082         0.74653           
Total           47    4.4327                  1.00000           
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##linked resource webpage: https://stat.ethz.ch/pipermail/r-sig-ecology/2018-November/005830.html

------------
##diagram for BC similarity change by time, bean drought (planted and unplanted)
metadata<- read.csv("bean-ww-planted-unplanted.csv")
#subset.metadata<- subset(metadata, planted=="planted")
keep.samples <- as.vector(metadata$sample_id)
keep.samples

phyloseq_merged_clean_decontam_final <- prune_samples(keep.samples, phyloseq_merged_clean_decontam)
phyloseq_merged_clean_decontam_final ##same sample numbers for drought and well watered == 60

otu_table()   OTU Table:         [ 78160 taxa and 60 samples ]
sample_data() Sample Data:       [ 60 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 78160 taxa by 7 taxonomic ranks ]


##do the same steps as above to check reads and rarefy, then generate BC distance object

phyloseq_rarefy<- rarefy_even_depth(phyloseq_merged_clean_decontam_final, sample.size = 25000, rngseed = TRUE, trimOTUs=FALSE)
## rarefaction removed 7 samples for drought dataset because they contained fewer reads than `sample.size`.
#Up to first five removed samples are: DNA_129DNA_146DNA_154DNA_155DNA_163, 164, 168


##rarefaction removed 2 samples from well watered dataset 
##2 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: DNA_129DNA_245



sample_data(phyloseq_rarefy)
phyloseq_rarefy ##drought
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 78160 taxa and 53 samples ]
sample_data() Sample Data:       [ 53 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 78160 taxa by 7 taxonomic ranks ]


phyloseq_rarefy ##wellwatered
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 78160 taxa and 58 samples ]
sample_data() Sample Data:       [ 58 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 78160 taxa by 7 taxonomic ranks ]



set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_rarefy, method = "bray")
phyloseq_bray

install.packages("usedist")
library(usedist)## can be useful in this case for functions like dist_subset and dist_groups

install.packages("dendextend")
library(dendextend)
dist_long<-dist_long(phyloseq_bray)
dist_long
write.csv(dist_long, "dist_long_wellwatered.csv")
dist_long<-read.csv("dist_long_wellwatered.csv")##change dist_long column names to iso1 and iso2 for formatting
metadata_new<-data.frame(sample_data(phyloseq_rarefy))
View(metadata_new)
write.csv(metadata_new, "metadata_dist_long_wellwatered.csv") 
metadata_new<- read.csv("metadata_dist_long_wellwatered_edited.csv")#change column names in metadata to group and sample-id
install.packages("harrietr", dependencies = TRUE)
install.packages('devtools')
setRepositories(ind = 1:2)
devtools::install_github("andersgs/harrietr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

library(harrietr)
library(ggtree)
update.packages("ggtree")
dist.df<- join_metadata(dist_long, metadata_new, isolate = 'sample_id', group = 'group', remove_ind=FALSE)
write.csv(dist.df,"dist_df_wellwatered.csv")

##bean planted BC versus time 

distdfbean<- read.csv("dist_df_wellwatered.csv")
distdfbean.planted <- distdfbean  %>%
  filter(planted_x=="planted"&
          planted_y== "planted")

View(distdfbean.planted)

write.csv(distdfbean.planted, "distdfbean_planted_wellwatered.csv")
beanpredr.planted<- distdfbean.planted  %>%
  filter(group_2=="pre-drought_planted" &
          group_1!="pre-drought_planted")

View(beanpredr.planted)
write.csv(beanpredr.planted, "beanpredr_planted_wellwatered.csv")
ggplot(beanpredr.planted, aes(x= group,y = distance)) + 
  geom_point(aes(), size = 4) +
  labs(y = "BC sim", x = "group") + 
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#bean unplanted BC versus time

distdfbean<- read.csv("dist_df_wellwatered.csv")
distdfbean.unplanted <- distdfbean  %>%
  filter(planted_x=="unplanted"&
           planted_y== "unplanted")

View(distdfbean.unplanted)

write.csv(distdfbean.unplanted, "distdfbean_unplanted_wellwatered.csv")
beanpredr.unplanted<- distdfbean.unplanted  %>%
  filter(group_2=="pre-drought_unplanted" &
           group_1!="pre-drought_unplanted")

View(beanpredr.unplanted)
write.csv(beanpredr.unplanted, "beanpredr_unplanted_wellwatered.csv")
ggplot(beanpredr.unplanted, aes(x= group,y = distance)) + 
  geom_point(aes(), size = 4) +
  labs(y = "BC sim", x = "group") + 
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  

##plottting BC sim with time, planted unplanted, drought in bean
bean_dr_BC_time<- read.csv("bean-wellwatered-BC-time-planted-unplanted.csv")
plot<- ggplot(bean_dr_BC_time, aes(x= group,y=distance, color=planted_x)) + 
  geom_point(aes()) + geom_smooth(method = lm, se=FALSE) +
  labs(y = "BC similarity", x = "time(days)") + theme(axis.text.x = element_text(face = "bold",colour = "black", size = 12,angle = 90, vjust = 0.5, hjust=1))+
  theme( axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =12, face = "bold", colour = "black"),
         legend.text = element_text(size = 10, face = "bold", colour = "black")) 
  #scale_colour_continuous(high = "navy", low = "salmon")
  
  #scale_colour_continuous(high = "navy", low = "salmon")
plot
ggsave(filename = "bean_ww_BC_time.tiff", plot = plot,
       width = 18,
       height = 18, units = c("cm"),
       dpi = 300)

