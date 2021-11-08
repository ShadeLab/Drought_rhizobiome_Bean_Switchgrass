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
metadata = read.csv("sample-metadata-decontam.csv", sep=",", row.names=1) ##contains true samples & negative controls only
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged



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



