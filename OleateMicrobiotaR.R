## reading Mothur files into Phyloseq
## install phyloseq this way
install.packages("devtools")
install.packages(Rtools.33)
source("https://bioconductor.org/biocLite.R")
biocLite("BiocInstaller")
library("devtools")
install_github("phyloseq", "joey711")
install_github("guokai8/microbial") # Install the package
library(microbial)
## others can be installed using tools -> install packages on R studio
## If run into errors such as below the try this command:
## setInternet2(TRUE) 
## "Cannot open compressed file 'knitr/DESCRIPTION', probable reason 'No such file or directory'"
library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library("grid")
library("directlabels")
library("knitr")
library("clustsig")
library("ellipse")
## import data 
sharedFile = read.table("oleate.opti_mcc.shared") ## this takes a while (~20mins on big computer)
sharedFile = t(sharedFile) ## transform data
rownames(sharedFile) = sharedFile[,1]  ## define rowname
colnames(sharedFile) = sharedFile[2,]  ## define column names
sharedFile
sharedFile <-as.matrix(sharedFile)
sharedFile = sharedFile[,2:57]       ## say which columns you want (remember order [row,colum]) 2 -> 57 for me
## as I had 57 samples
sharedFile = sharedFile[4:819,]     ## and which rows you want. I have 819 OTUs so I changed this number to 819
# original => sharedFile = sharedFile[4:819,]
class(sharedFile) <- "numeric"
head(sharedFile)
dim(sharedFile)
## look at head first 7 lines to see it's made correclty
## Import subsampled otu matrix (26880 seqs)
sharedsubFile = read.table('oleate.opti_mcc.0.03.subsample.shared')
sharedsubFile = t(sharedsubFile)
rownames(sharedsubFile) = sharedsubFile[,1]
colnames(sharedsubFile) = sharedsubFile[2,]
sharedsubFile = sharedsubFile[,2:57]
sharedsubFile = sharedsubFile[4:494,]
class(sharedsubFile) <- "numeric"
head(sharedsubFile)
dim(sharedsubFile)
sharedsubFile <-as.matrix(sharedsubFile)
## Import taxonomy file 
## As it is from mothur, there are not column headers for order, family genus etc
## in the cons.taxonomy file and this must be fixed first. This is how I did it:
## read cons.taxonomy file into excel and choose semicolon as a separator so each tax level
## is therefore in its own column. then delete header 'taxonomy' and put in appropriate
## header for each column (order or family or genus etc. 
## Copy and paste into notepad and save. Then carry on: 
taxFile = read.table('oleate.taxonomy', header=T, sep='\t')
rownames(taxFile) = taxFile[,1]
taxFile
taxFile = taxFile[,3:8]
taxFile = as.matrix(taxFile)
head(taxFile)
## import metadata file
metaFile = read.table('oleate.metadata', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:4]
head(metaFile)
## Create phyloseq object
OTU = otu_table(sharedFile, taxa_are_rows = TRUE)
OTUsub = otu_table(sharedsubFile, taxa_are_rows = TRUE)
TAX = tax_table(taxFile)
META = sample_data(metaFile)
physeq = phyloseq(OTU, TAX, META)
physeqSub = phyloseq(OTUsub, TAX, META)
physeqSub
## Get rid of any OTUs not present in any samples and get relative abundance
microSub <- prune_taxa(taxa_sums(physeqSub) > 0, physeqSub)
microSubRel = transform_sample_counts(microSub, function(x) x / sum(x) )
microSubRelFilt = filter_taxa(microSubRel, function(x) mean(x) > 1e-5, TRUE)
# for subsampled shared file
sharedSubRelAubd = transform_sample_counts(sharedsubFile, function(x) x / sum(x) )
# create a tree and a new phyloseq object
library("ape")
random_tree = rtree(ntaxa(microSub), rooted=TRUE, tip.label=taxa_names(microSub))
plot(random_tree)
microSubt = merge_phyloseq(microSub, random_tree)
# colors
# Box plot
bp + scale_fill_manual(values=c("brown3", "steelblue","grey50"))
# Scatter plot
sp + scale_color_manual(values=c("brown3", "steelblue","grey50"))
# Plotting rarefaction
rarecurve(t(otu_table(microSubt)), step=50, cex=0.5)
# Plotting taxonomy
plot_bar(microSub, fill="Phylum") + facet_wrap(~condition, scales = "free_x", nrow = 1)+theme(text = element_text(size=16))+geom_bar(stat="identity")+
  scale_fill_brewer(type = "seq", palette = "Spectral")+
  theme(
    legend.background = element_rect(
      fill = "lemonchiffon", 
      colour = "grey50", 
      size = 1
    )
  )
plot_bar(microSub, fill="Phylum") + facet_wrap(~site, scales = "free_x", nrow = 1)+theme(text = element_text(size=16))+geom_bar(stat="identity")+
  scale_fill_brewer(type = "seq", palette = "Spectral")+
  theme(
    legend.background = element_rect(
      fill = "lemonchiffon", 
      colour = "grey50", 
      size = 1
    )
  ) 
plot_bar(microSub, fill="Phylum")+facet_grid(condition~site, scales = "free")+theme(text = element_text(size=16))+geom_bar(stat="identity")+
  scale_fill_brewer(type = "seq", palette = "Spectral")+
  theme(
    legend.background = element_rect(
      fill = "lemonchiffon", 
      colour = "grey50", 
      size = 1
    )
  )
# Plotting PCoA
my.PCoA2 <- ordinate(microSub, "PCoA", "bray")
plot_ordination(microSub, my.PCoA2, type = "group", color = "condition", shape = "site")+theme(text = element_text(size=20))+ geom_point(size=5) + scale_color_manual(values=c("brown3", "steelblue","grey50"))
plot_ordination(microSub, my.PCoA2, type = "group", color = "condition")+ facet_wrap(~site, scales = "free_x", nrow = 1) + geom_point(size=5)+theme(text = element_text(size=20)) + scale_color_manual(values=c("brown3", "steelblue","grey50"))
plot_ordination(microSub, my.PCoA2, type = "group", color = "condition", shape = "site")+theme(text = element_text(size=20))+ geom_point(size=5)+ 
  geom_line() + scale_color_manual(values=c("brown3", "steelblue","grey50"))
plot_ordination(microSub, my.PCoA2, type = "group", color = "condition")+theme(text = element_text(size=20))+ geom_point(size=5)+ stat_ellipse(level=0.9)  + scale_color_manual(values=c("brown3", "steelblue","grey50"))
# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = phyloseq::distance(microSubt, method="unifrac", weighted=F)
ordination = ordinate(microSubt, method="PCoA", distance=wunifrac_dist)
plot_ordination(microSubt, ordination, color="condition") + theme(aspect.ratio=1)+
  ggtitle("PCoA: unweigthed Unifrac")+geom_point(size=5)  + scale_color_manual(values=c("brown3", "steelblue","grey50"))
#Plot PCoA using the weighted UniFrac as distance
wunifrac_dist = phyloseq::distance(microSubt, method="unifrac", weighted=T)
ordUF = ordinate(microSubt, method="PCoA", distance=wunifrac_dist)
plot_ordination(microSubt, ordUF, color = "condition") + 
  ggtitle("PCoA: weigthed Unifrac")+geom_point(size=5) + scale_color_manual(values=c("brown3", "steelblue","grey50"))
# Normalization
phy <- normalize(microSubt, method = "relative")
# Plot taxonomy by group
plotbar(phy,level="Phylum")+facet_wrap(~sample_data(microSub)$condition, scales = "free_x", nrow = 1) +scale_fill_manual(values=c("brown3", "steelblue","grey50"))
plotbar(phy,level="Phylum", group="condition") + theme(axis.text.x = element_text(size = 16))
plotbar(phy,level="Phylum", group="site") + theme(axis.text.x = element_text(size = 16))
# Plot alpha diversity
plotalpha(microSub, group = "condition", color = c("brown3","grey50", "steelblue"))
plotalpha(microSub, group = "site")
# Test for significance between groups
beta_condition <-betatest(phy,group="condition")
beta_site <-betatest(phy,group="site")
# Biomarkers discovery and plotting
bio <- biomarker(physeq,group="condition")
plotmarker(bio,level="Genus")
# LefSe testing and plotting
lda <- ldamarker(physeq,group="condition")
lda
write.csv(lda, "lda.csv", row.names = FALSE)
lda2 <- read.csv(file = 'lda.csv')
plotLDA(lda2,group=c("10%Fat","60%Fat"),lda=5,pvalue=0.05)+theme(axis.text.y = element_text(size = 16))+theme(text = element_text(size=20))
plotLDA(lda2,group=c("10%Fat","Oleate"),lda=5,pvalue=0.05)+theme(axis.text.y = element_text(size = 16))+theme(text = element_text(size=20))
plotLDA(lda2,group=c("Oleate","60%Fat"),lda=5,pvalue=0.05) +theme(axis.text.y = element_text(size = 16))+theme(text = element_text(size=20))
# Test for significant OTUs using DESeq2
BiocManager::install("DESeq2")
library(DESeq2)
sample_data(microSub)$condition <- as.factor(sample_data(microSub)$condition)
ds = phyloseq_to_deseq2(microSub, ~ condition)
ds = DESeq(ds)
alpha = 0.01
res = results(ds, contrast=c("condition", "Oleate", "60%Fat"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig
res_sig = cbind(as(res_sig, "data.frame"), as(tax_table(microSub)[rownames(res_sig), ], "matrix"))
ggplot(res_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
  geom_jitter(size=5, width = 0.2) + theme(text = element_text(size=16)) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
res2 = results(ds, contrast=c("condition", "10%Fat", "60%Fat"), alpha=alpha)
res2 = res2[order(res2$padj, na.last=NA), ]
res_sig2 = res2[(res2$padj < alpha), ]
res_sig2
res_sig2 = cbind(as(res_sig2, "data.frame"), as(tax_table(microSub)[rownames(res_sig2), ], "matrix"))
ggplot(res_sig2, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
  geom_jitter(size=5, width = 0.2) + theme(text = element_text(size=16)) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
# Plot heatmap
devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)
library(viridis)
#> Loading required package: viridisLite
library(RColorBrewer)
plot_taxa_heatmap(microSub,
                  subset.top = 20,
                  VariableA = "condition",
                  heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                  transformation = "log10"
)
# Venn Diagram
install.packages("VennDiagram")
library(VennDiagram)
install.packages("eulerr") # If not installed
library(eulerr)
library(microbiome)
#devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities)
pseq<-microSubt
table(meta(pseq)$condition, useNA = "always")
# convert to relative abundances
pseq.rel <- microbiome::transform(pseq, "compositional")
disease_states <- unique(as.character(meta(pseq.rel)$condition))
print(disease_states)
#Write a for loop to go through each of the disease_states one by one and combine identified core taxa into a list.
list_core <- c() # an empty object to store information

for (n in disease_states){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, condition == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
print(list_core)
# Specify colors and plot venn
# supplying colors in the order they appear in list_core
mycols <- c("brown3", "grey50", "steelblue")
venn.diagram(list_core,fill = mycols, filename = '#14_venn_diagramm.png',
             output=TRUE)
# Plot rarefaction
# set seed
set.seed(1)
subsamples <- seq(0, 5000, by=100)[-1]
#subsamples = c(10, 5000, 10000, 20000, 30000)
p <- plot_alpha_rcurve(microSub, index="observed", 
                       subsamples = subsamples,
                       lower.conf = 0.025, 
                       upper.conf = 0.975,
                       group="condition",
                       label.color = "brown3",
                       label.size = 3,
                       label.min = TRUE) 
# change color of line 
mycols <- c("brown3", "steelblue","grey50")

p <- p + scale_color_manual(values = mycols) + 
  scale_fill_manual(values = mycols)
print(p)
# Modified difftest function from microbial package
difftest2<-function(physeq,group,pvalue=0.05,padj=NULL,log2FC=0,gm_mean=TRUE,fitType="local",quiet=FALSE){
  if(!taxa_are_rows(physeq)){
    physeq<-t(physeq)
  }
  otu <- as(otu_table(physeq),"matrix")
  tax <- as.data.frame(as.matrix(tax_table(physeq)))
  colData<-as(sample_data(physeq),"data.frame")
  colData$condition<-colData[,group]
  contrasts<-levels(factor(unique(colData$condition)))[1:2]
  if(isTRUE(gm_mean)){
    countData<-round(otu, digits = 0)
  }else{
    countData<-round(otu, digits = 0)+1
  }
  dds <- DESeqDataSetFromMatrix(countData, colData, as.formula(~condition))
  if(isTRUE(gm_mean)){
    geoMeans = apply(counts(dds), 1, gm_mean)
    dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  }
  dds <- DESeq(dds, fitType=fitType)
  res <- results(dds,contrast=c("condition",contrasts),cooksCutoff = FALSE)
  res_tax = cbind(as.data.frame(res), as.matrix(countData[rownames(res), ]))
  if(!is.null(padj)){
    pval<-padj
    sig <- rownames(subset(res,padj<pval&abs(log2FoldChange)>log2FC))
  }else{
    pval<-pvalue
    sig <- rownames(subset(res,pvalue<pval&abs(log2FoldChange)>log2FC))
  }
  res_tax$Significant<- "No"
  res_tax$Significant <- ifelse(rownames(res_tax) %in% sig, "Yes", "No")
  res_tax <- cbind(res_tax, tax[rownames(res),])
  return(as.data.frame(res_tax))
}
# Test for significant bacteria using microbial package
mic_res <- difftest2(microSub,group="condition", gm_mean=FALSE)
plotdiff(mic_res,level="Genus",padj=0.001,log2FC = 3,fontsize.y = 14)
# Plotting comparison between top 3 families in different conditions
library(microbiomeutilities)
library(RColorBrewer)
mycols <- c("brown3", "grey50", "steelblue")
?plot_taxa_boxplot
pn <- plot_taxa_boxplot(microSubt,
                        taxonomic.level = "Family",
                        top.otu = 3, 
                        group = "condition",
                        add.violin= FALSE,
                        title = "Top three family", 
                        keep.other = FALSE,
                        group.order = c("10%Fat","Oleate","60%Fat"),
                        group.colors = mycols,
                        dot.size = 1)

print(pn + theme_biome_utils())
library(ggpubr)
p0.f <- format_to_besthit(microSubt)
top_taxa(p0.f, 5)
select.taxa <- c("Otu004:Lactobacillus(100)", "Otu005:Lachnospiraceae_unclassified(100)")
mycols <- c("brown3", "grey50", "steelblue")
p <- plot_listed_taxa(p0.f, select.taxa, 
                      group= "condition",
                      group.order = c("10%Fat","Oleate","60%Fat"),
                      group.colors = mycols,
                      add.violin = TRUE,
                      violin.opacity = 0.3,
                      dot.opacity = 0.25,
                      box.opacity = 0.25,
                      panel.arrange= "grid")
print(p + ylab("Relative abundance") + scale_y_continuous(labels = scales::percent))
# If more than two variables
comps <- make_pairs(sample_data(p0.f)$condition)
print(comps)
p <- p + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test")
p + scale_y_continuous(labels = scales::percent)
# Plotting top 4 genus
pa <- aggregate_taxa(microSubt, "Genus")
top_four <- top_taxa(pa, 4)
top_four
mycols <- c("brown3", "steelblue", "grey50")
p <- plot_listed_taxa(pa, top_four, 
                      group= "condition",
                      group.order = c("10%Fat","Oleate","60%Fat"),
                      group.colors = mycols,
                      add.violin = TRUE,
                      violin.opacity = 0.3,
                      dot.opacity = 0.25,
                      box.opacity = 0.25,
                      panel.arrange= "wrap")
comps <- make_pairs(sample_data(pa)$condition)
p <- p + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test")
print(p + ylab("Relative abundance") + scale_y_continuous(labels = scales::percent))
# Get taxa summary by group(s)
p0 <- microSubt
p0.gen <- aggregate_taxa(microSubt,"Genus")
x.d <- dominant_taxa(p0,level = "Genus", group="condition")
head(x.d$dominant_overview)
tx.sum1 <- taxa_summary(p0, "Phylum")


p0.rel <- transform(p0, "compositional")

grp_abund <- get_group_abundances(p0.rel, 
                                  level = "Phylum", 
                                  group="condition",
                                  transform = "compositional")
mycols <- c("brown3", "steelblue","grey50")

# clean names 
grp_abund$OTUID <- gsub("p__", "",grp_abund$OTUID)
grp_abund$OTUID <- ifelse(grp_abund$OTUID == "", 
                          "Unclassified", grp_abund$OTUID)


mean.plot <- grp_abund %>% # input data
  ggplot(aes(x= reorder(OTUID, mean_abundance), # reroder based on mean abundance
             y= mean_abundance,
             fill=condition)) + # x and y axis 
  geom_bar(     stat = "identity", 
                position=position_dodge()) + 
  scale_fill_manual("condition", values=mycols) + # manually specify colors
  theme_bw() + # add a widely used ggplot2 theme
  ylab("Mean Relative Abundance") + # label y axis
  xlab("Phylum") + # label x axis
  coord_flip() # rotate plot 
mean.plot
# Density of reads per sample
p0 <- microSubt
print_ps(p0)
kable(head(tax_table(p0)))
# reduce size for example
ps0 <- core(p0, detection = 10, prevalence = 20 / 100)
# Add a prefix to taxa labels
ps0.f2 <- format_to_besthit(ps0, prefix = "MyBug-")
kable(head(tax_table(ps0.f2))[3:6])
p <- plot_read_distribution(ps0.f2, groups = "condition", 
                            plot.type = "density")
print(p + theme_biome_utils())
ps0 <- core(p0, detection = 10, prevalence = 20 / 100)
pseq_df <- phy_to_ldf(ps0, transform.counts = NULL)
kable(head(pseq_df))
# Plot Density of phyla
pseq <- microSubt

# check 10%Fat
control_ps <- subset_samples(pseq, condition=="10%Fat")
p_hc <- taxa_distribution(control_ps) + 
  theme_biome_utils() + 
  labs(title = "10%Fat")

# check Oleate
oleate_ps <- subset_samples(pseq, condition=="Oleate")
p_oleate <- taxa_distribution(oleate_ps) + 
  theme_biome_utils() + 
  labs(title = "Oleate")

# check 60%Fat
HFD_ps <- subset_samples(pseq, condition=="60%Fat")
p_hfd <- taxa_distribution(HFD_ps) + 
  theme_biome_utils() + 
  labs(title = "60%Fat")

# harnessing the power of patchwork
p_hc / p_oleate / p_hfd + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "A")

