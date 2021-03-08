microSub

library(devtools) # Load the devtools package
install_github("guokai8/microbial") # Install the package

library(microbial)

#default normalize method is relative
phy <- normalize(microSub, method = "relative")

plotbar(phy,level="Phylum")+facet_wrap(~sample_data(microSub)$condition, scales = "free_x", nrow = 1)

plotbar(phy,level="Phylum", group="condition") + theme(axis.text.x = element_text(size = 16))  

plotbar(phy,level="Phylum", group="site") + theme(axis.text.x = element_text(size = 16))

plotalpha(microSub, group = "condition")
plotalpha(microSub, group = "site")

plotbeta(phy, group="condition")

beta <-betatest(phy,group="condition")

META
mic_res <- difftest2(microSub,group="condition", gm_mean=FALSE)
plotdiff(mic_res,level="Genus",padj=0.001,log2FC = 3,fontsize.y = 14)

bio <- biomarker(physeq,group="condition")
plotmarker(bio,level="Genus")

lda <- ldamarker(physeq,group="condition")
lda
write.csv(lda, "lda.csv", row.names = FALSE)
lda2 <- read.csv(file = 'lda.csv')
plotLDA(lda2,group=c("10%Fat","60%Fat"),lda=5,pvalue=0.05)+theme(axis.text.y = element_text(size = 16))
plotLDA(lda2,group=c("10%Fat","Oleate"),lda=5,pvalue=0.05)+theme(axis.text.y = element_text(size = 16))
plotLDA(lda2,group=c("Oleate","60%Fat"),lda=5,pvalue=0.05) +theme(axis.text.y = element_text(size = 16))

?plotLDA

theme(text = element_text(size=20)
      
      
      
      
      
      
      
      
      
      
plot_net(microSub, "bray", color = "site", laymeth = "svd")
plot_tree(microSub)                                                                      
gp.ch = subset_taxa(microSub, Phylum == "Bacteriodetes")
plot_bar(gp.ch, fill="Genus")


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



