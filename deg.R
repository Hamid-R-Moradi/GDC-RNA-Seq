setwd('../')
getwd()
list.files()

#install.packages('BiocManager')	
#BiocManager::install(c("biomaRt","DESeq2","apelgm","vsn","ggplot2","dplyr"))

library('biomaRt') # for Annotation
library('DESeq2') # for Differential Expression Analysis
#library('apeglm')
library('vsn') # for Normalization
library('ggplot2') # for Visualization
library('dplyr') # for pipeline

# Importing data

## count data
rawdata  <- read.delim2('countM.tsv',header = F) 
rawdata[1,] <- paste0('id',rawdata[1,]); colnames(rawdata) <- rawdata[1,]; rawdata <- rawdata[-1,]
genes  <- rawdata[,1] ; rawdata <- rawdata[,-1] # Added 1st col as rownames and deletedthe 1st col
rawdata[1,1:3] # checking the first 3 cols of countM

## pheno data
phendata <- read.delim2("gdc_sample_sheet.2024-04-27.tsv",header = T)  
phendata[,2] <- tolower(x=phendata[,2])
phendata[,2] <- paste0('id',phendata[,2])
phendata[1:3,2]

# reordering phenodata according to rawdata
row.names(phendata) <- phendata[,2]
all(colnames(rawdata) %in% row.names(phendata)) # checking if all sample names are correct
all(colnames(rawdata) == row.names(phendata)) # checking if all samples are in correct pos

phendata <- phendata[colnames(rawdata),] # sorting data according to colnames rawdata
all(colnames(rawdata) == row.names(phendata)) # checking if all samples are in correct pos

# fixing data to numeric, removing NA counts and adding geneIDs
rawdata <- sapply(rawdata,as.numeric)
rawdata[is.na(rawdata)] <- 0
rownames(rawdata) <- genes
rownames(rawdata)[1:10]

# constructing DESeq datastructure
dds <- DESeqDataSetFromMatrix(countData = rawdata,
                              colData = phendata,
                              design = ~ Sample.Type)
dds


# pre-filtering 
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
dds

unique(dds$Sample.Type)
dds$Sample.Type <- factor(dds$Sample.Type, levels= c("Solid Tissue Normal","Primary Tumor","Recurrent Tumor"))
# dds$Sample.Type <- relevel(dds$Sample.Type, ref = "Solid Tissue Normal")

# Data quality assessment
boxplot(log2(counts(dds)+1)[,10:15]) # data contain a median at approximately same range so they are comparable although I have not checked all of them

## variance stabalizing transformation 
dds.vst <- vst(dds,blind = F,nsub=5000) # since the conditions are expected to be to a large extent different blind option is set to False

par(mfrow = c(1,2))
plot(log2(counts(dds)+1)[,1:2],cex=.1,main='Normalized Read Counts') # assessing homoskedasticity, for less than 1024 (2^10) there is a large heteroskedasticity, so the variance is dependent on the mean
plot(assay(dds.vst)[,1:2],cex=.1,main='VST Transformed Read Counts')

dev.off()

##mean-sd # crashes my R - FIXED
par(mfrow = c(1,2))

msdplot <- meanSdPlot(log2(counts(dds)+1),ranks =F,plot=F)
msdplot$gg + ggtitle('sequencing depth log2 normalized(read counts)') +
	ylab('standard deviation')
msdplot.vst <- meanSdPlot(assay(dds.vst),ranks =F,plot=F)
msdplot.vst$gg + ggtitle('sequencing depth log2 normalized(read counts)') +
	ylab('standard deviation')

## clustering 
dist  <- as.dist(1- cor(assay(dds.vst),method='pearson'))
plot(hclust(dist),labels=phendata[,8],cex=.5,main="VST transfomed data \n distance: Pearson correlation")
dev.off()
##PCA
#pc <- plotPCA(dds.vst,intgroup='Sample.Type')
#pc <- pc + theme_bw() + ggtitle("VST transformed data")
#print(pc) # there are two p\suimary tumors samples which are clustered with Solid Tissue Normal samples

pc <- prcomp(t(x=assay(dds.vst)),scale=T)
var.explained <- pc$sdev^2/sum(pc$sdev^2)


pc$x %>% as.data.frame %>% mutate(Sample.Type=phendata$Sample.Type) %>% ggplot(aes(x=PC1,y=PC2,color=phendata$Sample.Type)) + geom_point(size=.8) + theme_bw() + 
	labs(x=paste0("PC1: ",round(var.explained[1]*100,1),"%"),
	     y=paste0("PC2: ",round(var.explained[2]*100,1),"%")) +
	guides(guide_legend(title = "Sample Type")) +
	theme(legend.position="top")

identify(x, y, labels = name, plot=TRUE)

#DEG Analysis
dds1 <- DESeq(dds)
res <- results(dds1,contrast=c("Sample.Type","Primary Tumor","Solid Tissue Normal"),alpha = 0.01)
res

# saving the DESeq analysis for future
saveRDS(dds1,file="fittedModel_LIHC.rds")
list.files()

reshfc <- lfcShrink(dds1,coef=2,type='apeglm')
reshfc
# --------------------- #
plotMA(res,ylim=c(-2,2))
plotMA(reshfc,ylim=c(-2,2))

plotCounts(dds, gene=which.min(res$padj), intgroup="Sample.Type")

mcols(res)$description
