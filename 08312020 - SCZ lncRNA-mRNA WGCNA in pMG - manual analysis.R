#==========================================================================================================#
#                                              required packages                                           #
#==========================================================================================================#

library(grid)
library(purrr)
library(WGCNA)
library(cluster)
library(DESeq2)
library(ggplot2)
library(limma)
library(sva)
library(variancePartition)
library(readr)
library(Hmisc)
library(ape)
library(RColorBrewer)
library(pheatmap)
library(vsn)
library(ggcorrplot)
library(biomaRt)
library(dplyr)
library(data.table)
library(DEGreport)
library(gprofiler2)
library(sna)
library(igraph)
library(GGally)
library(network)
library(ggplot2)
library(viridis) 



#==========================================================================================================#
#                               Analysis start: counts & covariate files                                   #
#==========================================================================================================#

#File directory
dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", "Internship", "Internship_Start", "Project 4 - Synthesis", "Analysis", "Counts")
setwd(dir)

#Reading in primary microglia data
pmg <- read.table("08242020.pmg.counts.complete.txt", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)

#Filter out missing values
pmg <- pmg[complete.cases(pmg), ]

#Grab covariate files
dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", "Internship", "Internship_Start", "Project 4 - Synthesis", "Analysis", "Metadata")
setwd(dir)

#Reading in covariates
all_covs <- read.table("metadata_allcovs.txt", header=T, row.names=1)[c(1:2,4:15)]
#all_covs_sub <- all_covs_sub[rownames(all_covs_sub) != "99_23",]
#Remove character variables
all_covs_sub <- all_covs_sub[c(1,4:14)]

#Making numerical variables readable
for (i in c("RIN", "Exonic.Rate", "Intergenic.Rate", "rRNA.rate", "Mean.Per.Base.Cov."))
   all_covs_sub[[i]] <- as.numeric(gsub(",",".",all_covs_sub[[i]]))
#Setting factor variables
for (i in c("batch", "sex", "status", "region_num"))
  all_covs_sub[[i]] <- as.factor(all_covs_sub[[i]])



#==========================================================================================================#
#                                    Filtering & outlier detection                                         #
#==========================================================================================================#

#Create a dds object for normalization and filtering with DESeq2 functions
dds <- DESeqDataSetFromMatrix(countData = pmg, colData = all_covs_sub, formula(~ sex  + batch + region_num + status))
keep <- rowSums(assay(dds)) > 1 #Minimal filtering: filter out all genes that are not expressed in the dataset
dds <- dds[keep,]
pmg <- data.frame(assay(dds))


#Detecing outliers based on hierarchically clustering sample correlation
d <- bicor(pmg, use = "p") #same results as with vst
hc <- hclust(dist(1 - d))
#plot(hc, xlab="", sub="", main = "Sample dissimilarity clustering",
#     labels = FALSE, hang = 0.04);
plot.phylo(as.phylo(hc), type = "p", edge.col = "blue", edge.width = 2, 
           show.node.label = TRUE, no.margin = TRUE)


#Checking sample clustering based on microglia core signature
dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", "Internship", "Internship_Start", "Project 4 - Synthesis", "Analysis", "Genelists")
setwd(dir)
coresign <- read_lines("PatirEtAl_coreMicrogliaGenes.csv")
matchedmicgenes <- read_lines("PatirEtAl_MatchedMicrogliaGenes.csv")
allmicgenes <- read_lines("PatirEtAl_AllMicGenes.csv")
dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", "Internship", "Internship_Start", "Project 4 - Synthesis", "Analysis")
setwd(dir)

#Annotation
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
IDannot_clean <- IDs(allmicgenes, t(pmg))
IDannot_clean.Expr <- pmg[rownames(pmg) %in% IDannot_clean[,1],]

#Cluster samples based on microglia gene signature expression
plotClusterTreeSamples(t(IDannot_clean.Expr))
hc_coreExpr <- hclust(dist(1 - t(IDannot_clean.Expr)))
plot.phylo(as.phylo(hc_coreExpr), type = "p", edge.col = "blue", edge.width = 2, 
           show.node.label = TRUE, no.margin = TRUE)


#Cluster samples based on microglia gene signature expression correlation
d_core <- bicor(IDannot_clean.Expr, use="p")
plotClusterTreeSamples(dist(1 - d_core))

hc_core <- hclust(dist(1 - d_core))
plot.phylo(as.phylo(hc_core), type = "p", edge.col = "blue", edge.width = 2, 
           show.node.label = TRUE, no.margin = TRUE)


#Statistical perspective: deviation from mean based on interquartile range
#VISUALIZING BOXPLOTS
barplot.df <- data.frame(Intersample.correlation=c(rowMeans(d), rowMedians(d)), Metric=c(rep("Correlation Means", 34), rep("Correlation Medians", 34)))
barplot.df %>%
  ggplot( aes(x=Metric, y=Intersample.correlation, fill=Metric)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
#  theme_ipsum() +
  theme(
    legend.position="Type",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Boxplots showing outliers based on inter-sample correlation means and medians before data correction") +
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))+
  xlab("")+
#geom_hline(yintercept=median(d_core)-3*mad(d_core),  linetype="dashed")+
  geom_segment(aes(x=1.7,xend=2.3,y=median(d)-3*mad(d),yend=median(d)-3*mad(d)),color="red", linetype="dashed")+
  geom_segment(aes(x=1.7,xend=2.3,y=median(d)+3*mad(d),yend=median(d)+3*mad(d)),color="red", linetype="dashed")+
  geom_segment(aes(x=0.7,xend=1.3,y=mean(d)-3*sd(d),yend=mean(d)-3*sd(d)),color="red", linetype="dashed")+
  geom_segment(aes(x=0.7,xend=1.3,y=mean(d)+3*sd(d),yend=mean(d)+3*sd(d)),color="red", linetype="dashed")



#==========================================================================================================#
#                        Data transformation/normalization and sample cluster QC                           #
#==========================================================================================================#

#Normalize
vsd <- vst(dds)

#Check for outliers after transformation:
plotClusterTreeSamples(t(assay(vsd)))
hc_Expr <- hclust(dist(1 - t(assay(vsd))))
plot.phylo(as.phylo(hc_Expr), type = "p", edge.col = "blue", edge.width = 2, 
           show.node.label = TRUE, no.margin = TRUE)

#Detecing outliers based on correlation
d <- bicor(assay(vsd), use = "p") #same results as with vst
plotClusterTreeSamples(dist(1 - d))
plotClusterTreeSamples(d)

hc <- hclust(dist(1 - d))
plot.phylo(as.phylo(hc), type = "p", edge.col = "blue", edge.width = 2, 
           show.node.label = TRUE, no.margin = TRUE)



#==========================================================================================================#
#                          PCA and Variance Partititioning before correction                               #
#==========================================================================================================#
#I decided to run the PCA on the top 10000 most connected
#Choosing a soft-threshold (using the WGCNA package to do so)
#Estimate SFT-criterion
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft.prefilter = pickSoftThreshold(
  t(assay(vsd)), 
  dataIsExpr = TRUE,
  weights = NULL,
  RsquaredCut = 0.8, 
  powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)), 
  removeFirst = FALSE, nBreaks = 20, blockSize = NULL, 
  corFnc = bicor, corOptions = list(use = 'p'), 
  networkType = "unsigned",
  moreNetworkConcepts = FALSE,
  gcInterval = NULL,
  verbose = 0, indent = 0)
sft.prefiler.fitind <- sft.min.prefilter$fitIndices

collectGarbage()

#Calculate connectivity to choose top ten thousand most connected genes
k.prefilter=softConnectivity(
  t(assay(vsd)),
  corFnc= "bicor",corOptions = "use = 'p',maxPOutliers = 0.05", #If using bicor, put maxPOutliers at 0.05 or 0.10
  weights = NULL,
  type = "unsigned",
  power = 6,
  blockSize = 20000,
  minNSamples = NULL,
  verbose = 0, indent = 0)

pre.filter <- t(assay(vsd))[, rank(-k.prefilter,ties.method="first")<=10000] 
pre.filter.subset <- vsd[rownames(assay(vsd)) %in% colnames(pre.filter)]


pca <- plotPCA.custom(pre.filter.subset, intgroup=c("status", "region", "batch", "age", "RIN", "Exonic.Rate"), ntop = 10000, returnData=TRUE)
#Grab proportion of variance explained by each PC:
PoV <- round(100 * attr(pca, "percentVar"))

#pcplot <- function(pca, aestocolor, title, PoV, colortitle, labels, colors){
#  ggplot(pca, aes(PC1, PC2, color=aestocolor)) +
#    geom_point(size=2) +
#    labs(title=title, x = paste0("PC1: ",PoV[1],"% variance"), y = paste0("PC2: ",PoV[2],"% variance"), color = colortitle) +
#    scale_color_manual(labels=labels, values = colors) +
#    geom_text(aes(label=name),vjust=2,check_overlap = FALSE,size = 4)+
#    theme_bw()+
#    theme(plot.title = element_text(hjust=0.5))
#}

#for (i in colnames(all.factors))
#  plot <- pcplot(pca, all.factors$i, title=i, PoV, "Status", c("CTL", "MDD"), c("turquoise", "indianred1"))
#  print(plot+plot)

#pcplot(pca, Diagnosis, title="DIAGNOSIS", PoV, "Status", c("CTL", "MDD"), c("turquoise", "indianred1"))

#Plotting
factor.covs <- data.frame(cbind("Diagnosis"=all_covs_sub$status, "Region"=as.factor(all_covs_sub$region), 
                                "Batch" = as.factor(all_covs_sub$batch)))
rownames(factor.covs) <- rownames(all_covs_sub)
num.covs <- data.frame(cbind("Age"=all_covs_sub$age, 
                          "RIN"=all_covs_sub$RIN, "EXONIC RATE"=all_covs_sub$Exonic.Rate))
rownames(num.covs) <- rownames(all_covs_sub)

plot.list <- vector("list", length(colnames(all.factors)))
names(plot.list) <- c("DIAGNOSIS", "REGION", "BATCH", "AGE", "RIN", "EXONIC RATE")
for (i in names(plot.list))

plot.list <- list("DIAGNOSIS"=c("CTL", "MDD"), "REGION"=c("GFM", "GTS"),
                      "BATCH"=c("1", "2", "3", "4"))

plot.list.clrs <- list("DIAGNOSIS"=c("turquoise", "indianred1"), "REGION"=c("turquoise", "indianred1"),
                  "BATCH"=c("turquoise", "indianred1", "blue", "magenta"))

plots <- NULL
for (i in names(plot.list))
  for (n in names(factor.covs))
    plots[i] <- ggplot(pca, aes(PC1, PC2, color=n)) +
  geom_point(size=2) +
  labs(title = i, x = paste0("PC1: ",PoV[1],"% variance"), y = paste0("PC2: ",PoV[2],"% variance"), color = i) +
  scale_color_manual(labels = plot.list[[i]], values = plot.list.clrs[[i]]) +
  geom_text(aes(label=name),vjust=2,check_overlap = FALSE,size = 4)+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))

plots2
for (i in c("AGE", "RIN", "EXONIC RATE"))
  for (n in names(factor.covs))
  plots2[i] <- ggplot(pca, aes(PC1, PC2, color=n)) +
  geom_point(size=2) +
  labs(title = i, x = paste0("PC1: ",PoV[1],"% variance"), y = paste0("PC2: ",PoV[2],"% variance"), color = i) +
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))

print(plots + plots2)



#==========================================================================================================#
#                                   Analysis start: covariate analysis                                     #
#==========================================================================================================#
#Covariation among technological & biological variables
COR = bicor(all_covs_sub, use="p")
pheatmap(COR, legend_breaks = c(-1,-0.8,-0.6,-0.4,-0.2, 0, 0.2 , 0.4, 0.6, 0.8, 1, max(COR)), legend = T, main = "Correlation matrix of measured covariates",
         legend_labels = c("-1","-0.8","-0.6","-0.4","-0.2","0", "0.2", "0.4", "0.6", "0.8", "1", "Bi-weight Midcorrelation\n\n"),
         color = colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(50), display_numbers = FALSE)

#Variance partition analysis
regform <- ~ batch + sex + age + status + region_num + RIN + Genes.Detected + Exonic.Rate + Intergenic.Rate + rRNA.rate + Mean.Per.Base.Cov.
varPart <- fitExtractVarPartModel(assay(vsd), regform, all_covs_sub)
vp <- sortCols(varPart)
plotVarPart(vp) +
  ggeasy::easy_rotate_x_labels(angle = 70, side =c("right"))

#Might as well check PC-covariate correlation
degCovariates(assay(vsd), all_covs_sub, fdr = 0.05, scale = TRUE, minPC = 5,
              correlation = "spearman", addCovDen = TRUE, legacy = FALSE,
              smart = TRUE, method = "lm", plot = TRUE)



#==========================================================================================================#
#                                       sva_network data correction                                        #
#==========================================================================================================#

vsd.corrected <- vsd
#Correction based on sv_network solely, only retaining regional effects

#SV identification and correction based on that
mat <- model.matrix(~ region_num, data = all_covs_sub)

n.svs <- num.sv(assay(vsd.corrected), mat, method="be", B = 20)
#n.svs <- num.sv(assay(vsd), as.matrix(all_covs_sub[c(1:8,10:12)]), method="be") #Number of mapped reads must be excluded due to colinearity
assay(vsd.corrected) <- data.frame(sva_network(assay(vsd.corrected), n.svs))

assay(vsd.corrected) <- assay(vsd.corrected) <- removeBatchEffect(assay(vsd.corrected), batch=all_covs_sub$batch)



#==========================================================================================================#
#                                Sample cluster QC after data correction                                   #
#==========================================================================================================#

#Let's check whether any outlier remain
plotClusterTreeSamples(t(assay(vsd.corrected)))
hc_Expr.corrected <- hclust(dist(t(assay(vsd.corrected))))
plot.phylo(as.phylo(hc_Expr.corrected), type = "p", edge.col = "blue", edge.width = 2, 
           show.node.label = TRUE, no.margin = TRUE)

#Detecing outliers based on correlation
d.corrected <- bicor(assay(vsd.corrected), use = "p") #same results as with vst
plotClusterTreeSamples(dist(1 - d.corrected),  main = "Sample dendrogram after normalization and PC-driven correction")

hc.corrected <- hclust(dist(1 - d.corrected))

plot.phylo(as.phylo(hc.corrected), type = "p", edge.col = "blue", edge.width = 2, 
           show.node.label = TRUE, no.margin = TRUE)



#==========================================================================================================#
#                          PCA and Variance Partitioning after correction                                  #
#==========================================================================================================#

#Check PC correlation with data subsequent to correction
degCovariates(assay(vsd.corrected), all_covs_sub, fdr = 0.05, scale = TRUE, minPC = 5,
              correlation = "spearman", addCovDen = TRUE, legacy = FALSE,
              smart = TRUE, method = "lm", plot = TRUE)

subset <- vsd.corrected[rownames(assay(vsd.corrected)) %in% colnames(Subtrans_norm)]
#PCA
pcs <- plotPCA.custom(subset , intgroup=c("RIN", "status", "sex", "age", "region", "Exonic.Rate"), ntop = 10000, returnData=TRUE)
#Grab proportion of variance explained by each PC:
PoV <- round(100 * attr(pcs, "percentVar"))

plots <- NULL
for (i in names(plot.list))
  for (n in names(factor.covs))
    plots[i] <- ggplot(pcs, aes(PC1, PC2, color=n)) +
  geom_point(size=2) +
  labs(title = i, x = paste0("PC1: ",PoV[1],"% variance"), y = paste0("PC2: ",PoV[2],"% variance"), color = i) +
  scale_color_manual(labels = plot.list[[i]], values = plot.list.clrs[[i]]) +
  geom_text(aes(label=name),vjust=2,check_overlap = FALSE,size = 4)+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))

plots2
for (i in c("AGE", "RIN", "EXONIC RATE"))
  for (n in names(factor.covs))
    plots2[i] <- ggplot(pcs, aes(PC1, PC2, color=n)) +
  geom_point(size=2) +
  labs(title = i, x = paste0("PC1: ",PoV[1],"% variance"), y = paste0("PC2: ",PoV[2],"% variance"), color = i) +
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))

print(plots + plots2)

#Plotting the VP after correction 
varPart2 <- fitExtractVarPartModel(assay(vsd.corrected), regform, all_covs_sub)
vp2 <- sortCols(varPart2)
plotVarPart(vp2) +
  ggeasy::easy_rotate_x_labels(angle = 70, side =c("right"))



#==========================================================================================================#
#                                       Data Transformation QC                                             #
#==========================================================================================================#

#Checking effect of transformation on sample distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd.corrected$condition,vsd.corrected$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, main="Sample distance matrix before correction")

#Checking variance distribution
meanSdPlot(assay(vsd))
meanSdPlot(as.matrix(assay(vsd.corrected)))

trans_norm <- data.frame(t(assay(vsd.corrected)))

trans_cor <- bicor(t(trans_norm), use = "p")
trans_cor.dat <- data.frame(Intersample.correlation=c(rowMeans(trans_cor), rowMedians(trans_cor)), Metric=c(rep("Correlation Means", 34), rep("Correlation Medians", 34)))

summary(rowMedians(trans_cor))
trans_cor.dat %>%
  ggplot( aes(x=Metric, y=Intersample.correlation, fill=Metric)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  #  theme_ipsum() +
  theme(
    legend.position="Type",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Boxplots showing outliers based on inter-sample correlation means and medians after data correction") +
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))+
  xlab("")+
  #geom_hline(yintercept=median(d_core)-3*mad(d_core),  linetype="dashed")+
  geom_segment(aes(x=1.7,xend=2.3,y=median(trans_cor)-3*mad(trans_cor),yend=median(trans_cor)-3*mad(trans_cor)),color="red", linetype="dashed")+
  geom_segment(aes(x=1.7,xend=2.3,y=median(trans_cor)+3*mad(trans_cor),yend=median(trans_cor)+3*mad(trans_cor)),color="red", linetype="dashed")+
  geom_segment(aes(x=0.7,xend=1.3,y=mean(trans_cor)-3*sd(trans_cor),yend=mean(trans_cor)-3*sd(trans_cor)),color="red", linetype="dashed")+
  geom_segment(aes(x=0.7,xend=1.3,y=mean(trans_cor)+3*sd(trans_cor),yend=mean(trans_cor)+3*sd(trans_cor)),color="red", linetype="dashed")



#==========================================================================================================#
#                       SFT criterion evaluation and dataframe subsetting                                  #
#==========================================================================================================#

#Choosing a soft-threshold
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft.min.prefilter = pickSoftThreshold(
  trans_norm, 
  dataIsExpr = TRUE,
  weights = NULL,
  RsquaredCut = 0.8, 
  powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)), 
  removeFirst = FALSE, nBreaks = 20, blockSize = NULL, 
  corFnc = bicor, corOptions = list(use = 'p'), 
  networkType = "unsigned",
  moreNetworkConcepts = FALSE,
  gcInterval = NULL,
  verbose = 0, indent = 0)
sft.minprefilter.fitIndices <- sft.min.prefilter$fitIndices

collectGarbage()

plot(sft.minprefilter.fitIndices[,1], -sign(sft.minprefilter.fitIndices[,3])*sft.minprefilter.fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology corrected Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft.minprefilter.fitIndices[,1], -sign(sft.minprefilter.fitIndices[,3])*sft.minprefilter.fitIndices[,2],
     labels=powers,cex=cex1,col="red");

#Setting an R^2 cutoff at height h
abline(h=0.8,col="red")

#Mean connectivity as a function of the soft-thresholding power
plot(sft.minprefilter.fitIndices[,1], sft.minprefilter.fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft.minprefilter.fitIndices[,1], sft.minprefilter.fitIndices[,5], labels=powers, cex=cex1,col="red")


#Calculate connectivity to choose top ten thousand most connected genes
k=softConnectivity(
  trans_norm,
  corFnc= "bicor",corOptions = "use = 'p',maxPOutliers = 0.05", #If using bicor, put maxPOutliers at 0.05 or 0.10
  weights = NULL,
  type = "unsigned",
  power = 6,
  blockSize = 20000,
  minNSamples = NULL,
  verbose = 0, indent = 0)

Subtrans_norm=trans_norm[, rank(-k,ties.method="first")<=10000] 

sft.top10k = pickSoftThreshold(
  Subtrans_norm, 
  dataIsExpr = TRUE,
  weights = NULL,
  RsquaredCut = 0.8, 
  powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)), 
  removeFirst = FALSE, nBreaks = 20, blockSize = NULL, 
  corFnc = bicor, corOptions = "use = 'p',maxPOutliers = 0.05",
  networkType = "unsigned",
  moreNetworkConcepts = TRUE, #Interesting to compare correcteds on
  gcInterval = NULL,
  verbose = 0, indent = 0)
sfttop10k.fitindices <- sft.top10k$fitIndices
#save(sft, sfttop10k,
#     file = "c.RData")
#When including top 10 thousand most connected genes the slope becomes very strange > use power of 6 to retain connectivity?

#Plotting the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

#Scale-free topology is reached at 14 or 15 (R=0.8); so should go for either, because k connectivity is at 60-80
plot(sft.top10k$fitIndices[,1], -sign(sft.top10k$fitIndices[,3])*sft.top10k$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology corrected Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft.top10k$fitIndices[,1], -sign(sft.top10k$fitIndices[,3])*sft.top10k$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

#Setting an R^2 cutoff at height h
abline(h=0.8,col="red")

#Mean connectivity as a function of the soft-thresholding power
plot(sft.top10k$fitIndices[,1], sft.top10k$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft.top10k$fitIndices[,1], sft.top10k$fitIndices[,5], labels=powers, cex=cex1,col="red")


###################################PCA after subsetting########################################################

#PCA with subsetted genes
vsd.corrected_sub <- vsd.corrected[colnames(Subtrans_norm),]
vsd_sub <- vsd[colnames(Subtrans_norm),]

pcs <- plotPCA.custom(vsd.corrected_sub, intgroup=c("RIN", "status", "sex", "age", "region"), ntop = 10000, returnData=TRUE)
#Grab proportion of variance explained by each PC:
PoV <- round(100 * attr(pcs, "percentVar"))

ggplot(pcs, aes(PC1, PC2, color=Diagnosis)) +
  geom_point(size=2) +
  labs(title = "DIAGNOSIS", x = paste0("PC1: ",PoV[1],"% variance"), y = paste0("PC2: ",PoV[2],"% variance"), color = "STATUS") +
  scale_color_manual(labels = c("CTL", "MDD"), values = c("turquoise", "indianred1")) +
  geom_text(aes(label=name),vjust=2,check_overlap = FALSE,size = 4)+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))



#==========================================================================================================#
#                                       WGCNA CONSTRUCTION                                                 #
#==========================================================================================================#

#Running one-step network construction and module detection
cor <- WGCNA::cor #Temporarily reassing correlation function since this will otherwise cause an error in blockwiseModules(); else use detach()
#Calculating the suitable threshold for module merging
nSamples = nrow(Subtrans_norm)
mergeCutHeight = dynamicMergeCut(nSamples, mergeCor = 0.9, Zquantile = 2.33) #Use Z = 1.96 for 95% confidence instead of 99% confidence at Z = 2.33
#I will set it manually higher, as right now it doesn't merge moduels that are too close

#Adjacency and TOM matrix construction
adjacency = adjacency(Subtrans_norm, power = 6, type="unsigned", corFnc = "bicor", 
                      corOptions = list(use = "p"));
TOM = TOMsimilarity(adjacency, suppressNegativeTOM = FALSE, TOMType ="signed", verbose = 0);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
#Plotting
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

#Module construction
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, #Much easier to play with the parameters manually
                           method="hybrid", deepSplit = 4, pamStage= TRUE, pamRespectsDendro = TRUE,
                           minClusterSize = 40, respectSmallClusters=TRUE, verbose=0);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods, zeroIsGrey = TRUE, naColor = "grey")
table(dynamicColors)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors")
MEList = moduleEigengenes(Subtrans_norm, colors = dynamicColors, softPower=6)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")+
  abline(h=MEDissThres, col = "red")
#Setting mergeCutHeight manually
mergeCutHeight = 0.25
MEDissThres = mergeCutHeight 
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function; iterated until no modules are merged
merge = mergeCloseModules(Subtrans_norm, dynamicColors, MEs = MEs, 
                          corFnc = bicor, corOptions = "use = 'p',maxPOutliers = 0.05",
                          cutHeight = MEDissThres, iterate= FALSE, verbose = 0,
                          unassdColor = if (is.numeric(colors)) 0 else "grey",
                          getNewMEs = TRUE, getNewUnassdME = TRUE)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = orderMEs(mergedMEs);
sort(table(moduleColors))

#Plotting
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#Multidimensional scaling plot of the first two PCs and the MEs
sizeGrWindow(12, 9)
oldMEs <- MEList$eigengenes
par(mar = c(6, 8.5, 3, 3));
distPC <- 1-abs(bicor(oldMEs,use="all.obs", maxPOutliers = 0.05))
distPC <- ifelse(is.na(distPC), 0, distPC)
MDS <- cmdscale(as.dist(distPC),2)
modNames = substring(names(oldMEs), 3)
col <- names(table(modNames))
plot(MDS, col=col, main="Multi-Dimensional Scaling (MDS) of Module Eigengenes (MEs) before Merging", cex=2, pch=19, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

sizeGrWindow(12, 9)
par(mar = c(6, 8.5, 3, 3));
distPC <- 1-abs(bicor(MEs,use="all.obs", maxPOutliers = 0.05))
distPC <- ifelse(is.na(distPC), 0, distPC)
MDS <- cmdscale(as.dist(distPC),2)
modNames = substring(names(MEs), 3)
col <- names(table(modNames))
plot(MDS, col=col, main="Multi-Dimensional Scaling (MDS) of Module Eigengenes (MEs) after Merging", cex=2, pch=19, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

#Plot which sample contributes to what module by plotting sample against module eigengenes in a heatmap
treelist = geneTree[[1]]
ordergenes = treelist[[3]]
plotMat((Subtrans_norm[ordergenes,]) , rlabels= moduleColors[ordergenes], clabels=
                 colnames(Subtrans_norm), rcols=moduleColors[ordergenes]) 

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
#For old MEs
annot.bfcut <- data.frame(colnames(oldMEs))
rownames(annot.bfcut) <- colnames(oldMEs)
colnames(annot.bfcut) <- "Module.Colors"

annot_clrs.bfcut <- list()
for (i in as.character(data.frame(table(dynamicColors))[,1]))
  annot_clrs.bfcut$Module.Colors[[i]] = i

names(annot_clrs.bfcut$Module.Colors) <- colnames(oldMEs)

pheatmap(bicor(oldMEs), cluster_rows = T, cluster_cols = T, annotation_legend = FALSE, 
                         annotation_col=annot.bfcut, annotation_colors = annot_clrs.bfcut,
                         color = colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(50), display_numbers = F, 
                         main="Module Eigengene Correlation Before Module Merging")

#Suppl. mat. figure; sample contribution
pheatmap(oldMEs, cluster_rows = T, cluster_cols = T, annotation_legend = FALSE, 
         annotation_col=annot.bfcut, annotation_colors = annot_clrs.bfcut,
         color = colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(50), display_numbers = F,
         main="Sample Contribution to Module Eigengene Before Module Merging")

pheatmap(MEs, cluster_rows = T, cluster_cols = T, annotation_legend = FALSE, 
         annotation_col=annot.bfcut, annotation_colors = annot_clrs.bfcut,
         color = colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(50), display_numbers = F,
         main="Sample Contribution to Module Eigengene After Module Merging")


#For new MEs
annot.afcut <- data.frame(colnames(MEs))
rownames(annot.afcut) <- colnames(MEs)
colnames(annot.afcut) <- "Module.Colors"

summarized.modClrs <- data.frame(table(mergedColors))[,1]
#data.frame(table(mergedColors))[!data.frame(table(mergedColors))[,1]=="grey",]

summarized.modClrs


annot_clrs.afcut <- list()
for (i in as.character(summarized.modClrs))
  annot_clrs.afcut$Module.Colors[[i]] = i

names(annot_clrs.afcut$Module.Colors) <- colnames(MEs)[order(colnames(MEs))]

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(2,2))
pheatmap(bicor(MEs), cluster_rows = T, cluster_cols = T, annotation_legend = FALSE, 
         legend_breaks = c(-1,-0.8,-0.6,-0.4,-0.2, 0, 0.2 , 0.4, 0.6, 0.8, 1, max(bicor(MEs))),
         legend_labels = c("-1","-0.8","-0.6","-0.4","-0.2", "0", "0.2", "0.4", "0.6", "0.8", "1", "Biweight midcorrelation \n \n"),
         annotation_col=annot.afcut, annotation_row=annot.afcut, annotation_colors = annot_clrs.afcut,
         color = colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(50), display_numbers = F, 
         main="Module Eigengene (ME) Correlation After Merging")

SampleAnnot <- cbind(data.frame("Batch" = all_covs_sub[,1]), data.frame("Region" = all_covs_sub[,4]), data.frame("Diagnosis" = all_covs_sub[,5]))
rownames(SampleAnnot) <- rownames(Subtrans_norm)
pheatmap(MEs, cluster_rows = T, cluster_cols = T, annotation_legend = T, annotation_row=SampleAnnot,
         annotation_col=annot.afcut, annotation_colors = annot_clrs.afcut, main="Sample Contribution to Module Eigengenes (ME)",
         color = colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(50), display_numbers = F)

#==========================================================================================================#
#                                        ME-trait QC check                                                 #
#==========================================================================================================#

#Before going any further with the analysis: correlate MEs with traits
METraitCor <- rcorr(as.matrix(all_covs_sub), as.matrix(MEs), type="spearman")
METraitCor <- bicor(as.matrix(all_covs_sub), as.matrix(MEs), use="p")
#MECorres <- METraitCor$r[13:45,13:45]
#TraitCor <- METraitCor$r[c(1:12), c(1:12)]
#METraitCorRes <- METraitCor$r[13:45,1:12]
#METraitCorResP <- METraitCor$P[13:45,1:12]
#METraitCorResFDR <- p.adjust(METraitCorResP, method="hochberg", n=length(METraitCorResP)) #none of which, except SV4, would survive correction

MECorres <- bicor(as.matrix(MEs), use="p")
TraitCor <- bicor(as.matrix(all_covs_sub), use="p")


#moduleTraitCor = bicor(MEs, select, "use = 'p'");
#moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

textMatrix = paste(signif(METraitCorRes, 2), "\n(",
#                   signif(METraitCorResP, 1), ")", "\n(",
                   signif(METraitCorResFDR, 0), ")",sep = "");
dim(textMatrix) = dim(METraitCorRes)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = METraitCorRes,
               xLabels = colnames(METraitCorRes),
               yLabels = rownames(METraitCorRes),
               ySymbols = rownames(METraitCorRes),
               textMatrix=textMatrix,
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x = 0.9,
               cex.lab.y = 0.9,
               zlim = c(-1,1),
               main = paste("Module eigengene (ME) correlation with study covariates"))

pheatmap(t(METraitCor), cluster_rows = T, cluster_cols = T, annotation_legend = FALSE, 
         annotation_row=annot.bfcut, annotation_colors = annot_clrs.bfcut, 
         color = colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(50), display_numbers = T, main="Merged Module Eigengene (ME) - Trait Correlation")

diag(dissTOM) <- NA
TOMplot(dissTOM^6, geneTree, as.character(dynamicColors))
#==========================================================================================================#
#                            Microglia gene enrichment in modules                                          #
#==========================================================================================================#

#Check method validity by checking which modules contain which microglia-specific genes;
coreIDs <- IDs(values=coresign, object=trans_norm)
table(names(Subtrans_norm) %in% coreIDs[,1])["TRUE"] 
table(moduleColors[colnames(Subtrans_norm) %in% coreIDs[,1]])

broadCoreIDs <- IDs(values=matchedmicgenes, object=trans_norm) 
table(names(Subtrans_norm) %in% broadCoreIDs[,1])["TRUE"] 
table(moduleColors[colnames(Subtrans_norm) %in% broadCoreIDs[,1]])

AllIDs <- IDs(values=allmicgenes, object=trans_norm)
table(names(Subtrans_norm) %in% AllIDs[,1])["TRUE"] 
pheatmap(data.frame(table(moduleColors[colnames(Subtrans_norm) %in% AllIDs[,1]])))



#==========================================================================================================#
#                                 Starting off the core analysis                                           #
#==========================================================================================================#

#Starting the analysis plan off:
#0) Check expression
#a) identifying which modules contain these, b) identifying their connectivity k, 
#c) identify their relevance (i.e., kME) to these modules, and 
#d) identify the module functionality, and e) identify their direct linking partners (high correlation?)
#- a) I achieve by looking at k, b) I achieve by correlating the lncRNAs contained within the module to the module's first PC,
#- c) I achieve by testing overlap (should I do Fisher exact test?), 
#- d) I achieve through functional enrichment, and e) I achieve by exporting it to the networks

#0)
#Testing expression of lncRNAs in microglia
#Grab ncRNAs
#Grab covariate files
dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", "Internship", "Internship_Start", "Project 4 - Synthesis", "Analysis", "Genelists")
setwd(dir)

SCZ_ncRNAs <- unique(read_lines("ncRNAgenes_Gandal.csv"))
SZ_ncRNAs.annot <- annotate(SCZ_ncRNAs)
mic_specific_ncRNAs <- unique(read_lines("microglia-specific ncRNAs.csv"))
mic_specific_ncRNAs.annot <- annotate(mic_specific_ncRNAs)

dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", "Internship", "Internship_Start", "Project 4 - Synthesis", "Analysis")
setwd(dir)

#Test for expression of SCZ-related ncRNAs in microglia
dds_select <- dds[rownames(dds) %in% SCZ_ncRNAs, ]
vsd_select <- assay(vsd[rownames(vsd) %in% SCZ_ncRNAs, ]) #Assay transformed values
select_sub <- Subtrans_norm[names(Subtrans_norm) %in% SCZ_ncRNAs]
Expr.SZ_ncRNAs.annot <- annotate(colnames(select_sub))

#Expr.SZ_ncRNAs.annot <- c(Expr.SZ_ncRNAs.annot[,c(1,4)], colnames(select_sub)[!colnames(select_sub) %in% Expr.SZ_ncRNAs.annot[,1]])
#Expr.SZ_ncRNAs.annot <- c(Expr.SZ_ncRNAs.annot[,c(1,4)])

scaledrows_sub <- t(scale(select_sub))
SampleAnnot <- cbind(data.frame("Batch" = all_covs_sub[,1]), data.frame("Region" = all_covs_sub[,4]), data.frame("Diagnosis" = all_covs_sub[,5]))
rownames(SampleAnnot) <- rownames(Subtrans_norm)
pheatmap(scaledrows_sub, scale="none", show_cols=FALSE, show_rownames=TRUE, cluster_rows=TRUE, labels_row = Expr.SZ_ncRNAs.annot.t,
         cluster_cols=TRUE, annotation_col=SampleAnnot, fontsize = 9, main="Gandal et al. SZ-dysregulated ncRNA expression")



#==========================================================================================================#
#                       Noncoding RNA list creation and expression analysis                                #
#==========================================================================================================#

#Creating a select list of genes of interest (for noncoding RNAs)
allmatchedmicnc <- names(trans_norm) %in% mic_specific_ncRNAs
table(allmatchedmicnc)["TRUE"]
matchedmicnc <- names(Subtrans_norm) %in% mic_specific_ncRNAs 
table(matchedmicnc)["TRUE"]
namesMatchedMic <- names(Subtrans_norm[names(Subtrans_norm) %in% mic_specific_ncRNAs])
table(moduleColors[matchedmicnc]) 


SCZ_ncRNAs.Expr.sub <- colnames(hm.subset.df(Subtrans_norm, SCZ_ncRNAs, SampleAnnot, "test", get=TRUE))

#While plotting heatmaps, retrieve subsetted dataframe to check on correlation between lncRNAs (coregulation)
Top10kMicNC.Expr<- hm.subset.df(Subtrans_norm, mic_specific_ncRNAs, SampleAnnot, c("Lake et al. 'microglia-specific' ncRNA expression"), get=TRUE)
AllMicNC.Expr <- hm.subset.df(trans_norm, mic_specific_ncRNAs, SampleAnnot, c("Lake et al. 'microglia-specific' ncRNA expression"), get=TRUE)

Top10kMicNC.Cor <- t(data.frame(bicor(Top10kMicNC.Expr, use = "p")))
AllMicNC.Cor <- t(data.frame(bicor(AllMicNC.Expr, use = "p")))

ggcorrplot(Top10kMicNC.Cor, hc.order = TRUE, type = "upper",
           outline.col = "white", tl.cex = 7,
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"))

ggcorrplot(AllMicNC.Cor, hc.order = TRUE, type = "upper",
           outline.col = "white", tl.cex = 7,
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"))



#==========================================================================================================#
#                                     Noncoding RNA module enrichment                                      #
#==========================================================================================================#

#a)
#Figure what SCZ-related ncRNAs are expressed in what portion of the modules
selectNCSCZgenes <- Subtrans_norm[colnames(Subtrans_norm) %in% SCZ_ncRNAs] 
probes = names(Subtrans_norm)
MatchWithSelectNC <- is.finite(match(probes, SCZ_ncRNAs))
table(MatchWithSelectNC)["TRUE"] 
selectSCZNCMod = moduleColors[MatchWithSelectNC]
df <- data.frame(table(selectSCZNCMod))



#==========================================================================================================#
#                                       Connectivity extraction                                            #
#==========================================================================================================#

#b)
#Figure k value (connectivity) of these module-specific lncRNAs:
GeneNames.ByMod.List <- NULL
GeneNames.ByMod.List <- extract.GeneNames.PerMod(summarized.modClrs, Subtrans_norm, moduleColors, names(Subtrans_norm))


hubsfromadj <- intramodularConnectivity(adjacency, mergedColors, scaleByMax = FALSE)
#Grab module-specific k statistics
HubGenes <- hubsfromadj

####################Grab k-statistics####################

summarized.modClrs <- names(table(moduleColors))

k.per.mod <- list()
for (i in summarized.modClrs)
  k.per.mod[[i]] <- extract.k(HubGenes, GeneNames.ByMod.List[[i]])

SCZ_ncRNA.GeneNames.PerMod <- extract.GeneNames.PerMod(summarized.modClrs, Subtrans_norm, moduleColors, SCZ_ncRNAs)

SCZ_ncRNA.k <- extract.k.perMod(summarized.modClrs, k.per.mod, SCZ_ncRNA.GeneNames.PerMod)


#==========================================================================================================#
#                                 MM extraction  (for all + SCZ_ncRNAs)                                    #
#==========================================================================================================#

#b.1)
#More meaningful, however, is the intramodular inspection on modular hub genes instead of generally highly interconnected genes
#Similarly one may measure direct correlation of given lncRNAs of a module with a given module by using:
#Trait-module correlation; using kME to compare networks (kME = module membership)

#Use the "select" data frame containing lncRNA expression per sample to test
#its centrality to certain modules by correlating lncRNA expression to module eigengene values, 
#thereby identifying modules in which any of the given select genes is central.

MMofallGenes <- data.frame(bicor(Subtrans_norm, MEs, use="p"))

####################Grab MM####################


All.GeneNames.PerMod <- extract.GeneNames.PerMod(summarized.modClrs, Subtrans_norm, moduleColors, probes)

All.MM <- NULL
for (i in summarized.modClrs){
  All.MM[[i]] <- get.MM(MMofallGenes, All.GeneNames.PerMod[[i]], MEmod=paste("ME", i, sep=""), color=i)
}

#For schizophrenia non-coding RNAs
#MMofSCZNCGenes = data.frame(bicor(select_sub, MEs, use="p")) #function signedKME gives same output

SCZ_ncRNAs.MM <- NULL
for (i in summarized.modClrs){
  SCZ_ncRNAs.MM[[i]] <- get.MM(MMofallGenes, SCZ_ncRNA.GeneNames.PerMod[[i]], MEmod=paste("ME", i, sep=""), color=i)
}



#==========================================================================================================#
#                                    Intramodular analysis                                                 #
#==========================================================================================================#

#c.2)
#Intramodular analysis:

#Choosing top hub genes
hubs = chooseTopHubInEachModule(Subtrans_norm, moduleColors, power=6, type="unsigned", omitColors=NA)
data.frame(hubs)

hubs.k <- NULL
for (i in summarized.modClrs){
  hubs.k[[i]] <- HubGenes[rownames(HubGenes) %in% hubs[[i]],]
}

##################Grab expression by gene list###################

hubs.Expr <- NULL
for (i in summarized.modClrs){
  hubs.Expr[[i]] <- Subtrans_norm[colnames(Subtrans_norm) %in% rownames(hubs.k[[i]])]
}
  


#==========================================================================================================#
#                           Feeding back: summary of DE genes in SZ                                        #
#==========================================================================================================#

dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", "Internship", "Internship_Start", "Project 4 - Synthesis", "Analysis", "Genelists")
setwd(dir)

DEGenes.Gandal <- read.table("DEGenes.txt", header=T, stringsAsFactors=FALSE)
DE.pcGenes.Gandal <- DEGenes.Gandal[DEGenes.Gandal$gene_type == "protein_coding",]
DE.pcGenes.Gandal$SCZ.fdr <- gsub(",",".",DE.pcGenes.Gandal$SCZ.fdr)
sign.DE.pcGenes.Gandal <- DE.pcGenes.Gandal[DE.pcGenes.Gandal$SCZ.fdr < 0.05,]$ensembl_gene_id
sign.DE.pcGenes.Gandal.sub.Expr <- colnames(hm.subset.df(Subtrans_norm, sign.DE.pcGenes.Gandal, annotation=SampleAnnot, title="sign.DE.pcGenes.Gandal", get=TRUE))

dysregSZRNAs <- list("sign.DE.pcRNAs.Gandal"=sign.DE.pcGenes.Gandal.sub.Expr, "sign.DE.ncRNAs"=SCZ_ncRNAs.Expr.sub)

dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", "Internship", "Internship_Start", "Project 4 - Synthesis", "Analysis")
setwd(dir)

dysreg.SZncRNA.SZpcRNA <- getitall(dysregSZRNAs, summarized.modClrs, Subtrans_norm, moduleColors, k.per.mod, MMofallGenes, hubs.Expr) #Function can be found below

##########################
#Several more goals from here on:
#1. Overarching:
#  1. Fisher's exact test to figure which list gene enrichments are statistically meaningful (i.e., module contains significantly many ncRNA that would imply a functional role of this module in SCZ)
#	2. Test for STAT1
#2. Finish primary goal no. 1:
#	1. t-test on k-Out statistics on SCZ dysregulated ncRNAs (test whether different from mean to make a statement about statistical implication of kOut value; i.e., significantly relevant to stability of module)
#3. Finish primary goal no. 2: Grab additionally TWAS prioritized genes from Gandal et al. 
#	1. identify k-statistics of SCZ-associated SNPs within modules (intramodular analysis of SCZ-related genetic risk loci)
#	2. correlate SCZ-associated SNP ncRNAs with module hub gene + ME
# 3. Identify other non-coding RNA SNP genes within microglia and their values
#	4. identify module enrichment of SCZ-associated SNP protein-coding genes
#		1. identify relevance to module (module hub / ME correlation + k-statistics); don't know whether important right now?
#		2. identify ncRNA partners within module (see STable 6):
#			1. correlate SCZ-associated SNP genes with module-specific (SCZ-dysregulated) ncRNA genes
#				1. test whether any of these ncRNAs are SNP enriched genes
#4. Functional annotation of modules:
#	1. IPA / gprofileR
#	2. GO-enrichment
#5. Consensus analysis:
#	1. Verify results in a) 15 healthy subset samples, b) Gijsje's & Katia's dataset
#6. Visualization:
#	1. Graphia



#==========================================================================================================#
#                         Analysing modules on all SNP-associated lncRNAs                                  #
#==========================================================================================================#

#3) Repeating the same process as in a-d (focused on dysregulated ncRNAs) for all SNP-associated noncoding RNAs

#Grab all necessary files
#Annotated risk loci genes from SCZ CLOZUK study
dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", "Internship", "Internship_Start", "Project 4 - Synthesis", "Analysis", "FUMA_CLOZUK")
setwd(dir)

CLOZUK_SNPs <- read.table("genes.txt", header=T, stringsAsFactors=FALSE)[c(1:2)]
CLOZUK_SNPs_allinfo <- read.table("genes.txt", header=T, stringsAsFactors=FALSE)
filteredCLOZUK <- CLOZUK_SNPs_allinfo[CLOZUK_SNPs_allinfo$minGwasP < 5*10^-8,]

table(mic_specific_ncRNAs %in% CLOZUK_SNPs[,1]) #None enriched for SNPs

dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", "Internship", "Internship_Start", "Project 4 - Synthesis", "Analysis")
setwd(dir)



#==========================================================================================================#
#                               Network summary for SCZ_ncRNA.SNPs                                         #
#==========================================================================================================#

########################3.1) Identify module statistics for all noncoding SNP genes#######################
allnoncoding <- filteredCLOZUK[,1][filteredCLOZUK[7] != "protein_coding"]

ncSNP.Expr <- hm.subset.df(Subtrans_norm, allnoncoding, SampleAnnot, c("SNP-associated lncRNAs"), get=TRUE)

ModuleAnalysis.GeneEnrichment(Subtrans_norm, allnoncoding, moduleColors)



#==========================================================================================================#
#                                k and MM for  prioritized SCZ_ncRNA.SNPs                                  #
#==========================================================================================================#

#Grab the TWAS prioritized ones
dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", "Internship", "Internship_Start", "Project 4 - Synthesis", "Analysis", "Genelists")
setwd(dir)

TWAS_SNPs <- read.table("TWAS_Gandal.txt", header=T, stringsAsFactors=FALSE) #Filtered on TWAS FDR < 0.05
TWAS_SNPs_corrected <- read_lines("TWAS_signGenesAfterFilter.csv")
TWAS_SNPs_evid <- read.table("TWAS_Gandal_multievid.txt", header=T, stringsAsFactors=FALSE) #Evidenced by SMR + HEIDI in addition

dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", "Internship", "Internship_Start", "Project 4 - Synthesis", "Analysis")
setwd(dir)

#Grab the prioritized ones from this:
SCZ_ncRNA.PrioSNPs <- TWAS_SNPs[1][TWAS_SNPs[7] != "protein_coding"]
SCZ_ncRNA.PrioSNPs.Corrected <- SCZ_ncRNA.PrioSNPs[SCZ_ncRNA.PrioSNPs %in% TWAS_SNPs_corrected] #23


########################3.2.0) Check expression#######################
SCZ_ncRNA.PrioSNPs.Corrected.Expr <- hm.subset.df(Subtrans_norm, SCZ_ncRNA.PrioSNPs.Corrected, SampleAnnot, c("Prioritized SNP-associated lncRNAs (schizophrenia) from Gandal et al. 2018 TWAS"), get=TRUE)

########################3.2.1) Grab module-specific prioritized noncoding SNP gene names#######################
#Check module enrichment
ModuleAnalysis.GeneEnrichment(Subtrans_norm, SCZ_ncRNA.PrioSNPs, moduleColors)
ModuleAnalysis.GeneEnrichment(Subtrans_norm, SCZ_ncRNA.PrioSNPs.Corrected, moduleColors)

#Automating it for every gene list I have
allgls <- list("DE ncRNAs"=SCZ_ncRNAs, "SNP-assoc. ncRNAs"=allnoncoding, "TWAS-prio. ncRNAs"=SCZ_ncRNA.PrioSNPs.Corrected)


#Grab gene names by module for every list 
ModEnrichment.ByMod.ByGeneList <- NULL
ModEnrichment.ByMod.ByGeneList <- GeneNames.ByMod.PerGeneList(allgls, summarized.modClrs, Subtrans_norm, moduleColors)



#==========================================================================================================#
#                            k and MM for evidenced prioritized SCZ_ncRNA.SNPs                             #
#==========================================================================================================#

#3.3) For the evidenced prioritized SNPs (multiple methods indicating these noncoding genes to be expression loci effected by SNP)
#Just grab it from the prioritized lists:
EvidGenes <- TWAS_SNPs_evid[1][TWAS_SNPs_evid[9] < 0.05]
SCZ_ncRNA.EvidPrioSNPs.SNPs <- EvidGenes[EvidGenes %in% SCZ_ncRNA.PrioSNPs]
SCZ_ncRNA.EvidPrioSNPs.Corrected.Expr <- hm.subset.df(Subtrans_norm, SCZ_ncRNA.EvidPrioSNPs.SNPs, SampleAnnot, c("Evidenced prioritized SNP-associated lncRNAs (schizophrenia) from Gandal et al. 2018 TWAS"), get=TRUE)
SCZ_ncRNA.PrioSNPs.Corrected %in% SCZ_ncRNA.EvidPrioSNPs.Corrected.Expr
#Module analysis
ModuleAnalysis.GeneEnrichment(Subtrans_norm, SCZ_ncRNA.EvidPrioSNPs.SNPs, moduleColors)

#3.3.1) Grab module specific names of evidenced prioritized genes

SCZ_ncRNA.EvidPrioSNPs.GeneNames <- extract.GeneNames.PerMod(summarized.modClrs, Subtrans_norm, moduleColors, SCZ_ncRNA.EvidPrioSNPs.SNPs)

allnoncoding.Expr <- colnames(hm.subset.df(Subtrans_norm, allnoncoding, SampleAnnot, "core ID Expr.", get=TRUE))
SCZ_ncRNA.PrioSNPs.Expr <- colnames(hm.subset.df(Subtrans_norm, SCZ_ncRNA.PrioSNPs, SampleAnnot, "core ID Expr.", get=TRUE))
SCZ_ncRNA.PrioSNPs.Corrected.Expr <- colnames(hm.subset.df(Subtrans_norm, SCZ_ncRNA.PrioSNPs.Corrected, SampleAnnot, "core ID Expr.", get=TRUE))
SCZ_ncRNA.EvidPrioSNPs.SNPs.Expr <- colnames(hm.subset.df(Subtrans_norm, SCZ_ncRNA.EvidPrioSNPs.SNPs, SampleAnnot, "core ID Expr.", get=TRUE))


SZ_ncRNA_SNPs.gls <- list("SNP-assoc. ncRNAs"=allnoncoding.Expr, "TWAS-prio. ncRNAs"=SCZ_ncRNA.PrioSNPs.Corrected.Expr)

ModEnrichment.ByMod.ByGeneList <- GeneNames.ByMod.PerGeneList(allgls, summarized.modClrs, Subtrans_norm, moduleColors)

#Summary with getitall:
SZncRNA.SNPs.res <- getitall(SZ_ncRNA_SNPs.gls, summarized.modClrs, Subtrans_norm, moduleColors, k.per.mod, MMofallGenes, hubs.Expr) #Function can be found below



#==========================================================================================================#
#                   Analysis of SNP-associated SZ dysreg. ncRNAs (SNP-dysreg overlap)                      #
#==========================================================================================================#

#3.4) Get specific for those dysregulated noncoding RNAs that overlap with SNP noncoding genes (all, prioritized, evidenced)
#CLOZUK SNP genes overlapping with SCZ dysregulated ncRNA genes:
CLOZUKgenes <- filteredCLOZUK[,1]
SCZncSNPsAllGenes <- filteredCLOZUK[filteredCLOZUK[,1] %in% SCZ_ncRNAs,] #21 of the SNP loci are within the SCZ ncRNAs (with GWAS P at < 5*10^-8)

#Gandal's TWAS prioritized ones:
SZprioncRNAs <- TWAS_SNPs[TWAS_SNPs[,1] %in% SCZ_ncRNAs,] #7 ncRNAs from our dataset that are prioritized SNP genes in Gandal's TWAS; two of which are within the corrected list (also the two we identified)
SZprioncRNAs.corrected <- TWAS_SNPs_corrected[TWAS_SNPs_corrected %in% SCZ_ncRNAs] #2 ncRNAs from our dataset that are prioritized SNP genes in Gandal's TWAS; two of which are within the corrected list (also the two we identified)
#Evidenced by multiple methods in Gandal's TWAS:
SCZncSNPprio <- EvidGenes[EvidGenes %in% SCZ_ncRNAs] #1 ncRNA that is a prioritized SNP in Gandal's TWAS (evidenced by multiple methods); this is the LINC00634 mentioend by Gandal's paper

#3.4.1) noncoding SNP overlap with dysregulated ncRNA gene overlap
dysSCZ_ncRNA.SNPs <- list("DE SNP-assoc. ncRNAs"=SCZncSNPsAllGenes[,1],"DE TWAS-prio ncRNAs"=SZprioncRNAs[,1])
dysSCZ_ncRNA.SNPs.ModEnrichment <- NULL
#dysSCZ_ncRNA.SNPs.ModEnrichment <- setNames(dysSCZ_ncRNA.SNPs.ModEnrichment, paste(names(dysSCZ_ncRNA.SNPs)))
for (i in dysSCZ_ncRNA.SNPs){
  hm.subset.df(Subtrans_norm, i, SampleAnnot, title=paste(names(dysSCZ_ncRNA.SNPs)), get=FALSE)
  dysSCZ_ncRNA.SNPs.ModEnrichment[[length(dysSCZ_ncRNA.SNPs.ModEnrichment)+1]] <- ModuleAnalysis.GeneEnrichment(Subtrans_norm, i, moduleColors)
}

#For only the CLOZUK + uncorrected prioritized we find some back in our dataset to be relevant
SCZncSNPsAllGenes.Expr <- colnames(hm.subset.df(Subtrans_norm, SCZncSNPsAllGenes[,1], SampleAnnot, "Expressed DE ncRNA SNP genes", get=TRUE))
SZprioncRNAs.Expr <- colnames(hm.subset.df(Subtrans_norm, SZprioncRNAs[,1], SampleAnnot, "Expressed DE ncRNA prio SNP genes", get=TRUE))

dysSCZ_ncRNA.SNPs <- list("DE SNP-assoc. ncRNAs"=SCZncSNPsAllGenes.Expr,"DE TWAS-prio ncRNAs"=SZprioncRNAs.Expr)

#Results for these SNP-associated ncRNAs (k, MM, Hub correlation)
dysSCZ_ncRNA.SNPs.res <- getitall(dysSCZ_ncRNA.SNPs, summarized.modClrs, Subtrans_norm, moduleColors, k.per.mod, MMofallGenes, hubs.Expr)

#To clear this up: the CLOZUK SNP gene list and Gandal's TWAS give slightly different results



#==========================================================================================================#
#                              Analysis of SNP-associated protein-coding RNAs                              #
#==========================================================================================================#

#3.5) Identify position of protein-coding SNP genes
#Test for protein-coding SNP genes, since one of the aims is to figure by what ncRNAs these are regulated within modules
#3.5.0.1) Grab all protein-coding SNP genes
SCZ_pcRNA.SNPs <- filteredCLOZUK[,1][filteredCLOZUK[7] == "protein_coding"] #703



#3.5.0.2) Grab prioritized protein-coding SNP-enriched genes
SCZ_pcRNA.prioSNPs <- TWAS_SNPs[,1][TWAS_SNPs[7] == "protein_coding"]
SCZ_pcRNA.prioSNPs.corrected <- TWAS_SNPs_corrected[TWAS_SNPs_corrected %in% SCZ_pcRNA.prioSNPs]

SCZ_pcRNA.prioSNPs.corrected %in% dysregSZRNAs$sign.DE.pcRNAs.Gandal
SCZ_pcRNA.EvidprioSNPs %in% dysregSZRNAs$sign.DE.pcRNAs.Gandal
#3.5.0.3) Grab evidenced prioritized protein-coding SNP-enriched genes
SCZ_pcRNA.EvidprioSNPs <- EvidGenes[EvidGenes %in% SCZ_pcRNA.prioSNPs.corrected] #25 in the evidence prioritized SNP genes

#Test for expression:
SZ_pcRNA.SNPs <- list("SNP-assoc. mRNAs"=SCZ_pcRNA.SNPs,#"SZ prio. SNP pcRNAs"=SCZ_pcRNA.prioSNPs, #leaving out this one for now
                          "TWAS-prio mRNAs"=SCZ_pcRNA.prioSNPs.corrected)

SZ_pcRNA.SNPs.Names <- list()
SZ_pcRNA.SNPs.ModEnrichment <- list()
SZ_pcRNA.SNPs.Expr <- list()
for (f in names(SZ_pcRNA.SNPs)){
  SZ_pcRNA.SNPs.Expr[[f]] <- colnames(hm.subset.df(Subtrans_norm, SZ_pcRNA.SNPs[[f]], SampleAnnot, title=paste(names(SZ_pcRNA.SNPs)), get=TRUE))
  SZ_pcRNA.SNPs.ModEnrichment[[f]] <- ModuleAnalysis.GeneEnrichment(Subtrans_norm, SZ_pcRNA.SNPs[[f]], moduleColors)
  SZ_pcRNA.SNPs.Names[[f]] <- extract.GeneNames.PerMod(summarized.modClrs, Subtrans_norm, moduleColors, SZ_pcRNA.SNPs[[f]])
}


#3.5.1.) Identify module relevance (ME/hub correlation, k-statistics) ; not important for now, will leave for later
SZ_pcRNA.SNPs.res <- getitall(SZ_pcRNA.SNPs, summarized.modClrs, Subtrans_norm, moduleColors, k.per.mod, MMofallGenes, hubs.Expr)



#==========================================================================================================#
#                                                 Cleaning                                                 #
#==========================================================================================================#

dysreg.SZncRNA.SZpcRNA.res.clean <- list(list.clean(dysreg.SZncRNA.SZpcRNA$k.perMod.perGL, summarized.modClrs), list.clean(dysreg.SZncRNA.SZpcRNA$MM.perMod.perGL, summarized.modClrs), cleanup.l3(dysreg.SZncRNA.SZpcRNA$`Expr + HubCor.perMod.perGL`, summarized.modClrs))
dysreg.SZncRNA.SZpcRNA.res.clean  <- setNames(dysreg.SZncRNA.SZpcRNA.res.clean , names(dysreg.SZncRNA.SZpcRNA))

SZncRNA.SNPs.res.clean <- list(list.clean(SZncRNA.SNPs.res$k.perMod.perGL, summarized.modClrs), list.clean(SZncRNA.SNPs.res$MM.perMod.perGL, summarized.modClrs), cleanup.l3(SZncRNA.SNPs.res$`Expr + HubCor.perMod.perGL`, summarized.modClrs))
SZncRNA.SNPs.res.clean  <- setNames(SZncRNA.SNPs.res.clean , names(SZncRNA.SNPs.res))

dysSCZ_ncRNA.SNPs.res.clean <- list(list.clean(dysSCZ_ncRNA.SNPs.res$k.perMod.perGL, summarized.modClrs), list.clean(dysSCZ_ncRNA.SNPs.res$MM.perMod.perGL, summarized.modClrs), cleanup.l3(dysSCZ_ncRNA.SNPs.res$`Expr + HubCor.perMod.perGL`,summarized.modClrs))
dysSCZ_ncRNA.SNPs.res.clean  <- setNames(dysSCZ_ncRNA.SNPs.res.clean , names(dysSCZ_ncRNA.SNPs.res))

SZ_pcRNA.SNPs.res.clean <- list(list.clean(SZ_pcRNA.SNPs.res$k.perMod.perGL, summarized.modClrs), list.clean(SZ_pcRNA.SNPs.res$MM.perMod.perGL, summarized.modClrs), cleanup.l3(SZ_pcRNA.SNPs.res$`Expr + HubCor.perMod.perGL`, summarized.modClrs))
SZ_pcRNA.SNPs.res.clean  <- setNames(SZ_pcRNA.SNPs.res.clean , names(SZ_pcRNA.SNPs.res))



##################################################MODULE ANALYSIS##################################################################

#==========================================================================================================#
#                                  Gene set enrichment analysis (GSEA)                                     #
#==========================================================================================================#

#Enrichment analysis
exprSCZ_ncRNAs <- colnames(selectNCSCZgenes)
#SCZ_ncRNA.mod.enrich <- GSEA.byMod(GeneNames.ByMod.List, exprSCZ_ncRNAs, universe=20000, get=TRUE)

mic.core.mod.enrich <- GSEA.byMod(GeneNames.ByMod.List, coreIDs[,1], universe=20000, get=TRUE)

coreIDs.Expr <- colnames(hm.subset.df(Subtrans_norm, coreIDs[,1], SampleAnnot, "core ID Expr.", get=TRUE))
broadCoreIDs.Expr <- colnames(hm.subset.df(Subtrans_norm, broadCoreIDs[,1], SampleAnnot, "broadCoreIDs Expr.", get=TRUE))
AllIDs.Expr <- colnames(hm.subset.df(Subtrans_norm, AllIDs[,1], SampleAnnot, "AllIDs Expr.", get=TRUE))

allMicGenes <- list("Microglia marker genes"=AllIDs.Expr)

list <- list("Microglia marker genes"=AllIDs.Expr,
             "DE ncRNAs"=dysregSZRNAs$sign.DE.ncRNAs, 
             "SNP-assoc. ncRNAs"=SZ_ncRNA_SNPs.gls$`SNP-assoc. ncRNAs`, "TWAS prio. ncRNAs"=SZ_ncRNA_SNPs.gls$`TWAS-prio. ncRNAs`,
             #             "DE SNP-assoc. ncRNAs"=dysSCZ_ncRNA.SNPs$`DE SNP-assoc. ncRNAs`, "DE TWAS-prio. ncRNAs" = dysSCZ_ncRNA.SNPs$`DE TWAS-prio ncRNAs`,
             "DE mRNAs"=dysregSZRNAs$sign.DE.pcRNAs.Gandal, 
             "SNP-assoc. mRNAs"=SZ_pcRNA.SNPs.Expr$`SNP-assoc. mRNAs`, "TWAS-prio. mRNAs"=SZ_pcRNA.SNPs.Expr$`TWAS-prio mRNAs`)


GSEA.1 <- list()
for (i in names(list))
  GSEA.1[[i]] <- GSEA.byMod(GeneNames.ByMod.List, list[[i]], universe=20000, get=TRUE)


######FOR LIST OVERLAP########

names(list)
#To test for list overlap
GSEA.l <- list()
for (i in names(list))
  for (n in names(list))
    GSEA.l[[paste(i, "x", n)]] <- data.frame(GSEA.byList(list[[i]], list[[n]], get=TRUE, universe=20000))

GSEA.lists <- NULL
GSEA.lists <-  GSEA.byList(list[[i]], list[[n]], get=TRUE, universe=20000)
for (i in names(list))
  for (n in names(list))
    GSEA.lists <- rbind(GSEA.lists, GSEA.byList(list[[i]], list[[n]], get=TRUE, universe=20000))

GSEA.lists <- GSEA.lists[-1,]

rownames(GSEA.lists) <- names(GSEA.l)
#length <- seq(from = 1, to = length(rownames(GSEA.lists)), by = 6)
#GSEA.lists.overlap <- NULL
#for (i in length)
#  GSEA.lists.overlap <- cbind(GSEA.lists.overlap, GSEA.lists[length(GSEA.lists.overlap):i,][3])

GSEA.lists.overlap <- NULL
GSEA.lists.overlap <- GSEA.lists[1:7,][3]
GSEA.lists.overlap <- cbind(GSEA.lists.overlap, GSEA.lists[8:14,][3])
GSEA.lists.overlap <- cbind(GSEA.lists.overlap, GSEA.lists[15:21,][3])
GSEA.lists.overlap <- cbind(GSEA.lists.overlap, GSEA.lists[22:28,][3])
GSEA.lists.overlap <- cbind(GSEA.lists.overlap, GSEA.lists[29:35,][3])
GSEA.lists.overlap <- cbind(GSEA.lists.overlap, GSEA.lists[36:42,][3])
GSEA.lists.overlap <- cbind(GSEA.lists.overlap, GSEA.lists[43:49,][3])

rownames(GSEA.lists.overlap) <- names(list)
colnames(GSEA.lists.overlap) <- names(list)

GSEA.lists.p <- NULL
GSEA.lists.p <- GSEA.lists[1:7,][5]
GSEA.lists.p <- cbind(GSEA.lists.p, GSEA.lists[8:14,][5])
GSEA.lists.p <- cbind(GSEA.lists.p, GSEA.lists[15:21,][5])
GSEA.lists.p <- cbind(GSEA.lists.p, GSEA.lists[22:28,][5])
GSEA.lists.p <- cbind(GSEA.lists.p, GSEA.lists[29:35,][5])
GSEA.lists.p <- cbind(GSEA.lists.p, GSEA.lists[36:42,][5])
GSEA.lists.p <- cbind(GSEA.lists.p, GSEA.lists[43:49,][5])
rownames(GSEA.lists.p) <- names(list)
colnames(GSEA.lists.p) <- names(list)


pheatmap(GSEA.lists.p)

rowlabels <- NULL
for (i in names(list))
  rowlabels <- c(rowlabels, paste(i, "(", length(list[[i]]), ")", sep=""))


pheatmap(as.dist(GSEA.lists.p), cluster_rows = F, cluster_cols = F, annotation_legend = FALSE, 
         legend_breaks = c(0, 0.2 , 0.4, 0.6, 0.8, 1, max(GSEA.lists.p)), legend = T, main = "GSEA: cross gene list enrichment",
         labels_cols = rowlabels, labels_row = rowlabels, labels_col=rowlabels, legend_labels = c("0", "0.2", "0.4", "0.6", "0.8", "1", "p-val\n\n"),
         color = colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(50), display_numbers = T)

pheatmap(GSEA.lists.overlap, cluster_rows = F, cluster_cols = F, annotation_legend = FALSE, 
         legend_breaks = c(0, 100, 200, 600, 1000, max(GSEA.lists.overlap)), legend = T, main = "GSEA: cross gene list enrichment",
         labels_cols = rowlabels, labels_row = rowlabels, labels_col=rowlabels, legend_labels = c("0", "100", "200", "600", "1000", "overlap\n\n"),
         color = colorRampPalette(c("white", "mediumpurple1"))(50), display_numbers = T)



######################EXPERIMENTING#######################################
# VISUALIZE ENRICHMENT
GSEA <- data.frame("empty"=length(summarized.modClrs))
  
for (i in names(GSEA.1))
  GSEA <- cbind(GSEA, GSEA.1[[i]]$GSEA["overlap"])

GSEA <- GSEA[!colnames(GSEA) %in% "empty"]

colnames(GSEA) <- names(list)

rowAnnot <- NULL
for (i in summarized.modClrs)
  rowAnnot <- c(rowAnnot, paste(i, "", "(", length(GeneNames.ByMod.List[[i]]), ")", sep=""))
#rownames(GSEA) <- rowAnnot

#Grab pvalues
pval.overlap <- data.frame("empty"=length(summarized.modClrs))

for (i in names(GSEA.1))
  pval.overlap <- cbind(pval.overlap, GSEA.1[[i]]$GSEA["p.value"])

pval.overlap <- pval.overlap[!colnames(pval.overlap) %in% "empty"]

colnames(pval.overlap) <- names(list)
######################EXPERIMENTING#######################################
annot.save <- annot
mod.annot <- data.frame("Module Colors"=summarized.modClrs)

annot_clrs <- list()
for (i in summarized.modClrs)
  annot_clrs$Module.Colors[[i]] = i

rowAnnot <- data.frame("Module Colors"=rowAnnot)
rownames(GSEA) <- rowAnnot[,1]
rownames(mod.annot) <- rowAnnot[,1]
rownames(pval.overlap)<- rowAnnot[,1]

pheatmap(t(GSEA), cluster_rows = F, cluster_cols = F, annotation_legend = FALSE, 
         annotation_col=mod.annot, annotation_colors = annot_clrs, labels_row = rowlabels,
         color = colorRampPalette(c("white", "mediumpurple1"))(50), display_numbers = TRUE)

pheatmap(t(GSEA), cluster_rows = F, cluster_cols = F, annotation_legend = FALSE, 
         legend_breaks = c(0, 50, 100, 150, 200, max(GSEA)), legend = T, main = "GSEA: gene list enrichment by module",
         legend_labels = c("0", "50", "100", "150", "200", "overlap\n\n"),
         annotation_col=mod.annot, annotation_colors = annot_clrs, labels_row = rowlabels,
         color = colorRampPalette(c("white", "mediumpurple1"))(50), display_numbers = T)


grid.text(1:40, x=rep(seq(0.05, 0.91, length.out=10), 10), 
          y=rep(seq(0, 1, 0.1)+0.05, each=10))
qnorm(.95)

pheatmap(t(pval.overlap), cluster_rows = F, cluster_cols = F, annotation_legend = FALSE, 
         annotation_col=mod.annot, annotation_colors = annot_clrs, annotation_row=rowAnnot,
         color = colorRampPalette(c("white", "yellowgreen"))(50), display_numbers = TRUE)

pheatmap(t(pval.overlap), cluster_rows = F, cluster_cols = F, annotation_legend = FALSE, 
         legend_breaks = c(0, 0.2 , 0.4, 0.6, 0.8, 1, max(GSEA.lists.p)), legend = T, main = "GSEA module x genelist significance of overlap",
         legend_labels = c("0", "0.2", "0.4", "0.6", "0.8", "1", "p-val\n\n"),
         annotation_col=mod.annot, annotation_colors = annot_clrs, labels_row = rowlabels,
         color = colorRampPalette(c("white", "mediumpurple1"))(50), display_numbers = T)

#How about making a heatmap for each of the lists separately?
pheatmap(t(pval.overlap), cluster_rows = F, cluster_cols = F, annotation_legend = FALSE, 
         legend_breaks = c(0, 0.2 , 0.4, 0.6, 0.8, 1, max(GSEA.lists.p)), legend = T, main = "GSEA module x genelist significance of overlap",
         legend_labels = c("0", "0.2", "0.4", "0.6", "0.8", "1", "p-val\n\n"),
         annotation_col=mod.annot, annotation_colors = annot_clrs, labels_row = rowlabels,
         color = colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(50), display_numbers = T)


#Red, brown, blue might be of interest (prioritized SNP ncRNA); red and brown also contain prio. pcRNA SNPs (evidenced); 
#while these modules also contain dysregulated ncRNAs and considerable amounts of microglia core signature genes
#Let's annotate them hence!
#==========================================================================================================#
#                              gprofiler module GO-term enrichment                                         #
#==========================================================================================================#
#ANNOTATING RELEVANT MODULES
#Create gene annotation
genenames.allgenes <- colnames(Subtrans_norm)
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl") #Human genes (GRCh38.p13)
genemap.allgenes <- getBM(values = genenames.allgenes,
                          filters = "ensembl_gene_id",
                          mart = ensembl,
                          attributes = c("ensembl_gene_id", "entrezgene_id",
                                         "hgnc_symbol", "external_gene_name",
                                         "description", "chromosome_name",
                                         "strand"))



annotations <- list()
for (i in names(GeneNames.ByMod.List))
  annotations[[i]] <- genemap.allgenes[,4][genemap.allgenes[,1] %in% GeneNames.ByMod.List[[i]]]


#Using gprofiler
res.annot <- NULL
for (i in names(annotations))
  res.annot[[i]] <- gost(annotations[[i]], organism = "hsapiens", ordered_query = FALSE, 
                                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                           measure_underrepresentation = FALSE, evcodes = TRUE, 
                                           user_threshold = 0.05, correction_method = "fdr", 
                                           domain_scope = "annotated", custom_bg = NULL, 
                                           numeric_ns = "", sources = NULL, as_short_link = FALSE)

res <- list()
for (i in names(GeneNames.ByMod.List))
  res[[i]] <- gost(GeneNames.ByMod.List[[i]], organism = "hsapiens", ordered_query = FALSE, 
                               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                               measure_underrepresentation = FALSE, evcodes = TRUE, 
                               user_threshold = 0.05, correction_method = "fdr", 
                               domain_scope = "annotated", custom_bg = NULL, 
                               numeric_ns = "", sources = NULL, as_short_link = FALSE)

#Let's grab for red, brown, blue
#res <- list()
#for (i in c("red", "brown", "blue"))
#  res[[i]] <- gost(GeneNames.ByMod.List[[i]], organism = "hsapiens", ordered_query = FALSE, 
#                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
#                   measure_underrepresentation = FALSE, evcodes = TRUE, 
#                   user_threshold = 0.05, correction_method = "g_SCS", 
#                   domain_scope = "annotated", custom_bg = NULL, 
#                   numeric_ns = "", sources = NULL, as_short_link = FALSE)

tab.res <- list()
for (i in names(res)){
  tab.res[[i]] <- res[[i]]$result[,c("term_id", "term_name", "p_value", "intersection")]
}

for (i in res)
  plot <- gostplot(i, capped = TRUE, interactive = TRUE) 
  print(plot)


#Works perfectly fine
for (i in names(tab.res)){
  plot <- GO.barplot(tab.res[[i]])
  print(plot)
}

for (i in c("brown", "cyan", "darkgrey", "green", "paleturquoise", "midnightblue", "yellow"))
  write.xlsx(tab.res[[i]], paste("GO results for",i,"module.xlsx"), row.names=T)


#==========================================================================================================#
#                                           Gene selection                                                 #
#==========================================================================================================#

###########################################For kWithin / kTotal criterion:################################################
length(GeneNames.ByMod.List$darkgrey)
#SZ.evidprio.pcRNA.SNPs.relevant.Genes <- kCrit.name.extract(summarized.modClrs, SZ_pcRNA.SNPs.res.clean$k.perMod.perGL$`SZ evid. prio. SNP pcRNAs`, GeneNames.ByMod.List)

k.perMod.perGL.list <- list("dysreg.SZncRNA.SZpcRNA.res.clean$k.perMod.perGL"=dysreg.SZncRNA.SZpcRNA.res.clean$k.perMod.perGL,
                            "SZncRNA.SNPs.res$k.perMod.perGL"=SZncRNA.SNPs.res.clean$k.perMod.perGL,
                            "dysSCZ_ncRNA.SNPs.res.clean$k.perMod.perGL"=dysSCZ_ncRNA.SNPs.res.clean$k.perMod.perGL,
                            "SZ_pcRNA.SNPs.res.clean$k.perMod.perGL"=SZ_pcRNA.SNPs.res.clean$k.perMod.perGL)

kCrit.byList <- list()
for (i in names(k.perMod.perGL.list))
  for (n in names(k.perMod.perGL.list[[i]]))
    kCrit.byList[[n]] <- kCrit.name.extract(sum.modClrs=summarized.modClrs, list=k.perMod.perGL.list[[i]][[n]], mod.name.list=GeneNames.ByMod.List)

###########################################For MM criterion:################################################

MM.perMod.perGL.list <- list("dysreg.SZncRNA.SZpcRNA.res.clean$MM.perMod.perGL"=dysreg.SZncRNA.SZpcRNA.res.clean$MM.perMod.perGL,
                             "SZncRNA.SNPs.res$MM.perMod.perGL"=SZncRNA.SNPs.res.clean$MM.perMod.perGL,
                             "dysSCZ_ncRNA.SNPs.res.clean$MM.perMod.perGL"=dysSCZ_ncRNA.SNPs.res.clean$MM.perMod.perGL,
                             "SZ_pcRNA.SNPs.res.clean$MM.perMod.perGL"=SZ_pcRNA.SNPs.res.clean$MM.perMod.perGL)

MMCrit.byList <- list()
for (i in names(MM.perMod.perGL.list))
  for (n in names(MM.perMod.perGL.list[[i]]))
    MMCrit.byList[[n]] <- MMCrit.name.extract(summarized.modClrs, MM.perMod.perGL.list[[i]][[n]], GeneNames.ByMod.List)

#If interested in checking in on whether correct:
#table(MM.perMod.perGL.list$`dysreg.SZncRNA.SZpcRNA.res.clean$MM.perMod.perGL`$sign.DE.pcRNAs.Gandal$yellow$`yellow Module Membership (MM)` > 0.6)


###########################################For hub correlation criterion:################################################

HubCor.ByMod.byGL.list <- list("dysreg.SZncRNA.SZpcRNA.res.clean$`Expr + HubCor.perMod.perGL`$HubCor.ByMod.byGL"=dysreg.SZncRNA.SZpcRNA.res.clean$`Expr + HubCor.perMod.perGL`$HubCor.ByMod.byGL,
                               "SZncRNA.SNPs.res$`Expr + HubCor.perMod.perGL`$HubCor.ByMod.byGLS"=SZncRNA.SNPs.res.clean$`Expr + HubCor.perMod.perGL`$HubCor.ByMod.byGL,
                               "dysSCZ_ncRNA.SNPs.res.clean$`Expr + HubCor.perMod.perGL`$HubCor.ByMod.byGL"=dysSCZ_ncRNA.SNPs.res.clean$`Expr + HubCor.perMod.perGL`$HubCor.ByMod.byGL,
                               "SZ_pcRNA.SNPs.res.clean$`Expr + HubCor.perMod.perGL`$HubCor.ByMod.byGL"=SZ_pcRNA.SNPs.res.clean$`Expr + HubCor.perMod.perGL`$HubCor.ByMod.byGL)

HubCorCrit.byList <- list()
for (i in names(HubCor.ByMod.byGL.list))
  for (n in names(HubCor.ByMod.byGL.list[[i]]))
    HubCorCrit.byList[[n]] <- HubCorCrit.name.extract(summarized.modClrs, HubCor.ByMod.byGL.list[[i]][[n]], GeneNames.ByMod.List)


#This one works best as it grabs ALL gene names (and doesn't replace them each time when genes are in the same module)
stringent <- list()
for (i in names(HubCorCrit.byList))
  for (n in summarized.modClrs){
    for (z in HubCorCrit.byList[[i]][[n]])
      for (k in kCrit.byList[[i]]$`kWithin criterion`[[n]])
        for (l in MMCrit.byList[[i]][[n]])
          if (k %in% z && k %in% l)
              stringent[[i]] <- unique(c(stringent[[i]], k))}

#If wanting to check whether correct:
#names <- NULL
#for (i in kCrit.byList$sign.DE.pcRNAs.Gandal$`kWithin criterion`)
#  names <- c(names, i)

#names.MM <- NULL
#for (i in MMCrit.byList$sign.DE.pcRNAs.Gandal)
#  names.MM <- c(names.MM, i)

#names %in% names.MM

kCrit.byList$sign.DE.pcRNAs.Gandal$`kWithin criterion`
#Grabbing specific values for these:
k.hub.MM.stringent <- getitall(stringent, summarized.modClrs, Subtrans_norm, moduleColors, k.per.mod, MMofallGenes, hubs.Expr)

k.hub.MM.stringent.clean <- list(list.clean(k.hub.MM.stringent$k.perMod.perGL, summarized.modClrs), list.clean(k.hub.MM.stringent$MM.perMod.perGL, summarized.modClrs), cleanup.l3(k.hub.MM.stringent$`Expr + HubCor.perMod.perGL`, summarized.modClrs))
k.hub.MM.stringent.clean  <- setNames(k.hub.MM.stringent.clean , names(k.hub.MM.stringent))

k.stringent <- NULL
for (i in names(k.hub.MM.stringent.clean$k.perMod.perGL))
  for (n in names(k.hub.MM.stringent.clean$k.perMod.perGL[[i]]))
    for (k in rownames(k.hub.MM.stringent.clean$k.perMod.perGL[[i]][[n]])){
      if (!k %in% rownames(k.stringent))
        k.stringent <- rbind(k.stringent, k.hub.MM.stringent.clean$k.perMod.perGL[[i]][[n]])}
#      else k.stringent <- rbind(k.stringent, NA)}
#Remove all doublets that were ignored by R in the if statement
k.stringent <- k.stringent[!rownames(k.stringent) %in% gsub('\\b\\w{1,15}\\b','',rownames(k.stringent)),]
stringent.names <- rownames(k.stringent)

k.hub.MM.stringent.clean$k.perMod.perGL$`SZ SNP pcRNAs`$brown

#For Hub correlation criterion
mods.stringent <- NULL
for (i in names(k.hub.MM.stringent.clean$k.perMod.perGL))
  mods.stringent <- c(mods.stringent, names(k.hub.MM.stringent.clean$k.perMod.perGL[[i]]))

mods.stringent <- unique(mods.stringent)
hubs.stringent <- hubs[names(hubs) %in% mods.stringent]

stringent.hubcor <- bicor(Subtrans_norm[colnames(Subtrans_norm) %in% stringent.names], Subtrans_norm[colnames(Subtrans_norm) %in% paste(hubs.stringent)])


######################Heatmap construction############################
colours <- c("Hub gene" = "indianred3", 
             
             "Microglia signature gene" = "mediumpurple1", #Needs to be more subtle
             
             "Hub gene & Microglia signature gene" = "red2",
             
             
             #This one is cool "#C6E021FF"
             
             "DE ncRNA" = "#E8E419FF",
             "Hub gene & DE ncRNA" = "orange",
             
             "SNP-assoc. ncRNA"="#C6E021FF",
             "TWAS-prio. ncRNA" = "#DDE318FF",
             "SNP-assoc. ncRNA & Microglia signature gene" = "mediumvioletred",
             "Non-module-relevant SNP-assoc. ncRNA & Microglia signature gene" = "violetred",
             
             
             "DE mRNA" = "lightblue1",
             "SNP-assoc. mRNA" = "skyblue",
             "TWAS-prio. mRNA" = "dodgerblue2",
             
             "SNP-assoc. mRNA & DE mRNA" = "deepskyblue1",
             
             "Hub gene & DE mRNA" = "mediumorchid2",

             "Hub gene & SNP-assoc. mRNA" = "#3A538BFF",
             "Hub gene & TWAS-prio. mRNA" = "#404588FF",
             
             "DE mRNA & Microglia signature gene" = "purple",
             "SNP-assoc. mRNA & Microglia signature gene" = "purple3",
             "TWAS-prio. mRNA & Microglia signature gene" = "magenta",
             "Non-module-relevant TWAS-prio. mRNA & Microglia signature gene" = "darkorange",
             "Non-module-relevant DE mRNA & Microglia signature gene" = "lightyellow",
             "Non-module-relevant SNP-assoc. mRNA & Microglia signature gene" = "pink"
             )


#Generating an annotation
stringent.pheatmap.annot <- NULL
for (i in names(k.hub.MM.stringent.clean$k.perMod.perGL))
  for (n in names(k.hub.MM.stringent.clean$k.perMod.perGL[[i]]))
    for (z in stringent.names)
      if (z %in% rownames(k.hub.MM.stringent.clean$k.perMod.perGL[[i]][[n]])){
        if (z %in% stringent.pheatmap.annot[,1])
          stringent.pheatmap.annot[stringent.pheatmap.annot[,1] %in% z,] <- c(z, n)
        else
          stringent.pheatmap.annot <- rbind(stringent.pheatmap.annot, c(z, n))
      }

#stringent.pheatmap.annot <- as.factor(stringent.pheatmap.annot[1:length(stringent.names)]) #Shortening for doublets

stringent.pheatmap.annot <- data.frame("GeneIDs"=stringent.pheatmap.annot[,1],"Module.Colors"=stringent.pheatmap.annot[,2])
#colnames(stringent.pheatmap.annot) <- "Module.Colors"
#stringent.pheatmap.annot <- data.frame(module.Colors.all=as.factor(stringent.pheatmap.annot))
rownames(stringent.pheatmap.annot) <- stringent.names


stringent.pheatmap.annot.hubs <- NULL
for (i in names(k.hub.MM.stringent.clean$k.perMod.perGL))
  for (n in names(k.hub.MM.stringent.clean$k.perMod.perGL[[i]]))
    if (!n %in% stringent.pheatmap.annot.hubs)
      stringent.pheatmap.annot.hubs <- c(stringent.pheatmap.annot.hubs, n)

annotated <- NULL
for (i in stringent.pheatmap.annot.hubs)
  annotated <- rbind(annotated, c(i, hubs[[i]]))

hubs.annot <- annotate(paste(hubs))

hubs.annot2 <- hubs.annot[hubs.annot$ensembl_gene_id %in% paste(hubs),]
hubs.annot2 <- data.frame("GeneIDs"=hubs.annot2$ensembl_gene_id, "GeneNames"=hubs.annot2$external_gene_name)
top.hubs <- hubs.annot2[hubs.annot2$GeneIDs %in% paste(hubs.stringent),]

stringent.pheatmap.annot.hubs <- top.hubs
#colnames(stringent.pheatmap.annot.hubs) <- "Module.Colors.hubs"

missing <- hubs[hubs %in% paste(hubs[paste(hubs) %in% colnames(stringent.hubcor)])[!paste(hubs[paste(hubs) %in% colnames(stringent.hubcor)]) %in% stringent.pheatmap.annot.hubs$GeneIDs]]

rownames(stringent.pheatmap.annot.hubs.names)[is.na(rownames(stringent.pheatmap.annot.hubs.names))] <- names(missing)
rownames(stringent.pheatmap.annot.hubs.names)[rownames(stringent.pheatmap.annot.hubs.names) %in% c("blue", "red")] <- paste(missing)

rnames <- rownames(stringent.pheatmap.annot.hubs.names)
rnames[c(4)] <- paste("blue")
rnames[c(10)] <- paste("red")
rownames(stringent.pheatmap.annot.hubs.names) <- rnames

#sorted.hubs.annot <- stringent.pheatmap.annot.hubs.names[sort(rownames(stringent.pheatmap.annot.hubs.names)),]
#sorted.hubs.annot$GeneIDs <- paste(hubs)[names(hubs) %in% rownames(sorted.hubs.annot)]

names <- as.character(stringent.pheatmap.annot.hubs.names$GeneNames)
names[is.na(names)] <- paste(hubs[names(hubs) %in% c("blue", "red")])
stringent.pheatmap.annot.hubs.names$GeneNames <- names

stringent.pheatmap.annot.hubs.names.ordered <- as.character(sorted.hubs.annot$GeneNames)

stringent.pheatmap.annot.names <- data.frame(unique(stringent.pheatmap.annot.names.all[,4]))
rownames(stringent.pheatmap.annot.names) <- unique(stringent.pheatmap.annot.names.all[,1])
annot.fornow <- cbind(stringent.pheatmap.annot.names, unique(stringent.pheatmap.annot.names.all[,1]))
stringent.pheatmap.annot.names.ordered <- annot.fornow[match(rownames(stringent.hubcor), rownames(annot.fornow)),]
colnames(stringent.pheatmap.annot.names.ordered)[1] <- "Genenames"


stringent.pheatmap.annot.rest <- NULL
for (i in names(k.hub.MM.stringent.clean$k.perMod.perGL))
  for (n in names(k.hub.MM.stringent.clean$k.perMod.perGL[[i]]))
    for (z in rownames(stringent.pheatmap.annot))
      if (z %in% rownames(k.hub.MM.stringent.clean$k.perMod.perGL[[i]][[n]])){
        if (z %in% stringent.pheatmap.annot.rest[,2])
          stringent.pheatmap.annot.rest[stringent.pheatmap.annot.rest[,2] %in% z] <- c(i, z)
        else 
          stringent.pheatmap.annot.rest <- rbind(stringent.pheatmap.annot.rest, c(i, z))
      }
  



stringent.pheatmap.annot.type <- stringent.pheatmap.annot.rest[!duplicated(stringent.pheatmap.annot.rest[,2]),]
stringent.pheatmap.annot.type[,1][11:12] <- "SZ evid. prio. SNP pcRNAs"

#Using the output of the network visualization function (last part of the analysis below) for heatmap annotation purposes
accum[,1] <- as.character(accum[,1])#Can be found below; output from network visualization code part
rownames(accum) <- accum[,1]
colnames(accum) <- c("Gene.Name", "RNA.type")
pheat.annot <- data.frame(accum)[!colnames(accum) == "Gene.Name"]
#clrs <- paste(viridis(18))
#clrs.seq <- clrs[seq(1, length(clrs), 2)]
RNA.type <- list("RNA.type"=viridis(40)[c(1,10,19:23,31:40)])
names(RNA.type$RNA.type) <- paste(as.character(unique(paste(as.character(accum[,2])))))

colours <- colours[names(colours) %in% unique(paste(as.character(accum[,2])))] #Is from the network graph part of the code

RNA.type <- list("RNA.type"=colours)
names(RNA.type$RNA.type) <- as.character(names(RNA.type$RNA.type))


stringent.pheatmap.annot.rest.type <- data.frame(RNA.type=as.factor(stringent.pheatmap.annot.type[,1]))
rownames(stringent.pheatmap.annot.rest.type) <- stringent.pheatmap.annot.type[,2]


rownames(stringent.hubcor) <- stringent.pheatmap.annot.names.ordered$Genenames
colnames(stringent.hubcor) <- stringent.pheatmap.annot.hubs.names$GeneNames
#rownames(stringent.pheatmap.annot.hubs) <- stringent.pheatmap.annot.hubs.names.ordered
new.for.annot <- data.frame("GeneIDs"=stringent.pheatmap.annot.hubs.names$GeneIDs, "Module.Colors"=rownames(stringent.pheatmap.annot.hubs.names))
rownames(new.for.annot) <- stringent.pheatmap.annot.hubs.names$GeneNames
rowannots.stringent.pheatmap <- accum[rownames(accum) %in% rownames(stringent.hubcor),]
#Sort them for pheatmap
rowannots.stringent.pheatmap <- rowannots.stringent.pheatmap[match(rownames(stringent.pheatmap.annot), rownames(rowannots.stringent.pheatmap)),]
#The solution is to not match the gene names with colors, but their modules.
pheatmap(stringent.hubcor, cluster_rows =T, cluster_cols = T, annotation_legend = F, 
         legend_breaks = c(-1,-0.8,-0.6,-0.4,-0.2, 0, 0.2 , 0.4, 0.6, 0.8, 1, max(stringent.hubcor)), legend = T, main = "Correlation of stringently selected genes (vertical) with respective module hub genes (horizontal)",
         legend_labels = c("-1","-0.8","-0.6","-0.4","-0.2","0", "0.2", "0.4", "0.6", "0.8", "1", "Bi-Midweight correlation \n \n"),
         annotation_row=cbind(stringent.pheatmap.annot[!colnames(stringent.pheatmap.annot)=="GeneIDs"],
                              rowannots.stringent.pheatmap[colnames(rowannots.stringent.pheatmap) != "Gene.Name"]), 
         annotation_colors = c(MMs.stringent.annot.clrs, RNA.type), #using the color annotation from below; annot_clrs from above works,too
         annotation_col= new.for.annot[colnames(new.for.annot)!="GeneIDs"],
         color = colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(50), display_numbers = F,
         width=9, height=10, cellwidth=10, cellheight=10)

rownames(rowannots.stringent.pheatmap) %in% rownames(stringent.pheatmap.annot)
#########################################For MM-criterion###################################
for (i in stringent.pheatmap.annot.hubs)
  MEnames.stringent <- paste("ME", i, sep="")

MMs.stringent <- bicor(Subtrans_norm[colnames(Subtrans_norm) %in% stringent.names], MEs[colnames(MEs) %in% MEnames.stringent])

stringent.pheatmap.annot.MMs <- stringent.pheatmap.annot.hubs
colnames(stringent.pheatmap.annot.MMs) <- "ME.Colors"

rownames(stringent.pheatmap.annot.MMs) <- MEnames.stringent


###########Heatmap construction##############
MMs.stringent.annot <- NULL
for (i in stringent.pheatmap.annot$Module.Colors)
  MMs.stringent.annot <- c(MMs.stringent.annot, i)


MMs.stringent.annot <- data.frame("Module Colors"=MMs.stringent.annot)
rownames(MMs.stringent.annot) <- rownames(MMs.stringent)
annot <- data.frame("Module Colors"=summarized.modClrs)


MMs.stringent.annot.clrs <- list()
for (i in stringent.pheatmap.annot$Module.Colors)
  MMs.stringent.annot.clrs$Module.Colors[[i]] = i
#Column annotation just for spatial properties of the legend in the heatmap
MM.stringent.column.annot <- stringent.pheatmap.annot.hubs
colnames(MM.stringent.column.annot)
rownames(MM.stringent.column.annot) <- MEnames.stringent


rownames(MMs.stringent) <- stringent.pheatmap.annot.names.ordered$Genenames
stringent.pheatmap.annot.rest.type.check <- stringent.pheatmap.annot.rest.type
rownames(stringent.pheatmap.annot.rest.type.check) <- stringent.pheatmap.annot.names.ordered[match(rownames(stringent.pheatmap.annot.rest.type), rownames(stringent.pheatmap.annot.names.ordered)),][,1]
rownames(stringent.pheatmap.annot) <- stringent.pheatmap.annot.names.ordered[match(rownames(stringent.pheatmap.annot.rest.type), rownames(stringent.pheatmap.annot.names.ordered)),][,1]
stringent.pheatmap.annot.names.ordered[rownames(stringent.pheatmap.annot.rest.type) %in% rownames(stringent.pheatmap.annot.names.ordered),]


#The solution is to not match the gene names with colors, but their modules.
pheatmap(MMs.stringent, cluster_rows = T, cluster_cols = T, annotation_legend = T, 
         legend_breaks = c(-0.6,-0.4,-0.2, 0, 0.2 , 0.4, 0.6, 0.8,  max(MMs.stringent)), legend = T, main = "Correlation of stringently selected genes (vertical) with respective module eigengene (horizontal) showing module membership (kME)",
         legend_labels = c("-0.6","-0.4","-0.2","0", "0.2", "0.4", "0.6", "0.8", "Bi-Midweight correlation \n \n"),
         annotation_row=cbind(stringent.pheatmap.annot[!colnames(stringent.pheatmap.annot)=="GeneIDs"],
                              rowannots.stringent.pheatmap[colnames(rowannots.stringent.pheatmap) != "Gene.Name"]), 
         annotation_col=MM.stringent.column.annot,
         annotation_colors=c(MMs.stringent.annot.clrs, RNA.type), labels_row=as.character(stringent.pheatmap.annot.names.ordered$Genenames),
         color = colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(50), display_numbers = FALSE,
         width=9, height=10, cellwidth=10, cellheight=10)

write.xlsx(MMs.stringent, "MMs.stringent.xlsx", row.names = T)

#################################k-criterion###################################
#Making a datatable for k criterion
k.stringent.genenames <- annot.fornow[match(rownames(k.stringent), rownames(annot.fornow)),]
colnames(k.stringent.genenames)[1] <- "Genenames"

k.stringent.complete <- cbind(k.stringent, stringent.pheatmap.annot$Module.Colors, stringent.pheatmap.annot.rest.type, k.stringent.genenames$Genenames)

colnames(k.stringent.complete)[c(9,11)] <- c("module.Color","Gene.Name")

library(xlsx)
write.xlsx(k.stringent.complete, "k.Criterion_unsigned_all.xlsx", row.names = T) #Table is correct now


#Annotationby RNA type


for (i in names(list))
  str.annot.1[[i]] <- annotate(list[[i]])

str.annot.2 <- list()
for (i in names(str.annot.1))
  str.annot.2[[i]] <- str.annot.1[[i]][,4]

str.annot <- list()
for (i in names(str.annot.2))
  for (n in k.stringent.complete$Gene.Name){
    if (n %in% str.annot.2[[i]]){
      str.annot[[i]] <- c(str.annot[[i]], n)}}

names(str.annot) <- c("Microglia signature gene", "DE ncRNA", "SNP-assoc. ncRNA", "DE mRNA", "SNP-assoc. mRNA", "TWAS-prio. mRNA")



#==========================================================================================================#
#                  Analysis of non-coding - protein-coding co-expression patterns                          #
#==========================================================================================================#

#Grab for every protein-coding SNP-assocaited RNA x dysreg. SNP-associatedRNA interaction the correlation value
SZ_ncRNA.SNPs.GeneExpr.ByMod <- GrabExpr.ByMod.ByGLs(dysSCZ_ncRNA.SNPs, summarized.modClrs) #Only for the CLOZUK + uncorrected prioritized genes (only genes present)
pcRNA.pcRNA.CoReg <- interact.CoReg(SZ_pcRNA.SNPs.GeneExpr.ByMod, SZ_pcRNA.SNPs.GeneExpr.ByMod)

ncRNA.pcRNA.CoReg <- interact.CoReg(SZ_ncRNA.SNPs.GeneExpr.ByMod,SZ_ncRNA.SNPs.GeneExpr.ByMod)

dysregSZ_ncRNAs.GeneExpr.ByMod <- GrabExpr.ByMod.ByGLs(allgls, summarized.modClrs) #Only for the CLOZUK + uncorrected prioritized genes (only genes present)
dysreg.SZ_ncRNA.CoReg <- interact.CoReg(dysregSZ_ncRNAs.GeneExpr.ByMod, dysregSZ_ncRNAs.GeneExpr.ByMod)

micGene.GeneExpr.ByMod <- Expr.byMod.ByGL(summarized.modClrs, AllIDs[,1])


###################################FOR SELECTED GENES###################################
#For stringently selected gene list:
#########################################pcRNAs#########################################
relevant.SZ.pcRNAs <- list("relevant SZ evid. prio. SNP pcRNAs"=stringent$`SZ evid. prio. SNP pcRNAs`, "relevant sign.DE.pcRNAs.Gandal"=stringent$sign.DE.pcRNAs.Gandal)
relevant.SZ.pcRNAs.Expr <- list.clean(GrabExpr.ByMod.ByGLs(relevant.SZ.pcRNAs, summarized.modClrs), summarized.modClrs)
relevant.SZ.pcRNAs.micCor <- CoReg.byMod.byGL(relevant.SZ.pcRNAs.Expr$`relevant SZ evid. prio. SNP pcRNAs`, micGene.GeneExpr.ByMod) #You have to set the gene list with the fewer modules in front

relevant.SZ.SNP.MicGenes.Expr <- list()
for (i in names(relevant.SZ.pcRNAs.Expr$`relevant SZ evid. prio. SNP pcRNAs`))
  relevant.SZ.SNP.MicGenes.Expr[[i]] <- c(relevant.SZ.pcRNAs.Expr$`relevant SZ evid. prio. SNP pcRNAs`[[i]], micGene.GeneExpr.ByMod[[i]])

relevant.SZ.SNP.micCor.sqr <- CoReg.byMod.byGL(relevant.SZ.SNP.MicGenes.Expr, relevant.SZ.SNP.MicGenes.Expr) #You have to set the gene list with the fewer modules in front

#########################################ncRNAs#########################################
#Using the liberal SNP-associated genes
relevant.SZ.ncRNAs <- list("sign.DE.ncRNAs.Gandal"=stringent$sign.DE.ncRNAs, "SNP-assoc. ncRNAs"=stringent$`SNP-assoc. ncRNAs`)


#How about I use all: liberal, liberal and then highlight the stringent ones in the plot
relevant.SZ.pcRNAs.liberal <- list("sign.DE.pcRNAs.Gandal"=liberal$sign.DE.pcRNAs.Gandal, "SZ evid. prio. SNP pcRNAs"=liberal$`SZ evid. prio. SNP pcRNAs`)
relevant.SZ.ncRNAs.liberal <- list("sign.DE.ncRNAs"=liberal$sign.DE.ncRNAs, "SZ evid. prio. SNP ncRNAs"=liberal$`evid. prio. SNP-assoc. ncRNAs`)
#liberal.2 <- list("SZ.pcRNAs.liberal"=relevant.SZ.pcRNAs.liberal, "SZ.ncRNAs.liberal"=relevant.SZ.ncRNAs.liberal) #None of the pcRNAs fall into the same modules as ncRNAs

liberal.Expr.2 <- list.clean(GrabExpr.ByMod.ByGLs(liberal, summarized.modClrs), summarized.modClrs)
stringent.Expr.2 <- list.clean(GrabExpr.ByMod.ByGLs(stringent, summarized.modClrs), summarized.modClrs)

#Rather this way:
liberal.Expr.safe <- list()
for (i in names(liberal.Expr.2))
  for (n in names(liberal.Expr.2[[i]]))
    liberal.Expr.safe[[n]][[i]] <- liberal.Expr.2[[i]][[n]]

stringent.Expr.safe <- list()
for (i in names(stringent.Expr.2))
  for (n in names(stringent.Expr.2[[i]]))
    stringent.Expr.safe[[n]][[i]] <- stringent.Expr.2[[i]][[n]]


library(purrr)
reduce <- purrr::reduce
liberal.Expr <- list()
  for (i in names(liberal.Expr.safe))
    liberal.Expr[[i]] <- reduce(liberal.Expr.safe[[i]], cbind)

stringent.Expr <- list()
for (i in names(stringent.Expr.safe))
  stringent.Expr[[i]] <- reduce(stringent.Expr.safe[[i]], cbind)

#Remove duplicates since all lists come together conjoint
for (i in names(liberal.Expr))
  liberal.Expr[[i]] <- liberal.Expr[[i]][!duplicated(colnames(liberal.Expr[[i]]))]
liberalSZRNA.Mic.Hub.Expr <- list()
for (i in names(liberal.Expr))
   liberalSZRNA.Mic.Hub.Expr[[i]] <- cbind(liberal.Expr[[i]], micGene.GeneExpr.ByMod[[i]], hubs.Expr[[i]])

for (i in names(stringent.Expr))
  stringent.Expr[[i]] <- stringent.Expr[[i]][!duplicated(colnames(stringent.Expr[[i]]))]
stringentSZRNA.Mic.Hub.Expr <- list()
for (i in names(stringent.Expr))
  stringentSZRNA.Mic.Hub.Expr[[i]] <- cbind(stringent.Expr[[i]], micGene.GeneExpr.ByMod[[i]], hubs.Expr[[i]])


#Remove duplicates again
for (i in names(liberalSZRNA.Mic.Hub.Expr))
  liberalSZRNA.Mic.Hub.Expr[[i]] <- liberalSZRNA.Mic.Hub.Expr[[i]][!duplicated(colnames(liberalSZRNA.Mic.Hub.Expr[[i]]))]

for (i in names(stringentSZRNA.Mic.Hub.Expr))
  stringentSZRNA.Mic.Hub.Expr[[i]] <- stringentSZRNA.Mic.Hub.Expr[[i]][!duplicated(colnames(stringentSZRNA.Mic.Hub.Expr[[i]]))]


#This is rather correct now
liberalSZRNA.Hub.Mic.CoReg <- list()
for (i in names(liberalSZRNA.Mic.Hub.Expr))
    liberalSZRNA.Hub.Mic.CoReg[[i]] <- BiCor.GenesOfInterest(data.frame(liberalSZRNA.Mic.Hub.Expr[[i]]), data.frame(liberalSZRNA.Mic.Hub.Expr[[i]]))

liberalSZRNA.Hub.Mic.CoReg

stringentSZRNA.Hub.Mic.CoReg <- list()
for (i in names(stringentSZRNA.Mic.Hub.Expr))
  stringentSZRNA.Hub.Mic.CoReg[[i]] <- BiCor.GenesOfInterest(data.frame(stringentSZRNA.Mic.Hub.Expr[[i]]), data.frame(stringentSZRNA.Mic.Hub.Expr[[i]]))


#ANNOTATION FOR VISUALIZATION#

annot <- list()
annot2 <- list()
for (i in names(liberalSZRNA.Hub.Mic.CoReg)){
  annot[[i]] <- annotate(rownames(liberalSZRNA.Hub.Mic.CoReg[[i]]))
  for (n in rownames(liberalSZRNA.Hub.Mic.CoReg[[i]]))
    if (!n %in% annot[[i]]$ensembl_gene_id)
      annot2[[i]] <- c(annot2[[i]], n)
    else if (n %in% annot[[i]]$ensembl_gene_id)
      annot2[[i]] <- c(annot2[[i]], annot[[i]][annot[[i]]$ensembl_gene_id %in% n,]$external_gene_name )}

for (i in names(annot2)){
  if (length(rownames(liberalSZRNA.Hub.Mic.CoReg[[i]])) == length(annot2[[i]])){
    rownames(liberalSZRNA.Hub.Mic.CoReg[[i]]) <- annot2[[i]]
    colnames(liberalSZRNA.Hub.Mic.CoReg[[i]]) <- annot2[[i]]}
  else (annot2[[i]] <- annot2[[i]][1:length(rownames(liberalSZRNA.Hub.Mic.CoReg[[i]]))])
    rownames(liberalSZRNA.Hub.Mic.CoReg[[i]]) <- annot2[[i]]
    colnames(liberalSZRNA.Hub.Mic.CoReg[[i]]) <- annot2[[i]]
}


annot.stri <- list()
annot2.stri <- list()
for (i in names(stringentSZRNA.Hub.Mic.CoReg)){
  annot.stri[[i]] <- annotate(rownames(stringentSZRNA.Hub.Mic.CoReg[[i]]))
  annot.stri[[i]] <- annot.stri[[i]][!duplicated(annot.stri[[i]][,1]),]
  for (n in rownames(stringentSZRNA.Hub.Mic.CoReg[[i]]))
    if (!n %in% annot.stri[[i]]$ensembl_gene_id)
      annot2.stri[[i]] <- c(annot2.stri[[i]], n)
    else if (n %in% annot.stri[[i]]$ensembl_gene_id)
      annot2.stri[[i]] <- c(annot2.stri[[i]], annot.stri[[i]][annot.stri[[i]]$ensembl_gene_id %in% n,]$external_gene_name )}

for (i in names(annot2.stri)){
  if (length(rownames(stringentSZRNA.Hub.Mic.CoReg[[i]])) == length(annot2.stri[[i]])){
    rownames(stringentSZRNA.Hub.Mic.CoReg[[i]]) <- annot2.stri[[i]]
    colnames(stringentSZRNA.Hub.Mic.CoReg[[i]]) <- annot2.stri[[i]]}
  else (annot2.stri[[i]] <- annot2.stri[[i]][1:length(rownames(stringentSZRNA.Hub.Mic.CoReg[[i]]))])
  rownames(stringentSZRNA.Hub.Mic.CoReg[[i]]) <- annot2.stri[[i]]
  colnames(stringentSZRNA.Hub.Mic.CoReg[[i]]) <- annot2.stri[[i]]
}


#==========================================================================================================#
#                        Network visualization - Cluster analysis (correlation plot)                       #
#==========================================================================================================#
#These co-expression patterns then need to be visualized. That means: 
#1) correlation matrices, #DONE
#2) network graphs of clusters, 
#3) and term enrichment per cluster. 


#I think for now this one suffices, because I will show it directly next to the network
for (i in names(liberalSZRNA.Hub.Mic.CoReg))
  pheatmap(liberalSZRNA.Hub.Mic.CoReg[[i]], cluster_rows = F, cluster_cols = F, annotation_legend = FALSE, 
           legend_breaks = c(-1,-0.8,-0.6,-0.4,-0.2, 0, 0.2 , 0.4, 0.6, 0.8, 1, max(stringent.hubcor)), legend = T, main = paste("Correlation", i, "module genes"),
           legend_labels = c("-1","-0.8","-0.6","-0.4","-0.2","0", "0.2", "0.4", "0.6", "0.8", "1", "Biweight midcorrelation\n\n"),
           color = colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(50), display_numbers = FALSE)


#For stringent genes
#For stringent genes
for (i in names(stringentSZRNA.Hub.Mic.CoReg)){
  plot <- pheatmap(stringentSZRNA.Hub.Mic.CoReg[[i]], cluster_rows = T, cluster_cols = T, annotation_legend = T, 
                   legend_breaks = c(-1,-0.8,-0.6,-0.4,-0.2, 0, 0.2 , 0.4, 0.6, 0.8, 1, max(stringent.hubcor)), legend = T, main = paste("Correlation", i, "module genes"),
                   legend_labels = c("-1","-0.8","-0.6","-0.4","-0.2","0", "0.2", "0.4", "0.6", "0.8", "1", "Bi-Midweight correlation\n\n"), annotation_col = pheat.annot, annotation_row =pheat.annot,
                   color = colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(50), display_numbers = F, annotation_colors=RNA.type,
                   width=9, height=10, cellwidth=17, cellheight=17)
  print(plot)
}

dev.off() 



#==========================================================================================================#
#                                 Network visualization - Cluster network plot                              #
#==========================================================================================================#

type <- list()
for (i in names(str.annot))
    type[[i]]<- data.frame("ID"=unique(paste(str.annot[[i]])), "RNA type"=1:length(unique(paste(str.annot[[i]]))))

for (i in names(type))
    type[[i]][2] <- paste(i)

for (i in names(type))
    rownames(type[[i]]) <- unique(paste(str.annot[[i]]))

scale01 <- function(x){(x-min(x))/(max(x)-min(x))}

#for (z in names(stringentSZRNA.Hub.Mic.CoReg)){
#if (length(table(duplicated(rownames(stringentSZRNA.Hub.Mic.CoReg[[z]]))))==2){
#  for (i in rownames(stringentSZRNA.Hub.Mic.CoReg[[z]][duplicated(rownames(stringentSZRNA.Hub.Mic.CoReg[[z]])),]))
#    stringentSZRNA.Hub.Mic.CoReg[[z]][rownames(stringentSZRNA.Hub.Mic.CoReg[[z]]) %in% i,][1] <- paste(i, ".1", sep="")}
#if (length(table(rownames(stringentSZRNA.Hub.Mic.CoReg[[z]]) %in% "SAMSN1") == 2))
#  rownames(stringentSZRNA.Hub.Mic.CoReg[[z]])[rownames(stringentSZRNA.Hub.Mic.CoReg[[z]]) %in% "SAMSN1"] <- as.character(c("SAMSN1.1","SAMSN1.2"))
#  rownames(stringentSZRNA.Hub.Mic.CoReg[[z]])[rownames(stringentSZRNA.Hub.Mic.CoReg[[z]]) %in% "PAK6"] <- as.character(c("PAK6.1","PAK6.2"))
#  rownames(stringentSZRNA.Hub.Mic.CoReg[[z]])[rownames(stringentSZRNA.Hub.Mic.CoReg[[z]]) %in% "RPP21"] <- as.character(c("RPP21.1","RPP21.2"))
#}

######NOTE######
#ENSG00000179743 = AL450998.2 is denoted as FLJ37453 in the newest annotation
#Hence it is a module-relevant DE ncRNA
str.annot$`DE ncRNA` <- c(str.annot$`DE ncRNA`, "AL450998.2")

accum <- NULL
#    rownames(stringentSZRNA.Hub.Mic.CoReg[[z]][duplicated(rownames(stringentSZRNA.Hub.Mic.CoReg[[z]])),])
for (z in names(stringentSZRNA.Hub.Mic.CoReg)){
  
  mat <- as.matrix(as.dist(stringentSZRNA.Hub.Mic.CoReg[[z]]))
  mat[is.na(mat)] <- 0
  
  g <- graph.adjacency( #Generating a graph object of the genes within a module of interest
    as.matrix(mat),
    mode="undirected",
    weighted=TRUE,
    diag=FALSE)
  g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE) #Removing multiple and looped edges
  
  E(g)$weight <- abs(E(g)$weight)
  
  g <- delete_edges(g, E(g)[which(E(g)$weight<0.5)])
  
  mst <- mst(g, algorithm="prim")
  mst.edgeweights <- E(mst)$weight * 2.5


  expr <- cbind(micGene.GeneExpr.ByMod[[z]], stringent.Expr[[z]], hubs.Expr[[z]])
  expr <- expr[!duplicated(colnames(expr))]
  annot.expr <- annotate(colnames(expr))
  annot.expr.byorder <- annot.expr[match(colnames(expr), annot.expr[,1]),][,4]
  colnames(expr) <- annot.expr.byorder
  
  vSizes <- (scale01(apply(t(expr), 1, mean)) + 1.0) * 10 #Doing the same for node sizes so that these are relative to expression level
  
  
  df <- NULL
  for (n in names(type)){
    for (i in rownames(type[[n]]))
      for (k in rownames(stringentSZRNA.Hub.Mic.CoReg[[z]]))
        if(i %in% k){
          if (!i %in% rownames(df))
            df <- rbind(df, type[[n]][rownames(type[[n]]) %in% i,])
          if (i %in% rownames(df))
            df[rownames(df) %in% i,] <- type[[n]][rownames(type[[n]]) %in% i,]}}
  
  
  #  df <- NULL
  #  for (n in names(type)){
  #    if(is.null(rownames(type[[n]])) == FALSE)
  #      if(rownames(type[[n]]) %in% rownames(df) == FALSE)
  #        df <- rbind(df, type[[n]])}
  
  df2 <- NULL
  for (n in rownames(stringentSZRNA.Hub.Mic.CoReg[[z]]))
    if (n %in% rownames(df) == FALSE){
      df2 <- data.frame("ID"=rownames(stringentSZRNA.Hub.Mic.CoReg[[z]][(n %in% rownames(df) == FALSE), (n %in% rownames(df) == FALSE)]),
                        "RNA.type"=paste("Microglia signature gene"), stringsAsFactors = FALSE)
#      if (length(table(df2$ID %in% "SAMSN1")) == 2)
#        df2$ID[df2$ID %in% "SAMSN1"] <- c("SAMSN1", "SAMSN1.2")
      rownames(df2) <- df2$ID #else stringentSZRNA.Hub.Mic.CoReg[[z]][(n %in% rownames(df) == FALSE), (n %in% rownames(df) == FALSE)]
    }
  
  df3 <- rbind(df, df2[!rownames(df2) %in% rownames(df),])
  df3 <- df3[!duplicated(df3["ID"]),]
  
  for (i in rownames(df3))
    if (i %in% hubs.annot2$GeneNames)
      df3[rownames(df3) == i,]$RNA.type <- "Hub gene"
  
  
  for (i in rownames(df3))
    for (n in names(str.annot)){
      if (i %in% str.annot[[n]])
        df3[rownames(df3) == i,]$RNA.type <- n}
  
  
  for (i in rownames(df3))
    for (n in names(str.annot))
      for (k in names(str.annot)){
        if (i %in% str.annot[[n]] && i %in% str.annot[[k]] && k != n){
          if (df3[rownames(df3) == i,]$RNA.type != "TWAS-prio. mRNA & Microglia signature gene"){
            df3[rownames(df3) == i,]$RNA.type <- paste(n, "&", k)}
        }
        if (n != "Microglia signature gene" && i %in% str.annot.2[[paste(n,"s", sep="")]]&& i %in% str.annot.2$`Microglia marker genes`){
          if (df3[rownames(df3) == i,]$RNA.type != paste(n, "& Microglia signature gene")){
            df3[rownames(df3) == i,]$RNA.type <- paste("Non-module-relevant", n, "& Microglia signature gene")}}
        if (df3[rownames(df3) == i,]$RNA.type == "TWAS-prio. mRNA & SNP-assoc. mRNA"){
          df3[rownames(df3) == i,]$RNA.type <- "TWAS-prio. mRNA"}
        if (df3[rownames(df3) == i,]$RNA.type == "TWAS-prio. ncRNA & SNP-assoc. ncRNA"){
          df3[rownames(df3) == i,]$RNA.type <- "TWAS-prio. ncRNA"}
        }
      

  
  for (i in rownames(df3))
    for (n in names(str.annot))
      if (i %in% str.annot[[n]] && i %in% hubs.annot2$GeneNames)
        df3[rownames(df3) == i,]$RNA.type <- paste("Hub gene &", n)
  
  for (i in rownames(df3))
    if (i %in% str.annot$`Microglia signature gene` && i %in% as.character(hubs.annot2$GeneNames))
      df3[rownames(df3) == i,]$RNA.type <- paste("Hub gene &", "Microglia signature gene")
  

  accum <- rbind(accum, df3)
  
  #Just tersting out what colors one might use
  colours.sub <- colours[names(colours) %in% unique(df3$RNA.type)]
  
  node.color = df3[match(names(V(mst)), rownames(df3)),]$RNA.type

  net.plot <- 
    ggnet2(
      mst, label = TRUE, 
      mode = "fruchtermanreingold",
      layout.par = list(cell.jitter=0.75),
      label.size = 4, edge.size=mst.edgeweights,
      edge.label=round(E(mst)$weight, digits=2),
      size = vSizes, size.legend="Scaled degree of expression",
#      node.color = df3$RNA.type, color.legend="RNA type",
      node.color = node.color, color.legend="RNA type",
      legend.size=11, edge.label.size=4,
      palette = colours.sub)+
    guides(size = FALSE)+
    labs(title = paste(z,"module cluster analysis"))+
    theme(plot.title = element_text(hjust=0.5))
  print(net.plot)
}
###########################################################END##########################################################
