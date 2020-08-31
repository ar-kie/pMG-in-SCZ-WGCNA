#==========================================================================================================#
#                              pMG lncRNA-mRNA WGCNA - functions by section                                #
#==========================================================================================================#

#==========================================================================================================#
#                                    Filtering & outlier detection                                         #
#==========================================================================================================#
#Annotation of microglia signature genes
IDs <- function(values, object){
  genemap_micgenes <- getBM(values,
                            filters = "external_gene_name",
                            mart = ensembl,
                            attributes = c("ensembl_gene_id", "entrezgene_id",
                                           "hgnc_symbol", "external_gene_name",
                                           "description", "chromosome_name",
                                           "strand"))
  geneNames_data <- unique(data.frame(genemap_micgenes[c(1,4)]))
  present <- geneNames_data[geneNames_data[,1] %in% colnames(object),]
}



#==========================================================================================================#
#                          PCA and Variance Partititioning before correction                               #
#==========================================================================================================#

#Tweaking DESeq2 plotPCA.custom function
plotPCA.custom <- function (object, intgroup = "condition", ntop = 44099, returnData = TRUE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]), scale=TRUE)
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", 
                              color = "group")) + geom_point(size = 3) + xlab(paste0("PC1: ", 
                                                                                     round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2: ", 
                                                                                                                                              round(percentVar[2] * 100), "% variance")) + coord_fixed()
}


#Solely for annotation purposes of expressed genes
annotate <- function(list){
  genemap.allgenes <- getBM(list,
                            filters = "ensembl_gene_id",
                            mart = ensembl,
                            attributes = c("ensembl_gene_id", "entrezgene_id",
                                           "hgnc_symbol", "external_gene_name",
                                           "description", "chromosome_name",
                                           "strand"))
}



#==========================================================================================================#
#                                 Starting off the core analysis                                           #
#==========================================================================================================#

Expr.SZ_ncRNAs.annot.t <- NULL
for (i in colnames(select_sub)){
  if (!i %in% Expr.SZ_ncRNAs.annot$ensembl_gene_id)
    Expr.SZ_ncRNAs.annot.t <- c(Expr.SZ_ncRNAs.annot.t, i)
  else if(i %in% Expr.SZ_ncRNAs.annot$ensembl_gene_id)
    Expr.SZ_ncRNAs.annot.t <- c(Expr.SZ_ncRNAs.annot.t, Expr.SZ_ncRNAs.annot[Expr.SZ_ncRNAs.annot$ensembl_gene_id %in% i,]$external_gene_name)}




#==========================================================================================================#
#                       Noncoding RNA list creation and expression analysis                                #
#==========================================================================================================#

#Function for plotting heatmap on subsetted gene list dataframe, and if wanted returning subsetted dataframe
hm.subset.df <- function(df, g.list, annotation, title, get=FALSE){
  df.sub <- df[names(df) %in% g.list]
  df.sub.annot <- annotate(names(df.sub))
  Expr.annot.t <- NULL
  for (i in names(df.sub))
    if (!i %in% df.sub.annot$ensembl_gene_id)
      Expr.annot.t <- c(Expr.annot.t, i)
  else if(i %in% df.sub.annot$ensembl_gene_id)
    Expr.annot.t <- c(Expr.annot.t, df.sub.annot[df.sub.annot$ensembl_gene_id %in% i,]$external_gene_name)
  scaled.df <- t(scale(df.sub))
  pheatmap(scaled.df, scale="none", show_cols=FALSE, show_rownames=TRUE, cluster_rows=TRUE, labels_row=Expr.annot.t,
           cluster_cols=TRUE, annotation_col=annotation, fontsize = 9, main=title)
  if (get){
    return(df.sub)}
}



#==========================================================================================================#
#                                       Connectivity extraction                                            #
#==========================================================================================================#
####################Function to grab gene names and expression####################
extract.GeneNames.PerMod <- function(sum.ModClrs, df, modClrs, gl){
  ModGenes <- NULL
  getModIDs <- vector("list", length(sum.ModClrs))
  for (i in sum.ModClrs){
    ModGenes[[i]] <- names(df)[modClrs==i]
  }
  getModIDs <- list() #To avoid an error like "too many elements supplied" or even more basic "replacement null", make sure R knows that it is supposed to be a list, otherwise it will replace the list with NULL
  for (i in sum.ModClrs){
    getModIDs[[i]] <- ModGenes[[i]][ModGenes[[i]] %in% gl]
  }
  return(getModIDs)
}


extract.Expr. <- function(df,mdclrs,color){
  df[mdclrs==color]
}

####################Function to grab k-statistics####################
extract.k <- function(df.k, modGeneNames){
  rank <- cbind("kTotalRank"=rank(-df.k[,1]),"kOutRank.all"=rank(-df.k[,3]))
  total.k <- cbind(df.k, rank)
  rownames(total.k) <- rownames(df.k)
  select <- total.k[rownames(total.k) %in% modGeneNames,]
  rank.select <- cbind("kWithinRank"=rank(-select[,2]),"kOutRank.select"=rank(-select[,3]))
  rank.and.kstats <- cbind(select, rank.select)
  return(rank.and.kstats)
}

extract.k.perMod <- function(sum.modClrs, hubdf, genelist){
  k <- list()
  for (i in sum.modClrs)
    if (length(hubdf[[i]]) > 0)
      k[[i]] <- hubdf[[i]][rownames(hubdf[[i]]) %in% genelist[[i]],]
    else return(NULL)
    return(k)
}


####################Functions to grab hub gene correlations####################
Expr.byMod.ByGL <- function(sumModClrs, gl){
  ExprList <- NULL
  for (i in summarized.modClrs){
    ExprList[[i]] <- extractExpr.GL(i, gl)
  }
  #  ExprList <- ExprList[lapply(ExprList, length) > 0]
  return(ExprList)
}

ExprAndHubCor.byMod.byGL.all <- function(all.gls, sumModClrs, HubsExpr){
  Expr.byMod.ByGL.List <- list()
  for (f in names(all.gls)){
    Expr.byMod.ByGL.List[[f]] <- Expr.byMod.ByGL(sumModClrs, all.gls[[f]])
  }
  HubCor.byMod.byGL.List <- list()
  for (n in names(Expr.byMod.ByGL.List)){
    HubCor.byMod.byGL.List[[n]] <- HubCor.byMod.byGL(sumModClrs, HubsExpr, ExprList=Expr.byMod.ByGL.List[[n]])
  }
  ExprAndHubCor <- list("Expr.ByMod.ByGL"=Expr.byMod.ByGL.List,"HubCor.ByMod.byGL"=HubCor.byMod.byGL.List)
  return(ExprAndHubCor)
}



#==========================================================================================================#
#                                 MM extraction  (for all + SCZ_ncRNAs)                                    #
#==========================================================================================================#

####################Function to grab MM####################
get.MM <- function(MMs, modGenes, MEmod, color){
  all_MMs <- MMs[rownames(MMs) %in% modGenes,][MEmod]
  PandFDR.all_MMs <- cbind(data.frame(corPvalueStudent(as.matrix(all_MMs), nSamples)),
                           data.frame(p.adjust(as.matrix(corPvalueStudent(as.matrix(all_MMs), nSamples)), method="hochberg")))
  mod <- cbind(all_MMs, PandFDR.all_MMs)
  rownames(mod) <- modGenes
  colnames(mod) <- c(paste(color,"Module Membership (MM)"),paste(color,"MM Pval"), paste(color,"MM FDR"))
  return(mod)
}

####################Function to grab ME and hub statistics by module####################
get.MM.perMod <- function(sumClrs, MMList, gl){
  MMs <- NULL
  for (i in sumClrs){
    MMs[[i]] <- get.MM(MMList, gl[[i]], MEmod=paste("ME", i, sep=""), color=i)
  }
  return(MMs)
}


##################Functions to grab expression by gene list###################
extractExpr.GL <- function(color, gl){
  expr <- extract.Expr.(Subtrans_norm, moduleColors, color)
  select <- expr[colnames(expr) %in% gl]
  return(select)
}


##################Function to grab hub correlation with any gene of interest###################

HubNCCor <- function(expr, expr2){
  if (dim(expr2)[2] == 0)
    return(NULL)
  else
    HubCor <- t(data.frame(bicor(expr, expr2, use = "p", maxPOutliers = 0.05)))
  HubCorP <- data.frame(p.adjust(as.matrix(corPvalueStudent(as.matrix(HubCor), nSamples)), method="hochberg"))
  HubCorAndP <- cbind("BiMidweight Correlation" = HubCor, "FDR" = HubCorP)
  colnames(HubCorAndP)[2] <- "FDR"
  return(HubCorAndP)
}



#==========================================================================================================#
#                               Network summary for SCZ_ncRNA.SNPs                                         #
#==========================================================================================================#


#Function to figure into what modules the genes of interest fall into
ModuleAnalysis.GeneEnrichment <- function(df, gl, modClrs){ 
  matched <- names(df) %in% gl
  res <- NULL
  res <- append(res, table(matched)["TRUE"])
  res <- append(res, table(modClrs[matched]))
  return(res)
}



#==========================================================================================================#
#                                k and MM for  prioritized SCZ_ncRNA.SNPs                                  #
#==========================================================================================================#
#Grab gene names by module for a multiple list object
GeneNames.ByMod.PerGeneList <- function(all.gls, sum.modClrs, df, modClrs){
  list <-list()
  for (n in names(all.gls)){
    list[[n]] <- extract.GeneNames.PerMod(sum.modClrs, df, modClrs,all.gls[[n]])
  }
  list <- setNames(list, paste(names(all.gls)))
  return(list)
}



#==========================================================================================================#
#                                                 Cleaning                                                 #
#==========================================================================================================#
#######################CLEANUP FUNCTIONS#######################
#Cleanup for the getitall function
#Simple list cleaning; both work
cleanup <- function(obj, sum.modClrs){
  clean <- vector("list", length(obj))
  clean <- setNames(clean, names(obj))
  for (i in sum.modClrs)
    if (length(lapply(rownames(obj[[i]]), length)) > 0 && length(lapply(colnames(obj[[i]]), length)) > 0)
      clean[[i]] <- obj[[i]]
  else (clean[[i]] <- NULL)
  return(clean)
}

#Double list cleaning
list.clean <- function(obj, sum.modClrs){
  l <- vector("list", length(obj))
  l <- setNames(l, names(obj))
  names <- names(obj)
  for (n in names)
    l[[n]] <- cleanup(obj[[n]], sum.modClrs)
  return(l)
}

#Deeper list level cleaning
cleanup.l3 <- function(obj, sum.modClrs){
  l3 <- NULL
  for (n in obj)
    l3[[length(l3)+1]] <- list.clean(n, sum.modClrs)
  l3 <- setNames(l3, names(obj))
  return(l3)
}



#==========================================================================================================#
#                                  Gene set enrichment analysis (GSEA)                                     #
#==========================================================================================================#
####################GSEA functions################################
#Enrichment function
#create a setenrichment function: universe is generally between 10000 and 20000 (genetic background). To get really fance you can also set TPM > 1 in your datafile and set universe to that number
setEnrichment <- function(set1, set2, universe = 20000){
  a = sum(set1 %in% set2)
  c = length(set1) - a
  
  b = length(set2) - a
  d = universe - length(set2) - c
  
  contingency_table = matrix(c(a, c, b, d), nrow = 2)
  # one-tailed test for enrichment only
  fisher_results = fisher.test(contingency_table, alternative = "greater")
  # returns data.frame containing the lengths of the two sets, the overlap, the enrichment ratio (odds ratio) and P value
  df <- tibble::tibble( set1_length = length(set1), set2_length = length(set2), overlap = a, ratio = fisher_results$estimate, p.value = fisher_results$p.value)
  return(df)
}


GSEA.byMod <- function(background.gl, gl.of.int., universe=10000, get=FALSE){
  df <- NULL
  for (i in background.gl)
    df <- data.frame(rbind(df, setEnrichment(i, gl.of.int., universe)))
  rownames(df) <- names(background.gl)
  enriched.mods <- rownames(df[df$p.value < 0.1,])
  if (get == TRUE){
    l <- list("GSEA"=df, "Enriched Modules" = enriched.mods)
    return(l)}
  else return(enriched.mods)
}


#To test for list overlap
GSEA.byList <- function(gl.of.int.1, gl.of.int.2, universe=10000, get=FALSE){
  df <- NULL
  df <- data.frame(rbind(df, setEnrichment(gl.of.int.1, gl.of.int.2, universe)))
  enriched.list <- rownames(df[df$p.value < 0.1,])
  if (get == TRUE){
    return(df)}
  else return(enriched.mods)
}



#==========================================================================================================#
#                              gprofiler module GO-term enrichment                                         #
#==========================================================================================================#
###############GO term enrichment###################
#Horizontal bar plot of most significant GO terms
GO.barplot <- function(tab.res){
  #  sort(tab.res$p_value)
  ordered <- tab.res[order(tab.res$p_value),]
  colnames(ordered) <- c("TERM ID", "GO_TERM", "LOG_pVAL", "INTERSECTION")
  if (ordered$LOG_pVAL[1] < 0.0001){
    ordered$LOG_pVAL <- -log10(ordered$LOG_pVAL)
    ylabel <- "-LOG10 FDR"}
  else{
    ordered$LOG_pVAL <- -log10(ordered$LOG_pVAL)
    ylabel <- "-LOG10 FDR"}
  
  if (length(ordered$LOG_pVAL) > 14){
    df <- data.frame(ordered[c(2:3)])[0:15,]}
  else 
    df <- data.frame(ordered[c(2:3)])[0:length(ordered$LOG_pVAL),]
  ggplot(data=df, aes(x=reorder(GO_TERM, LOG_pVAL), y=LOG_pVAL)) +
    labs(title = paste("GO-term enrichment of", i, "module"), x = "GO TERMs", y = ylabel) +
    theme_bw()+
    #    ylim(c(0,0.1))+
    geom_bar(stat="identity")+
    coord_flip()
}



#==========================================================================================================#
#                                           Gene selection                                                 #
#==========================================================================================================#
##############GENE SELECTION###################
HubCor.byMod.byGL <- function(sumModClrs, Hubs.Expr, ExprList){
  HubCorList <- NULL
  for (i in names(ExprList)){
    HubCorList[[i]] <- HubNCCor(expr=Hubs.Expr[[i]], ExprList[[i]])
  }
  return(HubCorList)
}

###################### fully automated functions#################################
getitall <- function(all.gls, sumModClrs, df, modClrs, hubs, MMList, Hubs.Expr){
  GenesByModByGL <- GeneNames.ByMod.PerGeneList(all.gls, sumModClrs, df, modClrs)
  k.perMod.perGL <- list()
  MM.perMod.perGL <- list()
  ExprHubCor.perMod.perGL <- NULL
  for (n in names(GenesByModByGL)){
    k.perMod.perGL[[n]] <- extract.k.perMod(sumModClrs, hubs, GenesByModByGL[[n]])
    MM.perMod.perGL[[n]] <- get.MM.perMod(sumModClrs, MMList, GenesByModByGL[[n]])
  }
  ExprHubCor.perMod.perGL <- ExprAndHubCor.byMod.byGL.all(all.gls, sumModClrs, Hubs.Expr)
  res.summary <- list("k.perMod.perGL"=k.perMod.perGL, "MM.perMod.perGL"=MM.perMod.perGL, "Expr + HubCor.perMod.perGL"=ExprHubCor.perMod.perGL)
  return(res.summary)
}

#k-criterion
kCrit.name.extract <- function(sum.modClrs, list, mod.name.list){ #I removed the kOut criterion because genes with low kWithin were selected
  df.1 <- list()
  df.2 <- list()
  for (n in sum.modClrs){
    for (i in colnames(data.frame(t(data.frame((list[[n]])))))){
      if (t(data.frame(t(list[[n]]))[i])[5] < 1000 | t(data.frame(t(list[[n]]))[i])[6] > 9000){ #Selecting genes that are either among the 10% of total connectivity, network stability (kOut), or kWithin (need gene name list per module)
        df.1[[n]] <- c(df.1[[n]], i)}
      if (t(data.frame(t(list[[n]]))[i])[7] < (0.1*length(mod.name.list[[n]]))){ #Selecting genes based on kWithin connectivity and kOut criterion; the lower the rank in kOut (i.e., the higher the rank number), the more the effect of the gene is within
        df.2[[n]] <- c(df.2[[n]],i)}}
  }
  l <- list("kTotal criterion" = df.1, "kWithin criterion" = df.2)
  return(l)
}

#MM-criterion
MMCrit.name.extract <- function(sum.modClrs, list, mod.name.list){
  df.1 <- list()
  for (n in sum.modClrs){
    for (i in colnames(data.frame(t(data.frame(list[[n]]))))){
      if (t(data.frame(list[[n]])[i,])[1] > 0.6 && t(data.frame(list[[n]])[i,])[3] < 0.01){ #Selecting genes based on MM
        df.1[[n]] <- c(df.1[[n]], i)}
      if (t(data.frame(list[[n]])[i,])[1] < -0.6 && t(data.frame(list[[n]])[i,])[3] < 0.01){
        df.1[[n]] <- c(df.1[[n]],i)
      }
    }
  }
  return(df.1)
}

#Hub correlation criterion
HubCorCrit.name.extract <- function(sum.modClrs, list, mod.name.list){
  df.1 <- list()
  for (n in sum.modClrs){
    for (i in colnames(data.frame(t(data.frame(list[[n]]))))){
      if (t(data.frame(list[[n]])[i,])[1] > 0.6 && t(data.frame(list[[n]])[i,])[2] < 0.01){ #Selecting genes based on MM
        df.1[[n]] <- c(df.1[[n]], i)}
      if (t(data.frame(list[[n]])[i,])[1] < -0.6 && t(data.frame(list[[n]])[i,])[2] < 0.01){
        df.1[[n]] <- c(df.1[[n]], i)
      }
    }
  }
  return(df.1)
}



#==========================================================================================================#
#                  Analysis of non-coding - protein-coding co-expression patterns                          #
#==========================================================================================================#
##########Co-expression pattern analysis##########
GrabExpr.ByMod.ByGLs <- function(gls, sumModClrs){
  GeneExpr.ByMod <- NULL
  for (f in gls){
    GeneExpr.ByMod[[length(GeneExpr.ByMod)+1]] <- Expr.byMod.ByGL(sumModClrs, f)
  }
  GeneExpr.ByMod <- setNames(GeneExpr.ByMod, paste(names(gls)))
  return(GeneExpr.ByMod)
}

SZ_pcRNA.SNPs.GeneExpr.ByMod <- GrabExpr.ByMod.ByGLs(SZ_pcRNA.SNPs, summarized.modClrs)

BiCor.GenesOfInterest <- function(expr, expr2){
  if (dim(expr)[2] == 0 | dim(expr2)[2] == 0)
    return(NULL)
  else
    Cor <- t(data.frame(bicor(expr, expr2, use = "p", maxPOutliers = 0.05)))
  return(Cor)
}

CoReg.byMod.byGL <- function(gl, gl2){
  extract <- list()
  for (n in names(gl)){
    extract[[n]] <- BiCor.GenesOfInterest(data.frame(gl[[n]]), data.frame(gl2[[n]]))
  }
  return(extract)
}

#Function to grab co-expression (correlation) values between genes of interest
interact.CoReg <- function(gl1, gl2){
  interact <- NULL
  for (i in gl1)
    for (n in gl2)
      interact[[length(interact)+1]] <- CoReg.byMod.byGL(i, n)
    gl.names <- cross2(names(gl1), names(gl2))
    names(interact) <- lapply(gl.names, function(x) paste0(x[1], " x ", x[2:length(x)]))
    return(interact)
}
