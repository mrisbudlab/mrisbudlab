library(tidyverse)
library(knitr)
library(HGNChelper)
library(Seurat)
library(scCustomize)
library(dplyr)
library(patchwork)
library(Rcpp)
library(Matrix)
library(openxlsx)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library("AnnotationDbi")
library('org.Mm.eg.db')
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)
library(europepmc)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue = "Immune system"
gs_list = gene_sets_prepare(db_, tissue)

##Lines 29-___ outline how the scRNA-seq analysis was conducted for this study
##Lines ___-___ categorizes the scRNA-seq data according to the cell-labeling criteria used for a CyTOF experiment in this study

SMJ.dir <- "#path to folder containing barcodes.tsv, features.tsv, and matrix.mtx files for 12M and 18M SM/J"
B6.dir <- "#path to folder containing barcodes.tsv, features.tsv, and matrix.mtx files for 12M and 18M BL/6"

SMJ.12M <- Read10X_GEO(paste0(SMJ.dir,'/SMJ_12M')) #/SMJ_12M is a subfolder containing just the 12M SM/J files
SMJ.12M <- CreateSeuratObject(counts = SMJ.12M[[1]][["Gene Expression"]], project = "SMJ.12M",
                              min.cells = 3, min.features = 200)

SMJ.18M <- Read10X_GEO(paste0(SMJ.dir,'/SMJ_18M')) #/SMJ_18M is a subfolder containing just the 18M SM/J files
SMJ.18M <- CreateSeuratObject(counts = SMJ.18M[[1]][["Gene Expression"]], project = "SMJ.18M",
                              min.cells = 3, min.features = 200)

B6.12M <- Read10X_GEO(paste0(B6.dir,'/B6_12M')) #/B6_12M is a subfolder containing just the 12M BL/6 files
B6.12M <- CreateSeuratObject(counts = B6.12M[[1]][["Gene Expression"]], project = "B6.12M",
                             min.cells = 3, min.features = 200)

B6.18M <- Read10X_GEO(paste0(B6.dir,'/B6_18M')) #/B6_18M is a subfolder containing just the 18M BL/6 files
B6.18M <- CreateSeuratObject(counts = B6.18M[[1]][["Gene Expression"]], project = "B6.18M",
                             min.cells = 3, min.features = 200)

scObj <- merge(x = SMJ.12M, y = c(SMJ.18M,B6.12M,B6.18M))
scObj[["percent.mt"]] <- PercentageFeatureSet(scObj, pattern = "mt-")

VlnPlot(scObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

MTplot <- FeatureScatter(scObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
Featplot <- FeatureScatter(scObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
MTplot + Featplot

scObj_refined <- subset(scObj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

scObj_norm <- NormalizeData(scObj_refined, normalization.method = "LogNormalize", scale.factor = 10000)
scObj_norm <- FindVariableFeatures(scObj_norm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scObj_norm), 10)
plot1 <- VariableFeaturePlot(scObj_norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

scObj_norm <- ScaleData(scObj_norm,features = rownames(scObj_norm))
scObj_norm <- RunPCA(scObj_norm, features = VariableFeatures(object = scObj_norm))
DimPlot(scObj_norm, reduction = "pca")
ElbowPlot(scObj_norm)

scObj_norm <- FindNeighbors(scObj_norm, dims = 1:15)
scObj_norm <- FindClusters(scObj_norm)
scObj_norm <- RunUMAP(scObj_norm, dims = 1:15)
DimPlot(scObj_norm, reduction = "umap", label = TRUE)
OGIDscObj_normPlot <- DimPlot(scObj_norm, reduction = "umap", group.by = "orig.ident", label = FALSE) + ggtitle('B6 & SMJ PBMCs at 12 and 18 Months')

es.max = sctype_score(scRNAseqData = scObj_norm[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

## merge by cluster

cL_results = do.call("rbind", lapply(unique(scObj_norm@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(scObj_norm@meta.data[scObj_norm@meta.data$seurat_clusters==cl, ])]), decreasing = TRUE)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scObj_norm@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

scObj_norm@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  scObj_norm@meta.data$customclassif[scObj_norm@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(scObj_norm, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  
ANNOscObj_normPlot <- DimPlot(scObj_norm, reduction = "umap", group.by = 'customclassif', label = TRUE) + ggtitle('B6 & SMJ PBMCs at 12 and 18 Months')

#Look at just the 12M SMJ data
SMJ12M_subset <- subset(scObj_norm,  subset = orig.ident =="SMJ.12M")
SMJ12M_plot <- DimPlot(SMJ12M_subset,reduction = "umap",label=T)
SMJ12M_plot <- DimPlot(SMJ12M_subset,reduction = "umap",group.by ='customclassif',label=T) + ggtitle('SMJ 12M') + NoLegend()

#Look at just the 18M SMJ data
SMJ18M_subset <- subset(scObj_norm,  subset = orig.ident =="SMJ.18M")
SMJ18M_plot <- DimPlot(SMJ18M_subset,reduction = "umap",label=T)
SMJ18M_plot <- DimPlot(SMJ18M_subset,reduction = "umap",group.by ='customclassif',label=T) + ggtitle('SMJ 18M') + NoLegend()

#Look at just the 12M B6 data
B612M_subset <- subset(scObj_norm,  subset = orig.ident =="B6.12M")
B612M_plot <- DimPlot(B612M_subset,reduction = "umap",label=T)
B612M_plot <- DimPlot(B612M_subset,reduction = "umap",group.by ='customclassif',label=T) + ggtitle('B6 12M')  + NoLegend()

#Look at just the 18M B6 data
B618M_subset <- subset(scObj_norm,  subset = orig.ident =="B6.18M")
B618M_plot <- DimPlot(B618M_subset,reduction = "umap",label=T)
B618M_plot <- DimPlot(B618M_subset,reduction = "umap",group.by ='customclassif',label=T) + ggtitle('B6 18M') + NoLegend()

# View each subset of data, according to strain and timepoint
OGIDscObj_normPlot + ANNOscObj_normPlot
SMJ12M_plot + SMJ18M_plot + B612M_plot + B618M_plot

#query how many of each cell type are in the each group & the overall data 
#this example shows the number of Macrophages
#query how many Macrophages are in the each group & overall data 
length(which(SMJ12M_subset@meta.data$customclassif=="Macrophages")) # 36 
length(which(SMJ18M_subset@meta.data$customclassif=="Macrophages")) # 393 
length(which(B612M_subset@meta.data$customclassif=="Macrophages")) # 622 
length(which(B618M_subset@meta.data$customclassif=="Macrophages")) # 1203 
length(which(scObj_norm@meta.data$customclassif=="Macrophages")) # total=2254

#Remove Erythrocytes & platelets (ErP)
##Erythrocytes were meant to be removed with RBC lysis
##Platelets were meant to be removed with centrifugation
##Remaining ErP do not yield information of biological significance
NoP_subset <- subset(scObj_norm,  subset = customclassif != "Platelets")
NoErP_subset <- subset(NoP_subset,  subset = customclassif != "Erythroid-like and erythroid precursor cells")
DimPlot(NoErP_subset,reduction = "umap",group.by ='orig.ident',label=T)
NoErP_subset <- RunPCA(NoErP_subset, features = VariableFeatures(object = NoErP_subset))
DimPlot(NoErP_subset, reduction = "pca")
ElbowPlot(NoErP_subset)

NoErP_subset <- FindNeighbors(NoErP_subset, dims = 1:15)
NoErP_subset <- FindClusters(NoErP_subset)
NoErP_subset <- RunUMAP(NoErP_subset, dims = 1:15)
UnSup_Plot <- DimPlot(NoErP_subset, reduction = "umap", label = TRUE,repel = FALSE)
OGIDs_plot <- DimPlot(NoErP_subset, reduction = "umap", group.by = "orig.ident", label = FALSE, cols = brewer.pal(4,"Set2"),repel = TRUE) + ggtitle('B6 & SMJ PBMCs at 12 & 18 Months')
ANNO_Plot <- DimPlot(NoErP_subset, reduction = "umap", group.by = "customclassif", label = FALSE, cols = brewer.pal(11,"Paired"),repel = TRUE) + ggtitle('B6 & SMJ PBMCs at 12 & 18 Months')
UnSup_Plot + OGIDs_plot + ANNO_Plot

#Look at just the 12M SMJ data from the NoErP_subset
SMJ12M_subset <- subset(NoErP_subset,  subset = orig.ident =="SMJ.12M")
SMJ12M_plot <- DimPlot(SMJ12M_subset,reduction = "umap",label=T)
SMJ12M_plot <- DimPlot(SMJ12M_subset,reduction = "umap",group.by ='customclassif',label = TRUE, cols = brewer.pal(11,"Paired"),repel = TRUE) + ggtitle('SMJ 12M') + NoLegend()

#Look at just the 18M SMJ data
SMJ18M_subset <- subset(NoErP_subset,  subset = orig.ident =="SMJ.18M")
SMJ18M_plot <- DimPlot(SMJ18M_subset,reduction = "umap",label=T)
SMJ18M_plot <- DimPlot(SMJ18M_subset,reduction = "umap",group.by ='customclassif',label = TRUE, cols = brewer.pal(11,"Paired"),repel = TRUE) + ggtitle('SMJ 18M') + NoLegend()

#Look at just the 12M B6 data
B612M_subset <- subset(NoErP_subset,  subset = orig.ident =="B6.12M")
B612M_plot <- DimPlot(B612M_subset,reduction = "umap",label=T)
B612M_plot <- DimPlot(B612M_subset,reduction = "umap",group.by ='customclassif',label = TRUE, cols = brewer.pal(11,"Paired"),repel = TRUE) + ggtitle('B6 12M')  + NoLegend()

#Look at just the 18M B6 data
B618M_subset <- subset(NoErP_subset,  subset = orig.ident =="B6.18M")
B618M_plot <- DimPlot(B618M_subset,reduction = "umap",label=T)
B618M_plot <- DimPlot(B618M_subset,reduction = "umap",group.by ='customclassif',label = TRUE, cols = brewer.pal(11,"Paired"),repel = TRUE) + ggtitle('B6 18M') + NoLegend()

OGIDs_plot + ANNO_Plot
B612M_plot + B618M_plot + SMJ12M_plot + SMJ18M_plot 

#Look at two subsets of the data
#This example looks at the 12M B6 & SMJ data
B6SMJ12M_subset <- subset(NoErP_subset, subset = orig.ident != "SMJ.18M")
B6SMJ12M_subset <- subset(B6SMJ12M_subset, subset = orig.ident != "B6.18M")
B6SMJ12M_plot1 <- DimPlot(B6SMJ12M_subset,reduction = "umap",label=T, group.by = "orig.ident") + ggtitle('B6 and SMJ 12M') + NoLegend()
B6SMJ12M_plot2 <- DimPlot(B6SMJ12M_subset,reduction = "umap",group.by ='customclassif',label=T) + ggtitle('B6 and SMJ 12M') + NoLegend()
B6SMJ12M_plot1 + B6SMJ12M_plot2


#Analyze one particular cell type
#this example analyzes just the Macrophages
Macrophage_subset <- subset(NoErP_subset,  subset = customclassif =="Macrophages")
Macrophage_Plot <- DimPlot(Macrophage_subset,reduction = "umap",group.by ='orig.ident',label=T)
Macrophage_subset <- RunPCA(Macrophage_subset, features = VariableFeatures(object = Macrophage_subset))
DimPlot(Macrophage_subset, reduction = "pca")
ElbowPlot(Macrophage_subset)

Macrophage_subset <- FindNeighbors(Macrophage_subset, dims = 1:15)
Macrophage_subset <- FindClusters(Macrophage_subset)
Macrophage_subset <- RunUMAP(Macrophage_subset, dims = 1:15)
Macrophage_PlotCl <- DimPlot(Macrophage_subset, reduction = "umap", label = TRUE) + ggtitle('Macrophages')
Macrophage_PlotID <- DimPlot(Macrophage_subset, reduction = "umap", group.by = "orig.ident", label = FALSE, cols = brewer.pal(4,"Set2"),repel = TRUE) + ggtitle('Macrophages')

Macrophage_PlotCl + Macrophage_PlotID


#Determine DEGs w/in relevant subclusters within a cell type
#this example looks at different comparisons among Macrophages
SMJvB6_Macrophages <- FindMarkers(Macrophage_subset, ident.1 = c(2,7),
                                  ident.2 = c(0,1,3,4,5,6,8), min.pct = 0.25)
Macrophage_names <- list('SMJvB6_Macrophages' = SMJvB6_Macrophages) #'c27vc586' = c27vc586.Macrophage.markers, 'c27vc0134' = c27vc0134.Macrophage.markers,'c27vc0134586' = c27vc0134586.Macrophage.markers)
#save these lists of DEGs in .xlsx format
write.xlsx(Macrophage_names, rowNames = TRUE, file="#path of your choosing", #append = TRUE if you add multiple analyses to the same file)

##Within the same excel sheet, before moving to the next step, save the lists of DEGs you'll want to further analyze
###for example, if you wanted to separate upregulated and downregulated DEGs, you would do that before moving on

###the following examples show GSEA analysis for upregulated (Macrophage1+) and downregulated (Macrophage1-) DEGs identified in Macrophages
MacrophageUp_df = read.xlsx(xlsxFile = "#file folder path/Macrophage.xlsx", sheet = "Macrophage1+")
MacrophageUp_gene_list <- MacrophageUp_df$avg_log2FC
names(MacrophageUp_gene_list) <- mapIds(org.Mm.eg.db, keys = MacrophageUp_df$gene, 'ENTREZID', 'SYMBOL')
MacrophageUp_gene_list <- na.omit(MacrophageUp_gene_list)
MacrophageUp_gene_list <- sort(MacrophageUp_gene_list,decreasing = TRUE)

MacrophageUp_gse <- gseGO(geneList=MacrophageUp_gene_list, 
                          ont ="BP", 
                          keyType = "ENTREZID", 
                          nPerm = 10000, 
                          minGSSize = 3, 
                          maxGSSize = 800, 
                          pvalueCutoff = 0.05, 
                          verbose = TRUE, 
                          OrgDb = organism, 
                          pAdjustMethod = "none")

as.data.frame(MacrophageUp_gse)

require(DOSE)
MacrophageUp_dotplot <- dotplot(MacrophageUp_gse, showCategory=15, split=".sign", font.size=7) + facet_grid(.~.sign) + ggtitle("Macrophage Up GSEA SMJ vs. B6")

MacrophageUp_ridgeplot <- ridgeplot(MacrophageUp_gse) + labs(x = "enrichment distribution")

MacrophageUp_cnetplot <- cnetplot(MacrophageUp_gse, categorySize="pvalue", foldChange=MacrophageUp_gene_list, showCategory = 3)

MacrophageUp_gseaplot <- gseaplot(MacrophageUp_gse, by = "all", title = MacrophageUp_gse$Description[1], geneSetID = 1)

MacrophageDown_df = read.xlsx(xlsxFile = "#file folder path/Macrophage.xlsx", sheet = "Macrophage1-")
MacrophageDown_gene_list <- MacrophageDown_df$avg_log2FC
names(MacrophageDown_gene_list) <- mapIds(org.Mm.eg.db, keys = MacrophageDown_df$gene, 'ENTREZID', 'SYMBOL')
MacrophageDown_gene_list <- na.omit(MacrophageDown_gene_list)
MacrophageDown_gene_list <- sort(MacrophageDown_gene_list,decreasing = TRUE)

MacrophageDown_gse <- gseGO(geneList=MacrophageDown_gene_list, 
                            ont ="BP", 
                            keyType = "ENTREZID", 
                            nPerm = 10000, 
                            minGSSize = 3, 
                            maxGSSize = 800, 
                            pvalueCutoff = 0.05, 
                            verbose = TRUE, 
                            OrgDb = organism, 
                            pAdjustMethod = "none")

as.data.frame(MacrophageDown_gse)

require(DOSE)
MacrophageDown_dotplot <- dotplot(MacrophageDown_gse, showCategory=15, split=".sign", font.size=7) + facet_grid(.~.sign) + ggtitle("Macrophage Down GSEA SMJ vs. B6")

MacrophageDown_ridgeplot <- ridgeplot(MacrophageDown_gse) + labs(x = "enrichment distribution")

MacrophageDown_cnetplot <- cnetplot(MacrophageDown_gse, categorySize="pvalue", foldChange=MacrophageDown_gene_list, showCategory = 3)

MacrophageDown_gseaplot <- gseaplot(MacrophageDown_gse, by = "all", title = MacrophageDown_gse$Description[1], geneSetID = 1)

MacrophageUp_dotplot + MacrophageDown_dotplot


################################################################

#scRNA-seq data organized according to CyTOF labeling criteria

DoHeatmap(object = NoErP_subset, features = c("Cd3d","Cd19","Klrb1","Cd4","Cd8a","Cd44","Sell","Itgam","Itgax","Ptprc"),group.by = 'customclassif', size = 3, angle = 90) + scale_fill_gradientn(colors = c("purple4", "pink", "turquoise"))  + guides(color="none") + 
  theme(text = element_text(size = 15)) #+ ggtitle("CyToF Assignments")