library(remotes)
remotes::install_version("Seurat", "4.4.0")

devtools::install_version("dbplyr", version = "2.3.4")

remotes::install_version("matrixStats", version="1.1.0")

library(tidyverse)
library(knitr)
library(HGNChelper)
library(Seurat)
library(umap)
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
library(SingleR)
library(SingleCellExperiment)
library(matrixStats)
library(scRNAseq)
library(scuttle)
library(TabulaMurisData)
library(ExperimentHub)


SMJ.dir <- "Path to directory with 12M and 20M SM/J data"
B6.dir <- "Path to directory with 12M and 20M SM/J data"

SMJ.12M <- Read10X_GEO(paste0(SMJ.dir,'/SMJ_12M'))
SMJ.12M <- CreateSeuratObject(counts = SMJ.12M[[1]][["Gene Expression"]], project = "SMJ.12M",
                              min.cells = 3, min.features = 200)

SMJ.20M <- Read10X_GEO(paste0(SMJ.dir,'/SMJ_20M'))
SMJ.20M <- CreateSeuratObject(counts = SMJ.20M[[1]][["Gene Expression"]], project = "SMJ.20M",
                              min.cells = 3, min.features = 200)

B6.12M <- Read10X_GEO(paste0(B6.dir,'/B6_12M'))
B6.12M <- CreateSeuratObject(counts = B6.12M[[1]][["Gene Expression"]], project = "B6.12M",
                             min.cells = 3, min.features = 200)

B6.20M <- Read10X_GEO(paste0(B6.dir,'/B6_20M'))
B6.20M <- CreateSeuratObject(counts = B6.20M[[1]][["Gene Expression"]], project = "B6.20M",
                             min.cells = 3, min.features = 200)

scObj <- merge(x = SMJ.12M, y = c(SMJ.20M,B6.12M,B6.20M))
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
CluPlot <- DimPlot(scObj_norm, reduction = "umap", label = TRUE) + NoLegend() 
OGIDscObj_normPlot <- DimPlot(scObj_norm, reduction = "umap", group.by = "orig.ident", label = FALSE) + NoLegend() #+ ggtitle('B6 & SM/J PBMCs at 12 and 20 Months') 

OGIDscObj_normPlot

ref <- celldex::MouseRNAseqData()
class(ref)
table(ref$label.main)

#can use "fine" or "main"

scObj_normSingR <- SingleR(test = as.SingleCellExperiment(scObj_norm), ref = ref, labels = ref$label.main)

scObj_norm$singlr_labels <- scObj_normSingR$labels

scObjPlot_SingR <- DimPlot(scObj_norm, reduction = 'umap', group.by = 'singlr_labels', label = FALSE) + NoLegend() #+ ggtitle('SingleR "MouseRNAseqData" Annotation') 

scObjPlot_SingR

scObjPlot_SingR + OGIDscObj_normPlot+ CluPlot

#Check which cell populations are Ifng+
Ifng_FeatPlot <- FeaturePlot(scObj_norm,c("Ifng"))
Ifng_FeatPlot+ scObjPlot_SingR + OGIDscObj_normPlot

#Look at just the 12M SMJ data
SMJ12M_subset <- subset(scObj_norm,  subset = orig.ident =="SMJ.12M")
SMJ12M_plot <- DimPlot(SMJ12M_subset,reduction = "umap",label=T)
SMJ12M_plot <- DimPlot(SMJ12M_subset,reduction = "umap",group.by ='singlr_labels',label=T) + ggtitle('SMJ 12M') + NoLegend()

#Look at just the 20M SMJ data
SMJ20M_subset <- subset(scObj_norm,  subset = orig.ident =="SMJ.20M")
SMJ20M_plot <- DimPlot(SMJ20M_subset,reduction = "umap",label=T)
SMJ20M_plot <- DimPlot(SMJ20M_subset,reduction = "umap",group.by ='singlr_labels',label=T) + ggtitle('SMJ 20M') + NoLegend()

#Look at just the 12M B6 data
B612M_subset <- subset(scObj_norm,  subset = orig.ident =="B6.12M")
B612M_plot <- DimPlot(B612M_subset,reduction = "umap",label=T)
B612M_plot <- DimPlot(B612M_subset,reduction = "umap",group.by ='singlr_labels',label=T) + ggtitle('B6 12M')  + NoLegend()

#Look at just the 20M B6 data
B620M_subset <- subset(scObj_norm,  subset = orig.ident =="B6.20M")
B620M_plot <- DimPlot(B620M_subset,reduction = "umap",label=T)
B620M_plot <- DimPlot(B620M_subset,reduction = "umap",group.by ='singlr_labels',label=T) + ggtitle('B6 20M') + NoLegend()

B612M_plot + B620M_plot + SMJ12M_plot + SMJ20M_plot

#counting cell populations based on the cell types identified by SingR
#query how many Erythrocytes are in the overall data 
length(which(SMJ12M_subset@meta.data$singlr_labels=="Erythrocytes")) # 14 
length(which(SMJ20M_subset@meta.data$singlr_labels=="Erythrocytes")) # 5 
length(which(B612M_subset@meta.data$singlr_labels=="Erythrocytes")) # 3 
length(which(B620M_subset@meta.data$singlr_labels=="Erythrocytes")) # 2
length(which(scObj_norm@meta.data$singlr_labels=="Erythrocytes")) # total=24

#query how many Endothelial cells are in each group & the overall data 
length(which(SMJ12M_subset@meta.data$singlr_labels=="Endothelial cells")) # 1 
length(which(SMJ20M_subset@meta.data$singlr_labels=="Endothelial cells")) # 1 
length(which(B612M_subset@meta.data$singlr_labels=="Endothelial cells")) # 4 
length(which(B620M_subset@meta.data$singlr_labels=="Endothelial cells")) # 5 
length(which(scObj_norm@meta.data$singlr_labels=="Endothelial cells")) # total=11

#query how many Fibroblasts cells are in each group & the overall data 
length(which(SMJ12M_subset@meta.data$singlr_labels=="Fibroblasts")) # 0 
length(which(SMJ20M_subset@meta.data$singlr_labels=="Fibroblasts")) # 2 
length(which(B612M_subset@meta.data$singlr_labels=="Fibroblasts")) # 0 
length(which(B620M_subset@meta.data$singlr_labels=="Fibroblasts")) # 0 
length(which(scObj_norm@meta.data$singlr_labels=="Fibroblasts")) # total=2

#query how many B cells are in each group & the overall data 
length(which(SMJ12M_subset@meta.data$singlr_labels=="B cells")) # 1361 
length(which(SMJ20M_subset@meta.data$singlr_labels=="B cells")) # 1583 
length(which(B612M_subset@meta.data$singlr_labels=="B cells")) # 2035 
length(which(B620M_subset@meta.data$singlr_labels=="B cells")) # 1021 
length(which(NoEr_subset@meta.data$singlr_labels=="B cells")) # total=6000

#query how many Dendritic cells are in each group & the overall data 
length(which(SMJ12M_subset@meta.data$singlr_labels=="Dendritic cells")) # 3 
length(which(SMJ20M_subset@meta.data$singlr_labels=="Dendritic cells")) # 2 
length(which(B612M_subset@meta.data$singlr_labels=="Dendritic cells")) # 3 
length(which(B620M_subset@meta.data$singlr_labels=="Dendritic cells")) # 7 
length(which(NoEr_subset@meta.data$singlr_labels=="Dendritic cells")) # total=15

#query how many Granulocytes cells are in each group & the overall data 
length(which(SMJ12M_subset@meta.data$singlr_labels=="Granulocytes")) # 27 
length(which(SMJ20M_subset@meta.data$singlr_labels=="Granulocytes")) # 24 
length(which(B612M_subset@meta.data$singlr_labels=="Granulocytes")) # 102 
length(which(B620M_subset@meta.data$singlr_labels=="Granulocytes")) # 134 
length(which(NoEr_subset@meta.data$singlr_labels=="Granulocytes")) # total=287

#query how many Macrophages cells are in each group & the overall data 
length(which(SMJ12M_subset@meta.data$singlr_labels=="Macrophages")) # 65 
length(which(SMJ20M_subset@meta.data$singlr_labels=="Macrophages")) # 387 
length(which(B612M_subset@meta.data$singlr_labels=="Macrophages")) # 94 
length(which(B620M_subset@meta.data$singlr_labels=="Macrophages")) # 108 
length(which(NoEr_subset@meta.data$singlr_labels=="Macrophages")) # total=654

#query how many Monocytes cells are in each group & the overall data 
length(which(SMJ12M_subset@meta.data$singlr_labels=="Monocytes")) # 510 
length(which(SMJ20M_subset@meta.data$singlr_labels=="Monocytes")) # 287 
length(which(B612M_subset@meta.data$singlr_labels=="Monocytes")) # 628 
length(which(B620M_subset@meta.data$singlr_labels=="Monocytes")) # 1208 
length(which(NoEr_subset@meta.data$singlr_labels=="Monocytes")) # total=2633

#query how many NK cells are in each group & the overall data 
length(which(SMJ12M_subset@meta.data$singlr_labels=="NK cells")) # 68 
length(which(SMJ20M_subset@meta.data$singlr_labels=="NK cells")) # 57 
length(which(B612M_subset@meta.data$singlr_labels=="NK cells")) # 151 
length(which(B620M_subset@meta.data$singlr_labels=="NK cells")) # 96 
length(which(NoEr_subset@meta.data$singlr_labels=="NK cells")) # total=372

#query how many T cells are in each group & the overall data 
length(which(SMJ12M_subset@meta.data$singlr_labels=="T cells")) # 361 
length(which(SMJ20M_subset@meta.data$singlr_labels=="T cells")) # 233 
length(which(B612M_subset@meta.data$singlr_labels=="T cells")) # 408 
length(which(B620M_subset@meta.data$singlr_labels=="T cells")) # 528 
length(which(NoEr_subset@meta.data$singlr_labels=="T cells")) # total=1530

#Analyze just the T cells
Tcell_subset <- subset(scObj_norm,  subset = singlr_labels =="T cells")
Tcell_Plot <- DimPlot(Tcell_subset,reduction = "umap",group.by ='orig.ident',label=T)
Tcell_subset <- RunPCA(Tcell_subset, features = VariableFeatures(object = Tcell_subset))
DimPlot(Tcell_subset, reduction = "pca")
ElbowPlot(Tcell_subset)

Tcell_subset <- FindNeighbors(Tcell_subset, dims = 1:15)
Tcell_subset <- FindClusters(Tcell_subset)
Tcell_subset <- RunUMAP(Tcell_subset, dims = 1:15)
Tcell_PlotCl <- DimPlot(Tcell_subset, reduction = "umap", label = TRUE) + ggtitle('T Cells')
Tcell_PlotID <- DimPlot(Tcell_subset, reduction = "umap", group.by = "orig.ident", label = FALSE, cols = brewer.pal(4,"Set2"),repel = TRUE) + NoLegend() #+ ggtitle('T Cells')

Tcell_PlotCl + Tcell_PlotID

#Ifng+ T Cells
IfngTcell_subset <- subset(Tcell_subset,  subset = Ifng >= 1)
length(which(IfngTcell_subset@meta.data$singlr_labels=="T cells")) # total=133
length(which(IfngTcell_subset@meta.data$orig.ident =="B6.12M")) # total=56
length(which(IfngTcell_subset@meta.data$orig.ident =="B6.20M")) # total=44
length(which(IfngTcell_subset@meta.data$orig.ident =="SMJ.12M")) # total=17
length(which(IfngTcell_subset@meta.data$orig.ident =="SMJ.20M")) # total=16

Tcell_IfngFeatPlot <- FeaturePlot(Tcell_subset,c("Ifng"))#"Cd4","Cd8a",
Tcell_IfngFeatPlot + Tcell_PlotCl + Tcell_PlotID

#Cd4+ T Cells
Tcell_4FeatPlot <- FeaturePlot(Tcell_subset,c("Cd4")) + ggtitle("CD4+")
Cd4Tcell_subset <- subset(Tcell_subset,  subset = Cd4 >= 1)
length(which(Cd4Tcell_subset@meta.data$singlr_labels=="T cells")) # total=427
length(which(Cd4Tcell_subset@meta.data$orig.ident =="B6.12M")) # total=84
length(which(Cd4Tcell_subset@meta.data$orig.ident =="B6.20M")) # total=176
length(which(Cd4Tcell_subset@meta.data$orig.ident =="SMJ.12M")) # total=93
length(which(Cd4Tcell_subset@meta.data$orig.ident =="SMJ.20M")) # total=74

IfngCd4Tcell_subset <- subset(Cd4Tcell_subset,  subset = Ifng >= 1)
length(which(IfngCd4Tcell_subset@meta.data$singlr_labels=="T cells")) # total=35
length(which(IfngCd4Tcell_subset@meta.data$orig.ident =="B6.12M")) # total=10
length(which(IfngCd4Tcell_subset@meta.data$orig.ident =="B6.20M")) # total=21
length(which(IfngCd4Tcell_subset@meta.data$orig.ident =="SMJ.12M")) # total=2
length(which(IfngCd4Tcell_subset@meta.data$orig.ident =="SMJ.20M")) # total=2

Cd4T_subset <- subset(Tcell_subset,  subset = seurat_clusters == c("0","1","2","6"))
Cd4T_plot1 <- DimPlot(Cd4T_subset,reduction = "umap",label=T)
Cd4T_plot2 <- DimPlot(Cd4T_subset,reduction = "umap",group.by ='orig.ident',label=F) + ggtitle('CD4+ T Cells')

SMJvB6_Cd4T <- FindMarkers(Cd4T_subset, ident.1 = c(0),ident.2 = c(2), min.pct = 0.25) #0 is the B6 cluster; 2 is the SM/J cluster
Cd4T_names <- list('SMJvB6_Cd4T' = SMJvB6_Cd4T)
write.xlsx(Cd4T_names, rowNames = TRUE, file="Path to directory where excel output will be housed/Cd4T.xlsx", append = TRUE)

#Cd8a+ T Cells
Cd8aTcell_subset <- subset(Tcell_subset,  subset = Cd8a >= 1)
length(which(Cd8aTcell_subset@meta.data$singlr_labels=="T cells")) # total=232
length(which(Cd8aTcell_subset@meta.data$orig.ident =="B6.12M")) # total=108
length(which(Cd8aTcell_subset@meta.data$orig.ident =="B6.20M")) # total=44
length(which(Cd8aTcell_subset@meta.data$orig.ident =="SMJ.12M")) # total=56
length(which(Cd8aTcell_subset@meta.data$orig.ident =="SMJ.20M")) # total=24

IfngCd8aTcell_subset <- subset(Cd8aTcell_subset,  subset = Ifng >= 1)
length(which(IfngCd8aTcell_subset@meta.data$singlr_labels=="T cells")) # total=35
length(which(IfngCd8aTcell_subset@meta.data$orig.ident =="B6.12M")) # total=26
length(which(IfngCd8aTcell_subset@meta.data$orig.ident =="B6.20M")) # total=2
length(which(IfngCd8aTcell_subset@meta.data$orig.ident =="SMJ.12M")) # total=3
length(which(IfngCd8aTcell_subset@meta.data$orig.ident =="SMJ.20M")) # total=4

Tcell_8FeatPlot <- FeaturePlot(Tcell_subset,c("Cd8a")) + ggtitle("CD8a+")

Cd8aT_subset <- subset(Tcell_subset,  subset = seurat_clusters == c("3","4","5","8"))
Cd8aT_plot1 <- DimPlot(Cd8aT_subset,reduction = "umap",label=T)
Cd8aT_plot2 <- DimPlot(Cd8aT_subset,reduction = "umap",group.by ='orig.ident',label=F) + ggtitle('CD8a+ T Cells')

SMJvB6_Cd8aT <- FindMarkers(Cd8aT_subset, ident.1 = c(3),ident.2 = c(5), min.pct = 0.25) #3 is the B6 cluster; 5 is the SM/J cluster
Cd8aT_names <- list('SMJvB6_Cd8aT' = SMJvB6_Cd8aT)
write.xlsx(Cd8aT_names, rowNames = TRUE, file="Path to directory where excel output will be housed/Cd8aT.xlsx", append = TRUE)
FeaturePlot(Cd8aT_subset,c("Saa3","Gm10260","Ccl5")) + Cd8aT_plot1 + Cd8aT_plot2

#Analyze just the B cells
Bcell_subset <- subset(scObj_norm,  subset = singlr_labels =="B cells")
BceLL_Plot <- DimPlot(Bcell_subset,reduction = "umap",group.by ='orig.ident',label=T)
Bcell_subset <- RunPCA(Bcell_subset, features = VariableFeatures(object = Bcell_subset))
DimPlot(Bcell_subset, reduction = "pca")
ElbowPlot(Bcell_subset)

Bcell_subset <- FindNeighbors(Bcell_subset, dims = 1:15)
Bcell_subset <- FindClusters(Bcell_subset)
Bcell_subset <- RunUMAP(Bcell_subset, dims = 1:15)
Bcell_PlotCl <- DimPlot(Bcell_subset, reduction = "umap", label = TRUE) + ggtitle('B Cells') + NoLegend()
Bcell_PlotID <- DimPlot(Bcell_subset, reduction = "umap", group.by = "orig.ident", label = FALSE, cols = brewer.pal(4,"Set2"),repel = TRUE) + ggtitle('B Cells') + NoLegend()

Bcell_PlotCl + Bcell_PlotID

SMJvB6_BCell <- FindMarkers(Bcell_subset, ident.1 = c(0,1,4,6,10),ident.2 = c(2,3,5,7,8,11,12,14), min.pct = 0.25) #(0,1,4,6,10) is the B6 cluster; (2,3,5,7,8,12,14) is the SM/J cluster
ySMJvB6_BCell <- FindMarkers(Bcell_subset, ident.1 = c(0,1,4,6,10),ident.2 = c(2,7), min.pct = 0.25) #(0,1,4,6,10) is the B6 cluster; (2,7) is the 12M SM/J cluster
oSMJvB6_BCell <- FindMarkers(Bcell_subset, ident.1 = c(0,1,4,6,10),ident.2 = c(5,8,12,14), min.pct = 0.25) #(0,1,4,6,10) is the B6 cluster; (5,8,12,14) is the 20M SM/J cluster
agingSMJ_BCell <- FindMarkers(Bcell_subset, ident.1 = c(2,7),ident.2 = c(5,8,12,14), min.pct = 0.25) #(2,7) is the B6 cluster; (5,8,12,14) is the 20M SM/J cluster
BCell_names <- list('SMJvB6_BCell' = SMJvB6_BCell,'ySMJvB6_BCell' = ySMJvB6_BCell,'oSMJvB6_BCell' = oSMJvB6_BCell,'agingSMJ_BCell' = agingSMJ_BCell)
write.xlsx(BCell_names, rowNames = TRUE, file="Path to directory where excel output will be housed/BCell.xlsx", append = TRUE)
FeaturePlot(Cd8aT_subset,c("Saa3","Gm10260","Ccl5")) + Cd8aT_plot1 + Cd8aT_plot2

SMJPain_BCell <- FindMarkers(Bcell_subset, ident.1 = c(0,1,4,6,10,2,7),ident.2 = c(5,8,12,14), min.pct = 0.25) #(0,1,4,6,10+2,7) is the B6+12MSMJ cluster; (5,8,12,14) is the 20M SM/J cluster
PainBCell_names <- list('SMJPain_BCell' = SMJPain_BCell)
write.xlsx(PainBCell_names, rowNames = TRUE, file="Path to directory where excel output will be housed/BCellPain.xlsx", append = TRUE)


#Analyze just the Granulocytes
Granulocytes_subset <- subset(scObj_norm,  subset = singlr_labels =="Granulocytes")
Granulocytes_Plot <- DimPlot(Granulocytes_subset,reduction = "umap",group.by ='orig.ident',label=T)
Granulocytes_subset <- RunPCA(Granulocytes_subset, features = VariableFeatures(object = Granulocytes_subset))
DimPlot(Granulocytes_subset, reduction = "pca")
ElbowPlot(Granulocytes_subset)

Granulocytes_subset <- FindNeighbors(Granulocytes_subset, dims = 1:15)
Granulocytes_subset <- FindClusters(Granulocytes_subset)
Granulocytes_subset <- RunUMAP(Granulocytes_subset, dims = 1:15)
Granulocytes_PlotCl <- DimPlot(Granulocytes_subset, reduction = "umap", label = TRUE) + ggtitle('Granulocytes') + NoLegend()
Granulocytes_PlotID <- DimPlot(Granulocytes_subset, reduction = "umap", group.by = "orig.ident", label = FALSE, cols = brewer.pal(4,"Set2"),repel = TRUE) + ggtitle('Granulocytes')+ NoLegend()

Granulocytes_PlotCl + Granulocytes_PlotID

Granulocyte_FeatPlot <- FeaturePlot(Granulocytes_subset,c("Ly6g","Fcer1a","Siglecf"))
Granulocyte_FeatPlot + Granulocytes_PlotCl + Granulocytes_PlotID

#Analyze just the Macrophages
Macrophages_subset <- subset(scObj_norm,  subset = singlr_labels =="Macrophages")
Macrophages_Plot <- DimPlot(Macrophages_subset,reduction = "umap",group.by ='orig.ident',label=T)
Macrophages_subset <- RunPCA(Macrophages_subset, features = VariableFeatures(object = Macrophages_subset))
DimPlot(Macrophages_subset, reduction = "pca")
ElbowPlot(Macrophages_subset)

Macrophages_subset <- FindNeighbors(Macrophages_subset, dims = 1:15)
Macrophages_subset <- FindClusters(Macrophages_subset)
Macrophages_subset <- RunUMAP(Macrophages_subset, dims = 1:15)
Macrophages_PlotCl <- DimPlot(Macrophages_subset, reduction = "umap", label = TRUE) + ggtitle('Macrophages') + NoLegend()
Macrophages_PlotID <- DimPlot(Macrophages_subset, reduction = "umap", group.by = "orig.ident", label = FALSE, cols = brewer.pal(4,"Set2"),repel = TRUE) + ggtitle('Macrophages') + NoLegend()

Macrophages_PlotCl + Macrophages_PlotID


SMJvB6_Macrophage <- FindMarkers(Macrophages_subset, ident.1 = c(1),ident.2 = c(0,2,3,5), min.pct = 0.25) #(1) is the B6 cluster; (0,2,3,5) is the SM/J cluster
ySMJvB6_Macrophage <- FindMarkers(Macrophages_subset, ident.1 = c(1),ident.2 = c(3), min.pct = 0.25) #(1) is the B6 cluster; (3) is the SM/J cluster
oSMJvB6_Macrophage <- FindMarkers(Macrophages_subset, ident.1 = c(1),ident.2 = c(0,2,5), min.pct = 0.25) #(1) is the B6 cluster; (0,2,5) is the SM/J cluster
agingSMJ_Macrophage <- FindMarkers(Macrophages_subset, ident.1 = c(3),ident.2 = c(0,2,5), min.pct = 0.25) #(3) is the 12M SM/J cluster; (0,2,5) is the 20M SM/J cluster
Macrophage_names <- list('SMJvB6_Macrophage' = SMJvB6_Macrophage,'ySMJvB6_Macrophage' = ySMJvB6_Macrophage,'oSMJvB6_Macrophage' = oSMJvB6_Macrophage,'agingSMJ_Macrophage' = agingSMJ_Macrophage)
write.xlsx(Macrophage_names, rowNames = TRUE, file="Path to directory where excel output will be housed/Macrophage.xlsx", append = TRUE)

SMJPain_Macrophage <- FindMarkers(Macrophages_subset, ident.1 = c(1,3),ident.2 = c(0,2,5), min.pct = 0.25) #(1+3) is the B6+SMJ cluster; (0,2,5) is the SM/J cluster
PainMacrophage_names <- list('SMJPain_Macrophage' = SMJPain_Macrophage)
write.xlsx(PainMacrophage_names, rowNames = TRUE, file="Path to directory where excel output will be housed/MacrophagePain.xlsx", append = TRUE)

#Analyze just the Monocytes
Monocytes_subset <- subset(scObj_norm,  subset = singlr_labels =="Monocytes")
Monocytes_Plot <- DimPlot(Monocytes_subset,reduction = "umap",group.by ='orig.ident',label=T)
Monocytes_subset <- RunPCA(Monocytes_subset, features = VariableFeatures(object = Monocytes_subset))
DimPlot(Monocytes_subset, reduction = "pca")
ElbowPlot(Monocytes_subset)

Monocytes_subset <- FindNeighbors(Monocytes_subset, dims = 1:15)
Monocytes_subset <- FindClusters(Monocytes_subset)
Monocytes_subset <- RunUMAP(Monocytes_subset, dims = 1:15)
Monocytes_PlotCl <- DimPlot(Monocytes_subset, reduction = "umap", label = TRUE) + ggtitle('Monocytes') + NoLegend()
Monocytes_PlotID <- DimPlot(Monocytes_subset, reduction = "umap", group.by = "orig.ident", label = FALSE, cols = brewer.pal(4,"Set2"),repel = TRUE) + ggtitle('Monocytes') + NoLegend()

Monocytes_PlotCl + Monocytes_PlotID

SMJvB6_Monocyte <- FindMarkers(Monocytes_subset, ident.1 = c(1, 2, 3, 4, 5, 7, 8, 9),ident.2 = c(0, 6, 10), min.pct = 0.25) #(1, 2, 3, 4, 5, 7, 8, 9) is the B6 cluster; (0, 6, 10) is the SM/J cluster
Monocyte_names <- list('SMJvB6_Monocyte' = SMJvB6_Monocyte)
write.xlsx(Monocyte_names, rowNames = TRUE, file="Path to directory where excel output will be housed/Monocyte.xlsx", append = TRUE)


#Analyze just the NK Cells
NK_subset <- subset(scObj_norm,  subset = singlr_labels =="NK cells")
NK_Plot <- DimPlot(NK_subset,reduction = "umap",group.by ='orig.ident',label=T)
NK_subset <- RunPCA(NK_subset, features = VariableFeatures(object = NK_subset))
DimPlot(NK_subset, reduction = "pca")
ElbowPlot(NK_subset)

NK_subset <- FindNeighbors(NK_subset, dims = 1:15)
NK_subset <- FindClusters(NK_subset)
NK_subset <- RunUMAP(NK_subset, dims = 1:15)
NK_PlotCl <- DimPlot(NK_subset, reduction = "umap", label = TRUE) + ggtitle('NK Cells') + NoLegend()
NK_PlotID <- DimPlot(NK_subset, reduction = "umap", group.by = "orig.ident", label = FALSE, cols = brewer.pal(4,"Set2"),repel = TRUE) + ggtitle('NK Cells') + NoLegend()

NK_PlotCl + NK_PlotID

SMJvB6_NK <- FindMarkers(NK_subset, ident.1 = c(1),ident.2 = c(0), min.pct = 0.25) #(1) is the B6 cluster; (0) is the SM/J cluster
NK_names <- list('SMJvB6_NK' = SMJvB6_NK)
write.xlsx(NK_names, rowNames = TRUE, file="Path to directory where excel output will be housed/NK.xlsx", append = TRUE)

NK_FeatPlot <- FeaturePlot(NK_subset,c("Cd3e","Cd19","Klrb1c","Nkg7")) 
NK_FeatPlot + NK_PlotCl + NK_PlotID

NKLymph_FeatPlot <- FeaturePlot(NK_subset,c("Klrb1c"))
NKLymph_FeatPlot + NK_PlotCl + NK_PlotID

NK_IfngFeatPlot <- FeaturePlot(NK_subset,c("Ifng")) 
NK_IfngFeatPlot + NK_PlotCl + NK_PlotID

IfngNKcell_subset <- subset(NK_subset,  subset = Ifng >= 1)
length(which(IfngNKcell_subset@meta.data$singlr_labels=="NK cells")) # total=119
length(which(IfngNKcell_subset@meta.data$orig.ident =="B6.12M")) # total=55
length(which(IfngNKcell_subset@meta.data$orig.ident =="B6.20M")) # total=28
length(which(IfngNKcell_subset@meta.data$orig.ident =="SMJ.12M")) # total=17
length(which(IfngNKcell_subset@meta.data$orig.ident =="SMJ.20M")) # total=19


#Analyze just the Dendritic Cells
DC_subset <- subset(scObj_norm,  subset = singlr_labels =="Dendritic cells")
DC_Plot <- DimPlot(DC_subset,reduction = "umap",group.by ='orig.ident',label=T)
DC_subset <- RunPCA(DC_subset, features = VariableFeatures(object = DC_subset))
DimPlot(DC_subset, reduction = "pca")
ElbowPlot(DC_subset)

DC_subset <- FindNeighbors(DC_subset, dims = 1:15)
DC_subset <- FindClusters(DC_subset)
DC_subset <- RunUMAP(DC_subset, dims = 1:15)
DC_PlotCl <- DimPlot(DC_subset, reduction = "umap", label = TRUE) + ggtitle('Dendritic Cells')
DC_PlotID <- DimPlot(DC_subset, reduction = "umap", group.by = "orig.ident", label = FALSE, cols = brewer.pal(4,"Set2"),repel = TRUE) + ggtitle('DC')

DC_PlotCl + DC_PlotID

################################################################

#scRNA-seq data organized in a heatmap according to CyTOF labeling criteria
DoHeatmap(object = scObj_norm, features = c("Cd3e","Cd19","Klrb1c","Cd4","Cd8a","Cd44","Sell","Itgam","Itgax","Ptprc"),group.by = 'singlr_labels', size = 3, angle = 90) + scale_fill_gradientn(colors = c("purple4", "pink", "turquoise"))  + guides(color="none") + 
  theme(text = element_text(size = 15)) #+ ggtitle("CyToF Assignments")

#CyTOF features
CyToF_FeatPlot <- FeaturePlot(scObj_norm,c("Cd3e","Cd19","Klrb1c","Cd4","Cd8a","Itgam","Itgax","Sell","Cd44","Ptprc"))
CyToF_FeatPlot + scObjPlot_SingR + OGIDscObj_normPlot

##Within the same excel sheet created in previous subset analyses, before moving to the next step, save the lists of DEGs you'll want to further analyze
###for example, if you wanted to separate upregulated and downregulated DEGs, you would do that before moving on

###the following examples show GSEA analysis for all (Monocyte1), upregulated Monocyte1+) and downregulated (Monocyte1-) DEGs identified in Monocytes
#Monocytes
Monocytes_df = read.xlsx(xlsxFile = "Path to directory where excel output is housed/Monocyte.xlsx", sheet = "SMJvB6_Monocyte1")
Monocytes_gene_list <- Monocytes_df$avg_log2FC
names(Monocytes_gene_list) <- mapIds(org.Mm.eg.db, keys = Monocytes_df$gene, 'ENTREZID', 'SYMBOL')
Monocytes_gene_list <- na.omit(Monocytes_gene_list)
Monocytes_gene_list <- sort(Monocytes_gene_list,decreasing = TRUE)

Monocytes_gse <- gseGO(geneList=Monocytes_gene_list, 
                       ont ="BP", 
                       keyType = "ENTREZID", 
                       nPerm = 10000, 
                       minGSSize = 3, 
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "none")

as.data.frame(Monocytes_gse)
require(DOSE)
Monocytes_dotplot <- dotplot(Monocytes_gse, showCategory=15, split=".sign", font.size=7) + facet_grid(.~.sign) + ggtitle("Monocytes GSEA SMJ vs. B6")
Monocytes_ridgeplot <- ridgeplot(Monocytes_gse) + labs(x = "enrichment distribution")
Monocytes_cnetplot <- cnetplot(Monocytes_gse, categorySize="pvalue", foldChange=Monocytes_gene_list, showCategory = 3)
Monocytes_gseaplot <- gseaplot(Monocytes_gse, by = "all", title = Monocytes_gse$Description[1], geneSetID = 1)
Monocytes_dotplot

MonocytesUp_df = read.xlsx(xlsxFile = "Path to directory where excel output is housed/Monocyte.xlsx", sheet = "SMJvB6_Monocyte1+")
MonocytesUp_gene_list <- MonocytesUp_df$avg_log2FC
names(MonocytesUp_gene_list) <- mapIds(org.Mm.eg.db, keys = MonocytesUp_df$gene, 'ENTREZID', 'SYMBOL')
MonocytesUp_gene_list <- na.omit(MonocytesUp_gene_list)
MonocytesUp_gene_list <- sort(MonocytesUp_gene_list,decreasing = TRUE)

MonocytesUp_gse <- gseGO(geneList=MonocytesUp_gene_list, 
                         ont ="BP", 
                         keyType = "ENTREZID", 
                         nPerm = 10000, 
                         minGSSize = 3, 
                         maxGSSize = 800, 
                         pvalueCutoff = 0.05, 
                         verbose = TRUE, 
                         OrgDb = organism, 
                         pAdjustMethod = "none")

as.data.frame(MonocytesUp_gse)
require(DOSE)
MonocytesUp_dotplot <- dotplot(MonocytesUp_gse, showCategory=15, split=".sign", font.size=7) + facet_grid(.~.sign) + ggtitle("Monocytes Up GSEA SMJ vs. B6")
MonocytesUp_ridgeplot <- ridgeplot(MonocytesUp_gse) + labs(x = "enrichment distribution")
MonocytesUp_cnetplot <- cnetplot(MonocytesUp_gse, categorySize="pvalue", foldChange=MonocytesUp_gene_list, showCategory = 3)
MonocytesUp_gseaplot <- gseaplot(MonocytesUp_gse, by = "all", title = MonocytesUp_gse$Description[1], geneSetID = 1)



MonocytesDown_df = read.xlsx(xlsxFile = "Path to directory where excel output is housed.xlsx", sheet = "SMJvB6_Monocyte1-")
MonocytesDown_gene_list <- MonocytesDown_df$avg_log2FC
names(MonocytesDown_gene_list) <- mapIds(org.Mm.eg.db, keys = MonocytesDown_df$gene, 'ENTREZID', 'SYMBOL')
MonocytesDown_gene_list <- na.omit(MonocytesDown_gene_list)
MonocytesDown_gene_list <- sort(MonocytesDown_gene_list,decreasing = TRUE)

MonocytesDown_gse <- gseGO(geneList=MonocytesDown_gene_list, 
                           ont ="BP", 
                           keyType = "ENTREZID", 
                           nPerm = 10000, 
                           minGSSize = 3, 
                           maxGSSize = 800, 
                           pvalueCutoff = 0.05, 
                           verbose = TRUE, 
                           OrgDb = organism, 
                           pAdjustMethod = "none")

as.data.frame(MonocytesDown_df)
require(DOSE)
MonocytesDown_dotplot <- dotplot(MonocytesDown_gse, showCategory=15, split=".sign", font.size=7) + facet_grid(.~.sign) + ggtitle("Monocytes Down GSEA SMJ vs. B6")
MonocytesDown_ridgeplot <- ridgeplot(MonocytesDown_gse) + labs(x = "enrichment distribution")
MonocytesDown_cnetplot <- cnetplot(MonocytesDown_gse, categorySize="pvalue", foldChange=MonocytesDown_gene_list, showCategory = 3)
MonocytesDown_gseaplot <- gseaplot(MonocytesDown_gse, by = "all", title = MonocytesDown_gse$Description[1], geneSetID = 1)
