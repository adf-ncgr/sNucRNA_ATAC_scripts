#!/usr/bin/env Rscript
#make_seurat $Revision: 1.4 $
library(Seurat)
library(cowplot)
WT1.data <-Read10X(data.dir = "WT1/filtered_feature_bc_matrix")
WT1.data <-WT1.data[grep('AT[1-5]G', rownames(WT1.data)),]
WT1 <- CreateSeuratObject(counts = WT1.data, project = "WT1", min.cells = 5)
doublets_WT1 <- try(read.table("WT1/doublets.txt"))
if (class(doublets_WT1) == 'try-error') {
doublets_WT1=NULL
} else {
WT1 <- subset(WT1, cells=WhichCells(WT1, doublets_WT1[,1]), invert=TRUE)
}
keepers_WT1 <- try(read.table("WT1/keepers_enum.txt"))
if (class(keepers_WT1) == 'try-error') {
keepers_WT1=NULL
} else {
WT1 <- subset(WT1, cells=WhichCells(WT1, keepers_WT1[,1]), invert=FALSE)
}
contaminants_WT1 <- try(read.table("WT1/contaminants_enum.txt"))
if (class(contaminants_WT1) == 'try-error') {
contaminants_WT1=NULL
} else {
WT1 <- subset(WT1, cells=WhichCells(WT1, contaminants_WT1[,1]), invert=TRUE)
}
WT1 <- subset(WT1, subset = nFeature_RNA  > 200)
WT1[["percent.mt"]] <- PercentageFeatureSet(WT1, pattern = "ATM")
WT1$protocol <- "WT1"
WT1 <- NormalizeData(WT1)
WT1 <- FindVariableFeatures(WT1, selection.method = "vst", nfeatures = 2000)
WT2.data <-Read10X(data.dir = "WT2/filtered_feature_bc_matrix")
WT2.data <-WT2.data[grep('AT[1-5]G', rownames(WT2.data)),]
WT2 <- CreateSeuratObject(counts = WT2.data, project = "WT2", min.cells = 5)
doublets_WT2 <- try(read.table("WT2/doublets.txt"))
if (class(doublets_WT2) == 'try-error') {
doublets_WT2=NULL
} else {
WT2 <- subset(WT2, cells=WhichCells(WT2, doublets_WT2[,1]), invert=TRUE)
}
keepers_WT2 <- try(read.table("WT2/keepers_enum.txt"))
if (class(keepers_WT2) == 'try-error') {
keepers_WT2=NULL
} else {
WT2 <- subset(WT2, cells=WhichCells(WT2, keepers_WT2[,1]), invert=FALSE)
}
contaminants_WT2 <- try(read.table("WT2/contaminants_enum.txt"))
if (class(contaminants_WT2) == 'try-error') {
contaminants_WT2=NULL
} else {
WT2 <- subset(WT2, cells=WhichCells(WT2, contaminants_WT2[,1]), invert=TRUE)
}
WT2 <- subset(WT2, subset = nFeature_RNA  > 200)
WT2[["percent.mt"]] <- PercentageFeatureSet(WT2, pattern = "ATM")
WT2$protocol <- "WT2"
WT2 <- NormalizeData(WT2)
WT2 <- FindVariableFeatures(WT2, selection.method = "vst", nfeatures = 2000)
WT3.data <-Read10X(data.dir = "WT3/filtered_feature_bc_matrix")
WT3.data <-WT3.data[grep('AT[1-5]G', rownames(WT3.data)),]
WT3 <- CreateSeuratObject(counts = WT3.data, project = "WT3", min.cells = 5)
doublets_WT3 <- try(read.table("WT3/doublets.txt"))
if (class(doublets_WT3) == 'try-error') {
doublets_WT3=NULL
} else {
WT3 <- subset(WT3, cells=WhichCells(WT3, doublets_WT3[,1]), invert=TRUE)
}
keepers_WT3 <- try(read.table("WT3/keepers_enum.txt"))
if (class(keepers_WT3) == 'try-error') {
keepers_WT3=NULL
} else {
WT3 <- subset(WT3, cells=WhichCells(WT3, keepers_WT3[,1]), invert=FALSE)
}
contaminants_WT3 <- try(read.table("WT3/contaminants_enum.txt"))
if (class(contaminants_WT3) == 'try-error') {
contaminants_WT3=NULL
} else {
WT3 <- subset(WT3, cells=WhichCells(WT3, contaminants_WT3[,1]), invert=TRUE)
}
WT3 <- subset(WT3, subset = nFeature_RNA  > 200)
WT3[["percent.mt"]] <- PercentageFeatureSet(WT3, pattern = "ATM")
WT3$protocol <- "WT3"
WT3 <- NormalizeData(WT3)
WT3 <- FindVariableFeatures(WT3, selection.method = "vst", nfeatures = 2000)
Arabidopsis_root.data <-Read10X(data.dir = "Arabidopsis_root/filtered_feature_bc_matrix")
Arabidopsis_root.data <-Arabidopsis_root.data[grep('AT[1-5]G', rownames(Arabidopsis_root.data)),]
Arabidopsis_root <- CreateSeuratObject(counts = Arabidopsis_root.data, project = "Arabidopsis_root", min.cells = 5)
doublets_Arabidopsis_root <- try(read.table("Arabidopsis_root/doublets.txt"))
if (class(doublets_Arabidopsis_root) == 'try-error') {
doublets_Arabidopsis_root=NULL
} else {
Arabidopsis_root <- subset(Arabidopsis_root, cells=WhichCells(Arabidopsis_root, doublets_Arabidopsis_root[,1]), invert=TRUE)
}
keepers_Arabidopsis_root <- try(read.table("Arabidopsis_root/keepers_enum.txt"))
if (class(keepers_Arabidopsis_root) == 'try-error') {
keepers_Arabidopsis_root=NULL
} else {
Arabidopsis_root <- subset(Arabidopsis_root, cells=WhichCells(Arabidopsis_root, keepers_Arabidopsis_root[,1]), invert=FALSE)
}
contaminants_Arabidopsis_root <- try(read.table("Arabidopsis_root/contaminants_enum.txt"))
if (class(contaminants_Arabidopsis_root) == 'try-error') {
contaminants_Arabidopsis_root=NULL
} else {
Arabidopsis_root <- subset(Arabidopsis_root, cells=WhichCells(Arabidopsis_root, contaminants_Arabidopsis_root[,1]), invert=TRUE)
}
Arabidopsis_root <- subset(Arabidopsis_root, subset = nFeature_RNA  > 200)
Arabidopsis_root[["percent.mt"]] <- PercentageFeatureSet(Arabidopsis_root, pattern = "ATM")
Arabidopsis_root$protocol <- "Arabidopsis_root"
Arabidopsis_root <- NormalizeData(Arabidopsis_root)
Arabidopsis_root <- FindVariableFeatures(Arabidopsis_root, selection.method = "vst", nfeatures = 2000)
Araroot.data <-Read10X(data.dir = "Araroot/filtered_feature_bc_matrix")
Araroot.data <-Araroot.data[grep('AT[1-5]G', rownames(Araroot.data)),]
Araroot <- CreateSeuratObject(counts = Araroot.data, project = "Araroot", min.cells = 5)
doublets_Araroot <- try(read.table("Araroot/doublets.txt"))
if (class(doublets_Araroot) == 'try-error') {
doublets_Araroot=NULL
} else {
Araroot <- subset(Araroot, cells=WhichCells(Araroot, doublets_Araroot[,1]), invert=TRUE)
}
keepers_Araroot <- try(read.table("Araroot/keepers_enum.txt"))
if (class(keepers_Araroot) == 'try-error') {
keepers_Araroot=NULL
} else {
Araroot <- subset(Araroot, cells=WhichCells(Araroot, keepers_Araroot[,1]), invert=FALSE)
}
contaminants_Araroot <- try(read.table("Araroot/contaminants_enum.txt"))
if (class(contaminants_Araroot) == 'try-error') {
contaminants_Araroot=NULL
} else {
Araroot <- subset(Araroot, cells=WhichCells(Araroot, contaminants_Araroot[,1]), invert=TRUE)
}
Araroot <- subset(Araroot, subset = nFeature_RNA  > 200)
Araroot[["percent.mt"]] <- PercentageFeatureSet(Araroot, pattern = "ATM")
Araroot$protocol <- "Araroot"
Araroot <- NormalizeData(Araroot)
Araroot <- FindVariableFeatures(Araroot, selection.method = "vst", nfeatures = 2000)
snRNA_Arabidopsis.data <-Read10X(data.dir = "snRNA_Arabidopsis/filtered_feature_bc_matrix")
snRNA_Arabidopsis.data <-snRNA_Arabidopsis.data[grep('AT[1-5]G', rownames(snRNA_Arabidopsis.data)),]
snRNA_Arabidopsis <- CreateSeuratObject(counts = snRNA_Arabidopsis.data, project = "snRNA_Arabidopsis", min.cells = 5)
doublets_snRNA_Arabidopsis <- try(read.table("snRNA_Arabidopsis/doublets.txt"))
if (class(doublets_snRNA_Arabidopsis) == 'try-error') {
doublets_snRNA_Arabidopsis=NULL
} else {
snRNA_Arabidopsis <- subset(snRNA_Arabidopsis, cells=WhichCells(snRNA_Arabidopsis, doublets_snRNA_Arabidopsis[,1]), invert=TRUE)
}
keepers_snRNA_Arabidopsis <- try(read.table("snRNA_Arabidopsis/keepers_enum.txt"))
if (class(keepers_snRNA_Arabidopsis) == 'try-error') {
keepers_snRNA_Arabidopsis=NULL
} else {
snRNA_Arabidopsis <- subset(snRNA_Arabidopsis, cells=WhichCells(snRNA_Arabidopsis, keepers_snRNA_Arabidopsis[,1]), invert=FALSE)
}
contaminants_snRNA_Arabidopsis <- try(read.table("snRNA_Arabidopsis/contaminants_enum.txt"))
if (class(contaminants_snRNA_Arabidopsis) == 'try-error') {
contaminants_snRNA_Arabidopsis=NULL
} else {
snRNA_Arabidopsis <- subset(snRNA_Arabidopsis, cells=WhichCells(snRNA_Arabidopsis, contaminants_snRNA_Arabidopsis[,1]), invert=TRUE)
}
snRNA_Arabidopsis <- subset(snRNA_Arabidopsis, subset = nFeature_RNA  > 200)
snRNA_Arabidopsis[["percent.mt"]] <- PercentageFeatureSet(snRNA_Arabidopsis, pattern = "ATM")
snRNA_Arabidopsis$protocol <- "snRNA_Arabidopsis"
snRNA_Arabidopsis <- NormalizeData(snRNA_Arabidopsis)
snRNA_Arabidopsis <- FindVariableFeatures(snRNA_Arabidopsis, selection.method = "vst", nfeatures = 2000)
arabido02_snRNA.data <-Read10X(data.dir = "arabido02_snRNA/filtered_feature_bc_matrix")
arabido02_snRNA.data <-arabido02_snRNA.data[grep('AT[1-5]G', rownames(arabido02_snRNA.data)),]
arabido02_snRNA <- CreateSeuratObject(counts = arabido02_snRNA.data, project = "arabido02_snRNA", min.cells = 5)
doublets_arabido02_snRNA <- try(read.table("arabido02_snRNA/doublets.txt"))
if (class(doublets_arabido02_snRNA) == 'try-error') {
doublets_arabido02_snRNA=NULL
} else {
arabido02_snRNA <- subset(arabido02_snRNA, cells=WhichCells(arabido02_snRNA, doublets_arabido02_snRNA[,1]), invert=TRUE)
}
keepers_arabido02_snRNA <- try(read.table("arabido02_snRNA/keepers_enum.txt"))
if (class(keepers_arabido02_snRNA) == 'try-error') {
keepers_arabido02_snRNA=NULL
} else {
arabido02_snRNA <- subset(arabido02_snRNA, cells=WhichCells(arabido02_snRNA, keepers_arabido02_snRNA[,1]), invert=FALSE)
}
contaminants_arabido02_snRNA <- try(read.table("arabido02_snRNA/contaminants_enum.txt"))
if (class(contaminants_arabido02_snRNA) == 'try-error') {
contaminants_arabido02_snRNA=NULL
} else {
arabido02_snRNA <- subset(arabido02_snRNA, cells=WhichCells(arabido02_snRNA, contaminants_arabido02_snRNA[,1]), invert=TRUE)
}
arabido02_snRNA <- subset(arabido02_snRNA, subset = nFeature_RNA  > 200)
arabido02_snRNA[["percent.mt"]] <- PercentageFeatureSet(arabido02_snRNA, pattern = "ATM")
arabido02_snRNA$protocol <- "arabido02_snRNA"
arabido02_snRNA <- NormalizeData(arabido02_snRNA)
arabido02_snRNA <- FindVariableFeatures(arabido02_snRNA, selection.method = "vst", nfeatures = 2000)
arabido05_snRNA.data <-Read10X(data.dir = "arabido05_snRNA/filtered_feature_bc_matrix")
arabido05_snRNA.data <-arabido05_snRNA.data[grep('AT[1-5]G', rownames(arabido05_snRNA.data)),]
arabido05_snRNA <- CreateSeuratObject(counts = arabido05_snRNA.data, project = "arabido05_snRNA", min.cells = 5)
doublets_arabido05_snRNA <- try(read.table("arabido05_snRNA/doublets.txt"))
if (class(doublets_arabido05_snRNA) == 'try-error') {
doublets_arabido05_snRNA=NULL
} else {
arabido05_snRNA <- subset(arabido05_snRNA, cells=WhichCells(arabido05_snRNA, doublets_arabido05_snRNA[,1]), invert=TRUE)
}
keepers_arabido05_snRNA <- try(read.table("arabido05_snRNA/keepers_enum.txt"))
if (class(keepers_arabido05_snRNA) == 'try-error') {
keepers_arabido05_snRNA=NULL
} else {
arabido05_snRNA <- subset(arabido05_snRNA, cells=WhichCells(arabido05_snRNA, keepers_arabido05_snRNA[,1]), invert=FALSE)
}
contaminants_arabido05_snRNA <- try(read.table("arabido05_snRNA/contaminants_enum.txt"))
if (class(contaminants_arabido05_snRNA) == 'try-error') {
contaminants_arabido05_snRNA=NULL
} else {
arabido05_snRNA <- subset(arabido05_snRNA, cells=WhichCells(arabido05_snRNA, contaminants_arabido05_snRNA[,1]), invert=TRUE)
}
arabido05_snRNA <- subset(arabido05_snRNA, subset = nFeature_RNA  > 200)
arabido05_snRNA[["percent.mt"]] <- PercentageFeatureSet(arabido05_snRNA, pattern = "ATM")
arabido05_snRNA$protocol <- "arabido05_snRNA"
arabido05_snRNA <- NormalizeData(arabido05_snRNA)
arabido05_snRNA <- FindVariableFeatures(arabido05_snRNA, selection.method = "vst", nfeatures = 2000)
protocol.anchors <- FindIntegrationAnchors(object.list = list(WT1,WT2,WT3,Arabidopsis_root,Araroot,snRNA_Arabidopsis,arabido02_snRNA,arabido05_snRNA), dims = 1:50)
protocol.combined <- IntegrateData(anchorset = protocol.anchors, dims = 1:50)
DefaultAssay(protocol.combined) <- "integrated"
dir.create("results-WT1,WT2,WT3,Arabidopsis_root,Araroot,snRNA_Arabidopsis,arabido02_snRNA,arabido05_snRNA-1-1-1-1-5-200-50-0.5")
pdf("results-WT1,WT2,WT3,Arabidopsis_root,Araroot,snRNA_Arabidopsis,arabido02_snRNA,arabido05_snRNA-1-1-1-1-5-200-50-0.5/vlnPlot-Genes.pdf")
VlnPlot(protocol.combined, features=c("nFeature_RNA"), pt.size=0)
dev.off()
pdf("results-WT1,WT2,WT3,Arabidopsis_root,Araroot,snRNA_Arabidopsis,arabido02_snRNA,arabido05_snRNA-1-1-1-1-5-200-50-0.5/vlnPlot-UMI.pdf")
VlnPlot(protocol.combined, features=c("nCount_RNA"), pt.size=0)
dev.off()
pdf("results-WT1,WT2,WT3,Arabidopsis_root,Araroot,snRNA_Arabidopsis,arabido02_snRNA,arabido05_snRNA-1-1-1-1-5-200-50-0.5/vlnPlot-Mito.pdf")
VlnPlot(protocol.combined, features=c("percent.mt"), pt.size=0)
dev.off()
protocol.combined <- ScaleData(protocol.combined, verbose = FALSE)
protocol.combined <- RunPCA(protocol.combined, npcs = 50, verbose = FALSE)
protocol.combined <- RunUMAP(protocol.combined, reduction = "pca", dims = 1:50)
protocol.combined <- FindNeighbors(protocol.combined, reduction = "pca", dims = 1:50)
protocol.combined <- FindClusters(protocol.combined, resolution = 0.5)
source("results-WT1,WT2,WT3,Arabidopsis_root,Araroot,snRNA_Arabidopsis,arabido02_snRNA,arabido05_snRNA-1-1-1-1-5-200-50-0.5/renaming.R")
pdf("results-WT1,WT2,WT3,Arabidopsis_root,Araroot,snRNA_Arabidopsis,arabido02_snRNA,arabido05_snRNA-1-1-1-1-5-200-50-0.5/splitdimplot-umap.pdf")
DimPlot(protocol.combined, reduction = "umap", split.by = "protocol")
dev.off()
write.table(table(protocol.combined@meta.data$seurat_clusters, protocol.combined@meta.data$orig.ident), "results-WT1,WT2,WT3,Arabidopsis_root,Araroot,snRNA_Arabidopsis,arabido02_snRNA,arabido05_snRNA-1-1-1-1-5-200-50-0.5/cluster_counts.txt", sep="\t", quote=F)
write.table(prop.table(table(protocol.combined@meta.data$seurat_clusters, protocol.combined@meta.data$orig.ident),2), "results-WT1,WT2,WT3,Arabidopsis_root,Araroot,snRNA_Arabidopsis,arabido02_snRNA,arabido05_snRNA-1-1-1-1-5-200-50-0.5/cluster_proportions.by_sample.txt", sep="\t", quote=F)
saveRDS(protocol.combined, "results-WT1,WT2,WT3,Arabidopsis_root,Araroot,snRNA_Arabidopsis,arabido02_snRNA,arabido05_snRNA-1-1-1-1-5-200-50-0.5/protocol.combined.rds")
write.csv(FetchData(protocol.combined, "ident"), "results-WT1,WT2,WT3,Arabidopsis_root,Araroot,snRNA_Arabidopsis,arabido02_snRNA,arabido05_snRNA-1-1-1-1-5-200-50-0.5/seurat_bc_clustermap.csv", quote=F)
write.csv(protocol.combined@reductions$umap@cell.embeddings, file = "results-WT1,WT2,WT3,Arabidopsis_root,Araroot,snRNA_Arabidopsis,arabido02_snRNA,arabido05_snRNA-1-1-1-1-5-200-50-0.5/umap.csv", quote=F)
