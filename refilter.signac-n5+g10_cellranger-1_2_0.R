#!/usr/bin/env Rscript
library(Signac)
library(Seurat)
set.seed(1234)
counts <- Read10X_h5(filename = "N5_cellranger-1_2_0/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "N5_cellranger-1_2_0/singlecell.csv",
  header = TRUE,
  row.names = 1
)
n5 <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)
n5$lib <- 'n5'
fragment.path <- 'N5_cellranger-1_2_0/fragments.tsv.gz'
n5 <- SetFragments(
  object = n5,
  file = fragment.path
)
n5 <- NucleosomeSignal(object = n5, region='1-1-24000000')
n5$pct_reads_in_peaks <- n5$peak_region_fragments / n5$total * 100
n5 <- subset(n5, subset = peak_region_fragments > 1000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & nucleosome_signal < 10)
n5 <- RunTFIDF(n5)
n5 <- FindTopFeatures(n5, min.cutoff = 'q0')
n5 <- RunSVD(
  object = n5,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
n5 <- RunUMAP(object = n5, reduction = 'lsi', dims = 1:30)


counts <- Read10X_h5(filename = "G10_cellranger-1_2_0/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "G10_cellranger-1_2_0/singlecell.csv",
  header = TRUE,
  row.names = 1
)
g10 <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)
g10$lib <- 'g10'
fragment.path <- 'G10_cellranger-1_2_0/fragments.tsv.gz'
g10 <- SetFragments(
  object = g10,
  file = fragment.path
)
g10 <- NucleosomeSignal(object = g10, region='1-1-24000000')
g10$pct_reads_in_peaks <- g10$peak_region_fragments / g10$total * 100
g10 <- subset(g10, subset = peak_region_fragments > 1000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & nucleosome_signal < 10)
g10 <- RenameCells(object = g10, add.cell.id = "B")
g10 <- RunTFIDF(g10)
g10 <- FindTopFeatures(g10, min.cutoff = 'q0')
g10 <- RunSVD(
  object = g10,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
g10 <- RunUMAP(object = g10, reduction = 'lsi', dims = 1:30)

intersecting.regions <- GetIntersectingFeatures(
  object.1 = n5,
  object.2 = g10,
  sep.1 = c(":", "-"),
  sep.2 = c(":", "-")
)



#temporary punt- try Harmony integration
atac_combined <- MergeWithRegions(
  object.1 = n5,
  object.2 = g10,
  sep.1 = c(":", "-"),
  sep.2 = c(":", "-")
)
atac_combined <- FindTopFeatures(atac_combined, min.cutoff = 0)
atac_combined <- RunTFIDF(atac_combined)
atac_combined <- RunSVD(atac_combined, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
atac_combined <- RunUMAP(atac_combined, dims = 1:30, reduction = 'lsi')
p5 <- DimPlot(atac_combined, group.by = 'lib', pt.size = 0.1) + ggplot2::ggtitle("Unintegrated")
library("harmony")
atac_combined <- RunHarmony(
  object = atac_combined,
  group.by.vars = 'lib',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)
atac_combined <- RunUMAP(atac_combined, dims = 1:30, reduction = 'harmony')
p6 <- DimPlot(atac_combined, group.by = 'lib', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")
pdf("refilter.harmony_integration.pdf")
CombinePlots(plots = list(p5, p6), ncol = 2)
dev.off()

#for clustering only with ATAC, per request of reviewer
atac_combined <- FindNeighbors(object = atac_combined, reduction = 'lsi', dims = 1:30)
atac_combined <- FindClusters(object = atac_combined, verbose = FALSE, algorithm = 3)
saveRDS(atac_combined, "refilter.atac_combined.harmonized_no_RNAseq_integration.RDS")
write.table(FetchData(atac_combined, "ident"), "refilter.atac_harmonized_no_RNAseq_integration.clusters.csv", sep=",", quote=F)
pdf("refilter.atac_harmonized_no_RNAseq_integration.clusters.pdf")
DimPlot(
object=atac_combined,
 group.by = 'ident',
  label = TRUE,
  repel = TRUE) + ggplot2::ggtitle('scATAC-seq self-clustering')
dev.off()

#https://davetang.org/muse/2017/08/04/read-gtf-file-r/
library(refGenome)
ens <- ensemblGenome()
read.gtf(ens, "../libault-atac/Arabidopsis_thaliana.TAIR10.41.useLocusIDs.gtf")
library(GenomicRanges)
my_gene <- getGenePositions(ens)
my_gr <- with(my_gene, GRanges(seqid, IRanges(start, end), strand, id = gene_id))
genebodyandpromoter.coords <- Extend(x = my_gr, upstream = 2000, downstream = 0)


gene.activities <- FeatureMatrix(
  #fragments = fragment.path,
  fragments =  'aggr-N5_G10-cellranger-1_2_0/renamed.fragments.tsv.gz',
  features = genebodyandpromoter.coords,
  #cells = colnames(n5),
  cells = colnames(atac_combined),
  chunk = 10
)
gene.key <- genebodyandpromoter.coords$id
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]
atac_combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
atac_combined <- NormalizeData(
  object = atac_combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac_combined$nCount_RNA)
)
DefaultAssay(atac_combined)<-'RNA'
nuc <- readRDS("../results-WT1,WT2,WT3,Arabidopsis_root,Araroot,snRNA_Arabidopsis,arabido02_snRNA,arabido05_snRNA-1-1-1-1-5-200-50-0.5/protocol.combined.rds")
transfer.anchors <- FindTransferAnchors(
  reference = nuc,
  query = atac_combined,
reduction = 'cca'
)
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  #refdata = nuc$celltype,
  #refdata = nuc$seurat_clusters,
  refdata = nuc$seurat_clusters_renamed,
  weight.reduction = atac_combined[['lsi']]
)
#this ensures that the combined plot gives the same colors to the same clusters, otherwise numeric vs lexical sorting throws things off
#predicted.labels$predicted.id <- strtoi(predicted.labels$predicted.id)
atac_combined <- AddMetaData(object = atac_combined,  metadata = predicted.labels)
saveRDS(atac_combined, "refilter.atac_combined.RDS")

write.table(FetchData(atac_combined, "predicted.id"), "refilter.n5+g10_cellranger-1_2_0.pred_labels.csv")
library(ggplot2)
plot1 <- DimPlot(
  object = nuc,
  #group.by = 'seurat_clusters',
  group.by = 'seurat_clusters_renamed',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(
object=atac_combined,
 group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
pdf("refilter.n5+g10_cellranger-1_2_0-integrated.pdf")
CombinePlots(list(plot1,plot2), ncol = 2)
dev.off()
write.csv(atac_combined@reductions$umap@cell.embeddings, file="refilter.n5+g10_cellranger-1_2_0-umap.csv")

write.table(FetchData(atac_combined, "ident"), "refilter.n5+g10_cellranger-1_2_0-harmony_integrated_clustermap")


#sc<-subset(nuc, subset = protocol == 'WT1' | protocol == 'WT2' | protocol == 'WT3')
#sn<-subset(nuc, subset = protocol != 'WT1' & protocol != 'WT2' & protocol != 'WT3')
#sn.markers <- FindAllMarkers(sn)
#sc.markers <- FindAllMarkers(sc)
#write.table(sn.markers, "sn.markers.txt", sep="\t", quote=F)
#write.table(sc.markers, "sc.markers.txt", sep="\t", quote=F)

#correlation matrix per https://github.com/satijalab/seurat/issues/1552
Idents(nuc)<-"protocol"
nuc<-RenameIdents(object=nuc, 'WT1' = 'sc', 'WT2' = 'sc', 'WT3' = 'sc', 'Arabidopsis_root' = 'sn', 'Araroot'='sn', 'snRNA_Arabidopsis'='sn', 'arabido02_snRNA'='sn', 'arabido05_snRNA'='sn')
#av.exp<-AverageExpression(nuc, add.ident='seurat_clusters')$RNA
av.exp<-AverageExpression(nuc, add.ident='seurat_clusters_renamed')$RNA
write.table(av.exp, "refilter.av.exp", sep="\t", quote=F)

Idents(atac_combined)<-predicted.labels
av.peak<-AverageExpression(atac_combined)$peaks
write.table(av.peak, "refilter.av.peak", sep="\t",quote=F)
q()


#TSS plots
tss.ranges <- GRanges(seqnames=seqnames(my_gr), ranges = IRanges(start = start(my_gr), width = 2),strand = strand(my_gr))
#get rid of organelles
seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')
seqlevelsStyle(tss.ranges) <- 'Ensembl'

atac_combined<-SetFragments(object=atac_combined,file='aggr-N5_G10-cellranger-1_2_0/renamed.fragments.tsv.gz')
atac_combined<- TSSEnrichment(object = atac_combined, tss.positions = tss.ranges[seqnames(tss.ranges) == 2])
library(ggplot2)

pdf("tss_plot.chr2.pdf")
TSSPlot(atac_combined) + ggtitle("TSS enrichment score: chr2") + NoLegend()
dev.off()

#TES plots
library(Seurat)
library(Signac)
library(ggplot2)
library(refGenome)
atac_combined=readRDS("refilter.atac_combined.RDS")
atac_combined<-SetFragments(object=atac_combined,file='aggr-N5_G10-cellranger-1_2_0/renamed.fragments.tsv.gz')
ens <- ensemblGenome()
read.gtf(ens, "../libault-atac/Arabidopsis_thaliana.TAIR10.41.useLocusIDs.gtf")
library(GenomicRanges)
my_gene <- getGenePositions(ens)
my_gr <- with(my_gene, GRanges(seqid, IRanges(start, end), strand, id = gene_id))
tes.ranges <- GRanges(seqnames=seqnames(my_gr), ranges = IRanges(start = end(my_gr), width = 2),strand = strand(my_gr))
#get rid of organelles
seqlevelsStyle(tes.ranges) <- 'UCSC'
tes.ranges <- keepStandardChromosomes(tes.ranges, pruning.mode = 'coarse')
seqlevelsStyle(tes.ranges) <- 'Ensembl'

atac_combined<- TSSEnrichment(object = atac_combined, tss.positions = tes.ranges[seqnames(tes.ranges) == 1])
pdf("tes_plot.chr1.pdf")
TSSPlot(atac_combined) + ggtitle("TES enrichment score: chr1") + NoLegend()
dev.off()



pdf("atac_combined.TSS_Fragments.pdf")
FeaturePlot(atac_combined, features=c('TSS_fragments'), pt.size=0.5)
dev.off()

pdf("fragmenthist.all.predicted.id.pdf")
FragmentHistogram(object = atac_combined, group.by = 'predicted.id.pdf', region=c("1-1-","2-1-","3-1-",4-1-",5-1-"))
dev.off()

library(Seurat)
library(Signac)
library(ggplot2)
atac_combined<-readRDS("refilter.atac_combined.RDS")
atac_combined<-SetFragments(object=atac_combined,file='aggr-N5_G10-cellranger-1_2_0/renamed.fragments.tsv.gz')
atac_combined<-SetIdent(atac_combined, value='predicted.id')
atac_combined <- RenameIdents(object = atac_combined,
'a14_1' = '1',
'b5_2' = '2',
'c13_3' = '3',
'd12_4' = '4',
'e2_5' = '5',
'f8_6' = '6',
'g4_7' = '7',
'h6_8' = '8',
'i9_9' = '9',
'j20_10' = '10',
'k0_11' = '11',
'l10_12' = '12',
'm11_13' = '13',
'n17_14' = '14',
'o1_15' = '15',
'p16_16' = '16',
'q7_17' = '17',
'r3_18' = '18',
's18_19' = '19',
't19_20' = '20',
'u15_21' = '21'
)
atac_combined[['seurat_clusters_renamed_again']] <- Idents(object = atac_combined)
pdf("fragmenthist.all.predicted.id.colors.pdf",width=14)
#FragmentHistogram(object = atac_combined, group.by = 'seurat_clusters_renamed_again', region=c("1-1-30427671","2-1-19698289","3-1-23459830","4-1-18585056","5-1-26975502"))+scale_fill_manual(values=c('#5B9BD5', '#296AB5', '#09C8F7', '#10E858', '#6DBE13', '#96F125', '#008E06', '#003899', '#747474', '#4A4A4A', '#8E5D0C', '#5F5FE8', '#AC8CD7', '#F1A3B5', '#E6CFD4', '#FF33CC', '#B31C46', '#F4C8A1', '#FB8E01', '#BE9146', '#FA5A0F'))+NoLegend()
FragmentHistogram(object = atac_combined, group.by = 'seurat_clusters_renamed_again', region=c("1-1-30427671","2-1-19698289","3-1-23459830","4-1-18585056","5-1-26975502"))+scale_fill_manual(values=c('#5B9BD5', '#296AB5', '#09C8F7', '#10E858', '#6DBE13', '#96F125', '#008E06', '#747474', '#4A4A4A', '#8E5D0C', '#5F5FE8', '#AC8CD7', '#F1A3B5', '#E6CFD4', '#FF33CC', '#B31C46', '#F4C8A1', '#FB8E01', '#BE9146', '#FA5A0F', '#BB4918'))+NoLegend()
dev.off()
