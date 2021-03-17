#!/usr/bin/env perl
use strict;
use Getopt::Long;
my $min_cells=5;
my $min_features=500;
my $filter_doublets=1;
my $use_keeper_enum=0;
my $filter_contaminants_enum=0;
my $filter_organelle_genes=0;
my $dims=20;
my $resolution=0.5;
GetOptions(
	"min_cells=i" => \$min_cells,
	"min_features=i" => \$min_features,
	"filter_doublets!" => \$filter_doublets,
	"filter_organelle_genes!" => \$filter_organelle_genes,
	"use_keeper_enum" => \$use_keeper_enum,
	"filter_contaminants_enum" => \$filter_contaminants_enum,
	"dims=i" => \$dims,
	"resolution=f" => \$resolution,
);
my $params_descriptor=join(",",@ARGV)."-$filter_doublets-$use_keeper_enum-$filter_contaminants_enum-$filter_organelle_genes-$min_cells-$min_features-$dims-$resolution";
my $outfile = $params_descriptor.".R";
my $results_dir = "results-".$params_descriptor;
open(OF, ">$outfile") || die $!;
#prior to starting R, source ~/grant_writing/2018/RH_Libault_PGRP_Resubmission/seurat/env/bin/activate for umap python
#this would be nicer but not working currently
#library(reticulate)
#use_virtualenv("/erdos/adf/libault/paper1/venv")
my $version = '$Revision: 1.4 $';
print OF <<END
#!/usr/bin/env Rscript
#make_seurat $version
library(Seurat)
library(cowplot)
END
;

foreach my $dataset (@ARGV) {
print OF "$dataset.data <-Read10X(data.dir = \"$dataset/filtered_feature_bc_matrix\")\n";
if ($filter_organelle_genes) {
	print OF "$dataset.data <-$dataset.data[grep('AT[1-5]G', rownames($dataset.data)),]\n";
}
print OF "$dataset <- CreateSeuratObject(counts = $dataset.data, project = \"$dataset\", min.cells = $min_cells)\n";
if ($filter_doublets) {
	print OF "doublets_$dataset <- try(read.table(\"$dataset/doublets.txt\"))\n";
	print OF "if (class(doublets_$dataset) == 'try-error') {\n";
	print OF "doublets_$dataset=NULL\n";
	print OF "} else {\n";
	print OF "$dataset <- subset($dataset, cells=WhichCells($dataset, doublets_$dataset\[,1\]), invert=TRUE)\n";
	print OF "}\n";
}
if ($use_keeper_enum) {
	print OF "keepers_$dataset <- try(read.table(\"$dataset/keepers_enum.txt\"))\n";
	print OF "if (class(keepers_$dataset) == 'try-error') {\n";
	print OF "keepers_$dataset=NULL\n";
	print OF "} else {\n";
	print OF "$dataset <- subset($dataset, cells=WhichCells($dataset, keepers_$dataset\[,1\]), invert=FALSE)\n";
	print OF "}\n";
}
if ($filter_contaminants_enum) {
	print OF "contaminants_$dataset <- try(read.table(\"$dataset/contaminants_enum.txt\"))\n";
	print OF "if (class(contaminants_$dataset) == 'try-error') {\n";
	print OF "contaminants_$dataset=NULL\n";
	print OF "} else {\n";
	print OF "$dataset <- subset($dataset, cells=WhichCells($dataset, contaminants_$dataset\[,1\]), invert=TRUE)\n";
	print OF "}\n";
}
print OF "$dataset <- subset($dataset, subset = nFeature_RNA  > $min_features)\n";
print OF "$dataset\[\[\"percent.mt\"\]\] <- PercentageFeatureSet($dataset, pattern = \"ATM\")\n";
print OF "$dataset\$protocol <- \"$dataset\"\n";
print OF "$dataset <- NormalizeData($dataset)\n";
print OF "$dataset <- FindVariableFeatures($dataset, selection.method = \"vst\", nfeatures = 2000)\n";
}

print OF "protocol.anchors <- FindIntegrationAnchors(object.list = list(",join(",",@ARGV),"), dims = 1:$dims)\n";
print OF "protocol.combined <- IntegrateData(anchorset = protocol.anchors, dims = 1:$dims)\n";
print OF "DefaultAssay(protocol.combined) <- \"integrated\"\n";

print OF "dir.create(\"$results_dir\")\n";
#Do some combined plotting pre-clustering to avoid the split by cluster
print OF "pdf(\"$results_dir/vlnPlot-Genes.pdf\")\n";
print OF "VlnPlot(protocol.combined, features=c(\"nFeature_RNA\"), pt.size=0)\n";
print OF "dev.off()\n";
print OF "pdf(\"$results_dir/vlnPlot-UMI.pdf\")\n";
print OF "VlnPlot(protocol.combined, features=c(\"nCount_RNA\"), pt.size=0)\n";
print OF "dev.off()\n";
print OF "pdf(\"$results_dir/vlnPlot-Mito.pdf\")\n";
print OF "VlnPlot(protocol.combined, features=c(\"percent.mt\"), pt.size=0)\n";
print OF "dev.off()\n";



# Run the standard workflow for visualization and clustering
print OF "protocol.combined <- ScaleData(protocol.combined, verbose = FALSE)\n";
#note-to-self; before $dims introduced, this value was fixed at 30, while other places now using $dims were 20
print OF "protocol.combined <- RunPCA(protocol.combined, npcs = $dims, verbose = FALSE)\n";
# t-SNE and Clustering
print OF "protocol.combined <- RunUMAP(protocol.combined, reduction = \"pca\", dims = 1:$dims)\n";
print OF "protocol.combined <- FindNeighbors(protocol.combined, reduction = \"pca\", dims = 1:$dims)\n";
print OF "protocol.combined <- FindClusters(protocol.combined, resolution = $resolution)\n";
print OF "pdf(\"$results_dir/splitdimplot-umap.pdf\")\n";
print OF "DimPlot(protocol.combined, reduction = \"umap\", split.by = \"protocol\")\n";
print OF "dev.off()\n";
print OF "write.table(table(protocol.combined\@meta.data\$seurat_clusters, protocol.combined\@meta.data\$orig.ident), \"$results_dir/cluster_counts.txt\", sep=\"\\t\", quote=F)\n";
print OF "write.table(prop.table(table(protocol.combined\@meta.data\$seurat_clusters, protocol.combined\@meta.data\$orig.ident),2), \"$results_dir/cluster_proportions.by_sample.txt\", sep=\"\\t\", quote=F)\n";
print OF "saveRDS(protocol.combined, \"$results_dir/protocol.combined.rds\")\n";
print OF "write.csv(FetchData(protocol.combined, \"ident\"), \"$results_dir/seurat_bc_clustermap.csv\", quote=F)\n";
print OF "write.csv(protocol.combined\@reductions\$umap\@cell.embeddings, file = \"$results_dir/umap.csv\", quote=F)\n";
