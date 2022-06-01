## mouse_seurat_interated_v1.R

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

## integrating our visium spatial transcriptomics data from mouse placenta with Marsh, et al. (2020) snRNA-seq from mouse placenta.

## Analysis of Marsh & Blelloch, eLife, 2020 murine placental snRNA-seq
# GSE152248 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152248
## E9.5-E14.5: E9.5, E10.5, E12.5, E14.5

## 17500 nuclei, 10X genomics 3' v3 chemistry, 20-30k reads/nucleus
## Filters: >500 <4000 nFeatures >25% mito
## 27,326 nuclei after QC/ 26 clusters (0.6res), XIST for male embryonic tissue
## 5 broad cell types, further divided into subgroups, 16,386 trophoblasts reclustered



## Analyze data on AagaardLab3
#

set.seed(seed=1)
setwd("/home/ebarrozo/visium/results")
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
library(SeuratDisk)
library(dplyr)
library(patchwork)


setwd("/home/ebarrozo/visium/results/seurat_mouse_v1")
load("/home/ebarrozo/visium/results/seurat_mouse_v1/murine_spatial_data-integrated-annotated_v1.RData")
	## see ## visium_mouse_seurat_v1.R

load("/home/ebarrozo/Marsh/results/seurat-v1/marsh_annotated_v1.RData")
	## see # Marsh_snRNA-seq_v2.R
DefaultAssay(seurat.object) <- "SCT"
DefaultAssay(marsh.seurat.object) <- "SCT"


dir.create("integration-v1")
setwd("integration-v1")
setwd("/home/ebarrozo/visium/results/seurat_mouse_v1/integration-v1")

############################################################# ############################################################# #############################################################
############################################################# Integration: merge all data using CCA anchors #############################################################
############################################################# ############################################################# #############################################################
## combine lists of genes
all.genes.combined <- union(all.genes, all.genes.marsh)
	## 33112 combined variable features

## combine lists of top variable features for later DE analysis/clustering
var.combined <- union(var.combined, marsh.var)
	## 8128 combined variable features

## Merge objects
all_merged <- merge(x = seurat.object, y = c(marsh.seurat.object), merge.data = TRUE, project = "all_merged")
rm(marsh.seurat.object)
rm(seurat.object)

ncol(all_merged)
	# 9701 spots/cells
options(future.globals.maxSize = 300000000000)

all_merged <- SCTransform(all_merged, assay = "SCT", verbose = T)
all_merged.list <- SplitObject(all_merged, split.by = "orig.ident")

# Make a reference list
reference.list <- all_merged.list[c("Mock-enox-1", "Mock-enox-2", "ZIKV-sham", "ZIKV-enox", "Uninfected")]

rm("all_merged")
rm(all_merged.list)

# Select transcripts used for data integration
reference.features <- SelectIntegrationFeatures(object.list = reference.list, nfeatures = 3000)
	# No variable features found for object1 in the object.list. Running FindVariableFeatures ...
	# Error: SCT assay is comprised of multiple SCT models. To change the variable features, please set manually with VariableFeatures<-
		## re-run SCTransform on the all_merged object
write.table(reference.features, "all_merged.reference.features.txt", sep="\t")


reference.list <- PrepSCTIntegration(object.list = reference.list, anchor.features = reference.features, verbose = T)

# Find integration anchors, this takes about 15-20 mins
all_merged.anchors <- FindIntegrationAnchors(object.list = reference.list, normalization.method = "SCT", anchor.features = reference.features, verbose = T)
## If there is an object with few cells/spots (<50?), you will see the error message below and need to change the k.filter to 1 fewer than tne lowest number of cells/spots you have
	# You're computing too large a percentage of total singular values, use a standard svd instead.
	# Warning messages:
		# 1: In FilterAnchors(object = object.pair, assay = assay, slot = slot,  :
  		# Number of anchor cells is less than k.filter. Retaining all anchors.
			# all_merged.anchors <- FindIntegrationAnchors(object.list = reference.list, reduction="rpca", normalization.method = "SCT", anchor.features = reference.features, verbose = TRUE)
					## Projecting new data onto SVD
rm(reference.features)
rm(reference.list)

seurat.object <- IntegrateData(anchorset = all_merged.anchors, normalization.method = "SCT")
	# Merging dataset 3 into 1
	#Extracting anchors for merged samples
	#Finding integration vectors
	#Finding integration vector weights
	#Error in idx[i, ] <- res[[i]][[1]] : 
	#  number of items to replace is not a multiple of replacement length
		#all_merged.integrated <- IntegrateData(anchorset = all_merged.anchors, normalization.method = "SCT", k.weight = 9)
			# Usually the k.weight is set to 100 by default, but here the smallest dataset only has 46 cells. Based on searching with the error message I found this shouldn't be a problem https://github.com/satijalab/seurat/issues/3930 https://github.com/satijalab/seurat/issues/1472

rm("all_merged.anchors")

# Add all viral+human genes as a list of features to call upon during clustering or differential expression of all genes, which may be compuationally demanding
	# all.genes <- rownames(seurat.object@assays[["Spatial"]])
		# combined all.genes from both objects above
all.genes <- all.genes.combined

# Use the 'integrated' assay for clustering and the SCT assay for differential expression
DefaultAssay(seurat.object) <- "integrated"

# Add all integrated genes as a list of features to call upon during clustering (or quick- differential expression)
integrated.genes <- rownames(seurat.object)
write.table(integrated.genes, "integrated.genes.txt", sep="\t")

############################################################# ############################################################# #############################################################
############################################################# Integrated: Dimension Reduction and Visualization #############################################################
############################################################# ############################################################# #############################################################

seurat.object <- ScaleData(seurat.object, features = integrated.genes, assay = "integrated", vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"))
seurat.object<- RunPCA(seurat.object, features = integrated.genes, assay = "integrated", slot = "scale.data")

seurat.object <- FindNeighbors(seurat.object, features = "integrated.genes", dims = 1:30)
seurat.object <- FindClusters(seurat.object, resolution = 0.6)
seurat.object <- RunUMAP(seurat.object, dims=1:30)
	## clusters 0 thru 16
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident")
p1
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE")
p2
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type")
p3
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Subtype")
p4
umap.combined <- CombinePlots(plots = list(p1, p2, p3, p4))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("integrated_UMAP-combined.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("integrated_UMAP-combined-wide.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "orig.ident", label=TRUE)
ggsave("integrated_UMAP_splitby_libID.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", split.by = "orig.ident", label=TRUE)
ggsave("integrated_UMAP_splitby_libID_Type.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Subtype", split.by = "orig.ident", label=TRUE)
ggsave("integrated_UMAP_splitby_libID_Subtype.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "seurat_clusters", label=TRUE)
ggsave("integrated_UMAP_splitby_clusters.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
#examine UMAPs with qc metrics 
p5 <- FeaturePlot(seurat.object, features = 'nFeature_Spatial')
p6 <- FeaturePlot(seurat.object, features = 'nCount_Spatial')
p7 <- FeaturePlot(seurat.object, features = 'percent.mt')
p8 <- FeaturePlot(seurat.object, features = 'percent.ribo')
p9 <- FeaturePlot(seurat.object, features = 'nFeature_RNA')
p10 <- FeaturePlot(seurat.object, features = 'nCount_RNA')
umap.combined <- CombinePlots(plots = list(p5, p9, p6, p10, p7, p8), ncol=2)
ggsave("integrated_UMAP_QCmetricsFeaturePlots.pdf", plot = umap.combined, device = "pdf", width = 8, height = 10, units = "in")


###### UMAP + Spatial DimPlots

p2 <- SpatialDimPlot(seurat.object, cols=5)
ggsave("Test-spatialdimplot.pdf", plot = p2, device = "pdf", width = 8, height = 20, units = "in")

p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("slice1.1")
p2 <- SpatialDimPlot(seurat.object, images=image.list)
p2
image.list <- c("slice1.2.1.1")
p3 <- SpatialDimPlot(seurat.object, images=image.list)
p3
image.list <- c("slice1.5.2")
p4 <- SpatialDimPlot(seurat.object, images=image.list, pt.size.factor = 3)
p4
image.list <- c("slice1.3.3.3")
p5 <- SpatialDimPlot(seurat.object, images=image.list)
p5
	## Check to confirm they are in the correct order
p6 <- wrap_plots(p1, p4, p2, p3, p5, ncol=5)
ggsave("integrated_UMAP_SpatialDimPlots-cropped.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")

p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("slice1.1")
p2 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p2
image.list <- c("slice1.2.1.1")
p3 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p3
image.list <- c("slice1.5.2")
p4 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p4
image.list <- c("slice1.3.3.3")
p5 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p5
p6 <- wrap_plots(p1, p4, p2, p3, p5, ncol=5)
ggsave("integrated_UMAP_SpatialDimPlots-uncropped.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")

Idents(seurat.object) <- "Type"

p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("slice1.1")
p2 <- SpatialDimPlot(seurat.object, images=image.list)
p2
image.list <- c("slice1.2.1.1")
p3 <- SpatialDimPlot(seurat.object, images=image.list)
p3
image.list <- c("slice1.5.2")
p4 <- SpatialDimPlot(seurat.object, images=image.list, pt.size.factor = 3)
p4
image.list <- c("slice1.3.3.3")
p5 <- SpatialDimPlot(seurat.object, images=image.list)
p5
    ## Check to confirm they are in the correct order
p6 <- wrap_plots(p1, p4, p2, p3, p5, ncol=5)
ggsave("integrated_UMAP_SpatialDimPlots-cropped_type.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")

p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("slice1.1")
p2 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p2
image.list <- c("slice1.2.1.1")
p3 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p3
image.list <- c("slice1.5.2")
p4 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p4
image.list <- c("slice1.3.3.3")
p5 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p5
p6 <- wrap_plots(p1, p4, p2, p3, p5, ncol=5)
ggsave("integrated_UMAP_SpatialDimPlots-uncropped_type.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")


Idents(seurat.object) <- "seurat_clusters"


## Spatial DimPlots split.by clusters
image.list <- c("slice1.5.2")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S07.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.1")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S09.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.2.1.1")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S12.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.3.3.3")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S13.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


image.list <- c("slice1.5.2")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0,5,7,11,12,14)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S07_macs.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.1")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0,5,7,11,12,14)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S09_macs.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.2.1.1")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0,5,7,11,12,14)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S12_macs.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.3.3.3")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0,5,7,11,12,14)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S13_macs.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


image.list <- c("slice1.5.2")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0,1,3,4,5,6,8,9,11,13,14)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S07_NK.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.1")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0,1,3,4,5,6,8,9,11,13,14)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S09_NK.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.2.1.1")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0,1,3,4,5,6,8,9,11,13,14)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S12_NK.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.3.3.3")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0,1,3,4,5,6,8,9,11,13,14)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S13_NK.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")




rm(umap.combined)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(p5)
rm(p6)
rm(p7)
rm(p8)
rm(p9)
rm(p10)
rm(image.list)

############################################################# ############################################################# #############################################################
############################################################# DE Analysis #############################################################
############################################################# ############################################################# #############################################################
		## Unresolved error using FindAllMarkers. None of the assays are working. 
		## Unresolved error using FindAllMarkers. None of the assays are working. 
		## Unresolved error using FindAllMarkers. None of the assays are working. 
> de_markers <- FindAllMarkers(seurat.object, features = var.combined, assay = "RNA", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
Calculating cluster 0
Calculating cluster 1
Calculating cluster 2
Calculating cluster 3
Calculating cluster 4
Calculating cluster 5
Calculating cluster 6
Calculating cluster 7
Calculating cluster 8
Calculating cluster 9
Calculating cluster 10
Calculating cluster 11
Calculating cluster 12
Calculating cluster 13
Calculating cluster 14
Calculating cluster 15
Warning: No DE genes identified
Warning: The following tests were not performed: 
Warning: When testing 0 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 1 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 2 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 3 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 4 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 5 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 6 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 7 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 8 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 9 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 10 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 11 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 12 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 13 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 14 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 15 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
> 
		## Unresolved error using FindAllMarkers. None of the assays are working. 
		## Unresolved error using FindAllMarkers. None of the assays are working. 
		## Unresolved error using FindAllMarkers. None of the assays are working. 


seurat.object <- SCTransform(seurat.object, assay="SCT", new.assay.name = "DE", do.scale=TRUE, return.only.var.genes= FALSE, verbose = T)

DefaultAssay(seurat.object) <- "DE"
DE.var.combined <- seurat.object@assays[["DE"]]@var.features
##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = var.combined, assay="Spatial")

DefaultAssay(seurat.object) <- "RNA"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = var.combined, assay="RNA")

Idents(seurat.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object, features = var.combined, assay = "DE", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "Spatial_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>0.693)
		## all clusters have sig. genes
			# de_markers <- FindAllMarkers(seurat.object, features = var.combined, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
			# top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
			# View(top5)
			# write.table(de_markers, "integrated_DEGs_byclusters_pos-0.2lnFC.txt", sep="\t")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("integrated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("integrated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("integrated_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("integrated_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("integrated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)

image.list <- c("slice1.5")
p2 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, pt.size.factor = 3, ncol=5)
p2
ggsave("top.markers_SpatialFeaturePlots-cropped-S07.pdf", plot = p2, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.1")
p3 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p3
ggsave("top.markers_SpatialFeaturePlots-cropped-S09.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.2.1")
p4 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=7)
p4
ggsave("top.markers_SpatialFeaturePlots-cropped-S12.pdf", plot = p4, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.3.3")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=7)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S13.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)

image.list <- c("slice1.5")
p2 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, pt.size.factor = 3, ncol=5)
p2
ggsave("top.markers_SpatialPlots-cropped-S07.pdf", plot = p2, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.1")
p3 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p3
ggsave("top.markers_SpatialPlots-cropped-S09.pdf", plot = p3, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.2.1")
p4 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=7)
p4
ggsave("top.markers_SpatialPlots-cropped-S12.pdf", plot = p4, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.3.3")
p5 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=7)
p5
ggsave("top.markers_SpatialPlots-cropped-S13.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")



rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# Gene set analysis #############################################################
############################################################# ############################################################# #############################################################
Idents(seurat.object) <- "Type"
seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(seurat.object, features = lymphoid.lineage)
png(file=paste0("lymphoid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Lymphoid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

## canonical markers for myeloid cells https://docs.abcam.com/pdf/immunology/myeloid_cell.pdf
# myeloid.lineage <- c("CD34", "CD117", "CD271", "CD33", "CD71", "CD61", "CD123", "CD44", "CD15", "CD16", "CD11b", "CD14", "CD1c", "CD141")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
myeloid.lineage <- c("CD34", "KIT", "NGFR", "CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
feature.plot <- DotPlot(seurat.object, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(seurat.object, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

feature.plot <- DotPlot(seurat.object, features = ZIKV.genes)
png(file=paste0("ZIKV.genes-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="ZIKV Transcripts") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
library.averages.heatmap <- DoHeatmap(seurat.object, features = ZIKV.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("ZIKV.genes_seurat.clusters_heatmap_counts.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = ZIKV.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("ZIKV.genes_seurat.clusters_heatmap_logfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
rm(library.averages.heatmap)

Idents(seurat.object) <- "Type"
seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell_subtype.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(seurat.object, features = lymphoid.lineage)
png(file=paste0("lymphoid.lineage-DotPlot_subtype.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Lymphoid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

## canonical markers for myeloid cells https://docs.abcam.com/pdf/immunology/myeloid_cell_subtype.pdf
# myeloid.lineage <- c("CD34", "CD117", "CD271", "CD33", "CD71", "CD61", "CD123", "CD44", "CD15", "CD16", "CD11b", "CD14", "CD1c", "CD141")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
myeloid.lineage <- c("CD34", "KIT", "NGFR", "CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
feature.plot <- DotPlot(seurat.object, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot_subtype.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(seurat.object, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot_subtype.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

feature.plot <- DotPlot(seurat.object, features = ZIKV.genes)
png(file=paste0("ZIKV.genes-DotPlot_subtype.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="ZIKV Transcripts") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
library.averages.heatmap <- DoHeatmap(seurat.object, features = ZIKV.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("ZIKV.genes_seurat.clusters_heatmap_counts_subtype.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = ZIKV.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("ZIKV.genes_seurat.clusters_heatmap_logfc_subtype.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
rm(library.averages.heatmap)

############################################################# ############################################################# #############################################################
seurat.object2 <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
seurat.object <- seurat.object2
rm(seurat.object2)


#Save tables with cell counts per library and per cluster
dir.create("cell.counts")
setwd("cell.counts")
# Save a table with cell counts per library
Idents(seurat.object) <- "orig.ident"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "sample.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(seurat.object) <- "seurat_clusters"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "cluster.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(seurat.object) <- "Subtype"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "Subtype.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(seurat.object) <- "Type"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "Type.cellcounts.txt", sep="\t")


dir.create("Type.cluster")
setwd("Type.cluster")

# "Trophoblast", "Fibroblast", "Stromal", "Smooth_muscle", "Endothelial", "Endometrial"
# c("Endometrial", "Endothelial", "Fibroblast", "Stromal", "Villious_cytotrophoblast", "Syncytiotrophoblast", "NK", "Macrophage", "Erythroblast")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Trophoblast"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Trophoblast.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Stromal"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Stromal.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Fibroblast"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Fibroblast.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Endothelial"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Endothelial.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Smooth_muscle"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Smooth_muscle.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Endometrial"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Endometrial.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Villious_cytotrophoblast"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Villious_cytotrophoblast.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Syncytiotrophoblast"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Syncytiotrophoblast.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Macrophage"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Macrophage.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("NK"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "NK.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Erythroblast"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Erythroblast.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")
setwd("..")


# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("Type.cluster/"))
manifest<-manifest%>%mutate(name=gsub(".cluster.cellcounts.txt","",value))
manifest
		## These commands below are required, change the prefix to match your samples
df<-read.table("Type.cluster/Trophoblast.cluster.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("cluster","Trophoblast")
 df2<-df
df<-df%>%pivot_longer(cols = `Trophoblast`,names_to = "Type",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("Type.cluster/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("cluster",X)
  df<-df%>%pivot_longer(cols = X,names_to = "Type",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of clusters ## there are 12 Types
# df2%>%mutate(perc=percs_by_group(count,group = Type))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])
df7<-my_fxn(X = manifest$value[6])
df8<-my_fxn(X = manifest$value[7])
df9<-my_fxn(X = manifest$value[8])
df10<-my_fxn(X = manifest$value[9])
df11<-my_fxn(X = manifest$value[10])
df12<-my_fxn(X = manifest$value[11])

final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)
final<-full_join(final,df7)
final<-full_join(final,df8)
final<-full_join(final,df9)
final<-full_join(final,df10)
final<-full_join(final,df11)
final<-full_join(final,df12)


final<-final%>%
  mutate(Type=gsub(".cluster.cellcounts.txt","",Type))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = Type))
#  pivot_wider(names_from = Type,values_from = count)

write.table(final,"Type.CellCounts_merged.tsv",sep = "\t",row.names = F)

final
​
​
# df3<-sapply(manifest$value, my_fxn)
# 
# manifest<-manifest%>%
#   mutate(results=lapply(FUN = my_fxn(manifest$value)))
# 
#                             
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18

final$cluster<-as.factor(final$cluster)

final$Type<-as.factor(final$Type)
​
final
​
library(ggpubr)
ggbarplot(data = final,
          x = "Type",
          y = "Percent",
          color = "cluster",
          fill = "cluster",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)

png(file=paste0("Type.cluster_diversity.plot.png"),
                res=300, 
                width=1500, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "Type",
          y = "Percent",
          color = "cluster",
          fill = "cluster",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()

####################################################################################################

# reference.list <- all_merged.list[c("Mock-enox-1", "Mock-enox-2", "ZIKV-sham", "ZIKV-enox", "Uninfected")]


dir.create("orig.ident.Type")
setwd("orig.ident.Type")
Idents(seurat.object)<-"orig.ident"
seurat.object2 <- subset(seurat.object, idents = c("Mock-enox-1"))
Idents(seurat.object2)<-"Type"
unfiltered.count <- table(Idents(seurat.object2))
write.table(unfiltered.count, "Mock-enox-1.Type.cellcounts.txt", sep="\t")
Idents(seurat.object)<-"orig.ident"
seurat.object2 <- subset(seurat.object, idents = c("Mock-enox-2"))
Idents(seurat.object2)<-"Type"
unfiltered.count <- table(Idents(seurat.object2))
write.table(unfiltered.count, "Mock-enox-2.Type.cellcounts.txt", sep="\t")
Idents(seurat.object)<-"orig.ident"
seurat.object2 <- subset(seurat.object, idents = c("ZIKV-sham"))
Idents(seurat.object2)<-"Type"
unfiltered.count <- table(Idents(seurat.object2))
write.table(unfiltered.count, "ZIKV-sham.Type.cellcounts.txt", sep="\t")
Idents(seurat.object)<-"orig.ident"
seurat.object2 <- subset(seurat.object, idents = c("ZIKV-enox"))
Idents(seurat.object2)<-"Type"
unfiltered.count <- table(Idents(seurat.object2))
write.table(unfiltered.count, "ZIKV-enox.Type.cellcounts.txt", sep="\t")
Idents(seurat.object)<-"orig.ident"
seurat.object2 <- subset(seurat.object, idents = c("Uninfected"))
Idents(seurat.object2)<-"Type"
unfiltered.count <- table(Idents(seurat.object2))
write.table(unfiltered.count, "Uninfected.Type.cellcounts.txt", sep="\t")
rm(seurat.object2)
setwd("..")


library(tidyverse)
library(mosaic)
manifest<-as_tibble(list.files("orig.ident.Type/"))
manifest<-manifest%>%mutate(name=gsub(".cluster.cellcounts.txt","",value))
manifest

## manually change file names to make sure it works before running function

df<-read.table("orig.ident.Type/Mock-enox-1.Type.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type","AF-p1")
df2<-df
df<-df%>%pivot_longer(cols = `Mock-enox-1`,names_to = "orig.ident",values_to = "Count")
df
manifest


my_fxn<-function(X){
  df<-read.table(paste0("orig.ident.Type/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Type",X)
  df<-df%>%pivot_longer(cols = X,names_to = "orig.ident",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

## this is required too
df2<-df

## make sure this matches the number of clusters ## there are 22 SubTypes
# df2%>%mutate(perc=percs_by_group(count,group = SubType))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])


​
final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)



​final
final<-final%>%mutate(orig.ident=gsub(".Type.cellcounts.txt","",orig.ident))
final<-final%>%mutate(Percent=percs_by_group(Count,group = orig.ident))
#  pivot_wider(names_from = SubType,values_from = count)
​
write.table(final,"orig.ident.Type.CellCounts_merged.tsv",sep = "\t",row.names = F)

final
​
​
# df3<-sapply(manifest$value, my_fxn)
# 
# manifest<-manifest%>%
#   mutate(results=lapply(FUN = my_fxn(manifest$value)))
# 
#                             
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.

## make sure you have enough colors: here, 22
mypal<-get_palette("ucscgb",23)


# final$cluster<-as.factor(final$cluster)
final$orig.ident<-as.factor(final$orig.ident)
final$Type<-as.factor(final$Type)
​
final
​
library(ggpubr)
ggbarplot(data = final,
          x = "orig.ident",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)

png(file=paste0("orig.ident.Type_diversity.plot.png"),
                res=300, 
                width=3000, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "orig.ident",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()



dir.create("orig.ident.cluster")
setwd("orig.ident.cluster")
# reference.list <- all_merged.list[c("Mock-enox-1", "Mock-enox-2", "ZIKV-sham", "ZIKV-enox", "Uninfected")]


Idents(seurat.object)<-"orig.ident"
clusters.0.uninfected <- subset(seurat.object, idents = c("Mock-enox-1"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.unfilcount <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Mock-enox-1.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"orig.ident"
clusters.0.uninfected <- subset(seurat.object, idents = c("Mock-enox-2"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Mock-enox-2.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"orig.ident"
clusters.0.uninfected <- subset(seurat.object, idents = c("ZIKV-sham"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "ZIKV-sham.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"orig.ident"
clusters.0.uninfected <- subset(seurat.object, idents = c("ZIKV-enox"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "ZIKV-enox.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"orig.ident"
clusters.0.uninfected <- subset(seurat.object, idents = c("Uninfected"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Uninfected.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

setwd("..")


library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("orig.ident.cluster/"))
manifest<-manifest%>%mutate(name=gsub(".cluster.cellcounts.txt","",value))
manifest
        ## These commands below are required, change the prefix to match your samples
df<-read.table("orig.ident.cluster/Mock-enox-1.cluster.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("cluster","Mock-enox-1")
 df2<-df
df<-df%>%pivot_longer(cols = `Mock-enox-1`,names_to = "orig.ident",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("orig.ident.cluster/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("cluster",X)
  df<-df%>%pivot_longer(cols = X,names_to = "orig.ident",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of clusters ## there are 12 orig.idents
# df2%>%mutate(perc=percs_by_group(count,group = orig.ident))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])
final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)

final<-final%>%
  mutate(orig.ident=gsub(".cluster.cellcounts.txt","",orig.ident))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = orig.ident))
#  pivot_wider(names_from = orig.ident,values_from = count)

write.table(final,"orig.ident.CellCounts_merged.tsv",sep = "\t",row.names = F)

final
​
​
# df3<-sapply(manifest$value, my_fxn)
# 
# manifest<-manifest%>%
#   mutate(results=lapply(FUN = my_fxn(manifest$value)))
# 
#                             
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18

final$cluster<-as.factor(final$cluster)

final$orig.ident<-as.factor(final$orig.ident)
​
final
​
library(ggpubr)
ggbarplot(data = final,
          x = "orig.ident",
          y = "Percent",
          color = "cluster",
          fill = "cluster",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)

png(file=paste0("orig.ident.cluster_diversity.plot.png"),
                res=300, 
                width=1500, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "orig.ident",
          y = "Percent",
          color = "cluster",
          fill = "cluster",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()





dir.create("seurat_clusters.Type")
setwd("seurat_clusters.Type")

# "Trophoblast", "Fibroblast", "Stromal", "Smooth_muscle", "Endothelial", "Endometrial"
# c("Endometrial", "Endothelial", "Fibroblast", "Stromal", "Villious_cytotrophoblast", "Syncytiotrophoblast", "NK", "Macrophage", "Erythroblast")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("0"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "0.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("1"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "1.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("2"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "2.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("3"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "3.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("4"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "4.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("5"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "5.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("6"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "6.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("7"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "7.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("8"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "8.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("9"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "9.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("10"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "10.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("11"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "11.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("12"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "12.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("13"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "13.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("14"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "14.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("15"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "15.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")


setwd("..")
# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("seurat_clusters.Type/"))
manifest<-manifest%>%mutate(name=gsub(".Type.cellcounts.txt","",value))
manifest
df<-read.table("seurat_clusters.Type/1.Type.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type","1")
 df2<-df
df<-df%>%pivot_longer(cols = `1`,names_to = "seurat_clusters",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("seurat_clusters.Type/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Type",X)
  df<-df%>%pivot_longer(cols = X,names_to = "seurat_clusters",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of Types ## there are 12 seurat_clusterss
# df2%>%mutate(perc=percs_by_group(count,group = seurat_clusters))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])
df7<-my_fxn(X = manifest$value[6])
df11<-my_fxn(X = manifest$value[10])
df12<-my_fxn(X = manifest$value[11])
df13<-my_fxn(X = manifest$value[12])
df14<-my_fxn(X = manifest$value[13])
df15<-my_fxn(X = manifest$value[14])
df16<-my_fxn(X = manifest$value[15])
df17<-my_fxn(X = manifest$value[16])
final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)
final<-full_join(final,df7)
final<-full_join(final,df8)
final<-full_join(final,df9)
final<-full_join(final,df10)
final<-full_join(final,df11)
final<-full_join(final,df12)
final<-full_join(final,df13)
final<-full_join(final,df14)
final<-full_join(final,df15)
final<-full_join(final,df16)
final<-full_join(final,df17)


final<-final%>%
  mutate(seurat_clusters=gsub(".Type.cellcounts.txt","",seurat_clusters))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = seurat_clusters))
#  pivot_wider(names_from = seurat_clusters,values_from = count)

write.table(final,"seurat_clusters.CellCounts_merged.tsv",sep = "\t",row.names = F)

final
​
​
# df3<-sapply(manifest$value, my_fxn)
# 
# manifest<-manifest%>%
#   mutate(results=lapply(FUN = my_fxn(manifest$value)))
# 
#                             
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18

final$Type<-as.factor(final$Type)

final$seurat_clusters<-as.factor(final$seurat_clusters)
​
final
​
library(ggpubr)
ggbarplot(data = final,
          x = "seurat_clusters",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)

png(file=paste0("seurat_clusters.Type_diversity.plot.png"),
                res=300, 
                width=2500, 
                height=1500)
feature.plot <- ggbarplot(data = final,
          x = "seurat_clusters",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()

####################################################################################################

####### The clusters are not partitioning clearly by type. Try to increase the resolution and see if that helps?? 

dir.create("res.up")
setwd("res.up")

# seurat.object <- ScaleData(seurat.object, features = integrated.genes, assay = "integrated", vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"))
seurat.object<- RunPCA(seurat.object, features = integrated.genes, assay = "integrated", slot = "scale.data")

seurat.object <- FindNeighbors(seurat.object, features = "integrated.genes", dims = 1:30)
seurat.object <- FindClusters(seurat.object, resolution = 1.7)
seurat.object <- RunUMAP(seurat.object, dims=1:30)
	## with res 1.5 I got 25 clusters and cluster 21 was macrophages with some potentially mis-assigned spatial spots.
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident")
p1
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE")
p2
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type")
p3
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Subtype")
p4
umap.combined <- CombinePlots(plots = list(p1, p2, p3, p4))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("integrated_UMAP-combined.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("integrated_UMAP-combined-wide.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")

dir.create("seurat_clusters.Type")
setwd("seurat_clusters.Type")

# "Trophoblast", "Fibroblast", "Stromal", "Smooth_muscle", "Endothelial", "Endometrial"
# c("Endometrial", "Endothelial", "Fibroblast", "Stromal", "Villious_cytotrophoblast", "Syncytiotrophoblast", "NK", "Macrophage", "Erythroblast")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("0"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "0.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("1"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "1.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("2"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "2.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("3"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "3.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("4"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "4.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("5"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "5.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("6"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "6.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("7"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "7.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("8"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "8.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("9"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "9.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("10"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "10.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("11"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "11.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("12"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "12.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("13"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "13.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("14"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "14.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("15"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "15.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("16"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "16.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("17"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "17.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("18"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "18.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("19"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "19.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")


Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("20"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "20.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")


Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("21"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "21.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("22"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "22.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("23"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "23.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("24"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "24.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("25"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "25.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("26"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "26.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")



setwd("..")
# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("seurat_clusters.Type/"))
manifest<-manifest%>%mutate(name=gsub(".Type.cellcounts.txt","",value))
manifest
df<-read.table("seurat_clusters.Type/1.Type.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type","1")
 df2<-df
df<-df%>%pivot_longer(cols = `1`,names_to = "seurat_clusters",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("seurat_clusters.Type/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Type",X)
  df<-df%>%pivot_longer(cols = X,names_to = "seurat_clusters",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of Types ## there are 12 seurat_clusterss
# df2%>%mutate(perc=percs_by_group(count,group = seurat_clusters))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])
df7<-my_fxn(X = manifest$value[6])
df11<-my_fxn(X = manifest$value[10])
df12<-my_fxn(X = manifest$value[11])
df13<-my_fxn(X = manifest$value[12])
df14<-my_fxn(X = manifest$value[13])
df15<-my_fxn(X = manifest$value[14])
df16<-my_fxn(X = manifest$value[15])
df17<-my_fxn(X = manifest$value[16])
df18<-my_fxn(X = manifest$value[17])
df19<-my_fxn(X = manifest$value[18])
df20<-my_fxn(X = manifest$value[19])
df21<-my_fxn(X = manifest$value[20])
df22<-my_fxn(X = manifest$value[21])
df23<-my_fxn(X = manifest$value[22])
df24<-my_fxn(X = manifest$value[23])
df25<-my_fxn(X = manifest$value[24])
df26<-my_fxn(X = manifest$value[25])
df27<-my_fxn(X = manifest$value[26])
final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)
final<-full_join(final,df7)
final<-full_join(final,df8)
final<-full_join(final,df9)
final<-full_join(final,df10)
final<-full_join(final,df11)
final<-full_join(final,df12)
final<-full_join(final,df13)
final<-full_join(final,df14)
final<-full_join(final,df15)
final<-full_join(final,df16)
final<-full_join(final,df17)
final<-full_join(final,df18)
final<-full_join(final,df19)
final<-full_join(final,df20)
final<-full_join(final,df21)
final<-full_join(final,df22)
final<-full_join(final,df23)
final<-full_join(final,df24)
final<-full_join(final,df25)
final<-full_join(final,df26)
final<-full_join(final,df27)

final<-final%>%
  mutate(seurat_clusters=gsub(".Type.cellcounts.txt","",seurat_clusters))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = seurat_clusters))
#  pivot_wider(names_from = seurat_clusters,values_from = count)

write.table(final,"seurat_clusters.CellCounts_merged.tsv",sep = "\t",row.names = F)

final
​
​
# df3<-sapply(manifest$value, my_fxn)
# 
# manifest<-manifest%>%
#   mutate(results=lapply(FUN = my_fxn(manifest$value)))
# 
#                             
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18

final$Type<-as.factor(final$Type)

final$seurat_clusters<-as.factor(final$seurat_clusters)
​
final
​
library(ggpubr)
ggbarplot(data = final,
          x = "seurat_clusters",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)

png(file=paste0("seurat_clusters.Type_diversity.plot.png"),
                res=300, 
                width=2500, 
                height=1500)
feature.plot <- ggbarplot(data = final,
          x = "seurat_clusters",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()

####################################################################################################


