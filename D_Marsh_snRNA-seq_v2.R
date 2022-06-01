# Marsh_snRNA-seq_v2.R 		using their published seurat object (no metadata labels for time)
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
setwd("/home/ebarrozo/Marsh/results")
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
library(SeuratDisk)

dir.create("seurat-v1")
setwd("seurat-v1")

setwd("/home/ebarrozo/Marsh/results/seurat-v1")

## load the data, make all genes CAPS-case, and make a seurat object 
seurat.object.data <-read.delim("/home/ebarrozo/Marsh/data/AllStages_AllNuclei_datamatrix.txt")
row.names(seurat.object.data) <- toupper(row.names(seurat.object.data))
seurat.object <- CreateSeuratObject(seurat.object.data, min.cells = 0, min.genes = 200, project = "marsh.published")
rm(seurat.object.data)

ncol(seurat.object)
	# 27326 unfiltered cells
Idents(seurat.object)

## Change the names of each library, 
seurat.object$orig.ident <- "Uninfected"   

## Assign and examine unfiltered quality control (QC) metrics
seurat.object <- PercentageFeatureSet(seurat.object, pattern = "^MT-", col.name = "percent.mt")
seurat.object <- PercentageFeatureSet(seurat.object, pattern = "^RP[SL]", col.name = "percent.ribo")

qc.vlnplot <- VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("seurat.object_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

        ## most cells <20% mt, nFeature <2000, ncount <10000, percent.ribo <40
## Filter based on the unfiltered QC plot above 
### eb-stringent filters based on the qc for this library, an approach supported by the literature
seurat.object <- subset(seurat.object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 10000 & percent.mt > 0.0000001 & percent.mt < 10 & percent.ribo < 20)
ncol(seurat.object)
#4333 filtered uninf cells
qc.vlnplot <- VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 5, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("seurat.object_filtered_QC_vlnplot.png", plot = qc.vlnplot, device = "png", width = 16, height = 8, units = "in")

# Log normalize data and scale by 10000, required for cell cycle scoring
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)

# Run cell cycle scoring
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
## human genes converted to mouse cc genes https://github.com/satijalab/seurat/issues/462
mouse.s.genes <- c("Mcm4",  "Exo1",  "Slbp",  "Gmnn", "Cdc45", "Msh2",  "Mcm6",  "Rrm2",  "Pold3", "Blm", "Ubr7",  "Mcm5",  "Clspn", "Hells", "Nasp",  "Rpa2",  "Rad51ap1", "Tyms",  "Rrm1",  "Rfc2", "Prim1", "Brip1", "Usp1",  "Ung", "Pola1", "Mcm2",  "Fen1",  "Tipin", "Pcna", "Cdca7", "Uhrf1", "Casp8ap2", "Cdc6",  "Dscc1", "Wdr76", "E2f8",  "Dtl", "Ccne2", "Atad2", "Gins2","Chaf1b","Pcna-ps2")
mouse.s.genes <- toupper(mouse.s.genes)
mouse.g2m.genes <- c("Nuf2", "Psrc1", "Ncapd2",  "Ccnb2", "Smc4", "Lbr",  "Tacc3", "Cenpa", "Kif23", "Cdca2", "Anp32e", "G2e3", "Cdca3", "Anln", "Cenpe", "Gas2l3",  "Tubb4b",  "Cenpf", "Dlgap5",  "Hjurp", "Cks1brt", "Gtse1", "Bub1", "Birc5", "Ube2c", "Rangap1", "Hmmr", "Ect2", "Tpx2", "Ckap5", "Cbx5", "Nek2", "Ttk", "Cdca8", "Nusap1", "Ctcf", "Cdc20", "Cks2", "Mki67", "Tmpo", "Ckap2l", "Aurkb", "Kif2c", "Cdk1", "Kif20b", "Top2a", "Aurka", "Ckap2", "Hmgb2", "Cdc25c", "Ndc80", "Kif11")
mouse.g2m.genes <- toupper(mouse.g2m.genes)
seurat.object <- CellCycleScoring(seurat.object, s.features = mouse.s.genes, g2m.features = mouse.g2m.genes, set.ident = FALSE)
seurat.object$CC.Difference <- seurat.object$S.Score - seurat.object$G2M.Score

all.genes <- rownames(seurat.object@assays[["RNA"]])
write.table(all.genes, "all.genes.txt", sep="\t")
	## 21513 genes

## Run SCTransform to normalize by neg. binomial
seurat.object <- SCTransform(seurat.object, assay = "RNA", verbose = T)

## Dimensionality reduction, clustering, and visualization
seurat.object <- RunPCA(seurat.object, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(seurat.object)
ggsave("seurat.object_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

seurat.object.var <- seurat.object@assays[["SCT"]]@var.features

# Cluster with the top 3k SCT variable genes
seurat.object <- RunPCA(seurat.object, features = seurat.object.var, assay = "SCT", slot = "scale.data")
seurat.object <- FindNeighbors(seurat.object, features = "seurat.object.var", dims = 1:30)
seurat.object <- FindClusters(seurat.object, resolution = 0.6)
seurat.object <- RunUMAP(seurat.object, dims=1:30)
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident")
p1
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE")
p2
	#clusters 0 thru 21
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Phase")

umap.combined <- CombinePlots(plots = list(p1, p2, p3))
ggsave("SCT_UMAP.png", plot = umap.combined, device = "png", width = 17, height = 12, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "orig.ident", label=TRUE)
ggsave("SCT_UMAP_splitby_libID.png", plot = p2, device = "png", width = 20, height = 3, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "seurat_clusters", label=TRUE)
ggsave("SCT_UMAP_splitby_clusters.png", plot = p2, device = "png", width = 20, height = 3, units = "in")
#examine UMAPs with qc metrics 
p5 <- FeaturePlot(seurat.object, features = 'nFeature_RNA')
p6 <- FeaturePlot(seurat.object, features = 'nCount_RNA')
p7 <- FeaturePlot(seurat.object, features = 'percent.mt')
p8 <- FeaturePlot(seurat.object, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("UMAP_QCmetricsFeaturePlots.png", plot = umap.combined, device = "png", width = 20, height = 12, units = "in")

##### Perform differential expression between clusters
DefaultAssay(seurat.object) <- "RNA"
Idents(seurat.object) <- "seurat_clusters"

seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = seurat.object.var, assay="RNA")

Idents(seurat.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "RNA", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "Marsh_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>0.693)

top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("seurat.object_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("seurat.object_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
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


############################################################# User Annotation of UMAP Clusters #############################################################
############################################################# ############################################################# #############################################################
dir.create("annotated")
setwd("annotated")
## For cluster Subtypes, I took the [marsh_DEGs_byclusters_pos-0.693lnFC.txt] 
		# cluster-Subtypes_DEGs_byclusters_top5k-pos-0.693lnFC.xlsx
	# with adj. p<0.05 and log2FC>2 and upload top marker genes into https://placentacellenrich.gdcb.iastate.edu and annotate based on vento-tormo et al or suryawanshi datasets. 
	## Also examining the human protein atlas (https://www.proteinatlas.org), Vento-Tormo and Suryawanshi datasets (https://placentacellenrich.gdcb.iastate.edu), and the PangloaDB (https://panglaodb.se/search.html)
Idents(seurat.object) <- "seurat_clusters"
cluster.Subtypes<- c("NK_1", "Stromal", "Syncytiotrophoblast_1", "Endometrial", "Fibroblast_1", "Villious_cytotrophoblast", "Endothelial_1", "Endothelial_2", "Endothelial_3", "Erythroblast", "Fibroblast_2", "Syncytiotrophoblast_2", "Endothelial_4", "Fibroblast_3", "Fibroblast_4", "Fibroblast_5", "Macrophage", "Fibroblast_6", "Endothelial_5", "Endothelial_6", "Fibroblast_7", "NK_2")
names(cluster.Subtypes) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Subtypes)
## add cluster Subtypes as a metadata column
seurat.object$Subtype <- Idents(seurat.object)

Idents(seurat.object) <- "seurat_clusters"
cluster.Type <- c("NK", "Stromal", "Syncytiotrophoblast", "Endometrial", "Fibroblast", "Villious_cytotrophoblast", "Endothelial", "Endothelial", "Endothelial", "Erythroblast", "Fibroblast", "Syncytiotrophoblast", "Endothelial", "Fibroblast", "Fibroblast", "Fibroblast", "Macrophage", "Fibroblast", "Endothelial", "Endothelial", "Fibroblast", "NK")
names(cluster.Type) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Type)
## add cluster Types as a metadata column
seurat.object$Type <- Idents(seurat.object)

## Re-order the clusters to group similar clusters together
Idents(seurat.object) <- "Subtype"
levels(seurat.object)
seurat.object2 <- subset(x = seurat.object, idents = c("Endometrial", "Endothelial_1", "Endothelial_2", "Endothelial_3", "Endothelial_4", "Endothelial_5", "Endothelial_6", "Fibroblast_1", "Fibroblast_2", "Fibroblast_3", "Fibroblast_4", "Fibroblast_5", "Fibroblast_6", "Fibroblast_7", "Stromal", "Villious_cytotrophoblast", "Syncytiotrophoblast_1", "Syncytiotrophoblast_2", "NK_1", "NK_2", "Macrophage", "Erythroblast"))
ncol(seurat.object2)
		# make sure you didn't lose any cells
        # 5368 cells
levels(seurat.object2)
levels(x=seurat.object2) <- c("Endometrial", "Endothelial_1", "Endothelial_2", "Endothelial_3", "Endothelial_4", "Endothelial_5", "Endothelial_6", "Fibroblast_1", "Fibroblast_2", "Fibroblast_3", "Fibroblast_4", "Fibroblast_5", "Fibroblast_6", "Fibroblast_7", "Stromal", "Villious_cytotrophoblast", "Syncytiotrophoblast_1", "Syncytiotrophoblast_2", "NK_1", "NK_2", "Macrophage", "Erythroblast")
levels(seurat.object2)
ncol(seurat.object)
	# 5368
seurat.object <- seurat.object2
rm(seurat.object2)
levels(seurat.object)

## Re-order the clusters to group similar clusters together
Idents(seurat.object) <- "Type"
levels(seurat.object)
seurat.object2 <- subset(x = seurat.object, idents = c("Endometrial", "Endothelial", "Fibroblast", "Stromal", "Villious_cytotrophoblast", "Syncytiotrophoblast", "NK", "Macrophage", "Erythroblast"))
ncol(seurat.object2)
        # make sure you didn't lose any cells
        # 3621 cells
levels(seurat.object2)
levels(x=seurat.object2) <- c("Endometrial", "Endothelial", "Fibroblast", "Stromal", "Villious_cytotrophoblast", "Syncytiotrophoblast", "NK", "Macrophage", "Erythroblast")
levels(seurat.object2)
ncol(seurat.object)
    # 3621
seurat.object <- seurat.object2
rm(seurat.object2)
levels(seurat.object)

Idents(seurat.object) <- "seurat_clusters"
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE")
p1
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", label= "TRUE")
p2
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Subtype")
p3
umap.combined <- CombinePlots(plots = list(p1, p2, p3))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("marsh_UMAP-annotated.pdf", plot = umap.combined, device = "pdf", width = 12, height = 7, units = "in")

## make essential figs with the cluster Subtype
Idents(seurat.object) <- "Subtype"
p10 <- DimPlot(seurat.object, reduction = "umap", group.by = "Subtype", label = TRUE, repel=TRUE)
ggsave("Annotated-marsh_UMAP-labelled.pdf", plot = p10, device = "pdf", width = 7, height = 5, units = "in")
p11 <- DimPlot(seurat.object, reduction = "umap", group.by = "Subtype")
ggsave("Annotated-marsh_UMAP-non-labelled.pdf", plot = p11, device = "pdf", width = 7, height = 5, units = "in")
ggsave("Annotated-marsh_UMAP-non-labelled_xl.pdf", plot = p11, device = "pdf", width = 10, height = 5, units = "in")

p14 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type")
ggsave("Type-marsh_UMAP-non-labelled.pdf", plot = p14, device = "pdf", width = 7, height = 5, units = "in")
ggsave("Type-marsh_UMAP-non-labelled_xl.pdf", plot = p14, device = "pdf", width = 12, height = 5, units = "in")
p12 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", label = TRUE, repel=TRUE)
ggsave("Type-marsh_UMAP-labelled.pdf", plot = p12, device = "pdf", width = 7, height = 5, units = "in")
ggsave("Type-marsh_UMAP-labelled_xl.pdf", plot = p12, device = "pdf", width = 10, height = 5, units = "in")


##### Perform differential expression between Types
DefaultAssay(seurat.object) <- "RNA"
Idents(seurat.object) <- "Type"

seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = seurat.object.var, assay="RNA")

Idents(seurat.object) <- "Type"
de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "RNA", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "Marsh_DEGs_byType_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>0.693)

top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("seurat.object_top5_markers_logfc_heatmaptop3k_type.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("seurat.object_top5_markers_counts_heatmaptop3k_type.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top2.marker-DotPlot_type.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot_type.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top3.marker-DotPlot_type.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top5.marker-DotPlot_type.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")

##### Perform differential expression between Subtypes
DefaultAssay(seurat.object) <- "RNA"
Idents(seurat.object) <- "Subtype"

seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = seurat.object.var, assay="RNA")

Idents(seurat.object) <- "Subtype"
de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "RNA", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "Marsh_DEGs_bySubtype_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>0.693)

top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("seurat.object_top5_markers_logfc_heatmaptop3k_subtype.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("seurat.object_top5_markers_counts_heatmaptop3k_subtype.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top2.marker-DotPlot_subtype.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot_subtype.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top3.marker-DotPlot_subtype.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top5.marker-DotPlot_subtype.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()


Idents(seurat.object) <- "Type"
seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "RNA")
## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell_Type.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(seurat.object, features = lymphoid.lineage)
png(file=paste0("lymphoid.lineage-DotPlot_Type.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Lymphoid Lineage Markers") + scale_y_discrete(name ="Type")
dev.off()

## canonical markers for myeloid cells https://docs.abcam.com/pdf/immunology/myeloid_cell_Type.pdf
# myeloid.lineage <- c("CD34", "CD117", "CD271", "CD33", "CD71", "CD61", "CD123", "CD44", "CD15", "CD16", "CD11b", "CD14", "CD1c", "CD141")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
myeloid.lineage <- c("CD34", "KIT", "NGFR", "CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
feature.plot <- DotPlot(seurat.object, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot_Type.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="Type")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(seurat.object, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot_Type.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="Type")

Idents(seurat.object) <- "Subtype"
#seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "RNA")
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(seurat.object, features = lymphoid.lineage)
png(file=paste0("lymphoid.lineage-DotPlot_Subtype.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Lymphoid Lineage Markers") + scale_y_discrete(name ="Subtype")
dev.off()

myeloid.lineage <- c("CD34", "KIT", "NGFR", "CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
feature.plot <- DotPlot(seurat.object, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot_Subtype.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="Subtype")
dev.off()

immunophenoSubtype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenoSubtype.lineage <- unique(immunophenoSubtype.lineage)
feature.plot <- DotPlot(seurat.object, features = immunophenoSubtype.lineage)
png(file=paste0("immunophenoSubtype.lineage-DotPlot_Subtype.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="ImmunephenoSubtype Markers") + scale_y_discrete(name ="Subtype")
dev.off()





all.genes.marsh <- all.genes
marsh.var <- seurat.object.var

rm(all.genes)
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
rm(p11)
rm(p12)
rm(p11)
rm(p12)
rm(p10)
rm(p14)
rm(p5)
rm(p6)
rm(p7)
rm(p8)
rm(qc.vlnplot)
rm(umap.combined)
rm(cluster.Subtypes)
rm(cluster.Type)

## Slim down the seurat.objects to save space. Don't need to keep the scale.data
marsh.seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
rm(seurat.object)

setwd("/home/ebarrozo/Marsh/results/seurat-v1")
save.image("marsh_annotated_v1.RData")
#load("/home/ebarrozo/Marsh/results/seurat-v1/marsh_annotated_v1.RData")



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
# df<-read.table("Type.cluster/AF_Epithelial.cluster.cellcounts.txt",header = T,sep = "\t",row.names = 1)
# colnames(df)<-c("cluster","AF_Epithelial")
 df2<-df
# df<-df%>%pivot_longer(cols = `AF_Epithelial`,names_to = "Type",values_to = "Count")
# df
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

