## visium_mouse_seurat_v1.R

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

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
dir.create("seurat_mouse_v1")
setwd("seurat_mouse_v1")

setwd("/home/ebarrozo/visium/results/seurat_mouse_v1")

############################################################# ############################################################# #############################################################
############################################################# S07: Load data, Quality Control (QC) filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S07-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S07.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S07.manual)
new.names <- gsub ("MM10-----------(.*)", "\\1", old.names)
new.names <- gsub ("NC-012532.1.V2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S07.manual <- RenameGenesSeurat(obj = S07.manual, newnames = new.names)
row.names(S07.manual)

## Save a table of all.genes
all.genes <- row.names(S07.manual)
write.table(all.genes, "all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
ZIKV.genes <- c("5PUTR", "ANCHC", "M", "E", "NS1", "NS2A", "NS2B", "NS3", "NS4A", "2K", "NS4B", "NS5", "3PUTR")

## Assign and examine unfiltered quality control (QC) metrics
S07.manual <- PercentageFeatureSet(S07.manual, pattern = "^MT-", col.name = "percent.mt")
S07.manual <- PercentageFeatureSet(S07.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S07.manual <- PercentageFeatureSet(S07.manual, features = ZIKV.genes, col.name = "percent.viral")

## Log normalize data and multiply by a factor of 10000
	## required for cellcyclescoring
S07.manual<- NormalizeData(S07.manual, normalization.method = "LogNormalize", scale.factor = 10000)

## human genes converted to mouse cc genes https://github.com/satijalab/seurat/issues/462
mouse.s.genes <- c("Mcm4",  "Exo1",  "Slbp",  "Gmnn", "Cdc45", "Msh2",  "Mcm6",  "Rrm2",  "Pold3", "Blm", "Ubr7",  "Mcm5",  "Clspn", "Hells", "Nasp",  "Rpa2",  "Rad51ap1", "Tyms",  "Rrm1",  "Rfc2", "Prim1", "Brip1", "Usp1",  "Ung", "Pola1", "Mcm2",  "Fen1",  "Tipin", "Pcna", "Cdca7", "Uhrf1", "Casp8ap2", "Cdc6",  "Dscc1", "Wdr76", "E2f8",  "Dtl", "Ccne2", "Atad2", "Gins2","Chaf1b","Pcna-ps2")
mouse.s.genes <- toupper(mouse.s.genes)
mouse.g2m.genes <- c("Nuf2", "Psrc1", "Ncapd2",  "Ccnb2", "Smc4", "Lbr",  "Tacc3", "Cenpa", "Kif23", "Cdca2", "Anp32e", "G2e3", "Cdca3", "Anln", "Cenpe", "Gas2l3",  "Tubb4b",  "Cenpf", "Dlgap5",  "Hjurp", "Cks1brt", "Gtse1", "Bub1", "Birc5", "Ube2c", "Rangap1", "Hmmr", "Ect2", "Tpx2", "Ckap5", "Cbx5", "Nek2", "Ttk", "Cdca8", "Nusap1", "Ctcf", "Cdc20", "Cks2", "Mki67", "Tmpo", "Ckap2l", "Aurkb", "Kif2c", "Cdk1", "Kif20b", "Top2a", "Aurka", "Ckap2", "Hmgb2", "Cdc25c", "Ndc80", "Kif11")
mouse.g2m.genes <- toupper(mouse.s.genes)

# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S07.manual <- CellCycleScoring(S07.manual, s.features = mouse.s.genes, g2m.features = mouse.g2m.genes, set.ident = FALSE)
S07.manual$CC.Difference <- S07.manual$S.Score - S07.manual$G2M.Score

acute.viral <- subset(S07.manual, features = rownames(S07.manual)[rownames(S07.manual) %in% ZIKV.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_all-cells.txt", sep="\t")

	#### Confirmed manually there are no ZIKV transcripts ###########  =SUM(B32287:PN32287)

acute.viral <- subset(acute.viral, percent.viral > 0)
acute.viral <- subset(acute.viral, features = rownames(S07.manual)[rownames(S07.manual) %in% ZIKV.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

ncol(S07.manual)
	# 429 spots

S07.manual$orig.ident <- "ZIKV-sham"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S07.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S07.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
		# pt.size.factor = 1.6 default
		# alpha - minimum and maximum transparency. Default is c(1, 1).
		# Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
		## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S07.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S07.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S07.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S07.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S07.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S07.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S07.manual)
	#  429
S07.manual2 <- subset(S07.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 10 & percent.ribo < 10)
ncol(S07.manual2)
	# 	414
qc.vlnplot <- VlnPlot(S07.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S07.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S07.manual <- S07.manual2
rm(S07.manual2)
qc.vlnplot <- VlnPlot(S07.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S07.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
############################################################# ############################################################# #############################################################
############################################################# S09: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################
## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S09-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S09.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = TRUE)
## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S09.manual)
new.names <- gsub ("MM10-----------(.*)", "\\1", old.names)
new.names <- gsub ("NC-012532.1.V2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

S09.manual <- RenameGenesSeurat(obj = S09.manual, newnames = new.names)
row.names(S09.manual)

## Save a table of all.genes
all.genes <- row.names(S09.manual)
write.table(all.genes, "S09.manual.all.genes.txt", sep="\t")

## Assign and examine unfiltered quality control (QC) metrics
S09.manual <- PercentageFeatureSet(S09.manual, pattern = "^MT-", col.name = "percent.mt")
S09.manual <- PercentageFeatureSet(S09.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S09.manual <- PercentageFeatureSet(S09.manual, features = ZIKV.genes, col.name = "percent.viral")

## Log normalize data and multiply by a factor of 10000
	## required for cellcyclescoring
S09.manual<- NormalizeData(S09.manual, normalization.method = "LogNormalize", scale.factor = 10000)

# Run cell cycle scoring 
S09.manual <- CellCycleScoring(S09.manual, s.features = mouse.s.genes, g2m.features = mouse.g2m.genes, set.ident = FALSE)
S09.manual$CC.Difference <- S09.manual$S.Score - S09.manual$G2M.Score

acute.viral <- subset(S09.manual, features = rownames(S09.manual)[rownames(S09.manual) %in% ZIKV.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S09.manual_viral.counts_all-cells.txt", sep="\t")
acute.viral <- subset(acute.viral, percent.viral > 0)
	# Error: No cells found
#acute.viral <- subset(acute.viral, features = rownames(S09.manual)[rownames(S09.manual) %in% ZIKV.genes])
#write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

ncol(S09.manual)
	# 1625 spots

S09.manual$orig.ident <- "Mock-enox-1"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S09.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S09.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S09.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S09.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S09.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S09.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S09.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S09.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S09.manual)
	#  1625
S09.manual2 <- subset(S09.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 10 & percent.ribo < 10)
ncol(S09.manual2)
	# 	1534
qc.vlnplot <- VlnPlot(S09.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S09.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S09.manual <- S09.manual2
rm(S09.manual2)
qc.vlnplot <- VlnPlot(S09.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S09.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
############################################################# ############################################################# #############################################################
############################################################# S12: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################
## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S12-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S12.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = TRUE)
## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S12.manual)
new.names <- gsub ("MM10-----------(.*)", "\\1", old.names)
new.names <- gsub ("NC-012532.1.V2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

S12.manual <- RenameGenesSeurat(obj = S12.manual, newnames = new.names)
row.names(S12.manual)

## Save a table of all.genes
all.genes <- row.names(S12.manual)
write.table(all.genes, "S12.manual.all.genes.txt", sep="\t")

## Assign and examine unfiltered quality control (QC) metrics
S12.manual <- PercentageFeatureSet(S12.manual, pattern = "^MT-", col.name = "percent.mt")
S12.manual <- PercentageFeatureSet(S12.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S12.manual <- PercentageFeatureSet(S12.manual, features = ZIKV.genes, col.name = "percent.viral")

## Log normalize data and multiply by a factor of 10000
	## required for cellcyclescoring
S12.manual<- NormalizeData(S12.manual, normalization.method = "LogNormalize", scale.factor = 10000)

# Run cell cycle scoring 
S12.manual <- CellCycleScoring(S12.manual, s.features = mouse.s.genes, g2m.features = mouse.g2m.genes, set.ident = FALSE)
S12.manual$CC.Difference <- S12.manual$S.Score - S12.manual$G2M.Score

acute.viral <- subset(S12.manual, features = rownames(S12.manual)[rownames(S12.manual) %in% ZIKV.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S12.manual_viral.counts_all-cells.txt", sep="\t")
acute.viral <- subset(acute.viral, percent.viral > 0)
	# Error: No cells found
#acute.viral <- subset(acute.viral, features = rownames(S12.manual)[rownames(S12.manual) %in% ZIKV.genes])
#write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

ncol(S12.manual)
	# 2302 spots

S12.manual$orig.ident <- "Mock-enox-2"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S12.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S12.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S12.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S12.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S12.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S12.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S12.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S12.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S12.manual)
	#  1625
S12.manual2 <- subset(S12.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 10 & percent.ribo < 10)
ncol(S12.manual2)
	# 	2290
qc.vlnplot <- VlnPlot(S12.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S12.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S12.manual <- S12.manual2
rm(S12.manual2)
qc.vlnplot <- VlnPlot(S12.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S12.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
############################################################# ############################################################# #############################################################
############################################################# S13: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S13-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S13.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = TRUE)
## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S13.manual)
new.names <- gsub ("MM10-----------(.*)", "\\1", old.names)
new.names <- gsub ("NC-013532.1.V2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

S13.manual <- RenameGenesSeurat(obj = S13.manual, newnames = new.names)
row.names(S13.manual)

## Save a table of all.genes
all.genes <- row.names(S13.manual)
write.table(all.genes, "S13.manual.all.genes.txt", sep="\t")

## Assign and examine unfiltered quality control (QC) metrics
S13.manual <- PercentageFeatureSet(S13.manual, pattern = "^MT-", col.name = "percent.mt")
S13.manual <- PercentageFeatureSet(S13.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S13.manual <- PercentageFeatureSet(S13.manual, features = ZIKV.genes, col.name = "percent.viral")

## Log normalize data and multiply by a factor of 10000
	## required for cellcyclescoring
S13.manual<- NormalizeData(S13.manual, normalization.method = "LogNormalize", scale.factor = 10000)

# Run cell cycle scoring 
S13.manual <- CellCycleScoring(S13.manual, s.features = mouse.s.genes, g2m.features = mouse.g2m.genes, set.ident = FALSE)
S13.manual$CC.Difference <- S13.manual$S.Score - S13.manual$G2M.Score

acute.viral <- subset(S13.manual, features = rownames(S13.manual)[rownames(S13.manual) %in% ZIKV.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S13.manual_viral.counts_all-cells.txt", sep="\t")
acute.viral <- subset(acute.viral, percent.viral > 0)
	# Error: No cells found
#acute.viral <- subset(acute.viral, features = rownames(S13.manual)[rownames(S13.manual) %in% ZIKV.genes])
#write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

ncol(S13.manual)
	# 1182 spots

S13.manual$orig.ident <- "ZIKV-enox"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S13.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S13.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S13.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S13.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S13.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S13.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S13.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S13.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S13.manual)
	#  1625
S13.manual2 <- subset(S13.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 10 & percent.ribo < 10)
ncol(S13.manual2)
	# 	1130
qc.vlnplot <- VlnPlot(S13.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S13.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S13.manual <- S13.manual2
rm(S13.manual2)
qc.vlnplot <- VlnPlot(S13.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S13.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

rm(plot1)
rm(plot2)
rm(acute.viral)
rm(plot3)
rm(qc.vlnplot)
rm(new.names)
rm(old.names)
rm(data_dir)
rm(to.upper)

############################################################# ############################################################# #############################################################
############################################################# S07: Dimension Reduction, Visualization, and Differential Expression (DE) #############################################################
############################################################# ############################################################# #############################################################

## Run SCTransform to normalize by neg. binomial
S07.manual <- SCTransform(S07.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S07.manual <- RunPCA(S07.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S07.manual)
ggsave("S07.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S07.manual <- FindNeighbors(S07.manual, reduction = "pca", dims = 1:30)
S07.manual <- FindClusters(S07.manual, verbose = T)
S07.manual <- RunUMAP(S07.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S07.manual, reduction = "umap", label = TRUE)
p1
	# clusters 08
p2 <- SpatialDimPlot(S07.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S07.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
	# As there are many colors, it can be challenging to visualize which voxel belongs to which cluster. We have a few strategies to help with this. Setting the label parameter places a colored box at the median of each cluster (see the plot above).
p3 <- wrap_plots(p1, p2, p4)
ggsave("S07.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S07.manual, cells.highlight = CellsByIdentities(object = S07.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8)), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0)
ggsave("S07.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S07.manual <- FindSpatiallyVariableFeatures(S07.manual, assay = "SCT", features = VariableFeatures(S07.manual)[1:1000], selection.method = "markvariogram")
S07.manual.var <- S07.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S07.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S07.manual, features = top.features, ncol = 3, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S07.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S07.manual <- NormalizeData(S07.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S07.manual <- ScaleData(S07.manual, features = S07.manual.var, assay="Spatial")

Idents(S07.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S07.manual, features = S07.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S07_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>0.693)
top5.heatmap <- DoHeatmap(S07.manual, features = top5$gene, raster = FALSE)
ggsave("S07.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S07.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S07.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S07.manual, features = unique.top2)
png(file=paste0("S07.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S07.manual, features = unique.top2)
png(file=paste0("S07.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S07.manual, features = unique.top2)
png(file=paste0("S07.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
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
############################################################# ############################################################# #############################################################
############################################################# S09: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S09.manual <- SCTransform(S09.manual, assay = "Spatial", verbose = FALSE)

## Dimensionality reduction, clustering, and visualization
S09.manual <- RunPCA(S09.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S09.manual)
ggsave("S09.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S09.manual <- FindNeighbors(S09.manual, reduction = "pca", dims = 1:30)
S09.manual <- FindClusters(S09.manual, verbose = T)
S09.manual <- RunUMAP(S09.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S09.manual, reduction = "umap", label = TRUE)
p1
	# clusters 0-9
p2 <- SpatialDimPlot(S09.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S09.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S09.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
	## Change idents to match the number of clusters
p3<- SpatialDimPlot(S09.manual, cells.highlight = CellsByIdentities(object = S09.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8,9)), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0)
ggsave("S09.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S09.manual <- FindSpatiallyVariableFeatures(S09.manual, assay = "SCT", features = VariableFeatures(S09.manual)[1:1000], selection.method = "markvariogram")
S09.manual.var <- S09.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S09.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S09.manual, features = top.features, ncol = 3, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S09.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S09.manual <- NormalizeData(S09.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S09.manual <- ScaleData(S09.manual, features = S09.manual.var, assay="Spatial")

Idents(S09.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S09.manual, features = S09.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S09_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>0.693)
top5.heatmap <- DoHeatmap(S09.manual, features = top5$gene, raster = FALSE)
ggsave("S09.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S09.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S09.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S09.manual, features = unique.top2)
png(file=paste0("S09.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S09.manual, features = unique.top2)
png(file=paste0("S09.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S09.manual, features = unique.top2)
png(file=paste0("S09.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
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
############################################################# ############################################################# #############################################################
############################################################# S12: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S12.manual <- SCTransform(S12.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S12.manual <- RunPCA(S12.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S12.manual)
ggsave("S12.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S12.manual <- FindNeighbors(S12.manual, reduction = "pca", dims = 1:30)
S12.manual <- FindClusters(S12.manual, verbose = T)
S12.manual <- RunUMAP(S12.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S12.manual, reduction = "umap", label = TRUE)
p1
	# clusters 0-12
p2 <- SpatialDimPlot(S12.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S12.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S12.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S12.manual, cells.highlight = CellsByIdentities(object = S12.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0)
ggsave("S12.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S12.manual <- FindSpatiallyVariableFeatures(S12.manual, assay = "SCT", features = VariableFeatures(S12.manual)[1:1000], selection.method = "markvariogram")
S12.manual.var <- S12.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S12.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S12.manual, features = top.features, ncol = 3, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S12.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S12.manual <- NormalizeData(S12.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S12.manual <- ScaleData(S12.manual, features = S12.manual.var, assay="Spatial")

Idents(S12.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S12.manual, features = S12.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S12_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>0.693)
top5.heatmap <- DoHeatmap(S12.manual, features = top5$gene, raster = FALSE)
ggsave("S12.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S12.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S12.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S12.manual, features = unique.top2)
png(file=paste0("S12.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S12.manual, features = unique.top2)
png(file=paste0("S12.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S12.manual, features = unique.top2)
png(file=paste0("S12.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
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
############################################################# ############################################################# #############################################################
############################################################# S13: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S13.manual <- SCTransform(S13.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S13.manual <- RunPCA(S13.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S13.manual)
ggsave("S13.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S13.manual <- FindNeighbors(S13.manual, reduction = "pca", dims = 1:30)
S13.manual <- FindClusters(S13.manual, verbose = T)
S13.manual <- RunUMAP(S13.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S13.manual, reduction = "umap", label = TRUE)
p1
	##### clusters 0 thru 8
p2 <- SpatialDimPlot(S13.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S13.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S13.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S13.manual, cells.highlight = CellsByIdentities(object = S13.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8)), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0)
ggsave("S13.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S13.manual <- FindSpatiallyVariableFeatures(S13.manual, assay = "SCT", features = VariableFeatures(S13.manual)[1:1000], selection.method = "markvariogram")
S13.manual.var <- S13.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S13.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S13.manual, features = top.features, ncol = 3, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S13.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S13.manual <- NormalizeData(S13.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S13.manual <- ScaleData(S13.manual, features = S13.manual.var, assay="Spatial")

Idents(S13.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S13.manual, features = S13.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S13_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>0.693)
		## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S13.manual, features = S13.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S13.manual, features = top5$gene, raster = FALSE)
ggsave("S13.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S13.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S13.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S13.manual, features = unique.top2)
png(file=paste0("S13.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S13.manual, features = unique.top2)
png(file=paste0("S13.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S13.manual, features = unique.top2)
png(file=paste0("S13.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
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
############################################################# ############################################################# #############################################################
############################################################# Integration: merge all data using CCA anchors #############################################################
############################################################# ############################################################# #############################################################

## Slim down the seurat.objects to save space. Don't need to keep the scale.data
S07.manual2 <- DietSeurat(S07.manual, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
View(S07.manual2)
	## from 188MB to 166 MB
S07.manual <- S07.manual2
rm(S07.manual2)
S09.manual <- DietSeurat(S09.manual, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
S12.manual <- DietSeurat(S12.manual, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
S13.manual <- DietSeurat(S13.manual, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

## combine lists of top variable features for later DE analysis/clustering
var.combined <- union(S07.manual.var, S09.manual.var)
var.combined <- union(var.combined, S12.manual.var)
var.combined <- union(var.combined, S13.manual.var)
	## 6564 combined variable features

## Merge objects
all_merged <- merge(x = S07.manual, y = c(S09.manual, S12.manual, S13.manual), merge.data = TRUE, project = "all_merged")
## make combined QC plot 
qc.vlnplot <- VlnPlot(all_merged, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), group.by = "orig.ident", pt.size = 0.000001, ncol = 4) + NoLegend()
ggsave("all_merged_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 15, height = 5, units = "in")
rm(qc.vlnplot)
rm(S07.manual)
rm(S09.manual)
rm(S12.manual)
rm(S13.manual)
rm(S07.manual.var)
rm(S09.manual.var)
rm(S12.manual.var)
rm(S13.manual.var)

options(future.globals.maxSize = 300000000000)

all_merged <- SCTransform(all_merged, assay = "Spatial", verbose = T)
all_merged.list <- SplitObject(all_merged, split.by = "orig.ident")


# Make a reference list
reference.list <- all_merged.list[c("Mock-enox-1", "Mock-enox-2", "ZIKV-sham", "ZIKV-enox")]

# rm("all_merged")
# rm(all_merged.list)

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
rm("all_merged")
rm(all_merged.list)
rm("all_merged.anchors")

# Add all viral+human genes as a list of features to call upon during clustering or differential expression of all genes, which may be compuationally demanding
all.genes <- rownames(seurat.object@assays[["Spatial"]])

# Use the 'integrated' assay only for clustering
DefaultAssay(seurat.object) <- "integrated"
# Add all integrated genes as a list of features to call upon during clustering (or quick- differential expression)
integrated.genes <- rownames(seurat.object)
write.table(integrated.genes, "integrated.genes.txt", sep="\t")

## Use the 'Spatial' object for QC and differential expression analysis of allllll genes, which is computationally demanding
#DefaultAssay(seurat.object) <- "Spatial"

setwd("/home/ebarrozo/visium/results/seurat_mouse_v1")
save.image("murine_spatial_data-integrated_v1.RData")

############################################################# ############################################################# #############################################################
############################################################# Integrated: Dimension Reduction and Visualization #############################################################
############################################################# ############################################################# #############################################################

seurat.object <- ScaleData(seurat.object, features = integrated.genes, assay = "integrated", vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"))
seurat.object<- RunPCA(seurat.object, features = integrated.genes, assay = "integrated", slot = "scale.data")

seurat.object <- FindNeighbors(seurat.object, features = "integrated.genes", dims = 1:30)
seurat.object <- FindClusters(seurat.object, resolution = 0.6)
seurat.object <- RunUMAP(seurat.object, dims=1:30)
	## clusters 0 thru 12
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident")
p1
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE")
p2
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Phase")
p3
umap.combined <- CombinePlots(plots = list(p1, p2, p3))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("integrated_UMAP.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "orig.ident", label=TRUE)
ggsave("integrated_UMAP_splitby_libID.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "seurat_clusters", label=TRUE)
ggsave("integrated_UMAP_splitby_clusters.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
#examine UMAPs with qc metrics 
p5 <- FeaturePlot(seurat.object, features = 'nFeature_Spatial')
p6 <- FeaturePlot(seurat.object, features = 'nCount_Spatial')
p7 <- FeaturePlot(seurat.object, features = 'percent.mt')
p8 <- FeaturePlot(seurat.object, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")


###### UMAP + Spatial DimPlots
p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("slice1.5")
p2 <- SpatialDimPlot(seurat.object, images=image.list, pt.size.factor = 3)
p2
image.list <- c("slice1.1")
p3 <- SpatialDimPlot(seurat.object, images=image.list)
p3
image.list <- c("slice1.2.1")
p4 <- SpatialDimPlot(seurat.object, images=image.list)
p4
image.list <- c("slice1.3.3")
p5 <- SpatialDimPlot(seurat.object, images=image.list)
p5
p6 <- wrap_plots(p1, p2, p3, p4, p5, ncol=5)
ggsave("integrated_UMAP_SpatialDimPlots-cropped.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")

p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("slice1.5")
p2 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p2
image.list <- c("slice1.1")
p3 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p3
image.list <- c("slice1.2.1")
p4 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p4
image.list <- c("slice1.3.3")
p5 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p5
p6 <- wrap_plots(p1, p2, p3, p4, p5, ncol=5)
ggsave("integrated_UMAP_SpatialDimPlots-uncropped.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")


## Spatial DimPlots split.by clusters
image.list <- c("slice1.5")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S07.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.1")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S09.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.2.1")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S12.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.3.3")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S13.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")

rm(umap.combined)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(p5)
rm(p6)
rm(p7)
rm(p8)
rm(image.list)

############################################################# ############################################################# #############################################################
############################################################# DE Analysis #############################################################
############################################################# ############################################################# #############################################################

seurat.object <- SCTransform(seurat.object, assay = "Spatial", verbose = T)
DefaultAssay(seurat.object) <- "SCT"
seurat.object.var <- seurat.object@assays[["SCT"]]@var.features

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = seurat.object.var, assay="Spatial")

Idents(seurat.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "integrated_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>0.693)
		## all clusters have sig. genes
			# de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
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
p4 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p4
ggsave("top.markers_SpatialFeaturePlots-cropped-S12.pdf", plot = p4, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.3.3")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
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
p4 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p4
ggsave("top.markers_SpatialPlots-cropped-S12.pdf", plot = p4, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.3.3")
p5 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
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
############################################################# User Annotation of UMAP Clusters #############################################################
############################################################# ############################################################# #############################################################
dir.create("annotated")
setwd("annotated")
## For cluster Subtypes, I took the [integrated_DEGs_byclusters_pos-0.693lnFC.txt] 
		# cluster-Subtypes_DEGs_byclusters_top5k-pos-0.693lnFC.xlsx
	# with adj. p<0.05 and log2FC>2 and upload top marker genes into https://placentacellenrich.gdcb.iastate.edu and annotate based on vento-tormo et al or suryawanshi datasets. 
	## Also examining the human protein atlas (https://www.proteinatlas.org), Vento-Tormo and Suryawanshi datasets (https://placentacellenrich.gdcb.iastate.edu), and the PangloaDB (https://panglaodb.se/search.html)
Idents(seurat.object) <- "seurat_clusters"
cluster.Subtypes<- c("Trophoblast_1", "Stromal", "Trophoblast_2", "Fibroblast_1", "Endothelial_1", "Fibroblast_2", "Smooth_muscle_1", "Smooth_muscle_2", "Endothelial_2", "Trophoblast_3", "Endometrial_1", "Trophoblast_4", "Endometrial_2")
names(cluster.Subtypes) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Subtypes)
## add cluster Subtypes as a metadata column
seurat.object$Subtype <- Idents(seurat.object)

Idents(seurat.object) <- "seurat_clusters"
cluster.Type <- c("Trophoblast", "Stromal", "Trophoblast", "Fibroblast", "Endothelial", "Fibroblast", "Smooth_muscle", "Smooth_muscle", "Endothelial", "Trophoblast", "Endometrial", "Trophoblast", "Endometrial")
names(cluster.Type) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Type)
## add cluster Types as a metadata column
seurat.object$Type <- Idents(seurat.object)


cluster.Type <-c("Trophoblast", "Stromal", "Trophoblast", "Fibroblast", "Endothelial", "Fibroblast", "Smooth_muscle", "Smooth_muscle", "Endothelial", "Trophoblast", "Endometrial", "Trophoblast", "Endometrial")
names(cluster.Type) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Type)
## add cluster Types as a metadata column
seurat.object$Type <- Idents(seurat.object)

## Re-order the clusters to group similar clusters together
Idents(seurat.object) <- "Subtype"
levels(seurat.object)
seurat.object2 <- subset(x = seurat.object, idents = c("Trophoblast_1", "Trophoblast_2", "Trophoblast_3", "Trophoblast_4", "Fibroblast_1", "Fibroblast_2", "Stromal",  "Smooth_muscle_1", "Smooth_muscle_2", "Endothelial_1", "Endothelial_2", "Endometrial_1", "Endometrial_2"))
ncol(seurat.object2)
		# make sure you didn't lose any cells
        # 5368 cells
levels(seurat.object2)
levels(x=seurat.object2) <- c("Trophoblast_1", "Trophoblast_2", "Trophoblast_3", "Trophoblast_4", "Fibroblast_1", "Fibroblast_2", "Stromal",  "Smooth_muscle_1", "Smooth_muscle_2", "Endothelial_1", "Endothelial_2", "Endometrial_1", "Endometrial_2")
levels(seurat.object2)
ncol(seurat.object)
	# 5368
seurat.object <- seurat.object2
rm(seurat.object2)
levels(seurat.object)

## Re-order the clusters to group similar clusters together
Idents(seurat.object) <- "Type"
levels(seurat.object)
seurat.object2 <- subset(x = seurat.object, idents = c("Trophoblast", "Fibroblast", "Stromal", "Smooth_muscle", "Endothelial", "Endometrial"))
ncol(seurat.object2)
        # make sure you didn't lose any cells
        # 3621 cells
levels(seurat.object2)
levels(x=seurat.object2) <- c("Trophoblast", "Fibroblast", "Stromal", "Smooth_muscle", "Endothelial", "Endometrial")
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
ggsave("integrated_UMAP-annotated.pdf", plot = umap.combined, device = "pdf", width = 12, height = 7, units = "in")

Idents(seurat.object) <- "Type"
###### UMAP + Spatial DimPlots
p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("slice1.5")
p2 <- SpatialDimPlot(seurat.object, images=image.list, pt.size.factor = 3)
p2
image.list <- c("slice1.1")
p3 <- SpatialDimPlot(seurat.object, images=image.list)
p3
image.list <- c("slice1.2.1")
p4 <- SpatialDimPlot(seurat.object, images=image.list)
p4
image.list <- c("slice1.3.3")
p5 <- SpatialDimPlot(seurat.object, images=image.list)
p5
p6 <- wrap_plots(p1, p2, p3, p4, p5, ncol=5)
ggsave("Type_annotated_UMAP_SpatialDimPlots-cropped.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")

p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("slice1.5")
p2 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p2
image.list <- c("slice1.1")
p3 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p3
image.list <- c("slice1.2.1")
p4 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p4
image.list <- c("slice1.3.3")
p5 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p5
p6 <- wrap_plots(p1, p2, p3, p4, p5, ncol=5)
ggsave("Type_annotated_UMAP_SpatialDimPlots-uncropped.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")

Idents(seurat.object) <- "Subtype"
###### UMAP + Spatial DimPlots
p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("slice1.5")
p2 <- SpatialDimPlot(seurat.object, images=image.list, pt.size.factor = 3)
p2
image.list <- c("slice1.1")
p3 <- SpatialDimPlot(seurat.object, images=image.list)
p3
image.list <- c("slice1.2.1")
p4 <- SpatialDimPlot(seurat.object, images=image.list)
p4
image.list <- c("slice1.3.3")
p5 <- SpatialDimPlot(seurat.object, images=image.list)
p5
p6 <- wrap_plots(p1, p2, p3, p4, p5, ncol=5)
ggsave("Subtype_annotated_UMAP_SpatialDimPlots-cropped.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")

p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("slice1.5")
p2 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p2
image.list <- c("slice1.1")
p3 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p3
image.list <- c("slice1.2.1")
p4 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p4
image.list <- c("slice1.3.3")
p5 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p5
p6 <- wrap_plots(p1, p2, p3, p4, p5, ncol=5)
ggsave("Subtype_annotated_UMAP_SpatialDimPlots-uncropped.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")
##### Perform differential expression between Types using the raw data
subtypes <- seurat.object <- SCTransform(seurat.object, assay = "Spatial", verbose = T)
DefaultAssay(seurat.object) <- "SCT"
seurat.object.var <- seurat.object@assays[["SCT"]]@var.features
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = seurat.object.var, assay="Spatial")

Idents(seurat.object) <- "Type"
de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "Type_annotated_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("Type_annotated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("Type_annotated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Type_annotated_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Type_annotated_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Type_annotated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

##### Perform differential expression between Subtypes using the raw data
subtypes <- seurat.object <- SCTransform(seurat.object, assay = "Spatial", verbose = T)
DefaultAssay(seurat.object) <- "SCT"
seurat.object.var <- seurat.object@assays[["SCT"]]@var.features
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = seurat.object.var, assay="Spatial")

Idents(seurat.object) <- "Subtype"
de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "Subtype_annotated_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("Subtype_annotated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("Subtype_annotated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Subtype_annotated_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Subtype_annotated_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Subtype_annotated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()


setwd("..")
############################################################# ############################################################# #############################################################
############################################################# Gene set analysis #############################################################
############################################################# ############################################################# #############################################################

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
############################################################# ############################################################# #############################################################
############################################################# DE Analysis Btwn orig.ident and uninf v inf. #############################################################
############################################################# ############################################################# #############################################################

dir.create("library")
setwd("library")
Idents(seurat.object)
Idents(seurat.object) <- "orig.ident"

p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident")
ggsave("integrated_UMAP-orig.ident.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes, assay="Spatial")

Idents(seurat.object) <- "orig.ident"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "integrated_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("integrated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("integrated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("integrated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="orig.ident")
dev.off()

image.list <- c("slice1.5")
p2 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, pt.size.factor = 3, ncol=5)
p2
ggsave("top.markers_SpatialPlots-cropped-S07.pdf", plot = p2, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.1")
p3 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p3
ggsave("top.markers_SpatialPlots-cropped-S09.pdf", plot = p3, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.2.1")
p4 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p4
ggsave("top.markers_SpatialPlots-cropped-S12.pdf", plot = p4, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.3.3")
p5 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialPlots-cropped-S13.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

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





############ Treatment
Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
	# "Mock-enox-1" "Mock-enox-2" "ZIKV-sham"   "ZIKV-enox"  
new.metadata <- c("Mock", "Mock", "ZIKV-sham", "ZIKV-enox")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
## add cluster Subtypes as a metadata column
seurat.object$Treatment <- Idents(seurat.object)
Idents(seurat.object) <- "Treatment"
Idents(seurat.object)
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident")
ggsave("integrated_UMAP-Treatment.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "Treatment_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("Treatment_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("Treatment_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Treatment_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Treatment")
dev.off()
rm(new.metadata)
rm(feature.plot)
rm(top5.heatmap)

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)

image.list <- c("slice1.5")
p2 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, pt.size.factor = 3, ncol=5)
p2
ggsave("treatment_top.markers_SpatialPlots-cropped-S07.pdf", plot = p2, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.1")
p3 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p3
ggsave("treatment_top.markers_SpatialPlots-cropped-S09.pdf", plot = p3, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.2.1")
p4 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p4
ggsave("treatment_top.markers_SpatialPlots-cropped-S12.pdf", plot = p4, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.3.3")
p5 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("treatment_top.markers_SpatialPlots-cropped-S13.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


## SpatialPlots with all samples showing a particular gene
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=2)
p5
ggsave("C3_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=2, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("C3_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


## SpatialPlots with all samples showing a particular genes # DAF is CD55, inhibits complement
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
marker.list <- c("C1QA", "C1QC", "C1R", "C1S", "C2", "C3", "CD55", "CD59")
p5 <- SpatialPlot(seurat.object, features = marker.list, images=image.list, ncol=2)
p5
ggsave("complement.select.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")
Idents(seurat.object) <- "orig.ident"
top5.heatmap <- DoHeatmap(seurat.object, features = marker.list, raster = FALSE)
ggsave("complement.select.markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 8, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = marker.list, slot="counts", raster = FALSE)
ggsave("complement.select.markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 8, units = "in")



marker.list <- c("C1QA", "C1QC", "C2", "C3", "CD55", "CFH", "F3", "KLK1", "SERPING1", "TIMP2", "EHD1")
p5 <- SpatialPlot(seurat.object, features = marker.list, images=image.list, ncol=2)
p5
ggsave("complement.top.pathway.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")

# http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_COMPLEMENT.html
feature.plot <- DotPlot(seurat.object, features = marker.list)
png(file=paste0("complement.select-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Complement Markers") + scale_y_discrete(name ="Condition")
dev.off()
marker.list <- c("ACTN2", "ADAM9", "ADRA2B", "AKAP10", "ANG", "ANXA5", "APOA4", "APOBEC3F", "APOBEC3G", "APOC1", "ATOX1", "BRPF3", "C1QA", "C1QC", "C1R", "C1S", "C2", "C3", "C4BPB", "C9", "CA2", "CALM1", "CALM3", "CASP1", "CASP10", "CASP3", "CASP4", "CASP5", "CASP7", "CASP9", "CBLB", "CCL5", "CD36", "CD40LG", "CD46", "CD55", "CD59", "CDA", "CDH13", "CDK5R1", "CEBPB", "CFB", "CFH", "CLU", "COL4A2", "CP", "CPM", "CPQ", "CR1", "CR2", "CSRP1", "CTSB", "CTSC", "CTSD", "CTSH", "CTSL", "CTSO", "CTSS", "CTSV", "CXCL1", "DGKG", "DGKH", "DOCK10", "DOCK4", "DOCK9", "DPP4", "DUSP5", "DUSP6", "DYRK2", "EHD1", "ERAP2", "F10", "F2", "F3", "F5", "F7", "F8", "FCER1G", "FCN1", "FDX1", "FN1", "FYN", "GATA3", "GCA", "GMFB", "GNAI2", "GNAI3", "GNB2", "GNB4", "GNG2", "GNGT2", "GP1BA", "GP9", "GPD2", "GRB2", "GZMA", "GZMB", "GZMK", "HNF4A", "HPCAL4", "HSPA1A", "HSPA5", "IL6", "IRF1", "IRF2", "IRF7", "ITGAM", "ITIH1", "JAK2", "KCNIP2", "KCNIP3", "KIF2A", "KLK1", "KLKB1", "KYNU", "L3MBTL4", "LAMP2", "LAP3", "LCK", "LCP2", "LGALS3", "LGMN", "LIPA", "LRP1", "LTA4H", "LTF", "LYN", "MAFF", "ME1", "MMP12", "MMP13", "MMP14", "MMP15", "MMP8", "MSRB1", "MT3", "NOTCH4", "OLR1", "PCLO", "PCSK9", "PDGFB", "PDP1", "PFN1", "PHEX", "PIK3CA", "PIK3CG", "PIK3R5", "PIM1", "PLA2G4A", "PLA2G7", "PLAT", "PLAUR", "PLEK", "PLG", "PLSCR1", "PPP2CB", "PPP4C", "PRCP", "PRDM4", "PREP", "PRKCD", "PRSS3", "PRSS36", "PSEN1", "PSMB9", "RABIF", "RAF1", "RASGRP1", "RBSN", "RCE1", "RHOG", "RNF4", "S100A12", "S100A13", "S100A9", "SCG3", "SERPINA1", "SERPINB2", "SERPINC1", "SERPINE1", "SERPING1", "SH2B3", "SIRT6", "SPOCK2", "SRC", "STX4", "TFPI2", "TIMP1", "TIMP2", "TMPRSS6", "TNFAIP3", "USP14", "USP15", "USP16", "USP8", "VCPIP1", "WAS", "XPNPEP1", "ZEB1", "ZFPM2")

# Warning message:
# In FetchData(object = object, vars = features, cells = cells) :
#  The following requested variables were not found (16 ): STX4, SERPINA1, S100A12, FCN1, ERAP2, CTSV, CR1, CD59, CASP5, CASP10, CA2, C4BPB, C1S, C1R, APOBEC3G, APOBEC3F

feature.plot <- DotPlot(seurat.object, features = marker.list)
png(file=paste0("complement.pathway-DotPlot.png"),
                res=300, 
                width=7000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Complement Markers") + scale_y_discrete(name ="Condition")
dev.off()

top5.heatmap <- DoHeatmap(seurat.object, features = marker.list, raster = FALSE)
ggsave("Treatment_complement.pathway_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 5, height = 20, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = marker.list, slot="counts", raster = FALSE)
ggsave("Treatment_complement.pathway_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 5, height = 20, units = "in")



## SpatialPlots with all samples showing a particular gene
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
mock.markers <- c("TPBPA", "TPBPB", "CEACAM14", "PAPPA2", "PRL3A1")
p5 <- SpatialPlot(seurat.object, features = mock.markers, images=image.list, ncol=2)
p5
ggsave("Mock.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
ZIKV.SHAM.markers <- c("PRAP1", "KAP", "GUCA2B", "C3", "CTLA2A")
p5 <- SpatialPlot(seurat.object, features = ZIKV.SHAM.markers, images=image.list, ncol=2)
p5
ggsave("ZIKV.SHAM.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
ZIKV.enox.markers <- c("PRL2B1", "GHRH", "HBA-A1", "GJB2", "LEPR")
p5 <- SpatialPlot(seurat.object, features = ZIKV.enox.markers, images=image.list, ncol=2)
p5
ggsave("ZIKV.enox.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")


## SpatialPlots with all samples showing a particular gene
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
marker.list <- c("PRL3A1")
p5 <- SpatialPlot(seurat.object, features = marker.list, images=image.list, ncol=2)
p5
ggsave("JZ.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
marker.list <- c("CAR4", "CTS6")
p5 <- SpatialPlot(seurat.object, features = marker.list, images=image.list, ncol=2)
p5
ggsave("LZ.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
marker.list <- c("PRL2A1", "MT2", "CPE")
p5 <- SpatialPlot(seurat.object, features = marker.list, images=image.list, ncol=2)
p5
ggsave("Decidua.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")



##Markers from Simmons et al. BMC Genomics (2008) doi:10.1186/1471-2164-9-352
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
marker.list <- c("PRL3B1", "PRL2C", "PRL8A1", "PRL8A6", "PRL8A8", "PRL7A2", "PRL3A1", "PRL2B1", "PRL5A1")
p5 <- SpatialPlot(seurat.object, features = marker.list, images=image.list, ncol=2)
p5
ggsave("Simmons-Spongiotrophoblast.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")

##Markers from Simmons et al. BMC Genomics (2008) doi:10.1186/1471-2164-9-352
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
marker.list <- c("PCDH12", "PRL6A1", "PRL2A1", "PRL7B1", "PRL7C1")
p5 <- SpatialPlot(seurat.object, features = marker.list, images=image.list, ncol=2)
p5
ggsave("Simmons-GlycogenTrophoblast.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")

##Markers from Simmons et al. BMC Genomics (2008) doi:10.1186/1471-2164-9-352
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
marker.list <- c("TPBPA", "PRL4A1", "PRL8A9", "PRL7D1")
p5 <- SpatialPlot(seurat.object, features = marker.list, images=image.list, ncol=2)
p5
ggsave("Simmons-SpTandGlyT.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")


Idents(seurat.object) <- "orig.ident"

de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "integrated_DEGs_by.orig.ident_pos-0.2lnFC.txt", sep="\t")
marker.list <- c("SRGAP1", "NEAT1", "EVI5", "FN1", "KCNQ1OT1", "UBR4","FOS", "TIMM44", "RRM2B", "OSGIN1", "SCN7A", "VXN")

feature.plot <- DotPlot(seurat.object, features = marker.list)
png(file=paste0("miRNA.hits-DotPlot.png"),
                res=300, 
                width=7000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="miRNA Hits") + scale_y_discrete(name ="Condition")
dev.off()

top5.heatmap <- DoHeatmap(seurat.object, features = marker.list, raster = FALSE)
ggsave("Treatment_miRNA-interactions_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 5, height = 10, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = marker.list, slot="counts", raster = FALSE)
ggsave("Treatment_miRNA-interactions_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 5, height = 10, units = "in")


image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
marker.list <- c("SRGAP1")
p5 <- SpatialPlot(seurat.object, features = marker.list, images=image.list, ncol=2)
p5
ggsave("SRGAP1.miRNA.hits_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")

Idents(seurat.object) <- "Annotation"


setwd("..")
############################################################# ############################################################# #############################################################
seurat.object2 <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
seurat.object <- seurat.object2
rm(seurat.object2)
############################################################# ############################################################# #############################################################
############################################################# Subset Analysis: ZIKV-sham v. ZIKV-enox #############################################################
############################################################# ############################################################# #############################################################
dir.create("zikv.enox")
setwd("zikv.enox")
Idents(seurat.object)
Idents(seurat.object) <- "orig.ident"
seurat.object2 <- subset(x = seurat.object, idents = c("ZIKV-enox", "ZIKV-sham"))
ncol(seurat.object2)
	# 1544
seurat.object2 <- SCTransform(seurat.object2, assay = "Spatial", verbose = T)
DefaultAssay(seurat.object2) <- "SCT"
seurat.object2.var <- seurat.object2@assays[["SCT"]]@var.features

seurat.object2 <- ScaleData(seurat.object2, features = seurat.object2.var, assay = "SCT", vars.to.regress = c("CC.Difference", "nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"))

seurat.object2<- RunPCA(seurat.object2, features = seurat.object2.var, assay = "SCT", slot = "scale.data")
seurat.object2 <- FindNeighbors(seurat.object2, features = "seurat.object2.var", dims = 1:30)
seurat.object2 <- FindClusters(seurat.object2, resolution = 0.6)
seurat.object2 <- RunUMAP(seurat.object2, dims=1:30)
	## clusters 0 thru 12
p1 <- DimPlot(seurat.object2, reduction = "umap", group.by = "orig.ident")
p1
p2 <- DimPlot(seurat.object2, reduction = "umap", group.by = "seurat_clusters", label= "TRUE")
p2
p3 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Phase")
p3
umap.combined <- CombinePlots(plots = list(p1, p2, p3))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("integrated_UMAP.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object2 <- NormalizeData(seurat.object2, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object2 <- ScaleData(seurat.object2, features = seurat.object2.var, assay="Spatial")

Idents(seurat.object2) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object2, features = seurat.object2.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "integrated_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, raster = FALSE)
ggsave("integrated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("integrated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("integrated_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("integrated_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("integrated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

##### Perform differential expression between orig.ident
Idents(seurat.object2) <- "orig.ident"
de_markers <- FindAllMarkers(seurat.object2, features = seurat.object2.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "orig.ident_DEGs_pos-0.693lnFC.txt", sep="\t")

ggsave("orig.ident_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("orig.ident_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
View(top2)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("orig.ident_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="orig.ident")
dev.off()

rm(seurat.object2)

setwd("..")
############################################################# ############################################################# #############################################################
############################################################# 		#############################################################
############################################################# ############################################################# #############################################################

setwd("/home/ebarrozo/visium/results/seurat_mouse_v1")
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
save.image("murine_spatial_data-integrated-annotated_v1.RData")





