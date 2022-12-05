## visium_mouse_seurat_v2.1.R
	# v1.1, changed slice= slice1 to slice = {samplename}.slice, since after integration, there were 80 images and only 1 that works. 
	## v2- added S31- ZIKV-Enox-2, S32- Mock-Veh-1, S33- ZIKV-Veh-2, S34- ZIKV-Veh-3 datasets
		## note- spaceranger for these samples revealed manual image alignment worked best
				#note- seurat changed from ln(0.693) to log2fc(2) for DE filters
	##v2.1 - revised filters on S31-S34 to retain at least 50% of spots, especially percent.mt. See if lower quality transcriptomes cluster

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

## Analyze data on AagaardLab3
# http://10.16.5.106:8787/

set.seed(seed=1)
setwd("/home/ebarrozo/visium/results")
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
library(patchwork)
dir.create("seurat_mouse_v2.1")
setwd("seurat_mouse_v2.1")

setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1")

## This script will have several parts:
	# I. QC filters
	# II. Analysis of individual samples
	# III. Analysis of visium samples together
	# IV. Prediction of niches using mouse placenta snRNA-seq data (Marsh et al)
	# V. Analysis of visium and snRNA-seq data together (mouse placenta atlas)
	# VI. Annotation of visium spatial transcriptomics using data from parts III, IV, and V
	# VII. Subset analysis (eg ZIKV-veh only, macrophages, etc)


############################################################# ############################################################# #############################################################
############################################################# I. Load, QC, and filter each visium dataset #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# S07: Load data, Quality Control (QC) filtering #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S07-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S07.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S07.slice", filter.matrix = TRUE, to.upper = TRUE)

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
write.table(acute.viral@assays[["Spatial"]]@counts, "S07.viral.counts_all-cells.txt", sep="\t")
plot2 <- SpatialFeaturePlot(S07.manual, features = "percent.viral", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2
	#### Confirmed manually there are no ZIKV transcripts ########### 

acute.viral <- subset(acute.viral, percent.viral > 0)
	# Error: No cells found
#acute.viral <- subset(acute.viral, features = rownames(S07.manual)[rownames(S07.manual) %in% ZIKV.genes])
#write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

ncol(S07.manual)
	# 429 spots

S07.manual$orig.ident <- "ZIKV-Veh-1"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S07.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot1
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
qc.vlnplot
ggsave("S07.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S07.manual)
	#  429
S07.manual2 <- subset(S07.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt < 10 & percent.ribo < 15)
ncol(S07.manual2)
	# 	427
qc.vlnplot <- VlnPlot(S07.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S07.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S07.manual <- S07.manual2
rm(S07.manual2)
qc.vlnplot <- VlnPlot(S07.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S07.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## plot UMI counts/spot 
plot1 <- VlnPlot(S07.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot1
plot2 <- SpatialFeaturePlot(S07.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
		# pt.size.factor = 1.6 default
		# alpha - minimum and maximum transparency. Default is c(1, 1).
		# Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
		## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S07.manual_filtered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S07.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S07.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S07.manual_filtered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S09: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################
## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S09-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S09.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S07.slice", filter.matrix = TRUE, to.upper = TRUE)
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

S09.manual$orig.ident <- "Mock-Enox-1"   

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
	qc.vlnplot

S09.manual2 <- subset(S09.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000  & percent.mt < 10 & percent.ribo < 15)
ncol(S09.manual2)
	# 	1535
qc.vlnplot <- VlnPlot(S09.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S09.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S09.manual <- S09.manual2
rm(S09.manual2)
qc.vlnplot <- VlnPlot(S09.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S09.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
## plot UMI counts/spot 
plot1 <- VlnPlot(S09.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S09.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S09.manual_filtered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S09.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S09.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S09.manual_filtered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S12: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################
## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S12-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S12.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S12.slice", filter.matrix = TRUE, to.upper = TRUE)
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

S12.manual$orig.ident <- "Mock-Enox-2"   

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
	#  2302
	qc.vlnplot
S12.manual2 <- subset(S12.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000  & percent.mt < 10 & percent.ribo < 15)
ncol(S12.manual2)
	# 	2290
qc.vlnplot <- VlnPlot(S12.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S12.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S12.manual <- S12.manual2
rm(S12.manual2)
qc.vlnplot <- VlnPlot(S12.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S12.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## plot UMI counts/spot 
plot1 <- VlnPlot(S12.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S12.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S12.manual_filtered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S12.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S12.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S12.manual_filtered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")


############################################################# ############################################################# #############################################################
############################################################# S13: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S13-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S13.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S13.slice", filter.matrix = TRUE, to.upper = TRUE)
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

S13.manual$orig.ident <- "ZIKV-Enox-1"   

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
	#  1182
S13.manual2 <- subset(S13.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000  & percent.mt < 10 & percent.ribo < 15)
ncol(S13.manual2)
	# 	1135
qc.vlnplot <- VlnPlot(S13.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S13.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S13.manual <- S13.manual2
rm(S13.manual2)
qc.vlnplot <- VlnPlot(S13.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S13.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
## plot UMI counts/spot 
plot1 <- VlnPlot(S13.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S13.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S13.manual_filtered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S13.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S13.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S13.manual_filtered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

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
############################################################# S31: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S31-manual/outs'
list.files(data_dir)

#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S31.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S31.slice", filter.matrix = TRUE, to.upper = TRUE)
## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S31.manual)
head(old.names)
tail(old.names)
new.names <- gsub ("MM10-----------(.*)", "\\1", old.names)
new.names <- gsub ("NC-013532.1.V2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

S31.manual <- RenameGenesSeurat(obj = S31.manual, newnames = new.names)
row.names(S31.manual)

## Save a table of all.genes
all.genes <- row.names(S31.manual)
write.table(all.genes, "S31.manual.all.genes.txt", sep="\t")

## Assign and examine unfiltered quality control (QC) metrics
S31.manual <- PercentageFeatureSet(S31.manual, pattern = "^MT-", col.name = "percent.mt")
S31.manual <- PercentageFeatureSet(S31.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S31.manual <- PercentageFeatureSet(S31.manual, features = ZIKV.genes, col.name = "percent.viral")

## Log normalize data and multiply by a factor of 10000
	## required for cellcyclescoring
S31.manual<- NormalizeData(S31.manual, normalization.method = "LogNormalize", scale.factor = 10000)

# Run cell cycle scoring 
S31.manual <- CellCycleScoring(S31.manual, s.features = mouse.s.genes, g2m.features = mouse.g2m.genes, set.ident = FALSE)
S31.manual$CC.Difference <- S31.manual$S.Score - S31.manual$G2M.Score

acute.viral <- subset(S31.manual, features = rownames(S31.manual)[rownames(S31.manual) %in% ZIKV.genes])
# write.table(acute.viral@assays[["Spatial"]]@counts, "S31.manual_viral.counts_all-cells.txt", sep="\t")
acute.viral <- subset(acute.viral, percent.viral > 0)
	# Error: No cells found
#acute.viral <- subset(acute.viral, features = rownames(S31.manual)[rownames(S31.manual) %in% ZIKV.genes])
#write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

ncol(S31.manual)
	# 706 spots

S31.manual$orig.ident <- "ZIKV-Enox-2"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S31.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S31.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("S31.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S31.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S31.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot4 <- wrap_plots(plot1, plot2)
plot4
ggsave("S31.manual_unfiltered_nFeature_Spatial.pdf", plot = plot4, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S31.manual, features = "percent.mt", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S31.manual, features = "percent.mt") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot8 <- wrap_plots(plot1, plot2)
plot8
ggsave("S31.manual_unfiltered_percent.mt_Spatial.pdf", plot = plot6, device = "pdf", width = 8, height = 4, units = "in")


qc.vlnplot <- VlnPlot(S31.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S31.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S31.manual)
qc.vlnplot
	#  706
S31.manual2 <- subset(S31.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 150000  & percent.mt < 25 & percent.ribo < 10)
ncol(S31.manual2)
	# 	668
qc.vlnplot <- VlnPlot(S31.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
qc.vlnplot
ggsave("S31.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot2 <- SpatialFeaturePlot(S31.manual2, features = "percent.mt") + theme(legend.position = "right")


S31.manual <- S31.manual2
rm(S31.manual2)
qc.vlnplot <- VlnPlot(S31.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S31.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## plot UMI counts/spot 
plot1 <- VlnPlot(S31.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S31.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot5 <- wrap_plots(plot1, plot2)
plot5
ggsave("S31.manual_filtered_nCount_Spatial.pdf", plot = plot5, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S31.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S31.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot6 <- wrap_plots(plot1, plot2)
plot6
ggsave("S31.manual_filtered_nFeature_Spatial.pdf", plot = plot6, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S31.manual, features = "percent.mt", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S31.manual, features = "percent.mt") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot7 <- wrap_plots(plot1, plot2)
plot7
ggsave("S31.manual_filtered_percent.mt_Spatial.pdf", plot = plot6, device = "pdf", width = 8, height = 4, units = "in")

plot9 <- plot3 + plot4+ plot5+ plot6 + plot7 + plot8

ggsave("S31.manual_Spatial-filters.pdf", plot = plot9, device = "pdf", width = 8, height = 16, units = "in")


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
############################################################# S32: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S32-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S32.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S32.slice", filter.matrix = TRUE, to.upper = TRUE)
## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S32.manual)
new.names <- gsub ("MM10-----------(.*)", "\\1", old.names)
new.names <- gsub ("NC-013532.1.V2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

S32.manual <- RenameGenesSeurat(obj = S32.manual, newnames = new.names)
row.names(S32.manual)

## Save a table of all.genes
all.genes <- row.names(S32.manual)
write.table(all.genes, "S32.manual.all.genes.txt", sep="\t")

## Assign and examine unfiltered quality control (QC) metrics
S32.manual <- PercentageFeatureSet(S32.manual, pattern = "^MT-", col.name = "percent.mt")
S32.manual <- PercentageFeatureSet(S32.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S32.manual <- PercentageFeatureSet(S32.manual, features = ZIKV.genes, col.name = "percent.viral")

## Log normalize data and multiply by a factor of 10000
	## required for cellcyclescoring
S32.manual<- NormalizeData(S32.manual, normalization.method = "LogNormalize", scale.factor = 10000)

# Run cell cycle scoring 
S32.manual <- CellCycleScoring(S32.manual, s.features = mouse.s.genes, g2m.features = mouse.g2m.genes, set.ident = FALSE)
S32.manual$CC.Difference <- S32.manual$S.Score - S32.manual$G2M.Score

acute.viral <- subset(S32.manual, features = rownames(S32.manual)[rownames(S32.manual) %in% ZIKV.genes])
# write.table(acute.viral@assays[["Spatial"]]@counts, "S32.manual_viral.counts_all-cells.txt", sep="\t")
acute.viral <- subset(acute.viral, percent.viral > 0)
	# Error: No cells found
#acute.viral <- subset(acute.viral, features = rownames(S32.manual)[rownames(S32.manual) %in% ZIKV.genes])
#write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

ncol(S32.manual)
	# 875 spots

S32.manual$orig.ident <- "Mock-Veh-1"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S32.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S32.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("S32.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S32.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S32.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("S32.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S32.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S32.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S32.manual)
qc.vlnplot
	#  875
S32.manual2 <- subset(S32.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 150000  & percent.mt < 25 & percent.ribo < 10)
ncol(S32.manual2)
	# 	854
qc.vlnplot <- VlnPlot(S32.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
qc.vlnplot
ggsave("S32.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S32.manual <- S32.manual2
rm(S32.manual2)
qc.vlnplot <- VlnPlot(S32.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S32.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## plot UMI counts/spot 
plot1 <- VlnPlot(S32.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S32.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("S32.manual_filtered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S32.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S32.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("S32.manual_filtered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")


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
############################################################# S33: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S33-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S33.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S33.slice", filter.matrix = TRUE, to.upper = TRUE)
## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S33.manual)
new.names <- gsub ("MM10-----------(.*)", "\\1", old.names)
new.names <- gsub ("NC-013532.1.V2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

S33.manual <- RenameGenesSeurat(obj = S33.manual, newnames = new.names)
row.names(S33.manual)

## Save a table of all.genes
all.genes <- row.names(S33.manual)
write.table(all.genes, "S33.manual.all.genes.txt", sep="\t")

## Assign and examine unfiltered quality control (QC) metrics
S33.manual <- PercentageFeatureSet(S33.manual, pattern = "^MT-", col.name = "percent.mt")
S33.manual <- PercentageFeatureSet(S33.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S33.manual <- PercentageFeatureSet(S33.manual, features = ZIKV.genes, col.name = "percent.viral")

## Log normalize data and multiply by a factor of 10000
	## required for cellcyclescoring
S33.manual<- NormalizeData(S33.manual, normalization.method = "LogNormalize", scale.factor = 10000)

# Run cell cycle scoring 
S33.manual <- CellCycleScoring(S33.manual, s.features = mouse.s.genes, g2m.features = mouse.g2m.genes, set.ident = FALSE)
S33.manual$CC.Difference <- S33.manual$S.Score - S33.manual$G2M.Score

acute.viral <- subset(S33.manual, features = rownames(S33.manual)[rownames(S33.manual) %in% ZIKV.genes])
# write.table(acute.viral@assays[["Spatial"]]@counts, "S33.manual_viral.counts_all-cells.txt", sep="\t")
acute.viral <- subset(acute.viral, percent.viral > 0)
	# Error: No cells found
#acute.viral <- subset(acute.viral, features = rownames(S33.manual)[rownames(S33.manual) %in% ZIKV.genes])
#write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

ncol(S33.manual)
	# 1709 spots

S33.manual$orig.ident <- "ZIKV-Veh-2"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S33.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S33.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("S33.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S33.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S33.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("S33.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S33.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S33.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S33.manual)
qc.vlnplot
	#  1709
S33.manual2 <- subset(S33.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 150000  & percent.mt < 35 & percent.ribo < 10)
ncol(S33.manual2)
	# 	1698
qc.vlnplot <- VlnPlot(S33.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
qc.vlnplot
ggsave("S33.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot2 <- SpatialFeaturePlot(S33.manual2, features = "nFeature_Spatial") + theme(legend.position = "right")
plot2

S33.manual <- S33.manual2
rm(S33.manual2)
qc.vlnplot <- VlnPlot(S33.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S33.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- VlnPlot(S33.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S33.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("S33.manual_filtered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S33.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S33.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("S33.manual_filtered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")


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
############################################################# S34: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S34-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S34.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S34.slice", filter.matrix = TRUE, to.upper = TRUE)
## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S34.manual)
new.names <- gsub ("MM10-----------(.*)", "\\1", old.names)
new.names <- gsub ("NC-013532.1.V2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

S34.manual <- RenameGenesSeurat(obj = S34.manual, newnames = new.names)
row.names(S34.manual)

## Save a table of all.genes
all.genes <- row.names(S34.manual)
write.table(all.genes, "S34.manual.all.genes.txt", sep="\t")

## Assign and examine unfiltered quality control (QC) metrics
S34.manual <- PercentageFeatureSet(S34.manual, pattern = "^MT-", col.name = "percent.mt")
S34.manual <- PercentageFeatureSet(S34.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S34.manual <- PercentageFeatureSet(S34.manual, features = ZIKV.genes, col.name = "percent.viral")

## Log normalize data and multiply by a factor of 10000
	## required for cellcyclescoring
S34.manual<- NormalizeData(S34.manual, normalization.method = "LogNormalize", scale.factor = 10000)

# Run cell cycle scoring 
S34.manual <- CellCycleScoring(S34.manual, s.features = mouse.s.genes, g2m.features = mouse.g2m.genes, set.ident = FALSE)
S34.manual$CC.Difference <- S34.manual$S.Score - S34.manual$G2M.Score

acute.viral <- subset(S34.manual, features = rownames(S34.manual)[rownames(S34.manual) %in% ZIKV.genes])
# write.table(acute.viral@assays[["Spatial"]]@counts, "S34.manual_viral.counts_all-cells.txt", sep="\t")
acute.viral <- subset(acute.viral, percent.viral > 0)
	# Error: No cells found
#acute.viral <- subset(acute.viral, features = rownames(S34.manual)[rownames(S34.manual) %in% ZIKV.genes])
#write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

ncol(S34.manual)
	# 1159 spots

S34.manual$orig.ident <- "ZIKV-Veh-3"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S34.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S34.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("S34.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S34.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S34.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("S34.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S34.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S34.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S34.manual)
qc.vlnplot
	#  1159
S34.manual2 <- subset(S34.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 150000  & percent.mt < 35 & percent.ribo < 10)
ncol(S34.manual2)
	# 	1136
qc.vlnplot <- VlnPlot(S34.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
qc.vlnplot
ggsave("S34.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot2 <- SpatialFeaturePlot(S34.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
plot2

S34.manual <- S34.manual2
rm(S34.manual2)
qc.vlnplot <- VlnPlot(S34.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S34.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

## plot UMI counts/spot 
plot1 <- VlnPlot(S34.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S34.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("S34.manual_filtered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S34.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S34.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("S34.manual_filtered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")


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
############################################################# II. Individual Visium sample analysis #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# S07: Dimension Reduction, Visualization, and Differential Expression (DE) #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
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
p3 <- SpatialFeaturePlot(S07.manual, features = top.features, ncol = 1, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S07.manual_spatial.var.features_ncol1.pdf", plot = p3, device = "pdf", width = 4, height = 16, units = "in")


##### Perform differential expression between clusters
S07.manual <- NormalizeData(S07.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S07.manual <- ScaleData(S07.manual, features = S07.manual.var, assay="Spatial")

Idents(S07.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S07.manual, features = S07.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "S07_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>2)
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
p3 <- SpatialFeaturePlot(S09.manual, features = top.features, ncol = 1, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S09.manual_spatial.var.features_ncol1.pdf", plot = p3, device = "pdf", width = 4, height = 16, units = "in")


##### Perform differential expression between clusters
S09.manual <- NormalizeData(S09.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S09.manual <- ScaleData(S09.manual, features = S09.manual.var, assay="Spatial")

Idents(S09.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S09.manual, features = S09.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "S09_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>2)
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
p3 <- SpatialFeaturePlot(S12.manual, features = top.features, ncol = 1, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S12.manual_spatial.var.features_ncol1.pdf", plot = p3, device = "pdf", width = 4, height = 16, units = "in")


##### Perform differential expression between clusters
S12.manual <- NormalizeData(S12.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S12.manual <- ScaleData(S12.manual, features = S12.manual.var, assay="Spatial")

Idents(S12.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S12.manual, features = S12.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "S12_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>2)
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
p3 <- SpatialFeaturePlot(S13.manual, features = top.features, ncol = 1, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S13.manual_spatial.var.features_ncol1.pdf", plot = p3, device = "pdf", width = 4, height = 16, units = "in")


##### Perform differential expression between clusters
S13.manual <- NormalizeData(S13.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S13.manual <- ScaleData(S13.manual, features = S13.manual.var, assay="Spatial")

Idents(S13.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S13.manual, features = S13.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "S13_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>2)
		## cluster 0 does not have sig. genes


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
############################################################# S31: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S31.manual <- SCTransform(S31.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S31.manual <- RunPCA(S31.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S31.manual)
ggsave("S31.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S31.manual <- FindNeighbors(S31.manual, reduction = "pca", dims = 1:30)
S31.manual <- FindClusters(S31.manual, verbose = T)
S31.manual <- RunUMAP(S31.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S31.manual, reduction = "umap", label = TRUE)
p1
	##### clusters 0 thru 5
p2 <- SpatialDimPlot(S31.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S31.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S31.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters

p3<- SpatialDimPlot(S31.manual, cells.highlight = CellsByIdentities(object = S31.manual, 
	## edit by number of clusters
	idents = c(0, 1, 2, 3, 4, 5)), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0)
ggsave("S31.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S31.manual <- FindSpatiallyVariableFeatures(S31.manual, assay = "SCT", features = VariableFeatures(S31.manual)[1:1000], selection.method = "markvariogram")
S31.manual.var <- S31.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S31.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S31.manual, features = top.features, ncol = 3, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S31.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")
p3 <- SpatialFeaturePlot(S31.manual, features = top.features, ncol = 1, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S31.manual_spatial.var.features_ncol1.pdf", plot = p3, device = "pdf", width = 4, height = 16, units = "in")


##### Perform differential expression between clusters
S31.manual <- NormalizeData(S31.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S31.manual <- ScaleData(S31.manual, features = S31.manual.var, assay="Spatial")

Idents(S31.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S31.manual, features = S31.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "S31_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>2)


top5.heatmap <- DoHeatmap(S31.manual, features = top5$gene, raster = FALSE)
ggsave("S31.manual_top5_markers_log2fc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S31.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S31.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S31.manual, features = unique.top2)
png(file=paste0("S31.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S31.manual, features = unique.top2)
png(file=paste0("S31.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S31.manual, features = unique.top2)
png(file=paste0("S31.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
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
############################################################# S32: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S32.manual <- SCTransform(S32.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S32.manual <- RunPCA(S32.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S32.manual)
ggsave("S32.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S32.manual <- FindNeighbors(S32.manual, reduction = "pca", dims = 1:30)
S32.manual <- FindClusters(S32.manual, verbose = T)
S32.manual <- RunUMAP(S32.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S32.manual, reduction = "umap", label = TRUE)
p1
	##### clusters 0 thru 5
p2 <- SpatialDimPlot(S32.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S32.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S32.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters

p3<- SpatialDimPlot(S32.manual, cells.highlight = CellsByIdentities(object = S32.manual, 
	## edit by number of clusters
	idents = c(0, 1, 2, 3, 4, 5)), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0)
ggsave("S32.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S32.manual <- FindSpatiallyVariableFeatures(S32.manual, assay = "SCT", features = VariableFeatures(S32.manual)[1:1000], selection.method = "markvariogram")
S32.manual.var <- S32.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S32.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S32.manual, features = top.features, ncol = 3, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S32.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")
p3 <- SpatialFeaturePlot(S32.manual, features = top.features, ncol = 1, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S32.manual_spatial.var.features_ncol1.pdf", plot = p3, device = "pdf", width = 4, height = 16, units = "in")


##### Perform differential expression between clusters
S32.manual <- NormalizeData(S32.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S32.manual <- ScaleData(S32.manual, features = S32.manual.var, assay="Spatial")

Idents(S32.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S32.manual, features = S32.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "S32_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>2)


top5.heatmap <- DoHeatmap(S32.manual, features = top5$gene, raster = FALSE)
ggsave("S32.manual_top5_markers_log2fc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S32.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S32.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S32.manual, features = unique.top2)
png(file=paste0("S32.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S32.manual, features = unique.top2)
png(file=paste0("S32.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S32.manual, features = unique.top2)
png(file=paste0("S32.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
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
############################################################# S33: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S33.manual <- SCTransform(S33.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S33.manual <- RunPCA(S33.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S33.manual)
ggsave("S33.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S33.manual <- FindNeighbors(S33.manual, reduction = "pca", dims = 1:30)
S33.manual <- FindClusters(S33.manual, verbose = T)
S33.manual <- RunUMAP(S33.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S33.manual, reduction = "umap", label = TRUE)
p1
	##### clusters 0 thru 6
p2 <- SpatialDimPlot(S33.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S33.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S33.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters

p3<- SpatialDimPlot(S33.manual, cells.highlight = CellsByIdentities(object = S33.manual, 
	## edit by number of clusters
	idents = c(0, 1, 2, 3, 4, 5, 6)), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0)
ggsave("S33.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S33.manual <- FindSpatiallyVariableFeatures(S33.manual, assay = "SCT", features = VariableFeatures(S33.manual)[1:1000], selection.method = "markvariogram")
S33.manual.var <- S33.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S33.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S33.manual, features = top.features, ncol = 3, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S33.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")
p3 <- SpatialFeaturePlot(S33.manual, features = top.features, ncol = 1, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S33.manual_spatial.var.features_ncol1.pdf", plot = p3, device = "pdf", width = 4, height = 16, units = "in")


##### Perform differential expression between clusters
S33.manual <- NormalizeData(S33.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S33.manual <- ScaleData(S33.manual, features = S33.manual.var, assay="Spatial")

Idents(S33.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S33.manual, features = S33.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "S33_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>2)


top5.heatmap <- DoHeatmap(S33.manual, features = top5$gene, raster = FALSE)
ggsave("S33.manual_top5_markers_log2fc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S33.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S33.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S33.manual, features = unique.top2)
png(file=paste0("S33.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S33.manual, features = unique.top2)
png(file=paste0("S33.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S33.manual, features = unique.top2)
png(file=paste0("S33.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
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
############################################################# S34: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S34.manual <- SCTransform(S34.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S34.manual <- RunPCA(S34.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S34.manual)
ggsave("S34.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S34.manual <- FindNeighbors(S34.manual, reduction = "pca", dims = 1:30)
S34.manual <- FindClusters(S34.manual, verbose = T)
S34.manual <- RunUMAP(S34.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S34.manual, reduction = "umap", label = TRUE)
p1
	##### clusters 0 thru 4
p2 <- SpatialDimPlot(S34.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S34.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S34.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters

p3<- SpatialDimPlot(S34.manual, cells.highlight = CellsByIdentities(object = S34.manual, 
	## edit by number of clusters
	idents = c(0, 1, 2, 3, 4)), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0)
ggsave("S34.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S34.manual <- FindSpatiallyVariableFeatures(S34.manual, assay = "SCT", features = VariableFeatures(S34.manual)[1:1000], selection.method = "markvariogram")
S34.manual.var <- S34.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S34.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S34.manual, features = top.features, ncol = 3, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S34.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")
p3 <- SpatialFeaturePlot(S34.manual, features = top.features, ncol = 1, alpha = c(0.1, 1), pt.size.factor = 4.0)
ggsave("S34.manual_spatial.var.features_ncol1.pdf", plot = p3, device = "pdf", width = 4, height = 16, units = "in")


##### Perform differential expression between clusters
S34.manual <- NormalizeData(S34.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S34.manual <- ScaleData(S34.manual, features = S34.manual.var, assay="Spatial")

Idents(S34.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S34.manual, features = S34.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "S34_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>2)


top5.heatmap <- DoHeatmap(S34.manual, features = top5$gene, raster = FALSE)
ggsave("S34.manual_top5_markers_log2fc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S34.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S34.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S34.manual, features = unique.top2)
png(file=paste0("S34.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S34.manual, features = unique.top2)
png(file=paste0("S34.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S34.manual, features = unique.top2)
png(file=paste0("S34.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
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
############################################################# III. analysis of visium data together #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# Visium integration: merge all data using CCA anchors #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################

## Slim down the seurat.objects to save space. Don't need to keep the scale.data
S07.manual2 <- DietSeurat(S07.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
View(S07.manual2)
	## from 188MB to 166 MB
S07.manual <- S07.manual2
rm(S07.manual2)
S09.manual <- DietSeurat(S09.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S12.manual <- DietSeurat(S12.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S13.manual <- DietSeurat(S13.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S31.manual <- DietSeurat(S31.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S32.manual <- DietSeurat(S32.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S33.manual <- DietSeurat(S33.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S34.manual <- DietSeurat(S34.manual, counts=TRUE, data=TRUE, scale.data = FALSE)

## combine lists of top variable features for later DE analysis/clustering
var.combined <- union(S07.manual.var, S09.manual.var)
var.combined <- union(var.combined, S12.manual.var)
var.combined <- union(var.combined, S13.manual.var)
var.combined <- union(var.combined, S31.manual.var)
var.combined <- union(var.combined, S32.manual.var)
var.combined <- union(var.combined, S33.manual.var)
var.combined <- union(var.combined, S34.manual.var)
	## 9010 combined variable features

## Merge objects
all_merged <- merge(x = S07.manual, y = c(S09.manual, S12.manual, S13.manual, S31.manual, S32.manual, S33.manual, S34.manual), merge.data = TRUE, project = "all_merged")
plot2 <- SpatialFeaturePlot(all_merged, features = "percent.viral", ncol = 1) + patchwork::plot_layout(ncol = 1)
####### FOR SOME REASON, INTEGRATION ADDS IMAGES ARBITRARILY
ggsave("all_merged_SpatialDimPlot_percent.viral_TEST.pdf", plot = plot2, device = "pdf", width = 45, height = 11, units = "in")

Idents(all_merged) <- "orig.ident"
## make combined QC plot 
qc.vlnplot <- VlnPlot(all_merged, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), group.by = "orig.ident", pt.size = 0.000001, ncol = 4) + NoLegend()
ggsave("all_merged_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 15, height = 5, units = "in")


## plot UMI counts/spot 
plot1 <- VlnPlot(all_merged, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(all_merged, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("all_merged_filtered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(all_merged, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(all_merged, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave("all_merged_filtered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

rm(plot2)
rm(plot3)
rm(plot1)
rm(qc.vlnplot)
rm(S07.manual)
rm(S09.manual)
rm(S12.manual)
rm(S13.manual)
rm(S31.manual)
rm(S32.manual)
rm(S33.manual)
rm(S34.manual)

rm(S07.manual.var)
rm(S09.manual.var)
rm(S12.manual.var)
rm(S13.manual.var)
rm(S31.manual.var)
rm(S32.manual.var)
rm(S33.manual.var)
rm(S34.manual.var)

options(future.globals.maxSize = 300000000000)

# setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1")
# save.image("murine_spatial_data-filtered_v1.RData")
# load("/home/ebarrozo/visium/results/seurat_mouse_v2/murine_spatial_data-filtered_v1.RData")


all_merged <- SCTransform(all_merged, assay = "Spatial", verbose = T)
all_merged.list <- SplitObject(all_merged, split.by = "orig.ident")


# Make a reference list
reference.list <- all_merged.list[c("Mock-Enox-1", "Mock-Enox-2", "ZIKV-Veh-1", "ZIKV-Enox-1", "ZIKV-Enox-2", "Mock-Veh-1", "ZIKV-Veh-2", "ZIKV-Veh-3")]

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

setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1")
save.image("murine_spatial_data-integrated_v1.RData")
# load("murine_spatial_data-integrated_v1.RData")
############################################################# ############################################################# #############################################################
############################################################# Visium: Dimension Reduction and Visualization #############################################################
############################################################# ############################################################# #############################################################
library(ggpubr)
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

Idents(seurat.object) <- "orig.ident"
levels(seurat.object) # "Mock-Enox-1" "Mock-Enox-2" "ZIKV-Veh-1"  "ZIKV-Enox-1" "ZIKV-Enox-2" "Mock-Veh-1"  "ZIKV-Veh-2" "ZIKV-Veh-3"
ncol(seurat.object)	# 9743 spots
new.metadata <- c("Mock-Enox", "Mock-Enox", "ZIKV-Veh", "ZIKV-Enox", "ZIKV-Enox", "Mock-Veh", "ZIKV-Veh", "ZIKV-Veh")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Treatment <- Idents(seurat.object)
Idents(seurat.object) <- "Treatment"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object) # "Mock-Enox-1" "Mock-Enox-2" "ZIKV-Veh-1"  "ZIKV-Enox-1" "ZIKV-Enox-2" "Mock-Veh-1"  "ZIKV-Veh-2" "ZIKV-Veh-3"
ncol(seurat.object)	# 9743 spots
new.metadata <- c("1", "1", "1", "1", "2", "2", "2", "2")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Batch <- Idents(seurat.object)
Idents(seurat.object) <- "Batch"

seurat.object <- ScaleData(seurat.object, features = integrated.genes, assay = "integrated", vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"))
seurat.object<- RunPCA(seurat.object, features = integrated.genes, assay = "integrated", slot = "scale.data")

seurat.object <- FindNeighbors(seurat.object, features = "integrated.genes", dims = 1:30)
seurat.object <- FindClusters(seurat.object, resolution = 0.6)
seurat.object <- RunUMAP(seurat.object, dims=1:30)
	## clusters 0 thru 12
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "Treatment", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Treatment')
p1
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p2
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Sample')
p3
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Batch", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Batch')
p4
umap.combined <- p1 + p2 + p3 + p4
ggsave("integrated_UMAP.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("UMAP_XL.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")

Idents(seurat.object) <- "seurat_clusters"
t1<-table(Idents(seurat.object))
write.table(t1, "seurat_clusters.counts.txt", sep="\t")
rm(t1)

# setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1")
# seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, , dimreducs=c("umap","pca"))
# save.image("murine_spatial_data-clustered_v1.RData")
# load("murine_spatial_data-clustered_v1.RData")

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
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Phase", label= "TRUE", repel=TRUE, raster=T)

umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8,p4))
ggsave("integrated_UMAP_QCmetricsFeaturePlots.pdf", plot = umap.combined, device = "pdf", width = 10, height = 6, units = "in")


plot2 <- SpatialFeaturePlot(seurat.object, features = "percent.viral", ncol = 1) + patchwork::plot_layout(ncol = 1)
####### FOR SOME REASON, INTEGRATION ADDS IMAGES ARBITRARILY
ggsave("integrated_SpatialDimPlot_percent.viral_TEST.pdf", plot = plot2, device = "pdf", width = 8, height = 45, units = "in")


## Using the test images above, copy and paste the images that work below. 
    ## identified a pattern and confirmed these all worked. 

#S07.slice.1					S09
#S12.slice.1					S12	
#S07.slice.3		ZIKV-Veh-1 S07
#S13.slice.3					S13
#S31.slice.4		S31 and the rest are correct
#S32.slice.5
#S33.slice.6
#S34.slice.7


image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_percent.viral.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
plot2<- SpatialDimPlot(seurat.object, images=image.list) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_cropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialDimPlot(seurat.object, images=image.list, crop=FALSE) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")



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
############################################################# Visium: DE Analysis_clusters_var.genes #############################################################
############################################################# ############################################################# #############################################################
dir.create("DE_clusters")
setwd("DE_clusters")

seurat.object <- SCTransform(seurat.object, assay = "Spatial", verbose = T)
DefaultAssay(seurat.object) <- "SCT"
seurat.object.var <- seurat.object@assays[["SCT"]]@var.features

Idents(seurat.object) <- "orig.ident"
sex.genes <- c("SRY", "DDX3Y")
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=sex.genes)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = sex.genes, raster = FALSE, slot="counts", assay="SCT", size=3, angle = 90) + NoLegend()
ggsave("fetalsex_mc_heatmap_orig.ident.png", plot = library.averages.heatmap, device = "png", width = 11, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = sex.genes, raster = FALSE, assay="SCT", size=3, angle = 90) + NoLegend()
ggsave("fetalsex_FC_heatmap_orig.ident.png", plot = library.averages.heatmap, device = "png", width = 11, height = 5, units = "in")
write.csv(cluster.averages@assays[["SCT"]]@counts, file = "orig.ident_DDX3Y_meancounts.csv")


##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = seurat.object.var, assay="Spatial")

Idents(seurat.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>2)
		## all clusters have sig. genes
		## only clusters 3, 5, 6, 7, 8, and 10 have significant genes passing the filter
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

image.list <- c("S07.slice.3")
p2 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, pt.size.factor = 3, ncol=2)
p2
ggsave("top.markers_SpatialFeaturePlots-cropped-S07.pdf", plot = p2, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S07.slice.1")
p3 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=2)
p3
ggsave("top.markers_SpatialFeaturePlots-cropped-S09.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S12.slice.1")
p4 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=2)
p4
ggsave("top.markers_SpatialFeaturePlots-cropped-S12.pdf", plot = p4, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S13.slice.3")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=2)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S13.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S31.slice.4")
p2 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, pt.size.factor = 3, ncol=2)
p2
ggsave("top.markers_SpatialFeaturePlots-cropped-S31.pdf", plot = p2, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S32.slice.5")
p3 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=2)
p3
ggsave("top.markers_SpatialFeaturePlots-cropped-S32.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S33.slice.6")
p4 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=2)
p4
ggsave("top.markers_SpatialFeaturePlots-cropped-S33.pdf", plot = p4, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S34.slice.7")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=2)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S34.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")


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
setwd("..")
############################################################# ############################################################# #############################################################
############################################################# Visium: DE Analysis-all.genes #############################################################
############################################################# ############################################################# #############################################################
dir.create("DE_clusters-all.genes")
setwd("DE_clusters-all.genes")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")

Idents(seurat.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-log2fc.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>2)
		## all clusters have sig. genes
		## only clusters 3, 5, 6, 7, 8, and 10 have significant genes passing the filter
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

image.list <- c("S07.slice.3")
p2 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, pt.size.factor = 3, ncol=2)
p2
ggsave("top.markers_SpatialFeaturePlots-cropped-S07.pdf", plot = p2, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S07.slice.1")
p3 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=2)
p3
ggsave("top.markers_SpatialFeaturePlots-cropped-S09.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S12.slice.1")
p4 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=2)
p4
ggsave("top.markers_SpatialFeaturePlots-cropped-S12.pdf", plot = p4, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S13.slice.3")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=2)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S13.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S31.slice.4")
p2 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, pt.size.factor = 3, ncol=2)
p2
ggsave("top.markers_SpatialFeaturePlots-cropped-S31.pdf", plot = p2, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S32.slice.5")
p3 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=2)
p3
ggsave("top.markers_SpatialFeaturePlots-cropped-S32.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S33.slice.6")
p4 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=2)
p4
ggsave("top.markers_SpatialFeaturePlots-cropped-S33.pdf", plot = p4, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S34.slice.7")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=2)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S34.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")


Idents(seurat.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0, logfc.threshold =  0)
write.table(de_markers, "integrated_DEGs_byclusters_pos-nofilter.txt", sep="\t")


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

setwd("..")
############################################################# ############################################################# #############################################################
############################################################# IV. Predictions of spatial transcriptome niches using snRNA-seq data #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# Load mouse placenta snRNA-seq dataset ############################## #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################

set.seed(seed=1)
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
library(patchwork)
setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1")
load("murine_spatial_data-clustered_v1.RData")

setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1")
dir.create("d")
setwd("atlas")

# see Marsh_snRNA-seq_v2.R for creation of snRNA-seq object
## integrating our visium spatial transcriptomics data from mouse placenta with Marsh, et al. (2020) snRNA-seq from mouse placenta.

## Analysis of Marsh & Blelloch, eLife, 2020 murine placental snRNA-seq
# GSE152248 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152248
## E9.5-E14.5: E9.5, E10.5, E12.5, E14.5

## 17500 nuclei, 10X genomics 3' v3 chemistry, 20-30k reads/nucleus
## Filters: >500 <4000 nFeatures >25% mito
## 27,326 nuclei after QC/ 26 clusters (0.6res), XIST for male embryonic tissue
## 5 broad cell types, further divided into subgroups, 16,386 trophoblasts reclustered

## add some more metadata to call upon later
Idents(seurat.object) <- "orig.ident"
levels(seurat.object) # "Mock-Enox-1" "Mock-Enox-2" "ZIKV-Veh-1"  "ZIKV-Enox-1" "ZIKV-Enox-2" "Mock-Veh-1"  "ZIKV-Veh-2" "ZIKV-Veh-3"
ncol(seurat.object)	# 9743 spots
new.metadata <- c("Visium", "Visium", "Visium", "Visium", "Visium", "Visium", "Visium", "Visium")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Platform <- Idents(seurat.object)
Idents(seurat.object) <- "Platform"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object) # "Mock-Enox-1" "Mock-Enox-2" "ZIKV-Veh-1"  "ZIKV-Enox-1" "ZIKV-Enox-2" "Mock-Veh-1"  "ZIKV-Veh-2" "ZIKV-Veh-3"
ncol(seurat.object)	# 9743 spots
new.metadata <- c("Mock-Enox", "Mock-Enox", "ZIKV-Veh", "ZIKV-Enox", "ZIKV-Enox", "Mock-Veh", "ZIKV-Veh", "ZIKV-Veh")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Treatment <- Idents(seurat.object)
Idents(seurat.object) <- "Treatment"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object) # "Mock-Enox-1" "Mock-Enox-2" "ZIKV-Veh-1"  "ZIKV-Enox-1" "ZIKV-Enox-2" "Mock-Veh-1"  "ZIKV-Veh-2" "ZIKV-Veh-3"
ncol(seurat.object)	# 9743 spots
new.metadata <- c("Mock", "Mock", "ZIKV", "ZIKV", "ZIKV", "Mock", "ZIKV", "ZIKV")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Virus <- Idents(seurat.object)
Idents(seurat.object) <- "Virus"

## make a copy of the spatial data in case the Marsh object is saved with the same name
visium.object <-seurat.object


load("/home/ebarrozo/Marsh/results/seurat-v1/marsh_annotated_v1.RData")
	## see # Marsh_snRNA-seq_v2.R	## saved as marsh.seurat.object
DefaultAssay(seurat.object) <- "SCT"
DefaultAssay(marsh.seurat.object) <- "SCT"

Idents(marsh.seurat.object) <- "orig.ident"
levels(marsh.seurat.object) # Uninfected
ncol(marsh.seurat.object)	# 4333 spots
new.metadata <- c("snRNA-seq")
names(new.metadata) <- levels(marsh.seurat.object)
marsh.seurat.object <- RenameIdents(marsh.seurat.object, new.metadata)
marsh.seurat.object$Platform <- Idents(marsh.seurat.object)
Idents(marsh.seurat.object) <- "Platform"



############################################################# ############################################################# #############################################################
############################################################# Predictions of Visium Niches using snRNA-seq Data #############################################################
############################################################# ############################################################# #############################################################
setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/atlas")
dir.create("predictions_v1")
setwd("predictions_v1")
## Run predictions for spatial transcriptome identities based on single-nuclei transcriptome annotations
Idents(marsh.seurat.object) <- "Type.Profile"
levels(marsh.seurat.object)
 ## [1] "NK"                       "Stromal"                  "Syncytiotrophoblast"      "Endometrial"             
	#[5] "Fibroblast"               "Villious_cytotrophoblast" "Endothelial"              "Erythroblast"            
# [9] "Macrophage" 

# change to VCT and SYT
Idents(marsh.seurat.object) <- "Type.Profile"
levels(marsh.seurat.object) 
 ## [1] "NK"                       "Stromal"                  "Syncytiotrophoblast"      "Endometrial"             
	#[5] "Fibroblast"               "Villious_cytotrophoblast" "Endothelial"              "Erythroblast"            
# [9] "Macrophage" 
ncol(marsh.seurat.object)	# 9743 spots
new.metadata <- c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma")
names(new.metadata) <- levels(marsh.seurat.object)
marsh.seurat.object <- RenameIdents(marsh.seurat.object, new.metadata)
marsh.seurat.object$Type <- Idents(marsh.seurat.object)
Idents(marsh.seurat.object) <- "Type.Profile"


type.list <-  levels(marsh.seurat.object)

## Run SCTransform on marsh data and visuaize
# marsh.seurat.object <- SCTransform(marsh.seurat.object, ncells = 3000, return.only.var.genes=TRUE, conserve.memory = TRUE, verbose = TRUE) 
# marsh.seurat.object <- SCTransform(marsh.seurat.object, ncells = 4333,method = "glmGamPoi",return.only.var.genes=TRUE, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt", "Phase"), conserve.memory = F, verbose = TRUE) 

## remove old umap and pca 
marsh.seurat.object <- DietSeurat(marsh.seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE)
marsh.seurat.object <- SCTransform(marsh.seurat.object, ncells = 4333,method = "glmGamPoi",return.only.var.genes=TRUE, vars.to.regress = c("Phase"), conserve.memory = F, verbose = TRUE) 
marsh.seurat.object <- RunPCA(marsh.seurat.object, assay="SCT", verbose = TRUE)  
var.genes <- marsh.seurat.object@assays[["SCT"]]@var.features
marsh.seurat.object <- FindNeighbors(marsh.seurat.object, assay="SCT",  features = "var.genes", dims = 1:30)
marsh.seurat.object <- FindClusters(marsh.seurat.object, assay="SCT", resolution = 0.6)
marsh.seurat.object <- RunUMAP(marsh.seurat.object, dims = 1:30)
p4 <- DimPlot(marsh.seurat.object, reduction = "umap", group.by = "Type.Profile", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
p4 		## not changing after var.reg; try removing dim.reducs with diet ;; forgot marsh.seurat.object <- in front of RunUMAP ;
p2 <- DimPlot(marsh.seurat.object, reduction = "umap", group.by = "Phase", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Phase')
p2   ## looks like by default, Phase heavily influences default clustering. Try regressing out CC differences.;; phase still has an effect. ; try glamGlmPoi; better, let's try regressing out only Phase
p3 <- FeaturePlot(marsh.seurat.object, reduction = "umap", features = 'percent.mt', raster=FALSE, cols = mypal3) + labs(title = NULL, color='percent.mt')
p3
p1 <- DimPlot(marsh.seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p1 
#examine UMAPs with qc metrics 
p5 <- FeaturePlot(marsh.seurat.object, features = 'nFeature_RNA')
p6 <- FeaturePlot(marsh.seurat.object, features = 'nCount_RNA')
p7 <- FeaturePlot(marsh.seurat.object, features = 'percent.mt')
p8 <- FeaturePlot(marsh.seurat.object, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
umap.combined
ggsave("Marsh_UMAP_QCmetricsFeaturePlots-RNA.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")

umap.combined <- p1+p4+p2+p3
ggsave("Marsh_UMAP.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("Marsh_XL.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")

## revise Marsh cluster annotations
############################################################# ############################################################# #############################################################
dir.create("DE_clusters-marsh")
setwd("DE_clusters-marsh")

marsh.seurat.object <- SCTransform(marsh.seurat.object, assay = "RNA", verbose = T)
DefaultAssay(marsh.seurat.object) <- "SCT"
marsh.seurat.object.var <- marsh.seurat.object@assays[["SCT"]]@var.features

Idents(marsh.seurat.object) <- "orig.ident"
sex.genes <- c("SRY", "DDX3Y")
cluster.averages <- AverageExpression(marsh.seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=sex.genes)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = sex.genes, raster = FALSE, slot="counts", assay="SCT", size=3, angle = 90) + NoLegend()
ggsave("fetalsex_mc_heatmap_orig.ident.png", plot = library.averages.heatmap, device = "png", width = 11, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = sex.genes, raster = FALSE, assay="SCT", size=3, angle = 90) + NoLegend()
ggsave("fetalsex_FC_heatmap_orig.ident.png", plot = library.averages.heatmap, device = "png", width = 11, height = 5, units = "in")
write.csv(cluster.averages@assays[["SCT"]]@counts, file = "orig.ident_DDX3Y_meancounts.csv")



##### Perform differential expression between clusters using the raw data
DefaultAssay(marsh.seurat.object) <- "RNA"
marsh.seurat.object <- NormalizeData(marsh.seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
marsh.seurat.object <- ScaleData(marsh.seurat.object, features = marsh.seurat.object.var, assay="RNA")

Idents(marsh.seurat.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(marsh.seurat.object, features = marsh.seurat.object.var, assay = "RNA", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "marsh_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>2)
		## all clusters have sig. genes
		## only clusters 3, 5, 6, 7, 8, and 10 have significant genes passing the filter
			top5.heatmap <- DoHeatmap(marsh.seurat.object, features = top5$gene, raster = FALSE)
ggsave("marsh_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(marsh.seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("marsh_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(marsh.seurat.object, features = unique.top2)
png(file=paste0("marsh_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(marsh.seurat.object, features = unique.top2)
png(file=paste0("marsh_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(marsh.seurat.object, features = unique.top2)
png(file=paste0("marsh_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
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

## further analyze immune cell clusters
dir.create("immunology.clusters")
setwd("immunology.clusters")
Idents(marsh.seurat.object) <- "seurat_clusters"
## use RNA for all DE analysis and plots
DefaultAssay(marsh.seurat.object) <- "SCT"
# marsh.seurat.object <- ScaleData(marsh.seurat.object, features = all.genes, assay = "SCT")
marsh.seurat.object <- ScaleData(marsh.seurat.object, features = all.genes, assay="SCT")

## kernel density estimation 
# https://www.nature.com/articles/s41591-021-01329-2#Sec8
major.t <- c("CD4", "CD8A", "CCR7", "PTPRC", "SELL", "CD27", "CD38", "CD44", "CXCR5", "CD40LG", "CCR7", "FOXP3", "IKZF2")
gd.t<- c("TRGV9", "TRDV2")
MAIT.t <- c("TRAV10", "TRAJ18", "TRAV1-2")
NK.t <- c("CD3G", "FCGR3A", "NCAM1", "NCR1")
TH1.t <- c("IFNG", "TBX21", "TNF")
TH2.t <- c("GATA3", "IL4")
TH17.t <- c("RORC", "IL17A", "IL17F", "IL21")
## https://www.bioconductor.org/packages/release/bioc/vignettes/Nebulosa/inst/doc/nebulosa_seurat.html#1_Overview
  ## https://academic.oup.com/bioinformatics/article-abstract/37/16/2485/6103785
# BiocManager::install("Nebulosa")
library("Nebulosa")
Idents(marsh.seurat.object) <- "seurat_clusters"

# SARS.genes <- c("ORF1AB-1-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
p1 <- plot_density(marsh.seurat.object, "CD4")
p2 <- FeaturePlot(marsh.seurat.object, "CD8A")
p5 <- CombinePlots(plots = list(p1, p2))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("UMAP_KernelDensity_CD4-CD8A.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots
p3 <- plot_density(marsh.seurat.object, c("CD4", "CD8A"))
p5 <- p3 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_CD4-CD8A_Joint.pdf", plot = p5, device = "pdf", width = 5, height = 7, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(marsh.seurat.object, c("CD4", "CD8A"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_major.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(marsh.seurat.object, c("TRGV9", "TRDV2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_gd.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(marsh.seurat.object, c("TRAV10", "TRAJ18"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MAIT.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(marsh.seurat.object, c("CD3G", "FCGR3A","NCAM1", "NCR1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_NK.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(marsh.seurat.object, c("IFNG", "TBX21", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH1.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 16, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(marsh.seurat.object, c("GATA3", "IL4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH2.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(marsh.seurat.object, c("RORC", "IL17A", "IL17F", "IL21"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH17.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")

## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(marsh.seurat.object, features = lymphoid.lineage)
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
feature.plot <- DotPlot(marsh.seurat.object, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(marsh.seurat.object, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

macrophage.lineage <- c("CD14", "CD68", "CSF1R", "ADGRE1", "SIGLEC1", "CD38", "GPR18", "FPR2", "MYC", "EGR2")
feature.plot <- DotPlot(marsh.seurat.object, features = macrophage.lineage)
png(file=paste0("macrophage.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Macrophage Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()


# 65. Aagaard-Tillery KM, Silver R, Dalton J. Immunology of normal pregnancy. Semin Fetal Neonatal Med. Oct 2006;11(5):279-95. doi:10.1016/j.siny.2006.04.003
# 66. Ivashkiv LB. Epigenetic regulation of macrophage polarization and function. Trends Immunol. May 2013;34(5):216-23. doi:10.1016/j.it.2012.11.001
# 67. Murray PJ. Macrophage Polarization. Annu Rev Physiol. 02 10 2017;79:541-566. doi:10.1146/annurev-physiol-022516-034339
# 68. Yao Y, Xu XH, Jin L. Macrophage Polarization in Physiological and Pathological Pregnancy. Front Immunol. 2019;10:792. doi:10.3389/fimmu.2019.00792
# 69. Mues B, Langer D, Zwadlo G, Sorg C. Phenotypic characterization of macrophages in human term placenta. Immunology. Jul 1989;67(3):303-7. 
# 70. Bulmer JN, Johnson PM. Macrophage populations in the human placenta and amniochorion. Clin Exp Immunol. Aug 1984;57(2):393-403. 
# 71. Loegl J, Hiden U, Nussbaumer E, et al. Hofbauer cells of M2a, M2b and M2c polarization may regulate feto-placental angiogenesis. Reproduction. 2016;152(5):447-455. doi:10.1530/REP-16-0159
# 74. Schliefsteiner C, Ibesich S, Wadsack C. Placental Hofbauer Cell Polarization Resists Inflammatory Cues In Vitro. Int J Mol Sci. Jan 22 2020;21(3)doi:10.3390/ijms21030736
# 105.  Ben Amara A, Gorvel L, Baulan K, et al. Placental macrophages are impaired in chorioamnionitis, an infectious pathology of the placenta. J Immunol. Dec 01 2013;191(11):5501-14. doi:10.4049/jimmunol.1300988

macrophage.lineage.markers.2 <- c("CD14", "ITGAM", "CSF1R", "ADGRE1", "CD80", "CD38", "GPR18", "FPR2", "CD86", "TNF", "IL12A", "MYC", "EGR2", "CD163", "MRC1", "CD209", "TGFB1", "IL10", "VEGFA", "CD163", "HLA-DRA", "CD86", "IL6")
macrophage.lineage.markers.2 <- unique(macrophage.lineage.markers.2)

feature.plot <- DotPlot(marsh.seurat.object, features = macrophage.lineage.markers.2)
png(file=paste0("macrophage.lineage.markers2-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Macrophage Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

# https://www.nature.com/articles/s41591-021-01329-2#Sec8
t.lineages <- c("CD4", "CD8A", "CCR7", "PTPRC", "SELL", "CD27", "CD38", "CD44", "CXCR5", "CD40LG", "CCR7", "FOXP3", "IKZF2", "TRGV9", "TRDV2", "TRAV10", "TRAJ18", "TRAV1-2", "CD3G", "FCGR3A", "NCAM1", "NCR1", "IFNG", "TBX21", "TNF", "GATA3", "IL4", "RORC", "IL17A", "IL17F", "IL21")
t.lineages <- unique(t.lineages)

# https://www.nature.com/articles/s41591-021-01329-2#Sec8
major.t <- c("CD4", "CD8A", "CCR7", "PTPRC", "SELL", "CD27", "CD38", "CD44", "CXCR5", "CD40LG", "CCR7", "FOXP3", "IKZF2")
gd.t<- c("TRGV9", "TRDV2")
MAIT.t <- c("TRAV10", "TRAJ18", "TRAV1-2")
NK.t <- c("CD3G", "FCGR3A", "NCAM1", "NCR1")
TH1.t <- c("IFNG", "TBX21", "TNF")
TH2.t <- c("GATA3", "IL4")
TH17.t <- c("RORC", "IL17A", "IL17F", "IL21")

#marsh.seurat.object <- ScaleData(marsh.seurat.object, features = all.genes, assay = "SCT")
cluster.averages <- AverageExpression(marsh.seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=t.lineages)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, slot="counts", assay="SCT")
ggsave("t.lineages_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, assay="SCT")
ggsave("t.lineages_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(marsh.seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=macrophage.lineage.markers.2)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, slot="counts", assay="SCT")
ggsave("macrophage.lineage.markers.2_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, assay="SCT")
ggsave("macrophage.lineage.markers.2_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(marsh.seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=lymphoid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, slot="counts", assay="SCT")
ggsave("lymphoid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, assay="SCT")
ggsave("lymphoid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(marsh.seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=myeloid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, slot="counts", assay="SCT")
ggsave("myeloid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, assay="SCT")
ggsave("myeloid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


############## Activation gating strategy https://www.nature.com/articles/s41590-021-01049-2/figures/10
## 2 Density Plots and a Joint density plot
p4 <- plot_density(marsh.seurat.object, c("CD3D", "NCAM1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TorNK_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(marsh.seurat.object, c("CD4", "KLRB1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRB1-CD161.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(marsh.seurat.object, c("CD4", "CD38"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD38.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(marsh.seurat.object, c("CD4", "CD69"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD69.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("CD4", "KLRK1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRK1-NKG2D.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("CD4", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("FOXP3", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-FOXP3-CD4.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("KIT", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("TNF", "TNFRSF8"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.TNFRSF8.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p4 <- plot_density(marsh.seurat.object, c("CD14", "ITGAM"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("SIGLEC1", "FOLR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("SIGLEC1", "DDX3Y"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC-fetalsex_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
	## Placenta-associated maternal macrophages PAMMS (FOLR2 neg, HLA positive); HBCs FOLR2+HLA- ; Thomas et al. Naomi McGovern, JEM 2020
p4 <- plot_density(marsh.seurat.object, c("FOLR2", "HLA-DRA"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMMs_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM1_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM2_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("ITGAM", "CD14"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M0.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("CD80", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("CD11B", "MYC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("CD163", "MRC1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2a.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("HLA-DRA", "IL6"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("CD14", "CD163"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2c.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("CD14", "CD1C"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(marsh.seurat.object, c("CD14", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p1 <- plot_density(marsh.seurat.object, "THBD")
p2 <- FeaturePlot(marsh.seurat.object, "THBD")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(marsh.seurat.object, "ITGB3")
p2 <- FeaturePlot(marsh.seurat.object, "ITGB3")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_Megakaryocyte.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(marsh.seurat.object, "CD1C")
p2 <- FeaturePlot(marsh.seurat.object, "CD1C")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(marsh.seurat.object, "IL6")
p2 <- FeaturePlot(marsh.seurat.object, "IL6")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(marsh.seurat.object, "TNF")
p2 <- FeaturePlot(marsh.seurat.object, "TNF")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(marsh.seurat.object, "MYC")
p2 <- FeaturePlot(marsh.seurat.object, "MYC")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(marsh.seurat.object, "FOLR2")
p2 <- FeaturePlot(marsh.seurat.object, "FOLR2")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_HBCs-FOLR2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(marsh.seurat.object, "HLA-DRA")
p2 <- FeaturePlot(marsh.seurat.object, "HLA-DRA")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_PAMM-HLA-DRA.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(marsh.seurat.object, "CD33")
p2 <- FeaturePlot(marsh.seurat.object, "CD33")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeliodProgenitor-M0.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p4 <- plot_density(marsh.seurat.object, c("KIT", "CD33"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_MyeliodProgenitor-M0-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
# marsh.seurat.object <- DietSeurat(marsh.seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

## note HLA genes will not work. Consider finding homologous gene to HLA-DR;; consider H2-XX eg H2-Q2 genes https://pubmed.ncbi.nlm.nih.gov/14602227/  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC17889/ 
setwd("..")

rm(cluster.averages)
rm(library.averages.heatmap)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(p5)
rm(p6)
rm(p7)
rm(p8)
rm(plot2)
rm(umap.combined)
rm(feature.plot)
rm(de_markers)

dir.create("DE_clusters-all.genes")
setwd("DE_clusters-all.genes")
	## clusters 0, 2, 9, 14, and 15 did not have any genes that passed the previous threshold. 
		## Consider all.genes instead of the top 3k variable and see what genes define them. If still no genes, see if it is a QC metric or Phase
##### Perform differential expression between clusters using the raw data

marsh.all.genes <- rownames(marsh.seurat.object)

DefaultAssay(marsh.seurat.object) <- "RNA"
marsh.seurat.object <- NormalizeData(marsh.seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
marsh.seurat.object <- ScaleData(marsh.seurat.object, features = marsh.all.genes, assay="RNA")

Idents(marsh.seurat.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(marsh.seurat.object, features = marsh.all.genes, assay = "RNA", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "all.genes_DEGs_byclusters_pos-log2fc.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>2)
		## all clusters have sig. genes
		## only clusters 3, 5, 6, 7, 8, and 10 have significant genes passing the filter
			top5.heatmap <- DoHeatmap(marsh.seurat.object, features = top5$gene, raster = FALSE)
ggsave("integrated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(marsh.seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("integrated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(marsh.seurat.object, features = unique.top2)
png(file=paste0("integrated_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(marsh.seurat.object, features = unique.top2)
png(file=paste0("integrated_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(marsh.seurat.object, features = unique.top2)
png(file=paste0("integrated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
setwd("..")
############################################################# ############################################################# #############################################################
## revise Type annotations
############################################################# ############################################################# #############################################################
## see DE results from Marsh paper: elife-60266-fig1-data1-v2.csv https://elifesciences.org/articles/60266 for annotations
## what are clusters 0, 2, and 9; DE reveals no markers, but CC phase suggests proliferative
## subset and do DE analysis on them directly to find differences
dir.create("DE_prolif.clusters")
setwd("DE_prolif.clusters")
Idents(marsh.seurat.object)<-"seurat_clusters"
marsh.seurat.object$Marsh.Cluster <- Idents(marsh.seurat.object)
Idents(marsh.seurat.object) <- "Marsh.Cluster"
marsh.seurat.object.prolif <- subset(marsh.seurat.object, idents = c("0","2", "9"))
marsh.seurat.object.prolif <- DietSeurat(marsh.seurat.object.prolif, counts=TRUE, data=TRUE, scale.data = FALSE)
marsh.seurat.object.prolif <- SCTransform(marsh.seurat.object.prolif)
marsh.seurat.object.prolif <- RunPCA(marsh.seurat.object.prolif, assay="SCT", verbose = TRUE)  
var.genes <- marsh.seurat.object.prolif@assays[["SCT"]]@var.features
marsh.seurat.object.prolif <- FindNeighbors(marsh.seurat.object.prolif, assay="SCT",  features = "var.genes", dims = 1:30)
marsh.seurat.object.prolif <- FindClusters(marsh.seurat.object.prolif, assay="SCT", resolution = 0.6)
marsh.seurat.object.prolif <- RunUMAP(marsh.seurat.object.prolif, dims = 1:30)
p1 <- DimPlot(marsh.seurat.object.prolif, reduction = "umap", group.by = "Type.Profile", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
p1 
p4 <- DimPlot(marsh.seurat.object.prolif, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Prolif.Cluster')
p4 
p2 <- DimPlot(marsh.seurat.object.prolif, reduction = "umap", group.by = "Marsh.Cluster", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Marsh.Cluster')
p2 
marsh.all.genes <- rownames(marsh.seurat.object.prolif)

umap.combined <- p1+p4+p2
ggsave("Marsh_UMAP-prolif.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")


DefaultAssay(marsh.seurat.object.prolif) <- "RNA"
marsh.seurat.object.prolif <- NormalizeData(marsh.seurat.object.prolif, normalization.method = "LogNormalize", scale.factor = 10000)
marsh.seurat.object.prolif <- ScaleData(marsh.seurat.object.prolif, features = marsh.all.genes, assay="RNA")

Idents(marsh.seurat.object.prolif) <- "seurat_clusters"
de_markers <- FindAllMarkers(marsh.seurat.object.prolif, features = marsh.all.genes, assay = "RNA", only.pos = F, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "prolif_all.genes_DEGs_byclusters_pos-log2fc.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(marsh.seurat.object.prolif, features = unique.top2)
png(file=paste0("prolif_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(marsh.seurat.object.prolif, features = unique.top2)
png(file=paste0("prolif_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(marsh.seurat.object.prolif)
setwd("..")

# change to VCT and SYT
Idents(marsh.seurat.object) <- "seurat_clusters"
levels(marsh.seurat.object) 
 ##  [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18"
ncol(marsh.seurat.object)	# 9743 spots
new.metadata <- c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Trophoblast", "FetalMesenchyme", "Prolif_Fibro", "Trophoblast", "DecidualStroma", "Immune", "Trophoblast", "DecidualStroma", "Trophoblast", "Trophoblast", "FetalMesenchyme", "Immune", "Endothelial")
names(new.metadata) <- levels(marsh.seurat.object)
marsh.seurat.object <- RenameIdents(marsh.seurat.object, new.metadata)
marsh.seurat.object$Type <- Idents(marsh.seurat.object)
Idents(marsh.seurat.object) <- "Type.Profile"
# [1] "Prolif_Fibro"         "Trophoblast"          "Prolif_Stromal"       "Endothelial"         
# [5] "FetalMesenchyme"      "Endometrial_Stromal " "Erythrocyte"          "Immune"              
# [9] "DecidualStroma"
## note most annotations are from Marsh annotations in Marsh Fig 1 w/ DE genes found here


## revise Type annotations
# change to VCT and SYT
Idents(marsh.seurat.object) <- "seurat_clusters"
levels(marsh.seurat.object) 
 ## 0 thru 18 clusters
ncol(marsh.seurat.object)	# 9743 spots
new.metadata <- c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2")
names(new.metadata) <- levels(marsh.seurat.object)
marsh.seurat.object <- RenameIdents(marsh.seurat.object, new.metadata)
marsh.seurat.object$Subtype <- Idents(marsh.seurat.object)
Idents(marsh.seurat.object) <- "Subtype.Profile"
levels(marsh.seurat.object) 
#  [1] "Prolif_Fibro"      "SynTI"             "Prolif_Stromal"    "Endothelial_1"     "FetalMesenchyme_1"
#  [6] "GC"                "Erythrocyte"       "SpT_1"             "FetalMesenchyme_2" "SynTII"           
# [11] "Immune_1"          "S-TGC"             "SpT_2"             "SYT"               "VCT"              
# [16] "FetalMesenchyme_3" "Immune_2"          "Endothelial_2"


## note subtype annotations mostly from Marsh trohpoblast subset analysis
Prolif= Proliferative, defined by G2M/S cell cycle phase assignments
VCT= Villous cytotrophoblast
SYT= Syncytiotrophoblast
SpT= spongiotrophoblasts w/n junctional zone (Marsh et al., 2020)
S-TGC=Sinusoidal trophoblast giant cells (Simmons,et al. 2008b; Marsh et al., 2020)
SynTI=outer most SYT layer, w/ maternal blood contact (Georgiades et al., 2002; Woods et al., 2018; Marsh et al., 2020)
SynTII= beneath outer most SYT layer, intermediate between SYT layer and fetal endothelium (Georgiades et al., 2002; Woods et al., 2018; Marsh et al., 2020)
GC= Glycogen cell, resides w/n junctional zone (Marsh et al., 2020)

setwd("..")
p4 <- DimPlot(marsh.seurat.object, reduction = "umap", group.by = "Type.Profile", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
p4 		
p2 <- DimPlot(marsh.seurat.object, reduction = "umap", group.by = "Phase", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Phase')
p2  
p3 <- DimPlot(marsh.seurat.object, reduction = "umap", group.by = "Subtype.Profile", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Subtype')
p3 
p1 <- DimPlot(marsh.seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p1 
#examine UMAPs with qc metrics 

umap.combined <- p1+p4+p2+p3
ggsave("Marsh_UMAP-final.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("Marsh_XL-final.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")

############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################

## Make sure visium data has SCTransformed data with PCA
seurat.object <- SCTransform(seurat.object, ncells = 3000, conserve.memory = TRUE, assay="Spatial", verbose = TRUE) 
seurat.object <- RunPCA(seurat.object, assay="SCT", verbose = TRUE)  

DefaultAssay(marsh.seurat.object) <- "SCT"
DefaultAssay(seurat.object) <- "SCT"
marsh.seurat.object <- UpdateSeuratObject(marsh.seurat.object) # Since your reference is generated in Seurat3, you just need to update your object. https://github.com/satijalab/seurat/issues/5422
marsh.seurat.object <- UpdateSCTAssays(marsh.seurat.object)
seurat.object <- UpdateSeuratObject(seurat.object) # Since your reference is generated in Seurat3, you just need to update your object. https://github.com/satijalab/seurat/issues/5422
seurat.object <- UpdateSCTAssays(seurat.object)
anchors <- FindTransferAnchors(reference = marsh.seurat.object, query = seurat.object, normalization.method = "SCT", verbose=T)
  ## Error in slot(object = reference[[reference.assay]], name = "SCTModel.list") : 
  	# no slot of name "SCTModel.list" for this object of class "Assay"
  		## now addressed above
predictions.assay <- TransferData(anchorset = anchors, refdata = marsh.seurat.object$Type, prediction.assay = TRUE, weight.reduction = seurat.object[["pca"]], dims = 1:30)

seurat.object[["Predictions"]] <- predictions.assay

DefaultAssay(seurat.object) <- "Predictions"
levels(seurat.object)
Idents(seurat.object) <- "orig.ident"
levels(seurat.object)

# Now we get prediction scores for each spot for each class. Of particular interest in the placenta are the macrophages and syncytial trophoblasts (SYTs). 
# Here we can distinguish between distinct regions for these subtypes, for example:
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
image.list <- c("S07.slice.3")

setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/atlas/predictions_v1")

  ##  c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Trophoblast", "FetalMesenchyme", "Prolif_Fibro", "Trophoblast",
  #  "Immune", "Trophoblast", "DecidualStroma", "Trophoblast", "Trophoblast", "FetalMesenchyme", "Immune", "Endothelial")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S07.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S07.slice.1")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S09.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S12.slice.1")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S12.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S13.slice.3")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S13.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S31.slice.4")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S31.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S32.slice.5")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S32.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S33.slice.6")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S33.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S33.slice.6")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S34.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

## Try plotting all 
	## for subtypes later c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2")
image.list <- c("S07.slice.3")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S07.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.1")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S09.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S12.slice.1")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S12.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S13.slice.3")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S13.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S31.slice.4")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S31.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S32.slice.5")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S32.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S33.slice.6")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S33.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S33.slice.6")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S34.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")


p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S34.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")


p4 <- FeaturePlot(seurat.object, reduction = "umap", features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), raster=T) + labs(color='Predictions')
p4 		
	# The following requested variables were not found: Prolif_Fibro, Prolif_Stromal, Endometrial_Stromal
ggsave("Visium-Marsh-Type-Predictions_UMAP.pdf", plot = p4, device = "pdf", width = 10, height = 12, units = "in")
############### Use this for cluster annotation
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
	
	## Try again with subtypes for Predictions


setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/atlas")
dir.create("predictions-Subtype_v1")
setwd("predictions-Subtype_v1")

## Make sure visium data has SCTransformed data with PCA
DefaultAssay(marsh.seurat.object) <- "SCT"
DefaultAssay(seurat.object) <- "SCT"
marsh.seurat.object <- UpdateSeuratObject(marsh.seurat.object) # Since your reference is generated in Seurat3, you just need to update your object. https://github.com/satijalab/seurat/issues/5422
marsh.seurat.object <- UpdateSCTAssays(marsh.seurat.object)
seurat.object <- UpdateSeuratObject(seurat.object) # Since your reference is generated in Seurat3, you just need to update your object. https://github.com/satijalab/seurat/issues/5422
seurat.object <- UpdateSCTAssays(seurat.object)
anchors <- FindTransferAnchors(reference = marsh.seurat.object, query = seurat.object, normalization.method = "SCT", verbose=T)
  ## Error in slot(object = reference[[reference.assay]], name = "SCTModel.list") : 
  	# no slot of name "SCTModel.list" for this object of class "Assay"
  		## now addressed above
predictions.assay.subtype <- TransferData(anchorset = anchors, refdata = marsh.seurat.object$Subtype, prediction.assay = TRUE, weight.reduction = seurat.object[["pca"]], dims = 1:30)

seurat.object[["Predictions.subtype"]] <- predictions.assay.subtype

DefaultAssay(seurat.object) <- "Predictions.subtype"
levels(seurat.object)
Idents(seurat.object) <- "orig.ident"
levels(seurat.object)

p4 <- FeaturePlot(seurat.object, reduction = "umap", features = c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2"), raster=T) + labs(color='Predictions')
p4 		
	# The following requested variables were not found: Prolif_Fibro, Prolif_Stromal, Endometrial_Stromal
ggsave("Visium-Marsh-Subtype-Predictions_UMAP.pdf", plot = p4, device = "pdf", width = 10, height = 12, units = "in")

# Now we get prediction scores for each spot for each class. Of particular interest in the placenta are the macrophages and syncytial trophoblasts (SYTs). 
# Here we can distinguish between distinct regions for these subSubtypes, for example:
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
image.list <- c("S07.slice.3")

## for subtypes  c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", 
  #  "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2")

p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("SYT", "SpT_1","SpT_2", "GC", "Immune_1", "Immune_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S07.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S07.slice.1")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("SYT", "SpT_1","SpT_2", "GC", "Immune_1", "Immune_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S09.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S12.slice.1")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("SYT", "SpT_1","SpT_2", "GC", "Immune_1", "Immune_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S12.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S13.slice.3")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("SYT", "SpT_1","SpT_2", "GC", "Immune_1", "Immune_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S13.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S31.slice.4")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("SYT", "SpT_1","SpT_2", "GC", "Immune_1", "Immune_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S31.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S32.slice.5")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("SYT", "SpT_1","SpT_2", "GC", "Immune_1", "Immune_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S32.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S33.slice.6")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("SYT", "SpT_1","SpT_2", "GC", "Immune_1", "Immune_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S33.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S33.slice.6")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("SYT", "SpT_1","SpT_2", "GC", "Immune_1", "Immune_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S34.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

## Try plotting all 
	## for subSubtypes later c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2")
image.list <- c("S07.slice.3")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.Subtypes_S07.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.1")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.Subtypes_S09.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S12.slice.1")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.Subtypes_S12.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S13.slice.3")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.Subtypes_S13.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S31.slice.4")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.Subtypes_S31.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S32.slice.5")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.Subtypes_S32.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S33.slice.6")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.Subtypes_S33.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S33.slice.6")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.Subtypes_S34.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")


p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.Subtypes_S34.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")


p4 <- FeaturePlot(seurat.object, reduction = "umap", features = c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2"), raster=T) + labs(color='Predictions')
p4 
# Warning messages:
#1: In FetchData.Seurat(object = object, vars = c(dims, "ident", features),  :
#  The following requested variables were not found (10 out of 11 shown): Prolif_Fibro, Prolif_Stromal, Endothelial_1, FetalMesenchyme_1, SpT_1, FetalMesenchyme_2, Immune_1, SpT_2, FetalMesenchyme_3, Immune_2
#2: In FeaturePlot(seurat.object, reduction = "umap", features = c("Prolif_Fibro",  :
#  All cells have the same value (0) of SYT.
########### only have: SynTI, GC, Erythrocyte, SynTII, S-TGC, and VCT
p4 <- FeaturePlot(seurat.object, reduction = "umap", features = c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2",  "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2"), raster=T) + labs(color='Predictions')
p4 		
	# The following requested variables were not found: Prolif_Fibro, Prolif_Stromal, Endometrial_Stromal
ggsave("Visium-Marsh-Subtype-Predictions_UMAP_final.pdf", plot = p4, device = "pdf", width = 10, height = 12, units = "in")

save.image("predictions_v1.RData")
setwd("..")
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################ ############################################################# #############################################################
############################################################# Annotation of Visium UMAP Clusters #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################

	## Based on the predictions above, overlap with snRNA-seq clusters, and DE markers from Visium UMAP clusters, we will rename the clusters into subtypes and types below
# see /Users/enricobarrozo/Library/CloudStorage/Box-Box/AagaardLab/Visium/SpaceRanger/results/seurat_mouse_v2.1/atlas/predictions_v1/DE_clusters-marsh/annotation-marsh_DEGs_byclusters_pos-log2FC.xlsx

## Annotate the Visium clusters using the Marsh Type and Subtype Predictions
## mostly Visium-Marsh-Subtype-Predictions_UMAP_final.pdf and Visium-Marsh-Type-Predictions_UMAP_final.pdf compared to integrated_UMAP.pdf seurat_clusters
# Visium.Type.Profile c("Trophoblast", "Trophoblast", "Trophoblast", "Trophoblast", "Trophoblast", "DecidualStroma", "Trophoblast", "FetalMesenchyme/Immune/Endothelial", "FetalMesenchyme", "Trophoblast", "Trophoblast", "FetalMesenchyme/Trophoblast")
# Visium.Subtype.Profile c("VCT", "SynTI", "SynTII", "Trophoblast", "GC", "DecidualStroma", "SynTI", "FetalMesenchyme/Immune/Endothelial", "FetalMesenchyme", "SynTII/S-TGC", "S-TGC", "VCT/SynTII/SynTI")
##

setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1")
dir.create("annotated")
setwd("annotated")
## For cluster Subtypes, I took the [integrated_DEGs_byclusters_pos-log2FC.txt] 
		# cluster-Subtypes_DEGs_byclusters_top5k-pos-log2FC.xlsx
	# with adj. p<0.05 and log2FC>2 and upload top marker genes into https://placentacellenrich.gdcb.iastate.edu and annotate based on vento-tormo et al or suryawanshi datasets. 
	## Also examining the human protein atlas (https://www.proteinatlas.org), Vento-Tormo and Suryawanshi datasets (https://placentacellenrich.gdcb.iastate.edu), and the PangloaDB (https://panglaodb.se/search.html)
Idents(seurat.object) <- "seurat_clusters"
cluster.Subtypes<- c("Trophoblast", "Trophoblast", "Trophoblast", "Trophoblast", "Trophoblast", "DecidualStroma", "Trophoblast", "FetalMesenchyme/Immune/Endothelial", "FetalMesenchyme", "Trophoblast", "Trophoblast", "FetalMesenchyme/Trophoblast")
names(cluster.Subtypes) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Subtypes)
## add cluster Subtypes as a metadata column
seurat.object$Type.Profile <- Idents(seurat.object)
levels(seurat.object) # [1] "Trophoblast"                        "DecidualStroma"                     "FetalMesenchyme/Immune/Endothelial"
# [4] "FetalMesenchyme"                    "FetalMesenchyme/Trophoblast"  

Idents(seurat.object) <- "seurat_clusters"
cluster.Type <- c("VCT", "SynTI", "SynTII", "Trophoblast", "GC", "DecidualStroma", "SynTI", "FetalMesenchyme/Immune/Endothelial", "FetalMesenchyme", "SynTII/S-TGC", "S-TGC", "VCT/SynTII/SynTI")
names(cluster.Type) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Type)
## add cluster Types as a metadata column
seurat.object$Subtype.Profile <- Idents(seurat.object)

# setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/annotated")
# seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs=c("umap"))
# save.image("murine_spatial_data-annotated_v1.RData")
# load("murine_spatial_data-annotated_v1.RData")


p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "Treatment", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Treatment')
p1
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p2
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Sample')
p3
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Batch", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Batch')
p4
p5 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type.Profile",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type.Profile')
p5
p6 <- DimPlot(seurat.object, reduction = "umap", group.by = "Subtype.Profile",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Subtype.Profile')
p6
p7 <- DimPlot(seurat.object, reduction = "umap", group.by = "Virus", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Virus')
p7
umap.combined <- p1 + p2 + p5 + p6
ggsave("Visium_UMAP-annotated_1.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("Visium_UMAP-annotated_XL_1.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")
ggsave("Visium_UMAP-annotated_M_1.pdf", plot = umap.combined, device = "pdf", width = 10, height = 6, units = "in")

p9 <- DimPlot(seurat.object, reduction = "umap", group.by = "Phase", label= "TRUE", repel=TRUE, raster=T)

umap.combined <- p7+p1+p3+p9
ggsave("Visium_UMAP-annotated_2.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("Visium_UMAP_annotated_XL_2.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")
ggsave("Visium_UMAP_annotated_M_2.pdf", plot = umap.combined, device = "pdf", width = 10, height = 6, units = "in")

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
umap.combined <- p4+p5+p6+p7
ggsave("Visium_UMAP-annotated_3.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("Visium_UMAP_annotated_XL_3.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")
ggsave("Visium_UMAP-annotated_M_3.pdf", plot = umap.combined, device = "pdf", width = 10, height = 6, units = "in")


Idents(seurat.object) <- "seurat_clusters"
t1<-table(Idents(seurat.object))
write.table(t1, "Visium_seurat_clusters.counts.txt", sep="\t")
rm(t1)


plot2 <- SpatialFeaturePlot(seurat.object, features = "percent.viral", ncol = 1) + patchwork::plot_layout(ncol = 1)
####### FOR SOME REASON, INTEGRATION ADDS IMAGES ARBITRARILY
ggsave("integrated_SpatialDimPlot_percent.viral_TEST.pdf", plot = plot2, device = "pdf", width = 8, height = 45, units = "in")



image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
Idents(seurat.object) <- "Type.Profile"
plot2 <- SpatialDimPlot(seurat.object,images=image.list, alpha=1, ncol=2)
plot2
ggsave("SpatialDimPlot_annotated-Type.Profile_all.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2 <- SpatialDimPlot(seurat.object, images=image.list, alpha=0.1, ncol=2) 
plot2
ggsave("SpatialDimPlot_annotated-Type.Profile_all-alpha.1.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
Idents(seurat.object) <- "Subtype.Profile"
plot2 <- SpatialDimPlot(seurat.object,images=image.list, alpha=1, ncol=2) 
ggsave("SpatialDimPlot_annotated-Subtype.Profile_all.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2 <- SpatialDimPlot(seurat.object, images=image.list, alpha=0.1, ncol=2) 
ggsave("SpatialDimPlot_annotated-Subtype.Profile_all-alpha.1.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")


image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
Idents(seurat.object) <- "Type.Profile"
plot2 <- SpatialDimPlot(seurat.object,images=image.list, alpha=1, ncol=2, crop=F) 
plot2
ggsave("SpatialDimPlot_annotated-Type.Profile_all-uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2 <- SpatialDimPlot(seurat.object, images=image.list, alpha=0.1, ncol=2, crop=F) 
ggsave("SpatialDimPlot_annotated-Type.Profile_all-alpha.1-uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
Idents(seurat.object) <- "Subtype.Profile"
plot2 <- SpatialDimPlot(seurat.object,images=image.list, alpha=1, ncol=2, crop=F) 
ggsave("SpatialDimPlot_annotated-Subtype.Profile_all-uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2 <- SpatialDimPlot(seurat.object, images=image.list, alpha=0.1, ncol=2, crop=F) 
ggsave("SpatialDimPlot_annotated-Subtype.Profile_all-alpha.1-uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")



dir.create("DE_Type.Profile")
setwd("DE_Type.Profile")
##### Perform differential expression between Types using the raw data
seurat.object <- SCTransform(seurat.object, assay = "Spatial", verbose = T)
DefaultAssay(seurat.object) <- "SCT"
seurat.object.var <- seurat.object@assays[["SCT"]]@var.features
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")

Idents(seurat.object) <- "Type.Profile"
de_markers <- FindAllMarkers(seurat.object, features = all.genes.visium, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "Type_annotated_DEGs_byclusters_pos-log2FC.txt", sep="\t")
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
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type.Profile")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Type_annotated_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type.Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Type_annotated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type.Profile")
dev.off()

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes.visium)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("unique.top.markers_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("unique.top.markers_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

setwd("..")

dir.create("DE_Subtype.Profile")
setwd("DE_Subtype.Profile")
##### Perform differential expression between Subtypes using the raw data
seurat.object <- SCTransform(seurat.object, assay = "Spatial", verbose = T)
DefaultAssay(seurat.object) <- "SCT"
seurat.object.var <- seurat.object@assays[["SCT"]]@var.features
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")

Idents(seurat.object) <- "Subtype.Profile"
de_markers <- FindAllMarkers(seurat.object, features = all.genes.visium, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "Subtype_annotated_DEGs_byclusters_pos-log2FC.txt", sep="\t")
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
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Subtype.Profile")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Subtype_annotated_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Subtype.Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Subtype_annotated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Subtype.Profile")
dev.off()

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes.visium)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("unique.top.markers_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("unique.top.markers_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

setwd("..")

dir.create("DE_orig.ident")
setwd("DE_orig.ident")
Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident")
ggsave("integrated_UMAP-orig.ident.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")

Idents(seurat.object) <- "orig.ident"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-log2FC.txt", sep="\t")
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
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="SampleID")
dev.off()

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes.visium)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("unique.top.markers_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("unique.top.markers_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

setwd("..")


dir.create("DE_Treatment")
setwd("DE_Treatment")
Idents(seurat.object) <- "Treatment"
levels(seurat.object)
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "Treatment")
ggsave("integrated_UMAP-Treatment.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")

Idents(seurat.object) <- "Treatment"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "DEGs_byTreatment_pos-log2FC.txt", sep="\t")
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
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Treatment")
dev.off()

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes.visium)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("unique.top.markers_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("unique.top.markers_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

setwd("..")

setwd("..")
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
	## FindSpatiallyVariableFeatures 
setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/annotated")
dir.create("DE_SpatiallyVariableFeatures_Treatment")
setwd("DE_SpatiallyVariableFeatures_Treatment")
seurat.object <- SCTransform(seurat.object, ncells = 3000, conserve.memory = TRUE, assay="Spatial", verbose = TRUE) 
############# DE based on predictions
# image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")? 
image.list <- c("S07.slice.3")
Idents(seurat.object) <- "Treatment"
all.genes <- rownames(seurat.object)
DefaultAssay(seurat.object) <- "SCT"
seurat.object <- FindSpatiallyVariableFeatures(seurat.object, assay = "SCT", features = VariableFeatures(seurat.object)[1:1000],
    selection.method = "markvariogram")
seurat.object <- FindSpatiallyVariableFeatures(seurat.object, assay = "SCT", selection.method = "markvariogram", features = VariableFeatures(seurat.object)[1:1000], only.pos=T, verbose=T)
top.clusters <- head(SpatiallyVariableFeatures(seurat.object), 4)
  # "SmoothMuscle" "SYT"          "VCT"          "Fibroblast" 
p1 <- SpatialPlot(object = seurat.object, features = top.clusters, ncol = 2, images = image.list)
p1
ggsave("Predictions_top.clusters_S07.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

top.clusters <- SpatiallyVariableFeatures(seurat.object)
 # [1] "SmoothMuscle"           "SYT"                    "VCT"                   
 #[4] "Fibroblast"             "DC2"                    "B"                     
 #[7] "RBC"                    "max"                    "VEC"                   
#[10] "EVT"                    "Stromal"                "PVC"                   
# [13] "Macrophage"             "Endometrial-epithelium" "DC"

# rm max
top.clusters.1 <- top.clusters[-8]

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
image.list <- c("S07.slice.3")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S07.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S07.slice.1")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S09.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S12.slice.1")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S12.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S13.slice.3")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S13.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")

image.list <- c("S31.slice.4")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S31.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S32.slice.5")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S32.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S33.slice.6")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S33.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S34.slice.7")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S34.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

rownames((predictions.assay))
# [1] "EVT"                    "Stromal"                "CD8T"                  
# [4] "NKT                     "Macrophage"             "VCT"                   
 #[7] "Fibroblast"             "B"                      "SmoothMuscle"          
#[10] "DC2"                    "SYT"                    "Megakaryocyte"         
#[13] "VEC"                    "PVC"                    "Endometrial-epithelium"
#[16] "DC"                     "RBC"                    "T"                     
#[19] "Endothelial"            "max" 

Idents(seurat.object) <- "Predictions"
levels(seurat.object)

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
plot2<- SpatialDimPlot(seurat.object, images=image.list) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_cropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialDimPlot(seurat.object, images=image.list, crop=FALSE) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
 

# https://learn.gencore.bio.nyu.edu/seurat-integration-and-label-transfer/
seurat.object <- AddMetaData(object = seurat.object, metadata = Predictions)
seurat.object$prediction.match <- seurat.object$predicted.id == seurat.object$Type
table(seurat.object$prediction.match)
table(seurat.object$predicted.id)

Idents(seurat.object) <- "prediction.match"
levels(seurat.object) # "Mock-Enox-1" "Mock-Enox-2" "ZIKV-Veh-1"  "ZIKV-Enox-1" "ZIKV-Enox-2" "Mock-Veh-1"  "ZIKV-Veh-2" "ZIKV-Veh-3"


Idents(seurat.object) <- "predicted.id"
levels(seurat.object) # "Mock-Enox-1" "Mock-Enox-2" "ZIKV-Veh-1"  "ZIKV-Enox-1" "ZIKV-Enox-2" "Mock-Veh-1"  "ZIKV-Veh-2" "ZIKV-Veh-3"


ncol(seurat.object)	# 9743 spots
new.metadata <- c("Mock-Enox", "Mock-Enox", "ZIKV-Veh", "ZIKV-Enox", "ZIKV-Enox", "Mock-Veh", "ZIKV-Veh", "ZIKV-Veh")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Treatment <- Idents(seurat.object)
Idents(seurat.object) <- "Treatment"


# save.image("predictions_v1.RData")
  # load('predictions_v1.RData')  ## 249 gb w/o removing atlas and predicitons.assay
setwd("..")
############################################################# ############################################################# #############################################################
############################################################# V. Analysis of visium and snRNA-seq data together (mouse placenta atlas) #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# Atlas Integration: merge all data using CCA anchors #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/annotated")
# load("murine_spatial_data-annotated_v1.RData")

all.genes.marsh <- rownames(marsh.seurat.object)

## combine lists of genes
all.genes.combined <- union(all.genes, all.genes.marsh)
	## 32856 combined genes

## combine lists of top variable features for later DE analysis/clustering
var.combined <- union(var.combined, marsh.seurat.object.var)
	## 10367 combined variable features


Idents(marsh.seurat.object) <- "Type"
levels(marsh.seurat.object) # [1] "Prolif_Fibro"         "Trophoblast"          "Prolif_Stromal"       "Endothelial"         [5] "FetalMesenchyme"    
#   "Endometrial_Stromal " "Erythrocyte"          "Immune"              [9] "DecidualStroma"
#names(new.metadata) <- levels(marsh.seurat.object)
#marsh.seurat.object <- RenameIdents(marsh.seurat.object, new.metadata)
marsh.seurat.object$Type.Profile <- Idents(marsh.seurat.object)
Idents(marsh.seurat.object) <- "Type.Profile"
Idents(marsh.seurat.object) <- "Subtype"
marsh.seurat.object$Subtype.Profile <- Idents(marsh.seurat.object)
Idents(marsh.seurat.object) <- "Subtype.Profile"

DefaultAssay(seurat.object) <- "SCT"
DefaultAssay(marsh.seurat.object) <- "SCT"
Idents(seurat.object) <- "Assay"
Idents(marsh.seurat.object) <- "Assay"
seurat.object <- SCTransform(seurat.object, assay = "Spatial", features=var.combined, verbose = T)
marsh.seurat.object <- SCTransform(marsh.seurat.object, assay = "RNA", features=var.combined, verbose = T)
seurat.object<- RunPCA(seurat.object, assay = "SCT", slot = "scale.data")
marsh.seurat.object<- RunPCA(marsh.seurat.object, assay = "SCT", slot = "scale.data")
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = T, dimreducs=c("umap", "pca"))
marsh.seurat.object <- DietSeurat(marsh.seurat.object, counts=TRUE, data=TRUE, scale.data = T, dimreducs=c("umap", "pca"))


## Merge objects
all_merged <- merge(x = seurat.object, y = c(marsh.seurat.object), merge.data = TRUE, project = "all_merged")

ncol(all_merged)
	# 9701 spots/cells
options(future.globals.maxSize = 300000000000)
 
all_merged <- SCTransform(all_merged, assay = "SCT", verbose = T)
all_merged<- RunPCA(all_merged, assay = "SCT", slot = "scale.data")
DefaultAssay(all_merged) <- "SCT"


all_merged <- DietSeurat(all_merged, counts=TRUE, data=TRUE, scale.data = T, dimreducs=c("pca"))
Idents(all_merged) <- "Type.Profile"
levels(all_merged)

all_merged.list <- SplitObject(all_merged, split.by = "Type.Profile")
levels(all_merged.list)
# Make a reference list

reference.list <- all_merged.list[c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma", "FetalMesenchyme/Trophoblast", "FetalMesenchyme/Immune/Endothelial")]

reference.features <- SelectIntegrationFeatures(object.list = reference.list)
write.table(reference.features, "all_merged.reference.features.txt", sep="\t")
reference.list <- PrepSCTIntegration(object.list = reference.list, anchor.features = reference.features, verbose = T)
all_merged.anchors <- FindIntegrationAnchors(object.list = reference.list, normalization.method = "SCT", reduction = "rpca", dims = 1:50, anchor.features = reference.features, verbose = T)
atlas.object <- IntegrateData(anchorset = all_merged.anchors, dims = 1:50)
rm(reference.list)
rm(all_merged.anchors)

# Use the 'integrated' assay for clustering and the SCT assay for differential expression
DefaultAssay(atlas.object) <- "integrated"
# Add all integrated genes as a list of features to call upon during clustering
integrated.genes <- rownames(atlas.object)
write.table(integrated.genes, "integrated.genes.txt", sep="\t")

atlas.object <- DietSeurat(atlas.object, counts=TRUE, data=TRUE, scale.data = FALSE)

rm(all_merged)
rm(all_merged.anchors)
rm(all_merged.list)
rm(reference.list)
rm(marsh.seurat.object)
rm(seurat.object)
all.genes.visium <- all.genes
all.genes <- all.genes.combined
setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/atlas")
# save.image("murine_atlas_data-integrated_v1.RData")

############################################################# ############################################################# #############################################################
############################################################# Atlas: Dimension Reduction and Visualization #############################################################
############################################################# ############################################################# #############################################################

set.seed(seed=1)
setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/atlas")
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
library(patchwork)

# load("murine_atlas_data-integrated_v1.RData")

library(ggpubr)
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

DefaultAssay(atlas.object) <- "integrated"
atlas.object <- ScaleData(atlas.object, features = integrated.genes, assay = "integrated", vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo", "nFeature_RNA", "nCount_RNA", "nFeature_Spatial", "nCount_Spatial"))
atlas.object<- RunPCA(atlas.object, features = integrated.genes, assay = "integrated", slot = "scale.data")
atlas.object <- FindNeighbors(atlas.object, features = "integrated.genes", dims = 1:30)
atlas.object <- FindClusters(atlas.object, resolution = 0.6)
atlas.object <- RunUMAP(atlas.object, dims=1:30)
	### 26 clusters
p2 <- DimPlot(atlas.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p2
px <- DimPlot(atlas.object, reduction = "umap", group.by = "Platform",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Platform')
px


atlas.object <- DietSeurat(atlas.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/atlas")
DefaultAssay(atlas.object) <- "SCT"
save.image("murine_atlas_data-integrated_v1.RData")

Idents(atlas.object) <- "seurat_clusters"
t1<-table(Idents(atlas.object))
write.table(t1, "seurat_clusters.counts.txt", sep="\t")
rm(t1)

p1 <- DimPlot(atlas.object, reduction = "umap", group.by = "Treatment", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Treatment')
p1
p2 <- DimPlot(atlas.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p2
p3 <- DimPlot(atlas.object, reduction = "umap", group.by = "orig.ident", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Sample')
p3
p4 <- DimPlot(atlas.object, reduction = "umap", group.by = "Batch", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Batch')
p4
p5 <- DimPlot(atlas.object, reduction = "umap", group.by = "Type.Profile",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type.Profile')
p5
p6 <- DimPlot(atlas.object, reduction = "umap", group.by = "Subtype.Profile",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Subtype.Profile')
p6
p7 <- DimPlot(atlas.object, reduction = "umap", group.by = "Virus", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Virus')
p7
umap.combined <- p1 + p2 + p5 + p6
ggsave("Visium_UMAP-annotated_1.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("Visium_UMAP-annotated_XL_1.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")
ggsave("Visium_UMAP-annotated_M_1.pdf", plot = umap.combined, device = "pdf", width = 10, height = 6, units = "in")

umap.combined <- px + p1 + p5 + p6
ggsave("Visium_UMAP-annotated_1a.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("Visium_UMAP-annotated_XL_1a.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")
ggsave("Visium_UMAP-annotated_M_1a.pdf", plot = umap.combined, device = "pdf", width = 10, height = 6, units = "in")


p9 <- DimPlot(atlas.object, reduction = "umap", group.by = "Phase", label= "TRUE", repel=TRUE, raster=T)

umap.combined <- p7+p1+p3+p9
ggsave("Visium_UMAP-annotated_2.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("Visium_UMAP_annotated_XL_2.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")
ggsave("Visium_UMAP_annotated_M_2.pdf", plot = umap.combined, device = "pdf", width = 10, height = 6, units = "in")

#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("integrated_UMAP.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
p2 <- DimPlot(atlas.object, reduction = "umap", split.by = "orig.ident", label=TRUE)
ggsave("integrated_UMAP_splitby_libID.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
p2 <- DimPlot(atlas.object, reduction = "umap", split.by = "seurat_clusters", label=TRUE)
ggsave("integrated_UMAP_splitby_clusters.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
#examine UMAPs with qc metrics 
p5 <- FeaturePlot(atlas.object, features = 'nFeature_Spatial')
p6 <- FeaturePlot(atlas.object, features = 'nCount_Spatial')
p7 <- FeaturePlot(atlas.object, features = 'percent.mt')
p8 <- FeaturePlot(atlas.object, features = 'percent.ribo')
umap.combined <- p4+p5+p6+p7
ggsave("Visium_UMAP-annotated_3.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("Visium_UMAP_annotated_XL_3.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")
ggsave("Visium_UMAP-annotated_M_3.pdf", plot = umap.combined, device = "pdf", width = 10, height = 6, units = "in")

Idents(atlas.object) <- "orig.ident"

p4 <- DimPlot(atlas.object, reduction = "umap", group.by = "FetalSex", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='FetalSex')
ggsave("UMAP-FetalSex.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")
p4 <- FeaturePlot(atlas.object, features="DDX3Y", reduction = "umap", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='DDX3Y')
ggsave("UMAP-DDX3Y.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")

Idents(atlas.object) <- "Platform"

#examine UMAPs with qc metrics 
p5 <- FeaturePlot(atlas.object, features = 'nFeature_RNA')
p6 <- FeaturePlot(atlas.object, features = 'nCount_RNA')
p7 <- FeaturePlot(atlas.object, features = 'percent.mt')
p8 <- FeaturePlot(atlas.object, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots-RNA.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")
#examine UMAPs with qc metrics 
p5 <- FeaturePlot(atlas.object, features = 'nFeature_Spatial'
p6 <- FeaturePlot(atlas.object, features = 'nCount_Spatial')
p7 <- FeaturePlot(atlas.object, features = 'percent.mt')
p8 <- FeaturePlot(atlas.object, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots-Spatial.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")

Idents(atlas.object) <- "Platform"
p2 <- DimPlot(atlas.object, reduction = "umap", split.by = "Platform", label=TRUE)
ggsave("integrated_UMAP_splitby_platform-platform.pdf", plot = p2, device = "pdf", width = 9, height = 3, units = "in")
Idents(atlas.object) <- "Type.Profile"
p2 <- DimPlot(atlas.object, reduction = "umap", split.by = "Platform", label=TRUE)
ggsave("integrated_UMAP_splitby_Platform-Type.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
Idents(atlas.object) <- "Subtype.Profile"
p2 <- DimPlot(atlas.object, reduction = "umap", split.by = "Platform", label=TRUE)
ggsave("integrated_UMAP_splitby_Platform-Subtype.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")


plot2 <- SpatialFeaturePlot(atlas.object, features = "percent.viral", ncol = 1) + patchwork::plot_layout(ncol = 1)
####### FOR SOME REASON, INTEGRATION ADDS IMAGES ARBITRARILY
ggsave("integrated_SpatialDimPlot_percent.viral_TEST.pdf", plot = plot2, device = "pdf", width = 8, height = 45, units = "in")


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
############################################################# Atlas: DE Analysis #############################################################
############################################################# ############################################################# #############################################################

dir.create("immunology.clusters")
setwd("immunology.clusters")
Idents(atlas.object) <- "seurat_clusters"
## use RNA for all DE analysis and plots
DefaultAssay(atlas.object) <- "SCT"
# atlas.object <- ScaleData(atlas.object, features = all.genes, assay = "SCT")
atlas.object <- ScaleData(atlas.object, features = all.genes, assay="SCT")

## kernel density estimation 
# https://www.nature.com/articles/s41591-021-01329-2#Sec8
major.t <- c("CD4", "CD8A", "CCR7", "PTPRC", "SELL", "CD27", "CD38", "CD44", "CXCR5", "CD40LG", "CCR7", "FOXP3", "IKZF2")
gd.t<- c("TRGV9", "TRDV2")
MAIT.t <- c("TRAV10", "TRAJ18", "TRAV1-2")
NK.t <- c("CD3G", "FCGR3A", "NCAM1", "NCR1")
TH1.t <- c("IFNG", "TBX21", "TNF")
TH2.t <- c("GATA3", "IL4")
TH17.t <- c("RORC", "IL17A", "IL17F", "IL21")
## https://www.bioconductor.org/packages/release/bioc/vignettes/Nebulosa/inst/doc/nebulosa_seurat.html#1_Overview
  ## https://academic.oup.com/bioinformatics/article-abstract/37/16/2485/6103785
# BiocManager::install("Nebulosa")
library("Nebulosa")
Idents(atlas.object) <- "seurat_clusters"

# SARS.genes <- c("ORF1AB-1-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
p1 <- plot_density(atlas.object, "CD4")
p2 <- FeaturePlot(atlas.object, "CD8A")
p5 <- CombinePlots(plots = list(p1, p2))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("UMAP_KernelDensity_CD4-CD8A.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots
p3 <- plot_density(atlas.object, c("CD4", "CD8A"))
p5 <- p3 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_CD4-CD8A_Joint.pdf", plot = p5, device = "pdf", width = 5, height = 7, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(atlas.object, c("CD4", "CD8A"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_major.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(atlas.object, c("TRGV9", "TRDV2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_gd.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(atlas.object, c("TRAV10", "TRAJ18"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MAIT.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(atlas.object, c("CD3G", "FCGR3A","NCAM1", "NCR1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_NK.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(atlas.object, c("IFNG", "TBX21", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH1.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 16, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(atlas.object, c("GATA3", "IL4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH2.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(atlas.object, c("RORC", "IL17A", "IL17F", "IL21"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH17.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")

## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(atlas.object, features = lymphoid.lineage)
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
feature.plot <- DotPlot(atlas.object, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(atlas.object, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

macrophage.lineage <- c("CD14", "CD68", "CSF1R", "ADGRE1", "SIGLEC1", "CD38", "GPR18", "FPR2", "MYC", "EGR2")
feature.plot <- DotPlot(atlas.object, features = macrophage.lineage)
png(file=paste0("macrophage.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Macrophage Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()


# 65. Aagaard-Tillery KM, Silver R, Dalton J. Immunology of normal pregnancy. Semin Fetal Neonatal Med. Oct 2006;11(5):279-95. doi:10.1016/j.siny.2006.04.003
# 66. Ivashkiv LB. Epigenetic regulation of macrophage polarization and function. Trends Immunol. May 2013;34(5):216-23. doi:10.1016/j.it.2012.11.001
# 67. Murray PJ. Macrophage Polarization. Annu Rev Physiol. 02 10 2017;79:541-566. doi:10.1146/annurev-physiol-022516-034339
# 68. Yao Y, Xu XH, Jin L. Macrophage Polarization in Physiological and Pathological Pregnancy. Front Immunol. 2019;10:792. doi:10.3389/fimmu.2019.00792
# 69. Mues B, Langer D, Zwadlo G, Sorg C. Phenotypic characterization of macrophages in human term placenta. Immunology. Jul 1989;67(3):303-7. 
# 70. Bulmer JN, Johnson PM. Macrophage populations in the human placenta and amniochorion. Clin Exp Immunol. Aug 1984;57(2):393-403. 
# 71. Loegl J, Hiden U, Nussbaumer E, et al. Hofbauer cells of M2a, M2b and M2c polarization may regulate feto-placental angiogenesis. Reproduction. 2016;152(5):447-455. doi:10.1530/REP-16-0159
# 74. Schliefsteiner C, Ibesich S, Wadsack C. Placental Hofbauer Cell Polarization Resists Inflammatory Cues In Vitro. Int J Mol Sci. Jan 22 2020;21(3)doi:10.3390/ijms21030736
# 105.  Ben Amara A, Gorvel L, Baulan K, et al. Placental macrophages are impaired in chorioamnionitis, an infectious pathology of the placenta. J Immunol. Dec 01 2013;191(11):5501-14. doi:10.4049/jimmunol.1300988

macrophage.lineage.markers.2 <- c("CD14", "ITGAM", "CSF1R", "ADGRE1", "CD80", "CD38", "GPR18", "FPR2", "CD86", "TNF", "IL12A", "MYC", "EGR2", "CD163", "MRC1", "CD209", "TGFB1", "IL10", "VEGFA", "CD163", "HLA-DRA", "CD86", "IL6")
macrophage.lineage.markers.2 <- unique(macrophage.lineage.markers.2)

feature.plot <- DotPlot(atlas.object, features = macrophage.lineage.markers.2)
png(file=paste0("macrophage.lineage.markers2-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Macrophage Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

# https://www.nature.com/articles/s41591-021-01329-2#Sec8
t.lineages <- c("CD4", "CD8A", "CCR7", "PTPRC", "SELL", "CD27", "CD38", "CD44", "CXCR5", "CD40LG", "CCR7", "FOXP3", "IKZF2", "TRGV9", "TRDV2", "TRAV10", "TRAJ18", "TRAV1-2", "CD3G", "FCGR3A", "NCAM1", "NCR1", "IFNG", "TBX21", "TNF", "GATA3", "IL4", "RORC", "IL17A", "IL17F", "IL21")
t.lineages <- unique(t.lineages)

# https://www.nature.com/articles/s41591-021-01329-2#Sec8
major.t <- c("CD4", "CD8A", "CCR7", "PTPRC", "SELL", "CD27", "CD38", "CD44", "CXCR5", "CD40LG", "CCR7", "FOXP3", "IKZF2")
gd.t<- c("TRGV9", "TRDV2")
MAIT.t <- c("TRAV10", "TRAJ18", "TRAV1-2")
NK.t <- c("CD3G", "FCGR3A", "NCAM1", "NCR1")
TH1.t <- c("IFNG", "TBX21", "TNF")
TH2.t <- c("GATA3", "IL4")
TH17.t <- c("RORC", "IL17A", "IL17F", "IL21")

#atlas.object <- ScaleData(atlas.object, features = all.genes, assay = "SCT")
cluster.averages <- AverageExpression(atlas.object, assays="SCT", slot="counts", return.seurat = TRUE, features=t.lineages)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, slot="counts", assay="SCT")
ggsave("t.lineages_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, assay="SCT")
ggsave("t.lineages_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(atlas.object, assays="SCT", slot="counts", return.seurat = TRUE, features=macrophage.lineage.markers.2)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, slot="counts", assay="SCT")
ggsave("macrophage.lineage.markers.2_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, assay="SCT")
ggsave("macrophage.lineage.markers.2_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(atlas.object, assays="SCT", slot="counts", return.seurat = TRUE, features=lymphoid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, slot="counts", assay="SCT")
ggsave("lymphoid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, assay="SCT")
ggsave("lymphoid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(atlas.object, assays="SCT", slot="counts", return.seurat = TRUE, features=myeloid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, slot="counts", assay="SCT")
ggsave("myeloid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, assay="SCT")
ggsave("myeloid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


############## Activation gating strategy https://www.nature.com/articles/s41590-021-01049-2/figures/10
## 2 Density Plots and a Joint density plot
p4 <- plot_density(atlas.object, c("CD3D", "NCAM1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TorNK_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(atlas.object, c("CD4", "KLRB1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRB1-CD161.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(atlas.object, c("CD4", "CD38"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD38.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(atlas.object, c("CD4", "CD69"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD69.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("CD4", "KLRK1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRK1-NKG2D.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("CD4", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("FOXP3", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-FOXP3-CD4.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("KIT", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("TNF", "TNFRSF8"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.TNFRSF8.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p4 <- plot_density(atlas.object, c("CD14", "ITGAM"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("SIGLEC1", "FOLR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("SIGLEC1", "DDX3Y"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC-fetalsex_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
	## Placenta-associated maternal macrophages PAMMS (FOLR2 neg, HLA positive); HBCs FOLR2+HLA- ; Thomas et al. Naomi McGovern, JEM 2020
p4 <- plot_density(atlas.object, c("FOLR2", "HLA-DRA"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMMs_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM1_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM2_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("ITGAM", "CD14"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M0.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("CD80", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("CD11B", "MYC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("CD163", "MRC1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2a.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("HLA-DRA", "IL6"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("CD14", "CD163"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2c.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("CD14", "CD1C"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(atlas.object, c("CD14", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p1 <- plot_density(atlas.object, "THBD")
p2 <- FeaturePlot(atlas.object, "THBD")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(atlas.object, "ITGB3")
p2 <- FeaturePlot(atlas.object, "ITGB3")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_Megakaryocyte.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(atlas.object, "CD1C")
p2 <- FeaturePlot(atlas.object, "CD1C")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(atlas.object, "IL6")
p2 <- FeaturePlot(atlas.object, "IL6")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(atlas.object, "TNF")
p2 <- FeaturePlot(atlas.object, "TNF")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(atlas.object, "MYC")
p2 <- FeaturePlot(atlas.object, "MYC")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(atlas.object, "FOLR2")
p2 <- FeaturePlot(atlas.object, "FOLR2")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_HBCs-FOLR2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(atlas.object, "HLA-DRA")
p2 <- FeaturePlot(atlas.object, "HLA-DRA")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_PAMM-HLA-DRA.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(atlas.object, "CD33")
p2 <- FeaturePlot(atlas.object, "CD33")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeliodProgenitor-M0.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p4 <- plot_density(atlas.object, c("KIT", "CD33"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_MyeliodProgenitor-M0-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
# atlas.object <- DietSeurat(atlas.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
setwd("..")

rm(cluster.averages)
rm(library.averages.heatmap)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(p5)
rm(p6)
rm(p7)
rm(p8)
rm(plot2)
rm(umap.combined)
rm(feature.plot)
rm(de_markers)

##### Perform differential expression between seurat_clusters
dir.create("DE_seurat_clusters")
setwd("DE_seurat_clusters")
	## all_merged <- SCTransform(all_merged, method = "glmGamPoi", assay = "SCT", conserve.memory=TRUE, verbose = T)
# atlas.object <- SCTransform(atlas.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes = FALSE, verbose = T)
DefaultAssay(atlas.object) <- "SCT"
# atlas.object <- NormalizeData(atlas.object, normalization.method = "LogNormalize", scale.factor = 10000)
# atlas.object <- ScaleData(atlas.object, features = all.genes, assay="SCT")
atlas.object <- SCTransform(atlas.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes=TRUE, conserve.memory=TRUE, verbose = T)

var.genes <- atlas.object@assays[["SCT"]]@var.features
combined.var.genes <- union(var.genes, SARS.genes.receptors.hb.fetalsex.mac)
atlas.object <- NormalizeData(atlas.object, features = combined.var.genes, assay="SCT")
atlas.object <- ScaleData(atlas.object, features = combined.var.genes, assay="SCT")
Idents(atlas.object) <- "seurat_clusters"
#de_markers <- FindAllMarkers(atlas.object, features = var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
#write.table(de_markers, "DEGs_by_seurat_clusters_pos-0.log2FC.txt", sep="\t")
# de_markers <- FindAllMarkers(atlas.object, features = combined.var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
de_markers <- FindAllMarkers(atlas.object, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_seurat_clusters_pos-0.693lnFC.txt", sep="\t")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(atlas.object, features = unique.top2)
png(file=paste0("seurat_clusters_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

cluster.averages <- AverageExpression(atlas.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="SCT")
ggsave("top2_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="SCT")
ggsave("top2_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(atlas.object, features = unique.top2)
png(file=paste0("seurat_clusters_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(atlas.object, features = unique.top2)
png(file=paste0("seurat_clusters_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()


## Slim down the atlas.objects to save space. Don't need to keep the scale.data
# atlas.object <- DietSeurat(atlas.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
setwd("..")
dir.create("DE_Type.Profile")
setwd("DE_Type.Profile")
	## all_merged <- SCTransform(all_merged, method = "glmGamPoi", assay = "SCT", conserve.memory=TRUE, verbose = T)
# atlas.object <- SCTransform(atlas.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes = FALSE, verbose = T)
DefaultAssay(atlas.object) <- "SCT"
# atlas.object <- NormalizeData(atlas.object, normalization.method = "LogNormalize", scale.factor = 10000)
# atlas.object <- ScaleData(atlas.object, features = all.genes, assay="SCT")
atlas.object <- SCTransform(atlas.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes=TRUE, conserve.memory=TRUE, verbose = T)

var.genes <- atlas.object@assays[["SCT"]]@var.features
combined.var.genes <- union(var.genes, SARS.genes.receptors.hb.fetalsex.mac)
atlas.object <- NormalizeData(atlas.object, features = combined.var.genes, assay="SCT")
atlas.object <- ScaleData(atlas.object, features = combined.var.genes, assay="SCT")
Idents(atlas.object) <- "Type.Profile"
#de_markers <- FindAllMarkers(atlas.object, features = var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
#write.table(de_markers, "DEGs_by_Type.Profile_pos-0.log2FC.txt", sep="\t")
# de_markers <- FindAllMarkers(atlas.object, features = combined.var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
de_markers <- FindAllMarkers(atlas.object, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_Type.Profile_pos-0.693lnFC.txt", sep="\t")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(atlas.object, features = unique.top2)
png(file=paste0("Type.Profile_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

cluster.averages <- AverageExpression(atlas.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="SCT")
ggsave("top2_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="SCT")
ggsave("top2_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(atlas.object, features = unique.top2)
png(file=paste0("Type.Profile_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(atlas.object, features = unique.top2)
png(file=paste0("Type.Profile_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()


## Slim down the atlas.objects to save space. Don't need to keep the scale.data
# atlas.object <- DietSeurat(atlas.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
setwd("..")
dir.create("DE_Subtype.Profile")
setwd("DE_Subtype.Profile")
	## all_merged <- SCTransform(all_merged, method = "glmGamPoi", assay = "SCT", conserve.memory=TRUE, verbose = T)
# atlas.object <- SCTransform(atlas.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes = FALSE, verbose = T)
DefaultAssay(atlas.object) <- "SCT"
# atlas.object <- NormalizeData(atlas.object, normalization.method = "LogNormalize", scale.factor = 10000)
# atlas.object <- ScaleData(atlas.object, features = all.genes, assay="SCT")
atlas.object <- SCTransform(atlas.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes=TRUE, conserve.memory=TRUE, verbose = T)

var.genes <- atlas.object@assays[["SCT"]]@var.features
combined.var.genes <- union(var.genes, SARS.genes.receptors.hb.fetalsex.mac)
atlas.object <- NormalizeData(atlas.object, features = combined.var.genes, assay="SCT")
atlas.object <- ScaleData(atlas.object, features = combined.var.genes, assay="SCT")
Idents(atlas.object) <- "Subtype.Profile"
#de_markers <- FindAllMarkers(atlas.object, features = var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
#write.table(de_markers, "DEGs_by_Subtype.Profile_pos-0.log2FC.txt", sep="\t")
# de_markers <- FindAllMarkers(atlas.object, features = combined.var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
de_markers <- FindAllMarkers(atlas.object, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_Subtype.Profile_pos-0.693lnFC.txt", sep="\t")



top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(atlas.object, features = unique.top2)
png(file=paste0("Subtype.Profile_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

cluster.averages <- AverageExpression(atlas.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="SCT")
ggsave("top2_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="SCT")
ggsave("top2_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(atlas.object, features = unique.top2)
png(file=paste0("Subtype.Profile_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(atlas.object, features = unique.top2)
png(file=paste0("Subtype.Profile_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()


## Slim down the atlas.objects to save space. Don't need to keep the scale.data
# atlas.object <- DietSeurat(atlas.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
setwd("..")
############################################################# ############################################################# #############################################################
############################################################# Atlas Cell Counts  #############################################################
############################################################# ############################################################# #############################################################

#Save tables with cell counts per library and per cluster
dir.create("cell.counts")
setwd("cell.counts")
# Save a table with cell counts per library
Idents(atlas.object) <- "orig.ident"
unfiltered.count <- table(Idents(atlas.object))
write.table(unfiltered.count, "sample.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(atlas.object) <- "seurat_clusters"
unfiltered.count <- table(Idents(atlas.object))
write.table(unfiltered.count, "cluster.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(atlas.object) <- "Subtype.Profile"
unfiltered.count <- table(Idents(atlas.object))
write.table(unfiltered.count, "Subtype.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(atlas.object) <- "Type.Profile"
unfiltered.count <- table(Idents(atlas.object))
write.table(unfiltered.count, "Type.cellcounts.txt", sep="\t")

############################################################# ############################################################# #############################################################
############################################################# Platform.Type.Profile #############################################################
############################################################# ############################################################# #############################################################


dir.create("Platform.Type.Profile")
setwd("Platform.Type.Profile")

# "0", "2", "1", "4", "3", "5"
# reference.list <- all_merged.list[c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", 
# "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma", "FetalMesenchyme/Trophoblast", "FetalMesenchyme/Immune/Endothelial")]


Idents(atlas.object)<-"Platform"
Type.Profile.0.uninfected <- subset(atlas.object, idents = c("Visium"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "Visium.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

Idents(atlas.object)<-"Platform"
Type.Profile.0.uninfected <- subset(atlas.object, idents = c("snRNA-seq"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "snRNA-seq.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

setwd("..")


# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("Platform.Type.Profile/"))
manifest<-manifest%>%mutate(name=gsub(".Type.Profile.cellcounts.txt","",value))
manifest
		## These commands below are required, change the prefix to match your samples
df<-read.table("Platform.Type.Profile/Visium.Type.Profile.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type.Profile","Visium")
 df2<-df
df<-df%>%pivot_longer(cols = `Visium`,names_to = "Platform",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("Platform.Type.Profile/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Type.Profile",X)
  df<-df%>%pivot_longer(cols = X,names_to = "Platform",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of Type.Profiles ## there are 12 Platforms
# df2%>%mutate(perc=percs_by_group(count,group = Platform))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])


final<-full_join(df2,df3)

final<-final%>%
  mutate(Platform=gsub(".Type.Profile.cellcounts.txt","",Platform))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = Platform))
write.table(final,"Platform-Type.Profile.CellCounts_merged.tsv",sep = "\t",row.names = F)
final
​                          
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18
final$Type.Profile<-as.factor(final$Type.Profile)
final$Platform<-as.factor(final$Platform)
final
​
library(ggpubr)
ggbarplot(data = final,
          x = "Platform",
          y = "Percent",
          color = "Type.Profile",
          fill = "Type.Profile",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal3)

png(file=paste0("Platform.Type.Profile_diversity.plot.png"),
                res=300, 
                width=1500, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "Platform",
          y = "Percent",
          color = "Type.Profile",
          fill = "Type.Profile",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal3)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()

############################################################# ############################################################# #############################################################
############################################################# Platform.Subtype.Profile #############################################################
############################################################# ############################################################# #############################################################


dir.create("Platform.Subtype.Profile")
setwd("Platform.Subtype.Profile")

# "0", "2", "1", "4", "3", "5"
# reference.list <- all_merged.list[c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", 
# "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma", "FetalMesenchyme/Trophoblast", "FetalMesenchyme/Immune/Endothelial")]


Idents(atlas.object)<-"Platform"
Subtype.Profile.0.uninfected <- subset(atlas.object, idents = c("Visium"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "Visium.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

Idents(atlas.object)<-"Platform"
Subtype.Profile.0.uninfected <- subset(atlas.object, idents = c("snRNA-seq"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "snRNA-seq.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

setwd("..")


# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("Platform.Subtype.Profile/"))
manifest<-manifest%>%mutate(name=gsub(".Subtype.Profile.cellcounts.txt","",value))
manifest
		## These commands below are required, change the prefix to match your samples
df<-read.table("Platform.Subtype.Profile/Visium.Subtype.Profile.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Subtype.Profile","Visium")
 df2<-df
df<-df%>%pivot_longer(cols = `Visium`,names_to = "Platform",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("Platform.Subtype.Profile/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Subtype.Profile",X)
  df<-df%>%pivot_longer(cols = X,names_to = "Platform",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of Subtype.Profiles ## there are 12 Platforms
# df2%>%mutate(perc=percs_by_group(count,group = Platform))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])


final<-full_join(df2,df3)

final<-final%>%
  mutate(Platform=gsub(".Subtype.Profile.cellcounts.txt","",Platform))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = Platform))
#  pivot_wider(names_from = Platform,values_from = count)

write.table(final,"Platform.Subtype.Profile.CellCounts_merged.tsv",sep = "\t",row.names = F)

final
​                           
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18
final$Subtype.Profile<-as.factor(final$Subtype.Profile)
final$Platform<-as.factor(final$Platform)
final
​
library(ggpubr)
ggbarplot(data = final,
          x = "Platform",
          y = "Percent",
          color = "Subtype.Profile",
          fill = "Subtype.Profile",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal3)

png(file=paste0("Platform.Subtype.Profile_diversity.plot.png"),
                res=300, 
                width=1500, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "Platform",
          y = "Percent",
          color = "Subtype.Profile",
          fill = "Subtype.Profile",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal3)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()


############################################################# ############################################################# #############################################################
############################################################# orig.ident.Type.Profile #############################################################
############################################################# ############################################################# #############################################################
setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/annotated")
# load("murine_spatial_data-annotated_v1.RData")

dir.create("cell.counts")
setwd("cell.counts")
# Save a table with cell counts per library
Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "sample.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(seurat.object) <- "seurat_clusters"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "cluster.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(seurat.object) <- "Subtype.Profile"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "Subtype.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(seurat.object) <- "Type.Profile"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "Type.cellcounts.txt", sep="\t")
# Save a table with cell counts per cluster
Idents(seurat.object) <- "Treatment"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "Treatment.cellcounts.txt", sep="\t")


dir.create("Treatment.Type.Profile")
setwd("Treatment.Type.Profile")
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

#Save tables with cell counts per library and per cluster
# "Mock-Enox-1" "Mock-Enox-2" "ZIKV-Veh-1"  "ZIKV-Enox-1" "ZIKV-Enox-2" "Mock-Veh-1"  "ZIKV-Veh-2" "ZIKV-Veh-3" 


Idents(seurat.object)<-"Treatment"
Type.Profile.0.uninfected <- subset(seurat.object, idents = c("Mock-Veh"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "Mock-Veh.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

Idents(seurat.object)<-"Treatment"
Type.Profile.0.uninfected <- subset(seurat.object, idents = c("Mock-Enox"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "Mock-Enox.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

Idents(seurat.object)<-"Treatment"
Type.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Veh"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Veh.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

Idents(seurat.object)<-"Treatment"
Type.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Enox"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Enox.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

setwd("..")


# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("Treatment.Type.Profile/"))
manifest<-manifest%>%mutate(name=gsub(".Type.Profile.cellcounts.txt","",value))
manifest
		## These commands below are required, change the prefix to match your samples
df<-read.table("Treatment.Type.Profile/Mock-Veh.Type.Profile.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type.Profile","Mock-Veh")
 df2<-df
df<-df%>%pivot_longer(cols = `Mock-Veh`,names_to = "Treatment",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("Treatment.Type.Profile/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Type.Profile",X)
  df<-df%>%pivot_longer(cols = X,names_to = "Treatment",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of Type.Profiles ## there are 12 Treatments
# df2%>%mutate(perc=percs_by_group(count,group = Treatment))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)

final<-final%>%
  mutate(Treatment=gsub(".Type.Profile.cellcounts.txt","",Treatment))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = Treatment))
write.table(final,"Treatment-Type.Profile.CellCounts_merged.tsv",sep = "\t",row.names = F)
final
​                          
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18
final$Type.Profile<-as.factor(final$Type.Profile)
final$Treatment<-as.factor(final$Treatment)
final
​
library(ggpubr)
ggbarplot(data = final,
          x = "Treatment",
          y = "Percent",
          color = "Type.Profile",
          fill = "Type.Profile",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal3)

png(file=paste0("Treatment.Type.Profile_diversity.plot.png"),
                res=300, 
                width=1500, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "Treatment",
          y = "Percent",
          color = "Type.Profile",
          fill = "Type.Profile",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal3)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()

############################################################# ############################################################# #############################################################
dir.create("Treatment.Subtype.Profile")
setwd("Treatment.Subtype.Profile")
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

#Save tables with cell counts per library and per cluster
# "Mock-Enox-1" "Mock-Enox-2" "ZIKV-Veh-1"  "ZIKV-Enox-1" "ZIKV-Enox-2" "Mock-Veh-1"  "ZIKV-Veh-2" "ZIKV-Veh-3" 


Idents(seurat.object)<-"Treatment"
Subtype.Profile.0.uninfected <- subset(seurat.object, idents = c("Mock-Veh"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "Mock-Veh.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

Idents(seurat.object)<-"Treatment"
Subtype.Profile.0.uninfected <- subset(seurat.object, idents = c("Mock-Enox"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "Mock-Enox.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

Idents(seurat.object)<-"Treatment"
Subtype.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Veh"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Veh.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

Idents(seurat.object)<-"Treatment"
Subtype.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Enox"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Enox.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

setwd("..")


# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("Treatment.Subtype.Profile/"))
manifest<-manifest%>%mutate(name=gsub(".Subtype.Profile.cellcounts.txt","",value))
manifest
		## These commands below are required, change the prefix to match your samples
df<-read.table("Treatment.Subtype.Profile/Mock-Veh.Subtype.Profile.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Subtype.Profile","Mock-Veh")
 df2<-df
df<-df%>%pivot_longer(cols = `Mock-Veh`,names_to = "Treatment",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("Treatment.Subtype.Profile/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Subtype.Profile",X)
  df<-df%>%pivot_longer(cols = X,names_to = "Treatment",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of Subtype.Profiles ## there are 12 Treatments
# df2%>%mutate(perc=percs_by_group(count,group = Treatment))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)

final<-final%>%
  mutate(Treatment=gsub(".Subtype.Profile.cellcounts.txt","",Treatment))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = Treatment))
write.table(final,"Treatment-Subtype.Profile.CellCounts_merged.tsv",sep = "\t",row.names = F)
final
​                          
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18
final$Subtype.Profile<-as.factor(final$Subtype.Profile)
final$Treatment<-as.factor(final$Treatment)
final
​
library(ggpubr)
ggbarplot(data = final,
          x = "Treatment",
          y = "Percent",
          color = "Subtype.Profile",
          fill = "Subtype.Profile",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal3)

png(file=paste0("Treatment.Subtype.Profile_diversity.plot.png"),
                res=300, 
                width=1500, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "Treatment",
          y = "Percent",
          color = "Subtype.Profile",
          fill = "Subtype.Profile",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal3)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()
############################################################# ############################################################# #############################################################
#

dir.create("orig.ident.Type.Profile")
setwd("orig.ident.Type.Profile")
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

#Save tables with cell counts per library and per cluster
# "Mock-Enox-1" "Mock-Enox-2" "ZIKV-Veh-1"  "ZIKV-Enox-1" "ZIKV-Enox-2" "Mock-Veh-1"  "ZIKV-Veh-2" "ZIKV-Veh-3" 


Idents(seurat.object)<-"orig.ident"
Type.Profile.0.uninfected <- subset(seurat.object, idents = c("Mock-Veh-1"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "Mock-Veh-1.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Type.Profile.0.uninfected <- subset(seurat.object, idents = c("Mock-Enox-1"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "Mock-Enox-1.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Type.Profile.0.uninfected <- subset(seurat.object, idents = c("Mock-Enox-2"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "Mock-Enox-2.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Type.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Veh-1"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Veh-1.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Type.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Veh-2"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Veh-2.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Type.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Veh-3"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Veh-3.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Type.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Enox-1"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Enox-1.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Type.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Enox-2"))
Idents(Type.Profile.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Type.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Enox-2.Type.Profile.cellcounts.txt", sep="\t")
rm("Type.Profile.0.uninfected")

setwd("..")


# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("orig.ident.Type.Profile/"))
manifest<-manifest%>%mutate(name=gsub(".Type.Profile.cellcounts.txt","",value))
manifest
		## These commands below are required, change the prefix to match your samples
df<-read.table("orig.ident.Type.Profile/Mock-Veh-1.Type.Profile.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type.Profile","Mock-Veh-1")
 df2<-df
df<-df%>%pivot_longer(cols = `Mock-Veh-1`,names_to = "orig.ident",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("orig.ident.Type.Profile/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Type.Profile",X)
  df<-df%>%pivot_longer(cols = X,names_to = "orig.ident",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of Type.Profiles ## there are 12 orig.idents
# df2%>%mutate(perc=percs_by_group(count,group = orig.ident))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])
df7<-my_fxn(X = manifest$value[6])
df8<-my_fxn(X = manifest$value[7])
df9<-my_fxn(X = manifest$value[8])
final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)
final<-full_join(final,df7)
final<-full_join(final,df8)
final<-full_join(final,df9)

final<-final%>%
  mutate(orig.ident=gsub(".Type.Profile.cellcounts.txt","",orig.ident))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = orig.ident))
write.table(final,"orig.ident-Type.Profile.CellCounts_merged.tsv",sep = "\t",row.names = F)
final
​                          
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18
final$Type.Profile<-as.factor(final$Type.Profile)
final$orig.ident<-as.factor(final$orig.ident)
final
​
library(ggpubr)
ggbarplot(data = final,
          x = "orig.ident",
          y = "Percent",
          color = "Type.Profile",
          fill = "Type.Profile",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal3)

png(file=paste0("orig.ident.Type.Profile_diversity.plot.png"),
                res=300, 
                width=1500, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "orig.ident",
          y = "Percent",
          color = "Type.Profile",
          fill = "Type.Profile",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal3)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()

############################################################# ############################################################# #############################################################
dir.create("orig.ident.Subtype.Profile")
setwd("orig.ident.Subtype.Profile")
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

#Save tables with cell counts per library and per cluster
# "Mock-Enox-1" "Mock-Enox-2" "ZIKV-Veh-1"  "ZIKV-Enox-1" "ZIKV-Enox-2" "Mock-Veh-1-1"  "ZIKV-Veh-2" "ZIKV-Veh-3" 


Idents(seurat.object)<-"orig.ident"
Subtype.Profile.0.uninfected <- subset(seurat.object, idents = c("Mock-Veh-1"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "Mock-Veh-1.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Subtype.Profile.0.uninfected <- subset(seurat.object, idents = c("Mock-Enox-1"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "Mock-Enox-1.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Subtype.Profile.0.uninfected <- subset(seurat.object, idents = c("Mock-Enox-2"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "Mock-Enox-2.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Subtype.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Veh-1"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Veh-1.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Subtype.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Veh-2"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Veh-2.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Subtype.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Veh-3"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Veh-3.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Subtype.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Enox-1"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Enox-1.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

Idents(seurat.object)<-"orig.ident"
Subtype.Profile.0.uninfected <- subset(seurat.object, idents = c("ZIKV-Enox-2"))
Idents(Subtype.Profile.0.uninfected) <- "Subtype.Profile"
unfiltered.count <- table(Idents(Subtype.Profile.0.uninfected))
write.table(unfiltered.count, "ZIKV-Enox-2.Subtype.Profile.cellcounts.txt", sep="\t")
rm("Subtype.Profile.0.uninfected")

setwd("..")


# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("orig.ident.Subtype.Profile/"))
manifest<-manifest%>%mutate(name=gsub(".Subtype.Profile.cellcounts.txt","",value))
manifest
		## These commands below are required, change the prefix to match your samples
df<-read.table("orig.ident.Subtype.Profile/Mock-Veh-1.Subtype.Profile.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Subtype.Profile","Mock-Veh-1")
 df2<-df
df<-df%>%pivot_longer(cols = `Mock-Veh-1`,names_to = "orig.ident",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("orig.ident.Subtype.Profile/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Subtype.Profile",X)
  df<-df%>%pivot_longer(cols = X,names_to = "orig.ident",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of Subtype.Profiles ## there are 12 orig.idents
# df2%>%mutate(perc=percs_by_group(count,group = orig.ident))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])
df7<-my_fxn(X = manifest$value[6])
df8<-my_fxn(X = manifest$value[7])
df9<-my_fxn(X = manifest$value[8])
final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)
final<-full_join(final,df7)
final<-full_join(final,df8)
final<-full_join(final,df9)

final<-final%>%
  mutate(orig.ident=gsub(".Subtype.Profile.cellcounts.txt","",orig.ident))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = orig.ident))
write.table(final,"orig.ident-Subtype.Profile.CellCounts_merged.tsv",sep = "\t",row.names = F)
final
​                          
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18
final$Subtype.Profile<-as.factor(final$Subtype.Profile)
final$orig.ident<-as.factor(final$orig.ident)
final
​
library(ggpubr)
ggbarplot(data = final,
          x = "orig.ident",
          y = "Percent",
          color = "Subtype.Profile",
          fill = "Subtype.Profile",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal3)

png(file=paste0("orig.ident.Subtype.Profile_diversity.plot.png"),
                res=300, 
                width=1500, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "orig.ident",
          y = "Percent",
          color = "Subtype.Profile",
          fill = "Subtype.Profile",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal3)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()
############################################################# ############################################################# #############################################################
#

############################################################# ############################################################# #############################################################
############################################################# End work with Visum + snRNA-seq atlas #############################################################
############################################################# ############################################################# #############################################################

############################################################# ############################################################# #############################################################
############################################################# VI. analysis of visium data using annotations #############################################################
############################################################# ############################################################# #############################################################

############################################################# ############################################################# #############################################################
############################################################# Gene set analysis #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/annotated")
# load("murine_spatial_data-annotated_v1.RData")

DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")

Idents(seurat.object) <- "orig.ident"
gene.list <- c("DDX3Y")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "DDX3Y", images=image.list, ncol=4)
p5
ggsave("DDX3Y_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "DDX3Y", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 3)
p5
ggsave("DDX3Y_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "SRY", images=image.list, ncol=4)
p5
ggsave("SRY_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "SRY", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 3)
p5
ggsave("SRY_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "CYB561D2", images=image.list, ncol=4)
p5
ggsave("CYB561D2_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "CYB561D2", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 2)
p5
ggsave("CYB561D2_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "TMSB15B2", images=image.list, ncol=4)
p5
ggsave("TMSB15B2_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "TMSB15B2", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 2)
p5
ggsave("TMSB15B2_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "TMSB15B1", images=image.list, ncol=4)
p5
ggsave("TMSB15B1_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "TMSB15B1", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 2)
p5
ggsave("TMSB15B1_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "CD14", images=image.list, ncol=4)
p5
ggsave("CD14_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "CD14", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 2)
p5
ggsave("CD14_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "SDK1", images=image.list, ncol=4)
p5
ggsave("SDK1_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "SDK1", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 2)
p5
ggsave("SDK1_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "SLC12A8", images=image.list, ncol=4)
p5
ggsave("SLC12A8_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "SLC12A8", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 2)
p5
ggsave("SLC12A8_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "VLDLR", images=image.list, ncol=4)
p5
ggsave("VLDLR_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "VLDLR", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 2)
p5
ggsave("VLDLR_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "PSTPIP1", images=image.list, ncol=4)
p5
ggsave("PSTPIP1_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "PSTPIP1", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 2)
p5
ggsave("PSTPIP1_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "IFFO2", images=image.list, ncol=4)
p5
ggsave("IFFO2_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "IFFO2", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 2)
p5
ggsave("IFFO2_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "TNK2", images=image.list, ncol=4)
p5
ggsave("TNK2_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "TNK2", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 2)
p5
ggsave("TNK2_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "PPIF", images=image.list, ncol=4)
p5
ggsave("PPIF_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "PPIF", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 2)
p5
ggsave("PPIF_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "SPATA2L", images=image.list, ncol=4)
p5
ggsave("SPATA2L_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "SPATA2L", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 2)
p5
ggsave("SPATA2L_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "HIST1H2BJ", images=image.list, ncol=4)
p5
ggsave("HIST1H2BJ_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "HIST1H2BJ", images=image.list, ncol=4, min.cutoff = 0.000001, max.cutoff = 2)
p5
ggsave("HIST1H2BJ_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")




Idents(seurat.object) <- "Treatment"
gene.list <- c("C1R", "C1S", "C1QA", "C1QB", "C2", "C3", "C4A", "C4B", "C5", "C6", "C7", "C8", "C9")


## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=4)
p5
ggsave("C3_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("C3_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object, features = "C3",)
p5
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = "C3")
p5
ggsave("C3_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = gene.list)
p5
ggsave("complement.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")


image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("C3_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")



gene.list <- c("C1RA", "C1RB", "C1QA", "C1QB", "C2", "C3", "C4A", "C4B", "C5", "C6", "C7", "C8A","C8B","C8G", "C9")

Idents(seurat.object) <- "Treatment"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("complement-c9-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="Treatment")
dev.off()

seurat.object <- ScaleData(seurat.object, features = gene.list, assay="Spatial")


top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("complement-c9_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("complement-c9_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("complement-c9_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("complement-c9_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "orig.ident"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("orig.ident_complement-c9-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("orig.ident_complement-c9_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("orig.ident_complement-c9_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("orig.ident_complement-c9_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("orig.ident_complement-c9_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "Type.Profile"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("Type.Profile_complement-c9-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("Type.Profile_complement-c9_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Type.Profile_complement-c9_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Type.Profile_complement-c9_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Type.Profile_complement-c9_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "Subtype.Profile"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("Subtype.Profile_complement-c9-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("Subtype.Profile_complement-c9_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Subtype.Profile_complement-c9_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Subtype.Profile_complement-c9_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Subtype.Profile_complement-c9_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
############################################################# ############################################################# #############################################################
dir.create("ZIKV-Veh_complement_.gene.list")
setwd("ZIKV-Veh_complement_.gene.list")
## ZIKV-Veh Treatment Gene List
Idents(seurat.object) <- "Treatment"
gene.list <- c("CTLA2A", "KAP","PRAP1", "C1R", "C1S", "C1QA", "C1QB", "C2", "C3", "C4A", "C4B", "C5", "C6", "C7", "C8", "C9")
seurat.object <- ScaleData(seurat.object, features = gene.list, assay="Spatial")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("ZIKV-Veh.gene.list_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("ZIKV-Veh.gene.list_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "orig.ident"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("orig.ident_ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("orig.ident_ZIKV-Veh.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("orig.ident_ZIKV-Veh.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("orig.ident_ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("orig.ident_ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "Type.Profile"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("Type.Profile_ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("Type.Profile_ZIKV-Veh.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Type.Profile_ZIKV-Veh.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Type.Profile_ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Type.Profile_ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "Subtype.Profile"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("Subtype.Profile_ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
setwd("..")

dir.create("miRNA.gene.list")
setwd("miRNA.gene.list")
Idents(seurat.object) <- "Treatment"

gene.list <- c("SLC12A8", "SDK1", "VLDLR", "PSTPIP1", "DNAJB1", "IFFO2", "TNK2", "CYB561D2", "PPIF", "SPATA2L", "HERC2","HIST1H2BJ","SLC44A4","FAM189A2","PARD3B", "CYB561D2", "TMSB15B1","TMSB15B2")
gene.list <- c("SLC12A8", "LHFPL2", "DENND1A", "PAX8-AS1", "XPO4", "NR3C2", "HECTD1", "LPIN2", "FER", "ACCS", "AFF4", "VLDLR", "LMTK2", "ARHGEF11", "DOCK9", "SDK1", "DDX46")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("mirna.gene.list_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("mirna.gene.list_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("mirna.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("mirna.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "orig.ident"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("orig.ident_mirna.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("orig.ident_mirna.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("orig.ident_mirna.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("orig.ident_mirna.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("orig.ident_mirna.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "Type.Profile"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("Type.Profile_mirna.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("Type.Profile_mirna.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Type.Profile_mirna.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Type.Profile_mirna.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Type.Profile_mirna.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "Subtype.Profile"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("Subtype.Profile_mirna.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("Subtype.Profile_mirna.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Subtype.Profile_mirna.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Subtype.Profile_mirna.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Subtype.Profile_mirna.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
setwd("..")

dir.create("immune.ZIKV-Veh.gene.list")
setwd("immune.ZIKV-Veh.gene.list")
Idents(seurat.object) <- "Treatment"
## ZIKV-Veh Treatment Gene List
Idents(seurat.object) <- "Treatment"
gene.list <- c("CTLA2A", "KAP", "PRAP1", "PLA1A", "CCL21A",  "C3", "LYZ2", "CCL8")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")

p5 <- SpatialPlot(seurat.object, features = "CTLA2A", images=image.list, ncol=4)
p5
ggsave("CTLA2A_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "CTLA2A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CTLA2A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object, features = "CTLA2A", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object, features = "CTLA2A",)
p5
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = "CTLA2A")
p5
ggsave("CTLA2A_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "CTLA2A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CTLA2A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


p5 <- SpatialPlot(seurat.object, features = "PLA1A", images=image.list, ncol=4)
p5
ggsave("PLA1A_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "PLA1A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("PLA1A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object, features = "PLA1A", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object, features = "PLA1A",)
p5
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = "PLA1A")
p5
ggsave("PLA1A_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "PLA1A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("PLA1A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


p5 <- SpatialPlot(seurat.object, features = "CCL21A", images=image.list, ncol=4)
p5
ggsave("CCL21A_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "CCL21A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CCL21A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object, features = "CCL21A", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object, features = "CCL21A",)
p5
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = "CCL21A")
p5
ggsave("CCL21A_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "CCL21A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CCL21A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=4)
p5
ggsave("C3_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("C3_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object, features = "C3",)
p5
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = "C3")
p5
ggsave("C3_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("C3_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


p5 <- SpatialPlot(seurat.object, features = "LYZ2", images=image.list, ncol=4)
p5
ggsave("LYZ2_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "LYZ2", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("LYZ2_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object, features = "LYZ2", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object, features = "LYZ2",)
p5
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = "LYZ2")
p5
ggsave("LYZ2_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "LYZ2", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("LYZ2_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


p5 <- SpatialPlot(seurat.object, features = "CCL8", images=image.list, ncol=4)
p5
ggsave("CCL8_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "CCL8", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CCL8_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object, features = "CCL8", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object, features = "CCL8",)
p5
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = "CCL8")
p5
ggsave("CCL8_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "CCL8", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CCL8_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


p5 <- SpatialPlot(seurat.object, features = "PRAP1", images=image.list, ncol=4)
p5
ggsave("PRAP1_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "PRAP1", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("PRAP1_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object, features = "PRAP1", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object, features = "PRAP1",)
p5
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = "PRAP1")
p5
ggsave("PRAP1_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "PRAP1", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("PRAP1_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

p5 <- SpatialPlot(seurat.object, features = "KAP", images=image.list, ncol=4)
p5
ggsave("KAP_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "KAP", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("KAP_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object, features = "KAP", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object, features = "KAP",)
p5
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = "KAP")
p5
ggsave("KAP_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "KAP", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("KAP_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")





Idents(seurat.object) <- "Treatment"

gene.list <- c("CTLA2A", "PRAP1", "CCL21A", "KAP", "SRGN", "PENK", "IGFBP6", "HTRA3", "AEBP1", "GPX3", "TAGLN", "FGL2", "DES", "LBP", "SLPI", "CNN1", "GUCA2B", "LY6C1", "C3", "PLA1A", "CTSK", "IFI27L2A", "POSTN", "LYZ2", "SFRP4", "GSN", "TIMP2", "ACTG2", "PLTP", "IGFBP2", "DCN", "SPARCL1", "CCL8")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("gene.list_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("gene.list_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "orig.ident"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("orig.ident_gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("orig.ident_gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("orig.ident_gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("orig.ident_gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("orig.ident_gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "Type.Profile"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("Type.Profile_gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("Type.Profile_gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Type.Profile_gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Type.Profile_gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Type.Profile_gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "Subtype.Profile"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("Subtype.Profile_gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("Subtype.Profile_gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Subtype.Profile_gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Subtype.Profile_gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Subtype.Profile_gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
setwd("..")
############################################################# ############################################################# #############################################################

dir.create("immunology.clusters")
setwd("immunology.clusters")
Idents(seurat.object) <- "seurat_clusters"
## use RNA for all DE analysis and plots
DefaultAssay(seurat.object) <- "Spatial"
# seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")

## kernel density estimation 
# https://www.nature.com/articles/s41591-021-01329-2#Sec8
major.t <- c("CD4", "CD8A", "CCR7", "PTPRC", "SELL", "CD27", "CD38", "CD44", "CXCR5", "CD40LG", "CCR7", "FOXP3", "IKZF2")
gd.t<- c("TRGV9", "TRDV2")
MAIT.t <- c("TRAV10", "TRAJ18", "TRAV1-2")
NK.t <- c("CD3G", "FCGR3A", "NCAM1", "NCR1")
TH1.t <- c("IFNG", "TBX21", "TNF")
TH2.t <- c("GATA3", "IL4")
TH17.t <- c("RORC", "IL17A", "IL17F", "IL21")
## https://www.bioconductor.org/packages/release/bioc/vignettes/Nebulosa/inst/doc/nebulosa_seurat.html#1_Overview
  ## https://academic.oup.com/bioinformatics/article-abstract/37/16/2485/6103785
# BiocManager::install("Nebulosa")
library("Nebulosa")
Idents(seurat.object) <- "seurat_clusters"

# SARS.genes <- c("ORF1AB-1-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
p1 <- plot_density(seurat.object, "CD4")
p2 <- FeaturePlot(seurat.object, "CD8A")
p5 <- CombinePlots(plots = list(p1, p2))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("UMAP_KernelDensity_CD4-CD8A.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots
p3 <- plot_density(seurat.object, c("CD4", "CD8A"))
p5 <- p3 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_CD4-CD8A_Joint.pdf", plot = p5, device = "pdf", width = 5, height = 7, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD4", "CD8A"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_major.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("TRGV9", "TRDV2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_gd.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("TRAV10", "TRAJ18"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MAIT.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD3G", "FCGR3A","NCAM1", "NCR1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_NK.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("IFNG", "TBX21", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH1.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 16, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("GATA3", "IL4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH2.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("RORC", "IL17A", "IL17F", "IL21"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH17.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")

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

macrophage.lineage <- c("CD14", "CD68", "CSF1R", "ADGRE1", "SIGLEC1", "CD38", "GPR18", "FPR2", "MYC", "EGR2")
feature.plot <- DotPlot(seurat.object, features = macrophage.lineage)
png(file=paste0("macrophage.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Macrophage Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()


# 65. Aagaard-Tillery KM, Silver R, Dalton J. Immunology of normal pregnancy. Semin Fetal Neonatal Med. Oct 2006;11(5):279-95. doi:10.1016/j.siny.2006.04.003
# 66. Ivashkiv LB. Epigenetic regulation of macrophage polarization and function. Trends Immunol. May 2013;34(5):216-23. doi:10.1016/j.it.2012.11.001
# 67. Murray PJ. Macrophage Polarization. Annu Rev Physiol. 02 10 2017;79:541-566. doi:10.1146/annurev-physiol-022516-034339
# 68. Yao Y, Xu XH, Jin L. Macrophage Polarization in Physiological and Pathological Pregnancy. Front Immunol. 2019;10:792. doi:10.3389/fimmu.2019.00792
# 69. Mues B, Langer D, Zwadlo G, Sorg C. Phenotypic characterization of macrophages in human term placenta. Immunology. Jul 1989;67(3):303-7. 
# 70. Bulmer JN, Johnson PM. Macrophage populations in the human placenta and amniochorion. Clin Exp Immunol. Aug 1984;57(2):393-403. 
# 71. Loegl J, Hiden U, Nussbaumer E, et al. Hofbauer cells of M2a, M2b and M2c polarization may regulate feto-placental angiogenesis. Reproduction. 2016;152(5):447-455. doi:10.1530/REP-16-0159
# 74. Schliefsteiner C, Ibesich S, Wadsack C. Placental Hofbauer Cell Polarization Resists Inflammatory Cues In Vitro. Int J Mol Sci. Jan 22 2020;21(3)doi:10.3390/ijms21030736
# 105.  Ben Amara A, Gorvel L, Baulan K, et al. Placental macrophages are impaired in chorioamnionitis, an infectious pathology of the placenta. J Immunol. Dec 01 2013;191(11):5501-14. doi:10.4049/jimmunol.1300988

macrophage.lineage.markers.2 <- c("CD14", "ITGAM", "CSF1R", "ADGRE1", "CD80", "CD38", "GPR18", "FPR2", "CD86", "TNF", "IL12A", "MYC", "EGR2", "CD163", "MRC1", "CD209", "TGFB1", "IL10", "VEGFA", "CD163", "HLA-DRA", "CD86", "IL6")
macrophage.lineage.markers.2 <- unique(macrophage.lineage.markers.2)

feature.plot <- DotPlot(seurat.object, features = macrophage.lineage.markers.2)
png(file=paste0("macrophage.lineage.markers2-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Macrophage Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

# https://www.nature.com/articles/s41591-021-01329-2#Sec8
t.lineages <- c("CD4", "CD8A", "CCR7", "PTPRC", "SELL", "CD27", "CD38", "CD44", "CXCR5", "CD40LG", "CCR7", "FOXP3", "IKZF2", "TRGV9", "TRDV2", "TRAV10", "TRAJ18", "TRAV1-2", "CD3G", "FCGR3A", "NCAM1", "NCR1", "IFNG", "TBX21", "TNF", "GATA3", "IL4", "RORC", "IL17A", "IL17F", "IL21")
t.lineages <- unique(t.lineages)

# https://www.nature.com/articles/s41591-021-01329-2#Sec8
major.t <- c("CD4", "CD8A", "CCR7", "PTPRC", "SELL", "CD27", "CD38", "CD44", "CXCR5", "CD40LG", "CCR7", "FOXP3", "IKZF2")
gd.t<- c("TRGV9", "TRDV2")
MAIT.t <- c("TRAV10", "TRAJ18", "TRAV1-2")
NK.t <- c("CD3G", "FCGR3A", "NCAM1", "NCR1")
TH1.t <- c("IFNG", "TBX21", "TNF")
TH2.t <- c("GATA3", "IL4")
TH17.t <- c("RORC", "IL17A", "IL17F", "IL21")

#seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=t.lineages)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, slot="counts", assay="Spatial")
ggsave("t.lineages_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, assay="Spatial")
ggsave("t.lineages_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=macrophage.lineage.markers.2)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("macrophage.lineage.markers.2_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, assay="Spatial")
ggsave("macrophage.lineage.markers.2_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=lymphoid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("lymphoid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, assay="Spatial")
ggsave("lymphoid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=myeloid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("myeloid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, assay="Spatial")
ggsave("myeloid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


############## Activation gating strategy https://www.nature.com/articles/s41590-021-01049-2/figures/10
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD3D", "NCAM1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TorNK_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD4", "KLRB1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRB1-CD161.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD4", "CD38"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD38.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD4", "CD69"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD69.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD4", "KLRK1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRK1-NKG2D.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD4", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("FOXP3", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-FOXP3-CD4.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("KIT", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("TNF", "TNFRSF8"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.TNFRSF8.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p4 <- plot_density(seurat.object, c("CD14", "ITGAM"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("SIGLEC1", "FOLR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("SIGLEC1", "DDX3Y"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC-fetalsex_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
	## Placenta-associated maternal macrophages PAMMS (FOLR2 neg, HLA positive); HBCs FOLR2+HLA- ; Thomas et al. Naomi McGovern, JEM 2020
p4 <- plot_density(seurat.object, c("FOLR2", "HLA-DRA"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMMs_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM1_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM2_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("ITGAM", "CD14"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M0.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD80", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD11B", "MYC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD163", "MRC1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2a.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("HLA-DRA", "IL6"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD14", "CD163"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2c.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD14", "CD1C"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD14", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p1 <- plot_density(seurat.object, "THBD")
p2 <- FeaturePlot(seurat.object, "THBD")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "ITGB3")
p2 <- FeaturePlot(seurat.object, "ITGB3")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_Megakaryocyte.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "CD1C")
p2 <- FeaturePlot(seurat.object, "CD1C")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "IL6")
p2 <- FeaturePlot(seurat.object, "IL6")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "TNF")
p2 <- FeaturePlot(seurat.object, "TNF")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "MYC")
p2 <- FeaturePlot(seurat.object, "MYC")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "FOLR2")
p2 <- FeaturePlot(seurat.object, "FOLR2")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_HBCs-FOLR2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "HLA-DRA")
p2 <- FeaturePlot(seurat.object, "HLA-DRA")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_PAMM-HLA-DRA.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "CD33")
p2 <- FeaturePlot(seurat.object, "CD33")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeliodProgenitor-M0.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p4 <- plot_density(seurat.object, c("KIT", "CD33"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_MyeliodProgenitor-M0-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
# seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
setwd("..")
############################################################# ############################################################# #############################################################

dir.create("ZIKV-Veh.gene.list")
setwd("ZIKV-Veh.gene.list")
## ZIKV-Veh Treatment Gene List
Idents(seurat.object) <- "Treatment"
gene.list <- c("CTLA2A", "KAP","PRAP1", "MALAT1", "MT-ND2", "MT-ATP6", "MT-CO1")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "CTLA2A", images=image.list, ncol=4)
p5
ggsave("CTLA2A_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "CTLA2A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CTLA2A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "PRAP1", images=image.list, ncol=4)
p5
ggsave("PRAP1_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "PRAP1", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("PRAP1_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "KAP", images=image.list, ncol=4)
p5
ggsave("KAP_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "KAP", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("KAP_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


p5 <- SpatialPlot(seurat.object, features = "CTLA2A", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object, features = "CTLA2A",)
p5
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = "CTLA2A")
p5
ggsave("CTLA2A_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object) <- "Treatment"
p5 <- VlnPlot(seurat.object, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")


image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "CTLA2A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CTLA2A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

Idents(seurat.object) <- "Treatment"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="Treatment")
dev.off()

seurat.object <- ScaleData(seurat.object, features = gene.list, assay="Spatial")


top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("ZIKV-Veh.gene.list_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("ZIKV-Veh.gene.list_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "orig.ident"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("orig.ident_ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("orig.ident_ZIKV-Veh.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("orig.ident_ZIKV-Veh.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("orig.ident_ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("orig.ident_ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "Type.Profile"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("Type.Profile_ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("Type.Profile_ZIKV-Veh.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Type.Profile_ZIKV-Veh.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Type.Profile_ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Type.Profile_ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "Subtype.Profile"
feature.plot <- DotPlot(seurat.object, features = gene.list)
png(file=paste0("Subtype.Profile_ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, raster = FALSE, size=5) 
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
setwd("..")

############################################################# ############################################################# #############################################################
############################################################# DE Analysis Btwn Treatment #############################################################
############################################################# ############################################################# #############################################################
setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/annotated")



dir.create("Treatment")
setwd("Treatment")
Idents(seurat.object) <- "Treatment"
levels(seurat.object)
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "Treatment")
ggsave("integrated_UMAP-Treatment.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")

Idents(seurat.object) <- "Treatment"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-log2FC.txt", sep="\t")
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
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Treatment")
dev.off()
setwd("..")
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# DE Analysis Btwn Virus #############################################################
############################################################# ############################################################# #############################################################
setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/annotated")

dir.create("Virus")
setwd("Virus")
Idents(seurat.object) <- "Virus"
levels(seurat.object)
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "Virus")
ggsave("integrated_UMAP-Virus.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")

Idents(seurat.object) <- "Virus"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("integrated_top15_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("integrated_top15_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("integrated_top30.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Virus")
dev.off()
setwd("..")
############################################################# ############################################################# #############################################################

############################################################# ############################################################# #############################################################
############################################################# DE Analysis Btwn orig.ident #############################################################
############################################################# ############################################################# #############################################################
setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/annotated")



dir.create("orig.ident")
setwd("orig.ident")
Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident")
ggsave("integrated_UMAP-orig.ident.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")

Idents(seurat.object) <- "orig.ident"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-log2FC.txt", sep="\t")
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
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="SampleID")
dev.off()
setwd("..")
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################

dir.create("GSEA_Treatment")
setwd("GSEA_Treatment")
Idents(seurat.object) <- "Treatment"

## use RNA for all DE analysis and plots
DefaultAssay(seurat.object) <- "Spatial"
# seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")

hypoxia.list <- c("ACKR3", "ADM", "ADORA2B", "AK4", "AKAP12", "ALDOA", "ALDOB", "ALDOC", "AMPD3", "ANGPTL4", "ANKZF1", "ANXA2", "ATF3", "ATP7A", "B3GALT6", "B4GALNT2", "BCAN", "BCL2", "BGN", "BHLHE40", "BNIP3L", "BRS3", "BTG1", "CA12", "CASP6", "CAV1", "CAVIN1", "CAVIN3", "CCN1", "CCN2", "CCN5", "CCNG2", "CDKN1A", "CDKN1B", "CDKN1C", "CHST2", "CHST3", "CITED2", "COL5A1", "CP", "CSRP2", "CXCR4", "DCN", "DDIT3", "DDIT4", "DPYSL4", "DTNA", "DUSP1", "EDN2", "EFNA1", "EFNA3", "EGFR", "ENO1", "ENO2", "ENO3", "ERO1A", "ERRFI1", "ETS1", "EXT1", "F3", "FAM162A", "FBP1", "FOS", "FOSL2", "FOXO3", "GAA", "GALK1", "GAPDH", "GAPDHS", "GBE1", "GCK", "GCNT2", "GLRX", "GPC1", "GPC3", "GPC4", "GPI", "GRHPR", "GYS1", "HAS1", "HDLBP", "HEXA", "HK1", "HK2", "HMOX1", "HOXB9", "HS3ST1", "HSPA5", "IDS", "IER3", "IGFBP1", "IGFBP3", "IL6", "ILVBL", "INHA", "IRS2", "ISG20", "JMJD6", "JUN", "KDELR3", "KDM3A", "KIF5A", "KLF6", "KLF7", "KLHL24", "LALBA", "LARGE1", "LDHA", "LDHC", "LOX", "LXN", "MAFF", "MAP3K1", "MIF", "MT1E", "MT2A", "MXI1", "MYH9", "NAGK", "NCAN", "NDRG1", "NDST1", "NDST2", "NEDD4L", "NFIL3", "NOCT", "NR3C1", "P4HA1", "P4HA2", "PAM", "PCK1", "PDGFB", "PDK1", "PDK3", "PFKFB3", "PFKL", "PFKP", "PGAM2", "PGF", "PGK1", "PGM1", "PGM2", "PHKG1", "PIM1", "PKLR", "PKP1", "PLAC8", "PLAUR", "PLIN2", "PNRC1", "PPARGC1A", "PPFIA4", "PPP1R15A", "PPP1R3C", "PRDX5", "PRKCA", "PYGM", "RBPJ", "RORA", "RRAGD", "S100A4", "SAP30", "SCARB1", "SDC2", "SDC3", "SDC4", "SELENBP1", "SERPINE1", "SIAH2", "SLC25A1", "SLC2A1", "SLC2A3", "SLC2A5", "SLC37A4", "SLC6A6", "SRPX", "STBD1", "STC1", "STC2", "SULT2B1", "TES", "TGFB3", "TGFBI", "TGM2", "TIPARP", "TKTL1", "TMEM45A", "TNFAIP3", "TPBG", "TPD52", "TPI1", "TPST2", "UGP2", "VEGFA", "VHL", "VLDLR", "WSB1", "XPNPEP1", "ZFP36", "ZNF292")
hypoxia.list <- intersect(all.genes, hypoxia.list)

seurat.object <- AddModuleScore(object = seurat.object, features = hypoxia.list, ctrl = 5, name = 'hypoxia.score')
p1 <- FeaturePlot(seurat.object, features = 'hypoxia.score1')
ggsave("FeaturePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'hypoxia.score1')
ggsave("RidgePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "hypoxia.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "hypoxia.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "hypoxia.score1", images=image.list, ncol=4)
p5
ggsave("hypoxia.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


complement.list <- c("ACTN2", "ADAM9", "ADRA2B", "AKAP10", "ANG", "ANXA5", "APOA4", "APOBEC3F", "APOBEC3G", "APOC1", "ATOX1", "BRPF3", "C1QA", "C1QC", "C1R", "C1S", "C2", "C3", "C4BPB", "C9", "CA2", "CALM1", "CALM3", "CASP1", "CASP10", "CASP3", "CASP4", "CASP5", "CASP7", "CASP9", "CBLB", "CCL5", "CD36", "CD40LG", "CD46", "CD55", "CD59", "CDA", "CDH13", "CDK5R1", "CEBPB", "CFB", "CFH", "CLU", "COL4A2", "CP", "CPM", "CPQ", "CR1", "CR2", "CSRP1", "CTSB", "CTSC", "CTSD", "CTSH", "CTSL", "CTSO", "CTSS", "CTSV", "CXCL1", "DGKG", "DGKH", "DOCK10", "DOCK4", "DOCK9", "DPP4", "DUSP5", "DUSP6", "DYRK2", "EHD1", "ERAP2", "F10", "F2", "F3", "F5", "F7", "F8", "FCER1G", "FCN1", "FDX1", "FN1", "FYN", "GATA3", "GCA", "GMFB", "GNAI2", "GNAI3", "GNB2", "GNB4", "GNG2", "GNGT2", "GP1BA", "GP9", "GPD2", "GRB2", "GZMA", "GZMB", "GZMK", "HNF4A", "HPCAL4", "HSPA1A", "HSPA5", "IL6", "IRF1", "IRF2", "IRF7", "ITGAM", "ITIH1", "JAK2", "KCNIP2", "KCNIP3", "KIF2A", "KLK1", "KLKB1", "KYNU", "L3MBTL4", "LAMP2", "LAP3", "LCK", "LCP2", "LGALS3", "LGMN", "LIPA", "LRP1", "LTA4H", "LTF", "LYN", "MAFF", "ME1", "MMP12", "MMP13", "MMP14", "MMP15", "MMP8", "MSRB1", "MT3", "NOTCH4", "OLR1", "PCLO", "PCSK9", "PDGFB", "PDP1", "PFN1", "PHEX", "PIK3CA", "PIK3CG", "PIK3R5", "PIM1", "PLA2G4A", "PLA2G7", "PLAT", "PLAUR", "PLEK", "PLG", "PLSCR1", "PPP2CB", "PPP4C", "PRCP", "PRDM4", "PREP", "PRKCD", "PRSS3", "PRSS36", "PSEN1", "PSMB9", "RABIF", "RAF1", "RASGRP1", "RBSN", "RCE1", "RHOG", "RNF4", "S100A12", "S100A13", "S100A9", "SCG3", "SERPINA1", "SERPINB2", "SERPINC1", "SERPINE1", "SERPING1", "SH2B3", "SIRT6", "SPOCK2", "SRC", "STX4", "TFPI2", "TIMP1", "TIMP2", "TMPRSS6", "TNFAIP3", "USP14", "USP15", "USP16", "USP8", "VCPIP1", "WAS", "XPNPEP1", "ZEB1", "ZFPM2")
complement.list <- intersect(all.genes, complement.list)
seurat.object <- AddModuleScore(object = seurat.object, features = complement.list, ctrl = 5, name = 'complement.score')
p1 <- FeaturePlot(seurat.object, features = 'complement.score1')
ggsave("FeaturePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'complement.score1')
ggsave("RidgePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "complement.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "complement.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "complement.score1", images=image.list, ncol=4)
p5
ggsave("complement.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


# http://www.gsea-msigdb.org/gsea/msigdb/cards/COATES_MACROPHAGE_M1_VS_M2_UP.html
M1 <- c("ABCA9", "ABHD1", "ACSS1", "ACYP2", "ADA", "AIF1", "ALAD", "ARFGEF3", "ARSB", "ATP6V0E2", "AXL", "BCKDHB", "BLOC1S6", "CADM1", "CAP1", "CAPN5", "CBX6", "CD59", "CFH", "CLBA1", "CNRIP1", "COLEC12", "COMT", "CRIM1", "CXCL14", "CXCR4", "DST", "DYNLT1", "EMC1", "ENO2", "FAM124A", "FAM135A", "FAM9B", "FGD2", "FILIP1L", "GALNT11", "GATM", "GDA", "GJA1", "GLO1", "GNB4", "HAUS2", "HDDC3", "HLA-DQA1", "HMGN3", "KCNJ10", "LAMA3", "LCORL", "LYPLAL1", "MAF", "MALAT1", "MARCKSL1", "MARCO", "MSR1", "NAT8L", "NRCAM", "OCEL1", "OGFRL1", "P2RY13", "PIANP", "PIK3AP1", "PLAAT3", "PLBD1", "PLXDC2", "PPP2R5C", "PTGER3", "RAB10", "RAPSN", "RASAL2", "RCBTB2", "RCN1", "RFX3", "RPL14", "SFI1", "SLC35A1", "SLC7A7", "SLCO2B1", "SRD5A3", "TGFBI", "TIFAB", "TM7SF3", "TOR3A", "TTC3", "TUBB2B", "TXNIP", "ZNF727")
M1 <- intersect(all.genes, M1)

seurat.object <- AddModuleScore(object = seurat.object, features = M1, ctrl = 5, name = 'M1.score')
p1 <- FeaturePlot(seurat.object, features = 'M1.score1')
ggsave("FeaturePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'M1.score1')
ggsave("RidgePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "M1.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "M1.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "M1.score1", images=image.list, ncol=4)
p5
ggsave("M1.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

M2<- c("ADAM8", "AK4", "AKR1E2", "ALDOC", "ANG", "ANKRD37", "ANOS1", "ARG1", "ARSK", "ATG4D", "ATP6V0D2", "BNIP3", "BST1", "C1QB", "C1QC", "C1QTNF12", "C5orf34", "CBLB", "CD24", "CD300C", "CD300LD", "CD5L", "CHIA", "COL20A1", "CRIPT", "CTSK", "CTSL", "CYTIP", "EFHD2", "EGLN3", "ERO1A", "F10", "F7", "FAM177A1", "FAM199X", "FAM241A", "FBLIM1", "FLRT2", "GASK1B", "GDF15", "GPRC5B", "HILPDA", "HLA-A", "HLA-B", "HS6ST1", "HSPA1A", "HSPA1B", "IGF2R", "IGHM", "IL18BP", "ITGB3", "LBP", "LIN7C", "LRRC27", "MCOLN3", "MFGE8", "MGST2", "MYO1F", "PADI4", "PDXDC1", "PINK1", "PKDCC", "PRDX2", "PROCR", "PTGER2", "QRICH1", "RGCC", "RIMBP2", "RPGRIP1", "RRAS2", "SCD", "SH2B2", "SLC2A1", "SNHG6", "SOAT1", "THBS1", "TJP2", "TLR1", "TMEM267", "TULP4", "UCHL1", "VEGFA", "XDH")
M2 <- intersect(all.genes, M2)
seurat.object <- AddModuleScore(object = seurat.object, features = M2, ctrl = 5, name = 'M2.score')
p1 <- FeaturePlot(seurat.object, features = 'M2.score1')
ggsave("FeaturePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'M2.score1')
ggsave("RidgePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "M2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "M2.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "M2.score1", images=image.list, ncol=4)
p5
ggsave("M2.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


## http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_PI3K_AKT_MTOR_SIGNALING.html
pi3k <- c("ACACA", "ACTR2", "ACTR3", "ADCY2", "AKT1", "AKT1S1", "AP2M1", "ARF1", "ARHGDIA", "ARPC3", "ATF1", "CAB39", "CAB39L", "CALR", "CAMK4", "CDK1", "CDK2", "CDK4", "CDKN1A", "CDKN1B", "CFL1", "CLTC", "CSNK2B", "CXCR4", "DAPP1", "DDIT3", "DUSP3", "E2F1", "ECSIT", "EGFR", "EIF4E", "FASLG", "FGF17", "FGF22", "FGF6", "GNA14", "GNGT1", "GRB2", "GRK2", "GSK3B", "HRAS", "HSP90B1", "IL2RG", "IL4", "IRAK4", "ITPR2", "LCK", "MAP2K3", "MAP2K6", "MAP3K7", "MAPK1", "MAPK10", "MAPK8", "MAPK9", "MAPKAP1", "MKNK1", "MKNK2", "MYD88", "NCK1", "NFKBIB", "NGF", "NOD1", "PAK4", "PDK1", "PFN1", "PIK3R3", "PIKFYVE", "PIN1", "PITX2", "PLA2G12A", "PLCB1", "PLCG1", "PPP1CA", "PPP2R1B", "PRKAA2", "PRKAG1", "PRKAR2A", "PRKCB", "PTEN", "PTPN11", "RAC1", "RAF1", "RALB", "RIPK1", "RIT1", "RPS6KA1", "RPS6KA3", "RPTOR", "SFN", "SLA", "SLC2A1", "SMAD2", "SQSTM1", "STAT2", "TBK1", "THEM4", "TIAM1", "TNFRSF1A", "TRAF2", "TRIB3", "TSC2", "UBE2D3", "UBE2N", "VAV3", "YWHAB"))
seurat.object <- AddModuleScore(object = seurat.object, features = pi3k, ctrl = 5, name = 'pi3k.score')
p1 <- FeaturePlot(seurat.object, features = 'pi3k.score1')
ggsave("FeaturePlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'pi3k.score1')
ggsave("RidgePlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "pi3k.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "pi3k.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

inflammation <- c("ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "ADM", "ADORA2B", "ADRM1", "AHR", "APLNR", "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1", "BEST1", "BST2", "BTG2", "C3AR1", "C5AR1", "CALCRL", "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7", "CCRL2", "CD14", "CD40", "CD48", "CD55", "CD69", "CD70", "CD82", "CDKN1A", "CHST2", "CLEC5A", "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1", "CXCL10", "CXCL11", "CXCL6", "CXCL8", "CXCL9", "CXCR6", "CYBB", "DCBLD2", "EBI3", "EDN1", "EIF2AK2", "EMP3", "EREG", "F3", "FFAR2", "FPR1", "FZD5", "GABBR1", "GCH1", "GNA15", "GNAI3", "GP1BA", "GPC3", "GPR132", "GPR183", "HAS2", "HBEGF", "HIF1A", "HPN", "HRH1", "ICAM1", "ICAM4", "ICOSLG", "IFITM1", "IFNAR1", "IFNGR2", "IL10", "IL10RA", "IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5", "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2", "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV", "MEP1A", "MET", "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA", "NLRP3", "NMI", "NMUR1", "NOD2", "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2", "PCDH7", "PDE4B", "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2", "PTGER4", "PTGIR", "PTPRE", "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16", "RHOG", "RIPK2", "RNF144B", "ROS1", "RTP4", "SCARF1", "SCN1B", "SELE", "SELENOS", "SELL", "SEMA4D", "SERPINE1", "SGMS2", "SLAMF1", "SLC11A2", "SLC1A2", "SLC28A2", "SLC31A1", "SLC31A2", "SLC4A4", "SLC7A1", "SLC7A2", "SPHK1", "SRI", "STAB1", "TACR1", "TACR3", "TAPBP", "TIMP1", "TLR1", "TLR2", "TLR3", "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10", "TNFSF15", "TNFSF9", "TPBG", "VIP"))
inflammation <- intersect(all.genes, inflammation)

seurat.object <- AddModuleScore(object = seurat.object, features = inflammation, ctrl = 5, name = 'inflammation.score')
p1 <- FeaturePlot(seurat.object, features = 'inflammation.score1')
ggsave("FeaturePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'inflammation.score1')
ggsave("RidgePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "inflammation.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "inflammation.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p5 <- SpatialPlot(seurat.object, features = "inflammation.score1", images=image.list, ncol=4)
p5
ggsave("inflammation.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

ifn.alpha <- c("ADAR", "B2M", "BATF2", "BST2", "C1S", "CASP1", "CASP8", "CCRL2", "CD47", "CD74", "CMPK2", "CMTR1", "CNP", "CSF1", "CXCL10", "CXCL11", "DDX60", "DHX58", "EIF2AK2", "ELF1", "EPSTI1", "GBP2", "GBP4", "GMPR", "HELZ2", "HERC6", "HLA-C", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3", "IL15", "IL4R", "IL7", "IRF1", "IRF2", "IRF7", "IRF9", "ISG15", "ISG20", "LAMP3", "LAP3", "LGALS3BP", "LPAR6", "LY6E", "MOV10", "MVB12A", "MX1", "NCOA7", "NMI", "NUB1", "OAS1", "OASL", "OGFR", "PARP12", "PARP14", "PARP9", "PLSCR1", "PNPT1", "PROCR", "PSMA3", "PSMB8", "PSMB9", "PSME1", "PSME2", "RIPK2", "RNF31", "RSAD2", "RTP4", "SAMD9", "SAMD9L", "SELL", "SLC25A28", "SP110", "STAT2", "TAP1", "TDRD7", "TENT5A", "TMEM140", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TRIM5", "TXNIP", "UBA7", "UBE2L6", "USP18", "WARS1"))
seurat.object <- AddModuleScore(object = seurat.object, features = ifn.alpha, ctrl = 5, name = 'ifn.alpha.score')
p1 <- FeaturePlot(seurat.object, features = 'ifn.alpha.score1')
ggsave("FeaturePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'ifn.alpha.score1')
ggsave("RidgePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "ifn.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "ifn.alpha.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

ifn.beta <- c("ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM", "BST2", "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1", "CMPK2", "CMTR1", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX58", "DDX60", "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "GBP4", "GBP6", "GCH1", "GPR18", "GZMA", "HELZ2", "HERC6", "HIF1A", "HLA-A", "HLA-B", "HLA-DMA", "HLA-DQA1", "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IFNAR2", "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1", "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7", "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1", "METTL7B", "MT2A", "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12", "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP", "PNPT1", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1", "RIPK2", "RNF213", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP", "SERPING1", "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2", "SP110", "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4", "STAT1", "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP", "TDRD7", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8", "VCAM1", "WARS1", "XAF1", "XCL1", "ZBP1", "ZNFX1"))
seurat.object <- AddModuleScore(object = seurat.object, features = ifn.beta, ctrl = 5, name = 'ifn.beta.score')
p1 <- FeaturePlot(seurat.object, features = 'ifn.beta.score1')
ggsave("FeaturePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'ifn.beta.score1')
ggsave("RidgePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "ifn.beta.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "ifn.beta.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tgfb <- c("ACVR1", "APC", "ARID4B", "BCAR3", "BMP2", "BMPR1A", "BMPR2", "CDH1", "CDK9", "CDKN1C", "CTNNB1", "ENG", "FKBP1A", "FNTA", "FURIN", "HDAC1", "HIPK2", "ID1", "ID2", "ID3", "IFNGR2", "JUNB", "KLF10", "LEFTY2", "LTBP2", "MAP3K7", "NCOR2", "NOG", "PMEPA1", "PPM1A", "PPP1CA", "PPP1R15A", "RAB31", "RHOA", "SERPINE1", "SKI", "SKIL", "SLC20A1", "SMAD1", "SMAD3", "SMAD6", "SMAD7", "SMURF1", "SMURF2", "SPTBN1", "TGFB1", "TGFBR1", "TGIF1", "THBS1", "TJP1", "TRIM33", "UBE2D3", "WWTR1", "XIAP"))
seurat.object <- AddModuleScore(object = seurat.object, features = tgfb, ctrl = 5, name = 'tgfb.score')
p1 <- FeaturePlot(seurat.object, features = 'tgfb.score1')
ggsave("FeaturePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'tgfb.score1')
ggsave("RidgePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "tgfb.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "tgfb.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

wnt <- c("ADAM17", "AXIN1", "AXIN2", "CCND2", "CSNK1E", "CTNNB1", "CUL1", "DKK1", "DKK4", "DLL1", "DVL2", "FRAT1", "FZD1", "FZD8", "GNAI1", "HDAC11", "HDAC2", "HDAC5", "HEY1", "HEY2", "JAG1", "JAG2", "KAT2A", "LEF1", "MAML1", "MYC", "NCOR2", "NCSTN", "NKD1", "NOTCH1", "NOTCH4", "NUMB", "PPARD", "PSEN2", "PTCH1", "RBPJ", "SKP2", "TCF7", "TP53", "WNT1", "WNT5B", "WNT6"))
seurat.object <- AddModuleScore(object = seurat.object, features = wnt, ctrl = 5, name = 'wnt.score')
p1 <- FeaturePlot(seurat.object, features = 'wnt.score1')
ggsave("FeaturePlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'wnt.score1')
ggsave("RidgePlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "wnt.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "wnt.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

il2 <- c("ABCB1", "ADAM19", "AGER", "AHCY", "AHNAK", "AHR", "ALCAM", "AMACR", "ANXA4", "APLP1", "ARL4A", "BATF", "BATF3", "BCL2", "BCL2L1", "BHLHE40", "BMP2", "BMPR2", "CA2", "CAPG", "CAPN3", "CASP3", "CCND2", "CCND3", "CCNE1", "CCR4", "CD44", "CD48", "CD79B", "CD81", "CD83", "CD86", "CDC42SE2", "CDC6", "CDCP1", "CDKN1C", "CISH", "CKAP4", "COCH", "COL6A1", "CSF1", "CSF2", "CST7", "CTLA4", "CTSZ", "CXCL10", "CYFIP1", "DCPS", "DENND5A", "DHRS3", "DRC1", "ECM1", "EEF1AKMT1", "EMP1", "ENO3", "ENPP1", "EOMES", "ETFBKMT", "ETV4", "F2RL2", "FAH", "FAM126B", "FGL2", "FLT3LG", "FURIN", "GABARAPL1", "GADD45B", "GALM", "GATA1", "GBP4", "GLIPR2", "GPR65", "GPR83", "GPX4", "GSTO1", "GUCY1B1", "HIPK2", "HK2", "HOPX", "HUWE1", "ICOS", "IFITM3", "IFNGR1", "IGF1R", "IGF2R", "IKZF2", "IKZF4", "IL10", "IL10RA", "IL13", "IL18R1", "IL1R2", "IL1RL1", "IL2RA", "IL2RB", "IL3RA", "IL4R", "IRF4", "IRF6", "IRF8", "ITGA6", "ITGAE", "ITGAV", "ITIH5", "KLF6", "LCLAT1", "LIF", "LRIG1", "LRRC8C", "LTB", "MAFF", "MAP3K8", "MAP6", "MAPKAPK2", "MUC1", "MXD1", "MYC", "MYO1C", "MYO1E", "NCOA3", "NCS1", "NDRG1", "NFIL3", "NFKBIZ", "NOP2", "NRP1", "NT5E", "ODC1", "P2RX4", "P4HA1", "PDCD2L", "PENK", "PHLDA1", "PHTF2", "PIM1", "PLAGL1", "PLEC", "PLIN2", "PLPP1", "PLSCR1", "PNP", "POU2F1", "PRAF2", "PRKCH", "PRNP", "PTCH1", "PTGER2", "PTH1R", "PTRH2", "PUS1", "RABGAP1L", "RGS16", "RHOB", "RHOH", "RNH1", "RORA", "RRAGD", "S100A1", "SCN9A", "SELL", "SELP", "SERPINB6", "SERPINC1", "SH3BGRL2", "SHE", "SLC1A5", "SLC29A2", "SLC2A3", "SLC39A8", "SMPDL3A", "SNX14", "SNX9", "SOCS1", "SOCS2", "SPP1", "SPRED2", "SPRY4", "ST3GAL4", "SWAP70", "SYNGR2", "SYT11", "TGM2", "TIAM1", "TLR7", "TNFRSF18", "TNFRSF1B", "TNFRSF21", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFSF10", "TNFSF11", "TRAF1", "TTC39B", "TWSG1", "UCK2", "UMPS", "WLS", "XBP1"))
seurat.object <- AddModuleScore(object = seurat.object, features = il2, ctrl = 5, name = 'il2.score')
p1 <- FeaturePlot(seurat.object, features = 'il2.score1')
ggsave("FeaturePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'il2.score1')
ggsave("RidgePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "il2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "il2.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

il6 <- c("A2M", "ACVR1B", "ACVRL1", "BAK1", "CBL", "CCL7", "CCR1", "CD14", "CD36", "CD38", "CD44", "CD9", "CNTFR", "CRLF2", "CSF1", "CSF2", "CSF2RA", "CSF2RB", "CSF3R", "CXCL1", "CXCL10", "CXCL11", "CXCL13", "CXCL3", "CXCL9", "DNTT", "EBI3", "FAS", "GRB2", "HAX1", "HMOX1", "IFNAR1", "IFNGR1", "IFNGR2", "IL10RB", "IL12RB1", "IL13RA1", "IL15RA", "IL17RA", "IL17RB", "IL18R1", "IL1B", "IL1R1", "IL1R2", "IL2RA", "IL2RG", "IL3RA", "IL4R", "IL6", "IL6ST", "IL7", "IL9R", "INHBE", "IRF1", "IRF9", "ITGA4", "ITGB3", "JUN", "LEPR", "LTB", "LTBR", "MAP3K8", "MYD88", "OSMR", "PDGFC", "PF4", "PIK3R5", "PIM1", "PLA2G2A", "PTPN1", "PTPN11", "PTPN2", "REG1A", "SOCS1", "SOCS3", "STAM2", "STAT1", "STAT2", "STAT3", "TGFB1", "TLR2", "TNF", "TNFRSF12A", "TNFRSF1A", "TNFRSF1B", "TNFRSF21", "TYK2"))
seurat.object <- AddModuleScore(object = seurat.object, features = il6, ctrl = 5, name = 'il6.score')
p1 <- FeaturePlot(seurat.object, features = 'il6.score1')
ggsave("FeaturePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'il6.score1')
ggsave("RidgePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "il6.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "il6.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

apoptosis <- c("ADD1", "AIFM3", "ANKH", "ANXA1", "APP", "ATF3", "AVPR1A", "BAX", "BCAP31", "BCL10", "BCL2L1", "BCL2L10", "BCL2L11", "BCL2L2", "BGN", "BID", "BIK", "BIRC3", "BMF", "BMP2", "BNIP3L", "BRCA1", "BTG2", "BTG3", "CASP1", "CASP2", "CASP3", "CASP4", "CASP6", "CASP7", "CASP8", "CASP9", "CAV1", "CCNA1", "CCND1", "CCND2", "CD14", "CD2", "CD38", "CD44", "CD69", "CDC25B", "CDK2", "CDKN1A", "CDKN1B", "CFLAR", "CLU", "CREBBP", "CTH", "CTNNB1", "CYLD", "DAP", "DAP3", "DCN", "DDIT3", "DFFA", "DIABLO", "DNAJA1", "DNAJC3", "DNM1L", "DPYD", "EBP", "EGR3", "EMP1", "ENO2", "ERBB2", "ERBB3", "EREG", "ETF1", "F2", "F2R", "FAS", "FASLG", "FDXR", "FEZ1", "GADD45A", "GADD45B", "GCH1", "GNA15", "GPX1", "GPX3", "GPX4", "GSN", "GSR", "GSTM1", "GUCY2D", "H1-0", "HGF", "HMGB2", "HMOX1", "HSPB1", "IER3", "IFITM3", "IFNB1", "IFNGR1", "IGF2R", "IGFBP6", "IL18", "IL1A", "IL1B", "IL6", "IRF1", "ISG20", "JUN", "KRT18", "LEF1", "LGALS3", "LMNA", "LUM", "MADD", "MCL1", "MGMT", "MMP2", "NEDD9", "NEFH", "PAK1", "PDCD4", "PDGFRB", "PEA15", "PLAT", "PLCB2", "PLPPR4", "PMAIP1", "PPP2R5B", "PPP3R1", "PPT1", "PRF1", "PSEN1", "PSEN2", "PTK2", "RARA", "RELA", "RETSAT", "RHOB", "RHOT2", "RNASEL", "ROCK1", "SAT1", "SATB1", "SC5D", "SLC20A1", "SMAD7", "SOD1", "SOD2", "SPTAN1", "SQSTM1", "TAP1", "TGFB2", "TGFBR3", "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF12A", "TNFSF10", "TOP2A", "TSPO", "TXNIP", "VDAC2", "WEE1", "XIAP"))
seurat.object <- AddModuleScore(object = seurat.object, features = apoptosis, ctrl = 5, name = 'apoptosis.score')
p1 <- FeaturePlot(seurat.object, features = 'apoptosis.score1')
ggsave("FeaturePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'apoptosis.score1')
ggsave("RidgePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "apoptosis.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "apoptosis.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tnf.alpha <- c("ABCA1", "ACKR3", "AREG", "ATF3", "ATP2B1", "B4GALT1", "B4GALT5", "BCL2A1", "BCL3", "BCL6", "BHLHE40", "BIRC2", "BIRC3", "BMP2", "BTG1", "BTG2", "BTG3", "CCL2", "CCL20", "CCL4", "CCL5", "CCN1", "CCND1", "CCNL1", "CCRL2", "CD44", "CD69", "CD80", "CD83", "CDKN1A", "CEBPB", "CEBPD", "CFLAR", "CLCF1", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL11", "CXCL2", "CXCL3", "CXCL6", "DDX58", "DENND5A", "DNAJB4", "DRAM1", "DUSP1", "DUSP2", "DUSP4", "DUSP5", "EDN1", "EFNA1", "EGR1", "EGR2", "EGR3", "EHD1", "EIF1", "ETS2", "F2RL1", "F3", "FJX1", "FOS", "FOSB", "FOSL1", "FOSL2", "FUT4", "G0S2", "GADD45A", "GADD45B", "GCH1", "GEM", "GFPT2", "GPR183", "HBEGF", "HES1", "ICAM1", "ICOSLG", "ID2", "IER2", "IER3", "IER5", "IFIH1", "IFIT2", "IFNGR2", "IL12B", "IL15RA", "IL18", "IL1A", "IL1B", "IL23A", "IL6", "IL6ST", "IL7R", "INHBA", "IRF1", "IRS2", "JAG1", "JUN", "JUNB", "KDM6B", "KLF10", "KLF2", "KLF4", "KLF6", "KLF9", "KYNU", "LAMB3", "LDLR", "LIF", "LITAF", "MAFF", "MAP2K3", "MAP3K8", "MARCKS", "MCL1", "MSC", "MXD1", "MYC", "NAMPT", "NFAT5", "NFE2L2", "NFIL3", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NINJ1", "NR4A1", "NR4A2", "NR4A3", "OLR1", "PANX1", "PDE4B", "PDLIM5", "PER1", "PFKFB3", "PHLDA1", "PHLDA2", "PLAU", "PLAUR", "PLEK", "PLK2", "PLPP3", "PMEPA1", "PNRC1", "PPP1R15A", "PTGER4", "PTGS2", "PTPRE", "PTX3", "RCAN1", "REL", "RELA", "RELB", "RHOB", "RIPK2", "RNF19B", "SAT1", "SDC4", "SERPINB2", "SERPINB8", "SERPINE1", "SGK1", "SIK1", "SLC16A6", "SLC2A3", "SLC2A6", "SMAD3", "SNN", "SOCS3", "SOD2", "SPHK1", "SPSB1", "SQSTM1", "STAT5A", "TANK", "TAP1", "TGIF1", "TIPARP", "TLR2", "TNC", "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFRSF9", "TNFSF9", "TNIP1", "TNIP2", "TRAF1", "TRIB1", "TRIP10", "TSC22D1", "TUBB2A", "VEGFA", "YRDC", "ZBTB10", "ZC3H12A", "ZFP36"))
seurat.object <- AddModuleScore(object = seurat.object, features = tnf.alpha, ctrl = 5, name = 'tnf.alpha.score')
p1 <- FeaturePlot(seurat.object, features = 'tnf.alpha.score1')
ggsave("FeaturePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'tnf.alpha.score1')
ggsave("RidgePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "tnf.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "tnf.alpha.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

setwd("..")


dir.create("GSEA_Subtype.Profile")
setwd("GSEA_Subtype.Profile")
Idents(seurat.object) <- "Subtype.Profile"

## use RNA for all DE analysis and plots
DefaultAssay(seurat.object) <- "Spatial"
# seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")

hypoxia.list <- c("ACKR3", "ADM", "ADORA2B", "AK4", "AKAP12", "ALDOA", "ALDOB", "ALDOC", "AMPD3", "ANGPTL4", "ANKZF1", "ANXA2", "ATF3", "ATP7A", "B3GALT6", "B4GALNT2", "BCAN", "BCL2", "BGN", "BHLHE40", "BNIP3L", "BRS3", "BTG1", "CA12", "CASP6", "CAV1", "CAVIN1", "CAVIN3", "CCN1", "CCN2", "CCN5", "CCNG2", "CDKN1A", "CDKN1B", "CDKN1C", "CHST2", "CHST3", "CITED2", "COL5A1", "CP", "CSRP2", "CXCR4", "DCN", "DDIT3", "DDIT4", "DPYSL4", "DTNA", "DUSP1", "EDN2", "EFNA1", "EFNA3", "EGFR", "ENO1", "ENO2", "ENO3", "ERO1A", "ERRFI1", "ETS1", "EXT1", "F3", "FAM162A", "FBP1", "FOS", "FOSL2", "FOXO3", "GAA", "GALK1", "GAPDH", "GAPDHS", "GBE1", "GCK", "GCNT2", "GLRX", "GPC1", "GPC3", "GPC4", "GPI", "GRHPR", "GYS1", "HAS1", "HDLBP", "HEXA", "HK1", "HK2", "HMOX1", "HOXB9", "HS3ST1", "HSPA5", "IDS", "IER3", "IGFBP1", "IGFBP3", "IL6", "ILVBL", "INHA", "IRS2", "ISG20", "JMJD6", "JUN", "KDELR3", "KDM3A", "KIF5A", "KLF6", "KLF7", "KLHL24", "LALBA", "LARGE1", "LDHA", "LDHC", "LOX", "LXN", "MAFF", "MAP3K1", "MIF", "MT1E", "MT2A", "MXI1", "MYH9", "NAGK", "NCAN", "NDRG1", "NDST1", "NDST2", "NEDD4L", "NFIL3", "NOCT", "NR3C1", "P4HA1", "P4HA2", "PAM", "PCK1", "PDGFB", "PDK1", "PDK3", "PFKFB3", "PFKL", "PFKP", "PGAM2", "PGF", "PGK1", "PGM1", "PGM2", "PHKG1", "PIM1", "PKLR", "PKP1", "PLAC8", "PLAUR", "PLIN2", "PNRC1", "PPARGC1A", "PPFIA4", "PPP1R15A", "PPP1R3C", "PRDX5", "PRKCA", "PYGM", "RBPJ", "RORA", "RRAGD", "S100A4", "SAP30", "SCARB1", "SDC2", "SDC3", "SDC4", "SELENBP1", "SERPINE1", "SIAH2", "SLC25A1", "SLC2A1", "SLC2A3", "SLC2A5", "SLC37A4", "SLC6A6", "SRPX", "STBD1", "STC1", "STC2", "SULT2B1", "TES", "TGFB3", "TGFBI", "TGM2", "TIPARP", "TKTL1", "TMEM45A", "TNFAIP3", "TPBG", "TPD52", "TPI1", "TPST2", "UGP2", "VEGFA", "VHL", "VLDLR", "WSB1", "XPNPEP1", "ZFP36", "ZNF292")
hypoxia.list <- intersect(all.genes, hypoxia.list)

seurat.object <- AddModuleScore(object = seurat.object, features = hypoxia.list, ctrl = 5, name = 'hypoxia.score')
p1 <- FeaturePlot(seurat.object, features = 'hypoxia.score1')
ggsave("FeaturePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'hypoxia.score1')
ggsave("RidgePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "hypoxia.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "hypoxia.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "hypoxia.score1", images=image.list, ncol=4)
p5
ggsave("hypoxia.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


complement.list <- c("ACTN2", "ADAM9", "ADRA2B", "AKAP10", "ANG", "ANXA5", "APOA4", "APOBEC3F", "APOBEC3G", "APOC1", "ATOX1", "BRPF3", "C1QA", "C1QC", "C1R", "C1S", "C2", "C3", "C4BPB", "C9", "CA2", "CALM1", "CALM3", "CASP1", "CASP10", "CASP3", "CASP4", "CASP5", "CASP7", "CASP9", "CBLB", "CCL5", "CD36", "CD40LG", "CD46", "CD55", "CD59", "CDA", "CDH13", "CDK5R1", "CEBPB", "CFB", "CFH", "CLU", "COL4A2", "CP", "CPM", "CPQ", "CR1", "CR2", "CSRP1", "CTSB", "CTSC", "CTSD", "CTSH", "CTSL", "CTSO", "CTSS", "CTSV", "CXCL1", "DGKG", "DGKH", "DOCK10", "DOCK4", "DOCK9", "DPP4", "DUSP5", "DUSP6", "DYRK2", "EHD1", "ERAP2", "F10", "F2", "F3", "F5", "F7", "F8", "FCER1G", "FCN1", "FDX1", "FN1", "FYN", "GATA3", "GCA", "GMFB", "GNAI2", "GNAI3", "GNB2", "GNB4", "GNG2", "GNGT2", "GP1BA", "GP9", "GPD2", "GRB2", "GZMA", "GZMB", "GZMK", "HNF4A", "HPCAL4", "HSPA1A", "HSPA5", "IL6", "IRF1", "IRF2", "IRF7", "ITGAM", "ITIH1", "JAK2", "KCNIP2", "KCNIP3", "KIF2A", "KLK1", "KLKB1", "KYNU", "L3MBTL4", "LAMP2", "LAP3", "LCK", "LCP2", "LGALS3", "LGMN", "LIPA", "LRP1", "LTA4H", "LTF", "LYN", "MAFF", "ME1", "MMP12", "MMP13", "MMP14", "MMP15", "MMP8", "MSRB1", "MT3", "NOTCH4", "OLR1", "PCLO", "PCSK9", "PDGFB", "PDP1", "PFN1", "PHEX", "PIK3CA", "PIK3CG", "PIK3R5", "PIM1", "PLA2G4A", "PLA2G7", "PLAT", "PLAUR", "PLEK", "PLG", "PLSCR1", "PPP2CB", "PPP4C", "PRCP", "PRDM4", "PREP", "PRKCD", "PRSS3", "PRSS36", "PSEN1", "PSMB9", "RABIF", "RAF1", "RASGRP1", "RBSN", "RCE1", "RHOG", "RNF4", "S100A12", "S100A13", "S100A9", "SCG3", "SERPINA1", "SERPINB2", "SERPINC1", "SERPINE1", "SERPING1", "SH2B3", "SIRT6", "SPOCK2", "SRC", "STX4", "TFPI2", "TIMP1", "TIMP2", "TMPRSS6", "TNFAIP3", "USP14", "USP15", "USP16", "USP8", "VCPIP1", "WAS", "XPNPEP1", "ZEB1", "ZFPM2")
complement.list <- intersect(all.genes, complement.list)
seurat.object <- AddModuleScore(object = seurat.object, features = complement.list, ctrl = 5, name = 'complement.score')
p1 <- FeaturePlot(seurat.object, features = 'complement.score1')
ggsave("FeaturePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'complement.score1')
ggsave("RidgePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "complement.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "complement.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "complement.score1", images=image.list, ncol=4)
p5
ggsave("complement.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


# http://www.gsea-msigdb.org/gsea/msigdb/cards/COATES_MACROPHAGE_M1_VS_M2_UP.html
M1 <- c("ABCA9", "ABHD1", "ACSS1", "ACYP2", "ADA", "AIF1", "ALAD", "ARFGEF3", "ARSB", "ATP6V0E2", "AXL", "BCKDHB", "BLOC1S6", "CADM1", "CAP1", "CAPN5", "CBX6", "CD59", "CFH", "CLBA1", "CNRIP1", "COLEC12", "COMT", "CRIM1", "CXCL14", "CXCR4", "DST", "DYNLT1", "EMC1", "ENO2", "FAM124A", "FAM135A", "FAM9B", "FGD2", "FILIP1L", "GALNT11", "GATM", "GDA", "GJA1", "GLO1", "GNB4", "HAUS2", "HDDC3", "HLA-DQA1", "HMGN3", "KCNJ10", "LAMA3", "LCORL", "LYPLAL1", "MAF", "MALAT1", "MARCKSL1", "MARCO", "MSR1", "NAT8L", "NRCAM", "OCEL1", "OGFRL1", "P2RY13", "PIANP", "PIK3AP1", "PLAAT3", "PLBD1", "PLXDC2", "PPP2R5C", "PTGER3", "RAB10", "RAPSN", "RASAL2", "RCBTB2", "RCN1", "RFX3", "RPL14", "SFI1", "SLC35A1", "SLC7A7", "SLCO2B1", "SRD5A3", "TGFBI", "TIFAB", "TM7SF3", "TOR3A", "TTC3", "TUBB2B", "TXNIP", "ZNF727")
M1 <- intersect(all.genes, M1)

seurat.object <- AddModuleScore(object = seurat.object, features = M1, ctrl = 5, name = 'M1.score')
p1 <- FeaturePlot(seurat.object, features = 'M1.score1')
ggsave("FeaturePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'M1.score1')
ggsave("RidgePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "M1.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "M1.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "M1.score1", images=image.list, ncol=4)
p5
ggsave("M1.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

M2<- c("ADAM8", "AK4", "AKR1E2", "ALDOC", "ANG", "ANKRD37", "ANOS1", "ARG1", "ARSK", "ATG4D", "ATP6V0D2", "BNIP3", "BST1", "C1QB", "C1QC", "C1QTNF12", "C5orf34", "CBLB", "CD24", "CD300C", "CD300LD", "CD5L", "CHIA", "COL20A1", "CRIPT", "CTSK", "CTSL", "CYTIP", "EFHD2", "EGLN3", "ERO1A", "F10", "F7", "FAM177A1", "FAM199X", "FAM241A", "FBLIM1", "FLRT2", "GASK1B", "GDF15", "GPRC5B", "HILPDA", "HLA-A", "HLA-B", "HS6ST1", "HSPA1A", "HSPA1B", "IGF2R", "IGHM", "IL18BP", "ITGB3", "LBP", "LIN7C", "LRRC27", "MCOLN3", "MFGE8", "MGST2", "MYO1F", "PADI4", "PDXDC1", "PINK1", "PKDCC", "PRDX2", "PROCR", "PTGER2", "QRICH1", "RGCC", "RIMBP2", "RPGRIP1", "RRAS2", "SCD", "SH2B2", "SLC2A1", "SNHG6", "SOAT1", "THBS1", "TJP2", "TLR1", "TMEM267", "TULP4", "UCHL1", "VEGFA", "XDH")
M2 <- intersect(all.genes, M2)
seurat.object <- AddModuleScore(object = seurat.object, features = M2, ctrl = 5, name = 'M2.score')
p1 <- FeaturePlot(seurat.object, features = 'M2.score1')
ggsave("FeaturePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'M2.score1')
ggsave("RidgePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "M2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "M2.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "M2.score1", images=image.list, ncol=4)
p5
ggsave("M2.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


## http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_PI3K_AKT_MTOR_SIGNALING.html
pi3k <- c("ACACA", "ACTR2", "ACTR3", "ADCY2", "AKT1", "AKT1S1", "AP2M1", "ARF1", "ARHGDIA", "ARPC3", "ATF1", "CAB39", "CAB39L", "CALR", "CAMK4", "CDK1", "CDK2", "CDK4", "CDKN1A", "CDKN1B", "CFL1", "CLTC", "CSNK2B", "CXCR4", "DAPP1", "DDIT3", "DUSP3", "E2F1", "ECSIT", "EGFR", "EIF4E", "FASLG", "FGF17", "FGF22", "FGF6", "GNA14", "GNGT1", "GRB2", "GRK2", "GSK3B", "HRAS", "HSP90B1", "IL2RG", "IL4", "IRAK4", "ITPR2", "LCK", "MAP2K3", "MAP2K6", "MAP3K7", "MAPK1", "MAPK10", "MAPK8", "MAPK9", "MAPKAP1", "MKNK1", "MKNK2", "MYD88", "NCK1", "NFKBIB", "NGF", "NOD1", "PAK4", "PDK1", "PFN1", "PIK3R3", "PIKFYVE", "PIN1", "PITX2", "PLA2G12A", "PLCB1", "PLCG1", "PPP1CA", "PPP2R1B", "PRKAA2", "PRKAG1", "PRKAR2A", "PRKCB", "PTEN", "PTPN11", "RAC1", "RAF1", "RALB", "RIPK1", "RIT1", "RPS6KA1", "RPS6KA3", "RPTOR", "SFN", "SLA", "SLC2A1", "SMAD2", "SQSTM1", "STAT2", "TBK1", "THEM4", "TIAM1", "TNFRSF1A", "TRAF2", "TRIB3", "TSC2", "UBE2D3", "UBE2N", "VAV3", "YWHAB"))
seurat.object <- AddModuleScore(object = seurat.object, features = pi3k, ctrl = 5, name = 'pi3k.score')
p1 <- FeaturePlot(seurat.object, features = 'pi3k.score1')
ggsave("FeaturePlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'pi3k.score1')
ggsave("RidgePlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "pi3k.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "pi3k.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

inflammation <- c("ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "ADM", "ADORA2B", "ADRM1", "AHR", "APLNR", "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1", "BEST1", "BST2", "BTG2", "C3AR1", "C5AR1", "CALCRL", "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7", "CCRL2", "CD14", "CD40", "CD48", "CD55", "CD69", "CD70", "CD82", "CDKN1A", "CHST2", "CLEC5A", "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1", "CXCL10", "CXCL11", "CXCL6", "CXCL8", "CXCL9", "CXCR6", "CYBB", "DCBLD2", "EBI3", "EDN1", "EIF2AK2", "EMP3", "EREG", "F3", "FFAR2", "FPR1", "FZD5", "GABBR1", "GCH1", "GNA15", "GNAI3", "GP1BA", "GPC3", "GPR132", "GPR183", "HAS2", "HBEGF", "HIF1A", "HPN", "HRH1", "ICAM1", "ICAM4", "ICOSLG", "IFITM1", "IFNAR1", "IFNGR2", "IL10", "IL10RA", "IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5", "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2", "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV", "MEP1A", "MET", "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA", "NLRP3", "NMI", "NMUR1", "NOD2", "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2", "PCDH7", "PDE4B", "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2", "PTGER4", "PTGIR", "PTPRE", "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16", "RHOG", "RIPK2", "RNF144B", "ROS1", "RTP4", "SCARF1", "SCN1B", "SELE", "SELENOS", "SELL", "SEMA4D", "SERPINE1", "SGMS2", "SLAMF1", "SLC11A2", "SLC1A2", "SLC28A2", "SLC31A1", "SLC31A2", "SLC4A4", "SLC7A1", "SLC7A2", "SPHK1", "SRI", "STAB1", "TACR1", "TACR3", "TAPBP", "TIMP1", "TLR1", "TLR2", "TLR3", "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10", "TNFSF15", "TNFSF9", "TPBG", "VIP"))
inflammation <- intersect(all.genes, inflammation)

seurat.object <- AddModuleScore(object = seurat.object, features = inflammation, ctrl = 5, name = 'inflammation.score')
p1 <- FeaturePlot(seurat.object, features = 'inflammation.score1')
ggsave("FeaturePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'inflammation.score1')
ggsave("RidgePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "inflammation.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "inflammation.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p5 <- SpatialPlot(seurat.object, features = "inflammation.score1", images=image.list, ncol=4)
p5
ggsave("inflammation.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

ifn.alpha <- c("ADAR", "B2M", "BATF2", "BST2", "C1S", "CASP1", "CASP8", "CCRL2", "CD47", "CD74", "CMPK2", "CMTR1", "CNP", "CSF1", "CXCL10", "CXCL11", "DDX60", "DHX58", "EIF2AK2", "ELF1", "EPSTI1", "GBP2", "GBP4", "GMPR", "HELZ2", "HERC6", "HLA-C", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3", "IL15", "IL4R", "IL7", "IRF1", "IRF2", "IRF7", "IRF9", "ISG15", "ISG20", "LAMP3", "LAP3", "LGALS3BP", "LPAR6", "LY6E", "MOV10", "MVB12A", "MX1", "NCOA7", "NMI", "NUB1", "OAS1", "OASL", "OGFR", "PARP12", "PARP14", "PARP9", "PLSCR1", "PNPT1", "PROCR", "PSMA3", "PSMB8", "PSMB9", "PSME1", "PSME2", "RIPK2", "RNF31", "RSAD2", "RTP4", "SAMD9", "SAMD9L", "SELL", "SLC25A28", "SP110", "STAT2", "TAP1", "TDRD7", "TENT5A", "TMEM140", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TRIM5", "TXNIP", "UBA7", "UBE2L6", "USP18", "WARS1"))
seurat.object <- AddModuleScore(object = seurat.object, features = ifn.alpha, ctrl = 5, name = 'ifn.alpha.score')
p1 <- FeaturePlot(seurat.object, features = 'ifn.alpha.score1')
ggsave("FeaturePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'ifn.alpha.score1')
ggsave("RidgePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "ifn.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "ifn.alpha.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

ifn.beta <- c("ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM", "BST2", "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1", "CMPK2", "CMTR1", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX58", "DDX60", "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "GBP4", "GBP6", "GCH1", "GPR18", "GZMA", "HELZ2", "HERC6", "HIF1A", "HLA-A", "HLA-B", "HLA-DMA", "HLA-DQA1", "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IFNAR2", "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1", "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7", "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1", "METTL7B", "MT2A", "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12", "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP", "PNPT1", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1", "RIPK2", "RNF213", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP", "SERPING1", "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2", "SP110", "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4", "STAT1", "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP", "TDRD7", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8", "VCAM1", "WARS1", "XAF1", "XCL1", "ZBP1", "ZNFX1"))
seurat.object <- AddModuleScore(object = seurat.object, features = ifn.beta, ctrl = 5, name = 'ifn.beta.score')
p1 <- FeaturePlot(seurat.object, features = 'ifn.beta.score1')
ggsave("FeaturePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'ifn.beta.score1')
ggsave("RidgePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "ifn.beta.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "ifn.beta.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tgfb <- c("ACVR1", "APC", "ARID4B", "BCAR3", "BMP2", "BMPR1A", "BMPR2", "CDH1", "CDK9", "CDKN1C", "CTNNB1", "ENG", "FKBP1A", "FNTA", "FURIN", "HDAC1", "HIPK2", "ID1", "ID2", "ID3", "IFNGR2", "JUNB", "KLF10", "LEFTY2", "LTBP2", "MAP3K7", "NCOR2", "NOG", "PMEPA1", "PPM1A", "PPP1CA", "PPP1R15A", "RAB31", "RHOA", "SERPINE1", "SKI", "SKIL", "SLC20A1", "SMAD1", "SMAD3", "SMAD6", "SMAD7", "SMURF1", "SMURF2", "SPTBN1", "TGFB1", "TGFBR1", "TGIF1", "THBS1", "TJP1", "TRIM33", "UBE2D3", "WWTR1", "XIAP"))
seurat.object <- AddModuleScore(object = seurat.object, features = tgfb, ctrl = 5, name = 'tgfb.score')
p1 <- FeaturePlot(seurat.object, features = 'tgfb.score1')
ggsave("FeaturePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'tgfb.score1')
ggsave("RidgePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "tgfb.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "tgfb.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

wnt <- c("ADAM17", "AXIN1", "AXIN2", "CCND2", "CSNK1E", "CTNNB1", "CUL1", "DKK1", "DKK4", "DLL1", "DVL2", "FRAT1", "FZD1", "FZD8", "GNAI1", "HDAC11", "HDAC2", "HDAC5", "HEY1", "HEY2", "JAG1", "JAG2", "KAT2A", "LEF1", "MAML1", "MYC", "NCOR2", "NCSTN", "NKD1", "NOTCH1", "NOTCH4", "NUMB", "PPARD", "PSEN2", "PTCH1", "RBPJ", "SKP2", "TCF7", "TP53", "WNT1", "WNT5B", "WNT6"))
seurat.object <- AddModuleScore(object = seurat.object, features = wnt, ctrl = 5, name = 'wnt.score')
p1 <- FeaturePlot(seurat.object, features = 'wnt.score1')
ggsave("FeaturePlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'wnt.score1')
ggsave("RidgePlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "wnt.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "wnt.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

il2 <- c("ABCB1", "ADAM19", "AGER", "AHCY", "AHNAK", "AHR", "ALCAM", "AMACR", "ANXA4", "APLP1", "ARL4A", "BATF", "BATF3", "BCL2", "BCL2L1", "BHLHE40", "BMP2", "BMPR2", "CA2", "CAPG", "CAPN3", "CASP3", "CCND2", "CCND3", "CCNE1", "CCR4", "CD44", "CD48", "CD79B", "CD81", "CD83", "CD86", "CDC42SE2", "CDC6", "CDCP1", "CDKN1C", "CISH", "CKAP4", "COCH", "COL6A1", "CSF1", "CSF2", "CST7", "CTLA4", "CTSZ", "CXCL10", "CYFIP1", "DCPS", "DENND5A", "DHRS3", "DRC1", "ECM1", "EEF1AKMT1", "EMP1", "ENO3", "ENPP1", "EOMES", "ETFBKMT", "ETV4", "F2RL2", "FAH", "FAM126B", "FGL2", "FLT3LG", "FURIN", "GABARAPL1", "GADD45B", "GALM", "GATA1", "GBP4", "GLIPR2", "GPR65", "GPR83", "GPX4", "GSTO1", "GUCY1B1", "HIPK2", "HK2", "HOPX", "HUWE1", "ICOS", "IFITM3", "IFNGR1", "IGF1R", "IGF2R", "IKZF2", "IKZF4", "IL10", "IL10RA", "IL13", "IL18R1", "IL1R2", "IL1RL1", "IL2RA", "IL2RB", "IL3RA", "IL4R", "IRF4", "IRF6", "IRF8", "ITGA6", "ITGAE", "ITGAV", "ITIH5", "KLF6", "LCLAT1", "LIF", "LRIG1", "LRRC8C", "LTB", "MAFF", "MAP3K8", "MAP6", "MAPKAPK2", "MUC1", "MXD1", "MYC", "MYO1C", "MYO1E", "NCOA3", "NCS1", "NDRG1", "NFIL3", "NFKBIZ", "NOP2", "NRP1", "NT5E", "ODC1", "P2RX4", "P4HA1", "PDCD2L", "PENK", "PHLDA1", "PHTF2", "PIM1", "PLAGL1", "PLEC", "PLIN2", "PLPP1", "PLSCR1", "PNP", "POU2F1", "PRAF2", "PRKCH", "PRNP", "PTCH1", "PTGER2", "PTH1R", "PTRH2", "PUS1", "RABGAP1L", "RGS16", "RHOB", "RHOH", "RNH1", "RORA", "RRAGD", "S100A1", "SCN9A", "SELL", "SELP", "SERPINB6", "SERPINC1", "SH3BGRL2", "SHE", "SLC1A5", "SLC29A2", "SLC2A3", "SLC39A8", "SMPDL3A", "SNX14", "SNX9", "SOCS1", "SOCS2", "SPP1", "SPRED2", "SPRY4", "ST3GAL4", "SWAP70", "SYNGR2", "SYT11", "TGM2", "TIAM1", "TLR7", "TNFRSF18", "TNFRSF1B", "TNFRSF21", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFSF10", "TNFSF11", "TRAF1", "TTC39B", "TWSG1", "UCK2", "UMPS", "WLS", "XBP1"))
seurat.object <- AddModuleScore(object = seurat.object, features = il2, ctrl = 5, name = 'il2.score')
p1 <- FeaturePlot(seurat.object, features = 'il2.score1')
ggsave("FeaturePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'il2.score1')
ggsave("RidgePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "il2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "il2.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

il6 <- c("A2M", "ACVR1B", "ACVRL1", "BAK1", "CBL", "CCL7", "CCR1", "CD14", "CD36", "CD38", "CD44", "CD9", "CNTFR", "CRLF2", "CSF1", "CSF2", "CSF2RA", "CSF2RB", "CSF3R", "CXCL1", "CXCL10", "CXCL11", "CXCL13", "CXCL3", "CXCL9", "DNTT", "EBI3", "FAS", "GRB2", "HAX1", "HMOX1", "IFNAR1", "IFNGR1", "IFNGR2", "IL10RB", "IL12RB1", "IL13RA1", "IL15RA", "IL17RA", "IL17RB", "IL18R1", "IL1B", "IL1R1", "IL1R2", "IL2RA", "IL2RG", "IL3RA", "IL4R", "IL6", "IL6ST", "IL7", "IL9R", "INHBE", "IRF1", "IRF9", "ITGA4", "ITGB3", "JUN", "LEPR", "LTB", "LTBR", "MAP3K8", "MYD88", "OSMR", "PDGFC", "PF4", "PIK3R5", "PIM1", "PLA2G2A", "PTPN1", "PTPN11", "PTPN2", "REG1A", "SOCS1", "SOCS3", "STAM2", "STAT1", "STAT2", "STAT3", "TGFB1", "TLR2", "TNF", "TNFRSF12A", "TNFRSF1A", "TNFRSF1B", "TNFRSF21", "TYK2"))
seurat.object <- AddModuleScore(object = seurat.object, features = il6, ctrl = 5, name = 'il6.score')
p1 <- FeaturePlot(seurat.object, features = 'il6.score1')
ggsave("FeaturePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'il6.score1')
ggsave("RidgePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "il6.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "il6.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

apoptosis <- c("ADD1", "AIFM3", "ANKH", "ANXA1", "APP", "ATF3", "AVPR1A", "BAX", "BCAP31", "BCL10", "BCL2L1", "BCL2L10", "BCL2L11", "BCL2L2", "BGN", "BID", "BIK", "BIRC3", "BMF", "BMP2", "BNIP3L", "BRCA1", "BTG2", "BTG3", "CASP1", "CASP2", "CASP3", "CASP4", "CASP6", "CASP7", "CASP8", "CASP9", "CAV1", "CCNA1", "CCND1", "CCND2", "CD14", "CD2", "CD38", "CD44", "CD69", "CDC25B", "CDK2", "CDKN1A", "CDKN1B", "CFLAR", "CLU", "CREBBP", "CTH", "CTNNB1", "CYLD", "DAP", "DAP3", "DCN", "DDIT3", "DFFA", "DIABLO", "DNAJA1", "DNAJC3", "DNM1L", "DPYD", "EBP", "EGR3", "EMP1", "ENO2", "ERBB2", "ERBB3", "EREG", "ETF1", "F2", "F2R", "FAS", "FASLG", "FDXR", "FEZ1", "GADD45A", "GADD45B", "GCH1", "GNA15", "GPX1", "GPX3", "GPX4", "GSN", "GSR", "GSTM1", "GUCY2D", "H1-0", "HGF", "HMGB2", "HMOX1", "HSPB1", "IER3", "IFITM3", "IFNB1", "IFNGR1", "IGF2R", "IGFBP6", "IL18", "IL1A", "IL1B", "IL6", "IRF1", "ISG20", "JUN", "KRT18", "LEF1", "LGALS3", "LMNA", "LUM", "MADD", "MCL1", "MGMT", "MMP2", "NEDD9", "NEFH", "PAK1", "PDCD4", "PDGFRB", "PEA15", "PLAT", "PLCB2", "PLPPR4", "PMAIP1", "PPP2R5B", "PPP3R1", "PPT1", "PRF1", "PSEN1", "PSEN2", "PTK2", "RARA", "RELA", "RETSAT", "RHOB", "RHOT2", "RNASEL", "ROCK1", "SAT1", "SATB1", "SC5D", "SLC20A1", "SMAD7", "SOD1", "SOD2", "SPTAN1", "SQSTM1", "TAP1", "TGFB2", "TGFBR3", "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF12A", "TNFSF10", "TOP2A", "TSPO", "TXNIP", "VDAC2", "WEE1", "XIAP"))
seurat.object <- AddModuleScore(object = seurat.object, features = apoptosis, ctrl = 5, name = 'apoptosis.score')
p1 <- FeaturePlot(seurat.object, features = 'apoptosis.score1')
ggsave("FeaturePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'apoptosis.score1')
ggsave("RidgePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "apoptosis.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "apoptosis.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tnf.alpha <- c("ABCA1", "ACKR3", "AREG", "ATF3", "ATP2B1", "B4GALT1", "B4GALT5", "BCL2A1", "BCL3", "BCL6", "BHLHE40", "BIRC2", "BIRC3", "BMP2", "BTG1", "BTG2", "BTG3", "CCL2", "CCL20", "CCL4", "CCL5", "CCN1", "CCND1", "CCNL1", "CCRL2", "CD44", "CD69", "CD80", "CD83", "CDKN1A", "CEBPB", "CEBPD", "CFLAR", "CLCF1", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL11", "CXCL2", "CXCL3", "CXCL6", "DDX58", "DENND5A", "DNAJB4", "DRAM1", "DUSP1", "DUSP2", "DUSP4", "DUSP5", "EDN1", "EFNA1", "EGR1", "EGR2", "EGR3", "EHD1", "EIF1", "ETS2", "F2RL1", "F3", "FJX1", "FOS", "FOSB", "FOSL1", "FOSL2", "FUT4", "G0S2", "GADD45A", "GADD45B", "GCH1", "GEM", "GFPT2", "GPR183", "HBEGF", "HES1", "ICAM1", "ICOSLG", "ID2", "IER2", "IER3", "IER5", "IFIH1", "IFIT2", "IFNGR2", "IL12B", "IL15RA", "IL18", "IL1A", "IL1B", "IL23A", "IL6", "IL6ST", "IL7R", "INHBA", "IRF1", "IRS2", "JAG1", "JUN", "JUNB", "KDM6B", "KLF10", "KLF2", "KLF4", "KLF6", "KLF9", "KYNU", "LAMB3", "LDLR", "LIF", "LITAF", "MAFF", "MAP2K3", "MAP3K8", "MARCKS", "MCL1", "MSC", "MXD1", "MYC", "NAMPT", "NFAT5", "NFE2L2", "NFIL3", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NINJ1", "NR4A1", "NR4A2", "NR4A3", "OLR1", "PANX1", "PDE4B", "PDLIM5", "PER1", "PFKFB3", "PHLDA1", "PHLDA2", "PLAU", "PLAUR", "PLEK", "PLK2", "PLPP3", "PMEPA1", "PNRC1", "PPP1R15A", "PTGER4", "PTGS2", "PTPRE", "PTX3", "RCAN1", "REL", "RELA", "RELB", "RHOB", "RIPK2", "RNF19B", "SAT1", "SDC4", "SERPINB2", "SERPINB8", "SERPINE1", "SGK1", "SIK1", "SLC16A6", "SLC2A3", "SLC2A6", "SMAD3", "SNN", "SOCS3", "SOD2", "SPHK1", "SPSB1", "SQSTM1", "STAT5A", "TANK", "TAP1", "TGIF1", "TIPARP", "TLR2", "TNC", "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFRSF9", "TNFSF9", "TNIP1", "TNIP2", "TRAF1", "TRIB1", "TRIP10", "TSC22D1", "TUBB2A", "VEGFA", "YRDC", "ZBTB10", "ZC3H12A", "ZFP36"))
seurat.object <- AddModuleScore(object = seurat.object, features = tnf.alpha, ctrl = 5, name = 'tnf.alpha.score')
p1 <- FeaturePlot(seurat.object, features = 'tnf.alpha.score1')
ggsave("FeaturePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'tnf.alpha.score1')
ggsave("RidgePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "tnf.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "tnf.alpha.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

setwd("..")

dir.create("GSEA_Type.Profile")
setwd("GSEA_Type.Profile")
Idents(seurat.object) <- "Type.Profile"

## use RNA for all DE analysis and plots
DefaultAssay(seurat.object) <- "Spatial"
# seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")
hypoxia.list <- c("ACKR3", "ADM", "ADORA2B", "AK4", "AKAP12", "ALDOA", "ALDOB", "ALDOC", "AMPD3", "ANGPTL4", "ANKZF1", "ANXA2", "ATF3", "ATP7A", "B3GALT6", "B4GALNT2", "BCAN", "BCL2", "BGN", "BHLHE40", "BNIP3L", "BRS3", "BTG1", "CA12", "CASP6", "CAV1", "CAVIN1", "CAVIN3", "CCN1", "CCN2", "CCN5", "CCNG2", "CDKN1A", "CDKN1B", "CDKN1C", "CHST2", "CHST3", "CITED2", "COL5A1", "CP", "CSRP2", "CXCR4", "DCN", "DDIT3", "DDIT4", "DPYSL4", "DTNA", "DUSP1", "EDN2", "EFNA1", "EFNA3", "EGFR", "ENO1", "ENO2", "ENO3", "ERO1A", "ERRFI1", "ETS1", "EXT1", "F3", "FAM162A", "FBP1", "FOS", "FOSL2", "FOXO3", "GAA", "GALK1", "GAPDH", "GAPDHS", "GBE1", "GCK", "GCNT2", "GLRX", "GPC1", "GPC3", "GPC4", "GPI", "GRHPR", "GYS1", "HAS1", "HDLBP", "HEXA", "HK1", "HK2", "HMOX1", "HOXB9", "HS3ST1", "HSPA5", "IDS", "IER3", "IGFBP1", "IGFBP3", "IL6", "ILVBL", "INHA", "IRS2", "ISG20", "JMJD6", "JUN", "KDELR3", "KDM3A", "KIF5A", "KLF6", "KLF7", "KLHL24", "LALBA", "LARGE1", "LDHA", "LDHC", "LOX", "LXN", "MAFF", "MAP3K1", "MIF", "MT1E", "MT2A", "MXI1", "MYH9", "NAGK", "NCAN", "NDRG1", "NDST1", "NDST2", "NEDD4L", "NFIL3", "NOCT", "NR3C1", "P4HA1", "P4HA2", "PAM", "PCK1", "PDGFB", "PDK1", "PDK3", "PFKFB3", "PFKL", "PFKP", "PGAM2", "PGF", "PGK1", "PGM1", "PGM2", "PHKG1", "PIM1", "PKLR", "PKP1", "PLAC8", "PLAUR", "PLIN2", "PNRC1", "PPARGC1A", "PPFIA4", "PPP1R15A", "PPP1R3C", "PRDX5", "PRKCA", "PYGM", "RBPJ", "RORA", "RRAGD", "S100A4", "SAP30", "SCARB1", "SDC2", "SDC3", "SDC4", "SELENBP1", "SERPINE1", "SIAH2", "SLC25A1", "SLC2A1", "SLC2A3", "SLC2A5", "SLC37A4", "SLC6A6", "SRPX", "STBD1", "STC1", "STC2", "SULT2B1", "TES", "TGFB3", "TGFBI", "TGM2", "TIPARP", "TKTL1", "TMEM45A", "TNFAIP3", "TPBG", "TPD52", "TPI1", "TPST2", "UGP2", "VEGFA", "VHL", "VLDLR", "WSB1", "XPNPEP1", "ZFP36", "ZNF292")
hypoxia.list <- intersect(all.genes, hypoxia.list)

seurat.object <- AddModuleScore(object = seurat.object, features = hypoxia.list, ctrl = 5, name = 'hypoxia.score')
p1 <- FeaturePlot(seurat.object, features = 'hypoxia.score1')
ggsave("FeaturePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'hypoxia.score1')
ggsave("RidgePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "hypoxia.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "hypoxia.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "hypoxia.score1", images=image.list, ncol=4)
p5
ggsave("hypoxia.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


complement.list <- c("ACTN2", "ADAM9", "ADRA2B", "AKAP10", "ANG", "ANXA5", "APOA4", "APOBEC3F", "APOBEC3G", "APOC1", "ATOX1", "BRPF3", "C1QA", "C1QC", "C1R", "C1S", "C2", "C3", "C4BPB", "C9", "CA2", "CALM1", "CALM3", "CASP1", "CASP10", "CASP3", "CASP4", "CASP5", "CASP7", "CASP9", "CBLB", "CCL5", "CD36", "CD40LG", "CD46", "CD55", "CD59", "CDA", "CDH13", "CDK5R1", "CEBPB", "CFB", "CFH", "CLU", "COL4A2", "CP", "CPM", "CPQ", "CR1", "CR2", "CSRP1", "CTSB", "CTSC", "CTSD", "CTSH", "CTSL", "CTSO", "CTSS", "CTSV", "CXCL1", "DGKG", "DGKH", "DOCK10", "DOCK4", "DOCK9", "DPP4", "DUSP5", "DUSP6", "DYRK2", "EHD1", "ERAP2", "F10", "F2", "F3", "F5", "F7", "F8", "FCER1G", "FCN1", "FDX1", "FN1", "FYN", "GATA3", "GCA", "GMFB", "GNAI2", "GNAI3", "GNB2", "GNB4", "GNG2", "GNGT2", "GP1BA", "GP9", "GPD2", "GRB2", "GZMA", "GZMB", "GZMK", "HNF4A", "HPCAL4", "HSPA1A", "HSPA5", "IL6", "IRF1", "IRF2", "IRF7", "ITGAM", "ITIH1", "JAK2", "KCNIP2", "KCNIP3", "KIF2A", "KLK1", "KLKB1", "KYNU", "L3MBTL4", "LAMP2", "LAP3", "LCK", "LCP2", "LGALS3", "LGMN", "LIPA", "LRP1", "LTA4H", "LTF", "LYN", "MAFF", "ME1", "MMP12", "MMP13", "MMP14", "MMP15", "MMP8", "MSRB1", "MT3", "NOTCH4", "OLR1", "PCLO", "PCSK9", "PDGFB", "PDP1", "PFN1", "PHEX", "PIK3CA", "PIK3CG", "PIK3R5", "PIM1", "PLA2G4A", "PLA2G7", "PLAT", "PLAUR", "PLEK", "PLG", "PLSCR1", "PPP2CB", "PPP4C", "PRCP", "PRDM4", "PREP", "PRKCD", "PRSS3", "PRSS36", "PSEN1", "PSMB9", "RABIF", "RAF1", "RASGRP1", "RBSN", "RCE1", "RHOG", "RNF4", "S100A12", "S100A13", "S100A9", "SCG3", "SERPINA1", "SERPINB2", "SERPINC1", "SERPINE1", "SERPING1", "SH2B3", "SIRT6", "SPOCK2", "SRC", "STX4", "TFPI2", "TIMP1", "TIMP2", "TMPRSS6", "TNFAIP3", "USP14", "USP15", "USP16", "USP8", "VCPIP1", "WAS", "XPNPEP1", "ZEB1", "ZFPM2")
complement.list <- intersect(all.genes, complement.list)
seurat.object <- AddModuleScore(object = seurat.object, features = complement.list, ctrl = 5, name = 'complement.score')
p1 <- FeaturePlot(seurat.object, features = 'complement.score1')
ggsave("FeaturePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'complement.score1')
ggsave("RidgePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "complement.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "complement.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "complement.score1", images=image.list, ncol=4)
p5
ggsave("complement.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


# http://www.gsea-msigdb.org/gsea/msigdb/cards/COATES_MACROPHAGE_M1_VS_M2_UP.html
M1 <- c("ABCA9", "ABHD1", "ACSS1", "ACYP2", "ADA", "AIF1", "ALAD", "ARFGEF3", "ARSB", "ATP6V0E2", "AXL", "BCKDHB", "BLOC1S6", "CADM1", "CAP1", "CAPN5", "CBX6", "CD59", "CFH", "CLBA1", "CNRIP1", "COLEC12", "COMT", "CRIM1", "CXCL14", "CXCR4", "DST", "DYNLT1", "EMC1", "ENO2", "FAM124A", "FAM135A", "FAM9B", "FGD2", "FILIP1L", "GALNT11", "GATM", "GDA", "GJA1", "GLO1", "GNB4", "HAUS2", "HDDC3", "HLA-DQA1", "HMGN3", "KCNJ10", "LAMA3", "LCORL", "LYPLAL1", "MAF", "MALAT1", "MARCKSL1", "MARCO", "MSR1", "NAT8L", "NRCAM", "OCEL1", "OGFRL1", "P2RY13", "PIANP", "PIK3AP1", "PLAAT3", "PLBD1", "PLXDC2", "PPP2R5C", "PTGER3", "RAB10", "RAPSN", "RASAL2", "RCBTB2", "RCN1", "RFX3", "RPL14", "SFI1", "SLC35A1", "SLC7A7", "SLCO2B1", "SRD5A3", "TGFBI", "TIFAB", "TM7SF3", "TOR3A", "TTC3", "TUBB2B", "TXNIP", "ZNF727")
M1 <- intersect(all.genes, M1)

seurat.object <- AddModuleScore(object = seurat.object, features = M1, ctrl = 5, name = 'M1.score')
p1 <- FeaturePlot(seurat.object, features = 'M1.score1')
ggsave("FeaturePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'M1.score1')
ggsave("RidgePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "M1.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "M1.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "M1.score1", images=image.list, ncol=4)
p5
ggsave("M1.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

M2<- c("ADAM8", "AK4", "AKR1E2", "ALDOC", "ANG", "ANKRD37", "ANOS1", "ARG1", "ARSK", "ATG4D", "ATP6V0D2", "BNIP3", "BST1", "C1QB", "C1QC", "C1QTNF12", "C5orf34", "CBLB", "CD24", "CD300C", "CD300LD", "CD5L", "CHIA", "COL20A1", "CRIPT", "CTSK", "CTSL", "CYTIP", "EFHD2", "EGLN3", "ERO1A", "F10", "F7", "FAM177A1", "FAM199X", "FAM241A", "FBLIM1", "FLRT2", "GASK1B", "GDF15", "GPRC5B", "HILPDA", "HLA-A", "HLA-B", "HS6ST1", "HSPA1A", "HSPA1B", "IGF2R", "IGHM", "IL18BP", "ITGB3", "LBP", "LIN7C", "LRRC27", "MCOLN3", "MFGE8", "MGST2", "MYO1F", "PADI4", "PDXDC1", "PINK1", "PKDCC", "PRDX2", "PROCR", "PTGER2", "QRICH1", "RGCC", "RIMBP2", "RPGRIP1", "RRAS2", "SCD", "SH2B2", "SLC2A1", "SNHG6", "SOAT1", "THBS1", "TJP2", "TLR1", "TMEM267", "TULP4", "UCHL1", "VEGFA", "XDH")
M2 <- intersect(all.genes, M2)
seurat.object <- AddModuleScore(object = seurat.object, features = M2, ctrl = 5, name = 'M2.score')
p1 <- FeaturePlot(seurat.object, features = 'M2.score1')
ggsave("FeaturePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'M2.score1')
ggsave("RidgePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "M2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "M2.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "M2.score1", images=image.list, ncol=4)
p5
ggsave("M2.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

## http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_PI3K_AKT_MTOR_SIGNALING.html
pi3k <- c("ACACA", "ACTR2", "ACTR3", "ADCY2", "AKT1", "AKT1S1", "AP2M1", "ARF1", "ARHGDIA", "ARPC3", "ATF1", "CAB39", "CAB39L", "CALR", "CAMK4", "CDK1", "CDK2", "CDK4", "CDKN1A", "CDKN1B", "CFL1", "CLTC", "CSNK2B", "CXCR4", "DAPP1", "DDIT3", "DUSP3", "E2F1", "ECSIT", "EGFR", "EIF4E", "FASLG", "FGF17", "FGF22", "FGF6", "GNA14", "GNGT1", "GRB2", "GRK2", "GSK3B", "HRAS", "HSP90B1", "IL2RG", "IL4", "IRAK4", "ITPR2", "LCK", "MAP2K3", "MAP2K6", "MAP3K7", "MAPK1", "MAPK10", "MAPK8", "MAPK9", "MAPKAP1", "MKNK1", "MKNK2", "MYD88", "NCK1", "NFKBIB", "NGF", "NOD1", "PAK4", "PDK1", "PFN1", "PIK3R3", "PIKFYVE", "PIN1", "PITX2", "PLA2G12A", "PLCB1", "PLCG1", "PPP1CA", "PPP2R1B", "PRKAA2", "PRKAG1", "PRKAR2A", "PRKCB", "PTEN", "PTPN11", "RAC1", "RAF1", "RALB", "RIPK1", "RIT1", "RPS6KA1", "RPS6KA3", "RPTOR", "SFN", "SLA", "SLC2A1", "SMAD2", "SQSTM1", "STAT2", "TBK1", "THEM4", "TIAM1", "TNFRSF1A", "TRAF2", "TRIB3", "TSC2", "UBE2D3", "UBE2N", "VAV3", "YWHAB"))
seurat.object <- AddModuleScore(object = seurat.object, features = pi3k, ctrl = 5, name = 'pi3k.score')
p1 <- FeaturePlot(seurat.object, features = 'pi3k.score1')
ggsave("FeaturePlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'pi3k.score1')
ggsave("RidgePlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "pi3k.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "pi3k.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

inflammation <- c("ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "ADM", "ADORA2B", "ADRM1", "AHR", "APLNR", "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1", "BEST1", "BST2", "BTG2", "C3AR1", "C5AR1", "CALCRL", "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7", "CCRL2", "CD14", "CD40", "CD48", "CD55", "CD69", "CD70", "CD82", "CDKN1A", "CHST2", "CLEC5A", "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1", "CXCL10", "CXCL11", "CXCL6", "CXCL8", "CXCL9", "CXCR6", "CYBB", "DCBLD2", "EBI3", "EDN1", "EIF2AK2", "EMP3", "EREG", "F3", "FFAR2", "FPR1", "FZD5", "GABBR1", "GCH1", "GNA15", "GNAI3", "GP1BA", "GPC3", "GPR132", "GPR183", "HAS2", "HBEGF", "HIF1A", "HPN", "HRH1", "ICAM1", "ICAM4", "ICOSLG", "IFITM1", "IFNAR1", "IFNGR2", "IL10", "IL10RA", "IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5", "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2", "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV", "MEP1A", "MET", "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA", "NLRP3", "NMI", "NMUR1", "NOD2", "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2", "PCDH7", "PDE4B", "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2", "PTGER4", "PTGIR", "PTPRE", "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16", "RHOG", "RIPK2", "RNF144B", "ROS1", "RTP4", "SCARF1", "SCN1B", "SELE", "SELENOS", "SELL", "SEMA4D", "SERPINE1", "SGMS2", "SLAMF1", "SLC11A2", "SLC1A2", "SLC28A2", "SLC31A1", "SLC31A2", "SLC4A4", "SLC7A1", "SLC7A2", "SPHK1", "SRI", "STAB1", "TACR1", "TACR3", "TAPBP", "TIMP1", "TLR1", "TLR2", "TLR3", "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10", "TNFSF15", "TNFSF9", "TPBG", "VIP"))
inflammation <- intersect(inflammation, all.genes)
seurat.object <- AddModuleScore(object = seurat.object, features = inflammation, ctrl = 5, name = 'inflammation.score')
p1 <- FeaturePlot(seurat.object, features = 'inflammation.score1')
ggsave("FeaturePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'inflammation.score1')
ggsave("RidgePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "inflammation.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "inflammation.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p5 <- SpatialPlot(seurat.object, features = "inflammation.score1", images=image.list, ncol=4)
p5
ggsave("inflammation.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

ifn.alpha <- c("ADAR", "B2M", "BATF2", "BST2", "C1S", "CASP1", "CASP8", "CCRL2", "CD47", "CD74", "CMPK2", "CMTR1", "CNP", "CSF1", "CXCL10", "CXCL11", "DDX60", "DHX58", "EIF2AK2", "ELF1", "EPSTI1", "GBP2", "GBP4", "GMPR", "HELZ2", "HERC6", "HLA-C", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3", "IL15", "IL4R", "IL7", "IRF1", "IRF2", "IRF7", "IRF9", "ISG15", "ISG20", "LAMP3", "LAP3", "LGALS3BP", "LPAR6", "LY6E", "MOV10", "MVB12A", "MX1", "NCOA7", "NMI", "NUB1", "OAS1", "OASL", "OGFR", "PARP12", "PARP14", "PARP9", "PLSCR1", "PNPT1", "PROCR", "PSMA3", "PSMB8", "PSMB9", "PSME1", "PSME2", "RIPK2", "RNF31", "RSAD2", "RTP4", "SAMD9", "SAMD9L", "SELL", "SLC25A28", "SP110", "STAT2", "TAP1", "TDRD7", "TENT5A", "TMEM140", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TRIM5", "TXNIP", "UBA7", "UBE2L6", "USP18", "WARS1"))
seurat.object <- AddModuleScore(object = seurat.object, features = ifn.alpha, ctrl = 5, name = 'ifn.alpha.score')
p1 <- FeaturePlot(seurat.object, features = 'ifn.alpha.score1')
ggsave("FeaturePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'ifn.alpha.score1')
ggsave("RidgePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "ifn.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "ifn.alpha.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

ifn.beta <- c("ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM", "BST2", "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1", "CMPK2", "CMTR1", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX58", "DDX60", "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "GBP4", "GBP6", "GCH1", "GPR18", "GZMA", "HELZ2", "HERC6", "HIF1A", "HLA-A", "HLA-B", "HLA-DMA", "HLA-DQA1", "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IFNAR2", "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1", "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7", "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1", "METTL7B", "MT2A", "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12", "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP", "PNPT1", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1", "RIPK2", "RNF213", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP", "SERPING1", "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2", "SP110", "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4", "STAT1", "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP", "TDRD7", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8", "VCAM1", "WARS1", "XAF1", "XCL1", "ZBP1", "ZNFX1"))
seurat.object <- AddModuleScore(object = seurat.object, features = ifn.beta, ctrl = 5, name = 'ifn.beta.score')
p1 <- FeaturePlot(seurat.object, features = 'ifn.beta.score1')
ggsave("FeaturePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'ifn.beta.score1')
ggsave("RidgePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "ifn.beta.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "ifn.beta.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tgfb <- c("ACVR1", "APC", "ARID4B", "BCAR3", "BMP2", "BMPR1A", "BMPR2", "CDH1", "CDK9", "CDKN1C", "CTNNB1", "ENG", "FKBP1A", "FNTA", "FURIN", "HDAC1", "HIPK2", "ID1", "ID2", "ID3", "IFNGR2", "JUNB", "KLF10", "LEFTY2", "LTBP2", "MAP3K7", "NCOR2", "NOG", "PMEPA1", "PPM1A", "PPP1CA", "PPP1R15A", "RAB31", "RHOA", "SERPINE1", "SKI", "SKIL", "SLC20A1", "SMAD1", "SMAD3", "SMAD6", "SMAD7", "SMURF1", "SMURF2", "SPTBN1", "TGFB1", "TGFBR1", "TGIF1", "THBS1", "TJP1", "TRIM33", "UBE2D3", "WWTR1", "XIAP"))
seurat.object <- AddModuleScore(object = seurat.object, features = tgfb, ctrl = 5, name = 'tgfb.score')
p1 <- FeaturePlot(seurat.object, features = 'tgfb.score1')
ggsave("FeaturePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'tgfb.score1')
ggsave("RidgePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "tgfb.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "tgfb.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

wnt <- c("ADAM17", "AXIN1", "AXIN2", "CCND2", "CSNK1E", "CTNNB1", "CUL1", "DKK1", "DKK4", "DLL1", "DVL2", "FRAT1", "FZD1", "FZD8", "GNAI1", "HDAC11", "HDAC2", "HDAC5", "HEY1", "HEY2", "JAG1", "JAG2", "KAT2A", "LEF1", "MAML1", "MYC", "NCOR2", "NCSTN", "NKD1", "NOTCH1", "NOTCH4", "NUMB", "PPARD", "PSEN2", "PTCH1", "RBPJ", "SKP2", "TCF7", "TP53", "WNT1", "WNT5B", "WNT6"))
seurat.object <- AddModuleScore(object = seurat.object, features = wnt, ctrl = 5, name = 'wnt.score')
p1 <- FeaturePlot(seurat.object, features = 'wnt.score1')
ggsave("FeaturePlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'wnt.score1')
ggsave("RidgePlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "wnt.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "wnt.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

il2 <- c("ABCB1", "ADAM19", "AGER", "AHCY", "AHNAK", "AHR", "ALCAM", "AMACR", "ANXA4", "APLP1", "ARL4A", "BATF", "BATF3", "BCL2", "BCL2L1", "BHLHE40", "BMP2", "BMPR2", "CA2", "CAPG", "CAPN3", "CASP3", "CCND2", "CCND3", "CCNE1", "CCR4", "CD44", "CD48", "CD79B", "CD81", "CD83", "CD86", "CDC42SE2", "CDC6", "CDCP1", "CDKN1C", "CISH", "CKAP4", "COCH", "COL6A1", "CSF1", "CSF2", "CST7", "CTLA4", "CTSZ", "CXCL10", "CYFIP1", "DCPS", "DENND5A", "DHRS3", "DRC1", "ECM1", "EEF1AKMT1", "EMP1", "ENO3", "ENPP1", "EOMES", "ETFBKMT", "ETV4", "F2RL2", "FAH", "FAM126B", "FGL2", "FLT3LG", "FURIN", "GABARAPL1", "GADD45B", "GALM", "GATA1", "GBP4", "GLIPR2", "GPR65", "GPR83", "GPX4", "GSTO1", "GUCY1B1", "HIPK2", "HK2", "HOPX", "HUWE1", "ICOS", "IFITM3", "IFNGR1", "IGF1R", "IGF2R", "IKZF2", "IKZF4", "IL10", "IL10RA", "IL13", "IL18R1", "IL1R2", "IL1RL1", "IL2RA", "IL2RB", "IL3RA", "IL4R", "IRF4", "IRF6", "IRF8", "ITGA6", "ITGAE", "ITGAV", "ITIH5", "KLF6", "LCLAT1", "LIF", "LRIG1", "LRRC8C", "LTB", "MAFF", "MAP3K8", "MAP6", "MAPKAPK2", "MUC1", "MXD1", "MYC", "MYO1C", "MYO1E", "NCOA3", "NCS1", "NDRG1", "NFIL3", "NFKBIZ", "NOP2", "NRP1", "NT5E", "ODC1", "P2RX4", "P4HA1", "PDCD2L", "PENK", "PHLDA1", "PHTF2", "PIM1", "PLAGL1", "PLEC", "PLIN2", "PLPP1", "PLSCR1", "PNP", "POU2F1", "PRAF2", "PRKCH", "PRNP", "PTCH1", "PTGER2", "PTH1R", "PTRH2", "PUS1", "RABGAP1L", "RGS16", "RHOB", "RHOH", "RNH1", "RORA", "RRAGD", "S100A1", "SCN9A", "SELL", "SELP", "SERPINB6", "SERPINC1", "SH3BGRL2", "SHE", "SLC1A5", "SLC29A2", "SLC2A3", "SLC39A8", "SMPDL3A", "SNX14", "SNX9", "SOCS1", "SOCS2", "SPP1", "SPRED2", "SPRY4", "ST3GAL4", "SWAP70", "SYNGR2", "SYT11", "TGM2", "TIAM1", "TLR7", "TNFRSF18", "TNFRSF1B", "TNFRSF21", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFSF10", "TNFSF11", "TRAF1", "TTC39B", "TWSG1", "UCK2", "UMPS", "WLS", "XBP1"))
seurat.object <- AddModuleScore(object = seurat.object, features = il2, ctrl = 5, name = 'il2.score')
p1 <- FeaturePlot(seurat.object, features = 'il2.score1')
ggsave("FeaturePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'il2.score1')
ggsave("RidgePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "il2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "il2.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

il6 <- c("A2M", "ACVR1B", "ACVRL1", "BAK1", "CBL", "CCL7", "CCR1", "CD14", "CD36", "CD38", "CD44", "CD9", "CNTFR", "CRLF2", "CSF1", "CSF2", "CSF2RA", "CSF2RB", "CSF3R", "CXCL1", "CXCL10", "CXCL11", "CXCL13", "CXCL3", "CXCL9", "DNTT", "EBI3", "FAS", "GRB2", "HAX1", "HMOX1", "IFNAR1", "IFNGR1", "IFNGR2", "IL10RB", "IL12RB1", "IL13RA1", "IL15RA", "IL17RA", "IL17RB", "IL18R1", "IL1B", "IL1R1", "IL1R2", "IL2RA", "IL2RG", "IL3RA", "IL4R", "IL6", "IL6ST", "IL7", "IL9R", "INHBE", "IRF1", "IRF9", "ITGA4", "ITGB3", "JUN", "LEPR", "LTB", "LTBR", "MAP3K8", "MYD88", "OSMR", "PDGFC", "PF4", "PIK3R5", "PIM1", "PLA2G2A", "PTPN1", "PTPN11", "PTPN2", "REG1A", "SOCS1", "SOCS3", "STAM2", "STAT1", "STAT2", "STAT3", "TGFB1", "TLR2", "TNF", "TNFRSF12A", "TNFRSF1A", "TNFRSF1B", "TNFRSF21", "TYK2"))
seurat.object <- AddModuleScore(object = seurat.object, features = il6, ctrl = 5, name = 'il6.score')
p1 <- FeaturePlot(seurat.object, features = 'il6.score1')
ggsave("FeaturePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'il6.score1')
ggsave("RidgePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "il6.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "il6.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

apoptosis <- c("ADD1", "AIFM3", "ANKH", "ANXA1", "APP", "ATF3", "AVPR1A", "BAX", "BCAP31", "BCL10", "BCL2L1", "BCL2L10", "BCL2L11", "BCL2L2", "BGN", "BID", "BIK", "BIRC3", "BMF", "BMP2", "BNIP3L", "BRCA1", "BTG2", "BTG3", "CASP1", "CASP2", "CASP3", "CASP4", "CASP6", "CASP7", "CASP8", "CASP9", "CAV1", "CCNA1", "CCND1", "CCND2", "CD14", "CD2", "CD38", "CD44", "CD69", "CDC25B", "CDK2", "CDKN1A", "CDKN1B", "CFLAR", "CLU", "CREBBP", "CTH", "CTNNB1", "CYLD", "DAP", "DAP3", "DCN", "DDIT3", "DFFA", "DIABLO", "DNAJA1", "DNAJC3", "DNM1L", "DPYD", "EBP", "EGR3", "EMP1", "ENO2", "ERBB2", "ERBB3", "EREG", "ETF1", "F2", "F2R", "FAS", "FASLG", "FDXR", "FEZ1", "GADD45A", "GADD45B", "GCH1", "GNA15", "GPX1", "GPX3", "GPX4", "GSN", "GSR", "GSTM1", "GUCY2D", "H1-0", "HGF", "HMGB2", "HMOX1", "HSPB1", "IER3", "IFITM3", "IFNB1", "IFNGR1", "IGF2R", "IGFBP6", "IL18", "IL1A", "IL1B", "IL6", "IRF1", "ISG20", "JUN", "KRT18", "LEF1", "LGALS3", "LMNA", "LUM", "MADD", "MCL1", "MGMT", "MMP2", "NEDD9", "NEFH", "PAK1", "PDCD4", "PDGFRB", "PEA15", "PLAT", "PLCB2", "PLPPR4", "PMAIP1", "PPP2R5B", "PPP3R1", "PPT1", "PRF1", "PSEN1", "PSEN2", "PTK2", "RARA", "RELA", "RETSAT", "RHOB", "RHOT2", "RNASEL", "ROCK1", "SAT1", "SATB1", "SC5D", "SLC20A1", "SMAD7", "SOD1", "SOD2", "SPTAN1", "SQSTM1", "TAP1", "TGFB2", "TGFBR3", "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF12A", "TNFSF10", "TOP2A", "TSPO", "TXNIP", "VDAC2", "WEE1", "XIAP"))
seurat.object <- AddModuleScore(object = seurat.object, features = apoptosis, ctrl = 5, name = 'apoptosis.score')
p1 <- FeaturePlot(seurat.object, features = 'apoptosis.score1')
ggsave("FeaturePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'apoptosis.score1')
ggsave("RidgePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "apoptosis.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "apoptosis.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tnf.alpha <- c("ABCA1", "ACKR3", "AREG", "ATF3", "ATP2B1", "B4GALT1", "B4GALT5", "BCL2A1", "BCL3", "BCL6", "BHLHE40", "BIRC2", "BIRC3", "BMP2", "BTG1", "BTG2", "BTG3", "CCL2", "CCL20", "CCL4", "CCL5", "CCN1", "CCND1", "CCNL1", "CCRL2", "CD44", "CD69", "CD80", "CD83", "CDKN1A", "CEBPB", "CEBPD", "CFLAR", "CLCF1", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL11", "CXCL2", "CXCL3", "CXCL6", "DDX58", "DENND5A", "DNAJB4", "DRAM1", "DUSP1", "DUSP2", "DUSP4", "DUSP5", "EDN1", "EFNA1", "EGR1", "EGR2", "EGR3", "EHD1", "EIF1", "ETS2", "F2RL1", "F3", "FJX1", "FOS", "FOSB", "FOSL1", "FOSL2", "FUT4", "G0S2", "GADD45A", "GADD45B", "GCH1", "GEM", "GFPT2", "GPR183", "HBEGF", "HES1", "ICAM1", "ICOSLG", "ID2", "IER2", "IER3", "IER5", "IFIH1", "IFIT2", "IFNGR2", "IL12B", "IL15RA", "IL18", "IL1A", "IL1B", "IL23A", "IL6", "IL6ST", "IL7R", "INHBA", "IRF1", "IRS2", "JAG1", "JUN", "JUNB", "KDM6B", "KLF10", "KLF2", "KLF4", "KLF6", "KLF9", "KYNU", "LAMB3", "LDLR", "LIF", "LITAF", "MAFF", "MAP2K3", "MAP3K8", "MARCKS", "MCL1", "MSC", "MXD1", "MYC", "NAMPT", "NFAT5", "NFE2L2", "NFIL3", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NINJ1", "NR4A1", "NR4A2", "NR4A3", "OLR1", "PANX1", "PDE4B", "PDLIM5", "PER1", "PFKFB3", "PHLDA1", "PHLDA2", "PLAU", "PLAUR", "PLEK", "PLK2", "PLPP3", "PMEPA1", "PNRC1", "PPP1R15A", "PTGER4", "PTGS2", "PTPRE", "PTX3", "RCAN1", "REL", "RELA", "RELB", "RHOB", "RIPK2", "RNF19B", "SAT1", "SDC4", "SERPINB2", "SERPINB8", "SERPINE1", "SGK1", "SIK1", "SLC16A6", "SLC2A3", "SLC2A6", "SMAD3", "SNN", "SOCS3", "SOD2", "SPHK1", "SPSB1", "SQSTM1", "STAT5A", "TANK", "TAP1", "TGIF1", "TIPARP", "TLR2", "TNC", "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFRSF9", "TNFSF9", "TNIP1", "TNIP2", "TRAF1", "TRIB1", "TRIP10", "TSC22D1", "TUBB2A", "VEGFA", "YRDC", "ZBTB10", "ZC3H12A", "ZFP36"))
seurat.object <- AddModuleScore(object = seurat.object, features = tnf.alpha, ctrl = 5, name = 'tnf.alpha.score')
p1 <- FeaturePlot(seurat.object, features = 'tnf.alpha.score1')
ggsave("FeaturePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'tnf.alpha.score1')
ggsave("RidgePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "tnf.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "tnf.alpha.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

setwd("..")



############################################################# ############################################################# #############################################################
############################################################# VII. Subset analysis #############################################################
############################################################# ############################################################# #############################################################
#
############################################################# ############################################################# #############################################################
############################################################# Subset Analysis: ZIKV-Veh-1 v. ZIKV-Enox-1 #############################################################
############################################################# ############################################################# #############################################################
dir.create("zikv.enox")
setwd("zikv.enox")
Idents(seurat.object)
Idents(seurat.object) <- "Treatment"
seurat.object2 <- subset(x = seurat.object, idents = c("ZIKV-Veh", "ZIKV-Enox"))
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
de_markers <- FindAllMarkers(seurat.object2, features = seurat.object2.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, raster = FALSE)
ggsave("integrated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("integrated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

Idents(seurat.object2) <- "Type.Profile"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("integrated_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type.Profile")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("integrated_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type.Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("integrated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type.Profile")
dev.off()

##### Perform differential expression between orig.ident
Idents(seurat.object2) <- "Treatment"
de_markers <- FindAllMarkers(seurat.object2, features = seurat.object2.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "Treatment_DEGs_pos-log2FC.txt", sep="\t")

ggsave("Treatment_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("Treatment_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
View(top2)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("orig.ident_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Treatment")
dev.off()

inflammation <- c("ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "ADM", "ADORA2B", "ADRM1", "AHR", "APLNR", "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1", "BEST1", "BST2", "BTG2", "C3AR1", "C5AR1", "CALCRL", "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7", "CCRL2", "CD14", "CD40", "CD48", "CD55", "CD69", "CD70", "CD82", "CDKN1A", "CHST2", "CLEC5A", "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1", "CXCL10", "CXCL11", "CXCL6", "CXCL8", "CXCL9", "CXCR6", "CYBB", "DCBLD2", "EBI3", "EDN1", "EIF2AK2", "EMP3", "EREG", "F3", "FFAR2", "FPR1", "FZD5", "GABBR1", "GCH1", "GNA15", "GNAI3", "GP1BA", "GPC3", "GPR132", "GPR183", "HAS2", "HBEGF", "HIF1A", "HPN", "HRH1", "ICAM1", "ICAM4", "ICOSLG", "IFITM1", "IFNAR1", "IFNGR2", "IL10", "IL10RA", "IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5", "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2", "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV", "MEP1A", "MET", "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA", "NLRP3", "NMI", "NMUR1", "NOD2", "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2", "PCDH7", "PDE4B", "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2", "PTGER4", "PTGIR", "PTPRE", "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16", "RHOG", "RIPK2", "RNF144B", "ROS1", "RTP4", "SCARF1", "SCN1B", "SELE", "SELENOS", "SELL", "SEMA4D", "SERPINE1", "SGMS2", "SLAMF1", "SLC11A2", "SLC1A2", "SLC28A2", "SLC31A1", "SLC31A2", "SLC4A4", "SLC7A1", "SLC7A2", "SPHK1", "SRI", "STAB1", "TACR1", "TACR3", "TAPBP", "TIMP1", "TLR1", "TLR2", "TLR3", "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10", "TNFSF15", "TNFSF9", "TPBG", "VIP"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = inflammation, ctrl = 5, name = 'inflammation.score')
p1 <- FeaturePlot(seurat.object2, features = 'inflammation.score1')
ggsave("FeaturePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'inflammation.score1')
ggsave("RidgePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "inflammation.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object2, feature1 = "inflammation.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")



rm(seurat.object2)

setwd("..")
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# Subset Analysis: ZIKV-Veh-1 v. ZIKV-Enox-1 #############################################################
############################################################# ############################################################# #############################################################
dir.create("Mock.ZIKV")
setwd("Mock.ZIKV")
Idents(seurat.object) <- "Treatment"
seurat.object2 <- subset(x = seurat.object, idents = c("ZIKV-Veh", "Mock-Veh", "Mock-Enox"))
ncol(seurat.object2)
	# 7940
Idents(seurat.object2) <- "Virus"
levels(seurat.object2)
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
p3 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Virus")
p3
p4 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Type.Profile")
p4
umap.combined <- CombinePlots(plots = list(p1, p2, p3))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("integrated_UMAP.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object2 <- NormalizeData(seurat.object2, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object2 <- ScaleData(seurat.object2, features = seurat.object2.var, assay="Spatial")

Idents(seurat.object2) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object2, features = seurat.object2.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, raster = FALSE)
ggsave("seurat_clusters_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("seurat_clusters_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

Idents(seurat.object2) <- "Type.Profile"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("Type.Profile_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type.Profile")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("Type.Profile_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type.Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("Type.Profile_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type.Profile")
dev.off()

##### Perform differential expression between orig.ident
Idents(seurat.object2) <- "Treatment"
de_markers <- FindAllMarkers(seurat.object2, features = seurat.object2.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "Treatment_DEGs_pos-log2FC.txt", sep="\t")

ggsave("Treatment_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("Treatment_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
View(top2)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("Treatment_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Treatment")
dev.off()

##### Perform differential expression between orig.ident
Idents(seurat.object2) <- "Virus"
de_markers <- FindAllMarkers(seurat.object2, features = seurat.object2.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "Virus_DEGs_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
View(top5)
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("Virus_top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("Virus_top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
View(top2)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("Virus_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Treatment")
dev.off()

rm(seurat.object2)

setwd("..")
############################################################# ############################################################# #############################################################
FetalMesenchyme/Immune/Endothelial

############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# Subset Analysis: ZIKV-Veh-1 v. ZIKV-Enox-1 #############################################################
############################################################# ############################################################# #############################################################
dir.create("immune")
setwd("immune")
Idents(seurat.object) <- "Subtype.Profile"
seurat.object2 <- subset(x = seurat.object, idents = c("FetalMesenchyme/Immune/Endothelial"))
ncol(seurat.object2)
	# 667
Idents(seurat.object2) <- "Treatment"
levels(seurat.object2)
seurat.object2 <- SCTransform(seurat.object2, assay = "Spatial", verbose = T)
DefaultAssay(seurat.object2) <- "SCT"
seurat.object2.var <- seurat.object2@assays[["SCT"]]@var.features

seurat.object2 <- ScaleData(seurat.object2, features = seurat.object2.var, assay = "SCT", vars.to.regress = c("CC.Difference", "nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"))

seurat.object2<- RunPCA(seurat.object2, features = seurat.object2.var, assay = "SCT", slot = "scale.data")
seurat.object2 <- FindNeighbors(seurat.object2, features = "seurat.object2.var", dims = 1:30)
seurat.object2 <- FindClusters(seurat.object2, resolution = 0.6)
seurat.object2 <- RunUMAP(seurat.object2, dims=1:30)
	## clusters 0 thru 12

Idents(seurat.object2) <- "Treatment"
t1<-table(Idents(seurat.object2))
write.table(t1, "Treatment.counts.txt", sep="\t")
rm(t1)

p1 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Treatment", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Treatment')
p1
p2 <- DimPlot(seurat.object2, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p2
p3 <- DimPlot(seurat.object2, reduction = "umap", group.by = "orig.ident", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Sample')
p3
p4 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Batch", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Batch')
p4
p5 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Type.Profile",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type.Profile')
p5
p6 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Subtype.Profile",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Subtype.Profile')
p6
p7 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Virus", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Virus')
p7
umap.combined <- p1 + p2 + p6 + p7
ggsave("Visium_UMAP-annotated_1.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("Visium_UMAP-annotated_XL_1.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")
ggsave("Visium_UMAP-annotated_M_1.pdf", plot = umap.combined, device = "pdf", width = 10, height = 6, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object2 <- NormalizeData(seurat.object2, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object2 <- ScaleData(seurat.object2, features = seurat.object2.var, assay="Spatial")

Idents(seurat.object2) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object2, features = seurat.object2.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, raster = FALSE)
ggsave("seurat_clusters_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("seurat_clusters_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

############################################################# ############################################################# #############################################################


dir.create("GSEA_Subtype.Profile")
setwd("GSEA_Subtype.Profile")
Idents(seurat.object2) <- "Subtype.Profile"

## use RNA for all DE analysis and plots
DefaultAssay(seurat.object2) <- "Spatial"
# seurat.object2 <- ScaleData(seurat.object2, features = all.genes, assay = "Spatial")
seurat.object2 <- ScaleData(seurat.object2, features = all.genes.visium, assay="Spatial")

hypoxia.list <- c("ACKR3", "ADM", "ADORA2B", "AK4", "AKAP12", "ALDOA", "ALDOB", "ALDOC", "AMPD3", "ANGPTL4", "ANKZF1", "ANXA2", "ATF3", "ATP7A", "B3GALT6", "B4GALNT2", "BCAN", "BCL2", "BGN", "BHLHE40", "BNIP3L", "BRS3", "BTG1", "CA12", "CASP6", "CAV1", "CAVIN1", "CAVIN3", "CCN1", "CCN2", "CCN5", "CCNG2", "CDKN1A", "CDKN1B", "CDKN1C", "CHST2", "CHST3", "CITED2", "COL5A1", "CP", "CSRP2", "CXCR4", "DCN", "DDIT3", "DDIT4", "DPYSL4", "DTNA", "DUSP1", "EDN2", "EFNA1", "EFNA3", "EGFR", "ENO1", "ENO2", "ENO3", "ERO1A", "ERRFI1", "ETS1", "EXT1", "F3", "FAM162A", "FBP1", "FOS", "FOSL2", "FOXO3", "GAA", "GALK1", "GAPDH", "GAPDHS", "GBE1", "GCK", "GCNT2", "GLRX", "GPC1", "GPC3", "GPC4", "GPI", "GRHPR", "GYS1", "HAS1", "HDLBP", "HEXA", "HK1", "HK2", "HMOX1", "HOXB9", "HS3ST1", "HSPA5", "IDS", "IER3", "IGFBP1", "IGFBP3", "IL6", "ILVBL", "INHA", "IRS2", "ISG20", "JMJD6", "JUN", "KDELR3", "KDM3A", "KIF5A", "KLF6", "KLF7", "KLHL24", "LALBA", "LARGE1", "LDHA", "LDHC", "LOX", "LXN", "MAFF", "MAP3K1", "MIF", "MT1E", "MT2A", "MXI1", "MYH9", "NAGK", "NCAN", "NDRG1", "NDST1", "NDST2", "NEDD4L", "NFIL3", "NOCT", "NR3C1", "P4HA1", "P4HA2", "PAM", "PCK1", "PDGFB", "PDK1", "PDK3", "PFKFB3", "PFKL", "PFKP", "PGAM2", "PGF", "PGK1", "PGM1", "PGM2", "PHKG1", "PIM1", "PKLR", "PKP1", "PLAC8", "PLAUR", "PLIN2", "PNRC1", "PPARGC1A", "PPFIA4", "PPP1R15A", "PPP1R3C", "PRDX5", "PRKCA", "PYGM", "RBPJ", "RORA", "RRAGD", "S100A4", "SAP30", "SCARB1", "SDC2", "SDC3", "SDC4", "SELENBP1", "SERPINE1", "SIAH2", "SLC25A1", "SLC2A1", "SLC2A3", "SLC2A5", "SLC37A4", "SLC6A6", "SRPX", "STBD1", "STC1", "STC2", "SULT2B1", "TES", "TGFB3", "TGFBI", "TGM2", "TIPARP", "TKTL1", "TMEM45A", "TNFAIP3", "TPBG", "TPD52", "TPI1", "TPST2", "UGP2", "VEGFA", "VHL", "VLDLR", "WSB1", "XPNPEP1", "ZFP36", "ZNF292")
hypoxia.list <- intersect(all.genes, hypoxia.list)

seurat.object2 <- AddModuleScore(object = seurat.object2, features = hypoxia.list, ctrl = 5, name = 'hypoxia.score')
# p1 <- FeaturePlot(seurat.object2, features = 'hypoxia.score1')
# ggsave("FeaturePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'hypoxia.score1')
ggsave("RidgePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "hypoxia.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "hypoxia.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "hypoxia.score1", images=image.list, ncol=4)
p5
ggsave("hypoxia.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


complement.list <- c("ACTN2", "ADAM9", "ADRA2B", "AKAP10", "ANG", "ANXA5", "APOA4", "APOBEC3F", "APOBEC3G", "APOC1", "ATOX1", "BRPF3", "C1QA", "C1QC", "C1R", "C1S", "C2", "C3", "C4BPB", "C9", "CA2", "CALM1", "CALM3", "CASP1", "CASP10", "CASP3", "CASP4", "CASP5", "CASP7", "CASP9", "CBLB", "CCL5", "CD36", "CD40LG", "CD46", "CD55", "CD59", "CDA", "CDH13", "CDK5R1", "CEBPB", "CFB", "CFH", "CLU", "COL4A2", "CP", "CPM", "CPQ", "CR1", "CR2", "CSRP1", "CTSB", "CTSC", "CTSD", "CTSH", "CTSL", "CTSO", "CTSS", "CTSV", "CXCL1", "DGKG", "DGKH", "DOCK10", "DOCK4", "DOCK9", "DPP4", "DUSP5", "DUSP6", "DYRK2", "EHD1", "ERAP2", "F10", "F2", "F3", "F5", "F7", "F8", "FCER1G", "FCN1", "FDX1", "FN1", "FYN", "GATA3", "GCA", "GMFB", "GNAI2", "GNAI3", "GNB2", "GNB4", "GNG2", "GNGT2", "GP1BA", "GP9", "GPD2", "GRB2", "GZMA", "GZMB", "GZMK", "HNF4A", "HPCAL4", "HSPA1A", "HSPA5", "IL6", "IRF1", "IRF2", "IRF7", "ITGAM", "ITIH1", "JAK2", "KCNIP2", "KCNIP3", "KIF2A", "KLK1", "KLKB1", "KYNU", "L3MBTL4", "LAMP2", "LAP3", "LCK", "LCP2", "LGALS3", "LGMN", "LIPA", "LRP1", "LTA4H", "LTF", "LYN", "MAFF", "ME1", "MMP12", "MMP13", "MMP14", "MMP15", "MMP8", "MSRB1", "MT3", "NOTCH4", "OLR1", "PCLO", "PCSK9", "PDGFB", "PDP1", "PFN1", "PHEX", "PIK3CA", "PIK3CG", "PIK3R5", "PIM1", "PLA2G4A", "PLA2G7", "PLAT", "PLAUR", "PLEK", "PLG", "PLSCR1", "PPP2CB", "PPP4C", "PRCP", "PRDM4", "PREP", "PRKCD", "PRSS3", "PRSS36", "PSEN1", "PSMB9", "RABIF", "RAF1", "RASGRP1", "RBSN", "RCE1", "RHOG", "RNF4", "S100A12", "S100A13", "S100A9", "SCG3", "SERPINA1", "SERPINB2", "SERPINC1", "SERPINE1", "SERPING1", "SH2B3", "SIRT6", "SPOCK2", "SRC", "STX4", "TFPI2", "TIMP1", "TIMP2", "TMPRSS6", "TNFAIP3", "USP14", "USP15", "USP16", "USP8", "VCPIP1", "WAS", "XPNPEP1", "ZEB1", "ZFPM2")
complement.list <- intersect(all.genes, complement.list)
seurat.object2 <- AddModuleScore(object = seurat.object2, features = complement.list, ctrl = 5, name = 'complement.score')
# p1 <- FeaturePlot(seurat.object2, features = 'complement.score1')
# ggsave("FeaturePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'complement.score1')
ggsave("RidgePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "complement.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "complement.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "complement.score1", images=image.list, ncol=4)
p5
ggsave("complement.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


# http://www.gsea-msigdb.org/gsea/msigdb/cards/COATES_MACROPHAGE_M1_VS_M2_UP.html
M1 <- c("ABCA9", "ABHD1", "ACSS1", "ACYP2", "ADA", "AIF1", "ALAD", "ARFGEF3", "ARSB", "ATP6V0E2", "AXL", "BCKDHB", "BLOC1S6", "CADM1", "CAP1", "CAPN5", "CBX6", "CD59", "CFH", "CLBA1", "CNRIP1", "COLEC12", "COMT", "CRIM1", "CXCL14", "CXCR4", "DST", "DYNLT1", "EMC1", "ENO2", "FAM124A", "FAM135A", "FAM9B", "FGD2", "FILIP1L", "GALNT11", "GATM", "GDA", "GJA1", "GLO1", "GNB4", "HAUS2", "HDDC3", "HLA-DQA1", "HMGN3", "KCNJ10", "LAMA3", "LCORL", "LYPLAL1", "MAF", "MALAT1", "MARCKSL1", "MARCO", "MSR1", "NAT8L", "NRCAM", "OCEL1", "OGFRL1", "P2RY13", "PIANP", "PIK3AP1", "PLAAT3", "PLBD1", "PLXDC2", "PPP2R5C", "PTGER3", "RAB10", "RAPSN", "RASAL2", "RCBTB2", "RCN1", "RFX3", "RPL14", "SFI1", "SLC35A1", "SLC7A7", "SLCO2B1", "SRD5A3", "TGFBI", "TIFAB", "TM7SF3", "TOR3A", "TTC3", "TUBB2B", "TXNIP", "ZNF727")
M1 <- intersect(all.genes, M1)

seurat.object2 <- AddModuleScore(object = seurat.object2, features = M1, ctrl = 5, name = 'M1.score')
# p1 <- FeaturePlot(seurat.object2, features = 'M1.score1')
# ggsave("FeaturePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'M1.score1')
ggsave("RidgePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "M1.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "M1.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "M1.score1", images=image.list, ncol=4)
p5
ggsave("M1.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

M2<- c("ADAM8", "AK4", "AKR1E2", "ALDOC", "ANG", "ANKRD37", "ANOS1", "ARG1", "ARSK", "ATG4D", "ATP6V0D2", "BNIP3", "BST1", "C1QB", "C1QC", "C1QTNF12", "C5orf34", "CBLB", "CD24", "CD300C", "CD300LD", "CD5L", "CHIA", "COL20A1", "CRIPT", "CTSK", "CTSL", "CYTIP", "EFHD2", "EGLN3", "ERO1A", "F10", "F7", "FAM177A1", "FAM199X", "FAM241A", "FBLIM1", "FLRT2", "GASK1B", "GDF15", "GPRC5B", "HILPDA", "HLA-A", "HLA-B", "HS6ST1", "HSPA1A", "HSPA1B", "IGF2R", "IGHM", "IL18BP", "ITGB3", "LBP", "LIN7C", "LRRC27", "MCOLN3", "MFGE8", "MGST2", "MYO1F", "PADI4", "PDXDC1", "PINK1", "PKDCC", "PRDX2", "PROCR", "PTGER2", "QRICH1", "RGCC", "RIMBP2", "RPGRIP1", "RRAS2", "SCD", "SH2B2", "SLC2A1", "SNHG6", "SOAT1", "THBS1", "TJP2", "TLR1", "TMEM267", "TULP4", "UCHL1", "VEGFA", "XDH")
M2 <- intersect(all.genes, M2)
seurat.object2 <- AddModuleScore(object = seurat.object2, features = M2, ctrl = 5, name = 'M2.score')
# p1 <- FeaturePlot(seurat.object2, features = 'M2.score1')
# ggsave("FeaturePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'M2.score1')
ggsave("RidgePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "M2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "M2.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "M2.score1", images=image.list, ncol=4)
p5
ggsave("M2.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


## http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_PI3K_AKT_MTOR_SIGNALING.html
pi3k <- c("ACACA", "ACTR2", "ACTR3", "ADCY2", "AKT1", "AKT1S1", "AP2M1", "ARF1", "ARHGDIA", "ARPC3", "ATF1", "CAB39", "CAB39L", "CALR", "CAMK4", "CDK1", "CDK2", "CDK4", "CDKN1A", "CDKN1B", "CFL1", "CLTC", "CSNK2B", "CXCR4", "DAPP1", "DDIT3", "DUSP3", "E2F1", "ECSIT", "EGFR", "EIF4E", "FASLG", "FGF17", "FGF22", "FGF6", "GNA14", "GNGT1", "GRB2", "GRK2", "GSK3B", "HRAS", "HSP90B1", "IL2RG", "IL4", "IRAK4", "ITPR2", "LCK", "MAP2K3", "MAP2K6", "MAP3K7", "MAPK1", "MAPK10", "MAPK8", "MAPK9", "MAPKAP1", "MKNK1", "MKNK2", "MYD88", "NCK1", "NFKBIB", "NGF", "NOD1", "PAK4", "PDK1", "PFN1", "PIK3R3", "PIKFYVE", "PIN1", "PITX2", "PLA2G12A", "PLCB1", "PLCG1", "PPP1CA", "PPP2R1B", "PRKAA2", "PRKAG1", "PRKAR2A", "PRKCB", "PTEN", "PTPN11", "RAC1", "RAF1", "RALB", "RIPK1", "RIT1", "RPS6KA1", "RPS6KA3", "RPTOR", "SFN", "SLA", "SLC2A1", "SMAD2", "SQSTM1", "STAT2", "TBK1", "THEM4", "TIAM1", "TNFRSF1A", "TRAF2", "TRIB3", "TSC2", "UBE2D3", "UBE2N", "VAV3", "YWHAB"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = pi3k, ctrl = 5, name = 'pi3k.score')
# p1 <- FeaturePlot(seurat.object2, features = 'pi3k.score1')
# ggsave("FeaturePlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'pi3k.score1')
ggsave("RidgePlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "pi3k.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "pi3k.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

inflammation <- c("ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "ADM", "ADORA2B", "ADRM1", "AHR", "APLNR", "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1", "BEST1", "BST2", "BTG2", "C3AR1", "C5AR1", "CALCRL", "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7", "CCRL2", "CD14", "CD40", "CD48", "CD55", "CD69", "CD70", "CD82", "CDKN1A", "CHST2", "CLEC5A", "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1", "CXCL10", "CXCL11", "CXCL6", "CXCL8", "CXCL9", "CXCR6", "CYBB", "DCBLD2", "EBI3", "EDN1", "EIF2AK2", "EMP3", "EREG", "F3", "FFAR2", "FPR1", "FZD5", "GABBR1", "GCH1", "GNA15", "GNAI3", "GP1BA", "GPC3", "GPR132", "GPR183", "HAS2", "HBEGF", "HIF1A", "HPN", "HRH1", "ICAM1", "ICAM4", "ICOSLG", "IFITM1", "IFNAR1", "IFNGR2", "IL10", "IL10RA", "IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5", "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2", "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV", "MEP1A", "MET", "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA", "NLRP3", "NMI", "NMUR1", "NOD2", "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2", "PCDH7", "PDE4B", "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2", "PTGER4", "PTGIR", "PTPRE", "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16", "RHOG", "RIPK2", "RNF144B", "ROS1", "RTP4", "SCARF1", "SCN1B", "SELE", "SELENOS", "SELL", "SEMA4D", "SERPINE1", "SGMS2", "SLAMF1", "SLC11A2", "SLC1A2", "SLC28A2", "SLC31A1", "SLC31A2", "SLC4A4", "SLC7A1", "SLC7A2", "SPHK1", "SRI", "STAB1", "TACR1", "TACR3", "TAPBP", "TIMP1", "TLR1", "TLR2", "TLR3", "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10", "TNFSF15", "TNFSF9", "TPBG", "VIP"))
inflammation <- intersect(all.genes, inflammation)

seurat.object2 <- AddModuleScore(object = seurat.object2, features = inflammation, ctrl = 5, name = 'inflammation.score')
# p1 <- FeaturePlot(seurat.object2, features = 'inflammation.score1')
# ggsave("FeaturePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'inflammation.score1')
ggsave("RidgePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "inflammation.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "inflammation.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "inflammation.score1", images=image.list, ncol=4)
p5
ggsave("inflammation.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

ifn.alpha <- c("ADAR", "B2M", "BATF2", "BST2", "C1S", "CASP1", "CASP8", "CCRL2", "CD47", "CD74", "CMPK2", "CMTR1", "CNP", "CSF1", "CXCL10", "CXCL11", "DDX60", "DHX58", "EIF2AK2", "ELF1", "EPSTI1", "GBP2", "GBP4", "GMPR", "HELZ2", "HERC6", "HLA-C", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3", "IL15", "IL4R", "IL7", "IRF1", "IRF2", "IRF7", "IRF9", "ISG15", "ISG20", "LAMP3", "LAP3", "LGALS3BP", "LPAR6", "LY6E", "MOV10", "MVB12A", "MX1", "NCOA7", "NMI", "NUB1", "OAS1", "OASL", "OGFR", "PARP12", "PARP14", "PARP9", "PLSCR1", "PNPT1", "PROCR", "PSMA3", "PSMB8", "PSMB9", "PSME1", "PSME2", "RIPK2", "RNF31", "RSAD2", "RTP4", "SAMD9", "SAMD9L", "SELL", "SLC25A28", "SP110", "STAT2", "TAP1", "TDRD7", "TENT5A", "TMEM140", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TRIM5", "TXNIP", "UBA7", "UBE2L6", "USP18", "WARS1"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = ifn.alpha, ctrl = 5, name = 'ifn.alpha.score')
# p1 <- FeaturePlot(seurat.object2, features = 'ifn.alpha.score1')
# ggsave("FeaturePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'ifn.alpha.score1')
ggsave("RidgePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "ifn.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "ifn.alpha.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

ifn.beta <- c("ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM", "BST2", "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1", "CMPK2", "CMTR1", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX58", "DDX60", "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "GBP4", "GBP6", "GCH1", "GPR18", "GZMA", "HELZ2", "HERC6", "HIF1A", "HLA-A", "HLA-B", "HLA-DMA", "HLA-DQA1", "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IFNAR2", "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1", "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7", "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1", "METTL7B", "MT2A", "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12", "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP", "PNPT1", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1", "RIPK2", "RNF213", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP", "SERPING1", "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2", "SP110", "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4", "STAT1", "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP", "TDRD7", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8", "VCAM1", "WARS1", "XAF1", "XCL1", "ZBP1", "ZNFX1"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = ifn.beta, ctrl = 5, name = 'ifn.beta.score')
# p1 <- FeaturePlot(seurat.object2, features = 'ifn.beta.score1')
# ggsave("FeaturePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'ifn.beta.score1')
ggsave("RidgePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "ifn.beta.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "ifn.beta.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tgfb <- c("ACVR1", "APC", "ARID4B", "BCAR3", "BMP2", "BMPR1A", "BMPR2", "CDH1", "CDK9", "CDKN1C", "CTNNB1", "ENG", "FKBP1A", "FNTA", "FURIN", "HDAC1", "HIPK2", "ID1", "ID2", "ID3", "IFNGR2", "JUNB", "KLF10", "LEFTY2", "LTBP2", "MAP3K7", "NCOR2", "NOG", "PMEPA1", "PPM1A", "PPP1CA", "PPP1R15A", "RAB31", "RHOA", "SERPINE1", "SKI", "SKIL", "SLC20A1", "SMAD1", "SMAD3", "SMAD6", "SMAD7", "SMURF1", "SMURF2", "SPTBN1", "TGFB1", "TGFBR1", "TGIF1", "THBS1", "TJP1", "TRIM33", "UBE2D3", "WWTR1", "XIAP"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = tgfb, ctrl = 5, name = 'tgfb.score')
# p1 <- FeaturePlot(seurat.object2, features = 'tgfb.score1')
# ggsave("FeaturePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'tgfb.score1')
ggsave("RidgePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "tgfb.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "tgfb.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

wnt <- c("ADAM17", "AXIN1", "AXIN2", "CCND2", "CSNK1E", "CTNNB1", "CUL1", "DKK1", "DKK4", "DLL1", "DVL2", "FRAT1", "FZD1", "FZD8", "GNAI1", "HDAC11", "HDAC2", "HDAC5", "HEY1", "HEY2", "JAG1", "JAG2", "KAT2A", "LEF1", "MAML1", "MYC", "NCOR2", "NCSTN", "NKD1", "NOTCH1", "NOTCH4", "NUMB", "PPARD", "PSEN2", "PTCH1", "RBPJ", "SKP2", "TCF7", "TP53", "WNT1", "WNT5B", "WNT6"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = wnt, ctrl = 5, name = 'wnt.score')
# p1 <- FeaturePlot(seurat.object2, features = 'wnt.score1')
# ggsave("FeaturePlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'wnt.score1')
ggsave("RidgePlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "wnt.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "wnt.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

il2 <- c("ABCB1", "ADAM19", "AGER", "AHCY", "AHNAK", "AHR", "ALCAM", "AMACR", "ANXA4", "APLP1", "ARL4A", "BATF", "BATF3", "BCL2", "BCL2L1", "BHLHE40", "BMP2", "BMPR2", "CA2", "CAPG", "CAPN3", "CASP3", "CCND2", "CCND3", "CCNE1", "CCR4", "CD44", "CD48", "CD79B", "CD81", "CD83", "CD86", "CDC42SE2", "CDC6", "CDCP1", "CDKN1C", "CISH", "CKAP4", "COCH", "COL6A1", "CSF1", "CSF2", "CST7", "CTLA4", "CTSZ", "CXCL10", "CYFIP1", "DCPS", "DENND5A", "DHRS3", "DRC1", "ECM1", "EEF1AKMT1", "EMP1", "ENO3", "ENPP1", "EOMES", "ETFBKMT", "ETV4", "F2RL2", "FAH", "FAM126B", "FGL2", "FLT3LG", "FURIN", "GABARAPL1", "GADD45B", "GALM", "GATA1", "GBP4", "GLIPR2", "GPR65", "GPR83", "GPX4", "GSTO1", "GUCY1B1", "HIPK2", "HK2", "HOPX", "HUWE1", "ICOS", "IFITM3", "IFNGR1", "IGF1R", "IGF2R", "IKZF2", "IKZF4", "IL10", "IL10RA", "IL13", "IL18R1", "IL1R2", "IL1RL1", "IL2RA", "IL2RB", "IL3RA", "IL4R", "IRF4", "IRF6", "IRF8", "ITGA6", "ITGAE", "ITGAV", "ITIH5", "KLF6", "LCLAT1", "LIF", "LRIG1", "LRRC8C", "LTB", "MAFF", "MAP3K8", "MAP6", "MAPKAPK2", "MUC1", "MXD1", "MYC", "MYO1C", "MYO1E", "NCOA3", "NCS1", "NDRG1", "NFIL3", "NFKBIZ", "NOP2", "NRP1", "NT5E", "ODC1", "P2RX4", "P4HA1", "PDCD2L", "PENK", "PHLDA1", "PHTF2", "PIM1", "PLAGL1", "PLEC", "PLIN2", "PLPP1", "PLSCR1", "PNP", "POU2F1", "PRAF2", "PRKCH", "PRNP", "PTCH1", "PTGER2", "PTH1R", "PTRH2", "PUS1", "RABGAP1L", "RGS16", "RHOB", "RHOH", "RNH1", "RORA", "RRAGD", "S100A1", "SCN9A", "SELL", "SELP", "SERPINB6", "SERPINC1", "SH3BGRL2", "SHE", "SLC1A5", "SLC29A2", "SLC2A3", "SLC39A8", "SMPDL3A", "SNX14", "SNX9", "SOCS1", "SOCS2", "SPP1", "SPRED2", "SPRY4", "ST3GAL4", "SWAP70", "SYNGR2", "SYT11", "TGM2", "TIAM1", "TLR7", "TNFRSF18", "TNFRSF1B", "TNFRSF21", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFSF10", "TNFSF11", "TRAF1", "TTC39B", "TWSG1", "UCK2", "UMPS", "WLS", "XBP1"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = il2, ctrl = 5, name = 'il2.score')
# p1 <- FeaturePlot(seurat.object2, features = 'il2.score1')
# ggsave("FeaturePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'il2.score1')
ggsave("RidgePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "il2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "il2.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

il6 <- c("A2M", "ACVR1B", "ACVRL1", "BAK1", "CBL", "CCL7", "CCR1", "CD14", "CD36", "CD38", "CD44", "CD9", "CNTFR", "CRLF2", "CSF1", "CSF2", "CSF2RA", "CSF2RB", "CSF3R", "CXCL1", "CXCL10", "CXCL11", "CXCL13", "CXCL3", "CXCL9", "DNTT", "EBI3", "FAS", "GRB2", "HAX1", "HMOX1", "IFNAR1", "IFNGR1", "IFNGR2", "IL10RB", "IL12RB1", "IL13RA1", "IL15RA", "IL17RA", "IL17RB", "IL18R1", "IL1B", "IL1R1", "IL1R2", "IL2RA", "IL2RG", "IL3RA", "IL4R", "IL6", "IL6ST", "IL7", "IL9R", "INHBE", "IRF1", "IRF9", "ITGA4", "ITGB3", "JUN", "LEPR", "LTB", "LTBR", "MAP3K8", "MYD88", "OSMR", "PDGFC", "PF4", "PIK3R5", "PIM1", "PLA2G2A", "PTPN1", "PTPN11", "PTPN2", "REG1A", "SOCS1", "SOCS3", "STAM2", "STAT1", "STAT2", "STAT3", "TGFB1", "TLR2", "TNF", "TNFRSF12A", "TNFRSF1A", "TNFRSF1B", "TNFRSF21", "TYK2"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = il6, ctrl = 5, name = 'il6.score')
# p1 <- FeaturePlot(seurat.object2, features = 'il6.score1')
# ggsave("FeaturePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'il6.score1')
ggsave("RidgePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "il6.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "il6.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

apoptosis <- c("ADD1", "AIFM3", "ANKH", "ANXA1", "APP", "ATF3", "AVPR1A", "BAX", "BCAP31", "BCL10", "BCL2L1", "BCL2L10", "BCL2L11", "BCL2L2", "BGN", "BID", "BIK", "BIRC3", "BMF", "BMP2", "BNIP3L", "BRCA1", "BTG2", "BTG3", "CASP1", "CASP2", "CASP3", "CASP4", "CASP6", "CASP7", "CASP8", "CASP9", "CAV1", "CCNA1", "CCND1", "CCND2", "CD14", "CD2", "CD38", "CD44", "CD69", "CDC25B", "CDK2", "CDKN1A", "CDKN1B", "CFLAR", "CLU", "CREBBP", "CTH", "CTNNB1", "CYLD", "DAP", "DAP3", "DCN", "DDIT3", "DFFA", "DIABLO", "DNAJA1", "DNAJC3", "DNM1L", "DPYD", "EBP", "EGR3", "EMP1", "ENO2", "ERBB2", "ERBB3", "EREG", "ETF1", "F2", "F2R", "FAS", "FASLG", "FDXR", "FEZ1", "GADD45A", "GADD45B", "GCH1", "GNA15", "GPX1", "GPX3", "GPX4", "GSN", "GSR", "GSTM1", "GUCY2D", "H1-0", "HGF", "HMGB2", "HMOX1", "HSPB1", "IER3", "IFITM3", "IFNB1", "IFNGR1", "IGF2R", "IGFBP6", "IL18", "IL1A", "IL1B", "IL6", "IRF1", "ISG20", "JUN", "KRT18", "LEF1", "LGALS3", "LMNA", "LUM", "MADD", "MCL1", "MGMT", "MMP2", "NEDD9", "NEFH", "PAK1", "PDCD4", "PDGFRB", "PEA15", "PLAT", "PLCB2", "PLPPR4", "PMAIP1", "PPP2R5B", "PPP3R1", "PPT1", "PRF1", "PSEN1", "PSEN2", "PTK2", "RARA", "RELA", "RETSAT", "RHOB", "RHOT2", "RNASEL", "ROCK1", "SAT1", "SATB1", "SC5D", "SLC20A1", "SMAD7", "SOD1", "SOD2", "SPTAN1", "SQSTM1", "TAP1", "TGFB2", "TGFBR3", "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF12A", "TNFSF10", "TOP2A", "TSPO", "TXNIP", "VDAC2", "WEE1", "XIAP"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = apoptosis, ctrl = 5, name = 'apoptosis.score')
# p1 <- FeaturePlot(seurat.object2, features = 'apoptosis.score1')
# ggsave("FeaturePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'apoptosis.score1')
ggsave("RidgePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "apoptosis.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "apoptosis.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tnf.alpha <- c("ABCA1", "ACKR3", "AREG", "ATF3", "ATP2B1", "B4GALT1", "B4GALT5", "BCL2A1", "BCL3", "BCL6", "BHLHE40", "BIRC2", "BIRC3", "BMP2", "BTG1", "BTG2", "BTG3", "CCL2", "CCL20", "CCL4", "CCL5", "CCN1", "CCND1", "CCNL1", "CCRL2", "CD44", "CD69", "CD80", "CD83", "CDKN1A", "CEBPB", "CEBPD", "CFLAR", "CLCF1", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL11", "CXCL2", "CXCL3", "CXCL6", "DDX58", "DENND5A", "DNAJB4", "DRAM1", "DUSP1", "DUSP2", "DUSP4", "DUSP5", "EDN1", "EFNA1", "EGR1", "EGR2", "EGR3", "EHD1", "EIF1", "ETS2", "F2RL1", "F3", "FJX1", "FOS", "FOSB", "FOSL1", "FOSL2", "FUT4", "G0S2", "GADD45A", "GADD45B", "GCH1", "GEM", "GFPT2", "GPR183", "HBEGF", "HES1", "ICAM1", "ICOSLG", "ID2", "IER2", "IER3", "IER5", "IFIH1", "IFIT2", "IFNGR2", "IL12B", "IL15RA", "IL18", "IL1A", "IL1B", "IL23A", "IL6", "IL6ST", "IL7R", "INHBA", "IRF1", "IRS2", "JAG1", "JUN", "JUNB", "KDM6B", "KLF10", "KLF2", "KLF4", "KLF6", "KLF9", "KYNU", "LAMB3", "LDLR", "LIF", "LITAF", "MAFF", "MAP2K3", "MAP3K8", "MARCKS", "MCL1", "MSC", "MXD1", "MYC", "NAMPT", "NFAT5", "NFE2L2", "NFIL3", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NINJ1", "NR4A1", "NR4A2", "NR4A3", "OLR1", "PANX1", "PDE4B", "PDLIM5", "PER1", "PFKFB3", "PHLDA1", "PHLDA2", "PLAU", "PLAUR", "PLEK", "PLK2", "PLPP3", "PMEPA1", "PNRC1", "PPP1R15A", "PTGER4", "PTGS2", "PTPRE", "PTX3", "RCAN1", "REL", "RELA", "RELB", "RHOB", "RIPK2", "RNF19B", "SAT1", "SDC4", "SERPINB2", "SERPINB8", "SERPINE1", "SGK1", "SIK1", "SLC16A6", "SLC2A3", "SLC2A6", "SMAD3", "SNN", "SOCS3", "SOD2", "SPHK1", "SPSB1", "SQSTM1", "STAT5A", "TANK", "TAP1", "TGIF1", "TIPARP", "TLR2", "TNC", "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFRSF9", "TNFSF9", "TNIP1", "TNIP2", "TRAF1", "TRIB1", "TRIP10", "TSC22D1", "TUBB2A", "VEGFA", "YRDC", "ZBTB10", "ZC3H12A", "ZFP36"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = tnf.alpha, ctrl = 5, name = 'tnf.alpha.score')
# p1 <- FeaturePlot(seurat.object2, features = 'tnf.alpha.score1')
# ggsave("FeaturePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'tnf.alpha.score1')
ggsave("RidgePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "tnf.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "tnf.alpha.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

setwd("..")

dir.create("GSEA_Type.Profile")
setwd("GSEA_Type.Profile")
Idents(seurat.object2) <- "Type.Profile"

## use RNA for all DE analysis and plots
DefaultAssay(seurat.object2) <- "Spatial"
# seurat.object2 <- ScaleData(seurat.object2, features = all.genes, assay = "Spatial")
seurat.object2 <- ScaleData(seurat.object2, features = all.genes.visium, assay="Spatial")
hypoxia.list <- c("ACKR3", "ADM", "ADORA2B", "AK4", "AKAP12", "ALDOA", "ALDOB", "ALDOC", "AMPD3", "ANGPTL4", "ANKZF1", "ANXA2", "ATF3", "ATP7A", "B3GALT6", "B4GALNT2", "BCAN", "BCL2", "BGN", "BHLHE40", "BNIP3L", "BRS3", "BTG1", "CA12", "CASP6", "CAV1", "CAVIN1", "CAVIN3", "CCN1", "CCN2", "CCN5", "CCNG2", "CDKN1A", "CDKN1B", "CDKN1C", "CHST2", "CHST3", "CITED2", "COL5A1", "CP", "CSRP2", "CXCR4", "DCN", "DDIT3", "DDIT4", "DPYSL4", "DTNA", "DUSP1", "EDN2", "EFNA1", "EFNA3", "EGFR", "ENO1", "ENO2", "ENO3", "ERO1A", "ERRFI1", "ETS1", "EXT1", "F3", "FAM162A", "FBP1", "FOS", "FOSL2", "FOXO3", "GAA", "GALK1", "GAPDH", "GAPDHS", "GBE1", "GCK", "GCNT2", "GLRX", "GPC1", "GPC3", "GPC4", "GPI", "GRHPR", "GYS1", "HAS1", "HDLBP", "HEXA", "HK1", "HK2", "HMOX1", "HOXB9", "HS3ST1", "HSPA5", "IDS", "IER3", "IGFBP1", "IGFBP3", "IL6", "ILVBL", "INHA", "IRS2", "ISG20", "JMJD6", "JUN", "KDELR3", "KDM3A", "KIF5A", "KLF6", "KLF7", "KLHL24", "LALBA", "LARGE1", "LDHA", "LDHC", "LOX", "LXN", "MAFF", "MAP3K1", "MIF", "MT1E", "MT2A", "MXI1", "MYH9", "NAGK", "NCAN", "NDRG1", "NDST1", "NDST2", "NEDD4L", "NFIL3", "NOCT", "NR3C1", "P4HA1", "P4HA2", "PAM", "PCK1", "PDGFB", "PDK1", "PDK3", "PFKFB3", "PFKL", "PFKP", "PGAM2", "PGF", "PGK1", "PGM1", "PGM2", "PHKG1", "PIM1", "PKLR", "PKP1", "PLAC8", "PLAUR", "PLIN2", "PNRC1", "PPARGC1A", "PPFIA4", "PPP1R15A", "PPP1R3C", "PRDX5", "PRKCA", "PYGM", "RBPJ", "RORA", "RRAGD", "S100A4", "SAP30", "SCARB1", "SDC2", "SDC3", "SDC4", "SELENBP1", "SERPINE1", "SIAH2", "SLC25A1", "SLC2A1", "SLC2A3", "SLC2A5", "SLC37A4", "SLC6A6", "SRPX", "STBD1", "STC1", "STC2", "SULT2B1", "TES", "TGFB3", "TGFBI", "TGM2", "TIPARP", "TKTL1", "TMEM45A", "TNFAIP3", "TPBG", "TPD52", "TPI1", "TPST2", "UGP2", "VEGFA", "VHL", "VLDLR", "WSB1", "XPNPEP1", "ZFP36", "ZNF292")
hypoxia.list <- intersect(all.genes, hypoxia.list)

seurat.object2 <- AddModuleScore(object = seurat.object2, features = hypoxia.list, ctrl = 5, name = 'hypoxia.score')
# p1 <- FeaturePlot(seurat.object2, features = 'hypoxia.score1')
# ggsave("FeaturePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'hypoxia.score1')
ggsave("RidgePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "hypoxia.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "hypoxia.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "hypoxia.score1", images=image.list, ncol=4)
p5
ggsave("hypoxia.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


complement.list <- c("ACTN2", "ADAM9", "ADRA2B", "AKAP10", "ANG", "ANXA5", "APOA4", "APOBEC3F", "APOBEC3G", "APOC1", "ATOX1", "BRPF3", "C1QA", "C1QC", "C1R", "C1S", "C2", "C3", "C4BPB", "C9", "CA2", "CALM1", "CALM3", "CASP1", "CASP10", "CASP3", "CASP4", "CASP5", "CASP7", "CASP9", "CBLB", "CCL5", "CD36", "CD40LG", "CD46", "CD55", "CD59", "CDA", "CDH13", "CDK5R1", "CEBPB", "CFB", "CFH", "CLU", "COL4A2", "CP", "CPM", "CPQ", "CR1", "CR2", "CSRP1", "CTSB", "CTSC", "CTSD", "CTSH", "CTSL", "CTSO", "CTSS", "CTSV", "CXCL1", "DGKG", "DGKH", "DOCK10", "DOCK4", "DOCK9", "DPP4", "DUSP5", "DUSP6", "DYRK2", "EHD1", "ERAP2", "F10", "F2", "F3", "F5", "F7", "F8", "FCER1G", "FCN1", "FDX1", "FN1", "FYN", "GATA3", "GCA", "GMFB", "GNAI2", "GNAI3", "GNB2", "GNB4", "GNG2", "GNGT2", "GP1BA", "GP9", "GPD2", "GRB2", "GZMA", "GZMB", "GZMK", "HNF4A", "HPCAL4", "HSPA1A", "HSPA5", "IL6", "IRF1", "IRF2", "IRF7", "ITGAM", "ITIH1", "JAK2", "KCNIP2", "KCNIP3", "KIF2A", "KLK1", "KLKB1", "KYNU", "L3MBTL4", "LAMP2", "LAP3", "LCK", "LCP2", "LGALS3", "LGMN", "LIPA", "LRP1", "LTA4H", "LTF", "LYN", "MAFF", "ME1", "MMP12", "MMP13", "MMP14", "MMP15", "MMP8", "MSRB1", "MT3", "NOTCH4", "OLR1", "PCLO", "PCSK9", "PDGFB", "PDP1", "PFN1", "PHEX", "PIK3CA", "PIK3CG", "PIK3R5", "PIM1", "PLA2G4A", "PLA2G7", "PLAT", "PLAUR", "PLEK", "PLG", "PLSCR1", "PPP2CB", "PPP4C", "PRCP", "PRDM4", "PREP", "PRKCD", "PRSS3", "PRSS36", "PSEN1", "PSMB9", "RABIF", "RAF1", "RASGRP1", "RBSN", "RCE1", "RHOG", "RNF4", "S100A12", "S100A13", "S100A9", "SCG3", "SERPINA1", "SERPINB2", "SERPINC1", "SERPINE1", "SERPING1", "SH2B3", "SIRT6", "SPOCK2", "SRC", "STX4", "TFPI2", "TIMP1", "TIMP2", "TMPRSS6", "TNFAIP3", "USP14", "USP15", "USP16", "USP8", "VCPIP1", "WAS", "XPNPEP1", "ZEB1", "ZFPM2")
complement.list <- intersect(all.genes, complement.list)
seurat.object2 <- AddModuleScore(object = seurat.object2, features = complement.list, ctrl = 5, name = 'complement.score')
# p1 <- FeaturePlot(seurat.object2, features = 'complement.score1')
# ggsave("FeaturePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'complement.score1')
ggsave("RidgePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "complement.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "complement.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "complement.score1", images=image.list, ncol=4)
p5
ggsave("complement.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


# http://www.gsea-msigdb.org/gsea/msigdb/cards/COATES_MACROPHAGE_M1_VS_M2_UP.html
M1 <- c("ABCA9", "ABHD1", "ACSS1", "ACYP2", "ADA", "AIF1", "ALAD", "ARFGEF3", "ARSB", "ATP6V0E2", "AXL", "BCKDHB", "BLOC1S6", "CADM1", "CAP1", "CAPN5", "CBX6", "CD59", "CFH", "CLBA1", "CNRIP1", "COLEC12", "COMT", "CRIM1", "CXCL14", "CXCR4", "DST", "DYNLT1", "EMC1", "ENO2", "FAM124A", "FAM135A", "FAM9B", "FGD2", "FILIP1L", "GALNT11", "GATM", "GDA", "GJA1", "GLO1", "GNB4", "HAUS2", "HDDC3", "HLA-DQA1", "HMGN3", "KCNJ10", "LAMA3", "LCORL", "LYPLAL1", "MAF", "MALAT1", "MARCKSL1", "MARCO", "MSR1", "NAT8L", "NRCAM", "OCEL1", "OGFRL1", "P2RY13", "PIANP", "PIK3AP1", "PLAAT3", "PLBD1", "PLXDC2", "PPP2R5C", "PTGER3", "RAB10", "RAPSN", "RASAL2", "RCBTB2", "RCN1", "RFX3", "RPL14", "SFI1", "SLC35A1", "SLC7A7", "SLCO2B1", "SRD5A3", "TGFBI", "TIFAB", "TM7SF3", "TOR3A", "TTC3", "TUBB2B", "TXNIP", "ZNF727")
M1 <- intersect(all.genes, M1)

seurat.object2 <- AddModuleScore(object = seurat.object2, features = M1, ctrl = 5, name = 'M1.score')
# p1 <- FeaturePlot(seurat.object2, features = 'M1.score1')
# ggsave("FeaturePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'M1.score1')
ggsave("RidgePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "M1.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "M1.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "M1.score1", images=image.list, ncol=4)
p5
ggsave("M1.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

M2<- c("ADAM8", "AK4", "AKR1E2", "ALDOC", "ANG", "ANKRD37", "ANOS1", "ARG1", "ARSK", "ATG4D", "ATP6V0D2", "BNIP3", "BST1", "C1QB", "C1QC", "C1QTNF12", "C5orf34", "CBLB", "CD24", "CD300C", "CD300LD", "CD5L", "CHIA", "COL20A1", "CRIPT", "CTSK", "CTSL", "CYTIP", "EFHD2", "EGLN3", "ERO1A", "F10", "F7", "FAM177A1", "FAM199X", "FAM241A", "FBLIM1", "FLRT2", "GASK1B", "GDF15", "GPRC5B", "HILPDA", "HLA-A", "HLA-B", "HS6ST1", "HSPA1A", "HSPA1B", "IGF2R", "IGHM", "IL18BP", "ITGB3", "LBP", "LIN7C", "LRRC27", "MCOLN3", "MFGE8", "MGST2", "MYO1F", "PADI4", "PDXDC1", "PINK1", "PKDCC", "PRDX2", "PROCR", "PTGER2", "QRICH1", "RGCC", "RIMBP2", "RPGRIP1", "RRAS2", "SCD", "SH2B2", "SLC2A1", "SNHG6", "SOAT1", "THBS1", "TJP2", "TLR1", "TMEM267", "TULP4", "UCHL1", "VEGFA", "XDH")
M2 <- intersect(all.genes, M2)
seurat.object2 <- AddModuleScore(object = seurat.object2, features = M2, ctrl = 5, name = 'M2.score')
# p1 <- FeaturePlot(seurat.object2, features = 'M2.score1')
# ggsave("FeaturePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'M2.score1')
ggsave("RidgePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "M2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "M2.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "M2.score1", images=image.list, ncol=4)
p5
ggsave("M2.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

## http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_PI3K_AKT_MTOR_SIGNALING.html
pi3k <- c("ACACA", "ACTR2", "ACTR3", "ADCY2", "AKT1", "AKT1S1", "AP2M1", "ARF1", "ARHGDIA", "ARPC3", "ATF1", "CAB39", "CAB39L", "CALR", "CAMK4", "CDK1", "CDK2", "CDK4", "CDKN1A", "CDKN1B", "CFL1", "CLTC", "CSNK2B", "CXCR4", "DAPP1", "DDIT3", "DUSP3", "E2F1", "ECSIT", "EGFR", "EIF4E", "FASLG", "FGF17", "FGF22", "FGF6", "GNA14", "GNGT1", "GRB2", "GRK2", "GSK3B", "HRAS", "HSP90B1", "IL2RG", "IL4", "IRAK4", "ITPR2", "LCK", "MAP2K3", "MAP2K6", "MAP3K7", "MAPK1", "MAPK10", "MAPK8", "MAPK9", "MAPKAP1", "MKNK1", "MKNK2", "MYD88", "NCK1", "NFKBIB", "NGF", "NOD1", "PAK4", "PDK1", "PFN1", "PIK3R3", "PIKFYVE", "PIN1", "PITX2", "PLA2G12A", "PLCB1", "PLCG1", "PPP1CA", "PPP2R1B", "PRKAA2", "PRKAG1", "PRKAR2A", "PRKCB", "PTEN", "PTPN11", "RAC1", "RAF1", "RALB", "RIPK1", "RIT1", "RPS6KA1", "RPS6KA3", "RPTOR", "SFN", "SLA", "SLC2A1", "SMAD2", "SQSTM1", "STAT2", "TBK1", "THEM4", "TIAM1", "TNFRSF1A", "TRAF2", "TRIB3", "TSC2", "UBE2D3", "UBE2N", "VAV3", "YWHAB"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = pi3k, ctrl = 5, name = 'pi3k.score')
# p1 <- FeaturePlot(seurat.object2, features = 'pi3k.score1')
# ggsave("FeaturePlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'pi3k.score1')
ggsave("RidgePlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "pi3k.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "pi3k.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_pi3k.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

inflammation <- c("ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "ADM", "ADORA2B", "ADRM1", "AHR", "APLNR", "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1", "BEST1", "BST2", "BTG2", "C3AR1", "C5AR1", "CALCRL", "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7", "CCRL2", "CD14", "CD40", "CD48", "CD55", "CD69", "CD70", "CD82", "CDKN1A", "CHST2", "CLEC5A", "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1", "CXCL10", "CXCL11", "CXCL6", "CXCL8", "CXCL9", "CXCR6", "CYBB", "DCBLD2", "EBI3", "EDN1", "EIF2AK2", "EMP3", "EREG", "F3", "FFAR2", "FPR1", "FZD5", "GABBR1", "GCH1", "GNA15", "GNAI3", "GP1BA", "GPC3", "GPR132", "GPR183", "HAS2", "HBEGF", "HIF1A", "HPN", "HRH1", "ICAM1", "ICAM4", "ICOSLG", "IFITM1", "IFNAR1", "IFNGR2", "IL10", "IL10RA", "IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5", "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2", "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV", "MEP1A", "MET", "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA", "NLRP3", "NMI", "NMUR1", "NOD2", "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2", "PCDH7", "PDE4B", "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2", "PTGER4", "PTGIR", "PTPRE", "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16", "RHOG", "RIPK2", "RNF144B", "ROS1", "RTP4", "SCARF1", "SCN1B", "SELE", "SELENOS", "SELL", "SEMA4D", "SERPINE1", "SGMS2", "SLAMF1", "SLC11A2", "SLC1A2", "SLC28A2", "SLC31A1", "SLC31A2", "SLC4A4", "SLC7A1", "SLC7A2", "SPHK1", "SRI", "STAB1", "TACR1", "TACR3", "TAPBP", "TIMP1", "TLR1", "TLR2", "TLR3", "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10", "TNFSF15", "TNFSF9", "TPBG", "VIP"))
inflammation <- intersect(inflammation, all.genes)
seurat.object2 <- AddModuleScore(object = seurat.object2, features = inflammation, ctrl = 5, name = 'inflammation.score')
# p1 <- FeaturePlot(seurat.object2, features = 'inflammation.score1')
# ggsave("FeaturePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'inflammation.score1')
ggsave("RidgePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "inflammation.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "inflammation.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "inflammation.score1", images=image.list, ncol=4)
p5
ggsave("inflammation.score1SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

ifn.alpha <- c("ADAR", "B2M", "BATF2", "BST2", "C1S", "CASP1", "CASP8", "CCRL2", "CD47", "CD74", "CMPK2", "CMTR1", "CNP", "CSF1", "CXCL10", "CXCL11", "DDX60", "DHX58", "EIF2AK2", "ELF1", "EPSTI1", "GBP2", "GBP4", "GMPR", "HELZ2", "HERC6", "HLA-C", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3", "IL15", "IL4R", "IL7", "IRF1", "IRF2", "IRF7", "IRF9", "ISG15", "ISG20", "LAMP3", "LAP3", "LGALS3BP", "LPAR6", "LY6E", "MOV10", "MVB12A", "MX1", "NCOA7", "NMI", "NUB1", "OAS1", "OASL", "OGFR", "PARP12", "PARP14", "PARP9", "PLSCR1", "PNPT1", "PROCR", "PSMA3", "PSMB8", "PSMB9", "PSME1", "PSME2", "RIPK2", "RNF31", "RSAD2", "RTP4", "SAMD9", "SAMD9L", "SELL", "SLC25A28", "SP110", "STAT2", "TAP1", "TDRD7", "TENT5A", "TMEM140", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TRIM5", "TXNIP", "UBA7", "UBE2L6", "USP18", "WARS1"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = ifn.alpha, ctrl = 5, name = 'ifn.alpha.score')
# p1 <- FeaturePlot(seurat.object2, features = 'ifn.alpha.score1')
# ggsave("FeaturePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'ifn.alpha.score1')
ggsave("RidgePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "ifn.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "ifn.alpha.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

ifn.beta <- c("ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM", "BST2", "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1", "CMPK2", "CMTR1", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX58", "DDX60", "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "GBP4", "GBP6", "GCH1", "GPR18", "GZMA", "HELZ2", "HERC6", "HIF1A", "HLA-A", "HLA-B", "HLA-DMA", "HLA-DQA1", "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IFNAR2", "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1", "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7", "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1", "METTL7B", "MT2A", "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12", "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP", "PNPT1", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1", "RIPK2", "RNF213", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP", "SERPING1", "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2", "SP110", "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4", "STAT1", "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP", "TDRD7", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8", "VCAM1", "WARS1", "XAF1", "XCL1", "ZBP1", "ZNFX1"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = ifn.beta, ctrl = 5, name = 'ifn.beta.score')
# p1 <- FeaturePlot(seurat.object2, features = 'ifn.beta.score1')
# ggsave("FeaturePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'ifn.beta.score1')
ggsave("RidgePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "ifn.beta.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "ifn.beta.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tgfb <- c("ACVR1", "APC", "ARID4B", "BCAR3", "BMP2", "BMPR1A", "BMPR2", "CDH1", "CDK9", "CDKN1C", "CTNNB1", "ENG", "FKBP1A", "FNTA", "FURIN", "HDAC1", "HIPK2", "ID1", "ID2", "ID3", "IFNGR2", "JUNB", "KLF10", "LEFTY2", "LTBP2", "MAP3K7", "NCOR2", "NOG", "PMEPA1", "PPM1A", "PPP1CA", "PPP1R15A", "RAB31", "RHOA", "SERPINE1", "SKI", "SKIL", "SLC20A1", "SMAD1", "SMAD3", "SMAD6", "SMAD7", "SMURF1", "SMURF2", "SPTBN1", "TGFB1", "TGFBR1", "TGIF1", "THBS1", "TJP1", "TRIM33", "UBE2D3", "WWTR1", "XIAP"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = tgfb, ctrl = 5, name = 'tgfb.score')
# p1 <- FeaturePlot(seurat.object2, features = 'tgfb.score1')
# ggsave("FeaturePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'tgfb.score1')
ggsave("RidgePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "tgfb.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "tgfb.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

wnt <- c("ADAM17", "AXIN1", "AXIN2", "CCND2", "CSNK1E", "CTNNB1", "CUL1", "DKK1", "DKK4", "DLL1", "DVL2", "FRAT1", "FZD1", "FZD8", "GNAI1", "HDAC11", "HDAC2", "HDAC5", "HEY1", "HEY2", "JAG1", "JAG2", "KAT2A", "LEF1", "MAML1", "MYC", "NCOR2", "NCSTN", "NKD1", "NOTCH1", "NOTCH4", "NUMB", "PPARD", "PSEN2", "PTCH1", "RBPJ", "SKP2", "TCF7", "TP53", "WNT1", "WNT5B", "WNT6"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = wnt, ctrl = 5, name = 'wnt.score')
# p1 <- FeaturePlot(seurat.object2, features = 'wnt.score1')
# ggsave("FeaturePlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'wnt.score1')
ggsave("RidgePlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "wnt.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "wnt.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_wnt.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

il2 <- c("ABCB1", "ADAM19", "AGER", "AHCY", "AHNAK", "AHR", "ALCAM", "AMACR", "ANXA4", "APLP1", "ARL4A", "BATF", "BATF3", "BCL2", "BCL2L1", "BHLHE40", "BMP2", "BMPR2", "CA2", "CAPG", "CAPN3", "CASP3", "CCND2", "CCND3", "CCNE1", "CCR4", "CD44", "CD48", "CD79B", "CD81", "CD83", "CD86", "CDC42SE2", "CDC6", "CDCP1", "CDKN1C", "CISH", "CKAP4", "COCH", "COL6A1", "CSF1", "CSF2", "CST7", "CTLA4", "CTSZ", "CXCL10", "CYFIP1", "DCPS", "DENND5A", "DHRS3", "DRC1", "ECM1", "EEF1AKMT1", "EMP1", "ENO3", "ENPP1", "EOMES", "ETFBKMT", "ETV4", "F2RL2", "FAH", "FAM126B", "FGL2", "FLT3LG", "FURIN", "GABARAPL1", "GADD45B", "GALM", "GATA1", "GBP4", "GLIPR2", "GPR65", "GPR83", "GPX4", "GSTO1", "GUCY1B1", "HIPK2", "HK2", "HOPX", "HUWE1", "ICOS", "IFITM3", "IFNGR1", "IGF1R", "IGF2R", "IKZF2", "IKZF4", "IL10", "IL10RA", "IL13", "IL18R1", "IL1R2", "IL1RL1", "IL2RA", "IL2RB", "IL3RA", "IL4R", "IRF4", "IRF6", "IRF8", "ITGA6", "ITGAE", "ITGAV", "ITIH5", "KLF6", "LCLAT1", "LIF", "LRIG1", "LRRC8C", "LTB", "MAFF", "MAP3K8", "MAP6", "MAPKAPK2", "MUC1", "MXD1", "MYC", "MYO1C", "MYO1E", "NCOA3", "NCS1", "NDRG1", "NFIL3", "NFKBIZ", "NOP2", "NRP1", "NT5E", "ODC1", "P2RX4", "P4HA1", "PDCD2L", "PENK", "PHLDA1", "PHTF2", "PIM1", "PLAGL1", "PLEC", "PLIN2", "PLPP1", "PLSCR1", "PNP", "POU2F1", "PRAF2", "PRKCH", "PRNP", "PTCH1", "PTGER2", "PTH1R", "PTRH2", "PUS1", "RABGAP1L", "RGS16", "RHOB", "RHOH", "RNH1", "RORA", "RRAGD", "S100A1", "SCN9A", "SELL", "SELP", "SERPINB6", "SERPINC1", "SH3BGRL2", "SHE", "SLC1A5", "SLC29A2", "SLC2A3", "SLC39A8", "SMPDL3A", "SNX14", "SNX9", "SOCS1", "SOCS2", "SPP1", "SPRED2", "SPRY4", "ST3GAL4", "SWAP70", "SYNGR2", "SYT11", "TGM2", "TIAM1", "TLR7", "TNFRSF18", "TNFRSF1B", "TNFRSF21", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFSF10", "TNFSF11", "TRAF1", "TTC39B", "TWSG1", "UCK2", "UMPS", "WLS", "XBP1"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = il2, ctrl = 5, name = 'il2.score')
# p1 <- FeaturePlot(seurat.object2, features = 'il2.score1')
# ggsave("FeaturePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'il2.score1')
ggsave("RidgePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "il2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "il2.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

il6 <- c("A2M", "ACVR1B", "ACVRL1", "BAK1", "CBL", "CCL7", "CCR1", "CD14", "CD36", "CD38", "CD44", "CD9", "CNTFR", "CRLF2", "CSF1", "CSF2", "CSF2RA", "CSF2RB", "CSF3R", "CXCL1", "CXCL10", "CXCL11", "CXCL13", "CXCL3", "CXCL9", "DNTT", "EBI3", "FAS", "GRB2", "HAX1", "HMOX1", "IFNAR1", "IFNGR1", "IFNGR2", "IL10RB", "IL12RB1", "IL13RA1", "IL15RA", "IL17RA", "IL17RB", "IL18R1", "IL1B", "IL1R1", "IL1R2", "IL2RA", "IL2RG", "IL3RA", "IL4R", "IL6", "IL6ST", "IL7", "IL9R", "INHBE", "IRF1", "IRF9", "ITGA4", "ITGB3", "JUN", "LEPR", "LTB", "LTBR", "MAP3K8", "MYD88", "OSMR", "PDGFC", "PF4", "PIK3R5", "PIM1", "PLA2G2A", "PTPN1", "PTPN11", "PTPN2", "REG1A", "SOCS1", "SOCS3", "STAM2", "STAT1", "STAT2", "STAT3", "TGFB1", "TLR2", "TNF", "TNFRSF12A", "TNFRSF1A", "TNFRSF1B", "TNFRSF21", "TYK2"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = il6, ctrl = 5, name = 'il6.score')
# p1 <- FeaturePlot(seurat.object2, features = 'il6.score1')
# ggsave("FeaturePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'il6.score1')
ggsave("RidgePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "il6.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "il6.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

apoptosis <- c("ADD1", "AIFM3", "ANKH", "ANXA1", "APP", "ATF3", "AVPR1A", "BAX", "BCAP31", "BCL10", "BCL2L1", "BCL2L10", "BCL2L11", "BCL2L2", "BGN", "BID", "BIK", "BIRC3", "BMF", "BMP2", "BNIP3L", "BRCA1", "BTG2", "BTG3", "CASP1", "CASP2", "CASP3", "CASP4", "CASP6", "CASP7", "CASP8", "CASP9", "CAV1", "CCNA1", "CCND1", "CCND2", "CD14", "CD2", "CD38", "CD44", "CD69", "CDC25B", "CDK2", "CDKN1A", "CDKN1B", "CFLAR", "CLU", "CREBBP", "CTH", "CTNNB1", "CYLD", "DAP", "DAP3", "DCN", "DDIT3", "DFFA", "DIABLO", "DNAJA1", "DNAJC3", "DNM1L", "DPYD", "EBP", "EGR3", "EMP1", "ENO2", "ERBB2", "ERBB3", "EREG", "ETF1", "F2", "F2R", "FAS", "FASLG", "FDXR", "FEZ1", "GADD45A", "GADD45B", "GCH1", "GNA15", "GPX1", "GPX3", "GPX4", "GSN", "GSR", "GSTM1", "GUCY2D", "H1-0", "HGF", "HMGB2", "HMOX1", "HSPB1", "IER3", "IFITM3", "IFNB1", "IFNGR1", "IGF2R", "IGFBP6", "IL18", "IL1A", "IL1B", "IL6", "IRF1", "ISG20", "JUN", "KRT18", "LEF1", "LGALS3", "LMNA", "LUM", "MADD", "MCL1", "MGMT", "MMP2", "NEDD9", "NEFH", "PAK1", "PDCD4", "PDGFRB", "PEA15", "PLAT", "PLCB2", "PLPPR4", "PMAIP1", "PPP2R5B", "PPP3R1", "PPT1", "PRF1", "PSEN1", "PSEN2", "PTK2", "RARA", "RELA", "RETSAT", "RHOB", "RHOT2", "RNASEL", "ROCK1", "SAT1", "SATB1", "SC5D", "SLC20A1", "SMAD7", "SOD1", "SOD2", "SPTAN1", "SQSTM1", "TAP1", "TGFB2", "TGFBR3", "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF12A", "TNFSF10", "TOP2A", "TSPO", "TXNIP", "VDAC2", "WEE1", "XIAP"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = apoptosis, ctrl = 5, name = 'apoptosis.score')
# p1 <- FeaturePlot(seurat.object2, features = 'apoptosis.score1')
# ggsave("FeaturePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'apoptosis.score1')
ggsave("RidgePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "apoptosis.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "apoptosis.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tnf.alpha <- c("ABCA1", "ACKR3", "AREG", "ATF3", "ATP2B1", "B4GALT1", "B4GALT5", "BCL2A1", "BCL3", "BCL6", "BHLHE40", "BIRC2", "BIRC3", "BMP2", "BTG1", "BTG2", "BTG3", "CCL2", "CCL20", "CCL4", "CCL5", "CCN1", "CCND1", "CCNL1", "CCRL2", "CD44", "CD69", "CD80", "CD83", "CDKN1A", "CEBPB", "CEBPD", "CFLAR", "CLCF1", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL11", "CXCL2", "CXCL3", "CXCL6", "DDX58", "DENND5A", "DNAJB4", "DRAM1", "DUSP1", "DUSP2", "DUSP4", "DUSP5", "EDN1", "EFNA1", "EGR1", "EGR2", "EGR3", "EHD1", "EIF1", "ETS2", "F2RL1", "F3", "FJX1", "FOS", "FOSB", "FOSL1", "FOSL2", "FUT4", "G0S2", "GADD45A", "GADD45B", "GCH1", "GEM", "GFPT2", "GPR183", "HBEGF", "HES1", "ICAM1", "ICOSLG", "ID2", "IER2", "IER3", "IER5", "IFIH1", "IFIT2", "IFNGR2", "IL12B", "IL15RA", "IL18", "IL1A", "IL1B", "IL23A", "IL6", "IL6ST", "IL7R", "INHBA", "IRF1", "IRS2", "JAG1", "JUN", "JUNB", "KDM6B", "KLF10", "KLF2", "KLF4", "KLF6", "KLF9", "KYNU", "LAMB3", "LDLR", "LIF", "LITAF", "MAFF", "MAP2K3", "MAP3K8", "MARCKS", "MCL1", "MSC", "MXD1", "MYC", "NAMPT", "NFAT5", "NFE2L2", "NFIL3", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NINJ1", "NR4A1", "NR4A2", "NR4A3", "OLR1", "PANX1", "PDE4B", "PDLIM5", "PER1", "PFKFB3", "PHLDA1", "PHLDA2", "PLAU", "PLAUR", "PLEK", "PLK2", "PLPP3", "PMEPA1", "PNRC1", "PPP1R15A", "PTGER4", "PTGS2", "PTPRE", "PTX3", "RCAN1", "REL", "RELA", "RELB", "RHOB", "RIPK2", "RNF19B", "SAT1", "SDC4", "SERPINB2", "SERPINB8", "SERPINE1", "SGK1", "SIK1", "SLC16A6", "SLC2A3", "SLC2A6", "SMAD3", "SNN", "SOCS3", "SOD2", "SPHK1", "SPSB1", "SQSTM1", "STAT5A", "TANK", "TAP1", "TGIF1", "TIPARP", "TLR2", "TNC", "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFRSF9", "TNFSF9", "TNIP1", "TNIP2", "TRAF1", "TRIB1", "TRIP10", "TSC22D1", "TUBB2A", "VEGFA", "YRDC", "ZBTB10", "ZC3H12A", "ZFP36"))
seurat.object2 <- AddModuleScore(object = seurat.object2, features = tnf.alpha, ctrl = 5, name = 'tnf.alpha.score')
# p1 <- FeaturePlot(seurat.object2, features = 'tnf.alpha.score1')
# ggsave("FeaturePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object2, features = 'tnf.alpha.score1')
ggsave("RidgePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "tnf.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
# p1 <- FeatureScatter(seurat.object2, feature1 = "tnf.alpha.score1", feature2 = "percent.viral")
# ggsave("FeatureScatter_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

setwd("..")



############################################################# ############################################################# #############################################################


############################################################# ############################################################# #############################################################
############################################################# DE Analysis Btwn Treatment #############################################################
############################################################# ############################################################# #############################################################

dir.create("Treatment")
setwd("Treatment")
Idents(seurat.object2) <- "Treatment"
levels(seurat.object2)
p1 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Treatment")
ggsave("integrated_UMAP-Treatment.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object2) <- "Spatial"
seurat.object2 <- NormalizeData(seurat.object2, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object2 <- ScaleData(seurat.object2, features = all.genes.visium, assay="Spatial")

Idents(seurat.object2) <- "Treatment"
de_markers <- FindAllMarkers(seurat.object2, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, raster = FALSE)
ggsave("integrated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("integrated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("integrated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Treatment")
dev.off()
setwd("..")
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
############################################################# DE Analysis Btwn Virus #############################################################
############################################################# ############################################################# #############################################################

dir.create("Virus")
setwd("Virus")
Idents(seurat.object2) <- "Virus"
levels(seurat.object2)
p1 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Virus")
ggsave("integrated_UMAP-Virus.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object2) <- "Spatial"
seurat.object2 <- NormalizeData(seurat.object2, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object2 <- ScaleData(seurat.object2, features = all.genes.visium, assay="Spatial")

Idents(seurat.object2) <- "Virus"
de_markers <- FindAllMarkers(seurat.object2, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, raster = FALSE)
ggsave("integrated_top15_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("integrated_top15_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("integrated_top30.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Virus")
dev.off()
setwd("..")
############################################################# ############################################################# #############################################################

############################################################# ############################################################# #############################################################
############################################################# DE Analysis Btwn orig.ident #############################################################
############################################################# ############################################################# #############################################################

dir.create("orig.ident")
setwd("orig.ident")
Idents(seurat.object2) <- "orig.ident"
levels(seurat.object2)
p1 <- DimPlot(seurat.object2, reduction = "umap", group.by = "orig.ident")
ggsave("integrated_UMAP-orig.ident.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object2) <- "Spatial"
seurat.object2 <- NormalizeData(seurat.object2, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object2 <- ScaleData(seurat.object2, features = all.genes.visium, assay="Spatial")

Idents(seurat.object2) <- "orig.ident"
de_markers <- FindAllMarkers(seurat.object2, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, raster = FALSE)
ggsave("integrated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("integrated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("integrated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="SampleID")
dev.off()
setwd("..")
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################

############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################


dir.create("ZIKV-Veh.gene.list")
setwd("ZIKV-Veh.gene.list")
## ZIKV-Veh Treatment Gene List
Idents(seurat.object2) <- "Treatment"
gene.list <- c("CTLA2A", "KAP","PRAP1", "MALAT1", "MT-ND2", "MT-ATP6", "MT-CO1")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "CTLA2A", images=image.list, ncol=4)
p5
ggsave("CTLA2A_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "CTLA2A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CTLA2A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "PRAP1", images=image.list, ncol=4)
p5
ggsave("PRAP1_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "PRAP1", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("PRAP1_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "KAP", images=image.list, ncol=4)
p5
ggsave("KAP_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "KAP", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("KAP_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


p5 <- SpatialPlot(seurat.object2, features = "CTLA2A", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object2, features = "CTLA2A",)
p5
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = "CTLA2A")
p5
ggsave("CTLA2A_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")


image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "CTLA2A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CTLA2A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

Idents(seurat.object2) <- "Treatment"
feature.plot <- DotPlot(seurat.object2, features = gene.list)
png(file=paste0("ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="Treatment")
dev.off()

seurat.object2 <- ScaleData(seurat.object2, features = gene.list, assay="Spatial")


top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("ZIKV-Veh.gene.list_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("ZIKV-Veh.gene.list_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object2) <- "orig.ident"
feature.plot <- DotPlot(seurat.object2, features = gene.list)
png(file=paste0("orig.ident_ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("orig.ident_ZIKV-Veh.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("orig.ident_ZIKV-Veh.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("orig.ident_ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("orig.ident_ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object2) <- "Type.Profile"
feature.plot <- DotPlot(seurat.object2, features = gene.list)
png(file=paste0("Type.Profile_ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("Type.Profile_ZIKV-Veh.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Type.Profile_ZIKV-Veh.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Type.Profile_ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Type.Profile_ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object2) <- "Subtype.Profile"
feature.plot <- DotPlot(seurat.object2, features = gene.list)
png(file=paste0("Subtype.Profile_ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
setwd("..")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
Idents(seurat.object2) <- "Type.Profile"
plot2 <- SpatialDimPlot(seurat.object2,images=image.list, alpha=1, ncol=2, crop=F) 
plot2
ggsave("SpatialDimPlot_annotated-Type.Profile_all-uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2 <- SpatialDimPlot(seurat.object2, images=image.list, alpha=0.1, ncol=2, crop=F) 
ggsave("SpatialDimPlot_annotated-Type.Profile_all-alpha.1-uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
Idents(seurat.object2) <- "Subtype.Profile"
plot2 <- SpatialDimPlot(seurat.object2,images=image.list, alpha=1, ncol=2, crop=F) 
ggsave("SpatialDimPlot_annotated-Subtype.Profile_all-uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

plot2 <- SpatialDimPlot(seurat.object2, images=image.list, alpha=0.1, ncol=2, crop=F) 
ggsave("SpatialDimPlot_annotated-Subtype.Profile_all-alpha.1-uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
Idents(seurat.object2) <- "orig.ident"
plot2 <- SpatialDimPlot(seurat.object2,images=image.list, alpha=10, ncol=4, crop=F) + NoLegend()

p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = "FOLR2", pt.size.factor = 3, ncol = 4, crop = F)
p1
ggsave("Hofbauer-immune.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = "CCR2", pt.size.factor = 3, ncol = 4, crop = F)
p1
ggsave("PAMM-immune.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")


p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S07.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S07.slice.1")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S09.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S12.slice.1")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S12.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S13.slice.3")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S13.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S31.slice.4")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S31.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S32.slice.5")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S32.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S33.slice.6")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S33.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S33.slice.6")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Trophoblast", "DecidualStroma", "Immune"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_Tropho.Immune.Decidua_S34.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

## Try plotting all 
	## for subtypes later c("Prolif_Fibro", "SynTI", "Prolif_Stromal", "Endothelial_1", "FetalMesenchyme_1", "GC", "Erythrocyte", "SpT_1", "FetalMesenchyme_2", "Prolif_Fibro", "SynTII", "Immune_1", "S-TGC", "SpT_2", "SYT", "VCT", "FetalMesenchyme_3", "Immune_2", "Endothelial_2")
image.list <- c("S07.slice.3")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S07.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.1")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S09.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S12.slice.1")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S12.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S13.slice.3")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S13.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S31.slice.4")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S31.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S32.slice.5")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S32.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S33.slice.6")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S33.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S33.slice.6")
p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S34.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")


p1 <- SpatialFeaturePlot(seurat.object2, images = image.list, features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
p1
ggsave("Predictions_cell.types_S34.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")


p4 <- FeaturePlot(seurat.object2, reduction = "umap", features = c("Prolif_Fibro", "Trophoblast", "Prolif_Stromal", "Endothelial", "FetalMesenchyme", "Endometrial_Stromal ", "Erythrocyte", "Immune", "DecidualStroma"), raster=T) + labs(color='Predictions')
p4 		
	# The following requested variables were not found: Prolif_Fibro, Prolif_Stromal, Endometrial_Stromal
ggsave("Visium-Marsh-Type-Predictions_UMAP.pdf", plot = p4, device = "pdf", width = 10, height = 12, units = "in")

dir.create("ZIKV-Veh_complement_.gene.list")
setwd("ZIKV-Veh_complement_.gene.list")
DefaultAssay(seurat.object2) <- "Spatial"

## ZIKV-Veh Treatment Gene List
Idents(seurat.object2) <- "Treatment"
seurat.object2 <- ScaleData(seurat.object2, features = gene.list, assay="Spatial")

gene.list <- c("CTLA2A", "KAP","PRAP1", "C1R", "C1S", "C1QA", "C1QB", "C2", "C3", "C4A", "C4B", "C5", "C6", "C7", "C8", "C9")
seurat.object2 <- ScaleData(seurat.object2, features = gene.list, assay="Spatial")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("ZIKV-Veh.gene.list_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("ZIKV-Veh.gene.list_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object2) <- "orig.ident"
feature.plot <- DotPlot(seurat.object2, features = gene.list)
png(file=paste0("orig.ident_ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("orig.ident_ZIKV-Veh.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("orig.ident_ZIKV-Veh.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("orig.ident_ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("orig.ident_ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object2) <- "Type.Profile"
feature.plot <- DotPlot(seurat.object2, features = gene.list)
png(file=paste0("Type.Profile_ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("Type.Profile_ZIKV-Veh.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Type.Profile_ZIKV-Veh.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Type.Profile_ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Type.Profile_ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object2) <- "Subtype.Profile"
feature.plot <- DotPlot(seurat.object2, features = gene.list)
png(file=paste0("Subtype.Profile_ZIKV-Veh.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Subtype.Profile_ZIKV-Veh.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
setwd("..")

dir.create("miRNA.gene.list")
setwd("miRNA.gene.list")
Idents(seurat.object2) <- "Treatment"

gene.list <- c("SLC12A8", "SDK1", "VLDLR", "PSTPIP1", "DNAJB1", "IFFO2", "TNK2", "CYB561D2", "PPIF", "SPATA2L", "HERC2","HIST1H2BJ","SLC44A4","FAM189A2","PARD3B", "CYB561D2", "TMSB15B1","TMSB15B2")
gene.list <- c("SLC12A8", "LHFPL2", "DENND1A", "PAX8-AS1", "XPO4", "NR3C2", "HECTD1", "LPIN2", "FER", "ACCS", "AFF4", "VLDLR", "LMTK2", "ARHGEF11", "DOCK9", "SDK1", "DDX46")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("mirna.gene.list_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("mirna.gene.list_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("mirna.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("mirna.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object2) <- "orig.ident"
feature.plot <- DotPlot(seurat.object2, features = gene.list)
png(file=paste0("orig.ident_mirna.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("orig.ident_mirna.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("orig.ident_mirna.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("orig.ident_mirna.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("orig.ident_mirna.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object2) <- "Type.Profile"
feature.plot <- DotPlot(seurat.object2, features = gene.list)
png(file=paste0("Type.Profile_mirna.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("Type.Profile_mirna.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Type.Profile_mirna.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Type.Profile_mirna.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Type.Profile_mirna.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object2) <- "Subtype.Profile"
feature.plot <- DotPlot(seurat.object2, features = gene.list)
png(file=paste0("Subtype.Profile_mirna.gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("Subtype.Profile_mirna.gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Subtype.Profile_mirna.gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Subtype.Profile_mirna.gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Subtype.Profile_mirna.gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
setwd("..")

dir.create("immune.ZIKV-Veh.gene.list")
setwd("immune.ZIKV-Veh.gene.list")
Idents(seurat.object2) <- "Treatment"

setwd("immune")
setwd("immune.ZIKV-Veh.gene.list")
Idents(seurat.object2) <- "Treatment"
## ZIKV-Veh Treatment Gene List
Idents(seurat.object2) <- "Treatment"
gene.list <- c("CTLA2A", "KAP", "PRAP1", "PLA1A", "CCL21A",  "C3", "LYZ2", "CCL8")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")

p5 <- SpatialPlot(seurat.object2, features = "CTLA2A", images=image.list, ncol=4)
p5
ggsave("CTLA2A_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "CTLA2A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CTLA2A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "CTLA2A", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object2, features = "CTLA2A",)
p5
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = "CTLA2A")
p5
ggsave("CTLA2A_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "CTLA2A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CTLA2A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


p5 <- SpatialPlot(seurat.object2, features = "PLA1A", images=image.list, ncol=4)
p5
ggsave("PLA1A_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "PLA1A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("PLA1A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "PLA1A", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object2, features = "PLA1A",)
p5
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = "PLA1A")
p5
ggsave("PLA1A_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "PLA1A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("PLA1A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


p5 <- SpatialPlot(seurat.object2, features = "CCL21A", images=image.list, ncol=4)
p5
ggsave("CCL21A_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "CCL21A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CCL21A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "CCL21A", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object2, features = "CCL21A",)
p5
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = "CCL21A")
p5
ggsave("CCL21A_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "CCL21A", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CCL21A_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

p5 <- SpatialPlot(seurat.object2, features = "C3", images=image.list, ncol=4)
p5
ggsave("C3_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "C3", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("C3_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "C3", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object2, features = "C3",)
p5
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = "C3")
p5
ggsave("C3_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "C3", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("C3_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


p5 <- SpatialPlot(seurat.object2, features = "LYZ2", images=image.list, ncol=4)
p5
ggsave("LYZ2_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "LYZ2", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("LYZ2_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "LYZ2", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object2, features = "LYZ2",)
p5
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = "LYZ2")
p5
ggsave("LYZ2_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "LYZ2", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("LYZ2_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


p5 <- SpatialPlot(seurat.object2, features = "CCL8", images=image.list, ncol=4)
p5
ggsave("CCL8_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "CCL8", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CCL8_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "CCL8", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object2, features = "CCL8",)
p5
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = "CCL8")
p5
ggsave("CCL8_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "CCL8", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("CCL8_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


p5 <- SpatialPlot(seurat.object2, features = "PRAP1", images=image.list, ncol=4)
p5
ggsave("PRAP1_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "PRAP1", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("PRAP1_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "PRAP1", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object2, features = "PRAP1",)
p5
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = "PRAP1")
p5
ggsave("PRAP1_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "PRAP1", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("PRAP1_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

p5 <- SpatialPlot(seurat.object2, features = "KAP", images=image.list, ncol=4)
p5
ggsave("KAP_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "KAP", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("KAP_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "KAP", images=image.list, ncol=4)
p5
p5 <- FeaturePlot(seurat.object2, features = "KAP",)
p5
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = "KAP")
p5
ggsave("KAP_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")
Idents(seurat.object2) <- "Treatment"
p5 <- VlnPlot(seurat.object2, features = gene.list)
p5
ggsave("ZIKV-Veh.gene.list.list_vlnplot_treatment.pdf", plot = p5, device = "pdf", width = 8, height = 8, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "KAP", images=image.list, ncol=4, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("KAP_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")




gene.list <- c("CTLA2A", "PRAP1", "CCL21A", "KAP", "SRGN", "PENK", "IGFBP6", "HTRA3", "AEBP1", "GPX3", "TAGLN", "FGL2", "DES", "LBP", "SLPI", "CNN1", "GUCA2B", "LY6C1", "C3", "PLA1A", "CTSK", "IFI27L2A", "POSTN", "LYZ2", "SFRP4", "GSN", "TIMP2", "ACTG2", "PLTP", "IGFBP2", "DCN", "SPARCL1", "CCL8")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("gene.list_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("gene.list_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object2) <- "orig.ident"
feature.plot <- DotPlot(seurat.object2, features = gene.list)
png(file=paste0("orig.ident_gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("orig.ident_gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("orig.ident_gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("orig.ident_gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("orig.ident_gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object2) <- "Type.Profile"
feature.plot <- DotPlot(seurat.object2, features = gene.list)
png(file=paste0("Type.Profile_gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("Type.Profile_gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Type.Profile_gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Type.Profile_gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Type.Profile_gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object2) <- "Subtype.Profile"
feature.plot <- DotPlot(seurat.object2, features = gene.list)
png(file=paste0("Subtype.Profile_gene.list-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="RNA Expression") + scale_y_discrete(name ="SampleID")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, raster = FALSE, size=5) 
ggsave("Subtype.Profile_gene.list_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = gene.list, slot="counts", raster = FALSE, size=5) 
ggsave("Subtype.Profile_gene.list_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 11, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=gene.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("Subtype.Profile_gene.list_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial", size=5) 
ggsave("Subtype.Profile_gene.list_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
setwd("..")

rm(seurat.object2)

setwd("..")
############################################################# ############################################################# #############################################################


############################################################# 		#############################################################
############################################################# ############################################################# #############################################################

setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1")
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
save.image("murine_spatial_data-integrated-annotated_v1.1.RData")
load("murine_spatial_data-integrated-annotated_v1.1.RData")


# seurat.object <- spatial.seurat.object

## use this plot below to determine names of images retained after integration. 
plot2 <- SpatialFeaturePlot(seurat.object, features = "DDX3Y", ncol = 1) + patchwork::plot_layout(ncol = 1)
plot2
####### FOR SOME REASON, INTEGRATION ADDS IMAGES ARBITRARILY
## ggsave("integrated_SpatialDimPlot_percent.viral_TEST.pdf", plot = plot2, device = "pdf", width = 45, height = 11, units = "in")


## 8.11.22 determine fetal sex
Idents(seurat.object) <- "orig.ident"
levels(seurat.object)	# [1] "Mock-Enox-1" "Mock-Enox-2" "ZIKV-Veh-1"   "ZIKV-Enox-1"

library.averages.heatmap <- VlnPlot(seurat.object, features = "SRY", assay="Spatial")
library.averages.heatmap
ggsave("vlnplot_SRY_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object) <- "orig.ident"
library.averages.heatmap <- VlnPlot(seurat.object, features = "XIST", assay="Spatial")
library.averages.heatmap
ggsave("vlnplot_XIST_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object) <- "orig.ident"
library.averages.heatmap <- VlnPlot(seurat.object, features = "HBB", assay="Spatial")
library.averages.heatmap
ggsave("vlnplot_HBB_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object) <- "orig.ident"
library.averages.heatmap <- VlnPlot(seurat.object, features = "DDX3X", assay="Spatial")
library.averages.heatmap
ggsave("vlnplot_DDX3X_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object) <- "orig.ident"
library.averages.heatmap <- VlnPlot(seurat.object, features = "DDX3Y", assay="Spatial")
library.averages.heatmap	### [1] "Mock-Enox-1" + "Mock-Enox-2" + "ZIKV-Veh-1"  -  "ZIKV-Enox-1" +
ggsave("vlnplot_DDX3Y_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")


image.list <- c("S07.slice.1")
p2 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p2
image.list <- c("S07.slice.1.2")
p3 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p3
image.list <- c("S12.slice.1")
p4 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p4
image.list <- c("S13.slice.3")
p5 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
image.list <- c("S07.slice.1")
plot2<- SpatialFeaturePlot(seurat.object, features = "DDX3Y", images=image.list, slot = "counts") + patchwork::plot_layout(ncol = 1)
plot2
ggsave("S07_SpatialFeaturePlot_percent.DDX3Y.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")

image.list <- c("S07.slice.3")
plot2<- SpatialFeaturePlot(seurat.object, features = "DDX3Y", images=image.list, slot = "counts") + patchwork::plot_layout(ncol = 1)
plot2
ggsave("S09_SpatialFeaturePlot_percent.DDX3Y.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")

image.list <- c("S12.slice.1")
plot2<- SpatialFeaturePlot(seurat.object, features = "DDX3Y", images=image.list, slot = "counts") + patchwork::plot_layout(ncol = 1)
plot2
ggsave("S12_SpatialFeaturePlot_percent.DDX3Y.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")

image.list <- c("S13.slice.3")
plot2<- SpatialFeaturePlot(seurat.object, features = "DDX3Y", images=image.list, slot = "counts") + patchwork::plot_layout(ncol = 1)
plot2
ggsave("S13_SpatialFeaturePlot_percent.DDX3Y.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")






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
Idents(seurat.object) <- "Subtype.Profile"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "Subtype.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(seurat.object) <- "Type.Profile"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "Type.cellcounts.txt", sep="\t")



dir.create("Type.cluster")
setwd("Type.cluster")

Idents(seurat.object) <- "orig.ident"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "sample.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(seurat.object) <- "seurat_clusters"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "cluster.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(seurat.object) <- "Subtype.Profile"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "Subtype.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(seurat.object) <- "Type.Profile"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "Type.cellcounts.txt", sep="\t")

# "Trophoblast", "Fibroblast", "Stromal", "Smooth_muscle", "Endothelial", "Endometrial"
# c("Endometrial", "Endothelial", "Fibroblast", "Stromal", "VCT", "SYT", "NKT, "Macrophage", "Erythroblast")

Idents(seurat.object)<-"Type.Profile"
Idents(seurat.object)
clusters.0.uninfected <- subset(seurat.object, idents = c("Trophoblast"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Trophoblast.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type.Profile"
clusters.0.uninfected <- subset(seurat.object, idents = c("Stromal"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Stromal.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type.Profile"
clusters.0.uninfected <- subset(seurat.object, idents = c("Fibroblast"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Fibroblast.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type.Profile"
clusters.0.uninfected <- subset(seurat.object, idents = c("Endothelial"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Endothelial.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type.Profile"
clusters.0.uninfected <- subset(seurat.object, idents = c("Smooth_muscle"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Smooth_muscle.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type.Profile"
clusters.0.uninfected <- subset(seurat.object, idents = c("Endometrial"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Endometrial.cluster.cellcounts.txt", sep="\t")
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
df<-read.table("Type.cluster/Endometrial.cluster.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("cluster","Endometrial")
 df2<-df
df<-df%>%pivot_longer(cols = `Endometrial`,names_to = "Type.Profile",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("Type.cluster/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("cluster",X)
  df<-df%>%pivot_longer(cols = X,names_to = "Type.Profile",values_to = "Count")
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

final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)
final<-full_join(final,df7)


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
          x = "Type.Profile",
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
          x = "Type.Profile",
          y = "Percent",
          color = "cluster",
          fill = "cluster",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()

####################################################################################################




dir.create("seurat_clusters.Type")
setwd("seurat_clusters.Type")

# "Trophoblast", "Fibroblast", "Stromal", "Smooth_muscle", "Endothelial", "Endometrial"
# c("Endometrial", "Endothelial", "Fibroblast", "Stromal", "VCT", "SYT", "NKT, "Macrophage", "Erythroblast")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("1"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "1.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("2"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "2.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("3"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "3.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("4"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "4.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("5"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "5.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("6"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "6.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("7"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "7.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("8"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "8.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("9"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "9.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("10"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "10.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("11"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "11.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("12"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "12.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("13"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "13.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("14"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "14.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("15"))
Idents(Types.0.uninfected) <- "Type.Profile"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "15.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")


setwd("..")


