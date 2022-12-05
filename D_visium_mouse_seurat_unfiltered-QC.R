## visium_mouse_seurat_unfiltered-QC.R
# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine



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


############################################################# ############################################################# #############################################################
############################################################# S31: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S31-manual/outs'
list.files(data_dir)

#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S31.manual <- Load10X_Spatial(data_dir, filename = "S31_filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S31.slice", filter.matrix = TRUE, to.upper = TRUE)
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
S32.manual <- Load10X_Spatial(data_dir, filename = "S32_filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S32.slice", filter.matrix = TRUE, to.upper = TRUE)
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
S33.manual <- Load10X_Spatial(data_dir, filename = "S33_filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S33.slice", filter.matrix = TRUE, to.upper = TRUE)
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
S34.manual <- Load10X_Spatial(data_dir, filename = "S34_filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S34.slice", filter.matrix = TRUE, to.upper = TRUE)
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


rm(plot1)
rm(plot2)
rm(acute.viral)
rm(plot3)
rm(qc.vlnplot)
rm(new.names)
rm(old.names)
rm(data_dir)
rm(to.upper)
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
ggsave("all_merged_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 15, height = 5, units = "in")


## plot UMI counts/spot 
plot1 <- SpatialFeaturePlot(S34.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S34.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
plot3<-wrap_plots(plot1, plot2)
plot11<-wrap_plots(plot1, plot2)
ggsave("S34.manual_unfiltered_nCount-nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")
## plot UMI counts/spot 
plot1 <- SpatialFeaturePlot(S33.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S33.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
plot3<-wrap_plots(plot1, plot2)
plot4<-wrap_plots(plot1, plot2)
ggsave("S33.manual_unfiltered_nCount-nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")
## plot UMI counts/spot 
plot1 <- SpatialFeaturePlot(S32.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S32.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
plot3<-wrap_plots(plot1, plot2)
plot5<-wrap_plots(plot1, plot2)
ggsave("S32.manual_unfiltered_nCount-nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")
## plot UMI counts/spot 
plot1 <- SpatialFeaturePlot(S31.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S31.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
plot3<-wrap_plots(plot1, plot2)
plot6<-wrap_plots(plot1, plot2)
ggsave("S31.manual_unfiltered_nCount-nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")
## plot UMI counts/spot 
plot1 <- SpatialFeaturePlot(S07.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S07.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
plot3<-wrap_plots(plot1, plot2)
plot7<-wrap_plots(plot1, plot2)
ggsave("S07.manual_unfiltered_nCount-nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")
## plot UMI counts/spot 
plot1 <- SpatialFeaturePlot(S09.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S09.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
plot3<-wrap_plots(plot1, plot2)
plot8<-wrap_plots(plot1, plot2)
ggsave("S09.manual_unfiltered_nCount-nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")
## plot UMI counts/spot 
plot1 <- SpatialFeaturePlot(S12.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S12.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
plot3<-wrap_plots(plot1, plot2)
plot9<-wrap_plots(plot1, plot2)
ggsave("S12.manual_unfiltered_nCount-nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")
## plot UMI counts/spot 
plot1 <- SpatialFeaturePlot(S13.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S13.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
plot3<-wrap_plots(plot1, plot2)
plot10<-wrap_plots(plot1, plot2)
ggsave("S13.manual_unfiltered_nCount-nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot12<-wrap_plots(plot7, plot8, plot9, plot10, plot6, plot5, plot4, plot11) + ncol=4



