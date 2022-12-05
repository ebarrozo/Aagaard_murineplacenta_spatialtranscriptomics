# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine
## Pearson_mouse-visium_tgfb_v1.R

Idents(seurat.object) <- "Treatment"
## use spatial for all DE analysis and plots
DefaultAssay(seurat.object) <- "Spatial"
# seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
seurat.object <- ScaleData(seurat.object, features = all.genes.visium, assay="Spatial")

setwd("/home/ebarrozo/visium/results/seurat_mouse_v2.1/annotated")
# http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_TGF_BETA_SIGNALING.html
tgfb <- c("ACVR1", "APC", "ARID4B", "BCAR3", "BMP2", "BMPR1A", "BMPR2", "CDH1", "CDK9", "CDKN1C", "CTNNB1", "ENG", "FKBP1A", "FNTA", "FURIN", "HDAC1", "HIPK2", "ID1", "ID2", "ID3", "IFNGR2", "JUNB", "KLF10", "LEFTY2", "LTBP2", "MAP3K7", "NCOR2", "NOG", "PMEPA1", "PPM1A", "PPP1CA", "PPP1R15A", "RAB31", "RHOA", "SERPINE1", "SKI", "SKIL", "SLC20A1", "SMAD1", "SMAD3", "SMAD6", "SMAD7", "SMURF1", "SMURF2", "SPTBN1", "TGFB1", "TGFBR1", "TGIF1", "THBS1", "TJP1", "TRIM33", "UBE2D3", "WWTR1", "XIAP")
# http://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_TGF_BETA_SIGNALING.html
tgfb.mouse  <-c("PPP1CA","NOG","CDKN1C","FNTA","SKIL","XIAP","IFNGR2","HDAC1","SLC20A1","SMAD1","BMPR2","RHOA","SMAD7","KLF10","TGIF1","SMAD3","HIPK2","CDK9","SMAD6","NCOR2","BMPR1A","MAP3K7","BCAR3","UBE2D3","SMURF2","RAB31","WWTR1","SMURF1","PPP1R15A","PMEPA1","TRIM33","ARID4B","ACVR1","APC","BMP2","CTNNB1","CDH1","ENG","FKBP1A","ID1","ID2","ID3","JUNB","FURIN","SERPINE1","SKI","SPTBN1","TGFB1","TGFBR1","THBS1","TJP1","LTBP2","PPM1A")
tgfb <- intersect(tgfb.mouse, all.genes.visium)
    # 53 genes
####################################################################################################
######################### Perform Pearson's Correlation Analysis of tgfb genes with all spatial transcriptomes   ########################################
########################################################################################################################
dir.create("Pearson_all_tgfb.mouse")
setwd("Pearson_all_tgfb.mouse")
Idents(seurat.object) <- "Treatment"
seurat.object <- AddModuleScore(object = seurat.object, features = tgfb.mouse, name = 'tgfb.mouse.scores')
p1<-RidgePlot(seurat.object, features = 'tgfb.mouse.scores1')
ggsave("RidgePlot_tgfb_Treatment.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "tgfb.mouse.scores1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tgfb_Treatment.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = tgfb.mouse, raster = FALSE)
ggsave("tgfb.mouse_heatmap_logfc_Treatment.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = tgfb.mouse, slot="counts", raster = FALSE)
ggsave("tgfb.mouse_heatmap_counts_Treatment.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=tgfb.mouse)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = tgfb.mouse, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("tgfb.mouse_mc_heatmap_Treatment.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = tgfb.mouse, raster = FALSE, assay="Spatial", size=5) 
ggsave("tgfb.mouse_FC_heatmap_Treatment.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

Idents(seurat.object) <- "Type.Profile"
top5.heatmap <- DoHeatmap(seurat.object, features = tgfb.mouse, raster = FALSE)
ggsave("tgfb.mouse_heatmap_logfc_Type.Profile.pdf", plot = top5.heatmap, device = "pdf", width = 17, height = 7, units = "in")

Idents(seurat.object) <- "Subtype.Profile"
top5.heatmap <- DoHeatmap(seurat.object, features = tgfb.mouse, raster = FALSE, size=5) + NoLegend()
ggsave("tgfb.mouse_heatmap_logfc_Subtype.Profile.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 7, units = "in")

Idents(seurat.object) <- "Treatment"


image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object, features = "tgfb.mouse.scores1", images=image.list, ncol=4)
p5
ggsave("tgfb.score-SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object, features = "TGFBR1", images=image.list, ncol=4)
p5
ggsave("TGFBR1-SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object, features = "TGFB1", images=image.list, ncol=4)
p5
ggsave("TGFB1-SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object, features = "SMAD2", images=image.list, ncol=4)
p5
ggsave("SMAD2-SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

################################################################################################
library(Seurat)
library(corrplot)
library(Hmisc)
library(ggcorrplot)
library(data.table)

## What identity
Idents(seurat.object) <- "Treatment"


## What identity is being subset?
    ## nothing here
    seurat.object3 <- seurat.object
# seurat.object3 <- subset(seurat.object, idents = "ZIKV-Veh")

ncol(seurat.object3)	# 9743

## What genes are were using for correlations?
features.input <- tgfb
features <- features.input


# set default assay from integrated back to Spatial. 
    ## Pearson's script originally written by Michael Mariani, Ph.D.
DefaultAssay(seurat.object3) <- "Spatial"
rows<-length(features)
cols<-length(features)    
x <- rep(NA, rows*cols)
cor.mat <- matrix(x, nrow=rows, ncol=cols)
scatter.list <- list()
count <- 1
for(i in 1:length(features))
{
  for(j in 1:length(features))     
  {
    scatter.list[[count]] <- Seurat::FeatureScatter(seurat.object3,
                                            as.character(features[i]),
                                            as.character(features[j]))
    cor.mat[i,j] <- scatter.list[[count]]$labels$title
    count <- count +1 
    print(paste0(features[i],",",features[j]))
  }
}

rownames(cor.mat) <- features
colnames(cor.mat) <- features


png(file=paste0("all_pearsons_heatmaps-tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         tl.cex=1)
dev.off()

png(file=paste0("all_pearsons_heatmaps-original_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 90,
         tl.cex=1)
dev.off()

hmisc.cor.out <- Hmisc::rcorr(cor.mat)
hmisc.cor.out$P
write.csv(hmisc.cor.out$P, file = "all_tgfb_pearson.csv")
##Compute a correlation matrix
corr <- round(cor(cor.mat), 1)
##Compute a matrix of correlation p-values
p.mat <- cor_pmat(cor.mat)
##Correlation matrix visualization
##Visualize the correlation matrix
##--------------------------------
##method = "square" (default)
png(file=paste0("all_pearsons_heatmaps_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         tl.cex=1)
dev.off()
##Correlation matrix visualization
png(file=paste0("all_pearsons_heatmaps_sig_origorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         tl.cex=1)

dev.off()

png(file=paste0("all_pearsons_upper_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("circle"),
         type = "lower", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 45,
         p.mat=p.mat,
         sig.level = 0.05,
         #insig="label_sig",
         insig="blank",
         addrect = 4,
         tl.cex=1)
dev.off()
## any significant? 
png(file=paste0("immune-ZIKV-Veh_pearsons_shade_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("shade"),
         type = "lower", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         #addrect = 12,
         tl.cex=1)
dev.off()
setwd("..")
####################################################################################################
######################### Perform Pearson's Correlation Analysis of tgfb genes with immune spatial transcriptomes   ########################################
########################################################################################################################
dir.create("Pearson_immune_tgfb")
setwd("Pearson_immune_tgfb")
Idents(seurat.object2) <- "Treatment"
seurat.object2 <- AddModuleScore(object = seurat.object2, features = tgfb.mouse, name = 'tgfb.mouse.scores')
p1<-RidgePlot(seurat.object2, features = 'tgfb.mouse.scores1')
ggsave("RidgePlot_tgfb_Treatment.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object2, features = "tgfb.mouse.scores1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tgfb_Treatment.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = tgfb.mouse, raster = FALSE)
ggsave("tgfb.mouse_heatmap_logfc_Treatment.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = tgfb.mouse, slot="counts", raster = FALSE)
ggsave("tgfb.mouse_heatmap_counts_Treatment.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=tgfb.mouse)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = tgfb.mouse, raster = FALSE, slot="counts", assay="Spatial", size=5)
ggsave("tgfb.mouse_mc_heatmap_Treatment.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = tgfb.mouse, raster = FALSE, assay="Spatial", size=5) 
ggsave("tgfb.mouse_FC_heatmap_Treatment.png", plot = library.averages.heatmap, device = "png", width = 8, height = 5, units = "in")

image.list <- c("S07.slice.3", "S07.slice.1", "S12.slice.1", "S13.slice.3", "S31.slice.4", "S32.slice.5", "S33.slice.6", "S34.slice.7")
p5 <- SpatialPlot(seurat.object2, features = "tgfb.mouse.scores1", images=image.list, ncol=4)
ggsave("tgfb.score-SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "TGFBR1", images=image.list, ncol=4)
ggsave("TGFBR1-SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "TGFB1", images=image.list, ncol=4)
ggsave("TGFB1-SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "SMAD2", images=image.list, ncol=4)
ggsave("SMAD2-SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "JUNB", images=image.list, ncol=4)
ggsave("JUNB-SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "RAR31", images=image.list, ncol=4)
ggsave("RAR31-SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "IFNGR2", images=image.list, ncol=4)
ggsave("IFNGR2-SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "SMAD3", images=image.list, ncol=4)
ggsave("SMAD3-SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
p5 <- SpatialPlot(seurat.object2, features = "NOG", images=image.list, ncol=4)
ggsave("NOG-SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

################################################################################################

################################################################################################s
################################################################################################
library(Seurat)
library(corrplot)
library(Hmisc)
library(ggcorrplot)
library(data.table)

## What identity
Idents(seurat.object2) <- "Treatment"


## What identity is being subset?
    ## nothing here
    seurat.object3 <- seurat.object2
# seurat.object3 <- subset(seurat.object2, idents = "ZIKV-Veh")

ncol(seurat.object3)    # 37

## What genes are were using for correlations?
features.input <- tgfb
features <- features.input


# set default assay from integrated back to Spatial. 
DefaultAssay(seurat.object3) <- "Spatial"
rows<-length(features)
cols<-length(features)    
x <- rep(NA, rows*cols)
cor.mat <- matrix(x, nrow=rows, ncol=cols)
scatter.list <- list()
count <- 1
for(i in 1:length(features))
{
  for(j in 1:length(features))     
  {
    scatter.list[[count]] <- Seurat::FeatureScatter(seurat.object3,
                                            as.character(features[i]),
                                            as.character(features[j]))
    cor.mat[i,j] <- scatter.list[[count]]$labels$title
    count <- count +1 
    print(paste0(features[i],",",features[j]))
  }
}

rownames(cor.mat) <- features
colnames(cor.mat) <- features
png(file=paste0("immune_pearsons_heatmaps-tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         tl.cex=1)
dev.off()

png(file=paste0("immune_pearsons_heatmaps-original_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 90,
         tl.cex=1)
dev.off()

hmisc.cor.out <- Hmisc::rcorr(cor.mat)
hmisc.cor.out$P
write.csv(hmisc.cor.out$P, file = "immune_tgfb_pearson.csv")
##Compute a correlation matrix
corr <- round(cor(cor.mat), 1)
##Compute a matrix of correlation p-values
p.mat <- cor_pmat(cor.mat)
##Correlation matrix visualization
##Visualize the correlation matrix
##--------------------------------
##method = "square" (default)
png(file=paste0("immune_pearsons_heatmaps_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         tl.cex=1)
dev.off()
##Correlation matrix visualization
png(file=paste0("immune_pearsons_heatmaps_sig_origorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         tl.cex=1)

dev.off()

png(file=paste0("all_pearsons_upper_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("circle"),
         type = "lower", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 45,
         p.mat=p.mat,
         sig.level = 0.05,
         #insig="label_sig",
         insig="blank",
         addrect = 4,
         tl.cex=1)
dev.off()
## any significant? 
png(file=paste0("immune-ZIKV-Veh_pearsons_shade_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("shade"),
         type = "lower", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         #addrect = 12,
         tl.cex=1)
dev.off()
setwd("..")

################################################################################################################################################################################################################################################################################################
####################################################################################################
######################### Perform Pearson's Correlation Analysis of tgfb genes with immune ZIKV-Veh spatial transcriptomes   ########################################
########################################################################################################################
dir.create("Pearson_immune-ZIKV-Veh_tgfb")
setwd("Pearson_immune-ZIKV-Veh_tgfb")
################################################################################################
library(Seurat)
library(corrplot)
library(Hmisc)
library(ggcorrplot)
library(data.table)

## What identity
Idents(seurat.object2) <- "Treatment"


## What identity is being subset?
    ## nothing here
    # seurat.object3 <- seurat.object2
 seurat.object3 <- subset(seurat.object2, idents = "ZIKV-Veh")

ncol(seurat.object3)    # 37

## What genes are were using for correlations?
features.input <- tgfb
features <- features.input


# set default assay from integrated back to Spatial. 
DefaultAssay(seurat.object3) <- "Spatial"
rows<-length(features)
cols<-length(features)    
x <- rep(NA, rows*cols)
cor.mat <- matrix(x, nrow=rows, ncol=cols)
scatter.list <- list()
count <- 1
for(i in 1:length(features))
{
  for(j in 1:length(features))     
  {
    scatter.list[[count]] <- Seurat::FeatureScatter(seurat.object3,
                                            as.character(features[i]),
                                            as.character(features[j]))
    cor.mat[i,j] <- scatter.list[[count]]$labels$title
    count <- count +1 
    print(paste0(features[i],",",features[j]))
  }
}

rownames(cor.mat) <- features
colnames(cor.mat) <- features
png(file=paste0("immune-ZIKV-Veh_pearsons_heatmaps-tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         tl.cex=1)
dev.off()

png(file=paste0("immune-ZIKV-Veh_pearsons_heatmaps-original_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 90,
         tl.cex=1)
dev.off()

hmisc.cor.out <- Hmisc::rcorr(cor.mat)
hmisc.cor.out$P
write.csv(hmisc.cor.out$P, file = "immune-ZIKV-Veh_tgfb_pearson.csv")
##Compute a correlation matrix
corr <- round(cor(cor.mat), 1)
##Compute a matrix of correlation p-values
p.mat <- cor_pmat(cor.mat)
##Correlation matrix visualization
##Visualize the correlation matrix
##--------------------------------
##method = "square" (default)
png(file=paste0("immune-ZIKV-Veh_pearsons_heatmaps_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         #addrect = 12,
         tl.cex=1)
dev.off()
png(file=paste0("immune-ZIKV-Veh_pearsons_shade_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("shade"),
         type = "lower", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         #addrect = 12,
         tl.cex=1)
dev.off()

##Correlation matrix visualization
png(file=paste0("immune-ZIKV-Veh_pearsons_heatmaps_sig_origorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         tl.cex=1)

dev.off()

png(file=paste0("all_pearsons_upper_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("circle"),
         type = "lower", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 45,
         p.mat=p.mat,
         sig.level = 0.05,
         #insig="label_sig",
         insig="blank",
         addrect = 4,
         tl.cex=1)
dev.off()
## any significant? 
png(file=paste0("immune-ZIKV-Veh_pearsons_shade_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("shade"),
         type = "lower", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         #addrect = 12,
         tl.cex=1)
dev.off()
setwd("..")

####################################################################################################
######################### Perform Pearson's Correlation Analysis of tgfb genes with all ZIKV-Veh spatial transcriptomes   ########################################
########################################################################################################################
dir.create("Pearson_all-ZIKV-Veh_tgfb")
setwd("Pearson_all-ZIKV-Veh_tgfb")
################################################################################################
library(Seurat)
library(corrplot)
library(Hmisc)
library(ggcorrplot)
library(data.table)

## What identity
Idents(seurat.object) <- "Treatment"


## What identity is being subset?
    ## nothing here
    # seurat.object3 <- seurat.object
 seurat.object3 <- subset(seurat.object, idents = "ZIKV-Veh")

ncol(seurat.object3)    # 37

## What genes are were using for correlations?
features.input <- tgfb
features <- features.input


# set default assay from integrated back to Spatial. 
DefaultAssay(seurat.object3) <- "Spatial"
rows<-length(features)
cols<-length(features)    
x <- rep(NA, rows*cols)
cor.mat <- matrix(x, nrow=rows, ncol=cols)
scatter.list <- list()
count <- 1
for(i in 1:length(features))
{
  for(j in 1:length(features))     
  {
    scatter.list[[count]] <- Seurat::FeatureScatter(seurat.object3,
                                            as.character(features[i]),
                                            as.character(features[j]))
    cor.mat[i,j] <- scatter.list[[count]]$labels$title
    count <- count +1 
    print(paste0(features[i],",",features[j]))
  }
}

rownames(cor.mat) <- features
colnames(cor.mat) <- features
png(file=paste0("all-ZIKV-Veh_pearsons_heatmaps-tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         tl.cex=1)
dev.off()

png(file=paste0("all-ZIKV-Veh_pearsons_heatmaps-original_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 90,
         tl.cex=1)
dev.off()

hmisc.cor.out <- Hmisc::rcorr(cor.mat)
hmisc.cor.out$P
write.csv(hmisc.cor.out$P, file = "all-ZIKV-Veh_tgfb_pearson.csv")
##Compute a correlation matrix
corr <- round(cor(cor.mat), 1)
##Compute a matrix of correlation p-values
p.mat <- cor_pmat(cor.mat)
##Correlation matrix visualization
##Visualize the correlation matrix
##--------------------------------
##method = "square" (default)
png(file=paste0("all-ZIKV-Veh_pearsons_heatmaps_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         tl.cex=1)
dev.off()
##Correlation matrix visualization
png(file=paste0("all-ZIKV-Veh_pearsons_heatmaps_sig_origorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         tl.cex=1)

dev.off()

png(file=paste0("all_pearsons_upper_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("circle"),
         type = "lower", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 45,
         p.mat=p.mat,
         sig.level = 0.05,
         #insig="label_sig",
         insig="blank",
         addrect = 4,
         tl.cex=1)
dev.off()
## any significant? 
png(file=paste0("immune-ZIKV-Veh_pearsons_shade_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("shade"),
         type = "lower", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         #addrect = 12,
         tl.cex=1)
dev.off()
################################################################################################
################################################################################################
setwd("..")



tgfb.sig <- c("ID1", "ID2", "ID3", "FKBP1A", "JUNB", "XIAP", "LTBP2", "PMEPA1", "SMURF2", "SKI", "BMPR2", "BMPR1A", "WWTR1", "SMAD7", "FURIN", "SERPINE1", "SPTBN1", "CDKN1C", "CDH1", "SLC20A1", "BCAR3")
######################### Perform Pearson's Correlation Analysis of tgfb genes with immune ZIKV-Veh spatial transcriptomes   ########################################
########################################################################################################################
dir.create("Pearson_immune-ZIKV-Veh_tgfb-sig")
setwd("Pearson_immune-ZIKV-Veh_tgfb-sig")
################################################################################################
library(Seurat)
library(corrplot)
library(Hmisc)
library(ggcorrplot)
library(data.table)

## What identity
Idents(seurat.object2) <- "Treatment"


## What identity is being subset?
    ## nothing here
    # seurat.object3 <- seurat.object2
 seurat.object3 <- subset(seurat.object2, idents = "ZIKV-Veh")

ncol(seurat.object3)    # 37

## What genes are were using for correlations?
features.input <- tgfb.sig
features <- features.input


# set default assay from integrated back to Spatial. 
DefaultAssay(seurat.object3) <- "Spatial"
rows<-length(features)
cols<-length(features)    
x <- rep(NA, rows*cols)
cor.mat <- matrix(x, nrow=rows, ncol=cols)
scatter.list <- list()
count <- 1
for(i in 1:length(features))
{
  for(j in 1:length(features))     
  {
    scatter.list[[count]] <- Seurat::FeatureScatter(seurat.object3,
                                            as.character(features[i]),
                                            as.character(features[j]))
    cor.mat[i,j] <- scatter.list[[count]]$labels$title
    count <- count +1 
    print(paste0(features[i],",",features[j]))
  }
}

rownames(cor.mat) <- features
colnames(cor.mat) <- features
png(file=paste0("immune-ZIKV-Veh_pearsons_heatmaps-tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         tl.cex=1)
dev.off()

png(file=paste0("immune-ZIKV-Veh_pearsons_heatmaps-original_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 90,
         tl.cex=1)
dev.off()

hmisc.cor.out <- Hmisc::rcorr(cor.mat)
hmisc.cor.out$P
write.csv(hmisc.cor.out$P, file = "immune-ZIKV-Veh_tgfb_pearson.csv")
##Compute a correlation matrix
corr <- round(cor(cor.mat), 1)
##Compute a matrix of correlation p-values
p.mat <- cor_pmat(cor.mat)
##Correlation matrix visualization
##Visualize the correlation matrix
##--------------------------------
##method = "square" (default)
png(file=paste0("immune-ZIKV-Veh_pearsons_heatmaps_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         #addrect = 12,
         tl.cex=1)
dev.off()
png(file=paste0("immune-ZIKV-Veh_pearsons_shade_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("shade"),
         type = "lower", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         #addrect = 12,
         tl.cex=1)
dev.off()

##Correlation matrix visualization
png(file=paste0("immune-ZIKV-Veh_pearsons_heatmaps_sig_origorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         tl.cex=1)

dev.off()

png(file=paste0("all_pearsons_upper_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("circle"),
         type = "lower", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 45,
         p.mat=p.mat,
         sig.level = 0.05,
         #insig="label_sig",
         insig="blank",
         addrect = 4,
         tl.cex=1)
dev.off()
## any significant? 
png(file=paste0("immune-ZIKV-Veh_pearsons_shade_sig_hclustorder_tgfb.png"),
                res=300, 
                width=2500, 
                height=2500)
corrplot(cor.mat, 
         method=c("shade"),
         type = "lower", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         #addrect = 12,
         tl.cex=1)
dev.off()
setwd("..")