# Aagaard_murineplacenta_spatialtranscriptomics
Gnotobiotic congenital ZIKV mouse model. Following timed mattings of four-to-six-week-old Swiss/webster mice from the Baylor College of Medicine Gnotobiotics Core, plugged dams were weighed and transferred to separate sterile insulators. Daily intraperitoneal injections of enoxacin (Sigma-Aldrich, Cat. 557304) at 10 mg/kg were previously shown to be safe and effective in a model of obesity32 and were done from embryonic days E1.5-E18.5. Dams were either received mock United States Pharmacopeia (USP)-grade phosphate-buffered saline (PBS) (VWR, Cat. 76081-514) subdermal injections or ZIKV infection with HN16 strain (first passage in Vero-E6 cells) at an MOI of 1x105 PFU on embryonic days E7.5, E8.5, and E9.5. On day 18.5, dams were humanely sacrificed, and uterine horns were dissected. Images of the uterine horns were taken to record resorption events. Individual fetuses were weighed. Placental tissues were stored at -80C in Trizol. All mouse experiments were conducted under the supervision of Baylor College of Medicine Institutional Animal Care and Use Committee approved protocols.

Spatial transcriptomics in murine placentae. At E18.5 placentas were dissected and preserved in RNA-later at -25C for 18 hours. Tissues were dehydrated by incubation in 30% sucrose at 4oC for 6 hours. At the Baylor College of Medicine Pathology and Histology Core placentas were fresh-frozen in optimal cutting temperature compound (OCT) with a liquid nitrogen and isopentane bath. Tissue blocks were cryosectioned at 10 microns for quality control H&E staining and 10x Genomics Visium Tissue Optimization. The Tissue Optimization slide revealed optimal permeabilization was observed after 12 minutes. Finally, the tissue blocks were cryosectioned and placed on Visium Gene Expression slides (v2) for H&E staining and Visium spatial transcriptomics library preparation following the manufacturer’s instructions. Visium libraries were sequenced (Illumina, NovaSeq 6000) at the Baylor College of Medicine Genomic and RNA Profiling Core.

Spatial transcriptomics analysis. Spatial transcriptomes were aligned to a custom human and ZIKV transcriptome ( and transcripts counts were quantified using Space Ranger (v1.3.0) and bash scripts. Downstream analyses were done in R using the Seurat scRNA-seq analysis package 34. H&E slides were annotated by a trained pathologist and transferred to spatial barcodes using the Loupe Browser (v5) and spatial transcriptomes were annotated based on PCA dimension reduction and clustering analysis to identify broad cell-types. The spatial transcriptomes were then integrated with published murine single-nuclei RNA-seq (snRNA-seq) dataset to assess cell-type annotations (GEO accession GSE152248). The Marsh & Blelloch, eLife. (2020) snRNA-seq dataset included murine placentas sampled at E9.5-E14.5: E9.5, E10.5, E12.5, E14.5 across 2,933 independently filtered (nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 10000 & percent.mt < 10 & percent.ribo < 20) high-quality transcriptomes.

S07 ZIKV-Vehicle
S09 Mock
S12 Mock
S13 ZIKV-Enoxacin

# Repository contents
A) Images necessary for alignment of spatial transcriptomes

B) Bash scripts for creating a custom ZIKV+Mouse transcriptome and alignment using Spaceranger

C) R script for QC filtering and initial clustering of spatial transcriptomics data using Seurat

D) R script for clustering of previously published placenta snRNA-seq data (Marsh et al) using Seurat

E) R script for integration of placenta spatial transcriptomics and snRNA-seq data using Seurat
