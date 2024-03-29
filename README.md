# Aagaard_murineplacenta_spatialtranscriptomics
Gnotobiotic congenital ZIKV mouse model. Following timed mattings of Swiss/webster mice (>8-weeks old) drawn from colonies of the Baylor College of Medicine Gnotobiotics Core, plugged dams were weighed and transferred to separate sterile insulators. Daily intraperitoneal injections of enoxacin (Sigma-Aldrich, 557304) at 10 mg/kg were previously shown to be safe and effective in a model of obesity(31) and were done from embryonic days E1.5-E18.5. Dams either received mock United States Pharmacopeia (USP)-grade PBS (VWR, 76081-514) subdermal injections or ZIKV infection with HN16 strain (first passage in Vero-E6 cells) at a dose of 1x105 PFU on embryonic days E7.5, E8.5, and E9.5. On day 18.5, dams were humanely sacrificed, and uterine horns were dissected. Images of the uterine horns were taken to record resorption events. Individual fetuses were weighed. Placental tissues were stored at -80°C in Trizol. All mouse experiments were conducted under the supervision of Baylor College of Medicine Institutional Animal Care and Use Committee approved protocols.

Spatial transcriptomics in murine placentae. At E18.5 placentas were dissected and preserved in RNA-later at -25°C for 18 hours. Tissues were dehydrated by incubation in 30% sucrose at 4oC for 6 hours. At the Baylor College of Medicine Pathology and Histology Core placentas were fresh-frozen in optimal cutting temperature compound (OCT) with a liquid nitrogen and isopentane bath. Tissue blocks were cryosectioned at 10 microns for quality control H&E staining and 10x Genomics Visium Tissue Optimization. The Tissue Optimization slide revealed optimal permeabilization was observed after 12 minutes. Finally, the tissue blocks were cryosectioned and placed on Visium Gene Expression slides (v2) for H&E staining and Visium spatial transcriptomics library preparation following the manufacturer’s instructions. Visium libraries were sequenced (Illumina, NovaSeq 6000) at the Baylor College of Medicine Genomic and RNA Profiling Core. 

Spatial transcriptomics analysis. Spatial transcriptomes were aligned to a custom human and ZIKV transcriptome and transcripts counts were quantified using SpaceRanger (v1.3.0) and bash scripts. Downstream analyses were done in R (v4.0.5) using the Seurat (v4.0.6) scRNA-seq analysis package (33). H&E slides were annotated by a trained pathologist and transferred to spatial barcodes using the Loupe Browser (v5), followed by annotation based on PCA dimension reduction and clustering analysis to identify broad transcriptome niches. Spatial transcriptomes were then integrated with a published murine single-nuclei RNA-seq (snRNA-seq) dataset to assess transcriptome annotations (GEO accession GSE152248). The Marsh & Blelloch, eLife. (2020) snRNA-seq dataset included murine placentas sampled at E9.5, E10.5, E12.5, and E14.5 across 2,933 independently filtered (nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 10000 & percent.mt < 10 & percent.ribo < 20) high-quality transcriptomes. Full DE results are in Table S2 along with Y chromosome expression. 

SampleIDs: 
S07-ZIKV-Veh-1
S09-Mock-Enox-1
S12-Mock-Enox-2
S13-ZIKV-Enox-1
S31-ZIKV-Enox-2
S32-Mock-Veh-1
S33-ZIKV-Veh-2
S34-ZIKV-Veh-3

GEO Contents: Spatial transcriptomics data are available at GEO Accession GSE205631. Use the 10x Genomics bam2fastq package for manual re-alignments. Processed .h5 matricies are also deposited. The script to generate the custom mouse + ZIKV reference (mm10_and_NC_012532.1.v2) is provided here.

Downloading from GEO/SRA: For each sample, we have uploaded the 10x BAM file (possorted_genome_bam.bam) and processed count matrix (filtered_feature_bc_matrix.h5). Do not use fastq-dump or sam-dump to pull the data as the 10x Bam is different frum regular BAM files. The 10x BAM is different from regular BAM files, so the sam-dump option does not work. Please download the BAM file with the link in the 'Data access' --> 'Original format' section (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR23094132&display=data-access). Please see the 10x Genomics Support Page on the topic: https://kb.10xgenomics.com/hc/en-us/articles/360035250172-Why-do-I-not-see-any-FASTQs-or-incomplete-FASTQs-for-the-SRA-of-my-interest- . Then run bamtofastq (https://support.10xgenomics.com/docs/bamtofastq) in the CellRanger package to extract all the fastq files in the format necessary for CellRanger alignments. For example, 1 10x BAM produces 12 total fastq files.

# Repository contents
A) Images necessary for alignment of spatial transcriptomes

B) Bash scripts for creating a custom ZIKV+Mouse transcriptome and alignment using Spaceranger

C) R script for QC filtering and initial clustering of spatial transcriptomics data using Seurat, clustering of previously published placenta snRNA-seq data (Marsh et al) using Seurat, and integration of placenta spatial transcriptomics and snRNA-seq data using Seurat

D) R script for comparing visium data before and after filtering

E) R script for Pearson's correlation analysis
