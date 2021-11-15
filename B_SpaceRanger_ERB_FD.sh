## SpaceRanger_ERB_v1.sh
# Enrico Barrozo, BCM, Aagaard Lab, 2021


## download latest spaceranger and mouse reference
## https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/installation
curl -o spaceranger-1.3.0.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-1.3.0.tar.gz?Expires=1634802546&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvc3BhdGlhbC1leHAvc3BhY2VyYW5nZXItMS4zLjAudGFyLmd6IiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNjM0ODAyNTQ2fX19XX0_&Signature=TFCSeNH0mg4Ym9ugg3dC~a41SRlGCP--XLaxFc~1tXFYA-TMaGVZgsMT1OKS9RLlxrkftF0IVSVPID4VElOt0hdapliDLGIGZazBueFfyWAlElOWEmkJdIvJUBKlN8cB04YQ9Gg~73Dt8LdPMUdPb4IFGKqG6rlZvOeJDqwKrKuZ1Yfo0GXC~eoIVXwNdSU9CBI7DGhyE1VGEOff6Hl9YVh~crYLFRTrZJ808aLD0UuT~HkdULp1P9dsrnxdmXUBNaKdmCBw5sBqyiB95adw5Gnyxyh-a2zmY6qlp01vlOcI0R4xs6lb4dTnq-Inxpr2sgldZHt1Q6pD3QiwBR-7uw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

curl -O https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-mm10-2020-A.tar.gz

cd
tar -xzvf spaceranger-1.3.0.tar.gz
tar -xzvf refdata-gex-mm10-2020-A.tar.gz

################################################# Make Custom References ###########################################################
## Download the latest and greatest references for each virus (.gtf and .fasta)
	## SARS-CoV-2 is NC_045512.2
	## ZIKV clinical strain we are using is NC_045512.2 ; revised gtf annotation details can be found at https://github.com/ebarrozo/Aagaard_ZIKV_AGOHITSCLIP/blob/main/A_zikv-hg38_eb_v3.sh
cd /home/ebarrozo
mkdir visium
cd visium
mkdir data
mkdir docs
mkdir scripts
mkdir results


########### mouse+ZIKV ref
cd /home/ebarrozo/dang/docs/ref
gffread GRCh38_latest_genomic.gff -T -o GRCh38_latest_genomic.gtf
gffread NC_012532.1.gff3 -T -o /home/ebarrozo/dang/docs/ref/NC_012532.1.gtf

cd
export PATH=/home/ebarrozo/cellranger-5.0.1:$PATH

cd /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs
spaceranger mkref --genome=mm10 \
--fasta=/home/ebarrozo/refdata-gex-mm10-2020-A/fasta/genome.fa \
--genes=/home/ebarrozo/refdata-gex-mm10-2020-A/genes/genes.gtf \
--nthreads=32 --genome=NC_012532.1.v2 --fasta=/home/ebarrozo/dang/docs/ref/NC_012532.1.fasta \
--genes=/home/ebarrozo/dang/docs/ref/NC_012532.1.v2.filtered.gtf
## mouse + ZIKV can now  be found at /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/mm10_and_NC_012532.1.v2



################################################# Alignment in SpaceRanger ###########################################################
cd
export PATH=/home/ebarrozo/cellranger-6.1.1:$PATH
export PATH=/home/ebarrozo/spaceranger-1.3.0:$PATH

cd /home/ebarrozo/visium/docs
	## make sure the .tif files of each slide from the pathology core are transfered here

cd /home/ebarrozo/visium/data
	## make sure are all fastq files are named correctly
			# Consistent with bcl2fastq/mkfastq, e.g. "MySample_S1_L001_R1_001.fastq.gz" "MySample_S1_L001_R2_001.fastq.gz"


cd /home/ebarrozo/visium/results
## mkfastq done by GARP/multiomics core
	# spaceranger mkfastq

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
cd /home/ebarrozo/visium/results


################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
cd /home/ebarrozo/visium/results
spaceranger count --id=S07-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/mm10_and_NC_012532.1.v2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25752_07_MPL \
                   --sample=KA_25752 \
                   --image=/home/ebarrozo/visium/docs/7.tif \
                   --slide=V10S21-369 \
                   --area=D1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-369-D1.json \
                   --localcores=48 \
                   --localmem=300

## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S07-auto \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/mm10_and_NC_012532.1.v2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25752_07_MPL \
                   --sample=KA_25752 \
                   --image=/home/ebarrozo/visium/docs/7.tif \
                   --slide=V10S21-369 \
                   --area=D1 \
                   --localcores=48 \
                   --localmem=300
################################################# 

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S09-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/mm10_and_NC_012532.1.v2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25753_09_MPL \
                   --sample=KA_25753 \
                   --image=/home/ebarrozo/visium/docs/9.tif \
                   --slide=V10S21-368 \
                   --area=A1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-368-A1.json \
                   --localcores=48 \
                   --localmem=300

## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S09-auto \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/mm10_and_NC_012532.1.v2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25753_09_MPL \
                   --sample=KA_25753 \
                   --image=/home/ebarrozo/visium/docs/9.tif \
                   --slide=V10S21-368 \
                   --area=A1 \
                   --localcores=48 \
                   --localmem=300
################################################# 

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S12-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/mm10_and_NC_012532.1.v2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25754_12_MPL \
                   --sample=KA_25754 \
                   --image=/home/ebarrozo/visium/docs/12.tif \
                   --slide=V10S21-368 \
                   --area=B1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-368-B1.json \
                   --localcores=48 \
                   --localmem=300

## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S12-auto \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/mm10_and_NC_012532.1.v2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25754_12_MPL \
                   --sample=KA_25754 \
                   --image=/home/ebarrozo/visium/docs/12.tif \
                   --slide=V10S21-368 \
                   --area=B1 \
                   --localcores=48 \
                   --localmem=300
#################################################

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S13-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/mm10_and_NC_012532.1.v2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25755_13_MPL \
                   --sample=KA_25755 \
                   --image=/home/ebarrozo/visium/docs/13.tif \
                   --slide=V10S21-368 \
                   --area=C1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-368-C1.json \
                   --localcores=48 \
                   --localmem=300

## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S13-auto \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/mm10_and_NC_012532.1.v2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25755_13_MPL \
                   --sample=KA_25755 \
                   --image=/home/ebarrozo/visium/docs/13.tif \
                   --slide=V10S21-368 \
                   --area=C1 \
                   --localcores=48 \
                   --localmem=300
################################################# 

A successful spaceranger count run concludes with a message similar to this:
2016-11-10 16:10:09 [runtime] (join_complete)   ID.sample345.SPATIAL_RNA_COUNTER_CS.SPATIAL_RNA_COUNTER_CS.SUMMARIZE_REPORTS
 
Outputs:
- Run summary HTML:                         /opt/sample345/outs/web_summary.html
- Outputs of spatial pipeline:              /opt/sample345/outs/spatial
- Run summary CSV:                          /opt/sample345/outs/metrics_summary.csv
- BAM:                                      /opt/sample345/outs/possorted_genome_bam.bam
- BAM index:                                /opt/sample345/outs/possorted_genome_bam.bam.bai
- Filtered feature-barcode matrices MEX:    /opt/sample345/outs/filtered_feature_bc_matrix
- Filtered feature-barcode matrices HDF5:   /opt/sample345/outs/filtered_feature_bc_matrix.h5
- Unfiltered feature-barcode matrices MEX:  /opt/sample345/outs/raw_feature_bc_matrix
- Unfiltered feature-barcode matrices HDF5: /opt/sample345/outs/raw_feature_bc_matrix.h5
- Secondary analysis output CSV:            /opt/sample345/outs/analysis
- Per-molecule read information:            /opt/sample345/outs/molecule_info.h5
- Loupe Browser file:                       /opt/sample345/outs/cloupe.cloupe
- Spatial Enrichment using Moran's I file:  /opt/sample345/outs/spatial_enrichment.csv
 
Pipestance completed successfully!