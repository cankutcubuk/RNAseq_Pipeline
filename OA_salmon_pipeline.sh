#!/bin/bash

####################################################
#         Salmon pipeline for bulk RNASeq		   #
#         cankutcubuk [at] {gmail} [dot] {com}     #       
#                     2025				       	   #      
#          https://www.qmul.ac.uk/whri/emr/	       #       
#                @ EMR-QMUL, London, UK		       #    
#################################################### 

outdir=$1
threads=$2
fastaPath=$3
sampleList=$4

# example usage: ./OA_salmon_pipeline.sh "/OA_RNAseq/OA_preprocessed_realignment2025/"  50 "/OA_RNAseq/40-37xxxxx/" "/OA_RNAseq/40-37xxxxx/OA_RNAseq_files.txt"
# OA_RNAseq_files.txt has one column. In each line, there is one sampleID. In the code below, "_R1_001.fastq.gz and _R2_001.fastq.gz" need to be modified according to suffixes used in the fastq files.

mkdir $outdir
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz -O $outdir"/salmon-1.10.0_linux_x86_64.tar.gz"
mkdir "/home/"$USER"/tmp/"
tar -xvzf $outdir"/salmon-1.10.0_linux_x86_64.tar.gz" -C /home/$USER/tmp/
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O $outdir"Homo_sapiens.Ensembl_v113.GRCh38.cdna.all.fa.gz"
mkdir $outdir"/salmon_index_ensembl_v113_GRCh38_cdna_all/"
/home/$USER/tmp/salmon-latest_linux_x86_64/bin/salmon index -t $outdir"Homo_sapiens.Ensembl_v113.GRCh38.cdna.all.fa.gz" -i $outdir"salmon_index_ensembl_v113_GRCh38_cdna_all/" -k 31 -p $threads

while read filename; do

/home/$USER/tmp/salmon-latest_linux_x86_64/bin/salmon quant -p $threads -i $outdir"/salmon_index_ensembl_v113_GRCh38_cdna_all/" \
 --gcBias --validateMappings --useVBOpt --libType A -1 $fastaPath/$filename"_R1_001.fastq.gz" -2 $fastaPath/$filename"_R2_001.fastq.gz" -o  $outdir/$filename/"salmon"

done < $sampleList

echo >&2 '
************
*** DONE ***
************'

