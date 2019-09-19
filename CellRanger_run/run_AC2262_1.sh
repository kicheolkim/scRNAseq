#!/bin/bash
#
#$ -S /bin/bash
#$ -o /scrapp/kkim1/log
#$ -e /scrapp/kkim1/log/stderr
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=64G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=500G
#$ -l h_rt=24:00:00
date
hostname
cellranger count --id=AC2262_1 --transcriptome=/scrapp/kkim1/refdata-cellranger-mm10-1.2.0 --fastqs=/scrapp/kkim1/scrna/fastq/gi_AC2262 --localcores=48 --sample=gi_AC2262_SI-GA-D3_1,gi_AC2262_SI-GA-D3_2,gi_AC2262_SI-GA-D3_3,gi_AC2262_SI-GA-D3_4
date