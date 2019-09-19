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
cellranger aggr --id=AC2262 --csv=/scrapp/kkim1/scrna/AC2262_aggr.csv --normalize=mapped 
date