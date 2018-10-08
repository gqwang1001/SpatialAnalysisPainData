#!/bin/csh
#$ -cwd
#$ -pe local 20
#$ -l mem_free=10G
#$ -l h_vmem=11G
#$ -M gqwang1001@gmail.com
#$ -be

R CMD BATCH main.R
