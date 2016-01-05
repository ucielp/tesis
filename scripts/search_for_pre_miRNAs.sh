#!/bin/bash



for locus_id in $(awk {'print $2'} data/miRNA_locus.csv) ; 
do
    grep -A 1 $locus_id /home/uciel/lab/Datos/Procesamiento/precursores_especies/db/TAIR10_cdna/TAIR10_cdna_20101214_updated

done    
    
