#!/bin/bash

#~ for i in $(ls *.csv) ; 
#~ do
    #~ filename=$(basename "$i")
    #~ extension="${filename##*.}"
    #~ filename="${filename%.*}"
    #~ echo $filename
	#~ sort -u -t$'\t' -k1,1 $i > $filename.txt 
    #~ awk -F '\t' '{print ">"$1"\n"$2}' $filename.txt > $filename.tmp
    #~ #Hago solo las sustituciones en el fasta no en el nombre
    #~ sed 's/T/U/g' $filename.tmp > $filename.tmp2
    #~ sed 's/>U/>T/g' $filename.tmp2 > $filename.infile
#~ done
#~ 
#~ 
#~ for i in $(ls *.infile) ; 
#~ do
    #~ filename=$(basename "$i")
    #~ extension="${filename##*.}"
    #~ filename="${filename%.*}"
    #~ t_coffee -in=$filename.infile -mode=regular -method=slow_pair -output=score_html score_pdf fasta_aln -run_name=$filename
    #~ t_coffee -other_pg seq_reformat -in $filename.fasta_aln -in3 $filename.fasta_aln -action +3evaluate idmat -output=color_pdf -out=$filename.pdf
#~ done
#~ 
#~ rm *.score_html
#~ rm *.dnd 
#~ rm *.fasta_aln
#~ rm *.txt
#~ rm *.score_pdf
#~ rm *.tmp
#~ rm *.tmp2

#~ # Hago el RNAhybrid
# Ojo con el miR que lo tengo que ir cambiando 
# miR396 UUCCACAGCUUUCUUGAACUG
# miR408 AUGCACUGCCUCUUCCCUGGC
# miR159 UUUGGAUUGAAGGGAGCUCUA
# miR167 UGAAGCUGCCAGCAUGAUCUA
# miR398 UGUGUUCUCAGGUCACCCCUU
 
for i in $(ls *.infile) ; 
do
    echo '>miR398' > query_file 
    echo 'UGUGUUCUCAGGUCACCCCUU' >> query_file
    filename=$(basename "$i")
    extension="${filename##*.}"
    filename="${filename%.*}"
    grep -A 1 $filename /home/uciel/lab/programas/patmatch_2013/databases/TAIR/TAIR10_cdna_20101214_updated > target_file
    RNAhybrid -d 0 -m 12000 -t target_file -q query_file > $filename.rnahybrid
done


#~ # Ignoro GRF5 AT3G13960
#~ for i in $(ls *.csv) ; 
#~ do
    #~ filename=$(basename "$i")
    #~ extension="${filename##*.}"
    #~ filename="${filename%.*}"
    #~ echo $filename
    #~ echo grep "$(grep Athaliana $i | awk '{print $2}')" /home/uciel/lab/programas/patmatch_2013/databases/TAIR/TAIR10_cdna_20101214_updated 
#~ done
 #~ 
