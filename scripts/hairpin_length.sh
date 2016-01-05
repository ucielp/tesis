# 
awk 'NR%2==0''{print length($0)}' /home/uciel/lab/plant_databases/miRBASE/hairpin_Viridiplantae.fa.prepared > out/hairpin_Viridiplantae_lengths.csv
# En R

x <- read.table("hairpin_Viridiplantae_lengths.csv",header=T,sep=',')
x_no_outliers <- x[x<1000]


svg(filename="hairpin_Viridiplantae_lengths_no_outliers.svg", 
    width=5, 
    height=5, 
    pointsize=10)
    boxplot(x_no_outliers)
dev.off()
