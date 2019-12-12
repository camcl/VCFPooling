#!/bin/sh
# Description: convert GL to GT, also adding the relevant ##FORMAT-line for GT in the VCF-header
# Usage: $1, $2 are command line parameters representing: VCF file with GL data to convert, VCF file with GT data to output
# Ex. $ sh path/to/script/awk_gl22gt.sh IMP.chr20.pooled.imputed.gl_default.chunk1000.vcf.gz IMP.chr20.pooled.gt.imputed.gl_default.chunk1000.vcf.gz

headln=$(bcftools view -h "$1" | wc -l)
samplesln=$(bcftools view -h "$1" | tail -1)
bcftools view -h "$1" | head -"$((headln - 1))" > temp.vcf
echo '##FORMAT=<ID=GL,Number=G,Type=Float,Description="three log10-scaled likelihoods for RR,RA,AA genotypes">\n'"$samplesln""\n" | tee -a temp.vcf  
bcftools view -H "$1" | awk -v fileout=temp.body.vcf 'BEGIN{ FS=OFS="\t" } {gsub(/1,0,0/, "0\|0"); gsub(/0,1,0/, "1\|0"); gsub(/0,0,1/, "1\|1"); gsub(/GL/, "GT", $9); print | "tee " fileout}' 
cat temp.body.vcf >> temp.vcf
bcftools view -Oz -o "$2" temp.vcf
bcftools index -f "$2"
rm temp.body.vcf temp.vcf






