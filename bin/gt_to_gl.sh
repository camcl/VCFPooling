#!/bin/bash

# Converts a VCF file written with GT to a VCF file with NOT log GL
# Usage example: 
# $ bash ~/PoolImpHuman/bin/bash-scripts/gt_to_gl.sh IMP.chr20.snps.gt.vcf.gz IMP.chr20.snps.gl.vcf.gz


# Updates the metadata in the header
# $ ~/PoolImpHuman/data/20200615$ bcftools view -h IMP.chr20.snps.gt.chunk10000.vcf.gz | sed 's/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">/##FORMAT=<ID=GL,Number=G,Type=Float,Description="Estimated Genotype Probability">/' | grep '##FORMAT'
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Estimated Genotype Probability">

# Replaces GT FORMAT field with GL for each variant
# $ bcftools view -H IMP.chr20.snps.gt.chunk10000.vcf.gz | sed 's/GT/GL/g'

# Converts the GT values to GL ones
# $ bcftools view -H IMP.chr20.snps.gt.chunk10000.vcf.gz | sed -e 's/1|1/0.0,0.0,1.0/g' -e 's/0|0/1.0,0.0,0.0/g' -e 's/1|0/0.0,1.0,0.0/g' -e 's/0|1/0.0,1.0,0.0/g'

fin=$1
fout=$2
fbname=$(basename "$fout" .vcf.gz)

bcftools view $fin | sed 's/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">/##FORMAT=<ID=GL,Number=G,Type=Float,Description="Estimated Genotype Probability">/' | sed 's/GT/GL/g' | sed -e 's/1|1/0.0,0.0,1.0/g' -e 's/0|0/1.0,0.0,0.0/g' -e 's/1|0/0.0,1.0,0.0/g' -e 's/0|1/0.0,1.0,0.0/g' > $fbname.vcf
bcftools view -Oz -o $fout $fbname.vcf
bcftools index -f $fout
rm $fbname.vcf


