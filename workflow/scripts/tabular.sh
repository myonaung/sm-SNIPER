#!/bin/bash

### This is to modify tabular vcf file for downstream manipulation####
longshot_out=$argv$1
pepper_out=$argv$2


path=out/tabular
cd ${path}


for f in *_longshot_amplicon.txt; do
    name=$(echo $f| cut -d_ -f1)
    sed '/^CHROM/ s/'${name}.'//g' ${f} | awk 'BEGIN{FS=OFS="\t"} $3=="."{$3="longshot"} 1' > ${longshot_out}
done 

for f in *_PEPPER_amplicon.txt; do
    name=$(echo $f| cut -d_ -f1)
    sed '/^CHROM/ s/'${name}.'//g' ${f} | awk 'BEGIN{FS=OFS="\t"} $3=="."{$3="PEPPER"} 1' > ${pepper_out}
done 
