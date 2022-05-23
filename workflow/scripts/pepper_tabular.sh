#!/bin/bash

file=$argv$1
name=$argv$2
sed '/^CHROM/ s/'${name}.'//g' ${file} | awk 'BEGIN{FS=OFS="\t"} $3=="."{$3="PEPPER"} 1'

