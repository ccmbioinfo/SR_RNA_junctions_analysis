#!/bin/bash
module load bedtools/2.29.2

bam=$1
chrom_file=$2
bed_file=$3

coverage_stat=`dirname $bam`"/"`basename $bam mdup.bam`"coverage.stat"

samtools view -b -@ 8 -q 255 $bam|coverageBed -S -split -hist -sorted -g $chrom_file -a $bed_file -b -|awk '$1=="all"'|awk 'BEGIN{sum=0;oX=0;tX=0;fX=0;}{sum=sum+$2*$3;if ($2>=50){fX=fX+$3;}if ($2>=10){tX=tX+$3;}if($2>=1){oX=oX+$3;}}END{printf "%0.2f\t%0.2f\t%0.2f\t%0.2f\n",sum/$4,100*oX/$4,100*tX/$4, 100*fX/$4}' > $coverage_stat
