# this is the alignment tool

# firsty index the reference
bwa index data/chrX.fa.gz

bwa mem data/chrX.fa.gz data/tu.r1.fq.gz data/tu.r2.fq.gz > results/alignment.bam

# samtools view -b results/alignment.bam chrX:20000000-40000000 > results/view.bam

echo "DONE!"