samtools sort results/alignment.bam -o results/sorted.bam

samtools depth results/sorted.bam > output.txt