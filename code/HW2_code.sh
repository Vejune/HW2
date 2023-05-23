#!/bin/sh

#duomenų atsisiuntimas
wget -O ./inputs/SRR18214264.fastq.gz "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR18214264"

wget -O ./inputs/ERR204044.fastq.gz "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=ERR204044"

wget -O ./inputs/SRR15131330.fastq.gz "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR15131330"

# Data QC
for i in inputs/*.gz
do fastqc $i
done

multiqc inputs/*fastqc.zip -o inputs/

#trim
for i in inputs/*.fastq.gz
do
    trim_galore -j 6 $i -o trim/
done

#trim QC
for i in trim/*.gz
do fastqc $i
done

multiqc trim/*fastqc.zip -o trim/

#Genome assembly
#Assemble your genomes using spades program
spades.py -t 4 --phred-offset 33 -s ./trim/SRR15131330_trimmed.fq.gz -o outputs/