!/bin/sh

#duomenų atsisiuntimas
wget -O ./inputs/SRR18214264.fastq.gz "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR18214264"

wget -O ./inputs/ERR204044.fastq.gz "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=ERR204044"

wget -O ./inputs/SRR15131330.fastq.gz "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR15131330"

efetch -db nucleotide -id CP015498  -format fasta > ./ref/CP015498.fasta

# Data QC
for i in inputs/*.gz
do fastqc $i
done

#In your code write a short comment about your data quality.

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

#In your code write a short comment about your results.


#Create a MultiQC plot that would include both fastqc results (raw and trimmed data). Upload this report to your git repository.
multiqc trim/*fastqc.zip -o trim/

#Genome assembly
#Assemble your genomes using spades program
#SRR15131330
spades.py -t 4 --phred-offset 33 -s ./trim/SRR15131330_trimmed.fq.gz -o SRR15131330_assembly/spades

abyss-pe -j 6 name=SRR15131330 k=48 in='trim/SRR15131330_trimmed.fq.gz'


#ERR204044
spades.py -t 4 --phred-offset 33 -s ./trim/ERR204044_trimmed.fq.gz -o ERR204044_assembly/spades

abyss-pe -j 6 name=ERR204044 k=48 in='trim/ERR204044_trimmed.fq.gz'

#SRR18214264
spades.py -t 4 --phred-offset 33 -s ./trim/SRR18214264_trimmed.fq.gz -o SRR18214264_assembly/spades

abyss-pe -j 6 name=SRR18214264 k=48 in='trim/SRR18214264_trimmed.fq.gz'

#labai pana6i kokyb4 abiem programom, spades \šiek tiek geresne, tai pasirinta toliau naudoti su šia programa surinktus genomus, juos perkeliu i gonome direktoriją

#Orientate your contigs (make scaffolds) using ragtag program 
ragtag.py correct ./ref/CP015498.fasta ./genomes/SRR15131330.fasta
ragtag.py correct ./ref/CP015498.fasta ./genomes/ERR204044.fasta
ragtag.py correct ./ref/CP015498.fasta ./genomes/SRR18214264.fasta

for i in ./trim/*_trimmed.fq.gz
do
    base=$(basename $i _trimmed.fq.gz)
    bwa index ./genomes/${base}_correct.fasta
    bwa mem -t 6 ./genomes/${base}_correct.fasta ./trim/${base}_trimmed.fq.gz 2> ./map_test/${base}_bwa_log.txt |\
    samtools view -bS -@ 6 | samtools sort -@ 6 -o ./map_test/${base}.bam
    samtools stats -in ./map_test/${base}.bam > ./map_test/map_stats_${base}.txt
done
