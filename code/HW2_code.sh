!/bin/sh

#duomenų atsisiuntimas
wget -O ./inputs/SRR18214264.fastq.gz "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR18214264"

wget -O ./inputs/ERR204044.fastq.gz "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=ERR204044"

wget -O ./inputs/SRR15131330.fastq.gz "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR15131330"

efetch -db nucleotide -id CP015498  -format fasta > ./ref/CP015498.fasta

##Data QA/QC
# Data QC
for i in inputs/*.gz
do fastqc $i
done

#duomenys atrodo gana padoriai pirmi 10 nt  pasiskirto netolygiai, 
#bet nusprendŽiau juos palikti, kad būtų kuo ilgesnė sekos genomo surinkimui
#yra duplikaciju

#Trim your data

for i in inputs/*.fastq.gz
do
    trim_galore -j 6 $i -o trim/
done

#Repeat fastqc and evaluate if there were any changes.
for i in trim/*.gz
do fastqc $i
done

#šiek tiek sumažėjo sekų, iš esmės nelabai pakito

#Create a MultiQC plot that would include both fastqc results (raw and trimmed data). Upload this report to your git repository.
multiqc inputs/*fastqc.zip -o inputs/
multiqc trim/*fastqc.zip -o trim/

##Genome assembly
#Assemble your genomes using spades program
#Create an alternative assembly. Choose a program by yourself.
#SRR15131330
spades.py -t 4 --phred-offset 33 -s ./trim/SRR15131330_trimmed.fq.gz -o SRR15131330_assembly/spades

abyss-pe -j 6 name=SRR15131330 k=48 in='trim/SRR15131330_trimmed.fq.gz'


#ERR204044
spades.py -t 4 --phred-offset 33 -s ./trim/ERR204044_trimmed.fq.gz -o ERR204044_assembly/spades

abyss-pe -j 6 name=ERR204044 k=48 in='trim/ERR204044_trimmed.fq.gz'

#SRR18214264
spades.py -t 4 --phred-offset 33 -s ./trim/SRR18214264_trimmed.fq.gz -o SRR18214264_assembly/spades

abyss-pe -j 6 name=SRR18214264 k=48 in='trim/SRR18214264_trimmed.fq.gz'

#Evaluate your assemblies using Quast. Try to use CP015498 (NCBI) as a reference.
#In your code describe your results: do all assemblies look the same? What are the main results?
#Using quast results select the best/better assembly for each sample. In your code explain why you chose a specific assembly.


#komentaras:
#labai panaši kokybė su abiem programom, spades šiek tiek geresnė, 
#tai pasirinkau toliau naudoti su šia programa surinktus genomus, juos perkeliu į gonome direktoriją

#Orientate your contigs (make scaffolds) using ragtag program 
ragtag.py correct ./ref/CP015498.fasta ./genomes/SRR15131330.fasta
ragtag.py correct ./ref/CP015498.fasta ./genomes/ERR204044.fasta
ragtag.py correct ./ref/CP015498.fasta ./genomes/SRR18214264.fasta

#Using appropriate mapper, map original reads to you assemblies. 
#Evaluate mapping fraction as well as genome coverage from mapped reads 
#(and as in other questions, provide a comment on your results).

for i in ./trim/*_trimmed.fq.gz
do
    base=$(basename $i _trimmed.fq.gz)
    bwa index ./genomes/${base}_correct.fasta
    bwa mem -t 6 ./genomes/${base}_correct.fasta ./trim/${base}_trimmed.fq.gz 2> ./map_test/${base}_bwa_log.txt |\
    samtools view -bS -@ 6 | samtools sort -@ 6 -o ./map_test/${base}.bam
    samtools stats -in ./map_test/${base}.bam > ./map_test/map_stats_${base}.txt
    plot-bamstats ./map_test/map_stats_${base}.txt -p ./map_test/plots/${base}
done

#Susimappino virš 99% nuskaitymu

##Genome analysis and annotation

#Using Gepard tool create dotplots to show similarities/dissimilarities between your samples. 
#Describe, your results (in your code). In the last question you will have to upload dotplots, so save them.

#Using BUSCO analysis tool, evaluate your assemblies. Provide a short comment on BUSCO results.

#Using GeneMarkS-2 tool (http://exon.gatech.edu/genemark/genemarks2.cgi) predict genes in your assembled genomes.

#Using RAST genome annotation server, predict and annotate genes in your assemblies.

#Using CP015498 genes and proteins models as well as local blast, predict genes in your assemblies.

#Ranka atsisiunčiau CP015498 cDNA ir baltymų sekas ir 
#įkėliau į ref direktoriją (CP015498_cdna.fasta ir CP015498_prot.fasta)

makeblastdb -in ./ref/CP015498_prot.fasta -dbtype prot -parse_seqids

makeblastdb -in ./ref/CP015498_cdna.fasta -dbtype nucl -parse_seqids

for i in ./genomes/*_correct.fasta
do
    base=$(basename $i _correct.fasta)
    blastx -db ./ref/CP015498_prot.fasta -query ./genomes/${base}_correct.fasta > ./Blast/${base}_blastp.txt
done

for i in ./genomes/*_correct.fasta
do
    base=$(basename $i _correct.fasta)
    blastn -db ./ref/CP015498_cdna.fasta -query ./genomes/${base}_correct.fasta -outfmt "6 qseqid qstart qend" > ./Blast/${base}_blastn_3.txt
done

#At the moment you should have 3 gene predictions. Compare and describe them 
#(you should compare number of predicted genes and genes overlap. You don't have to include functional annotations).

library(dplyr)
# ERR204044 genų sk palyginimas:
ERR204044_blast <- read.delim(file="ERR204044_blastn_3.txt", header = FALSE, sep = "\t")
ERR204044_blast <- ERR204044_blast %>% 
  rename("seqnames" = "V1",
         "start" = "V2",
         "end" = "V3")

ERR204044_rast <- rtracklayer::import('ERR204044_rast.gtf')
ERR204044_rast = as.data.frame(ERR204044_rast)

ERR204044_GeneMarkS <- rtracklayer::import('ERR204044_GeneMarkS.gtf')
ERR204044_GeneMarkS = as.data.frame(ERR204044_GeneMarkS)

ERR204044_blast %>% distinct() %>% nrow()
#2720

nrow(ERR204044_rast)
#2601

nrow(ERR204044_GeneMarkS)
#2360

inner_join(ERR204044_blast, ERR204044_rast, by=c("seqnames", "start", "end")) %>% nrow()
#1478

inner_join(ERR204044_blast, ERR204044_GeneMarkS, by=c("seqnames", "start", "end")) %>% nrow()
#1478

inner_join(ERR204044_rast, ERR204044_GeneMarkS, by=c("seqnames", "start", "end")) %>% nrow()
#2097

inner_join(ERR204044_rast, ERR204044_GeneMarkS, by=c("seqnames", "start", "end")) %>% inner_join(ERR204044_blast, by=c("seqnames", "start", "end")) %>% nrow() 
#1410

# SRR15131330 genų sk palyginimas:

SRR15131330_blast <- read.delim(file="SRR15131330_blastn_3.txt", header = FALSE, sep = "\t")
SRR15131330_blast <- SRR15131330_blast %>% 
  rename("seqnames" = "V1",
         "start" = "V2",
         "end" = "V3")

SRR15131330_rast <- rtracklayer::import('SRR15131330_rast.gtf')
SRR15131330_rast = as.data.frame(SRR15131330_rast)

SRR15131330_GeneMarkS <- rtracklayer::import('SRR15131330_GeneMarkS.gtf')
SRR15131330_GeneMarkS = as.data.frame(SRR15131330_GeneMarkS)

SRR15131330_blast %>% distinct() %>% nrow()
#3755
nrow(SRR15131330_rast)
#2733
nrow(SRR15131330_GeneMarkS)
#2538

inner_join(SRR15131330_blast, SRR15131330_rast, by=c("seqnames", "start", "end")) %>% nrow()
#2340

inner_join(SRR15131330_blast, SRR15131330_GeneMarkS, by=c("seqnames", "start", "end")) %>% nrow()
#2144

inner_join(SRR15131330_rast, SRR15131330_GeneMarkS, by=c("seqnames", "start", "end")) %>% nrow()
#2155

inner_join(SRR15131330_rast, SRR15131330_GeneMarkS, by=c("seqnames", "start", "end")) %>% inner_join(SRR15131330_blast, by=c("seqnames", "start", "end")) %>% nrow() 
#2037

# SRR18214264 genų sk palyginimas:

SRR18214264_blast <- read.delim(file="SRR18214264_blastn_3.txt", header = FALSE, sep = "\t")
SRR18214264_blast <- SRR18214264_blast %>% 
  rename("seqnames" = "V1",
         "start" = "V2",
         "end" = "V3")

SRR18214264_rast <- rtracklayer::import('SRR18214264_rast.gtf')
SRR18214264_rast = as.data.frame(SRR18214264_rast)

SRR18214264_GeneMarkS <- rtracklayer::import('SRR18214264_GeneMarkS.gtf')
SRR18214264_GeneMarkS = as.data.frame(SRR18214264_GeneMarkS)

SRR18214264_blast %>% distinct() %>% nrow()
#2568
nrow(SRR18214264_rast)
#2541
nrow(SRR18214264_GeneMarkS)
#2320

inner_join(SRR18214264_blast, SRR18214264_rast, by=c("seqnames", "start", "end")) %>% nrow()
#1595

inner_join(SRR18214264_blast, SRR18214264_GeneMarkS, by=c("seqnames", "start", "end")) %>% nrow()
#1545

inner_join(SRR18214264_rast, SRR18214264_GeneMarkS, by=c("seqnames", "start", "end")) %>% nrow()
#2036

inner_join(SRR18214264_rast, SRR18214264_GeneMarkS, by=c("seqnames", "start", "end")) %>% inner_join(SRR18214264_blast, by=c("seqnames", "start", "end")) %>% nrow() 
#1496

#Create a phylogenetic tree from 16S sequences. 
#In your tree include reference genome 16S sequences and an outgroup (you can use Staphylococcus as an outgroup)

#Išsitraukiau mano surinktų genomų 16S sekas iš RAST, o ref ir outgrout iš NCBI.
#Susikėliau visas sekas į UGENE ir sulygiavau su MUSCE algoritmu (palyginimo failas ~/HW2/palyginiai/16S.aln). 
#Su UGENE nubraižiau filogenetinė medį (medis: ~/HW2/palyginiai/16S.pgn)

#From your/BUSCO predictions, select five proteins and create a multi-gene tree. 
#Again, include a reference genome as well as an outgroup (organisms in this and the previous question should match).

5 baltymai:
NAD-dependent oxidoreductase
DNA topoisomerase IV subunit A
DNA primase
DNA polymerase III subunits gamma and tau
LSU ribosomal maturation GTPase RbgA

sekas paėmiau iš rast, o ref ir outgtoup iš NCBI

for i in ./palyginiai/*.fasta
do
    base=$(basename $i .fasta)
    muscle -in $i -out ./palyginiai/${base}.aln
done

Sulygiuotus baltymus sudėjau į vieną failą: ~/HW2/palyginiai/5_baltymaai.aln
Su UGENE nubraižiau filogenetinė medį (medis: ~/HW2/palyginiai/5_baltymai.pgn)