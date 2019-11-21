#!/bin/bash -l
## Grid Engine Example Job Script  

# -- Begin SGE embedded arguments --
#$ -V
#Pass all environment variables to job
#$ -cwd
#Use current working directory

#Grab the patient number from the filename
#long_file_name=$(pwd)

#$ -N LASV_Dedup
# Name of script

#$ -j y
#Combine standard error and output files.

#$-q highmem.q
#Use the all.q queue, and not any other queue.

#$-pe smp 4-16
#Ask for a parallel environment for multi-threading.

# -- End SGE embedded arguments --

module load bwa/0.7.7
module load samtools/1.2
module load picard/1.128
module load BEDTools/2.17.0
module load cutadapt/1.8.3
module load prinseq/0.20.3

# create temp directory for work on /scratch
mkdir -p /scratch/evk3/NIPV/
scratch='/scratch/evk3/NIPV'


SEEDFILE=./file_names.txt
file_num=$(awk "NR==$SGE_TASK_ID" $SEEDFILE)

echo "Patient number" $file_num
echo "SGE Value " $SGE_TASK_ID

sample_num=$(sed -r 's/_.*$//g' <<< "$file_num")
echo "Sample number is: " $sample_num

L1_READ1=$file_num
L1_READ2=$(awk "NR==($SGE_TASK_ID + 1)" $SEEDFILE)

echo $L1_READ1
echo $L1_READ2


#L1_READ1='./LASV-Togo-ILM-pool_S2_L001_R1_001.fastq.gz'
#L1_READ2='./LASV-Togo-ILM-pool_S2_L001_R2_001.fastq.gz'
#sample_num="LASV-Togo-ILM-pool"

cp /gb_JN808863_Organism_Nipah_virus_NIVBGD2008RAJBARI.fasta "$scratch"/


#Remove Illumina TruSeq adaptors from PE reads:
echo "Starting cutadapt"
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
         -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
         -o "$scratch"/"$sample_num"_R1_cutadapt.fastq.gz -p "$scratch"/"$sample_num"_R2_cutadapt.fastq.gz \
         $L1_READ1 $L1_READ2 

echo "Gunzipping now!"
gunzip -c "$scratch"/"$sample_num"_R1_cutadapt.fastq.gz > "$scratch"/"$sample_num"_R1_cutadapt.fastq
gunzip -c "$scratch"/"$sample_num"_R2_cutadapt.fastq.gz > "$scratch"/"$sample_num"_R2_cutadapt.fastq

#Remove low quality reads:
echo "starting printseq-lite"
prinseq-lite -fastq "$scratch"/"$sample_num"_R1_cutadapt.fastq -fastq2 "$scratch"/"$sample_num"_R2_cutadapt.fastq -min_qual_mean 25 -trim_qual_right 20 -min_len 50 -out_good "$scratch"/"$sample_num"_trimmed

echo "Indexing reference sequence using BWA"
bwa index "$scratch"/gb_JN808863_Organism_Nipah_virus_NIVBGD2008RAJBARI.fasta -p "$scratch"/gb_JN808863_Organism_Nipah_virus_NIVBGD2008RAJBARI.fasta

echo "Mapping reads to reference genome"
bwa mem -t $NSLOTS "$scratch"/gb_JN808863_Organism_Nipah_virus_NIVBGD2008RAJBARI.fasta "$scratch"/"$sample_num"_trimmed_1.fastq "$scratch"/"$sample_num"_trimmed_2.fastq > "$scratch"/"$sample_num"_S.sam

# Previously generated samtools index of reference genome.  Generates *.fai file and only need to do 1X.
# This step may not be necessary with the samtool pipeline below:
echo "Indexing Ebo genome with samtools"
samtools faidx "$scratch"/gb_JN808863_Organism_Nipah_virus_NIVBGD2008RAJBARI.fasta

echo "Starting samtools - convert SAM to BAM"
samtools view -S -b -o "$scratch"/"$sample_num"_S.bam "$scratch"/"$sample_num"_S.sam

echo "Starting samtools sort BAM file"
samtools sort -@ $NSLOTS "$scratch"/"$sample_num"_S.bam "$scratch"/"$sample_num"_S.sorted

echo "Starting samtools index BAM file"
samtools index "$scratch"/"$sample_num"_S.sorted.bam

JAVA_OPTS='-Xmx50g'
TMP_DIR=/tmp

echo "Removing duplicates with Picard."
picard MarkDuplicates INPUT="$scratch"/"$sample_num"_S.sorted.bam REMOVE_DUPLICATES=true METRICS_FILE="$scratch"/"$sample_num"_L1_RM_S.metrics.txt OUTPUT="$scratch"/"$sample_num"_S.dedup.sorted.bam

#Collect info about bases:
picard CollectQualityYieldMetrics INPUT="$scratch"/"$sample_num"_S.dedup.sorted.bam OUTPUT="$scratch"/"$sample_num"_S.dedup.sorted.bam.metrics


bedtools genomecov -d -split -ibam "$scratch"/"$sample_num"_S.sorted.bam  > "$scratch"/"$sample_num"_S.sorted.bam.txt

bedtools genomecov -d -split -ibam "$scratch"/"$sample_num"_S.dedup.sorted.bam  > "$scratch"/"$sample_num"_S.dedup.sorted.bam.txt

#copy node /scratch/evk3/ebo back to home dir
#copy results from node /scratch/evk3/ebo back to home dir
cp -R "$scratch"/"$sample_num"_* /Nipah/Enriched_clinical_n26

#clean up scratch directory
rm -rf /scratch/evk3

module unload bwa/0.7.7
module unload samtools/1.2
module unload picard/1.128
module unload BEDTools/2.17.0
module unload cutadapt/1.8.3
module unload prinseq/0.20.3

echo "Script finish"
