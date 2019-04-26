#!/bin/bash

SAMPLE=$1 

THREADS=6

READ_DIR=/stor/home/russd/scratch/10x_snps/HG3_Exome/raw_reads

REF_FASTA=/stor/home/russd/work/references/GRCh38/Homo_sapiens_assembly38.fasta
REF_BWA=/stor/home/russd/work/references/GRCh38/bwa/Homo_sapiens_assembly38.fasta
PICARD_JAR=/home/russd/software/picard-tools-2.1.1/picard.jar
GATK_JAR=/home/russd/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar


## Let's Do This 

echo 'Aligning Reads for '$SAMPLE
bwa mem -M -R "@RG\tID:$SAMPLE\tLB:$SAMPLE\tSM:$SAMPLE\tPL:illumina\tPM:nextseq\t" \
-t $THREADS \
$REF_BWA \
$READ_DIR/$SAMPLE\_1.fastq.gz \
$READ_DIR/$SAMPLE\_2.fastq.gz \
> ${SAMPLE}_aligned_reads.sam


echo 'Sorting Sam File > Bam'
java -jar $PICARD_JAR \
SortSam \
INPUT=${SAMPLE}_aligned_reads.sam \
OUTPUT=${SAMPLE}_sorted_reads.bam \
SORT_ORDER=coordinate

echo 'Collecting Alignment Metrics '
java -jar $PICARD_JAR \
CollectAlignmentSummaryMetrics \
R=$REF_FASTA \
I=${SAMPLE}_sorted_reads.bam \
O=${SAMPLE}_alignment_metrics.txt

java -jar $PICARD_JAR \
CollectInsertSizeMetrics \
INPUT=${SAMPLE}_sorted_reads.bam \
OUTPUT=${SAMPLE}_insert_metrics.txt \
HISTOGRAM_FILE=${SAMPLE}_insert_size_histogram.pdf


echo 'Marking Duplicates'
java -jar $PICARD_JAR \
MarkDuplicates \
INPUT=${SAMPLE}_sorted_reads.bam \
OUTPUT=${SAMPLE}_dedup_reads.bam \
REMOVE_DUPLICATES=true \
METRICS_FILE=${SAMPLE}_metrics.txt

echo 'Indexing Bam'
java -jar $PICARD_JAR \
BuildBamIndex \
INPUT=${SAMPLE}_dedup_reads.bam

echo 'Creating Realignment Targets'
java -jar $GATK_JAR \
-T RealignerTargetCreator \
-R $REF_FASTA \
-I ${SAMPLE}_dedup_reads.bam \
-o ${SAMPLE}_realignment_targets.list

echo 'Realigning Targets'
java -jar $GATK_JAR \
-T IndelRealigner \
-R $REF_FASTA \
-I ${SAMPLE}_dedup_reads.bam \
-targetIntervals ${SAMPLE}_realignment_targets.list \
-o ${SAMPLE}_realigned_reads.bam


echo 'Calculating Base Recalibration'
java -jar $GATK_JAR \
-nct $THREADS \
-T BaseRecalibrator \
-R $REF_FASTA \
-I ${SAMPLE}_realigned_reads.bam \
-knownSites /stor/home/russd/work/references/GRCh38/dbsnp_146.hg38.vcf.gz \
-knownSites /stor/home/russd/work/references/GRCh38/hapmap_3.3_grch38_pop_stratified_af.vcf.gz \
-knownSites /stor/home/russd/work/references/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-o ${SAMPLE}_BQSR.recal.table


echo 'Recalibrating Bases'
java -jar $GATK_JAR \
-nct $THREADS \
-T PrintReads \
-R $REF_FASTA \
-I ${SAMPLE}_realigned_reads.bam \
-BQSR ${SAMPLE}_BQSR.recal.table \
-o ${SAMPLE}_recal_reads.bam

