# 10X output to GATK-ready 

INPUT_BAM=$1 
echo 'preparing bam file:' $INPUT_BAM

INPUT_PREFIX=${INPUT_BAM/.bam/}
echo 'input prefix / sample name: ' $INPUT_PREFIX


THREADS=6
REF_FASTA=~/work/references/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa



java -jar ~/software/picard.jar AddOrReplaceReadGroups \
  I=$INPUT_BAM \
  O=$INPUT_PREFIX.rg.bam \
  RGID=1 \
  RGLB=lib1 \
  RGPL=illumina \
  RGPU=unit1 \
  RGSM=$INPUT_PREFIX
 

# MARK DUPLICATES
java -jar ~/software/picard.jar MarkDuplicates INPUT=$INPUT_PREFIX.rg.bam OUTPUT=$INPUT_PREFIX.dups.bam M=dup.metrics   # CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics 
# # mark duplicates using UMI info: 
# #java -jar ~/software/picard.jar UmiAwareMarkDuplicatesWithMateCigar INPUT=possorted_genome_bam.bam OUTPUT=possorted_genome_bam.dupsumi.bam M=dupumi.metrics  UMI_METRICS_FILE=dupumi.umi.metrics BARCODE_TAG=BC UMI_TAG_NAME=UR

# # sort and index  
samtools sort $INPUT_PREFIX.dups.bam > $INPUT_PREFIX.dups.sorted.bam
samtools index $INPUT_PREFIX.dups.sorted.bam $INPUT_PREFIX.dups.sorted.bam.bai 


# # indel realignment target creation 
echo 'realigning around indels - creating targets '
java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  -R $REF_FASTA \
  -known /stor/home/russd/work/references/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -I $INPUT_PREFIX.dups.sorted.bam \
  -o $INPUT_PREFIX.dups.sorted.intervals

# # indel realignment 
echo 'realigning around indels - performing realignment'
java -Xmx8G -Djava.io.tmpdir=/tmp -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  -R $REF_FASTA \
  -targetIntervals $INPUT_PREFIX.dups.sorted.intervals \
  --maxReadsForRealignment 50000 \
  -known /stor/home/russd/work/references/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -I $INPUT_PREFIX.dups.sorted.bam \
  -o $INPUT_PREFIX.dups.sorted.indelrealigned.bam



# # BASE QUALITY RECALIBRATION 
# # calculate sequencing behavior that nees calibration: 
echo 'base quality recalibration - creating recal table'
java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -nct $THREADS \
  -U ALLOW_N_CIGAR_READS \
  -R $REF_FASTA \
  -I $INPUT_PREFIX.dups.sorted.indelrealigned.bam \
  -o $INPUT_PREFIX.dups.sorted.indelrealigned.BQSR.recal.table \
  -knownSites /stor/home/russd/work/references/GRCh38/dbsnp_146.hg38.vcf.gz \
  -knownSites /stor/home/russd/work/references/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -knownSites /stor/home/russd/work/references/GRCh38/hapmap_3.3_grch38_pop_stratified_af.vcf.gz


# # export reads with new calibration 
echo 'base quality recalibration - performing recalibration'
java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
  -T PrintReads \
  -nct $THREADS \
  -R $REF_FASTA \
  -I $INPUT_PREFIX.dups.sorted.indelrealigned.bam \
  -BQSR $INPUT_PREFIX.dups.sorted.indelrealigned.BQSR.recal.table \
  -o $INPUT_PREFIX.dups.sorted.indelrealigned.BQSR.bam 

  # -U ALLOW_N_CIGAR_READS \














