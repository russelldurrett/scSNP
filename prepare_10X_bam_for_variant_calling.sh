# 10X output to GATK-ready 

THREADS=12
REF_FASTA=~/work/references/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa



# 10X output has some duplicates already marked.... 
# MARK DUPLICATES
#java -jar ~/software/picard.jar MarkDuplicates INPUT=possorted_genome_bam.bam OUTPUT=possorted_genome_bam.dups.bam M=dup.metrics   # CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics 
# mark duplicates using UMI info: 
#java -jar ~/software/picard.jar UmiAwareMarkDuplicatesWithMateCigar INPUT=possorted_genome_bam.bam OUTPUT=possorted_genome_bam.dupsumi.bam M=dupumi.metrics  UMI_METRICS_FILE=dupumi.umi.metrics BARCODE_TAG=BC UMI_TAG_NAME=UR

# DUP MARKING OUTPUT COMPARISON 
# marking duplicates with MarkDuplicates (not UMI/BC aware): 255963577 records, Marking 183084121 records as duplicates. 71.5% duplicate records 


# sorted already 



# split N' trim 
echo 'splitting and trimming N cigar strings'
java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
  -T SplitNCigarReads \
  -R $REF_FASTA  \
  -I possorted_genome_bam.bam \
  -o possorted_genome_bam.splitncigar.bam \
  -U ALLOW_N_CIGAR_READS \
  -rf ReassignOneMappingQuality \
  -RMQF 255 -RMQT 60




# INDEL REALIGNMENT (skippable, for now - focusing on SNPs)
# indel realignment target creation 
# indel realignment 

echo 'realigning around indels - creating targets '
java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  -R $REF_FASTA \
  -known /stor/home/russd/work/references/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -I possorted_genome_bam.splitncigar.bam \
  -o possorted_genome_bam.splitncigar.intervals

echo 'realigning around indels - performing realignment'
java -Xmx8G -Djava.io.tmpdir=/tmp -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R $REF_FASTA \
    -targetIntervals possorted_genome_bam.splitncigar.intervals \
    -known /stor/home/russd/work/references/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -I possorted_genome_bam.splitncigar.bam \
    -o possorted_genome_bam.splitncigar.indelrealigned.bam



# BASE QUALITY RECALIBRATION 
# calculate sequencing behavior that nees calibration: 
echo 'base quality recalibration - creating recal table'
java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -nct $THREADS \
  -U ALLOW_N_CIGAR_READS \
  -R $REF_FASTA \
  -I possorted_genome_bam.splitncigar.indelrealigned.bam \
  -o BQSR.recal.table \
  -knownSites /stor/home/russd/work/references/GRCh38/dbsnp_146.hg38.vcf.gz \
  -knownSites /stor/home/russd/work/references/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -knownSites /stor/home/russd/work/references/GRCh38/hapmap_3.3_grch38_pop_stratified_af.vcf.gz


# export reads with new calibration 
echo 'base quality recalibration - performing recalibration'
java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
  -T PrintReads \
  -nct $THREADS \
  -U ALLOW_N_CIGAR_READS \
  -R $REF_FASTA \
  -I possorted_genome_bam.splitncigar.indelrealigned.bam \
  -BQSR BQSR.recal.table \
  -o possorted_genome_bam.splitncigar.indelrealigned.BQSR.bam















