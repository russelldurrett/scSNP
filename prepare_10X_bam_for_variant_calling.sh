# 10X output to GATK-ready 

THREADS=6
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
java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R $REF_FASTA  -I possorted_genome_bam.bam -o possorted_genome_bam.splitncigar.bam -U ALLOW_N_CIGAR_READS -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60




# INDEL REALIGNMENT (skippable, for now - focusing on SNPs)
# indel realignment target creation 
# indel realignment 



# BASE QUALITY RECALIBRATION 
# calculate sequencing behavior that nees calibration: 
java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T BaseRecalibrator -nct $THREADS -U ALLOW_N_CIGAR_READS -R $REF_FASTA -I possorted_genome_bam.splitncigar.bam -o BQSR.recal.table \
  -knownSites /stor/home/russd/work/references/GRCh38/dbsnp_146.hg38.vcf.gz \
  -knownSites /stor/home/russd/work/references/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -knownSites /stor/home/russd/work/references/GRCh38/hapmap_3.3_grch38_pop_stratified_af.vcf.gz


# export reads with new calibration 
java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T PrintReads -nct $THREADS -U ALLOW_N_CIGAR_READS -R $REF_FASTA -I possorted_genome_bam.splitncigar.bam -BQSR BQSR.recal.table -o possorted_genome_bam.splitncigar.BQSR.bam















