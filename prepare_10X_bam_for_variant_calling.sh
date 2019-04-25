# 10X output to GATK-ready 

# mark duplicates
java -jar ~/software/picard.jar MarkDuplicates INPUT=possorted_genome_bam.bam OUTPUT=possorted_genome_bam.dups.bam M=dup.metrics     
    # java -jar picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics 
java -jar ~/software/picard.jar UmiAwareMarkDuplicatesWithMateCigar INPUT=possorted_genome_bam.bam OUTPUT=possorted_genome_bam.dupsumi.bam M=dupumi.metrics  UMI_METRICS_FILE=dupumi.umi.metrics BARCODE_TAG=BC UMI_TAG_NAME=UR

# marking duplicates with MarkDuplicates (not UMI/BC aware): 255963577 records, Marking 183084121 records as duplicates. 71.5% duplicate records 


# sorted?

# split n trim 
java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R ~/work/references/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -I possorted_genome_bam.bam -o possorted_genome_bam.splitncigar.bam -U ALLOW_N_CIGAR_READS -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60

# indel realignment target creation 
#java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R ~/work/references/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -I possorted_genome_bam.bam -o possorted_genome_bam.splitncigar.bam -U ALLOW_N_CIGAR_READS -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60
# indel realignment 
#java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R ~/work/references/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -I possorted_genome_bam.bam -o possorted_genome_bam.splitncigar.bam -U ALLOW_N_CIGAR_READS -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60

# base recalibration 
# PrintReads does on the fly base assesment and then calibration (?) 
java -jar GenomeAnalysisTK.jar -T PrintReads -R reference.fasta -I input.bam -BQSR recalibration_report.grp -o output.bam




#### Analysis-Ready RNAseq Reads 

# Variant Calling = HC in RNAseq Mode 
java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/work/references/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -I possorted_genome_bam.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o output.vcf
# -nct <threads> 
--maxNumHaplotypesInPopulation <128>
--sample_ploidy <2>


# Variant Calling = bcftools call 

 ~/software/bcftools-1.9/bcftools mpileup -Ou --max-depth 10000 -f ~/work/references/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa possorted_genome_bam.bam | ~/software/bcftools-1.9/bcftools call --threads 6  --multiallelic-caller --variants-only --skip-variants indels -o bcfcalls.vcf 
 # -Ob to .bcf        

# add AF and other tags to VCF 
export BCFTOOLS_PLUGINS=/home/russd/software/bcftools-1.9/plugins/
~/software/bcftools-1.9/bcftools +fill-tags tmp.vcf    # just the AF    bcftools +fill-tags test.vcf  -- -t AF

~/software/bcftools-1.9/bcftools +fill-tags tmp.vcf | ~/software/bcftools-1.9/bcftools query -i'QUAL>20 && DP>10' -f'%CHROM %POS %QUAL %DP\n' -





# Raw Variant Filtering (separate by variant type) 
# java -jar GenomeAnalysisTK.jar -T VariantFiltration -R hg_19.fasta -V input.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o output.vcf 
	# -window 35 -cluster 3 = removes all snps if more than 3 occur in 35bp window 








