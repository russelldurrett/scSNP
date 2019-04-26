# call variants on high-quality exome bam 

SAMPLE=$1 

THREADS=6

BAM_DIR=/stor/home/russd/scratch/10x_snps/HG3_Exome/alignment

REF_FASTA=/stor/home/russd/work/references/GRCh38/Homo_sapiens_assembly38.fasta
REF_BWA=/stor/home/russd/work/references/GRCh38/bwa/Homo_sapiens_assembly38.fasta
PICARD_JAR=/home/russd/software/picard-tools-2.1.1/picard.jar
GATK_JAR=/home/russd/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar


# Variant Calling = using bcftools call 

# ensure only use non-duplicate reads, skip indel calling, call only targets suggested by 10x data 
samtools view -bF 1024 $BAM_DIR/$SAMPLE\_recal_reads.bam | ~/software/bcftools-1.9/bcftools mpileup -Ou --max-depth 5000 -f $REF_FASTA - | ~/software/bcftools-1.9/bcftools call --threads $THREADS  --multiallelic-caller --skip-variants indels -T /stor/home/russd/scratch/10x_snps/FM1-4.bcfcalls.targets.txt > $SAMPLE\_bcfcalls.vcf 

# # use bcftools plugin to add AF and other tags to VCF (AF only works for multisample vcf...)
# export BCFTOOLS_PLUGINS=/home/russd/software/bcftools-1.9/plugins/
# ~/software/bcftools-1.9/bcftools +fill-tags $SAMPLE\_bcfcalls.vcf > $SAMPLE\_bcfcalls.tagged.vcf    # just the AF    bcftools +fill-tags test.vcf  -- -t AF


# FORMAT VCF TO ALLELE READ DEPTH DF 
# filter 1K read depth, at least 1% high quality variant reads 
# CHROM, POS, TOTAL_DEPTH, HQ_REF_DEPTH, HQ_ALT_DEPTH 
#~/software/bcftools-1.9/bcftools query -f'%CHROM %POS %DP %DP4\n' bcfcalls.vcf | tr ',' ' ' | awk '{print $1,$2,$3,$4+$5,$6+$7}' > bcfcalls.txt






