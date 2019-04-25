# call variants on high-quality bam 

THREADS=4
REF_FASTA=~/work/references/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa


#### Analysis-Ready RNAseq Reads 

# Variant Calling = HC in RNAseq Mode  - too stringent.... 
# java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/work/references/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -I possorted_genome_bam.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o output.vcf
# # -nct <threads> 
# --maxNumHaplotypesInPopulation <128>
# --sample_ploidy <2>


# Variant Calling = using bcftools call 

# ensure only use non-duplicate reads, skip indel calling, filter @ min 1K depth, min 1% reference and alternate reads (ensuring base heterogeneity) 
samtools view -bF 1024 possorted_genome_bam.splitncigar.BQSR.bam | ~/software/bcftools-1.9/bcftools mpileup -Ou --max-depth 20000 -f $REF_FASTA - | ~/software/bcftools-1.9/bcftools call --threads $THREADS  --multiallelic-caller --skip-variants indels | ~/software/bcftools-1.9/bcftools filter  -i'DP>1000 & ((DP4[2]+DP4[3])/(DP[0]+DP[1]+DP[2]+DP[3]))>0.01 & ((DP4[0]+DP4[1])/(DP[0]+DP[1]+DP[2]+DP[3]))>0.01' > bcfcalls.vcf 

# Variant Filter - window of 35 bases containing 3+ SNPs (probably not great quality) 
# only works with variant-only VCF... 
# java -jar ~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_FASTA -V bcfcalls.vcf -window 35 -cluster 3 -o bcfcalls.window35clust3.vcf 

# larger filter example 
# bcftools filter -sLowQual -g3 -G10  -e'%QUAL<10 || (RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || %MAX(DV)<=3 || %MAX(DV)/%MAX(DP)<=0.3' calls.vcf.gz


# # use bcftools plugin to add AF and other tags to VCF (AF only works for multisample vcf...)
# export BCFTOOLS_PLUGINS=/home/russd/software/bcftools-1.9/plugins/
# ~/software/bcftools-1.9/bcftools +fill-tags tmp.vcf    # just the AF    bcftools +fill-tags test.vcf  -- -t AF


# FORMAT VCF TO ALLELE READ DEPTH DF 
# filter 1K read depth, at least 1% high quality variant reads 
# CHROM, POS, TOTAL_DEPTH, HQ_REF_DEPTH, HQ_ALT_DEPTH 
~/software/bcftools-1.9/bcftools query -f'%CHROM %POS %DP %DP4\n' bcfcalls.vcf | tr ',' ' ' | awk '{print $1,$2,$3,$4+$5,$6+$7}' > bcfcalls.txt






