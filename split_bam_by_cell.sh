#split bam file into sams by cell barcode 

#samtools view -F 1024 possorted_genome_bam.splitncigar.BQSR.bam | cut -f 12- | tr "\t" "\n"  | grep  "^CB:Z:"  | cut -d ':' -f 3 | sort | uniq > cell_list.txt 
# instead get filtered list of cell barcodes from cellranger out 


cat cell_list.txt | while read CELL_BARCODE; do samtools view -F 1024 -h possorted_genome_bam.splitncigar.BQSR.bam |  awk -v tag="CB:Z:$CELL_BARCODE" '($0 ~ /^@/ || index($0,tag)>0)' > split_by_cell/${CELL_BARCODE}.sam ; done

