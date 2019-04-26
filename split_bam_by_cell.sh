#split bam file into sams by cell barcode 

#samtools view -F 1024 possorted_genome_bam.splitncigar.BQSR.bam | cut -f 12- | tr "\t" "\n"  | grep  "^CB:Z:"  | cut -d ':' -f 3 | sort | uniq > cell_list.txt 
# instead get filtered list of cell barcodes from cellranger out 

# too slow: (loops through bam file over and over filtering out a specific cell)
#cat cell_list.txt | while read CELL_BARCODE; do samtools view -F 1024 -h possorted_genome_bam.splitncigar.BQSR.bam |  awk -v tag="CB:Z:$CELL_BARCODE" '($0 ~ /^@/ || index($0,tag)>0)' > split_by_cell/${CELL_BARCODE}.sam ; done

# MUCH faster: (awks the cell barcode identity per line and prints to respective file. loops thorugh bam only once)
# still fucked up though... 
# samtools view ../possorted_genome_bam.splitncigar.BQSR.bam | head | \
# while read line; do
#         cell=$(awk '{split($12,cbs,":"); print cbs[3]}' <<< $line);
#   echo $line > $cell
# done

# faster and actually works! 
samtools view possorted_genome_bam.splitncigar.BQSR.bam | \
while read line; do
  cell=$(awk '{ for(i=1;i<=NF;++i){ if ($i ~ /CB:Z:/) {cell=split($i,cbs,":"); print cbs[3] } } }' <<< $line);
  if [ -z "$var" ] 
  then 
  	echo "Cell is not set! Cell: ", $cell, " |  ", $line
  else 
  	echo $line >> split_by_cell/$cell
  fi 
done



