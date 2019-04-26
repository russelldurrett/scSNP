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
# over 600K cell barcodes found though, need to filter only those that pass cellranger filters 

for filtered_cell in `cat cell_list.txt` 
do 
	touch split_by_cell/$filtered_cell 
done 


samtools view possorted_genome_bam.splitncigar.BQSR.bam | \
while read line; do
  cell=$(awk '{ for(i=1;i<=NF;++i){ if ($i ~ /CB:Z:/) {cell=split($i,cbs,":"); print cbs[3] } } }' <<< $line);
  if [ -z "$cell" ] 
  then 
  	do_nothing='cell is empty for some reason'
  else 
  	if [ -f split_by_cell/$cell ]
  	then 
  		echo $line >> split_by_cell/$cell
  	else 
  		do_nothing='cell barcode was not in filtered list from cellranger'
  	fi
  fi 
done
