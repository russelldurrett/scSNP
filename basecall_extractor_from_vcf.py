import collections 
od = collections.OrderedDict
import pysam
import pandas as pd 
import sys 

bamfile_handle = sys.argv[1]
vcf_handle = sys.argv[2]

print('Querying {} off of positions in VCF {}'.format(bamfile_handle, vcf_handle))


def afs(array_of_base_calls): 
	array_of_base_calls = [g.upper() for g in array_of_base_calls]
	d = od([(i,array_of_base_calls.count(i)/float(len(array_of_base_calls))) for i in set(array_of_base_calls)])
	return d 

def maf(array_of_base_calls): 
	afs_list = afs(array_of_base_calls).values()
	if len(afs_list)==1:
		return 0.0
	else:
		return sorted(afs_list)[-2]

bases = ['','G', 'A','T','C']


df = pd.DataFrame()

samfile = pysam.AlignmentFile(bamfile_handle, "rb" )
vcffile = open(vcf_handle, 'rb') 

for line in vcffile.readlines(): 
	line=str(line.decode('utf-8')).rstrip()
	if line[0]=='#': 
		line_is_comment=True
	else: 
		vcf_line=line.split('\t')
		query_chrom=vcf_line[0]
		query_position=vcf_line[1]
		query_ref_allele=vcf_line[3]
		for pileupcolumn in samfile.pileup(str(query_chrom), int(query_position)-1, int(query_position)+1):
			try: 
				pileupcolumn.get_query_sequences()
			except AssertionError: 
				do_nothing=True
				#print('Skipping', pileupcolumn.reference_name, pileupcolumn.pos, 'because of pileupcolumn.get_query_sequences AssertionError problem')
			else: 
				if pileupcolumn.pos==int(query_position):
					bd = od()
					bd['chr']=pileupcolumn.reference_name
					bd['pos']=pileupcolumn.pos 
					bd['ref']=query_ref_allele
					query_sequences = [g.upper() for g in pileupcolumn.get_query_sequences()]
					query_base_counts = [query_sequences.count(b) for b in set(bases+list(set(query_sequences)))]
					print('extracing calls on base', bd['chr'], bd['pos'], '\t', '\t'.join([str(b)+":"+str(c) for b,c in zip(bases,query_base_counts)]))
					for pileupread in pileupcolumn.pileups: 
						if pileupread.query_position!=None:
							rd = bd.copy()
							pileupalignment = pileupread.alignment
							pileupalignmenttags = od(pileupalignment.tags)
							if 'CB' in pileupalignmenttags.keys():
								rd['cell'] = pileupalignmenttags['CB'].split('-')[0]
							else:
								rd['cell'] = pileupalignmenttags['CR'].split('-')[0]
							if 'UB' in pileupalignmenttags.keys():
								rd['umi'] = pileupalignmenttags['UB'] 
							else:
								rd['umi'] = pileupalignmenttags['UR']
							rd['base_call'] = pileupalignment.query_sequence[pileupread.query_position]
							df=df.append(rd, ignore_index=True)
					print(len(df), 'base calls written thus far')
					df.to_csv('basecalls_by_cell_umi.working.tsv', sep='\t', index=False)

df.to_csv('basecalls_by_cell_umi.final.tsv', sep='\t', index=False)
samfile.close()
