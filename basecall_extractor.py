import collections 
od = collections.OrderedDict
import pysam
import pandas as pd 
import sys 

bamfile_handle = sys.argv[1]
query_position = sys.argv[2].split(':')

print('Querying {}, position {}'.format(bamfile_handle, query_position))


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

samfile = pysam.AlignmentFile(bamfile_handle, "rb" )


df = pd.DataFrame()

for pileupcolumn in samfile.pileup(str(query_position[0]), int(query_position[1])-1, int(query_position[1])+1):
	try: 
		pileupcolumn.get_query_sequences()
	except AssertionError: 
		print('Skipping', pileupcolumn.reference_name, pileupcolumn.pos, 'because of pileupcolumn.get_query_sequences AssertionError problem')
	else: 
		if pileupcolumn.pos==int(query_position[1]):
			bd = od()
			bd['chr']=pileupcolumn.reference_name
			bd['pos']=pileupcolumn.pos 
			query_sequences = [g.upper() for g in pileupcolumn.get_query_sequences()]
			query_base_counts = [query_sequences.count(b) for b in set(bases+list(set(query_sequences)))]
			print('Base', bd['chr'], bd['pos'], '\t', '\t'.join([str(b)+":"+str(c) for b,c in zip(bases,query_base_counts)]))
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
			# print(len(df), 'base calls written thus far')

# df.to_csv('mtDNA_basecalls.final.tsv', sep='\t', index=False)
print(df) 
samfile.close()
