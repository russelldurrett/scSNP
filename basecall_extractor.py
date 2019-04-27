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
		pileupcolumn.get_query_sequences() # gets wrong query sequences! off by one, must be zero-base error. ugh. 
	except AssertionError: 
		print('Skipping', pileupcolumn.reference_name, pileupcolumn.pos, 'because of pileupcolumn.get_query_sequences AssertionError problem')
	else: 
		if pileupcolumn.pos==int(query_position[1]):
			bd = od()
			bd['chr']=pileupcolumn.reference_name
			bd['pos']=pileupcolumn.pos 
			approximate_depth=pileupcolumn.n 
			print('Base', bd['chr'], bd['pos'], '\t', 'Approximate Depth: ', str(approximate_depth)) 
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
					rd['base_call'] = pileupalignment.query_sequence[pileupread.query_position-1]
					df=df.append(rd, ignore_index=True)
			# print(len(df), 'base calls written thus far')

print(df.head(10))
if len(df)>0: 
	df.to_csv('basecalls_by_cell_umi.{}_{}.tsv'.format(query_position[0],query_position[1]), sep='\t', index=False)	
	print(df.base_call.value_counts()) 
samfile.close()
