#!/usr/bin/python 

# script to concatenate single basecall files for each of the five samples 

import sys 
import os 
import pandas as pd 

def is_non_zero_file(fpath):  
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 1 else False


df = pd.DataFrame()
samples = ['TP0','FM1','FM7','FMV1','FMV7']
for sample in samples: 
	os.chdir('/stor/home/russd/scratch/10x_snps/{}/basecalls'.format(sample))
	basecall_files=os.listdir('/Users/red/Desktop/10x_snps/ex_basecalls')
	for basecall_file in basecall_files: 
		if is_non_zero_file(basecall_file):
			print('working on ', basecall_file)
			dftmp = pd.read_table(basecall_file)
			dftmp['sample']=sample
			# print(dftmp.head(2))
			df=pd.concat([df,dftmp], axis=0, ignore_index=True)
			df.to_csv("/stor/home/russd/scratch/10x_snps/concatenated_basecalls.working...tsv", sep='\t', index=False)
		else: 
			print('skipping empty file ', basecall_file)

print(df.head(5))
df.to_csv("/stor/home/russd/scratch/10x_snps/concatenated_basecalls.tsv", sep='\t', index=False)

