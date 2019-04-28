#!/usr/bin/python 

# script to process concatenated table 

import pandas as pd 
import scipy.stats as stats 

df = pd.read_table("/stor/home/russd/scratch/10x_snps/variants/concatenated_basecalls.dirty.tsv")
df = df[['chr', 'pos', 'sample', 'cell', 'umi', 'base_call']].dropna()

print('starting number of basecalls: ',len(df))
print(df.head())

# 14,644,537 basecalls 


print('removing unfiltered cells')
filtered_cells = pd.read_table('10x_filtered.cell_barcode_sample.tsv', sep='\t', names=['sample', 'cell']).apply(lambda d: '-'.join(d),1)
df = df[df[['sample','cell']].apply(lambda d: '-'.join(d),1).isin(filtered_cells)].reset_index(drop=True)
print('basecalls in filtered cells: ', len(df))
print(df.head())

# 8,987,536 basecalls


print('collapsing umis')
df=df.groupby(['chr', 'pos', 'sample', 'cell', 'umi'])['base_call'].apply(lambda d: stats.mode(d)[0][0])
print('after umi collapsing number of basecalls: ',len(df))
print(df.head())


# 7,209,889 basecalls 



df=df.reset_index()
df=df.astype({'pos':float, 'chr':str}).astype({'pos':int})
df['chr:pos']=df.apply(lambda r: "{}:{}".format(str(r['chr']).rstrip(),str(r['pos']).rstrip()), axis=1)
df.drop(['chr', 'pos'], 1, inplace=True)
print(df.head())




print('collapsing base_calls to counts per base')
df_counts=df.groupby(['chr:pos', 'cell', 'sample'])['base_call'].value_counts()
df_counts=df_counts.rename_axis(['chr:pos', 'cell', 'sample', 'base_call_gt']).astype(int).reset_index()
print(df_counts.head())




cell_counts=df_counts.pivot_table(index=['chr:pos', 'cell', 'sample'], columns='base_call_gt', values='base_call', aggfunc='sum').drop('N',1).fillna(0).astype(int)
sample_counts=df_counts.pivot_table(index=['chr:pos', 'sample'], columns='base_call_gt', values='base_call', aggfunc='sum').drop('N',1).fillna(0).astype(int)

print('saving intermediate files')
print('cell_counts:')
print(cell_counts.head())
cell_counts.to_csv('/stor/home/russd/scratch/10x_snps/concatenated_basecalls.cell_counts.intermediate.tsv', sep='\t', index=True)
print('sample_counts:')
print(sample_counts.head())
sample_counts.to_csv('/stor/home/russd/scratch/10x_snps/concatenated_basecalls.sample_counts.intermediate.tsv', sep='\t', index=True)





print('calculating outliers')
# get rid of outliers (seq error)
# if base is never sequenced more than 10 times in any particular sample, call it outlier 
base_call_minimum_nonoutliers = sample_counts.stack()[sample_counts.stack()>10].reset_index()[['chr:pos','base_call_gt']].drop_duplicates()
# if a base is never sequenced more than 5% of maximum base in any particular sample, call it outlier 
base_proportion_to_max_per_sample = sample_counts.astype(float).div(sample_counts.max(1).astype(float),axis=0).stack() 
base_proportion_to_max_nonoutliers = base_proportion_to_max_per_sample[base_proportion_to_max_per_sample>0.05].reset_index()[['chr:pos','base_call_gt']].drop_duplicates()

nonoutliers_passing_filters = pd.concat([base_call_minimum_nonoutliers, base_proportion_to_max_nonoutliers]) 
nonoutliers_passing_filters = nonoutliers_passing_filters[nonoutliers_passing_filters.duplicated()].reset_index().set_index(['chr:pos','base_call_gt'])

sample_counts = sample_counts.stack().reset_index().set_index(['chr:pos', 'base_call_gt']).loc[nonoutliers_passing_filters.index].reset_index().set_index(['chr:pos','sample', 'base_call_gt']).unstack().fillna(0).astype(int) 
sample_counts.columns = sample_counts.columns.droplevel()
cell_counts = cell_counts.stack().reset_index().set_index(['chr:pos', 'base_call_gt']).loc[nonoutliers_passing_filters.index].reset_index().set_index(['chr:pos','sample', 'cell', 'base_call_gt']).unstack().fillna(0).astype(int) 
cell_counts.columns = cell_counts.columns.droplevel()

print('done filtering outliers, writing to files')
print('cell_counts:')
print(cell_counts.head())
cell_counts.reset_index().to_csv('/stor/home/russd/scratch/10x_snps/concatenated_basecalls.cell_counts.final.tsv', sep='\t', index=False)
print('sample_counts:')
print(sample_counts.head())
sample_counts.reset_index().to_csv('/stor/home/russd/scratch/10x_snps/concatenated_basecalls.sample_counts.final.tsv', sep='\t', index=False)



# calculate reference base = we'll use maximum base in TP0 cells 
sample_max = sample_counts.idxmax(1) 
posrefs = sample_max.reset_index()[sample_max.reset_index()['sample']=='TP0'].reset_index(drop=True).drop('sample',1).set_index('chr:pos')
posrefs.columns=['base_call_gt']
posrefs[0]=0
posrefs = posrefs.reset_index().set_index(['chr:pos','base_call_gt'])


cell_mutants = cell_counts.copy().stack().reset_index().set_index(['chr:pos','base_call_gt'])
cell_mutants.update(posrefs)
cell_mutants = cell_mutants.reset_index().set_index(['chr:pos','sample','cell','base_call_gt']).unstack().astype(int)

posrefs.head()
cell_mutants.head()


sample_dict={'TP0':'1','FM1':'2','FM7':'3','FMV1':'4','FMV7':'5'}
cell_mutants = cell_mutants.reset_index()
cell_mutants['tmp'] = cell_mutants['sample'].apply(lambda d: sample_dict[d])
cell_mutants['tmp2']=cell_mutants['cell']+'-'+cell_mutants['tmp']
cell_mutants.drop(['tmp', 'cell', 'sample'],1,inplace=True)
cell_mutants=cell_mutants.rename(columns={'tmp2':'cell'}).set_index(['chr:pos','cell'])

cell_mutant_matrix = cell_mutants.sum(1).astype(bool).unstack().fillna(False).astype(int)   




### plotting data 

# # number of reads per cell/position 
# cell_counts.sum(1).value_counts()  

# # number of genotypes for each cell/position 
# cell_counts.astype(bool).sum(1).value_counts()   












