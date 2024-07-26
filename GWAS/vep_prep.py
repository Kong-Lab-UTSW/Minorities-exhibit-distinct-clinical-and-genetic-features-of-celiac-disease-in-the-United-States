import pandas as pd
import sys

variants=sys.argv[1]

vars= pd.read_csv(variants)

vars[['#CHROM','POS']]=vars['locus'].str.split(':',expand=True)
vars[['REF','ALT']]=vars['alleles'].str.replace('[','').str.replace(']','').str.replace("'","").str.split(', ',expand=True)
vars['ID']=vars['locus']+vars['alleles']
vars_anno=vars[['#CHROM', 'POS','ID','beta','standard_error','z_stat','p_value','fit']]
vars_anno.to_csv(variants.replace('csv','tsv'),index=False,sep='\t')


