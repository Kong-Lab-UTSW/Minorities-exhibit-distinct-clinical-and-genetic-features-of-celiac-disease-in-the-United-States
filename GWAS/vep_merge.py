import pandas as pd
import sys
import io

#open vcf as dataframe
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'chrom' })

counts=sys.argv[2]
count = read_vcf(counts)
#count['ID']=count['locus']+count['alleles']
anno=sys.argv[1]
gene2=pd.read_csv(anno,sep='\t')
total=pd.merge(gene2,count,left_on='ID',right_on='#Uploaded_variation')
total.to_csv(anno.replace('.tsv','_full.vep.csv'),index=False)
