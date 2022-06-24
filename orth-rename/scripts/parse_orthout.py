"""
parse output from ortholog inference
(Broccoli, GeneSpace, OrthoFinder),
output orthologs from species of interest into parsable dataframe
"""

import pandas as pd, numpy as np, argparse, sys

parser=argparse.ArgumentParser()
parser.add_argument('--orthologs', help='''path to orthology inference output:
                                           Broccoli: dir_step3/table_OGs_protein_names.txt
                                           GeneSpace: results/gffWithOgs.txt.gz
                                           Orthofinder: Orthogroups/Orthogroups.tsv''')
parser.add_argument('--method', choices=['Broccoli','OrthoFinder','GeneSpace'],
                    help='orthology inference method')
parser.add_argument('--s1', help='species 1 identifier')
parser.add_argument('--s2', help='species 2 identifier')
parser.add_argument('--out', help='outpath')
args=parser.parse_args()

def parse_genespace(df, s1, s2):
    # get rows from relevant species as seperate dataframes
    df1 = df[df.genome == s1]
    df2 = df[df.genome == s2]
    df1 = df1[['id', 'og']].reset_index(drop=True)
    df2 = df2[['id', 'og']].reset_index(drop=True)
    df1.rename(columns={'id':'s1'}, inplace=True)
    df2.rename(columns={'id':'s2'}, inplace=True)
    
    # collapse orthogroups to one row each
    df1 = df1.groupby('og')['s1'].apply(list).reset_index(name='s1')
    df2 = df2.groupby('og')['s2'].apply(list).reset_index(name='s2')

    # merge on common orthogroups
    df = pd.merge(df1, df2, on='og', how='left')

    return df    

def parse_broccoli(df, s1, s2):
    # get columns from relevant species as seperate dataframes
    df1 = df[['#OG_name', s1]]
    df2 = df[['#OG_name', s2]]

    # merge on common orthogroups
    df = pd.merge(df1, df2, on='#OG_name')

    # format cell values as lists
    df.loc[~pd.isna(df[s1]), s1] = np.asarray([i.split(sep=' ') for i in df.loc[~pd.isna(df[s1]), s1]], dtype='object')
    df.loc[~pd.isna(df[s2]), s2] = np.asarray([i.split(sep=' ') for i in df.loc[~pd.isna(df[s2]), s2]], dtype='object')

    return df

def parse_orthofinder(df, s1, s2):
    # get columns from relevant species as separate dataframes
    df1 = df[['Orthogroup', s1]]
    df2 = df[['Orthogroup', s2]]

    # merge on common orthogroups
    df = pd.merge(df1, df2, on='Orthogroup')

    # format cell values as lists
    df.loc[~pd.isna(df[s1]), s1] = np.asarray([i.split(sep=',') for i in df.loc[~pd.isna(df[s1]), s1]], dtype='object')
    df.loc[~pd.isna(df[s2]), s2] = np.asarray([i.split(sep=',') for i in df.loc[~pd.isna(df[s2]), s2]], dtype='object')

    return df

# read orthogroup dataframe
og = pd.read_csv(args.orthologs, sep='\t')

# find relevant column IDs
s1_ID = None
s2_ID = None
for c in og.columns:
    if args.s1 in c:
        s1_ID = c
    elif args.s2 in c:
        s2_ID = c
if not (s1_ID and s2_ID):
    sys.exit('Species identifier not found! Exiting...')

# format dataframe
if args.method=='Broccoli':
    ogdf = parse_broccoli(og, s1_ID, s2_ID)
    ogdf.rename(columns={s1_ID:args.s1, s2_ID:args.s2}, inplace=True)
elif args.method=='GeneSpace':
    ogdf = parse_genespace(og, s1_ID, s2_ID)
    ogdf.rename(columns={'s1':args.s1, 's2':args.s2}, inplace=True)
elif args.method=='OrthoFinder':
    ogdf = parse_orthofinder(og, s1_ID, s2_ID)
    ogdf.rename(columns={s1_ID:args.s1, s2_ID:args.s2}, inplace=True)

og.to_csv(args.out, sep='\t')
