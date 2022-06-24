"""
script for replacing CellRanger gene names with named orthologs.

Inputs:
    CellRanger output directory: "filtered_feature_bc_matrix/"
    Ortholog inference output: TSV file with columns named by species identifier. Each row contains one orthogroup, cell values are comma-separated lists of NCBI gene names (or NA if no orthologs identified). (parse_orthout.py output)

Output:
    Modified CellRanger output directory with features.tsv.gz, matrix.mtx.gz, and barcodes.tsv.gz
"""

import pandas as pd, numpy as np
from scipy import sparse
from scipy.io import mmread, mmwrite
import argparse, os

parser=argparse.ArgumentParser()
parser.add_argument('--cellranger', help='Path to CellRanger output directory (filtered_feature_bc_matrix/)')
#parser.add_argument('--features',help='Path to CellRanger output features file')
#parser.add_argument('--matrix', help='Path to CellRanger output matrix file')
parser.add_argument('--orthologs',help='Path to ortholog inference output')
parser.add_argument('--sp',help='Identifier for source species (the one the scRNAseq data came from)')
parser.add_argument('--st',help='Identifier for target species (the one whose gene names you want to use)')
parser.add_argument('--keep_name',action='store_true',help='Keep original gene names (rather than renaming to ortholog)?')
parser.add_argument('--strict',default=1,type=float,choices=[0,1,2],help=
                    '''How strict to be when choosing what to rename, drop.
                        0 = rename 1:1 orthologs, keep 1:many and 1:0 orthologs
                        1 = rename 1:1 orthologs, drop 1:many orthologs, keep 1:0 orthologs
                        2 = rename 1:1 orthologs, drop 1:many orthologs, drop 1:0 orthologs''')
parser.add_argument('--out',help='directory to output new "features.tsv.gz" and "matrix.mtx.gz" files to')
args=parser.parse_args()

# read in orthologs data and do some formatting
orth=pd.read_csv(args.orthologs, sep='\t')

# get column identifiers based on species identifiers
spcol=None
stcol=None
for c in orth.columns:
    if args.sp in c:
        spcol=c
    elif args.st in c:
        stcol=c
if not (spcol and stcol):
    sys.exit('Species identifier not found! Exiting...')

# columns as lists of strings, not just a long weird string
orth.loc[~pd.isna(orth[spcol]),spcol]=orth.loc[~pd.isna(orth[spcol]),spcol].apply(eval)
orth.loc[~pd.isna(orth[stcol]),stcol]=orth.loc[~pd.isna(orth[stcol]),stcol].apply(eval)
orth.loc[pd.isna(orth[spcol]),spcol] = [np.array(np.nan) for i in range(sum(pd.isna(orth[spcol])))]
orth.loc[pd.isna(orth[stcol]),stcol] = [np.array(np.nan) for i in range(sum(pd.isna(orth[stcol])))]

# read in cellranger data
feat=pd.read_csv(os.path.join(args.CellRanger, 'features.tsv.gz'), sep='\t', header=None)
matrix = mmread(os.path.join(args.CellRanger, 'matrix.mtx.gz'))
matrix = matrix.toarray()

### replace feature gene IDs with orthologs:

sp_orth=np.array(orth[spcol]) # source species orthologs as array
st_orth=np.array(orth[stcol]) # target species orthologs as array
ind_remove=[]                 # indices in features dataframe to drop

# loop thru each feature (slow, ugh)
for f in range(len(feat)):
    # find source species feature in ortholog dataframe
    spfeature = feat.loc[f,1]
    vec_func = np.vectorize(lambda arr: spfeature in arr)
    index = np.where(vec_func(sp_orth))[0]

    # if it has an ortholog
    if len(index) > 0:

        # get its orthogroup
        sporth=np.array(sp_orth[index[0]])
        storth=np.array(st_orth[index[0]])

        # if 1:1...
        if sporth.size == 1 and storth.size == 1:
            # (deal with 1:0 orthologs)
            if pd.isna(storth):
                if args.strict == 2:
                    ind_remove.append(f)
            # otherwise rename as target species ortholog
            else:
                if args.keep_name:
                    pass
                else:
                    feat.loc[f,0] = storth[0]
                    feat.loc[f,1] = storth[0]

    # (deal with 1:many, many:many orthologs)
        elif args.strict > 0:
            ind_remove.append(f)
    elif args.strict == 2:
        ind_remove.append(f)

# remove features to drop
feat.drop(ind_remove, inplace=True)
feat = feat.reset_index(drop=True)

# remove dropped features from matrix
mask = np.ones(matrix.shape[0], bool)
mask[ind_remove] = 0
matrix = matrix[mask]

# save data
outdir = args.out
if outdir[-1] != '/':
    outdir = outdir+'/'
feat.to_csv(os.path.join(args.out, 'features.tsv'), sep='\t', index=False, header=False)
matrix=sparse.csr_matrix(matrix)
mmwrite(os.path.join(args.out, 'matrix.mtx'), matrix)

# formatting
os.system('cp '+os.path.join(args.CellRanger, 'barcodes.tsv.gz')+' '+os.path.join(args.out, 'barcodes.tsv.gz'))
os.system('gzip '+os.path.join(args.out, 'features.tsv'))
os.system('gzip '+os.path.join(args.out, 'matrix.mtx'))
