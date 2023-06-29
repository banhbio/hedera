#!/usr/bin/env python
# coding: utf-8

import sys
import os
import pandas as pd
import glob

#input directory
path = "/aptmp/ideas2/ivyyue/tara_ocean_polar_protists_PRJEB9691/PRJEB9691/28_0.8_samples/coassembly/binning/delineage/"
tetra_files = glob.glob(os.path.join(path,"*tetra_count.tsv"))

for f in tetra_files:
    basename = f.split(".")[0]+"."+f.split(".")[1]+"."+f.split(".")[2] #for the case that bin fa files has two dots in filename
    df = pd.read_csv(f, sep='\t')
    cols = df.columns[df.columns.str.contains('A|T|C|G')]
    df[cols]  = df[cols].div(df[cols].sum(axis=1), axis=0)
    df.to_csv(f'{basename}_normalized.csv', sep=',', index=False, header=False)