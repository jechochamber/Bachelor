import numpy as np
import pandas as pd
import ViennaRNA as vrna
import polygraph.input, polygraph.sequence, polygraph.visualize

seqs=polygraph.input.read_seqs('/fast/AG_Ohler/jdemoli/bachelorgit/data/full_dataset.txt',incl_ids=True)
singleseqs= seqs.loc[seqs['Group']=='Random_gc50']

secstructs=[]
mfes=[]

for seq in singleseqs["Sequence"]:
    (ss,mfe)=vrna.fold(seq)
    secstructs.append(ss)
    mfes.append(mfe)

singleseqs.loc[:,'mfe']=mfes
singleseqs.loc[:,'ss']=secstructs

singleseqs.to_csv('/fast/AG_Ohler/jdemoli/bachelorgit/data/Randomgc50_mfestats.csv')