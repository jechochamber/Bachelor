import numpy as np
import pandas as pd
import ViennaRNA as vrna
import polygraph.input, polygraph.sequence, polygraph.visualize

seqs=polygraph.input.read_seqs('/fast/AG_Ohler/jdemoli/bachelorgit/data/full_dataset.txt',incl_ids=True)
seqs['Sequence Length'] = seqs.Sequence.apply(len)
humanseqs=seqs.loc[seqs['Group']=='Human']

human160seqs=humanseqs[humanseqs["Sequence Length"] >= 160]

secstructs=[]
mfes=[]

for i in range(len(human160seqs['Sequence'])):
    seq = human160seqs.iloc[i, 0]
    human160seqs.iloc[i, 0] = seq[len(seq) - 160: len(seq)]
human160seqs.loc[:, 'Sequence Length'] = human160seqs['Sequence'].apply(len)

for seq in human160seqs["Sequence"]:
    (ss,mfe)=vrna.fold(seq)
    secstructs.append(ss)
    mfes.append(mfe)

human160seqs.loc[:,'MFE']=mfes
human160seqs.loc[:,'ss']=secstructs

human160seqs.to_csv('/fast/AG_Ohler/jdemoli/bachelorgit/data/MFE_human160_stats.csv')