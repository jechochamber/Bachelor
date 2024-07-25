import numpy as np
import pandas as pd
import ViennaRNA as vrna
import polygraph.input, polygraph.sequence, polygraph.visualize

seqs=polygraph.input.read_seqs('/fast/AG_Ohler/jdemoli/bachelorgit/data/full_dataset.txt',incl_ids=True)


secstructs=[]
mfes=[]
for seq in seqs["Sequence"]:
    (ss,mfe)=vrna.fold(seq)
    secstructs.append(ss)
    mfes.append(mfe)
seqs['MFE']=mfes
seqs['ss']=secstructs

seqs.to_csv('/fast/AG_Ohler/jdemoli/bachelorgit/data/MFE_stats.csv')