import numpy as np
import pandas as pd
import ViennaRNA as vrna
import polygraph.input, polygraph.sequence, polygraph.visualize

seqs=polygraph.input.read_seqs('/fast/AG_Ohler/jdemoli/bachelorgit/data/full_dataset.txt',incl_ids=True)
humanseqs=seqs.loc[seqs['Group']=='Human']
randomdlseqs=seqs.loc[seqs['Group']=='Random_dl_gc_adjusted']
lengthdistributed=pd.concat([humanseqs,randomdlseqs])
lengthdistributed['Sequence Length']=lengthdistributed.Sequence.apply(len)

secstructs=[]
mfes=[]

for seq in lengthdistributed["Sequence"]:
    (ss,mfe)=vrna.fold(seq)
    secstructs.append(ss)
    mfes.append(mfe)

lengthdistributed['MFE']=mfes
lengthdistributed['ss']=secstructs

lengthdistributed.to_csv('/fast/AG_Ohler/jdemoli/bachelorgit/data/MFE_dl_stats.csv')