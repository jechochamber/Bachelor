import numpy as np
import pandas as pd
bases = np.array(["A","T","G","C"])
rng = np.random.default_rng()
seqs = rng.choice(bases, size=(1000,100))

df=pd.DataFrame(seqs)
df.to_csv("toy_dataset.csv",header=False,index=False)