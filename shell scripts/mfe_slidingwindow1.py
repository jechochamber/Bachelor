import pandas as pd
import numpy as np
import ViennaRNA as vrna
import pandas as pd
import numpy as np
import ViennaRNA as vrna
import sys
sys.path.append("..")
import src.rna_analysis

dataset = pd.read_csv('/fast/AG_Ohler/jdemoli/bachelorgit/data/small_dataset.csv')

dataset.loc[:,'sliding window mfe']=src.rna_analysis.sliding_window_mfe(dataset)

dataset.to_csv('/fast/AG_Ohler/jdemoli/bachelorgit/data/small_dataset.csv')

