
import pyranges as pr
from pprint import pprint
import numpy as np

import sys
sys.path.append("..")
import src.human_sequences

gtf= pr.read_gtf("../data/Homo_sapiens.GRCh38.113.gtf")
gtf_df= gtf.df
fputrs=gtf_df[gtf_df['Feature']=='five_prime_utr']
print(fputrs)

