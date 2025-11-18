import pandas as pd
import numpy as np
import ViennaRNA as vrna
import pandas as pd
import numpy as np
import ViennaRNA as vrna

dataset = pd.read_csv('../data/new_fish_uORFstats.csv')


mfes = []


# divides a string of DNA bases into windows of size window_size with a stride,
# only las window is allowed to be shorter than window_size (check while loop)
# DNA bases are converted to RNA bases in the process
def get_windows(seq, window_size=100, stride=50):
    seq_string = str(seq)
    #seq_string = seq_string.replace('T', 'U')

    windows = []
    start = 0
    end = window_size

    if len(seq_string) <= window_size:
        return [seq_string]

    while end <= len(seq_string):
        windows.append(seq_string[start:end])
        start += stride
        end += stride
    windows.append(seq_string[len(seq_string)-window_size:len(seq_string)])

    return windows


# gets mfe for each row in features dataframe, keeping mfe for window with lowest mfe per 5'UTR
x=0
for fputr in dataset['Sequence']:
    x+=1
    wins = get_windows(fputr)
    if x<=10:
        print(x,wins)
        lens=[]
        for win in wins:
            lens.append(len(win))
        print(x,lens)

    lowest_mfe = float('inf')

    for win in wins:
        string = str(win)
        (ss, mfe) = vrna.fold(string)
        mfe = mfe / len(string)

        if mfe < lowest_mfe:
            lowest_mfe = mfe

    mfes.append(lowest_mfe)

dataset.loc[:,'sliding window mfe'] = mfes

dataset.to_csv('../data/new_fish_mfe.csv')



