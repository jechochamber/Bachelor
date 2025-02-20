from __future__ import annotations
import numpy as np
import pandas as pd
import sys

sys.path.append("..")

#Needed functions
BASES = ("A", "C", "G", "U")


def substringcount(ini_str: str, sub_str: str) -> int:
    """We have to write this function ourselves because the.count() method for the markov_matrix function because alot of edge cases are not counted"""
    res = sum([1 for i in range(len(ini_str) - len(sub_str) + 1) if ini_str[i:i + len(sub_str)] == sub_str])

    return res


def markov_matrix(seqs: pd.DataFrame | str,seq_column : str,a : int|float = 0) -> tuple[dict, dict, np.ndarray]:
    tot_length = 0
    if type(seqs) == pd.DataFrame:
        for seq in seqs[seq_column]:
            tot_length += len(seq)
    if type(seqs) == str:
        tot_length = len(seqs)

    P_single = {}

    if type(seqs) == pd.DataFrame:
        for base in BASES:
            counter = 0
            for seq in seqs[seq_column]:
                counter += seq.count(base)
            P_single[base] = counter / tot_length
        tot_length -= len(seqs)  #prep for coming calculations, reasoning expanded upon in paper/thesis
    if type(seqs) == str:
        for base in BASES:
            counter = seqs.count(base)
            P_single[base] = (counter+a) / (tot_length+a*4)
        tot_length -= 1

    P_duplett = {}
    markov_matrix = np.zeros([4, 4])
    i, j = 0, 0

    if type(seqs) == pd.DataFrame:
        for base1 in BASES:
            for base2 in BASES:
                counter = 0
                for seq in seqs[seq_column]:
                    counter += substringcount(seq, base1 + base2)
                P_duplett[base1 + base2] = (counter+a) / (tot_length+a*16)
                markov_matrix[i, j] = P_duplett[base1 + base2] / P_single[base1]
                i += 1
            j += 1
            i = 0
    if type(seqs) == str:
        for base1 in BASES:
            for base2 in BASES:
                P_duplett[base1 + base2] = substringcount(seqs, base1 + base2) / tot_length
                if P_single[base1] > 0:
                    markov_matrix[i, j] = P_duplett[base1 + base2] / P_single[base1]
                else:
                    markov_matrix[i, j] = 0
                i += 1
            j += 1
            i = 0

    for i in range(len(markov_matrix.T)):
        sum_col = sum(markov_matrix.T[i, :])
        for j in range(len(markov_matrix.T[0])):
            if sum_col > 0:
                markov_matrix[j, i] = markov_matrix[j, i] / sum_col
    sanitycheck=[]
    for i in markov_matrix.T:
        sanitycheck.append(i.sum())
    #print(sanitycheck)
    return P_duplett, P_single, markov_matrix


# P(X|Nukleotid) →
# P(Nukleotid) ↓

def kld(matrix1: np.ndarray, matrix2: np.ndarray) -> float:
    dist1=matrix1.T.flatten()
    dist2 =  matrix2.T.flatten()
    kld=0
    for i,j in zip(dist1, dist2):
        if i==0 or j==0:
            kld+=0
        else:
            kld+=i*np.log(i/j)  #Add pseudocounts!!!!(+normalize)
                                #how to test pseudocounts?
    return kld

def main_loop(reference_data : pd.DataFrame,p_data : pd.DataFrame ,reference_seq_column : str,p_seq_column : str,a : int|float =0) -> list:
    """
    This loop calculates the Kullback-Leibner divergence between single Sequences(p_data) and reference sequence dataset(reference_data)
    The seq_header parameter is the name of column in which the sequence string is saved
    """
    reference_matrix = markov_matrix(reference_data,reference_seq_column)[2]
    print(reference_matrix)
    kldlist = []

    for seq in p_data[p_seq_column]:

        seq_matrix = markov_matrix(seq,p_seq_column,a=a)[2]
        kldlist.append(kld(reference_matrix,seq_matrix))



    return kldlist



reference_data = pd.read_csv("../data/new_dataset.csv")
reference_seqs = reference_data.loc[reference_data['Group'] == 'Human']
p_data = pd.read_csv("../data/random_train_pc.csv")
for i in p_data.index:
    p_data.loc[i, 'utr'] = p_data.loc[i, 'utr'].replace('T', 'U')

#print(p_data.head())

p_data["KLD2"]=main_loop(reference_seqs,p_data,"Sequence","utr",a=2)
print(p_data["KLD2"])
p_data.to_csv("../data/random_train_pc.csv",index=False)


#possible name for cutoff function: "kld sampling" oder "dataloader"