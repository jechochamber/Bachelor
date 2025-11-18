from __future__ import annotations
import numpy as np
import pandas as pd
import sys
import pickle
import time

sys.path.append("..")

#Needed functions
BASES = ("A", "C", "G", "U")


def substringcount(ini_str: str, sub_str: str) -> int:
    """We have to write this function ourselves because the.count() method for the markov_matrix function because alot of edge cases are not counted"""
    res = sum([1 for i in range(len(ini_str) - len(sub_str) + 1) if ini_str[i:i + len(sub_str)] == sub_str])

    return res


def markov_matrix(seqs: pd.DataFrame | str,seq_column : str = None,a : int|float = 0) -> tuple[dict, dict, np.ndarray]:
    """This function is rewritten here to be more functional in further loops"""
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
                P_duplett[base1 + base2] = counter / tot_length
                markov_matrix[i, j] = P_duplett[base1 + base2] / P_single[base1]
                i += 1
            j += 1
            i = 0
    if type(seqs) == str:
        for base1 in BASES:
            for base2 in BASES:
                P_duplett[base1 + base2] = (substringcount(seqs, base1 + base2) +a) / (tot_length+a*16)
                if P_single[base1] > 0:
                    markov_matrix[i, j] = P_duplett[base1 + base2] / P_single[base1]

                i += 1
            j += 1
            i = 0

    for i in range(len(markov_matrix.T)):
        sum_col = sum(markov_matrix.T[i, :])
        for j in range(len(markov_matrix.T[0])):
            if sum_col > 0:
                markov_matrix[j, i] = markov_matrix[j, i] / sum_col

    return P_duplett, P_single, markov_matrix


# P(X|Nukleotid) →
# P(Nukleotid) ↓

def kld(p_matrix: np.ndarray | float, q_matrix: np.ndarray | float) -> float:
    """Returns Kullback-Leibler-distance between two matrices"""
    p_data=p_matrix.T.flatten()
    q_data =  q_matrix.T.flatten()
    kld=0
    for p,q in zip(p_data, q_data):
        if p==0 or q==0:
            kld+=0
        else:
            kld+=p*np.log(p/q)
    return kld
def jsd(p_matrix: np.ndarray, q_matrix: np.ndarray) -> float:
    """Returns Jensen-Shannon-distance between two matrices"""
    p_data=p_matrix.T.flatten()
    q_data = q_matrix.T.flatten()
    jsd=0
    i=0
    for p,q in zip(p_data, q_data):
        i+=1
        if p==0 or q==0:
            jsd+=0
        else:
            m=0.5*(p+q)
            D=lambda a,b : a*np.log(a/b)#kld as a lambda function, to not overwrite the kld function that already exists
            jsd+=0.5*D(p,m)+0.5*D(q,m)
    try:
        return np.sqrt(jsd)
    except:
        print('following value couldn`t be calculated:', jsd)
        return np.nan



def jsd_loop(p_data, p_seq_column, q_matrix, a=1):
    """loop to calculate all JSDs for a dataset"""
    jsd_list = []

    for seq in p_data[p_seq_column]:
        p_matrix = markov_matrix(seq, a=a)[2]
        jsd_list.append(jsd(p_matrix, q_matrix))
    return jsd_list


def main_loop(reference_data : pd.DataFrame,q_data : pd.DataFrame ,reference_seq_column : str,q_seq_column : str,a : int|float =0) -> list:
    """
    This loop calculates the Kullback-Leibner divergence between single Sequences(p_data) and reference sequence dataset(reference_data)
    The seq_header parameter is the name of column in which the sequence string is saved
    """
    reference_matrix = markov_matrix(reference_data,reference_seq_column)[2]
    print(reference_matrix)
    kldlist = []

    for seq in q_data[q_seq_column]:

        seq_matrix = markov_matrix(seq,q_seq_column,a=a)[2]
        kldlist.append(kld(reference_matrix,seq_matrix))




    return kldlist

#WIP
def distance_calculation(reference_data : pd.DataFrame,q_data : pd.DataFrame ,reference_seq_column : str,q_seq_column : str,seqnum) -> tuple:
    """Function specified in the thesis to calculate a more accurate distance between sequences"""
    reference_seqs=np.random.choice(np.array(reference_data[reference_seq_column]),seqnum)
    q_seqs=np.array(q_data[q_seq_column])
    human_dist=[]
    for seq in reference_seqs:
        for second_seq in reference_seqs:
            human_dist.append(jsd(markov_matrix(seq)[2],markov_matrix(second_seq)[2]))
    list_of_dists=[]
    for seq in q_seqs:
        random_dist = []
        for second_seq in reference_seqs:
            jsd_of_seqs=jsd(markov_matrix(seq,a=1)[2],markov_matrix(second_seq,a=1)[2])
            random_dist.append(jsd_of_seqs)

        list_of_dists.append(random_dist)
    return human_dist,list_of_dists
