import numpy as np
import pandas as pd
import src.constants
from Bio import SeqIO
from scipy.optimize import curve_fit
from collections import Counter

"""This module is used to generate random DNA and RNA sequences."""
def randomseqs(seqnum, seqlength, group_name, gc_content=0.5, seed=None, BASES=src.constants.RNABASES):
    """generates a random set of DNA sequences for a group
    BASES: list of bases
    seqnum: number of sequences
    seqlength: length of sequence
    group_name: name of group
    seq_final= output der Funktion, numpy array mit seqnum Sequenzen die seqlength lang sind
    """
    BASES = np.array(BASES)
    rng = np.random.default_rng(seed=seed)
    probabilities = [(1 - gc_content) / 2, gc_content / 2, gc_content / 2, (1 - gc_content) / 2]
    seqs_raw = rng.choice(BASES, size=(seqnum, seqlength), p=probabilities)
    seperator = ""
    seqs_final = []
    for seq in seqs_raw:
        seqs_final.append(seperator.join(seq))
    seqs_final = np.array(seqs_final)
    group = np.array([group_name] * seqnum)
    seqs_final = np.vstack((seqs_final, group)).T

    return seqs_final

#We will define different statistical distributions that can be used in later methods
def lognormalfunc(x,myu,sigma,a):
    return a*np.exp(-1*(np.log(x)-myu)**2/(2*sigma)**2)/(x*sigma*np.sqrt(2*np.pi))

def get_params(data,group,func=lognormalfunc):
    seqs = data.loc[data['Group'] == group]
    lengthdata = list(seqs['Sequence Length'])
    lengthdata = Counter(lengthdata)
    x, y = [], []
    for i, j in lengthdata.items():
        x.append(i)
        y.append(j)
    y = np.array(y)

    params, cov_matrix = curve_fit(func, x, y)
    return params

def randomlengthseqs(params,maxlength,seqnum, group_name,func=lognormalfunc, gc_content=0.5, seed=None, BASES=src.constants.RNABASES):
    x=np.arange(1,maxlength)
    distribution=func(x,*params)
    distribution=distribution/sum(distribution)
    gcprobabilities = [(1 - gc_content) / 2, gc_content / 2, gc_content / 2, (1 - gc_content) / 2]
    rng=np.random.default_rng(seed=seed)
    lengths=rng.choice(x, seqnum , p=distribution)
    seqs_raw=[]
    for i in lengths:
        seqs_raw.append(rng.choice(BASES,size=(1,i),p=gcprobabilities))
    seperator = ""
    seqs_final = []
    for seq in seqs_raw:
        seqs_final.append(seperator.join(seq[0]))
    seqs_final = np.array(seqs_final)
    group = np.array([group_name] * seqnum)
    seqs_final = np.vstack((seqs_final, group)).T
    return seqs_final
