import numpy as np
import pandas as pd
import src.constants
from Bio import SeqIO

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
