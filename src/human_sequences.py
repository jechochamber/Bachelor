import numpy as np
import pandas as pd
import src.constants
from Bio import SeqIO

"""This module contains functions to read and tweak human sequences."""

def readFASTA(filename,groupname, BASES=src.constants.DNABASES):
    """We read FASTA file and return a numpy array of the sequences."""
    seqs = []
    seqIDs = []
    for seq_record in SeqIO.parse(filename, "fasta"):
        seqs.append(str(seq_record.seq))
        seqIDs.append(seq_record.id)
    seqs = np.vstack((np.array(seqIDs), np.array(seqs))).T
    group = np.array([[groupname] * len(seqs)])
    seqs = np.concatenate((seqs, group.T), axis=1)
    if BASES == src.constants.DNABASES:
        for i in range(len(seqs)):
            seqs[i, 1] = seqs[i, 1].replace('T', 'U')
    return seqs
def delDupes(seqs):
    """Deletes any Sequence that occurs twice or more in the dataset"""
    seen = set()
    dupeIDs = []
    for x in range(len(seqs)):
        if seqs[x, 1] in seen:
            seen.add(seqs[x, 1])
            dupeIDs.append(x)
        else:
            seen.add(seqs[x, 1])
    return np.delete(seqs, dupeIDs, 0)