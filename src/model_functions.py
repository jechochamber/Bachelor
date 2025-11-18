import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import torch.nn.functional as F
import torch.optim as optim
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr
import pickle
import numpy as np
import pandas as pd
import sys

# Collection of general functions used in models
def one_hot(seqs):
    """One-hot encode a sequence dataset"""
    conversion_dict = {
        'A': np.array([1.0, 0.0, 0.0, 0.0]),
        'C': np.array([0.0, 1.0, 0.0, 0.0]),
        'G': np.array([0.0, 0.0, 1.0, 0.0]),
        'U': np.array([0.0, 0.0, 0.0, 1.0]),
        'T': np.array([0.0, 0.0, 0.0, 1.0]),
    }
    enc_seqs = []
    for seq in seqs:
        enc_arr = conversion_dict[seq[0]]
        for i in seq[1:]:
            enc_arr = np.vstack((enc_arr, conversion_dict[i]))
        # enc_arr=enc_arr.T.reshape((1,4,50))
        enc_arr = torch.tensor(enc_arr.T, dtype=torch.float32)
        enc_seqs.append(enc_arr)
    enc_seqs = torch.tensor(np.array(enc_seqs), dtype=torch.float32)

    return enc_seqs

def earlystopper(val_loss,patience=10,epsilon=1e-7):
    """Earlystopping function"""
    global es_stopcounter
    global es_min_loss
    if patience >= es_stopcounter:
        if val_loss <= es_min_loss - epsilon:
            es_min_loss = val_loss
            es_stopcounter = 0
        else:
            es_stopcounter += 1
        return False
    else:
        return True