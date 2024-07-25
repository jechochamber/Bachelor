import numpy as np
import pandas as pd
import src.constants

"""This module contains functions to calculate different statistics of uORFs in RNA sequences."""


def codons(seqs):
    """This method splits sequences into codons and their respective reading frames"""
    codons = {
        "frame1": [],
        "frame2": [],
        "frame3": []
    }
    '''We have the option of either reading a pandas DataFrame of sequence or a single sequence string'''
    if type(seqs) == pd.DataFrame:
        '''The algorithm here is similar to the one for a single sequence string.'''
        seqstrings = seqs["Sequence"].values
        for seqstring in seqstrings:
            gnirtsqes = seqstring[
                        ::-1]  # We reverse the sequence first to get our reading frames in relation to the CDS(frame 1 is the same reading frame as the CDS
            frame1, frame2, frame3 = [], [], []

            for i in range(len(gnirtsqes)):
                if i + 3 <= len(gnirtsqes):
                    if i % 3 == 0:
                        frame1.append(gnirtsqes[i] + gnirtsqes[i + 1] + gnirtsqes[i + 2])
                    if i % 3 == 1:
                        frame2.append(gnirtsqes[i] + gnirtsqes[i + 1] + gnirtsqes[i + 2])
                    if i % 3 == 2:
                        frame3.append(gnirtsqes[i] + gnirtsqes[i + 1] + gnirtsqes[i + 2])
                else:
                    break

            frame1, frame2, frame3 = frame1[::-1], frame2[::-1], frame3[
                                                                 ::-1]  # we reverse our list of codons to be in the correct order
            for i in range(len(frame1)):
                frame1[i] = frame1[i][::-1]
            for i in range(len(frame2)):
                frame2[i] = frame2[i][::-1]
            for i in range(len(frame3)):
                frame3[i] = frame3[i][::-1]
            '''we then reversed our codons themselves'''
            codons['frame1'].append(frame1)
            codons['frame2'].append(frame2)
            codons['frame3'].append(frame3)

    if type(seqs) == str:
        seqstring = seqs
        gnirtsqes = seqstring[::-1]

        frame1, frame2, frame3 = [], [], []

        for i in range(len(gnirtsqes)):
            if i + 3 <= len(gnirtsqes):
                if i % 3 == 0:
                    frame1.append(gnirtsqes[i] + gnirtsqes[i + 1] + gnirtsqes[i + 2])
                if i % 3 == 1:
                    frame2.append(gnirtsqes[i] + gnirtsqes[i + 1] + gnirtsqes[i + 2])
                if i % 3 == 2:
                    frame3.append(gnirtsqes[i] + gnirtsqes[i + 1] + gnirtsqes[i + 2])
            else:
                break

        frame1, frame2, frame3 = frame1[::-1], frame2[::-1], frame3[::-1]
        for i in range(len(frame1)):
            frame1[i] = frame1[i][::-1]
        for i in range(len(frame2)):
            frame2[i] = frame2[i][::-1]
        for i in range(len(frame3)):
            frame3[i] = frame3[i][::-1]
        codons['frame1'].append(frame1)
        codons['frame2'].append(frame2)
        codons['frame3'].append(frame3)
    return pd.DataFrame(codons)


def uORFs(seqs, startcodon):
    """This method calculates uORFs in RNA sequences."""
    frames = codons(seqs)
    template = {
        "uORFs": list(np.zeros(len(frames))),
        "ouORFs": list(np.zeros(len(frames))),
        "mean_uORF_length": list(np.zeros(len(frames))),
        "max_uORF_length": list(np.zeros(len(frames))),
    }
    framedict = {}
    for frame in range(3):
        framedict[f"frame{frame + 1}"] = pd.DataFrame(template)
        for index in range(len(frames[f"frame{frame + 1}"])):
            seq = frames[f"frame{frame + 1}"][index]
            lengths = []
            for codon in range(len(seq)):
                if seq[codon] == startcodon:
                    uORF_length = 0
                    ouORF_bool = True
                    for nextcodon in range(len(seq) - codon - 1):
                        uORF_length += 3
                        if seq[nextcodon + codon + 1] in src.constants.stopcodons:
                            ouORF_bool = False
                            uORF_length +=3
                            framedict[f"frame{frame + 1}"].loc[index, "uORFs"] += 1
                            lengths.append(uORF_length)
                            break
                    if ouORF_bool:
                        framedict[f"frame{frame + 1}"].loc[index, "ouORFs"] += 1
                        break

            if len(lengths) > 0:
                framedict[f"frame{frame + 1}"].loc[index, "mean_uORF_length"] = sum(lengths) / len(lengths)
                framedict[f"frame{frame + 1}"].loc[index, "max_uORF_length"] = max(lengths)

    counts = pd.DataFrame()

    for key in framedict:
        for column in framedict[key].columns:
            counts[f"{key}_{column}"] = framedict[key].loc[:, column]

    counts["uORFs"] = counts["frame1_uORFs"] + counts["frame2_uORFs"] + counts["frame3_uORFs"]
    counts["ouORFs"] = counts["frame1_ouORFs"] + counts["frame2_ouORFs"] + counts["frame3_ouORFs"]

    counts = counts.replace(0, np.nan)
    counts["mean_uORF_length"] = counts.loc[:, ["frame1_mean_uORF_length", "frame2_mean_uORF_length",
                                                "frame3_mean_uORF_length"]].mean(axis=1)

    counts["max_uORF_length"] = counts.loc[:,
                                ["frame1_max_uORF_length", "frame2_max_uORF_length", "frame3_max_uORF_length"]].max(
        axis=1)
    counts = counts.replace(np.nan, 0)

    cols_in_front = ["uORFs", "ouORFs", "mean_uORF_length", "max_uORF_length"]
    counts = counts.loc[:,
             [c for c in cols_in_front if c in counts.columns] + [c for c in counts if c not in cols_in_front]]

    if type(seqs) == pd.core.frame.DataFrame:
        counts.insert(loc=0, column='SeqID', value=seqs.index)
        counts.set_index('SeqID', inplace=True)

    return counts


def counts_concat_all(seqs, AUG_counts, CUG_counts, GUG_counts, UUG_counts,ACG_counts):
    """this function summarizes all of the count matrices created with the previous function"""
    counts = pd.DataFrame()
    countdict = {"AUG": AUG_counts, "CUG": CUG_counts, "GUG": GUG_counts,"UUG": UUG_counts,"ACG": ACG_counts}
    for key in countdict:
        for column in countdict[key].columns:
            counts[f"{key}_{column}"] = countdict[key].loc[:, column]

    counts["all uORFs"] = counts["AUG_uORFs"] + counts["CUG_uORFs"] + counts["GUG_uORFs"]+ counts["UUG_uORFs"] + counts["ACG_uORFs"]
    counts["all ouORFs"]=counts["AUG_ouORFs"] + counts["CUG_ouORFs"] + counts["GUG_ouORFs"]+ counts["UUG_ouORFs"] + counts["ACG_ouORFs"]
    counts = counts.replace(0, np.nan)
    counts["all mean uORF lengths"] = counts.loc[:, ["AUG_mean_uORF_length", "CUG_mean_uORF_length","GUG_mean_uORF_length","UUG_mean_uORF_length","ACG_mean_uORF_length"]].mean(axis=1)
    counts["all max uORF lengths"] = counts.loc[:, ["AUG_max_uORF_length", "CUG_max_uORF_length","GUG_max_uORF_length","UUG_max_uORF_length","ACG_max_uORF_length"]].max(axis=1)
    counts = counts.replace(np.nan, 0)
    cols_in_front = ["all uORFs", "all ouORFs", "all mean uORF lengths", "all max uORF lengths", "AUG_uORFs", "CUG_uORFs",
                     "GUG_uORFs","UUG_uORFs", "ACG_uORFs" , "AUG_ouORFs", "CUG_ouORFs", "GUG_ouORFs","UUG_ouORFs","ACG_ouORFs", "AUG_mean_uORF_length", "CUG_mean_uORF_length",
                     "GUG_mean_uORF_length","UUG_mean_uORF_length","ACG_mean_uORF_length", "AUG_max_uORF_length", "CUG_max_uORF_length", "GUG_max_uORF_length","UUG_max_uORF_length", "ACG_max_uORF_length"]
    counts = counts.loc[:,
             [c for c in cols_in_front if c in counts.columns] + [c for c in counts if c not in cols_in_front]]

    if type(seqs) == pd.core.frame.DataFrame:
        counts.insert(loc=0, column='SeqID', value=seqs.index)
        counts.set_index('SeqID', inplace=True)
    return counts

def GC_dinucleotides(seqs):
    """WIP"""
    pass


def fast_uORfs(seqs, startcodon):
    """This function is still WIP and is not faster(yet)"""
    seq_pos = {"frame1": [], "frame2": [], "frame3": []}
    if type(seqs) == pd.DataFrame:
        seqstrings = seqs["Sequence"].values
    if type(seqs) == str:
        seqstrings = [seqs]
    len_loop = len(seqstrings)
    template = {
        "uORFs": list(np.zeros(len_loop)),
        "ouORFs": list(np.zeros(len_loop)),
        "mean_uORF_length": list(np.zeros(len_loop)),
        "max_uORF_length": list(np.zeros(len_loop)),

        "frame1_uORFs": list(np.zeros(len_loop)),
        "frame1_ouORFs": list(np.zeros(len_loop)),
        "frame1_mean_uORF_length": list(np.zeros(len_loop)),
        "frame1_max_uORF_length": list(np.zeros(len_loop)),

        "frame2_uORFs": list(np.zeros(len_loop)),
        "frame2_ouORFs": list(np.zeros(len_loop)),
        "frame2_mean_uORF_length": list(np.zeros(len_loop)),
        "frame2_max_uORF_length": list(np.zeros(len_loop)),

        "frame3_uORFs": list(np.zeros(len_loop)),
        "frame3_ouORFs": list(np.zeros(len_loop)),
        "frame3_mean_uORF_length": list(np.zeros(len_loop)),
        "frame3_max_uORF_length": list(np.zeros(len_loop)),
    }
    counts = pd.DataFrame(template)

    posdict = {"frame1": [[] for _ in range(len(seqstrings))], "frame2": [[] for _ in range(len(seqstrings))],
               "frame3": [[] for _ in range(len(seqstrings))]}
    for index in range(len(seqstrings)):
        seq = seqstrings[index]
        for i in range(len(seq)):
            rev_i = len(seq) - i
            frame = (rev_i % 3)
            if seq[i:i + 3] == startcodon:
                posdict[f'frame{frame + 1}'][index].append(i + 1)
            if seq[i:i + 3] in src.constants.stopcodons:
                posdict[f'frame{frame + 1}'][index].append((i + 1) * -1)

    for frame in range(3):
        poslist = posdict[f'frame{frame + 1}']
        for index in range(len_loop):
            pos = poslist[index]
            if len(pos) > 1:
                lengths = []
                if pos[-1] > 0:
                    counts.loc[index, "ouORFs"] += 1
                    counts.loc[index, f"frame{frame + 1}_ouORFs"] += 1

                for i in range(len(pos)):
                    if pos[i] > 0:

                        for j in range(len(pos) - i - 1):
                            if pos[i + j + 1] < 0:
                                counts.loc[index, "uORFs"] += 1
                                counts.loc[index, f"frame{frame + 1}_uORFs"] += 1
                                lengths.append((pos[i] + pos[i + j + 1]) * -1 + 3)
                                break
                if len(lengths) > 0:
                    counts.loc[index, f"frame{frame + 1}_mean_uORF_length"] = sum(lengths) / len(lengths)
                    counts.loc[index, f"frame{frame + 1}_max_uORF_length"] = max(lengths)
                if len(pos) == 1 and pos[0] > 0:
                    counts.loc[index, "ouORFs"] += 1
                    counts.loc[index, f"frame{frame + 1}_ouORFs"] += 1

    counts["mean_uORF_length"] = counts.loc[:, ["frame1_mean_uORF_length", "frame2_mean_uORF_length",
                                                "frame3_mean_uORF_length"]].mean(axis=1)
    counts["max_uORF_length"] = counts.loc[:,
                                ["frame1_max_uORF_length", "frame2_max_uORF_length", "frame3_max_uORF_length"]].max(
        axis=1)

    return counts

