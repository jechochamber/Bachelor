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
            gnirtsqes = seqstring[::-1]  # We reverse the sequence first to get our reading frames in relation to the CDS(frame 1 is the same reading frame as the CDS
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


def uORFs(seqs):
    frames = codons(seqs)
    template = {
        "uORFs": list(np.zeros(len(frames))),
        "ouORFs": list(np.zeros(len(frames))),
        "mean_uORF_length": list(np.zeros(len(frames))),

        "CUG_uORFs": list(np.zeros(len(frames))),
        "CUG_ouORFs": list(np.zeros(len(frames))),
        "mean_CUG_uORF_length": list(np.zeros(len(frames))),

        "GUG_uORFs": list(np.zeros(len(frames))),
        "GUG_ouORFs": list(np.zeros(len(frames))),
        "mean_GUG_uORF_length": list(np.zeros(len(frames))),
    }
    framedict = {}
    for frame in range(3):
        framedict[f"frame{frame + 1}"] = pd.DataFrame(template)
        for index in range(len(frames[f"frame{frame + 1}"])):
            seq = frames[f"frame{frame + 1}"][index]
            lengths = []
            CUG_lengths = []
            GUG_lengths = []
            for codon in range(len(seq)):
                if seq[codon] == 'AUG':
                    uORF_length = 0
                    ouORF_bool = True
                    for nextcodon in range(len(seq) - codon - 1):
                        uORF_length += 3
                        if seq[nextcodon + codon + 1] in src.constants.stopcodons:
                            ouORF_bool = False
                            framedict[f"frame{frame + 1}"].loc[index, "uORFs"] += 1
                            lengths.append(uORF_length)
                            break
                    if ouORF_bool:
                        framedict[f"frame{frame + 1}"].loc[index, "ouORFs"] += 1
                        break

                if seq[codon] == 'CUG':
                    uORF_length = 0
                    ouORF_bool = True
                    for nextcodon in range(len(seq) - codon - 1):
                        uORF_length += 3
                        if seq[nextcodon + codon + 1] in src.constants.stopcodons:
                            ouORF_bool = False
                            framedict[f"frame{frame + 1}"].loc[index, "CUG_uORFs"] += 1
                            CUG_lengths.append(uORF_length)
                            break
                    if ouORF_bool:
                        framedict[f"frame{frame + 1}"].loc[index, "CUG_ouORFs"] += 1
                        break

                if seq[codon] == 'GUG':
                    uORF_length = 0
                    ouORF_bool = True
                    for nextcodon in range(len(seq) - codon - 1):
                        uORF_length += 3
                        if seq[nextcodon + codon + 1] in src.constants.stopcodons:
                            ouORF_bool = False
                            framedict[f"frame{frame + 1}"].loc[index, "GUG_uORFs"] += 1
                            GUG_lengths.append(uORF_length)
                            break
                    if ouORF_bool:
                        framedict[f"frame{frame + 1}"].loc[index, "GUG_ouORFs"] += 1
                        break

            if len(lengths) > 0:
                framedict[f"frame{frame + 1}"].loc[index, "mean_uORF_length"] = sum(lengths) / len(lengths)

            if len(CUG_lengths) > 0:
                framedict[f"frame{frame + 1}"].loc[index, "mean_CUG_uORF_length"] = sum(CUG_lengths) / len(CUG_lengths)
            if len(GUG_lengths) > 0:
                framedict[f"frame{frame + 1}"].loc[index, "mean_GUG_uORF_length"] = sum(GUG_lengths) / len(GUG_lengths)

    counts = pd.DataFrame()

    for key in framedict:
        for column in framedict[key].columns:
            counts[f"{key}_{column}"] = framedict[key].loc[:, column]

    counts["uORFs"] = counts["frame1_uORFs"] + counts["frame2_uORFs"] + counts["frame3_uORFs"]
    counts["CUG_uORFs"] = counts["frame1_CUG_uORFs"] + counts["frame2_CUG_uORFs"] + counts["frame3_CUG_uORFs"]
    counts["GUG_uORFs"] = counts["frame1_GUG_uORFs"] + counts["frame2_GUG_uORFs"] + counts["frame3_GUG_uORFs"]
    counts["ouORFs"] = counts["frame1_ouORFs"] + counts["frame2_ouORFs"] + counts["frame3_ouORFs"]
    counts["CUG_ouORFs"] = counts["frame1_CUG_ouORFs"] + counts["frame2_CUG_ouORFs"] + counts["frame3_CUG_ouORFs"]
    counts["GUG_ouORFs"] = counts["frame1_GUG_ouORFs"] + counts["frame2_GUG_ouORFs"] + counts["frame3_CUG_ouORFs"]
    counts["all uORFs"] = counts["uORFs"] + counts["CUG_uORFs"] + counts["GUG_uORFs"]
    counts["all ouORFs"] = counts["ouORFs"] + counts["CUG_ouORFs"] + counts["GUG_ouORFs"]

    counts = counts.replace(0, np.nan)
    counts["mean_uORF_length"] = counts.loc[:, ["frame1_mean_uORF_length", "frame2_mean_uORF_length",
                                                "frame3_mean_uORF_length"]].mean(axis=1)
    counts["mean_CUG_uORF_length"] = counts.loc[:, ["frame1_mean_CUG_uORF_length", "frame2_mean_CUG_uORF_length",
                                                    "frame3_mean_CUG_uORF_length"]].mean(axis=1)
    counts["mean_GUG_uORF_length"] = counts.loc[:, ["frame1_mean_GUG_uORF_length", "frame2_mean_GUG_uORF_length",
                                                    "frame3_mean_GUG_uORF_length"]].mean(axis=1)
    counts["all mean lengths"] = counts.loc[:,
                                 ["mean_uORF_length", "mean_CUG_uORF_length", "mean_GUG_uORF_length"]].mean(axis=1)
    counts = counts.replace(np.nan, 0)

    cols_in_front = ["all uORFs", "all ouORFs", "all mean lengths", "uORFs", "CUG_uORFs", "GUG_uORFs", "ouORFs",
                     "CUG_ouORFs", "GUG_ouORFs", "mean_uORF_length", "mean_CUG_uORF_length", "mean_GUG_uORF_length"]
    counts = counts.loc[:,
             [c for c in cols_in_front if c in counts.columns] + [c for c in counts if c not in cols_in_front]]

    if type(seqs) == pd.core.frame.DataFrame:
        counts.insert(loc=0, column='SeqID', value=seqs.index)
        counts.set_index('SeqID', inplace=True)

    return counts

def GC_dinucleotides(seqs):
    pass
