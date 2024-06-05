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

def uORFs(seqs):
    """This Method takes a pandas Dataframe of sequences or a string
     of one sequence outputs a pandas DataFrame countaining the number
     of uORFs, ouORFs in each sequence and the mean length of uORFs per Sequence"""
    frames = codons(seqs)
    counts = {
        "uORF_countssum": list(np.zeros(len(frames))),
        "ouORF_countssum": list(np.zeros(len(frames))),

        "uORF_counts1": list(np.zeros(len(frames))),
        "uORF_counts2": list(np.zeros(len(frames))),
        "uORF_counts3": list(np.zeros(len(frames))),

        "uORF_meanlength1": list(np.zeros(len(frames))),
        "uORF_meanlength2": list(np.zeros(len(frames))),
        "uORF_meanlength3": list(np.zeros(len(frames))),

        "ouORF_counts1": list(np.zeros(len(frames))),
        "ouORF_counts2": list(np.zeros(len(frames))),
        "ouORF_counts3": list(np.zeros(len(frames))),
    }
    counts = pd.DataFrame(counts)

    '''We created the Dataframe in which we will store our data,
    in the following code we iterate through the 3 different reading frames 
    and search for start and stopcodons and depending on what we find, save our counts.'''

    for index in range(len(frames["frame1"])):
        seq = frames["frame1"][index]
        for codon in range(len(seq)):
            lengths = []
            if seq[codon] == "AUG":
                uORF_length = 0
                ouORF_bool = True
                for nextcodon in range(len(seq)):
                    uORF_length += 3
                    if seq[nextcodon] in src.constants.stopcodons:
                        ouORF_bool = False
                        counts.loc[index, "uORF_counts1"] += 1
                        counts.loc[index, "uORF_countssum"] += 1
                        lengths.append(uORF_length)
                        break
                if ouORF_bool:
                    counts.loc[index, "ouORF_countssum"] += 1
                    counts.loc[index, "ouORF_counts1"] += 1
                    break

    for index in range(len(frames["frame2"])):
        seq = frames["frame2"][index]
        for codon in range(len(seq)):
            lengths = []
            if seq[codon] == "AUG":
                uORF_length = 0
                ouORF_bool = True
                for nextcodon in range(len(seq)):
                    uORF_length += 3
                    if seq[nextcodon] in src.constants.stopcodons:
                        ouORF_bool = False
                        counts.loc[index, "uORF_counts2"] += 1
                        counts.loc[index, "uORF_countssum"] += 1
                        lengths.append(uORF_length)
                        break
                if ouORF_bool:
                    counts.loc[index, "ouORF_countssum"] += 1
                    counts.loc[index, "ouORF_counts2"] += 1
                    break

    for index in range(len(frames["frame3"])):
        seq = frames["frame3"][index]
        for codon in range(len(seq)):
            lengths = []
            if seq[codon] == "AUG":
                uORF_length = 3
                ouORF_bool = True
                for nextcodon in range(len(seq)):
                    uORF_length += 3
                    if seq[nextcodon] in src.constants.stopcodons:
                        ouORF_bool = False
                        counts.loc[index, "uORF_counts3"] += 1
                        counts.loc[index, "uORF_countssum"] += 1
                        lengths.append(uORF_length)
                        break
                if ouORF_bool:
                    counts.loc[index, "ouORF_countssum"] += 1
                    counts.loc[index, "ouORF_counts3"] += 1
                    break

    '''Here we add the last columns to our DataFrame and clean it up for the final output'''

    counts.iloc[:, 2:] = counts.iloc[:, 2:].replace(0, np.nan)
    counts["uORF_meanlength"] = counts.loc[:, ["uORF_meanlength1", "uORF_meanlength2", "uORF_meanlength3"]].mean(
        axis=1)
    counts = counts.replace(np.nan, 0)
    cols_in_front = ["uORF_countssum", "uORF_meanlength", "ouORF_countssum"]
    counts = counts.loc[:,[c for c in cols_in_front if c in counts.columns] + [c for c in counts if c not in cols_in_front]]
    counts.insert(loc=0,column='SeqID',value=seqs.index)
    counts.set_index('SeqID', inplace=True)

    return counts