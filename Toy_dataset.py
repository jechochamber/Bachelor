import numpy as np
import pandas as pd

class RanSequence:
    BASES = np.array(["A", "T", "G", "C"])
    def __init__(self, seqnum, seqlength, group_name="Random"):
        self.seqnum = seqnum
        self.seqlength = seqlength
        self.group_name = group_name

    def __str__(self):
        print(self.seqs)

    def save(self,fname="toy_dataset.txt",delimiter="\t",newline="\n"):
        np.savetxt(fname, self.seqs_final, delimiter=delimiter, newline=newline, fmt='%s')

    @staticmethod
    def seqs():

        rng = np.random.default_rng()
        seqs_raw = rng.choice(BASES, size=(self.n, self.nb))
        seperator = ""
        seqs_final = []
        for seq in seqs_raw:
            seqs_final.append(seperator.join(seq))
        seqs_final = np.array(seqs_final)
        group = np.array([self.group_name] * n)
        self.seqs_final = np.vstack((seqs_final, group_name)).T
        return seqs_final



rng = np.random.default_rng()
seqs_raw = rng.choice(bases, size=(n,10))
seperator=""
seqs_final=[]
for seq in seqs_raw:
    seqs_final.append(seperator.join(seq))
seqs_final = np.array(seqs_final)
group=np.array(["Random"]*n)
seqs_final=np.vstack((seqs_final,group)).T
print(seqs_final)
np.savetxt("toy_dataset.txt",seqs_final, delimiter="\t", newline="\n", fmt='%s')

#df=pd.DataFrame(seqs_final)
