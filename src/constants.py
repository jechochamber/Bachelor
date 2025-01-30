"""This is a collection of constants that will be used across the package."""
import numpy as np
DNABASES = ["A", "C", "G", "T"]
RNABASES = ["A", "C", "G", "U"]
startcodons = ["AUG", "CUG", "GUG"]
stopcodons = ["UAA", "UGA", "UAG"]
btv={"A":np.array([1,0,0,0]),
     "C":np.array([0,1,0,0]),
     "G":np.array([0,0,1,0]),
     "U":np.array([0,0,0,1]),}# BaseToVector -> used in markov chain sequences