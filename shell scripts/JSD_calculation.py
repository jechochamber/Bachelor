from src.KLD_calculation import *

reference_seqs = pd.read_csv("../data/human_train.csv")

for i in reference_seqs.index:
    reference_seqs.loc[i, 'utr'] = reference_seqs.loc[i, 'utr'].replace('T', 'U')
a=1

jsd_list=[]

for seq in reference_seqs["utr"]:
    p_matrix = markov_matrix(seq, a=a)[2]
    jsd_list.append(p_matrix[0][0])


pickle.dump(jsd_list, open("../data/human_jsd.pkl", "wb"))
jsd_arr=np.array(jsd_list)
np.save("../data/human_jsd.npy", jsd_arr)
print("done")

reference_data=pd.read_csv("../data/human_sample.csv")
q_data=pd.read_csv("../data/random_train.csv")
results=distance_calculation(reference_data,q_data,'utr','utr',50)

pickle.dump(results, open("../data/human_distance.pkl", "wb"))
