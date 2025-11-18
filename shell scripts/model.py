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
from src.model_functions import *
sys.path.append("..")
from src.KLD_calculation import jsd_loop,markov_matrix

# set device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
scaler = StandardScaler()
# load data
testset1=pd.read_csv('../data/random_test.csv')
testset2=pd.read_csv('../data/human_test.csv')
testseqs1=one_hot(list(testset1["utr"]))
testseqs2=one_hot(list(testset2["utr"]))

#Daten für loop
sorted_random_seqs=pd.read_csv('../data/JSD_random_seqs.csv')#<-Speicherort
human_train_seqs=pd.read_csv('../data/human_train.csv')#<-Speicherort
full_dataset=pd.concat([sorted_random_seqs,human_train_seqs])
fullmrl = full_dataset['rl'].values.reshape([-1,1])
scaler.fit(fullmrl)
correct_random_mrl=np.array(testset1["rl"])
correct_human_mrl=np.array(testset2["rl"])
batch_size=128

#Validierungsset
human_dataset=pd.read_csv('../data/human_val.csv')
humanseqs1=one_hot(list(human_dataset["utr"]))
print(humanseqs1.shape)
humanmrl1=human_dataset['rl'].values.reshape([-1,1])
scaledhumanmrl1=torch.tensor(scaler.transform(humanmrl1).reshape([-1]),dtype = torch.float32)
human_val=list(zip(humanseqs1,scaledhumanmrl1))

valloader=torch.utils.data.DataLoader(human_val,batch_size=batch_size,shuffle=True)

# neural network model
class Model(nn.Module):

    def __init__(self):
        super(Model, self).__init__()
        self.conv1 = nn.Conv1d(4, 120, kernel_size=8, padding='same')
        self.conv2 = nn.Conv1d(120, 120, kernel_size=8, padding='same')
        self.conv3 = nn.Conv1d(120, 120, kernel_size=8, padding='same')
        self.flat1 = nn.Flatten()
        self.fc1 = nn.Linear(6000, 40)
        self.drop1 = nn.Dropout(0.2)
        self.out = nn.Linear(40, 1)

    def forward(self, x):
        x = F.relu(self.conv1(x))
        #print(x.shape)
        x = F.relu(self.conv2(x))
        #print(x.shape)
        x = F.relu(self.conv3(x))
        #print(x.shape)
        x = self.flat1(x)
        #x= torch.transpose(x, 1, 2)
        x = F.relu(self.fc1(x))
        x = self.drop1(x)
        x = self.out(x)
        return x


net = Model()
net.to(device)

criterion = nn.MSELoss()  # (aL-y)^2
optimizer = optim.Adam(net.parameters(), lr=0.001, betas=(0.9, 0.999))  # epsilon ist standradmäßig bei 1e-8

def train_one_epoch(epoch, trainloader, valloader):
    net.train(True)
    global loss_plot
    global val_plot
    global batch_plot
    global val_batch_plot
    global es_min_loss
    global min_batch_plot
    global stopcondition

    running_loss = 0.0
    val_running_loss = 0.0
    for batch_index, data in enumerate(trainloader):
        inputs, correct_mrl = data[0].to(device), data[1].to(device)
        optimizer.zero_grad()  # flush gradient
        outputs = net(inputs)  # shape:
        outputs = torch.reshape(outputs, (-1,))
        loss = criterion(outputs, correct_mrl)
        loss_plot.append(loss.item())
        batch_plot.append(batch_index + 1 + epoch * len(trainloader))
        loss.backward()
        optimizer.step()

        for val_batch_index, val_data in enumerate(valloader):
            val_inputs, val_correct_mrl = val_data[0].to(device), val_data[1].to(device)
            val_outputs = net(val_inputs)
            val_outputs = torch.reshape(val_outputs, (-1,))
            val_loss = criterion(val_outputs, val_correct_mrl)
            val_running_loss += val_loss.item()
        val_mean_loss = val_running_loss / len(valloader)
        val_plot.append(val_mean_loss)
        val_batch_plot.append(batch_index + 1 + epoch * len(trainloader))
        val_running_loss = 0.0
        tmp_min = es_min_loss
        if earlystopper(val_mean_loss, patience=len(trainloader), epsilon=1e-5):
            print(f'Stopping at epoch {epoch + 1}, batch {batch_index + 1}')
            stopcondition = False
            break
        if es_min_loss != tmp_min:
            min_batch_plot.append(batch_index + 1 + epoch * len(trainloader))
            minimum_plot.append(es_min_loss)


r2_rt_plot = []
r2_ht_plot = []

# New training loop
for i in np.arange(0, 1.1, 0.1):
    r2_rt_set = []
    r2_ht_set = []
    for j in range(10):
        net = Model().to(device)
        net.load_state_dict(torch.load("../models/untrained.pt", weights_only=False))
        optimizer = optim.Adam(net.parameters(), lr=0.001, betas=(0.9, 0.999))  # fresh optimizer
        if i <= 0.0:
            nextstep_dataset = human_train_seqs
        elif i == 1.0:
            nextstep_dataset = pd.concat([human_train_seqs, sorted_random_seqs.loc[0:len(sorted_random_seqs) - 1]],
                                         ignore_index=True)
        else:
            nextstep_dataset = pd.concat(
                [human_train_seqs, sorted_random_seqs.loc[0:int(i * len(sorted_random_seqs) - 1)]], ignore_index=True)
        nextstep_seqs = one_hot(list(nextstep_dataset['utr']))
        print(f'shape of trainingset: {nextstep_seqs.shape}')
        nexstep_mrl = nextstep_dataset['rl'].values.reshape([-1, 1])
        nextstep_scaled_mrl = torch.tensor(scaler.transform(nexstep_mrl).reshape([-1]), dtype=torch.float32)
        nextstep_train = list(zip(nextstep_seqs, nextstep_scaled_mrl))
        trainloader = torch.utils.data.DataLoader(nextstep_train, batch_size=batch_size, shuffle=True)
        # Actual Training
        batch_plot = []
        loss_plot = []
        val_plot = []
        val_batch_plot = []
        es_stopcounter = 0
        es_min_loss = np.inf
        minimum_plot = []
        min_batch_plot = []
        stopcondition = True
        epoch = 0
        print(f'Dataset = Human(20000) + {int(i * 100)}% JSD sampled({int(i * len(sorted_random_seqs))})')
        while stopcondition:
            train_one_epoch(epoch, trainloader, valloader)
            epoch += 1

        torch.save(net.state_dict(), f'../models/models_run04/JSD_trained_h{int(i * 100)}.pt')  # <-Speicherort

        # Getting the outputs of network
        trained_scaled_random_mrl = net(testseqs1.to(device)).cpu().detach().numpy().reshape(-1, )
        trained_scaled_human_mrl = net(testseqs2.to(device)).cpu().detach().numpy().reshape(-1, )
        trained_random_mrl = scaler.inverse_transform(trained_scaled_random_mrl.reshape([-1, 1])).reshape(-1, )
        trained_human_mrl = scaler.inverse_transform(trained_scaled_human_mrl.reshape([-1, 1])).reshape(-1, )

        # recording R^2
        hr2 = pearsonr(correct_human_mrl, trained_human_mrl)[0] ** 2
        print('R^2 of random seqs', hr2)
        r2_ht_set.append(hr2)
        rr2 = pearsonr(correct_random_mrl, trained_random_mrl)[0] ** 2
        print('R^2 of human seqs', rr2)
        r2_rt_set.append(rr2)
        if hr2 == np.nan or rr2 == np.nan or hr2 <= 0.2 or rr2 <= 0.2:
            pickle.dump([batch_plot, val_batch_plot, loss_plot, val_plot, minimum_plot, min_batch_plot],
                        open(f'../models/brokentraining_h{int(i * 100)}.pkl', 'wb'))  # <-Speicherort
            print(f'Model broke at {int(i * 100)}% of sampled seqs')

    r2_ht_plot.append(np.array(r2_ht_set))
    r2_rt_plot.append(np.array(r2_rt_set))
pickle.dump([r2_ht_plot,r2_rt_plot], open('../models/models_run04/r2_plot_sampled_data.pkl', 'wb'))