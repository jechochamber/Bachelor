import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import torch.nn.functional as F
import torch.optim as optim
import pickle
import numpy as np
import pandas as pd
import sys

sys.path.append("..")


# encoding of data
def one_hot(seqs):
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


# set device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# load data
dataset1 = pd.read_csv('../data/random_train_pc.csv')
trainseqs1 = one_hot(list(dataset1["utr"]))
trainmrl1 = torch.tensor(np.array(dataset1["rl"]), dtype=torch.float32)
dataset1_reshaped = list(zip(trainseqs1, trainmrl1))
batch_size = 128

trainset, valset = torch.utils.data.random_split(dataset1_reshaped, [int(len(dataset1_reshaped) * 0.975),
                                                                     int(len(dataset1_reshaped) * 0.025)])

trainloader = torch.utils.data.DataLoader(trainset, batch_size=batch_size, shuffle=True)
valloader = torch.utils.data.DataLoader(valset, batch_size=batch_size, shuffle=False)


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
scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=5, factor=0.5)

def train_one_epoch(epoch):
    net.train(True)
    global loss_plot
    global val_plot
    global batch_plot
    global val_batch_plot

    running_loss = 0.0
    val_running_loss = 0.0
    for batch_index, data in enumerate(trainloader):
        inputs, correct_mrl = data[0].to(device), data[1].to(device)
        optimizer.zero_grad()  # flush gradient
        outputs = net(inputs)  # shape:
        outputs = torch.reshape(outputs, (-1,))
        loss = criterion(outputs, correct_mrl)
        loss_plot.append(loss.item())
        batch_plot.append(batch_index+1)
        loss.backward()
        optimizer.step()

        if (batch_index+1) % 3 == 2:
            for val_batch_index, val_data in enumerate(valloader):
                val_inputs, val_correct_mrl = val_data[0].to(device), val_data[1].to(device)
                val_outputs = net(val_inputs)
                val_outputs = torch.reshape(val_outputs, (-1,))
                val_loss = criterion(val_outputs, val_correct_mrl)
                val_running_loss += val_loss.item()
            val_plot.append(val_running_loss/len(valloader))
            val_batch_plot.append(batch_index+1)
            scheduler.step(val_running_loss / len(valloader))
            val_running_loss = 0.0






batch_plot=[]
loss_plot=[]
val_plot=[]
val_batch_plot=[]
for epoch_index in range(9):
    print(f'Epoch {epoch_index+1} of 9')
    train_one_epoch(epoch_index)
pickle.dump([batch_plot,val_batch_plot,loss_plot,val_plot], open('../data/plot_vals.pkl', 'wb'))
torch.save(net.state_dict(),'../data/random_unval.pt')