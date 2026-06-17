# Policy and value function (network) models

import torch, torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
import numpy as np



def init(module, weight_init, bias_init, gain=1):
    weight_init(module.weight.data, gain=gain)
    bias_init(module.bias.data)
    return module

def hidden_init(layer):
    """
    Used for parameter initialization
    """
    fan_in = layer.weight.data.size()[0]
    lim = 1. / np.sqrt(fan_in)
    return (-lim, lim)

class TwoLayerFCNet(nn.Module):
    def __init__(self, n_in=4, n_hidden=512, n_out=2):
        super().__init__()
        init_ = lambda m: init(m, nn.init.orthogonal_, lambda x: nn.init.constant_(x, 0), np.sqrt(2))
        
        ### <<< Your Code Here
        self.fc1 = init_(nn.Linear(n_in, 2*n_hidden))  #change with 2*n_in, n_hidden, n_hidden/2 remember to change weight initialization
        #self.fc2 = init_(nn.Linear(2*n_hidden, 2*n_hidden))
        self.fc2 = init_(nn.Linear(2*n_hidden, 2*n_hidden))
        self.dropout = nn.Dropout2d(p=0.4)
        self.fc3 = init_(nn.Linear(2*n_hidden, n_out))

        # self.fc1.weight.data.uniform_(hidden_init(self.fc1))
        # self.fc2.weight.data.uniform_(2*hidden_init(self.fc2))
        # self.fc3.weight.data.uniform_(-3e-3, 3e-3)
        ### Your Code Ends >>>

    def forward(self, x):
#        print(x)
        x = x.float()
        ### <<< Your Code Here
        x = F.relu(self.fc1(x))#torch.tanh(self.fc1(x))
        #x = F.relu(self.fc2(x))
        #x = self.dropout(x)
        x = F.relu(self.fc2(x))#torch.tanh(self.fc2(x))
        ### Your Code Ends >>>

        return self.fc3(x)
        
class TwoLayerLSTMandFCNet(nn.Module):

    def __init__(self, n_in=4, n_hidden=512, n_layers = 1, n_out=1, device=None):
        super().__init__()
        init_ = lambda m: init(m, nn.init.orthogonal_, lambda x: nn.init.constant_(x, 0), np.sqrt(2))

        self.n_layers = n_layers
        self.n_hidden = n_hidden
        self.device = device
        
        print('n_in',n_in)
        self.LSTM = nn.LSTM(input_size=n_in, hidden_size = n_hidden, num_layers = n_layers)
        self.bn_lstm= nn.BatchNorm1d(n_hidden)
        # self.dropout = nn.Dropout(0.2)
        # self.fc1 = init_(nn.Linear(n_hidden, n_hidden))  #change wuth 2*n_in, n_hidden, n_hidden/2 remember to change weight initialization
        # self.fc3 = init_(nn.Linear(n_hidden, n_out))
        self.fc = init_(nn.Linear(n_hidden, n_out))
            
    
    def _tensor_creation(self, x):
        temp = []
        
        for i in range(x.size(2)):
            _x = x[:,:,i].tolist()
            
            temp.append(_x)
            
        return torch.Tensor(temp)

    def forward(self, x):
        x = x.float()
        no_of_timesteps = x.shape[-1]
#        x = self._tensor_creation(x)
        
        x = x.permute(2,0,1)

        h_0 = Variable(torch.zeros(
            self.n_layers, x.size(1), self.n_hidden,device=self.device))
        
        c_0 = Variable(torch.zeros(
            self.n_layers, x.size(1), self.n_hidden,device=self.device))

        # Propagate input through LSTM
        ula, (h_out, c_0) = self.LSTM(x, (h_0, c_0))
        x = h_out[-1]
        # x = self.dropout(x)
        x = self.fc(x)
        return x
    
    
        

class SimpleCNN(nn.Module):
    def __init__(self, n_in=[4,84,84], conv_channels=[32, 64, 64],
                 conv_kernels=[8, 4, 3], conv_strides=[4, 2, 1], n_fc=[256], n_out=4):
        super().__init__()

        self.conv_layers = []
        c0 = n_in[0]
        h0 = n_in[1]
        assert n_in[1] == n_in[2], 'input must be square image'
        for c, k, s in zip(conv_channels, conv_kernels, conv_strides):
            
            ### <<< Your Code Here
            # append nn.Conv2d with kernel size k, stride s
            self.conv_layers.append(nn.Conv2d(c0,c,k,stride=s))
            # append nn.ReLU layer
            self.conv_layers.append(nn.ReLU())
            ### Your Code Ends >>>

            h0 = int(float(h0-k) / s + 1)
            c0 = c
        self.conv_layers = nn.Sequential(*self.conv_layers)

        self.fc_layers = []
        h0 = h0 * h0 * conv_channels[-1]
        for i, h in enumerate(n_fc):
            # append Linear and ReLU layers
            ### <<< Your Code Here:
            self.fc_layers.append((nn.Linear(h0,h)))
            ### Your Code Ends >>>
            
            self.fc_layers.append( nn.ReLU() )
            h0 = h
        if type(n_out) is list:
            self.out_layers = nn.ModuleList([ nn.Linear(h, o) for o in n_out ])
        else:
            self.fc_layers.append( nn.Linear(h, n_out) )
        self.fc_layers = nn.Sequential(*self.fc_layers)


    def forward(self, x, head=None):
        x = x.float() / 256

        ### <<< Your Code Here:
        # feed x into the self.conv_layers
        x = self.conv_layers(x)
        # (flatten) reshape x into a batch of vectors
        x = x.view(x.size(0), -1) 
        # feed x into the self.fc_layers
        x = self.fc_layers(x)
        ### Your Code Ends >>>
        if head is not None:
            x = self.out_layers[head](x)
        return x

class DuelQNet(nn.Module):
    def __init__(self):
        pass
