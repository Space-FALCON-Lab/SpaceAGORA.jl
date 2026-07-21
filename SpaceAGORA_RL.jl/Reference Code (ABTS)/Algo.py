# RL algorithms: DQN and simplified Actor-Critic

import torch
import torch.nn.functional as F
from torch.distributions import Categorical, Normal
import Model, Replay
import random
from misc.utils import AddBias, init
import numpy as np
import torch.nn.functional as Fnn
from collections import namedtuple

class DQN:
    def __init__(self, obs_space, act_space, lr=1e-4, replay_size=1000000, batch_size=32,
                 discount=0.99, target_update=2500, eps_decay=500000, train_start = 50000, device=None, path=None):
                 
        print('obs_space',obs_space.shape[0])
        self.obs_space, self.act_space = obs_space, act_space
        print(act_space)
        self.batch_size, self.discount, self.target_update = batch_size, discount, target_update
        self.device = device
        self.train_start = train_start
    

        if len(obs_space.shape) == 1:
            self.q_func = Model.TwoLayerFCNet(n_in=obs_space.shape[0], n_out=act_space.n).to(device)
            self.target_q_func = Model.TwoLayerFCNet(n_in=obs_space.shape[0], n_out=act_space.n).to(device)
            self.state_dtype = torch.float32
        elif len(obs_space.shape) > 1:
            self.q_func =Model.TwoLayerLSTMandFCNet(n_in=obs_space.shape[0], n_out=act_space.n, device=self.device).to(device)
            self.target_q_func =Model.TwoLayerLSTMandFCNet(n_in=obs_space.shape[0], n_out=act_space.n, device=self.device).to(device)
            self.state_dtype = torch.uint8
        else:
            print ("observation shape not supported:", obs_space.shape)
            raise
        self.q_func.train()
        self.target_q_func.train()

        print ('parameters to optimize:',
            [(name, p.shape, p.requires_grad) for name,p in self.q_func.named_parameters()],
            '\n')
        # self.optimizer = torch.optim.RMSprop(self.q_func.parameters(), lr=lr)
        self.optimizer = torch.optim.Adam(self.q_func.parameters(), lr=lr, betas=(0.9,0.999), eps=1e-6)

        # number of action steps done
        self.num_act_steps = 0
        #self.eps_start, self.eps_end, self.eps_decay = 1.0, 0.01, eps_decay
        self.eps_decay = eps_decay
        self.num_train_steps = 0
        self.double_q = True
        self.loss = 0

        self.replay = Replay.NaiveReplay(replay_size)

        if path:
            self.load(path)

    def compute_epsilon(self):

        ### <<< Your Code Here
        # linearly anneal from 1 to 0.01 in eps_decay steps
        eps = 1 + (0.01 - 1) * (self.num_act_steps-self.train_start) / self.eps_decay
        eps = max(0.01, eps)
        ### Your Code Ends >>>

        return eps

    def act(self, obses, test=False):
        # print(obses)
        if len(self.obs_space.shape) == 1:
            maxrate = obses[0][-1]*(0.5 - 0) + 0
        elif len(self.obs_space.shape) > 1:
            maxrate = max(obses[0][-1])*(0.5 - 0) + 0
        obses = torch.as_tensor(obses, device=self.device) #, dtype=self.state_dtype

        # if self.num_act_steps < self.eps_decay // 4:
        #     eps = 1.0 + (0.1 - 1.0) * self.num_act_steps / (self.eps_decay // 4)
        # elif self.num_act_steps < self.eps_decay:
        #     eps = 0.1 + (0.01 - 0.1) * self.num_act_steps / self.eps_decay
        # else:
        #     eps = 0.01
        eps = self.compute_epsilon()
        self.num_act_steps += 1

        with torch.no_grad():
            greedy_actions = self.q_func(obses).max(1)[1].tolist()
        actions = []
        for i in range(len(obses)):
            if random.random() < eps and not test and self.num_act_steps < self.train_start:
                #if maxrate<0.05:
                #    a = random.randrange(0,int((self.act_space.n-1)/2))
                #elif maxrate>0.25:
                #    a = random.randrange(int((self.act_space.n+1)/2),self.act_space.n)
                #else:
                a = random.randrange(0,self.act_space.n)
                # print('action random',a,'heatrate',maxrate)
            else:
                a = greedy_actions[i]
            actions.append(a)
        return actions

    def observe(self, obses, actions, transitions):
        for s,a,(sn,r,t,_) in zip(obses, actions, transitions):
            if t: # (s,a) leads to a terminal state
                sn = None
            self.replay.add((s,a,r,sn))
        if len(self.replay) > self.batch_size and self.replay.cur_batch is None:
            self.replay.sample_torch(self.batch_size, self.device)

    def train(self):
        #self.replay.get_current_batch(self.batch_size, self.device)
        state_batch, action_batch, reward_batch, non_terminal_mask, non_terminal_next_states = self.replay.cur_batch
        self.replay.cur_batch = None

        q_values = self.q_func(state_batch.float())[torch.arange(self.batch_size), action_batch]

        next_state_values = torch.zeros(self.batch_size, device=self.device)
        next_q_values = self.target_q_func(non_terminal_next_states.float()).detach()

        ### <<< Your Code Here
        
        if self.double_q:
            # Hint: try using self.q_func(non_terminal_next_states).argmax(1)
            next_state_values[non_terminal_mask] = next_q_values[
                torch.arange(len(next_q_values)), self.q_func(non_terminal_next_states.float()).argmax(1)]
        else:
            next_state_values[non_terminal_mask] = next_q_values.max(1)[0]
            
        target_q_values = next_state_values * self.discount + reward_batch



        self.loss = ((target_q_values - q_values)*(target_q_values - q_values)).mean()  # Hint: try using q_values, target_q_values
        
        ### Your Code Ends >>>

        self.optimizer.zero_grad()
        self.loss.backward()
        # torch.nn.utils.clip_grad_value_(self.q_func.parameters(), 1)
        torch.nn.utils.clip_grad_norm_(self.q_func.parameters(), 10)
        self.optimizer.step()
        # The ugly patch for using Adam on BlueWaters...
#        for group in self.optimizer.param_groups:
#            for p in group['params']:
#                state = self.optimizer.state[p]
#                if state['step'] >= 1022:
#                    state['step'] = 1022

        self.num_train_steps += 1
        if self.num_train_steps % self.target_update == 0:
            self.target_q_func.load_state_dict(self.q_func.state_dict())
        return self.loss.item()

    def test(self):
        #self.replay.get_current_batch(self.batch_size, self.device)
        state_batch, action_batch, reward_batch, non_terminal_mask, non_terminal_next_states = self.replay.cur_batch
        self.replay.cur_batch = None

        q_values = self.q_func(state_batch.float())[torch.arange(self.batch_size), action_batch]

        next_state_values = torch.zeros(self.batch_size, device=self.device)
        next_q_values = self.target_q_func(non_terminal_next_states.float()).detach()

        ### <<< Your Code Here
        
        if self.double_q:
            # Hint: try using self.q_func(non_terminal_next_states).argmax(1)
            next_state_values[non_terminal_mask] = next_q_values[
                torch.arange(len(next_q_values)), self.q_func(non_terminal_next_states.float()).argmax(1)]
        else:
            next_state_values[non_terminal_mask] = next_q_values.max(1)[0]
            
        target_q_values = next_state_values * self.discount + reward_batch



        self.loss = ((target_q_values - q_values)*(target_q_values - q_values)).mean()  # Hint: try using q_values, target_q_values
        
        ### Your Code Ends >>>
        return self.loss.item()

#    def save(self, path):
#        torch.save([self.q_func.state_dict(), self.target_q_func.state_dict(), self.optimizer.state_dict()], path)
    def save(self, path):
        torch.save({
            'num_act_steps': self.num_act_steps,
            'q_func_dict': self.q_func.state_dict(),
            'target_q_func': self.target_q_func.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'loss': self.loss,
            'replay': self.replay}, path)

#    def load(self, path):
#        s1, s2, s3 = torch.load(path, map_location=self.device)
#        self.q_func.load_state_dict(s1)
#        self.target_q_func.load_state_dict(s2)
#        self.optimizer.load_state_dict(s3)
    def load(self, path):
        print('in the load',path)
        checkpoint = torch.load(path, map_location=self.device)
        self.num_act_steps = checkpoint['num_act_steps']
        self.q_func.load_state_dict(checkpoint['q_func_dict'])
        self.target_q_func.load_state_dict(checkpoint['target_q_func'])
        self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        self.loss = checkpoint['loss']
        self.replay = checkpoint['replay']



# Buffer = namedtuple('Buffer', ('state', 'action', 'reward', 'terminal', 'last_state'))

class ActorCritic:
    def __init__(self, obs_space, act_space, lr=0.001, nproc=1, seg_len=16, discount=0.99,
                 entropy_coef=0.01, value_coef=0.5, device=None, shared_net=False):
        self.obs_space, self.act_space = obs_space, act_space
#        print(act_space)
        self.discount, self.entropy_coef, self.value_coef = discount, entropy_coef, value_coef
        self.nproc, self.seg_len = nproc, seg_len
        self.device = device

        if self.act_space.__class__.__name__ == "Discrete":
            num_outputs = act_space.n
            self.action_dtype = torch.uint8
        elif self.act_space.__class__.__name__ == "Box":
            num_outputs = act_space.shape[0]
            self.action_dtype = torch.float32

        self.logstd = AddBias(torch.zeros(num_outputs))
        # self.logstd = torch.nn.Parameter(torch.full((num_outputs,), 0.1))

        if len(obs_space.shape) == 1:
            shared_net = False
            self.actor = Model.TwoLayerFCNet(n_in=obs_space.shape[0], n_out=num_outputs).to(device)
            self.critic = Model.TwoLayerFCNet(n_in=obs_space.shape[0], n_out=1).to(device)
            self.state_dtype = torch.float32
        elif len(obs_space.shape) > 1:
#            print(obs_space.shape)
            if shared_net:
                self.actor_and_critic = Model.TwoLayerLSTMandFCNet(n_in=obs_space.shape[0], n_out=[num_outputs, 1], device=self.device).to(device)
                self.actor = lambda obses: self.actor_and_critic(obses, head=0)
                self.critic = lambda obses: self.actor_and_critic(obses, head=1)
            else:
                self.actor = Model.TwoLayerLSTMandFCNet(n_in=obs_space.shape[0], n_out=num_outputs, device=self.device).to(device)
                self.critic = Model.TwoLayerLSTMandFCNet(n_in=obs_space.shape[0], n_out=1, device=self.device).to(device)
            self.state_dtype = torch.uint8
        else:
            print ("observation shape not supported:", obs_space.shape)
            raise

        self.shared_net = shared_net
        # Other people like to use RMSprop for policy gradient... but neither of them is perfect
        # if shared_net:
        #     self.optimizer = torch.optim.RMSprop(self.actor_and_critic.parameters(), lr=lr, alpha=0.99, eps=1e-5)
        # else:
        #     self.optimizer = torch.optim.RMSprop(list(self.actor.parameters()) + list(self.critic.parameters()),
        #                                       lr=lr, alpha=0.99, eps=1e-5)

        if shared_net:
            self.optimizer = torch.optim.Adam(self.actor_and_critic.parameters(), lr=lr,
                betas=(0.9,0.99), eps=1e-6)
            print ('shared net = True, parameters to optimize:',
                [(name, p.shape, p.requires_grad) for name,p in self.actor_and_critic.named_parameters()],
                '\n')
        else:
            self.optimizer = torch.optim.Adam(list(self.actor.parameters()) + list(self.critic.parameters()), lr=lr,
                betas=(0.9,0.99), eps=1e-6)
            print ('shared net = False, parameters to optimize:',
                [(name, p.shape, p.requires_grad) for name,p in list(self.actor.named_parameters()) + list(self.critic.named_parameters())],
                '\n')


        self.states = torch.empty((nproc, seg_len) + obs_space.shape, dtype=self.state_dtype)
        self.actions = torch.empty((nproc, seg_len) + act_space.shape, dtype=self.action_dtype)
        self.rewards = torch.empty((nproc, seg_len), dtype=torch.float32)
        # self.advantages = torch.empty((nproc, seg_len), dtype=torch.float32)
        self.step = 0
        self.num_act_steps = 0
        self.last_reset = [0] * nproc
        self.loss = 0

    def act(self, obses, test=False):
        obses = torch.as_tensor(obses, device=self.device, dtype=self.state_dtype)
        if self.act_space.__class__.__name__ == "Discrete":
            dist = Categorical(logits=self.actor(obses))
            actions = dist.sample()
        elif self.act_space.__class__.__name__ == "Box":
            action_mean = self.actor(obses) #torch.tanh
            print(action_mean)
            zeros = torch.zeros(action_mean.size())
            action_logstd = self.logstd(zeros)

            dist = Normal(action_mean,action_logstd.exp())
            # stds = torch.clamp(self.logstd.exp(), 1e-3, 10)
            # dist = Normal(action_mean,stds)
            actions = dist.sample()[0]
            actions = np.clip(actions, self.act_space.low.min(), self.act_space.high.max())

        #dist = Categorical(self.actor(obses))
        # actions = dist.sample()
        # print(actions)
        self.num_act_steps += 1
        return actions.tolist()

    def observe(self, obses, actions, transitions):
        # print('actions in observe', actions)
        nproc = len(obses)
        # self.states[:, self.step] = obses
        # self.actions[:, self.step] = actions
        for i,(sn,r,terminal,_) in enumerate(transitions):
            # print('transitions',transitions)
            self.states[i, self.step] = torch.as_tensor(obses[i], dtype=self.state_dtype)
            self.actions[i, self.step] = torch.as_tensor(actions[i], dtype=self.action_dtype)
            # print('action in loop',self.actions[i, self.step],actions[i])
            self.rewards[i, self.step] = r
            R = 0.0
            if self.step == self.seg_len-1 and not terminal:
                R = self.critic(torch.as_tensor([sn], dtype=self.state_dtype, device=self.device)).item()
#                exit()
                terminal = True
            if terminal:
                # compute discounted reward backward along the trajectory
                for t in range(self.step, self.last_reset[i]-1, -1):
                    R = self.rewards[i, t] + self.discount * R
                    self.rewards[i, t] = R
                self.last_reset[i] = (self.step + 1) % self.seg_len
        self.step = (self.step + 1) % self.seg_len

    def train(self):
        torch.set_printoptions(precision=10)
        assert self.step == 0
        # print(self.actions)
        states = self.states.view(-1, *self.obs_space.shape).to(self.device)
        actions = self.actions.view(-1, *self.act_space.shape).to(self.device)

        # print('action in train',actions)
        rewards = self.rewards.view(-1).to(self.device)

        values = self.critic(states).squeeze(-1)
        advantage = rewards - values
        # advantage_detach = torch.clamp(advantage.detach(), -1, 1)
        # advantage = advantage.clamp(-2, 2).detach() + advantage - advantage.detach()
        # advantage = (advantage - advantage.mean().detach())/ (1.0 + advantage.std().detach())
        advantage = advantage / (1.0 + advantage.std().detach())
        advantage_detach = advantage.detach() #* self.advantage_normalizer

        # compute losses
        # probs = Fnn.softmax(self.actor(states))
        if self.act_space.__class__.__name__ == "Discrete":
            dist = Categorical(logits=self.actor(states))
        elif self.act_space.__class__.__name__ == "Box":
            action_mean = self.actor(states)
            zeros = torch.zeros(action_mean.size())
            action_logstd = self.logstd(zeros)
            dist = Normal(action_mean,action_logstd.exp())
        # dist = Categorical(probs)

        ### <<< Your Code Here
        # use 'advantage_detach', 'dist' and 'actions' to compute this; double-check the sign of your expression!
        policy_loss =(dist.log_prob(actions)*advantage_detach).mean()
       
        # use the entropy function of 'dist' to compute this
        entropy_loss = dist.entropy().mean()

        # value loss is the squared loss on variable 'advantage'
        value_loss = ((advantage)*(advantage)).mean()
        ### Your Code Ends >>>

        self.loss = -policy_loss - self.entropy_coef * entropy_loss + self.value_coef * value_loss

        # print (policy_loss.item(), entropy_loss.item(), value_loss.item())

        # do gradient descent
        self.optimizer.zero_grad()
        self.loss.backward()
        if self.shared_net:
            torch.nn.utils.clip_grad_norm_(self.actor_and_critic.parameters(), 0.5)
        else:
            torch.nn.utils.clip_grad_norm_(self.actor.parameters(), 0.5)
            torch.nn.utils.clip_grad_norm_(self.critic.parameters(), 0.5)
        self.optimizer.step()
#        print(self.optimizer.state)
#        # The ugly patch for using Adam on BlueWaters...
#        for group in self.optimizer.param_groups:
#            print(self.optimizer.state)
#            for p in group['params']:
#                print(self.optimizer.state[p])
#                state = self.optimizer.state[p]
#                if state['step'] >= 1022:
#                    state['step'] = 1022

        return self.loss.item()
    
    def test(self):
        states = self.states.view(-1, *self.obs_space.shape).to(self.device)
        actions = self.actions.view(-1, *self.act_space.shape).to(self.device)
#        print(actions)
        rewards = self.rewards.view(-1).to(self.device)

        values = self.critic(states).squeeze(-1)
        advantage = rewards - values
        # advantage_detach = torch.clamp(advantage.detach(), -1, 1)
        # advantage = advantage.clamp(-2, 2).detach() + advantage - advantage.detach()
        # advantage = (advantage - advantage.mean().detach())/ (1.0 + advantage.std().detach())
        advantage = advantage / (1.0 + advantage.std().detach())
        advantage_detach = advantage.detach() #* self.advantage_normalizer

        # compute losses
        # probs = Fnn.softmax(self.actor(states))
        if self.act_space.__class__.__name__ == "Discrete":
            dist = Categorical(logits=self.actor(states))
        elif self.act_space.__class__.__name__ == "Box":
            action_mean = self.actor(states)
            zeros = torch.zeros(action_mean.size())
            action_logstd = self.logstd(zeros)
            dist = Normal(action_mean,action_logstd.exp())

        ### <<< Your Code Here
        # use 'advantage_detach', 'dist' and 'actions' to compute this; double-check the sign of your expression!
        policy_loss =(dist.log_prob(actions)*advantage_detach).mean()
       
        # use the entropy function of 'dist' to compute this
        entropy_loss = dist.entropy().mean()

        # value loss is the squared loss on variable 'advantage'
        value_loss = ((advantage)*(advantage)).mean()
        ### Your Code Ends >>>

        self.loss = -policy_loss - self.entropy_coef * entropy_loss + self.value_coef * value_loss

        # print (policy_loss.item(), entropy_loss.item(), value_loss.item())

        
        return self.loss.item()


    def save(self, path):
        if self.shared_net:
            torch.save({'num_act_steps': self.num_act_steps,
                'actor_and_critic_dict': self.actor_and_critic.state_dict(),
                'optimizer_state_dict': self.optimizer.state_dict(),
                'loss': self.loss}, path)
#            torch.save([self.actor_and_critic.state_dict(), self.optimizer.state_dict()], path)
        else:
            torch.save({'num_act_steps': self.num_act_steps,
                'actor_dict': self.actor.state_dict(),
                'critic_dict': self.critic.state_dict(),
                'optimizer_state_dict': self.optimizer.state_dict(),
                'loss': self.loss}, path)
#            torch.save([self.actor.state_dict(), self.critic.state_dict(), self.optimizer.state_dict()], path)

    def load(self, path):
        checkpoint = torch.load(path, map_location=self.device)
        self.num_act_steps = checkpoint('num_act_steps')
        self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        self.loss = checkpoint['loss']
        if self.shared_net:
            self.actor_and_critic.state_dict(checkpoint['actor_and_critic_dict'])
        else:
            self.actor.load_state_dict(checkpoint['actor_dict'])
            self.critic.load_state_dict(checkpoint['critic_dict'])
            
#        states = torch.load(path, map_location=self.device)
#        if self.shared_net:
#            self.actor_and_critic.load_state_dict(states[0])
#        else:
#            self.actor.load_state_dict(states[0])
#            self.critic.load_state_dict(states[1])
#        self.optimizer.load_state_dict(states[-1])
