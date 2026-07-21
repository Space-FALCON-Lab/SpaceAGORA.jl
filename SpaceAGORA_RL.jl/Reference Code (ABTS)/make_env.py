"""Aerobraking task"""
import numpy as np
from subprocess import Popen, PIPE
# from gym import core, spaces
import sys

# directory = '/home/gfalcon2/scratch/Aerobraking/ABTS/'
# directory2 = '/home/gfalcon2/scratch/DRLAerobraking/DRLModel/'

directory = '/Users/Josephine/Research/Aerobraking/ABTS/'
directory2 = '/Users/Josephine/Research/DRLAerobraking/DRLModel/'
# directory = '/content/drive/Othercomputers/My MacBook Pro/Research/Aerobraking/ABTS/'
# directory2 = '/content/drive/Othercomputers/My MacBook Pro/Research/DRLAerobraking/DRLModel/'

sys.path.insert(1, directory)

from utils.Reference_system import *
import config
from physical_models.Density_models import density_exp as dm
import spaces
import pandas as pd
import math
from physics import *
from datetime import *
from misc import seeding

import os
import time

from call_aerobraking import *

class AerobrakingEnv():

    """
    Aerobraking is a mission maneuver that decrease the orbital period through 
    several passages into the atmosphere of a planet. The goal is to reach a 
    target orbit without incurring into thermal violations (heat rate > 0.25 W/cm2
    and heat rate < 0.05 W/cm^2). This can be acieved through a periapsis location
    control. The periapsis location is modified through the ABM maneuver which is 
    an impulsive maneuver performed at the apoapsis.
    **STATE:**
    The state consists of the orbital elements (periapsis altitude, apoapsis radius,
    INC, AOP, RAN), drag passage time of flight, maximum heat rate and atmospheric
    density.
    **ACTIONS:**
    The action is either performing a periapsis raise maneuver (+DeltaV) or a
    periapsis lowering maneuver (-DeltaV). The action space is [-1.5, -1, -0.5,
    -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 1.5].
    .. note::
        The drag passage dynamics is performed through the ABTS simulator, while
        the in-space orbit dynamics and the ABM is calculated here.
    **REFERENCE:**
    .. seealso::
        R. Sutton: Generalization in Reinforcement Learning:
        Successful Examples Using Sparse Coarse Coding (NIPS 1996)
    .. seealso::
        R. Sutton and A. G. Barto:
        Reinforcement learning: An introduction.
        Cambridge: MIT press, 1998.
    """



    def __init__(self,args):
        self.max_ra = 30000*1e3 #m
        self.min_ra = 3396*1e3 #m
        self.max_hp = 140*1e3 #m
        self.min_hp = 80*1e3 #m
        self.max_hr = np.inf
        self.min_hr = 0.0
        self.max_hl = np.full((100),np.inf)
        self.min_hl =  np.zeros(100)
        self.min_density = np.zeros(100)
        self.min_velocity = np.zeros(100)
        self.min_position = np.zeros(100)
        self.max_density = np.full((100),np.inf)
        self.max_velocity =  np.full((100),np.inf)
        self.max_position =  np.full((100),np.inf)
        self.max_angle = np.full((100), np.pi*2)
        self.min_angle = np.zeros(100)
        self.min_heatrate =  np.zeros(100)
        self.max_heatrate =  np.full((100),np.inf)
        self.min_sma =  np.zeros(100)
        self.max_sma =  np.full((100),np.inf)
        self.min_ecc =  np.zeros(100)
        self.max_ecc  =  np.full((100),np.inf)
        self.max_time = np.inf
        self.min_time = 0.0
        self.max_deltav = 1.5
        self.min_deltav = -1.5
        self.min_target = 3500*1e3
        self.max_target = 7000*1e3
        self.count = 0

        self.phase = args.phase
        if self.phase == 'Walkout':
            self.finalstate = 3900 * 1e3
            self.r_norm = 5000*10**3
        elif self.phase == 'Endgame':
            self.finalstate = 4906 *1e3
            self.r_norm = 10100*10**3
        elif self.phase == 'Main':
            self.finalstate = 3900 *1e3
            self.r_norm = 10100*10**3
        elif self.phase == 'Campaign':
            self.finalstate = 3900 *1e3
            self.r_norm = 30000*10**3

        class p():
            def __init__(self):
                p.Rp_e = 3.3962e6
                p.J2 = 1.96045e-3
                p.mu = 4.2828e13  # gravitational parameter, m^3/s^2
                p.rho_ref = 8.748923102971180e-07  # 3.8*10**-8#7.3*10**-8#0.02  # kg/m^3
                p.h_ref = 90 * 1e3  # 103*10**3#109 * 10 ** 3
                p.H = 6.308278108290950e+03  # 10.6 * 10 ** 3  # m
                p.R = 188.92  # J/KgK
                p.gamma = 1.33  # wrong check
                p.T = 150  # K constant # wrong anchor # the script is not considering this number but T =150
                p.g_ref = 3.71  # m/s^2
        self.p = p()
        
        self.seed()
        self.args = args
        self.output = args.output
        self.directory = directory2+args.data_dir
        self.nominal = args.nominal


        
#        low = np.concatenate([[self.min_target],[self.min_hr],[self.min_hr],[self.min_ra],[self.min_hp]])
        # high=np.concatenate([[self.max_target], self.max_hl.flat, self.max_sma.flat, self.max_ecc.flat, self.max_density.flat, self.max_velocity.flat, self.max_position.flat, self.max_heatrate.flat],dtype=object)
        # low = np.concatenate([[self.min_target], self.min_hl.flat,  self.min_sma.flat, self.min_ecc.flat, self.min_density.flat, self.min_velocity.flat,self.min_position.flat, self.min_heatrate.flat],dtype=object)

        if self.output == 'NN':
            high = np.array([2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0],dtype=np.float32,)

            low = np.array([0.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],dtype=np.float32,)
            self.observation_space = spaces.Box(low=low, high=high, dtype=np.float32)#np.float32)
        else:
            high = np.array([self.max_position,self.max_position,self.max_velocity,self.max_angle,self.max_angle,self.max_angle,self.max_angle, self.max_density, self.max_heatrate])
            low = np.array([self.min_position,self.max_position,self.min_velocity,self.min_angle, self.min_angle, self.min_angle, self.min_angle,self.min_density, self.min_heatrate])
            self.observation_space = spaces.Box(low=low, high=high)
        # self.action_space = spaces.Box(low=self.min_deltav,high=self.max_deltav,shape=(21,),dtype='float32')
        if self.args.algo == 'dqn':
            self.action_space = spaces.Discrete(13)
        else:
            self.action_space = spaces.Box(
                low=self.min_deltav, high=self.max_deltav, shape=(1,), dtype=np.float32)


    def seed(self, seed=None):
        self.np_random, self.seed = seeding.np_random(seed)
        return [self.seed]


    def get_initialstate(self,pid,sim):
        if self.phase == 'Endgame':
            self.ra = 6000* 1e3+self.np_random.uniform(low=(-2.5)* 1e3, high=(2.5)*1e3)#self.ra = self.np_random.uniform(low=(6023.95)* 1e3, high=(5000)*1e3)
        elif self.phase == 'Walkout':
            self.ra = 4906*1e3+self.np_random.uniform(low=(-2.5)* 1e3, high=(2.5)*1e3)
        elif self.phase == 'Main':
            self.ra = 10038*1e3+self.np_random.uniform(low=(-2.5)* 1e3, high=(2.5)*1e3) #10038
        elif self.phase == 'Campaign':
            self.ra = 28533*1e3+self.np_random.uniform(low=(-2.5)* 1e3, high=(2.5)*1e3)
        self.rai = self.ra
        self.raf = self.finalstate#self.np_random.uniform(low=4806 * 1e3, high=5000 * 1e3)
        self.vi = math.pi
        if self.nominal == 0:
            self.incl = self.np_random.uniform(low=88.6, high=98.6)
            if self.phase == 'Walkout':
                self.hp = 85 * 1e3 + self.np_random.uniform(low=(-1.5) * 1e3, high=(1.5) * 1e3)  # self.np_random.uniform(low=108 * 1e3, high=120 * 1e3)
                self.hp_prev = self.hp
                self.OMEGA = self.np_random.uniform(low=110, high=120)
                self.omega = self.np_random.uniform(low=60, high=89)
                self.year = 2001#self.np_random.randint(low=2001, high=2021)
                self.month = 12#self.np_random.randint(low=1, high=12)
            elif self.phase == 'Campaign':
                self.hp = 89 * 1e3 + self.np_random.uniform(low=(-1.5) * 1e3, high=(1.5) * 1e3)  # self.np_random.uniform(low=108 * 1e3, high=120 * 1e3)
                self.hp_prev = self.hp
                self.OMEGA = self.np_random.uniform(low=110, high=120)
                self.omega = self.np_random.uniform(low=91, high=130)
                self.year = 2001#self.np_random.randint(low=2001, high=2021)
                self.month = self.np_random.randint(low=10, high=12)
            else:
                self.hp = 92 * 1e3 + self.np_random.uniform(low=(-1.5) * 1e3, high=(1.5) * 1e3)  # self.np_random.uniform(low=108 * 1e3, high=120 * 1e3)
                self.hp_prev = self.hp
                self.OMEGA = self.np_random.uniform(low=110, high=120)
                self.omega = self.np_random.uniform(low=60, high=89)
                self.year = 2001#self.np_random.randint(low=2001, high=2021)
                self.month = 12#self.np_random.randint(low=1, high=12)
            self.day = self.np_random.randint(low=1, high=28)
            self.hour = self.np_random.randint(low=0, high=23)
            self.mins = self.np_random.randint(low=0, high=59)
            self.sec = self.np_random.randint(low=0, high=59)
            self.montecarlo = self.np_random.randint(low=0, high=1000)
            self.montecarlo_on = 1
        else:
            self.hp = 86 * 1e3 + self.np_random.uniform(low=(-1.) * 1e3, high=(1.) * 1e3)  # self.np_random.uniform(low=108 * 1e3, high=120 * 1e3)
            self.hp_prev = self.hp
            self.incl = 93.6+self.np_random.uniform(low=-0.25, high=0.25)
            self.OMEGA = 115+self.np_random.uniform(low=-0.25, high=0.25)
            self.omega = 66+self.np_random.uniform(low=-0.25, high=0.25)
            if self.phase == 'Endgame':
                self.year = 2001
                self.month = 12
                self.day = 26
            elif self.phase == 'Walkout':
                self.year = 2002
                self.month = 1
                self.day = 3
            elif self.phase == 'Main':
                self.year = 2001
                self.month = 12
                self.day = 18
                self.OMEGA = 115 + self.np_random.uniform(low=-0.25, high=0.25)
                self.omega = 89 + self.np_random.uniform(low=-0.25, high=0.25) #88
                self.hp = 92 * 1e3 + self.np_random.uniform(low=(-1.) * 1e3, high=(1.) * 1e3)
            elif self.phase == 'Campaign':
                self.year = 2001
                self.month = 11
                self.day = 6
                self.OMEGA = 114 + self.np_random.uniform(low=-0.25, high=0.25)
                self.omega = 109 + self.np_random.uniform(low=-0.25, high=0.25)
                self.hp = 87.5 * 1e3 + self.np_random.uniform(low=(-1.) * 1e3, high=(1.) * 1e3)  # self.np_random.uniform(low=108 * 1e3, high=120 * 1e3)
                self.hp_prev = self.hp
            self.hour = 0
            self.mins = 0
            self.sec = 0
            self.montecarlo = self.np_random.randint(low=0, high=1000)
            self.montecarlo_on = 1
        self.pid=pid
        self.time = 0.0
        self.deltav = 0.0
        print('initial hp', self.hp/1e3)
        if not sim:
            self.sim = 'Orbits'
        else:
            self.sim = 'Drag Passage'
        self.initialstate = [self.ra, self.hp, self.incl, self.OMEGA, self.omega, self.year, self.month, self.day, self.hour, self.mins,
                             self.sec,self.raf,self.montecarlo, self.montecarlo_on,self.pid,self.sim, self.directory]  # define initial conditions

        self.state = self.initialstate
        return

    def reset(self,sim,pid, rep=0):
        # print('Start New Episode')
        print('rep',rep)
        self.count = 0
        if rep == 0:
            self.get_initialstate(pid, sim)
        elif rep == 1:
            self.state = self.initialstate
            self.ra, self.hp, self.incl, self.OMEGA, self.omega, self.year, self.month, self.day, self.hour, self.mins, self.sec, self.raf, self.montecarlo, self.montecarlo_on, self.pid, self.sim, self.directory = self.initialstate
            self.hp_prev = self.hp
            self.rai, self.raf, self.vi = self.ra, self.finalstate, math.pi
            self.time, self.deltav  = 0.0, 0.0

        if self.output == 'NN':
            ra_i =  (self.ra - 3396200) / (self.r_norm - 3396200)
            obs = [0.0, ra_i, 1.0,0.0,0.0,0.0,0.0,0.0, 0.0]
        else:
            obs = np.array([np.full((100),0),np.full((100),1),np.full((100),1),np.full((100),0),np.full((100),0),np.full((100),0),np.full((100),0),np.full((100),0),np.full((100),0)])
        return obs#(observations, reward, terminal, {})#self._get_ob()

    def step(self, a, sim):
        # print(a)
        # if a != 0:
        #     print(a,np.clip(a, self.min_deltav, self.max_deltav)[0])
        self.count += 1
        self._set_action(action=a)

        self._env_setup(sim)
        self._get_newstate(sim)

        observations = self._get_ob(sim)
        terminal = self._terminal(observations,sim)
        reward = self._reward(observations,sim, self.action)
        observations = self._norm_obs(observations)
        # print(observations)

        return (observations, reward, terminal, {})

    def _norm_obs(self,obs):
        if self.output == 'NN':
            t = (np.array([(obs[0] - 400.0) / (1200.0 - 400.0)]))  # time normalized over between 600 and 800 secs
            time = (np.array([(obs[6] - 730486) / (731216 - 730486)]))
            ra = (np.array([(obs[1] - 3396200) / (self.r_norm - 3396200)]))
            hp = (np.array([(obs[2] - 85000) / (145000 - 85000)]))
            OMEGA = (np.array([(obs[4] - 1.57) / (3.14 - 1.57)]))
            inc = (np.array([(obs[5] - 1.39) / (1.75 - 1.39)]))
            omega = (np.array([(obs[3] - 0) / (3.14 - 0)]))
            rho = (np.array([(obs[7] - 0) / (1.5e-07 - 0)]))
            heat_rate = (np.array([(obs[8] - 0) / (0.5 - 0)]))
            print(np.concatenate((t, ra, hp, omega, OMEGA, inc, time, rho, heat_rate)))
            return np.concatenate((t, ra, hp, omega, OMEGA, inc,time, rho, heat_rate))
        elif self.output == 'LSTM':
            t = (np.array([(value - 400.0)/(1200.0-400.0) for value in obs[0]])) # time normalized over 60 days
            # time = (np.array([(obs[6] - 730486) / (731216 - 730486)]))
            ra = (np.array([(value - 3396200)/(self.r_norm-3396200) for value in obs[1]]))
            hp = (np.array([(value - 85000)/(145000-85000) for value in obs[2]]))
            omega =  (np.array([(value - 0)/(3.14-0) for value in obs[3]]))
            OMEGA =  (np.array([(value - 1.57)/(3.14-1.57) for value in obs[4]]))
            inc =  (np.array([(value - 1.39)/(1.75-1.39) for value in obs[5]]))
            vi =  (np.array([(value - 0)/(6.28-0) for value in obs[6]]))
            rho =  (np.array([(value - 0)/(4.5e-07-0) for value in obs[7]]))
            heat_rate =  (np.array([(value - 0)/(1-0) for value in obs[8]]))

            return np.array([t, ra,hp,omega, OMEGA, inc, vi,rho,heat_rate])

    def _terminal(self,observations,sim):
        if self.output == 'NN':
            ra = observations[1]
        elif self.output == 'LSTM':
            ra = observations[1][-1]
        if not self.args.test:
            temp = bool(ra<=(self.finalstate) or self._is_success(observations) or self._impact(observations) or self._outdragpassage(observations) or not self._under_thermal_loads(observations,sim))
        else:
            temp = bool(ra<=(self.finalstate) or self._is_success(observations) or self._impact(observations) or self._outdragpassage(observations))
        return temp

    def _time_mission(self):
        return self.time/(60*60*24)

    def _reward(self, observations, sim, action):
        a = action
        r = 0.0
        # reward goal
        if self._terminal(observations, sim):
            # penalties impact
            if self._impact(observations):
                r -= 6
            # penalties out of drag passage
            elif self._outdragpassage(observations):
                r -= 6
            elif self._is_success(observations):
                if self.output == 'NN':
                    ra = observations[1]
                elif self.output == 'LSTM':
                    ra = observations[1][-1]
                dist = (abs(ra - self.finalstate)/10**3)//1
                r += (-dist*2/5+10)
            elif self._noreached_goal(observations):
                r -= 4
                # if self._under_thermal_loads(observations, sim):


        # reward thermal loads
        if not self._under_thermal_loads(observations, sim):
            if self.output == 'NN':
                max_rate = observations[-1]
            elif self.output == 'LSTM':
                max_rate = max(observations[-1])
            if max_rate >0.3 and max_rate<=0.45:
                r -= 5
            elif max_rate >0.45:
                r -= 6
            else:
                r -= 4

        else:
            if self.output == 'NN':
                ra = observations[1]
                hr_cap = (observations[-1]-0.05)/(0.25-0.05)  ## for LSTM needs to be changed
            elif self.output == 'LSTM':
                ra = observations[1][-1]
                hr_cap = (max(observations[-1])-0.05)/(0.25-0.05)  ## for LSTM needs to be changed


            d_cap = (ra-self.finalstate)/(self.r_norm-self.finalstate)
            deltav_cap = abs(a[0])


            if hr_cap < 0:
                hr_cap = 0 # as it is in the middle of the corridor
                deltav_cap = abs(a[0]-self.max_deltav)

            if d_cap >= 0:
                r += (1-d_cap**(0.4))* (0.025/(0.1*(2*np.pi)**0.5)) *np.exp(-0.5*((hr_cap-0.5)/0.1)**2)+ (-0.1-deltav_cap*0.1)
                # r += (1 - d_cap ** (0.4)) * 0.05 / (0.5 + 25 * hr_cap ** 2) + (0 - deltav_cap * 0.1) - 0.1
                # r += 0.095*(1-d_cap**(0.4))* (0.08-1.77* hr_cap + 57.9 * hr_cap**2 -615.5 * hr_cap**3+3053.6 * hr_cap**4 - 7460.4 * hr_cap**5 + 9429.1 * hr_cap**6
                #        - 5950.6 * hr_cap**7 + 1487.6 * hr_cap**8)+ (0-deltav_cap*0.1)-0.1
            else:
                r += 0.0
            print('d_cap',d_cap,'hr_cap',hr_cap,'deltav_cap',deltav_cap,'r',r)

        return r

    def _get_ob(self, sim):
        """Returns the observation.
        """
        if self.output == 'NN':
            if sim != 2:
                load = np.array([float(max(self.data['heat_load'].tolist()))])
                heat_rate = np.array([float(max(self.data['heat_rate'].tolist()))])
                ra = np.array([float(self.data['a'].tolist()[-1]) *(1+float(self.data['e'].tolist()[-1]))])
                ra0 = np.array([float(self.data['a'].tolist()[0]) *(1+float(self.data['e'].tolist()[0]))])
                hp = np.array([float(min(self.data['alt'].tolist()))])
                temp = (self.data['rho'].tolist())
                self.time += float(self.data['time'].tolist()[-1])
                dragpassage_time = np.array([float(self.data['time'].tolist()[-1])])
                idx = np.round(np.linspace(0, len(temp) - 1, 100)).astype(int)

                rho = (np.array([float(max(self.data['rho'].tolist()))]))
                omega = (np.array([float(self.data['omega'].tolist()[-1])]))
                OMEGA = (np.array([float(self.data['OMEGA'].tolist()[-1])]))
                inc =  (np.array([float(self.data['i'].tolist()[-1])]))
                year = int(self.data['year'].tolist()[-1])
                month = int(self.data['month'].tolist()[-1])
                day = int(self.data['day'].tolist()[-1])
                time = (np.array([float(date(year,month,day).toordinal())]))
                sma =  (np.array(self.data['a'].tolist()))
                ecc = (np.array(self.data['e'].tolist()))
                vel = (np.array(self.data['vel_ii_mag'].tolist())[idx.astype(int)])
                pos = (np.array(self.data['pos_ii_mag'].tolist()))
                # heat_rate = (np.array(self.data['heat_rate'].tolist())[idx.astype(int)])
                time_man = np.array([self.time/(60*60*24)])
                deltav = np.array([self.deltav])
                target = np.array([self.raf])
                if sim ==0:
                    h = pos-self.p.Rp_e
                    # print('h',h)
                    # h = (self.data['alt'].tolist())
                    vi = (self.data['vi'].tolist())
                    # print('vi',vi)
                    # print('t',self.data['time'].tolist())
                    idx = (np.abs(np.asarray(vi) - 0)).argmin()
                    # print('idx',idx)
                    idx_h1 = (np.abs(np.asarray(h[0:idx])- 160*10**3)).argmin()
                    idx_h2 = (np.abs(np.asarray(h[idx:]) - 160*10**3)).argmin()
                    # print(idx_h1,idx_h2)
                    t = abs(float(self.data['time'].tolist()[idx+idx_h2])-float(self.data['time'].tolist()[idx_h1]))
                    dragpassage_time = np.array([t])
                    # print(t,dragpassage_time)

            else:
                target = np.array([self.raf])
                heat_rate = np.array([max(self.data['heat_rate'])])
                load = np.array([self.data['heat_load']])
                ra = np.multiply(self.data['a'],(1+self.data['e']))
                ra =np.array([ra[-1]])
                ra0 = np.array([ra[0]])
                omega = np.array([self.data['omega'][-1]])
                OMEGA = np.array([self.data['OMEGA'][-1]])
                inc = np.array([self.data['incl'][-1]])
                rho = np.array([max(self.data['rho'])])
                hp = np.multiply(self.data['a'],(1-self.data['e'])) - self.p.Rp_e
                hp = np.array([hp[-1]])
                self.time += float(self.data['time'][-1])
                dragpassage_time = np.array([float(self.data['time'][-1])])
                time_man = np.array([self.time/(60*60*24)])
            print('Passage',self.count, 'Ra0', round(ra0[0]*1e-3,3),'Ra', round(ra[0]*1e-3,3),
                  'Raf', round(self.raf*1e-3,3),'hp0',round(self.hp_prev*1e-3,3),'hp', round(hp[-1]*1e-3,3),
                  'action',self.action,'maxrate',round(heat_rate[0],3),
                  'dragpassagetime', round(dragpassage_time[0]), 'omega', round(math.degrees(omega[0]),3),'OMEGA',round(math.degrees(OMEGA[0]),3))

            self.hp_prev = hp[-1]
            return np.concatenate((dragpassage_time, ra, hp, omega, OMEGA, inc,time, rho, heat_rate))
                # (time_man,ra,hp,target,maxrate))
            
        if self.output == 'LSTM':
            if sim != 2:

                temp = (self.data['alt'].tolist()) # find where alt<135, idx following a gaussian distribution not linspace
                index = [i for i, v in enumerate(temp) if v <= 145*1e3]
                idx = np.round(np.linspace(index[0], index[-1], 100)).astype(int)
                time = (np.array(self.data['time'].tolist())[idx.astype(int)])
                time = [(t+self.time)/(60*60*24) for t in time]
                dragpassage_time = (np.array(self.data['time'].tolist())[idx.astype(int)])
                a = (np.array(self.data['a'].tolist())[idx.astype(int)])
                e = (np.array(self.data['e'].tolist())[idx.astype(int)])
                ra = np.array(np.multiply(a,(1+e)))
                hp = np.array(self.data['alt'].tolist())[idx.astype(int)]
                incl = (np.array(self.data['i'].tolist())[idx.astype(int)])
                OMEGA = (np.array(self.data['OMEGA'].tolist())[idx.astype(int)])
                omega = (np.array(self.data['omega'].tolist())[idx.astype(int)])
                vi = (np.array(self.data['vi'].tolist())[idx.astype(int)])
                heat_rate = (np.array(self.data['heat_rate'].tolist())[idx.astype(int)])
                rho = (np.array(self.data['rho'].tolist())[idx.astype(int)])

                # substitute max of the short list with the real maximum
                max_index = np.argmax(heat_rate)
                heat_rate[max_index] = max(self.data['heat_rate'].tolist())
                max_index = np.argmax(rho)
                rho[max_index] = max(self.data['rho'].tolist())

                self.time += float(self.data['time'].tolist()[-1])

                print('Ra0', round(ra[0] * 1e-3,3), 'Ra', round(ra[-1] * 1e-3,3), 'Raf', self.raf * 1e-3, 'hp0', round(self.hp_prev * 1e-3,4),
                      'hp', round(min(hp) * 1e-3,4) ,'OMEGA',round(math.degrees(OMEGA[-1]),3),'INC',round(math.degrees(incl[-1]),3),
                      'omega',round(math.degrees(omega[-1]),3),'action',self.action, 'maxrate',round(heat_rate[max_index],4))
                # print('len',len(np.array([time, ra,hp,vi,rho,heat_rate])))
                self.hp_prev = min(hp)
                return np.array([dragpassage_time, ra,hp,omega, OMEGA, incl, vi,rho,heat_rate])
                
            else:
                temp1 = self.data['a']
                temp2 = self.data['e']
                temp3 = self.data['vi']
                temp = np.divide((np.multiply(temp1, (1 - np.square(temp2)))),(1+np.multiply(temp2,np.cos(temp3)))) - self.p.Rp_e
                try:
                    index = [i for i, v in enumerate(temp) if v <= 145*1e3]
                    print('index', index[0], index[-1])
                except:
                    index = [0,len(temp)-1]
                    np.set_printoptions(threshold=sys.maxsize)
                    print('temp',temp)
                    # exit()

                idx = np.round(np.linspace(index[0], index[-1], 100)).astype(int)
                rho = (np.array(self.data['rho'])[idx.astype(int)])
                time = (np.array(self.data['time'])[idx.astype(int)])
                dragpassage_time = (np.array(self.data['time'])[idx.astype(int)])
                time = [(t+self.time)/(60*60*24) for t in time]
                a = (np.array(self.data['a'])[idx.astype(int)])
                e = (np.array(self.data['e'])[idx.astype(int)])
                ra = np.array(np.multiply(a,(1+e)))
                vi = (np.array(self.data['vi'])[idx.astype(int)])
                hp = np.divide((np.multiply(a, (1 - np.square(e)))),(1+np.multiply(e,np.cos(vi)))) - self.p.Rp_e
                incl = (np.array(self.data['incl'])[idx.astype(int)])
                OMEGA = (np.array(self.data['OMEGA'])[idx.astype(int)])
                omega = (np.array(self.data['omega'])[idx.astype(int)])
                heat_rate = (np.array(self.data['heat_rate'])[idx.astype(int)])
                self.time += float(self.data['time'][-1])

                # substitute max of the short list with the real maximum
                max_index = np.argmax(heat_rate)
                heat_rate[max_index] = max(self.data['heat_rate'])
                max_index = np.argmax(rho)
                rho[max_index] = max(self.data['rho'])

                return np.array([dragpassage_time, ra,hp,omega,OMEGA,incl, vi,rho,heat_rate])

    def _get_newstate(self, sim):
        """Returns new state.
        """

        if sim != 2:
            self.ra = float(self.data['a'].tolist()[-1]) * (1 + float(self.data['e'].tolist()[-1]))
            self.hp = float(self.data['a'].tolist()[-1]) * (1 - float(self.data['e'].tolist()[-1])) - self.p.Rp_e
            self.incl = math.degrees(float(self.data['i'].tolist()[-1]))
            self.OMEGA = math.degrees(float(self.data['OMEGA'].tolist()[-1]))
            self.omega = math.degrees(float(self.data['omega'].tolist()[-1]))
            self.year = int(self.data['year'].tolist()[-1])
            self.month = int(self.data['month'].tolist()[-1])
            self.day = int(self.data['day'].tolist()[-1])
            self.hour = int(self.data['hour'].tolist()[-1])
            self.mins = int(self.data['min'].tolist()[-1])
            self.sec = int(self.data['sec'].tolist()[-1])
        else:
            self.ra = (self.data['a'][-1]) * (1 + (self.data['e'][-1]))
            self.hp = (self.data['a'][-1]) * (1 - (self.data['e'][-1])) - self.p.Rp_e
            self.incl = self.data['incl'][-1]
            self.OMEGA = self.data['OMEGA'][-1]
            self.omega = self.data['omega'][-1]
        self.montecarlo_on = 1
        if sim != 0 and self.hp>1351*1e3:
            self._getstate_outerspace(ri =(self.p.Rp_e + 160 * 1e3), rf = self.ra)
        self.state = [self.ra, self.hp, self.incl, self.OMEGA, self.omega, self.year, self.month, self.day,
                             self.hour, self.mins,
                             self.sec, self.raf, self.montecarlo,self.montecarlo_on, self.pid, self.sim,self.directory]  # define initial conditions

        return

    def _set_action(self, action, i=0):
        """Applies the given action to the simulation.
        """
        # if a == 0:
        if self.args.algo == 'dqn':
            AVAIL_DELTAV = [-1, -0.5, -0.3, -0.2, -0.1,-0.05, 0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 1,]  # [-1,-0.1,-0.01, 0.0, 0.01, 0.1, 1]  #

            # if self.phase  == 'Walkout':
            #     AVAIL_DELTAV =[-1.5,-1, -0.5, -0.3, -0.2,-0.1, 0.0, 0.1, 0.2, 0.3, 0.5, 1, 1.5] # [-1,-0.1,-0.01, 0.0, 0.01, 0.1, 1]  #
            # elif self.phase == 'Endgame':
            #     AVAIL_DELTAV = [-1.5,-1, -0.5, -0.3, -0.2,-0.1, 0.0, 0.1, 0.2, 0.3, 0.5,1,1.5]  # [-1,-0.1,-0.01, 0.0, 0.01, 0.1, 1]  #

            AVAIL_PHI = [0, 180]

            deltav = abs(AVAIL_DELTAV[action])
            # print('action',AVAIL_DELTAV[action])
            if AVAIL_DELTAV[action]<=0:
                phi = AVAIL_PHI[0]
                # print("env LOWER MANEUVER!")
            else:
                phi = AVAIL_PHI[1]
                # print("env RAISE MANEUVER!")
        else:
            # print(action)
            velocity = min(max(action
                               , self.min_deltav), self.max_deltav)
            deltav = abs(velocity)
            if velocity >= 0:
                phi = 180
            elif velocity<0:
                phi = 0

        self.action = [deltav, phi]
        # print(self.action)
        if i == 0:
            self.deltav += deltav
            return self.deltav
        else:
            return deltav

    def _is_success(self,observations):
        """Indicates whether or not the achieved ct successfully achieved the desired goal.
        """
        if self.output == 'NN':
            ra = observations[1]
        elif self.output == 'LSTM':
            ra = observations[1][-1]
        if bool(abs(ra - self.finalstate)<=10*1e3):
            print('Success!',ra*1e-3,self.finalstate)

        return bool(abs(ra - self.finalstate)<=10*1e3)

    def _impact(self,observations):
        """Indicates whether or not impact.
        """
        hp = observations[2]
        if self.output == 'LSTM':
            hp = min(hp)

        return bool(hp<85*1e3)

    def _outdragpassage(self,observations):
        """Indicates whether or not the achieved goal successfully achieved the desired goal.
        """
        hp = observations[2]
        if self.output == 'LSTM':
            hp = min(hp)
        # print('hp',hp)
        return bool(hp > 145 * 1e3)

    def _noreached_goal(self,observations):
        if self.output == 'NN':
            ra = observations[1]
        elif self.output == 'LSTM':
            ra = observations[1][-1]
        return bool((ra - self.finalstate)<-10*1e3)

    def _under_thermal_loads(self,observations,sim):
        """Indicates whether or not achieved thermal violation.
        """
        if self.output == 'NN':
            ra = observations[1]
        elif self.output == 'LSTM':
            ra = observations[1][-1]

        if sim !=2:
            # load = np.array([float(max(self.data['heat_load'].tolist()))])
            if self.output == 'NN':
                maxrate = observations[-1]
            elif self.output == 'LSTM':
                maxrate = max(observations[-1])
            # maxrate = np.array([float(max(self.data['heat_rate'].tolist()))])
        else:
            # maxrate = np.array([max(self.data['heat_rate'])])
            if self.output == 'NN':
                maxrate = observations[-1]
            elif self.output == 'LSTM':
                maxrate = max(observations[-1])
            # load = np.array([self.data['heat_load']])
        if abs(ra-self.finalstate)>100*1e3:
            return bool(maxrate<= 0.25 and maxrate>= 0.05)
        else:
            return bool(maxrate<= 0.25)




    def _getstate_outerspace(self,ri,rf):
        a = (self.ra + self.hp + self.p.Rp_e) / 2
        e = (self.ra - (self.hp + self.p.Rp_e)) / (self.ra + (self.hp + self.p.Rp_e))
        i = math.radians(self.incl)
        OMEGA = math.radians(self.OMEGA)
        omega = math.radians(self.omega)
        vi = self.vi

        # evaluate new time
        OMEGA_rate = (- (1.5 * ((self.p.mu) ** 0.5 * self.p.J2 * self.p.Rp_e ** 2) / ((1 - e ** 2) * a ** (7.0 / 2.0))) * math.cos(
            i))  # rad/s
        omega_rate = (OMEGA_rate * (((5.0 / 2.0) * (math.sin(i) ** 2) - 2.0) / math.cos(i)))  # rad/s

        # define initial and final state
        if ri == self.ra:
            initial_state_angle = math.pi
        else:
            initial_state_angle = math.acos(1 / e * (a * (1 - e ** 2) / ri - 1))#math.pi
        if rf == self.ra:
            final_state_angle = math.pi
        else:
            final_state_angle = - math.acos(1 / e * (a * (1 - e ** 2) / rf - 1))

        E_initialstate = 2.0 * math.atan(
            ((1 - e) / (1 + e)) ** 0.5 * math.tan((initial_state_angle) * 0.5))  # eccentric anomaly
        E_finalstate = 2.0 * math.atan(
            ((1 - e) / (1 + e)) ** 0.5 * math.tan((final_state_angle) * 0.5))  # eccentric anomaly

        # evaluate time to reach next state
        delta_t = abs((a ** 3.0 / self.p.mu) ** 0.5 * (
                (E_finalstate - e * math.sin(E_finalstate)) - (E_initialstate - e * math.sin(E_initialstate))))
        # print(delta_t)
        self.time += delta_t
        self.omega = math.degrees(omega + omega_rate * delta_t)
        self.OMEGA = math.degrees(OMEGA + OMEGA_rate * delta_t)
        self.vi = final_state_angle
        date = datetime(year=int(self.year), month=int(self.month), day=int(self.day), hour=int(self.hour), minute=int(self.mins),
                        second=int(self.sec)) + timedelta(seconds=int(delta_t))
        self.year, self.month, self.day, self.hour, self.mins, self.sec = date.year, date.month, date.day, date.hour, date.minute, date.second

        self.state = [self.ra, self.hp, self.incl, self.OMEGA, self.omega, self.vi, self.year, self.month, self.day,
                      self.hour, self.mins,
                      self.sec, self.raf, self.montecarlo, self.montecarlo_on, self.pid, self.sim, self.directory]  # define initial conditions


    def _getstate_dragpassage(self):
        a = (self.ra + (self.hp + self.p.Rp_e)) / 2
        e = (self.ra - (self.hp + self.p.Rp_e)) / (self.ra + (self.hp + self.p.Rp_e))
        i = self.incl
        OMEGA = math.radians(self.OMEGA)
        omega = math.radians(self.omega)
        vi = self.vi
        OE = [a,e,i,OMEGA,omega,vi]
        [pos_ii_org, vel_ii_org] = orbitalelemtorv(OE,self.p)
        pos_ii = pos_ii_org
        vel_ii = vel_ii_org
        r0 = np.linalg.norm(pos_ii)  # Inertial position magnitude
        v0 = np.linalg.norm(vel_ii)  # Inertial velocity magnitude
        h0 = r0 - self.p.Rp_e

        h_ii = np.cross(pos_ii, vel_ii)
        arg = np.median([-1, 1, np.linalg.norm(h_ii) / (r0 * v0)])  # limit to[-1, 1]
        gamma0 = math.acos(arg)
        if np.inner(pos_ii, vel_ii) < 0:
            gamma0 = -gamma0

        initial_state_angle = OE[5]
        e = OE[1]
        a = OE[0]
        final_state_angle = -initial_state_angle
        E_initialstate = 2 * math.atan(
            ((1 - e) / (1 + e)) ** 0.5 * math.tan((initial_state_angle) / 2))  # eccentric anomaly
        E_finalstate = 2 * math.atan(
            ((1 - e) / (1 + e)) ** 0.5 * math.tan((final_state_angle) / 2))  # eccentric anomaly

        # evaluate time to reach next state
        delta_t = (a ** 3 / self.p.mu) ** 0.5 * (
                    (E_finalstate - e * math.sin(E_finalstate)) - (E_initialstate - e * math.sin(E_initialstate)))
        t_p = delta_t / 2

        step_time = delta_t/0.1
        t_cf = np.linspace(0, delta_t, int(step_time))

        cost_3 = v0 * gamma0

        h_cf = h0 + cost_3 * (t_cf - (np.square(t_cf) / (2 * t_p)))
        rho = dm(h=h_cf, p=self.p)[0]

        T = self.p.T
        RT = T * self.p.R
        S = (v0 / (2 * RT) ** 0.5)
        Area_SA = 7.26
        Area_SC = 2.2*1.7
        Area_tot = Area_SA+Area_SC
        CL90, CD90 = am(self, Area_SA, Area_SC, S)
        CL0, CD0 = am(self, Area_SA, Area_SC, S)
        mass = 461
        Rp = self.p.Rp_e

        aoa_profile = [math.pi/2] * len(t_cf)

        CD_t = np.add(CD0, np.divide(np.multiply(aoa_profile, (CD90 - CD0)), math.pi / 2))
        CL_t = np.add(CL0, np.divide(np.multiply(aoa_profile, (CL90 - CL0)), math.pi / 2))
        cost_1 = np.divide(rho * CD_t * Area_tot, 2 * mass)
        cost_2 = np.divide(rho * CL_t * Area_tot, 2 * mass)

        a0 = 0.0016  # 0.0034
        c0 = 5e-06  # 2.9648e-06
        mean_a = 3.38  # 3.5839
        mean_c = 2.6  # 2.5413
        mean_b = -8.25  # -14.5169
        mean_d = -0.001  # 2.5413

        f1 = -0.005 * v0 + 27.87
        f2 = (a0 * (mean_a ** (2 * abs(math.degrees(gamma0) + 3)) * math.exp(mean_b * (v0 / 1000 - 3.7))) + c0 * (
                mean_c ** (2 * abs(math.degrees(gamma0) + 3)) * math.exp(mean_d * (v0 / 1000 - 3.7)))) * (t_cf) / (
                         2 * t_p)

        f2_solar_panels = f2 * math.pi/2 * Area_SA / Area_tot
        f2_spacecraft = f2 * math.pi / 2 * Area_SC / Area_tot

        epsilon = f1 + f2_solar_panels + f2_spacecraft

        k1 = (cost_2 + np.divide(1, (Rp + h_cf)))
        k2 = np.multiply(cost_1 * cost_3, (1 - t_cf / t_p))
        k3 = -self.p.g_ref - epsilon
        cost = v0 - (k2[0] / k1[0] - ((k2[0] / k1[0]) ** 2 - 4 * (k3[0] / k1[0])) ** 0.5) / 2
        v_cf = (np.divide(k2, k1) - np.sqrt(
            np.square(np.divide(k2, k1)) - 4 * np.divide(k3, k1))) / 2 + cost

        gamma_cf = cost_3 * np.divide((1 - t_cf / t_p), v_cf)
        heat_rate = ht(self,v_cf,T,rho)
        heat_load = sum(heat_rate) * t_cf[-1] / len(t_cf)
        ang_mom = np.multiply(np.multiply((self.p.Rp_e+h_cf),v_cf),np.cos(gamma_cf))
        energy = np.square(v_cf)/2 - self.p.mu/(self.p.Rp_e+h_cf)


        a = -self.p.mu/(2*energy)

        semilatus = np.square(ang_mom)/self.p.mu
        e = np.sqrt(1-np.divide(semilatus,a))

        # evaluate new time
        OMEGA_rate = (- (1.5 * ((self.p.mu) ** 0.5 * self.p.J2 * self.p.Rp_e ** 2) / ((1 - e ** 2) * a ** (7.0 / 2.0))) * math.cos(
            i))  # rad/s
        omega_rate = (OMEGA_rate * (((5.0 / 2.0) * (math.sin(i) ** 2) - 2.0) / math.cos(i)))  # rad/s
        omega = np.degrees(omega + omega_rate * t_cf)
        OMEGA = np.degrees(OMEGA + OMEGA_rate * t_cf)

        incl = [self.incl]* len(t_cf)
        vi = np.arccos(np.multiply(np.divide(1, e), (np.divide(np.multiply(a, (1 - np.square(e))), (self.p.Rp_e+h_cf)) - 1)))

        vi = [-vi[i] if gamma_cf[i] >0 else vi[i] for i in range(0,len(vi))]
        t_cf = np.array([item + self.time for item in t_cf])

        date = datetime(year=int(self.year), month=int(self.month), day=int(self.day), hour=int(self.hour), minute=int(self.mins),
                        second=int(self.sec)) + timedelta(seconds=int(delta_t))
        self.year, self.month, self.day, self.hour, self.mins, self.sec = date.year, date.month, date.day, date.hour, date.minute, date.second
        data = {}
        data['a'] = a
        data['e'] = e
        data['incl'] = incl  # degrees
        data['omega'] = omega # degrees
        data['OMEGA'] = OMEGA # degrees
        data['vi'] = vi
        data['time'] = t_cf
        data['heat_rate'] = heat_rate
        data['heat_load'] = heat_load
        data['rho'] = rho
        return data



    def _env_setup(self,sim):
        """Initial configuration of the environment. Can be used to configure initial state
        and extract information from the simulation.
        """
        temp =0
        while temp ==0:
            retval = os.getcwd()
            os.chdir(directory)

            if sim ==0:
                os.chdir(directory)
                # os.chdir('/Users/Josephine/Dropbox/Aerobraking/ABTS')
                # if all orbit (apoapsis to apoapsis)
                call_ABTS(self)

            elif sim == 1:
                os.chdir(directory)
                # os.chdir('/Users/Josephine/Dropbox/Aerobraking/ABTS')
                # if only drag passage (160 to 160)
                self.get_newperiapsis()
                self._getstate_outerspace(ri =self.ra, rf =(self.p.Rp_e + 160 * 1e3))
                call_ABTS(self)

            elif sim == 2:
                # if only drag passage (160 to 160)
                self.get_newperiapsis()
                self._getstate_outerspace(ri =self.ra, rf =(self.p.Rp_e + 160 * 1e3))
                data = self._getstate_dragpassage()
                temp = 1


            if sim != 2:
                os.chdir(retval)
                filename = 'Sim'+str(self.pid)

                ii = 0

                try:
                    # print(self.directory + filename + '.csv')
                    data = pd.read_csv(self.directory + filename + '.csv')
                    temp =1
                except:
                    ii += 1
                    time.sleep(0.1)
                    # print(self.directory + filename + '.csv')
                    if ii > 100:
                        break
                data.head()
                os.remove(self.directory + filename + '.csv')
        self.data = data
        return

    def get_newperiapsis(self):
        a = (self.ra + (self.hp + self.p.Rp_e)) / 2
        e = (self.ra - (self.hp + self.p.Rp_e)) / (self.ra + (self.hp + self.p.Rp_e))
        i = self.incl
        OMEGA = self.OMEGA
        omega = self.omega
        vi = math.pi
        OE = [a,e,i,OMEGA,omega,vi]

        rp = a*(1-e) # delete
        v_prova = (self.p.mu*(1+e)/rp)**0.5  # delete
        # evaluate new periapsis
        [_,v] = orbitalelemtorv(OE,self.p)

        v = np.linalg.norm(v)
        if round(self.action[1]) == 0:
            v -= self.action[0]
            v_prova -= self.action[0]  # delete
        elif round(self.action[1]) == 180:
            v += self.action[0]
            v_prova += self.action[0]  # delete

        Energy = (v ** 2) * 0.5 - self.p.mu / (self.ra)
        a = - self.p.mu / (2 * Energy)
        # self.hp_prev = self.hp
        self.hp = 2*a-(self.ra) - self.p.Rp_e
        hp_prova = self.p.mu*(1+e)/v_prova**2  - self.p.Rp_e
        print('v',v,'v_prova',v_prova,'self.hp',self.hp,hp_prova)
#        print('Prev hp',prev_hp*1e-3,'New hp',self.hp*1e-3,'deltav',self.action[0],'phi',self.action[1])
        return

    def close(self):
        pass
