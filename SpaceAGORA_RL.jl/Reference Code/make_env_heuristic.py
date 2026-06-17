"""Aerobraking task"""
import numpy as np
from subprocess import Popen, PIPE
# from gym import core, spaces
import sys
from scipy import optimize
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


class AerobrakingEnvHeuristic():
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

    def __init__(self, args):
        self.max_ra = 30000 * 1e3  # m
        self.min_ra = 3396 * 1e3  # m
        self.max_hp = 140 * 1e3  # m
        self.min_hp = 80 * 1e3  # m
        self.max_hr = np.inf
        self.min_hr = 0.0
        self.max_hl = np.full((100), np.inf)
        self.min_hl = np.zeros(100)
        self.min_density = np.zeros(100)
        self.min_velocity = np.zeros(100)
        self.min_position = np.zeros(100)
        self.max_density = np.full((100), np.inf)
        self.max_velocity = np.full((100), np.inf)
        self.max_position = np.full((100), np.inf)
        self.max_angle = np.full((100), np.pi * 2)
        self.min_angle = np.zeros(100)
        self.min_heatrate = np.zeros(100)
        self.max_heatrate = np.full((100), np.inf)
        self.min_sma = np.zeros(100)
        self.max_sma = np.full((100), np.inf)
        self.min_ecc = np.zeros(100)
        self.max_ecc = np.full((100), np.inf)
        self.max_time = np.inf
        self.min_time = 0.0
        self.max_deltav = 1.5
        self.min_deltav = -1.5
        self.min_target = 3500 * 1e3
        self.max_target = 7000 * 1e3
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
        self.directory = directory2 + args.data_dir
        self.nominal = args.nominal
        self.delta_t = 0

        #        low = np.concatenate([[self.min_target],[self.min_hr],[self.min_hr],[self.min_ra],[self.min_hp]])
        # high=np.concatenate([[self.max_target], self.max_hl.flat, self.max_sma.flat, self.max_ecc.flat, self.max_density.flat, self.max_velocity.flat, self.max_position.flat, self.max_heatrate.flat],dtype=object)
        # low = np.concatenate([[self.min_target], self.min_hl.flat,  self.min_sma.flat, self.min_ecc.flat, self.min_density.flat, self.min_velocity.flat,self.min_position.flat, self.min_heatrate.flat],dtype=object)

        if self.output == 'NN':
            high = np.array([2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0], dtype=np.float32, )

            low = np.array([0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float32, )
            self.observation_space = spaces.Box(low=low, high=high, dtype=np.float32)  # np.float32)
        else:
            high = np.array([self.max_position, self.max_position, self.max_velocity, self.max_angle, self.max_angle,
                             self.max_angle, self.max_angle, self.max_density, self.max_heatrate])
            low = np.array([self.min_position, self.max_position, self.min_velocity, self.min_angle, self.min_angle,
                            self.min_angle, self.min_angle, self.min_density, self.min_heatrate])
            self.observation_space = spaces.Box(low=low, high=high)
        # self.action_space = spaces.Box(low=self.min_deltav,high=self.max_deltav,shape=(21,),dtype='float32')

    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        print('seed', seed)
        return [seed]

    def get_initialstate(self,pid,sim):
        if self.phase == 'Endgame':
            self.ra = 6000* 1e3+self.np_random.uniform(low=(-2.5)* 1e3, high=(2.5)*1e3)#self.ra = self.np_random.uniform(low=(6023.95)* 1e3, high=(5000)*1e3)
        elif self.phase == 'Walkout':
            self.ra = 4906*1e3+self.np_random.uniform(low=(-2.5)* 1e3, high=(2.5)*1e3)
        elif self.phase == 'Main':
            self.ra = 10038*1e3+self.np_random.uniform(low=(-2.5)* 1e3, high=(2.5)*1e3)
        elif self.phase == 'Campaign':
            self.ra = 28533*1e3+self.np_random.uniform(low=(-2.5)* 1e3, high=(2.5)*1e3)
        self.rai = self.ra
        self.raf = self.finalstate#self.np_random.uniform(low=4806 * 1e3, high=5000 * 1e3)
        self.vi = math.pi
        if self.nominal == 0:
            self.incl = self.np_random.uniform(low=88.6, high=98.6)
            if not self.phase == 'Campaign':
                self.hp = 92 * 1e3 + self.np_random.uniform(low=(-1.5) * 1e3, high=(1.5) * 1e3)  # self.np_random.uniform(low=108 * 1e3, high=120 * 1e3)
                self.hp_prev = self.hp
                self.OMEGA = self.np_random.uniform(low=110, high=120)
                self.omega = self.np_random.uniform(low=60, high=89)
                self.year = 2001#self.np_random.randint(low=2001, high=2021)
                self.month = 12#self.np_random.randint(low=1, high=12)
            else:
                self.hp = 92 * 1e3 + self.np_random.uniform(low=(-1.5) * 1e3, high=(1.5) * 1e3)  # self.np_random.uniform(low=108 * 1e3, high=120 * 1e3)
                self.hp_prev = self.hp
                self.OMEGA = self.np_random.uniform(low=110, high=120)
                self.omega = self.np_random.uniform(low=91, high=130)
                self.year = 2001#self.np_random.randint(low=2001, high=2021)
                self.month = self.np_random.randint(low=10, high=12)
            self.day = self.np_random.randint(low=1, high=28)
            self.hour = self.np_random.randint(low=0, high=23)
            self.mins = self.np_random.randint(low=0, high=59)
            self.sec = self.np_random.randint(low=0, high=59)
            self.montecarlo_real = self.np_random.randint(low=0, high=1000)
            self.montecarlo_predict = self.np_random.randint(low=0, high=1000)
            self.montecarlo = self.montecarlo_real
            self.montecarlo_on = 1
        else:
            self.hp = 86 * 1e3 + self.np_random.uniform(low=(-1.) * 1e3, high=(
                                                                                  1.) * 1e3)  # self.np_random.uniform(low=108 * 1e3, high=120 * 1e3)
            self.hp_prev = self.hp
            self.incl = 93.6 + self.np_random.uniform(low=-0.25, high=0.25)
            self.OMEGA = 115 + self.np_random.uniform(low=-0.25, high=0.25)
            self.omega = 66 + self.np_random.uniform(low=-0.25, high=0.25)
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
                self.omega = 89 + self.np_random.uniform(low=-0.25, high=0.25)
                self.hp = 92 * 1e3 + self.np_random.uniform(low=(-1.) * 1e3, high=(1.) * 1e3)
            elif self.phase == 'Campaign':
                self.year = 2001
                self.month = 11
                self.day = 6
                self.OMEGA = 114 + self.np_random.uniform(low=-0.25, high=0.25)
                self.omega = 109 + self.np_random.uniform(low=-0.25, high=0.25)
                self.hp = 84.6 * 1e3 + self.np_random.uniform(low=(-1.) * 1e3, high=(
                                                                                        1.) * 1e3)  # self.np_random.uniform(low=108 * 1e3, high=120 * 1e3)
                self.hp_prev = self.hp
            self.hour = 0
            self.mins = 0
            self.sec = 0
            self.montecarlo_real = self.np_random.randint(low=0, high=1000)
            self.montecarlo_predict = self.np_random.randint(low=0, high=1000)
            self.montecarlo = self.montecarlo_real
            self.montecarlo_on = 1
        self.pid = pid
        self.time = 0.0
        self.deltav = 0.0
        print('initial hp', self.hp / 1e3)
        if not sim:
            self.sim = 'Orbits'
        else:
            self.sim = 'Drag Passage'
        self.initialstate = [self.ra, self.hp, self.incl, self.OMEGA, self.omega, self.year, self.month, self.day,
                             self.hour, self.mins,
                             self.sec, self.raf, self.montecarlo, self.montecarlo_on, self.pid, self.sim,
                             self.directory]  # define initial conditions
        self.state = self.initialstate
        return

    def reset(self,sim,pid):
        # print('Start New Episode')
        self.count = 0
        self.get_initialstate(pid, sim)
        if self.output == 'NN':
            ra_i =  (self.ra - 3396200) / (10100000 - 3396200)
            obs = [0.0, ra_i, 1.0,0.0,0.0,0.0,0.0,0.0, 0.0]
        else:
            obs = np.array([np.full((100),0),np.full((100),1),np.full((100),1),np.full((100),0),np.full((100),0),np.full((100),0),np.full((100),0),np.full((100),0),np.full((100),0)])
        return obs#(observations, reward, terminal, {})#self._get_ob()

    def step(self, sim):
        self.count += 1
        self._predict(sim)
        self._env_setup(sim)
        self._get_newstate(sim)
        observations = self._get_ob(sim)
        terminal = self._terminal(observations,sim)
        reward = self._reward(observations,sim, self.action)
        if self.action[1] == 0:
            action = -self.action[0]
        elif self.action[1] == 180:
            action = self.action[0]
        return (observations, reward, action, terminal, {})


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
                    print(self.directory + filename + '.csv')
                    data = pd.read_csv(self.directory + filename + '.csv')
                    temp =1
                except:
                    ii += 1
                    time.sleep(0.1)
                    print(self.directory + filename + '.csv')
                    if ii > 100:
                        break
                data.head()
                os.remove(self.directory + filename + '.csv')
        self.data = data
        return


    def _get_ob(self, sim):
        """Returns the observation.
        """
        if self.output == 'NN':
            if sim != 2:
                load = np.array([float(max(self.data['heat_load'].tolist()))])
                heat_rate = np.array([float(max(self.data['heat_rate'].tolist()))])
                ra = np.array([float(self.data['a'].tolist()[-1]) * (1 + float(self.data['e'].tolist()[-1]))])
                ra0 = np.array([float(self.data['a'].tolist()[0]) * (1 + float(self.data['e'].tolist()[0]))])
                hp = np.array([float(min(self.data['alt'].tolist()))])
                temp = (self.data['rho'].tolist())
                self.time += float(self.data['time'].tolist()[-1])
                dragpassage_time = np.array([float(self.data['time'].tolist()[-1])])
                idx = np.round(np.linspace(0, len(temp) - 1, 100)).astype(int)

                rho = (np.array([float(max(self.data['rho'].tolist()))]))
                omega = (np.array([float(self.data['omega'].tolist()[-1])]))
                OMEGA = (np.array([float(self.data['OMEGA'].tolist()[-1])]))
                inc = (np.array([float(self.data['i'].tolist()[-1])]))
                year = int(self.data['year'].tolist()[-1])
                month = int(self.data['month'].tolist()[-1])
                day = int(self.data['day'].tolist()[-1])
                time = (np.array([float(date(year, month, day).toordinal())]))
                sma = (np.array(self.data['a'].tolist()))
                ecc = (np.array(self.data['e'].tolist()))
                vel = (np.array(self.data['vel_ii_mag'].tolist())[idx.astype(int)])
                pos = (np.array(self.data['pos_ii_mag'].tolist()))
                # heat_rate = (np.array(self.data['heat_rate'].tolist())[idx.astype(int)])
                time_man = np.array([self.time / (60 * 60 * 24)])
                deltav = np.array([self.deltav])
                target = np.array([self.raf])
                if sim == 0:
                    h = pos - self.p.Rp_e
                    # print('h',h)
                    # h = (self.data['alt'].tolist())
                    vi = (self.data['vi'].tolist())
                    # print('vi',vi)
                    # print('t',self.data['time'].tolist())
                    idx = (np.abs(np.asarray(vi) - 0)).argmin()
                    # print('idx',idx)
                    idx_h1 = (np.abs(np.asarray(h[0:idx]) - 160 * 10 ** 3)).argmin()
                    idx_h2 = (np.abs(np.asarray(h[idx:]) - 160 * 10 ** 3)).argmin()
                    # print(idx_h1,idx_h2)
                    t = abs(float(self.data['time'].tolist()[idx + idx_h2]) - float(self.data['time'].tolist()[idx_h1]))
                    dragpassage_time = np.array([t])
                    self.delta_t = float(self.data['time'].tolist()[-1])
                    # print(t,dragpassage_time)


            else:
                target = np.array([self.raf])
                heat_rate = np.array([max(self.data['heat_rate'])])
                load = np.array([self.data['heat_load']])
                ra = np.multiply(self.data['a'], (1 + self.data['e']))
                ra = np.array([ra[-1]])
                ra0 = np.array([ra[0]])
                omega = np.array([self.data['omega'][-1]])
                OMEGA = np.array([self.data['OMEGA'][-1]])
                inc = np.array([self.data['incl'][-1]])
                rho = np.array([max(self.data['rho'])])
                hp = np.multiply(self.data['a'], (1 - self.data['e'])) - self.p.Rp_e
                hp = np.array([hp[-1]])
                self.time += float(self.data['time'][-1])
                dragpassage_time = np.array([float(self.data['time'][-1])])
                time_man = np.array([self.time / (60 * 60 * 24)])
            print('Passage',self.count,'Ra0', round(ra0[0]*1e-3,3),'Ra', round(ra[-1]*1e-3,3),'Raf', round(self.raf*1e-3,3),
                  'hp0',round(self.hp_prev*1e-3,3),'hp', round(hp[-1]*1e-3,3),'action',self.action,'maxrate',round(heat_rate[0],3),
                  'dragpassagetime', round(dragpassage_time[0],2), 'omega', round(math.degrees(omega[0]),2),
                  'OMEGA',round(math.degrees(OMEGA[0]),2))

            # print('Ra0', ra0[0] * 1e-3, 'Ra', ra[-1] * 1e-3, 'Raf', self.raf * 1e-3, 'hp0', self.hp_prev * 1e-3, 'hp',
            #       hp[-1] * 1e-3, 'action', self.action, 'maxrate', heat_rate, 'dragpassagetime', dragpassage_time)
            self.hp_prev = hp[-1]
            return np.concatenate((dragpassage_time, ra, hp, omega, OMEGA, inc, time, rho, heat_rate))
            # (time_man,ra,hp,target,maxrate))

        if self.output == 'LSTM':
            if sim != 2:

                temp = (self.data['alt'].tolist())  # find where alt<135, idx following a gaussian distribution not linspace
                index = [i for i, v in enumerate(temp) if v <= 145 * 1e3]
                idx = np.round(np.linspace(index[0], index[-1], 100)).astype(int)
                time = (np.array(self.data['time'].tolist())[idx.astype(int)])
                time = [(t + self.time) / (60 * 60 * 24) for t in time]
                dragpassage_time = (np.array(self.data['time'].tolist())[idx.astype(int)])
                a = (np.array(self.data['a'].tolist())[idx.astype(int)])
                e = (np.array(self.data['e'].tolist())[idx.astype(int)])
                ra = np.array(np.multiply(a, (1 + e)))
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

                print('Ra0', round(ra[0] * 1e-3, 3), 'Ra', round(ra[-1] * 1e-3, 3), 'Raf', self.raf * 1e-3, 'hp0',
                      round(self.hp_prev * 1e-3, 4),
                      'hp', round(min(hp) * 1e-3, 4), 'OMEGA', round(math.degrees(OMEGA[-1]), 3), 'INC',
                      round(math.degrees(incl[-1]), 3),
                      'omega', round(math.degrees(omega[-1]), 3), 'action', self.action, 'maxrate',
                      round(heat_rate[max_index], 4))
                # print('len',len(np.array([time, ra,hp,vi,rho,heat_rate])))
                self.hp_prev = min(hp)
                return np.array([dragpassage_time, ra, hp, omega, OMEGA, incl, vi, rho, heat_rate])

            else:
                temp1 = self.data['a']
                temp2 = self.data['e']
                temp3 = self.data['vi']
                temp = np.divide((np.multiply(temp1, (1 - np.square(temp2)))),
                                 (1 + np.multiply(temp2, np.cos(temp3)))) - self.p.Rp_e
                try:
                    index = [i for i, v in enumerate(temp) if v <= 145 * 1e3]
                    print('index', index[0], index[-1])
                except:
                    index = [0, len(temp) - 1]
                    np.set_printoptions(threshold=sys.maxsize)
                    print('temp', temp)
                    # exit()

                idx = np.round(np.linspace(index[0], index[-1], 100)).astype(int)
                rho = (np.array(self.data['rho'])[idx.astype(int)])
                time = (np.array(self.data['time'])[idx.astype(int)])
                dragpassage_time = (np.array(self.data['time'])[idx.astype(int)])
                time = [(t + self.time) / (60 * 60 * 24) for t in time]
                a = (np.array(self.data['a'])[idx.astype(int)])
                e = (np.array(self.data['e'])[idx.astype(int)])
                ra = np.array(np.multiply(a, (1 + e)))
                vi = (np.array(self.data['vi'])[idx.astype(int)])
                hp = np.divide((np.multiply(a, (1 - np.square(e)))), (1 + np.multiply(e, np.cos(vi)))) - self.p.Rp_e
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

                return np.array([dragpassage_time, ra, hp, omega, OMEGA, incl, vi, rho, heat_rate])


    def _predict(self,sim):
        """Returns heat rate prediction"""
        self.montecarlo = self.montecarlo_predict
        self.action = [0.0,0]
        self._env_setup(sim)
        obs = self._get_ob(sim)
        # heat_rate = obs[-1]
        if sim != 2:
            maxrate = obs[-1]
        else:
            maxrate = max(obs[-1])


        if maxrate>=0.05 and maxrate<=0.25:
            self.action = [0.0, 180]
        elif maxrate<0.05:
            def func(action):
                self._reinitialize_state()
                self.action = [action, 0]
                self._env_setup(sim)
                obs = self._get_ob(sim)
                maxrate = obs[-1]
                return 10 * (maxrate - 0.15)

            try:
                action = optimize.bisect(func, 0.0, 1., xtol=5e-2)
            except:
                try:
                    action = optimize.bisect(func, 0.0, 6, xtol=5e-2)
                except:
                    action = optimize.bisect(func, 0.0, 3, xtol=5e-2)
            self.action = [action, 0]
            print(self.action)
        elif maxrate>0.25:
            def func(action):
                self._reinitialize_state()
                self.action = [action, 180]
                self._env_setup(sim)
                obs = self._get_ob(sim)
                maxrate = obs[-1]
                return 10 * (maxrate - 0.15)
            try:
                action = optimize.bisect(func, 1., 0., xtol=5e-2)
            except:
                try:
                    action = optimize.bisect(func, 6, 0., xtol=5e-2)
                except:
                    action = optimize.bisect(func, 3, 0., xtol=5e-2)
            self.action = [action, 180]
        self._reinitialize_state()
        self.montecarlo = self.montecarlo_real
        return

    def _reinitialize_state(self):
        self.ra = self.state[0]
        self.hp = self.state[1]
        self.incl = self.state[2]
        self.OMEGA = self.state[3]
        self.omega = self.state[4]
        self.year = self.state[5]
        self.month = self.state[6]
        self.day = self.state[7]
        self.hour = self.state[8]
        self.mins = self.state[9]
        self.sec = self.state[10]
        self.time -= self.delta_t


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
        if sim != 0 and self.hp > 1351 * 1e3:
            self._getstate_outerspace(ri=(self.p.Rp_e + 160 * 1e3), rf=self.ra)
        self.state = [self.ra, self.hp, self.incl, self.OMEGA, self.omega, self.year, self.month, self.day,
                      self.hour, self.mins,
                      self.sec, self.raf, self.montecarlo, self.montecarlo_on, self.pid, self.sim,
                      self.directory]  # define initial conditions

        return


    def get_newperiapsis(self):
        a = (self.ra + (self.hp + self.p.Rp_e)) / 2
        e = (self.ra - (self.hp + self.p.Rp_e)) / (self.ra + (self.hp + self.p.Rp_e))
        i = self.incl
        OMEGA = self.OMEGA
        omega = self.omega
        vi = math.pi
        OE = [a,e,i,OMEGA,omega,vi]

        # evaluate new periapsis
        [_,v] = orbitalelemtorv(OE,self.p)

        v = np.linalg.norm(v)
        if round(self.action[1]) == 0:
            v -= self.action[0]
            # print('here')
        elif round(self.action[1]) == 180:
            v += self.action[0]
        Energy = (v ** 2) * 0.5 - self.p.mu / (self.ra)
        a = - self.p.mu / (2 * Energy)
        # self.hp_prev = self.hp
        self.hp = 2*a-(self.ra) - self.p.Rp_e
        #        print('Prev hp',prev_hp*1e-3,'New hp',self.hp*1e-3,'deltav',self.action[0],'phi',self.action[1])
        return

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
        self.delta_t = abs((a ** 3.0 / self.p.mu) ** 0.5 * (
                (E_finalstate - e * math.sin(E_finalstate)) - (E_initialstate - e * math.sin(E_initialstate))))
        self.time += self.delta_t
        self.omega = math.degrees(omega + omega_rate * self.delta_t)
        self.OMEGA = math.degrees(OMEGA + OMEGA_rate * self.delta_t)
        self.vi = final_state_angle
        date = datetime(year=int(self.year), month=int(self.month), day=int(self.day), hour=int(self.hour), minute=int(self.mins),
                        second=int(self.sec)) + timedelta(seconds=int(self.delta_t))
        self.year, self.month, self.day, self.hour, self.mins, self.sec = date.year, date.month, date.day, date.hour, date.minute, date.second

        self.state_outerspace = [self.ra, self.hp, self.incl, self.OMEGA, self.omega, self.vi, self.year, self.month, self.day,
                      self.hour, self.mins,
                      self.sec, self.raf, self.montecarlo, self.montecarlo_on, self.pid, self.sim, self.directory]  # define initial conditions

    def _terminal(self,observations,sim):
        if self.output == 'NN':
            ra = observations[1]
        elif self.output == 'LSTM':
            ra = observations[1][-1]
        temp = bool(ra<=(self.finalstate+20*10**3) or self._is_success(observations) or self._impact(observations) or self._outdragpassage(observations))
        return temp

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
                hr_cap = 0.5 # as it is in the middle of the corridor
                deltav_cap = abs(a[0]-self.max_deltav)

            if d_cap >= 0:
                r += (1 - d_cap ** (0.4)) * (0.025 / (0.1 * (2 * np.pi) ** 0.5)) * np.exp(
                    -0.5 * ((hr_cap - 0.5) / 0.1) ** 2) + (-0.1 - deltav_cap * 0.1)

#                r += 0.095*(1-d_cap**(0.4))* (0.08-1.77* hr_cap + 57.9 * hr_cap**2 -615.5 * hr_cap**3+3053.6 * hr_cap**4 - 7460.4 * hr_cap**5 + 9429.1 * hr_cap**6
#                       - 5950.6 * hr_cap**7 + 1487.6 * hr_cap**8)+ (0-deltav_cap*0.1)-0.1
#                r += 0.3*(1-d_cap**(0.4))* (0.08-1.77* hr_cap + 57.9 * hr_cap**2 -615.5 * hr_cap**3+3053.6 * hr_cap**4 - 7460.4 * hr_cap**5 + 9429.1 * hr_cap**6
#                       - 5950.6 * hr_cap**7 + 1487.6 * hr_cap**8)* (1.5-deltav_cap)
                # r +=  0.35*(1-d_cap**(0.4))*(0.01+1.24* hr_cap - 15.3 * hr_cap**2 +56.5 * hr_cap**3-64.3 * hr_cap**4 + 21.8 * hr_cap**5)* (1.5-deltav_cap)
                #0.5*(1-d_cap**(0.4))*(0.01+1.24* hr_cap - 15.3 * hr_cap**2 +56.5 * hr_cap**3-64.3 * hr_cap**4 + 21.8 * hr_cap**5)* (1-deltav_cap)
                #r += 0.5*(1-d_cap**(0.4))* (2.2* hr_cap - 5.17 * hr_cap**2 + 8.04 * hr_cap**3 - 4.05*hr_cap**4)* (1-deltav_cap)
                # r += (1-d_cap**(0.4))* (-0.005821928 + 2.190743* hr_cap - 5.164967 * hr_cap**2 + 8.045404 * hr_cap**3 - 4.126949*hr_cap**4)* (1-deltav_cap)
            else:
                r += 0.0
            print('d_cap',d_cap,'hr_cap',hr_cap,'deltav_cap',deltav_cap,'r',r)

        return r


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

        return bool(hp<80*1e3)

    def _outdragpassage(self,observations):
        """Indicates whether or not the achieved goal successfully achieved the desired goal.
        """
        hp = observations[2]
        if self.output == 'LSTM':
            hp = min(hp)
        # print('hp',hp)
        return bool(hp > 135 * 1e3)
        
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

    def _time_mission(self):
        return self.time/(60*60*24)

