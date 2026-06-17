from misc.restart import restart

import torch
import Algo
import make_env
import numpy as np
import argparse, os
import time
import sys
import math

import subprocess
# import multiprocessing as mp
parser = argparse.ArgumentParser()
parser.add_argument('--niter', type=int, default=2000, help='# training iterations')  # 100000
parser.add_argument('--ntest', type=int, default=100, help='# n of testing episode')  # 100
parser.add_argument('--test', type=int, default=0, help='# test mode')
parser.add_argument('--env', type=str, default='Aerobraking')
parser.add_argument('--algo', type=str, default='a2c', choices=['dqn', 'a2c'])
parser.add_argument('--nproc', type=int, default=1, help='# parallel processes')
parser.add_argument('--lr', type=float, default=1e-5, help='learning rate')  # 0.0001
parser.add_argument('--train_freq', type=int, default=20, help='train every train_freq iterations')  # 16
parser.add_argument('--train_start', type=int, default=0, help='train start iteration')  # 5000
parser.add_argument('--batch_size', type=int, default=16, help='SGD batch size')  # 128
parser.add_argument('--discount', type=float, default=0.996, help='discount factor (or gamma)')
parser.add_argument('--replay_size', type=int, default=20000, help='for DQN, replay buffer size')
parser.add_argument('--target_update', type=int, default=1000,
                    help='for DQN, update target net every target_update grad steps')
# Difference between -v0, -v4, Deterministic, NoFrameSkip, etc. on https://github.com/openai/gym/issues/1280
parser.add_argument('--eps_decay', type=int, default=3000,
                    help='for DQN, epsilon-greedy exploration over eps_decay iterations')  # 3000
parser.add_argument('--ent_coef', type=float, default=0.1, help='for A2C, coefficient of entropy')
parser.add_argument('--value_coef', type=float, default=0.5, help='for A2C, coefficient of value function loss')
parser.add_argument('--print_freq', type=int, default=200)  # 500
parser.add_argument('--checkpoint_freq', type=int, default=5000)  # 500000
parser.add_argument('--save_dir', type=str, default='breakout/')
parser.add_argument('--log', type=str, default='log_testing.txt', help='the log file name, appended to save_dir')
parser.add_argument('--load', type=str, default='')
parser.add_argument('--parallel_env', type=int, default=1)
parser.add_argument('--data_dir', type=str, default='data_eval/')
parser.add_argument('--sim', type=int, default=0,
                    help='simulate from apoapsis to apoapsis (sim =0) or from AI to AE with ABTS(sim=1) or from AI to AE with closedform solution(sim=2) ')
parser.add_argument('--phase',type=str,default='Campaign',choices=['Walkout','Endgame','Main','Campaign'])
parser.add_argument('--nominal', type=int, default = 1)
parser.add_argument('--output', type=str, default='NN',
                    choices=['NN', 'LSTM'])
args_DRL = parser.parse_args()


# args_DRL.load = load


 # set up SIGTERM and SIGKILL handler to save the model
import sys, signal


def exp_moving_average(x, new):
    # x != x means x is nan
    if x == x:
        return 0.9 * x + 0.1 * new
    else:
        return new




## Start Testing
#
### Testing
args_DRL.sim = 0
args_DRL.test = 1
if args_DRL.nominal ==1:
    args_DRL.log = 'log_testing_Ody_duststorm.txt'


# import mpi4py
# mpi4py.rc.initialize = False  # do not initialize MPI automatically
# mpi4py.rc.finalize = False  # do not finalize MPI automatically
from mpi4py import MPI

## establish connection
# MPI.Init()
cmd = "/sbin/ifconfig"
out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE).communicate()
# ip = str(out).split("inet addr:")[1].split()[0]
ip = str(out).split("inet")[1].split()[0]

name = MPI.Get_processor_name()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

num_nodes = int(comm.Get_size())

ip = comm.gather(ip)

if rank != 0:
    ip = None

ip = comm.bcast(ip, root=0)

os.environ['MASTER_ADDR'] = ip[0]
os.environ['MASTER_PORT'] = '2222'

backend = 'mpi'
rank = comm.Get_rank()

## Initialize Env and Alg

# args_DRL.env = 'CartPole-v1'
args_DRL.env = 'Aerobraking'
if args_DRL.env == 'Aerobraking':
    env = make_env.AerobrakingEnv(args_DRL)
else:
    import gym
    env = gym.make(args_DRL.env)

if args_DRL.phase == 'Walkout':
    r_norm = 5000000
elif args_DRL.phase == 'Endgame':
    r_norm = 10100000
elif args_DRL.phase == 'Main':
    r_norm = 10100000
elif args_DRL.phase == 'Campaign':
    r_norm = 30000000
endloop = 0
if rank == 0:

    ## Logprint
    if args_DRL.save_dir[-1] != '/':
        args_DRL.save_dir += '/'
    try:
        os.makedirs(args_DRL.save_dir)
    except:
        pass

    if len(args_DRL.log):
        # custom print function, line buffering, bad coding habit but works...
        log_file = open(args_DRL.save_dir + args_DRL.log, 'w', buffering=1)


        def logprint(*args, **kwargs):
            print(*args, **kwargs)
            print(*args, **kwargs, file=log_file)

    else:
        logprint = print
    logprint('RANK', rank)
    logprint('N Nodes', num_nodes)
    logprint(args_DRL)


    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logprint('running on device', device)
    if args_DRL.algo == 'dqn':
        algo = Algo.DQN(
            env.observation_space, env.action_space,
            args_DRL.lr, args_DRL.replay_size, args_DRL.batch_size, args_DRL.discount, args_DRL.target_update,
            args_DRL.eps_decay,args_DRL.train_start,
            device, args_DRL.load)
    elif args_DRL.algo == 'a2c':
        #            args_DRL.train_freq = 64
        assert (num_nodes - 1) * args_DRL.train_freq == args_DRL.batch_size
        algo = Algo.ActorCritic(
            env.observation_space, env.action_space,
            args_DRL.lr, args_DRL.nproc, args_DRL.train_freq, args_DRL.discount, args_DRL.ent_coef,
            args_DRL.value_coef,
            device, args_DRL.load)
    args_DRL.train_start = 1




    logprint('-- Start Testing --')
    # initialize some statistics
    avg_episode_len = float('nan')
    avg_episode_reward = float('nan')
    avg_raw_episode_reward = float('nan')
    avg_episode_maneuvern = float('nan')
    avg_episode_thermalviolations = float('nan')
    n_episodes_completed = 0
    n_reached_goal = 0
    n_not_reached_goal = 0
    n_impact = 0
    n_outdragpassage = 0
    episode_thermalviolations = 0

    loss = float('nan')
    time_start = time.time()

    # initialize state and requests
    obses = env.observation_space
    req = [None] * num_nodes
    for ii in range(1, num_nodes):
        req[ii] = comm.irecv(source=ii)

    list = [*range(1, num_nodes, 1)]
    while not endloop and len(list) != 0:
        for ii in list:
            data = req[ii].test()
            if data[0]:
                totdata = data[1]
                results = totdata[1]
                obses = totdata[0]
                action = totdata[2]
                new_obses = totdata[3]
                first_action = totdata[4]

                new_action = algo.act(new_obses, args_DRL.test)
                temp = [new_action, endloop]
                seq = comm.isend(temp, dest=ii)
                seq.wait()


                if first_action:
                    algo.observe(obses, action, results)

                    # evaluate loss
                    if algo.num_act_steps % args_DRL.train_freq == 0:
                        loss = algo.test()
                # algo.observe(obses, action, results)

                # algo.test()

                # update statistics
                if len(totdata) > 5:
                    # print('Update Statistics')
                    n_episodes_completed += 1
                    stat = totdata[5]
                    episode_lens = stat[0]
                    episode_rewards = stat[1]
                    episode_maneuvern = stat[3]
                    n_reached_goal += stat[2]
                    n_impact += stat[4]
                    n_outdragpassage += stat[5]
                    if not bool(stat[2] or stat[4] or stat[5]):
                        n_not_reached_goal +=1 # successful but not reached goal
                    episode_thermalviolations = stat[6]
                    distance_goal = stat[7]
                    time_mission = stat[8]
                    deltav = stat[9]
                    episode_maxrate = stat[10]
                    episode_actions = stat[11]
                    episode_ra = stat[12]
                    episode_hp = stat[13]
                    episode_omega = stat[14]
                    corridor_man = stat[15]
                    episode_reward = stat[16]
                    avg_episode_len = exp_moving_average(avg_episode_len, episode_lens)
                    avg_episode_reward = exp_moving_average(avg_episode_reward, episode_rewards)
                    avg_episode_maneuvern = exp_moving_average(avg_episode_maneuvern, episode_maneuvern)
                    avg_episode_thermalviolations = exp_moving_average(avg_episode_thermalviolations,
                                                                       episode_thermalviolations)

                    logprint(
                        "n_ep %5d |loss %6.3f |ep_len %6.1f |ep_rew %6.2f |maneuvers_number %6.1f |thermalv_number %5d |reached_goal %5d |impact_number %5d |out_of_dragpassage %5d | no_reached_goal %5d | distance_goal %6.2f |days_mission %6.2f |deltav %6.2f |heatrate %s |actions %s |ra %s |hp %s |omega %s |r %s |corridor maneuver %5d" %
                        (n_episodes_completed, loss, episode_lens, episode_rewards, episode_maneuvern,
                         episode_thermalviolations, n_reached_goal, n_impact, n_outdragpassage, n_not_reached_goal, distance_goal*1e-3, time_mission,deltav, episode_maxrate, episode_actions, episode_ra, episode_hp, episode_omega,episode_reward,corridor_man))
                req[ii] = comm.irecv(source=ii)
                # print('ask new request')

                if n_episodes_completed >= args_DRL.ntest:
                    list.remove(ii)
                    if len(list) == 0:
                        break
            else:
                pass

    if n_episodes_completed >= args_DRL.ntest:
        endloop = True
        print('end loop')
        for ii in range(1, num_nodes):
            results = req[ii].wait()
            action = algo.act(results[0])
            temp = [action, endloop]
            seq = comm.isend(temp, dest=ii)
            seq.wait()
    logprint('-- Complete Training --')

elif rank != 0:
    if args_DRL.phase == 'Walkout':
        raf = 3900 * 1e3
    elif args_DRL.phase == 'Endgame':
        raf = 4906 * 1e3
    elif args_DRL.phase == 'Main':
        raf = 3900 * 1e3
    elif args_DRL.phase == 'Campaign':
        raf = 3900 * 1e3
    # print('here')
    while not endloop:
        episode_maxrate = []
        episode_ra = []
        episode_action = []
        episode_hp = []
        episode_reward = []
        episode_omega = []
        episode_OMEGA = []
        episode_inc = []
        corridor_man = 0
        if args_DRL.env == 'Aerobraking':
            if args_DRL.env == 'Aerobraking':
                # initialization
                first_action = 0
                action_zero = 6 if args_DRL.algo == 'dqn' else 0
                obses = [env.reset(args_DRL.sim, rank +0)]
                ra = obses[0][1] * (r_norm - 3396200) + 3396200
                action = [action_zero]
                episode_ra.append(round(ra / 10 ** 3, 3))

                # first drag passage with action 0
                results = [env.step(action[0], args_DRL.sim)]
                if args_DRL.output == 'NN':
                    ra = results[0][0][1]*(r_norm - 3396200) + 3396200
                    hp = results[0][0][2]*(145000 - 85000) + 85000
                    omega = results[0][0][3]*(3.14 - 0) + 0
                    OMEGA = results[0][0][4]*(3.14 - 1.57) + 1.57
                    inc = results[0][0][5]*(1.75 - 1.39) + 1.39
                    maxrate = results[0][0][-1]*(0.5 - 0) + 0
                else:
                    maxrate = max(results[0][0][-1])*(0.5 - 0) + 0
                    ra = results[0][0][1][-1]*(r_norm - 3396200) + 3396200
                    hp = results[0][0][2][-1]*(145000 - 85000) + 85000
                    omega = results[0][0][3][-1]*(3.14 - 0) + 0
                    OMEGA = results[0][0][4][-1]*(3.14 - 1.57) + 1.57
                    inc = results[0][0][5][-1]*(1.75 - 1.39) + 1.39
                episode_maxrate.append(round(maxrate, 3))
                episode_action.append(action_zero)
                episode_ra.append(round(ra / 10 ** 3, 3))
                episode_hp.append(round(hp / 10 ** 3, 3))
                episode_omega.append(round(math.degrees(omega), 2))
                episode_OMEGA.append(round(math.degrees(OMEGA), 2))
                episode_inc.append(round(math.degrees(inc), 2))

                # Second step only if really far from corridor
                # if far out the drag passage, perform maximum action to get reachable
                if maxrate > 0.4 or maxrate < 0.025:
                    if maxrate > 0.4:
                        action = [12]
                    elif maxrate < 0.025:
                        action = [0]
                    results = [env.step(action[0], args_DRL.sim)]
                    if args_DRL.output == 'NN':
                        maxrate = results[0][0][-1] * (0.5 - 0) + 0
                        ra = results[0][0][1] * (r_norm - 3396200) + 3396200
                        hp = results[0][0][2] * (145000 - 85000) + 85000
                        omega = results[0][0][3] * (3.14 - 0) + 0
                        OMEGA = results[0][0][4] * (3.14 - 1.57) + 1.57
                        inc = results[0][0][5] * (1.75 - 1.39) + 1.39
                    else:
                        maxrate = max(results[0][0][-1]) * (0.5 - 0) + 0
                        ra = results[0][0][1][-1] * (r_norm - 3396200) + 3396200
                        hp = min(results[0][0][2]) * (145000 - 85000) + 85000
                        omega = results[0][0][3][-1] * (3.14 - 0) + 0
                        OMEGA = results[0][0][4][-1] * (3.14 - 1.57) + 1.57
                        inc = results[0][0][5][-1] * (1.75 - 1.39) + 1.39
                    episode_maxrate.append(round(maxrate, 3))
                    episode_action.append(action[0])
                    episode_ra.append(round(ra / 10 ** 3, 3))
                    episode_hp.append(round(hp / 10 ** 3, 3))
                    episode_omega.append(round(math.degrees(omega), 2))
                    episode_OMEGA.append(round(math.degrees(OMEGA), 2))
                    episode_inc.append(round(math.degrees(inc), 2))
                    corridor_man +=1
        else:
            obses = [env.reset()]  # to delete #prev step
            action = [0]
            results = [env.step(action[0])]  # to delete #next step

        new_obses = next(zip(*results))
        terminal = 0
        episode_lens = 0
        episode_rewards = 0
        episode_maneuvern = 0
        episode_impact = 0
        episode_outdragpassage = 0
        episode_thermalviolations = 0
        i = 1
        distance_goal = 6000*1e3
        # action = algo.act(obses)
        # results = [env.step(action[0], args_DRL.sim)]
        data = [obses, results, action, new_obses, first_action]
        while not terminal:
            # print('before seq rank', rank)
            seq = comm.isend(data, dest=0)
            seq.wait()
            req = comm.irecv(source=0)
            temp = req.wait()
            new_action = temp[0]
            # emergency maneuvers
            # print(distance_goal,maxrate)
            # if distance_goal < 150*1e3 and new_action[0] <2:
            #     new_action = [12]
            # if distance_goal > 100*1e3 and maxrate < 0.03:
            #     new_action = [0]
            # if distance_goal > 100*1e3 and maxrate > 0.4:
            #     new_action = [12]

            endloop = temp[1]
            # print('rank', rank, 'endloop', endloop)
            if endloop:
                # print('rank', rank, 'endloop', endloop)
                env.close()
                break

            if args_DRL.env == 'Aerobraking':
                results = [env.step(new_action[0], args_DRL.sim)]
            else:
                results = [env.step(new_action[0])]  # to delete

            r = results[0][1]
            terminal = results[0][2]
            i += 1
            episode_lens += 1
            episode_rewards += r
            first_action = 1

            obses, action = new_obses, new_action
            new_obses = next(zip(*results))

            if args_DRL.env == 'Aerobraking':
                if args_DRL.output == 'NN':
                    ra = results[0][0][1]*(r_norm - 3396200) + 3396200
                    hp = results[0][0][2]*(145000 - 85000) + 85000
                    maxrate = results[0][0][-1]*(0.5 - 0) + 0
                    omega = results[0][0][3] * (3.14 - 0) + 0
                    OMEGA = results[0][0][4] * (3.14 - 1.57) + 1.57
                    inc = results[0][0][5] * (1.75 - 1.39) + 1.39
                    time_mission = env._time_mission()
                    distance_goal = ra - raf
                    # print('ra',ra*1e-3,'hp',hp*1e-3,'max rate',maxrate, 'time',time_mission/(60*60*24))
                else:
                    time_mission = env._time_mission()
                    ra = results[0][0][1][-1]*(r_norm - 3396200) + 3396200
                    hp = results[0][0][2][-1]*(145000 - 85000) + 85000
                    maxrate =max(results[0][0][-1])*(0.5 - 0) + 0
                    omega = results[0][0][3][-1] * (3.14 - 0) + 0
                    OMEGA = results[0][0][4][-1] * (3.14 - 1.57) + 1.57
                    inc = results[0][0][5][-1] * (1.75 - 1.39) + 1.39

                if maxrate > 0.25:
                    episode_thermalviolations += 1
                episode_maxrate.append(round(maxrate,3))
                episode_action.append(action[0])
                episode_ra.append(round(ra/10**3, 3))
                episode_hp.append(round(hp/10**3, 3))
                episode_reward.append(round(r,3))
                episode_omega.append(round(math.degrees(omega), 2))
                episode_OMEGA.append(round(math.degrees(OMEGA), 2))
                episode_inc.append(round(math.degrees(inc), 2))

            if abs(action[0] - action_zero)>=0.001:
                episode_maneuvern += 1
            if terminal:
                if args_DRL.env == 'Aerobraking':
                    delta_v = env._set_action(action_zero)
                    episode_reachedgoal = 0
                    # print('hp',hp*1e-3)
                    distance_goal = ra - raf
                    if abs(distance_goal) <=  10*1e3:
                        episode_reachedgoal = 1
                    elif distance_goal<-10*1e3 and maxrate <0.5:
                        episode_reachedgoal = 0

                    elif (distance_goal<-10*1e3 and maxrate >0.5) or hp<85*1e3:
                        episode_impact = 1
                        episode_reachedgoal = 0
                    elif (maxrate < 0.01 and distance_goal>100*1e3) or hp>145*1e3:
                        episode_outdragpassage = 1
                        episode_reachedgoal = 0
                else:
                    episode_reachedgoal = 0  # to delete
                data = [obses, results, action, new_obses,first_action,
                        [episode_lens, episode_rewards, episode_reachedgoal, episode_maneuvern, episode_impact,
                         episode_outdragpassage, episode_thermalviolations, distance_goal, time_mission,delta_v,
                         episode_maxrate, episode_action, episode_ra, episode_hp,episode_omega,corridor_man,episode_reward]]

                seq = comm.isend(data, dest=0)
                seq.wait()
                req = comm.irecv(source=0)
                temp = req.wait()
                break
            else:
                data = [obses, results, action, new_obses, first_action]
# MPI.Finalize()
# env.close()
# quit()
