# Main file: train an agent with specific algorithm on specific environment
# if __name__ == '__main__':

import torch
import Algo
import make_env
import numpy as np
import argparse, os
import time
from mpi4py import MPI
import subprocess
import math

# import multiprocessing as mp

parser = argparse.ArgumentParser()
parser.add_argument('--niter', type=int, default=1000000, help='# training iterations')  # 100000
parser.add_argument('--ntest', type=int, default=100, help='# n of testing episode')  # 100
parser.add_argument('--test', type=int, default=0, help='# test mode')
parser.add_argument('--env', type=str, default='Aerobraking')
parser.add_argument('--algo', type=str, default='dqn', choices=['dqn', 'a2c'])
parser.add_argument('--nproc', type=int, default=1, help='# parallel processes')
parser.add_argument('--lr', type=float, default=1e-4, help='learning rate') #0.0001
parser.add_argument('--train_freq', type=int, default=5, help='train every train_freq iterations')#16
parser.add_argument('--train_start', type=int, default=10000, help='train start iteration')  # 5000
parser.add_argument('--batch_size', type=int, default=256, help='SGD batch size')  # 128
parser.add_argument('--discount', type=float, default=0.95, help='discount factor (or gamma)')
parser.add_argument('--replay_size', type=int, default=50000, help='for DQN, replay buffer size')
parser.add_argument('--target_update', type=int, default=10000,
                    help='for DQN, update target net every target_update grad steps')
# Difference between -v0, -v4, Deterministic, NoFrameSkip, etc. on https://github.com/openai/gym/issues/1280
parser.add_argument('--eps_decay', type=int, default=5000,
                    help='for DQN, epsilon-greedy exploration over eps_decay iterations') #3000
parser.add_argument('--ent_coef', type=float, default=0.1, help='for A2C, coefficient of entropy')
parser.add_argument('--value_coef', type=float, default=0.5, help='for A2C, coefficient of value function loss')
parser.add_argument('--print_freq', type=int, default=200)  # 500
parser.add_argument('--checkpoint_freq', type=int, default=1000)  # 500000
parser.add_argument('--save_dir', type=str, default='breakout/')
parser.add_argument('--log', type=str, default='log.txt', help='the log file name, appended to save_dir')
parser.add_argument('--load', type=str, default='')
parser.add_argument('--parallel_env', type=int, default=1)
parser.add_argument('--data_dir', type=str, default='data/')
parser.add_argument('--sim', type=int, default=0,
                    help='simulate from apoapsis to apoapsis (sim =0) or from AI to AE with ABTS(sim=1) or from AI to AE with closedform solution(sim=2) ')
parser.add_argument('--phase',type=str,default='Main',choices=['Walkout','Endgame','Main','Campaign'])
parser.add_argument('--nominal', type=int, default=0)
parser.add_argument('--continuetraining', type=int, default=0)
parser.add_argument('--replaycreation', type=int, default=0, help = '0: no creation, 1:creation of the observation txt file, 2:read of obs file and creation of replay memory')
parser.add_argument('--output', type=str, default='NN',
                    choices=['NN','LSTM'])
parser.add_argument('--initialrank',type=int,default=0)

# args_DRL = parser.parse_args()
args_DRL , unknown = parser.parse_known_args()
## establish connection

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
comm = MPI.COMM_WORLD


## Initialize Env and Alg

# args_DRL.env = 'CartPole-v1'
# args_DRL.env = 'Pendulum-v0'
args_DRL.env = 'Aerobraking'
if args_DRL.env == 'Aerobraking':
    # if args_DRL.replaycreation == 2:
        # args_DRL.train_start  = 1
    env = make_env.AerobrakingEnv(args_DRL)
else:
    import gym
    env = gym.make(args_DRL.env)


episode_stop = 100000
endloop = False
test = False

if args_DRL.phase == 'Walkout':
    r_norm = 5000*10**3
elif args_DRL.phase == 'Endgame':
    r_norm = 10100*10**3
elif args_DRL.phase == 'Main':
    r_norm = 10100*10**3
elif args_DRL.phase == 'Campaign':
    r_norm = 30000*10**3

## Training
if rank == 0:
    # dir = '/content/drive/Othercomputers/My MacBook Pro/Research/DRLAerobraking/breakout'
    dir = args_DRL.save_dir
    if not args_DRL.continuetraining:
        for f in os.listdir(dir):
            if f != 'replaymemlogLSTM.csv' and f != 'replaymemlogNN.csv':
                os.remove(os.path.join(dir, f))
    else:
        file = []
        all_files = os.listdir(dir)
        for filename in all_files:
            if filename.lower().endswith('.pth'):
                file.append(int(filename.partition('.')[0]))
        file.sort()
        args_DRL.load = args_DRL.save_dir + str(file[-1]) + '.pth'

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
        print(log_file)


        def logprint(*args, **kwargs):
            print(*args, **kwargs)
            print(*args, **kwargs, file=log_file)

    else:
        logprint = print

    if args_DRL.replaycreation==1:
        # print('here')
        replaymemlog = open(args_DRL.save_dir + 'replaymemlogLSTM.csv', 'a', buffering=1)
        # def replayprint(*args, **kwargs):
        #     import csv
        #     csv.writer(replaymemlog).writerows()
            # print(*args, **kwargs, file=replaymemlog)

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
            device)
    elif args_DRL.algo == 'a2c':
        args_DRL.train_freq = 10
        args_DRL.batch_size = int(args_DRL.train_freq*(num_nodes-1))
        logprint('batch size',args_DRL.batch_size)
        assert (num_nodes-1) * args_DRL.train_freq == args_DRL.batch_size
        algo = Algo.ActorCritic(
            env.observation_space, env.action_space,
            args_DRL.lr, args_DRL.nproc, args_DRL.train_freq, args_DRL.discount, args_DRL.ent_coef, args_DRL.value_coef, args_DRL.train_start,
            device)
        args_DRL.train_start = 1

    if len(args_DRL.load):
        algo.load(args_DRL.load)
        logprint("load pretrained model from", args_DRL.load)

    ## if replay memory file populated
    if args_DRL.replaycreation == 2:
        import csv
        replaymemlog = open(args_DRL.save_dir + 'replaymemlogLSTM.csv', 'r')
        reader = csv.reader(replaymemlog,delimiter='\t')
        for row in reader:
            obses = row[0].split()
            action = row[1].split()

            if args_DRL.output == 'NN':
                results = row[2].split('[')
                results = results[1].split(']')
                results = results[0].split()
                obses = [float(i) for i in obses[1:-1]]
                action = int(action[0])
                results = [float(i) for i in results]
            elif args_DRL.output == 'LSTM':
                obses = row[0].translate({ ord("["): None })
                obses = obses.translate({ord("]"): None})
                obses = obses.split()
                ob = [float(i) for i in obses]

                results = row[2].translate({ ord("["): None })
                results = results.translate({ord("]"): None})
                results = results.split()
                tr = [float(i) for i in results]

                obses = [[],[],[],[],[],[]]
                results = [[],[],[],[],[],[]]
                count = 0
                for ii in range(0,6):
                    for jj in range(0,100):
                        obses[ii].append(ob[count])
                        results[ii].append(tr[count])
                        count += 1
                action = int(action[0])

            t = env._terminal(results, args_DRL.sim)
            a = env._set_action(action, i=1)
            r = env._reward(results, args_DRL.sim, [a])
            ## convert in format for algorithm
            trans = (np.array(obses), r, t, {})
            algo.observe(np.array([obses]), [action], [trans])
    # exit()
    logprint('-- Start Training --')
    # initialize some statistics
    n_reached_goal = 0
    n_not_reached_goal = 0
    n_episodes_completed = 0
    n_impact = 0
    n_outdragpassage = 0
    episode_thermalviolations = 0
    loss = float('nan')
    time_start = time.time()

    # set up SIGTERM and SIGKILL handler to save the model
    import sys, signal

    sig_reenter = False


    def sigterm_handler(signal, frame):
        global sig_reenter
        if not sig_reenter:
            sig_reenter = True
            algo.save(args_DRL.save_dir + str(algo.num_act_steps) + '.pth')
            logprint('sigterm received, save checkpoint to', args_DRL.save_dir + str(it) + '.pth')
            env.close()
            if len(args_DRL.log): log_file.close()
        else:
            logprint('forced exit')
        sys.exit(0)


    signal.signal(signal.SIGTERM, sigterm_handler)
    signal.signal(signal.SIGINT, sigterm_handler)


    def exp_moving_average(x, new):
        # x != x means x is nan
        if x == x:
            return 0.9 * x + 0.1 * new
        else:
            return new


    # buf = bytearray(1 << 30)
    obses = env.observation_space
    req = [None] * num_nodes
    for ii in range(1, num_nodes):
        req[ii] = comm.irecv(source=ii)

    n_episodes_completed_pre =0
    n_episodes_completed_test = 0
    n_episodes_completed_pre_test = 0
    while not endloop:
        for ii in range(1, num_nodes):
            data = req[ii].test()
            if data[0]:
                totdata = data[1]
                obses = totdata[0]
                results = totdata[1]
                action = totdata[2]
                new_obses= totdata[3]
                first_action = totdata[4]

                new_action = algo.act(new_obses)
                temp = [new_action, endloop]
                seq = comm.isend(temp, dest=ii)
                seq.wait()

                distance_goal = np.nan
                time_mission = np.nan
                deltav = np.nan
                # update statistics
                if len(totdata) > 5:
                    stat = totdata[5]
                    #print('Update Statistics')
                    n_episodes_completed += 1
                    n_reached_goal += stat[2]
                    n_impact += stat[4]
                    n_outdragpassage += stat[5]
                    if not bool(stat[2] or stat[4] or stat[5]):
                        n_not_reached_goal += 1

                    episode_lens = stat[0]
                    episode_rewards = stat[1]
                    episode_manevuern = stat[3]
                    episode_thermalviolations = stat[6]
                    distance_goal = stat[7]
                    time_mission = stat[8]
                    deltav = stat[9]
                    episode_maxrate = stat[10]
                    episode_action = stat[11]
                    episode_ra = stat[12]
                    episode_hp = stat[13]
                    episode_omega = stat[14]
                    # episode_OMEGA = stat[15]
                    # episode_inc = stat[16]
                    corridor_man = stat[15]
                    # episode_reward = stat[17]
                    rank = stat[-1]

                req[ii] = comm.irecv(source=ii)


                if first_action:
                    if args_DRL.replaycreation == 1:  # store all observation, action, results, to create replaymemory
                        import csv

                        temp1 = str(obses[0])
                        temp2 = str(results[0][0])
                        csv.writer(replaymemlog, delimiter='\t').writerow([temp1, str(action[0]), temp2])

                    algo.observe(obses, action, results)

                    # evaluate loss
                    if algo.num_act_steps % args_DRL.train_freq == 0 and algo.num_act_steps > args_DRL.train_start:
                        loss = algo.train()


                it = algo.num_act_steps
                if n_episodes_completed > n_episodes_completed_pre or algo.num_act_steps >= args_DRL.niter:#n_episodes_completed >= args_DRL.niter: #algo.num_act_steps % args_DRL.print_freq == 0
                    n_episodes_completed_pre = n_episodes_completed
                    n_episodes = n_episodes_completed
                    n_goal = n_reached_goal
                    n_imp = n_impact
                    n_out = n_outdragpassage
                    n_notgoal = n_not_reached_goal
                    time_now = time.time()
                    time_passed = time_now - time_start
                    time_remain = time_passed * (args_DRL.niter - 1 - it) / it
                    time_last = time_now
                    time_passed = '%02d:%02d' % (time_passed // (60*60), time_passed % (60*60))
                    time_remain = '%02d:%02d' % (time_remain // (60*60), time_remain % (60*60))
                    # print( distance_goal*1e-3, time_mission/(60*60*24))
                    logprint(
                            "iter %6d |rank %5d |loss %6.3f |n_ep %5d |ep_len %6.1f |ep_rew %6.2f |num_manevuers %6.1f |num_thermalviolations %5d | num_reachedgoal %5d | num_impact %5d | num_outdragpassage %5d |no_reached_goal %5d | distance_goal %6.2f |days_mission %6.2f  |deltav %6.2f " \
                            "|time %s rem %s |heatrate %s |actions %s |ra %s |hp %s |omega %s |corridor maneuver %5d" %
                            (it, rank, loss, n_episodes, episode_lens, episode_rewards, episode_manevuern, episode_thermalviolations,
                             n_goal, n_imp,n_out, n_notgoal, distance_goal*1e-3, time_mission, deltav, time_passed, time_remain,
                             episode_maxrate,episode_action,episode_ra,episode_hp,episode_omega,corridor_man))
                    # else:
                    #     logprint(
                    #         "TEST: iter %6d |loss %6.3f |n_ep %5d |ep_len %6.1f |ep_rew %6.2f |num_manevuers %6.1f |num_thermalviolations %5d | num_reachedgoal %5d | num_impact %5d | num_outdragpassage %5d | no_reached_goal %5d | distance_goal %6.2f |days_mission %6.2f |deltav %6.2f " \
                    #         "|time %s rem %s |heatrate %s |actions %s" %
                    #         (it, loss, n_episodes, episode_lens, episode_rewards, episode_manevuern,
                    #          episode_thermalviolations,
                    #          n_goal, n_imp, n_out, n_notgoal, distance_goal*1e-3, time_mission, deltav, time_passed, time_remain,episode_maxrate,episode_action))
                    #     if n_episodes_completed_test == 10:
                    #         print('END TESTING')
                    #         test = False
                    #         n_episodes_completed_pre_test = 0
                    #         n_episodes_completed_test = 0

                if it % args_DRL.checkpoint_freq == 0:
                    # logprint('save checkpoint to', args_DRL.save_dir + str(it) + '.pth')
                    algo.save(args_DRL.save_dir + str(it) + '.pth')
                    files_in_directory = os.listdir(args_DRL.save_dir)
                    filtered_files = [file for file in files_in_directory if
                                      file.endswith(".pth") and not file == str(algo.num_act_steps) + '.pth' and int(file.split('.')[0]) % 5000 != 0]
                    for file in filtered_files:
                        path_to_file = os.path.join(args_DRL.save_dir, file)
                        os.remove(path_to_file)


            else:
                pass

        if algo.num_act_steps == args_DRL.niter:
            episode_stop = n_episodes_completed
        if n_episodes_completed > episode_stop:
            endloop = True
            print('end loop')
            logprint('save checkpoint to', args_DRL.save_dir + str(algo.num_act_steps) + '.pth')
            algo.save(args_DRL.save_dir + str(algo.num_act_steps) + '.pth')
            files_in_directory = os.listdir(args_DRL.save_dir)
            filtered_files = [file for file in files_in_directory if file.endswith(".pth") and not file == str(algo.num_act_steps) + '.pth' and int(file.split('.')[0]) % 5000 != 0]
            for file in filtered_files:
                path_to_file = os.path.join(args_DRL.save_dir, file)
                os.remove(path_to_file)
            for ii in range(1, num_nodes):
                results = req[ii].wait()
                action = algo.act(results[0])
                temp = [action, endloop]
                seq = comm.isend(temp, dest=ii)
                seq.wait()
    logprint('-- Complete Training --')

elif rank != 0:
    it = 0
    rep = 0
    count_repetition = 0
    while not endloop:
        terminal = 0
        episode_lens = 1
        episode_rewards = 0
        episode_maneuvern = 0
        episode_impact = 0
        episode_outdragpassage = 0
        episode_thermalviolations = 0
        episode_maxrate = []
        episode_action = []
        episode_ra = []
        episode_hp = []
        episode_reward = []
        episode_omega = []
        episode_OMEGA = []
        episode_inc = []
        corridor_man = 0
        i = 1
        if args_DRL.env == 'Aerobraking':
            # initialization
            first_action = 0
            action_zero = 6 if args_DRL.algo == 'dqn' else 0
            obses = [env.reset(args_DRL.sim, rank+args_DRL.initialrank,rep)]
            ra = obses[0][1] * (r_norm - 3396200) + 3396200
            action = [action_zero]
            episode_ra.append(round(ra / 10 ** 3, 3))

            # first drag passage with action 0
            results = [env.step(action[0],args_DRL.sim)]
            if args_DRL.output == 'NN':
                ra = results[0][0][1] * (r_norm - 3396200) + 3396200
                hp = results[0][0][2] * (145000 - 85000) + 85000
                omega = results[0][0][3] * (3.14 - 0) + 0
                OMEGA = results[0][0][4] * (3.14 - 1.57) + 1.57
                inc = results[0][0][5] * (1.75 - 1.39) + 1.39
                maxrate = results[0][0][-1]*(0.5 - 0) + 0
            else:
                maxrate = max(results[0][0][-1])*(0.5 - 0) + 0
                ra = results[0][0][1][-1] * (r_norm - 3396200) + 3396200
                hp = min(results[0][0][2]) * (145000 - 85000) + 85000
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
                if maxrate >0.4:
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
                corridor_man += 1
        else:
            obses = [env.reset()] # to delete #prev step
            action = [0]
            results = [env.step(action)]#[env.step(action[0])] # to delete #next step

        new_obses = next(zip(*results))
        data = [obses,results,action,new_obses,first_action]
        while not terminal:
            seq = comm.isend(data, dest=0)
            seq.wait()
            req = comm.irecv(source=0)
            temp = req.wait()
            new_action = temp[0]


            endloop = temp[1]
            if endloop:
                env.close()
                break


            if args_DRL.env == 'Aerobraking':
                results = [env.step(new_action[0], args_DRL.sim)]
            else:
                results = [env.step(new_action)]#[env.step(new_action[0])] # to delete
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
                    # print('ra',ra*1e-3,'hp',hp*1e-3,'max rate',maxrate, 'time',time_mission/(60*60*24))
                else:
                    time_mission = env._time_mission()
                    ra = results[0][0][1][-1]*(r_norm - 3396200) + 3396200
                    hp = results[0][0][2][-1]*(145000 - 85000) + 85000
                    maxrate =max(results[0][0][-1])*(0.5 - 0) + 0
                    omega = results[0][0][3][-1] * (3.14 - 0) + 0
                    OMEGA = results[0][0][4][-1] * (3.14 - 1.57) + 1.57
                    inc = results[0][0][5][-1] * (1.75 - 1.39) + 1.39
                    # print('ra',ra*1e-3,'hp',hp*1e-3,'max rate',maxrate, 'time',time_mission/(60*60*24))
                if maxrate > 0.25:
                    episode_thermalviolations += 1
                episode_maxrate.append(round(maxrate,3))
                episode_action.append(action[0])
                episode_ra.append(round(ra / 10 ** 3, 3))
                episode_hp.append(round(hp / 10 ** 3, 3))
                episode_reward.append(round(r,3))
                episode_omega.append(round(math.degrees(omega), 2))
                episode_OMEGA.append(round(math.degrees(OMEGA), 2))
                episode_inc.append(round(math.degrees(inc), 2))

            if args_DRL.env == 'Aerobraking': 
              if abs(action[0] - action_zero)>=0.001:
                  episode_maneuvern += 1
            if terminal:

                print('Terminal Rank', rank, 'episode length', episode_lens)
                if args_DRL.env == 'Aerobraking':
                    if args_DRL.phase == 'Walkout':
                        raf = 3900 * 1e3
                    elif args_DRL.phase == 'Endgame':
                        raf = 4906 * 1e3
                    elif args_DRL.phase == 'Main':
                        raf = 3900 * 1e3
                    elif args_DRL.phase == 'Campaign':
                        raf = 3900 * 1e3
                    distance_goal = ra - raf
                    print('ra',ra,'raf',raf,'distance goal',distance_goal)
                    # print('ra', ra*1e-3, 'raf', raf*1e-3, abs(ra - raf)*1e-3, 'hp', hp/1e3)
                    if abs(distance_goal) <= 10*1e3:
                        episode_reachedgoal = 1
                    else:
                        episode_reachedgoal = 0

                    if maxrate > 0.25:
                        episode_impact = 1
                    elif maxrate < 0.05 and abs(distance_goal)>100:
                        episode_outdragpassage = 1
                else:
                    episode_reachedgoal = 0 # to delete
                    distance_goal = 0
                    delta_v = 0
                    time_mission = 0
                if args_DRL.env == 'Aerobraking':
                    delta_v = env._set_action(action_zero)
                    if episode_reachedgoal == 1 and count_repetition == 0:
                        rep = 1
                        count_repetition += 1
                        print('Firts if count rep',count_repetition,'rep',rep)
                    elif rep ==1 and count_repetition<9:
                        count_repetition += 1
                        print('Second if count rep',count_repetition,'rep',rep)
                    elif rep == 1 and count_repetition>=9:
                        rep = 0
                        count_repetition = 0
                        print('Third if count rep',count_repetition,'rep',rep)
                data = [obses, results, action, new_obses, first_action,
                        [episode_lens, episode_rewards, episode_reachedgoal, episode_maneuvern, episode_impact,
                         episode_outdragpassage, episode_thermalviolations, distance_goal, time_mission,delta_v,
                         episode_maxrate, episode_action, episode_ra, episode_hp,episode_omega,corridor_man, rank]]
                seq = comm.isend(data, dest=0)
                seq.wait()
                req = comm.irecv(source=0)
                temp = req.wait()
                break
            else:
                data = [obses,results,action,new_obses, first_action]

env.close()
