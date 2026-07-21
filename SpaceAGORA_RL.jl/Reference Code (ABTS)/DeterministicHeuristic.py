import make_env_heuristic
import numpy as np
import argparse, os
import time
import math
parser = argparse.ArgumentParser()
parser.add_argument('--log', type=str, default='logHeuristic4.txt', help='the log file name, appended to save_dir')
parser.add_argument('--save_dir', type=str, default='breakout/')
parser.add_argument('--nominal', type=int, default=1)
parser.add_argument('--output', type=str, default='NN',
                    choices=['NN','LSTM'])
parser.add_argument('--data_dir', type=str, default='data/')
parser.add_argument('--phase',type=str,default='Main',choices=['Walkout','Endgame','Main','Campaign'])
parser.add_argument('--sim', type=int, default=0,
                    help='simulate from apoapsis to apoapsis (sim =0) or from AI to AE with ABTS(sim=1) or from AI to AE with closedform solution(sim=2) ')
args_DRL , unknown = parser.parse_known_args()
env = make_env_heuristic.AerobrakingEnvHeuristic(args_DRL)

if args_DRL.nominal ==1:
    args_DRL.log = 'log_testing_Ody_Heuristic3.txt'
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
time_start = time.time()
episode_reachedgoal = 0
episode_noreachedgoal =0
for n_episodes in range(1,41):
    deltav = 0
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
    episode_omega = []
    episode_OMEGA = []
    episode_inc = []
    episode_reward = []
    i = 1

    obses = [env.reset(args_DRL.sim, 2)]
    ra = obses[0][1]
    episode_ra.append(round(ra / 10 ** 3, 3))

    terminal = 0
    while not terminal:
        results = [env.step(args_DRL.sim)]
        maxrate = results[0][0][-1]
        episode_maxrate.append(round(maxrate, 3))
        r = results[0][1]
        action = results[0][2]
        terminal = results[0][3]
        episode_action.append(round(action,3))
        deltav += abs(action)
        i += 1
        episode_lens += 1
        episode_rewards += r
        ra = results[0][0][1]
        hp = results[0][0][2]
        omega = results[0][0][3]
        OMEGA = results[0][0][4]
        inc = results[0][0][5]
        episode_ra.append(round(ra/10**3,3))
        episode_hp.append(round(hp / 10 ** 3, 3))
        episode_reward.append(round(r, 3))
        episode_omega.append(round(math.degrees(omega),2))
        episode_OMEGA.append(round(math.degrees(OMEGA),2))
        episode_inc.append(round(math.degrees(inc),2))
        time_mission = env._time_mission()
        if maxrate > 0.25 or maxrate < 0.05:
            episode_thermalviolations += 1
        if abs(action) >= 0.001:
            episode_maneuvern += 1
        if terminal:
            time_now = time.time()
            time_passed = time_now - time_start
            time_remain = time_passed * (41 - 1 - n_episodes) / n_episodes
            time_last = time_now
            time_passed = '%02d:%02d' % (time_passed // 60, time_passed % 60)
            time_remain = '%02d:%02d' % (time_remain // 60, time_remain % 60)
            print('episode length', episode_lens)
            if args_DRL.phase == 'Walkout':
                raf = 3900 * 1e3
            elif args_DRL.phase == 'Endgame':
                raf = 4906 * 1e3
            elif args_DRL.phase == 'Main':
                raf = 3900 * 1e3
            elif args_DRL.phase == 'Campaign':
                raf = 3900 * 1e3
            distance_goal = ra - raf
            n_imp = 0
            n_out = 0

            # print('ra', ra*1e-3, 'raf', raf*1e-3, abs(ra - raf)*1e-3, 'hp', hp/1e3)
            if abs(distance_goal) <= 10 * 1e3:
                episode_reachedgoal += 1
                episode_noreachedgoal += 0
            else:
                episode_reachedgoal += 0
                episode_noreachedgoal += 1
            logprint(
                "n_ep %5d |ep_len %6.1f |ep_rew %6.2f |num_manevuers %6.1f |num_thermalviolations %5d | num_reachedgoal %5d | num_impact %5d | num_outdragpassage %5d |no_reached_goal %5d | distance_goal %6.2f |days_mission %6.2f  |deltav %6.2f " \
                "|time %s rem %s |heatrate %s |actions %s |ra %s |hp %s |omega %s |OMEGA %s |inc %s |r %s" %
                (n_episodes, episode_lens, episode_rewards, episode_maneuvern, episode_thermalviolations,
                 episode_reachedgoal, n_imp, n_out, episode_noreachedgoal, distance_goal * 1e-3, time_mission, deltav, time_passed, time_remain,
                 episode_maxrate, episode_action, episode_ra, episode_hp, episode_omega, episode_OMEGA, episode_inc, episode_reward))
