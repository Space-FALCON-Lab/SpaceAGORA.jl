# Main file: train an agent with specific algorithm on specific environment
# if __name__ == '__main__':

import os

import subprocess
import time

directory = 'breakout/'
all_files = os.listdir(directory)
file = []
for filename in all_files:
    if filename == 'log.txt':
        args_file = filename
    if filename.lower().endswith('.pth'):
        file.append(int(filename.partition('.')[0]))
file.sort()

# args definition
lines = open(directory+args_file, 'r').readlines()[2].split('(')[1]
par = lines.split(', ')
key = []
value = []
for comb in par:
    key.append(comb.split('=')[0])
    value.append(comb.split('=')[1])
print(key)
for i in range(len(key)):
    if key[i] =='algo':
        algo = value[i][1:-1]
    elif key[i] =='phase':
        phase = value[i][1:-1]
    elif key[i] =='env':
        env = value[i]
    elif key[i] == 'output':
        output = value[i][1:-1]
    elif key[i] == 'batch_size':
        batch_size = value[i]

for f in file:
    load = 'breakout/' + str(f) + '.pth'
    log = 'log_testing_'+str(f)+'.txt'
    print(load)
    subprocess.call(["mpiexec", "-n", "7", "python3", "evaluation_function.py", '--load', load,'--log',log, '--algo', algo, '--env', env, '--output',output, '--batch_size',batch_size, '--phase',phase])
    # subprocess.Popen(["mpiexec", "-n", "5", "python3", "evaluation_function.py", '--load', load,'--log',log], stdout=subprocess.PIPE, shell=False).wait()
    # p.wait()

