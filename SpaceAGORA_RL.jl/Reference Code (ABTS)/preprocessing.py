import csv
network = 'LSTM'
replaymemlog = open('breakout/replaymemlog'+network+'.csv', 'r')
reader = csv.reader(replaymemlog,delimiter='\t')
results_max = [[], [], [], [], [], []]
results_min = [[], [], [], [], [], []]
for row in reader:
    obses = row[0].split()
    action = row[1].split()

    if network== 'NN':
        results = row[2].split('[')
        results = results[1].split(']')
        results = results[0].split()
        obses = [float(i) for i in obses[1:-1]]
        action = int(action[0])
        results = [float(i) for i in results]
    elif network == 'LSTM':
        obses = row[0].translate({ ord("["): None })
        obses = obses.translate({ord("]"): None})
        obses = obses.split()
        ob = [float(i) for i in obses]

        results = row[2].translate({ ord("["): None })
        results = results.translate({ord("]"): None})
        results = results.split()
        tr = [float(i) for i in results]

        # obses = [[],[],[],[],[],[],[],[],[]]
        # results = [[],[],[],[],[],[],[],[],[]]


        count = 0
        # print(tr[0:100])
        # exit()
        for ii in range(0,6):
            results_min[ii].append(min(tr[ii*100:(ii+1)*100]))
            results_max[ii].append(max(tr[ii*100:(ii+1)*100]))
            # results_max[ii].append(max(tr[ii]))
            count += 1


## Normalize_data
if network == network == 'LSTM':
    # print(results_min)

    min_data = [min(results_min[0]), min(results_min[1]), min(results_min[2]), min(results_min[3]), min(results_min[4]),
                        min(results_min[5])]
    max_data = [max(results_max[0]), max(results_max[1]), max(results_max[2]), max(results_max[3]), max(results_max[4]),
                    max(results_max[5])]

print(min_data,max_data)

        # min_data = [min(results[0]),min(results[1]),min(results[2]),min(results[3]), min(results[4]),
        #             min(results[5]),min(results[6]),min(results[7]),min(results[8])]
        # max_data = [max(results[0]),max(results[1]),max(results[2]),max(results[3]),max(results[4]),
        # #             max(results[5]),max(results[6]),max(results[7]),max(results[8])]
        # action = int(action[0])
        #
        # count = 0
        # for ii in range(0,6):
        #     for jj in range(0,100):
        #         obses[ii].append((ob[count]-min_data[ii])/(max_data[ii]-min_data[ii]))
        #         results[ii].append(tr[count]-min_data[ii])/(max_data[ii]-min_data[ii])
        #         count += 1

