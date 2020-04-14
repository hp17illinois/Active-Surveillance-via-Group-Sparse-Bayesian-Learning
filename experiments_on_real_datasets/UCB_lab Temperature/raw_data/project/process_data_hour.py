import numpy as np
import matplotlib.pyplot as plt
import csv

f = open('temperature_data.txt','r')
count = 0

sensor_dict = {}
date_hour_min_log = {}
for line in f.readlines():
    try:
        date,time,epoch,moteid,temperature,humidity,light,voltage = line.strip().split(' ')
        temperature = float(temperature)
        moteid = int(moteid)
    except:
        continue
    md = date.split('-',1)[1]
    hour = time.split(':',1)[0]
    date_hour_id = md + ' ' + hour

    if temperature > 50 or temperature < 0:
        continue
    if moteid in [5,15,18]: # delete 3 sensors, 54-3=51
        continue

    if date_hour_id not in date_hour_min_log:
        date_hour_min_log[date_hour_id] = 50

    if moteid not in sensor_dict:
        sensor_dict[moteid] = {}
    if date_hour_id not in sensor_dict[moteid]:
        sensor_dict[moteid][date_hour_id] = []
    sensor_dict[moteid][date_hour_id].append(temperature)
    #count += 1
    #if count > 30:
    #    break
f.close()

sensor_ave_t_dict = {}
for s_id in sensor_dict:
    if s_id not in sensor_ave_t_dict:
        sensor_ave_t_dict[s_id] = {}
    for date_hour_id in sensor_dict[s_id]:
        if date_hour_id not in sensor_ave_t_dict:
            sensor_ave_t_dict[s_id][date_hour_id] = {'log_num':0,'ave_t':-1}
        t_list = sensor_dict[s_id][date_hour_id]
        sensor_ave_t_dict[s_id][date_hour_id]['log_num'] = len(t_list)
        sensor_ave_t_dict[s_id][date_hour_id]['ave_t'] = sum(t_list)/len(t_list)
        if sensor_ave_t_dict[s_id][date_hour_id]['log_num'] < date_hour_min_log[date_hour_id]:
            date_hour_min_log[date_hour_id] = sensor_ave_t_dict[s_id][date_hour_id]['log_num']

start_id = 160
end_id = 265

'''
# find the best data range
date_hour_list = []
date_hour_log_list = []
for date_hour_id in sorted (date_hour_min_log.keys()):
    date_hour_list.append(date_hour_id)
    date_hour_log_list.append(date_hour_min_log[date_hour_id])
    
print(date_hour_list[start_id:end_id])
print(date_hour_log_list[start_id:end_id])
plt.figure()
plt.plot(date_hour_list[start_id:end_id],date_hour_log_list[start_id:end_id])
plt.show()
'''

sensor_data_time_dict = {}
for s_id in sensor_ave_t_dict:
    for date_hour_id in sensor_ave_t_dict[s_id]:
        if date_hour_id > '03-05 16' and date_hour_id < '03-10 00':
            if date_hour_id not in sensor_data_time_dict:
                sensor_data_time_dict[date_hour_id] = {}
            if s_id not in sensor_data_time_dict[date_hour_id]:
                sensor_data_time_dict[date_hour_id][s_id] = sensor_ave_t_dict[s_id][date_hour_id]['ave_t']


with open('Data_lab.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile)
    for date_hour_id in sorted (sensor_data_time_dict.keys()):
        print(date_hour_id)
        t_list = []
        for s_id in sorted (sensor_data_time_dict[date_hour_id].keys()):
            print(s_id,)
            t_list.append(sensor_data_time_dict[date_hour_id][s_id])
        spamwriter.writerow(t_list)


















