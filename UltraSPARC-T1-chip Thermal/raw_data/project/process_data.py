import numpy as np
import matplotlib.pyplot as plt
import csv

f = open('output2.txt','r')
count = 0

heat_map = []
heat_vec = []
for line in f.readlines():
    token = line.strip().split('  ')
    if token == ['']:
        heat_map.append(heat_vec)
        heat_vec = []
    else:
        token = list(map(float,token))
        heat_vec.extend(token)

    #count += 1
    #if count > 30:
    #    break
f.close()


with open('Data_chip.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile)
    for row in heat_map:
        spamwriter.writerow(row)


















