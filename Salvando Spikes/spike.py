import matplotlib.pyplot as plt
import pandas as pd
import json

file1 = pd.read_csv('1neuronI280.csv', sep=",", names=['time','v','m','n','h','spike'])
file2 = pd.read_csv('1neuronI300.csv', sep=",", names=['time','v','m','n','h','spike'])

size = len(file1['time'])

plt.xlabel("t(ms)")
plt.ylabel("Potencial (V)")
plt.title("Período de atividade para um neurônio")

"""
spikes: contains the instant t in which the spike occurs
period: contains the period between two consecutive spikes
period list is calculated from spikes list
"""

spikes1 = []
spikes2 = []
period1 = []
period2 = []

for i in range (0,size):
    if file1['spike'][i] == 1:
        plt.scatter(file1['time'][i-1], file1['v'][i-1], c="lightblue")
        spikes1.append(file1['time'][i-1])
    if file2['spike'][i] == 1:
        plt.scatter(file2['time'][i-1], file2['v'][i-1], c="lightcoral")
        spikes2.append(file2['time'][i-1]) 

plt.plot(file1['time'] ,file1['v'], "cornflowerblue", label = "280pA")
plt.plot(file2['time'] ,file2['v'], "brown", label = "300pA")
plt.legend(loc = "upper right")
#plt.savefig("spike.pdf")
plt.show()

for i in range (1,len(spikes1)):
    period1.append(spikes1[i] - spikes1[i-1])

for i in range (1,len(spikes2)):
    period2.append(spikes2[i] - spikes2[i-1])

#calculando o delay entre as duas correntes
size_delay = min([len(spikes1), len(spikes2)])
delay = []
for i in range (0,size_delay):
    delay.append(spikes1[i] - spikes2[i])

print('Current 280pA')
print(spikes1)
print()
print(period1)
print()
print('Current 300pA')
print(spikes2)
print()
print(period2)
print()
print('Delay')
print(delay)

resultados = {
    "spikes 1": spikes1,
    "period 1": period1,
    "spikes 2": spikes2,
    "period 2": period2,
    "delay": delay 
}

with open('./resultados.json', 'w') as fp:
    json.dump(resultados, fp, indent=6)