import matplotlib.pyplot as plt
import pandas as pd

file = pd.read_csv('seno.csv', sep=",", names=['x','seno','spike'])

plt.xlabel("x")
plt.ylabel("sen(x)")
plt.title("Período de uma série temporal")

"""
spikes: contains the instant t in which the spike occurs
period: contains the period between two consecutive spikes
period list is calculated from spikes list
"""

spikes = []
period = []

for i in range (0,len(file['spike'])):
    if file['spike'][i] != 0:
        plt.scatter(file['x'][i-1], file['seno'][i-1], c="k")
        spikes.append(file['x'][i-1])   

plt.plot(file['x'] ,file['seno'], "grey")
#plt.legend(loc="upper right")
plt.savefig("pico.png")

for i in range (1,len(spikes)):
    period.append(spikes[i] - spikes[i-1])

print(spikes)
print(period)