import matplotlib.pyplot as plt
import pandas as pd

file1 = pd.read_csv('1neuronI300.csv', sep=',', names = ['time','V','m','n','h'])
file2 = pd.read_csv('1neuronI280.csv', sep=',', names = ['time','V','m','n','h'])

limite = [0,2000]
lim_mnh = [0,2000]

plt.subplot(2,2,1)
plt.xlim(limite)
plt.xlabel('t(ms)')
plt.title('I = 280pA')
plt.plot(file1['V'], 'm', label = 'Cur = 280pA')
legend = plt.legend(loc = 'upper right')

plt.subplot(2,2,2)
plt.xlim(limite)
plt.xlabel('t(ms)')
plt.title('I = 300pA')
plt.plot(file2['V'], 'coral', label = 'Cur = 300pA')
legend = plt.legend(loc = 'upper right')

plt.subplot(2,2,3)
plt.xlim(lim_mnh)
plt.xlabel('t(ms)')
plt.plot(file1['m'], 'lime', label = 'm')
plt.plot(file1['n'], 'lightseagreen', label = 'n')
plt.plot(file1['h'], 'olive', label = 'h')
legend = plt.legend(loc = 'upper right')

plt.subplot(2,2,4)
plt.xlim(lim_mnh)
plt.xlabel('t(ms)')
plt.plot(file2['m'], 'r', label = 'm')
plt.plot(file2['n'], 'b', label = 'n')
plt.plot(file2['h'], 'yellow', label = 'h')
legend = plt.legend(loc = 'upper right')

plt.show()