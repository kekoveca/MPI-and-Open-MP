import numpy as np
import matplotlib.pyplot as plt

times = np.loadtxt('./4-life/times.txt')
x = np.arange(1, 29)
times = times[0] / times

plt.figure(figsize=(10,6),dpi=300)
plt.scatter(x, times, c = 'b')
plt.plot(x, times, c = 'b', linestyle = '--', alpha = 0.3)
plt.grid(visible=True, which='major', color='black', linestyle='-', alpha = 0.3)
plt.grid(visible=True, which='minor', color='black', linestyle='--', alpha = 0.3)
plt.xlabel('N, кол-во исполнителей')
plt.ylabel('$T_1/T_N$') 
plt.title('График ускорения параллельного алгоритма')
plt.minorticks_on()
plt.savefig('./4-life/life_accel.png')