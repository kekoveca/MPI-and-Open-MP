import numpy as np
import matplotlib.pyplot as plt

times = np.loadtxt('times.txt')
x = np.arange(1, 13)
times = times[0] / times

plt.figure(figsize=(10,6),dpi=300)
plt.scatter(x, times, c = 'b')
plt.plot(x, times, c = 'b', linestyle = '--', alpha = 0.3)
plt.grid(visible=True, which='major', color='black', linestyle='-', alpha = 0.3)
plt.grid(visible=True, which='minor', color='black', linestyle='--', alpha = 0.3)
plt.xlabel('N, кол-во исполнителей')
plt.ylabel('$T_1/T_N$') 
plt.title('График ускорения параллельного алгоритма\n6 ядер, 12 потоков')
plt.minorticks_on()
plt.savefig('life_accel.png')