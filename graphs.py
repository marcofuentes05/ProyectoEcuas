import math
import numpy as np
import matplotlib.pyplot as plt

m = 0.2;
l = 0.3;
g = 9.8;
k = 0.9;

def theta_1(t,c_1,c_2):
	return c_1*math.sin(math.sqrt(g/l)*t)+c_2*math.cos(math.sqrt(g/l)*t)

def theta_2(t,c_1,c_2):
	return c_1*math.cos(math.sqrt((g/l)+(2*k/(m)))*t + c_2)
t = np.arange(0,5*3.14,0.01)
t_1 = []
t_2 = []
t_3 = []
for i in t: t_1.append(theta_1(i,0.7,0.7))

for j in t: t_2.append(theta_2(j,1,0))

for k in range(len(t)): t_3.append(t_1[k]+t_2[k])
plt.title('Modos normales de oscilación del sistema')
plt.plot(t,t_1,'b:',label='Modo 1')
plt.plot(t,t_2,'r:',label='Modo 2')
plt.plot(t,t_3,'darkcyan',label='Combinacion lineal de modos normales')
plt.ylabel('Posicion angular (rad)')
plt.xlabel('Tiempo (s)')
plt.legend(loc='upper right')
plt.ylim(-3,3.5)
plt.show()