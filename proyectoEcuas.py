"""
Intento de usar RK4 para resolver el sistema de EDOS resultante del proyecto
de Ecuaciones Diferenciales 1 - Pendulos acoplados con un resorte
Segundo semestre de 2019
"""

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
m = 0.2;
l = 0.5;
g = 9.8;
k = 0.5;

def thetadot(t,x,y,theta,phi):
	return y

def ydotR(t,x,y,theta,phi):
	return -(g/l)*math.sin(theta) - (k/m)*cos(theta)*(math.sin(theta) - math.sin(phi))

def ydot(t,x,y,theta,phi):
	return theta*(-(g/l) - k/m) +phi*k/m

def phidot(t,x,y,theta,phi):
	return x

def xdotR(t,x,y,theta,phi):
	return -(g/l)*math.sin(phi) -(k/m)*math.sin(phi)*math.cos(phi) + (k/m)*math.sin(theta)*math.cos(theta)
def xdot(t,x,y,theta,phi):
	return phi*(-(g/l) - k/m) + theta*k/m

def RungeKutta(thetadot, ydot, phidot, xdot,tI, thetai, yi, phii, xi, tF, h):
	#Como es una aproximación numérica, el resultado se guarda en arreglos de datos tipo numérico.
	# El primer elemento de cada arreglo representa las condiciones iniciales.
	t,theta,y,phi,x = [tI] ,[thetai] ,[yi] ,[phii], [xi]
	i = 0
	while t[i] <= tF:
		#Se calculan los valores de K respectivos asociados a cada ecuación
		k1theta = thetadot	(t[i],x[i],y[i],theta[i],phi[i])
		k1y = ydot			(t[i],x[i],y[i],theta[i],phi[i])
		k1phi = phidot		(t[i],x[i],y[i],theta[i],phi[i])
		k1x = xdot			(t[i],x[i],y[i],theta[i],phi[i])

		k2theta = thetadot	(t[i]+h/2,	x[i]+(h/2)*k1x	,y[i]+(h/2)*k1y	,theta[i]+(h/2)*k1theta	,phi[i]+(h/2)*k1phi)
		k2y = ydot			(t[i]+h/2,	x[i]+(h/2)*k1x	,y[i]+(h/2)*k1y	,theta[i]+(h/2)*k1theta	,phi[i]+(h/2)*k1phi)
		k2phi = phidot		(t[i]+h/2,	x[i]+(h/2)*k1x	,y[i]+(h/2)*k1y	,theta[i]+(h/2)*k1theta	,phi[i]+(h/2)*k1phi)
		k2x = xdot			(t[i]+h/2,	x[i]+(h/2)*k1x	,y[i]+(h/2)*k1y	,theta[i]+(h/2)*k1theta	,phi[i]+(h/2)*k1phi)

		k3theta = thetadot	(t[i]+h/2,	x[i]+(h/2)*k2x	,y[i]+(h/2)*k2y	,theta[i]+(h/2)*k2theta	,phi[i]+(h/2)*k2phi)
		k3y = ydot			(t[i]+h/2,	x[i]+(h/2)*k2x	,y[i]+(h/2)*k2y	,theta[i]+(h/2)*k2theta	,phi[i]+(h/2)*k2phi)
		k3phi = phidot		(t[i]+h/2,	x[i]+(h/2)*k2x	,y[i]+(h/2)*k2y	,theta[i]+(h/2)*k2theta	,phi[i]+(h/2)*k2phi)
		k3x = xdot			(t[i]+h/2,	x[i]+(h/2)*k2x	,y[i]+(h/2)*k2y	,theta[i]+(h/2)*k2theta	,phi[i]+(h/2)*k2phi)

		k4theta = thetadot	(t[i]+h/2,	x[i]+(h/2)*k3x	,y[i]+(h/2)*k3y	,theta[i]+(h/2)*k3theta	,phi[i]+(h/2)*k3phi)
		k4y = ydot			(t[i]+h/2,	x[i]+(h/2)*k3x	,y[i]+(h/2)*k3y	,theta[i]+(h/2)*k3theta	,phi[i]+(h/2)*k3phi)
		k4phi = phidot		(t[i]+h/2,	x[i]+(h/2)*k3x	,y[i]+(h/2)*k3y	,theta[i]+(h/2)*k3theta	,phi[i]+(h/2)*k3phi)
		k4x = xdot			(t[i]+h/2,	x[i]+(h/2)*k3x	,y[i]+(h/2)*k3y	,theta[i]+(h/2)*k3theta	,phi[i]+(h/2)*k3phi)

		#Los valores calculados se agregan
		x.append(x[i] + h*(k1x + 2*k2x + 2*k3x + k4x)/6)
		y.append(y[i] + h*(k1y + 2*k2y + 2*k3y + k4y)/6)
		theta.append(theta[i] + h*(k1theta + 2*k2theta + 2*k3theta + k4theta)/6)
		phi.append(phi[i] + h*(k1phi + 2*k2phi + 2*k3phi + k4phi)/6)
		t.append(t[i] + h)
		i+= 1
	#Solo retornamos theta y phi porque x,y son variables mudas
	return t,theta,phi
"""
t,theta,phi = RungeKutta(thetadot, ydot, phidot,xdot,0,3.14/2,0,0,0,4*3.14,0.0001)

tst = []
for i in range(len(t)):
    tst.append(theta[i]+phi[i])

t,theta,phi = RungeKutta(thetadot, ydotR, phidot,xdotR,0,3.14/2,0,0,0,4*3.14,0.0001)

#Se muestran resultados
plt. plot(t, theta,'r:',label='$\\theta(t)$')
plt.plot(t, phi,'b--',label='$\\phi(t)$')
plt.ylabel('Posicion angular (rad)')
plt.xlabel('Tiempo (s)')
plt.title('Posicion de los pendulos, $\\theta(0) = \pi / 2$')
plt.legend(loc='upper right')
plt.show()


"""
#plt.plot(t, tst,'g--')
"""
plt. plot(t, theta,'r--',label='Posición de la masa 1')
plt.plot(t, phi,'b',label='Posición de la masa 2')
plt.ylabel('Posicion angular (rad)')
plt.xlabel('Tiempo (s)')
plt.title('Posicion de los pendulos')
plt.legend(loc='upper right')
#plt.ylim(-0.8,0.8)
plt.show()
"""
#Se hace una llamada al metodo y se almacenan los resultados en los arreglos t, theta, y phi
t,theta,phi = RungeKutta(thetadot, ydot, phidot,xdot,0,3.14/6,0,0,0,7*3.14,0.0001)
tr,thetar,phir = RungeKutta(thetadot, ydotR, phidot,xdotR,0,3.14/6,0,0,0,7*3.14,0.0001)

theta0 = []
for i in range(len(t)):
	theta0.append(thetar[i] - theta[i])

suma = 0
for j in range(len(theta0)):
	suma +=theta0[j]

resultado = suma/len(theta0)
print(resultado)
"""tst = []
for i in range(len(t)):
    tst.append(theta[i]+phi[i])

#Se muestran resultados
plt. plot(t, theta,'c',label='$\\theta(t)$')
plt.plot(t, phi,'b',label='$\\phi(t)$')
#plt.plot(t, tst,'g--')
plt.ylabel('Posición angular (rad)')
plt.xlabel('Tiempo (s)')
plt.title('Posicion de los pendulos')
plt.legend(loc='upper right')
plt.ylim(-.8,0.85)
plt.show()"""