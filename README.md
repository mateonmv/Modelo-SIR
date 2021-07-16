# Modelo-SIR

#EulerExplicito
import numpy as np
import sympy as sp
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button
from sympy import Symbol, solve, Eq
from math import *
import matplotlib.pyplot as plt

def fi(s, y):    #Población Infectada
    return b*(s)*y-k*y

i, f = 0, 180  # intervalo de integración
N = 200
h = (f-i)/N
b = 0.5
k = 0.333
s_0 = 1 #Poblacion total=7900000
i_0 = 1.27*(10**-6) #Infectados dia1

ts = np.linspace(i, f, N+1)
s = [s_0]
i = [i_0]
r = [0]

for t in ts[:-1]:
    i.append(i[-1]+h*fi(s[-1], i[-1]))
    s.append(s[-1]-(0.5*i[-1])*h*s[-1]) 
    r.append(r[-1]+i[-1]*k*h)

plt.plot(ts,s,color="red",label="Susceptible")
plt.plot(ts,i,color="b",label="Infectado")
plt.plot(ts,r,color="y",label="Recuperado")
plt.legend()
plt.title('EULER EXPLICITO')
plt.grid()
plt.xlabel('Dias')
plt.ylabel('Poblacion')
plt.show()

#EulerImplicito

import numpy as np
import sympy as sp
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button
from sympy import Symbol, solve, Eq
from math import *
import matplotlib.pyplot as plt

def fi(s, y):    #Población Infectada
    return b*(s)*y-k*y

i, f = 0, 180  # intervalo de integración
N = 200
h = (f-i)/N
b = 0.5
k = 0.333
s_0 = 1 #Poblacion total=7900000
i_0 = 1.27*(10**-6) #Infectados dia1

ts = np.linspace(i, f, N+1)
s = [s_0]
i = [i_0]
r = [0]
u_n = Symbol('u_n')

for t in ts[1:]:
    res = solve(Eq(i[-1] + h*fi(s[-1], u_n), u_n), u_n)
    i.append(res[0])
    s.append(s[-1]-(0.5*i[-1])*h*s[-1]) 
    r.append(r[-1]+0.333*i[-1]*h)
    

plt.plot(ts,s,color="red",label="Susceptible")
plt.plot(ts,i,color="b",label="Infectado")
plt.plot(ts,r,color="y",label="Recuperado")
plt.legend()
plt.title('EULER IMPLICITO')
plt.grid()
plt.xlabel('Dias')
plt.ylabel('Poblacion')
plt.show()

#MetodoTrapecios

import numpy as np
import sympy as sp
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button
from sympy import Symbol, solve, Eq
from math import *
import matplotlib.pyplot as plt
def fi(s, y):    #Población Infectada
    return b*(s)*y-k*y

i, f = 0, 180  # intervalo de integración
N = 200
h = (f-i)/N
b = 0.5
k = 0.333
s_0 = 1 #Poblacion total=7900000
i_0 = 1.27*(10**-6) #Infectados dia1

ts = np.linspace(i, f, N+1)
s = [s_0]
u_n = Symbol('u_n')
i = [i_0]
r = [0]
for t in ts[1:]:
    res = solve(Eq(i[-1] + (h/2)*(fi(s[-1], i[-1])+fi(s[-1], u_n)), u_n), u_n)
    i.append(res[0])
    s.append(s[-1]-0.5*i[-1]*s[-1]*h) 
    r.append(r[-1]+0.3333*i[-1]*h)

plt.plot(ts,s,color="red",label="Susceptible")
plt.plot(ts,i,color="b",label="Infectado")
plt.plot(ts,r,color="y",label="Recuperado")
plt.legend()
plt.title('Metodo Trapecios')
plt.grid()
plt.xlabel('Dias')
plt.ylabel('Poblacion')
plt.show()

#Heun

import numpy as np
from matplotlib import pyplot as plt
from sympy import Symbol, Eq, solve
from sympy.solvers.solvers import nsolve

def fi(s, y):    #Población Infectada
    return b*s*y-k*y

i, f = 0, 365 # intervalo de integración
N = 400
h = 0.2
b = 0.5
k = 0.333
s_0 = 1 #Poblacion total=7900000
i_0 = 1.27*(10**-6) #Infectados dia1

ts = np.linspace(i, f, N+1)
s = [s_0]
i = [i_0]
r = [0]

for t in ts[:-1]:
    k1=fi(s[-1],i[-1])
    k2=fi(s[-1]+(2/3)*h,i[-1]+(2/3)*k1*h)
    i_n=(i[-1]+(h/4)*(k1+3*k2))
    i.append(i_n+(0.5*s[-1]*i_n-0.333*i_n)*h)
    s.append(s[-1]-(0.5*i[-1])*h*s[-1]) 
    r.append(r[-1]+i[-1]*k*h)


plt.plot(ts,s,color="red",label="Susceptible")
plt.plot(ts,i,color="b",label="Infectado")
plt.plot(ts,r,color="y",label="Recuperado")
plt.legend()
plt.title('RK HEUN')
plt.grid()
plt.xlabel('DIAS')
plt.ylabel('POBLACION')
plt.show()



#EulerModificado

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button
from sympy import Symbol, Eq, solve

def fi(s, y):    #Población Infectada
    return b*s*y-k*y

i, f = 0, 365  # intervalo de integración
N = 400
h = 0.2
b = 0.5
k = 0.333
s_0 = 1 #n=7900000
i_0 = 1.27*(10**-6) #10 infectados inicialmente

ts = np.linspace(i, f, N+1)
s = [s_0]
i = [i_0]
r = [0]
u_n=Symbol('u_n')

for t in ts[1:]:
    k1=fi(s[-1],i[-1])
    k2=fi(s[-1]+h,i[-1]+k1*h)
    i_n=i[-1]+(h/2)*(k1+k2)
    i.append(i_n+(0.5*s[-1]*i_n-0.333*i_n)*h)
    s.append(s[-1]-0.5*i[-1]*s[-1]*h) 
    r.append(r[-1]+0.3333*i[-1]*h)
    

plt.plot(ts,s,color="red",label="Susceptible")
plt.plot(ts,i,color="b",label="Infectado")
plt.plot(ts,r,color="y",label="Recuperado")
plt.legend()
plt.title('Euler Modificado')
plt.grid()
plt.xlabel('Dias')
plt.ylabel('Poblacion')
plt.show()

#PuntoMedio

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button
from sympy import Symbol, Eq, solve

def fi(s, y):    #Población Infectada
    return b*s*y-k*y

i, f = 0, 365  # intervalo de integración
N = 400
h = 0.2
b = 0.5
k = 0.333
s_0 = 1 #n=7900000
i_0 = 1.27*(10**-6) #10 infectados inicialmente

ts = np.linspace(i, f, N+1)
s = [s_0]
i = [i_0]
r = [0]
u_n=Symbol('u_n')

for t in ts[1:]:
    k1=fi(s[-1],i[-1])
    k2=fi(s[-1]+0.5*h,i[-1]+0.5*k1*h)
    i_n= i[-1]+k2*h
    i.append(i_n+(0.5*s[-1]*i_n-0.333*i_n)*h)
    s.append(s[-1]-0.5*i[-1]*s[-1]*h) 
    r.append(r[-1]+0.333*i[-1]*h)
    

plt.plot(ts,s,color="red",label="Susceptible")
plt.plot(ts,i,color="b",label="Infectado")
plt.plot(ts,r,color="y",label="Recuperado")
plt.legend()
plt.title('Punto Medio')
plt.grid()
plt.xlabel('Dias')
plt.ylabel('Poblacion')
plt.show()

#RK 4toOrden

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button
from sympy import Symbol, Eq, solve

def fi(s, y):    #Población Infectada
    return b*(s)*y-k*y

i, f = 0, 365  # intervalo de integración
N = 400
h = 0.2
b = 0.5
k = 0.333
s_0 = 1 #Poblacion total=7900000
i_0 = 1.27*(10**-6) #Infectados dia1

ts = np.linspace(i, f, N+1)
s = [s_0]
i = [i_0]
r = [0]
u_n=Symbol('u_n')

for t in ts[1:]:
    k1=fi(s[-1], i[-1])
    k2=fi(s[-1]+(h/2),i[-1]+(h/2)*k1)
    k3=fi(s[-1]+(h/2),i[-1]+(h/2)*k2)
    k4=fi(s[-1]+h, i[-1]+k3*h)
    i_n=i[-1] + (h/6)*(k1+2*k2+2*k3+k4)
    i.append(i_n+(0.5*s[-1]*i_n-0.333*i_n)*h)
    s.append(s[-1]-(0.5*i[-1])*h*s[-1]) 
    r.append(r[-1]+0.333*i[-1]*h)


plt.plot(ts,s,color="red",label="Susceptible")
plt.plot(ts,i,color="b",label="Infectado")
plt.plot(ts,r,color="y",label="Recuperado")
plt.legend()
plt.title('RK 4TO ORDEN')
plt.grid()
plt.xlabel('Dias')
plt.ylabel('Poblacion')
plt.show()

