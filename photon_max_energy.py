import numpy as np
import matplotlib.pyplot as plt 

h = 6.626e-34
gelen_foton = 1

E = []
X = []

for i in range (0,180):
    X.append(i)
    i = np.deg2rad(i)
    ayrilan_foton = gelen_foton / (1+(gelen_foton/ (0.511)) *(1- np.cos(i)))
    E.append(gelen_foton - ayrilan_foton)
    
    
plt.plot(E,X)
plt.show()


