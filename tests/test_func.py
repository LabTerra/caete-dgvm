import numpy as np
import matplotlib.pyplot as plt
import test1 as f

water = np.linspace(1,500,100)
out = np.zeros(100,)
i = 0

for w in water:
    out[i] = f.photo.water_stress_modifier(w, 1.1, 110, 0.9)
    i +=1 
plt.plot(water, out, 'r-'); plt.show()
