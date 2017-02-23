import numpy as np
import matplotlib.pyplot as plt

#Read data for the extrapolation from file
data = np.loadtxt('extra')

x = data[:,0]
y = data[:,1]

#Perform polynomial fit of the data (linear in this case)
coeff = np.polyfit(x, y, 1)
polynomial = np.poly1d(coeff)

print "Fitted function :", polynomial

#Calculates R squared value
ybar = np.sum(y)/len(y)
ys = polynomial(x)
ssreg = np.sum((ys-ybar)**2) 
ssres = np.sum((y - ys)**2)
sstot = np.sum((y - ybar)**2)
R2 = 1.0 - (ssres/sstot)

print "R squared :", R2

#Calculates slope and intercept of the line (for the linear case)
print "Slope :", coeff[0], "Intercept :", coeff[1]

plt.plot(x, y, 'o')
plt.plot(x, ys)
plt.ylabel('y')
plt.xlabel('x')
plt.show()


