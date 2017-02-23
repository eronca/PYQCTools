import numpy as np
import matplotlib.pyplot as plt

def parabola(coeff, x):
    ypar = coeff[0]*pow(x,2.0)+coeff[1]*x+coeff[2]
    return ypar 

#Read values of the 3 points to fit from file
data = np.loadtxt('max_data.txt')

amat = np.zeros((3,3))

amat[:,0] = pow(data[:,0],2.0)
amat[:,1] = data[:,0]
amat[:,2] = np.ones(3)

bvec = data[:,1]

#Solve an equations system to calculate the parameter of the parabola
coeff = np.linalg.solve(amat, bvec)

#Check whether the solution is correct
assert(np.allclose(np.dot(amat, coeff), bvec))

print "Coefficients of the Parabola : ", coeff[0], coeff[1], coeff[2]

#Calculate the vertex of the parabola
vert_coords = np.zeros(2)
Delta = pow(coeff[1],2.0)-4.0*coeff[0]*coeff[2]
vert_coords[0] = -coeff[1]/(2.0*coeff[0])
vert_coords[1] = -Delta/(4.0*coeff[0])

print "Vertex Coordinates : ", vert_coords
print "Peak position (eV) : ", vert_coords[0]*27.2113

x = np.arange(data.T[0][2], data.T[0][0], 0.0001)
y = parabola(coeff,x)

plt.plot(x, y, 'r', linewidth=2.5, label = '')
plt.plot(data[:,0], data[:,1], 'ko', label = '')
plt.plot(vert_coords[0], vert_coords[1], color = 'b', marker = 'd', label = '')
plt.show()
plt.close()
