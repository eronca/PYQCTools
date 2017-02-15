import sys, os
import numpy as np

#VARIABLES DEFINITION
inpfile = 'in.xyz'
outfile = 'out.xyz'
#Define atom to set as the origin of axes
centatom = 0
#Atom to set along the x axis
atom1 = 1
#Atom to define the x-y plane of the molecule 
atom2 = 2 

#READ COORDINATES
with open(inpfile, 'r') as fin:
   natoms = int(fin.readline())
   fin.readline()
   names = []
   coords = np.zeros((natoms,3))
   for i in range(natoms):
       line = fin.readline()
       names.append(line.split()[0])
       coords[i][0] = line.split()[1]
       coords[i][1] = line.split()[2]
       coords[i][2] = line.split()[3]
fin.close()

#TRANSLATE COORDINATED TO THE ORIGIN
tcoords = np.zeros(coords.shape)
for i in range(natoms):
    tcoords[i] = coords[i] - coords[centatom]

#DEFINE UNIT VECTORS
ex = (tcoords[atom1] - tcoords[centatom])/np.linalg.norm(tcoords[atom1] - tcoords[centatom])
e2 = (tcoords[atom2] - tcoords[centatom])/np.linalg.norm(tcoords[atom2] - tcoords[centatom])
ez = np.cross(ex,e2)/np.linalg.norm(np.cross(ex,e2))
ey = np.cross(ex,ez)/np.linalg.norm(np.cross(ex,ez))

#CALCULATE THE ROTATION MATRIX
rot_mat = np.zeros((3,3))
rot_mat.T[0] = ex
rot_mat.T[1] = ey
rot_mat.T[2] = ez

#ORIENT COORDINATES
rcoords = np.dot(tcoords,rot_mat)

#WRITE NEW COORDINATES
with open(outfile, 'w') as fout:
   fout.write(str(natoms)+'\n')
   fout.write('\n')
   for i in range(natoms):
       fout.write(names[i]+'    '+str(rcoords[i][0])+'    '+str(rcoords[i][1])+'    '+str(rcoords[i][2])+'\n')
fout.close()
