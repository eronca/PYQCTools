#!/usr/bin/env python
#
# Author: Enrico Ronca <enrico.r8729@gmail.com>
#                                                                                                                                    \


import os
from math import sqrt
import numpy

def write_coupling(U,f,nsites):
    for i in range(1,nsites+1):
      print >>f, '{0:.12e}'.format(U), i,i,i,i

def write_hopping(t,f,nsites):
    print >>f, '{0:.12e}'.format(t), 1,2,0,0
    for i in range(2,nsites):
      print >>f, '{0:.12e}'.format(t), i,i-1,0,0
      print >>f, '{0:.12e}'.format(t), i,i+1,0,0
    print >>f, '{0:.12e}'.format(t), nsites,nsites-1,0,0

def run(nsites, t, U, fname):

    nelec = nsites

    f = open(fname,'w')
    print >>f, '&FCI NORB=',nsites,', NELEC= ', nelec, ' ,MS2= ',0,','
    f.write('ORBSYM=')
    for i in range(nsites):
        f.write("%d," % 1)
    f.write('\n')
    print >>f, '&END'
    #dumping coupling integrals
    write_coupling(U,f,nsites)
    #dumping of hopping integrals
    write_hopping(t,f,nsites)
    print >> f, '{0:.12e}'.format(0),0,0,0,0
    f.close()


if __name__=="__main__":
    import sys
    run(int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4])
