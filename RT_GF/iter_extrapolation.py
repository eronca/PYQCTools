#!/usr/bin/python                                                                                                                                             \


import os
import numpy as np
from scipy import fft, arange
from distutils.util import strtobool

def run_extrapolation(real_RTGF, imag_RTGF, npoints):
    time_array = real_RTGF.transpose()[0]
    delta_t = time_array[1]-time_array[0]

    nfit = npoints/2

    x_imag = imag_RTGF.transpose()[1][-npoints:]
    x_real = real_RTGF.transpose()[1][-npoints:] 

    x = x_real + 1j*x_imag

    R = np.zeros((nfit,nfit), dtype = np.complex) 
    r = np.zeros((nfit), dtype = np.complex) 
    for j in range(1,nfit+1):
     for i in range(1,nfit+1):
      r[j-1] = 0
      for n in range(nfit,npoints):
        r[j-1] += np.conjugate(x[n-j])*x[n]
        R[j-1][i-1] += np.conjugate(x[n-j])*x[n-i]

    invR = np.linalg.inv(R)
    a = - np.dot(invR,r)

    A = np.eye(nfit,k=-1,dtype = np.complex)
    A[0] = -a

    eig, R = np.linalg.eig(A)
    eig1, L = np.linalg.eig(np.transpose(A))
    X = x[nfit:]
    X = X[::-1]
    #print np.dot(R,np.dot(np.diag(eig),np.dot(np.linalg.inv(R),X)))
    for i in range(A.shape[0]):
      if (np.abs(eig[i]) > 1.0):
          eig[i] /= np.abs(eig[i])
          #eig[i] = 1/np.conjugate(eig[i])

    A = np.dot(R,np.dot(np.diag(eig),np.linalg.inv(R)))
    #print np.dot(R,np.dot(np.diag(eig),np.dot(np.linalg.inv(R),X)))

    new_vals = np.zeros((4), dtype = np.complex)
    X = x[nfit:]
    X = X[::-1]
    for n_new in range(4):
      X = np.dot(A,X)
      new_vals[n_new]  = X[0]

    xfull = real_RTGF.transpose()[1] + 1j*imag_RTGF.transpose()[1] 
    expanded = np.concatenate([xfull,new_vals])
    time = time_array[0]

    with open('new_full_real.out', 'w') as fout:
       fout.write('#     Time          A(t)\n')
       for i in range(xfull.shape[0]+4):
           fout.write('%12.8f  %12.8f\n' % (time, expanded[i].real))
           time += delta_t

    time = time_array[0]
    with open('new_full_imag.out', 'w') as fout:
       fout.write('#     Time          A(t)\n')
       for i in range(xfull.shape[0]+4):
           fout.write('%12.8f  %12.8f\n' % (time, expanded[i].imag))
           time += delta_t


if __name__=="__main__":
    import sys, math

    tot_points = int(sys.argv[1])
   
    real_RTGF = np.loadtxt("rt_real.txt") 
    imag_RTGF = np.loadtxt("rt_imag.txt")
    npoints = real_RTGF.transpose()[0].shape[0]

    init_npoints = npoints
   
    run_extrapolation(real_RTGF, imag_RTGF, init_npoints)

    while(npoints+4 < tot_points):
        real_RTGF = np.loadtxt("new_full_real.out")
        imag_RTGF = np.loadtxt("new_full_imag.out")
        npoints = real_RTGF.transpose()[0].shape[0]

        run_extrapolation(real_RTGF, imag_RTGF, init_npoints)
