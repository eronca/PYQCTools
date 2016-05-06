#!/usr/bin/python                                                                                                                                             \


import os
import numpy as np
from scipy import fft, arange
from distutils.util import strtobool

def run_extrapolation(real_RTGF, imag_RTGF, fullrange):
    time_array = real_RTGF.transpose()[0]
    npoints = time_array.shape[0]
    delta_t = time_array[1]-time_array[0]

    nfit = npoints/2

    x_imag = imag_RTGF.transpose()[1] 
    x_real = real_RTGF.transpose()[1] 

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

    new_vals = np.zeros((x.shape[0]), dtype = np.complex)
    X = x[nfit:]
    X = X[::-1]
    for n_new in range(x.shape[0]):
      X = np.dot(A,X)
      new_vals[n_new]  = X[0]

    if (fullrange == True):
       expanded = np.concatenate([x,new_vals])
       time = time_array[0]
       filename = 'new_full_data.out'
    else :
       expanded = new_vals
       time = time_array[npoints-1]+delta_t
       filename = 'predicted.out'

    with open(filename, 'w') as fout:
       fout.write('#     Time          A(t)\n')
       if (fullrange == True):
          for i in range(2*npoints):
              fout.write('%12.8f  %12.8f\n' % (time, expanded[i].real))
              time += delta_t
       else :
          for i in range(npoints):
              fout.write('%12.8f  %12.8f\n' % (time, expanded[i].real))
              time += delta_t


if __name__=="__main__":
    import sys, math

    fullrange = strtobool(sys.argv[1])
   
    real_RTGF = np.loadtxt("rt_real.txt") 
    imag_RTGF = np.loadtxt("rt_imag.txt")
   
    run_extrapolation(real_RTGF, imag_RTGF, fullrange)
