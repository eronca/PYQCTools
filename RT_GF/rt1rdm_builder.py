#!/usr/bin/python     

import os, sys, math
import numpy as np

def input_name(nt):
    inputname = 'green.'+str(nt)+'.'+str(nt)+'.txt'
    return inputname

def output_name(nt,component,delta_t):
    outputname = component+'_1rdm_'+str(nt*delta_t)+'.txt'
    return outputname

def run(tot_t, delta_t, scratch):

  #Define variables
  scratch = scratch+'/node0/'
  tot_t = float(tot_t)
  delta_t = float(delta_t)
  nsteps = int(tot_t/delta_t) 

  #Calculate and save 1RDM at t = 0
  nt = 0
  Re_1rdm = np.loadtxt(scratch+input_name(nt))
  Im_1rdm = np.zeros((Re_1rdm.shape))

  np.savetxt(output_name(nt,'Re',delta_t),Re_1rdm,fmt='%10.5f')
  np.savetxt(output_name(nt,'Im',delta_t),Im_1rdm,fmt='%10.5f')

  #Calculate and save 1RDM at the following times
  for nt in range(1,nsteps+1):
    ReRe_1rdm = np.loadtxt(scratch+input_name(nt))
    ImRe_1rdm = np.loadtxt(scratch+input_name(nt+100000))
    ImIm_1rdm = np.loadtxt(scratch+input_name(nt+200000))

    Re_1rdm = ReRe_1rdm - ImIm_1rdm
    Im_1rdm = -2*ImRe_1rdm

    np.savetxt(output_name(nt,'Re',delta_t),Re_1rdm,fmt='%10.5f')
    np.savetxt(output_name(nt,'Im',delta_t),Im_1rdm,fmt='%10.5f')


if __name__=="__main__":
  import sys, math
  
  tot_t = sys.argv[1]
  delta_t = sys.argv[2]
  scratch = sys.argv[3]

  run(tot_t, delta_t, scratch)
