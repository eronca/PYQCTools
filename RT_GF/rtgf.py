#!/usr/bin/python                                                                                                                                             \


import os
from math import sqrt
import numpy

def run(GFfile, time, real):
  GF = numpy.loadtxt(GFfile)
  GFtrace = numpy.trace(GF)

  if (real) :
      with open('rt_real.txt', 'a') as fout:
          fout.write('%8.4f  %12.8f\n' % (time, GFtrace))
  else: 
      with open('rt_imag.txt', 'a') as fout:
          fout.write('%8.4f  %12.8f\n' % (time, GFtrace))


if __name__=="__main__":
    import sys, math

    prop_time = float(sys.argv[1])
    time_step = float(sys.argv[2])
    green_dir = sys.argv[3]

    steps_num = int(prop_time/time_step)

    with open('rt_real.txt', 'a') as fout:
         fout.write('#     Time          A(Time)\n')

    real_files = green_dir+'/green.%d.%d.txt' % (0,0)
    run(real_files, 0.0, True)

    with open('rt_imag.txt', 'a') as fout:
         fout.write('#     Time          A(Time)\n')
         fout.write('%8.4f  %12.8f\n' % (0.000, 0.000))

    time = time_step
    for itime in range(1, steps_num+1):
        real_files = green_dir+'/green.%d.%d.txt' % (itime,itime)
        imag_files = green_dir+'/green.%d.%d.txt' % (30000+itime,30000+itime)

        real = True
        run(real_files, time, real)
        real = False
        run(imag_files, time, real)

        time += time_step
