#!/usr/bin/env python
#
# Author: Enrico Ronca <enrico.r8729@gmail.com>
#                                                                                                                                    \


import os
from math import sqrt
import numpy

def run(GFfile, omega):
 
  #int1

  GF = numpy.loadtxt(GFfile)
  GFtrace = numpy.trace(GF)
  print omega, GFtrace

if __name__=="__main__":
    import sys
    run(sys.argv[1], sys.argv[2])
