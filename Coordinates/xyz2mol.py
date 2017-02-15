import sys, os
import numpy as np

with open('mol', 'w') as fout:
  with open('inp.xyz', 'r') as fin:
    n = int(fin.readline())
    fin.readline()
    for i in range(n):
        line = fin.readline()
        string = '[\''+line.split()[0]+'\',  ('+line.split()[1]+', '+line.split()[2]+', '+line.split()[3]+')],\n'
        fout.write(string)
