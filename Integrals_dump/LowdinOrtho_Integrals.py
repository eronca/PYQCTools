#!/usr/bin/env python
#
# Author: Enrico Ronca <enrico.r8729@gmail.com>
#

from pyscf import gto,ao2mo,fci,lo,scf
import scipy as sp
import numpy as np
import pyscf.tools.fcidump as fcid

if __name__ == '__main__':

  #Build molecule
  mol = gto.Mole()
  mol.verbose = 0

  mol.unit = 'bohr'
  mol.atom = [
   ['H', ( 0., 0.    ,  i*1.80000  )] for i in range(10)
  ]

  #mol.symmetry = True
  mol.basis = {'H': 'sto-6g'}
  mol.build()

  m = scf.RHF(mol)
  m.verbose = 4
  ehf = m.scf()
  print "ORB ENERGY : ", m.mo_energy

  #1-e Integrals
  hcore = mol.intor('cint1e_nuc_sph') + mol.intor('cint1e_kin_sph')
  ovlp = mol.intor('cint1e_ovlp_sph')

  Sihalf = lo.orth.lowdin(ovlp)
  h1e = np.einsum('pi,pq,qj->ij',Sihalf,hcore,Sihalf)
  nmo = np.shape(h1e)[0]

  #2-e Integrals
  eri_4fold = ao2mo.kernel(mol,Sihalf)

  #FCI-DUMP
  fname = 'FCIDUMP'
  nelec = mol.nelectron
  energy_nuc = mol.energy_nuc()
  fcid.from_integrals(fname,h1e,eri_4fold,nmo,nelec,energy_nuc)

  #Energy Test
  e = fci.direct_spin0.kernel(h1e,eri_4fold,nmo,nelec)
  print e[0]+energy_nuc
