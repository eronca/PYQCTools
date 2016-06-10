#!/usr/bin/env python
#
# Author: Enrico Ronca <enrico.r8729@gmail.com>
#

import numpy

IRREP_MAP = {'D2h': (1,         # Ag
                     4,         # B1g
                     6,         # B2g
                     7,         # B3g
                     8,         # Au
                     5,         # B1u
                     3,         # B2u
                     2),        # B3u
             'C2v': (1,         # A1
                     4,         # A2
                     2,         # B1
                     3),        # B2
             'C2h': (1,         # Ag
                     4,         # Bg
                     2,         # Au
                     3),        # Bu
             'D2' : (1,         # A
                     4,         # B1
                     3,         # B2
                     2),        # B3
             'Cs' : (1,         # A'
                     2),        # A"
             'C2' : (1,         # A
                     2),        # B
             'Ci' : (1,         # Ag
                     2),        # Au
             'C1' : (1,)}

def write_dip_sym(dip,f,reorder,tol):
    for i in xrange(0,dip.shape[0]):
        for j in xrange(0,i+1):
            if (abs(dip[i,j]) > tol):
                print >>f, '{0:.12e}'.format(dip[i,j]), reorder[0][i]+1,reorder[1][j]+1,0,0


def fcidump(mc, orbsym):
    from pyscf import fci
    from pyscf import mcscf

    mo_cas = mc.mo_coeff[:,mc.ncore:mc.ncore+mc.ncas]

    AOdipInt = mol.intor('cint1e_r_sph', comp=3)

    MOdipInt = numpy.zeros([3,mo_cas.shape[1],mo_cas.shape[1]])

    for i in range(AOdipInt.shape[0]):
        MOdipInt[i] = reduce(numpy.dot, (mo_cas.T, AOdipInt[i], mo_cas))

    nelec = mc.nelecas[0] + mc.nelecas[1]

    ncore = mc.ncore
    ncas = mc.ncas
    sz = mc.ncas
    tol = float(1e-15)
    reorder = range(mc.mo_coeff.shape[0])

    reorder_orbsym = orbsym[:]
    reorder_orbsym = [0]*len(orbsym)

    for i in xrange(mc.ncas):
        reorder[i] = i
        reorder_orbsym[i] = orbsym[i]

    cas_order = reorder[:mc.ncas]

    f = open('FCIDUMP_DIP_X','w')
    print >>f, '&FCI NORB=',ncas,', NELEC= ', nelec, ' ,MS2= ',0,',' 
    f.write('ORBSYM=')
    for i in range(ncas):
        f.write("%d," % reorder_orbsym[i])
    f.write('\n')
    print >>f, '&END'
    write_dip_sym(MOdipInt[0],f,(cas_order,cas_order),tol)
    f.close()

    f = open('FCIDUMP_DIP_Y','w')
    print >>f, '&FCI NORB=',ncas,', NELEC= ', nelec, ' ,MS2= ',0,','
    f.write('ORBSYM=')
    for i in range(ncas):
        f.write("%d," % reorder_orbsym[i])
    f.write('\n')
    print >>f, '&END'
    write_dip_sym(MOdipInt[1],f,(cas_order,cas_order),tol)
    f.close()

    f = open('FCIDUMP_DIP_Z','w')
    print >>f, '&FCI NORB=',ncas,', NELEC= ', nelec, ' ,MS2= ',0,','
    f.write('ORBSYM=')
    for i in range(ncas):
        f.write("%d," % reorder_orbsym[i])
    f.write('\n')
    print >>f, '&END'
    write_dip_sym(MOdipInt[2],f,(cas_order,cas_order),tol)
    f.close()










    

    

if __name__ == '__main__':
    from functools import reduce
    from pyscf import gto
    from pyscf import scf
    from pyscf import ao2mo
    from pyscf import fci
    from pyscf import mcscf
    from pyscf import symm

    mol = gto.Mole()
    mol.verbose = 0
    mol.atom = [
         ['O',    (0.000000,    0.000000,    1.210000)],
         ['C',    (0.000000,    0.000000,    0.000000)],
         ['H',    (0.000000,    0.943435,   -0.599878)],
         ['H',    (0.000000,   -0.943435,   -0.599878)],
     ]
    mol.basis = {'H': 'cc-pVTZ', 'C': 'cc-pVTZ', 'O': 'cc-pVTZ'}
    mol.symmetry = 'c1'
    mol.build()

    m = scf.RHF(mol)

    m.kernel()
    irrep_name = m.mol.irrep_id
    sym = symm.label_orb_symm(m.mol, irrep_name,m.mol.symm_orb,m.mo_coeff)
    orbsym = [IRREP_MAP['C1'][i % 10] for i in sym]

    m.verbose = 4
    ehf = m.scf()
    mc=mcscf.CASCI(m,16,16)
    mc.verbose = 4
#    ci_e = mc.kernel()[0]
#    print 'CI energy', ci_e
#    mc.casci()

    mo_cas_coeff = mc.mo_coeff[:,mc.ncore:mc.ncore+mc.ncas]
    AOdipInt = mol.intor('cint1e_r_sph', comp=3)

    MOdipInt = numpy.zeros([3,mo_cas_coeff.shape[1],mo_cas_coeff.shape[1]])
   
    for i in range(AOdipInt.shape[0]):
        MOdipInt[i] = reduce(numpy.dot, (mo_cas_coeff.T, AOdipInt[i], mo_cas_coeff)) 

    print MOdipInt

    fcidump(mc,orbsym)
