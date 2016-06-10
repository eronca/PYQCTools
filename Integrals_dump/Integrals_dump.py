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

def writeh2e_sym(h2e,f,reorder,tol):
    for i in xrange(0,h2e.shape[0]):
        for j in xrange(0,i+1):
            for k in xrange(0,h2e.shape[2]):
                for l in xrange(0,k+1):
                    if (abs(h2e[i,j,k,l]) > tol and i*h2e.shape[0]+j >= k*h2e.shape[2]+l ):
                        print >>f, '{0:.12e}'.format(h2e[i,j,k,l]), reorder[0][i]+1,reorder[1][j]+1,reorder[2][k]+1,reorder[3][l]+1

def writeh1e_sym(h1e,f,reorder,tol):
    for i in xrange(0,h1e.shape[0]):
        for j in xrange(0,i+1):
            if (abs(h1e[i,j]) > tol):
                print >>f, '{0:.12e}'.format(h1e[i,j]), reorder[0][i]+1,reorder[1][j]+1,0,0


def fcidump(mc, orbsym):
    from pyscf import fci
    from pyscf import mcscf

    mo_core = mc.mo_coeff[:,:mc.ncore]
    mo_cas = mc.mo_coeff[:,mc.ncore:mc.ncore+mc.ncas]
    mo_virt = mc.mo_coeff[:,mc.ncore+mc.ncas:]
    h2e = ao2mo.incore.general(mc._scf._eri,[mo_cas,mo_cas,mo_cas,mo_cas],compact=False)
    h2e = h2e.reshape(mc.ncas,mc.ncas,mc.ncas,mc.ncas)
    h1e,energy_core = mc.h1e_for_cas()
    energy_core += mc.mol.energy_nuc()

    nelec = mc.nelecas[0] + mc.nelecas[1]

    core_dm = numpy.dot(mo_core,mo_core.T) *2
    core_vhf = mc.get_veff(mc.mol,core_dm)

    ncore = mc.ncore
    ncas = mc.ncas
    nvirt = mo_virt.shape[1]
    sz = mc.ncas
    tol = float(1e-15)
    reorder = range(mc.mo_coeff.shape[0])

    reorder_orbsym = orbsym[:]
    reorder_orbsym = [0]*len(orbsym)

    for i in xrange(mc.ncas):
        reorder[i] = i
        reorder_orbsym[i] = orbsym[i]

    cas_order = reorder[:mc.ncas]

    f = open('FCIDUMP','w')
    print >>f, '&FCI NORB=',ncas,', NELEC= ', nelec, ' ,MS2= ',0,',' 
    f.write('ORBSYM=')
    for i in range(ncas):
        f.write("%d," % reorder_orbsym[i])
    f.write('\n')
    print >>f, '&END'
    #h2e in active space
    writeh2e_sym(h2e,f,(cas_order,cas_order,cas_order,cas_order),tol)
    #h1e in active space
    writeh1e_sym(h1e,f,(cas_order,cas_order),tol)
    print >> f, energy_core,0,0,0,0
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
    mol.output = 'INT_new.out'
    mol.atom = [
        ['N', ( 0., 0.    , -1.090105133/2    )],
        ['N', ( 0., 0.    ,  1.090105133/2   )],
    ]
    mol.symmetry = True
    mol.basis = {'N': "cc-pVDZ"}
    mol.build()

    m = scf.RHF(mol)

    m.kernel()
    irrep_name = m.mol.irrep_id
    sym = symm.label_orb_symm(m.mol, irrep_name,m.mol.symm_orb,m.mo_coeff)
    orbsym = [IRREP_MAP['D2h'][i % 10] for i in sym]

    print "Orbital Energies"
    print m.mo_energy

    m.verbose = 4
    ehf = m.scf()
    mc=mcscf.CASCI(m,28,14)
    mc.verbose = 4
    #ci_e = mc.kernel()[0]
    #print 'CI energy', ci_e
    #mc.casci()

    fcidump(mc,orbsym)


