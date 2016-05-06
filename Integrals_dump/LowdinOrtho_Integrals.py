#!/usr/bin/env python

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

if __name__ == '__main__':
    from functools import reduce
    from pyscf import gto
    from pyscf import scf
    from pyscf import ao2mo
    from pyscf import fci
    from pyscf import mcscf
    from pyscf import symm 
    from pyscf import lo 
    import scipy as sp
    import numpy as np
    import pyscf.tools.fcidump as fcid

    mol = gto.Mole()
    mol.verbose = 0
    mol.unit = 'au'
    mol.atom = [
        ['H', ( 0., 0.    ,  0.0   )],
        ['H', ( 0., 0.    ,  1.40112   )],
        ['H', ( 0., 0.    ,  3.90112   )],
        ['H', ( 0., 0.    ,  5.30224   )],
        ['H', ( 0., 0.    ,  7.80224   )],
        ['H', ( 0., 0.    ,  9.20336   )],
        ['H', ( 0., 0.    ,  11.70336   )],
        ['H', ( 0., 0.    ,  13.10448   )],
        ['H', ( 0., 0.    ,  15.60448   )],
        ['H', ( 0., 0.    ,  17.00560   )],
    ]
    mol.symmetry = True
    mol.basis = {'H': 'sto-6g'}
    mol.build()

    m = scf.RHF(mol)

    m.kernel()
    irrep_name = m.mol.irrep_id
    sym = symm.label_orb_symm(m.mol, irrep_name,m.mol.symm_orb,m.mo_coeff)
    orbsym = [IRREP_MAP['D2h'][i % 10] for i in sym]

    m.verbose = 4
    ehf = m.scf()
    print "ORB ENERGY : ", m.mo_energy

    #1-e Integrals
    hcore = mol.intor('cint1e_nuc_sph') + mol.intor('cint1e_kin_sph')
    ovlp = mol.intor('cint1e_ovlp_sph')

    Sihalf = lo.orth.vec_lowdin(m.mo_coeff,ovlp)
    h1e = np.einsum('pi,pq,qj->ij',Sihalf,hcore,Sihalf)
    energy_nuc = mol.energy_nuc()
    nmo = np.shape(h1e)[0]
    nelec = mol.nelectron

    #2-e Integrals
    eri_4fold = ao2mo.kernel(mol,Sihalf)

    #FCI-DUMP
    fcid.from_integrals('FCIDUMP',h1e,eri_4fold,nmo,nelec,energy_nuc,0,orbsym)
