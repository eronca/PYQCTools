#!/usr/bin/env python
import numpy
from pyscf.mrpt.nevpt2 import sc_nevpt

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


def writeh2e(h2e,f,reorder,tol):
    for i in xrange(0,h2e.shape[0]):
        for j in xrange(0,h2e.shape[1]):
            for k in xrange(0,h2e.shape[2]):
                for l in xrange(0,h2e.shape[3]):
                    if (abs(h2e[i,j,k,l]) > tol):
                        #if ( j==k or j == l) :
                        #if (j==k and k==l) :
                            print >>f, ' % .12e  %3d  %3d  %3d  %3d' % (h2e[i,j,k,l], reorder[0][i]+1, reorder[1][j]+1, reorder[2][k]+1, reorder[3][l]+1)


def writeh1e(h1e,f,reorder,tol):
    for i in xrange(0,h1e.shape[0]):
        for j in xrange(0,h1e.shape[1]):
            if (abs(h1e[i,j]) > tol):
                print >>f, ' % .12e  %3d  %3d  %3d  %3d' % (h1e[i,j], reorder[0][i]+1, reorder[1][j]+1, 0, 0)

def writeh2e_sym(h2e,f,reorder,tol):
    for i in xrange(0,h2e.shape[0]):
        for j in xrange(0,i+1):
            for k in xrange(0,h2e.shape[2]):
                for l in xrange(0,k+1):
                    if (abs(h2e[i,j,k,l]) > tol and i*h2e.shape[0]+j >= k*h2e.shape[2]+l ):
                        print >>f, ' % .12e  %3d  %3d  %3d  %3d' % (h2e[i,j,k,l], reorder[0][i]+1, reorder[1][j]+1, reorder[2][k]+1, reorder[3][l]+1)

def writeh1e_sym(h1e,f,reorder,tol):
    for i in xrange(0,h1e.shape[0]):
        for j in xrange(0,i+1):
            if (abs(h1e[i,j]) > tol):
                print >>f, ' % .12e  %3d  %3d  %3d  %3d' % (h1e[i,j], reorder[0][i]+1, reorder[1][j]+1, 0, 0)


def dyall_perturb(mc, orbsym, nevpt_class = None):
    from pyscf import fci
    from pyscf import mcscf

    #Calculate 1 and 2-electron integrals for every nevpt2 class as well as the core energy
    mo_core = mc.mo_coeff[:,:mc.ncore]
    mo_cas = mc.mo_coeff[:,mc.ncore:mc.ncore+mc.ncas]
    mo_virt = mc.mo_coeff[:,mc.ncore+mc.ncas:]
    h2e = ao2mo.incore.general(mc._scf._eri,[mo_cas,mo_cas,mo_cas,mo_cas],compact=False)
    h2e = h2e.reshape(mc.ncas,mc.ncas,mc.ncas,mc.ncas)
    h1e,energy_core = mc.h1e_for_cas()
    energy_core += mc.mol.energy_nuc()

    #Calcolate two-electron integrals for the (0) class
    h2e_Sijrs = ao2mo.incore.general(mc._scf._eri,[mo_virt,mo_core,mo_virt,mo_core],compact=False)
    h2e_Sijrs = h2e_Sijrs.reshape(mo_virt.shape[1],mc.ncore,mo_virt.shape[1],mc.ncore)

    #Calcolate two-electron integrals for the (+1) class
    h2e_Sijr = ao2mo.incore.general(mc._scf._eri,[mo_virt,mo_core,mo_cas,mo_core],compact=False)
    h2e_Sijr = h2e_Sijr.reshape(mo_virt.shape[1],mc.ncore,mc.ncas,mc.ncore)

    #Calcolate two-electron integrals for the (-1) class
    h2e_Srsi = ao2mo.incore.general(mc._scf._eri,[mo_virt,mo_core,mo_virt,mo_cas],compact=False)
    h2e_Srsi = h2e_Srsi.reshape(mo_virt.shape[1],mc.ncore,mo_virt.shape[1],mc.ncas)

    #Calcolate two-electron integrals for the (+2) class
    h2e_Sij = ao2mo.incore.general(mc._scf._eri,[mo_cas,mo_core,mo_cas,mo_core],compact=False)
    h2e_Sij = h2e_Sij.reshape(mc.ncas,mc.ncore,mc.ncas,mc.ncore)

    #Calcolate two-electron integrals for the (-2) class
    h2e_Srs = ao2mo.incore.general(mc._scf._eri,[mo_virt,mo_cas,mo_virt,mo_cas],compact=False)
    h2e_Srs = h2e_Srs.reshape(mo_virt.shape[1],mc.ncas,mo_virt.shape[1],mc.ncas)

    #Calcolate two-electron integrals for the (0') class
    h2e_Sir1 = ao2mo.incore.general(mc._scf._eri,[mo_virt,mo_core,mo_cas,mo_cas],compact=False)
    h2e_Sir1 = h2e_Sir1.reshape(mo_virt.shape[1],mc.ncore,mc.ncas,mc.ncas)
    h2e_Sir2 = ao2mo.incore.general(mc._scf._eri,[mo_virt,mo_cas,mo_cas,mo_core],compact=False)
    h2e_Sir2 = h2e_Sir2.reshape(mo_virt.shape[1],mc.ncas,mc.ncas,mc.ncore)

    #Calcolate two-electron integrals for the (+1') class
    h2e_Si = ao2mo.incore.general(mc._scf._eri,[mo_cas,mo_core,mo_cas,mo_cas],compact=False)
    h2e_Si = h2e_Si.reshape(mc.ncas,mc.ncore,mc.ncas,mc.ncas)

    #Calcolate two-electron integrals for the (-1') class
    h2e_Sr = ao2mo.incore.general(mc._scf._eri,[mo_virt,mo_cas,mo_cas,mo_cas],compact=False)
    h2e_Sr = h2e_Sr.reshape(mo_virt.shape[1],mc.ncas,mc.ncas,mc.ncas)

    nelec = mc.nelecas[0] + mc.nelecas[1]
   
    core_dm = numpy.dot(mo_core,mo_core.T) *2
    core_vhf = mc.get_veff(mc.mol,core_dm)

    #Calcolate one-electron integrals for the (0') class
    h1e_Sir = reduce(numpy.dot, (mo_virt.T,mc.get_hcore()+core_vhf,mo_core))

    #Calcolate one-electron integrals for the (+1') class
    h1e_Si =  reduce(numpy.dot, (mo_cas.T, mc.get_hcore()+core_vhf , mo_core))

    #Calcolate one-electron integrals for the (-1') class
    h1e_Sr =  reduce(numpy.dot, (mo_virt.T,mc.get_hcore()+core_vhf , mo_cas))

    #Calculate orbital energies in the outer space to build the Dyall Hamiltonian
    #FIXME
    #need onepdm to get refined orbital energy of non-active orbitals.
    dm1 = mcscf.make_rdm1(mc)
    total_vhf = mc.get_veff(mc.mol,dm1)
    fock = reduce(numpy.dot, (mc.mo_coeff.T,mc.get_hcore()+total_vhf , mc.mo_coeff))

    orbe = fock.diagonal()

    #Calculate the C constant of the Dyall Hamiltonian
    C = energy_core

    for i in range(mc.ncore):
        C -= 2.0*orbe[i]

    ncore = mc.ncore
    ncas = mc.ncas
    nvirt = mo_virt.shape[1]
    sz = mc.ncas
    tol = float(1e-15)

    #Take care of the orbitals order
    reorder = range(mc.mo_coeff.shape[0])
    core_order = reorder[:mc.ncore]
    cas_order = reorder[mc.ncore:mc.ncore+mc.ncas]
    virt_order = reorder[mc.ncore+mc.ncas:]

    cas_order_init = [x - ncore for x in cas_order]

    #Set permutational symmetry
    perm_symm = True    
    if (nevpt_class != None):
       perm_symm = False

    #Dump the Total FCIDUMP files as a perturber for all the NEVPT2 Classes
    h2e_all = ao2mo.incore.general(mc._scf._eri,[mc.mo_coeff,mc.mo_coeff,mc.mo_coeff,mc.mo_coeff],compact=False)
    h2e_all = h2e_all.reshape(ncore+ncas+nvirt,ncore+ncas+nvirt,ncore+ncas+nvirt,ncore+ncas+nvirt)
    h1e_all = reduce(numpy.dot, (mc.mo_coeff.T,mc.get_hcore(), mc.mo_coeff))

    if (nevpt_class == None):
        f = open('PERTURB_TOT','w')
        print >>f, ' &FCI NORB=',ncore+ncas+nvirt,', NELEC= ', mc.mol.nelectron, ' ,MS2= ',0,','
        f.write(' ORBSYM=')
        for i in xrange(ncore+ncas+nvirt):
            f.write("%d," % orbsym[i])
        f.write('\n')
        print >>f, ' &END'
        #h2e in active space
        if (perm_symm == False):
           writeh2e(h2e_all,f,(reorder,reorder,reorder,reorder),tol)
        else:
           writeh2e_sym(h2e_all,f,(reorder,reorder,reorder,reorder),tol)
        #h1e in active space
        if (perm_symm == False):
           writeh1e(h1e_all,f,(reorder,reorder),tol)
        else:
           writeh1e_sym(h1e_all,f,(reorder,reorder),tol)
        f.close()

    #Dump the integrals for DYALL hamiltonian
    f = open('DYALL','w')
    print >>f, ' &FCI NORB=',ncore+ncas+nvirt,', NELEC= ', mc.mol.nelectron, ' ,MS2= ',0,',' 
    f.write(' ORBSYM=')
    for i in xrange(ncore+ncas+nvirt):
        f.write("%d," % orbsym[i])
    f.write('\n')
    print >>f, ' &END'
    #h2e in active space
    if (perm_symm == False):
       writeh2e(h2e,f,(cas_order,cas_order,cas_order,cas_order),tol)
    else:
       writeh2e_sym(h2e,f,(cas_order,cas_order,cas_order,cas_order),tol)
    #h1e in active space
    for i in xrange(ncore):
        print >>f, ' % .12e  %3d  %3d  %3d  %3d' % (orbe[i], i+1, i+1, 0, 0)

    if (perm_symm == False):
       writeh1e(h1e,f,(cas_order,cas_order),tol)
    else:
       writeh1e_sym(h1e,f,(cas_order,cas_order),tol)

    for i in xrange(nvirt):
        print >>f, ' % .12e  %3d  %3d  %3d  %3d' % (orbe[ncore+ncas+i], ncore+ncas+i+1, ncore+ncas+i+1, 0, 0)

    print >> f, C,0,0,0,0
    f.close()

    if (nevpt_class != None):
        f = open('PERTURB','w')
        print >>f, ' &FCI NORB=',ncore+ncas+nvirt,', NELEC= ', mc.mol.nelectron, ' ,MS2= ',0,','
        f.write(' ORBSYM=')
        for i in xrange(ncore+ncas+nvirt):
            f.write("%d," % orbsym[i])
        f.write('\n')
        print >>f, ' &END'
        if (nevpt_class == "0"): 
           writeh2e(h2e_Sijrs,f,(virt_order,core_order,virt_order,core_order),tol)
 
        if (nevpt_class == "+1"): 
           writeh2e(h2e_Sijr,f,(virt_order,core_order,cas_order,core_order),tol)
 
        if (nevpt_class == "-1"): 
           writeh2e(h2e_Srsi,f,(virt_order,core_order,virt_order,cas_order),tol)
 
        if (nevpt_class == "+2"): 
           writeh2e(h2e_Sij,f,(cas_order,core_order,cas_order,core_order),tol)
 
        if (nevpt_class == "-2"): 
           writeh2e(h2e_Srs,f,(virt_order,cas_order,virt_order,cas_order),tol)
 
        if (nevpt_class == "0p"): 
           writeh2e(h2e_Sir1,f,(virt_order,core_order,cas_order,cas_order),tol)
           writeh2e(h2e_Sir2,f,(virt_order,cas_order,cas_order,core_order),tol)
 
        if (nevpt_class == "+1p"): 
           writeh2e(h2e_Si,f,(cas_order,core_order,cas_order,cas_order),tol)
 
        if (nevpt_class == "-1p"): 
           writeh2e(h2e_Sr,f,(virt_order,cas_order,cas_order,cas_order),tol)
 
        if (nevpt_class == "0p"): 
           writeh1e(h1e_Sir,f,(virt_order,core_order),tol)
 
        if (nevpt_class == "+1p"): 
           writeh1e(h1e_Si,f,(cas_order,core_order),tol)
 
        if (nevpt_class == "-1p"): 
           writeh1e(h1e_Sr,f,(virt_order,cas_order),tol)

        print "Remember to switch off permutational symmetry and screening in Block" 
 
        f.close()









    

    

if __name__ == '__main__':
    from functools import reduce
    from pyscf import gto
    from pyscf import scf
    from pyscf import ao2mo
    from pyscf import fci
    from pyscf import mcscf
    from pyscf import symm
    from pyscf.mrpt import nevpt2

    mol = gto.Mole()
    mol.verbose = 5
    mol.atom = [
        ['N', ( 0., 0.    , -.9  )],
        ['N', ( 0., 0.    ,  .9  )],
    ]
    mol.symmetry = 'D2h'
    mol.basis = {'N': 'cc-pVDZ'}
    mol.build()

    m = scf.RHF(mol)

    m.kernel()
    irrep_name = m.mol.irrep_id
    sym = symm.label_orb_symm(m.mol, irrep_name,m.mol.symm_orb,m.mo_coeff)
    orbsym = [IRREP_MAP['D2h'][i % 10] for i in sym]

    mc=mcscf.CASCI(m,8,10)
    #mc=mcscf.CASSCF(m,8,10)

    #ci_e = mc.kernel()[0]
    #print 'CI energy', ci_e
    mc.casci()

    nevpt_class = '-1p'

    dyall_perturb(mc, orbsym, nevpt_class)
    #nevpt2.NEVPT(mc).kernel()
