import numpy as np
import scipy.linalg
import pyscf.gto
from pyscf.gto import ft_ao
import pyscf.scf
import pyscf.mcscf
import pyscf.ao2mo
import pyscf.lib.logger as logger
import math
from pyscf import fcigf,lo,fci
from pyscf.fci import cistring,fci_slow

################################################################################
#This script calculates the complex polarization functions as described in 
#Ref: PRB,84,245117 (2011). Some approximations due to resolutions 
#of the identity introduced during the evaluation of the exponential and 
#in particular to Open Boundary Conditions (OBC) of the system have been used.
#
#A slower increase of the function to 1 due to OBC is usually observed.
################################################################################

def run(r, N):
    mol = pyscf.gto.Mole()
    mol.verbose = 0
    mol.unit = 'au'
    mol.atom = [
        ['H', ( 0., 0.    ,  r*i   )] for  i in range(N)
    ]
    mol.basis = {'H': 'sto-6g'}
    mol.build()
    
    mf = pyscf.scf.RHF(mol)
    mf.conv_threshold = 1e-10
    ehf = mf.scf()
    dm = mf.make_rdm1()
    #print "SCF energy:", ehf, "\n"

    ovlp = mol.intor('cint1e_ovlp_sph')
    LoOrb = lo.orth.lowdin(ovlp)

    #e, c = scipy.linalg.eigh(ovlp)
    #V = np.dot(c*1./(np.sqrt(e)), c.T)

    #s = (2*np.pi/N)*np.arange(-N/2+0.5,+N/2-0.4)
    #X = np.dot(ovlp, V)
    #M = np.dot(X, np.dot(np.diag(np.exp(+1j*s)), X.T))

    L = (N*r)

    #ov = scipy.linalg.det(np.dot(mo_occ.T, np.dot(M, mo_occ)))
    #mo_occ = mf.mo_coeff[:,mf.mo_occ>1]
    #mod_comppol = np.abs(ov*ov)

    mo_occ = mf.mo_occ>1

    #Calculate the z component (parallel to the chain) of the dipole moment integrals
    #in the  AO basis
    dipz = mol.intor_symmetric('cint1e_r_sph', comp=3)[2]

    #The dipole moment integrals are transformed in the MO basis.
    mo_dipz = reduce(np.dot,(mf.mo_coeff.T, dipz, mf.mo_coeff))

    #The dipole moment integrals (MO basis) matrix is diagonalized.
    e, c = scipy.linalg.eigh(mo_dipz)

    #The exponential is calculated.
    comppol_ints = reduce(np.dot, (c, np.diag(np.exp(1j*(2*np.pi/L)*e)), c.T))
    U = comppol_ints[np.ix_(mo_occ,mo_occ)]
    comppol = np.linalg.det(U)
    comppol = comppol*comppol

    #It takes the absolute value of the complex polarization function.
    mod_comppol = np.absolute(comppol)

    print r, mod_comppol
    

    
    

start = 1
stop = 100

for i in range(start, stop, -1):
    r = 0.1*i
    N = 10
    #print "R = ", r, " ang"
    run(r, N)

