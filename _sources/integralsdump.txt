Tools for Integrals Dumping
==============================================
                     
:attr:`Integrals_dump.py`: It dumps 1 and 2-elctron integrals in the MO basis inside a CASCI space in FCIDUMP format.
                           The PySCF input to calculate the integrals is already included in the script.

:attr:`DipoleIntegrals_dump.py`: It dumps dipole integrals in the MO basis in a CASCI space in FCIDUMP format.
                                 The PySCF input to calculate the integrals is already included in the script.

:attr:`LowdinOrtho_Integrals.py`: It dumps 1 and 2-electron integrals in the Localized basis obtained by Lowdin Orthogonalization in FCIDUMP format.
                                  The PySCF input to calculate the integrals is already included in the script.

:attr:`hubbard_1d`: It dumps 1 and 2-electron integrals got the 1D Hubbard model.
                     
                    **Example**::
                   
                      from PYQCTools.Integrals_dump import hubbard_1d		

                      hubbard_1d.run(nsites, t, U, output, pbc)		

                    :attr:`nsites`: Number of sites, :attr:`t`: Hopping constant, :attr:`U`: Coupling constant, :attr:`output`: Output file name, :attr:`pbc`: 'True' of 'False' respectively if periodic boudary conditions need to be included or not.

:attr:`MPSPT_integrals.py`: It dumps the DYALL and PERTURB files to run NEVPT2 calculations by MPS-PT using the BLOCK Code for DMRG.
                            Integrals are dumped in FCIDUMP format.
                            The PySCF input to calculate the integrals is already included in the script.
