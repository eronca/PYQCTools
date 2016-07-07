Tools for Omega Space Green's Functions
==============================================

:attr:`gf_trace.py`: Calculate the Density of States (DOS) value from an :math:`\omega`-dependent Green's Function. 
                     It makes the trace of the Green's Function associated with a certain frequency value. 

                     **Example**::

                       from PYQCTools.Omega_GF import gf_trace           

                       gf_trace.run(green.txt, omega_value)
                 
                     :attr:`green.txt`: formatted text file containing the Green's function, :attr:`omega_value`: double frequency value.
