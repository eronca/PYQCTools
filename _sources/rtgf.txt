Tools for Real-Time Green's Functions
==============================================

:attr:`rtgf.py`: Calculate the Density of States (DOS) values during a time propagation.
                 It makes the trace of time-dependent Green's Functions calculated 
                 along a time propagation and return the DOS values as a function of time both for the real and for the imaginary part of the Green's Function.  

                     **Example :**

                     >>> python rtgf.py prop_time time_step /PATH/scratch

                     :attr:`prop_time`: double value of the full propagation time (period).

                     :attr:`time_step`: double value of the time-step.

                     :attr:`scratch`: directory containing text files of the real (green.$t.$t.txt) and imaginary (green.30000+$t.30000+$t.txt) Green's Functions, where $t indicate the specific time-step. 

                 The script save two output files, :attr:`rt_real.txt` and :attr:`rt_imag.txt`, containing the real and imaginary part of the time-dependent DOS respectively. 
