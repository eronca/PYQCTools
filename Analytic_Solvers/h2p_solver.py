import numpy as np
import scipy
from scipy.integrate import odeint
from scipy.optimize import newton, leastsq
import matplotlib.pyplot as plt
import matplotlib as mpl

##################################################################################
## The analytic solution of H2+ has been obtained following the procedure
## described in Grivet, J. Chem. Edu., 79, 127, (2002). 
##################################################################################

font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 30,
        }

mpl.rcParams['xtick.labelsize'] = 30
mpl.rcParams['ytick.labelsize'] = 30

mpl.rc('font', family='serif')
mpl.rc('font', serif='Times')
mpl.rc('text', usetex='true')
mpl.rcParams.update({'font.size': 30})

def f(z,x,c2,L,m):
   return [z[1], (2.*(m+1)*x*z[1]+(m*(m+1)+c2*x*x-L)*z[0])/(1-x*x)]

def g(z,x,c2,L,D,m):
   return [z[1], (-2.*(m+1)*x*z[1]-(2.*D*x+m*(m+1)-L+c2*x*x)*z[0])/(x*x-1)]

def S_eta(L_vec, c2_val):

   c2 = c2_val

   traj = np.zeros([x_ang.shape[0]*2+1,L_vec.shape[0]])
   n=0
   for L in L_vec:
   
      slope = (m*(m+1)+c2-L)/(m+1)/2.
      z0 = [1+step*slope,slope]
   
      z = odeint(f,z0,x_ang,args=(c2,L,m))
   
      temp=1-pow(x_ang,2.0)
      temp=pow(temp,m/2.)
   
      zz=temp*z[:,0]
   
      first_zz = np.array([1])
      zz=np.append(first_zz, zz)
   
      sloper = -(m*(m+1)+c2-L)/(m+1)/2.
      z0r = [1-step*sloper,sloper]
   
      zr = odeint(f,z0r,x_angr,args=(c2,L,m))
   
      zzr=temp*zr[:,0]
      zzr=np.append(first_zz, zzr)
   
      traj[:,n] = np.append(zz,zzr[::-1][1:])
   
      n += 1
   
   x_tot = np.arange(-1,1.+step,step)
   
   figure = plt.figure(figsize=(12, 11))
   plt.plot(x_tot,traj, linewidth=2.0, label = '')
   plt.ylabel('S($\\eta$)')#, fontdict=font)
   plt.xlabel('$\\eta$')#, fontdict=font)
   plt.xlim(-1.0,1.0)
   plt.ylim(0.3,1.0)
   plt.locator_params(axis='x', nbins=10)
   plt.locator_params(axis='y', nbins=10)
   plt.tick_params(axis='x', pad=15)
   #plt.legend(loc=1)
   
   plt.show()
   plt.close()


def newton_ang_func(L_val,c2,m):

   L = L_val

   slope = (m*(m+1)+c2-L)/(m+1)/2.
   z0 = [1+step*slope,slope]

   z = odeint(f,z0,x_ang,args=(c2,L,m))

   temp=1-pow(x_ang,2.0)
   temp=pow(temp,m/2.)

   zz=temp*z[:,0]

   first_zz = np.array([1])
   zz=np.append(first_zz, zz)

   sloper = -(m*(m+1)+c2-L)/(m+1)/2.
   z0r = [1-step*sloper,sloper]

   zr = odeint(f,z0r,x_angr,args=(c2,L,m))

   zzr=temp*zr[:,0]
   zzr=np.append(first_zz, zzr)
 
   return  z[:,1][-1]

def R_xi(E_vec,L0):
    
   traj = np.zeros([x_rad.shape[0]+1,E_vec.shape[0]])

   n=0
   for E in E_vec:

      c2 = D*D*E/4.

      L = newton(newton_ang_func,L0,args=(c2,m),tol=1e-8,maxiter=200)

      slope = -(-L+m*(m+1.)+2.*D+c2)/(2.*(m+1.))
      z0 = [1+step*slope,slope]

      z = odeint(g,z0,x_rad,args=(c2,L,D,m))

      temp=pow(x_rad,2.0)-1.
      temp=pow(temp,m/2.)

      zz=temp*z[:,0]

      first_zz = np.array([1])
      zz=np.append(first_zz, zz)

      traj[:,n] = zz

      n += 1

   xt = np.append(np.array([1]),x_rad)

   figure = plt.figure(figsize=(12, 11))
   plt.plot(xt,traj, linewidth=2.0, label = '')
   plt.ylabel('R($\\xi$)')#, fontdict=font)
   plt.xlabel('$\\xi$')#, fontdict=font)
   plt.xlim(1.0,10.0)
   plt.ylim(-1.0,1.0)
   plt.locator_params(axis='x', nbins=10)
   plt.locator_params(axis='y', nbins=10)
   plt.tick_params(axis='x', pad=15)
   #plt.legend(loc=1)

   plt.show()
   plt.close()

def newton_rad_func(E_val,D,m,L0):

   E = E_val

   c2 = D*D*E/4.

   L = newton(newton_ang_func,L0,args=(c2,m),tol=1e-8,maxiter=500)

   slope = -(-L+m*(m+1.)+2.*D+c2)/(2.*(m+1.))

   z0 = [1+step*slope,slope]

   z = odeint(g,z0,x_rad,args=(c2,L,D,m))

   temp=pow(x_rad,2.0)-1.
   temp=pow(temp,m/2.)

   zz=temp*z[:,0]

   first_zz = np.array([1])
   zz=np.append(first_zz, zz)

   return z[:,1][-1]

def L_calc(E_val,D,m,L0):
   E = E_val
   c2 = D*D*E/4.
   L = newton(newton_ang_func,L0,args=(c2,m),tol=1e-8,maxiter=500)
   return L

def pes(D_vec,m,L0,E_start):

   PES = np.zeros((D_vec.shape[0],D_vec.shape[0]))
   i = 0
   for D in D_vec:    
     E_elec = newton(newton_rad_func,E_start,args=(D,m,L0),tol=1e-8,maxiter=200)
     L0 = L_calc(E_elec,D,m,L0)
     E_nuc = 2./D
     E_start = E_elec
     PES[i][0] = D
     PES[i][1] = E_elec+E_nuc
     print PES[i][0], PES[i][1]
     i += 1
   
   return PES


if __name__ == '__main__':

   #D = 1.4
   m=0

   step = 0.001
   x_ang=np.arange(-1.+step,0.+step,step)
   x_angr=np.arange(1.-step,0.-step,-step)
   x_rad = np.arange(1.+step,10.+step,step)

   ##Angular Part S-eta plot generator
   #L_vec = np.array([-1.8, -1.7, -1.6, -1.59449, -1.5, -1.4])
   #S_eta(L_vec, -4)

   ##Radial Part R-xi plot generator
   #E_vec = np.array([-2.3,-2.5,-2.6,-2.7])
   #R_xi(E_vec,-1.8)

   ##Calculation of the electronic energy
   #E_elec = newton(newton_rad_func,-2.5,args=(D,m,-1.8),tol=1e-8,maxiter=200)

   ##Calculation of the nuclear energy
   #E_nuc = 2./D

   ##Calculation of the total energy
   #E_tot = E_elec+E_nuc
   #print "The Total Energy is:", E_tot, "Ryd"

   ##Calculate PES 
   D_vec = np.arange(1.0,40.0, 0.1)
   PES = pes(D_vec,m,-1.8,-2.5)
   #print PES
