#Eddington distribution function for the truncated density profile

import numpy as np
import sympy as sp
from scipy.interpolate import interp1d
from scipy.special import hyp2f1
import matplotlib.pylab as pl
from scipy.integrate import quad

G_N = 4.302e-3 #(pc/solar mass) (km/s)^2
c_light = 2.99792458e5 # km/s

func = None
Menc_fun = None
psi_fun = None
r_sp = 1.0
M_IMBH = 0.0
gamma = 0.0
rho_sp = 1.0

alpha = 0.0
r_t = 0.0

rmin = 0.0
psi_interp = None


def loadDistribution(M_IMBH_in, rho0_in, gamma_in, pure=False):
    global func, r_sp, rho_sp, M_IMBH, gamma, alpha, r_t, rmin, psi_interp
    
    M_IMBH = 1.0*M_IMBH_in
    rho_sp = 1.0*rho0_in
    gamma = 1.0*gamma_in
    
    r_sp = ((3-gamma)*(0.2**(3.0-gamma))*M_IMBH/(2*np.pi*rho_sp))**(1.0/3.0)
    print("r_sp [pc]:", r_sp)
    fid = "distributions/distribution_M=" + str(int(M_IMBH)) + "_rho0=" + "{0:.2f}".format(rho_sp) + "_gamma=" + "{0:.2f}".format(gamma)
    if (pure):
        fname = fid + "_pure.dat"
    else:
        fname = fid + ".dat"

    #fname = "distributions/distribution_M=" + str(int(M_IMBH)) + "_rho0=" + "{0:.2f}".format(rho_sp) + "_gamma=" + "{0:.2f}".format(gamma) + ".dat"
    
    fdat = np.loadtxt(fname)
    alpha = fdat[0,0]
    r_t = fdat[0,1]
    print("alpha = ", alpha)
    func = interp1d(fdat[1:,0], fdat[1:,1], kind='linear', bounds_error=False, fill_value=0.0)
    
    rmin = 6*M_IMBH*G_N/c_light**2
    
    print("Calculating potential...")
    rlist = np.logspace(np.log10(rmin), np.log10(10*r_sp),1000)
    
    psilist = np.asarray([psi_1(r) for r in rlist])
    psi_interp = interp1d(rlist, psilist, bounds_error=False, fill_value=0.0)
    
    #Redefine some sympy stuff
    #Menc_fun = r_tr**3*4*np.pi*sp.integrate(rhoDM(x1)*x1**2, (x1, 0, x2))
    #psi_fun = -r_tr**2*G_N*sp.integrate((M_PBH + Menc_fun)/x2**2, (x2, 1000, x3))
    
# Density
def rho_DM(r):
    return rho_sp*(r/r_sp)**-gamma/(1+r/r_t)**alpha
   
   
def rho_DM_true(r):
    return rho_sp*(r/r_sp)**-gamma
    
#Enclosed mass
def Menc_fun(r):
    return (4*np.pi*rho_sp/(3-gamma))*r**3*(r/r_sp)**-gamma*hyp2f1(alpha, 3-gamma, 4-gamma, -r/r_t)

def Menc(r):
    return M_IMBH + np.clip(Menc_fun(r) - Menc_fun(rmin), 0, 1e30)
    
#Potential
def psi_1(r):
    #print(r)
    integ = lambda x: Menc(x)/x**2
    
    if (r < r_sp):
        points = np.logspace(np.log10(10*r_sp), np.log10(r), 20)
    else: 
        points = None
    
    return -G_N*quad(integ, 100*r_sp, r, points=points)[0]

psi = np.vectorize(psi_1)


#Maximum speed at a given radius x = r/r_tr
def vmax(r):
    return np.sqrt(2.0*psi_interp(r))

    
#Speed distribution f(v) at a given radius r
def f_scalar(r, v):
    if (v >= vmax(r)):
        return 0.0
    else:
        return 4.0*np.pi*(v**2)*func(psi_interp(r) - 0.5*v**2)/rho_DM(r)

f = np.vectorize(f_scalar)


