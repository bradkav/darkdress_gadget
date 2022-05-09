import numpy as np
import pygadgetic

import units
import math

from tqdm import tqdm

import os.path
import os
from scipy.integrate import quad, cumtrapz
from scipy.interpolate import interp1d

from scipy.optimize import brenth

from matplotlib import pylab as pl



#import mpi4py.MPI

#-------------
#You should just be able to import whichever eddington module
#you're interested in and everything should work from here on...
import eddington_EMRI as edd

#-------------

L_sim = 1.0*units.L_M #Simulation units in pc   

#--------------------------------------
#------------------EMRI----------------
#--------------------------------------
def AddDressedPBH_EMRI( x0, v0, r_soft, M_IMBH, rho0, gamma, nDM = 100, verbose=False, haloID="nothing", pure=False):
    """Add a dressed PBH to the initial conditions...
    
    Parameters:
        body - the pygadgetic 'body' object (usually called my_body)
        DMinds - indices specifying where to put the DM particles
                in the list (the DM particles usually come before the PBH
                particles)
        PBHind - single index specifying where to place the PBH
                in the list (close to the end usually, so -1 or -2)
        nDM    - number of DM particles around this PBH
        x0     - initial position of PBH (in pc)
        v0     - initial velocity of PBH+DM halo (in km/s)
        r_soft - softening length in parsec
        a      - Semi major axis in parsec (to 3 decimal places). 
                 Tabulated values are a = [0.001, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08] 
        haloID - string identifying a text file to load the halo from (in folder /halos)
                 If file not found, a new halo is generated and saved in /halos.
                 Set haloID = "nothing" to ignore this option.
        pure   - Determine whether to use a pure power law with index gamma
                 (or a soothly truncated one). pure=True uses a pure power law.
    """

    fix_CoM = True
    #TO BE ADDED AS PARAMETERS
    #N_inner = 100
    #delta_Rm = 50
    #PBH mass and truncation radius imported from the eddington file for
    #self-consistency
    edd.loadDistribution(M_IMBH, rho0, gamma, pure)
    
    r_min = 2.8*r_soft
    r_max = 1e-6 #What should I set this to...?

    r_isco = 6*M_IMBH*units.G_N/units.C_LIGHT**2

        

    print(r_isco, r_min)

    if (r_min < r_isco):
        r_min = r_isco
    
    #print mHalo
    r_sp = edd.r_sp
    mHalo = edd.Menc(r_max) - M_IMBH

    
    
    #Initialise the masses, positions and velocities
    m_vals = np.zeros(nDM + 1)
    pos_vals = np.zeros((nDM + 1,3))
    vel_vals = np.zeros((nDM + 1,3))
    
    #Set up the central black hole
    m_vals[0] = M_IMBH
    pos_vals[0,:] = np.zeros(3)
    vel_vals[0,:] = np.zeros(3)
    
    
    #PBH position and velocity (before CoM velocity is subtracted...)
    xPBH=np.array([0.,0.,0.])
    vPBH=np.array([0.,0.,0.])
    
    halofile = "halos/" + haloID + ".txt"
    
    #Check to see whether a halo file already exists...
    if (haloID != "nothing" and os.path.isfile("halos/" + haloID + ".txt")):
        print(("    Loading halo from file. HaloID:", haloID))
        
        #Load DM phase space coordinates from file
        xvals, yvals, zvals, vxvals, vyvals, vzvals, mvals = np.loadtxt(halofile, unpack=True)
        body.mass[DMinds] = mvals
    else:
        if (haloID != "nothing"):
            print(("   Halo file <" + halofile+"> not found. Generating from scratch..."))
    
        #Generate the mass profile
        print("   Generating mass profile...")
    
        #rlist = np.logspace(np.log10(r_min), np.log10(r_max), 500)
        rlist = np.logspace(np.log10(r_isco), np.log10(r_max), 1000)
        Menc = 0.0*rlist

        #Radial distribution function for the DM halo
        P_r_1 = lambda r: 4.0*np.pi*r**2*edd.rhoDM(r)
        P_r = np.vectorize(P_r_1)

        for i in range(len(rlist)):
            Menc[i] = edd.Menc(rlist[i]) - M_IMBH
            #Menc[i] = quad(P_r, r_min, rlist[i])[0]

        #print("Menc[0] = ", Menc[0])
        #Menc -= Menc[0]
        M_max = Menc[-1]

        
        Minterp = interp1d(Menc/M_max, rlist, kind='linear')
        
        

        #DM positions
        #Cut off a fraction p of the radial distribution
        #near the PBH (fraction inside r_soft)
        p_cut = (edd.Menc(r_min) - M_IMBH)/M_max
        print(p_cut)
        rvals = Minterp((p_cut + np.random.rand(nDM))/(1+p_cut))
        
            
        print("   Minimum particle radius [pc]", np.min(rvals))
        print("   Number of particles below r  = 1e-7 pc:", np.sum(rvals < 1e-7))
        print("   Number of particles below r  = 1e-8 pc:", np.sum(rvals < 1e-8))
        print("   Number of particles below r  = 1e-9 pc:", np.sum(rvals < 1e-9))
        print("  ")
    
        #Generate some random directions for setting particle positions
        ctvals = 2.0*np.random.rand(nDM) - 1.0
        thetavals = np.arccos(ctvals)
        phivals = 2*np.pi*np.random.rand(nDM)

        xvals = rvals*np.cos(phivals)*np.sin(thetavals)
        yvals = rvals*np.sin(phivals)*np.sin(thetavals)
        zvals = rvals*np.cos(thetavals)

        pos_vals[1:,:] = np.array([xvals, yvals, zvals]).T
        m_vals[1:] = mHalo/nDM
        
        #Deal with C.o.M.    
        totmass = np.sum(m_vals)
        x_CoM = np.zeros(3)
        x_CoM[0] = np.sum(pos_vals[:,0]*m_vals)
        x_CoM[1] = np.sum(pos_vals[:,1]*m_vals)
        x_CoM[2] = np.sum(pos_vals[:,2]*m_vals)
        x_CoM /= totmass
        
        if (fix_CoM):                
            pos_vals -= x_CoM

            #Now rotate so that IMBH is in x-y plane
            x_BH = pos_vals[0,:]/L_sim
            print(x_BH)
            theta_rot = np.arctan(-x_BH[2]/x_BH[1])
            print(theta_rot)
            pos_old = 1.0*pos_vals
            pos_vals[:,0] = pos_old[:,0] #Keep x coordinates the same
            pos_vals[:,1] = np.cos(theta_rot)*pos_old[:,1] - np.sin(theta_rot)*pos_old[:,2]
            pos_vals[:,2] = np.sin(theta_rot)*pos_old[:,1] + np.cos(theta_rot)*pos_old[:,2]

            x_BH = pos_vals[0,:]/L_sim
            print(x_BH)
        
        print("    ")
        print("    Total CoM position:", np.sqrt(np.sum(x_CoM**2)))
        print("    Total Halo mass:", totmass - M_IMBH)
        print("    ")
            
            
        #-------------
        disp_env = os.environ.get('DISPLAY')        
        #cols = ['r','g','b','c']
        
        if (disp_env is not None):
            pl.figure()
            ax = pl.gca()
            ax.set_xscale("log")
            pl.hist(rvals, bins=np.logspace(-15, -0, 151), alpha = 0.5)
        
            pl.axvline(r_soft, linestyle='--', color='k')
            pl.axvline(r_max, linestyle='--', color='k')
        
            #pl.axvline(np.log10(r_outer), linestyle=':', color='k')
            pl.show()
        #---------
            
        vvals = np.zeros(nDM)
        for ind in tqdm(range(nDM), desc="Sampling velocities..."):
            r = rvals[ind]
            #Now sample f(v) at given r to get the speed v
            found = 0
    
            while (found == 0):
                
                v = np.random.rand(1)*edd.vmax(r)
                #Use 5/vmax as the 'maximum' values of f(v)
                #but in some cases it might not be enough...
                if (np.random.rand(1)*(5.0/edd.vmax(r)) < edd.f(r, v)):
                    #pl.show()
                    found = 1
                    vvals[ind] = v

        #Get a new set of random directions for the velocities
        ctvals2 = 2.0*np.random.rand(nDM) - 1.0
        thetavals2 = np.arccos(ctvals2)
        phivals2 = 2*np.pi*np.random.rand(nDM)

        vxvals = vvals*np.cos(phivals2)*np.sin(thetavals2)
        vyvals = vvals*np.sin(phivals2)*np.sin(thetavals2)
        vzvals = vvals*np.cos(thetavals2)

        
        
        vel_vals[1:,:] = np.array([vxvals, vyvals, vzvals]).T
        
    p_CoM = np.zeros(3)
    p_CoM[0] = np.sum(vel_vals[:,0]*m_vals)
    p_CoM[1] = np.sum(vel_vals[:,1]*m_vals)
    p_CoM[2] = np.sum(vel_vals[:,2]*m_vals)

    if (fix_CoM):
        vel_vals -= p_CoM/totmass
    
    print(("   v_PBH [pc/kyr]:", np.sqrt(np.sum(vel_vals[0,:]**2))*3.24078e-14*(3600*24.0*365*1000)))
     
    
    #Save the output to a halo file if needed
    if (haloID != "nothing"):
        headertxt = "Number of DM particles: " + str(nDM) + ". Softening length [pc]: " + str(r_soft)
        headertxt += "\nColumns: x [pc], y [pc], z [pc], vx [km/s], vy [km/s], vz [km/s], m [M_solar]"
        np.savetxt("halos/" + haloID + ".txt", zip(pos_vals[1:,0],pos_vals[1:,1],pos_vals[1:,2],vel_vals[1:,0],vel_vals[1:,1],vel_vals[1:,2], m_vals[1:]), header=headertxt)
    
    
    
    #Add on the CoM position and velocity
    pos_vals += np.asarray(x0)
    vel_vals += v0
    
    pos_vals /= L_sim

    #print("Mtot = ", np.sum(m_vals))
    #print("CoM positions [M]:", np.sum(pos_vals*np.atleast_2d(m_vals).T, axis=0)/np.sum(m_vals))
    
    return m_vals, pos_vals, vel_vals
