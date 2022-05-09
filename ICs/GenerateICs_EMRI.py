### this is an example of how to generate an initial condition file

import numpy as np
import pygadgetic
#-------------

import eddington_EMRI as edd
from scipy.interpolate import interp1d
import argparse

import units
from collections import Counter

#This is the module that add the PBHs to the initial condition file
import PBH


#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-lN','--log2N', help='Log2 of number of DM particle per PBH', type=int, default=15)
parser.add_argument('-rs','--r_soft', help='Softening length in km', type=float, default=1e-3)

args = parser.parse_args()
#e = args.ecc
lN = args.log2N
r_soft = args.r_soft

M_PBH = 1000



rho0 = 226.0
gamma = 7.0/3.0

edd.loadDistribution(M_PBH, rho0, gamma)

r_sp = edd.r_sp

#Number of DM particles per halo
nDM = 2**lN

print("  Generating ICs with %d DM pseudo-particles per PBH..."%nDM)


#L_sim = 3.24078e-14 #simulation units in pc
L_sim = 1.0*self.L_M

#Calculate the truncation radius (just for information

delta_Rm = 1.0

N_inner = nDM

print("  Softening length (km):", r_soft)
print("  Softening length (pc):", r_soft*L_sim)
print("  ")
print("  PBH mass (M_solar):", M_PBH)
print("  Spike radius r_tr (km):", r_sp/L_sim)
print("  Spike radius r_tr (pc):", r_sp)
print("  ")
#print "  Apoapsis:", (1+e)*a
#print "  Periapsis:",(1-e)*a

"""

try:
    f = open('../run/ICs.txt','w')
    IDstr = "bin_N" + str(lN) + "_s" + str(x_soft/1e-2)+"_a"+str(a) + "_e"+str(e)
    f.write('ID    '+IDstr + '\n')
    f.write('TYPE    bin\n')
    f.write('lN    %d\n'%(lN,))
    f.write('x_soft    %f\n'%(x_soft,))
    f.write('a    %f\n'%(a,))
    f.write('e    %f\n'%(e,))
    f.close()
except IOError:
    print "File '../run/ICs.txt' not found - continuing anyway..."
"""
print("  ")





#mu in units of pc (km/s)^2
#a in units of pc
#T = 2*np.pi*np.sqrt(a**3*(3.086e13)**2/mu)/3.15576e13 
#print "  Orbital period (Myr):", T
#print "  Orbital period (simulation time):",T/(0.019139*scaling)

#vapo = 0.5*np.sqrt((1.0-e)*mu/((1.0+e)*a))
#print "  Initial speed:", vapo


#Orbital parameters (for circular orbit)
d0 = 1.084e-8
print("   Initial BH separation (pc):",edd.Menc(d0))

print("   DM density at d0 (M_sun/pc^3):", edd.rho_DM(d0))
print("   True DM density at d0 (M_sun/pc^3):", edd.rho_DM_true(d0))

v0 = np.sqrt(edd.G_N*(1.0 + edd.Menc(d0))/d0)

haloID = "halo_M=" + str(int(M_PBH)) + "_rho0=" + "{0:.2f}".format(rho0) + "_gamma=" + "{0:.2f}".format(gamma) + "_rs=" + str(r_soft)

mlist, xlist, vlist = PBH.AddDressedPBH_EMRI([d0*1.0/(M_PBH+1.0), 0, 0],[0, v0*1.0/(M_PBH+1.0), 0],  r_soft*L_sim, M_PBH, rho0, gamma, nDM, haloID="nothing", verbose=True)
#PBH.AddDressedPBH_seg(my_body,np.arange(0,nDM),-1, nDM, [0, 0, 0],[0, 0, 0], r_soft, a, r_seg = np.sqrt(r_soft*r_tr), N_ratio = N_ratio, haloID=haloID1, verbose=True)



inds = (mlist.argsort())[::-1]
xlist_sorted = xlist[inds,:]
vlist_sorted = vlist[inds,:]
mlist_sorted = mlist[inds]


#Add the other BH
mlist_sorted = np.append(mlist_sorted, 1.0)
xlist_sorted = np.append(xlist_sorted, np.atleast_2d([-(d0/L_sim)*M_PBH/(M_PBH + 1.0), 0, 0]), axis=0)
vlist_sorted = np.append(vlist_sorted, np.atleast_2d([0, -v0*M_PBH/(M_PBH + 1.0), 0]), axis=0)

cnt = Counter(mlist_sorted)

mvals = np.array([k for k, v in cnt.items()])
n = np.array([v for k, v in cnt.items()])

n_species = len(mvals)
n_particles = np.zeros(n_species, dtype='int')
n_particles[0] = 1
n_particles[1] = nDM
n_particles[2] = 1
#n_particles = n[mvals.argsort()[::-1]]
if (np.sum(n) != nDM + 2):
    raise ValueError("Number of particles doesn't seem to match...")

##define number of particles
npart = np.zeros(6, dtype='int')
npart[1:(n_species+1)] = n_particles

print(npart)

#for i in range(N_shell):
#    npart[i+2] = nDM_shell[i]

total_number_of_particles=np.sum(npart) #total number of particles

##create objects
my_header=pygadgetic.Header()
my_body=pygadgetic.Body(npart)

my_body.pos = xlist_sorted
my_body.mass = mlist_sorted
my_body.vel = vlist_sorted

print(my_body.mass[0])
print(my_body.mass[1:])

#Checking CoM properties
print("   CoM position [pc]:", np.sum(np.atleast_2d(my_body.mass).T*my_body.pos*1e-5, axis=0)/np.sum(my_body.mass))
print("   CoM velocity [pc/kyr]:", np.sum(np.atleast_2d(my_body.mass).T*my_body.vel, axis=0)*3.24078e-14*(3600*24.0*365*1000)/np.sum(my_body.mass))

print("   ")
print("   v_rms [km/s]:", np.mean(np.sqrt(np.sum(my_body.vel**2, axis=1))))

#PBH.AddDressedPBH(my_body,np.arange(0,nDM),-1, nDM, [0, 0, 0],[0, 0, 0], r_soft, a, haloID=haloID1, verbose=True)
#print "  Second Halo:"
#PBH.AddDressedPBH(my_body,np.arange(nDM,2*nDM),-1, nDM, [-apo/2.0, 0, 0],[0, -vapo, 0], x_soft,a,  haloID=haloID2,verbose=True)

##fill in the header
my_header.NumPart_ThisFile = np.array(npart)
my_header.NumPart_Total = np.array(npart)

#id
my_body.id[:]=np.arange(0,total_number_of_particles) #generate an array from 0 to total_number_of_particles

print("  Printing to file...")
##now writes the initial condition file
#try:
my_name="../run/EMRI1.dat"
pygadgetic.dump_ic(my_header,my_body,my_name)
#except IOError:
#    my_name="run/PBH1.dat"
#    pygadgetic.dump_ic(my_header,my_body,my_name)
