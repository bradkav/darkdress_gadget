## darkdress_gadget



#### Requirements

* **Gadget2** - https://wwwmpa.mpa-garching.mpg.de/gadget/
* **pyGadgetIC** - https://github.com/ldocao/pygadgetic (put the `pygadgetic/` folder from `pyGadgetIC` into the `ICs/` folder)


#### How-to

This code allows you to generate initial conditions to feed into Gadget2, in order to simulate DM spikes around black holes, as well as BH binaries in the presence of DM spikes. 

**Initial Conditions:** You can run `ICs/GenerateICs_single.py` to generate initial conditions for a single BH surrounded by a DM spike or `ICs/GenerateICs_EMRI.py` to generate initial conditions for an IMRI/EMRI, with a DM spike around the larger central black hole. These files will output an initial conditions file called `EMRI1.dat` into the `run/` folder.

The initial conditions files rely on tabulated distribution functions for the DM spike, found in `ICs/distributions/`. I can share the notebooks I used to generate these, if necessary (they need a bit of cleaning up). 

**Units:** The IC files are set up to have lengths in units of `L = 1000*G_N/c^2 ~ 4.78e-11 pc` and masses in Solar Mass. The Gadget2 param files use these same units.

**Gadget2 param file:** The Gadget2 param files are in `run/` (take a look for example at `EMRI_HPC.param`). This specifies timestep parameters, softening lengths etc. Gagdet2 allows the possibility to 'label' different kinds of simulation particles. In our case, "Halo" is the central BH, "Disk" are the DM particles in the spike, and "Bulge" is the orbiting BH. The param file also lets you specify the initial conditions file and the output directory. 

**Running Gadget2:** You can run gadget2 (once installed) using 
```
mpiexec Gadget2 EMRI.param > output
```
or something similar. There are a number of helper scripts, such as `SubmitJob.sh` which will take a set of template param files and submit scripts, copy them into a new (target) directory and then submit the Gadget2 job to a cluster. 

**Reading the results:** Reading and intrepreting the gadget2 snapshots can be done in a few ways. You can use [`yt`](https://yt-project.org), or [`pyGadgetReader`](https://ascl.net/1411.001). Some examples of how to use `pyGadgetReader` in this context can be found [here](https://github.com/bradkav/BlackHolesDarkDress/blob/master/Nbody/PlotSeparation_single.py).

