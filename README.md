# README #

##################
# Documentation  #
##################

The underlying principle of this code, as well as some studies, are documented in https://arxiv.org/abs/1707.00334  and https://journals.aps.org/prd/abstract/10.1103/PhysRevD.97.023021 . If you used this code for a publication, please cite the corresponding paper https://journals.aps.org/prd/abstract/10.1103/PhysRevD.97.023021 . A bug fix ( https://github.com/harmscho/NuTauSim/commit/ddacd1aeebab7cdfa6c155464b79d549e8b9b02e ) lead to an erratum which corresponds to the current version of the this code. The erratum you can find here https://journals.aps.org/prd/abstract/10.1103/PhysRevD.99.069902 

If you want to contribute to this project, please use the forking workflow as explained in:

https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow



##################
# HOW TO COMPILE #
##################

Method 1.
----------
In the command line run
> ./com

Method 2.
----------
In the command line type 
> make

##############
# HOW TO USE #
##############


Step 1
----------
The first step is to produce a tau-exit probability table. 
The is done with the Simu_elost executables.

Simu_elost [Energy eV] [Angle] [Number of Events] [Nu Cross Section Model] [Tau Energy Loss Model] [Water Layer Thickness] [Water Layer Density] [Output Tag] [Output Directory] [Directory of the Executable]  

Example: 
./Simu_elost 1.e20 91.0 1e2 0 0 4.0 0.92 test_tmp ./

This injects tau neutrinos with energy 10^20 eV at a zenith angle of 91 degrees (incidence angle 89 degrees, exit angle 89 degrees, emergence angle 1 degree). 100 (1e2) instances of injected tau neutrinos are simulated. The cross-section model is the middle (a.k.a. standard) curve (0) and the tau energy loss rate model is ALLM (0). The water layer thickness is 4.0 km with density 0.92 g/cm^3 (ice). The output files have tag "test_tmp" on the and they will be written to the local directory ./

One can give the energy a value of 0, in this case tau-energies are generated (randomly) uniformly in log-space in the range E=1e15 eV and E=1e21. Only in this case the tau-neutrino energies are stored in the last column in the output file. 
 
Step 2
----------
We produce a look-up table of tau-exit probabilities. 
For example generate_inputs_v3.py generates the inputs and pbs scripts to run the tau-exit probability results.
Once the runs are complete, the look-up tables are compiled with make_LUT.py and placed in the LUTs subdirectory.
NOTE: the user will have to edit the script to specify their relevant outdir and exe directories. Keep an eye out for any directories written in by hand.
The script also sets up pbs  

Step 3
----------
For a cluster running the pbs queue system, you submit via the command line
> qsub pbs_[your_script_name_here].sh

NOTE: the scripts included here are not made for general use. It is up to the user to edit them to be compatible with the computing environment.


