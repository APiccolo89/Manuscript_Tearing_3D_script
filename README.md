Initial setup and post processing routines: 
1. [Install LaMEM from the folder here]
2. Initial setup: 
a. Run LaMEM and save the grid: 
mpiexec -n128 -ParamFile Initial_setup.dat -mode save_grid
b. Run the Matlab file that is generated
c. => After the creation of the marker: 
mpiexec -n128 -ParamFile Initial_setup.dat
3. Post Processing: 
a. First I run all the experiments: 
a1: I use the script contained in Python3D to postprocess and extract the relevant data:
-> I divided the experiments in three different chunk as a function of the tau_lim 
b. I used Python_Visual_Data_Base2 to merge the data bases, and access the data of each of the test. 
