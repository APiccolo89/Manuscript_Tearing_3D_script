1) Install LaMEM from the folder within this repository
2) Run Experiments 
Step 1.: Run LaMEM opt using mode save_grid:
mpiexec -n128 -ParamFile Initial_setup.dat -mode save_grid
Step 2.: Run the Initial_setup_generator.m [Takes a lot, it is a matlab code]
        a. Markers and Topography
Step 3.: Run LaMEM again: 
mpiexec -n128 -ParamFile Initial_setup.dat 
in case you need to restart the run: 
mpiexec -n128 -ParamFile Initial_setup.dat -mode restart
I changed the data with a script with an innocent script
The sinking velocity is changed directly in the matlab script: 
For running again all the experiments: 
Phase 5 => Change the Vn of Phase 5
Phase 6 => Change the cohesion value to the desired value.
3) After running the experiments: Run the following script: 
a. Python3D/Slab_Break_off3D.py


b. After producing all the numerical experiments, you can use the script in Python_Visual_Data_Base2 (or if you want to cut to the chase, you can use the avaiable Database in the zenodo repository): 
"python3 Data_Base_Reader.py  "../../Data_Bases" "../../New_Results" False False False"
1st argument: (here is "../../Data_Bases") -> Location of the database 
2nd argument (here is "../../New_Results") -> Location of the folder where do you want to save the database
3rd argument: (Here is False) -> save_h5 file: it should produced a smaller database with only surface data. 
4th argument: (Here is False) -> print a txt files of the free surface data per each simulation {if you want to use for petrell or other software}
5th argument: (Here is False) -> merge database. The previous step of post-processing produces 3 different database. These database needs to be merged for using the other functionality of Data_Base_Reader. If you do for the first time the python command should be: 
"python3 Data_Base_Reader.py  "path_2_Database" "path_2_save" False False True"