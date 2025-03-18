# Instructions for Running the Experiments

## 1) Initial Setup  
In the initial setup folder, you will find the `Initial_setup.dat` file and the MATLAB script used to generate the initial setup.

To modify the velocity of convergence (i.e., the sinking velocity in the manuscript), update the following lines in the MATLAB script:  

```matlab
Thermal_TypeTrench = Thermal_Type;
Thermal_TypeTrench.Age = convert_Age_velocity(30,1); 
Thermal_TypeTrench.vel = {[convert_Age_velocity(5,2), convert_Age_velocity(5,2)], 'none'}; 
Thermal_TypeTrench.Type = 'McKenzie';
```

The MATLAB script will also generate an `.h5` file containing all metadata related to the geometry. These data are later used by Python scripts for computing relevant scaling and generating vectors needed for referencing the tearing process.

---

## 2) Running Experiments  

### Step 1: Install LaMEM  
Install LaMEM from the folder within this repository.

### Step 2: Run Experiments  
1. **Run LaMEM in `save_grid` mode**:  
   ```sh
   mpiexec -n 128 -ParamFile Initial_setup.dat -mode save_grid
   ```
2. **Run the MATLAB script** (`Initial_setup_generator.m`):  
   - This step generates markers and topography.  
3. **Run LaMEM again**:  
   ```sh
   mpiexec -n 128 -ParamFile Initial_setup.dat
   ```
   - If you need to restart the run:  
     ```sh
     mpiexec -n 128 -ParamFile Initial_setup.dat -mode restart
     ```

---

## 3) Modifying Experimental Parameters  

- The sinking velocity is modified directly in the MATLAB script.  
- To rerun all experiments, adjust the following parameters:  
  - **Phase 5** → Change the `Vn` of Phase 5.  
  - **Phase 6** → Change the cohesion value as needed.  

### Experimental Categorization  
Experiments are grouped by `tau_lim` values:  

| Experiment Name | `tau_lim` Value |
|---------------|--------------|
| **PR_200**   | `200 MPa` |
| **PR_r**     | `400 MPa` |
| **PR_600**   | `600 MPa` |

Within each group, tests are categorized based on **sinking velocity (`v_s`)**:  

| Sinking Velocity (`v_s`) | Identifier |
|-----------------|------------|
| `5 cm/yr`  | **TSD2** |
| `10 cm/yr` | **TSD3** |
| `2.5 cm/yr` | **TSD4** |

**Important:** Run the MATLAB script each time after modifying these values. Additionally, update the activation volume of **Phase 5** as required.

### Naming Convention for Tests  
Each test should follow this naming format:  

```
TSDx_Vy/_PR/_PR2
```

Where:  
- `TSDx` → Temperature-based classification  
- `Vy` → Sinking velocity  
- `_PR` → `tau_lim = 200 MPa`  
- `_PR2` → `tau_lim = 600 MPa`  
- Absence of `_PR` → `tau_lim = 400 MPa`  

#### Example Test Names  
- `TSD2_V10_PR` → Test with `v_s = 5 cm/yr`, activation volume of `10 cm³/mol`, and `tau_lim = 200 MPa`.  
- `TSD2_V10_PR2` → Same as above but with `tau_lim = 600 MPa`.  
- `TSD2_V10` → `tau_lim = 400 MPa`.  

### Folder Organization  
1. Create three main folders: `PR_r`, `PR_200`, `PR_600`.  
2. Within each, create subfolders for each experiment following the naming convention.  
3. Run the MATLAB script.  
4. Modify `Initial_data.dat` accordingly.  
5. Run the experiments.

---

## 4) Post-Processing the Results  

### Step 1: Run the `Slab_Break_off3D.py` script  
```sh
python3 Slab_Break_off3D.py "folder PR" TestName "folder to save"
```
- **1st parameter** → Folder containing tests for a specific `tau_lim` group.  
- **2nd parameter** → Name of the test to process.  
- **3rd parameter** → Path to save results.  

### Step 2: Generate the Database  
You can either generate a new database or use the available one from the **Zenodo repository**.  

```sh
python3 Data_Base_Reader.py "../../Data_Bases" "../../New_Results" False False False
```

- **1st parameter** → Path to the database (e.g., `"../../Data_Bases"`).  
- **2nd parameter** → Folder where the processed database should be saved.  
- **3rd parameter (`False`)** → Option to save an `.h5` file (creates a smaller database with only surface data).  
- **4th parameter (`False`)** → Option to export free surface data in `.txt` format for use in external software like **Petrel**.  
- **5th parameter (`False`)** → Option to merge the database.  

If merging databases for the first time, run:  
```sh
python3 Data_Base_Reader.py "path_2_Database" "path_2_save" False False True
```

---

## Additional Notes  
- Some figures are assembled using **Inkscape**. The script provides the required parts for assembly.  
- The database includes additional test cases, such as:  
  - Simulations **without activation volume**.  
  - Simulations **without a stress limiter**.  
  - The effect of varying **convergence rates** (higher solver convergence reduces random errors on the free surface but does not significantly affect the average).  

---