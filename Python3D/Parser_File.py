# -*- coding: utf-8 -*-
"""
Let' try to use classes superclasses and whatsoever: 
    Phase_S --> Upper_crust = Phase_S(6,bla/SSB2.dat)
       short tip: i was try to evoke locals()[] to dynamically change the value of 
    a variable. However, it seems that it is considering only the local scope
    of the function. The quick hack to solve the issue was to call the class before 
    hand with the default values, using its local dictionary. But i'm suspecting 
    that there must be a clever and elegant solution that is more pythonic to
    do so.
    However, it seems the case that exec do the job, ok, of course it is a blasphemy
    in term of programming, but i need to dive yet into python before doing 
    fancy implementation
    Additional remark and for the future: it should be " Xval + space" otherwise
    it will look for different variable that feature a similar ending (i.e. eta0)
   
Organisation of the work flow: 
    a) Reference viscosity, standard condition, and reference stress. 
    -> eta0 = B tau ^(1-n) exp (H)
       if B,H are defined -> compute the eta0 for the given tau0 [spell out this approach]
    -> Reference temperature and pressure (T = base of the lithosphere, P = base of the lithosphere)
    b) If a real rheology is used: 
    -> given the reference tau
    -> given the reference T and P 
    -> compute eta0 
    -> effective compliance (i.e. B*Exp(H))
"""
import re
import numpy as np


        

class Phase_S(): 
    def __init__(self,Phase_ID,Phase_data):
        #############################################
        # Argument parsing 
        #############################################
        try:
            PD = Phase_data[Phase_ID][0]

            self.Density = self._parse_data_density(PD)
            self.Rheology = self._parse_Rheology(PD)
            self.Thermal  = self._parse_thermal(PD)
            self.Exist    = 1
        
        except:
            print("phase %d does not exist, and default value are set for each field" %Phase_ID)
            self.Density = Density()
            self.Thermal = Thermal()
            self.Rheology = Rheology()
            self.Exist    = 0 

    def _parse_data_density(self,Phase_data):
        """
        
        Parameters
        ----------
        Phase_data : str. Piece of input file that belongs to Phase = Phase_ID
            
        Read=> Parse => Organize 

        Returns Rheology,Density,Thermal 
        -------
        
        """
        # Dictionary Density 
        rx_dict = {
            'rho': re.compile(r'rho '),
            'alpha'  : re.compile(r'alpha '),
            'beta'   : re.compile(r'beta ')
            }
        key = ["rho","alpha","beta"]
        Df = Density()
        self._parse_properties(Phase_data,key,rx_dict,Df)
        return Df 


    def _parse_Rheology(self,Phase_data):
        import re
        
        """
        

        Parameters
        ----------
        Phase_data : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        Df = Rheology()
        
        
        Df.Diffusion = self._parse_Diffusion_Creep(Phase_data)
        Df.Dislocation = self._parse_Dislocation_Creep(Phase_data)
        Df.Plastic     = self._parse_Plasticity(Phase_data)
        Df.Elastic  = self._parse_Elasticity(Phase_data)
        Df.Peirls      = Peirls()
        return Df
       # Rheo = Rheology(Elastic, Diffusion, Dislocation, Plastic, WeakeningP, WeakeningS)
        
    
    def _parse_Diffusion_Creep(self,Phase_data):
        
        rx_dict = {
            'eta': re.compile(r'eta '),
            'B'  : re.compile(r'Bd '),
            'E'   : re.compile(r'Ed '),
            'V'   : re.compile(r'Vd '),
            'd' : re.compile(r'dd '),
            'diff_prof' : re.compile(r'diff_prof ')
            }
        key = ['diff_prof','eta','B','E','V']
       
        
        Df = Diffusion()
        Df=self._parse_properties(Phase_data,key,rx_dict,Df)        
        if(Df.eta != -10e23): 
            Df.B = 1/(2*Df.eta)
        return Df 
    
    def _parse_Dislocation_Creep(self,Phase_data):
        
        rx_dict = {
            'eta0': re.compile(r'eta0 '),
            'B'  : re.compile(r'Bn '),
            'E'   : re.compile(r'En '),
            'V'   : re.compile(r'Vn '),
            'n' : re.compile(r'n  '),
            'tau0': re.compile(r'tau0 '),
            'disl_prof': re.compile(r'disl_prof ')
            }
        key = ['disl_prof','eta0','B','E','V','tau0','n']
       
        
        Df = Dislocation()
        
        Df=self._parse_properties(Phase_data,key,rx_dict,Df)
        
        if(Df.eta0 != -10e23): 
            Df.Bn = (Df.tau0*1e6)**(1-Df.n)/(2*Df.eta0)
        
        return Df 
    
    
    def _parse_Plasticity(self,Phase_data):
        
        rx_dict = {
            'ch': re.compile(r'ch '),
            'fr'  : re.compile(r'fr '),
            'chSoftID' : re.compile(r'chSoftID '),
            'frSoftID' : re.compile(r'frSoftID')
            }
        key = ['ch','fr']
       
        
        Df = Plastic()
        
        Df=self._parse_properties(Phase_data,key,rx_dict,Df)
        
        # Make a separate function of it

        return Df 
    
    def _parse_thermal(self,Phase_data):
         
         rx_dict = {
             'k': re.compile(r'Cp '),
             'Cp'  : re.compile(r'k '),
             }
         key = ['k','Cp']
        
         
         Df = Thermal()
         
         Df= self._parse_properties(Phase_data,key,rx_dict,Df)
         
         # Make a separate function of it

         return Df 
     
    def _parse_Elasticity(self,Phase_data):
             
         rx_dict = {
             'G': re.compile(r'G '),
             }
         key = ['G']
            
             
         Df = Elastic()
             
         Df = self._parse_properties(Phase_data,key,rx_dict,Df)
             
        # Make a separate function of it

         return Df 
                  
                
        
    def _parse_properties(self, Phase_data,key,rx_dict,Df):
        """
        Parameters
        ----------
        Phase_data : String
            Phase data is a chunk of the LaMEM input file belonging to the current 
            phase that is analyzed
        key : String list
            For a given property there are a set of parameter that belongs to it, i.e., 
            diffusion rheology required: Ed, Vd, Bd or/and diff_profile. The key list together 
            with dictionary allow to connect key->input LaMEM structure-> Df structure
        rx_dict : Dictionary
        Df : Property class 
            Df is a sub class (Diffusion,Dislocation,Density, Thermal) containing the 
            information of a specific property (i.e. Df.Ed = X)

        Df
        -------

        """
        counter = 0  # Line Counter
        # Scan each line of the phase and find the properties
        # For each line check each of the key of the dictionary. 
        # If one of the key match with current line. Save in the respective field 
        # the number of line for the second part of the algorithm
        for line in Phase_data: 
            for k in key: 
                buf = _parse_line_input(line,k,rx_dict)
          
                if(buf != "nothing"): 
                    
                    exp  = "Df.%s = %d" %(k,counter)
                    exec(exp)
                    break
            counter +=1
        # Scan each of the property found and convert it to a number 
        # The first part of the alghoritm is simply collecting where to find the infromation
        # within the phase data base. Now, we repeat the loop over the key of the dictionary
        # to convert the information contained in the line into a data structure 
        for k in key: 
            if(eval(k,globals(),Df.__dict__)> 0): 
                line_no = (eval(k,globals(),Df.__dict__))
                if isinstance(line_no,int):
                    b_l = Phase_data[line_no]
                    if (k == 'diff_prof') | (k == 'disl_prof') | (k == 'peirl_prof'):
                        Df = self._parse_Rheological_flow_law(b_l,k,rx_dict,Df)                
                    else:
                        x=self._find_number(b_l,k)
                        exp = "Df.%s = %6f" %(k,x)
                        exec(exp)
            else: 
                x = -10e23  # If the field does not exist, place hold with nan inf
                exp = "Df.%s = %6f" %(k,x) # Execute the string
                exec(exp)
        return Df 
                
    def _find_number(self, line,k):
        import re
        split = line.split()
        for i in split :
            if re.findall(r'#',i):
                number = -10e23
                break    
            if (i != k):
                t = re.findall(r'[0-9]+',i)
            else: 
                t = []
            if len(t)>0:
                number = float(i)
                break
        return number 
    
    def _parse_Rheological_flow_law(self,line,k,rx_dict,Df):
        # To do for the next iteration: just create an other class where to store the dictionary in a reasonable way (most likely RB database)
        if k == 'diff_prof':
           
            rx_ = {
                'Diffusion_DryOlivine': re.compile(r' Dry_Olivine_diff_creep-Hirth_Kohlstedt_2003'),
                'Diffusion_WetPlagio'  : re.compile(r' Wet_Plagioclase_RybackiDresen_2000'),
                'Diffusion_WetOlivine' : re.compile(r' Wet_Olivine_diff_creep-Hirth_Kohlstedt_2003')
                }
        
        elif k == 'disl_prof':
            
            rx_ = {
                'Dislocation_DryOlivine': re.compile(r' Dry_Olivine_disl_creep-Hirth_Kohlstedt_2003'),
                'Dislocation_WetPlagio'  : re.compile(r' Wet_Plagioclase_RybackiDresen_2000'),
                'Dislocation_WetOlivine' : re.compile(r' Wet_Olivine_diff_creep-Hirth_Kohlstedt_2003')
                }
        elif k == 'peir_prof':
            
            rx_ = {
                'Peirl_creep': re.compile(r' Olivine_Peierls-Kameyama_1999'),
                }
        # Extract the key from the dictionary related to the current rheology 
        KEYS = list(rx_.keys())
        # Loop over the rheology and find which rheology law was used 
        for i in KEYS:
            
            profile = _parse_line_input(line,i,rx_)
            
            if profile != 'nothing':
                rheo_law = i 
                break
        # Create the name of the rheology to introduce in the rheology 
        rheo_name =rx_.get(rheo_law)
        
        rheo_name = rheo_name.pattern 
        
        # Extract the rheological database
        RB = Rheological_data_Base()
        RB = eval(i,globals(),RB.__dict__)
        # Extract the relevant property from the value used in the current phase
        KEYS_DB = list(Df.__dict__.keys())
        
        for ik in KEYS_DB:
            if ik == k: # If this attribute is the profile, set the name
                exp = "Df.%s = '%s'" %(ik,rheo_name)
                exec(exp) 
            if hasattr(RB, ik): # Does this attribute exist 
                
                if eval(ik,globals(),Df.__dict__) == -10e23: # If does not exist the value already (i.e. modify in the input script)
                    exp = "Df.%s = RB.%s" %(ik,ik) 
                    exec(exp)  
                    
        return Df 
        

class Density(): 
    def __init__(self):
        self.alpha = -10e23
        self.beta  = -10e23
        self.rho  = -10e23

class Rheology(): 
    def __init__(self):
        self.Elastic     = Elastic()
        self.Diffusion   = Diffusion()
        self.Dislocation = Dislocation()
        self.Plastic     = Plastic()
        self.Peirls      = Peirls()
        
        
class Elastic():
    def __init__(self):
        self.G = -10e23 
        
class Diffusion(): 
    def __init__(self):
        self.diff_prof = -10e23
        self.E    = -10e23 
        self.V     = -10e23 
        self.d = -10e23 
        self.B     = -10e23 
        self.eta    = -10e23 
        
class Dislocation(): 
    def __init__(self):
        self.disl_prof =-10e23
        self.E     = -10e23 
        self.V     = -10e23 
        self.n      = -10e23 
        self.B     = -10e23 
        self.tau0   = -10e23 
        self.eta0   = -10e23 
        
class Thermal(): 
    def __init__(self):
        self.Cp     = -10e23 
        self.k   = -10e23 
        
class Plastic(Rheology):
    def __init__(self):
        self.ch = -10e23
        self.fr = -10e23 
        self.chSoftID = 10e-23
        self.frSoftID = 10e-23


class Peirls():
    def __init__(self):
        self.peirl_prof = -10e23
        self.E = -10e23
        self.V = -10e23
        self.q  = -10e23
        self.gamma = -10e23
        self.B = -10e23        

        
def _parse_line_input(line,key,rx_dict):
    import re

    """
    Do a regex search against all defined regexes and
    return the key and match result of the first matching regex
    """
    rx = rx_dict.get(key) 
    match = rx.search(line)
    if match:
        # Sometimes it detects the commented line. Therefore it is better to remove trail character 
        # from the left, and check if it is a comment 
        line_comment = line.lstrip()
        if line_comment[0]!='#':
            return match
        else:
            return 'nothing'
    return "nothing"

def _parse_input_file(Input_File): 
    import re 
    """
    Parameters
    ----------
    input_path : "str"
    A) Read Input, find <Phase_Start> <Phase_End>
    B) Divide into N Phases 
    C) Store the chunk of information
    Return
    -------
    Class object Phases
    """
        
    rx_dict = {
        'Start_Phase': re.compile(r'<MaterialStart>'),
        'End_Phase'  : re.compile(r'<MaterialEnd>')
        }
    key2 = 'End_Phase'
    key1 = 'Start_Phase'


    counter = 0
    Phase_list = []
    d_s = -1
    d_e = -1
    for line in Input_File: 
        
        matchS = _parse_line_input(line,key1,rx_dict)
        if matchS != "nothing":
            d_s = counter
        matchE = _parse_line_input(line,key2,rx_dict)
        if matchE != "nothing":
            d_e = counter
        if ( ((d_s)>0)  & ((d_e)>0)):
            Phase_list.append([Input_File[d_s:d_e]])
            d_s = -1
            d_e = -1
            match_S = []
            match_E = []
        counter += 1 
            
    return Phase_list 
            
def _parse_geometry_slab(Input_File):
    import re 
    rx_dict = {
        'BoxS': re.compile(r'<BoxStart>'),
        'BoxE'  : re.compile(r'<BoxEnd>')
        }
    key2 = 'BoxE'
    key1 = 'BoxS'
    counter = 0
    Phase_list = []
    d_s = -1
    d_e = -1
    for line in Input_File: 
        
        matchS = _parse_line_input(line,key1,rx_dict)
        if matchS != "nothing":
            d_s = counter
        matchE = _parse_line_input(line,key2,rx_dict)
        if matchE != "nothing":
            d_e = counter
        if ( ((d_s)>0)  & ((d_e)>0)):
            Phase_list.append([Input_File[d_s:d_e]])
            d_s = -1
            d_e = -1
            match_S = []
            match_E = []
        counter += 1 
    
    rx_dict2 = {
        'bounds' : re.compile(r'bounds = ')
        }
    key = 'bounds'
    for line in Phase_list[0][0]: 
        
        matchS = _parse_line_input(line,key,rx_dict2)
        if matchS != 'nothing': 
            LL = line 
            break 
    
    spli = line.split()
    Bounds = []
    for i in spli: 
        t = re.findall(r'[0-9]+',i)
        if len(t)>0:
            Bounds.append(float(i))
    
    D0 = Bounds[1]-Bounds[0]
    L0 = Bounds[5]-Bounds[4]
    print('D0 is %6f and L0 is %6f km' %(D0,L0))
    return D0, L0 
        


class Rheological_flow_law():
    """
    Class that contains the rheological flow law parameters. 
    """
    def __init__(self,E,V,n,m,d0,B,F,MPa,r,water,q,gamma,taup):
        self.E = E
        self.V = V
        self.n = n
        self.m = m
        self.d = d0
        self.B = self._correction(B,F,n,m,MPa,d0,r,water)
        self.R = 8.314
        self.q  = q
        self.gamma = gamma
        self.taup = taup
    def _correction(self,B,F,n,m,MPa,d0,r,water):
        # Correct for accounting the typology of the experiment
        if F == 1: # Simple shear
            B = B*2**(n-1)
        if F == 2 : # Uniaxial
            B = B*(3**(n+1)/2)/2
        # Convert the unit of measure
        if MPa == 1:
            B = B*10**(-n*6); 
        # Implicitly account the water content and the grain size
        B = B*d0**(-m)*water**(r)
        return B 
        
    
    


class Rheological_data_Base():
    """
    Global data base of the rheology employed for the project 
    """
    def __init__(self):
        # Dislocation creep
        # Dry Olivine
        E = 530.0e3
        V = 15e-6
        n = 3.5 
        m = 0.0
        B = 1.1e5
        r = 1.0
        d0 = 1
        water = 1.0
        q   = -1e23
        taup = -1e23
        gamma = -1e23 
        self.Dislocation_DryOlivine = Rheological_flow_law(E,V,n,m,d0,B,1,1,r,water,q,gamma,taup)
        # Wet Olivine
        E = 520.0e3
        V = 22e-6
        n = 3.5 
        m = 0.0
        B = 1600
        r = 1.2
        d0 = 1
        water = 1000.0
        self.Dislocation_WetOlivine = Rheological_flow_law(E,V,n,m,d0,B,1,1,r,water,q,gamma,taup)
        # Wet Plagio
        E = 345.0e3
        V = 38e-6
        n = 3.0 
        m = 0.0
        B = 1.5849
        r = 1.0
        d0 = 1
        water = 158.4893
        self.Dislocation_WetPlagio  = Rheological_flow_law(E,V,n,m,d0,B,2,1,r,water,q,gamma,taup)
        # Diffusion creep 
        E = 375.0e3
        V = 5e-6
        n = 1.0 
        m = 3.0
        B = 1.5e9
        r = 1.0
        d0 = 10e3
        water = 1.0
        self.Diffusion_DryOlivine   = Rheological_flow_law(E,V,n,m,d0,B,1,1,r,water,q,gamma,taup)
        E = 375.0e3
        V = 10e-6
        n = 1.0 
        m = 3.0
        B = 2.5e7
        r = 0.8
        d0 = 10e3
        water = 1000
        self.Diffusion_WetOlivine   = Rheological_flow_law(E,V,n,m,d0,B,1,1,r,water,q,gamma,taup)
        E = 159.0e3
        V = 38e-6
        n = 1.0 
        m = 3.0
        B = 0.1995
        r = 1.0
        d0 = 100
        water = 158.4893
        self.Diffusion_WetPlagio    = Rheological_flow_law(E,V,n,m,d0,B,2,1,r,water,q,gamma,taup)

class Phase_Data_Base:
        
        def __init__(self,path):
                self.initial_input = self.__read_input(path)
                self.Phase_data    = _parse_input_file(self.initial_input)
                self.Rheo_Database: Rheological_data_Base()
                self.Phase_0_= Phase_S(0, self.Phase_data)
                self.Phase_1_= Phase_S(1, self.Phase_data)
                self.Phase_2_= Phase_S(2, self.Phase_data)
                self.Phase_3_= Phase_S(3, self.Phase_data)
                self.Phase_4_= Phase_S(4, self.Phase_data)
                self.Phase_5_= Phase_S(5, self.Phase_data)
                self.Phase_6_= Phase_S(6, self.Phase_data)
                self.Phase_7_= Phase_S(7, self.Phase_data)
                self.Phase_8_= Phase_S(8, self.Phase_data)
                self.Phase_9_= Phase_S(9, self.Phase_data)
                self.Phase_10_= Phase_S(10, self.Phase_data)
                self.Phase_11_= Phase_S(11, self.Phase_data)
                self.Phase_12_= Phase_S(12, self.Phase_data)
                self.Phase_13_= Phase_S(13, self.Phase_data)
                self.Phase_14_= Phase_S(14, self.Phase_data)
                self.Phase_15_= Phase_S(15, self.Phase_data)
                self.Phase_16_= Phase_S(16, self.Phase_data)
                self.Phase_17_= Phase_S(17, self.Phase_data)
                self.Phase_18_= Phase_S(18, self.Phase_data)
                self.Phase_19_= Phase_S(19, self.Phase_data)
                self.Phase_20_= Phase_S(20, self.Phase_data)
                self.Phase_21_= Phase_S(21, self.Phase_data)
                self.Phase_22_= Phase_S(22, self.Phase_data)
                self.Phase_23_= Phase_S(23, self.Phase_data)
                #Phase_24_: Phase_S(24, self.Phase_data)
                #Phase_25_: Phase_S(25, self.Phase_data)
                #Phase_26_: Phase_S(26, self.Phase_data)
                #Phase_27_: Phase_S(27, self.Phase_data)
                #Phase_28_: Phase_S(28, self.Phase_data)
                #Phase_29_: Phase_S(29, self.Phase_data)
                #Phase_30_: Phase_S(30, self.Phase_data)
                #Phase_31_: Phase_S(31, self.Phase_data)
                #Phase_32_: Phase_S(32, self.Phase_data)
                
            #Phase_aggregator_1_= self.Phase_Aggregator(0,self.Phase_Data)
            #Phase_aggregator_2_= self.Phase_Aggregator(1,self.Phase_Data)
            #Phase_aggregator_3_= self.Phase_Aggregator(2,self.Phase_Data)
            #Phase_aggregator_5_= self.Phase_Aggregator(3,self.Phase_Data)
            #Phase_aggregator_6_= self.Phase_Aggregator(4,self.Phase_Data)
            #WLaw_1= W_law(0,self.Phase_Data)
            #WLaw_2= W_law(1,self.Phase_Data)
            #WLaw_3= W_law(2,self.Phase_Data)

    #def W_law(id_num,)
     #place_holder
    #def Phase_Aggregator
            
        def __read_input(self,path):
            import os 
            name = os.path.join(path,"SSB_Nevena.dat")
            f = open(name, 'r')
            dd = f.readlines()
            f.close()
            return dd 
