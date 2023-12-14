
from Read_VTK_files_LAMEM import * #Read_VTK_files_LAMEM
from Read_VTK_files_LAMEM import  _file_list
from Parser_File import *
from Slab_detachment import *

def _write_h5_database(ptsave,Group,TestName,S:SLAB,IC:Initial_condition,IG:Initial_Geometry,C:Coordinate_System,Fs:FS,Ph_DB:Phase_Data_Base,time):

    import h5py
    data_base_name = os.path.join(ptsave,"Data_base_Slab_detachment_3D.hdf5")
    print(data_base_name)
    f = h5py.File(data_base_name, 'a')
    e = False
    node = "/"+Group
    if node in f.keys():
            f[node]
            e = True
    if e == False: 
        f.create_group(Group)
    node_test = node+"/"+TestName
    e = False
    if node in f.keys():
            f[node]
            e = True
    if e == False: 
        f.create_group(node_test)
    node_coordinate_system = node_test+"/"+"Coordinate_System"
    # function to save the initial coordinate system 
    f = save_test_data(f,node_coordinate_system,C)
    # function to save the initial geometry
    node_IG = node_test+"/"+"IG"
    f= save_test_data(f,node_IG,IG)
    # function to save the initial Phase database
    node_PB = node_test+"/"+"Phase_DB"
    f= save_test_data(f,node_PB,Ph_DB)
    # function to save the initial condition
    node_PB = node_test+"/"+"IC"
    f= save_test_data(f,node_PB,IC)
    # function to save the slab detachment 
    node_S = node_test+"/"+"Slab_Detachment"
    f= save_test_data(f,node_S,S)
    # function to save the free surface
    node_FS = node_test+"/"+"FS"
    f= save_test_data(f,node_FS,Fs)
    buf_name =node_test+"/METADATA"
    if buf_name in f.keys():
        del f[buf_name]      # load the data
        f.create_dataset(buf_name,data = Ph_DB.initial_input)
    else:
        f.create_dataset(buf_name,data = Ph_DB.initial_input)
    buf_name = node_test+"/time"
    if buf_name in f.keys():
        del f[buf_name]      # load the data
        f.create_dataset(buf_name,data = np.array(time))
    else:
        f.create_dataset(buf_name,data = np.array(time))
    
    f.close() 

    return 1 

def save_test_data(f,path_DB,Type_DATA):
    
    classes = ['Phase_Data_Base','Initial_Geometry']
    typ = np.nan
    if isinstance(Type_DATA,Phase_Data_Base):
        typ = 1
    elif isinstance(Type_DATA,Initial_Geometry):
        typ = 2
    else:
        typ = 0 
    
    
    keys_data=Type_DATA.__dict__.keys()
    
    for v in keys_data:
                buf_name = path_DB+"/"+v
                if typ>0:
                    if typ==1:
                        if (v != 'initial_input')& (v != 'Phase_data') & (v!= 'Rheo_Database'):
                            buf_cl  = eval(v,globals(),Type_DATA.__dict__)
                            key_level1 = buf_cl.__dict__.keys()
                            if buf_cl.Exist == 1:
                                for iv1 in key_level1:
                                    buf_name1 = buf_name+"/"+iv1
                                    buf_cl1  = eval(iv1,globals(),buf_cl.__dict__)
                                    if iv1 != 'Exist':
                                        key_level2 = buf_cl1.__dict__.keys() 
                                        if iv1 == 'Rheology':
                                            for iv2 in key_level2:
                                                buf_name2 = buf_name1+"/"+iv2
                                                buf_cl2 =  eval(iv2,globals(),buf_cl1.__dict__)
                                                key_level3 = buf_cl2.__dict__.keys()
                                                for iv3 in key_level3:
                                                    buf_name3 = buf_name2+"/"+iv3
                                                    x = eval(iv3,globals(),buf_cl2.__dict__)
                                                    if isinstance(x,str):
                                                        x = x
                                                    else:
                                                        x = np.array(x)

                                                    if buf_name3 in f.keys():
                                                        del f[buf_name3]
                                                        f.create_dataset(buf_name3,data=x)
                                                    else:
                                                        f.create_dataset(buf_name3,data=x)
                                        else:
                                            for iv2 in key_level2:
                                                buf_name2 = buf_name1+"/"+iv2
                                                x = np.array(eval(iv2,globals(),buf_cl1.__dict__),dtype=float)
                                                if buf_name2 in f.keys():
                                                    del f[buf_name2]
                                                    f.create_dataset(buf_name2,data=x)
                                                else:
                                                    f.create_dataset(buf_name2,data=x)
                    if typ==2: 
                        buf_cl  = eval(v,globals(),Type_DATA.__dict__)
                        key_level1 = buf_cl.__dict__.keys()
                        for iv1 in key_level1:
                            buf_name1 = buf_name+"/"+iv1
                            buf_cl1 = eval(iv1,globals(),buf_cl.__dict__)
                            if isinstance(buf_cl1,str):
                                x = buf_cl1
                            else:
                                if isinstance(buf_cl1,object)==False:
                                    x = np.array(buf_cl1)
                                else: 
                                    x = []
                                    for i in range(len(buf_cl1)):
                                        x = np.append(x,np.array(buf_cl1[i]))
                            if buf_name1 in f.keys():
                                del f[buf_name1]
                                f.create_dataset(buf_name1,data=x)
                                    
                            else:
                                f.create_dataset(buf_name1,data=x)
                else:  
                    buf_cl  = eval(v,globals(),Type_DATA.__dict__)
                    if (v!= 'LGV' )&(v!="Label") & (v!= 'Colormap' ) & (v!= 'Val' ) & (v!= 'CV'):
                        if buf_name in f.keys():
                            del f[buf_name]      # load the data
                            f.create_dataset(buf_name,data=buf_cl)
                        else:
                            f.create_dataset(buf_name,data=buf_cl)
    return f