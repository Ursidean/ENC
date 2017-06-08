"""
Rule evaluation using Kappa, Kappa Simulation and clumpiness
"""


import numpy as np
import csv        
from read_map import read_map
from kappa import kappa
from kappa import ksim
from clumpy_module import clumpiness_index
from run_metro import run_metro
from set_rand import set_rand
from set_NR import set_lp_rule

'''Data_specification--------------------------------------------------------'''
#Path to base
base_path = "C:\\ENC\\"
#Land-use class names
luc_names=["Natural areas","Arable land","Permanent crops", "Pastures", 
"Agricultural areas", "Residential", "Industry & commerce", "Recreation areas", 
"Forest", "Road & rail", "Seaports", "Airports", "Mine & dump sites", 
"Fresh water", "Marine water"]
#Maximum_neighbourhood_distance_considered      
dmax=5
#Number_of_land_use_classes
luc=len(luc_names)
#Number_of_passive_land_use_classes 
pas=1
#Number_of_feature_land_use_classes
fea=6
#Number_of_active_land_use_classes
act=luc-(fea+pas)
# Geonamica command line
geo_cmd="C:\\Program Files (x86)\\Geonamica\\Metronamica\\GeonamicaCmd.exe"
# Calibration iterations testing
base_seed = 1000
max_runs = 10

# Iterative testing
case_studies = ["Berlin","Budapest","Madrid","Rome"]
for x in range(0,len(case_studies)):
    case_study = case_studies[x]
    print "Testing case study:" + case_study
    #Map data path
    data_path=base_path+"Example_case_study_data\\"
    #Original map path (time slice 1)
    map1_path=data_path+case_study+"\\"+case_study.lower()+"_1990.asc"
    #Actual map path (time slice 2)
    map2_path=data_path+case_study+"\\"+case_study.lower()+"_2000.asc"
    #Actual map path (time slice 3)
    map3_path=data_path+case_study+"\\"+case_study.lower()+"_2006.asc"
    #Mask map path 
    mask_path=data_path+case_study+"\\"+case_study.lower()+"_mask.asc"
    # Output path 
    output_path = base_path+"Example_case_study_output\\"
    '''Store_data-----------------------------------------------------------'''
    #Read in maps and mask
    map1=read_map(map1_path)
    map2=read_map(map2_path)
    map3=read_map(map3_path)
    mask=read_map(mask_path)
    '''Process_inputs-------------------------------------------------------'''
    #Determine map properties
    mapshape=np.shape(map1)         #Find the dimensions of the data
    row=mapshape[0]                 #Number of rows
    column=mapshape[1]              #Number of columns
    #Module: Partition required distances into discrete rings
    from considered_distances import considered_distances
    x=considered_distances(dmax)
    cd=x[0]                         #List of considered distances [0.00, 1.00, 1.41, 2.00 etc.]
    cdl=x[1]                        #The number of distances considered
    #Generate list of maximum possible neighbourhood sizes based on dmax
    N_all=[1,8,12,16,32,28,40,40,20]
    #Partition to size of distances considered
    N=[]
    for c in range(0,dmax):
        N.append(N_all[c])    
    '''Input_rules----------------------------------------------------------'''
    rules = {}
    for i in range(0,act):
        for j in range(0,luc):
            key="from "+luc_names[j]+" to "+luc_names[i+pas]
            rules[key]=[0,0,0,5]
    # Load the attraction rules file
    att_rule_file = output_path + case_study + "\\Rules\\att_rules.txt"
    att_rules = np.loadtxt(att_rule_file)
    # Load the rules file
    #Read inputs from csv file
    rules_file = output_path + case_study + "\\Rules\\inertia_adjusted_rules.csv"
    with open (rules_file, 'rb') as f:
        readCSV=csv.reader(f)
        next(f)                 #This skips the header line
        for row in readCSV:
            i=row[0]
            j=row[1]
            key="from "+i+" to "+j 
            rules[key][0]=row[2]
            rules[key][1]=row[3]
            rules[key][2]=row[4]
            rules[key][3]=row[5]
    # Generate results
    run_count = 0
    # Metrics
    metrics={}
    # Testing initial rules
    working_directory_2000 = ("C:\\Users\\a1210607\\"
                            "Geonamica\\Metronamica\\" + case_study + "\\")
    working_directory_2006 = ("C:\\Users\\a1210607\\"
                            "Geonamica\\Metronamica\\" + case_study + "_2006\\")                         
    
    project_file_2000 = working_directory_2000 + case_study + ".geoproj"
    project_file_2006 = working_directory_2006 + case_study + "_2006.geoproj"
    
    log_file_2000 = base_path + "LogSettings.xml"
    log_file_2006 = base_path + "LogSettings2006.xml"
    
    smap_path_2000 = ("C:\\Users\\a1210607\\Geonamica\\"
                      "Metronamica\\" + case_study + "\\Log\\Land_use\\"
                      "Land use map_2000-Jan-01 00_00_00.rst")
    smap_path_2006 = ("C:\\Users\\a1210607\\Geonamica\\"
                    "Metronamica\\" + case_study + "_2006\\Log\\Land_use\\"
                    "Land use map_2006-Jan-01 00_00_00.rst")
    # Input the rules
    for i in range(0,luc):
        for j in range(0,act):
            key="from "+luc_names[i]+" to "+luc_names[j+pas]
            fu_elem=j
            lu_elem=i
            y0=rules[key][0]
            y1=rules[key][1]
            y2=rules[key][2]
            xe=rules[key][3]
            set_lp_rule(project_file_2000,fu_elem,lu_elem,y0,y1,y2,xe)
            set_lp_rule(project_file_2006,fu_elem,lu_elem,y0,y1,y2,xe)
    while run_count < max_runs:
        # Perform analysis for year 2000
        seed = base_seed + run_count
        set_rand(project_file_2000,seed)
        run_metro(project_file_2000,log_file_2000,working_directory_2000,geo_cmd)
        smap = read_map(smap_path_2000)
        it_kappa = kappa(map2,smap,mask)
        it_ksim = ksim(map1,map2,smap,mask)
        it_clu = clumpiness_index(smap,mask,luc)
        key1="kappa-2000-it-"+str(seed)
        metrics[key1]=it_kappa
        key2="ksim-2000-it-"+str(seed)
        metrics[key2]=it_ksim
        key3="clu-2000-it-"+str(seed)
        metrics[key3]=it_clu
        # Perform analysis for year 2006
        set_rand(project_file_2006,seed)
        run_metro(project_file_2006,log_file_2006,working_directory_2006,geo_cmd)
        smap = read_map(smap_path_2006)
        it_kappa = kappa(map3,smap,mask)
        it_ksim = ksim(map2,map3,smap,mask)
        it_clu = clumpiness_index(smap,mask,luc)
        key1="kappa-2006-it-"+str(seed)
        metrics[key1]=it_kappa
        key2="ksim-2006-it-"+str(seed)
        metrics[key2]=it_ksim
        key3="clu-2006-it-"+str(seed)
        metrics[key3]=it_clu
        # Add 1 to iterator
        run_count = run_count + 1
    
    # Write output to a .csv file
    metrics_output_file = output_path + case_study + "\\final_rule_metrics.csv"
    store = [0]*(3+luc)
    with open (metrics_output_file, "wb") as csv_file:
        writer=csv.writer(csv_file)
        values=["Key","Kappa","KSIM"]+luc_names
        writer.writerow(values)
        for x in range(0,max_runs):
            seed = base_seed + x
            key1="kappa-2000-it-"+str(seed)
            key2="ksim-2000-it-"+str(seed)
            key3="clu-2000-it-"+str(seed)
            store[0]="2000-"+str(x)
            store[1]=metrics[key1]
            store[2]=metrics[key2]
            for b in range(0,luc):
                store[3+b]=metrics[key3][b]
            writer.writerow(store)
        for x in range(0,max_runs):
            seed=base_seed+x
            key1="kappa-2006-it-"+str(seed)
            key2="ksim-2006-it-"+str(seed)
            key3="clu-2006-it-"+str(seed)
            store[0]="2006-"+str(x)
            store[1]=metrics[key1]
            store[2]=metrics[key2]
            for b in range(0,luc):
                store[3+b]=metrics[key3][b]
            writer.writerow(store)

'''
import winsound
freq=2500
dur=1000
winsound.Beep(freq,dur)
'''