#to insert in the main PhysiBoSS/PhysiCell folder
#make sure the output folder contains the sub output folder for each model
#the script creates a .csv file that can be studied easily through pandas/sumpy/matplotlib/seaborn

import pcDataLoader as pc
import pandas as pd
import os
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt

#setting general output folder in the main PhysiCell/PhysiBoSS folder
path = "./output"

#defining the dictionary with the cell cycle phases
phase_dict = {}

phase_dict[4] = "G0G1_phase"
phase_dict[10] = "S_phase"
phase_dict[11] = "G2M_phase"
phase_dict[100] = "apoptotic"


import argparse
 
parser = argparse.ArgumentParser(description="Building PhysiCell simulation for a particular ",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-f", "--folder", help="PhysiCell output folder", required=True)
args = parser.parse_args()

name_folder = args.folder
if os.path.isdir(name_folder):
    xml_files = []
    for folder, cose, files in os.walk(name_folder):
        for name in files:
            if name.endswith(".xml") and name.startswith("output"):
                xml_files.append(name)
                
    xml_files.sort()

    full_data = pd.DataFrame(columns=["time_step", "ID", "phase", "duration"])    

    for file in xml_files:
        
        mcds = pc.pyMCDS(file, name_folder)
        
        time = mcds.data["metadata"]["current_time"]
        for i in mcds.data["discrete_cells"]["data"]["ID"]:
            ID = i
            index = np.where(mcds.data["discrete_cells"]["data"]["ID"] == ID)
            phase = phase_dict[int(mcds.data["discrete_cells"]["data"]["current_phase"][index])]
            phase_duration = int(mcds.data["discrete_cells"]["data"]["elapsed_time_in_phase"][index][0])
            new_entry = pd.Series({"time_step":time, "ID":ID,
                        "phase":phase, "duration":phase_duration})
            full_data = pd.concat([full_data, new_entry.to_frame().T], ignore_index=True)

    print(full_data)
    full_data.to_csv(os.path.join(folder, "data.csv"), index=False)

