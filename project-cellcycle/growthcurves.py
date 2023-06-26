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
    pop_sizes = {}
    for file in xml_files:
        
        mcds = pc.pyMCDS(file, name_folder)
        time = mcds.data["metadata"]["current_time"]
        df = mcds.get_cell_df()
        # print(df)
        pop_sizes.update({time: df[df["current_phase"] != 100].shape[0]})
        
        # for cell in mcds.data["discrete_cells"]["data"]:
        #     print(cell)
        #     print(mcds.data["discrete_cells"]["data"]["current_phase"][cell])
        # pop_size = len([cell for cell in mcds.data["discrete_cells"] if cell["data"]["current_phase"] != 100])
    
        # pop_sizes.update({time: pop_size})
        
    lists = sorted(pop_sizes.items())
    x, y = zip(*lists)
    plt.plot(x,y)
plt.show()
