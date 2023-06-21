import os
import pandas as pd
import csv
import xml.etree.ElementTree as ET

physicell_path = "./PhysiCell"

import argparse
 
parser = argparse.ArgumentParser(description="Building PhysiCell simulation for a particular ",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-b", "--bnd_file", help="MaBoSS BND file", required=True)
parser.add_argument("-c", "--cfg_file", help="MaBoSS CFG file", required=True)
parser.add_argument("-o", "--physicell_output_folder", help="PhysiCell output folder", required=True)
parser.add_argument("-p", "--physicell_file", help="PhysiCell XML config file", required=True)
args = parser.parse_args()

error = False
if not os.path.exists(args.bnd_file):
    error = True
    print("BND file %s does not exist !" % args.bnd_file)

if not os.path.exists(args.cfg_file):
    error = True
    print("CFG file %s does not exist !" % args.cfg_file)
    
if not os.path.exists(args.physicell_file):
    error = True
    print("Physicell XML config file %s does not exist !" % args.physicell_file)
    
if not os.path.isdir(args.physicell_output_folder):
    error = True
    print("PhysiCell output folder does not exist !" % args.physicell_output_folder)

if error:
    print("Returning...")
    exit(1)

prefix = os.path.basename(args.bnd_file[:-4])

tree = ET.parse(args.physicell_file)
root = tree.getroot()

save = root.find("save")
output_folder = save.find("folder")
output_folder.text = os.path.relpath(os.path.join(args.physicell_output_folder, prefix), physicell_path)

os.makedirs(os.path.join(args.physicell_output_folder, prefix), exist_ok=True)

bnd_entry = root.find("./cell_definitions/cell_definition/phenotype/intracellular/bnd_filename")
bnd_entry.text = str(os.path.relpath(args.bnd_file, physicell_path)) 

cfg_entry = root.find("./cell_definitions/cell_definition/phenotype/intracellular/cfg_filename")
cfg_entry.text = str(os.path.relpath(args.cfg_file, physicell_path)) 

output_xml_file = os.path.join(physicell_path, "config", f"PhysiCell_settings_cell_cycle_{prefix}.xml") 
tree.write(output_xml_file)

print("New PhysiCell XML config file has been save as %s" % output_xml_file)



