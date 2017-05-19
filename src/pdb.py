'''
Created on May 5, 2017
@author: javi
'''
import os
import util
'''
REMOVE PDB HEADER, SAVE ONLY THE ATOMS
'''
def remove_header(pdb_path):
    pdb_temp_path = pdb_path + '_clean'
    with open(pdb_temp_path,'w') as new_file:
        with open(pdb_path) as old_file:
            for line in old_file:
                if(line.startswith("ATOM")):  
                    new_file.write(line)     
    old_file.close()
    new_file.close()
    os.remove(pdb_path)
    os.rename(pdb_temp_path, pdb_path)
     
