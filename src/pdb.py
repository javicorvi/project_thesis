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
     
# The MIT License
# 
# Copyright (c) 2010-2016 Anders S. Christensen
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import Bio.PDB

# Select what residues numbers you wish to align
# and put them in a list

def align_pdb():
    start_id = 1
    end_id   = 110
    atoms_to_be_aligned = range(start_id, end_id + 1)
    
    # Start the parser
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    
    # Get the structures
    ref_structure = pdb_parser.get_structure("reference", "2TRX.pdb")
    sample_structure = pdb_parser.get_structure("samle", "2H6X.pdb")
    
    # Use the first model in the pdb-files for alignment
    # Change the number 0 if you want to align to another structure
    ref_model    = ref_structure[0]
    sample_model = sample_structure[0]
    
    # Make a list of the atoms (in the structures) you wish to align.
    # In this case we use CA atoms whose index is in the specified range
    ref_atoms = []
    sample_atoms = []
    
    # Iterate of all chains in the model in order to find all residues
    for ref_chain in ref_model:
        # Iterate of all residues in each model in order to find proper atoms
        for ref_res in ref_chain:
        # Check if residue number ( .get_id() ) is in the list
            if ref_res.get_id()[1] in atoms_to_be_aligned:
                # Append CA atom to list
                ref_atoms.append(ref_res['CA'])
    
    # Do the same for the sample structure
    for sample_chain in sample_model:
        for sample_res in sample_chain:
            if sample_res.get_id()[1] in atoms_to_be_aligned:
                sample_atoms.append(sample_res['CA'])
    
    # Now we initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
    
    # Print RMSD:
    print super_imposer.rms
    
    # Save the aligned version of 1UBQ.pdb
    io = Bio.PDB.PDBIO()
    io.set_structure(sample_structure) 
    io.save("1UBQ_aligned.pdb")