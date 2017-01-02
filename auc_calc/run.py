import os
from subprocess import call
import time
def run(pdb_name,model_name,chain_name,contact_map_path,clustered_sequences_path,result_auc_file_name, result_zmip_path):
    start_time = time.time()
    print "auc_calculation_begins"
    #call(["julia", "auc_script.jl" ])
    call(["julia", "auc_calc/auc_script.jl",pdb_name,model_name,chain_name,contact_map_path,clustered_sequences_path,result_auc_file_name,result_zmip_path])
    print "auc_calculation_ends"
    print("--- %s seconds ---" % (time.time() - start_time))    
    
def msa_set_reference(fasta_path):
    start_time = time.time()
    print "msa_set_reference"
    #call(["julia", "auc_script.jl" ])
    call(["julia", "auc_calc/setProteinReferenceToMSA.jl",fasta_path])
    print "msa_set_reference"
    print("--- %s seconds ---" % (time.time() - start_time)) 

def auc(fasta_path,contact_map):
    start_time = time.time()
    print "auc"
    #call(["julia", "auc_script.jl" ])
    call(["julia", "auc_calc/auc.jl",fasta_path,contact_map])
    print "auc"
    print("--- %s seconds ---" % (time.time() - start_time)) 
    
def buslje09(fasta_path, zmip_result_path):
    start_time = time.time()
    print "buslje09"
    #call(["julia", "auc_script.jl" ])
    call(["julia", "auc_calc/buslje09.jl",fasta_path, zmip_result_path])
    print "buslje09"
    print("--- %s seconds ---" % (time.time() - start_time))   
        
def buslje09_(input_folder, zmip_result_path,pattern_array):
    start_time = time.time()
    print "buslje09_"
    for filename in os.listdir(input_folder):
        if filename.endswith(".cluster") & any(r in filename for r in pattern_array):
            buslje09(input_folder + filename , zmip_result_path + "zmip" + filename + ".dat")
    print "buslje09_"
    print("--- %s seconds ---" % (time.time() - start_time))
    
