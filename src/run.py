#from utils import generate_table as table
#from utils import clustering as clust


#from plot import auc as plot
#from auc_calc import run as auc
#from scpe import scpe 

import generate_table
import scpe2 as scpe
import urllib
import dataanalisys 
import msa 
import plot2 as plot
import util 

'''
Calculate the MI for the natural MSA putting the protein as the reference
'''
execute_natural_mi_msa=False
'''
PDBs to evolve. 
Take each of this structures and run the all process.
'''
structures=["2trx"]
'''
Execute the evolution of the protein with SCPE.
Generate several families; taking account de parameters beta, nsus and runs 
'''
execute_scpe = False
beta = ["1.00"]
run = ["5000"]
nsus = ["1.0"]
'''
Execute the clustering for the families generated with SCPE to avoid redundant sequences with high identity
'''
execute_clustering = False
'''
Execute the AUC for all the sequences genrated and clustered
'''
execute_auc_process = False
'''
Generate the table with the results as a txt for documentation
'''
execute_table_results=False
'''
Generate the plot for the diferents auc taking into account the beta, nsus and runs
'''
execute_auc_plot = False
'''
Calculates the natural AUC with the protein as reference, with the contact map generates by scpe
'''
execute_natural_auc=False






execute_spearman = True

execute_sincronize_natural_evol_msas=False



data_path="../data/"
mitos_scripts_path=""


#structures=["2trx","4wxt"]


#beta = ["0.05","0.15","0.40","0.60","0.80","1.0"]
#run = ["5000","10000","15000","20000"]
#nsus = ["1.0","2.0","3.0","4.0","5.0"]

#"secuencias-beta1.00-nsus3.00-runs20000.fasta_0.62.cluster",0.8477336520262113


pattern=["sequences"]
pdbs_folder=data_path+"pdbs/"
contact_map_path=data_path+"contact_map/"
clustered_sequences_path=data_path+"clustered_sequences/test_data/"
curated_sequences_path=data_path+"curated_sequences/test_data/" 
scpe_sequences=data_path+"scpe_sequences/test_data/"
result_auc_path=data_path+"results/"
mi_results_path=data_path+"mi_data/test_data/"
mi_results_plot_path=data_path+"mi_data/plots/test_data/"
zmip_natural_result_path = data_path+"natural/zmip_PF00085_THIO_ECOLI_reference.dat"


if(execute_natural_mi_msa):
    msa.setProteinReference(data_path+"natural/PF00085.fasta")
    dataanalisys.buslje09(data_path+"natural/PF00085_THIO_ECOLI_reference.fasta",zmip_natural_result_path)    
    
'''
Iterates over the structures, pdbs and execute the all process 
'''        
for s in structures:
    pdb_name=s
    pdb_data = urllib.urlopen('http://files.rcsb.org/download/'+pdb_name+'.pdb').read()
    #download de fasta sino esta , idem pdb
    pdb_complete_path=pdbs_folder+pdb_name+".pdb"
    pdb_file = open(pdb_complete_path, "w")
    pdb_file.write(pdb_data)
    pdb_file.close()
    model_name="1"
    chain_name="A"
    contact_map=contact_map_path+"contact_map_" + pdb_name + ".dat"
    result_auc_file_name=result_auc_path+"result_auc_"+pdb_name+".dat"
    #'results/auc.dat'
    scpe_sequences_file=scpe_sequences+"sequences_"+pdb_name
    if(execute_scpe):
        scpe.run(pdb_complete_path,beta,run,nsus,chain_name,scpe_sequences_file,contact_map)
    if(execute_clustering):
        msa.clustering("0.62",scpe_sequences, clustered_sequences_path,pattern)
    if(execute_auc_process):
        dataanalisys.auc_job(pdb_name,model_name,chain_name,contact_map,clustered_sequences_path,result_auc_file_name,result_auc_path)
    
    
    if(execute_table_results):
        generate_table.run_job(data_path+'results/auc.dat', data_path+'results/auc_table.dat', ["0.05","0.15","0.40","0.60","0.80","1.0"],["5000","10000","15000","20000"], ["1.0","2.0","3.0","4.0","5.0"])
    if(execute_auc_plot):
        plot.plot_auc(data_path+'results/auc.dat',data_path+'results', ["0.05","0.15","0.40","0.60","0.80","1.0"],  ["5000","10000","15000","20000"], ["1.0","2.0","3.0","4.0","5.0"])
    
    if(execute_natural_auc):
        util.sincronize_contact_map(contact_map,contact_map+"sync",2,106)#todo automatizar la sincronizacion
        dataanalisys.auc(data_path+"natural/PF00085_THIO_ECOLI_reference.fasta", contact_map+"sync")
    
    if(execute_spearman):
        if(execute_sincronize_natural_evol_msas):
            util.sincronize_natural_evol_msas(clustered_sequences_path, curated_sequences_path,pattern,2,-3)
            dataanalisys.buslje09_(curated_sequences_path,mi_results_path, pattern)
        util.sincronize_contact_map(contact_map,contact_map+"sync",2,106)    
        dataanalisys.run_analisys(zmip_natural_result_path, mi_results_path, pattern,contact_map+"sync",mi_results_plot_path)    
    
    #el auc del MSA evolucionado y sincronizao con la matriz de contacto sincronizado solo a lo natural 
    #tambien da 0,83 alto el auc.  
    #
    #auc.run(pdb_name,model_name,chain_name,contact_map+"sync",curated_sequences_path,result_auc_file_name,result_auc_path)
    
    '''
    natural_zmip = utils.load_zmip(mi_results_path + filename)
    utils.load_zmip(mi_results_path + filename)
    contact_map = utils.load_contact_map(contact_map+"sync")
    
    '''                                        