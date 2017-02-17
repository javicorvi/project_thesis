'''
Principal script for run the evolution of a protein bla ...


'''

import generate_table
import scpe
import urllib
import dataanalisys 
import msa 
import plot 
import util 
import os
#Ver que es ejecucion
#2   if(atom[j].sequential < thisresidue-1 || atom[j].sequential > thisresidue+1) w=4
#3  if(atom[j].sequential < thisresidue-3 || atom[j].sequential > thisresidue+3) w=4
#4  if(atom[j].sequential < thisresidue-3 || atom[j].sequential > thisresidue+3) w=0
#5 if(atom[j].sequential < thisresidue-1 || atom[j].sequential > thisresidue+1) w=0
execution_name="1thx_s3_w3"
window = 3
'''
Calculate the MI for the natural MSA putting the protein as the reference
'''
execute_natural_mi_msa=False
'''
PDBs to evolve. 
Take each of this structures and run the all process.
'''
structures=["1thx_edit"]
'''
Execute the evolution of the protein with SCPE.
Generate several families; taking account de parameters beta, nsus and runs 
'''
execute_scpe = False
beta = ["1.00"]
run = ["20000"]
nsus = ["3.0"]
'''
Execute the clustering for the families generated with SCPE to avoid redundant sequences with high identity
'''
execute_clustering = False
'''
Execute the AUC for all the sequences generated and clustered
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
'''
Run calculations of MI with contact maps. Top Ranks. Spearman and more.
All the plots where in the data_path .. /plots
'''
data_analisys = True

execute_sincronize_natural_evol_msas=False



data_path="../"+execution_name+"/"
mitos_scripts_path=""



#structures=["2trx","4wxt"]


#beta = ["0.05","0.15","0.40","0.60","0.80","1.0"]
#run = ["5000","10000","15000","20000"]
#nsus = ["1.0","2.0","3.0","4.0","5.0"]

#"secuencias-beta1.00-nsus3.00-runs20000.fasta_0.62.cluster",0.8477336520262113


pattern=["sequences"]
pdbs_folder=data_path+"pdbs/"
contact_map_path=data_path+"contact_map/"
clustered_sequences_path=data_path+"clustered_sequences/"
curated_sequences_path=data_path+"curated_sequences/" 
scpe_sequences=data_path+"scpe_sequences/"
result_auc_path=data_path+"results/"
mi_results_path=data_path+"mi_data/"
mi_results_plot_path=data_path+"mi_data/plots/"
zmip_natural_result_path = data_path+"natural/"
zmip_natural_result_file = zmip_natural_result_path + "zmip_PF00085_THIO_ECOLI_reference_1thx.dat"

if not os.path.exists(data_path):
    os.makedirs(data_path)
if not os.path.exists(pdbs_folder):
    os.makedirs(pdbs_folder)
if not os.path.exists(contact_map_path):
    os.makedirs(contact_map_path)
if not os.path.exists(clustered_sequences_path):
    os.makedirs(clustered_sequences_path)
if not os.path.exists(curated_sequences_path):
    os.makedirs(curated_sequences_path)
if not os.path.exists(scpe_sequences):
    os.makedirs(scpe_sequences)
if not os.path.exists(mi_results_path):
    os.makedirs(mi_results_path)
if not os.path.exists(mi_results_plot_path):
    os.makedirs(mi_results_plot_path)
if not os.path.exists(zmip_natural_result_path):
    os.makedirs(zmip_natural_result_path)    
if not os.path.exists(result_auc_path):
    os.makedirs(result_auc_path)  
    
if(execute_natural_mi_msa):
    #msa.setProteinReference("../data/natural/PF00085.fasta")
    dataanalisys.buslje09("../data/natural/PF00085_THIO_ECOLI_reference_1thx.fst",zmip_natural_result_file)    
    
'''
Iterates over the structures, pdbs and execute the all process 
'''        
for s in structures:
    pdb_name=s
    #pdb_data = urllib.urlopen('http://files.rcsb.org/download/'+pdb_name+'.pdb').read()
    file_pdb = open("../pdbs/1thx_edit.pdb")
    pdb_data = file_pdb.read()
    #download de fasta sino esta , idem pdb
    #pdb_complete_path=pdbs_folder+pdb_name+".pdb"
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
        dataanalisys.auc("../data/natural/PF00085_THIO_ECOLI_reference.fasta", contact_map+"sync")
    if(data_analisys):
        if(execute_sincronize_natural_evol_msas):
            util.sincronize_natural_evol_msas(clustered_sequences_path, curated_sequences_path,pattern,0,-3)
            dataanalisys.buslje09_(curated_sequences_path,mi_results_path, pattern)
        util.sincronize_contact_map(contact_map,contact_map+"sync",0,103)    
        dataanalisys.run_analisys(zmip_natural_result_file, mi_results_path, pattern,contact_map+"sync",mi_results_plot_path, window)    
    
    #msa.conservation(clustered_sequences_path)   
        
    #el auc del MSA evolucionado y sincronizao con la matriz de contacto sincronizado solo a lo natural 
    #tambien da 0,83 alto el auc.  
    #
    #auc.run(pdb_name,model_name,chain_name,contact_map+"sync",curated_sequences_path,result_auc_file_name,result_auc_path)
    
    #http://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html#sphx-glr-auto-examples-model-selection-plot-roc-py
    
    '''
    natural_zmip = utils.load_zmip(mi_results_path + filename)
    utils.load_zmip(mi_results_path + filename)
    contact_map = utils.load_contact_map(contact_map+"sync")
    
    '''                                        