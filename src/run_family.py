'''
Principal script for run the evolution of a protein bla ...


'''
import glob
import gzip
import generate_table
import scpe
import urllib
import dataanalisys 
import msa 
import plot 
import util 
import os
import time
import constants as const
import pdb

import logging
logging.basicConfig(filename=const.log_file_path,level=logging.DEBUG,format="%(asctime)s - %(levelname)s - %(filename)s - %(funcName)s - %(message)s ")
consoleHandler = logging.StreamHandler()
rootLogger = logging.getLogger()
rootLogger.addHandler(consoleHandler)
#Ver que es ejecucion
#2   if(atom[j].sequential < thisresidue-1 || atom[j].sequential > thisresidue+1) w=4
#3  if(atom[j].sequential < thisresidue-3 || atom[j].sequential > thisresidue+3) w=4
#4  if(atom[j].sequential < thisresidue-3 || atom[j].sequential > thisresidue+3) w=0
#5 if(atom[j].sequential < thisresidue-1 || atom[j].sequential > thisresidue+1) w=0
 
window = 1
'''
Family Evolution
'''
execute_family_evol=True

'''
Calculate the MI for the natural MSA putting the protein as the reference
'''
execute_natural_mi_msa=False

'''
Calculate the Conservation of the families natural MSA
'''
execute_msa_natural_information=False

'''
Execute the evolution of the protein with SCPE.
Generate several families; taking account de parameters beta, nsus and runs 
'''
execute_scpe = True
beta = ["1.00"]
runs = ["20000"]
nsus = ["3.0"]
'''
Execute the clustering for the families generated with SCPE to avoid redundant sequences with high identity
'''
execute_clustering = True
'''
Execute the analisys of the MSA: Seq Logo.
'''
execute_msa_information = True
'''
Execute the MI calcultation busjle09
'''
execute_mi = True
'''
Execute the analisys of the information
'''
execute_dataanalisys = False
'''
Execute the analisys of the information between all the PDBS and MSA generated. All together
'''
execute_joined_pdb_analisys = False

'''
Recorta el msa
'''
execute_cut_msa = True


execute_download_pdbs=False
'''
Pattern to execute process
'''
pattern=["sequences"]

'''
Iterates over the structures, pdbs and execute the scpe and the clusterization 
'''        
input_families_folder="../FAMILIES_3/"
def run_families_evol():
    logging.info('Begin of the execution process')
    start_time = time.time()
    for family_folder in os.listdir(input_families_folder):
        logging.info(family_folder)
        #todo manejar errores
        
        msa_gz_path=glob.glob(input_families_folder + family_folder+"/*.final.gz")
        if(len(msa_gz_path)==0):
            logging.error('No existe alineamiento de la familia ' + family_folder)
            return
        if(len(msa_gz_path)>1):
            logging.error('Existe mas de un alineamiento de la familia ' + family_folder)
            return
        
        msa_gz_path=msa_gz_path[0]
        aux_path=msa_gz_path.split('/')
        family_pdb_information = aux_path[0]+"/"+aux_path[1]+"/"+aux_path[2] +"/"+aux_path[2]+"_pdb_level.csv" 
        pdb_to_evol_df = util.find_pdb_to_evolve(family_pdb_information)
        
        if(execute_download_pdbs):
            download_pdbs(input_families_folder,family_folder,pdb_to_evol_df)
        
        if(execute_family_evol):
            family_evol(input_families_folder, family_folder, pdb_to_evol_df)
        
        pdb_to_evol_df = validate_pdb_evolution(input_families_folder + family_folder, pdb_to_evol_df)
        if(execute_joined_pdb_analisys):
            dataanalisys.comparative_conservation(input_families_folder + family_folder, family_folder, pdb_to_evol_df)
            dataanalisys.sum_contact_map(input_families_folder + family_folder, pdb_to_evol_df)
            dataanalisys.comparative_mi_information(input_families_folder + family_folder, family_folder , 2 ,window, pdb_to_evol_df)  
    logging.info('End of the execution process')
    logging.info('Time of Execution --- %s seconds ---' % (time.time() - start_time))

def validate_pdb_evolution(family_folder, pdb_to_evol_df):
    family_folder_pdb = family_folder+"/PDB/"
    for index,pdb_protein_to_evolve in pdb_to_evol_df.iterrows():
        pdb_folder = family_folder_pdb + pdb_protein_to_evolve['pdb_folder_name']
        if(os.path.isdir(pdb_folder)):
            contact_map = pdb_folder + "/contact_map.dat"
            cmap = util.load_contact_map(contact_map)
            if(cmap.shape[0]!=69):#fix load natura size
                pdb_to_evol_df=pdb_to_evol_df.drop(index)
    return pdb_to_evol_df           
def family_evol(input_families_folder, family_folder, pdb_to_evol_df):
    start_time = time.time()
    try:
        logging.info('Family Evol ' + family_folder)
        #todo manejar errores
        msa_gz_path=glob.glob(input_families_folder + family_folder+"/*.final.gz")
        if(len(msa_gz_path)==0):
            logging.error('No existe alineamiento de la familia ' + family_folder)
            return
        if(len(msa_gz_path)>1):
            logging.error('Existe mas de un alineamiento de la familia ' + family_folder)
            return
        
        msa_gz_path=msa_gz_path[0]
        aux_path=msa_gz_path.split('/')
        msa_file_name_fasta = aux_path[0]+"/"+aux_path[1]+"/"+aux_path[2] +"/"+aux_path[2]+".fasta"    
        zmip_natural_path = msa_file_name_fasta + "_zmip.dat"
        
        family_pdb_information = aux_path[0]+"/"+aux_path[1]+"/"+aux_path[2] +"/"+aux_path[2]+"_pdb_level.csv"    
        
        if(execute_natural_mi_msa):
            msa.natural_msa_mi(msa_gz_path, msa_file_name_fasta, zmip_natural_path)
        
        #Natural Conservation Information
        if(execute_msa_natural_information):
            msa.msa_information(msa_file_name_fasta, msa_file_name_fasta, aux_path[2])
        
        #web_logo.create_web_logo(msa_file_name_fasta, msa_file_name_fasta + "_logo_sh.png",msa_file_name_fasta + "_data_sh.csv", 'png', aux_path[2], logo_type='SHANNON')
        #web_logo.create_web_logo(msa_file_name_fasta, msa_file_name_fasta + "_logo_kl.png",msa_file_name_fasta + "_data_kl.csv", 'png', aux_path[2], logo_type='KL')
        
        #pdb_paths_files = input_families_folder +  family_folder  +"/PDB/*.pdb.gz"
        pdb_paths_files = input_families_folder +  family_folder  +"/PDB/"
        
        sufix="_superimposed.pdb.gz"
        #pdb_to_evol_df = util.find_pdb_to_evolve(family_pdb_information)
        logging.info('Begin of the PDBs Evolution ')
        for index,pdb_protein_to_evolve in pdb_to_evol_df.iterrows():
            pdb_name=pdb_protein_to_evolve['pdb']
            file_name_pdb =pdb_protein_to_evolve['seq'].replace("/","_").replace("-","_") + "_" + pdb_protein_to_evolve['pdb']+"_"+pdb_protein_to_evolve['chain'] + sufix
            complete_file_name_pdb = pdb_paths_files + file_name_pdb
            logging.info('Begin of the PDB ' + file_name_pdb)
            pdb_file_complete_filename_to_evolve =  pdb_paths_files + pdb_protein_to_evolve['pdb_folder_name'] + "/"+pdb_protein_to_evolve['pdb_folder_name']+"_complete.pdb"
            
            util.remove_header(pdb_file_complete_filename_to_evolve)
            
            
            
            #for pdb_gz in glob.glob(pdb_paths_files):
            #aca arrancar otro try
            #unzip pdb and move to pdb folder
            try:
                with gzip.open(complete_file_name_pdb, 'rb') as f:
                    aux_path=f.filename.split('/')
                    pdb_file_name=os.path.basename(f.filename[:-3])
                    pdb_folder=aux_path[0]+"/"+aux_path[1]+"/"+aux_path[2]+"/"+aux_path[3]+"/"+pdb_file_name[:-17]
                    if not os.path.exists(pdb_folder):
                        os.makedirs(pdb_folder)
                    pdb_complete_path=pdb_folder+"/"+pdb_file_name    
                    pdb_file = open(pdb_complete_path ,"w")
                    file_content = f.read()
                    pdb_file.write(file_content)
                    pdb_file.close()
                #chain name to evol
                chain_name = pdb_file_name[-18:-17]
                #the contact map will be created by scpe 
                contact_map=pdb_folder+"/contact_map.dat"
                
                contact_map_syncronized = pdb_folder+"/contact_map_sync.dat"
                #the folder to put de evol scpe sequences
                scpe_sequences=pdb_folder+"/scpe_sequences/"
                #the folder to put de evol clustered sequences
                clustered_sequences_path = pdb_folder + "/clustered_sequences/"
                #the folder to put de evol clustered sequences sincronized
                sincronized_evol_path = pdb_folder + "/sincronized_evol_path/"
                #MSA information
                msa_information_path = clustered_sequences_path + "/information/"
                mi_data = pdb_folder + "/mi_data/"
                mi_data_analisys = mi_data + "info/"
                if not os.path.exists(scpe_sequences):
                    os.makedirs(scpe_sequences)
                if not os.path.exists(clustered_sequences_path):
                    os.makedirs(clustered_sequences_path)
                if not os.path.exists(sincronized_evol_path):
                    os.makedirs(sincronized_evol_path)
                if not os.path.exists(msa_information_path):
                    os.makedirs(msa_information_path)    
                if not os.path.exists(mi_data):
                    os.makedirs(mi_data)
                if not os.path.exists(mi_data_analisys):
                    os.makedirs(mi_data_analisys)    
                
                scpe_sequences_file=scpe_sequences+"sequences_"+pdb_name
                
                #return_variables_optimizated = run_methaherustic_for_optimization_parameters(pdb_file_complete_filename_to_evolve, chain_name, scpe_sequences, clustered_sequences_path,contact_map)
                
                if(execute_scpe):
                    scpe.run(pdb_file_complete_filename_to_evolve,beta,runs,nsus,chain_name,scpe_sequences_file,contact_map)
                if(execute_clustering):
                    msa.clustering("0.62",scpe_sequences, clustered_sequences_path)
                    util.delete_files(scpe_sequences+'*')
                    util.delete_files(clustered_sequences_path+'*.clstr')
                if(execute_cut_msa):
                    util.synchronize_evol_with_cutted_pdb(pdb_file_complete_filename_to_evolve, pdb_complete_path, clustered_sequences_path, sincronized_evol_path, contact_map, contact_map_syncronized)
                #if(execute_msa_information):
                        #msa.msa_information_process(clustered_sequences_path, msa_information_path)
                if(execute_msa_information):
                    msa.msa_information_process(sincronized_evol_path, msa_information_path)        
                    #util.zip_files(clustered_sequences_path+'*.cluster')
                    #util.delete_files(clustered_sequences_path+'*.cluster')
                if(execute_mi):
                    dataanalisys.buslje09_(sincronized_evol_path,mi_data)
                
                if(execute_dataanalisys):
                    dataanalisys.run_analisys(zmip_natural_path, mi_data, pattern, contact_map, mi_data_analisys, window)
                    #create_web_logo('../2trx_s0_w0/clustered_sequences/sequences_2trx_edit-beta1.00-nsus3.00-runs20000.fasta_0.62.cluster', 'loco.png', 'png', 'Logo')    
                    
                logging.info('End of the PDB ' + file_name_pdb)
            except Exception as inst:
                print inst
                x = inst.args
                print x
                logging.error('The PDB was not evolutionated ' + file_name_pdb)
    except Exception as inst:
        print inst
        x = inst.args
        print x
        logging.error('The family was not evolutionated  ' + family_folder)
    logging.info('End family Evol  ' + family_folder)
    logging.info('--- %s seconds ---' % (time.time() - start_time))   
from math import *

def run_methaherustic_for_optimization_parameters(pdb_file_complete_filename_to_evolve, chain_name, scpe_sequences, clustered_sequences_path,contact_map):
    '''import scipy.optimize as optimize
    
    def f(beta,nsus,runs):
        return beta+nsus+runs
    result=optimize.minimize(f, args=(1,1,1), method='SLSQP', bounds=((1,4),(1,10),(1,1000)))
    print(result)
    #result = optimize.minimize(f, betas, nsus,runs)
    '''
    import pandas
    df = pandas.DataFrame(columns=['BETA','NSUS','RUN','AUC'])
    beta = ["1.00"]
    runs = ["1000"]
    nsus = ["3.0"]
    for b in beta:
        for sus in nsus:
            for r in runs: 
                value = run(pdb_file_complete_filename_to_evolve,b,r,sus,chain_name, scpe_sequences, clustered_sequences_path,contact_map )
    
    
def run(pdb_file_complete_filename_to_evolve,beta,runs,nsus,chain,scpe_sequences,clustered_sequences_path,contact_map_path):    
    scpe_sequence = scpe_sequences + 'sequences_beta'+beta+'_nsus'+nsus+'_runs'+runs
    scpe.run(pdb_file_complete_filename_to_evolve, [beta], [runs],[nsus],chain,scpe_sequences,contact_map_path)
    msa.clustering("0.62",scpe_sequences, clustered_sequences_path)
    
def download_pdbs(input_families_folder, family_folder, pdb_to_evol_df):
    pdb_paths_files = input_families_folder +  family_folder  +"/PDB/"
    for index,pdb_protein_to_evolve in pdb_to_evol_df.iterrows():
        pdb_folder = pdb_paths_files + pdb_protein_to_evolve["pdb_folder_name"]
        if not os.path.exists(pdb_folder):
            os.makedirs(pdb_folder)
        pdb_data = urllib.urlopen('http://files.rcsb.org/download/'+pdb_protein_to_evolve["pdb"]+'.pdb').read()
        pdb_complete_path=pdb_folder +"/"+pdb_protein_to_evolve["pdb_folder_name"]+"_complete.pdb"
        pdb_file = open(pdb_complete_path, "w")
        pdb_file.write(pdb_data)
        pdb_file.close()


    
    
    
'''    
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
            util.sincronize_natural_evol_msas(clustered_sequences_path, curated_sequences_path,pattern,2,-3)
            dataanalisys.buslje09_(curated_sequences_path,mi_results_path, pattern)
        util.sincronize_contact_map(contact_map,contact_map+"sync",2,106)    
        dataanalisys.run_analisys(zmip_natural_result_file, mi_results_path, pattern,contact_map+"sync",mi_results_plot_path, window)    
    
    msa.conservation(clustered_sequences_path)   
    '''    
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


run_families_evol()                                   