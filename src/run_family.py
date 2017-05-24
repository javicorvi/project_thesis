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
import pandas

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
execute_natural_mi_msa=True
'''
Calculate the Conservation of the families natural MSA
'''
execute_msa_natural_information=True

execute_download_pdbs=True

'''
Execute the evolution of the protein with SCPE.
Generate several families; taking account de parameters beta, nsus and runs 
'''
execute_scpe = False
'''
DefaultParameters
'''
beta = "1.00"
runs = "1000"
nsus = "3.0"
'''
Execute the clustering for the families generated with SCPE to avoid redundant sequences with high identity
'''
execute_clustering = False
'''
Recorta el msa
'''
execute_cut_msa = False
'''
Execute the analisys of the MSA: Seq Logo.
'''
execute_msa_information = False
'''
Execute the MI calcultation busjle09
'''
execute_mi = False
'''
Execute the analisys of the information
'''
execute_dataanalisys = False
'''
Execute the analisys of the information between all the PDBS and MSA generated. All together
'''
execute_joined_pdb_analisys = False

'''
Pattern to execute process
'''
pattern=["sequences"]

'''
Iterates over the structures, pdbs and execute the scpe and the clusterization 
'''        
input_families_folder="../FAMILIES_TO_EVOL/"
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
        family_pdb_evol_info_path = aux_path[0]+"/"+aux_path[1]+"/"+aux_path[2] +"/"+aux_path[2]+"_evol_info.csv"
        
        
        if(not os.path.isfile(family_pdb_evol_info_path)):
            pdb_to_evol_df = util.find_pdb_to_evolve(family_pdb_information)
            pdb_to_evol_df.to_csv(family_pdb_evol_info_path)
        else:
            pdb_to_evol_df = pandas.read_csv(family_pdb_evol_info_path,header=0)    
            
            
        if(execute_download_pdbs):
            download_pdbs(input_families_folder,family_folder,pdb_to_evol_df)
        
        if(execute_family_evol):
            family_evol(input_families_folder, family_folder, pdb_to_evol_df, family_pdb_evol_info_path)
        
        pdb_to_evol_df.to_csv(family_pdb_evol_info_path)
        

        #pdb_to_evol_df = validate_pdb_evolution(input_families_folder + family_folder, pdb_to_evol_df)

        
        if(execute_joined_pdb_analisys):
            dataanalisys.comparative_conservation(input_families_folder + family_folder, family_folder, pdb_to_evol_df)
            dataanalisys.sum_contact_map(input_families_folder + family_folder, pdb_to_evol_df)
            dataanalisys.comparative_mi_information(input_families_folder + family_folder, family_folder , 2 ,window, pdb_to_evol_df)  
        
        #compute_joined_msas(input_families_folder,family_folder,pdb_to_evol_df)
        
           
            
            
    logging.info('End of the execution process')
    logging.info('Time of Execution --- %s seconds ---' % (time.time() - start_time))

'''
Compute joined msas
'''
def compute_joined_msas(input_families_folder,family_folder,pdb_to_evol_df):
    joined_path = input_families_folder + family_folder + "/joined/"
    name = "joined_evol_msa"
    joined_fasta_path = joined_path + name +  ".fasta"
    if not os.path.exists(joined_path):
        os.makedirs(joined_path)
    folder_to_join = "sincronized_evol_path/"
    fasta_files=glob.glob(input_families_folder + family_folder+"/PDB/*/"+folder_to_join+"*.fasta")
    count = 0
    with open(joined_fasta_path, "w") as joined_fasta:
        for fasta in fasta_files:
            with open(fasta) as infile:
                for line in infile:
                    if('>' in line):
                        line = line.replace('\n','_'+str(count)+'\n')
                        count=count+1
                    joined_fasta.write(line)
            infile.close() 
    joined_fasta.close()
    
    
    mi_data_output_path = joined_path + name + ".csv"
    msa_conservation_path =  joined_path
    dataanalisys.evol_analisys(joined_fasta_path, mi_data_output_path, msa_conservation_path, name)
                  
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
def family_evol(input_families_folder, family_folder, pdb_to_evol_df, family_pdb_evol_info_path):
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
            pdb_file_complete =  pdb_paths_files + pdb_protein_to_evolve['pdb_folder_name'] + "/"+pdb_protein_to_evolve['pdb_folder_name']+"_complete.pdb"
            pdb_file_complete_filename_to_evolve = pdb_paths_files + pdb_protein_to_evolve['pdb_folder_name'] + "/"+pdb_protein_to_evolve['pdb_folder_name']+"_clean.pdb"
            #util.remove_header(pdb_file_complete_filename_to_evolve)
            
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
                    cutted_pdb_path=pdb_folder+"/"+pdb_file_name    
                    pdb_file = open(cutted_pdb_path ,"w")
                    file_content = f.read()
                    pdb_file.write(file_content)
                    pdb_file.close()
                #chain name to evol
                chain_name = pdb_file_name[-18:-17]
                
                
                
                #the contact map will be created by scpe 
                optimization_folder=input_families_folder +  family_folder + "/optimization_folder/"
                
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
                msa_information_path = sincronized_evol_path + "conservation/"
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
                if not os.path.exists(optimization_folder):
                    os.makedirs(optimization_folder)
                    
                util.clean_pdb(pdb_file_complete,pdb_file_complete_filename_to_evolve, chain_name)        
                
                if(index==1):
                    return_variables_optimizated = run_methaherustic_for_optimization_parameters(pdb_name,optimization_folder, pdb_file_complete_filename_to_evolve, cutted_pdb_path, chain_name)
                    beta=return_variables_optimizated[0]
                    nsus=return_variables_optimizated[1]
                    runs=return_variables_optimizated[2]
                
                  
                    #pdb_to_evol_df.set_value(index,"beta",beta)
                    #pdb_to_evol_df.set_value(index,"nsus",nsus)
                    #pdb_to_evol_df.set_value(index,"runs",runs)
                    #pdb_to_evol_df.to_csv(family_pdb_evol_info_path)
                
                #TODO 
                #Levantar maximo de optimizacion
                if(pdb_protein_to_evolve['status']!='okey'):
                    pdb_to_evol_df.set_value(index,"beta",beta)
                    pdb_to_evol_df.set_value(index,"nsus",nsus)
                    pdb_to_evol_df.set_value(index,"runs",runs)
                    evol_protein(pdb_to_evol_df, index,pdb_file_complete_filename_to_evolve, cutted_pdb_path, beta, runs, nsus, chain_name, scpe_sequences, clustered_sequences_path, sincronized_evol_path, zmip_natural_path, mi_data, mi_data_analisys ,msa_information_path, contact_map, contact_map_syncronized)

                 
                '''    
                if(execute_scpe):
                    sufix = "sequences-beta"+beta+"-nsus"+nsus+"-runs"+runs
                    file_name = sufix +".fasta"
                    output_msa_path = scpe_sequences + file_name
                    scpe.run_singular(pdb_file_complete_filename_to_evolve, beta, runs,nsus,chain_name,output_msa_path,contact_map)
                    #scpe.run(pdb_file_complete_filename_to_evolve,beta,runs,nsus,chain_name,scpe_sequences_file,contact_map)
                if(execute_clustering):
                    clustered_sequences_path = clustered_sequences_path + file_name + ".cluster"
                    clustered_tmp_sequences_path = clustered_sequences_path + file_name + ".clstr"
                    msa.clustering_singular("0.62",output_msa_path, clustered_sequences_path)
                    util.delete_files(output_msa_path)
                    util.delete_files(clustered_tmp_sequences_path)
                    #msa.clustering("0.62",scpe_sequences, clustered_sequences_path)
                    #util.delete_files(scpe_sequences+'*')
                    #util.delete_files(clustered_sequences_path+'*.clstr')
                if(execute_cut_msa):
                    sincronized_evol_path = sincronized_evol_path + file_name
                    util.synchronize_evol_with_cutted_pdb_singular(pdb_file_complete_filename_to_evolve, cutted_pdb_path, clustered_sequences_path, sincronized_evol_path, contact_map, contact_map_syncronized)
                    util.delete_files(clustered_sequences_path)
                    #util.synchronize_evol_with_cutted_pdb(pdb_file_complete_filename_to_evolve, cutted_pdb_path, clustered_sequences_path, sincronized_evol_path, contact_map, contact_map_syncronized)
                #if(execute_msa_information):
                        #msa.msa_information_process(clustered_sequences_path, msa_information_path)
                if(execute_mi):
                    #dataanalisys.buslje09_(sincronized_evol_path,mi_data)
                    mi_data_path = mi_data + "zmip_"+sufix+".csv"
                    dataanalisys.buslje09(sincronized_evol_path,mi_data_path)
                '''
                #TODO incorporar estas dos funciones para que sean singulares.
                #if(execute_msa_information):
                    #msa.msa_information_process(sincronized_evol_path, msa_information_path)        
                    #util.zip_files(clustered_sequences_path+'*.cluster')
                    #util.delete_files(clustered_sequences_path+'*.cluster')
                
                #if(execute_dataanalisys):
                    #dataanalisys.run_analisys(pdb_to_evol_df,index, zmip_natural_path, mi_data, pattern, contact_map_syncronized, mi_data_analisys, window)
                    #create_web_logo('../2trx_s0_w0/clustered_sequences/sequences_2trx_edit-beta1.00-nsus3.00-runs20000.fasta_0.62.cluster', 'loco.png', 'png', 'Logo')    
                pdb_to_evol_df.set_value(index, "status","okey")
                logging.info('End of the PDB ' + file_name_pdb)
            except Exception as inst:
                print inst
                x = inst.args
                print x
                logging.error('The PDB was not evolutionated ' + file_name_pdb)
                pdb_to_evol_df.set_value(index, "status","error")
            pdb_to_evol_df.to_csv(family_pdb_evol_info_path)    
    except Exception as inst:
        print inst
        x = inst.args
        print x
        logging.error('The family was not evolutionated  ' + family_folder)
    logging.info('End family Evol  ' + family_folder)
    logging.info('--- %s seconds ---' % (time.time() - start_time))   

def run_methaherustic_for_optimization_parameters(pdb_name,optimization_folder,pdb_file_complete_filename_to_evolve, cutted_pdb_path, chain_name):
    logging.info('Run Optimization For ' + pdb_name)
    start_time_total = time.time()
    columns=["pdb","beta","nsus","run","auc","auc_01","execution_time"]
    df = pandas.DataFrame(columns=columns)
    scpe_sequences = optimization_folder + "scpe_sequences/"
    clustered_sequences_path = optimization_folder + "clustered_sequences_path/"
    sincronized_evol_path = optimization_folder + "sincronized_evol_path/"
    mi_data_path = optimization_folder + "mi_data_path/"
    contact_map = optimization_folder + "contact_map.dat"
    contact_map_sync = optimization_folder + "contact_map_sync.dat" 
    if not os.path.exists(scpe_sequences):
        os.makedirs(scpe_sequences)
    if not os.path.exists(clustered_sequences_path):
        os.makedirs(clustered_sequences_path)
    if not os.path.exists(sincronized_evol_path):
        os.makedirs(sincronized_evol_path)
    if not os.path.exists(mi_data_path):
        os.makedirs(mi_data_path)    
    
    '''beta = ["0.5","1.00","1.5","2.00","2.5","3.00","3.5","4.0"]
    runs = ["1000","20000","30000"]
    nsus = ["1.0","2.0","3.0","4.0","5.0"]
    '''
    beta = ["0.5","1.00","1.5","2.00","2.5","3.00"]
    runs = ["20000"]
    #runs = ["1000","20000","30000"]
    nsus = ["1.0","2.0","3.0","4.0","5.0"]
    
    auc_max = 0
    index=1
    for b in beta:
        for sus in nsus:
            for r in runs:
                start_time = time.time()
                logging.info('Calculation of beta ' + b + ' nsus ' + sus + ' run ' + r)
                df.set_value(index, 'pdb', pdb_name)
                df.set_value(index, 'beta', b)
                df.set_value(index, 'nsus', sus)
                df.set_value(index, 'run', r)
                try: 
                    auc,auc01 = run(pdb_file_complete_filename_to_evolve, cutted_pdb_path, b,r,sus,chain_name, scpe_sequences, clustered_sequences_path,sincronized_evol_path, mi_data_path,contact_map,contact_map_sync )
                    df.set_value(index, 'auc', auc)
                    df.set_value(index, 'auc_01', auc01)
                    if(auc>auc_max):
                        parameters = (b,sus,r) 
                except Exception as inst:
                    print inst
                    x = inst.args
                    print x
                    df.set_value(index, 'auc', 'error')
                    df.set_value(index, 'auc_01', 'error')
                    logging.error('Error with beta ' + b + ' nsus ' + sus + ' run ' + r)        
                df.set_value(index, 'execution_time', time.time() - start_time)
                index=index+1
                df.to_csv(optimization_folder+"optimization.csv")    
    df.set_value(index, 'execution_time_optimization_total', time.time() - start_time_total)            
    df.to_csv(optimization_folder+"optimization2.csv")                
    return parameters

def evol_protein(data_frame_evol, index,pdb_file_complete_filename_to_evolve,cutted_pdb_path, beta,runs,nsus,chain,scpe_sequences,clustered_sequences_folder_path,sincronized_evol_path,zmip_natural_result_path,mi_data_path, mi_data_analisys,msa_conservation_path,contact_map_path, contact_map_sync):    
    start_time = time.time()
    sufix = "sequences-beta"+beta+"-nsus"+nsus+"-runs"+runs
    file_name = sufix +".fasta"
    output_msa_path = scpe_sequences + file_name
    clustered_sequences_path = clustered_sequences_folder_path + file_name + ".cluster"
    clustered_tmp_sequences_path = clustered_sequences_path + ".clstr"
    sincronized_evol_file_path = sincronized_evol_path + file_name 
    mi_data_path = mi_data_path + "zmip_"+sufix+".csv"
    scpe.run_singular(pdb_file_complete_filename_to_evolve, beta, runs,nsus,chain,output_msa_path,contact_map_path)
    msa.clustering_singular("0.62",output_msa_path, clustered_sequences_path)
    util.delete_files(output_msa_path)
    util.delete_files(clustered_tmp_sequences_path)
    #se realiza la optimizacion sobre el msa ya recortado
    util.synchronize_evol_with_cutted_pdb_singular(pdb_file_complete_filename_to_evolve, cutted_pdb_path, clustered_sequences_path, sincronized_evol_file_path, contact_map_path, contact_map_sync)
    util.delete_files(clustered_sequences_path)
    dataanalisys.evol_analisys(sincronized_evol_file_path, mi_data_path, msa_conservation_path + sufix, file_name)
    dataanalisys.run_analisys_singular(data_frame_evol, index, zmip_natural_result_path, mi_data_path, contact_map_sync, mi_data_analisys, window)
    data_frame_evol.set_value(index, 'execution_time', time.time() - start_time)
def run(pdb_file_complete_filename_to_evolve,cutted_pdb_path, beta,runs,nsus,chain,scpe_sequences,clustered_sequences_folder_path,sincronized_evol_path,mi_data_path,contact_map_path, contact_map_sync):    
    sufix = "sequences-beta"+beta+"-nsus"+nsus+"-runs"+runs
    file_name = sufix +".fasta"
    output_msa_path = scpe_sequences + file_name
    clustered_sequences_path = clustered_sequences_folder_path + file_name + ".cluster"
    clustered_tmp_sequences_path = clustered_sequences_folder_path + file_name + ".clstr"
    util.delete_files(sincronized_evol_path+"*")
    sincronized_evol_path = sincronized_evol_path + file_name 
    mi_data_path = mi_data_path + "zmip_"+sufix+".csv"
    scpe.run_singular(pdb_file_complete_filename_to_evolve, beta, runs,nsus,chain,output_msa_path,contact_map_path)
    msa.clustering_singular("0.62",output_msa_path, clustered_sequences_path)
    util.delete_files(output_msa_path)
    util.delete_files(clustered_tmp_sequences_path)
    #se realiza la optimizacion sobre el msa ya recortado
    #util.synchronize_evol_with_cutted_pdb_singular(pdb_file_complete_filename_to_evolve, cutted_pdb_path, clustered_sequences_path, sincronized_evol_path, contact_map_path, contact_map_sync)
    #util.delete_files(clustered_sequences_path)
    dataanalisys.buslje09(clustered_sequences_path,mi_data_path)
    
    #target,scores=dataanalisys.getTargetScores(mi_data_path,contact_map_sync,window)
    target,scores=dataanalisys.getTargetScores(mi_data_path,contact_map_path,window)
    auc,auc01 = util.getAUC(target,scores)
    return auc,auc01
    
def download_pdbs(input_families_folder, family_folder, pdb_to_evol_df):
    logging.info('download_pdbs inicio ' + family_folder)
    pdb_paths_files = input_families_folder +  family_folder  +"/PDB/"
    for index,pdb_protein_to_evolve in pdb_to_evol_df.iterrows():
        pdb_folder = pdb_paths_files + pdb_protein_to_evolve["pdb_folder_name"]
        if not os.path.exists(pdb_folder):
            os.makedirs(pdb_folder)
        logging.info('download_pdbs ' + pdb_protein_to_evolve["pdb"])    
        pdb_data = urllib.urlopen('http://files.rcsb.org/download/'+pdb_protein_to_evolve["pdb"]+'.pdb').read()
        pdb_complete_path=pdb_folder +"/"+pdb_protein_to_evolve["pdb_folder_name"]+"_complete.pdb"
        pdb_file = open(pdb_complete_path, "w")
        pdb_file.write(pdb_data)
        pdb_file.close()
    logging.info('download_pdbs fin ' + family_folder)    

    
    
    
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