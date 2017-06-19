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

optimized_scpe_variables = True

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
input_families_folder="../FAMILIES_FIVE/"
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
            pdb_to_evol_df = pandas.read_csv(family_pdb_evol_info_path,header=0,index_col='cluster' )    
        
        if(execute_download_pdbs):
            download_pdbs(input_families_folder,family_folder,pdb_to_evol_df)
        
        if(execute_family_evol):
            family_evol(input_families_folder, family_folder, pdb_to_evol_df, family_pdb_evol_info_path)
            pdb_to_evol_df.to_csv(family_pdb_evol_info_path)
        
        

        if(execute_joined_pdb_analisys):
            pdb_to_compare=pdb_to_evol_df.loc[pdb_to_evol_df['status'] == 'okey']
            
            #dataanalisys.comparative_conservation(input_families_folder + family_folder, family_folder, pdb_to_compare)
            
            #dataanalisys.sum_contact_map(input_families_folder + family_folder, pdb_to_compare)
            
            #dataanalisys.comparative_mi_information(input_families_folder + family_folder, family_folder , 1 ,window, pdb_to_compare)  
            
            dataanalisys.compute_joined_msas(input_families_folder + family_folder ,pdb_to_compare)
        
           
            
            
    logging.info('End of the execution process')
    logging.info('Time of Execution --- %s seconds ---' % (time.time() - start_time))

        
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
        
        
        
        with gzip.open(msa_gz_path, 'rb') as f:
            aux_path=f.filename.split('/')
            msa_filename=os.path.basename(f.filename)
            msa_complete_filename_stock=aux_path[0]+"/"+aux_path[1]+"/"+aux_path[2]+"/"+msa_filename[:-3]
            msa_file = open(msa_file_name_fasta ,"w")
            file_content = f.read()
            msa_file.write(file_content)
            msa_file.flush()
            msa_file.close()
        
        
        
        msa.convertMSAToFasta(msa_complete_filename_stock, msa_file_name_fasta)
        
        
        
        
        if(execute_natural_mi_msa):
            msa.natural_msa_mi(msa_complete_filename, msa_file_name_fasta, zmip_natural_path)
        #Natural Conservation Information
        if(execute_msa_natural_information):
            msa.msa_information(msa_file_name_fasta, msa_file_name_fasta, aux_path[2])
            
        
        
        
        #web_logo.create_web_logo(msa_file_name_fasta, msa_file_name_fasta + "_logo_sh.png",msa_file_name_fasta + "_data_sh.csv", 'png', aux_path[2], logo_type='SHANNON')
        #web_logo.create_web_logo(msa_file_name_fasta, msa_file_name_fasta + "_logo_kl.png",msa_file_name_fasta + "_data_kl.csv", 'png', aux_path[2], logo_type='KL')
        
        #pdb_paths_files = input_families_folder +  family_folder  +"/PDB/*.pdb.gz"
        pdb_paths_files = input_families_folder +  family_folder  +"/PDB/"
        
        
        #optimizacion
        optimized_family = False
        optimization_folder=input_families_folder +  family_folder + "/optimization_folder/"
        optimization_file_path = optimization_folder+"optimization.csv"
        #si existe el arhivo de optimizacion entonces no hay que optimizar
        if(os.path.isfile(optimization_file_path)):
            optimized_family = True
        
        
        
        
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
            protein = pdb_protein_to_evolve['seq']
           
            
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
                
                
                #start_residue = int(pdb_protein_to_evolve['seq'][pdb_protein_to_evolve['seq'].index("/")+1:pdb_protein_to_evolve['seq'].index("-")]) 
                #end_residue = int(pdb_protein_to_evolve['seq'][pdb_protein_to_evolve['seq'].index("-")+1:]) 
                start_residue, end_residue = util.find_pdb_start_end_for_protein(msa_complete_filename_stock, pdb_protein_to_evolve['seq'], pdb_name, chain_name)
                #pdb_to_evol_df.set_value(index,"start_residue",start_residue)
                #pdb_to_evol_df.set_value(index,"end_residue",end_residue)
                    
                util.clean_pdb(pdb_file_complete,pdb_file_complete_filename_to_evolve, chain_name, start_residue, end_residue)   
                
                if(optimized_scpe_variables and not optimized_family):
                    run_optimization_parameters(protein, pdb_name,optimization_folder,optimization_file_path, pdb_file_complete_filename_to_evolve, cutted_pdb_path, chain_name, start_residue, end_residue)
                    optimized_family=True
                    
                    
                optimization_df = pandas.read_csv(optimization_file_path)
                max_scpe_params = optimization_df.ix[optimization_df['auc'].idxmax()]
                beta=max_scpe_params['beta']
                nsus=max_scpe_params['nsus']
                runs=max_scpe_params['run']
                
                
                if(pdb_protein_to_evolve['status']!='okey'):
                    pdb_to_evol_df.set_value(index,"beta",beta)
                    pdb_to_evol_df.set_value(index,"nsus",nsus)
                    pdb_to_evol_df.set_value(index,"runs",runs)
                    pdb_to_evol_df.set_value(index,"start_pdb_residue",start_residue)
                    pdb_to_evol_df.set_value(index,"end_pdb_residue",end_residue)
                    evol_protein(pdb_to_evol_df, index,pdb_file_complete_filename_to_evolve, cutted_pdb_path, beta, runs, nsus, chain_name, scpe_sequences, clustered_sequences_path, sincronized_evol_path, zmip_natural_path, mi_data, mi_data_analisys ,msa_information_path, contact_map, contact_map_syncronized)
                
                pdb_to_evol_df.set_value(index, "status","okey")
                logging.info('End of the PDB ' + file_name_pdb)
            except Exception as inst:
                print inst
                x = inst.args
                print x
                logging.error('The PDB was not evolutionated ' + file_name_pdb)
                logging.error(inst)
                pdb_to_evol_df.set_value(index, "status","error")
            pdb_to_evol_df.to_csv(family_pdb_evol_info_path)    
    except Exception as inst:
        print inst
        x = inst.args
        print x
        logging.error('The family was not evolutionated  ' + family_folder)
    logging.info('End family Evol  ' + family_folder)
    logging.info('--- %s seconds ---' % (time.time() - start_time))   


def evol_2trx():
    chain_name = "A"
    #the contact map will be created by scpe 
    pdb_folder = "../THIO_ECOLI_4_107_2TRX_A"
    
    pdb_file_complete = pdb_folder + "/THIO_ECOLI_4_107_2TRX_A_complete.pdb"
    pdb_file_complete_filename_to_evolve = pdb_folder + "/THIO_ECOLI_4_107_2TRX_A_clean.pdb" 
    optimization_folder= pdb_folder + "/optimization_folder/"
    optimization_file_path = optimization_folder+"optimization.csv"
    
    zmip_natural_result_path = pdb_folder+"/natural/"
    if not os.path.exists(zmip_natural_result_path):
        os.makedirs(zmip_natural_result_path)
    
    cutted_pdb_path = "/cutted_pdb_path/"
    zmip_natural_result_file = zmip_natural_result_path + "zmip_PF00085_THIO_ECOLI_reference.csv"
    msa_file_name_fasta = pdb_folder +  "/PF00085_THIO_ECOLI_reference.fasta"
    #setea como referencia la 2trx a la familia
    dataanalisys.buslje09(msa_file_name_fasta, zmip_natural_result_file)
    msa.msa_information(msa_file_name_fasta, msa_file_name_fasta, "PF00085")
    util.clean_pdb(pdb_file_complete,pdb_file_complete_filename_to_evolve, chain_name)   
    #run_optimization_parameters("THIO_ECOLI/4-107", "2TRX",optimization_folder,optimization_file_path, pdb_file_complete_filename_to_evolve, cutted_pdb_path, chain_name)

def analisys_2trx_evol():
    chain_name = "A"
    #the contact map will be created by scpe 
    pdb_folder = "../THIO_ECOLI_4_107_2TRX_A"
    pdb_file_complete = pdb_folder + "/THIO_ECOLI_4_107_2TRX_A_complete.pdb"
    pdb_file_complete_filename_to_evolve = pdb_folder + "/THIO_ECOLI_4_107_2TRX_A_clean.pdb" 
    optimization_folder= pdb_folder + "/optimization_folder/"
    optimization_file_path = optimization_folder+"optimization.csv"
   
    
    zmip_natural_result_path = pdb_folder+"/natural/"
    if not os.path.exists(zmip_natural_result_path):
        os.makedirs(zmip_natural_result_path)
    
    cutted_pdb_path = "/cutted_pdb_path/"
    zmip_natural_result_file = zmip_natural_result_path + "zmip_PF00085_THIO_ECOLI_reference.csv"
    msa_file_name_fasta = pdb_folder +  "/PF00085_THIO_ECOLI_reference.fasta"
    #setea como referencia la 2trx a la familia
    
    #dataanalisys.buslje09(msa_file_name_fasta, zmip_natural_result_file)
    clustered_sequences_path = optimization_folder + "clustered_sequences_path/"
    curated_sequences_path = optimization_folder + "curated_sequences_path/"
    contact_map = optimization_folder + "contact_map.dat"
    contact_map_sync = optimization_folder + "contact_map_sync.dat" 
    mi_data_path_curated = curated_sequences_path + "mi_data_path/"
    if not os.path.exists(curated_sequences_path):
        os.makedirs(curated_sequences_path) 
    conservation_path = curated_sequences_path + "conservation/"
    if not os.path.exists(conservation_path):
        os.makedirs(conservation_path) 
    if not os.path.exists(conservation_path):
        os.makedirs(conservation_path)  
    if not os.path.exists(mi_data_path_curated):
        os.makedirs(mi_data_path_curated)  
    
    util.sincronize_natural_evol_msas(clustered_sequences_path, curated_sequences_path,pattern,2,-3)
    util.sincronize_contact_map(contact_map,contact_map_sync,2,106)    
    optimization_df = pandas.read_csv(optimization_file_path,header=0,index_col=0 )
    
    for index,row_optimization in optimization_df.iterrows():
        beta=str(row_optimization['beta'])
        nsus=str(row_optimization['nsus'])
        runs=str(int(row_optimization['run']))
        sufix = "sequences-beta"+beta+"-nsus"+nsus+"-runs"+runs
        curated_seq = curated_sequences_path +sufix+".fasta.cluster"
        conservation_file = conservation_path +sufix+".fasta.cluster" 
        mi_data_file = mi_data_path_curated + "zmip_"+sufix+".csv" 
        #if(row_optimization['analysis']!='okey' and row_optimization['auc_']>=0.80):
        if(row_optimization['auc_']>=0.80):
            try: 
                dataanalisys.evol_analisys(curated_seq, mi_data_file, conservation_file, sufix)
                dataanalisys.run_analisys_singular(optimization_df, index, zmip_natural_result_file, mi_data_file, contact_map_sync, mi_data_path_curated, window)
                optimization_df.set_value(index, 'analysis', 'okey')
            except Exception as inst:
                print inst
                x = inst.args
                print x
                optimization_df.set_value(index, 'analysis', 'error')
        optimization_df.to_csv(optimization_file_path)  
        
          
def run_optimization_parameters(protein, pdb_name,optimization_folder, optimization_file_path,pdb_file_complete_filename_to_evolve, cutted_pdb_path, chain, start_residue, end_residue ):
    logging.info('Run Optimization For ' + pdb_name)
    start_time_total = time.time()
    columns=["protein","pdb","chain","beta","nsus","run","auc","auc_01","execution_time", "start_residue","end_residue"]
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
    beta = ["0.5","1.0","2.0","3.0","5.0", "7.0", "10.0","15.0","20.0"]
    runs = ["1000","5000","10000","20000"]
    nsus = ["1.0","2.0","3.0","5.0","7.0","10.0","15.0"]
    
    
    runs = ["20000"]
    
    
    
     
    auc_max = 0
    index=1
    for b in beta:
        for sus in nsus:
            for r in runs:
                start_time = time.time()
                logging.info('Calculation of beta ' + b + ' nsus ' + sus + ' run ' + r)
                df.set_value(index, 'protein', protein)
                df.set_value(index, 'pdb', pdb_name)
                df.set_value(index, 'chain', chain)
                df.set_value(index, 'start_residue', start_residue)
                df.set_value(index, 'end_residue', end_residue)
                df.set_value(index, 'beta', b)
                df.set_value(index, 'nsus', sus)
                df.set_value(index, 'run', r)
                try: 
                    auc,auc01,count_seq_scpe,count_seq_cluster,seq_lenght = evol_optimization(pdb_file_complete_filename_to_evolve, cutted_pdb_path, b,r,sus,chain, scpe_sequences, clustered_sequences_path,sincronized_evol_path, mi_data_path,contact_map,contact_map_sync )
                    df.set_value(index, 'auc', auc)
                    df.set_value(index, 'auc_01', auc01)
                    df.set_value(index, 'count_seq_scpe', count_seq_scpe)
                    df.set_value(index, 'count_seq_cluster', count_seq_cluster)
                    df.set_value(index, 'seq_lenght', seq_lenght)
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
                df.to_csv(optimization_file_path)    
    df['execution_time_optimization_total']=time.time() - start_time_total            
    df.to_csv(optimization_file_path)                
    return parameters

def run_metaheuristic_optimization_parameters(protein, pdb_name,optimization_folder, optimization_file_path,pdb_file_complete_filename_to_evolve, cutted_pdb_path, chain, start_residue, end_residue ):
    logging.info('Run Meta Heuristic Optimization For ' + pdb_name)
    start_time_total = time.time()
    columns=["protein","pdb","chain","beta","nsus","run","auc","auc_01","execution_time", "start_residue","end_residue"]
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
    beta = ["0.5","1.0","2.0","3.0","5.0", "7.0", "10.0","15.0","20.0"]
    
    nsus = ["1.0","2.0","3.0","5.0","7.0","10.0","15.0"]
    
    runs = ["20000"]
    
    max_executions = 50
    
    initial_sol = {'beta': '5.0', 'nsus': '7.0', 'runs': '20000'}
    auc,auc01,count_seq_scpe,count_seq_cluster,seq_lenght = evol_optimization(pdb_file_complete_filename_to_evolve, cutted_pdb_path, initial_sol.beta,initial_sol.runs,initial_sol.nsus,chain, scpe_sequences, clustered_sequences_path,sincronized_evol_path, mi_data_path,contact_map,contact_map_sync )
    initial_sol_value=auc01
     
    auc_max = 0
    index=1
    for s in max_executions:
        start_time = time.time()
        logging.info('Calculation of beta ' + b + ' nsus ' + sus + ' run ' + r)
        df.set_value(index, 'protein', protein)
        df.set_value(index, 'pdb', pdb_name)
        df.set_value(index, 'chain', chain)
        df.set_value(index, 'start_residue', start_residue)
        df.set_value(index, 'end_residue', end_residue)
        df.set_value(index, 'beta', initial_sol.beta)
        df.set_value(index, 'nsus', initial_sol.nsus)
        df.set_value(index, 'run', initial_sol.runs)
        try:
            find_neightbord(initial_sol)  
            auc,auc01,count_seq_scpe,count_seq_cluster,seq_lenght = evol_optimization(pdb_file_complete_filename_to_evolve, cutted_pdb_path, b,r,sus,chain, scpe_sequences, clustered_sequences_path,sincronized_evol_path, mi_data_path,contact_map,contact_map_sync )
            df.set_value(index, 'auc', auc)
            df.set_value(index, 'auc_01', auc01)
            df.set_value(index, 'count_seq_scpe', count_seq_scpe)
            df.set_value(index, 'count_seq_cluster', count_seq_cluster)
            df.set_value(index, 'seq_lenght', seq_lenght)
            if(auc01>initial_sol_value):
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
            df.to_csv(optimization_file_path)    
    df['execution_time_optimization_total']=time.time() - start_time_total            
    df.to_csv(optimization_file_path)                
    return parameters

def find_neightbord(solution):
    beta = solution.beta + 1
    nsus = solution.nsus + 1
    return -1

def evol_protein(data_frame_evol, index,pdb_file_complete_filename_to_evolve,cutted_pdb_path, beta,runs,nsus,chain,scpe_sequences,clustered_sequences_folder_path,sincronized_evol_path,zmip_natural_result_path,mi_data_path, mi_data_analisys,msa_conservation_path,contact_map_path, contact_map_sync):    
    start_time = time.time()
    beta=str(beta)
    nsus=str(nsus)
    runs=str(runs)
    sufix = "sequences-beta"+beta+"-nsus"+nsus+"-runs"+runs
    file_name = sufix +".fasta"
    output_msa_path = scpe_sequences + file_name
    clustered_sequences_path = clustered_sequences_folder_path + file_name + ".cluster"
    clustered_tmp_sequences_path = clustered_sequences_path + ".clstr"
    sincronized_evol_file_path = sincronized_evol_path + file_name 
    mi_data_path = mi_data_path + "zmip_"+sufix+".csv"
    scpe.run_singular(pdb_file_complete_filename_to_evolve, beta, runs,nsus,chain,output_msa_path,contact_map_path)
    count_seq_scpe = msa.count_sequences(output_msa_path)
    data_frame_evol.set_value(index, 'count_seq_scpe', count_seq_scpe)
    msa.clustering_singular("0.62",output_msa_path, clustered_sequences_path)
    count_seq_cluster = msa.count_sequences(clustered_sequences_path)
    data_frame_evol.set_value(index, 'count_seq_cluster', count_seq_cluster)
    seq_lenght = msa.count_aminoacids(clustered_sequences_path)
    data_frame_evol.set_value(index, 'seq_lenght', seq_lenght)
    util.delete_files(output_msa_path)
    util.delete_files(clustered_tmp_sequences_path)
    #se realiza la optimizacion sobre el msa ya recortado
    util.synchronize_evol_with_cutted_pdb_singular(pdb_file_complete_filename_to_evolve, cutted_pdb_path, clustered_sequences_path, sincronized_evol_file_path, contact_map_path, contact_map_sync)
    seq_cutted_lenght = msa.count_aminoacids(sincronized_evol_file_path)
    data_frame_evol.set_value(index, 'seq_cutted_lenght', seq_cutted_lenght)
    #util.delete_files(clustered_sequences_path)
    dataanalisys.evol_analisys(sincronized_evol_file_path, mi_data_path, msa_conservation_path + sufix, file_name)
    dataanalisys.run_analisys_singular(data_frame_evol, index, zmip_natural_result_path, mi_data_path, contact_map_sync, mi_data_analisys, window)
    data_frame_evol.set_value(index, 'execution_time', time.time() - start_time)
def evol_optimization(pdb_file_complete_filename_to_evolve,cutted_pdb_path, beta,runs,nsus,chain,scpe_sequences,clustered_sequences_folder_path,sincronized_evol_path,mi_data_path,contact_map_path, contact_map_sync):    
    sufix = "sequences-beta"+beta+"-nsus"+nsus+"-runs"+runs
    file_name = sufix +".fasta"
    output_msa_path = scpe_sequences + file_name
    clustered_sequences_path = clustered_sequences_folder_path + file_name + ".cluster"
    clustered_tmp_sequences_path = clustered_sequences_path + ".clstr"
    util.delete_files(sincronized_evol_path+"*")
    sincronized_evol_path = sincronized_evol_path + file_name 
    mi_data_path = mi_data_path + "zmip_"+sufix+".csv"
    scpe.run_singular(pdb_file_complete_filename_to_evolve, beta, runs,nsus,chain,output_msa_path,contact_map_path)
    count_seq_scpe = msa.count_sequences(output_msa_path)
    msa.clustering_singular("0.62",output_msa_path, clustered_sequences_path)
    count_seq_cluster = msa.count_sequences(clustered_sequences_path)
    seq_lenght = msa.count_aminoacids(clustered_sequences_path)
    util.delete_files(output_msa_path)
    util.delete_files(clustered_tmp_sequences_path)
    #se realiza la optimizacion sobre el msa ya recortado
    #util.synchronize_evol_with_cutted_pdb_singular(pdb_file_complete_filename_to_evolve, cutted_pdb_path, clustered_sequences_path, sincronized_evol_path, contact_map_path, contact_map_sync)
    #util.delete_files(clustered_sequences_path)
    dataanalisys.buslje09(clustered_sequences_path,mi_data_path)
    
    #target,scores=dataanalisys.getTargetScores(mi_data_path,contact_map_sync,window)
    target,scores=dataanalisys.getTargetScores(mi_data_path,contact_map_path,window)
    auc,auc01 = util.getAUC(target,scores)
    return auc,auc01,count_seq_scpe,count_seq_cluster,seq_lenght
    
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