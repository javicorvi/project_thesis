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


#Ver que es ejecucion
#2   if(atom[j].sequential < thisresidue-1 || atom[j].sequential > thisresidue+1) w=4
#3  if(atom[j].sequential < thisresidue-3 || atom[j].sequential > thisresidue+3) w=4
#4  if(atom[j].sequential < thisresidue-3 || atom[j].sequential > thisresidue+3) w=0
#5 if(atom[j].sequential < thisresidue-1 || atom[j].sequential > thisresidue+1) w=0
 
window = 3
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
'''
PDBs to evolve. 
Take each of this structures and run the all process.
'''

'''
Execute the evolution of the protein with SCPE.
Generate several families; taking account de parameters beta, nsus and runs 
'''
execute_scpe = True
beta = ["1.00"]
run = ["20000"]
nsus = ["3.0"]
'''
Execute the clustering for the families generated with SCPE to avoid redundant sequences with high identity
'''
execute_clustering = True
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
Pattern to execute process
'''
pattern=["sequences"]

'''
Iterates over the structures, pdbs and execute the scpe and the clusterization 
'''        
input_families_folder="../FAMILIES_2/"
def run_families_evol():
    start_time = time.time()
    print "run_families_evol"
    for family_folder in os.listdir(input_families_folder):
        print (family_folder)
        if(execute_family_evol):
            family_evol(input_families_folder, family_folder)
        dataanalisys.comparative_conservation(input_families_folder + family_folder)
        dataanalisys.sum_contact_map(input_families_folder + family_folder)  
    print "run_families_evol"
    print("--- %s seconds ---" % (time.time() - start_time))     
       

def family_evol(input_families_folder, family_folder):
    start_time = time.time()
    try:
        print "family_evol " + family_folder
        #todo manejar errores
        msa_gz_path=glob.glob(input_families_folder + family_folder+"/*.final.gz")
        if(len(msa_gz_path)==0):
            print ("No existe alineamiento de la familia " + family_folder)
            return
        if(len(msa_gz_path)>1):
            print ("Existe mas de un alineamiento de la familia " + family_folder)
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
        pdb_to_evol_df = util.find_pdb_to_evolve(family_pdb_information)
        for index,pdb_protein_to_evolve in pdb_to_evol_df.iterrows():
            p=pdb_protein_to_evolve['pdb']
            file_name_pdb =pdb_protein_to_evolve['seq'].replace("/","_").replace("-","_") + "_" + pdb_protein_to_evolve['pdb']+"_"+pdb_protein_to_evolve['chain'] + sufix
            complete_file_name_pdb = pdb_paths_files + file_name_pdb
            print complete_file_name_pdb
            #for pdb_gz in glob.glob(pdb_paths_files):
            #aca arrancar otro try
            #unzip pdb and move to pdb folder
            try:
                with gzip.open(complete_file_name_pdb, 'rb') as f:
                    aux_path=f.filename.split('/')
                    pdb_file_name=os.path.basename(f.filename[:-3])
                    print (pdb_file_name)
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
                #the folder to put de evol scpe sequences
                scpe_sequences=pdb_folder+"/scpe_sequences/"
                #the folder to put de evol clustered sequences
                clustered_sequences_path = pdb_folder + "/clustered_sequences/"
                #MSA information
                msa_information_path = clustered_sequences_path + "/information/"
                mi_data = pdb_folder + "/mi_data/"
                mi_data_analisys = mi_data + "info/"
                if not os.path.exists(scpe_sequences):
                    os.makedirs(scpe_sequences)
                if not os.path.exists(clustered_sequences_path):
                    os.makedirs(clustered_sequences_path)
                if not os.path.exists(msa_information_path):
                    os.makedirs(msa_information_path)    
                if not os.path.exists(mi_data):
                    os.makedirs(mi_data)
                if not os.path.exists(mi_data_analisys):
                    os.makedirs(mi_data_analisys)    
                
                    
                scpe_sequences_file=scpe_sequences+"sequences_"+pdb_file_name
                
                if(execute_scpe):
                    scpe.run(pdb_complete_path,beta,run,nsus,chain_name,scpe_sequences_file,contact_map)
                if(execute_clustering):
                    msa.clustering("0.62",scpe_sequences, clustered_sequences_path)
                    util.delete_files(scpe_sequences+'*')
                    util.delete_files(clustered_sequences_path+'*.clstr')
                if(execute_msa_information):
                        msa.msa_information_process(clustered_sequences_path, msa_information_path)
                        #util.zip_files(clustered_sequences_path+'*.cluster')
                        #util.delete_files(clustered_sequences_path+'*.cluster')
                if(execute_mi):
                    dataanalisys.buslje09_(clustered_sequences_path,mi_data)
                if(execute_dataanalisys):
                    dataanalisys.run_analisys(zmip_natural_path, mi_data, pattern, contact_map, mi_data_analisys, window)
                    #create_web_logo('../2trx_s0_w0/clustered_sequences/sequences_2trx_edit-beta1.00-nsus3.00-runs20000.fasta_0.62.cluster', 'loco.png', 'png', 'Logo')    
            except Exception as inst:
                print inst
                x = inst.args
                print x
                print "The PDB was not evolutionated " + file_name_pdb 
    except Exception as inst:
        print inst
        x = inst.args
        print x
        print "The family was not evolutionated " + family_folder 
    print "end family_evol " + family_folder
    print("--- %s seconds ---" % (time.time() - start_time))   

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