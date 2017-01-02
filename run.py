from utils import generate_table as table
from utils import clustering as clust
from utils import information as info
from utils import utils as utils
from plot import auc as plot
from auc_calc import run as auc
from scpe import scpe 
import urllib


execute_scpe = False
execute_clustering = False
execute_auc = False
execute_plot = False
execute_spearman = True
execute_table_results=False
execute_sincronize_natural_evol_alignments=False
execute_natural_msa=False
execute_natural_auc=False

structures=["2trx"]
#structures=["2trx","4wxt"]

beta = ["1.00"]
run = ["20000"]
nsus = ["3.0"]
#beta = ["0.05","0.15","0.40","0.60","0.80","1.0"]
#run = ["5000","10000","15000","20000"]
#nsus = ["1.0","2.0","3.0","4.0","5.0"]

#"secuencias-beta1.00-nsus3.00-runs20000.fasta_0.62.cluster",0.8477336520262113


pattern=["sequences"]
pdbs_folder="pdbs/"
contact_map_path="contact_map/"
clustered_sequences_path="clustered_sequences/test_data/"
curated_sequences_path="curated_sequences/test_data/" 
scpe_sequences="scpe_sequences/test_data/"
result_auc_path="results/"
mi_results_path="mi_data/test_data/"
mi_results_plot_path="mi_data/plots/test_data/"
zmip_natural_result_path = "natural/zmip_PF00085_THIO_ECOLI_reference.dat"

if(execute_natural_msa):
    auc.msa_set_reference("natural/PF00085.fasta")
    auc.buslje09("natural/PF00085_THIO_ECOLI_reference.fasta",zmip_natural_result_path)        
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
        clust.run("0.62",scpe_sequences, clustered_sequences_path,pattern)
    if(execute_auc):
        auc.run(pdb_name,model_name,chain_name,contact_map,clustered_sequences_path,result_auc_file_name,result_auc_path)
    if(execute_table_results):
        table.generate_table('results/auc.dat', 'results/auc_table.dat', beta, run, nsus)
    if(execute_plot):
        plot.plot_auc('results/auc.dat','results', beta, run, nsus)
    if(execute_natural_auc):
        utils.sincronize_contact_map(contact_map,contact_map+"sync",2,106)
        auc.auc("natural/PF00085_THIO_ECOLI_reference.fasta", contact_map+"sync")
    if(execute_spearman):
        if(execute_sincronize_natural_evol_alignments):
            utils.sincronize_natural_evol_alignments(clustered_sequences_path, curated_sequences_path,pattern,2,-3)
            auc.buslje09_(curated_sequences_path,mi_results_path, pattern)
        utils.sincronize_contact_map(contact_map,contact_map+"sync",2,106)    
        info.spearman_zmip(zmip_natural_result_path, mi_results_path, pattern,contact_map+"sync",mi_results_plot_path)    
    
    #el auc del MSA evolucionado y sincronizao con la matriz de contacto sincronizado solo a lo natural 
    #tambien da 0,83 alto el auc.  
    #
    #auc.run(pdb_name,model_name,chain_name,contact_map+"sync",curated_sequences_path,result_auc_file_name,result_auc_path)
    
    '''
    natural_zmip = utils.load_zmip(mi_results_path + filename)
    utils.load_zmip(mi_results_path + filename)
    contact_map = utils.load_contact_map(contact_map+"sync")
    
    '''                                        