import dataanalisys
import util
data_paths=[
            "../2trx_s3_w3/mi_data/zmipsequences_2trx_edit-beta1.00-nsus3.00-runs20000.fasta_0.62.cluster.dat_order.csv",
            "../1thx_s3_w3/mi_data/zmipsequences_1thx_edit-beta1.00-nsus3.00-runs20000.fasta_0.62.cluster.dat_order.csv"]

mi=util.load_zmip(data_paths[0])
mi2=util.load_zmip(data_paths[1])
minatural=util.load_zmip("../2trx_s3_w3/natural/zmip_PF00085_THIO_ECOLI_reference_1thx.dat")



dataanalisys.top_rank_desa(minatural,mi,mi2,1,contact_map,outputpath+filename+'top_1percent_withcon.png',filename,result_file)