import web_logo
filename="sequences_2trx_edit-beta1.00-nsus3.00-runs20000.fasta_0.62.cluster"
input_folder="../2trx_s3_w3/curated_sequences/"
input_folder="../2trx_s3_w3/curated_sequences/"
web_logo.create_web_logo(input_folder + filename, input_folder + filename + "_logo_sh.png",input_folder + filename + "_data_sh.csv", 'png', filename, logo_type='SHANNON')
web_logo.create_web_logo(input_folder + filename, input_folder + filename + "_logo_kl.png",input_folder + filename + "_data_kl.csv", 'png', filename, logo_type='KL')
                