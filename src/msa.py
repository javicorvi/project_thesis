'''
Created on Jan 6, 2017

@author: javi
'''
import os
from subprocess import call
from Bio.Align.AlignInfo import PSSM
import glob
import gzip
import time
import re
import dataanalisys
import web_logo
import pandas
import logging
import constants as cons
'''
Global Variables
'''
exten=".fasta"
letters = {'A':0,'R':0,'N':0,'D':0,'B':0,'C':0,'E':0,'Q':0,'Z':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0}


'''
Set the protein as reference in the msa.
This process called MIToS script to acomplish this function 
'''
def setProteinReference(fasta_path):
    start_time = time.time()
    logging.info('Begin of the method')
    #call(["julia", "auc_script.jl" ])
    call([cons.julia_exe_path, "mitos/setProteinReferenceToMSA.jl",fasta_path])
    logging.info('Time of Execution --- %s seconds ---' % (time.time() - start_time)) 
    
def clustering(clust, input_folder, output_folder, pattern_array=[".fasta"]):
    start_time = time.time()
    print "clustering_begins"
    for filename in os.listdir(input_folder):
        if filename.endswith(exten) & any(r in filename for r in pattern_array):
            print(filename)
            filenameclust = filename + "_"+clust + ".cluster"
            print(filenameclust)
            try:
                call(["cdhit", "-i" , input_folder+"/"+filename ,"-o", output_folder+"/"+filenameclust,"-c",clust,"-n", "4", "-M", "6000"])
            except Exception:
                print "The clusterization  get an exception with de pdb file " + input_folder+"/"+filename
                raise Exception ("The clusterization  get an exception with de pdb file " + input_folder+"/"+filename)    

    print "clustering_ends"
    print("--- %s seconds ---" % (time.time() - start_time))

def clustering_singular(clust, input_file, output_file):
    start_time = time.time()
    print(input_file)
    try:
        call(["cdhit", "-i" , input_file ,"-o", output_file,"-c",clust,"-n", "4", "-M", "6000"])
    except Exception:
        print "The clusterization  get an exception with de pdb file " + input_file
        raise Exception ("The clusterization  get an exception with de pdb file " + input_file)    
    print "clustering_ends"
    print("--- %s seconds ---" % (time.time() - start_time))

'''
Crea los logo shannon y kl 
'''
def msa_information_process(input_folder, output_folder, pattern_array=[".cluster"]):
    start_time = time.time()
    print "msa_information_process"
    for filename in os.listdir(input_folder):
        if filename.endswith(".cluster") & any(r in filename for r in pattern_array):
            print(filename)
            try:
                web_logo.create_web_logo(input_folder + filename, output_folder + filename + "_logo_sh.png",output_folder + filename + "_data_sh.csv", 'png', filename, logo_type='SHANNON')
                web_logo.create_web_logo(input_folder + filename, output_folder + filename + "_logo_kl.png",output_folder + filename + "_data_kl.csv", 'png', filename, logo_type='KL')
                #web_logo.create_web_logo(input_folder + filename, output_folder + filename + "_logo_rob.png",output_folder + filename + "_data_rob.csv", 'png', filename, logo_type='ROB')
            except Exception as inst:
                print inst
                print "LOGO  get an exception with the msa  " + input_folder+"/"+filename
                raise Exception ("LOGO  get an exception with the msa  " + input_folder+"/"+filename)    

    print "msa_information_process"
    print("--- %s seconds ---" % (time.time() - start_time))

def msa_information(input_msa, output_msa, msa_name):
    start_time = time.time()
    print "msa_information"
    print(msa_name)
    try:
        web_logo.create_web_logo(input_msa, output_msa + "_logo_sh.png",output_msa + "_data_sh.csv", 'png', msa_name, logo_type='SHANNON')
        web_logo.create_web_logo(input_msa, output_msa + "_logo_kl.png",output_msa + "_data_kl.csv", 'png', msa_name, logo_type='KL')
    except Exception as inst:
        print inst
        print "LOGO  get an exception with the msa  " + input_msa
        raise Exception ("LOGO  get an exception with the msa  " + input_msa)    

    print "msa_information"
    print("--- %s seconds ---" % (time.time() - start_time))
    
'''
Read the conservation info stored in file_conservation_path  
'''
def read_conservation(file_conservation_path):
    start_time = time.time()
    print "read_conservation"
    print file_conservation_path
    fields=["#","Entropy"]
    df = pandas.read_csv(file_conservation_path, delim_whitespace=True,header=7,usecols=fields)
    #df.shape[0]
    print "msa_information"
    print("--- %s seconds ---" % (time.time() - start_time))
    return df
    
def conservation(msa_path):
    import numpy as np
    import scipy.stats as sc
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from Bio.Alphabet import IUPAC
    from Bio.SubsMat import FreqTable
    import Bio.Alphabet as Alphabet
    from Bio import motifs
    for filename in os.listdir(msa_path):
        if filename.endswith(".cluster"):
            alignment = AlignIO.read(msa_path+filename, "fasta", alphabet=Alphabet.ProteinAlphabet())
            columns_quantity = []
            columns_frequency = []
            #summary_align = AlignInfo.SummaryInfo(alignment)
            #pssm = summary_align.pos_specific_score_matrix()
            #print pssm
            for x in range(0, len(alignment[0].seq)-1):
                column = alignment[:, x]
                quantity=letters
                for f in column:
                    print(f)
                    quantity[f]+=1
                double = 20/len(alignment)
                print len(alignment)
                print (quantity)
                #frequency=list(map(lambda x: x/len(alignment), quantity)) 
                frequency = dict(map(lambda (k,v): (k, v/len(alignment)), quantity.iteritems()))
                print frequency
                columns_quantity.append(quantity)
                columns_frequency.append(frequency)
            print (columns_quantity)
               
            
            '''
            m = motifs.create(alignment,alphabet=Alphabet.ProteinAlphabet())
            print (m)
            
            alfa = summary_align.alignment._alphabet
            base_alpha = Alphabet._get_base_alphabet(alfa) 
            print(summary_align)
            print(alfa)
            print(base_alpha)
            data=summary_align.information_content(5,30)
            print(data)'''
            
    #n is the number of data points
    ''''n=10
    kld = np.zeros(n, n)
    for i in range(0, n):
        for j in range(0, n):
            if(i != j):
                kld[i, j] = sc.entropy(distributions[i, :], distributions[j, :])'''

#Conver from Stockholm MSA to Fasta MSA            
def convertMSAToFasta(msa_, new_msa):
    with open(new_msa,'w') as new_file:
        with open(msa_) as old_file:
            for line in old_file:
                if('#' not in line):
                    new_line = re.split(r'\t+', line.rstrip('\t'))
                    if(len(new_line)==2):
                        new_file.write(">"+new_line[0]+"\n")
                        new_file.write(new_line[1])
    old_file.close()
    new_file.close()

# Unzip file and the convert file to fasta format and save it.  Return the path of de msa fasta format 
def natural_msa_mi(msa_gz, msa_file_name_fasta, result_zmip_path):
    try:
        with gzip.open(msa_gz, 'rb') as f:
            aux_path=f.filename.split('/')
            msa_filename=os.path.basename(f.filename)
            msa_complete_filename=aux_path[0]+"/"+aux_path[1]+"/"+aux_path[2]+"/"+msa_filename[:-3]
            msa_file = open(msa_complete_filename ,"w")
            file_content = f.read()
            msa_file.write(file_content)
            msa_file.flush()
            msa_file.close()
        #convierto el msa a formato fasta
        convertMSAToFasta(msa_complete_filename, msa_file_name_fasta)
        dataanalisys.buslje09(msa_file_name_fasta, result_zmip_path)
    except BaseException as inst:
        logging.error('Error execution MI form the natural MSA ' )
        raise Exception('Error execution MI form the natural MSA')
    
def lettercount(pos):
    return {c: pos.count(c) for c in pos}

def frequency():
    sequences = ['AATC','GCCT','ATCA']
    f = zip(*sequences)
    counts = [{letter: column.count(letter) for letter in column} for column in f]
    print(counts)
    import csv
    with open('test', "wb") as f:
        writer = csv.writer(f)
        for row in counts:
            for l in ['A','T','G','C']:
                if(row.has_key(l)):
                    print l  + str(row[l])
                else:
                    print l + "0"    
        f.close()   
#frequency()        
        