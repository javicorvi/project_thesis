'''
Created on Jan 6, 2017

@author: javi
'''
import os
from subprocess import call
from Bio.Align.AlignInfo import PSSM
exten=".fasta"
import time

letters = {'A':0,'R':0,'N':0,'D':0,'B':0,'C':0,'E':0,'Q':0,'Z':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0}

def setProteinReference(fasta_path):
    start_time = time.time()
    print "msa_set_reference"
    #call(["julia", "auc_script.jl" ])
    call(["julia", "mitos/setProteinReferenceToMSA.jl",fasta_path])
    print "msa_set_reference"
    print("--- %s seconds ---" % (time.time() - start_time)) 
    
def clustering(clust, input_folder, output_folder, pattern_array=[".fasta"]):
    start_time = time.time()
    print "clustering_begins"
    for filename in os.listdir(input_folder):
        if filename.endswith(exten) & any(r in filename for r in pattern_array):
            print(filename)
            filenameclust = filename + "_"+clust + ".cluster"
            print(filenameclust)
            call(["cdhit", "-i" , input_folder+"/"+filename ,"-o", output_folder+"/"+filenameclust,"-c",clust,"-n", "4", "-M", "3000"])
    print "clustering_ends"
    print("--- %s seconds ---" % (time.time() - start_time))
    
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
            
