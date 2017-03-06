import weblogolib
import corebio
import os
import csv


## Using the Robinson Robinson Frequencies 
amino_acid_composition = dict( 
A = 0.087135727479, R = 0.040894944605, N = 0.040418278548, D = 0.046870296325, 
C = 0.033468612677, Q = 0.038269289735, E = 0.049525516559, G = 0.088606655336, 
H = 0.033621241997, I = 0.036899088289, L = 0.085361634465, K = 0.080483022246, 
M = 0.014743987313, F = 0.039767240243, P = 0.050677889818, S = 0.069597088795, 
T = 0.058530491824, W = 0.010489421950, Y = 0.029922503029, V = 0.064717068767 ) 
aa_composition = [amino_acid_composition[_k]  for _k in weblogolib.unambiguous_protein_alphabet] 

#https://www.ncbi.nlm.nih.gov/CBBresearch/Yu/logoddslogo/apidocs/logoddslogolib.logodds-pysrc.html#info1

def create_web_logo(sequences_path, logo_output_file=None, data_output_file=None ,out_format=None, title=None, units='bits',
                    alphabet=corebio.seq.unambiguous_protein_alphabet,logo_type='SHANNON'):
    if out_format is None:
        extension = os.path.splitext(logo_output_file)[1] if logo_output_file is not None else ''
        out_format = extension[1:] if extension else 'png'
    
    fin = open(sequences_path) 
    sequences = weblogolib.read_seq_data(fin)      
    seqs = corebio.seq.SeqList([corebio.seq.Seq(s, alphabet) for s in sequences], alphabet)
    seqs.alphabet = alphabet
    
    prior=None
    if(logo_type=='SHANNON'):
        prior=None
    elif(logo_type=='KL'):
        prior = weblogolib.parse_prior('auto',alphabet)
    elif(logo_type=='ROB'):
        prior = weblogolib.parse_prior(aa_composition,alphabet)
    data = weblogolib.LogoData.from_seqs(seqs, prior)
    
    #data = weblogolib.LogoData.from_seqs(seqs)
    
    #data = weblogolib.LogoData.from_seqs(seqs,prior=aa_composition)
    options = weblogolib.LogoOptions()
    if title is not None:
        options.logo_title = title
    options.unit_name = units
    options.show_fineprint = False
    
    if out_format == 'png':
        options.resolution = 400.0
    
    format = weblogolib.LogoFormat(data, options)
    
    formatters = {
        'png': weblogolib.png_formatter,
        'svg': weblogolib.svg_formatter
    }
    
    image = formatters[out_format](data, format)
    if logo_output_file is None:
        return image
    
    with open(logo_output_file, 'wb') as fout:
        fout.write(image)
        fout.close()
    with open(data_output_file, 'wb') as fout_:
        fout_.write(str(data))  
        fout_.close()  
        
        
#create_web_logo('../cap.fa', "../cap2.png", "cap.csv", 'png')        