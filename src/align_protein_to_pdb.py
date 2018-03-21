

from collections import defaultdict
import re, os, sys
from Bio import pairwise2
from Bio.Align.Applications import MuscleCommandline
import string
import random
from Bio import SeqIO
import multiprocessing as mp


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def searchCDsequence(cdstopath, seqid):
    with open(cdstopath, 'r') as f:
        for line in f:
            if re.search(seqid, line):
                actual_seqid, seq = line.split("\t")
                colmaps = []
                for i in range(1,len(seq)):
                    if seq[i-1] == "." or seq[i-1] == "-":
                        continue
                    colmaps.append(i)
                seq = re.sub(r'\.', "", seq)
                seq = re.sub("-", "", seq)
                # seq = re.sub(r'\s', "", seq)
                return(seq, colmaps)

def MusclePairAlign(id1, seq1, id2, seq2):
    tmpname = id_generator()
    if not os.path.exists(os.path.join("tmp")):
        os.makedirs(os.path.join("tmp"))
    with open(os.path.join("tmp", tmpname+".fasta"), 'w') as o:
        o.write(">"+id1+"\n")
        o.write(seq1+"\n")
        o.write(">"+id2+"\n")
        o.write(seq2+"\n")
    muscle_line = MuscleCommandline(input=os.path.join("tmp", tmpname+".fasta"), out=os.path.join("tmp", tmpname+".fasta.out"))
    stdout, stderr = muscle_line()
    handle = open(os.path.join("tmp", tmpname+".fasta.out"), "rU")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    # print str(records[0].id)
    # print str(records[0].seq)
    # print str(records[1].id)
    # print str(records[1].seq)
    return(records[0].seq,records[1].seq)


def mapCDtoPDB(line,base_dir):
    columns = line.split(",")
                
    cd = # nombre base del alineamiento
    cdpdbseq = # id de la secuencia en el alineamiento
    cdfastapath = os.path.join(base_dir, "stockholms", cd+".sto")
    cdsequence, colmap = searchCDsequence(cdfastapath, cdpdbseq)

    pdbid = "3QK3"
    chain = "B"
    # resolution = float(columns[4])
    desde = # desde donde del PDB
    hasta = # hasta donde del PDB

    if os.path.exists( pdbid+"_"+chain+".pdb"):     # si existe el PDB
        with open(pdbid+"_"+chain+".pdb", 'r') as f:    # abro el PDB
            pdbtext = f.read()
            
            ### aca usa BioPython para obtener 2 listas, una con los numeros de residuos [resid_list] 
            ### y otra con los residuos (en letras) [resname_list]

            ## del estilo
            ## resid_list = [23,24,25,26,26]
            ## resname_list = ["A","N","P","V","K"]

            status = True
            resid_list, resname_list = pdblib.getPDBSequence(pdbtext, chain)
            with open("mi_mapa_seq2pdb.cdmap"), 'w') as outmap:
                print "processing ..."
                pdbseqmap = dict(zip(range(0,len(resid_list)),resid_list))
                pdbseq = "".join(resname_list)

                aln1, aln2 = MusclePairAlign("pdbseq", pdbseq, "cdseq", cdsequence.upper())

                # aln1 is PDB
                # aln2 is cdseq
                print aln1
                print aln2

                # alignments = pairwise2.align.globalxx(pdbseq,cdsequence.upper())
                # top_aln = alignments[0]
                # aln1, aln2, score, begin, end = top_aln        #aln1 contains pdbseq, aln2 contains CD seq
                
                outmap.write("# pdbseq aln1 "+str(aln1)+"\n")
                outmap.write("# fasseq aln2 "+str(aln2)+"\n")
                outmap.write("# alignment_col pdbres seqres pdbseq seqaa "+"\n")
                cdoffset = 0
                cdabspos = 0
                pdboffset = 0
                coverage = 0
                for i in range(0, len(aln1)):
                    if aln1[i] != "-" and re.search(r'\w', aln1[i]): 
                        if aln2[i] != "-" and re.search(r'\w', aln2[i]):
                            #### print colmap[i-cdoffset],":",resid_list[i-pdboffset],aln1[i], desde+i-cdoffset, aln2[i], cdsequence[i-cdoffset]
                            # print str(colmap[i-cdoffset])+","+str(resid_list[i-pdboffset])+","+str(desde+i-cdoffset)+","+aln2[i]+","+ cdsequence[i-cdoffset]+"\n"
                            outmap.write(str(colmap[i-cdoffset])+","+str(resid_list[i-pdboffset])+","+str(desde+i-cdoffset)+","+aln1[i]+","+ cdsequence[i-cdoffset]+"\n")
                            coverage += 1
                        else:
                            cdoffset+=1
                    else:
                        pdboffset+=1
                rel = float(coverage)/float(len(cdsequence))
                outmap.write("# coverage "+str(rel)+"\n")
                return None
            else:
                print "Error with PDB "+pdbid

if __name__ == "__main__":
    mapCDtoPDB()