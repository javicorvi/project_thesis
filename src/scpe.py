from subprocess import call
def run(pdb_file,beta,runs,nsus,chain,result_file,contact_map_path):
    for b in beta:
        for sus in nsus:
            for r in runs:
                try: 
                    call(["scpe/scpe_program/scpe.exe", "-i" ,pdb_file,"-q","scpe/scpe_program/QijCodonesEmp" ,"-c", chain,"-b",b,"-A", sus, "-a", sus, "-R", r,"-1","1","-2", result_file,"-3","25.0","-4","1","-5",contact_map_path])
                except Exception:
                    print "The SCPE execution get an exception with de pdb file " + pdb_file
                    raise Exception ("The SCPE execution get an exception with de pdb file " + pdb_file)   
                
                

def run_singular(pdb_file,beta,run,nsus,chain,output_msa_path,contact_map_path):
    try: 
        call(["scpe/scpe_program/scpe.exe", "-i" ,pdb_file,"-q","scpe/scpe_program/QijCodonesEmp" ,"-c", chain,"-b",beta,"-A", nsus, "-a", nsus, "-R", run,"-1","1","-2", output_msa_path,"-3","25.0","-4","1","-5",contact_map_path])
    except Exception:
        print "The SCPE execution get an exception with de pdb file " + pdb_file
        raise Exception ("The SCPE execution get an exception with de pdb file " + pdb_file)                   



