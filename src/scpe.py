from subprocess import call
def run(pdb_file,beta,runs,nsus,chain,result_file,contact_map_path):
    for b in beta:
        for sus in nsus:
            for r in runs: 
                call(["scpe/scpe_program/scpe.exe", "-i" ,pdb_file,"-q","scpe/scpe_program/QijCodonesEmp" ,"-c", "A","-b",b,"-A", sus, "-a", sus, "-R", r,"-1","1","-2", result_file,"-3","25.0","-4","1","-5",contact_map_path])




