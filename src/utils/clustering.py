import os
from subprocess import call
exten=".fasta"
import time
def run(clust, input_folder, output_folder, pattern_array=[".fasta"]):
    start_time = time.time()
    print "clustering_begins"
    for filename in os.listdir(input_folder):
        if filename.endswith(exten) & any(r in filename for r in pattern_array):
            print(filename)
            filenameclust = filename + "_"+clust + ".cluster"
            print(filenameclust)
            call(["cd-hit", "-i" , input_folder+"/"+filename ,"-o", output_folder+"/"+filenameclust,"-c",clust,"-n", "4", "-M", "3000"])
    print "clustering_ends"
    print("--- %s seconds ---" % (time.time() - start_time))        