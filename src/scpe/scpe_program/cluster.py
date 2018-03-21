import os
from subprocess import call
clust="0.62"
exten=".fasta"
#runsToCheck = ["runs1","runs10","runs100","runs500"]
runsToCheck = ["beta0.05-nsus5.00-runs15000","beta0.15-nsus5.00-runs5000"]
for filename in os.listdir("."):
    if filename.endswith(exten) & any(r in filename for r in runsToCheck):
        print(filename)
        filenameclust = filename + "_"+clust + ".cluster"
        print(filenameclust)
        call(["cd-hit", "-i" , filename ,"-o", "clusters3/"+filenameclust,"-c","0.62","-n", "4", "-M", "3000"])
    