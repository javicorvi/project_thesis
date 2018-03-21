from subprocess import call
beta = ["0.05"]
run = ["1"]
nsus = ["1.0"]
for b in beta:
	for sus in nsus:
		for r in run: 
			call(["./scpe.exe", "-i" ,"./Example/2trx.pdb" ,"-c", "A","-b",b,"-A", sus, "-a", sus, "-R",r])



