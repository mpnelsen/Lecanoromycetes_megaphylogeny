import random
import subprocess
import os

#os.mkdir('~/path...')

#file path and name (original alignment)
f="concat.phy"
cons="family_constraints.tre"

#number of bootstrap replicates
n=100

#number of processors to use
proc=10

#original data set - make binary, get starting tree, find ml tree
#subprocess.call(args="~/ExaML/parser/parser -s {0} -q 19may2014_concat_partfile.txt -m DNA -n forexaml".format(f), shell=True);
#y=random.randint(1,1000000)
#subprocess.call(args="~/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -T {0} -y -m GTRCAT -q 19may2014_concat_partfile.txt -s {1} -n MLparstart.tre -p {2}".format(proc, f, y), shell=True);
#subprocess.call(args="mpirun.openmpi -np {0} ~/ExaML/examl/examl -s forexaml.binary -D -m PSR -t RAxML_parsimonyTree.MLparstart.tre -n ML_TREE".format(proc), shell=True);

#create bootstrap data sets
subprocess.call(args="raxmlHPC-PTHREADS -T {0} -# {1} -b 12345 -f j -q 19may2014_concat_partfile.txt -m GTRCAT -s {2} -n boot_data".format(proc,n,f), shell=True);

#create parsimony starting tree for each data set	
for x in range (0, n):
	y=random.randint(1,1000000)
	#subprocess.call(args="raxmlHPC-PTHREADS-SSE3 -T {0} -y -m GTRCAT -q 19may2014_concat_partfile.txt -s {1}.BS{2} -n {2}_parstart.tre -p {3} -g {4}".format(proc,f,x,y,cons), shell=True);
	subprocess.call(args="~/ExaML/parser/parse-examl -s {0}.BS{1} -q 19may2014_concat_partfile.txt -m DNA -n forexaml{1}".format(f,x), shell=True);
	subprocess.call(args="mpirun -np {0} ~/ExaML/examl/examl -s forexaml{1}.binary -D -m PSR -g {2} -p {3} -n {1}_OUT ".format(proc,x,cons,y), shell=True);

#concatenate bootstrap trees into one file
subprocess.call(args="cat ExaML_result.* > bootstrap_trees_combined.tre", shell=True);

#now...for constraint...
#	subprocess.call(args="~/ExaML/parser/parser -s {0} -q 19may2014_concat_partfile.txt -m DNA -n forexaml_family_cons1".format(f), shell=True);
#y=random.randint(1,1000000)
#subprocess.call(args="mpirun.openmpi -np {0} ~/ExaML/examl/examl -s forexaml_family_cons1.binary -D -m PSR -g {1} -p {2} -n subclass_family_constraint_ML_TREE".format(proc,cons,y), shell=True);
