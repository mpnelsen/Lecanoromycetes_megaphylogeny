import random
import subprocess
import os

#os.mkdir('~/path...')

#file path and name (original alignment)
f="concat.phy"
cons="family_constraints.tre"

#number of bootstrap replicates
#n=100

#number of processors to use
proc=10

#original data set - make binary, get starting tree, find ml tree
#subprocess.call(args="~/ExaML/parser/parser -s {0} -q 19may2014_concat_partfile.txt -m DNA -n forexaml".format(f), shell=True);
#y=random.randint(1,1000000)
#subprocess.call(args="~/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -T {0} -y -m GTRCAT -q 19may2014_concat_partfile.txt -s {1} -n MLparstart.tre -p {2}".format(proc, f, y), shell=True);
#subprocess.call(args="mpirun.openmpi -np {0} ~/ExaML/examl/examl -s forexaml.binary -D -m PSR -t RAxML_parsimonyTree.MLparstart.tre -n ML_TREE".format(proc), shell=True);

#now...for constraint...
subprocess.call(args="~/ExaML/parser/parser -s {0} -q 19may2014_concat_partfile.txt -m DNA -n forexaml_family_cons1".format(f), shell=True);
y=random.randint(1,1000000)
subprocess.call(args="mpirun.openmpi -np {0} ~/ExaML/examl/examl -s forexaml_family_cons1.binary -D -m PSR -g {1} -p {2} -n subclass_family_constraint_ML_TREE".format(proc,cons,y), shell=True);
