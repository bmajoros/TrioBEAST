#!/usr/bin/env python
#=========================================================================
# Copyright (C)William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import os
import ProgramName
from Pipe import Pipe
from SlurmWriter import SlurmWriter
from Rex import Rex
rex=Rex()

BASE="/hpc/group/majoroslab/trios"
MCMC_SAMPLES=1000
GENE_RANGE="0-999"
LAMBDA=1.2
DENOVO_RATE=0.0001
RECOMB_RATE=0.001
PRIOR_AFFECTED=0.5
JOB_NAME="TRIO"
MAX_PARALLEL=300
additional_SBATCH_lines=""

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <inputs-dir> <outputs-dir> <slurm-dir>\n")
(inDir,outDir,slurmDir)=sys.argv[1:]

slurm=SlurmWriter()
files=os.listdir(inDir)
for f in files:
    if(not rex.find("^phased-(.+).essex",f)): continue
    parms=rex[1]
    inFile=inDir+"/"+f
    outFile=outDir+"/trio-"+parms+".txt"
    cmd="cd "+BASE+";\n"
    cmd+="git/run1.py git/TripleHets "+inFile+" "+str(MCMC_SAMPLES)+\
        " "+GENE_RANGE+" "+str(LAMBDA)+" "+str(DENOVO_RATE)+" "+\
        str(RECOMB_RATE)+" "+str(PRIOR_AFFECTED)+" >  "+outFile
    slurm.addCommand(cmd)
slurm.mem(1500)
slurm.setQueue("scavenger")
slurm.writeArrayScript(slurmDir,JOB_NAME,MAX_PARALLEL,
                           additional_SBATCH_lines)
    


