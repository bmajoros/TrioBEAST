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
import os
import sys
import ProgramName
from SlurmWriter import SlurmWriter
from Rex import Rex
rex=Rex()

BASE="/hpc/group/majoroslab/trios"
MCMC_SAMPLES=1000
JOB_NAME="INDIV"
MAX_PARALLEL=300
additional_SBATCH_lines=""

def makeCommand(subdir,parmString,indivID,model,Lambda,outDir):
    data=subdir+"/phased-"+parmString+".essex"
    outfile=outDir+"/indiv-"+indivID+"-"+parmString+".txt"
    cmd="cd "+BASE+";\n"
    cmd+="git/run-phased-indiv.py "+model+" "+data+\
        " "+indivID+" "+str(MCMC_SAMPLES)+" "+str(Lambda)+" > "+outfile
    return cmd
    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <slurms-dir> <sims-dir> <model> <lambda> <out-dir>\n")
(slurmsDir,simsDir,model,Lambda,outDir)=sys.argv[1:]

slurm=SlurmWriter()
files=os.listdir(simsDir)
for filename in files:
    if(not rex.find("^phased-(.+).essex",filename)): continue
    parmString=rex[1]
    cmd1=makeCommand(simsDir,parmString,"mother",model,Lambda,outDir)
    cmd2=makeCommand(simsDir,parmString,"father",model,Lambda,outDir)
    cmd3=makeCommand(simsDir,parmString,"child",model,Lambda,outDir)
    slurm.addCommand(cmd1)
    slurm.addCommand(cmd2)
    slurm.addCommand(cmd3)
slurm.mem(1500)
slurm.setQueue("scavenger")
slurm.writeArrayScript(slurmsDir,JOB_NAME,MAX_PARALLEL,
                       additional_SBATCH_lines)

