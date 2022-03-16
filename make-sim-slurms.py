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
import ProgramName
from Pipe import Pipe

NUM_GENES=1000
SITES_PER_GENE=[2, 5, 10, 15]
READS_PER_SITE=[5, 10, 15, 20, 25, 30]
THETAS=[0.5, 0.4, 0.3, 0.2, 0.1]
RECOMB_RATE=0.001
BASE="/hpc/group/majoroslab/trios"


def sim(sites,reads,theta,vcfFile,motherID,fatherID,outDir):
    truthFile=outDir+"/truth-sites"+str(sites)+"-reads"+str(reads)+\
        "-theta"+str(theta)+".essex"
    dataFile=outDir+"/data-sites"+str(sites)+"-reads"+str(reads)+\
        "-theta"+str(theta)+".essex"
    cmd=BASE+"/git/sim1 "+vcfFile+" "+motherID+" "+fatherID+" "+\
        str(NUM_GENES)+" "+str(sites)+" "+str(reads)+" "+str(RECOMB_RATE)+\
        " "+str(theta)+" "+truthFile+" "+dataFile
    #print(cmd)
    Pipe.run(cmd)
    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <in.vcf.gz> <mother-ID> <father_ID> <out-subdir>\n")
(vcfFile,motherID,fatherID,outDir)=sys.argv[1:]

if(outDir[0]=="~" or outDir[0]=="/"):
    raise Exception("use relative path in subdirectory")

for sites in SITES_PER_GENE:
    for reads in READS_PER_SITE:
        for theta in THETAS:
            sim(sites,reads,theta,vcfFile,motherID,fatherID,outDir)


