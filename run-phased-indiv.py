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
from EssexParser import EssexParser
from Stan import Stan
from StanParser import StanParser
import numpy as np
import TempFilename

DEBUG=False
WARMUP=300
STDERR=TempFilename.generate(".stderr")
INPUT_FILE=TempFilename.generate(".staninputs")
INIT_FILE=TempFilename.generate(".staninit")
OUTPUT_TEMP=TempFilename.generate(".stanoutputs")

def getPair(sxNode):
    return [sxNode[0],sxNode[1]]

def isHet(gt):
    return gt[0]!=gt[1]

def parseGene(sxGene,indivID):
    geneID=sxGene[0]
    sites=sxGene.findChildren("site")
    countsArray=[]
    phasedArray=[]
    for sxSite in sites:
        sxGenotype=sxSite.findChild("genotypes").findChild(indivID)
        sxCount=sxSite.findChild("counts").findChild(indivID)
        gt=getPair(sxGenotype)
        if(not isHet(gt)): continue
        counts=getPair(sxCount)
        countsArray.append(counts)
        phased=int(sxSite.getAttribute("phased"))
        phasedArray.append(phased)
    return (geneID,countsArray,phasedArray)

def writeInitializationFile(filename):
    OUT=open(filename,"wt")
    print("theta <- 1",file=OUT)
    OUT.close()

def writeInputsFile(stan,counts,phased,filename):
    OUT=open(filename,"wt")

    # N_SITES
    N_SITES=len(counts)
    print("N_SITES <- ",N_SITES,file=OUT)

    # int<lower=0> count[N_SITES,2];
    #stanCounts=np.zeros((N_SITES,2),int)
    #for i in range(N_SITES): stanCounts[i]=counts[i]
    stan.writeTwoDimArray("count",counts,N_SITES,2,OUT)

    # int<lower=0,upper=1> isPhased[N_SITES]
    stan.writeOneDimArray("isPhased",phased,N_SITES,OUT)
    
    OUT.close()

def run(stan,counts,phased,NUM_SAMPLES,LAMBDA):
    # Write inputs file for STAN
    writeInputsFile(stan,counts,phased,INPUT_FILE)
    writeInitializationFile(INIT_FILE)

    # Run STAN model
    if(DEBUG):
        cmd=stan.getCmd(WARMUP,NUM_SAMPLES,INPUT_FILE,OUTPUT_TEMP,STDERR,
                        INIT_FILE)
        print(cmd)
        exit()
    else:
        stan.run(WARMUP,NUM_SAMPLES,INPUT_FILE,OUTPUT_TEMP,STDERR,INIT_FILE)
        #cmd=stan.getCmd(WARMUP,NUM_SAMPLES,INPUT_FILE,OUTPUT_TEMP,STDERR,
        #                INIT_FILE)
        #print(cmd)

    # Parse MCMC output
    parser=StanParser(OUTPUT_TEMP)

    # Get estimates
    (median,CI_left,CI_right)=parser.getMedianAndCI(0.95,"theta")
    left=parser.getLeftTail("theta",1.0/LAMBDA)
    right=parser.getRightTail("theta",LAMBDA)
    #P_alt=left+right
    P_alt=max(left,right)

    # Return estimates
    return (median,P_alt,CI_left,CI_right)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <model> <phased-data.essex> <mother|father|child> <#MCMC-samples> <lambda>\n")
(modelFile,dataFile,indivID,NUM_SAMPLES,LAMBDA)=sys.argv[1:]
NUM_SAMPLES=int(NUM_SAMPLES)
LAMBDA=float(LAMBDA)

stan=Stan(modelFile)
parser=EssexParser(dataFile)
print("gene\tmedian\tP(alt)\tCI_left\tCI_right")
while(True):
    sxGene=parser.nextElem()
    if(sxGene is None): break
    (geneID,counts,phased)=parseGene(sxGene,indivID)
    if(len(counts)==0):
        (median,P_alt,CI_left,CI_right)=(1,0,0,1000)
        print(geneID,median,P_alt,CI_left,CI_right,sep="\t")
        continue
    (median,P_alt,CI_left,CI_right)=\
        run(stan,counts,phased,NUM_SAMPLES,LAMBDA)
    print(geneID,median,P_alt,CI_left,CI_right,sep="\t")
