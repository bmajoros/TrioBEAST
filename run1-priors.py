#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import os
import math
import ProgramName
from Rex import Rex
rex=Rex()
import TempFilename
from Pipe import Pipe
from SummaryStats import SummaryStats
import getopt
from EssexParser import EssexParser
from Stan import Stan
from StanParser import StanParser
import numpy as np

DEBUG=False
WARMUP=300
STDERR=TempFilename.generate(".stderr")
INPUT_FILE=TempFilename.generate(".staninputs")
INIT_FILE=TempFilename.generate(".staninit")
OUTPUT_TEMP=TempFilename.generate(".stanoutputs")

MODES = [ "00 00 00 = all unaffected",
"00 00 10 = child has a de novo in the causal variant",
"00 00 01 = child has a de novo in the causal variant",
"01 00 00 = mother affected, child doesn't inherit",
"01 00 10 = mother affected and recombines, child inherits",
"00 01 00 = father affected, child doesn't inherit",
"00 01 00 = father affected and recombines, child inherits",
"10 00 10 = mother affected, child inherits",
"10 00 00 = mother affected and recombines, child doesn't inherit",
"00 10 01 = father affected, child inherits",
"00 10 00 = father affected and recombines, child doesn't inherit"
]
NUM_MODES=11

# Modes with child ASE: 1, 2, 4, 7, 9

#=========================================================================
class Mode:
    def __init__(self,index,posterior):
        self.index=index
        self.posterior=posterior
    def getDescription(self):
        return MODES[self.index]
    def isChildAffected(self):
        desc=self.getDescription()
        if(not rex.find("^\d\d \d\d (\d\d)",desc)):
            raise Exception("Cannot parse: "+desc)
        childGT=rex[1]
        return childGT[0]=="1" or childGT[1]=="1"
#=========================================================================
class Site:
    def __init__(self,ID,phased):
        self.ID=ID
        self.phased=phased
        self.counts=np.zeros((3,2),int) # [indiv][haplotype]
        self.het=[0]*3 # [indiv]
#=========================================================================
class Gene:
    def __init__(self,ID):
        self.ID=ID
        self.sites=[]
    def addSite(self,site):
        self.sites.append(site)
#=========================================================================
def parseGene(sxGene):
    ID=sxGene[0];
    gene=Gene(ID);
    numSites=sxGene.numElements()-1
    for i in range(numSites):
        sxSite=sxGene[i+1];
        ID=sxSite[0]
        sxGenotypes=sxSite.findChild("genotypes")
        sxCounts=sxSite.findChild("counts")
        isPhased=sxSite.getAttribute("phased")
        site=Site(ID,isPhased)
        site.het[0]=isHet(sxGenotypes,"mother")
        site.het[1]=isHet(sxGenotypes,"father")
        site.het[2]=isHet(sxGenotypes,"child")
        site.counts[0]=getCounts(sxCounts,"mother")
        site.counts[1]=getCounts(sxCounts,"father")
        site.counts[2]=getCounts(sxCounts,"child")
        gene.addSite(site)
    return gene

def getCounts(sxCounts,label):
    child=sxCounts.findChild(label)
    return np.array([child[0],child[1]])

def isHet(sxGenotypes,label):
    child=sxGenotypes.findChild(label)
    return child[0]!=child[1]

def printFields(fields,hFile):
    numFields=len(fields)
    for i in range(7,numFields):
        print(i-6,"=",fields[i],sep="",end="",file=hFile)
        if(i<numFields-1): print("\t",end="",file=hFile)
    print(file=hFile)

def getFieldIndex(label,fields):
    numFields=len(fields)
    index=None
    for i in range(7,numFields):
        if(fields[i]==label): index=i
    return index

def writeInitializationFile(filename):
    OUT=open(filename,"wt")
    #print("theta <- 1",file=OUT)
    print("probDenovo <- ",0.5,file=OUT)
    print("probRecomb <- ",0.5,file=OUT)
    OUT.close()

def writeInputsFile(stan,gene,probDenovo,probRecomb,probAffected,filename):
    OUT=open(filename,"wt")

    # N_SITES
    N_SITES=len(gene.sites)
    print("N_SITES <- ",N_SITES,file=OUT)

    # int<lower=0,upper=1> het[N_SITES,3]
    het=np.zeros((N_SITES,3),int)
    for i in range(N_SITES): het[i]=gene.sites[i].het
    stan.writeTwoDimArray("het",het,N_SITES,3,OUT)

    # int<lower=0> count[N_SITES,3,2]
    count=np.zeros((N_SITES,3,2),int)
    for i in range(N_SITES): count[i]=gene.sites[i].counts
    stan.writeThreeDimArray("count",count,N_SITES,3,2,OUT)

    # int<lower=0,upper=1> isPhased[N_SITES]
    phased=[]
    for site in gene.sites: phased.append(site.phased)
    stan.writeOneDimArray("isPhased",phased,N_SITES,OUT)

    # Probabilities
    print("probAffected <- ",probAffected,file=OUT)
    OUT.close()

def runGene(stan,gene,numSamples,probDenovo,probRecomb,probAffected,
            Lambda):
    # Write inputs file for STAN
    writeInputsFile(stan,gene,probDenovo,probRecomb,probAffected,INPUT_FILE)
    writeInitializationFile(INIT_FILE)
    
    # Run STAN model
    if(DEBUG):
        cmd=stan.getCmd(WARMUP,numSamples,INPUT_FILE,OUTPUT_TEMP,STDERR,
                        INIT_FILE)
        print(cmd)
        exit()
    else:
        stan.run(WARMUP,numSamples,INPUT_FILE,OUTPUT_TEMP,STDERR,INIT_FILE)

    # Parse MCMC output
    parser=StanParser(OUTPUT_TEMP)

    # Get estimates
    (median,CI_left,CI_right)=parser.getMedianAndCI(0.95,"theta")
    left=parser.getLeftTail("theta",1.0/Lambda)
    right=parser.getRightTail("theta",Lambda)
    #P_alt=left+right
    P_alt=max(left,right)

    # Estimate mode of inheritance
    modes=getInheritanceMode(parser)

    # Return estimates
    return (median,P_alt,CI_left,CI_right,modes)

def getInheritancePosterior(i,parser,denom):
    array=[]
    numer=parser.getVariable("numerator."+str(i+1))
    numSamples=len(numer)
    for j in range(numSamples):
        posterior=math.exp(numer[j]-denom[j]);
        array.append(posterior)
    return sum(array)/len(array)

def getInheritanceMode(parser):
    modes=[]
    denom=parser.getVariable("denominator");
    for i in range(NUM_MODES):
        posterior=getInheritancePosterior(i,parser,denom)
        modes.append(Mode(i,posterior))
    modes.sort(key=lambda x: 1-x.posterior)
    return modes

def probChildIsAffected(modes):
    total=0
    for mode in modes:
        if(mode.isChildAffected()): total+=mode.posterior
    return total

#=========================================================================
# main()
#=========================================================================
(options,args)=getopt.getopt(sys.argv[1:],"s:")
if(len(args)!=8):
    exit(ProgramName.get()+" [-s file] <model> <input.essex> <#MCMC-samples> <firstGene-lastGene> <lambda=1.2> <P(de novo)> <P(recomb)> <P(affected)>\n   -s = save raw STAN file\n   gene range is zero-based and inclusive\n")
(model,inputFile,numSamples,geneRange,Lambda,probDenovo,probRecomb,
 probAffected)=args
stanFile=None
for pair in options:
    (key,value)=pair
    if(key=="-s"): stanFile=value
if(not rex.find("(\d+)-(\d+)",geneRange)):
    exit(geneRange+": specify range of gene: first-last")
firstIndex=int(rex[1])
lastIndex=int(rex[2])
Lambda=float(Lambda)

#print(OUTPUT_TEMP) ###

# Process each gene
geneIndex=0
parser=EssexParser(inputFile)
print("Gene\tP(ASE)\tFoldChg\t95%CredIntv\tP(child_affect)")
while(True):
    elem=parser.nextElem()
    if(elem is None): break
    if(elem.getTag()!="gene"): raise Exception("Expecting 'gene' tag in eesex")
    if(geneIndex<firstIndex):
        geneIndex+=1
        continue
    elif(geneIndex>lastIndex): break
    gene=parseGene(elem)
    if(gene is None): continue
    stan=Stan(model)
    (median,P_alt,CI_left,CI_right,modes)=\
      runGene(stan,gene,numSamples,probDenovo,
              probRecomb,probAffected,Lambda)
    probChildAffected=probChildIsAffected(modes)
    P_alt=round(P_alt,3)
    median=round(median,3)
    CI_left=round(CI_left,3); CI_right=round(CI_right,3)
    print(gene.ID,P_alt,median,str(CI_left)+"-"+str(CI_right),
          round(probChildAffected,4),sep="\t")
    for mode in modes:
        if(mode.posterior>=0.01):
            print("\t",round(mode.posterior*100),"% : ",MODES[mode.index],
                  sep="")
    geneIndex+=1
os.remove(STDERR)
os.remove(INPUT_FILE)
if(stanFile is None):
    os.remove(OUTPUT_TEMP)
else:
    os.system("cp "+OUTPUT_TEMP+" "+stanFile)


