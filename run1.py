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

WARMUP=300
STDERR=TempFilename.generate(".stderr")
INPUT_FILE=TempFilename.generate(".staninputs")
INIT_FILE=TempFilename.generate(".staninit")
OUTPUT_TEMP=TempFilename.generate(".stanoutputs")

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
    return 

def isHet(sxGenotypes,label):
    child=sxGenotypes.findChild(label)
    return child[0]!=child[1]

OLD:
def parseCounts(sxCounts):
    if(sxCounts is None): raise Exception("no counts")
    n=sxCounts.numElements()
    counts=[]
    for i in range(n):
        c=int(sxCounts.getIthElem(i))
        counts.append(c)
    return counts

def parseBackground(sxBackground):
    if(sxBackground is None): raise Exception("no background")
    n=sxBackground.numElements()
    bg=[]
    for i in range(n):
        p=float(sxBackground.getIthElem(i))
        bg.append(p)
    return bg

def parseRep(sxRep):
    repNum=int(sxRep.getIthElem(0))
    background=parseBackground(sxRep.findChild("background"))
    counts=parseCounts(sxRep.findChild("counts"))
    rep=Replicate(repNum,background,counts)
    return rep

def parseGuide(sxGuide):
    guideID=sxGuide.getIthElem(0)
    sxReps=sxGuide.findChildren("rep")
    if(len(sxReps)==0): raise Exception("no reps found for guide "+guideID)
    reps=[]
    for sxRep in sxReps:
        rep=parseRep(sxRep)
        #if(not DISCARD_ZEROS or not rep.hasZeros()):
        if(rep.numZeros()==0 or
           not DISCARD_ZEROS and rep.numBins()-rep.numZeros()>=2):
            reps.append(rep)
    if(len(reps)==0): return None
    guide=Guide(guideID,reps)
    return guide

def parseGene(sxGene):
    geneID=sxGene.getIthElem(0)
    sxGuides=sxGene.findChildren("guide")
    if(len(sxGuides)==0): raise Exception("No guides found for gene "+geneID)
    guides=[]
    for sxGuide in sxGuides:
        guide=parseGuide(sxGuide)
        if(guide is None): 
            #print("XXX DROPPING GUIDE "+sxGuide[0])
            continue
        guides.append(guide)
    if(len(guides)==0): 
        #print("XXX DROPPING GENE "+geneID)
        return None
    gene=Gene(geneID,guides)
    return gene

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

def writeInitializationFile(filename,gene):
    OUT=open(filename,"wt")
    print("theta <- 1",file=OUT)
    OUT.close()

def writeInputsFile(gene,filename):
    numGuides=gene.numGuides()
    numReps=gene.maxReps()
    K=gene.maxBins()

    # NG, NREP, K
    OUT=open(filename,"wt")
    print("NG <-",numGuides,file=OUT) ###
    print("MAXREPS <- ",numReps,file=OUT)
    print("K <-",K,file=OUT)

    # NREP[NG]
    NREP=[]
    for guide in gene.guides:
        NREP.append(guide.numReps())
    writeOneDimArray("NREP",NREP,numGuides,OUT)

    # NUM_BINS[NG,NREP]
    NUM_BINS=initArray2D(numGuides,numReps,-1)
    i=0
    for guide in gene.guides:
        j=0
        for rep in guide.reps:
            NUM_BINS[i][j]=rep.numBins()
            j+=1
        i+=1
    writeTwoDimArray("NUM_BINS",NUM_BINS,numGuides,numReps,OUT)

    # vector[K] PI_NEG[NG,NREP]
    PI_NEG=initArray3D(numGuides,numReps,K,-1)
    i=0
    for guide in gene.guides:
        j=0
        for rep in guide.reps:
            for k in range(len(rep.background)):
                PI_NEG[i][j][k]=rep.background[k]
            j+=1
        i+=1
    writeThreeDimArray("DIR_NEG",PI_NEG,numGuides,numReps,K,OUT)

    # X[NG,NREP,K]
    X=initArray3D(numGuides,numReps,K,-1)
    i=0
    for guide in gene.guides:
        j=0
        for rep in guide.reps:
            for k in range(len(rep.background)):
                X[i][j][k]=rep.counts[k]
            j+=1
        i+=1
    writeThreeDimArray("X",X,numGuides,numReps,K,OUT)

    # Sigma
    print("sigma <- ",sigma,file=OUT)

    # dropout
    print("dropout <- ",DROPOUT,file=OUT)
    OUT.close()

def runGene(gene,model,numSamples,outfile,sigma,Lambda):
    # Write inputs file for STAN
    writeInputsFile(gene,INPUT_FILE)

    # Run STAN model
    init=" init="+INIT_FILE if SHOULD_INITIALIZE else ""
    cmd=model+" sample thin=1"+\
        " num_samples="+numSamples+\
        " num_warmup="+str(WARMUP)+\
        " data file="+INPUT_FILE+\
        init+\
        " output file="+OUTPUT_TEMP+" refresh=0 > "+STDERR
    #print(cmd); 
    #exit(cmd)
    os.system(cmd)

    # Parse MCMC output
    parser=StanParser(OUTPUT_TEMP)
    betas=betasFromW(parser,gene.numGuides()) ### 
    #betas=parser.getVariable("B")
    betas=[-5*x for x in betas] # Convert beta from 6-bin scale to 2-bin scale
    if(outfile!=""): os.system("cp "+OUTPUT_TEMP+" "+outfile)

    # Get posterior median
    betas.sort(key=lambda x: x)
    n=len(betas)
    median=None
    mid=int(n/2)
    if(n%2==0): median=(betas[mid]+betas[mid+1])/2.0
    else: median=betas[mid]
    (lower,upper)=getCI(betas,0.95)
    median=round(median,3)

    # Get P_reg
    greater=countRight(betas,math.log(Lambda))
    less=countLeft(betas,math.log(1.0/Lambda))
    Preg=max(float(greater)/float(n),float(less)/float(n))

    # Return estimates
    median=math.exp(median)
    lower=math.exp(lower)
    upper=math.exp(upper)
    return (Preg,median,lower,upper)

#=========================================================================
# main()
#=========================================================================
(options,args)=getopt.getopt(sys.argv[1:],"s:")
if(len(args)!=7):
    exit(ProgramName.get()+" [-s file] <model> <input.essex> <samples-dir-or-dot> <#MCMC-samples> <firstGene-lastGene> <beta-shrinkage-parm> <lambda=1.2>\n   -s = save raw STAN file\n   gene range is zero-based and inclusive\n   use . for samples-dir if no samples desired")
(model,inputFile,samplesDir,numSamples,geneRange,sigma,Lambda)=args
stanFile=None
for pair in options:
    (key,value)=pair
    if(key=="-s"): stanFile=value
if(not rex.find("(\d+)-(\d+)",geneRange)):
    exit(geneRange+": specify range of gene: first-last")
firstIndex=int(rex[1])
lastIndex=int(rex[2])
Lambda=float(Lambda)
#Pipe.run("rm -r "+samplesDir)
#Pipe.run("mkdir "+samplesDir)

# Process each gene
geneIndex=0
parser=EssexParser(inputFile)
print("Gene\tSignif\tP(reg)\tFoldChg\t95%CredIntv")
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
    if(ONLY_2BIN): gene.dropReps(6)
    elif(ONLY_6BIN): gene.dropReps(2)
    if(gene.numGuides()==0): 
        geneIndex+=1
        continue
    if(REVERSE_BINS): gene.reverseBins()
    if(IGNORE_MIDDLE_BINS): gene.dropMiddleBins()
    elif(ONLY_MIDDLE_BINS): gene.onlyMiddleBins()
    if(SHOULD_INITIALIZE):
        writeInitializationFile(INIT_FILE,gene)
    outfile="" if samplesDir=="." else samplesDir+"/"+gene.ID+".samples"
    (Preg,posterior,lowerCI,upperCI)=\
        runGene(gene,model,numSamples,outfile,sigma,Lambda)
    Preg=round(Preg,3)
    posterior=round(posterior,3)
    lowerCI=round(lowerCI,3); upperCI=round(upperCI,3)
    regulatory="no"
    if(Preg>=0.7): regulatory="maybe"
    if(Preg>=0.9): regulatory="yes"
    print(gene.ID,regulatory,Preg,posterior,str(lowerCI)+"-"+str(upperCI),
          sep="\t")
    geneIndex+=1
os.remove(STDERR)
os.remove(INPUT_FILE)
if(stanFile is None):
    os.remove(OUTPUT_TEMP)
else:
    os.system("cp "+OUTPUT_TEMP+" "+stanFile)


