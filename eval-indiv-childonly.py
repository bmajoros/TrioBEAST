#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2021 William H. Majoros <bmajoros@alumni.duke.edu>
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
from Rex import Rex
rex=Rex()

def getAffected(node,label):
    child=node.findChild(label)
    return int(child[0])>0 or int(child[1])>0

def hasHets(root,indiv):
    sites=root.findChildren("site")
    for site in sites:
        node=site.findChild("genotypes").findChild(indiv)
        if(node is None): raise Exception("Can't find node in hasHets")
        if(node[0]!=node[1]): return True
    return False

def loadPreds(predDir,indiv,parmString):
    preds=[]
    filename=predDir+"/indiv-"+indiv+"-"+parmString+".txt"
    with open(filename,"rt") as IN:
        header=IN.readline()
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=5): raise Exception("Can't parse: "+line)
            (geneID,median,Palt,left,right)=fields
            preds.append(float(Palt))
    return preds

def predict(Palt,threshold):
    return Palt>=threshold

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <truth.essex> <predictions-dir> <posterior-threshold>\n")
(truthFile,predDir,threshold)=sys.argv[1:]
threshold=float(threshold)

rex.findOrDie("truth-(.*).essex",truthFile)
parmString=rex[1]
childPreds=loadPreds(predDir,"child",parmString)

essexReader=EssexParser(truthFile)
numTrios=0; right=0; wrong=0
triosWithEvidence=0; rightWithEvidence=0; wrongWithEvidence=0
geneIndex=0
while(True):
    root=essexReader.nextElem()
    if(root is None): break
    numTrios+=1
    sxAffected=root.findChild("affected")
    trueChild=getAffected(sxAffected,"child")
    childPred=predict(childPreds[geneIndex],threshold)
    geneIndex+=1
    if(trueChild==childPred): right+=1
    else: wrong+=1
N=right+wrong
acc=float(right)/float(N)
print(round(acc*100,3),"% correct",sep="")



