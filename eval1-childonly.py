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

class BufferedReader:
    def __init__(self,filename):
        self.fh=open(filename,"rt")
        self.fh.readline() # discard the header line
    def nextPred(self):
        while(True):
            line=self.fh.readline()
            while(not rex.find("^GENE",line)): 
                line=self.fh.readline()
                if(line is None): return None
            fields=line.rstrip().split()
            if(len(fields)!=5): raise Exception("Can't parse line: "+line)
            (ID,Palt,theta,CI,Pchild)=fields
            return (ID,Palt,theta,CI,Pchild)
          
def getAffected(node,label):
    child=node.findChild(label)
    return child[0]=="1" or child[1]=="1"

def compact(mother,father,child):
    return mother[0]+mother[1]+" "+father[0]+father[1]+" "+child[0]+child[1]

def isAffected(pair):
    return int(pair[0])>0 or int(pair[1])>0

def hasHets(root,indiv):
    sites=root.findChildren("site")
    for site in sites:
        node=site.findChild("genotypes").findChild(indiv)
        if(node is None): raise Exception("Can't find node in hasHets")
        if(node[0]!=node[1]): return True
    return False

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <truth.essex> <model-output> <threshold>\n")
(truthFile,outputFile,threshold)=sys.argv[1:]
threshold=float(threshold)

essexReader=EssexParser(truthFile)
predReader=BufferedReader(outputFile)
numTrios=0; right=0; wrong=0
while(True):
    root=essexReader.nextElem()
    if(root is None): break
    prediction=predReader.nextPred()
    if(prediction is None): break
    (ID,Palt,theta,CI,Pchild)=prediction
    if(ID!=root[0]): raise Exception(ID+" != "+root[0]);
    numTrios+=1
    sxAffected=root.findChild("affected")
    trueChild=int(getAffected(sxAffected,"child"))
    childAffected=float(Pchild)>=threshold
    if(trueChild==childAffected): right+=1
    else: wrong+=1
N=right+wrong
acc=float(right)/float(N)
print(round(acc*100,3),"% correct",sep="")


