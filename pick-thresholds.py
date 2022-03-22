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
from Rex import Rex
rex=Rex()

def getAcc(filename,accMap):
    if(not rex.find("eval.*-([\d\.]+).txt",filename)):
        raise Exception("Cannot parse filename: "+filename)
    threshold=rex[1]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=4): continue
            if(fields[0]=="sites"): continue
            (sites,reads,theta,acc)=fields
            key=sites+"\t"+reads
            acc=float(acc)
            if(accMap.get(key) is None):
                accMap[key]={}
            if(accMap[key].get(threshold) is None):
                accMap[key][threshold]=[]
            accMap[key][threshold].append(acc)
    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)<2):
    exit(ProgramName.get()+" <accuracy-files>\n")
filenames=sys.argv[1:]

print("sites\treads\tthreshold\tacc")
accMap={}
for filename in filenames:
    getAcc(filename,accMap)
keys=accMap.keys()
for key in keys:
    bestThreshold=None
    bestAcc=None
    thresholds=accMap[key].keys()
    for threshold in thresholds:
        values=accMap[key][threshold]
        N=len(values)
        acc=sum(values)/N
        #print(key,threshold,ave,sep="\t")
        if(bestAcc is None or acc>bestAcc):
            bestAcc=acc
            bestThreshold=threshold
    print(key,bestThreshold,round(bestAcc,3),sep="\t")
    




