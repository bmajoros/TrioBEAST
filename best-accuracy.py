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

def getAcc(filename):
    if(not rex.find("eval.*-([\d\.]+).txt",filename)):
        raise Exception("Cannot parse filename: "+filename)
    threshold=rex[1]
    Sum=0; N=0
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=4): continue
            if(fields[0]=="sites"): continue
            acc=float(fields[3])
            Sum+=acc
            N+=1
    ave=Sum/N
    return (ave,threshold)
    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)<2):
    exit(ProgramName.get()+" <accuracy-files>\n")
filenames=sys.argv[1:]

bestAcc=None
bestThreshold=None
for filename in filenames:
    (acc,threshold)=getAcc(filename)
    if(bestAcc is None or acc>bestAcc):
        bestAcc=acc
        bestThreshold=threshold
print("best threshold=",bestThreshold,", best acc=",bestAcc)        


