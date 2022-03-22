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

def load(filename):
    table=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(fields[0]=="sites"): continue
            (sites,reads,theta,acc)=fields
            table.append([sites,reads,theta,acc])
    return table

def mismatch(row1,row2):
    for i in range(3):
        if(row1[i]!=row2[i]): return True
    return False
    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <trio-eval> <indiv-eval>\n")
(trioFile,indivFile)=sys.argv[1:]

print("sites\treads\ttheta\ttrio\tindiv")
trioTable=load(trioFile)
indivTable=load(indivFile)
N=len(trioTable)
if(len(indivTable)!=N): raise Exception("Unequal number of rows")
for i in range(N):
    trioRow=trioTable[i]
    indivRow=indivTable[i]
    if(mismatch(trioRow,indivRow)): raise Exception("Different row order")
    trioRow.append(indivRow[3])
    line="\t".join(trioRow)
    print(line)
    


