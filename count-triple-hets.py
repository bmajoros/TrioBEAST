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

def isHet(sxGenotypes,label):
    child=sxGenotypes.findChild(label)
    return child[0]!=child[1]

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <in.essex>\n")
(filename,)=sys.argv[1:]

siteCounts=[]
tripleHetCounts=[]
parser=EssexParser(filename)
while(True):
    sxGene=parser.nextElem()
    if(sxGene is None): break
    genotypes=sxGene.findDescendents("genotypes")
    sitesPerGene=0
    tripleHets=0
    for site in genotypes:
        sitesPerGene+=1
        if(isHet(site,"mother") and isHet(site,"father")
           and isHet(site,"child")):
            tripleHets+=1
    siteCounts.append(sitesPerGene)
    tripleHetCounts.append(tripleHets)
OUT=open("triple-hets.txt","wt")
for x in tripleHetCounts:
    print(x,file=OUT)
OUT.close()



