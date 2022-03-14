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
from GffTranscriptReader import GffTranscriptReader
from Rex import Rex
rex=Rex()

def getChromLine(filename):
    pipe=Pipe("cat "+filename+" | gunzip ")
    while(True):
        line=pipe.readline()
        if(line is None): break
        if(rex.find("^#CHROM",line)):
            pipe.close()
            return line
    return None

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <*.gtf> <indexed.vcf.gz>\n")
(gffFile,vcfFile)=sys.argv[1:]

chromLine=getChromLine(vcfFile)
print(chromLine)
reader=GffTranscriptReader()
genes=reader.loadGenes(gffFile)
for gene in genes:
    Chr=gene.getSubstrate()
    geneID=gene.getID()
    if(rex.find("chr(.+)",Chr)): Chr=rex[1]
    exons=gene.getMergedExons()
    for exon in exons:
        cmd="tabix "+vcfFile+" "+Chr+":"+str(exon.begin)+"-"+str(exon.end)
        pipe=Pipe(cmd)
        while(True):
            line=pipe.readline()
            if(line is None or line==""): break
            fields=line.rstrip().split()
            fields[2]=geneID+":"+fields[2]
            line="\t".join(fields)
            print(line)
            
            





