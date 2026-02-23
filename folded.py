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
import math
import ProgramName

def C(p,n,maxK):
    s=0
    q=1-p
    for k in range(maxK+1):
        s+=math.comb(n,k)*(p**k)*(q**(n-k))
    return s

def condMoment(m,p,n,maxK):
    s=0
    q=1-p
    for k in range(maxK+1):
        s+=(k**m)*math.comb(n,k)*(p**k)*(q**(n-k))
    return 1/C(p,n,maxK)*s

def MSE(n,maxK,p):
    return 1/(n**2)*condMoment(2,p,n,maxK)-2*p/n*condMoment(1,p,n,maxK)+p*p

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <trueP> <N> <maxK>\n")
(trueP,n,maxK)=sys.argv[1:]
trueP=float(trueP)
n=int(n)
maxK=int(maxK)

#mse=MSE(n,maxK,trueP)
#print(mse)

for p in [x/20 for x in range(1,10)]:
    print(p,round(MSE(n,maxK,p),5),sep="\t")





