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
import math

def gam(x):
    return math.gamma(x)

def oneSite(m,n,x,c):
    numer=(x+1)*(x+2)-2*c*(x+1)*(n+3)
    denom=(n+2)*(n+3)
    return float(numer)/float(denom)+c*c

def getH(n,m,x,y):
    numer=gam(n+m+2)
    denom=gam(x+y+1)*gam(n+m-x-y+1)+gam(x+m-y+1)*gam(n-x+y+1)
    return numer/denom

def twoSites(m,n,x,y,c):
    H=getH(n,m,x,y)
    term1=gam(x+y+3)/gam(n+m+4)*gam(n+m-x-y+1)
    term2=gam(x+m-y+3)/gam(n+m+4)*gam(n-x+y+1)
    term3=-2*c*gam(x+y+2)/gam(n+m+x+y+3)*gam(n+m+1)
    term4=-c*c*gam(x+m-y+2)/gam(n+m+3)*gam(n-x+y+1)
    total=H*(term1+term2+term3+term4)+c*c
    return total


#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <m> <n> <x> <y> <c>\n")
(m,n,x,y,c)=sys.argv[1:]
m=int(m); n=int(n); x=int(x); y=int(y); c=float(c)

diff=twoSites(m,n,x,y,c)-oneSite(m,n,x,c)
print(diff)




