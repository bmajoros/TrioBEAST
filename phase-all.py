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
import os
import ProgramName
from Pipe import Pipe
from Rex import Rex
rex=Rex()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <subdir>\n")
(subdir,)=sys.argv[1:]

files=os.listdir(subdir)
for f in files:
    if(not rex.find("^data-(.+)",f)): continue
    rest=rex[1]
    inFile=subdir+"/"+f
    outFile=subdir+"/phased-"+rest
    cmd="git/phase-trio "+inFile+" "+outFile
    #print(cmd+"\n")
    Pipe.run(cmd)
    


