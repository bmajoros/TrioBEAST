#!/usr/bin/env python
#=========================================================================
# Copyright (C)2021 William H. Majoros <bmajoros@alumni.duke.edu>
#
# Permission is hereby granted, free of charge, to any person obtaining 
# a copy of this software and associated documentation files (the 
# "Software"), to deal in the Software without restriction, including 
# without limitation the rights to use, copy, modify, merge, publish, 
# distribute, sublicense, and/or sell copies of the Software, and to 
# permit persons to whom the Software is furnished to do so, subject 
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
        self.buffer=None
    def nextPred(self):
        if(self.buffer is not None):
            temp=self.buffer
            self.buffer=None
            return temp
        while(True):
            line=self.fh.readline()
            while(not rex.find("^GENE",line)): 
                line=self.fh.readline()
                if(line is None): return None
            fields=line.rstrip().split()
            (ID,Palt,theta,CI)=fields
            line=self.fh.readline()
            if(not rex.find("(\d)% : (\d)(\d) (\d)(\d) (\d)(\d)",line)):
                raise Exception("Can't parse: "+line)
            posterior=rex[1]; mother=[rex[2],rex[3]];
            father=[rex[4],rex[5]]; child=[rex[6],rex[7]]
            return [ID,Palt,theta,posterior,mother,father,child]
          
def getAffected(node,label):
    child=node.findChild(label)
    return [child[0],child[1]]

def isCorrect(trueMother,trueFather,trueChild,predMother,predFather,predChild):
    return trueMother==predMother and trueFather==predFather and \
        trueChild==predChild

def compact(mother,father,child):
    return mother[0]+mother[1]+" "+father[0]+father[1]+" "+child[0]+child[1]

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <truth.essex> <model-output>\n")
(truthFile,outputFile)=sys.argv[1:]

essexReader=EssexParser(truthFile)
predReader=BufferedReader(outputFile)
right=0; wrong=0
while(True):
    root=essexReader.nextElem()
    if(root is None): break
    prediction=predReader.nextPred()
    if(prediction is None): break
    (ID,Palt,theta,posterior,mother,father,child)=prediction
    if(ID!=root[0]): raise Exception(ID+" != "+root[0]);
    sxAffected=root.findChild("affected")
    trueMother=getAffected(sxAffected,"mother")
    trueFather=getAffected(sxAffected,"father")
    trueChild=getAffected(sxAffected,"child")
    if(isCorrect(trueMother,trueFather,trueChild,mother,father,child)):
        right+=1
    else:
        wrong+=1
       #print(ID,compact(trueMother,trueFather,trueChild),"<=>",
       #      compact(mother,father,child))
N=right+wrong
acc=float(right)/float(N)
print(round(acc*100,3),"correct")



