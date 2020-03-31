#!/usr/bin/env python

import commands
import re
import os

import sys

def AddFiles( directory ):
    outSamples = re.split(r'\n',commands.getoutput("ls -1 "+directory))
    sampleColl = ""
    count = 0 
    for idx, sample in  enumerate(outSamples):
        sampleColl = sampleColl+" "+directory+"/"+sample
        if ((idx%100==0 and idx!=0) or idx==(len(outSamples)-1)):
            count = count+1
            hadd = "hadd -f  tree_"+directory+"_"+str(count)+".root "+sampleColl
            print hadd
            os.system(hadd)
            sampleColl = ""
    haddAll = "hadd -f  tree_"+directory+".root "+"tree_"+directory+"_*.root"
    print haddAll
    os.system(haddAll)
    remove = "rm tree_"+directory+"_*.root"
    os.system(remove)
#####################################################
#AddFiles("ttbar")
AddFiles("ttbar_test")
