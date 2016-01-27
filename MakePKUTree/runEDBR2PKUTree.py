import os
import random
import sys
import ctypes
import time

import ROOT

from optparse import OptionParser
parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
(options, args) = parser.parse_args()


ROOT.gSystem.Load("EDBR2PKUTree_C.so")
from ROOT import EDBR2PKUTree, TFile, TTree



if __name__ == '__main__': 
    #input_dir="./BG_signal_samples/"
    input_dir="./2016_lowtrigger_samples/"

    infile=open("filelist.txt")

    totalsamples=infile.readlines()

    for sample in totalsamples:
        if sample[0]!="#":

            print "Loading ", sample
            sampleconfig=sample.split(';')
            samplename=sampleconfig[0]
            sampleXS=sampleconfig[1]
            sampleNum=sampleconfig[2]

            samplename=samplename.replace(' ', '').replace('\n','').replace('/n','') 
            sampleXS  =float(sampleXS.replace(' ', '').replace('\n','').replace('/n','') )
            sampleNum =int(sampleNum.replace(' ', '').replace('\n','').replace('/n','') )

            if (sampleXS-1.)==0: IsData=0
            elif "WJetsToLNu_HT" in samplename: IsData=2 
            else: IsData=1

            print samplename
            print sampleXS
            print sampleNum 
            print IsData

            filein =TFile(input_dir+samplename);
            dir1 = filein.Get("treeDumper");
            treein = dir1.Get("EDBRCandidates")
            #print "Open Tree Done"

            ##channel="el"
            ##if "singleMuon" in samplename: 
            ##    print "ignore mu data in el channel"
            ##else: 
            ##    analyzer_el=EDBR2PKUTree(treein,channel+"_out_"+samplename); 
            ##    analyzer_el.Loop(channel, sampleXS, sampleNum, IsData); 
            ##    analyzer_el.endJob();

            channel="mu"
            if "singleEl" in samplename: 
                print "ignore el data in mu channel"
            else:
                analyzer_mu=EDBR2PKUTree(treein,channel+"_out_"+samplename);
                analyzer_mu.Loop(channel, sampleXS, sampleNum, IsData);
                analyzer_mu.endJob();
