#! /usr/bin/env python
import os
import glob
import math
import array
import sys
import time

from array import array

import ROOT
from ROOT import gROOT, gStyle, gSystem, TLatex, TH1D, TString, TPaveText, TGaxis
import subprocess
from subprocess import Popen
from optparse import OptionParser

############################################
#             Job steering                 #
############################################

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

### to make the full analysis fit + datacards
parser.add_option('--makeCards', action='store_true', dest='makeCards', default=False, help='make datacards plus whole analysis')

### to compute limits
parser.add_option('--computeLimits', action='store_true', dest='computeLimits', default=False, help='compute limits')

### to plot limits
parser.add_option('--plotLimits', action='store_true', dest='plotLimits', default=False, help='plot limits')
parser.add_option('--lnubbBR', action='store_true', dest='lnubbBR', default=False, help='computing with munu BR')

### to do just signal lineshape fits
parser.add_option('--fitSignal', action='store_true', dest='fitSignal', default=False, help='do signal lineshape fits')

### other options 
parser.add_option('--channel',action="store",type="string",dest="channel",default="mu")
parser.add_option('--massPoint',action="store",type="int",dest="massPoint",default=-1)
parser.add_option('--cPrime',action="store",type="int",dest="cPrime",default=-1)
parser.add_option('--odir',action="store",type="string",dest="odir",default=".")
parser.add_option('--category',action="store",type="string",dest="category",default="HP") #"HP")
parser.add_option('--closuretest', action='store',type="int", dest='closuretest', default=0, help='closure test; 0: no test; 1: A1->A2; 2: A->B')
parser.add_option('--batchMode', action='store_true', dest='batchMode', default=False, help='no X11 windows')
parser.add_option('--limitMode', action='store',type="int", dest='limitMode', default=0, help='limit Mode; 0: Asymptotic ; 1: ProfileLikelihood ; 2: FullCLs ; 3: MaxLikelihoodFit')
parser.add_option('--isReReco', action='store',type="int", dest='isReReco', default=1, help='limit Mode; 0: Old signal samples ; 1: New signal Samples')
parser.add_option('--noSys', action='store',type="int", dest='noSys', default=0,help='run limit without systematic')
parser.add_option('--plotPvalue', action='store',type="int", default=0, dest='plotPvalue', help='plot p value')
parser.add_option('--signalWidth', action='store',type="int", default=0, dest='signalWidth', help='analysis on non-narrow signals')


(options, args) = parser.parse_args()

#########################################
### Global Variables for running jobs ###
#########################################

### mass point for signal to be fitted
#mass       = [1000,2000,3000,4000] 
#mass_width = [1000,2000,3000,4000]  
mass       = [600,700,750,800,900,1000]
mass_width = [600,700,750,800,900,1000]

bw_width   = [0.001,0.05,0.1,0.15,0.20,0.25,0.30,0.35,0.40] # ignore, no use in wh 

### mass window for couting analysis
ccmlo = [500,600,650,700, 800, 900]   # 
ccmhi = [700,800,850,900,1000,1100] #
ccmlo_width = [500,600,650,700, 800, 900]   # 
ccmhi_width = [700,800,850,900,1000,1100] #

### jet mass range
mjlo = [ 40, 40, 40, 40, 40, 40]
#mjhi = [ 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130]
mjhi = [ 130, 130, 130, 130, 130, 130]

mjlo_width = [ 40, 40, 40, 40, 40, 40]
#mjhi_width = [ 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130]
mjhi_width = [ 130, 130, 130, 130, 130, 130]

### mlvj range min and max used when run with option --makeCards
#fit range
mlo = [  400, 400, 400, 600, 600, 600]
mhi = [ 1000,1000,1000,1400,1400,1400]

mlo_width = [  400, 400, 400, 600, 600, 600]
mhi_width = [ 1000,1000,1000,1400,1400,1400]

### mlvj range min and max used when run with option --fitSignal
mlo_sig = [  400, 400, 400, 600, 600, 600]
mhi_sig = [ 1000,1000,1000,1400,1400,1400]
mlo_sig_width = [  400, 400, 400, 600, 600, 600]
mhi_sig_width = [ 1000,1000,1000,1400,1400,1400]

### shape to be used for bkg when --makeCards
shape    = ["ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1","Exp","Exp","Exp"]
shapeAlt = ["ErfPow2_v1"  ,"ErfPow2_v1"  ,"ErfPow2_v1"  ,"Pow","Pow","Pow"]

### shape to be used for bkg when --fitSignal

shape_sig_width  = ["BWDoubleCB","BWDoubleCB", "BWDoubleCB", "BWDoubleCB"]

#shape_sig_narrow = ["CB_v1","CB_v1","CB_v1","CB_v1","CB_v1"]

shape_sig_narrow = ["DoubleCB_v1","DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1"]


### signal mass fraction for non narrow samples
mass_fraction = [0.15,0.05,0.3,0.3]

#BRnew  = [00];
#cprime = [10];

##################################################
#  cross-sections for HVT and LH model  #
##################################################

#xsDict2 = {800:0.0004931,900:0.00022506,
#          1000:0.00011035,1100:0.000056883,1200:0.000030626,1300:0.000017003,
#          1400:0.0000097456,1500:0.0000056979,1600:0.0000034149,1700:0.0000020677,
#          1800:0.00000127,1900:7.9677e-07,2000:5.0345e-07,2100:3.198e-07,2200:2.0502e-07,
#          2300:1.324e-07,2400:8.6099e-08,2500:5.6338e-08,
#          2600:1.324e-08,2700:8.6099e-09,2800:5.6338e-09,2900:1.324e-09,3000:8.6099e-09}
## HVT model xsec*Br(W'->WH)
## HVT model xsec*Br(W'->WH)
#xsDict2 = {1000:2.37,1500:0.235,2000:0.04797,2500:0.0094,3000:0.00292,4000:0.0002397}
#xsDict2_munubb = {1000:2.37*0.4444,1500:0.235*0.4444,2000:0.04797*0.4444,2500:0.0094*0.4444,3000:0.00292*0.4444,4000:0.0002397*0.4444}
## LH model xsec*Br(W'->WH)
##xsDict = {1000:2.37,1500:0.235,2000:0.04797,2500:0.0094,3000:0.00292,4000:0.0002397}
##xsDict_munubb = {1000:2.37*0.4444,1500:0.235*0.4444,2000:0.04797*0.4444,2500:0.0094*0.4444,3000:0.00292*0.4444,4000:0.0002397*0.4444}
xsDict =  {
        600:    7.147913*1e-3, 
        700:    2.900665*1e-3, 
        750:    1.950000*1e-3,    
        800:    1.336325*1e-3, 
        900:    0.671447*1e-3, 
        1000:   0.360174*1e-3         
        }
xsDict_munubb =  {
        600:    7.147913*1e-3*0.3, 
        700:    2.900665*1e-3*0.3, 
        750:    1.950000*1e-3*0.3,    
        800:    1.336325*1e-3*0.3, 
        900:    0.671447*1e-3*0.3, 
        1000:   0.360174*1e-3*0.3         
        }

xsDict2 =  {
        600:    7.147913*1e-3, 
        700:    2.900665*1e-3, 
        750:    1.950000*1e-3,    
        800:    1.336325*1e-3, 
        900:    0.671447*1e-3, 
        1000:   0.360174*1e-3         
        }
xsDict2_munubb =  {
        600:    7.147913*1e-3*0.3, 
        700:    2.900665*1e-3*0.3, 
        750:    1.950000*1e-3*0.3,    
        800:    1.336325*1e-3*0.3, 
        900:    0.671447*1e-3*0.3, 
        1000:   0.360174*1e-3*0.3         
        }



################################
## style setup for doUL plots ##
################################
def setStyle():

  gStyle.SetPadBorderMode(0);
  gStyle.SetFrameBorderMode(0);
  gStyle.SetPadBottomMargin(0.12);
  gStyle.SetPadLeftMargin(0.12);
  gStyle.SetCanvasColor(ROOT.kWhite);
  gStyle.SetCanvasDefH(600); #Height of canvas
  gStyle.SetCanvasDefW(600); #Width of canvas
  gStyle.SetCanvasDefX(0);   #POsition on screen
  gStyle.SetCanvasDefY(0);

  gStyle.SetPadTopMargin(0.05);
  gStyle.SetPadBottomMargin(0.15);#0.13);
  gStyle.SetPadLeftMargin(0.15);#0.16);
  gStyle.SetPadRightMargin(0.05);#0.02);



 # For the Pad:
  gStyle.SetPadBorderMode(0);
  # gStyle.SetPadBorderSize(Width_t size = 1);
  gStyle.SetPadColor(ROOT.kWhite);
  gStyle.SetPadGridX(ROOT.kFALSE);
  gStyle.SetPadGridY(ROOT.kFALSE);
  gStyle.SetGridColor(0);
  gStyle.SetGridStyle(3);
  gStyle.SetGridWidth(1);

  # For the frame:
  gStyle.SetFrameBorderMode(0);
  gStyle.SetFrameBorderSize(1);
  gStyle.SetFrameFillColor(0);
  gStyle.SetFrameFillStyle(0);
  gStyle.SetFrameLineColor(1);
  gStyle.SetFrameLineStyle(1);
  gStyle.SetFrameLineWidth(1);

  gStyle.SetAxisColor(1, "XYZ");
  gStyle.SetStripDecimals(ROOT.kTRUE);
  gStyle.SetTickLength(0.03, "XYZ");
  gStyle.SetNdivisions(505, "XYZ");
  gStyle.SetPadTickX(1);  # To get tick marks on the opposite side of the frame
  gStyle.SetPadTickY(1);
  gStyle.SetGridColor(0);
  gStyle.SetGridStyle(3);
  gStyle.SetGridWidth(1);


  gStyle.SetTitleColor(1, "XYZ");
  gStyle.SetTitleFont(42, "XYZ");
  gStyle.SetTitleSize(0.05, "XYZ");
  # gStyle.SetTitleXSize(Float_t size = 0.02); # Another way to set the size?
  # gStyle.SetTitleYSize(Float_t size = 0.02);
  gStyle.SetTitleXOffset(1.15);#0.9);
  gStyle.SetTitleYOffset(1.3); # => 1.15 if exponents
  gStyle.SetLabelColor(1, "XYZ");
  gStyle.SetLabelFont(42, "XYZ");
  gStyle.SetLabelOffset(0.007, "XYZ");
  gStyle.SetLabelSize(0.045, "XYZ");

  gStyle.SetPadBorderMode(0);
  gStyle.SetFrameBorderMode(0);
  gStyle.SetTitleTextColor(1);
  gStyle.SetTitleFillColor(10);
  gStyle.SetTitleFontSize(0.05);




###########################################################
### submit fit jobs on condor batch -> makeCards option ###
###########################################################

def submitBatchJob( command, fn ):

    currentDir = os.getcwd();

    # create a dummy bash/sh -> give to condor the instruction
    outScript=open(fn+".sh","w");

    outScript.write('#!/bin/bash');
    outScript.write("\n"+'date');
    outScript.write("\n"+'source /uscmst1/prod/sw/cms/bashrc prod');
    outScript.write("\n"+'echo "condor dir: " ${_CONDOR_SCRATCH_DIR}');
    outScript.write("\n"+'cd '+currentDir);
    outScript.write("\n"+'eval `scram runtime -sh`');
    outScript.write("\n"+'cd -');
    outScript.write("\n"+'export PATH=${PATH}:'+currentDir);
    outScript.write("\n"+'echo ${PATH}');
    outScript.write("\n"+'ls');
    outScript.write("\n"+command);
    outScript.write("\n"+'tar -cvzf outputFrom_'+fn+'.tar.gz *');  ## pruce a tar file as output
    outScript.close();

    #### link a condor script to your shell script
    condorScript=open("condor_"+fn,"w");
    condorScript.write('universe = vanilla')
    condorScript.write("\n"+"Executable = "+fn+".sh")
    condorScript.write("\n"+'Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000')
    condorScript.write("\n"+'Should_Transfer_Files = YES')
    condorScript.write("\n"+'Transfer_Input_Files = g1_exo_doFit_class.py')
    condorScript.write("\n"+'WhenToTransferOutput = ON_EXIT_OR_EVICT')
    condorScript.write("\n"+'Output = out_$(Cluster).stdout')
    condorScript.write("\n"+'Error = out_$(Cluster).stderr')
    condorScript.write("\n"+'Error = out_$(Cluster).stderr')
    condorScript.write("\n"+'Log = out_$(Cluster).log')
    condorScript.write("\n"+'Notification = Error')
    condorScript.write("\n"+'Queue 1')
    condorScript.close();

    os.system("condor_submit "+"condor_"+fn)


##################################################################
### submit combine jobs on condor batch -> computeLimit option ###
##################################################################


#def submitBatchJobCombine( command, fn, channel, mass, cprime, BRnew, purity, isReReco, datacard =""):
def submitBatchJobCombine( command, fn, channel, mass, purity, isReReco, datacard =""):


    currentDir = os.getcwd();
    #### create a dummy bash/csh
    outScript=open(fn+".sh","w");

    outScript.write('#!/bin/bash');
    outScript.write("\n"+'date');
    outScript.write("\n"+'source /uscmst1/prod/sw/cms/bashrc prod');
    outScript.write("\n"+'echo "condor dir: " ${_CONDOR_SCRATCH_DIR}');
    outScript.write("\n"+'cd '+currentDir);
    outScript.write("\n"+'eval `scram runtime -sh`');
    outScript.write("\n"+'cd -');
    outScript.write("\n"+'export PATH=${PATH}:'+currentDir);
    outScript.write("\n"+'echo ${PATH}');
    outScript.write("\n"+'ls');
    outScript.write("\n"+command);
    if not purity=="combo" :
     outScript.write("\n"+'tar -cvzf outputFrom_'+fn+'.tar.gz *'+str(mass)+'*'+str(channel)+'*'+str(purity)+'*');
    else:
     outScript.write("\n"+'tar -cvzf outputFrom_'+fn+'.tar.gz *'+str(mass)+'*'+str(purity)+'*');
        
    outScript.close();

    if datacard == "" :
     if not purity=="combo" :
        if isReReco == 0:
	 file1 = "wwlvj_MWp_%03d_bb_%s_%s_unbin.txt"%(mass,channel,purity);
	 file2 = "wwlvj_MWp_%03d_bb_%s_%s_workspace.root"%(mass,channel,purity);
        else: 
	 file1 = "wwlvj_MWp_%03d_bb_%s_%s_unbin.txt"%(mass,channel,purity);
	 file2 = "wwlvj_MWp_%03d_bb_%s_%s_workspace.root"%(mass,channel,purity);
     else:
       if isReReco == 0 : 
	  file1 = "wwlvj_MWp_%03d_bb_%s_unbin.txt"%(mass, purity);
	  file2 = "wwlvj_MWp_%03d_bb_mu_HP_workspace.root"%(mass);
	  file3 = "wwlvj_MWp_%03d_bb_el_HP_workspace.root"%(mass);
	  file4 = "wwlvj_MWp_%03d_bb_mu_LP_workspace.root"%(mass);
	  file5 = "wwlvj_MWp_%03d_bb_el_LP_workspace.root"%(mass);
       else:
	  file1 = "wwlvj_MWp_%03d_bb_%s_unbin.txt"%(mass, purity);
	  file2 = "wwlvj_MWp_%03d_bb_mu_HP_workspace.root"%(mass);
	  file3 = "wwlvj_MWp_%03d_bb_el_HP_workspace.root"%(mass);
	  file4 = "wwlvj_MWp_%03d_bb_mu_LP_workspace.root"%(mass);
	  file5 = "wwlvj_MWp_%03d_bb_el_LP_workspace.root"%(mass);

    else:
        if isReReco == 0:
          file1 = datacard;
          datacard.ReplaceAll("unbin","workspace");
          datacard.ReplaceAll("txt","root");
          file2 = datacard
        else:
          file1 = datacard;
          datacard.ReplaceAll("unbin","workspace");
          datacard.ReplaceAll("txt","root");
          file2 = datacard
                                                      
    # link a condor script to your shell script
    condorScript=open("subCondor_"+fn,"w");
    condorScript.write('universe = vanilla')
    condorScript.write("\n"+"Executable = "+fn+".sh")
    condorScript.write("\n"+'Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000')
    condorScript.write("\n"+'Should_Transfer_Files = YES')
    if not purity=="combo" : condorScript.write("\n"+'Transfer_Input_Files = '+file1+', '+file2)
    else : condorScript.write("\n"+'Transfer_Input_Files = '+file1+', '+file2+', '+file3+', '+file4+', '+file5)
    condorScript.write("\n"+'WhenToTransferOutput = ON_EXIT_OR_EVICT')
                                                                    
    condorScript.write("\n"+'Output = out_$(Cluster).stdout')
    condorScript.write("\n"+'Error = out_$(Cluster).stderr')
    condorScript.write("\n"+'Error = out_$(Cluster).stderr')
    condorScript.write("\n"+'Log = out_$(Cluster).log')
    condorScript.write("\n"+'Notification = Error')
    condorScript.write("\n"+'Queue 1')
    condorScript.close();


    os.system("condor_submit "+"subCondor_"+fn)



##################################################
### Get Limit Value from Combine -M Asymptotic ###
##################################################

def getAsymLimits(file):
    
    
    f = ROOT.TFile(file);
    t = f.Get("limit");
    entries = t.GetEntries();
    
    lims = [0,0,0,0,0,0];
    
    for i in range(entries):
        
        t.GetEntry(i);
        t_quantileExpected = t.quantileExpected;
        t_limit = t.limit;
        
        #print "limit: ", t_limit, ", quantileExpected: ",t_quantileExpected;
        
        if t_quantileExpected == -1.: lims[0] = t_limit;
        elif t_quantileExpected >= 0.024 and t_quantileExpected <= 0.026: lims[1] = t_limit;
        elif t_quantileExpected >= 0.15  and t_quantileExpected <= 0.17:  lims[2] = t_limit;
        elif t_quantileExpected == 0.5: lims[3] = t_limit;
        elif t_quantileExpected >= 0.83  and t_quantileExpected <= 0.85:  lims[4] = t_limit;
        elif t_quantileExpected >= 0.974 and t_quantileExpected <= 0.976: lims[5] = t_limit;
        else: print "Unknown quantile!"
        print "entries: ", entries;
        print "obs: ", lims[0], ", 2s: ", lims[1], lims[1], ", 1s: ", lims[2], lims[4], ", exp: ", lims[3];
    
    return lims;


####################################################
### Get PValue from combine -M ProfileLikelihood ###
####################################################

def getPValueFromCard( file ):

    f = ROOT.TFile(file);
    t = f.Get("limit");
    entries = t.GetEntries();
    
    lims = 1;
    
    for i in range(1):
        
        t.GetEntry(i);
        lims = t.limit
    
    return lims;


##########################################
### Make Limit Plot --> Brazilian Plot ###
##########################################

def doULPlot( suffix ):

    
    xbins     = array('d', [])
    xbins_env = array('d', [])
    ybins_exp = array('d', [])
    ybins_obs = array('d', [])
    ybins_1s  = array('d', [])
    ybins_2s  = array('d', [])
    ybins_xs_02 = array('d', [])
    ybins_xs_05 = array('d', [])
    
    for i in range(len(mass)):
        curFile = "higgsCombine_lim_%03d%s.Asymptotic.mH%03d.root"%(mass[i],suffix,mass[i]);
	print "curFile: %s"%curFile;
	if options.lnubbBR:
	  sf = xsDict_munubb[mass[i]]; #*0.1057*0.577; # BR(W->munu)*BR(H->bb)
	  sf2 = xsDict2_munubb[mass[i]]; #*0.1057*0.577;
	else:
          sf = xsDict[mass[i]];
          sf2 = xsDict2[mass[i]];
        sf=sf;#*20.;
        sf2=sf2;#*20.; 
        curAsymLimits = getAsymLimits(curFile);
        xbins.append( mass[i] );
        xbins_env.append( mass[i] );
        ybins_exp.append( curAsymLimits[3]*sf2 );
        ybins_obs.append( curAsymLimits[0]*sf2 );
        ybins_2s.append( curAsymLimits[1]*sf2 );
        ybins_1s.append( curAsymLimits[2]*sf2 );
        ybins_xs_02.append(sf);#/20.); #*0.25);
        ybins_xs_05.append(sf2);#/20.); #*0.25);
    
    for i in range( len(mass)-1, -1, -1 ):
        curFile = "higgsCombine_lim_%03d%s.Asymptotic.mH%03d.root"%(mass[i],suffix,mass[i]);
        print "curFile: %s"%curFile;
        if options.lnubbBR:
          sf = xsDict_munubb[mass[i]]; #*0.1057*0.577; # BR(W->munu)*BR(H->bb)
          sf2 = xsDict2_munubb[mass[i]]; #*0.1057*0.577;
        else:
          sf = xsDict[mass[i]];
          sf2 = xsDict2[mass[i]];
        sf=sf;#/20.;
        sf2=sf2;#/20.;
        curAsymLimits = getAsymLimits(curFile);
        xbins_env.append( mass[i] );
        ybins_2s.append( curAsymLimits[5]*sf2 );
        ybins_1s.append( curAsymLimits[4]*sf2 );
    
    curGraph_exp = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_exp);
    curGraph_obs = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_obs);
    curGraph_xs_02 = ROOT.TGraph(nPoints,xbins,ybins_xs_02);
    curGraph_xs_05 = ROOT.TGraph(nPoints,xbins,ybins_xs_05);
    curGraph_1s = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_1s);
    curGraph_2s = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_2s);
    
    curGraph_obs.SetMarkerStyle(20);
    curGraph_obs.SetLineWidth(3);
    curGraph_obs.SetLineStyle(1);
    curGraph_obs.SetMarkerSize(1.6);
    curGraph_exp.SetMarkerSize(1.3);
    curGraph_exp.SetMarkerColor(ROOT.kBlack);

    curGraph_exp.SetLineStyle(2);
    curGraph_exp.SetLineWidth(3);
    curGraph_exp.SetMarkerSize(2);
    curGraph_exp.SetMarkerStyle(24);
    curGraph_exp.SetMarkerColor(ROOT.kBlack);

    curGraph_xs_02.SetLineStyle(ROOT.kDashed);
    curGraph_xs_02.SetFillStyle(3344);
    curGraph_xs_02.SetLineWidth(2);
    curGraph_xs_02.SetMarkerSize(2);
    curGraph_xs_02.SetLineColor(ROOT.kBlue);

    curGraph_xs_05.SetLineStyle(ROOT.kSolid);
    curGraph_xs_05.SetFillStyle(3344);
    curGraph_xs_05.SetLineWidth(2);
    curGraph_xs_05.SetMarkerSize(2);
    curGraph_xs_05.SetLineColor(ROOT.kRed);


    curGraph_1s.SetFillColor(ROOT.kGreen);
    curGraph_1s.SetFillStyle(1001);
    curGraph_1s.SetLineStyle(ROOT.kDashed);
    curGraph_1s.SetLineWidth(3);

    curGraph_2s.SetFillColor(ROOT.kYellow);
    curGraph_2s.SetFillStyle(1001);
    curGraph_2s.SetLineStyle(ROOT.kDashed);
    curGraph_2s.SetLineWidth(3);
    
        
    oneLine = ROOT.TF1("oneLine","1",799,1001);
    oneLine.SetLineColor(ROOT.kRed);
    oneLine.SetLineWidth(3);

    setStyle();

    can_SM = ROOT.TCanvas("can_SM","can_SM",630,600);
    
    if options.lnubbBR:
      hrl_SM = can_SM.DrawFrame(750,1e-3, 2550, 1); #0.0005,2550,1);
      hrl_SM.GetYaxis().SetTitle("#sigma_{95%} (pp #rightarrow G_{RS} #rightarrow munubb) (pb)");
    else:
      #hrl_SM = can_SM.DrawFrame(750,1e-4, 4100, 100);
      hrl_SM = can_SM.DrawFrame(550,1e-4, 1050, 1e2);
      hrl_SM.GetYaxis().SetTitle("#sigma_{95%} (pp #rightarrow G_{Bulk} #rightarrow WW) (pb)");
    hrl_SM.GetYaxis().SetTitleOffset(1.35);
    hrl_SM.GetYaxis().SetTitleSize(0.045);
    hrl_SM.GetYaxis().SetTitleFont(42);

    hrl_SM.GetXaxis().SetTitle("M_{G} (GeV)");
    hrl_SM.GetXaxis().SetTitleSize(0.045);
    hrl_SM.GetXaxis().SetTitleFont(42);

    hrl_SM.GetYaxis().SetNdivisions(505);
    can_SM.SetGridx(1);
    can_SM.SetGridy(1);

    curGraph_2s.Draw("F");
    curGraph_1s.Draw("Fsame");
    #curGraph_obs.Draw("PCsame");
    curGraph_exp.Draw("Csame");
    curGraph_xs_02.Draw("Csame");
#    curGraph_xs_05.Draw("Csame");
       
    leg2 = ROOT.TLegend(0.36,0.78,0.8,0.92);

    leg2.SetFillColor(0);
    leg2.SetShadowColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.027);

    #leg2.AddEntry(curGraph_obs,"Asympt. CL_{S} Observed","LP")
    leg2.AddEntry(curGraph_1s,"Asympt. CL_{S}  Expected #pm 1#sigma","LF")
    leg2.AddEntry(curGraph_2s,"Asympt. CL_{S}  Expected #pm 2#sigma","LF")
    if options.lnubbBR:
#      leg2.AddEntry(curGraph_xs_05,"HVT B(gv=3):xsec_{W'}*BR(W'#rightarrowWH)*BR(WH#rightarrowmunubb)", "L");  ##sigma_{TH} x BR(W' #rightarrow WH)","L")
      #leg2.AddEntry(curGraph_xs_02,"RS model","L");
      leg2.AddEntry(curGraph_xs_02,"Bulk model","L");
    else:
#      leg2.AddEntry(curGraph_xs_05,"HVT B(gv=3):xsec_{W'} * BR(W' #rightarrow WH)", "L");  ##sigma_{TH} x BR(W' #rightarrow WH)","L")
      #leg2.AddEntry(curGraph_xs_02,"RS model","L");
      leg2.AddEntry(curGraph_xs_02,"Bulk model","L");

    #xaxis2 = TGaxis(790.0,0.005,810,0.005,699.0,701.0,100,"UI");
    #xaxis2.Draw();
    #t800 =  TPaveText(0.135,0.108,0.18,0.108,"brNDC");
    #t800.SetFillColor(ROOT.kWhite);
    #t800.SetTextSize(0.045);
    #t800.SetTextAlign(11);
    #t800.SetTextFont(42);
    #t800.SetBorderSize(0);
    #t800.AddText("800");
    #t800.Draw();
                    


    ROOT.gPad.SetLogy();

    can_SM.Update();
    can_SM.RedrawAxis();
    can_SM.RedrawAxis("g");
    can_SM.Update();

    leg2.Draw();

    banner = TLatex(0.95, 0.96, "2.2 fb^{-1} (13 TeV)");
    banner.SetNDC(); banner.SetTextSize(0.038); banner.SetTextFont(42); banner.SetTextAlign(31); banner.SetLineWidth(2); banner.Draw();
    CMStext = TLatex(0.15,0.96,"CMS");
    CMStext.SetNDC(); CMStext.SetTextSize(0.041); CMStext.SetTextFont(61); CMStext.SetTextAlign(11); CMStext.SetLineWidth(2); CMStext.Draw();
    if suffix =="_el_HP" :
        Extratext = TLatex(0.241, 0.96, "Preliminary W#rightarrow e#nu");
        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
    if suffix =="_mu_HP" :
        Extratext = TLatex(0.241, 0.96, "Preliminary W#rightarrow #mu#nu");
        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
    if suffix =="_combo" :
        Extratext = TLatex(0.241, 0.96, "Preliminary e+#mu combined");
        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
        
        #label_sqrt = TPaveText(0.4,0.953,0.96,0.975, "brNDC");
        #label_sqrt.SetFillColor(ROOT.kWhite);
        #label_sqrt.SetBorderSize(0);
        #label_sqrt.SetTextSize(0.038);
        #label_sqrt.SetTextFont(62);
        #label_sqrt.SetTextAlign(31); # align right
        #label_sqrt.AddText("L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
        #label_sqrt.Draw();        

    os.system("mkdir -p %s/LimitResult/"%(os.getcwd()));
    os.system("mkdir -p %s/LimitResult/Limit_ExpTail/"%(os.getcwd()));

    can_SM.SaveAs("./LimitResult/Limit_ExpTail/Lim%s.png"%(suffix));
    can_SM.SaveAs("./LimitResult/Limit_ExpTail/Lim%s.pdf"%(suffix));
    can_SM.SaveAs("./LimitResult/Limit_ExpTail/Lim%s.root"%(suffix));
    can_SM.SaveAs("./LimitResult/Limit_ExpTail/Lim%s.C"%(suffix));


#################
### Main Code ###    
#################
    
if __name__ == '__main__':

    
    CHAN = options.channel;
    
    moreCombineOpts = "";

    ### Set the working directory
    if options.computeLimits or options.plotLimits:
	os.chdir("cards_allCats");    

    ### put in functionality to test just one mass point or just one cprime

    nMasses = len(mass);
    mLo = 0;
    mHi = nMasses;

    nMasses_width = len(mass_width);
    mLo_width = 0;
    mHi_width = nMasses_width;

    if options.massPoint > 0 and not options.signalWidth:
        curIndex = mass.index(options.massPoint);
        mLo = curIndex;
        mHi = curIndex+1;

    elif options.massPoint > 0 and options.signalWidth:
        curIndex = mass_width.index(options.massPoint);
        mLo_width = curIndex;
        mHi_width = curIndex+1;
        
    cpLo = 1;
    cpHi = 2;

    ### Make cards option analysis
    if options.makeCards:
     if options.signalWidth == 0:
        for i in range(mLo,mHi):
            print "##################################################";
            print "##################################################";
            print "############# R U N N I N G F I T S ##############";
            print "mass = ",mass[i];
            print "###################################################";
            print "###################################################";
                
            time.sleep(0.3);
            if options.isReReco == 0 :
             #command = "python g1_exo_doFit_class.py %s MWp_%03d_bb %02d %02d %02d %02d %02d %02d %s %s -b -m %01d --inPath %s --category %s --closuretest %01d"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], 1, os.getcwd(), options.category,options.closuretest);
             command = "python g1_exo_doFit_class.py %s BulkGravWW%03d %02d %02d %02d %02d %02d %02d %s %s -b -m %01d --inPath %s --category %s --closuretest %01d"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], 1, os.getcwd(), options.category,options.closuretest);
            else : 
             #command = "python g1_exo_doFit_class.py %s MWp_%03d_bb %02d %02d %02d %02d %02d %02d %s %s -b -m %01d --inPath %s --category %s --closuretest %01d"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], 1, os.getcwd(), options.category,options.closuretest);
             command = "python g1_exo_doFit_class.py %s BulkGravWW%03d %02d %02d %02d %02d %02d %02d %s %s -b -m %01d --inPath %s --category %s --closuretest %01d"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], 1, os.getcwd(), options.category,options.closuretest);

            print command ;
#####################if not running at lpc, don't use batchMode
            if options.batchMode :
             fn = "fitScript_%s_%03d_%s"%(options.channel,mass[i],options.category);
             submitBatchJob(command,fn);

            if not options.batchMode:
             #print command;
             os.system(command);

     elif options.signalWidth == 1:
      for k in range(len(mass_fraction)):
        for i in range(mLo_width,mHi_width):
            print "###################################################";
            print "###################################################";
            print "####### R U N N I N G F I T S   W I D T H #########";
            print "mass = ",mass_width[i]," mass_fraction = ",mass_fraction[k];
            print "##################################################";
            print "##################################################";
                
            time.sleep(0.3);
            if options.isReReco == 1 :
             if int(mass_width[i]*mass_fraction[k]) < 100:
              command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_M%03d_W%02d %02d %02d %02d %02d %02d %02d %s %s -b -m --inPath %s --category %s --closuretest %01d"%(CHAN, mass_width[i], int(mass_width[i]*mass_fraction[k]),ccmlo_width[i], ccmhi_width[i], mjlo_width[i], mjhi_width[i], mlo_width[i], mhi_width[i], shape[i], shapeAlt[i], os.getcwd(), options.category,options.closuretest);
             elif int(mass_width[i]*mass_fraction[k]) >= 100 and int(mass_width[i]*mass_fraction[k]) < 1000:
              command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_M%03d_W%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --inPath %s --category %s --closuretest %01d"%(CHAN, mass_width[i], int(mass_width[i]*mass_fraction[k]),ccmlo_width[i], ccmhi_width[i], mjlo_width[i], mjhi_width[i], mlo_width[i], mhi_width[i], shape[i], shapeAlt[i], os.getcwd(), options.category,options.closuretest);
             else:
              command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_M%03d_W%04d %02d %02d %02d %02d %02d %02d %s %s -b -m --inPath %s --category %s --closuretest %01d"%(CHAN, mass_width[i], int(mass_width[i]*mass_fraction[k]),ccmlo_width[i], ccmhi_width[i], mjlo_width[i], mjhi_width[i], mlo_width[i], mhi_width[i], shape[i], shapeAlt[i], os.getcwd(), options.category,options.closuretest);

            print command ;

            if options.batchMode :
             fn = "fitScript_%s_W%s_%03d_%s"%(options.channel,mass_width[i],int(mass_width[i]*mass_fraction[k]),options.category);
             submitBatchJob(command,fn);
            if not options.batchMode:
             os.system(command);

    ### Make signal fit option
    if options.fitSignal:
     if options.signalWidth == 0:
        for i in range(mLo,mHi):
                
            print "##############################################";
            print "##############################################";
            print "############# FIT SIGNAL SHAPES ##############";
            print "mass = ",mass[i];
            print "##############################################";
            print "##############################################";
               
            time.sleep(0.3);
            if options.isReReco == 0 :
              command = "python g1_exo_doFit_class.py %s MWp_%03d_bb %02d %02d %02d %02d %02d %02d %s %s -b -m --inPath %s --category %s --closuretest %01d --fitSignal 1"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo_sig[i], mhi_sig[i], shape_sig_narrow[i], shape_sig_width[i], os.getcwd(), options.category,options.closuretest);
            elif options.isReReco == 1 : 
              command = "python g1_exo_doFit_class.py %s MWp_%03d_bb %02d %02d %02d %02d %02d %02d %s %s -b -m --inPath %s --category %s --closuretest %01d --fitSignal 1"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo_sig[i], mhi_sig[i], shape_sig_narrow[i], shape_sig_width[i], os.getcwd(), options.category,options.closuretest);
          
            print command ;

            if options.batchMode :
             fn = "fitScript_signal_%s_%03d_%s"%(options.channel,mass[i],options.category);
             submitBatchJob(command,fn);

            if not options.batchMode:
             os.system(command);


     elif options.signalWidth == 1:
      for k in range(len(mass_fraction)):
        for i in range(mLo_width,mHi_width):
                
            print "###################################";
            print "###################################";
            print "####### FIT SIGNAL SHAPES #########";
            print "mass = ",mass_width[i]," mass_fraction = ",mass_fraction[k];
            print "###################################";
            print "###################################";
                
            time.sleep(0.3);
            if options.isReReco == 1 :
             if int(mass_width[i]*mass_fraction[k]) < 100:
              command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_M%03d_W%02d %02d %02d %02d %02d %02d %02d %s %s -b -m --inPath %s --category %s --closuretest %01d --fitSignal 1"%(CHAN, mass_width[i], int(mass_width[i]*mass_fraction[k]),ccmlo_width[i], ccmhi_width[i], mjlo_width[i], mjhi_width[i], mlo_sig_width[i], mhi_sig_width[i], shape_sig_narrow[i], shape_sig_width[i], os.getcwd(), options.category,options.closuretest);
             elif int(mass_width[i]*mass_fraction[k]) >= 100 and int(mass_width[i]*mass_fraction[k]) <= 1000:
              command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_M%03d_W%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --inPath %s --category %s --closuretest %01d --fitSignal 1"%(CHAN, mass_width[i], int(mass_width[i]*mass_fraction[k]),ccmlo_width[i], ccmhi_width[i], mjlo_width[i], mjhi_width[i], mlo_sig_width[i], mhi_sig_width[i], shape_sig_narrow[i], shape_sig_width[i], os.getcwd(), options.category,options.closuretest);
             else :
              command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_M%03d_W%04d %02d %02d %02d %02d %02d %02d %s %s -b -m --inPath %s --category %s --closuretest %01d --fitSignal 1"%(CHAN, mass_width[i], int(mass_width[i]*mass_fraction[k]),ccmlo_width[i], ccmhi_width[i], mjlo_width[i], mjhi_width[i], mlo_sig_width[i], mhi_sig_width[i], shape_sig_narrow[i], shape_sig_width[i], os.getcwd(), options.category,options.closuretest);

            print command ;

            if options.batchMode :
             fn = "fitScript_%s_W%s_%03d_%s"%(options.channel,mass_width[i],int(mass_width[i]*mass_fraction[k]),options.category);
             submitBatchJob(command,fn);
            if not options.batchMode:
             os.system(command);

    #######################################################################################################
                  
    ### Compute Limits
    if options.computeLimits:

     if options.channel != "em" :
        for i in range(mLo,mHi):
            print "##################"+str(mass[i])+"#####################";
            time.sleep(0.3);
            ### Asymptotic Limit + profileLikelihood to have an hint of the r value
            ## mu HP
            if options.limitMode==0:
              ##if options.isReReco == 0 and not options.noSys:
              ## runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkGravWW%03d_%s_HP_unbin.txt -v 2"%(mass[i],mass[i],"mu",mass[i],"mu");
              ##elif options.isReReco == 0 and options.noSys:
              ## runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkGravWW%03d_%s_HP_unbin.txt -v 2"%(mass[i],mass[i],"mu",mass[i],"mu");
              ##elif options.isReReco == 1 and not options.noSys:
              ## runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkGravWW%03d_%s_HP_unbin.txt -v 2"%(mass[i],mass[i],"mu",mass[i],"mu");
              ##elif options.isReReco == 1 and options.noSys:
              ## runCmmd2 = "combine -M Asymptotic -S 0 --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkGravWW%03d_%s_HP_unbin.txt -v 2"%(mass[i],mass[i],"mu",mass[i],"mu");

              runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkGravWW%03d_%s_HP_unbin.txt -v 2 -S %d"%(mass[i],mass[i],"mu",mass[i],"mu", options.noSys);

              if options.batchMode:
               fn = "combineScript_Asymptotic_%s_%03d%s_%s"%("mu",mass[i],"","HP");
               submitBatchJobCombine(runCmmd2,fn,"mu",mass[i],"HP",options.isReReco);
              else:  
               print runCmmd2;
               #os.system(runCmmd2);

              time.sleep(0.1);

     if options.channel == "em" :
        for i in range(mLo_width,mHi_width):

              for iwidth in bw_width :               
               if options.isReReco == 0 and not options.noSys:
                datacard = TString();
                datacard.Form("wwlvj_MWp_%03d_bb_%s_HP_unbin_W%.3f.txt"%(mass[i],"em",iwidth));
                datacard.ReplaceAll("0.","_");
                datacard.ReplaceAll("_txt","txt");                
                runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_bb_%s_HP_W_%s -d %s -v 2"%(mass[i],mass[i],"em",iwidth,datacard);
               elif options.isReReco == 0 and options.noSys:
                datacard = TString();
                datacard.Form("wwlvj_MWp_%03d_bb_%s_HP_unbin_W%.3f.txt"%(mass[i],"em",iwidth));
                datacard.ReplaceAll("0.","_");
                datacard.ReplaceAll("_txt","txt");                
                runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_bb_%s_HP_W_%s -d %s -v 2"%(mass[i],mass[i],"em",iwidth,datacard);
               elif options.isReReco == 1 and not  options.noSys:
                datacard = TString();
                datacard.Form("wwlvj_MWp_%03d_bb_%s_HP_unbin_W%.3f.txt"%(mass[i],"em",iwidth));
                datacard.ReplaceAll("0.","_");
                datacard.ReplaceAll("_txt",".txt");                
                runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_bb_%s_HP_%s -d %s -v 2"%(mass[i],mass[i],"em",iwidth,datacard);
               elif options.isReReco == 1 and options.noSys:
                datacard = TString();
                datacard.Form("wwlvj_MWp_%03d_bb_%s_HP_unbin_W%.3f.txt"%(mass[i],"em",iwidth));
                datacard.ReplaceAll("0.","_");
                datacard.ReplaceAll("_txt",".txt");                
                runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_bb_%s_HP_%s -d %s -v 2"%(mass[i],mass[i],"em",iwidth,datacard);
              
               print runCmmd2;
               if options.batchMode:
                fn = "combineScript_Asymptotic_%s_%03d%s_%s_W_%s"%("em",mass[i],"","HP",datacard);
                submitBatchJobCombine(runCmmd2,fn,"em",mass[i],"HP",options.isReReco,datacard);
               else:  
                os.system(runCmmd2);

               time.sleep(0.1);

     if options.channel != "em" :
      if options.limitMode == 1 :
        for i in range(mLo,mHi):            
            ### pvalue evaluation
                
              ## mu HP
              print "##################### mu HP #####################";   
              if options.isReReco == 0:   
               runCmmd = "combine -M ProfileLikelihood --significance --pvalue -m %03d -n _pval_obs_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt\n"%(mass[i],mass[i],"mu",mass[i],"mu");
               runCmmd += "combine -M ProfileLikelihood --significance --pvalue -m %03d -n _pval_exp_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt -t -1 --expectSignal=1 --toysFreq"%(mass[i],mass[i],"mu",mass[i],"mu");
              else: 
               runCmmd = "combine -M ProfileLikelihood --significance --pvalue -m %03d -n _pval_obs_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt\n"%(mass[i],mass[i],"mu",mass[i],"mu");
               runCmmd += "combine -M ProfileLikelihood --significance --pvalue -m %03d -n _pval_exp_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"mu",mass[i],"mu");
              
              if options.batchMode:
               fn = "combineScript_HybridNew_%s_%03d%s_%s"%("mu",mass[i],"","HP");
               submitBatchJobCombine( runCmmd, fn, "mu",mass[i],"HP",options.isReReco);
              else:  
               os.system(runCmmd);

              time.sleep(0.1);

     if options.channel == "em" :
      if options.limitMode == 1:
        for i in range(mLo_width,mHi_width):

              ## semplified for el+mu
              runCmmd ="";
              if options.isReReco == 0 :
               runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_obs_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt\n"%(mass[i],mass[i],"em",mass[i],"em");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"em",mass[i],"em");
              else: 
               runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_obs_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt\n"%(mass[i],mass[i],"em",mass[i],"em");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"em",mass[i],"em");
              
              if options.batchMode:
               fn = "combineScript_ProfileLikelihood_%s_%03d%s_%s"%("em",mass[i],"","HP");
               #submitBatchJobCombine(runCmmd,fn,"em",mass[i],cprime[0], BRnew[0],"HP",options.isReReco);
               submitBatchJobCombine(runCmmd,fn,"em",mass[i],"HP",options.isReReco);
              else:  
               os.system(runCmmd);

     if options.channel != "em" :
        for i in range(mLo,mHi):


            ### Full CLs -> naive setup
            if options.limitMode ==2:

             if options.channel != "em":

              ## mu HP
              runCmmd3="";
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_bb_%s_HP -d wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 \n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_bb_%s_HP -d wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.5 \n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_bb_%s_HP -d wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.16\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_bb_%s_HP -d wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.84\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_bb_%s_HP -d wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.025\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_bb_%s_HP -d wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.975\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");

                                                                     
              if options.batchMode:
               fn = "combineScript_HybridNew_%s_%03d%s_%s"%("mu",mass[i],"","HP");
               submitBatchJobCombine( runCmmd3,fn, "mu",mass[i],"HP" );
              else:  
               os.system(runCmmd3);

              time.sleep(0.1);

             elif options.channel == "em":

              ## combined
              runCmmd3="";
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_bb_%s_HP -d wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 \n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"em",mass[i],"em");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_bb_%s_HP -d wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.5\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"em",mass[i],"em");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_bb_%s_HP -d wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.16\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"em",mass[i],"em");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_bb_%s_HP -d wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.84\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"em",mass[i],"em");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_bb_%s_HP -d wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.025\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"em",mass[i],"em");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_bb_%s_HP -d wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.975\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"em",mass[i],"em");

                                                                     
              if options.batchMode:
               fn = "combineScript_HybridNew_%s_%03d%s"%("em",mass[i],"");
               submitBatchJobCombine(runCmmd3,fn, "",mass[i],"em");
              else:  
               os.system(runCmmd3);

            ### maximum likelihood  fit plus toys for S+B -> get bias estimation
            if options.limitMode == 3 :

             tmp_rMin=0.01; tmp_rMax=1000;
             if mass[i]>=2000: tmp_rMin=1000; tmp_rMax=50000;
             if mass[i]>=1500 and mass[i]<2000: tmp_rMin=100; tmp_rMax=10000;
             if mass[i]>=700 and mass[i]<1500:  tmp_rMin=10; tmp_rMax=5000;
             if mass[i]<700: tmp_rMin=0.001; tmp_rMax=1000;

             if options.channel != "em":
                
              ## mu HP
              if options.isReReco == 0:
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt -t 500 --expectSignal=600"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"mu",mass[i],"mu");
              else :
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt -t 500 --expectSignal=600"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"mu",mass[i],"mu");
                 
              print runCmmd
              if options.batchMode:
               fn = "combineScript_MaxLikelihoodFit_%s_%03d%s_%s"%("mu",mass[i],"","HP");
               submitBatchJobCombine(runCmmd,fn,"mu",mass[i],"HP",options.isReReco);
              else:  
                os.system(runCmmd);

              time.sleep(0.1);

              ## mu LP
              if options.isReReco == 0:
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_bb_%s_LP wwlvj_MWp_%03d_bb_%s_LP_unbin.txt -t 500 --expectSignal=2000 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"mu",mass[i],"mu");
              else:
               #runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_%03d_bb_%s_LP wwlvj_MWp_%03d_bb_%s_LP_unbin.txt\n"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"mu",mass[i],"mu");
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_bb_%s_LP wwlvj_MWp_%03d_bb_%s_LP_unbin.txt -t 500 --expectSignal=2000 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"mu",mass[i],"mu");
                 
              if options.batchMode:
               fn = "combineScript_MaxLikelihoodFit_%s_%03d%s_%s"%("mu",mass[i],"","LP");
               submitBatchJobCombine(runCmmd,fn,"mu",mass[i],"LP",options.isReReco);
              else:  
               os.system(runCmmd);

              time.sleep(0.1);

              ## el HP
              if options.isReReco == 0:
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt -t 500 --expectSignal=600 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"el",mass[i],"el");
              else:
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt -t 500 --expectSignal=600 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"el",mass[i],"el");
                 
              if options.batchMode:
               fn = "combineScript_MaxLikelihoodFit_%s_%03d%s_%s"%("el",mass[i],"","HP");
               submitBatchJobCombine(runCmmd,fn,"el",mass[i],"HP",options.isReReco);
              else:  
                os.system(runCmmd);

              time.sleep(0.1);

              ## el LP
              if options.isReReco == 0:
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_bb_%s_LP wwlvj_MWp_%03d_bb_%s_LP_unbin.txt -t 500 --expectSignal=2000 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"el",mass[i],"el");
              else:
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_bb_%s_LP wwlvj_MWp_%03d_bb_%s_LP_unbin.txt -t 500 --expectSignal=2000 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"el",mass[i],"el");
                 
              if options.batchMode:
               fn = "combineScript_MaxLikelihoodFit_%s_%03d%s_%s"%("el",mass[i],"","LP");
               submitBatchJobCombine(runCmmd,fn,"el",mass[i],"LP",options.isReReco);
              else :  
               os.system(runCmmd);

              time.sleep(0.1);

              ## combined
              if options.isReReco == 0:
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s wwlvj_MWp_%03d_bb_%s_unbin.txt -t 500 --expectSignal=350 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"combo",mass[i],"combo");
              else:
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s wwlvj_MWp_%03d_bb_%s_unbin.txt -t 500 --expectSignal=350 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"combo",mass[i],"combo");
                 
              if options.batchMode:
               fn = "combineScript_MaxLikelihoodFit_%s_%03d%s"%("combo",mass[i],"");
               submitBatchJobCombine(runCmmd,fn,"combo",mass[i],"combo",options.isReReco);
              else:  
               os.system(runCmmd);


             elif options.channel == "em":

              ## combined
              if options.isReReco == 0:
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s wwlvj_MWp_%03d_bb_%s_unbin.txt -t 500 --expectSignal=350 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"em",mass[i],"em");
             else:
              runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s wwlvj_MWp_%03d_bb_%s_unbin.txt -t 500 --expectSignal=350 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"em",mass[i],"em");
                 
              if options.batchMode:
               fn = "combineScript_MaxLikelihoodFit_%s_%03d%s"%("em",mass[i],"");
               submitBatchJobCombine(runCmmd,fn,"em",mass[i],"em",options.isReReco);
              else:  
               os.system(runCmmd);


                 
    ### make the plots    
    if options.plotLimits:

        #nGraphs = nCprimes*2 + 2;

        tGraphs = [];
        nPoints = len(mass);
        
        xbins = array('d', [])
        ybins_mu_allp = array('d', [])
        ybins_combo = array('d', [])

        #xbins_exp       = array('d', [])
        #ybins_mu_allp_exp = array('d', [])

        if options.channel != "em":

         for i in range(mLo,mHi):
             print "############################"+str(mass[i])+"#############################";            
             print "getting card: higgsCombine_pval_obs_%03d_bb_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);

             orootname_mu_allp = "higgsCombine_pval_obs_%03d_bb_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);

             #print "getting card: higgsCombine_pval_exp_%03d_bb_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);

             #orootname_mu_allp_exp = "higgsCombine_pval_exp_%03d_bb_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);

             xbins.append( mass[i] );
             #xbins_exp.append( mass[i] );

             if options.plotPvalue :

              ybins_mu_allp.append( getPValueFromCard(orootname_mu_allp) );

              #ybins_mu_allp_exp.append( getPValueFromCard(orootname_mu_allp_exp) );

         if options.plotPvalue : 

          gr_mu_allp = ROOT.TGraph(nPoints,xbins,ybins_mu_allp);
          gr_mu_allp.SetLineColor( 1 ); gr_mu_allp.SetMarkerColor( 1 ); gr_mu_allp.SetMarkerStyle( 20 ); gr_mu_allp.SetLineWidth( 3 );gr_mu_allp.SetMarkerSize( 1.6 );

          #gr_mu_allp_exp = ROOT.TGraph(nPoints,xbins_exp,ybins_mu_allp_exp);
          #gr_mu_allp_exp.SetLineColor( 2 ); gr_mu_allp_exp.SetMarkerColor( 2 ); gr_mu_allp_exp.SetMarkerStyle( 20 ); gr_mu_allp_exp.SetLineWidth( 3 );gr_mu_allp_exp.SetMarkerSize( 1.6 );
          #gr_mu_allp_exp.SetLineStyle(8);

          oneSLine = ROOT.TF1("oneSLine","1.58655253931457074e-01",800,3000);
          oneSLine.SetLineColor(ROOT.kRed); oneSLine.SetLineWidth(2); oneSLine.SetLineStyle(2);
          twoSLine = ROOT.TF1("twoSLine","2.27501319481792155e-02",800,3000);
          twoSLine.SetLineColor(ROOT.kRed); twoSLine.SetLineWidth(2); twoSLine.SetLineStyle(2);
          threeSLine = ROOT.TF1("threeSLine","1.34989803163009588e-03",800,3000);
          threeSLine.SetLineColor(ROOT.kRed); threeSLine.SetLineWidth(2); threeSLine.SetLineStyle(2);
          fourSLine = ROOT.TF1("fourSLine","3.16712418331199785e-05",800,3000);
          fourSLine.SetLineColor(ROOT.kRed); fourSLine.SetLineWidth(2); fourSLine.SetLineStyle(2);
    
          #banner = TLatex(0.43,0.91,("CMS Preliminary, 19.7 fb^{-1} at #sqrt{s}=8TeV"));
          #banner.SetNDC(); banner.SetTextSize(0.028);
          banner = TLatex(0.90, 0.91, "3. fb^{-1} (13 TeV)");
          banner.SetNDC(); banner.SetTextSize(0.038); banner.SetTextFont(42); banner.SetTextAlign(31); banner.SetLineWidth(2);
          CMStext = TLatex(0.1,0.91,"CMS");
          CMStext.SetNDC(); CMStext.SetTextSize(0.041); CMStext.SetTextFont(61); CMStext.SetTextAlign(11); CMStext.SetLineWidth(2);
          Extratext = TLatex(0.186,0.91,"Preliminary");
          Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2);
    
          ban1s = TLatex(2400,1.58655253931457074e-01,("1 #sigma"));
          ban1s.SetTextSize(0.028); ban1s.SetTextColor(2)
          ban2s = TLatex(2400,2.27501319481792155e-02,("2 #sigma"));
          ban2s.SetTextSize(0.028); ban2s.SetTextColor(2)
          ban3s = TLatex(2400,1.34989803163009588e-03,("3 #sigma"));
          ban3s.SetTextSize(0.028); ban3s.SetTextColor(2);
          ban4s = TLatex(2400,3.16712418331199785e-05,("4 #sigma"));
          ban4s.SetTextSize(0.028); ban4s.SetTextColor(2)
    
          leg3 = ROOT.TLegend(0.2,0.2,0.5,0.35);
          leg3.SetFillStyle(0);
          leg3.SetBorderSize(1);
          leg3.AddEntry( gr_mu_allp, "obs signif, #mu (HP)", "pl" );
          #leg3.AddEntry( gr_mu_allp_exp, "exp signif, HVT", "pl" );

          can2 = ROOT.TCanvas("can2","can2",800,800);
          hrl2 = can2.DrawFrame(799,1e-5,2500,0.6);
          hrl2.GetYaxis().SetTitle("p-value");
          hrl2.GetXaxis().SetTitle("mass (GeV)");
          can2.SetGrid();
          ROOT.gPad.SetLogy();
          gr_mu_allp.Draw("PL");
          #gr_mu_allp_exp.Draw("PLsame");
          oneSLine.Draw("same");
          twoSLine.Draw("same");
          threeSLine.Draw("same");
          fourSLine.Draw("same");
          banner.Draw();
       	  CMStext.Draw();
	  Extratext.Draw();
          leg3.Draw();
          ban1s.Draw();
          ban2s.Draw();
          ban3s.Draw();
          ban4s.Draw();
          can2.SaveAs("./LimitResult/Limit_ExpTail/pvals_mu_HP.pdf","pdf");
          can2.SaveAs("./LimitResult/Limit_ExpTail/pvals_mu_HP.png","png");
          can2.SaveAs("./LimitResult/Limit_ExpTail/pvals_mu_HP.root","root");
          can2.SaveAs("./LimitResult/Limit_ExpTail/pvals_mu_HP.C","C");

         doULPlot("_mu_HP");

         ####################################################

        elif options.channel == "em":

         for i in range(mLo,mHi):
             print "############################"+str(mass[i])+"#############################";            

             print "getting card: higgsCombine_pval_obs_%03d_bb_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"em",mass[i]);

             orootname_em_hp = "higgsCombine_pval_obs_%03d_bb_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"em",mass[i]);

             print "getting card: higgsCombine_pval_exp_%03d_bb_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"em",mass[i]);

             orootname_em_hp_exp = "higgsCombine_pval_exp_%03d_bb_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"em",mass[i]);

             xbins.append( mass[i] );
             xbins_exp.append( mass[i] );
             if options.plotPvalue :
              ybins_em_hp.append( getPValueFromCard(orootname_em_hp) );
              ybins_em_hp_exp.append( getPValueFromCard(orootname_em_hp_exp) );

         if options.plotPvalue :
          gr_em_hp = ROOT.TGraph(nPoints,xbins,ybins_em_hp);
          gr_em_hp.SetLineColor( 1 ); gr_em_hp.SetMarkerColor( 1 ); gr_em_hp.SetMarkerStyle( 20 ); gr_em_hp.SetLineWidth( 3 );gr_em_hp.SetMarkerSize( 1.6 );

          gr_em_hp_exp = ROOT.TGraph(nPoints,xbins_exp,ybins_em_hp_exp);
          gr_em_hp_exp.SetLineColor( 2 ); gr_em_hp_exp.SetMarkerColor( 2 ); gr_em_hp_exp.SetMarkerStyle( 20 ); gr_em_hp_exp.SetLineWidth( 3 );gr_em_hp_exp.SetMarkerSize( 1.6 );

          oneSLine = ROOT.TF1("oneSLine","1.58655253931457074e-01",800,3000);
          oneSLine.SetLineColor(ROOT.kRed); oneSLine.SetLineWidth(2); oneSLine.SetLineStyle(2);
          twoSLine = ROOT.TF1("twoSLine","2.27501319481792155e-02",800,3000);
          twoSLine.SetLineColor(ROOT.kRed); twoSLine.SetLineWidth(2); twoSLine.SetLineStyle(2);
          threeSLine = ROOT.TF1("threeSLine","1.34989803163009588e-03",800,3000);
          threeSLine.SetLineColor(ROOT.kRed); threeSLine.SetLineWidth(2); threeSLine.SetLineStyle(2);
          fourSLine = ROOT.TF1("fourSLine","3.16712418331199785e-05",800,3000);
          fourSLine.SetLineColor(ROOT.kRed); fourSLine.SetLineWidth(2); fourSLine.SetLineStyle(2);
    
          banner = TLatex(0.43,0.91,("CMS Preliminary, 19.7 fb^{-1} at #sqrt{s}=8TeV"));
          banner.SetNDC(); banner.SetTextSize(0.028);
    
          ban1s = TLatex(2400,1.58655253931457074e-01,("1 #sigma"));
          ban1s.SetTextSize(0.028); ban1s.SetTextColor(2)
          ban2s = TLatex(2400,2.27501319481792155e-02,("2 #sigma"));
          ban2s.SetTextSize(0.028); ban2s.SetTextColor(2)
          ban3s = TLatex(2400,1.34989803163009588e-03,("3 #sigma"));
          ban3s.SetTextSize(0.028); ban3s.SetTextColor(2);
          ban4s = TLatex(2400,3.16712418331199785e-05,("4 #sigma"));
          ban4s.SetTextSize(0.028); ban4s.SetTextColor(2)
    
          leg2 = ROOT.TLegend(0.2,0.2,0.5,0.35);
          leg2.SetFillStyle(0);
          leg2.SetBorderSize(1);
          leg2.AddEntry( gr_em_hp, "obs signif, e+#mu (HP)", "pl" );
          leg2.AddEntry( gr_em_hp_exp, "exp signif, #tilde{k}=0.5", "pl" );

    
          can = ROOT.TCanvas("can","can",800,800);
          hrl = can.DrawFrame(799,1e-3,2500,0.6);
          hrl.GetYaxis().SetTitle("p-value");
          hrl.GetXaxis().SetTitle("mass (GeV)");
          can.SetGrid();
          ROOT.gPad.SetLogy();
          gr_em_hp.Draw("PL");
          gr_em_hp_exp.Draw("PLsame");
          oneSLine.Draw("same");
          twoSLine.Draw("same");
          threeSLine.Draw("same");
          fourSLine.Draw("same");
          banner.Draw();
          leg2.Draw();
          ban1s.Draw();
          ban2s.Draw();
          ban3s.Draw();
          ban4s.Draw();
          can.SaveAs("~/LimitResult/Limit_ExpTail/pvals_em_HP.pdf","pdf");
          can.SaveAs("~/LimitResult/Limit_ExpTail/pvals_em_HP.png","png");
          can.SaveAs("~/LimitResult/Limit_ExpTail/pvals_em_HP.root","root");
          can.SaveAs("~/LimitResult/Limit_ExpTail/pvals_em_HP.C","C");

         doULPlot("_em_HP");
