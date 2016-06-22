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

parser.add_option('-b', action='store_true', dest='noX', default=True, help='no X11 windows')

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
#parser.add_option('--cPrime',action="store",type="int",dest="cPrime",default=-1)
parser.add_option('--odir',action="store",type="string",dest="odir",default=".")
parser.add_option('--category',action="store",type="string",dest="category",default="HP") #"HP")
parser.add_option('--closuretest', action='store',type="int", dest='closuretest', default=0, help='closure test; 0: no test; 1: A1->A2; 2: A->B')
#parser.add_option('--batchMode', action='store_true', dest='batchMode', default=False, help='no X11 windows')
parser.add_option('--limitMode', action='store',type="int", dest='limitMode', default=0, help='limit Mode; 0: Asymptotic ; 1: ProfileLikelihood ; 2: FullCLs ; 3: MaxLikelihoodFit')
#parser.add_option('--isReReco', action='store',type="int", dest='isReReco', default=1, help='limit Mode; 0: Old signal samples ; 1: New signal Samples')
parser.add_option('--noSys', action='store',type="int", dest='noSys', default=0,help='run limit without systematic')
#parser.add_option('--plotPvalue', action='store',type="int", default=0, dest='plotPvalue', help='plot p value')
#parser.add_option('--signalWidth', action='store',type="int", default=0, dest='signalWidth', help='analysis on non-narrow signals')


(options, args) = parser.parse_args()

#########################################
### Global Variables for running jobs ###
#########################################

##### mass point for signal to be fitted
##mass  = [600,700,750,800,900,1000]
##### mass window for couting analysis
##ccmlo = [500,600,650,700, 800, 900] 
##ccmhi = [700,800,850,900,1000,1100]
##### jet mass range
##mjlo = [  40,  40,  40,  40,  40,  40]
##mjhi = [ 150, 150, 150, 150, 150, 150]
##### mlvj range min and max used when run with option --makeCards
###fit range
##mlo = [  600, 600, 600, 600, 600, 600]
##mhi = [ 1400,1400,1400,1400,1400,1400]
##### shape to be used for bkg when --makeCards
##shape    = ["Exp","Exp","Exp","Exp","Exp","Exp"]
##shapeAlt = ["Pow","Pow","Pow","Pow","Pow","Pow"]
##### shape to be used for bkg when --fitSignal
##shape_sig_width  = ["BWDoubleCB" ,"BWDoubleCB" , "BWDoubleCB" , "BWDoubleCB" , "BWDoubleCB" ,  "BWDoubleCB"]
##shape_sig_narrow = ["DoubleCB_v1","DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1"]


##### mass point for signal to be fitted
##mass  = [1000,2000,2500,3000,4000]
##### mass nalysis
##ccmlo = [ 900,1900,2400,2900,3900] 
##ccmhi = [1100,2100,2600,3100,4100]
##### jet mass range
##mjlo = [   40,  40,  40,  40,  40]
##mjhi = [  150, 150, 150, 150, 150]
##### mlvj range min and max used when run with option --makeCards
###fit range
##mlo = [  800, 800, 800, 800, 800]
##mhi = [ 5000,5000,5000,5000,5000]
##### shape to be used for bkg when --makeCards
###shape    = [   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN"]
###shapeAlt = ["ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail"]
##shape    = ["Exp", "Exp", "Exp", "Exp", "Exp"]
##shapeAlt = ["Pow", "Pow", "Pow", "Pow", "Pow"]
##### shape to be used for bkg when --fitSignal
##shape_sig_width  = [ "BWDoubleCB", "BWDoubleCB" , "BWDoubleCB" , "BWDoubleCB" ,  "BWDoubleCB"]
##shape_sig_narrow = ["DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1"]

### mass point for signal to be fitted
mass  = [600,700,750,800,900,1000,2000,2500,3000,3500,4000,4500]
### mass window for couting analysis
ccmlo = [500,600,650,700, 800, 900,1900,2400,2900,3400,3900,4400 ] 
ccmhi = [700,800,850,900,1000,1100,2100,2600,3100,3600,4100,4600 ]
### jet mass range
mjlo = [  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40]
mjhi = [ 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150]
### mlvj range min and max used when run with option --makeCards
#fit range
mlo = [  600, 600, 600, 600, 600, 800, 800, 800, 800, 800, 800, 800]
mhi = [ 1400,1400,1400,1400,1400,5000,5000,5000,5000,5000,5000,5000]
### shape to be used for bkg when --makeCards
#shape    = ["Exp","Exp","Exp","Exp","Exp",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN"]
#shapeAlt = ["Pow","Pow","Pow","Pow","Pow","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail"]
shape    = ["Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp"]
shapeAlt = ["Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow"]
### shape to be used for bkg when --fitSignal
shape_sig_width  = ["BWDoubleCB" ,"BWDoubleCB" , "BWDoubleCB" , "BWDoubleCB" , "BWDoubleCB" ,  "BWDoubleCB", "BWDoubleCB" ,"BWDoubleCB" , "BWDoubleCB" , "BWDoubleCB" , "BWDoubleCB" ,  "BWDoubleCB"]
shape_sig_narrow = ["DoubleCB_v1","DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1","DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1", "DoubleCB_v1"]



### signal mass fraction for non narrow samples
#mass_fraction = [0.15,0.05,0.3,0.3]

##################################################
#  cross-sections for HVT and LH model  #
##################################################
xsDict =  {
        600:    406.830851*1e-3, 
        700:    165.094354*1e-3,
        750:    110.500282*1e-3,
        800:    76.0582930*1e-3, 
        900:    38.2161000*1e-3, 
        1000:   20.4996930*1e-3,
        1200: 6.841404   *1e-3,  
        1400: 2.625789675*1e-3,  
        1800: 0.49976082 *1e-3,  
        2000: 0.239518815*1e-3,  
        2500: 0.044885148*1e-3,  
        3000: 0.009824256*1e-3,  
        3500: 0.002409958*1e-3,  
        4000: 0.000622786*1e-3,  
        4500: 0.000168503*1e-3
        }
xsDict_munubb =  {
        600:    406.830851*1e-3*0.3, 
        700:    165.094354*1e-3*0.3, 
        750:    110.500282*1e-3*0.3,    
        800:    76.0582930*1e-3*0.3, 
        900:    38.2161000*1e-3*0.3, 
        1000:   20.4996930*1e-3*0.3,
        1200: 6.841404   *1e-3*0.3,  
        1400: 2.625789675*1e-3*0.3,  
        1800: 0.49976082 *1e-3*0.3,  
        2000: 0.239518815*1e-3*0.3,  
        2500: 0.044885148*1e-3*0.3,  
        3000: 0.009824256*1e-3*0.3,  
        3500: 0.002409958*1e-3*0.3,  
        4000: 0.000622786*1e-3*0.3,  
        4500: 0.000168503*1e-3*0.3
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
    ybins_th = array('d', [])
    
    for i in range(len(mass)):
        curFile = "higgsCombine_lim_%03d%s.Asymptotic.mH%03d.root"%(mass[i],suffix,mass[i]);
        print "curFile: %s"%curFile;
        if options.lnubbBR:
            sf = xsDict_munubb[mass[i]]; #*0.1057*0.577; # BR(W->munu)*BR(H->bb)
        else:
            sf = xsDict[mass[i]];
        curAsymLimits = getAsymLimits(curFile);
        xbins.append( mass[i] );
        xbins_env.append( mass[i] );
        ybins_exp.append( curAsymLimits[3]*sf );
        ybins_obs.append( curAsymLimits[0]*sf );
        ybins_2s.append( curAsymLimits[1]*sf );
        ybins_1s.append( curAsymLimits[2]*sf );
        ybins_th.append(sf);#/20.); #*0.25);
    
    for i in range( len(mass)-1, -1, -1 ):
        curFile = "higgsCombine_lim_%03d%s.Asymptotic.mH%03d.root"%(mass[i],suffix,mass[i]);
        print "curFile: %s"%curFile;
        if options.lnubbBR:
          sf = xsDict_munubb[mass[i]]; #*0.1057*0.577; # BR(W->munu)*BR(H->bb)
        else:
          sf = xsDict[mass[i]];
        curAsymLimits = getAsymLimits(curFile);
        xbins_env.append( mass[i] );
        ybins_2s.append( curAsymLimits[5]*sf );
        ybins_1s.append( curAsymLimits[4]*sf );

    nPoints = len(mass);
    curGraph_exp = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_exp);
    curGraph_obs = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_obs);
    curGraph_th = ROOT.TGraph(nPoints,xbins,ybins_th);
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

    #curGraph_th.SetLineStyle(ROOT.kDashed);
    curGraph_th.SetFillStyle(3344);
    curGraph_th.SetLineWidth(2);
    curGraph_th.SetMarkerSize(2);
    curGraph_th.SetLineColor(ROOT.kRed);


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
      #hrl_SM = can_SM.DrawFrame(550,1e-4, 1050, 1e2);
      #hrl_SM = can_SM.DrawFrame(600,1e-3, 1000, 1e2);
      hrl_SM = can_SM.DrawFrame(mass[0],1e-4, mass[nPoints-1], 1e2);
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
    curGraph_th.Draw("Csame");
       
    leg2 = ROOT.TLegend(0.36,0.78,0.8,0.92);

    leg2.SetFillColor(0);
    leg2.SetShadowColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.027);

    leg2.AddEntry(curGraph_obs,"Asympt. CL_{S} Observed","LP")
    leg2.AddEntry(curGraph_1s,"Asympt. CL_{S}  Expected #pm 1#sigma","LF")
    leg2.AddEntry(curGraph_2s,"Asympt. CL_{S}  Expected #pm 2#sigma","LF")
    if options.lnubbBR:
      #leg2.AddEntry(curGraph_th,"RS model","L");
      leg2.AddEntry(curGraph_th,"Bulk model","L");
    else:
      #leg2.AddEntry(curGraph_th,"RS model","L");
      leg2.AddEntry(curGraph_th,"Bulk model","L");

    ROOT.gPad.SetLogy();

    can_SM.Update();
    can_SM.RedrawAxis();
    can_SM.RedrawAxis("g");
    can_SM.Update();

    leg2.Draw();

    banner = TLatex(0.95, 0.96, "2.1 fb^{-1} (13 TeV)");
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

    if options.massPoint > 0:
        curIndex = mass.index(options.massPoint);
        mLo = curIndex;
        mHi = curIndex+1;

    cpLo = 1;
    cpHi = 2;

    ### Make cards option analysis
    if options.makeCards:
        for i in range(mLo,mHi):
            print "##################################################";
            print "##################################################";
            print "############# R U N N I N G F I T S ##############";
            print "mass = ",mass[i];
            print "###################################################";
            print "###################################################";
                
            time.sleep(0.3);
            command_makeCards = "python g1_exo_doFit_class.py %s BulkGravWW%03d %02d %02d %02d %02d %02d %02d %s %s -b -m %01d --inPath %s --category %s --closuretest %01d  > /dev/null 2>&1"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], 1, os.getcwd(), options.category,options.closuretest);
            print command_makeCards ;
            os.system(command_makeCards);
                 
    ### Compute Limits
    if options.computeLimits:

        for i in range(mLo,mHi):
            print "##################"+str(mass[i])+"#####################";
            time.sleep(0.3);
            ### Asymptotic Limit + profileLikelihood to have an hint of the r value
            if options.limitMode==0:
              runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkGravWW%03d_%s_HP_unbin.txt -v 2 -S %d"%(mass[i],mass[i],options.channel ,mass[i],options.channel, options.noSys);
            elif options.limitMode==1:
              runCmmd2 = "combine -M ProfileLikelihood --significance --pvalue -m %03d -n _pval_obs_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt\n"%(mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd2 += "combine -M ProfileLikelihood --significance --pvalue -m %03d -n _pval_exp_%03d_bb_%s_HP wwlvj_MWp_%03d_bb_%s_HP_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"mu",mass[i],"mu");
            print runCmmd2;
            #os.system(runCmmd2);
            time.sleep(0.1);

    ### make the plots    
    if options.plotLimits:
        doULPlot("_%s_HP"%(options.channel));
