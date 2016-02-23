#! /usr/bin/env python
import os
import glob
import math
import array
import ROOT
import ntpath
import sys
import subprocess
from subprocess import Popen
from optparse import OptionParser


from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack, TGraph,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHist,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite


############################################
#              Job steering                #
############################################

parser = OptionParser()

parser.add_option('-a', '--additioninformation',action="store",type="string",dest="additioninformation",default="EXO")
#parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-b', action='store_true', dest='noX', default=True, help='no X11 windows')
parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="mu")
#parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="el")

parser.add_option('-p', '--psmodel',action="store",type="string",dest="psmodel",default="pythia")

parser.add_option('-s','--simple', action='store', dest='simple', default=True, help='pre-limit in simple mode')
parser.add_option('-m','--multi', action='store', dest='multi', default=False, help='pre-limit in multi mode')
parser.add_option('--fitwtagger', action='store_true', dest='fitwtagger', default=False, help='fit wtagger jet in ttbar control sample')
parser.add_option('--fitwtaggersim', action='store_true', dest='fitwtaggersim', default=False, help='fit wtagger jet in ttbar control sample with mu and el samples simultaneously')
parser.add_option('--check', action='store_true', dest='check', default=True, help='check the workspace for limit setting')
parser.add_option('--control', action='store_true', dest='control', default=False, help='control plot')
parser.add_option('--fitsignal', action='store',type="int", dest='fitsignal', default=0, help='fit only signal lineshape with a chosen model')
parser.add_option('--closuretest', action='store',type="int", dest='closuretest', default=0, help='closure test; 0: no test; 1: A1->A2; 2: A->B')


parser.add_option('--cprime', action="store",type="int",dest="cprime",default=10)
parser.add_option('--BRnew', action="store",type="int",dest="BRnew",default=0)

parser.add_option('--inPath', action="store",type="string",dest="inPath",default=".")

parser.add_option('--category', action="store",type="string",dest="category",default="HP")
#parser.add_option('--category', action="store",type="string",dest="category",default="ALLP")


(options, args) = parser.parse_args()

ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/Util_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")

from ROOT import draw_error_band, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error, RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf,  RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

###############################
## doFit Class Implemetation ##
###############################

class doFit_wj_and_wlvj:

    def __init__(self, in_channel,in_signal_sample, in_mlvj_signal_region_min=650, in_mlvj_signal_region_max=850, in_mj_min=40, in_mj_max=150, in_mlvj_min=600., in_mlvj_max=1500., fit_model="ErfExp_v1", fit_model_alter="ErfPow_v1", input_workspace=None):

        self.setTDRStyle();

        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9) ;
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9) ;

        ### set the channel type --> electron or muon
        self.channel=in_channel;
        self.leg = TLegend(); 
        self.MODEL_4_mlvj=fit_model;
        self.MODEL_4_mlvj_alter=fit_model_alter;
                
        print "########################################################################################"
        print "######## define class: binning, variables, cuts, files and nuisance parameters ########"
        print "########################################################################################"

        ### Set the mj binning for plots
        self.BinWidth_mj=5.;

        ### Set the binning for mlvj plots as a function of the model
        if options.fitsignal == 1:
            self.BinWidth_mlvj=20.;
        elif self.MODEL_4_mlvj=="ErfPowExp_v1" or self.MODEL_4_mlvj=="ErfPow2_v1" or self.MODEL_4_mlvj=="ErfExp_v1":
            self.BinWidth_mlvj=50.;
        else:
            self.BinWidth_mlvj=50.;
           

        #narrow the BinWidth_mj and BinWidth_mlvj by a factor of 5. Because Higgs-Combination-Tools will generate a binned sample, so need the bin width narrow. So, as a easy selution, we will increase the bin-width by a factor of 5 when ploting m_j m_WW
        self.narrow_factor=1.;

        ## correct the binning of mj 
        self.BinWidth_mj=self.BinWidth_mj/self.narrow_factor;
        nbins_mj=int((in_mj_max-in_mj_min)/self.BinWidth_mj);
        in_mj_max=in_mj_min+nbins_mj*self.BinWidth_mj;
                   
        ## correct the binning of mlvj 
        self.BinWidth_mlvj=self.BinWidth_mlvj/self.narrow_factor;
        nbins_mlvj=int((in_mlvj_max-in_mlvj_min)/self.BinWidth_mlvj);
        in_mlvj_max=in_mlvj_min+nbins_mlvj*self.BinWidth_mlvj;

        ## define jet mass variable
        rrv_mass_j = RooRealVar("rrv_mass_j","Pruned jet mass",(in_mj_min+in_mj_max)/2.,in_mj_min,in_mj_max,"GeV");
        rrv_mass_j.setBins(nbins_mj);

        ## define invariant mass WW variable
        rrv_mass_lvj= RooRealVar("rrv_mass_lvj","M_{WW}",(in_mlvj_min+in_mlvj_max)/2.,in_mlvj_min,in_mlvj_max,"GeV");
        rrv_mass_lvj.setBins(nbins_mlvj);

        ## set the model used for the background parametrization
        self.MODEL_4_mlvj=fit_model;
        self.MODEL_4_mlvj_alter=fit_model_alter;

        ## create the workspace and import them
        if input_workspace is None:
            self.workspace4fit_ = RooWorkspace("workspace4fit_","workspace4fit_");
        else:
            self.workspace4fit_ = input_workspace;
        getattr(self.workspace4fit_,"import")(rrv_mass_j);
        getattr(self.workspace4fit_,"import")(rrv_mass_lvj);

        #prepare workspace for unbin-Limit -> just for the stuff on which running the limit 
        self.workspace4limit_ = RooWorkspace("workspace4limit_","workspace4limit_");

        ## different code operation mode -> just normal analysis
        if options.closuretest ==0:
            self.mj_sideband_lo_min = in_mj_min;
            self.mj_sideband_lo_max = 65;#65;
            self.mj_signal_min = 65;#65;
            self.mj_signal_max = 95;#85 or 105;
            self.mj_sideband_hi_min = 135;#105;
            self.mj_sideband_hi_max = in_mj_max;
        if options.closuretest ==1: ##closure test A1->A2
            self.mj_sideband_lo_min = in_mj_min;
            self.mj_sideband_lo_max = 55;
            self.mj_signal_min = 55;
            self.mj_signal_max = 60;#65;
            self.mj_sideband_hi_min = 105;
            self.mj_sideband_hi_max = in_mj_max;
        if options.closuretest == 2: #closure test A->B
            self.mj_sideband_lo_min = in_mj_min;
            self.mj_sideband_lo_max = 60;#65;
            self.mj_signal_min = 100;
            self.mj_signal_max = 115;
            self.mj_sideband_hi_min = 115;
            self.mj_sideband_hi_max = in_mj_max;

        ## zone definition in the jet mass 
        rrv_mass_j.setRange("sb_lo",self.mj_sideband_lo_min,self.mj_sideband_lo_max);
        rrv_mass_j.setRange("signal_region",self.mj_signal_min,self.mj_signal_max);
        rrv_mass_j.setRange("sb_hi",self.mj_sideband_hi_min,self.mj_sideband_hi_max);
        rrv_mass_j.setRange("sblo_to_sbhi",self.mj_sideband_lo_min,self.mj_sideband_hi_max);

        ## signal region definition in the mlvj variable in case of counting limit
        self.mlvj_signal_min = in_mlvj_signal_region_min
        self.mlvj_signal_max = in_mlvj_signal_region_max
        rrv_mass_lvj.setRange("signal_region",self.mlvj_signal_min,self.mlvj_signal_max);
        rrv_mass_lvj.setRange("high_mass",2500,in_mlvj_max);

        #prepare the data and mc files --> set the working directory and the files name
#        if TString(in_signal_sample).Contains("inclusive"):
#         self.file_Directory="AnaSigTree_new/";
#        else:
        #self.file_Directory="data_Jan12_1/";
        self.file_Directory="data_Feb1/";
            

        self.PS_model= options.psmodel
        
        self.signal_sample=in_signal_sample;

        #if options.closuretest == 0:
        #    self.file_data = (self.channel+"_PKUTree_data_xww.root");#keep blind!!!!
        #else:
        #    self.file_data = (self.channel+"_PKUTree_data_xww.root");#keep blind!!!!
        
        if options.control==1:
            self.file_data = (self.channel+"_PKUTree_15D.root");
        else: 
            #self.file_data = (self.channel+"_PKUTree_pdata.root");#keep blind!!!!
            self.file_data = (self.channel+"_PKUTree_15D.root");

        #self.file_pseudodata = (self.channel+"_PKUTree_allBkg_xww.root");#fake data
        self.file_signal     = (self.channel+"_PKUTree_%s.root"%(self.signal_sample));
        self.file_WJets0_mc  = (self.channel+"_PKUTree_WJetsPt180_xww.root");
        self.file_VV_mc      = (self.channel+"_PKUTree_VV_xww.root");# WW+WZ
        self.file_TTbar_mc   = (self.channel+"_PKUTree_TTBARpowheg_xww.root");
        self.file_STop_mc    = (self.channel+"_PKUTree_SingleTop_xww.root");

        ## event categorization as a function of the purity and the applied selection
        self.wtagger_label = options.category;
        
        if self.wtagger_label=="HP" :
            self.wtagger_cut_max=0.45;
            self.wtagger_cut_min=0.0;

        if self.wtagger_label=="LP":
            self.wtagger_cut_max=0.75;
            self.wtagger_cut_min=0.45;

        #if self.wtagger_label=="nocut":
        if self.wtagger_label=="ALLP":
            self.wtagger_cut_max=10000;

        self.categoryID=0;
        #self.categoryID=-1;
        if self.wtagger_label=="ALLP" and self.channel=="el": self.categoryID=0;
        if self.wtagger_label=="ALLP" and self.channel=="mu": self.categoryID=1;
        if self.wtagger_label=="LP" and self.channel=="el": self.categoryID=-2
        if self.wtagger_label=="HP" and self.channel=="el": self.categoryID=-1
        if self.wtagger_label=="LP" and self.channel=="mu": self.categoryID=2
        if self.wtagger_label=="HP" and self.channel=="mu": self.categoryID=1
            
        #medium wtagger_eff reweight between data and mc #Wtagger_forV SF have be add to ntuple weight;
        if self.channel=="el" and self.wtagger_label=="ALLP":
            self.rrv_htagger_eff_reweight_forT=RooRealVar("rrv_htagger_eff_reweight_forT","rrv_htagger_eff_reweight_forT",1); #1.11); #0.968);
            self.rrv_htagger_eff_reweight_forT.setError(0.); #0.05); #0.214);
            self.rrv_htagger_eff_reweight_forV=RooRealVar("rrv_htagger_eff_reweight_forV","rrv_htagger_eff_reweight_forV",1.0); #0.891);
            self.rrv_htagger_eff_reweight_forV.setError(0.);

        if self.channel=="mu" and self.wtagger_label=="ALLP":
            self.rrv_htagger_eff_reweight_forT=RooRealVar("rrv_htagger_eff_reweight_forT","rrv_htagger_eff_reweight_forT",1); #1.11); #1.31);
            self.rrv_htagger_eff_reweight_forT.setError(0.); #0.05); #0.212);
            self.rrv_htagger_eff_reweight_forV=RooRealVar("rrv_htagger_eff_reweight_forV","rrv_htagger_eff_reweight_forV",1.0); #1.277);
            self.rrv_htagger_eff_reweight_forV.setError(0.);

        if self.channel=="mu" and self.wtagger_label=="HP":
            #self.rrv_htagger_eff_reweight_forT=RooRealVar("rrv_htagger_eff_reweight_forT","rrv_htagger_eff_reweight_forT", 0.975);
            self.rrv_htagger_eff_reweight_forT=RooRealVar("rrv_htagger_eff_reweight_forT","rrv_htagger_eff_reweight_forT", 1.00);
            #self.rrv_htagger_eff_reweight_forT=RooRealVar("rrv_htagger_eff_reweight_forT","rrv_htagger_eff_reweight_forT", 0.85);
            self.rrv_htagger_eff_reweight_forT.setError(0.02*self.rrv_htagger_eff_reweight_forT.getVal());
            #self.rrv_htagger_eff_reweight_forV=RooRealVar("rrv_htagger_eff_reweight_forV","rrv_htagger_eff_reweight_forV",0.891);
            self.rrv_htagger_eff_reweight_forV=RooRealVar("rrv_htagger_eff_reweight_forV","rrv_htagger_eff_reweight_forV",1);
            self.rrv_htagger_eff_reweight_forV.setError(0.0717773);

        if self.channel=="el" and self.wtagger_label=="HP":
            #self.rrv_htagger_eff_reweight_forT=RooRealVar("rrv_htagger_eff_reweight_forT","rrv_htagger_eff_reweight_forT",1.011);#0.968);
            self.rrv_htagger_eff_reweight_forT=RooRealVar("rrv_htagger_eff_reweight_forT","rrv_htagger_eff_reweight_forT",1.00);
            self.rrv_htagger_eff_reweight_forT.setError(0.03*self.rrv_htagger_eff_reweight_forT.getVal());
            #self.rrv_htagger_eff_reweight_forV=RooRealVar("rrv_htagger_eff_reweight_forV","rrv_htagger_eff_reweight_forV",1.005);#0.891);
            self.rrv_htagger_eff_reweight_forV=RooRealVar("rrv_htagger_eff_reweight_forV","rrv_htagger_eff_reweight_forV",1.0);#0.891);
            self.rrv_htagger_eff_reweight_forV.setError(0.0717773);

        if self.channel=="mu" and self.wtagger_label=="LP":
            self.rrv_htagger_eff_reweight_forT=RooRealVar("rrv_htagger_eff_reweight_forT","rrv_htagger_eff_reweight_forT",1.31);
            self.rrv_htagger_eff_reweight_forT.setError(0.048103*self.rrv_htagger_eff_reweight_forT.getVal());
            self.rrv_htagger_eff_reweight_forV=RooRealVar("rrv_htagger_eff_reweight_forV","rrv_htagger_eff_reweight_forV",1.277);
            self.rrv_htagger_eff_reweight_forV.setError(0.303);

        if self.channel=="el" and self.wtagger_label=="LP":
            self.rrv_htagger_eff_reweight_forT=RooRealVar("rrv_htagger_eff_reweight_forT","rrv_htagger_eff_reweight_forT",1.39);
            self.rrv_htagger_eff_reweight_forT.setError(0.08*self.rrv_htagger_eff_reweight_forT.getVal());
            self.rrv_htagger_eff_reweight_forV=RooRealVar("rrv_htagger_eff_reweight_forV","rrv_htagger_eff_reweight_forV",1.277);
            self.rrv_htagger_eff_reweight_forV.setError(0.303);

        if self.channel=="em" and self.wtagger_label=="HP":
            self.rrv_htagger_eff_reweight_forT=RooRealVar("rrv_htagger_eff_reweight_forT","rrv_htagger_eff_reweight_forT", 0.971);
            self.rrv_htagger_eff_reweight_forT.setError(0.02*self.rrv_htagger_eff_reweight_forT.getVal());
            self.rrv_htagger_eff_reweight_forV=RooRealVar("rrv_htagger_eff_reweight_forV","rrv_htagger_eff_reweight_forV",0.891);
            self.rrv_htagger_eff_reweight_forV.setError(0.0717773);

        if self.channel=="em" and self.wtagger_label=="LP":
            self.rrv_htagger_eff_reweight_forT=RooRealVar("rrv_htagger_eff_reweight_forT","rrv_htagger_eff_reweight_forT", 1.34);
            self.rrv_htagger_eff_reweight_forT.setError(0.02*self.rrv_htagger_eff_reweight_forT.getVal());
            self.rrv_htagger_eff_reweight_forV=RooRealVar("rrv_htagger_eff_reweight_forV","rrv_htagger_eff_reweight_forV",1.277);
            self.rrv_htagger_eff_reweight_forV.setError(0.303);


        print "wtagger efficiency correction for Top sample: %s +/- %s"%(self.rrv_htagger_eff_reweight_forT.getVal(), self.rrv_htagger_eff_reweight_forT.getError());
        print "wtagger efficiency correction for V sample: %s +/- %s"%(self.rrv_htagger_eff_reweight_forV.getVal(), self.rrv_htagger_eff_reweight_forV.getError());


        #correct the W-jet mass peak difference between data and MC
        self.mean_shift=1.36;
        self.sigma_scale=1.0;#1.102;
        print " mean correction for the W peak : ",self.mean_shift," Resolution correction : ",self.sigma_scale
        
        #result files: The event number, parameters and error write into a txt file. The dataset and pdfs write into a root file
        if not os.path.isdir("cards_%s_%s_%s_g1"%(options.additioninformation, self.channel,self.wtagger_label)):
            os.system("mkdir cards_%s_%s_%s_g1"%(options.additioninformation, self.channel,self.wtagger_label));
        #self.rlt_DIR="cards_%s_%s_%s_g1/"%(options.additioninformation, self.channel,self.wtagger_label)
        #self.rlt_DIR="from_Jordan/Feb11unblindedTTbarSysMassShift/cards_em/"
        self.rlt_DIR="cards_allCats/"

        ## extra text file
        self.file_rlt_txt = self.rlt_DIR+"other_wwlvj_%s_%s_%s.txt"%(self.signal_sample,self.channel,self.wtagger_label)
        ## workspace for limit
        self.file_rlt_root = self.rlt_DIR+"wwlvj_%s_%s_%s_workspace.root"%(self.signal_sample,self.channel,self.wtagger_label)
        self.file_rlt_root1 = self.rlt_DIR+"wwlvj_%s_mu_%s_workspace.root"%(self.signal_sample,self.wtagger_label)
        self.file_rlt_root2 = self.rlt_DIR+"wwlvj_%s_el_%s_workspace.root"%(self.signal_sample,self.wtagger_label)
        ## datacard for the ubninned limit
        #self.file_datacard_unbin = self.rlt_DIR+"wwlvj_%s_%s_%02d_%02d_%s_unbin.txt"%(self.signal_sample,self.channel,options.cprime,options.BRnew,self.wtagger_label)
        self.file_datacard_unbin = self.rlt_DIR+"wwlvj_%s_%s_%s_unbin.txt"%(self.signal_sample,self.channel,self.wtagger_label)
        ## workspace for the binned limit
        #self.file_datacard_counting = self.rlt_DIR+"wwlvj_%s_%s_%02d_%02d_%s_counting.txt"%(self.signal_sample,self.channel,options.cprime,options.BRnew,self.wtagger_label)
        self.file_datacard_counting = self.rlt_DIR+"wwlvj_%s_%s_%s_counting.txt"%(self.signal_sample,self.channel,self.wtagger_label)
        
        self.file_out=open(self.file_rlt_txt,"w");
        self.file_out.write("Welcome:\n");
        self.file_out.close()
        self.file_out=open(self.file_rlt_txt,"a+");

        ## color palette 
        self.color_palet={ #color palet
            'data' : 1,
            'WJets' : 2,
            'VV' : 4,
            'STop' : 7,
            'TTbar' : 210,
            'ggH' : 1,
            'vbfH' : 12,
            'Signal': 1,
            'Uncertainty' : kBlack,
            'Other_Backgrounds' : kBlue
        }

        ## for basic selection         
        self.vpt_cut   = 200;
        self.MET_cut = 40;
        self.lpt_cut   = 40;
        if self.channel=="el":
            self.MET_cut= 80; 
            self.lpt_cut = 45;
        self.deltaPhi_METj_cut =2.0;

        # parameters of data-driven method to get the WJets background event number.
        self.number_WJets_insideband=-1;
        self.datadriven_alpha_WJets_unbin=-1;
        self.datadriven_alpha_WJets_counting=-1;

        ### uncertainty for datacard
        self.lumi_uncertainty    = 0.05;
        self.XS_STop_uncertainty = 0.30 ;
        self.XS_VV_uncertainty   = 0.25 ;
        self.XS_TTbar_uncertainty   = 0.12; #0.20 ;#qun
        self.XS_TTbar_NLO_uncertainty = 0.063 ;# from AN-12-368 table8
        self.XS_STop_NLO_uncertainty  = 0.05 ;# from AN-12-368 table8
        self.XS_VV_NLO_uncertainty    = 0.10 ;# from AN-12-368 table8

        #el and mu trigger and eff uncertainty, AN2012_368_v5 12.3
        self.lep_trigger_uncertainty = 0.01;
        #self.lep_eff_uncertainty     = 0.02; #qun
        #b tag scale uncertainty
        self.btag_scale_uncertainty  = 0.025;
        self.signal_btag_uncertainty = 0.019;#0.002;#qun
      
        if self.channel == "mu":
         self.signal_lepton_energy_scale_uncertainty = 0.02;#0.007 ;#qun
         self.signal_lepton_energy_res_uncertainty   = 0.0025;#0.001 ;#qun
         self.signal_jet_energy_res_uncertainty      = 0.015;#0.003 ; #qun
         self.lep_eff_uncertainty     = 0.01;
        elif self.channel == "el":
         self.signal_lepton_energy_scale_uncertainty = 0.01;#0.002 ; #qun
         self.signal_lepton_energy_res_uncertainty   = 0.004;#0.001 ;#qun
         self.signal_jet_energy_res_uncertainty      = 0.015;#0.003 ;#qun
         self.lep_eff_uncertainty     = 0.03;
        elif self.channel == "em":
         self.signal_lepton_energy_scale_uncertainty = 0.045 ;
         self.signal_lepton_energy_res_uncertainty   = 0.001 ;
         self.signal_jet_energy_res_uncertainty      = 0.003 ;
         self.lep_eff_uncertainty     = 0.02;

        self.eff_vtag_model = 0. ;
        self.eff_vtag_mass_sf = 0.;

        label_tstring=TString(self.signal_sample);
        if label_tstring.Contains("600") and (not label_tstring.Contains("1600")):
            self.signal_jet_energy_scale_uncertainty = 0.01 ;
            self.xs_rescale = 0.32298/0.052087 ;            
        if label_tstring.Contains("700") and (not label_tstring.Contains("1700")):
            self.signal_jet_energy_scale_uncertainty = 0.011 ;
            self.xs_rescale = 0.11827/0.019006 ;            
        if label_tstring.Contains("750") and (not label_tstring.Contains("1750")):#tmp use 700 results
            self.signal_jet_energy_scale_uncertainty = 0.011 ;
            self.xs_rescale = 0.11827/0.019006 ;            
        if label_tstring.Contains("800") and (not label_tstring.Contains("1800")):
            self.signal_jet_energy_scale_uncertainty =  0.0068 ;#0.011 ;#qun
            self.xs_rescale = 0.04931/0.0079064 ;            
        if label_tstring.Contains("900") and (not label_tstring.Contains("1900")):
            self.signal_jet_energy_scale_uncertainty = 0.011 ;
            self.xs_rescale = 0.022506/0.0036364 ;            
        if label_tstring.Contains("1000"):
            self.signal_jet_energy_scale_uncertainty = 0.0214 ;#0.011 ;#qun
            self.xs_rescale = 0.011035/0.0017742;            
            self.signal_btag_uncertainty = 0.0;#qun
        if label_tstring.Contains("1100"):
            self.signal_jet_energy_scale_uncertainty = 0.014 ;
            self.xs_rescale = 0.0056883/0.00091785;            
        if label_tstring.Contains("1200"):
            self.signal_jet_energy_scale_uncertainty = 0.014 ;
            self.xs_rescale = 0.0030626/0.00049262;            
        if label_tstring.Contains("1300"):
            self.signal_jet_energy_scale_uncertainty = 0.014 ;
            self.xs_rescale = 0.0017003/0.00027418;            
        if label_tstring.Contains("1400"):
            self.signal_jet_energy_scale_uncertainty = 0.014 ;
            self.xs_rescale = 0.00097456/0.00015697;            
        if label_tstring.Contains("1500"):
            self.signal_jet_energy_scale_uncertainty = 0.015 ;
            self.xs_rescale = 0.00056979/9.2073e-05;            
        if label_tstring.Contains("1600"):
            self.signal_jet_energy_scale_uncertainty = 0.015 ;
            self.xs_rescale = 0.00034149/5.4715e-05;            
        if label_tstring.Contains("1700"):
            self.signal_jet_energy_scale_uncertainty = 0.016 ;
            self.xs_rescale = 0.00020677/3.3199e-05;            
        if label_tstring.Contains("1800"):
            self.signal_jet_energy_scale_uncertainty = 0.016 ;
            self.xs_rescale = 0.000127/2.0367e-05;            
        if label_tstring.Contains("1900"):
            self.signal_jet_energy_scale_uncertainty = 0.018 ;
            self.xs_rescale = 7.9677e-05/1.2723e-05;            
        if label_tstring.Contains("2000"):
            self.signal_jet_energy_scale_uncertainty = 0.0284 ;#0.018 ;#qun
            self.signal_btag_uncertainty = 0.0;
            self.xs_rescale = 5.0345e-05/8.0046e-06;            
        if label_tstring.Contains("2100"):
            self.signal_jet_energy_scale_uncertainty = 0.02 ;
            self.xs_rescale = 3.198e-05/5.0566e-06;            
        if label_tstring.Contains("2200"):
            self.signal_jet_energy_scale_uncertainty = 0.02 ;
            self.xs_rescale = 2.0502e-05/3.2608e-06;            
        if label_tstring.Contains("2300"):
            self.signal_jet_energy_scale_uncertainty = 0.023 ;
            self.xs_rescale = 1.324e-05/2.0938e-06;            
        if label_tstring.Contains("2400"):
            self.signal_jet_energy_scale_uncertainty = 0.026 ;
            self.xs_rescale = 8.6099e-06/1.3566e-06;            
        if label_tstring.Contains("2500"):
            self.signal_jet_energy_scale_uncertainty = 0.03 ;
            self.xs_rescale = 5.6338e-06/8.8518e-07;            
        if label_tstring.Contains("3000"):
            self.signal_jet_energy_scale_uncertainty = 0.0464 ;
            self.signal_btag_uncertainty = 0.07146;
        if label_tstring.Contains("4000"):
            self.signal_jet_energy_scale_uncertainty = 0.0464 ;
            self.signal_btag_uncertainty = 0.07146;

        #what is xs_rescale?
        self.xs_rescale = 1.;            

        #### sigma and mean signal systematic inflation
        self.mean_signal_uncertainty_jet_scale  = 0.005; #0.013 ;#qun
        self.mean_signal_uncertainty_lep_scale  = 0.005; #0.001 ;#qun
        self.sigma_signal_uncertainty_jet_scale = 0.0375; #0.033 ;#qun
        self.sigma_signal_uncertainty_jet_res   = 0.0325; #0.030 ;#qun
        self.sigma_signal_uncertainty_lep_scale = 0.0175; #0.005 ;#qun

        #### Set systematic on the Wjets shape   and TTbar due to PS, fitting function etc..
        self.shape_para_error_WJets0 = 1.4;
        self.shape_para_error_alpha  = 1.4;
        self.shape_para_error_TTbar = 2.0;
        self.shape_para_error_VV    = 1.;
        self.shape_para_error_STop  = 1.;
                                                                
        # shape parameter uncertainty
        self.FloatingParams=RooArgList("floatpara_list");
        
        self.nPV_min=0;
        self.nPV_max=200;

    ## Set basic TDR style for canvas, pad ..etc ..
    def setTDRStyle(self):
        self.tdrStyle =TStyle("tdrStyle","Style for P-TDR");
        #For the canvas:
        self.tdrStyle.SetCanvasBorderMode(0);
        self.tdrStyle.SetCanvasColor(kWhite);
        self.tdrStyle.SetCanvasDefH(600); #Height of canvas
        self.tdrStyle.SetCanvasDefW(600); #Width of canvas
        self.tdrStyle.SetCanvasDefX(0); #POsition on screen
        self.tdrStyle.SetCanvasDefY(0);
      
        #For the Pad:
        self.tdrStyle.SetPadBorderMode(0);
        self.tdrStyle.SetPadColor(kWhite);
        self.tdrStyle.SetPadGridX(False);
        self.tdrStyle.SetPadGridY(False);
        self.tdrStyle.SetGridColor(0);
        self.tdrStyle.SetGridStyle(3);
        self.tdrStyle.SetGridWidth(1);
      
        #For the frame:
        self.tdrStyle.SetFrameBorderMode(0);
        self.tdrStyle.SetFrameBorderSize(1);
        self.tdrStyle.SetFrameFillColor(0);
        self.tdrStyle.SetFrameFillStyle(0);
        self.tdrStyle.SetFrameLineColor(1);
        self.tdrStyle.SetFrameLineStyle(1);
        self.tdrStyle.SetFrameLineWidth(1);
      
        #For the histo:
        self.tdrStyle.SetHistLineColor(1);
        self.tdrStyle.SetHistLineStyle(0);
        self.tdrStyle.SetHistLineWidth(1);
        self.tdrStyle.SetEndErrorSize(2);
        self.tdrStyle.SetErrorX(0.);
        self.tdrStyle.SetMarkerStyle(20);
      
        #For the fit/function:
        self.tdrStyle.SetOptFit(1);
        self.tdrStyle.SetFitFormat("5.4g");
        self.tdrStyle.SetFuncColor(2);
        self.tdrStyle.SetFuncStyle(1);
        self.tdrStyle.SetFuncWidth(1);
      
        #For the date:
        self.tdrStyle.SetOptDate(0);
      
        #For the statistics box:
        self.tdrStyle.SetOptFile(0);
        self.tdrStyle.SetOptStat(0); #To display the mean and RMS:
        self.tdrStyle.SetStatColor(kWhite);
        self.tdrStyle.SetStatFont(42);
        self.tdrStyle.SetStatFontSize(0.025);
        self.tdrStyle.SetStatTextColor(1);
        self.tdrStyle.SetStatFormat("6.4g");
        self.tdrStyle.SetStatBorderSize(1);
        self.tdrStyle.SetStatH(0.1);
        self.tdrStyle.SetStatW(0.15);
      
        #Margins:
        self.tdrStyle.SetPadTopMargin(0.05);
        self.tdrStyle.SetPadBottomMargin(0.13);
        self.tdrStyle.SetPadLeftMargin(0.18);
        self.tdrStyle.SetPadRightMargin(0.06);
      
        #For the Global title:
        self.tdrStyle.SetOptTitle(0);
        self.tdrStyle.SetTitleFont(42);
        self.tdrStyle.SetTitleColor(1);
        self.tdrStyle.SetTitleTextColor(1);
        self.tdrStyle.SetTitleFillColor(10);
        self.tdrStyle.SetTitleFontSize(0.05);
      
        #For the axis titles:
        self.tdrStyle.SetTitleColor(1, "XYZ");
        self.tdrStyle.SetTitleFont(42, "XYZ");
        self.tdrStyle.SetTitleSize(0.03, "XYZ");
        self.tdrStyle.SetTitleXOffset(0.9);
        self.tdrStyle.SetTitleYOffset(1.5);
      
        #For the axis labels:
        self.tdrStyle.SetLabelColor(1, "XYZ");
        self.tdrStyle.SetLabelFont(42, "XYZ");
        self.tdrStyle.SetLabelOffset(0.007, "XYZ");
        self.tdrStyle.SetLabelSize(0.03, "XYZ");
      
        #For the axis:
        self.tdrStyle.SetAxisColor(1, "XYZ");
        self.tdrStyle.SetStripDecimals(kTRUE);
        self.tdrStyle.SetTickLength(0.03, "XYZ");
        self.tdrStyle.SetNdivisions(510, "XYZ");
        self.tdrStyle.SetPadTickX(1); #To get tick marks on the opposite side of the frame
        self.tdrStyle.SetPadTickY(1);
      
        #Change for log plots:
        self.tdrStyle.SetOptLogx(0);
        self.tdrStyle.SetOptLogy(0);
        self.tdrStyle.SetOptLogz(0);
      
        #Postscript options:
        self.tdrStyle.SetPaperSize(20.,20.);
        self.tdrStyle.cd();

    def read_2workspaces(self, logy=0):

        ### Taket the workspace for limits  
        file1 = TFile(self.file_rlt_root1) ;
        file2 = TFile(self.file_rlt_root2) ;
        workspace1 = file1.Get("workspace4limit_") ;
        workspace2 = file2.Get("workspace4limit_") ;
        workspace1.Print()
        workspace2.Print()

        ##### iterate on the workspace element parameters
        ##print "----------- Parameter Workspace -------------";
        ##parameters_workspace = workspace.allVars();
        ##par = parameters_workspace.createIterator();
        ##par.Reset();
        ##param = par.Next()
        ##while (param):
        ##    param.Print();
        ##    param=par.Next()
        ##print "---------------------------------------------";

        ##workspace1.data("data_obs_xww_mu_%s"%(self.wtagger_label)).Print()
        ##workspace2.data("data_obs_xww_el_%s"%(self.wtagger_label)).Print()

        ##print "----------- Pdf in the Workspace -------------";
        ##pdfs_workspace = workspace.allPdfs();
        ##par = pdfs_workspace.createIterator();
        ##par.Reset();
        ##param=par.Next()
        ##while (param):
        ##    param.Print();
        ##    param = par.Next()
        ##print "----------------------------------------------";

        rrv_x = workspace1.var("rrv_mass_lvj")
        data_obs  = workspace1.data("data_obs_xww_mu_%s"%(self.wtagger_label));
        data_obs2 = workspace2.data("data_obs_xww_el_%s"%(self.wtagger_label));
        data_obs.Print()
        data_obs2.Print()
        #data_obs  = RooDataSet#("data_obs_xww_em_%s"%(self.wtagger_label),"data_obs_xww_em_%s"%(self.wtagger_label), rrv_x,);
        #data_obs.append(data_obs1)
        #data_obs.append(data_obs1)
        data_obs.append(data_obs2)
        data_obs.Print()
        #data_obs.append(workspace1.data("data_obs_xww_mu_%s"%(self.wtagger_label)))
        #data_obs.append(workspace2.data("data_obs_xww_el_%s"%(self.wtagger_label)))
        if TString(self.signal_sample).Contains("Bulk_WW"):
            model_pdf_signal1 = workspace1.pdf("BulkWW_xww_mu_%s"%(self.wtagger_label));
            model_pdf_signal2 = workspace2.pdf("BulkWW_xww_el_%s"%(self.wtagger_label));
        else: 
            model_pdf_signal1 = workspace1.pdf("%s_xww_mu_%s"%(self.signal_sample,self.wtagger_label));
            model_pdf_signal2 = workspace2.pdf("%s_xww_el_%s"%(self.signal_sample,self.wtagger_label));
            
        model_pdf_WJets1  = workspace1.pdf("WJets_xww_mu_%s"%(self.wtagger_label));
        model_pdf_VV1     = workspace1.pdf("VV_xww_mu_%s"%(self.wtagger_label));
        model_pdf_TTbar1  = workspace1.pdf("TTbar_xww_mu_%s"%(self.wtagger_label));
        model_pdf_STop1   = workspace1.pdf("STop_xww_mu_%s"%(self.wtagger_label));

        model_pdf_WJets2  = workspace2.pdf("WJets_xww_el_%s"%(self.wtagger_label));
        model_pdf_VV2     = workspace2.pdf("VV_xww_el_%s"%(self.wtagger_label));
        model_pdf_TTbar2  = workspace2.pdf("TTbar_xww_el_%s"%(self.wtagger_label));
        model_pdf_STop2   = workspace2.pdf("STop_xww_el_%s"%(self.wtagger_label));

        model_pdf_signal1.Print();
        model_pdf_WJets1.Print();
        model_pdf_VV1.Print();
        model_pdf_TTbar1.Print();
        model_pdf_STop1.Print();

        model_pdf_signal2.Print();
        model_pdf_WJets2.Print();
        model_pdf_VV2.Print();
        model_pdf_TTbar2.Print();
        model_pdf_STop2.Print();

        if TString(self.signal_sample).Contains("Bulk_WW"): 
            rrv_number_signal1 = workspace1.var("rate_BulkWW_xww_for_unbin");
            rrv_number_signal2 = workspace2.var("rate_BulkWW_xww_for_unbin");
        else: 
            rrv_number_signal1 = workspace1.var("rate_%s_xww_for_unbin"%(self.signal_sample));
            rrv_number_signal2 = workspace2.var("rate_%s_xww_for_unbin"%(self.signal_sample));
            
         
        rrv_number_WJets1  = workspace1.var("rate_WJets_xww_for_unbin");
        rrv_number_VV1     = workspace1.var("rate_VV_xww_for_unbin");
        rrv_number_TTbar1  = workspace1.var("rate_TTbar_xww_for_unbin");
        rrv_number_STop1   = workspace1.var("rate_STop_xww_for_unbin");

        rrv_number_WJets2  = workspace2.var("rate_WJets_xww_for_unbin");
        rrv_number_VV2     = workspace2.var("rate_VV_xww_for_unbin");
        rrv_number_TTbar2  = workspace2.var("rate_TTbar_xww_for_unbin");
        rrv_number_STop2   = workspace2.var("rate_STop_xww_for_unbin");

        rrv_number_signal1.Print();
        rrv_number_WJets1.Print();
        rrv_number_VV1.Print();
        rrv_number_TTbar1.Print();
        rrv_number_STop1.Print();

        rrv_number_signal2.Print();
        rrv_number_WJets2.Print();
        rrv_number_VV2.Print();
        rrv_number_TTbar2.Print();
        rrv_number_STop2.Print();


        #### Prepare the final plot starting from total background 
        rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC_xww","rrv_number_Total_background_MC_xww",
                rrv_number_WJets1.getVal()+
                rrv_number_VV1.getVal()+
                rrv_number_TTbar1.getVal()+
                rrv_number_STop1.getVal()+
                rrv_number_WJets2.getVal()+
                rrv_number_VV2.getVal()+
                rrv_number_TTbar2.getVal()+
                rrv_number_STop2.getVal());

        rrv_number_Total_background_MC.setError(TMath.Sqrt(
                rrv_number_WJets1.getError()* rrv_number_WJets1.getError()+
                rrv_number_VV1.getError()* rrv_number_VV1.getError()+
                rrv_number_TTbar1.getError()* rrv_number_TTbar1.getError()+
                rrv_number_STop1.getError() *rrv_number_STop1.getError()+
                rrv_number_WJets2.getError()* rrv_number_WJets2.getError()+
                rrv_number_VV2.getError()* rrv_number_VV2.getError()+
                rrv_number_TTbar2.getError()* rrv_number_TTbar2.getError()+
                rrv_number_STop2.getError() *rrv_number_STop2.getError()
                ));
        data_obs.Print()
        rrv_number_Total_background_MC.Print()

        #### Total pdf 
        model_Total_background_MC = RooAddPdf("model_Total_background_MC_xww","model_Total_background_MC_xww",RooArgList(model_pdf_WJets1, model_pdf_WJets2,model_pdf_VV1,model_pdf_VV2,model_pdf_TTbar1,model_pdf_TTbar2,model_pdf_STop1,model_pdf_STop2),RooArgList(rrv_number_WJets1,rrv_number_WJets2,rrv_number_VV1,rrv_number_VV2,rrv_number_TTbar1,rrv_number_TTbar2,rrv_number_STop1,rrv_number_STop2));

        scale_number_signal = (rrv_number_signal1.getVal() + rrv_number_signal2.getVal() )/data_obs.sumEntries()
        #### scale factor in order to scale MC to data in the final plot -> in order to avoid the normalization to data which is done by default in rooFit
        scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()/data_obs.sumEntries()
                         
        #### create the frame
        mplot = rrv_x.frame(RooFit.Title("check_workspace"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));
        data_obs.plotOn(mplot , RooFit.Name("data_invisible"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0));

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets"), RooFit.Components("WJets_xww_*,VV_xww_*,TTbar_xww_*,STop_xww_*"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines());
        
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("VV"), RooFit.Components("VV_xww_*,TTbar_xww_*,STop_xww_*"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines());

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("TTbar"), RooFit.Components("TTbar_xww_*,STop_xww_*"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines());

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("STop"), RooFit.Components("STop_xww_*"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines());

        #solid line
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets_line_invisible"), RooFit.Components("WJets_xww_*,VV_xww_*,TTbar_xww_*,STop_xww_*"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("VV_line_invisible"), RooFit.Components("VV_xww_*,TTbar_xww_*,STop_xww_*"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("TTbar_line_invisible"), RooFit.Components("TTbar_xww_*,STop_xww_*"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("STop_line_invisible"), RooFit.Components("STop_xww_*"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());


        ### signal scale to be visible in the plots
        #label_tstring = TString(self.signal_sample);
        #if label_tstring.Contains("600") and (not label_tstring.Contains("1600")):
        #    signal_scale=20*self.xs_rescale;
        #elif label_tstring.Contains("700") and (not label_tstring.Contains("1700")):
        #    signal_scale=20*self.xs_rescale;
        #elif label_tstring.Contains("750") and (not label_tstring.Contains("1750")):
        #    signal_scale=20*self.xs_rescale;
        #elif label_tstring.Contains("800") and (not label_tstring.Contains("1800")):
        #    signal_scale=20*self.xs_rescale;
        #else:
        #    signal_scale=self.xs_rescale;
        signal_scale=10;


        model_pdf_signal = RooAddPdf("model_pdf_signal","model_pdf_signal",RooArgList(model_pdf_signal1, model_pdf_signal2),RooArgList(rrv_number_signal1,rrv_number_signal2));
        model_pdf_signal.plotOn(mplot,RooFit.Normalization(scale_number_signal*signal_scale),RooFit.Name("%s #times %s"%(self.signal_sample, signal_scale)),RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["Signal"]), RooFit.LineStyle(2), RooFit.VLines());

        #### plot the observed data using poissonian error bar
        self.getData_PoissonInterval(data_obs,mplot);
        
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Invisible());

        mplot_pull=self.get_pull(rrv_x,mplot);

        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        #draw_error_band(model_Total_background_MC, rrv_x.GetName(), rrv_number_Total_background_MC,self.FloatingParams,workspace1 ,mplot,self.color_palet["Uncertainty"],"F");

        mplot.Print();
        self.leg = self.legend4Plot(mplot,0,1,-0.01,-0.05,0.11,0.);
        self.leg.SetTextSize(0.036);
        mplot.addObject(self.leg);
        
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.8);
            

        parameters_list = RooArgList();
        self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s_%s_%s_g1/m_lvj_fitting/"%(options.additioninformation, "em",self.PS_model,self.wtagger_label),"check_workspace_for_limit","",0,1);
        
        #if workspace1.var("rrv_num_floatparameter_in_last_fitting"):
        #    self.nPar_float_in_fitTo = int(workspace1.var("rrv_num_floatparameter_in_last_fitting").getVal());
        #else:
        #    self.nPar_float_in_fitTo = self.FloatingParams.getSize();
        #nBinX = mplot.GetNbinsX();
        #ndof  = nBinX-self.nPar_float_in_fitTo;
        #print "nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare( self.nPar_float_in_fitTo )*ndof, ndof );

    ### in order to get the pull

    def read_workspace(self, logy=0):

        ### Taket the workspace for limits  
        file = TFile(self.file_rlt_root) ;
        workspace = file.Get("workspace4limit_") ;
        workspace.Print()

        ### iterate on the workspace element parameters
        print "----------- Parameter Workspace -------------";
        parameters_workspace = workspace.allVars();
        par = parameters_workspace.createIterator();
        par.Reset();
        param = par.Next()
        while (param):
            param.Print();
            param=par.Next()
        print "---------------------------------------------";

        workspace.data("data_obs_xww_%s_%s"%(self.channel,self.wtagger_label)).Print()

        print "----------- Pdf in the Workspace -------------";
        pdfs_workspace = workspace.allPdfs();
        par = pdfs_workspace.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            param.Print();
            param = par.Next()
        print "----------------------------------------------";

        rrv_x = workspace.var("rrv_mass_lvj")
        data_obs = workspace.data("data_obs_xww_%s_%s"%(self.channel,self.wtagger_label));
        if TString(self.signal_sample).Contains("Bulk_WW"):
            model_pdf_signal = workspace.pdf("BulkWW_xww_%s_%s"%(self.channel,self.wtagger_label));
        else: 
            model_pdf_signal = workspace.pdf("%s_xww_%s_%s"%(self.signal_sample,self.channel,self.wtagger_label));
            
        model_pdf_WJets  = workspace.pdf("WJets_xww_%s_%s"%(self.channel,self.wtagger_label));
        model_pdf_VV     = workspace.pdf("VV_xww_%s_%s"%(self.channel,self.wtagger_label));
        model_pdf_TTbar  = workspace.pdf("TTbar_xww_%s_%s"%(self.channel,self.wtagger_label));
        model_pdf_STop   = workspace.pdf("STop_xww_%s_%s"%(self.channel,self.wtagger_label));

        model_pdf_signal.Print();
        model_pdf_WJets.Print();
        model_pdf_VV.Print();
        model_pdf_TTbar.Print();
        model_pdf_STop.Print();

        if TString(self.signal_sample).Contains("Bulk_WW"): 
            rrv_number_signal = workspace.var("rate_BulkWW_xww_for_unbin");
        else: 
            rrv_number_signal = workspace.var("rate_%s_xww_for_unbin"%(self.signal_sample));
            
         
        rrv_number_WJets  = workspace.var("rate_WJets_xww_for_unbin");
        rrv_number_VV     = workspace.var("rate_VV_xww_for_unbin");
        rrv_number_TTbar  = workspace.var("rate_TTbar_xww_for_unbin");
        rrv_number_STop   = workspace.var("rate_STop_xww_for_unbin");

        rrv_number_signal.Print();
        rrv_number_WJets.Print();
        rrv_number_VV.Print();
        rrv_number_TTbar.Print();
        rrv_number_STop.Print();

        #### Prepare the final plot starting from total background 
        rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC_xww","rrv_number_Total_background_MC_xww",
                rrv_number_WJets.getVal()+
                rrv_number_VV.getVal()+
                rrv_number_TTbar.getVal()+
                rrv_number_STop.getVal());

        rrv_number_Total_background_MC.setError(TMath.Sqrt(
                rrv_number_WJets.getError()* rrv_number_WJets.getError()+
                rrv_number_VV.getError()* rrv_number_VV.getError()+
                rrv_number_TTbar.getError()* rrv_number_TTbar.getError()+
                rrv_number_STop.getError() *rrv_number_STop.getError()
                ));

        #### Total pdf 
        model_Total_background_MC = RooAddPdf("model_Total_background_MC_xww","model_Total_background_MC_xww",RooArgList(model_pdf_WJets,model_pdf_VV,model_pdf_TTbar,model_pdf_STop),RooArgList(rrv_number_WJets,rrv_number_VV,rrv_number_TTbar,rrv_number_STop));

        scale_number_signal = rrv_number_signal.getVal()/data_obs.sumEntries()
        #### scale factor in order to scale MC to data in the final plot -> in order to avoid the normalization to data which is done by default in rooFit
        scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()/data_obs.sumEntries()
                         
        #### create the frame
        mplot = rrv_x.frame(RooFit.Title("check_workspace"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));
        data_obs.plotOn(mplot , RooFit.Name("data_invisible"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0));

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets"), RooFit.Components("WJets_xww_%s_%s,VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines());
        
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("VV"), RooFit.Components("VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines());

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("TTbar"), RooFit.Components("TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines());

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("STop"), RooFit.Components("STop_xww_%s_%s"%(self.channel,self.wtagger_label)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines());

        #solid line
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets_line_invisible"), RooFit.Components("WJets_xww_%s_%s,VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("VV_line_invisible"), RooFit.Components("VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("TTbar_line_invisible"), RooFit.Components("TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("STop_line_invisible"), RooFit.Components("STop_xww_%s_%s"%(self.channel,self.wtagger_label)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

        ### signal scale to be visible in the plots
        #label_tstring = TString(self.signal_sample);
        #if label_tstring.Contains("600") and (not label_tstring.Contains("1600")):
        #    signal_scale=20*self.xs_rescale;
        #elif label_tstring.Contains("700") and (not label_tstring.Contains("1700")):
        #    signal_scale=20*self.xs_rescale;
        #elif label_tstring.Contains("750") and (not label_tstring.Contains("1750")):
        #    signal_scale=20*self.xs_rescale;
        #elif label_tstring.Contains("800") and (not label_tstring.Contains("1800")):
        #    signal_scale=20*self.xs_rescale;
        #else:
        #    signal_scale=self.xs_rescale;
        signal_scale=10;


        model_pdf_signal.plotOn(mplot,RooFit.Normalization(scale_number_signal*signal_scale),RooFit.Name("%s #times %s"%(self.signal_sample, signal_scale)),RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["Signal"]), RooFit.LineStyle(2), RooFit.VLines());

        #### plot the observed data using poissonian error bar
        self.getData_PoissonInterval(data_obs,mplot);
        
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Invisible());

        mplot_pull=self.get_pull(rrv_x,mplot);

        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        draw_error_band(model_Total_background_MC, rrv_x.GetName(), rrv_number_Total_background_MC,self.FloatingParams,workspace ,mplot,self.color_palet["Uncertainty"],"F");

        mplot.Print();
        self.leg = self.legend4Plot(mplot,0,1,-0.01,-0.05,0.11,0.);
        self.leg.SetTextSize(0.036);
        mplot.addObject(self.leg);
        
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.8);
            

        parameters_list = RooArgList();
        self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s_%s_%s_g1/m_lvj_fitting/"%(options.additioninformation, self.channel,self.PS_model,self.wtagger_label),"check_workspace_for_limit","",0,1);
        
        if workspace.var("rrv_num_floatparameter_in_last_fitting"):
            self.nPar_float_in_fitTo = int(workspace.var("rrv_num_floatparameter_in_last_fitting").getVal());
        else:
            self.nPar_float_in_fitTo = self.FloatingParams.getSize();
        nBinX = mplot.GetNbinsX();
        ndof  = nBinX-self.nPar_float_in_fitTo;
        print "nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare( self.nPar_float_in_fitTo )*ndof, ndof );

    def get_pull(self, rrv_x, mplot_orig):

        #print "############### draw the pull plot ########################"
        hpull = mplot_orig.pullHist();
        x = ROOT.Double(0.); y = ROOT.Double(0) ;
        for ipoint in range(0,hpull.GetN()):
          hpull.GetPoint(ipoint,x,y);
          if(y == 0):
           hpull.SetPoint(ipoint,x,10)
       
        mplot_pull = rrv_x.frame(RooFit.Title("Pull Distribution"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));
        medianLine = TLine(rrv_x.getMin(),0.,rrv_x.getMax(),0); medianLine.SetLineWidth(2); medianLine.SetLineColor(kRed);
        mplot_pull.addObject(medianLine);
        mplot_pull.addPlotable(hpull,"P");
        mplot_pull.SetTitle("");
        mplot_pull.GetXaxis().SetTitle("");
        mplot_pull.GetYaxis().SetRangeUser(-5,5);
        mplot_pull.GetYaxis().SetTitleSize(0.10);
        mplot_pull.GetYaxis().SetLabelSize(0.10);
        mplot_pull.GetXaxis().SetTitleSize(0.10);
        mplot_pull.GetXaxis().SetLabelSize(0.10);
        mplot_pull.GetYaxis().SetTitleOffset(0.40);
        mplot_pull.GetYaxis().SetTitle("#frac{Data-Fit}{#sigma_{data}}");
        mplot_pull.GetYaxis().CenterTitle();

        return mplot_pull;

    def getData_PoissonInterval(self,data_obs,mplot):
        rrv_x = self.workspace4fit_.var("rrv_mass_lvj");
        datahist   = data_obs.binnedClone(data_obs.GetName()+"_binnedClone",data_obs.GetName()+"_binnedClone");
        data_histo = datahist.createHistogram("histo_data",rrv_x) ;
        data_histo.SetName("data");
        data_plot  = RooHist(data_histo);
        data_plot.SetMarkerStyle(20);
        data_plot.SetMarkerSize(1.5);
        
        alpha = 1 - 0.6827;
        for iPoint  in range(data_plot.GetN()):
          N = data_plot.GetY()[iPoint];
          if N==0 : L = 0;
          else : L = (ROOT.Math.gamma_quantile(alpha/2,N,1.));
          U =  ROOT.Math.gamma_quantile_c(alpha/2,N+1,1);
          data_plot.SetPointEYlow(iPoint, N-L);
          data_plot.SetPointEYhigh(iPoint,U-N);
          data_plot.SetPointEXlow(iPoint,0);        
          data_plot.SetPointEXhigh(iPoint,0);        
        
        mplot.addPlotable(data_plot,"PE");


    #### in order to make the banner on the plots
    def banner4Plot(self, iswithpull=0):
      #print "############### draw the banner ########################"

      if iswithpull:
       if self.channel=="el":
        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 13TeV, W#rightarrow e #nu "%(self.GetLumi())));
#        banner = TLatex(0.18,0.96,"CMS                                             L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
       elif self.channel=="mu":
        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 13TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
#        banner = TLatex(0.18,0.96,"CMS                                             L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
       elif self.channel=="em":
        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 13TeV, W#rightarrow #mu,e #nu "%(self.GetLumi())));
#        banner = TLatex(0.18,0.96,"CMS                                             L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
       banner.SetNDC(); banner.SetTextSize(0.041);
      else:
       if self.channel=="el":
        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 13TeV, W#rightarrow e #nu "%(self.GetLumi())));
#        banner = TLatex(0.18,0.96,"CMS                                         L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
       if self.channel=="mu":
        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 13TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
        #banner = TLatex(0.18,0.96,"CMS                                         L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
       if self.channel=="em":
        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 13TeV, W#rightarrow #mu,e #nu "%(self.GetLumi())));
#        banner = TLatex(0.18,0.96,"CM                                          L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
       banner.SetNDC(); banner.SetTextSize(0.033);
                                                                                                         
      return banner;

    ### in order to make the legend
    def legend4Plot(self, plot, left=1, isFill=1, x_offset_low=0., y_offset_low=0., x_offset_high =0., y_offset_high =0., TwoCoulum =1.):
        #print "############### draw the legend ########################"
        if left==-1:
            theLeg = TLegend(0.65+x_offset_low, 0.58+y_offset_low, 0.93+x_offset_low, 0.87+y_offset_low, "", "NDC");
            theLeg.SetName("theLegend");
            theLeg.SetLineColor(0);
            theLeg.SetTextFont(42);
            theLeg.SetTextSize(.04);
        else:
            theLeg = TLegend(0.41+x_offset_low, 0.61+y_offset_low, 0.76+x_offset_high, 0.93+y_offset_high, "", "NDC");            
            theLeg.SetName("theLegend");
            if TwoCoulum :
                theLeg.SetNColumns(2);

        theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0);
        theLeg.SetBorderSize(0);
        theLeg.SetLineColor(0);
        theLeg.SetLineWidth(0);
        theLeg.SetLineStyle(0);
        theLeg.SetTextSize(0.040);
        theLeg.SetTextFont(42);

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        

        ##if   self.categoryID == 0: legHeader="(e#nu)"; #, 1JALLP)";
        ##elif self.categoryID == 1: legHeader="(#mu#nu)"; #, 1JALLP)";
        ##elif self.categoryID ==-2: legHeader="(e#nu, 1JLP)";
        ##elif self.categoryID ==-1: legHeader="(e#nu, 1JHP)";
        ##elif self.categoryID == 2: legHeader="(#mu#nu, 1JLP)";
        ##elif self.categoryID == 1: legHeader="(#mu#nu, 1JHP)";
        if   self.categoryID <0: legHeader="(e#nu)"
        elif self.categoryID >0: legHeader="(#mu#nu)";

        for obj in range(int(plot.numItems()) ):
          objName = plot.nameOf(obj);
          if objName == "errorband" : objName = "Uncertainty";
          print objName;
          if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
objName ==objName_before ):
            theObj = plot.getObject(obj);
            objTitle = objName;
            drawoption= plot.getDrawOptions(objName).Data()
            if drawoption=="P":drawoption="PE"
            if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):  objName_before=objName; continue ;
            elif TString(objName).Contains("Graph") :  objName_before=objName; continue ;
            elif TString(objName).Data()=="data" : theLeg.AddEntry(theObj, "CMS Data "+legHeader,"PE");  objName_before=objName;                 
            else: objName_before=objName; continue ;

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        
                   
        for obj in range(int(plot.numItems()) ):
          objName = plot.nameOf(obj);
          if objName == "errorband" : objName = "Uncertainty";
          print objName;
          if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
objName ==objName_before ):
            theObj = plot.getObject(obj);
            objTitle = objName;
            drawoption= plot.getDrawOptions(objName).Data()
            if drawoption=="P":drawoption="PE"
            if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):  objName_before=objName; continue ;
            elif TString(objName).Contains("Graph") :  objName_before=objName; continue ;
            elif TString(objName).Data()=="WJets" : theLeg.AddEntry(theObj, "W+jets","F");  objName_before=objName;                 
            else:  objName_before=objName; continue ;

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        


        for obj in range(int(plot.numItems()) ):
            objName = plot.nameOf(obj);
            if objName == "errorband" : objName = "Uncertainty";
            print objName;
            if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
objName ==objName_before ):
                theObj = plot.getObject(obj);
                objTitle = objName;
                drawoption= plot.getDrawOptions(objName).Data()
                if drawoption=="P":drawoption="PE"
                if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):
                    theLeg.AddEntry(theObj, objName,"F");
                elif TString(objName).Contains("Graph") :
                    if not (objName_before=="Graph" or objName_before=="Uncertainty"): theLeg.AddEntry(theObj, "Uncertainty","F");
                else:
                    if TString(objName).Data()=="STop" : theLeg.AddEntry(theObj, "Single Top","F");
                    elif TString(objName).Data()=="TTbar" : theLeg.AddEntry(theObj, "t#bar{t}","F");
                    elif TString(objName).Data()=="VV" : theLeg.AddEntry(theObj, "WW/WZ","F");
                    elif TString(objName).Data()=="data" :  objName_before=objName; entryCnt = entryCnt+1; continue ;
                    elif TString(objName).Data()=="WJets" : objName_before=objName; entryCnt = entryCnt+1; continue;
                    elif TString(objName).Contains("vbfH"): theLeg.AddEntry(theObj, (TString(objName).ReplaceAll("vbfH","qqH")).Data() ,"L");
                    elif TString(objName).Contains("Uncertainty"): theLeg.AddEntry(theObj, objTitle,drawoption);
                    elif TString(objName).Contains("Bulk"):
                        objName_signal_graviton = theObj
                        objNameLeg_signal_graviton = objName
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M600") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M600"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.6 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M700") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M700"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.7 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M800") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M800"):
                       #    objName_signal_graviton = theObj ; 
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.8 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M900") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M900"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.9 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1000") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1000"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1100") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1100"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.1 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1200") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1200"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.2 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1300") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1300"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.3 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1400") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1400"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.4 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1500") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1500"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.5 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1600") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1600"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.6 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1700") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1700"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.7 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1800") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1800"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.8 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1900") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1900"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.9 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2000") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2000"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=2 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2100") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2100"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.1 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2200") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2200"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.2 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2300") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2300"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.3 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2400") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2400"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.4 TeV #tilde{k}=0.5 (#times100)";
                       #if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2500") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2500"):
                       #    objName_signal_graviton = theObj ;
                       #    objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.5 TeV #tilde{k}=0.5 (#times100)";
                    else : theLeg.AddEntry(theObj, objTitle,drawoption);
                entryCnt=entryCnt+1;
            objName_before=objName;
        if objName_signal_graviton !="" :
           theLeg.AddEntry(objName_signal_graviton, TString(objNameLeg_signal_graviton).Data() ,"L");
        return theLeg;

    #### draw canvas with plots with pull
    def draw_canvas_with_pull(self, mplot, mplot_pull,parameters_list,in_directory, in_file_name, in_model_name="", show_constant_parameter=0, logy=0):# mplot + pull
        #print "############### draw the canvas with pull ########################"
        mplot.GetXaxis().SetTitleOffset(1.1);
        mplot.GetYaxis().SetTitleOffset(1.3);
        mplot.GetXaxis().SetTitleSize(0.055);
        mplot.GetYaxis().SetTitleSize(0.055);
        mplot.GetXaxis().SetLabelSize(0.045);
        mplot.GetYaxis().SetLabelSize(0.045);
        mplot_pull.GetXaxis().SetLabelSize(0.14);
        mplot_pull.GetYaxis().SetLabelSize(0.14);
        mplot_pull.GetYaxis().SetTitleSize(0.15);
        mplot_pull.GetYaxis().SetNdivisions(205);

                                                                          
        cMassFit = TCanvas("cMassFit","cMassFit", 600,600);
        # if parameters_list is empty, don't draw pad3
        par_first=parameters_list.createIterator();
        par_first.Reset();
        param_first=par_first.Next()
        doParameterPlot = 0 ;
        if param_first and doParameterPlot != 0:
         pad1=TPad("pad1","pad1",0.,0. ,0.8,0.24);
         pad2=TPad("pad2","pad2",0.,0.24,0.8,1. );
         pad3=TPad("pad3","pad3",0.8,0.,1,1);
         pad1.Draw();
         pad2.Draw();
         pad3.Draw();
        else:
         pad1=TPad("pad1","pad1",0.,0. ,0.99,0.24);
         pad2=TPad("pad2","pad2",0.,0.24,0.99,1. );
         pad1.Draw();
         pad2.Draw();
                                                                                                                                                                              
        pad2.cd();
        mplot.Draw();
        banner = self.banner4Plot(1);
        banner.Draw();

        pad1.cd();
        mplot_pull.Draw();

        if param_first and doParameterPlot != 0:

            pad3.cd();
            latex=TLatex();
            latex.SetTextSize(0.1);
            par=parameters_list.createIterator();
            par.Reset();
            param=par.Next()
            i=0;
            while param:
                if (not param.isConstant() ) or show_constant_parameter:
                    param.Print();
                    icolor=1;#if a paramenter is constant, color is 2
                    if param.isConstant(): icolor=2
                    latex.DrawLatex(0,0.9-i*0.04,"#color[%s]{%s}"%(icolor,param.GetName()) );
                    latex.DrawLatex(0,0.9-i*0.04-0.02," #color[%s]{%4.3e +/- %2.1e}"%(icolor,param.getVal(),param.getError()) );
                    i=i+1;
                param=par.Next();

        ## create the directory where store the plots
        Directory = TString(in_directory+self.signal_sample+"_%02d_%02d/"%(options.cprime,options.BRnew));
        if not Directory.EndsWith("/"):Directory = Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
             os.system("mkdir -p "+Directory.Data());

        rlt_file = TString(Directory.Data()+in_file_name);
        if rlt_file.EndsWith(".root"):
            TString(in_model_name).ReplaceAll(".root","");
            rlt_file.ReplaceAll(".root","_"+in_model_name+"_with_pull.png");
        else:
            TString(in_model_name).ReplaceAll(".root","");
            rlt_file.ReplaceAll(".root","");
            rlt_file=rlt_file.Append("_"+in_model_name+"_with_pull.png");

        cMassFit.SaveAs(rlt_file.Data());
        if logy:
            mplot.GetYaxis().SetRangeUser(1e-3,mplot.GetMaximum()*200);
            pad2.SetLogy() ;
            pad2.Update();
            cMassFit.Update();
            rlt_file.ReplaceAll(".png","_log.png");
            cMassFit.SaveAs(rlt_file.Data());


        #rlt_file.ReplaceAll(".png",".pdf");
        #cMassFit.SaveAs(rlt_file.Data());
        #
        #rlt_file.ReplaceAll(".pdf",".root");
        #cMassFit.SaveAs(rlt_file.Data());

        #string_file_name = TString(in_file_name);
        #if string_file_name.EndsWith(".root"):
        #    string_file_name.ReplaceAll(".root","_"+in_model_name);
        #else:
        #    string_file_name.ReplaceAll(".root","");
        #    string_file_name.Append("_"+in_model_name);

        #if logy:
        #    mplot.GetYaxis().SetRangeUser(1e-3,mplot.GetMaximum()*200);
        #    pad2.SetLogy() ;
        #    pad2.Update();
        #    cMassFit.Update();
        #    rlt_file.ReplaceAll(".root","_log.root");
        #    cMassFit.SaveAs(rlt_file.Data());
        #    rlt_file.ReplaceAll(".root",".pdf");
        #    cMassFit.SaveAs(rlt_file.Data());
        #    rlt_file.ReplaceAll(".pdf",".png");
        #    cMassFit.SaveAs(rlt_file.Data());

        ##self.draw_canvas(mplot,in_directory,string_file_name.Data(),0,logy,1);

    #### jusr drawing canvas with no pull
    def draw_canvas(self, in_obj,in_directory, in_file_name, is_range=0, logy=0, frompull=0):
        #print "############### draw the canvas without pull ########################"
        cMassFit = TCanvas("cMassFit","cMassFit", 600,600);

        if frompull and logy :
            in_obj.GetYaxis().SetRangeUser(1e-2,in_obj.GetMaximum()/200)
        elif not frompull and logy :
            in_obj.GetYaxis().SetRangeUser(0.00001,in_obj.GetMaximum())
            

        if is_range:
            h2=TH2D("h2","",100,400,1400,4,0.00001,4);
            h2.Draw();
            in_obj.Draw("same")
        else :
            in_obj.Draw()

        in_obj.GetXaxis().SetTitleSize(0.045);
        in_obj.GetXaxis().SetTitleOffset(1.15);
        in_obj.GetXaxis().SetLabelSize(0.04);

        in_obj.GetYaxis().SetTitleSize(0.055);
        in_obj.GetYaxis().SetTitleOffset(1.40);
        in_obj.GetYaxis().SetLabelSize(0.04);

        self.leg.SetTextSize(0.031); 

        banner = self.banner4Plot();
        banner.Draw();
        
        Directory=TString(in_directory+self.signal_sample+"_%02d_%02d/"%(options.cprime,options.BRnew));
        if not Directory.EndsWith("/"):Directory=Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
              os.system("mkdir -p "+Directory.Data());

        rlt_file=TString(Directory.Data()+in_file_name);
        if rlt_file.EndsWith(".root"):
            rlt_file.ReplaceAll(".root","_rlt_without_pull_and_paramters.png");
        else:
            rlt_file.ReplaceAll(".root","");
            rlt_file = rlt_file.Append(".png");

        cMassFit.SaveAs(rlt_file.Data());
        if logy:
            in_obj.GetYaxis().SetRangeUser(1e-3,in_obj.GetMaximum()*200);
            cMassFit.SetLogy() ;
            cMassFit.Update();
            rlt_file.ReplaceAll(".png","_log.png");
            cMassFit.SaveAs(rlt_file.Data());

        #rlt_file.ReplaceAll(".png",".pdf");
        #cMassFit.SaveAs(rlt_file.Data());

        #rlt_file.ReplaceAll(".pdf",".root");
        #cMassFit.SaveAs(rlt_file.Data());

        #if logy:
        #    in_obj.GetYaxis().SetRangeUser(1e-3,in_obj.GetMaximum()*200);
        #    cMassFit.SetLogy() ;
        #    cMassFit.Update();
        #    rlt_file.ReplaceAll(".root","_log.root");
        #    cMassFit.SaveAs(rlt_file.Data());
        #    rlt_file.ReplaceAll(".root",".pdf");
        #    cMassFit.SaveAs(rlt_file.Data());
        #    rlt_file.ReplaceAll(".pdf",".png");
        #    cMassFit.SaveAs(rlt_file.Data());
       

    ##### Get Lumi for banner title
    def GetLumi(self):
        if self.channel=="el":   return 2.18;
        elif self.channel=="mu": return 2.17;
        elif self.channel=="em": return 1.0;


### function to check the workspace once it has already created
def check_workspace(channel, higgs):
    boostedW_fitter = doFit_wj_and_wlvj("mu",higgs);
    boostedW_fitter.read_workspace()
    boostedW_fitter.read_2workspaces()

#### Main Code
if __name__ == '__main__':
    channel=options.channel;

    if options.check:
        print '################# check workspace for %s sample'%(channel);
        check_workspace(channel,"BulkGravWW750");

