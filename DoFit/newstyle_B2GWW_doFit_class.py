#! /usr/bin/env python
import os
from array import array
import ROOT
import ntpath
import sys
from optparse import OptionParser
import CMS_lumi, tdrstyle

############################################
#              Job steering                #
############################################

parser = OptionParser()

parser.add_option('-a', '--additioninformation', action = "store", type = "string", dest = "additioninformation", default = "B2GWW")
parser.add_option('-b', action = 'store_true', dest = 'noX', default = True, help = 'no X11 windows')
parser.add_option('-c', '--channel', action = "store", type = "string", dest = "opt_channel", default = "mu")
#parser.add_option('-c', '--channel', action = "store", type = "string", dest = "opt_channel", default = "el")
parser.add_option('-w', '--width', action = "store", type = "int", dest='width', default=1, help = "bin width")#1: 100; 0: 50

parser.add_option('-s', '--simple', action = 'store', dest = 'simple', default = True, help = 'pre-limit in simple mode')
parser.add_option('-m', '--multi', action = 'store', dest = 'multi', default = False, help = 'pre-limit in multi mode')
parser.add_option('--fitwtagger', action = 'store_true', dest = 'fitwtagger', default = False, help = 'fit wtagger jet in ttbar control sample')
parser.add_option('--fitwtaggersim', action = 'store_true', dest = 'fitwtaggersim', default = False, help = 'fit wtagger jet in ttbar control sample with mu and el samples simultaneously')
parser.add_option('--check', action = 'store_true', dest = 'check', default = False, help = 'check the workspace for limit setting')
parser.add_option('--combine', action = 'store_true', dest = 'combine', default = False, help = 'combine el and mu')
parser.add_option('--control', action = 'store_true', dest = 'control', default = False, help = 'control plot')
parser.add_option('--realdata', action='store', type = "int", dest='realdata', default=0, help='real data or pdata')
parser.add_option('--keepblind', action='store', type = "int", dest='keepblind', default=1, help='keep blind for real data')
parser.add_option('--fitsignal', action='store', type = "int", dest='fitsignal', default=0, help='fit only signal lineshape with a chosen model')
parser.add_option('--closuretest', action='store', type = "int", dest='closuretest', default=0, help='closure test; 0: no test; 1: A1->A2; 2: A->B')
parser.add_option('--MCStudy', action = 'store', dest = 'MCStudy', default = False, help = 'Do MCStudy for each fitTo')
#parser.add_option('--MCStudy', action = 'store', dest = 'MCStudy', default = True, help = 'Do MCStudy for each fitTo')

parser.add_option('--inPath', action = "store", type = "string", dest = "inPath", default = ".")
parser.add_option('--category', action = "store", type = "string", dest = "opt_wtagger_category", default = "HP")#Wtagger Category: HP, LP, ALLP
#parser.add_option('--category', action = "store", type = "string", dest = "category", default = "ALLP")

(options, args) = parser.parse_args()

ROOT.gROOT.ProcessLine(".L %s/PDFs/PdfDiagonalizer.cc+"%(options.inPath))
ROOT.gROOT.ProcessLine(".L %s/PDFs/HWWLVJRooPdfs.cxx+"%(options.inPath))
ROOT.gROOT.ProcessLine(".L %s/PDFs/Util.cxx+"%(options.inPath))
##ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
##ROOT.gSystem.Load(options.inPath+"/PDFs/Util_cxx.so")
##ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")

from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1, TFile, TLine, TLegend, TH1D, TH2D, THStack, TGraph, TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf, RooCBShape, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG, RooDataSet, RooExponential, RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar, RooFormulaVar, RooDataHist, RooHist, RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf, RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen, kOrange, kBlack, kBlue, kCyan, kMagenta, kWhite, draw_error_band, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error, RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf,  RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf, MCStudy

###############################
## doFit Class Implemetation ##
###############################

class DoFit:
    """Do WV->lvjet analysis"""

    def __init__(self, in_channel, in_signal_sample,
            in_mlvj_signal_region_min=500, in_mlvj_signal_region_max=700, 
            in_mj_min=40, in_mj_max=150, in_mlvj_min=600., in_mlvj_max=1500.,
            fit_model = "ExpN", fit_model_alter = "Pow",
            input_workspace=None):

        tdrstyle.setTDRStyle()
        TGaxis.SetMaxDigits(3)

        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9)
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9)

        ### set the channel type --> electron or muon
        self.channel = in_channel
        self.leg = TLegend()
        self.MODEL_4_mlvj = fit_model
        self.MODEL_4_mlvj_alter = fit_model_alter

        print "########################################################################################"
        print "######## define class: binning, variables, cuts, files and nuisance parameters ########"
        print "########################################################################################"

        ### Set the mj binning for plots
        self.BinWidth_mj = 5.

        ### Set the binning for mlvj plots as a function of the model
        if options.fitsignal == 1:
            self.BinWidth_mlvj = 20.
        else:
            if options.width == 1:
                self.BinWidth_mlvj = 100.
            else:
                self.BinWidth_mlvj = 50.


        #narrow the BinWidth_mj and BinWidth_mlvj by a factor of 5. Because Higgs-Combination-Tools will generate a binned sample, so need the bin width narrow. So, as a easy selution, we will increase the bin-width by a factor of 5 when ploting m_j m_WW
        self.binwidth_narrow_factor = 5.

        ## rounding the mj upper edge
        self.BinWidth_mj = self.BinWidth_mj/self.binwidth_narrow_factor
        nbins_mj = int((in_mj_max-in_mj_min)/self.BinWidth_mj)
        in_mj_max = in_mj_min+nbins_mj*self.BinWidth_mj

        ## rounding the mj upper edge
        self.BinWidth_mlvj = self.BinWidth_mlvj/self.binwidth_narrow_factor
        nbins_mlvj = int((in_mlvj_max-in_mlvj_min)/self.BinWidth_mlvj)
        in_mlvj_max = in_mlvj_min+nbins_mlvj*self.BinWidth_mlvj

        ## define jet mass variable
        rrv_mass_j = RooRealVar("rrv_mass_j", "m_{jet}", (in_mj_min+in_mj_max)/2., in_mj_min, in_mj_max, "GeV")
        rrv_mass_j.setBins(nbins_mj)

        ## define invariant mass WW variable
        rrv_mass_lvj = RooRealVar("rrv_mass_lvj", "m_{WW}", (in_mlvj_min+in_mlvj_max)/2., in_mlvj_min, in_mlvj_max, "GeV")
        rrv_mass_lvj.setBins(nbins_mlvj)
        rrv_mass_lvj.Print("v")

        ## set the model used for the background parametrization
        self.MODEL_4_mlvj = fit_model
        self.MODEL_4_mlvj_alter = fit_model_alter

        ## create the workspace and import them
        if input_workspace is None:
            self.workspace4fit_ = RooWorkspace("workspace4fit_", "workspace4fit_")
        else:
            self.workspace4fit_ = input_workspace
        getattr(self.workspace4fit_, "import")(rrv_mass_j)
        getattr(self.workspace4fit_, "import")(rrv_mass_lvj)

        #prepare workspace for unbin-Limit -> just for the stuff on which running the limit
        self.workspace4limit_ = RooWorkspace("workspace4limit_", "workspace4limit_")

        self.signal_sample = in_signal_sample
        self.signal_model =  filter(str.isalpha, in_signal_sample)
        self.signal_mass =  int(filter(str.isdigit, in_signal_sample))  #GeV
        if self.signal_model== "WprimeWZHVTA": 
            self.signal_model= "WprimeWZ-HVT-A"

        #print self.signal_sample
        #print self.signal_model
        #print self.signal_mass

        ## different code operation mode -> just normal analysis
        if options.closuretest == 0:
            self.mj_sideband_lo_min = in_mj_min
            self.mj_sideband_lo_max = 65
            if self.signal_model == "BulkGravWW":
                self.mj_signal_min = 65
                self.mj_signal_max = 95
            elif self.signal_model == "WprimeWZ":
                self.mj_signal_min = 75
                self.mj_signal_max = 105
            elif self.signal_model == "WprimeWZ-HVT-A":
                self.mj_signal_min = 75
                self.mj_signal_max = 105
            else:
                print self.signal_model
                raw_input("failed to find the signal!")
            self.mj_sideband_hi_min = 135
            self.mj_sideband_hi_max = in_mj_max
        if options.closuretest == 1: ##closure test A1->A2
            self.mj_sideband_lo_min = in_mj_min
            self.mj_sideband_lo_max = 55
            self.mj_signal_min = 55
            self.mj_signal_max = 65
            self.mj_sideband_hi_min = 135
            self.mj_sideband_hi_max = in_mj_max
        if options.closuretest == 2: #closure test A->B
            self.mj_sideband_lo_min = in_mj_min
            self.mj_sideband_lo_max = 65
            self.mj_signal_min = 135
            self.mj_signal_max = 140
            self.mj_sideband_hi_min = 140
            self.mj_sideband_hi_max = in_mj_max

        ## zone definition in the jet mass
        rrv_mass_j.setRange("sb_lo", self.mj_sideband_lo_min, self.mj_sideband_lo_max)
        rrv_mass_j.setRange("signal_region", self.mj_signal_min, self.mj_signal_max)
        rrv_mass_j.setRange("sb_hi", self.mj_sideband_hi_min, self.mj_sideband_hi_max)
        rrv_mass_j.setRange("sblo_to_sbhi", self.mj_sideband_lo_min, self.mj_sideband_hi_max)

        ## signal region definition in the mlvj variable in case of counting limit
        self.mlvj_signal_min = in_mlvj_signal_region_min
        self.mlvj_signal_max = in_mlvj_signal_region_max
        rrv_mass_lvj.setRange("signal_region", self.mlvj_signal_min, self.mlvj_signal_max)
        rrv_mass_lvj.setRange("high_mass", 2500, in_mlvj_max)


        #if self.channel == "el":   self.Lumi= 2.6
        #elif self.channel == "mu": self.Lumi= 2.6
        #elif self.channel == "em": self.Lumi= 2.6
        ##prepare the data and mc files --> set the working directory and the files name
        #self.file_Directory = "PKUTree_final_2p6fb_Jun30/"


        if self.channel == "el":
            self.Lumi = 12.9 
        elif self.channel == "mu":
            self.Lumi = 12.9
        elif self.channel == "em":
            self.Lumi = 12.9
        #self.file_Directory = "PKUTree_final_6p26fb_Jul18/"
        self.file_Directory = "data_12p9_63mb_afterapproval/"


        if options.realdata == 1:
            self.file_data = (self.channel+"_PKUTree_16B.root")#keep blind!!!!
        else:
            self.file_data = (self.channel+"_PKUTree_pdata_SF.root")#keep blind!!!!
        #raw_input(self.file_data)

        self.file_signal     = (self.channel+"_PKUTree_%s.root"%(self.signal_sample))
        self.file_WJets0_mc  = (self.channel+"_PKUTree_WJetsPt180_xww.root")
        self.file_VV_mc      = (self.channel+"_PKUTree_VV_xww.root")# WW+WZ
        self.file_TTbar_mc   = (self.channel+"_PKUTree_TTBARpowheg_xww.root")
        self.file_STop_mc    = (self.channel+"_PKUTree_SingleTop_xww.root")

        ## event categorization as a function of the purity and the applied selection
        self.wtagger_category = options.opt_wtagger_category

        self.categoryID = 0
        if   self.wtagger_category == "ALLP" and self.channel == "el":
            self.categoryID = 0
        elif self.wtagger_category == "ALLP" and self.channel == "mu":
            self.categoryID = 1
        elif self.wtagger_category == "LP"   and self.channel == "el":
            self.categoryID = -2
        elif self.wtagger_category == "LP"   and self.channel == "mu":
            self.categoryID = 2
        elif self.wtagger_category == "HP"   and self.channel == "el":
            self.categoryID = -1
        elif self.wtagger_category == "HP"   and self.channel == "mu":
            self.categoryID = 1
        elif self.channel == "em" and options.combine:
            print "for combine"
        else:
            raw_input("Fail to find correct categoryID. Please check your wtaggerI:%s and channel:%s"%(self.wtagger_category, self.channel))

        if options.width == 1: #high mass signal >1.0TeV
            self.tau21_cut = 0.6
        else:#for 600GeV - 1TeV
            self.tau21_cut = 0.45
        print "self.tau21_cut = %s"%(self.tau21_cut)
        #medium wtagger_eff reweight between data and mc #Wtagger_forV SF have be add to ntuple weight
        if   self.wtagger_category == "LP":
            self.rrv_wtagger_eff_reweight_forT = RooRealVar("rrv_wtagger_eff_reweight_forT", "rrv_wtagger_eff_reweight_forT", 0.748)
            self.rrv_wtagger_eff_reweight_forT.setError(0.038)
            self.rrv_wtagger_eff_reweight_forV = RooRealVar("rrv_wtagger_eff_reweight_forV", "rrv_wtagger_eff_reweight_forV", 0.858)
            self.rrv_wtagger_eff_reweight_forV.setError(0.550)
        elif self.wtagger_category == "HP":

            if self.tau21_cut == 0.6:
                self.rrv_wtagger_eff_reweight_forT = RooRealVar("rrv_wtagger_eff_reweight_forT", "rrv_wtagger_eff_reweight_forT", 0.775)
                self.rrv_wtagger_eff_reweight_forT.setError(0.012)
                self.rrv_wtagger_eff_reweight_forV = RooRealVar("rrv_wtagger_eff_reweight_forV", "rrv_wtagger_eff_reweight_forV", 1.002)
                self.rrv_wtagger_eff_reweight_forV.setError(0.019)
            elif self.tau21_cut == 0.45:
                self.rrv_wtagger_eff_reweight_forT = RooRealVar("rrv_wtagger_eff_reweight_forT", "rrv_wtagger_eff_reweight_forT", 0.753)
                self.rrv_wtagger_eff_reweight_forT.setError(0.014)
                self.rrv_wtagger_eff_reweight_forV = RooRealVar("rrv_wtagger_eff_reweight_forV", "rrv_wtagger_eff_reweight_forV", 0.977)
                self.rrv_wtagger_eff_reweight_forV.setError(0.047)
            else:
                print "self.tau21_cut = %s"%(self.tau21_cut)
                raw_input("wrong tau21 cut value")
        else:
            raw_input("Fail to find correct categoryID. Please check your wtaggerI:%s and channel:%s"%(self.wtagger_category, self.channel))

        if not options.realdata == 1:
            self.rrv_wtagger_eff_reweight_forT.setVal(1.0)
            self.rrv_wtagger_eff_reweight_forV.setVal(1.0)
            self.rrv_wtagger_eff_reweight_forT.setError(self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal())
            self.rrv_wtagger_eff_reweight_forV.setError(self.rrv_wtagger_eff_reweight_forV.getError()/self.rrv_wtagger_eff_reweight_forV.getVal())

        print "VTagger efficiency correction for Top sample: %s +/- %s"%(self.rrv_wtagger_eff_reweight_forT.getVal(), self.rrv_wtagger_eff_reweight_forT.getError())
        print "VTagger efficiency correction for V sample: %s +/- %s"%(self.rrv_wtagger_eff_reweight_forV.getVal(), self.rrv_wtagger_eff_reweight_forV.getError())


        #correct the W-jet mass peak difference between data and MC
        if not options.realdata == 1:
            self.mean_shift = 0.0# 0.0 means no correction
            self.sigma_scale = 1.0# 1.0 means no correction
        else:
            if self.tau21_cut == 0.6:
                self.mean_shift = 84.4-83.4
                self.sigma_scale = 8.6/7.91
            elif self.tau21_cut == 0.45:
                self.mean_shift = 84.9-83.8
                self.sigma_scale = 7.92/7.52
        print " mean correction for the W peak : ", self.mean_shift, " Resolution correction : ", self.sigma_scale

        #SF of WJets and TTBar for ControlPlots, and these SF also passed to the final limit calculation
        if self.channel == "mu":
            self.controlplot_WJets_scale = 1.01
            self.controlplot_TTbar_scale = 0.78
        if self.channel == "el":
            self.controlplot_WJets_scale = 1.00
            self.controlplot_TTbar_scale = 0.84
        self.controlplot_WJets_scale = 1.0
        self.controlplot_TTbar_scale = 1.0

        #result files: The event number, parameters and error write into a txt file. The dataset and pdfs write into a root file
        self.datacardsDir =  "cards_%s_closuretest%s_%s_%s"%(options.additioninformation, options.closuretest, self.wtagger_category, fit_model)
        if not os.path.isdir(self.datacardsDir):
            os.system("mkdir -p %s"%(self.datacardsDir))

        #self.plotsDir = "plots_%s_%s_%s_closuretest%s"%(options.additioninformation, self.channel, self.wtagger_category, options.closuretest)
        self.plotsDir = "doFit_plots_%s_%s/closuretest%s_%s_%s"%(options.additioninformation, fit_model, options.closuretest, self.channel, self.wtagger_category)
        if not os.path.isdir(self.plotsDir):
            os.system("mkdir -p %s"%(self.plotsDir))


        ## workspace for limit
        self.file_rlt_root = self.datacardsDir+"/wwlvj_%s_%s_%s_workspace.root"%(self.signal_sample, self.channel, self.wtagger_category)
        ## datacard for the ubninned limit
        self.file_datacard_unbin = self.datacardsDir+"/wwlvj_%s_%s_%s_unbin.txt"%(self.signal_sample, self.channel, self.wtagger_category)
        ## workspace for the binned limit
        self.file_datacard_counting = self.datacardsDir+"/wwlvj_%s_%s_%s_counting.txt"%(self.signal_sample, self.channel, self.wtagger_category)
        ## appendix text file
        self.file_rlt_txt = self.datacardsDir+"/other_wwlvj_%s_%s_%s.txt"%(self.signal_sample, self.channel, self.wtagger_category)

        self.file_out = open(self.file_rlt_txt, "w")
        self.file_out.write("Welcome:\n")
        self.file_out.close()
        self.file_out = open(self.file_rlt_txt, "a+")

        ## color palette
        self.color_palet = { #color palet
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
        ### signal scaled up to save to datacard, the signal will be scalue up by additional 20 in the plots
        if self.signal_model == "BulkGravWW":
            table_signalscale =  {
                    600: 1,
                    700: 1,
                    750: 1,
                    800: 1,
                    900: 1,
                    1000: 1,
                    1200: 1,
                    1400: 1,
                    1600: 3,
                    1800: 5,
                    2000: 10,
                    2500: 50,
                    3000: 300,
                    3500: 2000,
                    4000: 5000,
                    4500: 20000 }
        elif self.signal_model == "WprimeWZ":
            table_signalscale =  {
                    800: 1,
                    1000: 1,
                    1200: 1,
                    1400: 1,
                    1600: 1,
                    1800: 1,
                    2000: 1,
                    2500: 1,
                    3000: 3,
                    3500: 20,
                    4000: 50,
                    4500: 200 }
        elif self.signal_model == "WprimeWZ-HVT-A":
            table_signalscale =  {
                    600: 1,
                    800: 1,
                    1000: 1,
                    1200: 1,
                    1400: 1,
                    1600: 1,
                    1800: 1,
                    2000: 1,
                    2500: 1,
                    3000: 3,
                    3500: 20,
                    4000: 50,
                    4500: 200 }
        self.signal_scale     = table_signalscale[self.signal_mass] ## scale factor for datacard and plot
        if self.signal_model == "BulkGravWW":
            self.signal_scale_plot = 20 ## one more scale factor for plot
        elif self.signal_model == "WprimeWZ":
            self.signal_scale_plot = 4 ## one more scale factor for plot
        elif self.signal_model == "WprimeWZ-HVT-A":
            self.signal_scale_plot = 4 ## one more scale factor for plot

        # parameters of data-driven method to get the WJets background event number.
        self.number_WJets_insideband = -1
        self.datadriven_alpha_WJets_unbin = -1
        self.datadriven_alpha_WJets_counting = -1

        ### uncertainty for datacard
        self.lumi_uncertainty    = 0.062

        if self.signal_model == "BulkGravWW":
            table_signal_pdfuncertainty =  {
                    600 : 0.09,
                    700 : 0.10,
                    750 : 0.10,
                    800 : 0.11,
                    900 : 0.13,
                    1000: 0.13,
                    1200: 0.14,
                    1400: 0.17,
                    1600: 0.19,
                    1800: 0.22,
                    2000: 0.25,
                    2500: 0.34,
                    3000: 0.45,
                    3500: 0.59,
                    4000: 0.78,
                    4500: 1.0   }
            table_signal_scaleuncertainty =  {
                    600 : 0.08,
                    700 : 0.09,
                    750 : 0.10,
                    800 : 0.10,
                    900 : 0.11,
                    1000: 0.11,
                    1200: 0.12,
                    1400: 0.13,
                    1600: 0.14,
                    1800: 0.15,
                    2000: 0.16,
                    2500: 0.17,
                    3000: 0.19,
                    3500: 0.20,
                    4000: 0.22,
                    4500: 0.23  }
        elif self.signal_model == "WprimeWZ" or self.signal_model == "WprimeWZ-HVT-A":
            table_signal_pdfuncertainty =  {
                    600 : 0.04,
                    700 : 0.04,
                    750 : 0.04,
                    800 : 0.04,
                    900 : 0.04,
                    1000: 0.04,
                    1200: 0.04,
                    1400: 0.05,
                    1600: 0.05,
                    1800: 0.05,
                    2000: 0.05,
                    2500: 0.07,
                    3000: 0.08,
                    3500: 0.12,
                    4000: 0.18,
                    4500: 0.31  }
            table_signal_scaleuncertainty =  {
                    600 : 0.01,
                    700 : 0.02,
                    750 : 0.02,
                    800 : 0.03,
                    900 : 0.03,
                    1000: 0.04,
                    1200: 0.05,
                    1400: 0.06,
                    1600: 0.07,
                    1800: 0.07,
                    2000: 0.08,
                    2500: 0.10,
                    3000: 0.11,
                    3500: 0.13,
                    4000: 0.14,
                    4500: 0.15  }

        else:
            print self.signal_model
            raw_input("failed to find the signal model")

        signal_uncertainty_PDF  = table_signal_pdfuncertainty[self.signal_mass]
        signal_uncertainty_scale = table_signal_scaleuncertainty[self.signal_mass]

        self.XS_Signal_uncertainty = (signal_uncertainty_PDF**2+ signal_uncertainty_scale**2)**0.5 #pdf uncertainty 13% + scale uncertainty 11%
        self.XS_STop_uncertainty = 0.050
        self.XS_VV_uncertainty   = 0.20
        #self.XS_TTbar_uncertainty   = 0.12 #0.20 #qun
        #self.XS_TTbar_NLO_uncertainty = 0.063 # from AN-12-368 table8
        #self.XS_STop_NLO_uncertainty  = 0.05 # from AN-12-368 table8
        #self.XS_VV_NLO_uncertainty    = 0.10 # from AN-12-368 table8

        #el and mu trigger and eff uncertainty, B2G-16-4
        self.lep_trigger_uncertainty = 0.05
        self.btag_scale_ttbar_uncertainty  = 0.064
        self.btag_scale_singletop_uncertainty  = 0.051
        self.btag_scale_wjets_uncertainty  = 0.0 #need calculate according to the Norm of WJets, TTbar, and Single Top
        self.btag_scale_signal_uncertainty = 0.006

        if self.channel == "mu":
            self.signal_lepton_energy_scale_uncertainty = 0.007
            self.signal_lepton_energy_res_uncertainty   = 0.001
            self.lep_eff_uncertainty     = 0.05
        elif self.channel == "el":
            self.signal_lepton_energy_scale_uncertainty = 0.002
            self.signal_lepton_energy_res_uncertainty   = 0.001
            self.lep_eff_uncertainty     = 0.05
        elif self.channel == "em" and options.combine:
            print "for combine"
        else:
            raw_input("Fail to find correct channel. Please check your channel:%s"%(self.channel))

        self.eff_vtag_model = 0.

        self.signal_jet_energy_res_uncertainty      = 0.003
        self.signal_jet_energy_scale_uncertainty_low = 0.982
        self.signal_jet_energy_scale_uncertainty_high = 1.012
        self.signal_jet_mass_scale_uncertainty_low = 0.965
        self.signal_jet_mass_scale_uncertainty_high = 1.020
        self.signal_jet_mass_res_uncertainty_low = 1.006
        self.signal_jet_mass_res_uncertainty_high = 0.985

        #### sigma and mean signal systematic inflation
        self.mean_signal_uncertainty_jet_scale  = 0.005
        self.mean_signal_uncertainty_lep_scale  = 0.005
        self.sigma_signal_uncertainty_jet_scale = 0.0375
        self.sigma_signal_uncertainty_jet_res   = 0.0325
        self.sigma_signal_uncertainty_lep_scale = 0.0175

        #### Set systematic on the Wjets shape   and TTbar due to PS, fitting function etc..
        self.shape_para_error_WJets0 = 1.4
        self.shape_para_error_alpha  = 1.4
        self.shape_para_error_TTbar  = 2.0
        self.shape_para_error_VV     = 1.
        self.shape_para_error_STop   = 1.

        # shape parameter uncertainty
        self.FloatingParams = RooArgList("floatpara_list")

    #### Method to make a RooAbsPdf giving label, model name, spectrum, if it is mc or not and a constraint list for the parameters
    def make_Pdf(self, label, in_model_name, mass_spectrum= "_mj", ConstraintsList=[], ismc = 0):
        if TString(mass_spectrum).Contains("_mj"):
            rrv_x = self.workspace4fit_.var("rrv_mass_j")
        if TString(mass_spectrum).Contains("_mlvj"):
            rrv_x = self.workspace4fit_.var("rrv_mass_lvj")

        # W mass: 80.385
        if in_model_name == "Voig":
            print "########### Voigtian Pdf for mJ ############"
            rrv_mean_voig = RooRealVar("rrv_mean_voig"+label+"_"+self.channel, "rrv_mean_voig"+label+"_"+self.channel, 84, 78, 88)
            rrv_width_voig = RooRealVar("rrv_width_voig"+label+"_"+self.channel, "rrv_width_voig"+label+"_"+self.channel, 7., 1, 40)
            rrv_sigma_voig = RooRealVar("rrv_sigma_voig"+label+"_"+self.channel, "rrv_sigma_voig"+label+"_"+self.channel, 5, 0.01, 20)
            model_pdf = RooVoigtian("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_mean_voig, rrv_width_voig, rrv_sigma_voig)

        # Higgs mass 600-1000
        if in_model_name == "Voig_mlvj":
            print "########### Voigtian Pdf for Higgs mlvj ############"
            rrv_mean_voig = RooRealVar("rrv_mean_voig"+label+"_"+self.channel, "rrv_mean_voig"+label+"_"+self.channel, 650, 550, 1200)
            rrv_width_voig = RooRealVar("rrv_width_voig"+label+"_"+self.channel, "rrv_width_voig"+label+"_"+self.channel, 100., 10, 600)
            rrv_sigma_voig = RooRealVar("rrv_sigma_voig"+label+"_"+self.channel, "rrv_sigma_voig"+label+"_"+self.channel, 200, 10, 400)
            model_pdf = RooVoigtian("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_mean_voig, rrv_width_voig, rrv_sigma_voig)

        ## BW for the W mass peak
        if in_model_name == "BW":
            print "########### BW Pdf for mj fit ############"
            rrv_mean_BW = RooRealVar("rrv_mean_BW"+label+"_"+self.channel, "rrv_mean_BW"+label+"_"+self.channel, 84, 78, 88)
            rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.channel, "rrv_width_BW"+label+"_"+self.channel, 20, 1, 40)
            model_pdf = RooBreitWigner("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_mean_BW, rrv_width_BW)

        ##  Voig for W mass peak
        if in_model_name == "2Voig":

            print "########### Double Voigtian for mj fit ############"
            rrv_mean_voig    = RooRealVar("rrv_mean_voig"+label+"_"+self.channel, "rrv_mean_voig"+label+"_"+self.channel, 84, 78, 88)#W mass 80.385
            rrv_shift_2Voig  = RooRealVar("rrv_shift_2Voig"+label+"_"+self.channel, "rrv_shift_https://github.com/zixu/EXO2016/tree/master/DoFit2Voig"+label+"_"+self.channel, 10.8026)# Z mass: 91.1876 shift=91.1876-80.385=10.8026
            rrv_mean_shifted = RooFormulaVar("rrv_mean_voig2"+label+"_"+self.channel, "@0+@1", RooArgList(rrv_mean_voig, rrv_shift_2Voig))

            rrv_width_voig = RooRealVar("rrv_width_voig"+label+"_"+self.channel, "rrv_width_voig"+label+"_"+self.channel, 16., 6, 26)
            rrv_sigma_voig = RooRealVar("rrv_sigma_voig"+label+"_"+self.channel, "rrv_sigma_voig"+label+"_"+self.channel, 5., 0., 10.)

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel, "rrv_frac"+label+"_"+self.channel, 0.8, 0.5, 1.)

            model_voig1 = RooVoigtian("model_voig1"+label+"_"+self.channel+mass_spectrum, "model_voig1"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_mean_voig, rrv_width_voig, rrv_sigma_voig)

            model_voig2 = RooVoigtian("model_voig2"+label+"_"+self.channel+mass_spectrum, "model_voig2"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_mean_shifted, rrv_width_voig, rrv_sigma_voig)
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, RooArgList(model_voig1, model_voig2), RooArgList(rrv_frac))

        ## Gaus for the W peak
        if in_model_name == "Gaus":
            print "########### Gaus for W peak  ############"
            rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel, "rrv_mean_gaus"+label+"_"+self.channel, 84, 78, 88)
            rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel, "rrv_sigma_gaus"+label+"_"+self.channel, 7, 1, 15)
            model_pdf = RooGaussian("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_mean_gaus, rrv_sigma_gaus)

        ## Gaus for the signal lineshape
        if in_model_name == "Gaus_mlvj":
            print "########### Gaus for Higgs mlvj ############"
            rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel, "rrv_mean_gaus"+label+"_"+self.channel, 1000, 800, 1100)
            rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel, "rrv_sigma_gaus"+label+"_"+self.channel, 200, 50, 300)
            model_pdf = RooGaussian("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_mean_gaus, rrv_sigma_gaus)

        ## Crystal Ball for the W mass peak
        if in_model_name == "CB":
            print "########### Cystal Ball for mj fit ############"
            rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel, "rrv_mean_CB"+label+"_"+self.channel, 84, 78, 88)
            rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel, "rrv_sigma_CB"+label+"_"+self.channel, 7, 4, 10)
            rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel, "rrv_alpha_CB"+label+"_"+self.channel,-2,-4,-0.5)
            rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.channel, "rrv_n_CB"+label+"_"+self.channel, 2, 0., 4)
            model_pdf = RooCBShape("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_mean_CB, rrv_sigma_CB, rrv_alpha_CB, rrv_n_CB)

        ## Crystal  ball shape for Bulk GR samples and signal
        if in_model_name == "CB_mlvj":
            print "########### Crystal Ball for Higgs and  Bulk GR  mlvj ############"
            label_tstring = TString(label)
            rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 700, 550, 2500)
            rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 50, 20 , 120)
            rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 4, 1, 5)
            rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 20., 10, 40)
            model_pdf = RooCBShape("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_mean_CB, rrv_sigma_CB, rrv_alpha_CB, rrv_n_CB)


        ## Crystal  ball shape for Bulk GR samples and signal
        if in_model_name == "DoubleCB_mlvj":
            label_tstring = TString(label)
            print "########### Double CB for Bulk graviton mlvj ############"

            if label_tstring.Contains(self.signal_model+"600"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 621, 550, 680)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 47.9, 40, 60)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 31)#, 15, 60)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 5)#, 0.1, 10.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2.2, 0.01, 35)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1.94, 1.0, 4.0)
            elif label_tstring.Contains(self.signal_model+"700"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 720, 650, 780)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 55, 30 , 70)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 45.)#, 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2, 0.1, 4.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2.9, 1., 5.)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2, 0.1, 4.)
            elif label_tstring.Contains(self.signal_model+"750"):#because the M_ljv lower limit is 800GeV, so, just use right side
                #rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 777, 700, 820)
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 777, 700, 950)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 59.4, 50, 70)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 45.)#, 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1.66, 1., 3.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 4., 1, 10)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1.8, 0.5, 5.)
            elif label_tstring.Contains(self.signal_model+"800"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 826, 770, 900)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 61.2, 50, 80)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 45)#4., 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1.52, 1, 3.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 5.4, 1, 10)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1.7, 1.0, 3.)
            elif label_tstring.Contains(self.signal_model+"900"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 929, 870, 970)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 66.2, 50, 80)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 45)#., 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1.4, 1., 3.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 5., 1, 10)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1.8, 1., 3.)
            elif label_tstring.Contains(self.signal_model+"1000"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1029, 800, 1270)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 71.6, 60 , 80)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 45., 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1.36, 0.5, 3.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 3.1, 0.1, 8)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2.15, 0.5, 5.)

            elif label_tstring.Contains(self.signal_model+"1200"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1229, 970, 2070)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 71.6, 60 , 80)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 45., 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1.36, 0.5, 3.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 3.1, 0.1, 8)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2.15, 0.5, 5.)


            elif label_tstring.Contains(self.signal_model+"1400"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1429, 970, 2070)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 71.6, 60 , 80)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 45., 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1.36, 0.5, 3.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 3.1, 0.1, 8)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2.15, 0.5, 5.)


            elif label_tstring.Contains(self.signal_model+"1600"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1629, 970, 2070)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 71.6, 60 , 80)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 45., 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1.36, 0.5, 3.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 3.1, 0.1, 8)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2.15, 0.5, 5.)


            elif label_tstring.Contains(self.signal_model+"1800"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1829, 970, 2070)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 71.6, 60 , 80)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 45., 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1.36, 0.5, 3.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 3.1, 0.1, 8)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2.15, 0.5, 5.)

            elif label_tstring.Contains(self.signal_model+"2000"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2025, 1500, 3500)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 100, 50 , 150)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 4., 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 5, 0.1, 10.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 4., 0.01, 35)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 3., 0.5, 6.)

            elif label_tstring.Contains(self.signal_model+"2500"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2500, 2000, 3000)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 100, 20 , 300)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 4., 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 5, 0.1, 10.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2., 0.01, 35)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 3., 0.5, 6.)

            elif label_tstring.Contains(self.signal_model+"3000"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 3000, 2500, 3500)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 100, 20 , 300)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 4., 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 5, 0.1, 10.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2., 0.01, 35)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 3., 0.5, 6.)

            elif label_tstring.Contains(self.signal_model+"3500"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 3500, 2500, 4000)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 150, 100 , 300)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 4., 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 5, 0.1, 10.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2., 0.01, 35)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 3., 0.5, 6.)

            elif label_tstring.Contains(self.signal_model+"4000"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 4000, 3500, 4500)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 160, 100 , 300)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 4., 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 5, 0.1, 10.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2., 0.01, 35)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 3., 0.5, 6.)

            elif label_tstring.Contains(self.signal_model+"4500"):#because the M_ljv lower limit is 800GeV, so, just use right side
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 4545, 4000, 5000)
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 155, 100 , 300)
                rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 7.5, 0.01, 45)
                rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1, 0.1, 10.)
                rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2.7, 0.01, 35)
                rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2., 0.5, 6.)

            else :
                raw_input("Please check your signal or backgroud: "+label)
                #rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 1000, 500, 5500)
                #rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 50, 20 , 120)
                #rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 22., 0.01, 45)
                #rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 5, 0.1, 10.)
                #rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 2., 0.01, 35)
                #rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_category, 3., 0.5, 6.)

            if options.closuretest>0:
                #rrv_mean_CB.setConstant(kTRUE)
                rrv_sigma_CB.setConstant(kTRUE)
                rrv_n1_CB.setConstant(kTRUE)
                rrv_n2_CB.setConstant(kTRUE)
                rrv_alpha1_CB.setConstant(kTRUE)
                rrv_alpha2_CB.setConstant(kTRUE)

            rrv_mean_scale_p1 = RooRealVar("CMS_sig_p1_jes", "CMS_sig_p1_jes", 0)
            rrv_mean_scale_p1.setConstant(kTRUE)
            if self.channel == "mu" :
                rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_m", "CMS_sig_p1_scale_m", 0)
                rrv_mean_scale_p2.setConstant(kTRUE)
            elif self.channel == "el" :
                rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_e", "CMS_sig_p1_scale_e", 0)
                rrv_mean_scale_p2.setConstant(kTRUE)
            elif self.channel == "em":
                rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_em", "CMS_sig_p1_scale_em", 0)
                rrv_mean_scale_p2.setConstant(kTRUE)

            rrv_mean_scale_X1 = RooRealVar("rrv_mean_shift_scale_lep"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_shift_scale_lep"+label+"_"+self.channel+"_"+self.wtagger_category, float(self.mean_signal_uncertainty_lep_scale))
            rrv_mean_scale_X1.setConstant(kTRUE)
            rrv_mean_scale_X2 = RooRealVar("rrv_mean_shift_scale_jes"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_mean_shift_scale_jes"+label+"_"+self.channel+"_"+self.wtagger_category, float(self.mean_signal_uncertainty_jet_scale))
            rrv_mean_scale_X2.setConstant(kTRUE)

            rrv_total_mean_CB = RooFormulaVar("rrv_total_mean_CB"+label+"_"+self.channel, "@0*(1+@1*@2)*(1+@3*@4)", RooArgList(rrv_mean_CB, rrv_mean_scale_p1, rrv_mean_scale_X1, rrv_mean_scale_p2, rrv_mean_scale_X2))

            if self.channel == "mu":
                rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_m", "CMS_sig_p2_scale_m", 0)
                rrv_sigma_scale_p1.setConstant(kTRUE)
            elif self.channel == "el":
                rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_e", "CMS_sig_p2_scale_e", 0)
                rrv_sigma_scale_p1.setConstant(kTRUE)
            elif self.channel == "em":
                rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_em", "CMS_sig_p2_scale_em", 0)
                rrv_sigma_scale_p1.setConstant(kTRUE)

            rrv_sigma_scale_p2 = RooRealVar("CMS_sig_p2_jer", "CMS_sig_p2_jer", 0)
            rrv_sigma_scale_p3 = RooRealVar("CMS_sig_p2_jes", "CMS_sig_p2_jes", 0)
            rrv_sigma_scale_p2.setConstant(kTRUE)
            rrv_sigma_scale_p3.setConstant(kTRUE)

            rrv_mean_sigma_X1 = RooRealVar("rrv_sigma_shift_lep_scale"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_shift_scale"+label+"_"+self.channel+"_"+self.wtagger_category, float(self.sigma_signal_uncertainty_lep_scale))
            rrv_mean_sigma_X2 = RooRealVar("rrv_sigma_shift_jes"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_shift_scale"+label+"_"+self.channel+"_"+self.wtagger_category, float(self.sigma_signal_uncertainty_jet_scale))
            rrv_mean_sigma_X3 = RooRealVar("rrv_sigma_shift_res"+label+"_"+self.channel+"_"+self.wtagger_category, "rrv_sigma_shift_res"+label+"_"+self.channel+"_"+self.wtagger_category, float(self.sigma_signal_uncertainty_jet_res))
            rrv_mean_sigma_X1.setConstant(kTRUE)
            rrv_mean_sigma_X2.setConstant(kTRUE)
            rrv_mean_sigma_X3.setConstant(kTRUE)

            rrv_total_sigma_CB = RooFormulaVar("rrv_total_sigma_CB"+label+"_"+self.channel, "@0*(1+@1*@2)*(1+@3*@4)*(1+@5*@6)", RooArgList(rrv_sigma_CB, rrv_sigma_scale_p1, rrv_mean_sigma_X1, rrv_sigma_scale_p2, rrv_mean_sigma_X2, rrv_sigma_scale_p3, rrv_mean_sigma_X3))

            model_pdf = ROOT.RooDoubleCrystalBall("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_total_mean_CB, rrv_total_sigma_CB, rrv_alpha1_CB, rrv_n1_CB, rrv_alpha2_CB, rrv_n2_CB)

        ## ExpN pdf for W+jets bkg fit
        if in_model_name == "ExpN":

            print "########### ExpN funtion for W+jets mlvj ############"
            ##rrv_c_ExpN = RooRealVar("rrv_c_ExpN"+label+"_"+self.channel, "rrv_c_ExpN"+label+"_"+self.channel,-5.32e-3,-1e-1,-1e-6)
            ##rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel, "rrv_n_ExpN"+label+"_"+self.channel, -64.3, -1e4, 1e5)
            rrv_c_ExpN = RooRealVar("rrv_c_ExpN"+label+"_"+self.channel, "rrv_c_ExpN"+label+"_"+self.channel,-3e-3,-1e-1,-1e-6)
            rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel, "rrv_n_ExpN"+label+"_"+self.channel, 2000, -1e5, 1e5)
            #rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel, "rrv_n_ExpN"+label+"_"+self.channel, 2000, 0, 1e5)

            model_pdf = ROOT.RooExpNPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c_ExpN, rrv_n_ExpN)


        ## levelled exp for W+jets bkg fit
        if in_model_name == "ExpTail":
            print "########### ExpTai = levelled exp funtion for W+jets mlvj ############"
            label_tstring = TString(label)
            ##if self.wtagger_category == "LP":
            ##    rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel, "rrv_s_ExpTail"+label+"_"+self.channel, 250,-1.e6, 1e6)
            ##    rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel, "rrv_a_ExpTail"+label+"_"+self.channel, 1e-1,-1.e-2, 1e6)
            ##else:
            ##    rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel, "rrv_s_ExpTail"+label+"_"+self.channel, 100, 0., 500)
            ##    rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel, "rrv_a_ExpTail"+label+"_"+self.channel, 0.5,-1, 5)
            ##rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel, "rrv_s_ExpTail"+label+"_"+self.channel, 200, 0., 500)
            ##rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel, "rrv_a_ExpTail"+label+"_"+self.channel, 0.,-1, 1)
            ##model_pdf     = ROOT.RooExpTailPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_s_ExpTail, rrv_a_ExpTail)

            rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel, "rrv_s_ExpTail"+label+"_"+self.channel, 2,-10., 10)
            rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel, "rrv_a_ExpTail"+label+"_"+self.channel, 0.,-10, 10)
            rrv_sp_ExpTail      = RooFormulaVar("rrv_sp_ExpTail"+label+"_"+self.channel, "@0*100.", RooArgList(rrv_s_ExpTail))
            rrv_ap_ExpTail      = RooFormulaVar("rrv_ap_ExpTail"+label+"_"+self.channel, "@0/100.", RooArgList(rrv_a_ExpTail))
            model_pdf     = ROOT.RooExpTailPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_sp_ExpTail, rrv_ap_ExpTail)


        ## sum of two exponential
        if in_model_name == "2Exp":
            print "########### 2Exp = levelled exp funtion for W+jets mlvj ############"
            rrv_c0_2Exp   = RooRealVar("rrv_c0_2Exp"+label+"_"+self.channel, "rrv_c0_2Exp"+label+"_"+self.channel, -5e-3, -8e-3,-4e-3)
            rrv_c1_2Exp   = RooRealVar("rrv_c1_2Exp"+label+"_"+self.channel, "rrv_c1_2Exp"+label+"_"+self.channel, -1e-3, -4e-3,-1e-4)
            rrv_frac_2Exp = RooRealVar("rrv_frac_2Exp"+label+"_"+self.channel, "rrv_frac_2Exp"+label+"_"+self.channel, 0., 0., 1e-2)
            model_pdf  = ROOT.Roo2ExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c0_2Exp, rrv_c1_2Exp, rrv_frac_2Exp)

        ## sum of two exponential
        if in_model_name == "Exp" or in_model_name == "Exp_sr":
            print "########### Exp = levelled exp funtion for W+jets mlvj ############"
            rrv_c_Exp = RooRealVar("rrv_c_Exp"+label+"_"+self.channel, "rrv_c_Exp"+label+"_"+self.channel,-0.05,-0.2, 0.)
            model_pdf = ROOT.RooExponential("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c_Exp)

        ## Erf times for mj spectrum
        if in_model_name == "ErfExp" :
            print "########### Erf*Exp for mj fit  ############"
            ##rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel,-1.4e-2,-0.1,-1e-4)
            ##rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel, "rrv_offset_ErfExp"+label+"_"+self.channel, 60., 40., 80)
            ##rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, 20., 10., 40.)
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel,-4.4e-2,-0.1,-1e-4)
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel, "rrv_offset_ErfExp"+label+"_"+self.channel, 74., 50., 90)
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, 30., 15., 50.)
            model_pdf         = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c_ErfExp, rrv_offset_ErfExp, rrv_width_ErfExp)

        ## different initial values -> for mlvj
        if in_model_name == "ErfExp_v1" :
            print "########### Erf*Exp for mlvj fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel,-0.006,-0.1, 0.)
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel, "rrv_offset_ErfExp"+label+"_"+self.channel, 450., 400., 550.)
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, 70., 10, 100.)
            model_pdf         = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c_ErfExp, rrv_offset_ErfExp, rrv_width_ErfExp)

        ## different initial values -> for mlvj
        if in_model_name == "ErfExp_v2" :
            print "########### Erf*Exp for mlvj fit  ############"
            rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel,-0.005,-0.1, 0.)
            rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel, "rrv_offset_ErfExp"+label+"_"+self.channel, 450., 400., 500.)
            rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, 50., 10, 100.)
            rrv_residue_ErfExp = RooRealVar("rrv_residue_ErfExp"+label+"_"+self.channel, "rrv_residue_ErfExp"+label+"_"+self.channel, 0., 0., 1.)
            model_pdf = RooGenericPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, "(TMath::Exp(%s*%s) + %s)*(1.+TMath::Erf((%s-%s)/%s))/2. "%(rrv_c_ErfExp.GetName(), rrv_x.GetName(), rrv_residue_ErfExp.GetName(), rrv_x.GetName(), rrv_offset_ErfExp.GetName(), rrv_width_ErfExp.GetName()), RooArgList(rrv_x, rrv_c_ErfExp, rrv_offset_ErfExp, rrv_width_ErfExp, rrv_residue_ErfExp))

        ## different initial values -> for mlvj
        if in_model_name == "ErfExp_v3" : #different init-value and range
            print "########### Erf*Exp for mlvj fit  ############"
            rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel,-0.005,-0.1, 0.)
            rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel, "rrv_offset_ErfExp"+label+"_"+self.channel, 450., 400, 500.)
            rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, 50., 10, 100.)
            rrv_residue_ErfExp = RooRealVar("rrv_residue_ErfExp"+label+"_"+self.channel, "rrv_residue_ErfExp"+label+"_"+self.channel, 0., 0., 1.)
            rrv_high_ErfExp    = RooRealVar("rrv_high_ErfExp"+label+"_"+self.channel, "rrv_high_ErfExp"+label+"_"+self.channel, 1., 0., 400)
            rrv_high_ErfExp.setConstant(kTRUE)
            model_pdf = RooGenericPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, "(TMath::Exp(%s*%s) + %s)* TMath::Power(((1+TMath::Erf((%s-%s)/%s))/2.), %s)"%(rrv_c_ErfExp.GetName(), rrv_x.GetName(), rrv_residue_ErfExp.GetName(), rrv_x.GetName(), rrv_offset_ErfExp.GetName(), rrv_width_ErfExp.GetName(), rrv_high_ErfExp.GetName()), RooArgList(rrv_x, rrv_c_ErfExp, rrv_offset_ErfExp, rrv_high_ErfExp, rrv_width_ErfExp, rrv_residue_ErfExp))

        ## Exp+Gaus or mj spectrum
        if in_model_name == "ExpGaus":
            print "########### Exp + Gaus for mj  fit  ############"
            rrv_c_Exp       = RooRealVar("rrv_c_Exp"+label+"_"+self.channel, "rrv_c_Exp"+label+"_"+self.channel, 0.05,-0.2, 0.2)
            exp             = ROOT.RooExponential("exp"+label+"_"+self.channel, "exp"+label+"_"+self.channel, rrv_x, rrv_c_Exp)

            ##rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel, "rrv_mean1_gaus"+label+"_"+self.channel, 84, 78, 88)
            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel, "rrv_mean1_gaus"+label+"_"+self.channel, 84, 65, 88)
            rrv_sigma1_gaus = RooRealVar("rrv_smgma1_gaus"+label+"_"+self.channel, "rrv_sigma1_gaus"+label+"_"+self.channel, 7, 4, 10)
            rrv_high        = RooRealVar("rrv_high"+label+"_"+self.channel, "rrv_high"+label+"_"+self.channel, 0.5, 0., 1.)
            gaus            = RooGaussian("gaus"+label+"_"+self.channel, "gaus"+label+"_"+self.channel, rrv_x, rrv_mean1_gaus, rrv_sigma1_gaus)

            model_pdf       = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, RooArgList(exp, gaus), RooArgList(rrv_high))

        ## Erf*Exp + Gaus for mj spectrum
        if in_model_name == "ErfExpGaus":
            print "########### Erf*Exp + Gaus for mj  fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel,-0.05,-0.4, 0.)
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel, "rrv_offset_ErfExp"+label+"_"+self.channel, 100., 10., 300.)
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, 30., 10, 100.)

            erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel, "erfExp"+label+"_"+self.channel, rrv_x, rrv_c_ErfExp, rrv_offset_ErfExp, rrv_width_ErfExp)

            rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel, "rrv_mean_gaus"+label+"_"+self.channel, 82, 78, 87)
            rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel, "rrv_sigma_gaus"+label+"_"+self.channel, 7, 4, 10)
            gaus = RooGaussian("gaus"+label+"_"+self.channel, "gaus"+label+"_"+self.channel, rrv_x, rrv_mean_gaus, rrv_sigma_gaus)

            rrv_high  = RooRealVar("rrv_high"+label+"_"+self.channel, "rrv_high"+label+"_"+self.channel, 0.7, 0., 1.)
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, RooArgList(erfExp, gaus), RooArgList(rrv_high))

        ## Erf*Exp + Gaus for mj spectrum with offset == mean
        if in_model_name == "ErfExpGaus_sp":
            print "########### Erf*Exp + Gaus for mj  fit  ############"
            rrv_c_ErfExp     = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel,-0.05,-0.2, 0.)
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, 30., 10, 200.)
            erfExp           = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel, "erfExp"+label+"_"+self.channel, rrv_x, rrv_c_ErfExp, rrv_mean1_gaus, rrv_width_ErfExp)

            rrv_mean1_gaus   = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel, "rrv_mean1_gaus"+label+"_"+self.channel, 84, 78, 88)
            rrv_sigma1_gaus  = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel, "rrv_sigma1_gaus"+label+"_"+self.channel, 7, 4, 10)
            gaus             = RooGaussian("gaus"+label+"_"+self.channel, "gaus"+label+"_"+self.channel, rrv_x, rrv_mean1_gaus, rrv_sigma1_gaus)

            rrv_high  = RooRealVar("rrv_high"+label+"_"+self.channel, "rrv_high"+label+"_"+self.channel, 0.5, 0., 1.)
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, RooArgList(erfExp, gaus), RooArgList(rrv_high))

        ## Erf*Exp+Gaus or mj spectrum
        if in_model_name == "ErfExpGaus_v0":
            print "########### Erf*Exp + Gaus for mj  fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel,-0.05,-0.2, 0.)
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel, "rrv_offset_ErfExp"+label+"_"+self.channel, 100., 10., 140.)
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, 30., 10, 100.)
            erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel, "erfExp"+label+"_"+self.channel, rrv_x, rrv_c_ErfExp, rrv_offset_ErfExp, rrv_width_ErfExp)

            rrv_mean_gaus     = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel, "rrv_mean_gaus"+label+"_"+self.channel, 84, 78, 88)
            rrv_sigma_gaus    = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel, "rrv_sigma_gaus"+label+"_"+self.channel, 7, 4, 10)
            gaus              = RooGaussian("gaus"+label+"_"+self.channel, "gaus"+label+"_"+self.channel, rrv_x, rrv_mean_gaus, rrv_sigma_gaus)

            rrv_high   = RooRealVar("rrv_high"+label+"_"+self.channel, "rrv_high"+label+"_"+self.channel, 0.7, 0., 1.)
            model_pdf  = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, RooArgList(erfExp, gaus), RooArgList(rrv_high))

        ## Erf*Exp+Gaus or mj spectrum
        if in_model_name == "ErfExpGaus_v1":
            print "########### Erf*Exp + Gaus for mlvj fit  ############"
            rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel,-0.007,-0.1, 0.)
            rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel, "rrv_offset_ErfExp"+label+"_"+self.channel, 800., 10., 1400.)
            rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, 24., 10, 150.)
            erfExp             = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel, "erfExp"+label+"_"+self.channel, rrv_x, rrv_c_ErfExp, rrv_offset_ErfExp, rrv_width_ErfExp)

            rrv_mean_gaus   = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel, "rrv_mean_gaus"+label+"_"+self.channel, 700, 500, 1200)
            rrv_sigma_gaus  = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel, "rrv_sigma_gaus"+label+"_"+self.channel, 150, 10, 300)
            gaus            = RooGaussian("gaus"+label+"_"+self.channel, "gaus"+label+"_"+self.channel, rrv_x, rrv_mean_gaus, rrv_sigma_gaus)

            rrv_high  = RooRealVar("rrv_high"+label+"_"+self.channel, "rrv_high"+label+"_"+self.channel, 0.1, 0., 1.)
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, RooArgList(erfExp, gaus), RooArgList(rrv_high))

        ## Erf*Exp+Gaus or mj spectrum
        if in_model_name == "ErfExpGaus_sp_v1":
            print "########### Erf*Exp + Gaus for mlvj fit  ############"
            rrv_c_ErfExp     = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel,-0.007,-0.1, 0.)
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, 24., 10, 150.)
            rrv_mean_gaus    = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel, "rrv_mean_gaus"+label+"_"+self.channel, 900, 860, 1200)
            erfExp           = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel, "erfExp"+label+"_"+self.channel, rrv_x, rrv_c_ErfExp, rrv_mean_gaus, rrv_width_ErfExp)

            rrv_sigma_gaus   = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel, "rrv_sigma_gaus"+label+"_"+self.channel, 150, 10, 300)
            gaus = RooGaussian("gaus"+label+"_"+self.channel, "gaus"+label+"_"+self.channel, rrv_x, rrv_mean_gaus, rrv_sigma_gaus)

            rrv_high  = RooRealVar("rrv_high"+label+"_"+self.channel, "rrv_high"+label+"_"+self.channel, 0.1, 0., 1.)
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, RooArgList(erfExp, gaus), RooArgList(rrv_high))

        ## Erf*Exp+Gaus or mj spectrum
        if in_model_name == "ErfExpGaus_v2":
            print "########### Erf*Exp + Gaus for mj fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel,-0.05,-10., 0.)
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel, "rrv_offset_ErfExp"+label+"_"+self.channel, 100., 10., 140.)
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, 30., 10, 100.)
            rrv_mean_gaus     = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel, "rrv_mean_gaus"+label+"_"+self.channel, 84, 78, 88)
            rrv_sigma_gaus    = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel, "rrv_sigma_gaus"+label+"_"+self.channel, 7, 4, 10)
            rrv_high          = RooRealVar("rrv_high"+label+"_"+self.channel, "rrv_high"+label+"_"+self.channel, 200., 0., 1000.)
            model_pdf = ROOT.RooErfExp_Gaus_Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c_ErfExp, rrv_offset_ErfExp, rrv_width_ErfExp, rrv_mean_gaus, rrv_sigma_gaus, rrv_high)

        ## Erf*Exp + 2Gaus
        if in_model_name == "ErfExp2Gaus":
            print "########### Erf*Exp + 2Gaus for mj fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel,-0.05,-0.2, 0.)
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel, "rrv_offset_ErfExp"+label+"_"+self.channel, 100., 10., 140.)
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, 30., 10, 100.)
            erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel, "erfExp"+label+"_"+self.channel, rrv_x, rrv_c_ErfExp, rrv_offset_ErfExp, rrv_width_ErfExp)

            rrv_mean1_gaus   = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel, "rrv_mean1_gaus"+label+"_"+self.channel, 84, 78, 88)
            rrv_mean2_gaus   = RooRealVar("rrv_mean2_gaus"+label+"_"+self.channel, "rrv_mean2_gaus"+label+"_"+self.channel, 180, 170, 190)
            rrv_sigma1_gaus  = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel, "rrv_sigma1_gaus"+label+"_"+self.channel, 7, 4, 10)
            rrv_sigma2_gaus  = RooRealVar("rrv_sigma2_gaus"+label+"_"+self.channel, "rrv_sigma2_gaus"+label+"_"+self.channel, 10, 7, 15)
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel, "gaus1"+label+"_"+self.channel, rrv_x, rrv_mean1_gaus, rrv_sigma1_gaus)
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel, "gaus2"+label+"_"+self.channel, rrv_x, rrv_mean2_gaus, rrv_sigma2_gaus)

            rrv_high1 = RooRealVar("rrv_high1"+label+"_"+self.channel, "rrv_high1"+label+"_"+self.channel, 0.6, 0., 1.)
            rrv_high2 = RooRealVar("rrv_high2"+label+"_"+self.channel, "rrv_high2"+label+"_"+self.channel, 0.4, 0., 1.)
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, RooArgList(erfExp, gaus1, gaus2), RooArgList(rrv_high1, rrv_high2))

        ## Gaus + Gaus for mj spectrum
        if in_model_name == "2Gaus":
            print "########### 2Gaus for mj fit  ############"
            mean1_tmp      = 8.3141e+01
            mean1_tmp_err      = 1.63e-01
            deltamean_tmp  = 6.9129e+00
            deltamean_tmp_err  = 1.24e+00
            sigma1_tmp     = 7.5145e+00
            sigma1_tmp_err     = 1.99e-01
            scalesigma_tmp = 3.6819e+00
            scalesigma_tmp_err = 2.11e-01
            frac_tmp       = 6.7125e-01
            frac_tmp_err       = 2.09e-02

            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel, "rrv_mean1_gaus"+label+"_"+self.channel, mean1_tmp, 70, 90)
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel, "rrv_sigma1_gaus"+label+"_"+self.channel, sigma1_tmp, 7, 10)
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel, "gaus1"+label+"_"+self.channel, rrv_x, rrv_mean1_gaus, rrv_sigma1_gaus)

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel, "rrv_deltamean_gaus"+label+"_"+self.channel, 0.)
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel, "@0+@1", RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus))
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel, "rrv_scalesigma_gaus"+label+"_"+self.channel, scalesigma_tmp, 2.8, 4.4)
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel, "@0*@1", RooArgList(rrv_sigma1_gaus, rrv_scalesigma_gaus))
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel, "gaus2"+label+"_"+self.channel, rrv_x, rrv_mean2_gaus, rrv_sigma2_gaus)

            rrv_frac  = RooRealVar("rrv_frac"+label+"_"+self.channel, "rrv_frac"+label+"_"+self.channel, frac_tmp, 0.6, 0.75)
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, RooArgList(gaus1, gaus2), RooArgList(rrv_frac), 1)

        ## 2Gaus+2Gaus for VV mj spectrum -> WZ and WW
        if in_model_name == "2_2Gaus":

            print "########### 2Gaus +2Gaus for mj fit  ############"
            mean1_tmp      = 8.3141e+01 #mean1_tmp_err      = 1.63e-01
            deltamean_tmp  = 9.0e+00 #deltamean_tmp_err  = 1.24e+00
            sigma1_tmp     = 7.5145e+00 #sigma1_tmp_err     = 1.99e-01
            scalesigma_tmp = 3.6819e+00 #scalesigma_tmp_err = 2.11e-01
            frac_tmp       = 6.7125e-01 #frac_tmp_err       = 2.09e-02

            rrv_shift = RooRealVar("rrv_shift"+label+"_"+self.channel, "rrv_shift"+label+"_"+self.channel, 10.8026) # Z mass: 91.1876 shift=91.1876-80.385=10.8026

            rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel, "rrv_mean1_gaus"+label+"_"+self.channel, mean1_tmp, mean1_tmp-20, mean1_tmp+4)
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel, "rrv_sigma1_gaus"+label+"_"+self.channel, sigma1_tmp, 0, sigma1_tmp*2)
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel, "gaus1"+label+"_"+self.channel, rrv_x, rrv_mean1_gaus, rrv_sigma1_gaus)

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel, "rrv_deltamean_gaus"+label+"_"+self.channel, deltamean_tmp, deltamean_tmp*-1, deltamean_tmp*3)
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel, "@0+@1", RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus))
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel, "rrv_scalesigma_gaus"+label+"_"+self.channel, scalesigma_tmp, 0, scalesigma_tmp*2)
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel, "@0*@1", RooArgList(rrv_sigma1_gaus, rrv_scalesigma_gaus))
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel, "gaus2"+label+"_"+self.channel, rrv_x, rrv_mean2_gaus, rrv_sigma2_gaus)

            rrv_frac1 = RooRealVar("rrv_frac1"+label+"_"+self.channel, "rrv_frac1"+label+"_"+self.channel, frac_tmp, 0, frac_tmp*2)
            gausguas_1 = RooAddPdf("gausguas_1"+label+"_"+self.channel+mass_spectrum, "gausguas_1"+label+"_"+self.channel+mass_spectrum, RooArgList(gaus1, gaus2), RooArgList(rrv_frac1), 1)

            rrv_mean3_gaus = RooFormulaVar("rrv_mean3_gaus"+label+"_"+self.channel, "@0+@1", RooArgList(rrv_mean1_gaus, rrv_shift))
            rrv_mean4_gaus = RooFormulaVar("rrv_mean4_gaus"+label+"_"+self.channel, "@0+@1", RooArgList(rrv_mean2_gaus, rrv_shift))
            gaus3 = RooGaussian("gaus3"+label+"_"+self.channel, "gaus3"+label+"_"+self.channel, rrv_x, rrv_mean3_gaus, rrv_sigma1_gaus)
            gaus4 = RooGaussian("gaus4"+label+"_"+self.channel, "gaus4"+label+"_"+self.channel, rrv_x, rrv_mean4_gaus, rrv_sigma2_gaus)
            gausguas_2 = RooAddPdf("gausguas_2"+label+"_"+self.channel+mass_spectrum, "gausguas_2"+label+"_"+self.channel+mass_spectrum, RooArgList(gaus3, gaus4), RooArgList(rrv_frac1), 1)

            rrv_frac  = RooRealVar("rrv_frac"+label+"_"+self.channel, "rrv_frac"+label+"_"+self.channel, 0.74)#, 0.5, 1.0)
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, RooArgList(gausguas_1, gausguas_2), RooArgList(rrv_frac), 1)

        ## Erf*Exp + 2Gaus for mj spectrum
        if in_model_name == "2Gaus_ErfExp":

            print "########### 2Gaus + Erf*Exp for mj fit  ############"
            mean1_tmp      = 8.3141e+01
            mean1_tmp_err      = 1.63e-01
            deltamean_tmp  = 6.9129e+00
            deltamean_tmp_err  = 1.24e+00
            sigma1_tmp     = 7.5145e+00
            sigma1_tmp_err     = 1.99e-01
            scalesigma_tmp = 3.6819e+00
            scalesigma_tmp_err = 2.11e-01
            frac_tmp       = 6.7125e-01
            frac_tmp_err       = 2.09e-02

            ##rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel, "rrv_mean1_gaus"+label+"_"+self.channel, mean1_tmp, mean1_tmp-4, mean1_tmp+4)
            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel, "rrv_mean1_gaus"+label+"_"+self.channel, mean1_tmp, mean1_tmp-20, mean1_tmp+4)
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel, "rrv_sigma1_gaus"+label+"_"+self.channel, sigma1_tmp, sigma1_tmp-4, sigma1_tmp+4)
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel, "gaus1"+label+"_"+self.channel, rrv_x, rrv_mean1_gaus, rrv_sigma1_gaus)

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel, "rrv_deltamean_gaus"+label+"_"+self.channel, deltamean_tmp)#, deltamean_tmp, deltamean_tmp)
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel, "@0+@1", RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus))
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel, "rrv_scalesigma_gaus"+label+"_"+self.channel, scalesigma_tmp)#, scalesigma_tmp, scalesigma_tmp)
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel, "@0*@1", RooArgList(rrv_sigma1_gaus, rrv_scalesigma_gaus))
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel, "gaus2"+label+"_"+self.channel, rrv_x, rrv_mean2_gaus, rrv_sigma2_gaus)

            rrv_frac_2gaus = RooRealVar("rrv_frac_2gaus"+label+"_"+self.channel, "rrv_frac_2gaus"+label+"_"+self.channel, frac_tmp)#, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4)

            c0_tmp     = -2.9893e-02
            c0_tmp_err     = 6.83e-03
            offset_tmp = 7.9350e+01
            offset_tmp_err = 9.35e+00
            width_tmp  = 3.3083e+01
            width_tmp_err  = 2.97e+00

            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel, c0_tmp, c0_tmp-4e-2, c0_tmp+4e-2)
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel, "rrv_offset_ErfExp"+label+"_"+self.channel, offset_tmp)#, offset_tmp-offset_tmp_err*4, offset_tmp+offset_tmp_err*4)
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, width_tmp, width_tmp-10, width_tmp+10)
            erfexp = ROOT.RooErfExpPdf("erfexp"+label+"_"+self.channel+mass_spectrum, "erfexp"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c_ErfExp, rrv_offset_ErfExp, rrv_width_ErfExp)

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel, "rrv_frac"+label+"_"+self.channel, 0.5, 0., 1.)
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, RooArgList(erfexp, gaus1, gaus2), RooArgList(rrv_frac, rrv_frac_2gaus), 1)


        ## Erf*Exp+Voig+Gaus for mj spectrum
        if in_model_name == "ErfExpVoigGaus":
            print "########### Erf*Exp + Voig + Gaus for mj fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel, "rrv_c_ErfExp"+label+"_"+self.channel,-0.1,-10., 0.)
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel, "rrv_offset_ErfExp"+label+"_"+self.channel, 100., 10., 140.)
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel, "rrv_width_ErfExp"+label+"_"+self.channel, 30., 10, 100.)
            rrv_mean_voig     = RooRealVar("rrv_mean_voig"+label+"_"+self.channel, "rrv_mean_voig"+label+"_"+self.channel, 84, 78, 88)
            rrv_width_voig    = RooRealVar("rrv_width_voig"+label+"_"+self.channel, "rrv_width_voig"+label+"_"+self.channel, 7, 1, 20)
            rrv_sigma_voig    = RooRealVar("rrv_sigma_voig"+label+"_"+self.channel, "rrv_sigma_voig"+label+"_"+self.channel, 5, 1, 100)
            rrv_high1         = RooRealVar("rrv_high1"+label+"_"+self.channel, "rrv_high1"+label+"_"+self.channel, 1, 0., 200.)
            rrv_mean_gaus     = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel, "rrv_mean_gaus"+label+"_"+self.channel, 174)#, 160, 187)
            rrv_sigma_gaus    = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel, "rrv_sigma_gaus"+label+"_"+self.channel, 20)#, 0.1, 100)
            rrv_high2 = RooRealVar("rrv_high2"+label+"_"+self.channel, "rrv_high2"+label+"_"+self.channel, 0.)#, 0., 0.)
            model_pdf = ROOT.RooErfExp_Voig_Gaus_Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c_ErfExp, rrv_offset_ErfExp, rrv_width_ErfExp, rrv_mean_voig, rrv_width_voig, rrv_sigma_voig, rrv_high1, rrv_mean_gaus, rrv_sigma_gaus, rrv_high2)

        ## User1 function
        if in_model_name == "User1":
            print "########### User 1 Pdf  for mlvj fit ############"
            rrv_p0     = RooRealVar("rrv_p0_User1"+label+"_"+self.channel, "rrv_p0_User1"+label+"_"+self.channel, 30, 0, 60)
            if self.wtagger_category == "HP":
                rrv_p1 = RooRealVar("rrv_p1_User1"+label+"_"+self.channel, "rrv_p1_User1"+label+"_"+self.channel, -2, -9, 0)
            else:
                rrv_p1 = RooRealVar("rrv_p1_User1"+label+"_"+self.channel, "rrv_p1_User1"+label+"_"+self.channel, -2, -4, 0.)
            model_pdf = RooUser1Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_p0, rrv_p1)

        ## QCD pdf
        if in_model_name == "QCD":
            print "########### QCD Pdf  for mlvj fit ############"
            rrv_p0 = RooRealVar("rrv_p0_QCD"+label+"_"+self.channel, "rrv_p0_QCD"+label+"_"+self.channel, 0,-200, 200)
            rrv_p1 = RooRealVar("rrv_p1_QCD"+label+"_"+self.channel, "rrv_p1_QCD"+label+"_"+self.channel, 0,-200, 200)
            rrv_p2 = RooRealVar("rrv_p2_QCD"+label+"_"+self.channel, "rrv_p2_QCD"+label+"_"+self.channel, 0,-200, 200)
            model_pdf = RooQCDPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_p0, rrv_p1, rrv_p2)

        if in_model_name == "QCD_v2":#can replace exp
            print "########### QCD Pdf  for mlvj fit ############"
            rrv_p0 = RooRealVar("rrv_p0_QCD"+label+"_"+self.channel, "rrv_p0_QCD"+label+"_"+self.channel, -15,-50, 0)
            rrv_p1 = RooRealVar("rrv_p1_QCD"+label+"_"+self.channel, "rrv_p1_QCD"+label+"_"+self.channel, 20, 0, 250)
            rrv_p2 = RooRealVar("rrv_p2_QCD"+label+"_"+self.channel, "rrv_p2_QCD"+label+"_"+self.channel, 0,-20, 20)
            model_pdf = RooQCDPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_p0, rrv_p1, rrv_p2)

        ## For mlvj fit -> Pow function can replace exp
        if in_model_name == "Pow" or in_model_name == "Pow_sr" :
            print "########### Pow Pdf  for mlvj fit ############"
            rrv_c = RooRealVar("rrv_c_Pow"+label+"_"+self.channel, "rrv_c_Pow"+label+"_"+self.channel, -5, -20, 0)
            model_pdf = RooPowPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c)

        ## For mlvj fit -> Pow function can replace exp
        if in_model_name == "Pow2":
            print "########### Pow2 Pdf  for mlvj fit ############"
            rrv_c0 = RooRealVar("rrv_c0_Pow2"+label+"_"+self.channel, "rrv_c0_Pow2"+label+"_"+self.channel, 5, 0, 20)
            rrv_c1 = RooRealVar("rrv_c1_Pow2"+label+"_"+self.channel, "rrv_c1_Pow2"+label+"_"+self.channel, 0, -5 , 5)
            model_pdf = RooPow2Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c0, rrv_c1)

        ## For mlvj fit ->Erf*Pow can replace Erf*Exp
        if in_model_name == "ErfPow_v1":
            print "########### Erf*Pow Pdf  for mlvj fit ############"
            rrv_c      = RooRealVar("rrv_c_ErfPow"+label+"_"+self.channel, "rrv_c_ErfPow"+label+"_"+self.channel, -5,-10, 0)
            rrv_offset = RooRealVar("rrv_offset_ErfPow"+label+"_"+self.channel, "rrv_offset_ErfPow"+label+"_"+self.channel, 450, 350, 550)
            rrv_width  = RooRealVar("rrv_width_ErfPow"+label+"_"+self.channel, "rrv_width_ErfPow"+label+"_"+self.channel, 50, 20, 90)
            model_pdf  = RooErfPowPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c, rrv_offset, rrv_width)

        ## For mlvj fit ->Erf*Pow can replace Erf*Exp -> in the sideband
        if in_model_name == "ErfPow2_v1":
            print "########### Erf*Pow2 Pdf  for mlvj fit ############"
            rrv_c0 = RooRealVar("rrv_c0_ErfPow2"+label+"_"+self.channel, "rrv_c0_ErfPow2"+label+"_"+self.channel, 10, 1, 20)
            rrv_c1 = RooRealVar("rrv_c1_ErfPow2"+label+"_"+self.channel, "rrv_c1_ErfPow2"+label+"_"+self.channel, 5,-5, 10)
            rrv_offset = RooRealVar("rrv_offset_ErfPow2"+label+"_"+self.channel, "rrv_offset_ErfPow2"+label+"_"+self.channel, 450, 400, 520)
            rrv_width  = RooRealVar("rrv_width_ErfPow2"+label+"_"+self.channel, "rrv_width_ErfPow2"+label+"_"+self.channel, 30, 10, 80)
            model_pdf  = RooErfPow2Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c0, rrv_c1, rrv_offset, rrv_width)

        ## For mlvj fit ->Erf*Pow can replace Erf*Exp for sr
        if in_model_name == "ErfPow2_v1_sr":
            print "########### Erf*Pow2 Pdf  for mlvj fit in the SR  ############"
            rrv_c0 = RooRealVar("rrv_c0_ErfPow2"+label+"_"+self.channel, "rrv_c0_ErfPow2"+label+"_"+self.channel, 4, 2, 8)
            rrv_c1 = RooRealVar("rrv_c1_ErfPow2"+label+"_"+self.channel, "rrv_c1_ErfPow2"+label+"_"+self.channel, -0.5,-2, 0)
            rrv_offset = RooRealVar("rrv_offset_ErfPow2"+label+"_"+self.channel, "rrv_offset_ErfPow2"+label+"_"+self.channel, 490, 440, 520)
            rrv_width  = RooRealVar("rrv_width_ErfPow2"+label+"_"+self.channel, "rrv_width_ErfPow2"+label+"_"+self.channel, 50, 30, 80)
            model_pdf = RooErfPow2Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c0, rrv_c1, rrv_offset, rrv_width)

        ## For mlvj fit ->Erf*Pow*Exp can replace Erf*Exp
        if in_model_name == "ErfPowExp_v1":
            print "########### Erf*Pow*Exp Pdf  for mlvj fit   ############"
            if self.channel == "mu":
                rrv_c0 = RooRealVar("rrv_c0_ErfPowExp"+label+"_"+self.channel, "rrv_c0_ErfPowExp"+label+"_"+self.channel, 11, 5, 20)
                rrv_c1 = RooRealVar("rrv_c1_ErfPowExp"+label+"_"+self.channel, "rrv_c1_ErfPowExp"+label+"_"+self.channel, 0,-2, 2)
                rrv_offset = RooRealVar("rrv_offset_ErfPowExp"+label+"_"+self.channel, "rrv_offset_ErfPowExp"+label+"_"+self.channel, 470, 420, 520)
                rrv_width  = RooRealVar("rrv_width_ErfPowExp"+label+"_"+self.channel, "rrv_width_ErfPowExp"+label+"_"+self.channel, 40, 10, 80)
            else:
                rrv_c0 = RooRealVar("rrv_c0_ErfPowExp"+label+"_"+self.channel, "rrv_c0_ErfPowExp"+label+"_"+self.channel, 17, 0, 40)
                rrv_c1 = RooRealVar("rrv_c1_ErfPowExp"+label+"_"+self.channel, "rrv_c1_ErfPowExp"+label+"_"+self.channel, 1.7,-10, 10)
                rrv_offset = RooRealVar("rrv_offset_ErfPowExp"+label+"_"+self.channel, "rrv_offset_ErfPowExp"+label+"_"+self.channel, 520, 400, 650)
                rrv_width  = RooRealVar("rrv_width_ErfPowExp"+label+"_"+self.channel, "rrv_width_ErfPowExp"+label+"_"+self.channel, 40, 10, 80)
            model_pdf  = RooErfPowExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c0, rrv_c1, rrv_offset, rrv_width)

        ## For mlvj fit ->Erf*Pow*Exp can replace Erf*Exp
        if in_model_name == "ErfPowExp_v1_sr":
            print "########### Erf*Pow*Exp Pdf for mlvj fit in SR  ############"
            rrv_c0 = RooRealVar("rrv_c0_ErfPowExp"+label+"_"+self.channel, "rrv_c0_ErfPowExp"+label+"_"+self.channel, 6, 2, 15)
            rrv_c1 = RooRealVar("rrv_c1_ErfPowExp"+label+"_"+self.channel, "rrv_c1_ErfPowExp"+label+"_"+self.channel, -1,-3, 2)
            rrv_offset = RooRealVar("rrv_offset_ErfPowExp"+label+"_"+self.channel, "rrv_offset_ErfPowExp"+label+"_"+self.channel, 490, 440, 520)
            rrv_width  = RooRealVar("rrv_width_ErfPowExp"+label+"_"+self.channel, "rrv_width_ErfPowExp"+label+"_"+self.channel, 50, 30, 70)
            model_pdf = RooErfPowExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c0, rrv_c1, rrv_offset, rrv_width)

        ## For mlvj fit ->Erf*Pow*Exp can replace Erf*Exp
        if in_model_name == "ErfPowExp_v1_0":#difference inital value
            print "########### Erf*Pow*Exp Pdf for mlvj fit in SR  ############"
            rrv_c0 = RooRealVar("rrv_c0_ErfPowExp"+label+"_"+self.channel, "rrv_c0_ErfPowExp"+label+"_"+self.channel, 20, 15, 40)
            rrv_c1 = RooRealVar("rrv_c1_ErfPowExp"+label+"_"+self.channel, "rrv_c1_ErfPowExp"+label+"_"+self.channel, 1.6, 0.5, 5)
            rrv_offset = RooRealVar("rrv_offset_ErfPowExp"+label+"_"+self.channel, "rrv_offset_ErfPowExp"+label+"_"+self.channel, 470, 420, 520)
            rrv_width  = RooRealVar("rrv_width_ErfPowExp"+label+"_"+self.channel, "rrv_width_ErfPowExp"+label+"_"+self.channel, 47, 30, 60)
            model_pdf  = RooErfPowExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c0, rrv_c1, rrv_offset, rrv_width)

        ## Keys
        if in_model_name == "Keys":
            print "########### Erf*Pow*Exp Pdf for Keys  ############"
            rdataset = self.workspace4fit_.data("rdataset_%s_signal_region_mlvj"%(self.signal_sample))
            model_pdf = RooKeysPdf("model_pdf"+label+"_"+self.channel+mass_spectrum, "model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rdataset)

        ## return the pdf
        getattr(self.workspace4fit_, "import")(model_pdf)
        return self.workspace4fit_.pdf("model_pdf"+label+"_"+self.channel+mass_spectrum)

    ########### Gaussian contraint of a parameter of a pdf
    def addConstraint(self, rrv_x, x_mean, x_sigma, ConstraintsList):
        print "########### Add to Contraint List some parameters  ############"
        rrv_x_mean = RooRealVar(rrv_x.GetName()+"_mean", rrv_x.GetName()+"_mean", x_mean)
        rrv_x_sigma = RooRealVar(rrv_x.GetName()+"_sigma", rrv_x.GetName()+"_sigma", x_sigma)
        constrainpdf_x = RooGaussian("constrainpdf_"+rrv_x.GetName(), "constrainpdf_"+rrv_x.GetName(), rrv_x, rrv_x_mean, rrv_x_sigma)
        ## import in the workspace and save the name of constriant pdf
        getattr(self.workspace4fit_, "import")(constrainpdf_x)
        ConstraintsList.append(constrainpdf_x.GetName())

    ### take the dataset, the model , the parameters in order to fix them as constant --> for extended pdf
    def get_General_mj_Model(self, label):
        print "########### Fixing a general mj model  ############"
        rdataset_General_mj = self.workspace4fit_.data("rdataset%s_%s_mj"%(label, self.channel))
        model_General =  self.workspace4fit_.pdf("model"+label+"_"+self.channel+"_mj")
        rdataset_General_mj.Print()
        model_General.Print()
        ## get the parameters and cycle on them
        parameters_General = model_General.getParameters(rdataset_General_mj)
        par = parameters_General.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            if (TString(label).Contains("VV") or TString(label).Contains("STop") or TString(label).Contains("TTbar")):
                param.Print()
            param.setConstant(kTRUE)
            param = par.Next()
        ## return the pdf after having fixed the paramters
        return self.workspace4fit_.pdf("model%s_%s_mj"%(label, self.channel))

    ### fix only the ttbar component using the default label --> for extended pdf
    def get_TTbar_mj_Model(self, label = "_TTbar_xww"):
        print "########### Fixing only the TTbar mj Shape  ############"
        return self.get_General_mj_Model(label)

    ### fix only the stop component using the default label --> for extended pdf
    def get_STop_mj_Model(self, label = "_STop_xww"):
        print "########### Fixing only the Stop mj Shape  ############"
        return self.get_General_mj_Model(label)

    ### fix only the VV component using the default label --> for extended pdf
    def get_VV_mj_Model(self, label = "_VV_xww"):
        print "########### Fixing only the VV mj Shape  ############"
        return self.get_General_mj_Model(label)

    ### fix only the WJets model --> for extended pdf (just fix shape parameters of width, offset of ErfExp and p1 of User1 function
    def get_WJets_mj_Model(self, label):
        print "########### Fixing only the WJets mj Shape --> just the printed parameters  ############"
        rdataset_WJets_mj = self.workspace4fit_.data("rdataset%s_%s_mj"%(label, self.channel))
        model_WJets =  self.workspace4fit_.pdf("model"+label+"_"+self.channel+"_mj")
        rdataset_WJets_mj.Print()
        model_WJets.Print()
        parameters_WJets = model_WJets.getParameters(rdataset_WJets_mj)
        par = parameters_WJets.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            paraName = TString(param.GetName())
            if (paraName.Contains("rrv_width_ErfExp_WJets") or paraName.Contains("rrv_offset_ErfExp_WJets") or paraName.Contains("rrv_p1_User1_WJets")):
                param.setConstant(kTRUE)
                param.Print()
            else:
                param.setConstant(0)
            param = par.Next()
        return self.workspace4fit_.pdf("model%s_%s_mj"%(label, self.channel))

    ### fix a given model taking the label, and the region --> for extended pdf --> all the parameter of the pdf + normalization
    def fix_Model(self, label, mlvj_region = "_signal_region", mass_spectrum = "_mlvj"):
        print "########### Fixing an Extended Pdf for mlvj  ############"
        rdataset = self.workspace4fit_.data("rdataset%s%s_%s%s"%(label, mlvj_region, self.channel, mass_spectrum))
        model = self.get_mlvj_Model(label, mlvj_region)
        rdataset.Print()
        model.Print()
        parameters = model.getParameters(rdataset)
        par = parameters.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            param.setConstant(kTRUE)
            param = par.Next()



    #### get a generic mlvj model from the workspace
    def get_mlvj_Model(self, label, mlvj_region):
        return self.workspace4fit_.pdf("model"+label+mlvj_region+"_"+self.channel+"_mlvj")

    #### get a general mlvj model and fiz the paramters --> for extended pdf
    def get_General_mlvj_Model(self, label, mlvj_region = "_signal_region"):
        print "########### Fixing a general mlvj model  ############"
        rdataset_General_mlvj = self.workspace4fit_.data("rdataset%s%s_%s_mlvj"%(label, mlvj_region, self.channel))
        model_General = self.get_mlvj_Model(label, mlvj_region)
        rdataset_General_mlvj.Print()
        model_General.Print()
        parameters_General = model_General.getParameters(rdataset_General_mlvj)
        par = parameters_General.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            param.setConstant(kTRUE)
            param.Print()
            param = par.Next()
        return self.get_mlvj_Model(label, mlvj_region)

    ###### get TTbar model mlvj in a region
    def get_TTbar_mlvj_Model(self, mlvj_region = "_signal_region"):
        print "########### Fixing TTbar mlvj model for the region", mlvj_region, "  ############"
        return self.get_General_mlvj_Model("_TTbar_xww", mlvj_region)

    ###### get Single Top model mlvj in a region
    def get_STop_mlvj_Model(self, mlvj_region = "_signal_region"):
        print "########### Fixing Stop mlvj model for the region", mlvj_region, "  ############"
        return self.get_General_mlvj_Model("_STop_xww", mlvj_region)

    ###### get Signal model mlvj in a region
    def get_signal_mlvj_Model(self, mlvj_region = "_signal_region"):
        print "########### Fixing signal mlvj model for the region", mlvj_region, "  ############"
        return self.get_General_mlvj_Model("_%s_xww"%(self.signal_sample), mlvj_region)

    ###### get VV mlvj in a region
    def get_VV_mlvj_Model(self, mlvj_region = "_signal_region"):
        print "########### Fixing VV mlvj for the region", mlvj_region, "  ############"
        return self.get_General_mlvj_Model("_VV_xww", mlvj_region)

    ###### get W+jets mlvj in a region
    def get_WJets_mlvj_Model(self, mlvj_region = "_signal_region"):
        rdataset_WJets_mlvj = self.workspace4fit_.data("rdataset_WJets_xww_%s_mlvj"%(mlvj_region))
        model_WJets = self.get_mlvj_Model("_WJets0_xww", mlvj_region)
        print "######## get Wjet mlvj model for the region --> set constant just the normalization from mj fit", mlvj_region, " ########"
        rdataset_WJets_mlvj.Print()
        model_WJets.Print()
        parameters_WJets = model_WJets.getParameters(rdataset_WJets_mlvj)
        par = parameters_WJets.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            paraName = TString(param.GetName())
            param.Print()
            if paraName.Contains("rrv_number_WJets"): ## set the correct normalization for W+jets if we are inside the signal region and fix it as constant
                if self.workspace4fit_.var("rrv_number_WJets_xww_in_mj%s_from_fitting_%s"%(mlvj_region, self.channel)):
                    self.workspace4fit_.var("rrv_number_WJets_xww_in_mj%s_from_fitting_%s"%(mlvj_region, self.channel)).Print()
                    param.setVal(self.workspace4fit_.var("rrv_number_WJets_xww_in_mj%s_from_fitting_%s"%(mlvj_region, self.channel)).getVal())
                if mlvj_region == "_signal_region":
                    param.setConstant(kTRUE)
            param.Print()
            param = par.Next()
        return self.get_mlvj_Model("_WJets0_xww", mlvj_region)


    ### change a dataset to a histpdf roofit object
    def change_dataset_to_histpdf(self, x, dataset):
        print "######## change the dataset into a histpdf  ########"
        datahist = dataset.binnedClone(dataset.GetName()+"_binnedClone", dataset.GetName()+"_binnedClone")
        histpdf = RooHistPdf(dataset.GetName()+"_histpdf", dataset.GetName()+"_histpdf", RooArgSet(x), datahist)
        dataset.Print()
        histpdf.Print()
        getattr(self.workspace4fit_, "import")(histpdf)

    ### change from a dataset to a histogramm of Roofit
    def change_dataset_to_histogram(self, x, dataset, label = ""):
        print "######## change the dataset into a histogramm for mj distribution ########"
        datahist = dataset.binnedClone(dataset.GetName()+"_binnedClone", dataset.GetName()+"_binnedClone")
        nbin = int((x.getMax()-x.getMin())/self.BinWidth_mj)
        if label == "":
            return datahist.createHistogram("histo_%s"%(dataset.GetName()), x, RooFit.Binning(nbin , x.getMin(), x.getMax()))
        else:
            return datahist.createHistogram("histo_"+label, x, RooFit.Binning(nbin, x.getMin(), x.getMax()))


    ### Define the Extended Pdf for and mJ fit giving: label, fit model name, list constraint and ismc
    def make_Model(self, label, in_model_name, mass_spectrum = "_mj", ConstraintsList=[], ismc_wjet=0, area_init_value=500):
        ##### define an extended pdf from a standard Roofit One
        print " "
        print "###############################################"
        print "## Make model : ", label, " ", in_model_name, "##"
        print "###############################################"
        print " "

        rrv_number = RooRealVar("rrv_number"+label+"_"+self.channel+mass_spectrum, "rrv_number"+label+"_"+self.channel+mass_spectrum, area_init_value, 0., 1e7)
        ## call the make RooAbsPdf method
        model_pdf = self.make_Pdf(label, in_model_name, mass_spectrum,ConstraintsList, ismc_wjet)
        print "######## Model Pdf ########"
        model_pdf.Print()

        ## create the extended pdf
        model = RooExtendPdf("model"+label+"_"+self.channel+mass_spectrum, "model"+label+"_"+self.channel+mass_spectrum, model_pdf, rrv_number)
        print "######## Model Extended Pdf ########"

        #### put all the parameters and the shape in the workspace
        getattr(self.workspace4fit_, "import")(rrv_number)
        getattr(self.workspace4fit_, "import")(model)
        self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum).Print()
        ## return the total extended pdf
        return self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum)

    ### Method for a single MC fit of the mj spectra giving: file name, label, model name
    def fit_mj_single_MC(self, in_file_name, label, in_model_name, additioninformation = ""):

        print "############### Fit mj single MC sample", in_file_name, " ", label, "  ", in_model_name, " ##################"
        ## import variable and dataset
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
        rdataset_mj = self.workspace4fit_.data("rdataset4fit"+label+"_"+self.channel+"_mj")
        rdataset_mj.Print()

        ## make the extended model
        model = self.make_Model(label, in_model_name)
        rfresult = model.fitTo(rdataset_mj, RooFit.Save(1), RooFit.SumW2Error(kTRUE) , RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit"))
        rfresult = model.fitTo(rdataset_mj, RooFit.Save(1), RooFit.SumW2Error(kTRUE) , RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"))
        rfresult.Print()

        if options.MCStudy:
            MCStudy(model, rrv_mass_j, "%s/other/%s/"%(self.plotsDir, self.signal_sample)+"try_"+model.GetName()+".png")
        ## Plot the result
        mplot = rrv_mass_j.frame(RooFit.Title(label+" fitted by "+in_model_name), RooFit.Bins(int(rrv_mass_j.getBins()/self.binwidth_narrow_factor)))
        rdataset_mj.plotOn(mplot, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0))

        ## draw the error band for an extend pdf
        draw_error_band_extendPdf(rdataset_mj, model, rfresult, mplot, 2, "L")
        ## re-draw the dataset
        rdataset_mj.plotOn(mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0))
        ## draw the function
        model.plotOn(mplot)# remove RooFit.VLines() in order to get right pull in the 1st bin

        ## Get the pull
        mplot_pull = self.get_pull(rrv_mass_j, mplot)
        mplot.GetYaxis().SetRangeUser(1e-2, mplot.GetMaximum()*1.5)

        nPar = rfresult.floatParsFinal().getSize()
        nBinX = mplot.GetNbinsX()
        ndof  = nBinX-nPar
        print "#################### nPar = %s, nBinX = %s , chiSquare = %s/%s"%(nPar, nBinX, mplot.chiSquare(nPar)*ndof, ndof)

        datahist = rdataset_mj.binnedClone(rdataset_mj.GetName()+"_binnedClone", rdataset_mj.GetName()+"_binnedClone")
        parameters_list = model.getParameters(rdataset_mj)
        self.draw_canvas_with_pull1(rrv_mass_j, datahist, mplot, mplot_pull, ndof, parameters_list, "%s/m_j_fitting%s/"%(self.plotsDir, additioninformation), label+in_file_name, in_model_name)

        #normalize the number of total events to lumi --> correct the number to scale to the lumi
        self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal())
        self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal())

        ##if TString(label).Contains("ggH"):
        ##    self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal())
        ##    self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getError())
        self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").Print()

        ##### apply the correction of the mean and sigma from the ttbar control sample to the STop, TTbar and VV
        par = parameters_list.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            if (TString(label).Contains("VV") or TString(label).Contains("STop") or TString(label).Contains("TTbar")):
                #param.Print()
                if TString(param.GetName()).Contains("rrv_mean1_gaus"):
                    param.setRange(param.getMin()+self.mean_shift, param.getMax()+self.mean_shift)
                    param.setVal(param.getVal()+self.mean_shift)
                if TString(param.GetName()).Contains("rrv_deltamean_gaus"):
                    param.setRange(param.getMin()-self.mean_shift, param.getMax()-self.mean_shift)
                    param.setVal(param.getVal()-self.mean_shift)
                if TString(param.GetName()).Contains("rrv_sigma1_gaus"):
                    param.setVal(param.getVal()*self.sigma_scale)
                    param.setRange(param.getMin()*self.sigma_scale, param.getMax()*self.sigma_scale)
                if TString(param.GetName()).Contains("rrv_scalesigma_gaus"):
                    param.setRange(param.getMin()/self.sigma_scale, param.getMax()/self.sigma_scale)
                    param.setVal(param.getVal()/self.sigma_scale)
            param = par.Next()

    ######## ++++++++++++++
    def ControlPlots(self):
        ##event selection
        MET_cut = 40
        lpt_cut   = 50
        if self.channel == "el":
            MET_cut = 80
            lpt_cut = 55

        rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j")
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")

        #cut= "(CategoryID == 1 || CategoryID == -1 || CategoryID == 2 || CategoryID == -2|| CategoryID == 4 || CategoryID == -4) && m_lvj> 100 && m_lvj<3000 &&((massVhadJEC>40 && massVhadJEC<65)||(massVhadJEC>135&&massVhadJEC<150)) && l_pt> %s && MET_et>%s"%(lpt_cut, MET_cut)
        #self.Make_Controlplots(cut, "preselection")

        #cut= "(CategoryID == 3 || CategoryID == -3) && m_lvj> 100 && m_lvj<4000 && massVhadJEC>40 && massVhadJEC<150 && l_pt>%s && MET_et>%s"%(lpt_cut, MET_cut)
        #self.Make_Controlplots(cut, "TopControl", 1)

        cut= "(CategoryID == 1) && m_lvj> 600 && m_lvj<14000 && (massVhadJEC>65 && massVhadJEC<95) && l_pt> %s && MET_et>%s"%(lpt_cut, MET_cut)
        self.Make_Controlplots(cut, "fulselection")


    ######## ++++++++++++++
    def Make_Controlplots(self, cut, tag, TTBarControl = 0):
        self.make_controlplot("m_lvj", cut, tag, 26, 200, 1500, "mass(lvj)", "Events/(50 GeV)", 0, TTBarControl)
        self.make_controlplot("massVhadJEC", cut, tag, 23, 40, 155, "mass(j)", "Events/(5 GeV)", 0 , TTBarControl)
        self.make_controlplot("W_pt", cut, tag, 30, 200, 800, "W_pt", "Events/(20 GeV)", 0 , TTBarControl)
        self.make_controlplot("l_pt", cut, tag, 26, 0, 520, "l_pt", "Events/(20 GeV)", 0 , TTBarControl)
        self.make_controlplot("l_eta", cut, tag, 20,-2.5, 2.5, "l_eta", "Events(0.25)", 0 , TTBarControl)
        self.make_controlplot("MET_et", cut, tag, 30, 0, 600, "MET_et", "Events/(20 GeV)", 0 , TTBarControl)
        self.make_controlplot("nPV", cut, tag, 20, 0, 40, "nPV", "Events/(2)", 0 , TTBarControl)
        self.make_controlplot("tau21", cut, tag, 20, 0, 1, "tau21", "Events/(0.05)", 0 , TTBarControl)
        self.make_controlplot("nbtag", cut, tag, 5,-0.5, 4.5, "number of b-jets", "Events", 0 , TTBarControl)

    ######## ++++++++++++++
    def make_controlplot(self, variable, cut, tag, nbin, x_min, x_max, xtitle = "", ytitle = "", logy=0 , TTBarControl=0):
        tmp_lumi = self.GetLumi()
        tmp_signal_scale = 20
        weight_mc_forSignal = "weight*%s*%s"%(tmp_lumi, tmp_signal_scale)
        weight_mc_forV = "weight*%s*%s"%(tmp_lumi, self.rrv_wtagger_eff_reweight_forV.getVal())
        weight_mc_forT = "weight*%s*%s"%(tmp_lumi, self.rrv_wtagger_eff_reweight_forT.getVal())
        weight_mc_forG = "weight*%s"%(tmp_lumi) #General
        weight_mc_forWJets = "weight*%s*%s"%(tmp_lumi, self.controlplot_WJets_scale)
        weight_mc_forTTBar = "weight*%s*%s*%s"%(tmp_lumi, self.controlplot_TTbar_scale, self.rrv_wtagger_eff_reweight_forT.getVal())

        weightcut_mc_forSignal = "(%s)*(%s)"%(weight_mc_forSignal, cut)
        weightcut_mc_forV = "(%s)*(%s)"%(weight_mc_forV, cut)
        weightcut_mc_forT = "(%s)*(%s)"%(weight_mc_forT, cut)
        weightcut_mc_forG = "(%s)*(%s)"%(weight_mc_forG, cut)
        weightcut_mc_forWJets = "(%s)*(%s)"%(weight_mc_forWJets, cut)
        weightcut_mc_forTTBar = "(%s)*(%s)"%(weight_mc_forTTBar, cut)
        weightcut_data = "%s"%(cut)
        print "weightcut_mc_forV = "+weightcut_mc_forV
        print "weightcut_mc_forT = "+weightcut_mc_forT
        print "weightcut_mc_forG = "+weightcut_mc_forG
        print "weightcut_mc_forWJets = "+weightcut_mc_forWJets
        print "weightcut_mc_forTTBar = "+weightcut_mc_forTTBar
        hist_data = TH1D("hist_data", "hist_data"+";%s;%s"%(xtitle, ytitle), nbin, x_min, x_max)
        hist_Signal = TH1D("hist_Signal", "hist_Signal"+";%s;%s"%(xtitle, ytitle), nbin, x_min, x_max)
        hist_Signal.Sumw2()
        hist_WJets = TH1D("hist_WJets", "hist_WJets"+";%s;%s"%(xtitle, ytitle), nbin, x_min, x_max)
        hist_WJets.Sumw2()
        hist_TTbar = TH1D("hist_TTbar", "hist_TTbar"+";%s;%s"%(xtitle, ytitle), nbin, x_min, x_max)
        hist_TTbar.Sumw2()
        hist_STop = TH1D("hist_STop", "hist_STop"+";%s;%s"%(xtitle, ytitle), nbin, x_min, x_max)
        hist_STop.Sumw2()
        hist_VV   = TH1D("hist_VV", "hist_VV"+";%s;%s"%(xtitle, ytitle), nbin, x_min, x_max)
        hist_VV.Sumw2()
        hist_TotalMC = TH1D("hist_TotalMC", "hist_TotalMC"+";%s;%s"%(xtitle, ytitle), nbin, x_min, x_max)
        hist_TotalMC.Sumw2()


        hstack_TotalMC = THStack("hstack_TotalMC", "hstack_TotalMC"+";%s;%s"%(xtitle, ytitle))
        if TTBarControl == 0:
            hstack_TotalMC.Add(hist_STop)
            hstack_TotalMC.Add(hist_TTbar)
            hstack_TotalMC.Add(hist_VV)
            hstack_TotalMC.Add(hist_WJets)
        else:
            hstack_TotalMC.Add(hist_WJets)
            hstack_TotalMC.Add(hist_VV)
            hstack_TotalMC.Add(hist_STop)
            hstack_TotalMC.Add(hist_TTbar)

        hist_data.SetLineColor(self.color_palet["data"])
        hist_data.SetFillColor(self.color_palet["data"])
        hist_Signal.SetLineColor(self.color_palet["Signal"])
        hist_Signal.SetFillColor(self.color_palet["Signal"])
        hist_Signal.SetFillStyle(0)
        hist_Signal.SetLineWidth(2)
        hist_WJets.SetLineColor(kBlack)
        hist_WJets.SetFillColor(self.color_palet["WJets"])
        hist_TTbar.SetLineColor(kBlack)
        hist_TTbar.SetFillColor(self.color_palet["TTbar"])
        hist_STop.SetLineColor(kBlack)
        hist_STop.SetFillColor(self.color_palet["STop"])
        hist_VV.SetLineColor(kBlack)
        hist_VV.SetFillColor(self.color_palet["VV"])

        tree_data   = TChain("PKUTree")
        tree_data.Add(self.file_Directory+self.file_data)
        tree_Signal = TChain("PKUTree")
        tree_Signal.Add(self.file_Directory+self.file_signal)
        tree_WJets  = TChain("PKUTree")
        tree_WJets.Add(self.file_Directory+self.file_WJets0_mc)
        tree_TTbar  = TChain("PKUTree")
        tree_TTbar.Add(self.file_Directory+self.file_TTbar_mc)
        tree_STop   = TChain("PKUTree")
        tree_STop.Add(self.file_Directory+self.file_STop_mc)
        tree_VV     = TChain("PKUTree")
        tree_VV.Add(self.file_Directory+self.file_VV_mc)

        tree_data.Draw("%s >> hist_data"%(variable), weightcut_data)
        tree_Signal.Draw("%s >> hist_Signal"%(variable), weightcut_mc_forSignal)
        #tree_WJets.Draw("%s >> hist_WJets"%(variable), weightcut_mc_forG)
        tree_WJets.Draw("%s >> hist_WJets"%(variable), weightcut_mc_forWJets)
        #tree_TTbar.Draw("%s >> hist_TTbar"%(variable), weightcut_mc_forT)
        tree_TTbar.Draw("%s >> hist_TTbar"%(variable), weightcut_mc_forTTBar)
        tree_STop.Draw("%s >> hist_STop"%(variable), weightcut_mc_forT)
        tree_VV.Draw("%s >> hist_VV"%(variable), weightcut_mc_forV)

        hist_TotalMC.Add(hist_WJets)
        hist_TotalMC.Add(hist_TTbar)
        hist_TotalMC.Add(hist_STop)
        hist_TotalMC.Add(hist_VV)

        #canvas_controlplot = TCanvas("canvas_controlplot"+variable, "canvas_controlplot"+variable, 600, 600)
        canvas_controlplot = self.get_canvas("canvas_controlplot"+variable)

        canvas_controlplot.cd()
        hist_data.GetYaxis().SetRangeUser(1e-2, TMath.Max(hist_data.GetMaximum(), hist_TotalMC.GetMaximum())*1.8)
        hist_data.Draw("e")
        hstack_TotalMC.Draw("HIST same")
        hist_data.Draw("same e")
        hist_Signal.Draw("same HIST")

        #hist_TotalMC.Draw("same e")
        gr_MCStat = TGraph(hist_TotalMC.GetNbinsX()*4)
        gr_MCStat.SetName("MCStat")
        gr_MCStat.SetTitle("MC Stat")
        gr_MCStat.SetFillColor(1)
        gr_MCStat.SetFillStyle(3013)
        count = 0
        while (count<hist_TotalMC.GetNbinsX()):
            gr_MCStat.SetPoint(count*2  , hist_TotalMC.GetBinCenter(count+1)-0.5*hist_TotalMC.GetBinWidth(count+1), hist_TotalMC.GetBinContent(count+1)+hist_TotalMC.GetBinError(count+1))
            gr_MCStat.SetPoint(count*2+1, hist_TotalMC.GetBinCenter(count+1)+0.5*hist_TotalMC.GetBinWidth(count+1), hist_TotalMC.GetBinContent(count+1)+hist_TotalMC.GetBinError(count+1))
            count = count+1
        count_re = 0
        while (count>0):
            gr_MCStat.SetPoint(hist_TotalMC.GetNbinsX()*2+count_re*2  , hist_TotalMC.GetBinCenter(count)+0.5*hist_TotalMC.GetBinWidth(count), hist_TotalMC.GetBinContent(count)-hist_TotalMC.GetBinError(count))
            gr_MCStat.SetPoint(hist_TotalMC.GetNbinsX()*2+count_re*2+1, hist_TotalMC.GetBinCenter(count)-0.5*hist_TotalMC.GetBinWidth(count), hist_TotalMC.GetBinContent(count)-hist_TotalMC.GetBinError(count))
            count_re = count_re+1
            count = count-1
        print count, count_re

        gr_MCStat.Draw("F")

        hist_data.GetXaxis().SetTitleOffset(1.1)
        hist_data.GetYaxis().SetTitleOffset(1.3)
        hist_data.GetXaxis().SetTitleSize(0.03)
        hist_data.GetYaxis().SetTitleSize(0.03)
        hist_data.GetXaxis().SetLabelSize(0.03)
        hist_data.GetYaxis().SetLabelSize(0.03)

        banner = self.banner4Plot()
        banner.Draw()

        theLeg = TLegend(0.41, 0.61, 0.76, 0.93, "", "NDC")
        theLeg.SetName("theLegend")
        theLeg.SetBorderSize(0)
        theLeg.SetLineColor(0)
        theLeg.SetFillColor(0)
        theLeg.SetFillStyle(0)
        theLeg.SetLineWidth(0)
        theLeg.SetLineStyle(0)
        theLeg.SetTextFont(42)
        theLeg.SetTextSize(.04)
        theLeg.SetNColumns(2)

        theLeg.SetFillColor(0)
        theLeg.SetFillStyle(0)
        theLeg.SetBorderSize(0)
        theLeg.SetLineColor(0)
        theLeg.SetLineWidth(0)
        theLeg.SetLineStyle(0)
        theLeg.SetTextSize(0.040)
        theLeg.SetTextFont(42)

        theLeg.AddEntry(hist_data, "CMS Data", "ep")
        theLeg.AddEntry(hist_WJets, "WJets", "F")
        theLeg.AddEntry(hist_VV, "VV", "F")
        theLeg.AddEntry(hist_TTbar, "TTbar", "F")
        theLeg.AddEntry(hist_STop, "Single Top", "F")
        theLeg.AddEntry(gr_MCStat, "MC Stat", "F")
        theLeg.AddEntry(hist_Signal, self.signal_sample+" #times %s"%(tmp_signal_scale), "L")
        theLeg.SetY1NDC(0.9 - 0.05*6 - 0.005)
        theLeg.SetY1(theLeg.GetY1NDC())
        theLeg.Draw()

        #Directory = TString("plots_%s_%s_%s/controlplot_wtaggercut%s/"%(options.additioninformation, self.channel, self.wtagger_category, self.wtagger_category)+self.signal_sample+"/")
        Directory = TString(self.plotsDir+"/controlplot/"+self.signal_sample+"/")

        if not Directory.EndsWith("/"):
            Directory = Directory.Append("/")
        if not os.path.isdir(Directory.Data()):
            os.system("mkdir -p  "+Directory.Data())

        #only draw png
        rlt_file = TString(Directory.Data()+"controlplot_"+variable+"_"+tag+".png")
        canvas_controlplot.SaveAs(rlt_file.Data())
        if logy:
            canvas_controlplot.SetLogy()
            canvas_controlplot.Update()
            rlt_file.ReplaceAll(".png", "_log.png")
            canvas_controlplot.SaveAs(rlt_file.Data())

        #draw png and pdf
        #rlt_file = TString(Directory.Data()+"controlplot_"+variable+"_"+tag+".png")
        #canvas_controlplot.SaveAs(rlt_file.Data())
        #rlt_file.ReplaceAll(".png", ".pdf")
        #canvas_controlplot.SaveAs(rlt_file.Data())

        #if logy:
        #    canvas_controlplot.SetLogy()
        #    canvas_controlplot.Update()
        #    rlt_file.ReplaceAll(".pdf", "_log.pdf")
        #    canvas_controlplot.SaveAs(rlt_file.Data())
        #    rlt_file.ReplaceAll(".pdf", ".png")
        #    canvas_controlplot.SaveAs(rlt_file.Data())

    ### Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
    def fit_mlvj_model_single_MC(self, in_file_name, label, in_range, mlvj_model, deco = 0, show_constant_parameter = 0, logy = 1, ismc = 0):

        print "############### Fit mlvj single MC sample ", in_file_name, " ", label, "  ", mlvj_model, "  ", in_range, " ##################"
        ## import variable and dataset
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rdataset = self.workspace4fit_.data("rdataset4fit"+label+in_range+"_"+self.channel+"_mlvj")
        constrainslist = []

        ## make the extended pdf model
        model = self.make_Model(label+in_range, mlvj_model, "_mlvj", constrainslist, ismc)

        ## make the fit
        rfresult = model.fitTo(rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) , RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit"))
        rfresult = model.fitTo(rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) , RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"))
        rfresult.Print()

        if options.MCStudy:
            MCStudy(model, rrv_mass_lvj, "%s/other/%s/"%(self.plotsDir, self.signal_sample)+"try_"+model.GetName()+".png")
        ## set the name of the result of the fit and put it in the workspace
        rfresult.SetName("rfresult"+label+in_range+"_"+self.channel+"_mlvj")
        getattr(self.workspace4fit_, "import")(rfresult)

        ## plot the result
        mplot = rrv_mass_lvj.frame(RooFit.Title("M_{lvj"+in_range+"} fitted by "+mlvj_model), RooFit.Bins(int(rrv_mass_lvj.getBins()/self.binwidth_narrow_factor)))
        rdataset.plotOn(mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0))
        ## plot the error band but don't store the canvas (only plotted without -b option
        draw_error_band_extendPdf(rdataset, model, rfresult, mplot, 2, "L")
        rdataset.plotOn(mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0))
        model.plotOn(mplot)#, RooFit.VLines()) in order to have the right pull

        nPar = rfresult.floatParsFinal().getSize()
        nBinX = mplot.GetNbinsX()
        ndof  = nBinX-nPar
        print mplot.chiSquare()
        print "#################### JENchi2 nPar = %s, chiSquare = %s/%s"%(nPar , mplot.chiSquare(nPar)*ndof, ndof)
        datahist = rdataset.binnedClone(rdataset.GetName()+"_binnedClone", rdataset.GetName()+"_binnedClone")
        rdataset.Print()
        ## get the pull
        mplot_pull      = self.get_pull(rrv_mass_lvj, mplot)
        parameters_list = model.getParameters(rdataset)
        mplot.GetYaxis().SetRangeUser(1e-2, mplot.GetMaximum()*1.5)

        self.draw_canvas_with_pull1(rrv_mass_lvj, datahist, mplot, mplot_pull, ndof, parameters_list, "%s/m_lvj_fitting/"%(self.plotsDir), in_file_name, "m_lvj"+in_range+mlvj_model, show_constant_parameter, logy)


        ## if the shape parameters has to be decorrelated
        if deco :
            print "################### Decorrelated mlvj single mc shape ################"
            model_pdf = self.workspace4fit_.pdf("model_pdf%s%s_%s_mlvj"%(label, in_range, self.channel)) ## take the pdf from the workspace
            rfresult_pdf = model_pdf.fitTo(rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit"))
            rfresult_pdf = model_pdf.fitTo(rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"))
            rfresult_pdf.Print()
            if options.MCStudy:
                MCStudy(model, rrv_mass_lvj, "%s/other/%s/"%(self.plotsDir, self.signal_sample)+"try_"+model.GetName()+".png")

            ## temp workspace for the pdf diagonalizer
            wsfit_tmp = RooWorkspace("wsfit_tmp"+label+in_range+"_"+self.channel+"_mlvj")
            Deco      = PdfDiagonalizer("Deco"+label+in_range+"_"+self.channel+"_"+self.wtagger_category+"_mlvj", wsfit_tmp, rfresult_pdf) ## in order to have a good name
            print "##################### diagonalize "
            model_pdf_deco = Deco.diagonalize(model_pdf) ## diagonalize
            print "##################### workspace for decorrelation "
            wsfit_tmp.Print("v")
            print "##################### original  parameters "
            model_pdf.getParameters(rdataset).Print("v")
            print "##################### original  decorrelated parameters "
            model_pdf_deco.getParameters(rdataset).Print("v")
            print "##################### original  pdf "
            model_pdf.Print()
            print "##################### decorrelated pdf "
            model_pdf_deco.Print()

            ## import in the workspace and print the diagonalizerd pdf
            getattr(self.workspace4fit_, "import")(model_pdf_deco)

            ### define a frame for TTbar or other plots
            mplot_deco = rrv_mass_lvj.frame(RooFit.Bins(int(rrv_mass_lvj.getBins()/self.binwidth_narrow_factor)))

            if label == "_TTbar_xww" and in_range == "_signal_region":

                rdataset.plotOn(mplot_deco, RooFit.Name("Powheg Sample"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0))
                model_pdf_deco.plotOn(mplot_deco, RooFit.Name("TTbar_Powheg"), RooFit.LineColor(kBlack))

                mplot_deco.GetYaxis().SetRangeUser(1e-2, mplot_deco.GetMaximum()*1.5)

                rrv_number_dataset = RooRealVar("rrv_number_dataset", "rrv_number_dataset", rdataset.sumEntries())
                rrv_number_dataset.setError(0.)
                draw_error_band(rdataset, model_pdf, rrv_number_dataset, rfresult_pdf, mplot_deco, self.color_palet["Uncertainty"], "F") ## draw the error band with the area
                self.workspace4fit_.var("rrv_number_TTbar_xww_signal_region_%s_mlvj"%(self.channel)).Print()
            else:
                rdataset.plotOn(mplot_deco, RooFit.Name("Data"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0))
                model_pdf_deco.plotOn(mplot_deco, RooFit.Name(label), RooFit.LineColor(kBlack))

                mplot_deco.GetYaxis().SetRangeUser(1e-2, mplot_deco.GetMaximum()*1.5)

                rrv_number_dataset = RooRealVar("rrv_number_dataset", "rrv_number_dataset", rdataset.sumEntries())
                rrv_number_dataset.setError(0.)
                draw_error_band(rdataset, model_pdf, rrv_number_dataset, rfresult_pdf, mplot_deco, self.color_palet["Uncertainty"], "F") ## don't store the number in the workspace

            self.leg = self.legend4Plot(mplot_deco, 0) ## add the legend
            mplot_deco.addObject(self.leg)

            #self.draw_canvas(mplot_deco, "plots_%s_%s_%s/other/"%(options.additioninformation, self.channel, self.wtagger_category), "m_lvj"+label+in_range+in_range+mlvj_model+"_deco", 0, logy)
            self.draw_canvas(mplot_deco, self.plotsDir+"/other/", "m_lvj"+label+in_range+in_range+mlvj_model+"_deco", 0, logy)

        ### Number of the event in the dataset and lumi scale factor --> set the proper number for bkg extraction or for signal region
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").Print()
        self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).Print()
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").setVal(self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal())
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").setError(self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal())

        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").Print()


    #### method to fit the WJets normalization inside the mj signal region -> and write the jets mass sys if available
    def fit_WJetsNorm(self, scaleJetMass = 0): # to get the normalization of WJets in signal_region

        print "############### Fit mj Normalization ##################"
        ## fit the two version of pdf for Wjets shape if available
        self.fit_WJetsNormalization_in_Mj_signal_region("_WJets0_xww")
        self.fit_WJetsNormalization_in_Mj_signal_region("_WJets01_xww")

        ## in case fit also the scaled jet mass distributions in order to have the jet mass scale sys included
        if scaleJetMass :
            self.fit_WJetsNormalization_in_Mj_signal_region("_WJets0_xww_massup", "massup")
            self.fit_WJetsNormalization_in_Mj_signal_region("_WJets0_xww_massdn", "massdn")
            self.fit_WJetsNormalization_in_Mj_signal_region("_WJets1_xww")

        ## take the normalization numbers
        rrv_WJets0  = self.workspace4fit_.var("rrv_number_WJets0_xww_in_mj_signal_region_from_fitting_%s"%(self.channel))
        rrv_WJets01 = self.workspace4fit_.var("rrv_number_WJets01_xww_in_mj_signal_region_from_fitting_%s"%(self.channel))
        rrv_WJets0.Print()
        rrv_WJets01.Print()
        if scaleJetMass :
            rrv_WJets1 = self.workspace4fit_.var("rrv_number_WJets1_xww_in_mj_signal_region_from_fitting_%s"%(self.channel))
            rrv_WJets1.Print()
            #rrv_WJets0massup.Print()
            #rrv_WJets0massdn.Print()

        ### total uncertainty combining the result with two different shapes
        total_uncertainty = TMath.Sqrt(TMath.Power(rrv_WJets0.getError(), 2) + TMath.Power(rrv_WJets01.getVal()-rrv_WJets0.getVal(), 2))
        rrv_WJets0.setError(total_uncertainty)
        rrv_WJets0.Print()

        ##jet mass uncertainty on WJets normalization and the other bkg component
        if self.workspace4fit_.var("rrv_number_WJets0_xww_massup_in_mj_signal_region_from_fitting_%s"%(self.channel)) and self.workspace4fit_.var("rrv_number_WJets0_xww_massdn_in_mj_signal_region_from_fitting_%s"%(self.channel)):
            rrv_WJets0massup = self.workspace4fit_.var("rrv_number_WJets0_xww_massup_in_mj_signal_region_from_fitting_%s"%(self.channel))
            rrv_WJets0massdn = self.workspace4fit_.var("rrv_number_WJets0_xww_massdn_in_mj_signal_region_from_fitting_%s"%(self.channel))
            self.WJets_normlization_uncertainty_from_jet_mass = (TMath.Abs(rrv_WJets0massup.getVal()-rrv_WJets0.getVal())+TMath.Abs(rrv_WJets0massdn.getVal()-rrv_WJets0.getVal()))/2./rrv_WJets0.getVal()

        rrv_STop  = self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww__%s_mj"%(self.channel))

        if self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww_massup_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww_massdn_%s_mj"%(self.channel)) :
            rrv_STopmassup = self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww_massup_%s_mj"%(self.channel))
            rrv_STopmassdn = self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww_massdn_%s_mj"%(self.channel))
            self.STop_normlization_uncertainty_from_jet_mass = (TMath.Abs(rrv_STopmassup.getVal()-rrv_STop.getVal())+TMath.Abs(rrv_STopmassdn.getVal()-rrv_STop.getVal()))/2./rrv_STop.getVal()

        rrv_TTbar = self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww__%s_mj"%(self.channel))
        if self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww_massup_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww_massdn_%s_mj"%(self.channel)):
            rrv_TTbarmassup = self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww_massup_%s_mj"%(self.channel))
            rrv_TTbarmassdn = self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww_massdn_%s_mj"%(self.channel))
            self.TTbar_normlization_uncertainty_from_jet_mass = (TMath.Abs(rrv_TTbarmassup.getVal()-rrv_TTbar.getVal())+TMath.Abs(rrv_TTbarmassdn.getVal()-rrv_TTbar.getVal()))/2./rrv_TTbar.getVal()

        rrv_VV = self.workspace4fit_.var("rrv_number_dataset_signal_region_VV_xww_%s_mj"%(self.channel))
        if self.workspace4fit_.var("rrv_number_dataset_signal_region_VV_xww_massup_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_VV_xww_massdn_%s_mj"%(self.channel)):
            rrv_VVmassup = self.workspace4fit_.var("rrv_number_dataset_signal_region_VV_xww_massup_%s_mj"%(self.channel))
            rrv_VVmassdn = self.workspace4fit_.var("rrv_number_dataset_signal_region_VV_xww_massdn_%s_mj"%(self.channel))
            self.VV_normlization_uncertainty_from_jet_mass = (TMath.Abs(rrv_VVmassup.getVal()-rrv_VV.getVal())+TMath.Abs(rrv_VVmassdn.getVal()-rrv_VV.getVal()))/2./rrv_VV.getVal()

    #### make the mj sideband fit on data ti get the Wjets normaliztion
    def fit_WJetsNormalization_in_Mj_signal_region(self, label, massscale = ""):

        print "############### Fit mj Normalization: ", label, " ", massscale, " ##################"
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
        ## get real data in mj distribution --> mass up and down have only an effect on Wjets shape -> effect on the normalization -> evaluated in the MC and fit data
        rdataset_data_mj = self.workspace4fit_.data("rdataset_data_xww_%s_mj"%(self.channel))

        ### Fix TTbar, VV and STop
        model_TTbar = self.get_TTbar_mj_Model("_TTbar_xww"+massscale)
        model_STop  = self.get_STop_mj_Model("_STop_xww"+massscale)
        model_VV    = self.get_VV_mj_Model("_VV_xww"+massscale)
        ## only two parameters are fix, offset and width while the exp is floating , otherwise if shape different User1 or ErfExp everything is flaoting
        model_WJets = self.get_WJets_mj_Model(label)

        ## Total Pdf and fit only in sideband
        model_data = RooAddPdf("model_data_xww%s_%s_mj"%(massscale, self.channel), "model_data_xww%s_%s_mj"%(massscale, self.channel), RooArgList(model_WJets, model_VV, model_TTbar, model_STop))
        rfresult = model_data.fitTo(rdataset_data_mj, RooFit.Save(1) , RooFit.Range("sb_lo,sb_hi") , RooFit.Extended(kTRUE), RooFit.NumCPU(2), RooFit.Minimizer("Minuit"))
        rfresult = model_data.fitTo(rdataset_data_mj, RooFit.Save(1) , RooFit.Range("sb_lo,sb_hi") , RooFit.Extended(kTRUE), RooFit.NumCPU(2), RooFit.Minimizer("Minuit2"))
        rfresult.Print()
        if options.MCStudy:
            MCStudy(model_data, rrv_mass_j, "%s/other/%s/"%(self.plotsDir, self.signal_sample)+"try_"+model_data.GetName()+".png")
        rfresult.covarianceMatrix().Print()
        getattr(self.workspace4fit_, "import")(model_data)

        ## Total numver of event
        rrv_number_data_mj = RooRealVar("rrv_number_data_xww%s_%s_mj"%(massscale, self.channel), "rrv_number_data_xww%s_%s_mj"%(massscale, self.channel),
                self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale, self.channel)).getVal()+
                self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale, self.channel)).getVal()+
                self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale, self.channel)).getVal()+
                self.workspace4fit_.var("rrv_number%s_%s_mj"%(label, self.channel)).getVal())

        rrv_number_data_mj.setError(TMath.Sqrt(self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale, self.channel)).getError()*
            self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale, self.channel)).getError()+
            self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale, self.channel)).getError()*
            self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale, self.channel)).getError()+
            self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale, self.channel)).getError()*
            self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale, self.channel)).getError()+
            self.workspace4fit_.var("rrv_number%s_%s_mj"%(label, self.channel)).getError()*
            self.workspace4fit_.var("rrv_number%s_%s_mj"%(label, self.channel)).getError()))
        getattr(self.workspace4fit_, "import")(rrv_number_data_mj)

        ## if fit on Wjets default with the default shape
        if TString(label).Contains("_WJets0"):

            ## make the final plot
            mplot = rrv_mass_j.frame(RooFit.Title(""), RooFit.Bins(int(rrv_mass_j.getBins()/self.binwidth_narrow_factor)))
            rdataset_data_mj.plotOn(mplot, RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0))

            ## plot solid style
            model_data.plotOn(mplot, RooFit.Name("VV"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj,model_VV_xww_%s_mj"%(label, self.channel, self.channel, self.channel, self.channel)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())
            model_data.plotOn(mplot, RooFit.Name("TTbar"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj"%(label, self.channel, self.channel, self.channel)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())
            model_data.plotOn(mplot, RooFit.Name("STop"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj"%(label, self.channel, self.channel)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())
            model_data.plotOn(mplot, RooFit.Name("WJets"), RooFit.Components("model%s_%s_mj"%(label, self.channel)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())

            ## plot "dashed" style area
            model_data.plotOn(mplot, RooFit.Name("VV_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj,model_VV_xww_%s_mj"%(label, self.channel, self.channel, self.channel, self.channel)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.FillStyle(3002), RooFit.Range(rrv_mass_j.getMin(), rrv_mass_j.getMax()), RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())
            model_data.plotOn(mplot, RooFit.Name("TTbar_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj"%(label, self.channel, self.channel, self.channel)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.FillStyle(3002), RooFit.Range(rrv_mass_j.getMin(), rrv_mass_j.getMax()), RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())
            model_data.plotOn(mplot, RooFit.Name("STop_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj"%(label, self.channel, self.channel)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.FillStyle(3002), RooFit.Range(rrv_mass_j.getMin(), rrv_mass_j.getMax()), RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())
            model_data.plotOn(mplot, RooFit.Name("WJets_invisible"), RooFit.Components("model%s_%s_mj"%(label, self.channel)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.FillStyle(3002), RooFit.Range(rrv_mass_j.getMin(), rrv_mass_j.getMax()), RooFit.LineColor(kBlack), RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())

            ### solid line
            model_data.plotOn(mplot, RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj"%(label, self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2) , RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())
            model_data.plotOn(mplot, RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj"%(label, self.channel, self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())
            model_data.plotOn(mplot, RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj"%(label, self.channel, self.channel, self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2) , RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())
            model_data.plotOn(mplot, RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj,model_VV_%s_mj"%(label, self.channel, self.channel, self.channel, self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2) , RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())

            ### dash line
            model_data.plotOn(mplot, RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj"%(label, self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(), rrv_mass_j.getMax()), RooFit.LineStyle(kDashed) , RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())
            model_data.plotOn(mplot, RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj"%(label, self.channel, self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(), rrv_mass_j.getMax()), RooFit.LineStyle(kDashed) , RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())
            model_data.plotOn(mplot, RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj"%(label, self.channel, self.channel, self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(), rrv_mass_j.getMax()), RooFit.LineStyle(kDashed) , RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())
            model_data.plotOn(mplot, RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj,model_VV_xww_%s_mj"%(label, self.channel, self.channel, self.channel, self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(), rrv_mass_j.getMax()), RooFit.LineStyle(kDashed) , RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines())

            model_data.plotOn(mplot, RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj,model_VV_xww_%s_mj"%(label, self.channel, self.channel, self.channel, self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(), rrv_mass_j.getMax()), RooFit.LineStyle(kDashed) , RooFit.NormRange("sb_lo,sb_hi"))

            rdataset_data_mj.plotOn(mplot, RooFit.Name("data"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0))

            ### draw the error band using the sum of all the entries component MC + fit
            draw_error_band(rdataset_data_mj, model_data, rrv_number_data_mj, rfresult, mplot, self.color_palet["Uncertainty"], "F")
            rdataset_data_mj.plotOn(mplot, RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0))

            ### Get the pull and plot it
            mplot_pull = self.get_pull(rrv_mass_j, mplot)

            ### signal window zone with vertical lines
            lowerLine = TLine(self.mj_signal_min, 0., self.mj_signal_min, mplot.GetMaximum()*0.9)
            lowerLine.SetLineWidth(2)
            lowerLine.SetLineColor(kBlack)
            lowerLine.SetLineStyle(9)
            upperLine = TLine(self.mj_signal_max, 0., self.mj_signal_max, mplot.GetMaximum()*0.9)
            upperLine.SetLineWidth(2)
            upperLine.SetLineColor(kBlack)
            upperLine.SetLineStyle(9)
            #lowerLine = TLine(65, 0., 65, mplot.GetMaximum()*0.9); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kBlack); lowerLine.SetLineStyle(9)
            #middleLine = TLine(85, 0., 85, mplot.GetMaximum()*0.9); middleLine.SetLineWidth(2); middleLine.SetLineColor(kBlack); middleLine.SetLineStyle(9)
            #upperLine = TLine(95, 0., 95, mplot.GetMaximum()*0.9); upperLine.SetLineWidth(2); upperLine.SetLineColor(kBlack); upperLine.SetLineStyle(9)
            mplot.addObject(lowerLine)
            #mplot.addObject(middleLine)
            mplot.addObject(upperLine)

            pt = ROOT.TPaveText(0.3592965, 0.02847153, 0.50, 0.1008991, "NDC")
            pt.SetTextFont(42)
            pt.SetTextSize(0.04995005)
            #pt.SetTextAlign(12)
            pt.SetFillColor(0)
            pt.SetBorderSize(0)
            text = pt.AddText("#leftarrow signal region #rightarrow")
            text.SetTextFont(62)
            #mplot.addObject(pt)

            pt1 = ROOT.TPaveText(0.3555276, 0.1183816, 0.4535176, 0.1908092, "NDC")
            pt1.SetTextFont(42)
            pt1.SetTextSize(0.037)
            #pt.SetTextAlign(12)
            pt1.SetFillColor(0)
            pt1.SetFillStyle(0)
            pt1.SetBorderSize(0)
            text = pt1.AddText("WW")
            text.SetTextFont(62)
            text = pt1.AddText("category")
            text.SetTextFont(62)
            #mplot.addObject(pt1)

            pt2 = ROOT.TPaveText(0.4723618, 0.1183816, 0.5678392, 0.1908092, "NDC")
            pt2.SetTextFont(42)
            pt2.SetTextSize(0.037)
            #pt.SetTextAlign(12)
            pt2.SetFillColor(0)
            pt2.SetBorderSize(0)
            pt2.SetFillStyle(0)
            text = pt2.AddText("WZ")
            text.SetTextFont(62)
            text = pt2.AddText("category")
            text.SetTextFont(62)
            #mplot.addObject(pt2)

            ### legend of the plot
            self.leg = self.legend4Plot(mplot, 0, 1, 0., 0., 0.13, 0.02, 1, 0, 1)
            #self.leg = self.legend4Plot(mplot, 0, 1,-0.10,-0.01, 0.10, 0.01)
            mplot.addObject(self.leg)
            mplot.GetYaxis().SetRangeUser(1e-2, mplot.GetMaximum()*1.8)

            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize()
            nBinX = mplot.GetNbinsX()
            ndof  = nBinX-self.nPar_float_in_fitTo
            print mplot.chiSquare()
            print "#################### JENchi2 nPar = %s, chiSquare = %s/%s"%(self.nPar_float_in_fitTo , mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof)
            datahist = rdataset_data_mj.binnedClone(rdataset_data_mj.GetName()+"_binnedClone", rdataset_data_mj.GetName()+"_binnedClone")

            parameters_list = model_data.getParameters(rdataset_data_mj)
            self.draw_canvas_with_pull1(rrv_mass_j, datahist, mplot, mplot_pull, ndof, parameters_list, "%s/m_j_fitting/"%(self.plotsDir), "m_j_sideband%s"%(label), "", 1, 0, 1)

            ### call the function for getting the normalizatio in signal region for data, TTbar, STop, VV and W+jets = label -> store in a output txt file
            self.get_mj_normalization_insignalregion("_data_xww")
            self.get_mj_normalization_insignalregion("_TTbar_xww")
            self.get_mj_normalization_insignalregion("_STop_xww")
            self.get_mj_normalization_insignalregion("_VV_xww")
            self.get_mj_normalization_insignalregion(label)

        #### to calculate the WJets's normalization and error in M_J signal_region. The error must contain the shape error: model_WJets have new parameters fitting data
        fullInt   = model_WJets.createIntegral(RooArgSet(rrv_mass_j), RooArgSet(rrv_mass_j))
        signalInt = model_WJets.createIntegral(RooArgSet(rrv_mass_j), RooArgSet(rrv_mass_j), ("signal_region"))
        fullInt_val = fullInt.getVal()
        signalInt_val = signalInt.getVal()/fullInt_val
        ## take the value from the fit (normalization) and multiply it from the ratio of the integrals
        rrv_number_WJets_in_mj_signal_region_from_fitting = RooRealVar("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label, self.channel), "rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label, self.channel), self.workspace4fit_.var("rrv_number%s_%s_mj"%(label, self.channel)).getVal()*signalInt_val)

        #### Error on the normalization --> from a dedicated function taking into account shape uncertainty
        rrv_number_WJets_in_mj_signal_region_from_fitting.setError(Calc_error_extendPdf(rdataset_data_mj, model_WJets, rfresult, "signal_region"))
        print "########## error on the normalization due to shape + norm = %s"%(rrv_number_WJets_in_mj_signal_region_from_fitting.getError())
        getattr(self.workspace4fit_, "import")(rrv_number_WJets_in_mj_signal_region_from_fitting)
        rrv_number_WJets_in_mj_signal_region_from_fitting.Print()


    ##### Counting of the events of each component in the signal region taking the lavel for the model
    def get_mj_normalization_insignalregion(self, label):
        print "################## get mj normalization ", label, " ################## "
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
        model      = self.workspace4fit_.pdf("model"+label+"_"+self.channel+"_mj")

        fullInt   = model.createIntegral(RooArgSet(rrv_mass_j), RooArgSet(rrv_mass_j))
        sb_loInt  = model.createIntegral(RooArgSet(rrv_mass_j), RooArgSet(rrv_mass_j), ("sb_lo"))
        signalInt = model.createIntegral(RooArgSet(rrv_mass_j), RooArgSet(rrv_mass_j), ("signal_region"))
        sb_hiInt  = model.createIntegral(RooArgSet(rrv_mass_j), RooArgSet(rrv_mass_j), ("sb_hi"))

        fullInt_val   = fullInt.getVal()
        sb_loInt_val  = sb_loInt.getVal()/fullInt_val
        sb_hiInt_val  = sb_hiInt.getVal()/fullInt_val
        signalInt_val = signalInt.getVal()/fullInt_val

        print "########### Events Number in MC Dataset: #############"
        self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").Print()

        print "########### Events Number get from fit: ##############"
        rrv_tmp = self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj")
        rrv_tmp.Print()
        print "Events Number in sideband_low :%s"%(rrv_tmp.getVal()*sb_loInt_val)
        print "Events Number in Signal Region:%s"%(rrv_tmp.getVal()*signalInt_val)
        print "Events Number in sideband_high:%s"%(rrv_tmp.getVal()*sb_hiInt_val)
        print "Total Number in sidebands :%s"%(rrv_tmp.getVal()*(sb_loInt_val+sb_hiInt_val))
        print "Ratio signal_region/sidebands :%s"%(signalInt_val/(sb_loInt_val+sb_hiInt_val))

        ##### Save numbers in the output text file
        self.file_out.write("\n%s++++++++++++++++++++++++++++++++++++"%(label))
        self.file_out.write("\nEvents Number in sideband_low from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal()))
        self.file_out.write("\nEvents Number in Signal Region from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").getVal()))
        self.file_out.write("\nEvents Number in sideband_high from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal()))
        self.file_out.write("\nTotal Number in sidebands from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal()+ self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal()))
        self.file_out.write("\nRatio signal_region/sidebands from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").getVal()/(self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal()+ self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal())))

        self.file_out.write("\nEvents Number in sideband_low from fitting:%s"%(rrv_tmp.getVal()*sb_loInt_val))
        self.file_out.write("\nEvents Number in Signal Region from fitting:%s"%(rrv_tmp.getVal()*signalInt_val))
        self.file_out.write("\nEvents Number in sideband_high from fitting:%s"%(rrv_tmp.getVal()*sb_hiInt_val))
        self.file_out.write("\nTotal Number in sidebands from fitting:%s"%(rrv_tmp.getVal()*(sb_loInt_val+sb_hiInt_val)))
        self.file_out.write("\nRatio signal_region/sidebands from fitting:%s"%(signalInt_val/(sb_loInt_val+sb_hiInt_val)))

    ##### Method to fit data mlvj shape in the sideband -> first step for the background extraction of the shape
    def fit_mlvj_in_Mj_sideband(self, label, mlvj_region, mlvj_model, logy = 0):


        print "############### Fit mlvj in mj sideband: ", label, " ", mlvj_region, "  ", mlvj_model, " ##################"
        rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j")
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rdataset_data_mlvj = self.workspace4fit_.data("rdataset_data_xww%s_%s_mlvj"%(mlvj_region, self.channel))

        ## get the minor component shapes in the sb low
        model_VV_backgrounds    = self.get_VV_mlvj_Model("_sb_lo")
        number_VV_sb_lo_mlvj    = self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel))
        model_TTbar_backgrounds = self.get_TTbar_mlvj_Model("_sb_lo")
        number_TTbar_sb_lo_mlvj = self.workspace4fit_.var("rrv_number_TTbar_xww_sb_lo_%s_mlvj"%(self.channel))
        model_STop_backgrounds  = self.get_STop_mlvj_Model("_sb_lo")
        number_STop_sb_lo_mlvj  = self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel))

        self.workspace4fit_.var("rrv_number_TTbar_xww_sb_lo_%s_mlvj"%(self.channel)).Print()
        self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel)).Print()
        self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel)).Print()

        ### Make the Pdf for the WJets
        model_pdf_WJets = self.make_Pdf("%s_sb_lo_from_fitting"%(label), mlvj_model, "_mlvj")
        model_pdf_WJets.Print()
        ### inititalize the value to what was fitted with the mc in the sideband
        number_WJets_sb_lo = self.workspace4fit_.var("rrv_number%s_sb_lo_%s_mlvj"%(label, self.channel)).clone("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label, self.channel))
        model_WJets = RooExtendPdf("model%s_sb_lo_from_fitting_%s_mlvj"%(label, self.channel), "model%s_sb_lo_from_fitting_%s_mlvj"%(label, self.channel), model_pdf_WJets, number_WJets_sb_lo)
        model_pdf_WJets.Print()
        number_WJets_sb_lo.Print()

        ## Add the other bkg component fixed to the total model
        model_data = RooAddPdf("model_data%s%s_mlvj"%(label, mlvj_region), "model_data%s%s_mlvj"%(label, mlvj_region), RooArgList(model_WJets, model_VV_backgrounds, model_TTbar_backgrounds, model_STop_backgrounds))

        rfresult = model_data.fitTo(rdataset_data_mlvj, RooFit.Save(1) , RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit"))
        rfresult = model_data.fitTo(rdataset_data_mlvj, RooFit.Save(1) , RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"))
        rfresult.Print()
        if options.MCStudy:
            MCStudy(model_data, rrv_mass_lvj, "%s/other/%s/"%(self.plotsDir, self.signal_sample)+"MCStudy_"+model_data.GetName()+".png", 400, int(rdataset_data_mlvj.sumEntries()))
        rfresult.covarianceMatrix().Print()
        getattr(self.workspace4fit_, "import")(model_data)

        model_WJets.Print()
        model_WJets.getParameters(rdataset_data_mlvj).Print("v")
        self.workspace4fit_.pdf("model_pdf%s_sb_lo_%s_mlvj"%(label, self.channel)).getParameters(rdataset_data_mlvj).Print("v")

        ### data in the sideband plus error from fit
        rrv_number_data_sb_lo_mlvj = RooRealVar("rrv_number_data_xww_sb_lo_%s_mlvj"%(self.channel), "rrv_number_data_xww_sb_lo_%s_mlvj"%(self.channel),
                self.workspace4fit_.var("rrv_number_TTbar_xww_sb_lo_%s_mlvj"%(self.channel)).getVal()+
                self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel)).getVal()+
                self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel)).getVal()+
                self.workspace4fit_.var("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label, self.channel)).getVal())

        rrv_number_data_sb_lo_mlvj.setError(TMath.Sqrt(self.workspace4fit_.var("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label, self.channel)).getError()*
            self.workspace4fit_.var("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label, self.channel)).getError()+
            self.workspace4fit_.var("rrv_number_TTbar_xww_sb_lo_%s_mlvj"%(self.channel)).getError()*
            self.workspace4fit_.var("rrv_number_TTbar_xww_sb_lo_%s_mlvj"%(self.channel)).getError()+
            self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel)).getError()*
            self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel)).getError()+
            self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel)).getError()*
            self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel)).getError()))

        getattr(self.workspace4fit_, "import")(rrv_number_data_sb_lo_mlvj)

        ### plot for WJets default + default shape
        if TString(label).Contains("_WJets0"):

            mplot = rrv_mass_lvj.frame(RooFit.Title("M_lvj fitted in M_j sideband "), RooFit.Bins(int(rrv_mass_lvj.getBins()/self.binwidth_narrow_factor)))

            rdataset_data_mlvj.plotOn(mplot , RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0))

            model_data.plotOn(mplot, RooFit.Components("model%s_sb_lo_from_fitting_%s_mlvj,model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_%s_mlvj"%(label, self.channel, self.channel, self.channel, self.channel)), RooFit.Name("WJets"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_%s_mlvj"%(self.channel, self.channel, self.channel)), RooFit.Name("VV"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj"%(self.channel, self.channel)), RooFit.Name("TTbar"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines())

            model_data.plotOn(mplot, RooFit.Components("model_STop_xww_sb_lo_%s_mlvj"%(self.channel)), RooFit.Name("STop"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())

            #solid line
            model_data.plotOn(mplot, RooFit.Components("model%s_sb_lo_from_fitting_%s_mlvj,model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_%s_mlvj"%(label, self.channel, self.channel, self.channel, self.channel)), RooFit.Name("WJets_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_%s_mlvj"%(self.channel, self.channel, self.channel)), RooFit.Name("VV_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj"%(self.channel, self.channel)), RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

            model_data.plotOn(mplot, RooFit.Components("model_STop_xww_sb_lo_%s_mlvj"%(self.channel)), RooFit.Name("STop_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())


            ### draw the error band
            draw_error_band(rdataset_data_mlvj, model_data, self.workspace4fit_.var("rrv_number_data_xww_sb_lo_%s_mlvj"%(self.channel)) , rfresult, mplot, self.color_palet["Uncertainty"], "F")
            model_data.plotOn(mplot , RooFit.VLines(), RooFit.Invisible())
            model_data.plotOn(mplot , RooFit.Invisible())
            self.plot_data_with_poissoninterval(rrv_mass_lvj, rdataset_data_mlvj, mplot)

            mplot.GetYaxis().SetRangeUser(1e-2, mplot.GetMaximum()*1.5)

            ### Add the legend to the plot
            self.leg = self.legend4Plot(mplot, 0, 1, 0., 0.06, 0.16, 0.)
            mplot.addObject(self.leg)

            ### calculate the chi2
            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize()
            nBinX = mplot.GetNbinsX()
            ndof  = nBinX-self.nPar_float_in_fitTo
            print mplot.chiSquare()
            print "#################### nPar = %s, chiSquare = %s/%s"%(self.nPar_float_in_fitTo , mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof)
            ### write the result in the output
            self.file_out.write("\n fit_mlvj_in_Mj_sideband: nPar = %s, chiSquare = %s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof))

            ### get the pull plot and store the canvas
            mplot_pull = self.get_pull(rrv_mass_lvj, mplot)
            parameters_list = model_data.getParameters(rdataset_data_mlvj)
            datahist = rdataset_data_mlvj.binnedClone(rdataset_data_mlvj.GetName()+"_binnedClone", rdataset_data_mlvj.GetName()+"_binnedClone")

            self.draw_canvas_with_pull1(rrv_mass_lvj, datahist, mplot, mplot_pull, ndof, parameters_list, "%s/m_lvj_fitting/"%(self.plotsDir), "m_lvj_sb_lo%s"%(label), mlvj_model, 1, 1)

        #### Decorrelate the parameters in order to have a proper shape in the workspace
        wsfit_tmp = RooWorkspace("wsfit_tmp%s_sb_lo_from_fitting_mlvj"%(label))
        Deco      = PdfDiagonalizer("Deco%s_sb_lo_from_fitting_%s_%s_mlvj"%(label, self.channel, self.wtagger_category), wsfit_tmp, rfresult)
        print"#################### diagonalize data sideband fit "
        model_pdf_WJets_deco = Deco.diagonalize(model_pdf_WJets)
        print"#################### print parameters "
        model_pdf_WJets_deco.Print("v")
        model_pdf_WJets_deco.getParameters(rdataset_data_mlvj).Print("")
        getattr(self.workspace4fit_, "import")(model_pdf_WJets_deco)

        #### Call the alpha evaluation in automatic
        self.get_WJets_mlvj_correction_sb_lo_to_signal_region(label, mlvj_model)

        ### Fix the pdf of signal, TTbar, STop and VV in the signal region
        self.fix_Model("_%s_xww"%(self.signal_sample), "_signal_region", "_mlvj")
        self.fix_Model("_TTbar_xww", "_signal_region", "_mlvj")
        self.fix_Model("_STop_xww", "_signal_region", "_mlvj")
        self.fix_Model("_VV_xww", "_signal_region", "_mlvj")

        ### Call the evaluation of the normalization in the signal region for signal, TTbar, VV, STop, and WJets after the extrapolation via alpha
        self.get_mlvj_normalization_insignalregion("_%s_xww"%(self.signal_sample))
        self.get_mlvj_normalization_insignalregion("_TTbar_xww")
        self.get_mlvj_normalization_insignalregion("_STop_xww")
        self.get_mlvj_normalization_insignalregion("_VV_xww")
        self.get_mlvj_normalization_insignalregion(label, "model_pdf%s_signal_region_%s_after_correct_mlvj"%(label, self.channel))


    ##### Function that calculate the normalization inside the mlvj signal region (mass window around the resonance in order to fill datacards)
    def get_mlvj_normalization_insignalregion(self, label, model_name = ""):

        print "############### get mlvj normalization inside SR ", label, " ", model_name, " ##################"
        if model_name == "":
            model = self.workspace4fit_.pdf("model"+label+"_signal_region"+"_"+self.channel+"_mlvj")
        else:
            model = self.workspace4fit_.pdf(model_name)

        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")

        fullInt   = model.createIntegral(RooArgSet(rrv_mass_lvj), RooArgSet(rrv_mass_lvj))
        signalInt = model.createIntegral(RooArgSet(rrv_mass_lvj), RooArgSet(rrv_mass_lvj), ("signal_region"))
        highMassInt = model.createIntegral(RooArgSet(rrv_mass_lvj), RooArgSet(rrv_mass_lvj), ("high_mass"))

        fullInt_val = fullInt.getVal()
        signalInt_val = signalInt.getVal()/fullInt_val
        highMassInt_val = highMassInt.getVal()/fullInt_val

        ## integal in the signal region
        print "######### integral in SR: ", label+"signalInt = %s"%(signalInt_val)

        print "####### Events Number in MC Dataset:"
        self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj").Print()
        self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj").Print()

        print "########## Events Number get from fit:"
        rrv_tmp = self.workspace4fit_.var("rrv_number"+label+"_signal_region"+"_"+self.channel+"_mlvj")
        print "Events Number in Signal Region from fitting: %s"%(rrv_tmp.getVal()*signalInt_val)

        #### store the info in the output file
        self.file_out.write("\n%s++++++++++++++++++++++++++++++++++++"%(label))
        self.file_out.write("\nEvents Number in All Region from dataset : %s"%(self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj").getVal()))
        self.file_out.write("\nEvents Number in Signal Region from dataset: %s"%(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj").getVal()))
        self.file_out.write("\nRatio signal_region/all_range from dataset :%s"%(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj").getVal()/self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj").getVal()))
        self.file_out.write("\nEvents Number in All Region from fitting : %s\n"%(rrv_tmp.getVal()))
        self.file_out.write("\nEvents Number in Signal Region from fitting: %s\n"%(rrv_tmp.getVal()*signalInt_val))
        self.file_out.write("\nEvents Number in High Mass Region from fitting: %s\n"%(rrv_tmp.getVal()*highMassInt_val))
        self.file_out.write("\nRatio signal_region/all_range from fitting :%s"%(signalInt_val))

        if not self.workspace4fit_.var("rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj"):
            rrv_number_fitting_signal_region_mlvj = RooRealVar("rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj", "rrv_number_fitting_signal_region"+label+"_"+
                    self.channel+"_mlvj", rrv_tmp.getVal()*signalInt_val)
            getattr(self.workspace4fit_, "import")(rrv_number_fitting_signal_region_mlvj)
        else :
            self.workspace4fit_.var("rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj").setVal(rrv_tmp.getVal()*signalInt_val)

        self.workspace4fit_.var("rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj").Print()

    ### method to get the alpha function to extrapolate the wjets in the signal region
    def get_WJets_mlvj_correction_sb_lo_to_signal_region(self, label, mlvj_model):

        print" ############# get the extrapolation function alpha from MC : ", label, "   ", mlvj_model, " ###############"
        #tmp_Style = self.tdrStyle.Clone("tmp_Style")
        #tmp_Style.SetPadRightMargin(0.08)
        #tmp_Style.SetPadTickY(0)
        #tmp_Style.cd()

        ### take input var and datasets from 4fit collection --> mc not scaled to lumi --> just a shape here
        rrv_x = self.workspace4fit_.var("rrv_mass_lvj")
        rdataset_WJets_sb_lo_mlvj = self.workspace4fit_.data("rdataset4fit%s_sb_lo_%s_mlvj"%(label, self.channel))
        rdataset_WJets_signal_region_mlvj = self.workspace4fit_.data("rdataset4fit%s_signal_region_%s_mlvj"%(label, self.channel))
        rdataset_WJets_sb_lo_mlvj.Print("v")
        rdataset_WJets_signal_region_mlvj.Print("v")

        ### create a frame for the next plots
        mplot = rrv_x.frame(RooFit.Title("correlation_pdf"), RooFit.Bins(int(rrv_x.getBins()/self.binwidth_narrow_factor)))
        mplot.GetYaxis().SetTitle("arbitrary units")

        ### model used for Higgs analysis --> parameters in the SR has to be fitted, not yet done in order to take into account correlations between mj and mlvj
        alpha_constrains = RooArgSet()
        if mlvj_model == "ErfExp_v1":

            rrv_c_sb       = self.workspace4fit_.var("rrv_c_ErfExp%s_sb_lo_%s"%(label, self.channel))
            rrv_offset_sb  = self.workspace4fit_.var("rrv_offset_ErfExp%s_sb_lo_%s"%(label, self.channel))
            rrv_width_sb   = self.workspace4fit_.var("rrv_width_ErfExp%s_sb_lo_%s"%(label, self.channel))

            rrv_delta_c      = RooRealVar("rrv_delta_c_ErfExp%s_%s"%(label, self.channel), "rrv_delta_c_ErfExp%s_%s"%(label, self.channel), 0.,
                    -100*rrv_c_sb.getError(), 100*rrv_c_sb.getError())
            rrv_delta_offset = RooRealVar("rrv_delta_offset_ErfExp%s_%s"%(label, self.channel), "rrv_delta_offset_ErfExp%s_%s"%(label, self.channel), 0.,
                    -100*rrv_offset_sb.getError(), 100*rrv_offset_sb.getError())
            rrv_delta_width = RooRealVar("rrv_delta_width_ErfExp%s_%s"%(label, self.channel), "rrv_delta_width_ErfExp%s_%s"%(label, self.channel), 0.,
                    -100*rrv_width_sb.getError(), 100*rrv_width_sb.getError())

            rrv_c_sr      = RooFormulaVar("rrv_c_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_c_sb, rrv_delta_c))
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_offset_sb, rrv_delta_offset))
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_width_sb, rrv_delta_width))

            correct_factor_pdf = RooAlpha("correct_factor_pdf", "correct_factor_pdf", rrv_x, rrv_c_sr, rrv_offset_sr, rrv_width_sr, rrv_c_sb, rrv_offset_sb, rrv_width_sb, rrv_x.getMin(), rrv_x.getMax())

        if mlvj_model == "ErfPow_v1":

            rrv_c_sb      = self.workspace4fit_.var("rrv_c_ErfPow%s_sb_lo_%s"%(label, self.channel))
            rrv_offset_sb = self.workspace4fit_.var("rrv_offset_ErfPow%s_sb_lo_%s"%(label, self.channel))
            rrv_width_sb  = self.workspace4fit_.var("rrv_width_ErfPow%s_sb_lo_%s"%(label, self.channel))

            rrv_delta_c      = RooRealVar("rrv_delta_c_ErfPow%s_%s"%(label, self.channel), "rrv_delta_c_ErfPow%s_%s"%(label, self.channel), 0.,
                                          -100*rrv_c_sb.getError(), 100*rrv_c_sb.getError())
            rrv_delta_offset = RooRealVar("rrv_delta_offset_ErfPow%s_%s"%(label, self.channel), "rrv_delta_offset_ErfPow%s_%s"%(label, self.channel), 0.,
                                          -100*rrv_offset_sb.getError(), 100*rrv_offset_sb.getError())
            rrv_delta_width  = RooRealVar("rrv_delta_width_ErfPow%s_%s"%(label, self.channel), "rrv_delta_width_ErfPow%s_%s"%(label, self.channel), 0.,
                                          -100*rrv_width_sb.getError(), 100*rrv_width_sb.getError())

            rrv_c_sr      = RooFormulaVar("rrv_c_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_c_sb, rrv_delta_c))
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_offset_sb, rrv_delta_offset))
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_width_sb, rrv_delta_width))

            correct_factor_pdf = RooAlpha4ErfPowPdf("correct_factor_pdf", "correct_factor_pdf", rrv_x, rrv_c_sr, rrv_offset_sr, rrv_width_sr, rrv_c_sb, rrv_offset_sb, rrv_width_sb)

        if mlvj_model == "ErfPow2_v1":

            rrv_c0_sb     = self.workspace4fit_.var("rrv_c0_ErfPow2%s_sb_lo_%s"%(label, self.channel))
            rrv_c1_sb     = self.workspace4fit_.var("rrv_c1_ErfPow2%s_sb_lo_%s"%(label, self.channel))
            rrv_offset_sb = self.workspace4fit_.var("rrv_offset_ErfPow2%s_sb_lo_%s"%(label, self.channel))
            rrv_width_sb  = self.workspace4fit_.var("rrv_width_ErfPow2%s_sb_lo_%s"%(label, self.channel))

            c0_sb_constrain = RooGaussian ("c0_sb_constrain", "c0_sb_constrain", rrv_c0_sb, RooFit.RooConst(rrv_c0_sb.getVal()), RooFit.RooConst(rrv_c0_sb.getError()))
            alpha_constrains.add(c0_sb_constrain)
            c1_sb_constrain = RooGaussian ("c1_sb_constrain", "c1_sb_constrain", rrv_c1_sb, RooFit.RooConst(rrv_c1_sb.getVal()), RooFit.RooConst(rrv_c1_sb.getError()))
            alpha_constrains.add(c1_sb_constrain)
            offset_sb_constrain = RooGaussian ("offset_sb_constrain", "offset_sb_constrain", rrv_offset_sb, RooFit.RooConst(rrv_offset_sb.getVal()), RooFit.RooConst(rrv_offset_sb.getError()))
            alpha_constrains.add(offset_sb_constrain)
            width_sb_constrain = RooGaussian ("width_sb_constrain", "width_sb_constrain", rrv_width_sb, RooFit.RooConst(rrv_width_sb.getVal()), RooFit.RooConst(rrv_width_sb.getError()))
            alpha_constrains.add(width_sb_constrain)

            rrv_delta_c0  = RooRealVar("rrv_delta_c0_ErfPow2%s_%s"%(label, self.channel), "rrv_delta_c0_ErfPow2%s_%s"%(label, self.channel),
                    self.workspace4fit_.var("rrv_c0_ErfPow2%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c0_sb.getVal(),
                    self.workspace4fit_.var("rrv_c0_ErfPow2%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c0_sb.getVal()-4*rrv_c0_sb.getError(),
                    self.workspace4fit_.var("rrv_c0_ErfPow2%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c0_sb.getVal()+4*rrv_c0_sb.getError())

            rrv_delta_c1 = RooRealVar("rrv_delta_c1_ErfPow2%s_%s"%(label, self.channel), "rrv_delta_c1_ErfPow2%s_%s"%(label, self.channel),
                    self.workspace4fit_.var("rrv_c1_ErfPow2%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c1_sb.getVal(),
                    self.workspace4fit_.var("rrv_c1_ErfPow2%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c1_sb.getVal()-4*rrv_c1_sb.getError(),
                    self.workspace4fit_.var("rrv_c1_ErfPow2%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c1_sb.getVal()+4*rrv_c1_sb.getError())

            rrv_delta_offset = RooRealVar("rrv_delta_offset_ErfPow2%s_%s"%(label, self.channel), "rrv_delta_offset_ErfPow2%s_%s"%(label, self.channel),
                    self.workspace4fit_.var("rrv_offset_ErfPow2%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_offset_sb.getVal(),
                    self.workspace4fit_.var("rrv_offset_ErfPow2%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_offset_sb.getVal()-4*rrv_offset_sb.getError(),
                    self.workspace4fit_.var("rrv_offset_ErfPow2%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_offset_sb.getVal()+4*rrv_offset_sb.getError())

            rrv_delta_width = RooRealVar("rrv_delta_width_ErfPow2%s_%s"%(label, self.channel), "rrv_delta_width_ErfPow2%s_%s"%(label, self.channel),
                    self.workspace4fit_.var("rrv_width_ErfPow2%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_width_sb.getVal(),
                    self.workspace4fit_.var("rrv_width_ErfPow2%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_width_sb.getVal()-4*rrv_width_sb.getError(),
                    self.workspace4fit_.var("rrv_width_ErfPow2%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_width_sb.getVal()+4*rrv_width_sb.getError())


            rrv_c0_sr     = RooFormulaVar("rrv_c0_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_c0_sb, rrv_delta_c0))
            rrv_c1_sr     = RooFormulaVar("rrv_c1_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_c1_sb, rrv_delta_c1))
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_offset_sb, rrv_delta_offset))
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_width_sb, rrv_delta_width))

            correct_factor_pdf = RooAlpha4ErfPow2Pdf("correct_factor_pdf", "correct_factor_pdf", rrv_x, rrv_c0_sr, rrv_c1_sr, rrv_offset_sr, rrv_width_sr, rrv_c0_sb, rrv_c1_sb, rrv_offset_sb, rrv_width_sb)

        if mlvj_model == "ErfPowExp_v1": ## take initial value from what was already fitted in the SR

            rrv_c0_sb     = self.workspace4fit_.var("rrv_c0_ErfPowExp%s_sb_lo_%s"%(label, self.channel))
            rrv_c1_sb     = self.workspace4fit_.var("rrv_c1_ErfPowExp%s_sb_lo_%s"%(label, self.channel))
            rrv_offset_sb = self.workspace4fit_.var("rrv_offset_ErfPowExp%s_sb_lo_%s"%(label, self.channel))
            rrv_width_sb  = self.workspace4fit_.var("rrv_width_ErfPowExp%s_sb_lo_%s"%(label, self.channel))
            c0_sb_constrain = RooGaussian ("c0_sb_constrain", "c0_sb_constrain", rrv_c0_sb, RooFit.RooConst(rrv_c0_sb.getVal()), RooFit.RooConst(rrv_c0_sb.getError()))
            alpha_constrains.add(c0_sb_constrain)
            c1_sb_constrain = RooGaussian ("c1_sb_constrain", "c1_sb_constrain", rrv_c1_sb, RooFit.RooConst(rrv_c1_sb.getVal()), RooFit.RooConst(rrv_c1_sb.getError()))
            alpha_constrains.add(c1_sb_constrain)
            offset_sb_constrain = RooGaussian ("offset_sb_constrain", "offset_sb_constrain", rrv_offset_sb, RooFit.RooConst(rrv_offset_sb.getVal()), RooFit.RooConst(rrv_offset_sb.getError()))
            alpha_constrains.add(offset_sb_constrain)
            width_sb_constrain = RooGaussian ("width_sb_constrain", "width_sb_constrain", rrv_width_sb, RooFit.RooConst(rrv_width_sb.getVal()), RooFit.RooConst(rrv_width_sb.getError()))
            alpha_constrains.add(width_sb_constrain)

            tmp_rrv_c0_sr    = self.workspace4fit_.var("rrv_c0_ErfPowExp%s_signal_region_%s"%(label, self.channel))
            tmp_rrv_c1_sr    = self.workspace4fit_.var("rrv_c1_ErfPowExp%s_signal_region_%s"%(label, self.channel))
            tmp_rrv_offset_sr = self.workspace4fit_.var("rrv_offset_ErfPowExp%s_signal_region_%s"%(label, self.channel))
            tmp_rrv_width_sr = self.workspace4fit_.var("rrv_width_ErfPowExp%s_signal_region_%s"%(label, self.channel))

            rrv_delta_c0  = RooRealVar("rrv_delta_c0_ErfPowExp%s_%s"%(label, self.channel), "rrv_delta_c0_ErfPowExp%s_%s"%(label, self.channel),
                    tmp_rrv_c0_sr.getVal()-rrv_c0_sb.getVal(),
                    tmp_rrv_c0_sr.getVal()-rrv_c0_sb.getVal()-rrv_c0_sb.getError()-tmp_rrv_c0_sr.getError(),
                    tmp_rrv_c0_sr.getVal()-rrv_c0_sb.getVal()+rrv_c0_sb.getError()+tmp_rrv_c0_sr.getError())
            rrv_delta_c0.setError(TMath.Sqrt(rrv_c0_sb.getError()*rrv_c0_sb.getError() + tmp_rrv_c0_sr.getError()*tmp_rrv_c0_sr.getError()))
            #delta_c0_constrain = RooGaussian ("delta_c0_constrain", "delta_c0_constrain", rrv_delta_c0, RooFit.RooConst(rrv_delta_c0.getVal()), RooFit.RooConst(rrv_delta_c0.getError()))
            #alpha_constrains.add(delta_c0_constrain)

            rrv_delta_c1 = RooRealVar("rrv_delta_c1_ErfPowExp%s_%s"%(label, self.channel), "rrv_delta_c1_ErfPowExp%s_%s"%(label, self.channel),
                    tmp_rrv_c1_sr.getVal()-rrv_c1_sb.getVal(),
                    tmp_rrv_c1_sr.getVal()-rrv_c1_sb.getVal()-rrv_c1_sb.getError()-tmp_rrv_c1_sr.getError(),
                    tmp_rrv_c1_sr.getVal()-rrv_c1_sb.getVal()+rrv_c1_sb.getError()+tmp_rrv_c1_sr.getError())
            rrv_delta_c1.setError(TMath.Sqrt(rrv_c1_sb.getError()*rrv_c1_sb.getError() + tmp_rrv_c1_sr.getError()*tmp_rrv_c1_sr.getError()))
            #delta_c1_constrain = RooGaussian ("delta_c1_constrain", "delta_c1_constrain", rrv_delta_c1, RooFit.RooConst(rrv_delta_c1.getVal()), RooFit.RooConst(rrv_delta_c1.getError()))
            #alpha_constrains.add(delta_c1_constrain)


            rrv_delta_offset = RooRealVar("rrv_delta_offset_ErfPowExp%s_%s"%(label, self.channel), "rrv_delta_offset_ErfPowExp%s_%s"%(label, self.channel),
                    tmp_rrv_offset_sr.getVal()-rrv_offset_sb.getVal(),
                    tmp_rrv_offset_sr.getVal()-rrv_offset_sb.getVal()-rrv_offset_sb.getError()-tmp_rrv_offset_sr.getError(),
                    tmp_rrv_offset_sr.getVal()-rrv_offset_sb.getVal()+rrv_offset_sb.getError()+tmp_rrv_offset_sr.getError())
            rrv_delta_offset.setError(TMath.Sqrt(rrv_offset_sb.getError()*rrv_offset_sb.getError() + tmp_rrv_offset_sr.getError()*tmp_rrv_offset_sr.getError()))
            #delta_offset_constrain = RooGaussian ("delta_offset_constrain", "delta_offset_constrain", rrv_delta_offset, RooFit.RooConst(rrv_delta_offset.getVal()), RooFit.RooConst(rrv_delta_offset.getError()))
            #alpha_constrains.add(delta_offset_constrain)

            rrv_delta_width = RooRealVar("rrv_delta_width_ErfPowExp%s_%s"%(label, self.channel), "rrv_delta_width_ErfPowExp%s_%s"%(label, self.channel),
                    tmp_rrv_width_sr.getVal()-rrv_width_sb.getVal(),
                    tmp_rrv_width_sr.getVal()-rrv_width_sb.getVal()-rrv_width_sb.getError()-tmp_rrv_width_sr.getError(),
                    tmp_rrv_width_sr.getVal()-rrv_width_sb.getVal()+rrv_width_sb.getError()+tmp_rrv_width_sr.getError())
            rrv_delta_width.setError(TMath.Sqrt(rrv_width_sb.getError()*rrv_width_sb.getError() + tmp_rrv_width_sr.getError()*tmp_rrv_width_sr.getError()))
            #delta_width_constrain = RooGaussian ("delta_width_constrain", "delta_width_constrain", rrv_delta_width, RooFit.RooConst(rrv_delta_width.getVal()), RooFit.RooConst(rrv_delta_width.getError()))
            #alpha_constrains.add(delta_width_constrain)

            rrv_c0_sr     = RooFormulaVar("rrv_c0_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_c0_sb, rrv_delta_c0))
            rrv_c1_sr     = RooFormulaVar("rrv_c1_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_c1_sb, rrv_delta_c1))
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_offset_sb, rrv_delta_offset))
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_width_sb, rrv_delta_width))

            correct_factor_pdf = RooAlpha4ErfPowExpPdf("correct_factor_pdf", "correct_factor_pdf", rrv_x, rrv_c0_sr, rrv_c1_sr, rrv_offset_sr, rrv_width_sr, rrv_c0_sb, rrv_c1_sb, rrv_offset_sb, rrv_width_sb)

        if mlvj_model == "Exp":
            rrv_c_sb    = self.workspace4fit_.var("rrv_c_Exp%s_sb_lo_%s"%(label, self.channel))
            c_sb_constrain = RooGaussian ("c_sb_constrain", "c_sb_constrain", rrv_c_sb, RooFit.RooConst(rrv_c_sb.getVal()), RooFit.RooConst(rrv_c_sb.getError()))
            alpha_constrains.add(c_sb_constrain)

            rrv_delta_c = RooRealVar("rrv_delta_c_Exp%s_%s"%(label, self.channel), "rrv_delta_c_Exp%s_%s"%(label, self.channel),
                    self.workspace4fit_.var("rrv_c_Exp%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c_sb.getVal(),
                    self.workspace4fit_.var("rrv_c_Exp%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c_sb.getVal()-4*rrv_c_sb.getError(),
                    self.workspace4fit_.var("rrv_c_Exp%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c_sb.getVal()+4*rrv_c_sb.getError())

            correct_factor_pdf = RooExponential("correct_factor_pdf", "correct_factor_pdf", rrv_x, rrv_delta_c)

        if mlvj_model == "2Exp":
            rrv_c0_sb    = self.workspace4fit_.var("rrv_c0_2Exp%s_sb_lo_%s"%(label, self.channel))
            rrv_delta_c0 = RooRealVar("rrv_delta_c0_2Exp%s_%s"%(label, self.channel), "rrv_delta_c0_2Exp%s_%s"%(label, self.channel),
                    self.workspace4fit_.var("rrv_c0_2Exp%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c0_sb.getVal(),
                    self.workspace4fit_.var("rrv_c0_2Exp%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c0_sb.getVal()-4*rrv_c0_sb.getError(),
                    self.workspace4fit_.var("rrv_c0_2Exp%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c0_sb.getVal()+4*rrv_c0_sb.getError())
            rrv_c0_sr = RooFormulaVar("rrv_c0_2Exp_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_c0_sb, rrv_delta_c0))

            rrv_c1_sb = self.workspace4fit_.var("rrv_c1_2Exp%s_sb_lo_%s"%(label, self.channel))
            rrv_delta_c1 = RooRealVar("rrv_delta_c1_2Exp%s_%s"%(label, self.channel), "rrv_delta_c1_2Exp%s_%s"%(label, self.channel),
                    self.workspace4fit_.var("rrv_c1_2Exp%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c1_sb.getVal(),
                    self.workspace4fit_.var("rrv_c1_2Exp%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c1_sb.getVal()-4*rrv_c1_sb.getError(),
                    self.workspace4fit_.var("rrv_c1_2Exp%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c1_sb.getVal()+4*rrv_c1_sb.getError())
            rrv_c1_sr = RooFormulaVar("rrv_c1_2Exp_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_c1_sb, rrv_delta_c1))

            rrv_frac_sb    = self.workspace4fit_.var("rrv_frac_2Exp%s_sb_lo_%s"%(label, self.channel))
            rrv_delta_frac = RooRealVar("rrv_delta_frac_2Exp%s_%s"%(label, self.channel), "rrv_delta_frac_2Exp%s_%s"%(label, self.channel),
                    self.workspace4fit_.var("rrv_frac_2Exp%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_frac_sb.getVal(),
                    self.workspace4fit_.var("rrv_frac_2Exp%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_frac_sb.getVal()-4*rrv_frac_sb.getError(),
                    self.workspace4fit_.var("rrv_frac_2Exp%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_frac_sb.getVal()+4*rrv_frac_sb.getError())
            rrv_frac_sr = RooFormulaVar("rrv_frac_2Exp_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_frac_sb, rrv_delta_frac))

            correct_factor_pdf = RooAlpha42ExpPdf("correct_factor_pdf", "correct_factor_pdf", rrv_x, rrv_c0_sr, rrv_c1_sr, rrv_frac_sr, rrv_c0_sb, rrv_c1_sb, rrv_frac_sb)

        if mlvj_model == "Pow":

            rrv_c_sb    = self.workspace4fit_.var("rrv_c_Pow%s_sb_lo_%s"%(label, self.channel))
            c_sb_constrain = RooGaussian ("c_sb_constrain", "c_sb_constrain", rrv_c_sb, RooFit.RooConst(rrv_c_sb.getVal()), RooFit.RooConst(rrv_c_sb.getError()))
            alpha_constrains.add(c_sb_constrain)

            rrv_delta_c = RooRealVar("rrv_delta_c_Pow%s_%s"%(label, self.channel), "rrv_delta_c_Pow%s_%s"%(label, self.channel), 0., -100*rrv_c_sb.getError(), 100*rrv_c_sb.getError())
            correct_factor_pdf = RooPowPdf("correct_factor_pdf", "correct_factor_pdf", rrv_x, rrv_delta_c)

        if mlvj_model == "ExpN":
            rrv_c_sb  = self.workspace4fit_.var("rrv_c_ExpN%s_sb_lo_%s"%(label, self.channel))
            rrv_n_sb  = self.workspace4fit_.var("rrv_n_ExpN%s_sb_lo_%s"%(label, self.channel))
            rrv_delta_c = RooRealVar("rrv_delta_c_ExpN%s_%s"%(label, self.channel), "rrv_delta_c_ExpN%s_%s"%(label, self.channel),
                    self.workspace4fit_.var("rrv_c_ExpN%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c_sb.getVal(),
                    self.workspace4fit_.var("rrv_c_ExpN%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c_sb.getVal()-4*rrv_c_sb.getError(),
                    self.workspace4fit_.var("rrv_c_ExpN%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_c_sb.getVal()+4*rrv_c_sb.getError())
            rrv_delta_n = RooRealVar("rrv_delta_n_ExpN%s_%s"%(label, self.channel), "rrv_delta_n_ExpN%s_%s"%(label, self.channel),
                    self.workspace4fit_.var("rrv_n_ExpN%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_n_sb.getVal(),
                    self.workspace4fit_.var("rrv_n_ExpN%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_n_sb.getVal()-4*rrv_n_sb.getError(),
                    self.workspace4fit_.var("rrv_n_ExpN%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_n_sb.getVal()+4*rrv_n_sb.getError())

            c_sb_constrain = RooGaussian ("c_sb_constrain", "c_sb_constrain", rrv_c_sb, RooFit.RooConst(rrv_c_sb.getVal()), RooFit.RooConst(rrv_c_sb.getError()))
            n_sb_constrain = RooGaussian ("n_sb_constrain", "n_sb_constrain", rrv_n_sb, RooFit.RooConst(rrv_n_sb.getVal()), RooFit.RooConst(rrv_n_sb.getError()))
            alpha_constrains.add(c_sb_constrain)
            alpha_constrains.add(n_sb_constrain)

            correct_factor_pdf = RooExpNPdf("correct_factor_pdf", "correct_factor_pdf", rrv_x, rrv_delta_c, rrv_delta_n)

        if mlvj_model == "ExpTail":
            rrv_s_sb = self.workspace4fit_.var("rrv_s_ExpTail%s_sb_lo_%s"%(label, self.channel))
            rrv_a_sb = self.workspace4fit_.var("rrv_a_ExpTail%s_sb_lo_%s"%(label, self.channel))

            rrv_delta_s = RooRealVar("rrv_delta_s_ExpTail%s_%s"%(label, self.channel), "rrv_delta_s_ExpTail%s_%s"%(label, self.channel),
                    self.workspace4fit_.var("rrv_s_ExpTail%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_s_sb.getVal(),
                    self.workspace4fit_.var("rrv_s_ExpTail%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_s_sb.getVal()-4*rrv_s_sb.getError(),
                    self.workspace4fit_.var("rrv_s_ExpTail%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_s_sb.getVal()+4*rrv_s_sb.getError())
            rrv_delta_a = RooRealVar("rrv_delta_a_ExpTail%s_%s"%(label, self.channel), "rrv_delta_a_ExpTail%s_%s"%(label, self.channel),
                    self.workspace4fit_.var("rrv_a_ExpTail%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_a_sb.getVal(),
                    self.workspace4fit_.var("rrv_a_ExpTail%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_a_sb.getVal()-4*rrv_a_sb.getError(),
                    self.workspace4fit_.var("rrv_a_ExpTail%s_signal_region_%s"%(label, self.channel)).getVal()-rrv_a_sb.getVal()+4*rrv_a_sb.getError())
            rrv_a_sr = RooFormulaVar("rrv_a_ExpTail_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_a_sb, rrv_delta_a))
            rrv_s_sr = RooFormulaVar("rrv_s_ExpTail_sr%s_%s"%(label, self.channel), "@0+@1", RooArgList(rrv_s_sb, rrv_delta_s))

            s_sb_constrain = RooGaussian ("s_sb_constrain", "s_sb_constrain", rrv_s_sb, RooFit.RooConst(rrv_s_sb.getVal()), RooFit.RooConst(rrv_s_sb.getError()))
            a_sb_constrain = RooGaussian ("a_sb_constrain", "a_sb_constrain", rrv_a_sb, RooFit.RooConst(rrv_a_sb.getVal()), RooFit.RooConst(rrv_a_sb.getError()))
            alpha_constrains.add(s_sb_constrain)
            alpha_constrains.add(a_sb_constrain)

            correct_factor_pdf = RooAlpha4ExpTailPdf("correct_factor_pdf", "correct_factor_pdf", rrv_x, rrv_s_sr, rrv_a_sr, rrv_s_sb, rrv_a_sb)

        if mlvj_model == "Pow2":

            rrv_c0_sb    = self.workspace4fit_.var("rrv_c0_Pow2%s_sb_lo_%s"%(label, self.channel))
            rrv_c1_sb    = self.workspace4fit_.var("rrv_c1_Pow2%s_sb_lo_%s"%(label, self.channel))
            rrv_delta_c0 = RooRealVar("rrv_delta_c0_Pow2%s_%s"%(label, self.channel), "rrv_delta_c0_Pow2%s_%s"%(label, self.channel), 0., -100*rrv_c0_sb.getError(), 100*rrv_c0_sb.getError())
            rrv_delta_c1 = RooRealVar("rrv_delta_c1_Pow2%s_%s"%(label, self.channel), "rrv_delta_c1_Pow2%s_%s"%(label, self.channel), 0., -100*rrv_c1_sb.getError(), 100*rrv_c1_sb.getError())
            correct_factor_pdf = RooPow2Pdf("correct_factor_pdf", "correct_factor_pdf", rrv_x, rrv_delta_c0, rrv_delta_c1)

        ### define the category and do the simultaneous fit taking the combined dataset of events in mlvj sb and sr
        combData4fit = self.workspace4fit_.data("combData4fit%s_%s"%(label, self.channel))
        combData4fit.Print("v")
        correct_factor_pdf.Print("v")

        model_pdf_sb_lo_WJets         = self.workspace4fit_.pdf("model_pdf%s_sb_lo_%s_mlvj"%(label, self.channel))
        model_pdf_signal_region_WJets = RooEffProd("model_pdf%s_signal_region_%s_mlvj"%(label, self.channel), "model_pdf%s_signal_region_%s_mlvj"%(label, self.channel) , model_pdf_sb_lo_WJets, correct_factor_pdf)

        simPdf = RooSimultaneous("simPdf", "simPdf", self.data_category)
        simPdf.addPdf(model_pdf_sb_lo_WJets, "sideband")
        simPdf.addPdf(model_pdf_signal_region_WJets, "signal_region")

        simPdf.Print("v")
        model_pdf_sb_lo_WJets.getParameters(rdataset_WJets_sb_lo_mlvj).Print("v")
        model_pdf_signal_region_WJets.getParameters(rdataset_WJets_signal_region_mlvj).Print("v")

        #constrainslist is necessary
        #without the constrainslist, the parameters' error of sb function are much large than before and the alpha uncertainty band is much large
        #but the difference is very small when use 1sigma constrain or 3 sigma constrain for the parameters
        rfresult = simPdf.fitTo(combData4fit, RooFit.Save(kTRUE), RooFit.Extended(kFALSE), RooFit.SumW2Error(kTRUE),
                RooFit.Minimizer("Minuit"), RooFit.ExternalConstraints(alpha_constrains))
        rfresult = simPdf.fitTo(combData4fit, RooFit.Save(kTRUE), RooFit.Extended(kFALSE), RooFit.SumW2Error(kTRUE),
                RooFit.Minimizer("Minuit2"), RooFit.ExternalConstraints(alpha_constrains))
        rfresult.Print("v")
        rfresult.covarianceMatrix().Print("v")

        #### Decorrelate the parameters in the alpha shape
        wsfit_tmp = RooWorkspace("wsfit_tmp%s_sim_mlvj"%(label))
        print "############### diagonalizer alpha "
        Deco      = PdfDiagonalizer("Deco%s_sim_%s_%s_mlvj"%(label, self.channel, self.wtagger_category), wsfit_tmp, rfresult)
        correct_factor_pdf_deco = Deco.diagonalize(correct_factor_pdf)
        correct_factor_pdf_deco.Print()
        correct_factor_pdf_deco.getParameters(rdataset_WJets_signal_region_mlvj).Print("v")
        getattr(self.workspace4fit_, "import")(correct_factor_pdf_deco)

        ## in case of default Wjets with default shape
        if TString(label).Contains("_WJets0"):

            ### only mc plots in the SB region
            mplot_sb_lo = rrv_x.frame(RooFit.Title("WJets sb low"), RooFit.Bins(int(rrv_x.getBins()/self.binwidth_narrow_factor)))

            #rdataset_WJets_sb_lo_mlvj.plotOn(mplot_sb_lo, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0))
            #model_pdf_sb_lo_WJets.plotOn(mplot_sb_lo)
            combData4fit.plotOn(mplot_sb_lo, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0), RooFit.Cut("data_category == data_category::sideband"))
            simPdf.plotOn(mplot_sb_lo, RooFit.Slice(self.data_category, "sideband"), RooFit.ProjWData(RooArgSet(self.data_category), combData4fit))
            mplot_pull_sideband = self.get_pull(rrv_x, mplot_sb_lo)
            parameters_list     = model_pdf_sb_lo_WJets.getParameters(rdataset_WJets_sb_lo_mlvj)
            mplot_sb_lo.GetYaxis().SetRangeUser(1e-2, mplot_sb_lo.GetMaximum()*1.5)

            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize()
            nBinX = mplot_sb_lo.GetNbinsX()
            ndof  = nBinX-self.nPar_float_in_fitTo
            print mplot_sb_lo.chiSquare()
            print "#################### JENchi2 nPar = %s, chiSquare = %s/%s"%(self.nPar_float_in_fitTo , mplot_sb_lo.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof)
            datahist = rdataset_WJets_sb_lo_mlvj.binnedClone(rdataset_WJets_sb_lo_mlvj.GetName()+"_binnedClone", rdataset_WJets_sb_lo_mlvj.GetName()+"_binnedClone")

            self.draw_canvas_with_pull1(rrv_x, datahist, mplot_sb_lo, mplot_pull_sideband, ndof, parameters_list, "%s/other/"%(self.plotsDir), "m_lvj%s_sb_lo_sim"%(label), "", 1, 1)

            ### only mc plots in the SR region
            mplot_signal_region = rrv_x.frame(RooFit.Title("WJets sr"), RooFit.Bins(int(rrv_x.getBins()/self.binwidth_narrow_factor)))

            #rdataset_WJets_signal_region_mlvj.plotOn(mplot_signal_region, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0))
            #model_pdf_signal_region_WJets.plotOn(mplot_signal_region)
            combData4fit.plotOn(mplot_signal_region, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0), RooFit.Cut("data_category == data_category::signal_region"))
            simPdf.plotOn(mplot_signal_region, RooFit.Slice(self.data_category, "signal_region"), RooFit.ProjWData(RooArgSet(self.data_category), combData4fit))

            mplot_pull_signal_region = self.get_pull(rrv_x, mplot_signal_region)
            parameters_list = model_pdf_signal_region_WJets.getParameters(rdataset_WJets_signal_region_mlvj)
            mplot_signal_region.GetYaxis().SetRangeUser(1e-2, mplot_signal_region.GetMaximum()*1.5)

            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize()
            nBinX = mplot_signal_region.GetNbinsX()
            ndof  = nBinX-self.nPar_float_in_fitTo
            print mplot_signal_region.chiSquare()
            print "#################### JENchi2 nPar = %s, chiSquare = %s/%s"%(self.nPar_float_in_fitTo , mplot_signal_region.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof)

            self.draw_canvas_with_pull1(rrv_x, datahist, mplot_signal_region, mplot_pull_signal_region, ndof, parameters_list, "%s/other/"%(self.plotsDir), "m_lvj%s_signal_region_sim"%(label), "", 1, 1)

        ### Total plot shape in sb_lo, sr and alpha
        model_pdf_sb_lo_WJets.plotOn(mplot, RooFit.Name("Sideband"), RooFit.LineStyle(10))
        model_pdf_signal_region_WJets.plotOn(mplot, RooFit.LineColor(kRed) , RooFit.LineStyle(8), RooFit.Name("Signal Region"))
        correct_factor_pdf_deco.plotOn(mplot, RooFit.LineColor(kBlack), RooFit.Name("#alpha"))

        ### plot also what is get from other source if available : alternate PS and shape: 1 PS and 01 is shape or fitting function
        if TString(label).Contains("_WJets0"):
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_xww_sim_%s_%s_mlvj"%(self.channel, self.wtagger_category)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_xww_sim_%s_%s_mlvj"%(self.channel, self.wtagger_category)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3), RooFit.Name("#alpha: Alternate PS"))

            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_xww_sim_%s_%s_mlvj"%(self.channel, self.wtagger_category)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_xww_sim_%s_%s_mlvj"%(self.channel, self.wtagger_category)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7), RooFit.Name("#alpha: Alternate Function"))

        paras = RooArgList()
        ### Make a list of paramters as a function of the model after decorrelation
        if mlvj_model == "ErfExp_v1" or mlvj_model == "ErfPow_v1" or mlvj_model == "2Exp" :
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig0"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig1"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig2"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig3"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig4"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig5"%(label, self.channel, self.wtagger_category)))

        if mlvj_model == "ErfPow2_v1" or mlvj_model == "ErfPowExp_v1" :
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig0"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig1"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig2"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig3"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig4"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig5"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig6"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig7"%(label, self.channel, self.wtagger_category)))

        if mlvj_model == "Exp" or mlvj_model == "Pow":
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig0"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig1"%(label, self.channel, self.wtagger_category)))

        if mlvj_model == "ExpN" or mlvj_model == "ExpTail" or mlvj_model == "Pow2":
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig0"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig1"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig2"%(label, self.channel, self.wtagger_category)))
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig3"%(label, self.channel, self.wtagger_category)))

        if TString(label).Contains("_WJets0") or TString(label).Contains("_WJets1"): ### draw error band ar 1 and 2 sigma using the decorrelated shape
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj"%(label, self.channel, self.wtagger_category), "rrv_mass_lvj", paras, self.workspace4fit_, 1 , mplot, kGray+3, "F", 3001, "#alpha #pm", 20, 400)
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj"%(label, self.channel, self.wtagger_category), "rrv_mass_lvj", paras, self.workspace4fit_, 2 , mplot, kGreen+2, "F", 3002, "#alpha #pm", 20, 400)
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj"%(label, self.channel, self.wtagger_category), "rrv_mass_lvj", paras, self.workspace4fit_, 1 , mplot, kGray+3, "F", 3001, "#alpha_invisible #pm", 20, 400)

        ### plot on the same canvas
        correct_factor_pdf_deco.plotOn(mplot, RooFit.LineColor(kBlack), RooFit.Name("#alpha_invisible"))

        if TString(label).Contains("_WJets0") : ## add also the plot of alternate ps and function on the canvas
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_xww_sim_%s_%s_mlvj"%(self.channel, self.wtagger_category)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_xww_sim_%s_%s_mlvj"%(self.channel, self.wtagger_category)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3), RooFit.Name("#alpha_invisible: Alternate PS"))
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_xww_sim_%s_%s_mlvj"%(self.channel, self.wtagger_category)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_xww_sim_%s_%s_mlvj"%(self.channel, self.wtagger_category)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7), RooFit.Name("#alpha_invisible: Alternate Function"))

        elif TString(label).Contains("_WJets01") : ## add also the plot of alternate ps and function on the canvas
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj"%(self.channel, self.wtagger_category)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj"%(self.channel, self.wtagger_category)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3), RooFit.Name("#alpha_invisible: Alternate PS"))
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets0_sim_%s_%s_mlvj"%(self.channel, self.wtagger_category)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets0_sim_%s_%s_mlvj"%(self.channel, self.wtagger_category)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7), RooFit.Name("#alpha_invisible: Alternate Function"))

        ### Add the legend
        #self.leg = self.legend4Plot(mplot, 1, 0, -0.01, -0.14, 0.01, -0.06, 0.)
        self.leg = self.legend4Plot(mplot, 1, 0, 0.13, 0.0, 0.0, 0.0, 0.)
        mplot.addObject(self.leg)

        ## set the Y axis in arbitrary unit
        ##if self.signal_sample == "ggH600" or self.signal_sample == "ggH700":
        ##    tmp_y_max = 0.25
        ##else: tmp_y_max = 0.28
        ##if self.signal_mass <1200:
        ##    tmp_y_max = 0.30
        ##else:
        ##    tmp_y_max = 0.40
        tmp_y_max = 0.25

        mplot.GetYaxis().SetRangeUser(1e-3, tmp_y_max)

        #### Draw another axis with the real value of alpha
        model_pdf_sb_lo_WJets.getVal(RooArgSet(rrv_x))
        model_pdf_signal_region_WJets.getVal(RooArgSet(rrv_x))
        correct_factor_pdf_deco.getVal(RooArgSet(rrv_x))
        tmp_alpha_ratio = (model_pdf_signal_region_WJets.getVal(RooArgSet(rrv_x))/model_pdf_sb_lo_WJets.getVal(RooArgSet(rrv_x)))
        tmp_alpha_pdf   = correct_factor_pdf_deco.getVal(RooArgSet(rrv_x)) * mplot.getFitRangeBinW() ## value of the pdf in each point
        tmp_alpha_scale = tmp_alpha_ratio/tmp_alpha_pdf

        #add alpha scale axis
        axis_alpha = TGaxis(rrv_x.getMax(), 1e-3, rrv_x.getMax(), tmp_y_max, 1e-3*tmp_alpha_scale, tmp_y_max*tmp_alpha_scale, 510, "+L")
        axis_alpha.SetTitle("#alpha")
        axis_alpha.SetTitleOffset(0.65)
        axis_alpha.SetTitleSize(0.05)
        axis_alpha.SetLabelSize(0.045)
        axis_alpha.SetTitleFont(42)
        axis_alpha.SetLabelFont(42)
        #axis_alpha.RotateTitle(1)
        mplot.addObject(axis_alpha)

        #self.draw_canvas(mplot, "plots_%s_%s_%s/other/"%(options.additioninformation, self.channel, self.wtagger_category), "correction_pdf%s_%s_M_lvj_signal_region_to_sideband"%(label, mlvj_model), 0, 1)
        self.draw_canvas(mplot, self.plotsDir+"/other/", "correction_pdf%s_%s_M_lvj_signal_region_to_sideband"%(label, mlvj_model), 0, 1, 0, 1)

        correct_factor_pdf_deco.getParameters(rdataset_WJets_sb_lo_mlvj).Print("v")
        model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco = self.workspace4fit_.pdf("model_pdf%s_sb_lo_from_fitting_%s_mlvj_Deco%s_sb_lo_from_fitting_%s_%s_mlvj"%(label, self.channel, label, self.channel, self.wtagger_category))
        model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco.Print("v")

        ### Wjets shape in the SR correctedfunction * sb
        model_pdf_WJets_signal_region_after_correct_mlvj = RooProdPdf("model_pdf%s_signal_region_%s_after_correct_mlvj"%(label, self.channel), "model_pdf%s_signal_region_%s_after_correct_mlvj"%(label, self.channel), model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco, self.workspace4fit_.pdf("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj"%(label, self.channel, self.wtagger_category)))
        model_pdf_WJets_signal_region_after_correct_mlvj.Print()
        ### fix the parmaters and import in the workspace
        getattr(self.workspace4fit_, "import")(model_pdf_WJets_signal_region_after_correct_mlvj)

        ##### calculate the normalization and alpha for limit datacard
        self.workspace4fit_.var("rrv_number%s_signal_region_%s_mlvj"%(label, self.channel)).Print()
        self.workspace4fit_.var("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label, self.channel)).Print()
        self.workspace4fit_.var("rrv_number%s_signal_region_%s_mlvj"%(label, self.channel)).setVal(self.workspace4fit_.var("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label, self.channel)).getVal())
        self.workspace4fit_.var("rrv_number%s_signal_region_%s_mlvj"%(label, self.channel)).setError(self.workspace4fit_.var("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label, self.channel)).getError())

        self.workspace4fit_.var("rrv_number%s_signal_region_%s_mlvj"%(label, self.channel)).setConstant(kTRUE)


    ##### Method used to cycle on the events and for the dataset to be fitted
    def get_mj_and_mlvj_dataset(self, in_file_name, label, jet_mass):# to get the shape of m_lvj, jet_mass = "jet_mass_pr"

        print "################### get_mj_and_mlvj_dataset : ", in_file_name, "  ", label, "  ##################"

        fileIn_name = TString(options.inPath+"/"+self.file_Directory+in_file_name)
        #print fileIn_name
        fileIn = TFile(fileIn_name.Data())
        treeIn = fileIn.Get("PKUTree")

        rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j")
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rrv_weight   = RooRealVar("rrv_weight", "rrv_weight", 0. , 10000000.)

        ##### dataset of m_j -> scaleed and not scaled to lumi
        rdataset_mj     = RooDataSet("rdataset"+label+"_"+self.channel+"_mj", "rdataset"+label+"_"+self.channel+"_mj", RooArgSet(rrv_mass_j, rrv_weight), RooFit.WeightVar(rrv_weight))
        rdataset4fit_mj = RooDataSet("rdataset4fit"+label+"_"+self.channel+"_mj", "rdataset4fit"+label+"_"+self.channel+"_mj", RooArgSet(rrv_mass_j, rrv_weight), RooFit.WeightVar(rrv_weight))
        ##### dataset of m_lvj -> scaled and not scaled to lumi in different region
        rdataset_sb_lo_mlvj = RooDataSet("rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj", "rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj", RooArgSet(rrv_mass_lvj, rrv_weight), RooFit.WeightVar(rrv_weight))

        rdataset_signal_region_mlvj = RooDataSet("rdataset"+label+"_signal_region"+"_"+self.channel+"_mlvj", "rdataset"+label+"_signal_region"+"_"+self.channel+"_mlvj", RooArgSet(rrv_mass_lvj, rrv_weight), RooFit.WeightVar(rrv_weight))

        rdataset_sb_hi_mlvj = RooDataSet("rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj", "rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj", RooArgSet(rrv_mass_lvj, rrv_weight), RooFit.WeightVar(rrv_weight))

        rdataset4fit_sb_lo_mlvj = RooDataSet("rdataset4fit"+label+"_sb_lo"+"_"+self.channel+"_mlvj", "rdataset4fit"+label+"_sb_lo"+"_"+self.channel+"_mlvj", RooArgSet(rrv_mass_lvj, rrv_weight), RooFit.WeightVar(rrv_weight))

        rdataset4fit_signal_region_mlvj = RooDataSet("rdataset4fit"+label+"_signal_region"+"_"+self.channel+"_mlvj", "rdataset4fit"+label+"_signal_region"+"_"+self.channel+"_mlvj", RooArgSet(rrv_mass_lvj, rrv_weight), RooFit.WeightVar(rrv_weight))

        rdataset4fit_sb_hi_mlvj = RooDataSet("rdataset4fit"+label+"_sb_hi"+"_"+self.channel+"_mlvj", "rdataset4fit"+label+"_sb_hi"+"_"+self.channel+"_mlvj", RooArgSet(rrv_mass_lvj, rrv_weight), RooFit.WeightVar(rrv_weight))

        ### categorize the event in sideband and signal region --> combined dataset
        self.data_category = RooCategory("data_category", "data_category")
        self.data_category.defineType("sideband")
        self.data_category.defineType("signal_region")

        combData = RooDataSet("combData"+label+"_"+self.channel, "combData"+label+"_"+self.channel, RooArgSet(rrv_mass_lvj, self.data_category, rrv_weight), RooFit.WeightVar(rrv_weight))
        combData4fit = RooDataSet("combData4fit"+label+"_"+self.channel, "combData4fit"+label+"_"+self.channel, RooArgSet(rrv_mass_lvj, self.data_category, rrv_weight), RooFit.WeightVar(rrv_weight))

        print "###### N entries: ", treeIn.GetEntries()
        hnum_4region = TH1D("hnum_4region"+label+"_"+self.channel, "hnum_4region"+label+"_"+self.channel, 4,-1.5, 2.5)# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total
        hnum_2region = TH1D("hnum_2region"+label+"_"+self.channel, "hnum_2region"+label+"_"+self.channel, 2,-0.5, 1.5)# m_lvj 0: signal_region; 1: total

        tmp_lumi = self.GetLumi()
        for i in range(treeIn.GetEntries()):
            if i % 100000 == 0:
                print "iEntry: ", i
            treeIn.GetEntry(i)

            if i == 0:
                tmp_scale_to_lumi = treeIn.lumiWeight*tmp_lumi

            tmp_jet_mass = getattr(treeIn, jet_mass)

            self.isGoodEvent = 0
            ## event in the whole range
            if self.channel == "mu" or self.channel == "el":
                ##if treeIn.CategoryID == self.categoryID and treeIn.m_lvj> rrv_mass_lvj.getMin() and treeIn.m_lvj<rrv_mass_lvj.getMax() and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() :
                ## temporary using tau21 cut to replace categoryID
                if TMath.Abs(treeIn.CategoryID)< 3 and treeIn.tau21<=self.tau21_cut and treeIn.m_lvj> rrv_mass_lvj.getMin() and treeIn.m_lvj<rrv_mass_lvj.getMax() and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and treeIn.passFilter_HBHE>0 and treeIn.passFilter_GlobalHalo>0 and treeIn.passFilter_HBHEIso>0 and treeIn.passFilter_ECALDeadCell>0 and treeIn.passFilter_GoodVtx>0 and treeIn.passFilter_EEBadSc>0 and treeIn.passFilter_badMuon>0 and treeIn.passFilter_badChargedHadron>0:
                    self.isGoodEvent = 1

            if self.channel == "mu" and treeIn.mtVlepnew<50:
                self.isGoodEvent = 0

            if label == "_data" or label == "_data_xww" :
                if tmp_jet_mass>105 and tmp_jet_mass<135:
                    #keep blind for VH analysis
                    self.isGoodEvent = 0
                if options.keepblind == 1 and tmp_jet_mass>65 and tmp_jet_mass<135:
                    self.isGoodEvent = 0

                if self.channel == "mu" and treeIn.l_pt> 500 and TMath.Abs(treeIn.l_eta)>1.2:
                    self.isGoodEvent=0


            if self.isGoodEvent == 1:
                ## HLT weight for mu channel
                ##tmp_HLT_weight = 1.0
                ##if self.channel == "mu":
                ##    if   treeIn.l_pt<60 and treeIn.l_pt>=50 and TMath.Abs(treeIn.l_eta)< 0.9 and TMath.Abs(treeIn.l_eta)>= 0  : tmp_HLT_weight= 0.975
                ##    elif treeIn.l_pt<60 and treeIn.l_pt>=50 and TMath.Abs(treeIn.l_eta)< 1.2 and TMath.Abs(treeIn.l_eta)>= 0.9: tmp_HLT_weight= 0.921
                ##    elif treeIn.l_pt<60 and treeIn.l_pt>=50 and TMath.Abs(treeIn.l_eta)< 2.1 and TMath.Abs(treeIn.l_eta)>= 1.2: tmp_HLT_weight= 0.886
                ##    elif treeIn.l_pt<80 and treeIn.l_pt>=60 and TMath.Abs(treeIn.l_eta)< 0.9 and TMath.Abs(treeIn.l_eta)>= 0  : tmp_HLT_weight= 0.972
                ##    elif treeIn.l_pt<80 and treeIn.l_pt>=60 and TMath.Abs(treeIn.l_eta)< 1.2 and TMath.Abs(treeIn.l_eta)>= 0.9: tmp_HLT_weight= 0.935
                ##    elif treeIn.l_pt<80 and treeIn.l_pt>=60 and TMath.Abs(treeIn.l_eta)< 2.1 and TMath.Abs(treeIn.l_eta)>= 1.2: tmp_HLT_weight= 0.883
                ##    elif                    treeIn.l_pt>=80 and TMath.Abs(treeIn.l_eta)< 0.9 and TMath.Abs(treeIn.l_eta)>= 0  : tmp_HLT_weight= 0.966
                ##    elif                    treeIn.l_pt>=80 and TMath.Abs(treeIn.l_eta)< 1.2 and TMath.Abs(treeIn.l_eta)>= 0.9: tmp_HLT_weight= 0.902
                ##    elif                    treeIn.l_pt>=80 and TMath.Abs(treeIn.l_eta)< 2.1 and TMath.Abs(treeIn.l_eta)>= 1.2: tmp_HLT_weight= 0.850
                ##elif self.channel == "el":
                ##    tmp_HLT_weight = 0.92
                tmp_HLT_weight = 1.0
                #print "l_pt = %s, l_eta = %s, HLT = %s"%(treeIn.l_pt, treeIn.l_eta, tmp_HLT_weight)
                #raw_input("zixu")
                ### weigh MC events
                tmp_event_weight     = treeIn.weight*tmp_HLT_weight*tmp_lumi
                tmp_event_weight4fit = (treeIn.weight)*tmp_HLT_weight*tmp_lumi/tmp_scale_to_lumi
                #tmp_event_weight4fit = tmp_event_weight
                ###### wtagger_eff_reweight
                if label == "_data" or label == "_data_xww" :
                    tmp_event_weight = 1.
                    tmp_event_weight4fit = 1.

                else:
                    if TString(label).Contains("_TTbar") or TString(label).Contains("_STop") :
                        tmp_event_weight = tmp_event_weight*self.rrv_wtagger_eff_reweight_forT.getVal()
                    else:
                        tmp_event_weight = tmp_event_weight*self.rrv_wtagger_eff_reweight_forV.getVal()
                    ###scale factor from control plots
                    ##if TString(label).Contains("_TTbar"):
                    ##    tmp_event_weight = tmp_event_weight*self.controlplot_TTbar_scale
                    ##if TString(label).Contains("_WJets"):
                    ##    tmp_event_weight = tmp_event_weight*self.controlplot_WJets_scale

                rrv_mass_lvj.setVal(treeIn.m_lvj)

                if tmp_jet_mass >= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max:
                    rdataset_sb_lo_mlvj.add(RooArgSet(rrv_mass_lvj), tmp_event_weight)
                    rdataset4fit_sb_lo_mlvj.add(RooArgSet(rrv_mass_lvj), tmp_event_weight4fit)

                    self.data_category.setLabel("sideband")
                    combData.add(RooArgSet(rrv_mass_lvj, self.data_category), tmp_event_weight)
                    combData4fit.add(RooArgSet(rrv_mass_lvj, self.data_category), tmp_event_weight4fit)

                if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max:
                    rdataset_signal_region_mlvj.add(RooArgSet(rrv_mass_lvj), tmp_event_weight)
                    rdataset4fit_signal_region_mlvj.add(RooArgSet(rrv_mass_lvj), tmp_event_weight4fit)

                    self.data_category.setLabel("signal_region")
                    combData.add(RooArgSet(rrv_mass_lvj, self.data_category), tmp_event_weight)
                    combData4fit.add(RooArgSet(rrv_mass_lvj, self.data_category), tmp_event_weight4fit)
                    hnum_2region.Fill(1, tmp_event_weight)

                    if treeIn.m_lvj >=self.mlvj_signal_min and treeIn.m_lvj <self.mlvj_signal_max:
                        hnum_2region.Fill(0, tmp_event_weight)

                if tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max:
                    rdataset_sb_hi_mlvj.add(RooArgSet(rrv_mass_lvj), tmp_event_weight)
                    rdataset4fit_sb_hi_mlvj.add(RooArgSet(rrv_mass_lvj), tmp_event_weight4fit)

                rrv_mass_j.setVal(tmp_jet_mass)
                rdataset_mj.add(RooArgSet(rrv_mass_j), tmp_event_weight)
                rdataset4fit_mj.add(RooArgSet(rrv_mass_j), tmp_event_weight4fit)

                if tmp_jet_mass >=self.mj_sideband_lo_min and tmp_jet_mass <self.mj_sideband_lo_max:
                    hnum_4region.Fill(-1, tmp_event_weight)
                if tmp_jet_mass >=self.mj_signal_min and tmp_jet_mass <self.mj_signal_max :
                    hnum_4region.Fill(0, tmp_event_weight)
                if tmp_jet_mass >=self.mj_sideband_hi_min and tmp_jet_mass <self.mj_sideband_hi_max:
                    hnum_4region.Fill(1, tmp_event_weight)

                hnum_4region.Fill(2, tmp_event_weight)

        if not label== "_data" and not label == "_data_xww": ## correct also because events in 4fit dataset where not rescaled in the cycle
            if TString(label).Contains("_TTbar") or TString(label).Contains("_STop") :
                tmp_scale_to_lumi = tmp_scale_to_lumi*self.rrv_wtagger_eff_reweight_forT.getVal()
            else:
                tmp_scale_to_lumi = tmp_scale_to_lumi*self.rrv_wtagger_eff_reweight_forV.getVal()
            ##if TString(label).Contains("_TTbar") :
            ##    tmp_scale_to_lumi = tmp_scale_to_lumi*self.controlplot_TTbar_scale
            ##if TString(label).Contains("_WJets") :
            ##    tmp_scale_to_lumi = tmp_scale_to_lumi*self.controlplot_WJets_scale

        ### scaler to lumi for MC in 4fit datasets
        rrv_scale_to_lumi = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel, "rrv_scale_to_lumi"+label+"_"+self.channel, tmp_scale_to_lumi)
        rrv_scale_to_lumi.Print()
        getattr(self.workspace4fit_, "import")(rrv_scale_to_lumi)

        ### prepare m_lvj dataset to be compared with the fit results
        rrv_number_dataset_signal_region_mlvj = RooRealVar("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj", "rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj", hnum_2region.GetBinContent(1))
        rrv_number_dataset_AllRange_mlvj = RooRealVar("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj", "rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj", hnum_2region.GetBinContent(2))

        getattr(self.workspace4fit_, "import")(rrv_number_dataset_signal_region_mlvj)
        getattr(self.workspace4fit_, "import")(rrv_number_dataset_AllRange_mlvj)

        ### import the dataser
        getattr(self.workspace4fit_, "import")(rdataset_sb_lo_mlvj)
        getattr(self.workspace4fit_, "import")(rdataset_signal_region_mlvj)
        getattr(self.workspace4fit_, "import")(rdataset_sb_hi_mlvj)
        getattr(self.workspace4fit_, "import")(rdataset4fit_sb_lo_mlvj)
        getattr(self.workspace4fit_, "import")(rdataset4fit_signal_region_mlvj)
        getattr(self.workspace4fit_, "import")(rdataset4fit_sb_hi_mlvj)
        getattr(self.workspace4fit_, "import")(combData)
        getattr(self.workspace4fit_, "import")(combData4fit)

        ### write in the output
        self.file_out.write("\n%s events number in m_lvj from dataset: %s"%(label, rdataset_signal_region_mlvj.sumEntries()))


        ### prepare m_j dataset
        rrv_number_dataset_sb_lo_mj = RooRealVar("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj", "rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj", hnum_4region.GetBinContent(1))
        rrv_number_dataset_signal_region_mj = RooRealVar("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj", "rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj", hnum_4region.GetBinContent(2))
        rrv_number_dataset_sb_hi_mj = RooRealVar("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj", "rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj", hnum_4region.GetBinContent(3))
        getattr(self.workspace4fit_, "import")(rrv_number_dataset_sb_lo_mj)
        getattr(self.workspace4fit_, "import")(rrv_number_dataset_signal_region_mj)
        getattr(self.workspace4fit_, "import")(rrv_number_dataset_sb_hi_mj)

        getattr(self.workspace4fit_, "import")(rdataset_mj)
        getattr(self.workspace4fit_, "import")(rdataset4fit_mj)

        #### print everything

        rdataset_sb_lo_mlvj.Print()
        rdataset_signal_region_mlvj.Print()
        rdataset_sb_hi_mlvj.Print()
        rdataset_mj.Print()
        rdataset4fit_sb_lo_mlvj.Print()
        rdataset4fit_signal_region_mlvj.Print()
        rdataset4fit_sb_hi_mlvj.Print()
        rdataset4fit_mj.Print()
        rrv_number_dataset_signal_region_mlvj.Print()
        rrv_number_dataset_AllRange_mlvj.Print()
        rrv_number_dataset_sb_lo_mj.Print()
        rrv_number_dataset_signal_region_mj.Print()
        rrv_number_dataset_sb_hi_mj.Print()
        print rdataset_signal_region_mlvj.sumEntries()
        print rrv_number_dataset_signal_region_mlvj.getVal()
        print rrv_number_dataset_AllRange_mlvj.getVal()

    ##### Prepare the workspace for the limit and to store info to be printed in the datacard
    def prepare_limit(self, mode, isTTbarFloating = 0, isVVFloating = 0, isSTopFloating = 0):
        print "####################### prepare_limit for %s method ####################"%(mode)

        getattr(self.workspace4limit_, "import")(self.workspace4fit_.var("rrv_mass_lvj"))

        ### whole number of events from the considered signal sample, WJets, VV, TTbar, STop -> couting
        getattr(self.workspace4limit_, "import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_%s_xww_%s_mlvj"%(self.signal_sample, self.channel)).clone("rate_%s_xww_for_counting"%(self.signal_sample)))
        getattr(self.workspace4limit_, "import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_WJets0_xww_%s_mlvj"%(self.channel)).clone("rate_WJets_xww_for_counting"))
        getattr(self.workspace4limit_, "import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_VV_xww_%s_mlvj"%(self.channel)).clone("rate_VV_xww_for_counting"))
        getattr(self.workspace4limit_, "import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_TTbar_xww_%s_mlvj"%(self.channel)).clone("rate_TTbar_xww_for_counting"))
        getattr(self.workspace4limit_, "import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_STop_xww_%s_mlvj"%(self.channel)).clone("rate_STop_xww_for_counting"))

        ### number of signal, Wjets, VV, TTbar and STop --> unbin
        getattr(self.workspace4limit_, "import")(self.workspace4fit_.var("rrv_number_%s_xww_signal_region_%s_mlvj"%(self.signal_sample, self.channel)).clone("rate_%s_xww_for_unbin"%(self.signal_sample)))
        getattr(self.workspace4limit_, "import")(self.workspace4fit_.var("rrv_number_WJets0_xww_signal_region_%s_mlvj"%(self.channel)).clone("rate_WJets_xww_for_unbin"))
        getattr(self.workspace4limit_, "import")(self.workspace4fit_.var("rrv_number_VV_xww_signal_region_%s_mlvj"%(self.channel)).clone("rate_VV_xww_for_unbin"))
        getattr(self.workspace4limit_, "import")(self.workspace4fit_.var("rrv_number_TTbar_xww_signal_region_%s_mlvj"%(self.channel)).clone("rate_TTbar_xww_for_unbin"))
        getattr(self.workspace4limit_, "import")(self.workspace4fit_.var("rrv_number_STop_xww_signal_region_%s_mlvj"%(self.channel)).clone("rate_STop_xww_for_unbin"))

        ### Set the error properly -> taking into account lumi, Vtagger and theoretical uncertainty on XS -> for VV, TTbar and STop
        self.workspace4limit_.var("rate_VV_xww_for_unbin").setError(self.workspace4limit_.var("rate_VV_xww_for_unbin").getVal()*TMath.Sqrt(self.lumi_uncertainty*self.lumi_uncertainty + self.rrv_wtagger_eff_reweight_forV.getError()/self.rrv_wtagger_eff_reweight_forV.getVal()*self.rrv_wtagger_eff_reweight_forV.getError()/self.rrv_wtagger_eff_reweight_forV.getVal() +self.XS_VV_uncertainty*self.XS_VV_uncertainty))
        self.workspace4limit_.var("rate_STop_xww_for_unbin").setError(self.workspace4limit_.var("rate_STop_xww_for_unbin").getVal()*TMath.Sqrt(self.lumi_uncertainty*self.lumi_uncertainty + self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal()*self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal() +self.XS_STop_uncertainty*self.XS_STop_uncertainty))
        self.workspace4limit_.var("rate_TTbar_xww_for_unbin").setError(self.workspace4limit_.var("rate_TTbar_xww_for_unbin").getVal()*TMath.Sqrt(self.lumi_uncertainty*self.lumi_uncertainty + self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal()*self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal()))

        ### Get the dataset for data into the signal region
        getattr(self.workspace4limit_, "import")(self.workspace4fit_.data("rdataset_data_xww_signal_region_%s_mlvj"%(self.channel)).Clone("data_obs_xww_%s_%s"%(self.channel, self.wtagger_category)))
        ### Take the corrected pdf from the alpha method for the WJets
        if mode == "sideband_correction_method1":
            getattr(self.workspace4limit_, "import")(self.workspace4fit_.pdf("model_pdf_WJets0_xww_signal_region_%s_after_correct_mlvj"%(self.channel)).clone("WJets_xww_%s_%s"%(self.channel, self.wtagger_category)))

        if isTTbarFloating:
            getattr(self.workspace4limit_, "import")(self.workspace4fit_.pdf("model_pdf_TTbar_xww_signal_region_%s_mlvj_Deco_TTbar_xww_signal_region_%s_%s_mlvj"%(self.channel, self.channel, self.wtagger_category)).clone("TTbar_xww_%s_%s"%(self.channel, self.wtagger_category)))
        else :
            getattr(self.workspace4limit_, "import")(self.workspace4fit_.pdf("model_pdf_TTbar_xww_signal_region_%s_mlvj"%(self.channel)).clone("TTbar_xww_%s_%s"%(self.channel, self.wtagger_category)))

        if isSTopFloating :
            getattr(self.workspace4limit_, "import")(self.workspace4fit_.pdf("model_pdf_STop_xww_signal_region_%s_mlvj_Deco_xww_STop_signal_region_%s_%s_mlvj"%(self.channel, self.channel, self.wtagger_category)).clone("STop_xww_%s_%s"%(self.channel, self.wtagger_category)))
        else :
            getattr(self.workspace4limit_, "import")(self.workspace4fit_.pdf("model_pdf_STop_xww_signal_region_%s_mlvj"%(self.channel)).clone("STop_xww_%s_%s"%(self.channel, self.wtagger_category)))

        if isVVFloating :
            getattr(self.workspace4limit_, "import")(self.workspace4fit_.pdf("model_pdf_VV_xww_signal_region_%s_mlvj_Deco_VV_xww_signal_region_%s_%s_mlvj"%(self.channel, self.channel, self.wtagger_category)).clone("VV_%s_%s"%(self.channel, self.wtagger_category)))
        else:
            getattr(self.workspace4limit_, "import")(self.workspace4fit_.pdf("model_pdf_VV_xww_signal_region_%s_mlvj"%(self.channel)).clone("VV_xww_%s_%s"%(self.channel, self.wtagger_category)))

        getattr(self.workspace4limit_, "import")(self.workspace4fit_.pdf("model_pdf_%s_xww_signal_region_%s_mlvj"%(self.signal_sample, self.channel)).clone(self.signal_sample+"_xww_%s_%s"%(self.channel, self.wtagger_category)))

        ### Fix all the Pdf parameters
        rrv_x = self.workspace4limit_.var("rrv_mass_lvj")

        fix_Pdf(self.workspace4limit_.pdf("TTbar_xww_%s_%s"%(self.channel, self.wtagger_category)), RooArgSet(rrv_x))
        fix_Pdf(self.workspace4limit_.pdf("STop_xww_%s_%s"%(self.channel, self.wtagger_category)), RooArgSet(rrv_x))
        fix_Pdf(self.workspace4limit_.pdf("VV_xww_%s_%s"%(self.channel, self.wtagger_category)), RooArgSet(rrv_x))
        fix_Pdf(self.workspace4limit_.pdf("WJets_xww_%s_%s"%(self.channel, self.wtagger_category)), RooArgSet(rrv_x))

        ##if TString(self.signal_sample).Contains("BulkG_WW"):
        ##    fix_Pdf(self.workspace4limit_.pdf("BulkWW_xww_%s_%s"%(self.channel, self.wtagger_category)), RooArgSet(rrv_x))
        ##else:
        ##    fix_Pdf(self.workspace4limit_.pdf(self.signal_sample+"_xww_%s_%s"%(self.channel, self.wtagger_category)), RooArgSet(rrv_x))
        fix_Pdf(self.workspace4limit_.pdf(self.signal_sample+"_xww_%s_%s"%(self.channel, self.wtagger_category)), RooArgSet(rrv_x))

        print " ############## Workspace for limit "
        parameters_workspace = self.workspace4limit_.allVars()
        par = parameters_workspace.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            param.Print()
            param = par.Next()

        params_list = []
        ### main modality for the alpha function method
        if mode == "sideband_correction_method1":

            if self.MODEL_4_mlvj == "ErfExp_v1" or self.MODEL_4_mlvj == "ErfPow_v1" or self.MODEL_4_mlvj == "2Exp" :
                ### uncertainty inflation on the Wjets shape from fitting data in sb_lo
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)
                ### Add to the parameter list
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))

                ### Do the same for alpha paramter
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig4"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig5"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                ### Add to the parameter list
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig4"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig5"%(self.channel, self.wtagger_category)))

                ### Do the same for the TTbar
                if isTTbarFloating !=0 :
                    self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_TTbar)
                    self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_TTbar)
                    self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_TTbar)

                    params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                    params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                    params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))


            if self.MODEL_4_mlvj == "ErfPow2_v1" or self.MODEL_4_mlvj == "ErfPowExp_v1" :
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig4"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig4"%(self.channel, self.wtagger_category)))

                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig4"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig5"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig6"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig7"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig4"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig5"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig6"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig7"%(self.channel, self.wtagger_category)))


                if isTTbarFloating !=0 :
                    self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_TTbar)
                    self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_TTbar)
                    self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_TTbar)
                    self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_TTbar)

                    params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                    params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                    params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
                    params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))


            if self.MODEL_4_mlvj == "Exp" or self.MODEL_4_mlvj == "Pow" :

                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))


                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))

                if isTTbarFloating !=0 :
                    self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_TTbar)
                    params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))

            if self.MODEL_4_mlvj == "ExpN" or self.MODEL_4_mlvj == "ExpTail" or self.MODEL_4_mlvj == "Pow2" :

                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_WJets0)

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))

                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_alpha)

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))


                ### TTbar use exp
                if isTTbarFloating !=0:
                    print "##################### TTbar will float in the limit procedure + final plot ######################"
                    self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_TTbar)
                    params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))

                ### VV use ExpTail:
                if isVVFloating !=0:
                    print "##################### VV will float in the limit procedure + final plot ######################"
                    self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_VV)
                    self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_VV)
                    params_list.append(self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                    params_list.append(self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))

                ### STop use Exp:
                if isSTopFloating !=0:
                    print "##################### STop will float in the limit procedure + final plot ######################"
                    self.workspace4limit_.var("Deco_STop_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).setError(self.shape_para_error_STop)
                    params_list.append(self.workspace4limit_.var("Deco_STop_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))


        #### add signal shape parameters' uncertainty -> increase the uncertainty on the mean and the sigma since we are using a CB or a Double CB or a BWxDB or BWxCB
        if self.workspace4limit_.var("rrv_mean_CB_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category)):

            self.workspace4limit_.var("rrv_mean_shift_scale_lep_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category)).setError(self.mean_signal_uncertainty_lep_scale)
            self.workspace4limit_.var("rrv_mean_shift_scale_jes_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category)).setError(self.mean_signal_uncertainty_jet_scale)
            self.workspace4limit_.var("rrv_sigma_shift_lep_scale_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category)).setError(self.sigma_signal_uncertainty_lep_scale)
            self.workspace4limit_.var("rrv_sigma_shift_jes_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category)).setError(self.sigma_signal_uncertainty_jet_scale)
            self.workspace4limit_.var("rrv_sigma_shift_res_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category)).setError(self.sigma_signal_uncertainty_jet_res)

            if self.channel == "mu":
                self.workspace4limit_.var("CMS_sig_p1_scale_m").setError(1)
                self.workspace4limit_.var("CMS_sig_p2_scale_m").setError(1)
                params_list.append(self.workspace4limit_.var("CMS_sig_p1_scale_m"))
                params_list.append(self.workspace4limit_.var("CMS_sig_p2_scale_m"))
            elif self.channel == "el":
                self.workspace4limit_.var("CMS_sig_p1_scale_e").setError(1)
                self.workspace4limit_.var("CMS_sig_p2_scale_e").setError(1)
                params_list.append(self.workspace4limit_.var("CMS_sig_p1_scale_e"))
                params_list.append(self.workspace4limit_.var("CMS_sig_p2_scale_e"))
            elif self.channel == "em":
                self.workspace4limit_.var("CMS_sig_p1_scale_em").setError(1)
                self.workspace4limit_.var("CMS_sig_p2_scale_em").setError(1)
                params_list.append(self.workspace4limit_.var("CMS_sig_p1_scale_em"))
                params_list.append(self.workspace4limit_.var("CMS_sig_p2_scale_em"))

            self.workspace4limit_.var("CMS_sig_p1_jes").setError(1)
            self.workspace4limit_.var("CMS_sig_p2_jes").setError(1)
            self.workspace4limit_.var("CMS_sig_p2_jer").setError(1)

            params_list.append(self.workspace4limit_.var("CMS_sig_p1_jes"))
            params_list.append(self.workspace4limit_.var("CMS_sig_p2_jes"))
            params_list.append(self.workspace4limit_.var("CMS_sig_p2_jer"))


        ### calculate the shape uncertainty for cut-and-counting
        self.rrv_counting_uncertainty_from_shape_uncertainty = RooRealVar("rrv_counting_uncertainty_from_shape_uncertainty_%s"%(self.channel), "rrv_counting_uncertainty_from_shape_uncertainty_%s"%(self.channel), 0)
        self.rrv_counting_uncertainty_from_shape_uncertainty.setError(Calc_error("WJets_xww_%s_%s"%(self.channel, self.wtagger_category), "rrv_mass_lvj" , self.FloatingParams, self.workspace4limit_, "signal_region"))
        self.rrv_counting_uncertainty_from_shape_uncertainty.Print()

        print " param list ", params_list

        ### Print the datacard for unbin and couting analysis
        self.print_limit_datacard("unbin", params_list)
        self.print_limit_datacard("counting")

        self.setFloatingParams(self.workspace4limit_, self.FloatingParams, mode, isTTbarFloating, isVVFloating, isSTopFloating)
        ### Add the floating list to the combiner --> the pdf which are not fixed are floating by default
        getattr(self.workspace4limit_, "import")(self.FloatingParams)

        ### Save the workspace
        self.save_workspace_to_file()


    ##adding floating params to a list
    def setFloatingParams(self, workspace, param_list, mode, isTTbarFloating = 0, isVVFloating = 0, isSTopFloating = 0):
        if mode == "sideband_correction_method1":
            if self.MODEL_4_mlvj == "ErfExp_v1" or self.MODEL_4_mlvj == "ErfPow_v1" or self.MODEL_4_mlvj == "2Exp" :
                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))

                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig4"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig5"%(self.channel, self.wtagger_category)))

                if isTTbarFloating!=0:
                    param_list.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                    param_list.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                    param_list.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))

            elif self.MODEL_4_mlvj == "ErfPow2_v1" or self.MODEL_4_mlvj == "ErfPowExp_v1" :

                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig4"%(self.channel, self.wtagger_category)))

                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig4"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig5"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig6"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig7"%(self.channel, self.wtagger_category)))

                if isTTbarFloating!=0:
                    param_list.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                    param_list.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                    param_list.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
                    param_list.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))

            elif self.MODEL_4_mlvj == "Exp" or self.MODEL_4_mlvj == "Pow" :

                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))

                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))

                if isTTbarFloating!=0:
                    param_list.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))


            elif self.MODEL_4_mlvj == "ExpN" or self.MODEL_4_mlvj == "ExpTail" or self.MODEL_4_mlvj == "Pow2" :

                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))


                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
                param_list.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))

                if isTTbarFloating!=0:
                    param_list.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))

                if isVVFloating!=0:
                    param_list.add(workspace.var("Deco_VV_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
                    param_list.add(workspace.var("Deco_VV_xww_signal_region_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))

                if isSTopFloating!=0:
                    param_list.add(workspace.var("Deco_STop_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
            else:
                print self.MODEL_4_mlvj
                raw_input("wrong MODEL_4_mlvj")

            if workspace.var("rrv_mean_CB_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category)):
                if self.channel == "mu":
                    param_list.add(workspace.var("CMS_sig_p1_scale_m"))
                    param_list.add(workspace.var("CMS_sig_p2_scale_m"))
                elif self.channel == "el":
                    param_list.add(workspace.var("CMS_sig_p1_scale_e"))
                    param_list.add(workspace.var("CMS_sig_p2_scale_e"))
                elif self.channel == "em":
                    param_list.add(workspace.var("CMS_sig_p1_scale_em"))
                    param_list.add(workspace.var("CMS_sig_p2_scale_em"))

                param_list.add(workspace.var("CMS_sig_p1_jes"))
                param_list.add(workspace.var("CMS_sig_p2_jes"))
                param_list.add(workspace.var("CMS_sig_p2_jer"))
        else:
            raw_input("this mode:%s is not support!"%(mode))

    #### Method used in order to save the workspace in a output root file
    def save_workspace_to_file(self):
        self.workspace4limit_.writeToFile(self.file_rlt_root)
        self.file_out.close()

    #### Method used to print the general format of the datacard for both counting and unbinned analysis
    def print_limit_datacard(self, mode, params_list = []):
        print "############## print_limit_datacard for %s ################"%(mode)
        if not (mode == "unbin" or mode == "counting"):
            print "print_limit_datacard use wrong mode: %s"%(mode)
            raw_input("ENTER")

        ### open the datacard
        datacard_out = open(getattr(self, "file_datacard_%s"%(mode)), "w")

        ### start to print inside
        datacard_out.write("imax 1")
        datacard_out.write("\njmax 4")
        datacard_out.write("\nkmax *")
        datacard_out.write("\n--------------- ")

        if mode == "unbin":
            fnOnly = ntpath.basename(self.file_rlt_root) ## workspace for limit --> output file for the workspace
            datacard_out.write("\nshapes %s_xww  CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category, fnOnly, self.workspace4limit_.GetName(), self.channel, self.wtagger_category))
            datacard_out.write("\nshapes WJets_xww  CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel, self.wtagger_category, fnOnly, self.workspace4limit_.GetName(), self.channel, self.wtagger_category))
            datacard_out.write("\nshapes TTbar_xww  CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel, self.wtagger_category, fnOnly, self.workspace4limit_.GetName(), self.channel, self.wtagger_category))
            datacard_out.write("\nshapes STop_xww   CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel, self.wtagger_category, fnOnly, self.workspace4limit_.GetName(), self.channel, self.wtagger_category))
            datacard_out.write("\nshapes VV_xww     CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel, self.wtagger_category, fnOnly, self.workspace4limit_.GetName(), self.channel, self.wtagger_category))
            datacard_out.write("\nshapes data_obs   CMS_xww_%s1J%s  %s %s:$PROCESS_xww_%s_%s"%(self.channel, self.wtagger_category, fnOnly, self.workspace4limit_.GetName(), self.channel, self.wtagger_category))
            datacard_out.write("\n--------------- ")

        datacard_out.write("\nbin CMS_xww_%s1J%s "%(self.channel, self.wtagger_category))
        if mode == "unbin":
            datacard_out.write("\nobservation %0.2f "%(self.workspace4limit_.data("data_obs_xww_%s_%s"%(self.channel, self.wtagger_category)).sumEntries()))
        if mode == "counting":
            datacard_out.write("\nobservation %0.2f "%(self.workspace4limit_.var("observation_for_counting_xww").getVal()))

        datacard_out.write("\n------------------------------")
        datacard_out.write("\nbin CMS_xww_%s1J%s CMS_xww_%s1J%s CMS_xww_%s1J%s CMS_xww_%s1J%s CMS_xww_%s1J%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category))
        datacard_out.write("\nprocess %s_xww WJets_xww TTbar_xww STop_xww VV_xww "%(self.signal_sample)) ## just one signal sample
        datacard_out.write("\nprocess -1 1 2 3 4")

        ### rates for the different process
        if mode == "unbin":
            datacard_out.write("\nrate %0.7f %0.3f %0.3f %0.3f %0.3f "%(self.signal_scale*self.workspace4limit_.var("rate_%s_xww_for_unbin"%(self.signal_sample)).getVal(), self.workspace4limit_.var("rate_WJets_xww_for_unbin").getVal(), self.workspace4limit_.var("rate_TTbar_xww_for_unbin").getVal(), self.workspace4limit_.var("rate_STop_xww_for_unbin").getVal(), self.workspace4limit_.var("rate_VV_xww_for_unbin").getVal()))
        elif mode == "counting":
            datacard_out.write("\nrate %0.7f %0.3f %0.3f %0.3f %0.3f "%(self.signal_scale*self.workspace4limit_.var("rate_%s_xww_for_counting"%(self.signal_sample)).getVal(), self.workspace4limit_.var("rate_WJets_xww_for_counting").getVal(), self.workspace4limit_.var("rate_TTbar_xww_for_counting").getVal(), self.workspace4limit_.var("rate_STop_xww_for_counting").getVal(), self.workspace4limit_.var("rate_VV_xww_for_counting").getVal()))
        else: raw_input("wrong mode:"+mode)

        datacard_out.write("\n-------------------------------- ")
        ### luminosity nouisance
        datacard_out.write("\nlumi_13TeV lnN %0.3f - %0.3f %0.3f %0.3f"%(1.+self.lumi_uncertainty, 1.+self.lumi_uncertainty, 1.+self.lumi_uncertainty, 1.+self.lumi_uncertainty))
        ### Signal XS  nouisance in boosted regime
        datacard_out.write("\nCMS_xww_XS_Signal lnN %0.3f - - - -"%(1+self.XS_Signal_uncertainty))
        ### STop XS  nouisance in boosted regime
        datacard_out.write("\nCMS_xww_XS_STop lnN - - - %0.3f -"%(1+self.XS_STop_uncertainty))
        ### VV XS  nouisance in boosted regime
        datacard_out.write("\nCMS_xww_XS_VV lnN - - - - %0.3f"%(1+self.XS_VV_uncertainty))
        ### WJets Normalization from data fit -> data driven
        if self.number_WJets_insideband >0:
            datacard_out.write("\nCMS_xww_WJ_norm gmN %0.3f %0.3f - - -"%(self.number_WJets_insideband, getattr(self, "datadriven_alpha_WJets_xww_%s"%(mode))))
        else:
            datacard_out.write("\nCMS_xww_WJ_norm_%s_%s lnN - %0.3f - - -"%(self.channel, self.wtagger_category, 1+ self.workspace4limit_.var("rate_WJets_xww_for_unbin").getError()/self.workspace4limit_.var("rate_WJets_xww_for_unbin").getVal()))

        ### Top normalization due to SF in the ttbar CR
        ##datacard_out.write("\nCMS_xww_Top_norm_%s_%s lnN - - %0.3f/%0.3f %0.3f/%0.3f -"%(self.channel, self.wtagger_category, 1+self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1-self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1+self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1-self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal()))
        datacard_out.write("\nCMS_xww_Top_norm_%s_%s lnN - %0.3f/%0.3f %0.3f/%0.3f %0.3f/%0.3f -"%(self.channel, self.wtagger_category, 1-0.6*self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1+0.6*self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1+self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1-self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1+self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1-self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal()))



        ### V-Tagger SF nouisance
        if self.wtagger_category == "HP":
            datacard_out.write("\nCMS_eff_vtag_tau21_sf lnN %0.3f/%0.3f - - - %0.3f/%0.3f"%(1+self.rrv_wtagger_eff_reweight_forV.getError(), 1-self.rrv_wtagger_eff_reweight_forV.getError(), 1+self.rrv_wtagger_eff_reweight_forV.getError(), 1-self.rrv_wtagger_eff_reweight_forV.getError()))
        elif self.wtagger_category == "LP":
            datacard_out.write("\nCMS_eff_vtag_tau21_sf lnN %0.3f/%0.3f - - - %0.3f/%0.3f"%(1-self.rrv_wtagger_eff_reweight_forV.getError(), 1+self.rrv_wtagger_eff_reweight_forV.getError(), 1-self.rrv_wtagger_eff_reweight_forV.getError(), 1+self.rrv_wtagger_eff_reweight_forV.getError()))

        ### btag scale factor on the MC background
        #calculate wjets btag unc 
        self.btag_scale_wjets_uncertainty= -1*(self.btag_scale_ttbar_uncertainty* self.workspace4limit_.var("rate_TTbar_xww_for_unbin").getVal() + self.btag_scale_singletop_uncertainty* self.workspace4limit_.var("rate_STop_xww_for_unbin").getVal())/ self.workspace4limit_.var("rate_WJets_xww_for_unbin").getVal() 
        datacard_out.write("\nCMS_xww_btagger lnN - %0.3f %0.3f %0.3f -"%( 1+self.btag_scale_wjets_uncertainty, 1+self.btag_scale_ttbar_uncertainty, 1+self.btag_scale_singletop_uncertainty))

        ### btag scale factor on the MC background
        datacard_out.write("\n#CMS_eff_vtag_model lnN %0.3f - - - %0.3f"%(1+self.eff_vtag_model, 1+self.eff_vtag_model))

        ### jet Mass effect only if available -> shapes changing due to the jet mass uncertainty (JEC for CA8/AK7) -> affects also WJets
        if (self.workspace4fit_.var("rrv_number_WJets0_xww_massup_in_mj_signal_region_from_fitting_%s"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_WJets0_xww_massdn_in_mj_signal_region_from_fitting_%s"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww_massup_%s_mj"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww_massdown_%s_mj"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww_massup_%s_mj"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww_massdn_%s_mj"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_dataset_signal_region_VV_xww_massup_%s_mj"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_dataset_signal_region_VV_xww_massdn_%s_mj"%(self.channel))) :
            datacard_out.write("\nJetMass_%s lnN - %0.3f %0.3f %0.3f %0.3f"%(self.channel, 1+self.WJets_normlization_uncertainty_from_jet_mass, 1+self.TTbar_normlization_uncertainty_from_jet_mass, 1+self.STop_normlization_uncertainty_from_jet_mass, 1+self.VV_normlization_uncertainty_from_jet_mass))

        if self.channel == "mu":
            self.channel_short = "m"
        elif self.channel == "el":
            self.channel_short = "e"
        elif self.channel == "em":
            self.channel_short = "em"

        ### trigger efficiency
        datacard_out.write("\nCMS_xww_trigger_%s lnN %0.3f - %0.3f %0.3f %0.3f"%(self.channel_short, 1+self.lep_trigger_uncertainty, 1+self.lep_trigger_uncertainty, 1+self.lep_trigger_uncertainty, 1+self.lep_trigger_uncertainty))
        ### Lepton SF
        datacard_out.write("\nCMS_eff_%s lnN %0.3f - %0.3f %0.3f %0.3f"%(self.channel_short, 1+self.lep_eff_uncertainty, 1+self.lep_eff_uncertainty, 1+self.lep_eff_uncertainty, 1+self.lep_eff_uncertainty))

        ############ Evaluated just for signal, in principle also on all the backgrounds with the same topology
        ### Lepton Energy scale
        datacard_out.write("\nCMS_scale_%s lnN %0.3f - - - -"%(self.channel_short, 1+self.signal_lepton_energy_scale_uncertainty))
        ### Lepton Energy Resolution
        datacard_out.write("\nCMS_res_%s lnN %0.3f - - - -"%(self.channel_short, 1+self.signal_lepton_energy_res_uncertainty))
        ### fat jet energy scale
        datacard_out.write("\nCMS_scale_j  lnN %0.3f/%0.3f - - - -"%(self.signal_jet_energy_scale_uncertainty_low, self.signal_jet_energy_scale_uncertainty_high))
        ### fat jet mass scale
        datacard_out.write("\nCMS_mass_scale_j  lnN %0.3f/%0.3f - - - -"%(self.signal_jet_mass_scale_uncertainty_low, self.signal_jet_mass_scale_uncertainty_high))
        ### fat jet mass res
        datacard_out.write("\nCMS_mass_res_j  lnN %0.3f/%0.3f - - - -"%(self.signal_jet_mass_res_uncertainty_low, self.signal_jet_mass_res_uncertainty_high))
        ### fat jet energy resolution
        datacard_out.write("\nCMS_res_j  lnN %0.3f - - - -"%(1+self.signal_jet_energy_res_uncertainty))
        ### btag on the signal
        datacard_out.write("\nCMS_xww_btag_eff lnN %0.3f - - - -"%(1+self.btag_scale_signal_uncertainty))

        ### print shapes parameter to be taken int account
        if mode == "unbin":
            for ipar in params_list:
                print "Name %s", ipar.GetName()
                if TString(ipar.GetName()).Contains("Deco_TTbar_xww_signal_region"):
                    datacard_out.write("\n%s param %0.1f %0.1f "%(ipar.GetName(), ipar.getVal(), ipar.getError()))
                else:
                    datacard_out.write("\n%s param %0.1f %0.1f "%(ipar.GetName(), ipar.getVal(), ipar.getError()))
        if mode == "counting":
            datacard_out.write("\nShape_%s_%s lnN - - %0.3f - - -"%(self.channel, self.wtagger_category, 1+self.rrv_counting_uncertainty_from_shape_uncertainty.getError()))


    #### Read the final workspace and produce the latest plots
    def read_workspace(self, logy=1, mode = "sideband_correction_method1", isTTbarFloating=1, isVVFloating=0, isSTopFloating=0):
        print "--------------------- read_workspace -------------------------"

        ### Taket the workspace for limits
        ##file = TFile(self.file_rlt_root)
        ##workspace = file.Get("workspace4limit_")
        workspace = TFile(self.file_rlt_root).Get("workspace4limit_")
        workspace.Print()

        ### iterate on the workspace element parameters
        print "----------- Parameter Workspace -------------"
        parameters_workspace = workspace.allVars()
        par = parameters_workspace.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            param.Print()
            param = par.Next()
        print "---------------------------------------------"

        workspace.data("data_obs_xww_%s_%s"%(self.channel, self.wtagger_category)).Print()

        print "----------- Pdf in the Workspace -------------"
        pdfs_workspace = workspace.allPdfs()
        par = pdfs_workspace.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            param.Print()
            param = par.Next()
        print "----------------------------------------------"

        rrv_x = workspace.var("rrv_mass_lvj")
        data_obs = workspace.data("data_obs_xww_%s_%s"%(self.channel, self.wtagger_category))
        #if TString(self.signal_sample).Contains("RS1G_WW"):
        #    model_pdf_signal = workspace.pdf("RSWW_xww_%s_%s"%(self.channel, self.wtagger_category))
        #elif TString(self.signal_sample).Contains("BulkGraviton"):
        #    model_pdf_signal = workspace.pdf("BulkWW_xww_%s_%s"%(self.channel, self.wtagger_category))
        #elif TString(self.signal_sample).Contains("Wprime_WZ"):
        #    model_pdf_signal = workspace.pdf("WprimeWZ_xww_%s_%s"%(self.channel, self.wtagger_category))
        #else:
        #    model_pdf_signal = workspace.pdf("%s_xww_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category))
        model_pdf_signal = workspace.pdf("%s_xww_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category))

        model_pdf_WJets  = workspace.pdf("WJets_xww_%s_%s"%(self.channel, self.wtagger_category))
        model_pdf_VV     = workspace.pdf("VV_xww_%s_%s"%(self.channel, self.wtagger_category))
        model_pdf_TTbar  = workspace.pdf("TTbar_xww_%s_%s"%(self.channel, self.wtagger_category))
        model_pdf_STop   = workspace.pdf("STop_xww_%s_%s"%(self.channel, self.wtagger_category))

        model_pdf_signal.Print()
        model_pdf_WJets.Print()
        model_pdf_VV.Print()
        model_pdf_TTbar.Print()
        model_pdf_STop.Print()

        #if TString(self.signal_sample).Contains("RS1G_WW"):
        #    rrv_number_signal = workspace.var("rate_RSWW_xww_for_unbin")
        #elif TString(self.signal_sample).Contains("BulkGraviton"):
        #    rrv_number_signal = workspace.var("rate_BulkWW_xww_for_unbin")
        #elif TString(self.signal_sample).Contains("Wprime_WZ"):
        #    rrv_number_signal = workspace.var("rate_WprimeWZ_xww_for_unbin")
        #else:
        #    rrv_number_signal = workspace.var("rate_%s_xww_for_unbin"%(self.signal_sample))
        rrv_number_signal = workspace.var("rate_%s_xww_for_unbin"%(self.signal_sample))


        rrv_number_WJets  = workspace.var("rate_WJets_xww_for_unbin")
        rrv_number_VV     = workspace.var("rate_VV_xww_for_unbin")
        rrv_number_TTbar  = workspace.var("rate_TTbar_xww_for_unbin")
        rrv_number_STop   = workspace.var("rate_STop_xww_for_unbin")

        rrv_number_signal.Print()
        rrv_number_WJets.Print()
        rrv_number_VV.Print()
        rrv_number_TTbar.Print()
        rrv_number_STop.Print()

        #### Prepare the final plot starting from total background
        rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC_xww", "rrv_number_Total_background_MC_xww",
                rrv_number_WJets.getVal()+
                rrv_number_VV.getVal()+
                rrv_number_TTbar.getVal()+
                rrv_number_STop.getVal())

        rrv_number_Total_background_MC.setError(TMath.Sqrt(
            rrv_number_WJets.getError()* rrv_number_WJets.getError()+
            rrv_number_VV.getError()* rrv_number_VV.getError()+
            rrv_number_TTbar.getError()* rrv_number_TTbar.getError()+
            rrv_number_STop.getError() *rrv_number_STop.getError()))

        #### Total pdf
        model_Total_background_MC = RooAddPdf("model_Total_background_MC_xww", "model_Total_background_MC_xww", RooArgList(model_pdf_WJets, model_pdf_VV, model_pdf_TTbar, model_pdf_STop), RooArgList(rrv_number_WJets, rrv_number_VV, rrv_number_TTbar, rrv_number_STop))

        if data_obs.sumEntries() != 0:
            #### scale factor in order to scale MC to data in the final plot -> in order to avoid the normalization to data which is done by default in rooFit
            scale_number_signal = rrv_number_signal.getVal()/data_obs.sumEntries()
            scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()/data_obs.sumEntries()
        else:
            scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()
            scale_number_signal = rrv_number_signal.getVal()
        #### create the frame
        mplot = rrv_x.frame(RooFit.Title("check_workspace"), RooFit.Bins(int(rrv_x.getBins()/self.binwidth_narrow_factor)))
        mplotP = rrv_x.frame(RooFit.Title("check_workspaceP"), RooFit.Bins(int(rrv_x.getBins()/self.binwidth_narrow_factor)))
        data_obs.plotOn(mplot , RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0))
        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("WJets"), RooFit.Components("WJets_xww_%s_%s,VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("VV"), RooFit.Components("VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("TTbar"), RooFit.Components("TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("STop"), RooFit.Components("STop_xww_%s_%s"%(self.channel, self.wtagger_category)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())

        #solid line
        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("WJets_line_invisible"), RooFit.Components("WJets_xww_%s_%s,VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("VV_line_invisible"), RooFit.Components("VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("TTbar_line_invisible"), RooFit.Components("TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("STop_line_invisible"), RooFit.Components("STop_xww_%s_%s"%(self.channel, self.wtagger_category)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_pdf_signal.plotOn(mplot, RooFit.Normalization(scale_number_signal*self.signal_scale*self.signal_scale_plot), RooFit.Name("%s #times %s"%(self.signal_sample, self.signal_scale*self.signal_scale_plot)), RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["Signal"]), RooFit.LineStyle(2), RooFit.VLines())

        #### plot the observed data using poissonian error bar
        self.plot_data_with_poissoninterval(rrv_x, data_obs, mplot)

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Invisible())

        mplot_pull = self.get_pull(rrv_x, mplot)

        datahist = data_obs.binnedClone(data_obs.GetName()+"_binnedClone", data_obs.GetName()+"_binnedClone")
        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        if self.FloatingParams.getSize() == 0:
            self.setFloatingParams(workspace, self.FloatingParams, mode, isTTbarFloating, isVVFloating, isSTopFloating)
        self.FloatingParams.Print("v")
        #raw_input("zixu")

        hdata = datahist.createHistogram(rrv_x.GetName(), int(rrv_x.getBins()/self.binwidth_narrow_factor))
        draw_error_band(model_Total_background_MC, rrv_x.GetName(), rrv_number_Total_background_MC, self.FloatingParams,
                workspace, mplot, mplotP, hdata, self.color_palet["Uncertainty"], "F")

        mplot.Print()
        self.leg = self.legend4Plot(mplot, 0, 1,-0.01,-0.05, 0.11, 0.)
        #self.leg.SetTextSize(0.036)
        mplot.addObject(self.leg)
        #pt1 = ROOT.TPaveText(0.6180905, 0.4355644, 0.8291457, 0.507992, "NDC")
        #pt1.SetTextFont(42)
        #pt1.SetTextSize(0.05)
        #pt1.SetFillColor(0)
        #pt1.SetFillStyle(0)
        #pt1.SetBorderSize(0)
        #text = pt1.AddText("")
        #if options.category.find('Z') != -1: text = pt1.AddText("WZ category")
        #elif options.category.find('W') != -1: text = pt1.AddText("WW category")
        #text.SetTextFont(62)
        #mplot.addObject(pt1)

        mplot.GetYaxis().SetRangeUser(1e-2, mplot.GetMaximum()*1.2)

        if workspace.var("rrv_num_floatparameter_in_last_fitting"):
            self.nPar_float_in_fitTo = int(workspace.var("rrv_num_floatparameter_in_last_fitting").getVal())
        else:
            self.nPar_float_in_fitTo = self.FloatingParams.getSize()
        nBinX = mplot.GetNbinsX()
        ndof  = nBinX-self.nPar_float_in_fitTo
        print "nPar = %s, chiSquare = %s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof)

        parameters_list = RooArgList()
        self.draw_canvas_with_pull2(rrv_x, datahist, mplot, mplotP, mplot_pull, ndof, parameters_list, "%s/m_lvj_fitting/"%(self.plotsDir), "check_workspace_for_limit", self.MODEL_4_mlvj, 0, 1)


    #### Read the final workspace and produce the latest plots
    def read_workspace_postfit(self, logy=1, mode = "sideband_correction_method1", isTTbarFloating=1, isVVFloating=0, isSTopFloating=0):
        print "--------------------- read_workspace -------------------------"

        ### Taket the workspace for limits
        ##file = TFile(self.file_rlt_root)
        ##workspace = file.Get("workspace4limit_")
        workspace = TFile(self.file_rlt_root).Get("workspace4limit_")
        workspace.Print()

        ### iterate on the workspace element parameters
        print "----------- Parameter Workspace -------------"
        parameters_workspace = workspace.allVars()
        par = parameters_workspace.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            param.Print()
            param = par.Next()
        print "---------------------------------------------"

        workspace.data("data_obs_xww_%s_%s"%(self.channel, self.wtagger_category)).Print()

        print "----------- Pdf in the Workspace -------------"
        pdfs_workspace = workspace.allPdfs()
        par = pdfs_workspace.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            param.Print()
            param = par.Next()
        print "----------------------------------------------"

        rrv_x = workspace.var("rrv_mass_lvj")
        data_obs = workspace.data("data_obs_xww_%s_%s"%(self.channel, self.wtagger_category))
        #if TString(self.signal_sample).Contains("RS1G_WW"):
        #    model_pdf_signal = workspace.pdf("RSWW_xww_%s_%s"%(self.channel, self.wtagger_category))
        #elif TString(self.signal_sample).Contains("BulkGraviton"):
        #    model_pdf_signal = workspace.pdf("BulkWW_xww_%s_%s"%(self.channel, self.wtagger_category))
        #elif TString(self.signal_sample).Contains("Wprime_WZ"):
        #    model_pdf_signal = workspace.pdf("WprimeWZ_xww_%s_%s"%(self.channel, self.wtagger_category))
        #else:
        #    model_pdf_signal = workspace.pdf("%s_xww_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category))
        model_pdf_signal = workspace.pdf("%s_xww_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category))

        model_pdf_WJets  = workspace.pdf("WJets_xww_%s_%s"%(self.channel, self.wtagger_category))
        model_pdf_VV     = workspace.pdf("VV_xww_%s_%s"%(self.channel, self.wtagger_category))
        model_pdf_TTbar  = workspace.pdf("TTbar_xww_%s_%s"%(self.channel, self.wtagger_category))
        model_pdf_STop   = workspace.pdf("STop_xww_%s_%s"%(self.channel, self.wtagger_category))

        model_pdf_signal.Print()
        model_pdf_WJets.Print()
        model_pdf_VV.Print()
        model_pdf_TTbar.Print()
        model_pdf_STop.Print()

        #if TString(self.signal_sample).Contains("RS1G_WW"):
        #    rrv_number_signal = workspace.var("rate_RSWW_xww_for_unbin")
        #elif TString(self.signal_sample).Contains("BulkGraviton"):
        #    rrv_number_signal = workspace.var("rate_BulkWW_xww_for_unbin")
        #elif TString(self.signal_sample).Contains("Wprime_WZ"):
        #    rrv_number_signal = workspace.var("rate_WprimeWZ_xww_for_unbin")
        #else:
        #    rrv_number_signal = workspace.var("rate_%s_xww_for_unbin"%(self.signal_sample))
        rrv_number_signal = workspace.var("rate_%s_xww_for_unbin"%(self.signal_sample))


        rrv_number_WJets  = workspace.var("rate_WJets_xww_for_unbin")
        rrv_number_VV     = workspace.var("rate_VV_xww_for_unbin")
        rrv_number_TTbar  = workspace.var("rate_TTbar_xww_for_unbin")
        rrv_number_STop   = workspace.var("rate_STop_xww_for_unbin")

        rrv_number_signal.Print()
        rrv_number_WJets.Print()
        rrv_number_VV.Print()
        rrv_number_TTbar.Print()
        rrv_number_STop.Print()

        #### Prepare the final plot starting from total background
        rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC_xww", "rrv_number_Total_background_MC_xww",
                rrv_number_WJets.getVal()+
                rrv_number_VV.getVal()+
                rrv_number_TTbar.getVal()+
                rrv_number_STop.getVal())

        rrv_number_Total_background_MC.setError(TMath.Sqrt(
            rrv_number_WJets.getError()* rrv_number_WJets.getError()+
            rrv_number_VV.getError()* rrv_number_VV.getError()+
            rrv_number_TTbar.getError()* rrv_number_TTbar.getError()+
            rrv_number_STop.getError() *rrv_number_STop.getError()))

        #### Total pdf
        model_Total_background_MC = RooAddPdf("model_Total_background_MC_xww", "model_Total_background_MC_xww", RooArgList(model_pdf_WJets, model_pdf_VV, model_pdf_TTbar, model_pdf_STop), RooArgList(rrv_number_WJets, rrv_number_VV, rrv_number_TTbar, rrv_number_STop))

        if data_obs.sumEntries() != 0:
            #### scale factor in order to scale MC to data in the final plot -> in order to avoid the normalization to data which is done by default in rooFit
            scale_number_signal = rrv_number_signal.getVal()/data_obs.sumEntries()
            scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()/data_obs.sumEntries()
        else:
            scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()
            scale_number_signal = rrv_number_signal.getVal()
        #### create the frame
        mplot = rrv_x.frame(RooFit.Title("check_workspace"), RooFit.Bins(int(rrv_x.getBins()/self.binwidth_narrow_factor)))
        mplotP = rrv_x.frame(RooFit.Title("check_workspaceP"), RooFit.Bins(int(rrv_x.getBins()/self.binwidth_narrow_factor)))
        data_obs.plotOn(mplot , RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0))
        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("WJets"), RooFit.Components("WJets_xww_%s_%s,VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("VV"), RooFit.Components("VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("TTbar"), RooFit.Components("TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("STop"), RooFit.Components("STop_xww_%s_%s"%(self.channel, self.wtagger_category)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())

        #solid line
        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("WJets_line_invisible"), RooFit.Components("WJets_xww_%s_%s,VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("VV_line_invisible"), RooFit.Components("VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("TTbar_line_invisible"), RooFit.Components("TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("STop_line_invisible"), RooFit.Components("STop_xww_%s_%s"%(self.channel, self.wtagger_category)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_pdf_signal.plotOn(mplot, RooFit.Normalization(scale_number_signal*self.signal_scale*self.signal_scale_plot), RooFit.Name("%s #times %s"%(self.signal_sample, self.signal_scale*self.signal_scale_plot)), RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["Signal"]), RooFit.LineStyle(2), RooFit.VLines())

        #### plot the observed data using poissonian error bar
        self.plot_data_with_poissoninterval(rrv_x, data_obs, mplot)

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Invisible())

        mplot_pull = self.get_pull(rrv_x, mplot)

        datahist = data_obs.binnedClone(data_obs.GetName()+"_binnedClone", data_obs.GetName()+"_binnedClone")
        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        if self.FloatingParams.getSize() == 0:
            self.setFloatingParams(workspace, self.FloatingParams, mode, isTTbarFloating, isVVFloating, isSTopFloating)
        self.FloatingParams.Print("v")
        #raw_input("zixu")

        hdata = datahist.createHistogram(rrv_x.GetName(), int(rrv_x.getBins()/self.binwidth_narrow_factor))
        draw_error_band(model_Total_background_MC, rrv_x.GetName(), rrv_number_Total_background_MC, self.FloatingParams,
                workspace, mplot, mplotP, hdata, self.color_palet["Uncertainty"], "F")

        mplot.Print()
        self.leg = self.legend4Plot(mplot, 0, 1,-0.01,-0.05, 0.11, 0.)
        #self.leg.SetTextSize(0.036)
        mplot.addObject(self.leg)
        #pt1 = ROOT.TPaveText(0.6180905, 0.4355644, 0.8291457, 0.507992, "NDC")
        #pt1.SetTextFont(42)
        #pt1.SetTextSize(0.05)
        #pt1.SetFillColor(0)
        #pt1.SetFillStyle(0)
        #pt1.SetBorderSize(0)
        #text = pt1.AddText("")
        #if options.category.find('Z') != -1: text = pt1.AddText("WZ category")
        #elif options.category.find('W') != -1: text = pt1.AddText("WW category")
        #text.SetTextFont(62)
        #mplot.addObject(pt1)

        mplot.GetYaxis().SetRangeUser(1e-2, mplot.GetMaximum()*1.2)

        if workspace.var("rrv_num_floatparameter_in_last_fitting"):
            self.nPar_float_in_fitTo = int(workspace.var("rrv_num_floatparameter_in_last_fitting").getVal())
        else:
            self.nPar_float_in_fitTo = self.FloatingParams.getSize()
        nBinX = mplot.GetNbinsX()
        ndof  = nBinX-self.nPar_float_in_fitTo
        print "nPar = %s, chiSquare = %s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof)

        parameters_list = RooArgList()
        self.draw_canvas_with_pull2(rrv_x, datahist, mplot, mplotP, mplot_pull, ndof, parameters_list, "%s/m_lvj_fitting/"%(self.plotsDir), "check_workspace_for_limit", self.MODEL_4_mlvj, 0, 1)



##    def read_workspace_postfit(self, logy = 0, isTTbarFloating = 1, isVVFloating = 0, isSTopFloating = 0):
##        print "--------------------- read_workspace -------------------------"
##
##        ### Taket the workspace for limits
##        ##file = TFile(self.file_rlt_root)
##        ##workspace = file.Get("workspace4limit_")
##        workspace = TFile(self.file_rlt_root).Get("workspace4limit_")
##        workspace.Print()
##
##
##        self.file_rlt_root_postfit = self.datacardsDir+"/higgsCombine_lim_BulkGravWW750_mu_HP.Asymptotic.mH750.root"
##        file_postfit = TFile(self.file_rlt_root_postfit)
##        raw_input("zijun postfit 1")
##        workspace_postfit = file_postfit.Get("w")
##        raw_input("zijun postfit 1")
##        file_postfit.Print()
##        workspace_postfit.Print()
##        raw_input("zijun postfit 1")
##        #print workspace_postfit.loadSnapshot("clean")
##
##
##
##        ### iterate on the workspace element parameters
##        print "----------- Parameter Workspace -------------"
##        parameters_workspace = workspace.allVars()
##        par = parameters_workspace.createIterator()
##        par.Reset()
##        param = par.Next()
##        while (param):
##            param.Print()
##            param = par.Next()
##        print "---------------------------------------------"
##
##
##        ### iterate on the workspace element parameters
##        print "----------- Parameter postfit Workspace -------------"
##        parameters_workspace = workspace_postfit.allVars()
##        par = parameters_workspace.createIterator()
##        par.Reset()
##        param = par.Next()
##        while (param):
##            param.Print()
##            param = par.Next()
##        print "---------------------------------------------"
##        raw_input("zixu postfit")
##
##
##        workspace.data("data_obs_xww_%s_%s"%(self.channel, self.wtagger_category)).Print()
##
##        print "----------- Pdf in the Workspace -------------"
##        pdfs_workspace = workspace.allPdfs()
##        par = pdfs_workspace.createIterator()
##        par.Reset()
##        param = par.Next()
##        while (param):
##            param.Print()
##            param = par.Next()
##        print "----------------------------------------------"
##
##        rrv_x = workspace.var("rrv_mass_lvj")
##        data_obs = workspace.data("data_obs_xww_%s_%s"%(self.channel, self.wtagger_category))
##        ##if TString(self.signal_sample).Contains("RS1G_WW"):
##        ##    model_pdf_signal = workspace.pdf("RSWW_xww_%s_%s"%(self.channel, self.wtagger_category))
##        ##elif TString(self.signal_sample).Contains("BulkGraviton"):
##        ##    model_pdf_signal = workspace.pdf("BulkWW_xww_%s_%s"%(self.channel, self.wtagger_category))
##        ##elif TString(self.signal_sample).Contains("Wprime_WZ"):
##        ##    model_pdf_signal = workspace.pdf("WprimeWZ_xww_%s_%s"%(self.channel, self.wtagger_category))
##        ##else:
##        ##    model_pdf_signal = workspace.pdf("%s_xww_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category))
##        model_pdf_signal = workspace.pdf("%s_xww_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category))
##
##        model_pdf_WJets  = workspace.pdf("WJets_xww_%s_%s"%(self.channel, self.wtagger_category))
##        model_pdf_VV     = workspace.pdf("VV_xww_%s_%s"%(self.channel, self.wtagger_category))
##        model_pdf_TTbar  = workspace.pdf("TTbar_xww_%s_%s"%(self.channel, self.wtagger_category))
##        model_pdf_STop   = workspace.pdf("STop_xww_%s_%s"%(self.channel, self.wtagger_category))
##
##        model_pdf_signal.Print()
##        model_pdf_WJets.Print()
##        model_pdf_VV.Print()
##        model_pdf_TTbar.Print()
##        model_pdf_STop.Print()
##
##        ##if TString(self.signal_sample).Contains("RS1G_WW"):
##        ##    rrv_number_signal = workspace.var("rate_RSWW_xww_for_unbin")
##        ##elif TString(self.signal_sample).Contains("BulkGraviton"):
##        ##    rrv_number_signal = workspace.var("rate_BulkWW_xww_for_unbin")
##        ##elif TString(self.signal_sample).Contains("Wprime_WZ"):
##        ##    rrv_number_signal = workspace.var("rate_WprimeWZ_xww_for_unbin")
##        ##else:
##        ##    rrv_number_signal = workspace.var("rate_%s_xww_for_unbin"%(self.signal_sample))
##        rrv_number_signal = workspace.var("rate_%s_xww_for_unbin"%(self.signal_sample))
##
##
##        rrv_number_WJets  = workspace.var("rate_WJets_xww_for_unbin")
##        rrv_number_VV     = workspace.var("rate_VV_xww_for_unbin")
##        rrv_number_TTbar  = workspace.var("rate_TTbar_xww_for_unbin")
##        rrv_number_STop   = workspace.var("rate_STop_xww_for_unbin")
##
##        rrv_number_signal.Print()
##        rrv_number_WJets.Print()
##        rrv_number_VV.Print()
##        rrv_number_TTbar.Print()
##        rrv_number_STop.Print()
##
##        #### Prepare the final plot starting from total background
##        rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC_xww", "rrv_number_Total_background_MC_xww",
##                rrv_number_WJets.getVal()+
##                rrv_number_VV.getVal()+
##                rrv_number_TTbar.getVal()+
##                rrv_number_STop.getVal())
##
##        rrv_number_Total_background_MC.setError(TMath.Sqrt(
##            rrv_number_WJets.getError()* rrv_number_WJets.getError()+
##            rrv_number_VV.getError()* rrv_number_VV.getError()+
##            rrv_number_TTbar.getError()* rrv_number_TTbar.getError()+
##            rrv_number_STop.getError() *rrv_number_STop.getError()))
##
##        #### Total pdf
##        model_Total_background_MC = RooAddPdf("model_Total_background_MC_xww", "model_Total_background_MC_xww", RooArgList(model_pdf_WJets, model_pdf_VV, model_pdf_TTbar, model_pdf_STop), RooArgList(rrv_number_WJets, rrv_number_VV, rrv_number_TTbar, rrv_number_STop))
##
##        if data_obs.sumEntries() != 0:
##            #### scale factor in order to scale MC to data in the final plot -> in order to avoid the normalization to data which is done by default in rooFit
##            scale_number_signal = rrv_number_signal.getVal()/data_obs.sumEntries()
##            scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()/data_obs.sumEntries()
##        else:
##            scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()
##            scale_number_signal = rrv_number_signal.getVal()
##        #### create the frame
##        mplot = rrv_x.frame(RooFit.Title("check_workspace"), RooFit.Bins(int(rrv_x.getBins()/self.binwidth_narrow_factor)))
##        mplotP = rrv_x.frame(RooFit.Title("check_workspaceP"), RooFit.Bins(int(rrv_x.getBins()/self.binwidth_narrow_factor)))
##        data_obs.plotOn(mplot , RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0))
##        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("WJets"), RooFit.Components("WJets_xww_%s_%s,VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())
##
##        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("VV"), RooFit.Components("VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())
##
##        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("TTbar"), RooFit.Components("TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines())
##
##        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("STop"), RooFit.Components("STop_xww_%s_%s"%(self.channel, self.wtagger_category)), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())
##
##        #solid line
##        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("WJets_line_invisible"), RooFit.Components("WJets_xww_%s_%s,VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())
##
##        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("VV_line_invisible"), RooFit.Components("VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())
##
##        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("TTbar_line_invisible"), RooFit.Components("TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel, self.wtagger_category, self.channel, self.wtagger_category)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())
##
##        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("STop_line_invisible"), RooFit.Components("STop_xww_%s_%s"%(self.channel, self.wtagger_category)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())
##
##        model_pdf_signal.plotOn(mplot, RooFit.Normalization(scale_number_signal*self.signal_scale*self.signal_scale_plot), RooFit.Name("%s #times %s"%(self.signal_sample, self.signal_scale*self.signal_scale_plot)), RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["Signal"]), RooFit.LineStyle(2), RooFit.VLines())
##
##        #### plot the observed data using poissonian error bar
##        self.plot_data_with_poissoninterval(rrv_x, data_obs, mplot)
##
##        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Invisible())
##
##        mplot_pull = self.get_pull(rrv_x, mplot)
##
##        datahist = data_obs.binnedClone(data_obs.GetName()+"_binnedClone", data_obs.GetName()+"_binnedClone")
##        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
##        #draw_error_band(model_Total_background_MC, rrv_x.GetName(), rrv_number_Total_background_MC, self.FloatingParams, workspace , mplot, mplotP, datahist, self.color_palet["Uncertainty"], "F")
##        self.FloatingParams.Print()
##        print self.FloatingParams.getSize()
##        if self.FloatingParams.getSize() == 0:
##            if self.MODEL_4_mlvj == "ErfExp_v1" or self.MODEL_4_mlvj == "ErfPow_v1" or self.MODEL_4_mlvj == "2Exp" :
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))
##
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig4"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig5"%(self.channel, self.wtagger_category)))
##
##                if isTTbarFloating!=0:
##                    self.FloatingParams.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##                    self.FloatingParams.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
##                    self.FloatingParams.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
##
##            elif self.MODEL_4_mlvj == "ErfPow2_v1" or self.MODEL_4_mlvj == "ErfPowExp_v1" :
##
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig4"%(self.channel, self.wtagger_category)))
##
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig4"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig5"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig6"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig7"%(self.channel, self.wtagger_category)))
##
##                if isTTbarFloating!=0:
##                    self.FloatingParams.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##                    self.FloatingParams.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
##                    self.FloatingParams.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
##                    self.FloatingParams.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))
##
##            elif self.MODEL_4_mlvj == "Exp" or self.MODEL_4_mlvj == "Pow" :
##
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
##
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
##
##                if isTTbarFloating!=0:
##                    self.FloatingParams.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##
##            elif self.MODEL_4_mlvj == "ExpN" or self.MODEL_4_mlvj == "ExpTail" or self.MODEL_4_mlvj == "Pow2" :
##                print "Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)
##                workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)).Print()
##                raw_input("haha")
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
##
##
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig2"%(self.channel, self.wtagger_category)))
##                self.FloatingParams.add(workspace.var("Deco_WJets0_xww_sim_%s_%s_mlvj_eig3"%(self.channel, self.wtagger_category)))
##
##                if isTTbarFloating!=0:
##                    self.FloatingParams.add(workspace.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##
##                if isVVFloating!=0:
##                    self.FloatingParams.add(workspace.var("Deco_VV_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##                    self.FloatingParams.add(workspace.var("Deco_VV_xww_signal_region_%s_%s_mlvj_eig1"%(self.channel, self.wtagger_category)))
##
##                if isSTopFloating!=0:
##                    self.FloatingParams.add(workspace.var("Deco_STop_xww_signal_region_%s_%s_mlvj_eig0"%(self.channel, self.wtagger_category)))
##            else:
##                raw_input("this self.MODEL_4_mlvj:%s is not support!"%(self.MODEL_4_mlvj))
##            if workspace.var("rrv_mean_CB_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_category)):
##                if self.channel == "mu":
##                    self.FloatingParams.add(workspace.var("CMS_sig_p1_scale_m"))
##                    self.FloatingParams.add(workspace.var("CMS_sig_p2_scale_m"))
##                elif self.channel == "el":
##                    self.FloatingParams.add(workspace.var("CMS_sig_p1_scale_e"))
##                    self.FloatingParams.add(workspace.var("CMS_sig_p2_scale_e"))
##                elif self.channel == "em":
##                    self.FloatingParams.add(workspace.var("CMS_sig_p1_scale_em"))
##                    self.FloatingParams.add(workspace.var("CMS_sig_p2_scale_em"))
##
##                self.FloatingParams.add(workspace.var("CMS_sig_p1_jes"))
##                self.FloatingParams.add(workspace.var("CMS_sig_p2_jes"))
##                self.FloatingParams.add(workspace.var("CMS_sig_p2_jer"))
##        self.FloatingParams.Print()
##        #raw_input("zixu")
##        hdata = datahist.createHistogram(rrv_x.GetName(), int(rrv_x.getBins()/self.binwidth_narrow_factor))
##        draw_error_band(model_Total_background_MC, rrv_x.GetName(), rrv_number_Total_background_MC, self.FloatingParams, workspace , mplot, mplotP, hdata, self.color_palet["Uncertainty"], "F")
##
##        mplot.Print()
##        self.leg = self.legend4Plot(mplot, 0, 1,-0.01,-0.05, 0.11, 0.)
##        #self.leg.SetTextSize(0.036)
##        mplot.addObject(self.leg)
##        #pt1 = ROOT.TPaveText(0.6180905, 0.4355644, 0.8291457, 0.507992, "NDC")
##        #pt1.SetTextFont(42)
##        #pt1.SetTextSize(0.05)
##        #pt1.SetFillColor(0)
##        #pt1.SetFillStyle(0)
##        #pt1.SetBorderSize(0)
##        #text = pt1.AddText("")
##        #if options.category.find('Z') != -1: text = pt1.AddText("WZ category")
##        #elif options.category.find('W') != -1: text = pt1.AddText("WW category")
##        #text.SetTextFont(62)
##        #mplot.addObject(pt1)
##
##        mplot.GetYaxis().SetRangeUser(1e-2, mplot.GetMaximum()*1.2)
##
##        if workspace.var("rrv_num_floatparameter_in_last_fitting"):
##            self.nPar_float_in_fitTo = int(workspace.var("rrv_num_floatparameter_in_last_fitting").getVal())
##        else:
##            self.nPar_float_in_fitTo = self.FloatingParams.getSize()
##        nBinX = mplot.GetNbinsX()
##        ndof  = nBinX-self.nPar_float_in_fitTo
##        print "nPar = %s, chiSquare = %s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof)
##
##        parameters_list = RooArgList()
##        self.draw_canvas_with_pull2(rrv_x, datahist, mplot, mplotP, mplot_pull, ndof, parameters_list, "%s/m_lvj_fitting/"%(self.plotsDir), "check_postfit_workspace_for_limit", self.MODEL_4_mlvj, 0, 1)



    def combine_2workspaces(self, logy = 0):

        ### Taket the workspace for limits
        file1 = TFile(self.file_rlt_root.replace("_em_", "_mu_")) #mu
        file2 = TFile(self.file_rlt_root.replace("_em_", "_el_")) #el
        workspace1 = file1.Get("workspace4limit_")
        workspace2 = file2.Get("workspace4limit_")
        workspace1.Print()
        workspace2.Print()

        workspace_comb = RooWorkspace("workspace_comb", "workspace_comb")
        getattr(workspace_comb, "import")(workspace1.allPdfs())
        getattr(workspace_comb, "import")(workspace2.allPdfs())


        ##### iterate on the workspace element parameters
        ##print "----------- Parameter Workspace -------------"
        ##parameters_workspace = workspace.allVars()
        ##par = parameters_workspace.createIterator()
        ##par.Reset()
        ##param = par.Next()
        ##while (param):
        ##    param.Print()
        ##    param = par.Next()
        ##print "---------------------------------------------"

        ##workspace1.data("data_obs_xww_mu_%s"%(self.wtagger_category)).Print()
        ##workspace2.data("data_obs_xww_el_%s"%(self.wtagger_category)).Print()

        ##print "----------- Pdf in the Workspace -------------"
        ##pdfs_workspace = workspace.allPdfs()
        ##par = pdfs_workspace.createIterator()
        ##par.Reset()
        ##param = par.Next()
        ##while (param):
        ##    param.Print()
        ##    param = par.Next()
        ##print "----------------------------------------------"

        rrv_x = workspace_comb.var("rrv_mass_lvj")
        rrv_x.Print("v")
        data_obs  = workspace1.data("data_obs_xww_mu_%s"%(self.wtagger_category))
        data_obs2 = workspace2.data("data_obs_xww_el_%s"%(self.wtagger_category))
        data_obs.Print()
        data_obs2.Print()
        data_obs.append(data_obs2)
        data_obs.Print()
        #if TString(self.signal_sample).Contains("Bulk"):
        #    model_pdf_signal1 = workspace_comb.pdf("BulkWW_xww_mu_%s"%(self.wtagger_category))
        #    model_pdf_signal2 = workspace_comb.pdf("BulkWW_xww_el_%s"%(self.wtagger_category))
        #else:
        #    model_pdf_signal1 = workspace_comb.pdf("%s_xww_mu_%s"%(self.signal_sample, self.wtagger_category))
        #    model_pdf_signal2 = workspace_comb.pdf("%s_xww_el_%s"%(self.signal_sample, self.wtagger_category))
        model_pdf_signal1 = workspace_comb.pdf("%s_xww_mu_%s"%(self.signal_sample, self.wtagger_category))
        model_pdf_signal2 = workspace_comb.pdf("%s_xww_el_%s"%(self.signal_sample, self.wtagger_category))

        model_pdf_WJets1  = workspace_comb.pdf("WJets_xww_mu_%s"%(self.wtagger_category))
        model_pdf_VV1     = workspace_comb.pdf("VV_xww_mu_%s"%(self.wtagger_category))
        model_pdf_TTbar1  = workspace_comb.pdf("TTbar_xww_mu_%s"%(self.wtagger_category))
        model_pdf_STop1   = workspace_comb.pdf("STop_xww_mu_%s"%(self.wtagger_category))

        model_pdf_WJets2  = workspace_comb.pdf("WJets_xww_el_%s"%(self.wtagger_category))
        model_pdf_VV2     = workspace_comb.pdf("VV_xww_el_%s"%(self.wtagger_category))
        model_pdf_TTbar2  = workspace_comb.pdf("TTbar_xww_el_%s"%(self.wtagger_category))
        model_pdf_STop2   = workspace_comb.pdf("STop_xww_el_%s"%(self.wtagger_category))

        model_pdf_signal1.Print()
        model_pdf_WJets1.Print()
        model_pdf_VV1.Print()
        model_pdf_TTbar1.Print()
        model_pdf_STop1.Print()

        model_pdf_signal2.Print()
        model_pdf_WJets2.Print()
        model_pdf_VV2.Print()
        model_pdf_TTbar2.Print()
        model_pdf_STop2.Print()

        rrv_number_signal1 = workspace1.var("rate_%s_xww_for_unbin"%(self.signal_sample))
        rrv_number_signal2 = workspace2.var("rate_%s_xww_for_unbin"%(self.signal_sample))

        rrv_number_WJets1  = workspace1.var("rate_WJets_xww_for_unbin")
        rrv_number_VV1     = workspace1.var("rate_VV_xww_for_unbin")
        rrv_number_TTbar1  = workspace1.var("rate_TTbar_xww_for_unbin")
        rrv_number_STop1   = workspace1.var("rate_STop_xww_for_unbin")

        rrv_number_WJets2  = workspace2.var("rate_WJets_xww_for_unbin")
        rrv_number_VV2     = workspace2.var("rate_VV_xww_for_unbin")
        rrv_number_TTbar2  = workspace2.var("rate_TTbar_xww_for_unbin")
        rrv_number_STop2   = workspace2.var("rate_STop_xww_for_unbin")

        rrv_number_signal1.Print()
        rrv_number_WJets1.Print()
        rrv_number_VV1.Print()
        rrv_number_TTbar1.Print()
        rrv_number_STop1.Print()

        rrv_number_signal2.Print()
        rrv_number_WJets2.Print()
        rrv_number_VV2.Print()
        rrv_number_TTbar2.Print()
        rrv_number_STop2.Print()


        #### Prepare the final plot starting from total background
        rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC_xww", "rrv_number_Total_background_MC_xww",
                rrv_number_WJets1.getVal()+
                rrv_number_VV1.getVal()+
                rrv_number_TTbar1.getVal()+
                rrv_number_STop1.getVal()+
                rrv_number_WJets2.getVal()+
                rrv_number_VV2.getVal()+
                rrv_number_TTbar2.getVal()+
                rrv_number_STop2.getVal())

        rrv_number_Total_background_MC.setError(TMath.Sqrt(
            rrv_number_WJets1.getError()* rrv_number_WJets1.getError()+
            rrv_number_VV1.getError()* rrv_number_VV1.getError()+
            rrv_number_TTbar1.getError()* rrv_number_TTbar1.getError()+
            rrv_number_STop1.getError() *rrv_number_STop1.getError()+
            rrv_number_WJets2.getError()* rrv_number_WJets2.getError()+
            rrv_number_VV2.getError()* rrv_number_VV2.getError()+
            rrv_number_TTbar2.getError()* rrv_number_TTbar2.getError()+
            rrv_number_STop2.getError() *rrv_number_STop2.getError()))
        data_obs.Print()
        rrv_number_Total_background_MC.Print()

        #### Total pdf
        model_Total_background_MC = RooAddPdf("model_Total_background_MC_xww", "model_Total_background_MC_xww", RooArgList(model_pdf_WJets1, model_pdf_WJets2, model_pdf_VV1, model_pdf_VV2, model_pdf_TTbar1, model_pdf_TTbar2, model_pdf_STop1, model_pdf_STop2), RooArgList(rrv_number_WJets1, rrv_number_WJets2, rrv_number_VV1, rrv_number_VV2, rrv_number_TTbar1, rrv_number_TTbar2, rrv_number_STop1, rrv_number_STop2))

        scale_number_signal = (rrv_number_signal1.getVal() + rrv_number_signal2.getVal())/data_obs.sumEntries()
        #### scale factor in order to scale MC to data in the final plot -> in order to avoid the normalization to data which is done by default in RooFit
        scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()/data_obs.sumEntries()

        #### create the frame
        mplot = rrv_x.frame(RooFit.Title("check_workspace"), RooFit.Bins(int(rrv_x.getBins()/self.binwidth_narrow_factor)))
        mplotP = rrv_x.frame(RooFit.Title("check_workspaceP"), RooFit.Bins(int(rrv_x.getBins()/self.binwidth_narrow_factor)))

        data_obs.plotOn(mplot , RooFit.Name("data_invisible"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0))

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("WJets"), RooFit.Components("WJets_xww_*,VV_xww_*,TTbar_xww_*,STop_xww_*"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("VV"), RooFit.Components("VV_xww_*,TTbar_xww_*,STop_xww_*"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("TTbar"), RooFit.Components("TTbar_xww_*,STop_xww_*"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("STop"), RooFit.Components("STop_xww_*"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())

        #solid line
        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("WJets_line_invisible"), RooFit.Components("WJets_xww_*,VV_xww_*,TTbar_xww_*,STop_xww_*"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("VV_line_invisible"), RooFit.Components("VV_xww_*,TTbar_xww_*,STop_xww_*"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("TTbar_line_invisible"), RooFit.Components("TTbar_xww_*,STop_xww_*"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("STop_line_invisible"), RooFit.Components("STop_xww_*"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())


        model_pdf_signal = RooAddPdf("model_pdf_signal", "model_pdf_signal", RooArgList(model_pdf_signal1, model_pdf_signal2), RooArgList(rrv_number_signal1, rrv_number_signal2))
        model_pdf_signal.plotOn(mplot, RooFit.Normalization(scale_number_signal*self.signal_scale*self.signal_scale_plot), RooFit.Name("%s #times %s"%(self.signal_sample, self.signal_scale*self.signal_scale_plot)), RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["Signal"]), RooFit.LineStyle(2), RooFit.VLines())

        #### plot the observed data with poissonian error bar
        self.plot_data_with_poissoninterval(rrv_x, data_obs, mplot)

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Invisible())


        datahist = data_obs.binnedClone(data_obs.GetName()+"_binnedClone", data_obs.GetName()+"_binnedClone")
        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        ##workspace_comb = RooWorkspace("workspace_comb", "workspace_comb")
        ##getattr(workspace_comb, "import")(model_Total_background_MC)

        #for ExpN
        floatingparams = RooArgList("floatpara_list")
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sb_lo_from_fitting_el_%s_mlvj_eig0"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sb_lo_from_fitting_el_%s_mlvj_eig1"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sb_lo_from_fitting_el_%s_mlvj_eig2"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_el_%s_mlvj_eig0"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_el_%s_mlvj_eig1"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_el_%s_mlvj_eig2"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_el_%s_mlvj_eig3"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_TTbar_xww_signal_region_el_%s_mlvj_eig0"%(self.wtagger_category)))

        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sb_lo_from_fitting_mu_%s_mlvj_eig0"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sb_lo_from_fitting_mu_%s_mlvj_eig1"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sb_lo_from_fitting_mu_%s_mlvj_eig2"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_mu_%s_mlvj_eig0"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_mu_%s_mlvj_eig1"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_mu_%s_mlvj_eig2"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_mu_%s_mlvj_eig3"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_TTbar_xww_signal_region_mu_%s_mlvj_eig0"%(self.wtagger_category)))

        floatingparams.Print("v")
        if floatingparams.getSize() == 0:
            raw_input("zixu")

        hdata = datahist.createHistogram(rrv_x.GetName(), int(rrv_x.getBins()/self.binwidth_narrow_factor))
        draw_error_band(model_Total_background_MC, rrv_x.GetName(), rrv_number_Total_background_MC, floatingparams,
                workspace_comb , mplot, mplotP, hdata, self.color_palet["Uncertainty"], "F")


        mplot_pull = self.get_pull(rrv_x, mplot)

        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        #draw_error_band(model_Total_background_MC, rrv_x.GetName(), rrv_number_Total_background_MC, self.FloatingParams, workspace1 , mplot, self.color_palet["Uncertainty"], "F")

        mplot.Print()
        self.leg = self.legend4Plot(mplot, 0, 1,-0.01,-0.05, 0.11, 0.)
        #self.leg.SetTextSize(0.036)
        mplot.addObject(self.leg)

        mplot.GetYaxis().SetRangeUser(1e-2, mplot.GetMaximum()*1.2)


        parameters_list = RooArgList()

        self.nPar_float_in_fitTo = self.FloatingParams.getSize()
        #if workspace.var("rrv_num_floatparameter_in_last_fitting"):
        #    self.nPar_float_in_fitTo = int(workspace.var("rrv_num_floatparameter_in_last_fitting").getVal())
        #else:
        #    self.nPar_float_in_fitTo = self.FloatingParams.getSize()
        nBinX = mplot.GetNbinsX()
        ndof  = nBinX-self.nPar_float_in_fitTo
        print "nPar = %s, chiSquare = %s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof)


        #self.draw_canvas_with_pull2(rrv_x, datahist, mplot, mplotP, mplot_pull, ndof, parameters_list, "%s/m_lvj_fitting/"%(self.plotsDir), "check_workspace_for_limit", "", 0, 1)
        self.draw_canvas_with_pull2(rrv_x, datahist, mplot, mplotP, mplot_pull, ndof, parameters_list, "plot_em_prefit/m_lvj_fitting/", "check_workspace_for_limit", self.MODEL_4_mlvj, 0, 1)

        #if workspace1.var("rrv_num_floatparameter_in_last_fitting"):
        #    self.nPar_float_in_fitTo = int(workspace1.var("rrv_num_floatparameter_in_last_fitting").getVal())
        #else:
        #    self.nPar_float_in_fitTo = self.FloatingParams.getSize()
        #nBinX = mplot.GetNbinsX()
        #ndof  = nBinX-self.nPar_float_in_fitTo
        #print "nPar = %s, chiSquare = %s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof)

    ### in order to get the pull


    def read_postfit_workspaces(self, logy = 0):

        ### Taket the workspace for limits
        self.file_rlt_root_postfit = self.datacardsDir+"higgsCombinewwlvj_BulkGraviton_newxsec750_em_HP_unbin.Asymptotic.mH750.root"
        file_postfit = TFile(self.file_rlt_root_postfit)
        workspace_postfit = file_postfit.Get("w")
        file_postfit.Print()
        workspace_postfit.Print()
        #raw_input("zijun postfit 1")
        print workspace_postfit.loadSnapshot("clean")
        #raw_input("zijun9")
        file1 = TFile(self.file_rlt_root1)
        file2 = TFile(self.file_rlt_root2)
        workspace1 = file1.Get("workspace4limit_")
        workspace2 = file2.Get("workspace4limit_")

        #workspace_comb = RooWorkspace("workspace_comb", "workspace_comb")
        #getattr(workspace_comb, "import")(workspace1.allPdfs())
        #getattr(workspace_comb, "import")(workspace2.allPdfs())
        workspace_comb = workspace_postfit


        ### iterate on the workspace element parameters
        print "----------- Parameter Workspace postfit -------------"
        parameters_workspace = workspace_postfit.allVars()
        par = parameters_workspace.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            param.Print()
            param = par.Next()
        print "---------------------------------------------"


        print "----------- Parameter Workspace 1 -------------"
        parameters_workspace = workspace1.allVars()
        par = parameters_workspace.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            param.Print()
            param = par.Next()
        print "---------------------------------------------"
        print "----------- Parameter Workspace 2 -------------"
        parameters_workspace = workspace1.allVars()
        par = parameters_workspace.createIterator()
        par.Reset()
        param = par.Next()
        while (param):
            param.Print()
            param = par.Next()
        print "---------------------------------------------"
        #raw_input("zijun10")

        ##workspace1.data("data_obs_xww_mu_%s"%(self.wtagger_category)).Print()
        ##workspace2.data("data_obs_xww_el_%s"%(self.wtagger_category)).Print()

        ##print "----------- Pdf in the Workspace -------------"
        ##pdfs_workspace = workspace.allPdfs()
        ##par = pdfs_workspace.createIterator()
        ##par.Reset()
        ##param = par.Next()
        ##while (param):
        ##    param.Print()
        ##    param = par.Next()
        ##print "----------------------------------------------"

        rrv_x = workspace_comb.var("rrv_mass_lvj")
        #data_obs  = workspace_comb.data("data_obs")
        #data_obs.Print()


        data_obs  = workspace1.data("data_obs_xww_mu_%s"%(self.wtagger_category))
        data_obs2 = workspace2.data("data_obs_xww_el_%s"%(self.wtagger_category))
        data_obs.Print()
        data_obs2.Print()
        data_obs.append(data_obs2)
        data_obs.Print()


        #if TString(self.signal_sample).Contains("Bulk"):
        #    model_pdf_signal1 = workspace_comb.pdf("BulkWW_xww_mu_%s"%(self.wtagger_category))
        #    model_pdf_signal2 = workspace_comb.pdf("BulkWW_xww_el_%s"%(self.wtagger_category))
        #else:
        #    model_pdf_signal1 = workspace_comb.pdf("%s_xww_mu_%s"%(self.signal_sample, self.wtagger_category))
        #    model_pdf_signal2 = workspace_comb.pdf("%s_xww_el_%s"%(self.signal_sample, self.wtagger_category))
        model_pdf_signal1 = workspace_comb.pdf("shapeSig_BulkWW_xww_ch1")
        model_pdf_signal2 = workspace_comb.pdf("shapeSig_BulkWW_xww_ch2")

        model_pdf_WJets1  = workspace_comb.pdf("shapeBkg_WJets_xww_ch1")
        model_pdf_VV1     = workspace_comb.pdf("shapeBkg_VV_xww_ch1")
        model_pdf_TTbar1  = workspace_comb.pdf("shapeBkg_TTbar_xww_ch1")
        model_pdf_STop1   = workspace_comb.pdf("shapeBkg_STop_xww_ch1")
        model_pdf_WJets2  = workspace_comb.pdf("shapeBkg_WJets_xww_ch2")
        model_pdf_VV2     = workspace_comb.pdf("shapeBkg_VV_xww_ch2")
        model_pdf_TTbar2  = workspace_comb.pdf("shapeBkg_TTbar_xww_ch2")
        model_pdf_STop2   = workspace_comb.pdf("shapeBkg_STop_xww_ch2")

#        model_pdf_WJets1  = workspace_comb.pdf("WJets_xww_mu_%s"%(self.wtagger_category))
#        model_pdf_VV1     = workspace_comb.pdf("VV_xww_mu_%s"%(self.wtagger_category))
#        model_pdf_TTbar1  = workspace_comb.pdf("TTbar_xww_mu_%s"%(self.wtagger_category))
#        model_pdf_STop1   = workspace_comb.pdf("STop_xww_mu_%s"%(self.wtagger_category))
#
#        model_pdf_WJets2  = workspace_comb.pdf("WJets_xww_el_%s"%(self.wtagger_category))
#        model_pdf_VV2     = workspace_comb.pdf("VV_xww_el_%s"%(self.wtagger_category))
#        model_pdf_TTbar2  = workspace_comb.pdf("TTbar_xww_el_%s"%(self.wtagger_category))
#        model_pdf_STop2   = workspace_comb.pdf("STop_xww_el_%s"%(self.wtagger_category))

        model_pdf_signal1.Print()
        model_pdf_WJets1.Print()
        model_pdf_VV1.Print()
        model_pdf_TTbar1.Print()
        model_pdf_STop1.Print()

        model_pdf_signal2.Print()
        model_pdf_WJets2.Print()
        model_pdf_VV2.Print()
        model_pdf_TTbar2.Print()
        model_pdf_STop2.Print()
        #raw_input("PDF out")

        ##if TString(self.signal_sample).Contains("Bulk"):
        ##    rrv_number_signal1 = workspace1.var("rate_BulkWW_xww_for_unbin")
        ##    rrv_number_signal2 = workspace2.var("rate_BulkWW_xww_for_unbin")
        ##else:
        ##    rrv_number_signal1 = workspace1.var("rate_%s_xww_for_unbin"%(self.signal_sample))
        ##    rrv_number_signal2 = workspace2.var("rate_%s_xww_for_unbin"%(self.signal_sample))
        rrv_number_signal1 = workspace1.var("rate_%s_xww_for_unbin"%(self.signal_sample))
        rrv_number_signal2 = workspace2.var("rate_%s_xww_for_unbin"%(self.signal_sample))


        rrv_number_WJets1  = workspace1.var("rate_WJets_xww_for_unbin")
        rrv_number_VV1     = workspace1.var("rate_VV_xww_for_unbin")
        rrv_number_TTbar1  = workspace1.var("rate_TTbar_xww_for_unbin")
        rrv_number_STop1   = workspace1.var("rate_STop_xww_for_unbin")

        rrv_number_WJets2  = workspace2.var("rate_WJets_xww_for_unbin")
        rrv_number_VV2     = workspace2.var("rate_VV_xww_for_unbin")
        rrv_number_TTbar2  = workspace2.var("rate_TTbar_xww_for_unbin")
        rrv_number_STop2   = workspace2.var("rate_STop_xww_for_unbin")

        rrv_number_signal1.Print()
        rrv_number_WJets1.Print()
        rrv_number_VV1.Print()
        rrv_number_TTbar1.Print()
        rrv_number_STop1.Print()

        rrv_number_signal2.Print()
        rrv_number_WJets2.Print()
        rrv_number_VV2.Print()
        rrv_number_TTbar2.Print()
        rrv_number_STop2.Print()


        #### Prepare the final plot starting from total background
        rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC_xww", "rrv_number_Total_background_MC_xww",
                rrv_number_WJets1.getVal()+
                rrv_number_VV1.getVal()+
                rrv_number_TTbar1.getVal()+
                rrv_number_STop1.getVal()+
                rrv_number_WJets2.getVal()+
                rrv_number_VV2.getVal()+
                rrv_number_TTbar2.getVal()+
                rrv_number_STop2.getVal())

        rrv_number_Total_background_MC.setError(TMath.Sqrt(
            rrv_number_WJets1.getError()* rrv_number_WJets1.getError()+
            rrv_number_VV1.getError()* rrv_number_VV1.getError()+
            rrv_number_TTbar1.getError()* rrv_number_TTbar1.getError()+
            rrv_number_STop1.getError() *rrv_number_STop1.getError()+
            rrv_number_WJets2.getError()* rrv_number_WJets2.getError()+
            rrv_number_VV2.getError()* rrv_number_VV2.getError()+
            rrv_number_TTbar2.getError()* rrv_number_TTbar2.getError()+
            rrv_number_STop2.getError() *rrv_number_STop2.getError()))
        data_obs.Print()
        rrv_number_Total_background_MC.Print()

        #### Total pdf
        model_Total_background_MC = RooAddPdf("model_Total_background_MC_xww", "model_Total_background_MC_xww", RooArgList(model_pdf_WJets1, model_pdf_WJets2, model_pdf_VV1, model_pdf_VV2, model_pdf_TTbar1, model_pdf_TTbar2, model_pdf_STop1, model_pdf_STop2), RooArgList(rrv_number_WJets1, rrv_number_WJets2, rrv_number_VV1, rrv_number_VV2, rrv_number_TTbar1, rrv_number_TTbar2, rrv_number_STop1, rrv_number_STop2))

        scale_number_signal = (rrv_number_signal1.getVal() + rrv_number_signal2.getVal())/data_obs.sumEntries()
        #### scale factor in order to scale MC to data in the final plot -> in order to avoid the normalization to data which is done by default in rooFit
        scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()/data_obs.sumEntries()


        postfit_file = open("postfit_param_list.txt")
        while 1:
            line = postfit_file.readline()
            if not line:
                break
            param = []
            param_value = []
            param = line.split()

            param_name = param[0]
            param_value = float(param[2])
            param_error = float(param[4])

            workspace_comb.var(param_name).Print()

            workspace_comb.var(param_name).setVal(param_value)
            workspace_comb.var(param_name).setError(param_error)
            workspace_comb.var(param_name).Print()

        raw_input("zijun postfit")






        #### create the frame
        mplot = rrv_x.frame(RooFit.Title("check_workspace"), RooFit.Bins(int(rrv_x.getBins()/self.binwidth_narrow_factor)))
        mplotP = rrv_x.frame(RooFit.Title("check_workspaceP"), RooFit.Bins(int(rrv_x.getBins()/self.binwidth_narrow_factor)))

        data_obs.plotOn(mplot , RooFit.Name("data_invisible"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0))

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("WJets"), RooFit.Components("*WJets*,*VV*,*TTbar*,*STop*"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("VV"), RooFit.Components("*VV*,*TTbar*,*STop*"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("TTbar"), RooFit.Components("*TTbar*,*STop*"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("STop"), RooFit.Components("*STop*"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())

        #solid line
        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("WJets_line_invisible"), RooFit.Components("*WJets*,*VV*,*TTbar*,*STop*"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("VV_line_invisible"), RooFit.Components("*VV*,*TTbar*,*STop*"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("TTbar_line_invisible"), RooFit.Components("*TTbar*,*STop*"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Name("STop_line_invisible"), RooFit.Components("*STop*"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())


        model_pdf_signal = RooAddPdf("model_pdf_signal", "model_pdf_signal", RooArgList(model_pdf_signal1, model_pdf_signal2), RooArgList(rrv_number_signal1, rrv_number_signal2))
        model_pdf_signal.plotOn(mplot, RooFit.Normalization(scale_number_signal*self.signal_scale*self.signal_scale_plot), RooFit.Name("%s #times %s"%(self.signal_sample, self.signal_scale*self.signal_scale_plot)), RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["Signal"]), RooFit.LineStyle(2), RooFit.VLines())

        #### plot the observed data using poissonian error bar
        self.plot_data_with_poissoninterval(rrv_x, data_obs, mplot)

        model_Total_background_MC.plotOn(mplot, RooFit.Normalization(scale_number_Total_background_MC), RooFit.Invisible())


        datahist = data_obs.binnedClone(data_obs.GetName()+"_binnedClone", data_obs.GetName()+"_binnedClone")
        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        ##workspace_comb = RooWorkspace("workspace_comb", "workspace_comb")
        ##getattr(workspace_comb, "import")(model_Total_background_MC)


        floatingparams = RooArgList("floatpara_list")
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sb_lo_from_fitting_el_%s_mlvj_13TeV_eig0"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sb_lo_from_fitting_el_%s_mlvj_13TeV_eig1"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sb_lo_from_fitting_el_%s_mlvj_13TeV_eig2"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_el_%s_mlvj_13TeV_eig0"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_el_%s_mlvj_13TeV_eig1"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_el_%s_mlvj_13TeV_eig2"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_el_%s_mlvj_13TeV_eig3"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_TTbar_xww_signal_region_el_%s_mlvj_13TeV_eig0"%(self.wtagger_category)))

        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sb_lo_from_fitting_mu_%s_mlvj_13TeV_eig0"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sb_lo_from_fitting_mu_%s_mlvj_13TeV_eig1"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sb_lo_from_fitting_mu_%s_mlvj_13TeV_eig2"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_mu_%s_mlvj_13TeV_eig0"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_mu_%s_mlvj_13TeV_eig1"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_mu_%s_mlvj_13TeV_eig2"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_WJets0_xww_sim_mu_%s_mlvj_13TeV_eig3"%(self.wtagger_category)))
        floatingparams.add(workspace_comb.var("Deco_TTbar_xww_signal_region_mu_%s_mlvj_13TeV_eig0"%(self.wtagger_category)))

        floatingparams.Print("v")
        #raw_input("zixu")
        draw_error_band(model_Total_background_MC, rrv_x.GetName(), rrv_number_Total_background_MC, floatingparams, workspace_comb , mplot, mplotP, datahist, self.color_palet["Uncertainty"], "F")


        mplot_pull = self.get_pull(rrv_x, mplot)

        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        #draw_error_band(model_Total_background_MC, rrv_x.GetName(), rrv_number_Total_background_MC, self.FloatingParams, workspace1 , mplot, self.color_palet["Uncertainty"], "F")

        mplot.Print()
        self.leg = self.legend4Plot(mplot, 0, 1,-0.01,-0.05, 0.11, 0.)
        #self.leg.SetTextSize(0.036)
        mplot.addObject(self.leg)

        mplot.GetYaxis().SetRangeUser(1e-2, mplot.GetMaximum()*1.2)


        parameters_list = RooArgList()

        self.nPar_float_in_fitTo = self.FloatingParams.getSize()
        #if workspace.var("rrv_num_floatparameter_in_last_fitting"):
        #    self.nPar_float_in_fitTo = int(workspace.var("rrv_num_floatparameter_in_last_fitting").getVal())
        #else:
        #    self.nPar_float_in_fitTo = self.FloatingParams.getSize()
        nBinX = mplot.GetNbinsX()
        ndof  = nBinX-self.nPar_float_in_fitTo
        print "nPar = %s, chiSquare = %s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof)


        #self.draw_canvas_with_pull2(rrv_x, datahist, mplot, mplotP, mplot_pull, ndof, parameters_list, "%s/m_lvj_fitting/"%(self.plotsDir), "check_workspace_for_limit", "", 0, 1)
        self.draw_canvas_with_pull2(rrv_x, datahist, mplot, mplotP, mplot_pull, ndof, parameters_list, "plot_em_postfit/m_lvj_fitting/", "check_workspace_for_limit", self.MODEL_4_mlvj, 0, 1)

        #if workspace1.var("rrv_num_floatparameter_in_last_fitting"):
        #    self.nPar_float_in_fitTo = int(workspace1.var("rrv_num_floatparameter_in_last_fitting").getVal())
        #else:
        #    self.nPar_float_in_fitTo = self.FloatingParams.getSize()
        #nBinX = mplot.GetNbinsX()
        #ndof  = nBinX-self.nPar_float_in_fitTo
        #print "nPar = %s, chiSquare = %s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof)

    ### in order to get the pull
    def get_pull(self, rrv_x, mplot_orig):
        print "############### draw the pull plot ########################"
        hpull = mplot_orig.pullHist()
        x = ROOT.Double(0.)
        y = ROOT.Double(0.)
        pull_edge=4;
        for ipoint in range(0, hpull.GetN()):
            hpull.GetPoint(ipoint, x, y)
            if(y == 0):
                hpull.SetPoint(ipoint, x,-100)#remove from PULL plots
            if y>4 or y<-4:
                pull_edge=8
        gt = ROOT.TH1F("gt", "gt", int(rrv_x.getBins()/self.binwidth_narrow_factor), rrv_x.getMin(), rrv_x.getMax())
        gt.SetMinimum(-1*pull_edge+0.01)
        gt.SetMaximum(pull_edge-0.01)
        gt.SetDirectory(0)
        gt.SetStats(0)
        gt.SetLineStyle(0)
        gt.SetMarkerStyle(20)
        gt.GetXaxis().SetTitle(rrv_x.getTitle(1).Data())
        gt.GetXaxis().SetLabelFont(42)
        gt.GetXaxis().SetLabelOffset(0.02)
        gt.GetXaxis().SetLabelSize(0.15)
        gt.GetXaxis().SetTitleSize(0.15)
        gt.GetXaxis().SetTitleOffset(1.2)
        gt.GetXaxis().SetTitleFont(42)
        gt.GetYaxis().SetTitle("#frac{Data-Fit}{#sigma_{data}}")
        gt.GetYaxis().CenterTitle(True)
        gt.GetYaxis().SetNdivisions(205)
        gt.GetYaxis().SetLabelFont(42)
        gt.GetYaxis().SetLabelOffset(0.007)
        gt.GetYaxis().SetLabelSize(0.15)
        gt.GetYaxis().SetTitleSize(0.15)
        gt.GetYaxis().SetTitleOffset(0.35)
        gt.GetYaxis().SetTitleFont(42)
        hpull.SetHistogram(gt)
        return hpull

    def plot_data_with_poissoninterval(self, rrv_x, data_obs, mplot):
        """ adding obs data to the mplots, and poisson error"""
        datahist = data_obs.binnedClone(data_obs.GetName()+"_binnedClone", data_obs.GetName()+"_binnedClone")
        data_histo = datahist.createHistogram(rrv_x.GetName(), int(rrv_x.getBins()/self.binwidth_narrow_factor))
        data_histo.SetName("data")
        data4plot = RooHist(data_histo)
        data4plot.SetMarkerStyle(20)
        data4plot.SetMarkerSize(1.5)

        alpha = 0.3173 #1 - 0.6827
        for iPoint  in range(data4plot.GetN()):
            N = data4plot.GetY()[iPoint]
            if N == 0 : #not show mark for event=0 bins
                L = 0
                U = 0
            else:
                L = (ROOT.Math.gamma_quantile(alpha/2, N, 1.))
                U =  ROOT.Math.gamma_quantile_c(alpha/2, N+1, 1)
            data4plot.SetPointEYlow(iPoint, N-L)
            data4plot.SetPointEYhigh(iPoint, U-N)
            data4plot.SetPointEXlow(iPoint, 0)
            data4plot.SetPointEXhigh(iPoint, 0)
            #print L, N, U
        #raw_input("zixu")

        mplot.addPlotable(data4plot, "PE")

    def get_canvas(self, cname, isalpha = False):
        #tdrstyle.setTDRStyle()
        CMS_lumi.lumi_13TeV = "%s fb^{-1}"%(self.GetLumi())
        CMS_lumi.writeExtraText = 1
        CMS_lumi.extraText = "Preliminary"

        iPos = 11
        if(iPos == 0):
            CMS_lumi.relPosX = 0.15

        H_ref = 600
        W_ref = 800
        W = W_ref
        H  = H_ref

        T = 0.08*H_ref
        B = 0.12*H_ref
        L = 0.12*W_ref
        R = 0.06*W_ref

        canvas = ROOT.TCanvas(cname, cname,W,H)
        canvas.SetFillColor(0)
        canvas.SetBorderMode(0)
        canvas.SetFrameFillStyle(0)
        canvas.SetFrameBorderMode(0)
        canvas.SetLeftMargin(L/W)
        canvas.SetRightMargin(R/W)
        canvas.SetTopMargin(T/H)
        canvas.SetBottomMargin(B/H+0.03)
        canvas.SetTickx()
        canvas.SetTicky()
        if isalpha:
            canvas.SetTicky(0)

        return canvas


    #### in order to make the banner on the plots
    def banner4Plot(self, iswithpull = 0):
        #print "############### draw the banner ########################"
        if iswithpull:
            if self.channel == "el":
                banner = TLatex(0.3, 0.96, ("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 13TeV, W#rightarrow e #nu "%(self.GetLumi())))
            elif self.channel == "mu":
                banner = TLatex(0.3, 0.96, ("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 13TeV, W#rightarrow #mu #nu "%(self.GetLumi())))
            elif self.channel == "em":
                banner = TLatex(0.3, 0.96, ("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 13TeV, W#rightarrow #mu, e #nu "%(self.GetLumi())))
            banner.SetNDC()
            banner.SetTextSize(0.041)
        else:
            if self.channel == "el":
                banner = TLatex(0.22, 0.96, ("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 13TeV, W#rightarrow e #nu "%(self.GetLumi())))
            if self.channel == "mu":
                banner = TLatex(0.22, 0.96, ("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 13TeV, W#rightarrow #mu #nu "%(self.GetLumi())))
            if self.channel == "em":
                banner = TLatex(0.22, 0.96, ("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 13TeV, W#rightarrow #mu, e #nu "%(self.GetLumi())))
            banner.SetNDC()
            banner.SetTextSize(0.033)

        return banner

    ### in order to make the legend
    def legend4Plot(self, plot, left=1, isFill=1, x_offset_low=0., y_offset_low=0., x_offset_high = 0., y_offset_high = 0., TwoCoulum = 1., isalpha=False, ismj=False):
        print "############### draw the legend ########################"
        if left == -1:
            theLeg = TLegend(0.65+x_offset_low, 0.58+y_offset_low, 0.93+x_offset_low, 0.87+y_offset_low, "", "NDC")
            theLeg.SetName("theLegend")
            theLeg.SetLineColor(0)
            theLeg.SetTextFont(42)
            theLeg.SetTextSize(.04)
        else:
            theLeg = TLegend(0.37+x_offset_low, 0.50+y_offset_low, 0.72+x_offset_high, 0.82+y_offset_high, "", "NDC")
            #theLeg = TLegend(0.41+x_offset_low, 0.61+y_offset_low, 0.76+x_offset_high, 0.93+y_offset_high, "", "NDC")
            theLeg.SetName("theLegend")
            if ismj:
                theLeg = TLegend(0.3715365, 0.505, 0.8526448, 0.845, "", "NDC")
            if TwoCoulum :
                theLeg.SetNColumns(2)
            if isalpha:
                theLeg = TLegend(0.3944724, 0.4370629, 0.7650754, 0.8374126, "", "NDC")

        theLeg.SetFillColor(0)
        theLeg.SetFillStyle(0)
        theLeg.SetBorderSize(0)
        theLeg.SetLineColor(0)
        theLeg.SetLineWidth(0)
        theLeg.SetLineStyle(0)
        theLeg.SetTextSize(0.05)
        theLeg.SetTextFont(42)

        entryCnt = 0
        objName_before = ""
        objName_signal_graviton = ""
        objNameLeg_signal_graviton = ""

        if   self.channel == 'el':
            legHeader = "W#rightarrowe#nu"
        elif self.channel == 'mu':
            legHeader = "W#rightarrow#mu#nu"
        else:
            legHeader = "W#rightarrowl#nu"

        for obj in range(int(plot.numItems())):
            objName = plot.nameOf(obj)
            if objName == "errorband" :
                objName = "Uncertainty"
            #print objName
            if not (((plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty"))) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or objName == objName_before):
                theObj = plot.getObject(obj)
                objTitle = objName
                drawoption = plot.getDrawOptions(objName).Data()
                if drawoption == "P":
                    drawoption = "PE"
                if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):
                    objName_before = objName
                    continue
                elif TString(objName).Contains("Graph") :
                    objName_before = objName
                    continue
                elif TString(objName).Data() == "data" :
                    theLeg.AddEntry(theObj, "Data "+legHeader, "PE")
                    objName_before = objName
                else:
                    objName_before = objName
                    continue

        entryCnt = 0
        objName_before = ""
        objName_signal_graviton = ""
        objNameLeg_signal_graviton = ""

        for obj in range(int(plot.numItems())):
            objName = plot.nameOf(obj)
            if objName == "errorband" :
                objName = "Uncertainty"
            #print objName
            if not (((plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty"))) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or objName == objName_before):
                theObj = plot.getObject(obj)
                objTitle = objName
                drawoption = plot.getDrawOptions(objName).Data()
                if drawoption == "P":
                    drawoption = "PE"
                if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):
                    objName_before = objName
                    continue
                elif TString(objName).Contains("Graph") :
                    objName_before = objName
                    continue
                elif TString(objName).Data() == "WJets" :
                    theLeg.AddEntry(theObj, "W+jets", "F")
                    objName_before = objName
                else:
                    objName_before = objName
                    continue

        entryCnt = 0
        objName_before = ""
        objName_signal_graviton = ""
        objNameLeg_signal_graviton = ""


        for obj in range(int(plot.numItems())):
            objName = plot.nameOf(obj)
            if objName.find("TPave") != -1:
                continue
            if objName == "errorband" :
                objName = "Uncertainty"
            #print objName
            if not (((plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty"))) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or objName == objName_before):
                theObj = plot.getObject(obj)
                objTitle = objName
                drawoption = plot.getDrawOptions(objName).Data()
                if drawoption == "P":
                    drawoption = "PE"
                if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):
                    theLeg.AddEntry(theObj, objName, "F")
                elif TString(objName).Contains("Graph") :
                    if not (objName_before == "Graph" or objName_before == "Uncertainty"):
                        theLeg.AddEntry(theObj, "Uncertainty", "F")
                else:
                    if TString(objName).Data() == "STop" :
                        theLeg.AddEntry(theObj, "Single Top", "F")
                    elif TString(objName).Data() == "TTbar" :
                        theLeg.AddEntry(theObj, "t#bar{t}", "F")
                    elif TString(objName).Data() == "VV" :
                        theLeg.AddEntry(theObj, "WW/WZ", "F")
                    elif TString(objName).Data() == "data" :
                        objName_before = objName
                        entryCnt = entryCnt+1
                        continue
                    elif TString(objName).Data() == "WJets" :
                        objName_before = objName
                        entryCnt = entryCnt+1
                        continue
                    elif TString(objName).Contains("vbfH"):
                        theLeg.AddEntry(theObj, (TString(objName).ReplaceAll("vbfH", "qqH")).Data() , "L")
                    elif TString(objName).Contains("Uncertainty"):
                        theLeg.AddEntry(theObj, objTitle, drawoption)
                    elif TString(objName).Contains("Wprime"):
                        objName_signal_graviton = theObj
                        tmp_mass = int(filter(str.isdigit, (objName.split())[0]))/1000.0 #TeV, for example: "WprimeWZ800 #times 20" -> 800/1000 = 0.8
                        objNameLeg_signal_graviton = "W' M_{W}=%s TeV (#times%s)"%(tmp_mass, self.signal_scale*self.signal_scale_plot)
                    elif TString(objName).Contains("RS") or TString(objName).Contains("Bulk"):
                        prefix = ""
                        if TString(objName).Contains("RS"):
                            prefix = "RS"
                        elif TString(objName).Contains("Bulk"):
                            prefix = "Bulk"
                        objName_signal_graviton = theObj
                        tmp_mass = int(filter(str.isdigit, (objName.split())[0]))/1000.0 #TeV, for example: "WprimeWZ800 #times 20" -> 800/1000 = 0.8
                        objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=%s TeV (#times%s)"%(tmp_mass, self.signal_scale*self.signal_scale_plot)

                    else : theLeg.AddEntry(theObj, objTitle, drawoption)
                entryCnt = entryCnt+1
            objName_before = objName
        if objName_signal_graviton != "":
            theLeg.AddEntry(objName_signal_graviton, TString(objNameLeg_signal_graviton).Data() , "L")
        return theLeg

    #### jusr drawing canvas with no pull
    def draw_canvas(self, in_obj, in_directory, in_file_name, is_range = 0, logy = 0, frompull = 0, isalpha = 0):
        #print "############### draw the canvas without pull ########################"
        #cMassFit = TCanvas("cMassFit", "cMassFit", 600, 800)
        cMassFit = self.get_canvas("cMassFit", isalpha)

        if frompull and logy :
            in_obj.GetYaxis().SetRangeUser(1e-2, in_obj.GetMaximum()/200)
        elif not frompull and logy :
            in_obj.GetYaxis().SetRangeUser(0.00001, in_obj.GetMaximum())

        if is_range:
            h2 = TH2D("h2", "", 100, 400, 1400, 4, 0.00001, 4)
            h2.Draw()
            in_obj.Draw("same")
        else :
            in_obj.Draw()

        in_obj.GetXaxis().SetTitleSize(0.045)
        in_obj.GetXaxis().SetTitleOffset(1.15)
        in_obj.GetXaxis().SetLabelSize(0.04)

        in_obj.GetYaxis().SetTitleSize(0.055)
        if isalpha:
            in_obj.GetYaxis().SetTitleOffset(1.00)
        else:
            in_obj.GetYaxis().SetTitleOffset(1.40)
        in_obj.GetYaxis().SetLabelSize(0.04)

        self.leg.SetTextSize(0.031)

        banner = self.banner4Plot()
        banner.Draw()

        Directory = TString(in_directory+self.signal_sample+"/")
        if not Directory.EndsWith("/"):
            Directory = Directory.Append("/")
        if not os.path.isdir(Directory.Data()):
            os.system("mkdir -p "+Directory.Data())

        rlt_file = TString(Directory.Data()+in_file_name)
        if rlt_file.EndsWith(".root"):
            rlt_file.ReplaceAll(".root", "_rlt_without_pull_and_paramters.png")
        else:
            rlt_file.ReplaceAll(".root", "")
            rlt_file = rlt_file.Append(".png")

        cMassFit.SaveAs(rlt_file.Data())

        rlt_file.ReplaceAll(".png", ".pdf")
        cMassFit.SaveAs(rlt_file.Data())

        ##rlt_file.ReplaceAll(".pdf", ".root")
        ##cMassFit.SaveAs(rlt_file.Data())

        ##if logy:
        ##    in_obj.GetYaxis().SetRangeUser(1e-3, in_obj.GetMaximum()*200)
        ##    cMassFit.SetLogy()
        ##    cMassFit.Update()
        ##    rlt_file.ReplaceAll(".root", "_log.root")
        ##    cMassFit.SaveAs(rlt_file.Data())
        ##    rlt_file.ReplaceAll(".root", ".pdf")
        ##    cMassFit.SaveAs(rlt_file.Data())
        ##    rlt_file.ReplaceAll(".pdf", ".png")
        ##    cMassFit.SaveAs(rlt_file.Data())

        if logy:
            in_obj.GetYaxis().SetRangeUser(1e-3, in_obj.GetMaximum()*200)
            cMassFit.SetLogy()
            cMassFit.Update()
            rlt_file.ReplaceAll(".pdf", "_log.pdf")
            cMassFit.SaveAs(rlt_file.Data())
            rlt_file.ReplaceAll(".pdf", ".png")
            cMassFit.SaveAs(rlt_file.Data())

    #### draw canvas with plots with pull
    def draw_canvas_with_pull1(self, rrv_x, datahist, mplot, mplot_pull, ndof, parameters_list, in_directory, in_file_name, in_model_name = "", show_constant_parameter=0, logy=0, ismj=0, isPull=0):

        print "############### draw the canvas with pull ########################"
        mplot.GetXaxis().SetTitle("")
        mplot.GetYaxis().SetTitleSize(0.07)
        mplot.GetYaxis().SetTitleOffset(0.9)
        mplot.GetYaxis().SetLabelSize(0.06)
        mplot.GetXaxis().SetLabelSize(0)

        cMassFit = self.get_canvas("cMassFit")#TCanvas("cMassFit", "cMassFit", 600, 600)
        # if parameters_list is empty, don't draw pad3
        par_first = parameters_list.createIterator()
        par_first.Reset()
        param_first = par_first.Next()
        doParameterPlot = 0
        if param_first and doParameterPlot != 0:
            pad1 = TPad("pad1", "pad1", 0., 0. , 0.8, 0.24)
            pad2 = TPad("pad2", "pad2", 0., 0.24, 0.8, 1.)
            pad3 = TPad("pad3", "pad3", 0.8, 0., 1, 1)
            pad1.Draw()
            pad2.Draw()
            pad3.Draw()
        else:
            pad1 = TPad("pad1", "pad1", 0., 0. , 1, 0.30) #pad1 - pull
            pad2 = TPad("pad2", "pad2", 0., 0.3, 1., 1.) #pad0

        pad2.SetRightMargin(0.1)
        pad2.SetTopMargin(0.1)
        pad2.SetBottomMargin(0.0001)
        pad1.SetRightMargin(0.1)
        pad1.SetTopMargin(0)
        pad1.SetBottomMargin(0.4)
        pad1.Draw()
        pad2.Draw()
        pad2.cd()
        mplot.Draw()
        pad1.cd()
        mplot_pull.Draw("AP")

        mplot.Print("v")
        if mplot.FindObject("errorband") != 0:
            print " Pull  !!!! "
            mplot_pull.Print("v")

        medianLine = TLine(mplot.GetXaxis().GetXmin(), 0., mplot.GetXaxis().GetXmax(), 0)
        medianLine.SetLineWidth(2)
        medianLine.SetLineColor(kRed)
        medianLine.Draw()
        mplot_pull.Draw("Psame")

        if param_first and doParameterPlot != 0:
            pad3.cd()
            latex = TLatex()
            latex.SetTextSize(0.1)
            par = parameters_list.createIterator()
            par.Reset()
            param = par.Next()
            i = 0
            while param:
                if (not param.isConstant()) or show_constant_parameter:
                    param.Print()
                    icolor = 1#if a paramenter is constant, color is 2
                    if param.isConstant():
                        icolor = 2
                    latex.DrawLatex(0, 0.9-i*0.04, "#color[%s]{%s}"%(icolor, param.GetName()))
                    latex.DrawLatex(0, 0.9-i*0.04-0.02, " #color[%s]{%4.3e +/- %2.1e}"%(icolor, param.getVal(), param.getError()))
                    i = i+1
                param = par.Next()
        cMassFit.Update()
        pad2.cd()
        CMS_lumi.CMS_lumi(pad2, 4, 11)
        pad2.cd()
        pad2.Update()
        pad2.RedrawAxis()
        frame = pad2.GetFrame()
        frame.Draw()
        cMassFit.cd()
        cMassFit.Update()

        ## create the directory where store the plots
        Directory = TString(in_directory+self.signal_sample)
        if not Directory.EndsWith("/"):
            Directory = Directory.Append("/")
        if not os.path.isdir(Directory.Data()):
            os.system("mkdir -p "+Directory.Data())

        rlt_file = TString(Directory.Data()+in_file_name)
        if rlt_file.EndsWith(".root"):
            TString(in_model_name).ReplaceAll(".root", "")
            rlt_file.ReplaceAll(".root", "_"+in_model_name+"_with_pull.png")
        else:
            TString(in_model_name).ReplaceAll(".root", "")
            rlt_file.ReplaceAll(".root", "")
            rlt_file = rlt_file.Append("_"+in_model_name+"_with_pull.png")

        cMassFit.SaveAs(rlt_file.Data())

        rlt_file.ReplaceAll(".png", ".pdf")
        cMassFit.SaveAs(rlt_file.Data())

        ##rlt_file.ReplaceAll(".pdf", ".root")
        ##cMassFit.SaveAs(rlt_file.Data())

        ##string_file_name = TString(in_file_name)
        ##if string_file_name.EndsWith(".root"):
        ##    string_file_name.ReplaceAll(".root", "_"+in_model_name)
        ##else:
        ##    string_file_name.ReplaceAll(".root", "")
        ##    string_file_name.Append("_"+in_model_name)

        if logy:
            mplot.GetYaxis().SetRangeUser(0.002, mplot.GetMaximum()*200)
            pad2.SetLogy()
            pad2.Update()
            cMassFit.Update()
            ##rlt_file.ReplaceAll(".root", "_log.root")
            ##cMassFit.SaveAs(rlt_file.Data())
            ##rlt_file.ReplaceAll(".root", ".pdf")
            ##cMassFit.SaveAs(rlt_file.Data())
            rlt_file.ReplaceAll(".pdf", "_log.pdf")
            cMassFit.SaveAs(rlt_file.Data())
            rlt_file.ReplaceAll(".pdf", ".png")
            cMassFit.SaveAs(rlt_file.Data())

        #self.draw_canvas(mplot, in_directory, string_file_name.Data(), 0, logy, 1)

    #### draw canvas with plots with pull
    def draw_canvas_with_pull2(self, rrv_x, datahist, mplot, mplotP, mplot_pull, ndof, parameters_list, in_directory, in_file_name, in_model_name = "", show_constant_parameter=0, logy=0, ismj=0, isPull=0):
        # mplot + pull

        print "############### draw the canvas with pull ########################"
        chi2_ = calculate_chi2(datahist, rrv_x, mplot, ndof, ismj)
        mplot.GetXaxis().SetTitle("")
        mplot.GetYaxis().SetTitleSize(0.07)
        mplot.GetYaxis().SetTitleOffset(0.9)
        mplot.GetYaxis().SetLabelSize(0.06)
        mplot.GetXaxis().SetLabelSize(0)

        cMassFit = self.get_canvas("cMassFit")
        # if parameters_list is empty, don't draw pad3
        par_first = parameters_list.createIterator()
        par_first.Reset()
        param_first = par_first.Next()
        doParameterPlot = 0
        if param_first and doParameterPlot != 0:
            pad1 = TPad("pad1", "pad1", 0., 0. , 0.8, 0.24)
            pad2 = TPad("pad2", "pad2", 0., 0.24, 0.8, 1.)
            pad3 = TPad("pad3", "pad3", 0.8, 0., 1, 1)
            pad1.Draw()
            pad2.Draw()
            pad3.Draw()
        else:
            pad1 = TPad("pad1", "pad1", 0., 0. , 1, 0.30) #pad1 - pull
            pad2 = TPad("pad2", "pad2", 0., 0.3, 1., 1.) #pad0
            pad2.SetRightMargin(0.1)
            pad2.SetTopMargin(0.1)
            pad2.SetBottomMargin(0.0001)
            pad1.SetRightMargin(0.1)
            pad1.SetTopMargin(0)
            pad1.SetBottomMargin(0.4)
            pad1.Draw()
            pad2.Draw()

        #draw the main plot
        pad2.cd()

        if ismj:
            pt = ROOT.TPaveText(0.6243719, 0.4080919, 0.8756281, 0.547952, "NDC")
            pt.SetTextSize(0.03746254)
        else:
            pt = ROOT.TPaveText(0.5175879, 0.7152847, 0.8027638, 0.8551449, "NDC")
            pt.SetTextSize(0.054)

        pt.SetTextFont(62)
        pt.SetTextAlign(12)
        pt.SetFillColor(0)
        pt.SetBorderSize(0)
        pt.SetFillStyle(0)
        text = pt.AddText("#chi^2/d.o.f = %.2f/%i = %.2f" %(chi2_[0], chi2_[1], chi2_[0]/chi2_[1]))
        text.SetTextFont(62)
        mplot.Draw()
        #pt.Draw()
        #banner = self.banner4Plot(1)
        #banner.Draw()

        #Draw the PULL plot
        pad1.cd()
        mplot_pull.Draw("AP")
        mplotP.Draw("same")
        medianLine = TLine(mplot.GetXaxis().GetXmin(), 0., mplot.GetXaxis().GetXmax(), 0)
        medianLine.SetLineWidth(2)
        medianLine.SetLineColor(kRed)
        medianLine.Draw()
        mplot_pull.Draw("Psame")

        if param_first and doParameterPlot != 0:
            pad3.cd()
            latex = TLatex()
            latex.SetTextSize(0.1)
            par = parameters_list.createIterator()
            par.Reset()
            param = par.Next()
            i = 0
            while param:
                if (not param.isConstant()) or show_constant_parameter:
                    param.Print()
                    icolor = 1#if a paramenter is constant, color is 2
                    if param.isConstant():
                        icolor = 2
                    latex.DrawLatex(0, 0.9-i*0.04, "#color[%s]{%s}"%(icolor, param.GetName()))
                    latex.DrawLatex(0, 0.9-i*0.04-0.02, " #color[%s]{%4.3e +/- %2.1e}"%(icolor, param.getVal(), param.getError()))
                    i = i+1
                param = par.Next()

        cMassFit.Update()
        pad2.cd()
        CMS_lumi.CMS_lumi(pad2, 4, 11)
        pad2.cd()
        pad2.Update()
        pad2.RedrawAxis()
        frame = pad2.GetFrame()
        frame.Draw()
        cMassFit.cd()
        cMassFit.Update()
        ## create the directory where store the plots
        Directory = TString(in_directory+self.signal_sample)
        if not Directory.EndsWith("/"):
            Directory = Directory.Append("/")
        if not os.path.isdir(Directory.Data()):
            os.system("mkdir -p "+Directory.Data())

        rlt_file = TString(Directory.Data()+in_file_name)
        if rlt_file.EndsWith(".root"):
            TString(in_model_name).ReplaceAll(".root", "")
            rlt_file.ReplaceAll(".root", "_"+in_model_name+"_with_pull.png")
        else:
            TString(in_model_name).ReplaceAll(".root", "")
            rlt_file.ReplaceAll(".root", "")
            rlt_file = rlt_file.Append("_"+in_model_name+"_with_pull.png")

        cMassFit.SaveAs(rlt_file.Data())

        rlt_file.ReplaceAll(".png", ".pdf")
        cMassFit.SaveAs(rlt_file.Data())

        #rlt_file.ReplaceAll(".pdf", ".root")
        #cMassFit.SaveAs(rlt_file.Data())

        #string_file_name = TString(in_file_name)
        #if string_file_name.EndsWith(".root"):
        #    string_file_name.ReplaceAll(".root", "_"+in_model_name)
        #else:
        #    string_file_name.ReplaceAll(".root", "")
        #    string_file_name.Append("_"+in_model_name)

        if logy:
            mplot.GetYaxis().SetRangeUser(0.002, mplot.GetMaximum()*200)
            pad2.SetLogy()
            pad2.Update()
            cMassFit.Update()
            #rlt_file.ReplaceAll(".root", "_log.root")
            #cMassFit.SaveAs(rlt_file.Data())
            #rlt_file.ReplaceAll(".root", ".pdf")
            #cMassFit.SaveAs(rlt_file.Data())
            rlt_file.ReplaceAll(".pdf", "_log.pdf")
            cMassFit.SaveAs(rlt_file.Data())
            rlt_file.ReplaceAll(".pdf", ".png")
            cMassFit.SaveAs(rlt_file.Data())

        #self.draw_canvas(mplot, in_directory, string_file_name.Data(), 0, logy, 1)


    ##### Get Lumi for banner title
    def GetLumi(self):
        return self.Lumi
        #if self.channel == "el":   return 2.6
        #elif self.channel == "mu": return 2.6
        #elif self.channel == "em": return 2.6

    #### function to run the selection on data to build the datasets
    def get_data(self):
        print "############### get_data ########################"
        self.get_mj_and_mlvj_dataset(self.file_data, "_data_xww", "massVhadJEC")
        getattr(self.workspace4limit_, "import")(self.workspace4fit_.var("rrv_number_dataset_signal_region_data_xww_%s_mlvj"%(self.channel)).clone("observation_for_counting_xww"))

    #### Define the steps to fit signal distribution in the mj and mlvj spectra
    def fit_Signal(self, model_narrow = "DoubleCB_mlvj", model_width = "BWDoubleCB"):
        print "############# fit_Signal #################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_signal, "_%s_xww"%(self.signal_sample), "massVhadJEC")# to get the shape of m_lvj
        self.fit_mlvj_model_single_MC(self.file_signal, "_%s_xww"%(self.signal_sample), "_signal_region", model_narrow, 0, 0, 0, 1)
        #raw_input("signal shape zixu")
        print "________________________________________________________________________"

    ##### Define the steps to fit WJets MC in the mj and mlvj spectra
    def fit_WJets(self):
        print "######################### fit_WJets ########################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_WJets0_mc, "_WJets0_xww", "massVhadJEC")# to get the shape of m_lvj
        self.get_mj_and_mlvj_dataset(self.file_WJets0_mc, "_WJets01_xww", "massVhadJEC")# to get the shape of m_lvj

        ### Fit in mj depends on the mlvj lower limit -> fitting the turn on at low mass or not
        if self.MODEL_4_mlvj == "ErfPowExp_v1" or self.MODEL_4_mlvj == "ErfPow2_v1" or self.MODEL_4_mlvj == "ErfExp_v1" :
            self.fit_mj_single_MC(self.file_WJets0_mc, "_WJets0_xww", "ErfExp")
            self.fit_mj_single_MC(self.file_WJets0_mc, "_WJets01_xww", "User1")
        else:
            self.fit_mj_single_MC(self.file_WJets0_mc, "_WJets0_xww", "User1")
            self.fit_mj_single_MC(self.file_WJets0_mc, "_WJets01_xww", "ErfExp")

        #### Fit the mlvj in sb_lo, signal region using two different model as done in the mj
        self.fit_mlvj_model_single_MC(self.file_WJets0_mc, "_WJets0_xww", "_sb_lo", self.MODEL_4_mlvj, 0, 0, 1, 1)
        #raw_input("ENTER")
        self.fit_mlvj_model_single_MC(self.file_WJets0_mc, "_WJets0_xww", "_signal_region", self.MODEL_4_mlvj, 0, 0, 1, 1)
        #raw_input("ENTER")
        self.fit_mlvj_model_single_MC(self.file_WJets0_mc, "_WJets01_xww", "_sb_lo", self.MODEL_4_mlvj_alter, 0, 0, 1, 1)
        #raw_input("ENTER")
        self.fit_mlvj_model_single_MC(self.file_WJets0_mc, "_WJets01_xww", "_signal_region", self.MODEL_4_mlvj_alter, 0, 0, 1, 1)
        #raw_input("ENTER")

        print "________________________________________________________________________"


    ##### Define the steps to fit VV MC in the mj and mlvj spectra
    def fit_VV(self):
        print "############################# fit_VV ################################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_VV_mc, "_VV_xww", "massVhadJEC")

        self.fit_mj_single_MC(self.file_VV_mc, "_VV_xww", "2_2Gaus")

        if self.MODEL_4_mlvj == "ErfPowExp_v1" or self.MODEL_4_mlvj == "ErfPow2_v1" or self.MODEL_4_mlvj == "ErfExp_v1":
            self.fit_mlvj_model_single_MC(self.file_VV_mc, "_VV_xww", "_sb_lo", "ErfExp_v1", 0, 0, 1)
            self.fit_mlvj_model_single_MC(self.file_VV_mc, "_VV_xww", "_signal_region", self.MODEL_4_mlvj, 1, 0, 1)

        else:
            self.fit_mlvj_model_single_MC(self.file_VV_mc, "_VV_xww", "_sb_lo", self.MODEL_4_mlvj, 0, 0, 1)
            self.fit_mlvj_model_single_MC(self.file_VV_mc, "_VV_xww", "_signal_region", self.MODEL_4_mlvj, 1, 0, 1)

        print "________________________________________________________________________"

    ##### Define the steps to fit TTbar MC in the mj and mlvj spectra
    def fit_TTbar(self):
        print "################################ fit_TTbar #########################################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_TTbar_mc, "_TTbar_xww", "massVhadJEC")# to get the shape of m_lvj

        if self.MODEL_4_mlvj == "ErfPowExp_v1" or self.MODEL_4_mlvj == "ErfPow2_v1" or self.MODEL_4_mlvj == "ErfExp_v1":
            if self.wtagger_category == "LP":
                self.fit_mj_single_MC(self.file_TTbar_mc, "_TTbar_xww", "ExpGaus")
            else:
                self.fit_mj_single_MC(self.file_TTbar_mc, "_TTbar_xww", "2Gaus_ErfExp")
        else:
            if self.wtagger_category == "LP" :
                self.fit_mj_single_MC(self.file_TTbar_mc, "_TTbar_xww", "ExpGaus")
            else:
                self.fit_mj_single_MC(self.file_TTbar_mc, "_TTbar_xww", "2Gaus_ErfExp")

        if self.MODEL_4_mlvj == "ErfPowExp_v1" or self.MODEL_4_mlvj == "ErfPow2_v1" or self.MODEL_4_mlvj == "ErfExp_v1" :
            #self.fit_mlvj_model_single_MC(self.file_TTbar_mc, "_TTbar_xww", "_sb_lo", "ErfExp_v1", 0, 0, 1)
            #self.fit_mlvj_model_single_MC(self.file_TTbar_mc, "_TTbar_xww", "_sb_lo", self.MODEL_4_mlvj, 0, 0, 1)
            #self.fit_mlvj_model_single_MC(self.file_TTbar_mc, "_TTbar_xww", "_signal_region", self.MODEL_4_mlvj, 1, 0, 1)
            self.fit_mlvj_model_single_MC(self.file_TTbar_mc, "_TTbar_xww",        "_sb_lo", "ErfPow2_v1", 0, 0, 1)
            self.fit_mlvj_model_single_MC(self.file_TTbar_mc, "_TTbar_xww", "_signal_region", "ErfPow2_v1", 1, 0, 1)

        else:
            self.fit_mlvj_model_single_MC(self.file_TTbar_mc, "_TTbar_xww", "_sb_lo", self.MODEL_4_mlvj)
            self.fit_mlvj_model_single_MC(self.file_TTbar_mc, "_TTbar_xww", "_signal_region", self.MODEL_4_mlvj, 1, 0, 1)

        print "________________________________________________________________________"


    #### Define the steps to fit STop MC in the mj and mlvj spectra
    def fit_STop(self):
        print "############################## fit_STop  #################################"
        self.get_mj_and_mlvj_dataset(self.file_STop_mc, "_STop_xww", "massVhadJEC")
        #self.fit_mj_single_MC(self.file_STop_mc, "_STop_xww", "ExpGaus")
        self.fit_mj_single_MC(self.file_STop_mc, "_STop_xww", "2Gaus_ErfExp")

        if self.MODEL_4_mlvj == "ErfPowExp_v1" or self.MODEL_4_mlvj == "ErfPow2_v1" or self.MODEL_4_mlvj == "ErfExp_v1":
            self.fit_mlvj_model_single_MC(self.file_STop_mc, "_STop_xww", "_sb_lo", "ErfExp_v1", 0, 0, 1)
            self.fit_mlvj_model_single_MC(self.file_STop_mc, "_STop_xww", "_signal_region", "ErfExp_v1", 1, 0, 1)
        else:
            self.fit_mlvj_model_single_MC(self.file_STop_mc, "_STop_xww", "_sb_lo", "Exp", 0, 0, 1)
            self.fit_mlvj_model_single_MC(self.file_STop_mc, "_STop_xww", "_signal_region", self.MODEL_4_mlvj, 1, 0, 1)

        print "________________________________________________________________________"

    ##### Fit of all the MC in both mj and mlvj : Signal, TTbar, STop, VV and Wjets
    def fit_AllSamples_Mj_and_Mlvj(self):
        print "################### fit_AllSamples_Mj_and_Mlvj #####################"
        self.fit_Signal()
        self.fit_WJets()
        self.fit_TTbar()
        self.fit_VV()
        self.fit_STop()
        print "________________________________________________________________________"


    ##### Analysis with sideband alpha correction
    def analysis_sideband_correction_method1(self):
        print "##################### Start sideband correction full analysis ##############"
        ### Fit all MC components in both mj and mlvj
        self.fit_AllSamples_Mj_and_Mlvj()
        ### take the real data
        self.get_data()
        ### fit the WJets Normalization into the signal region -> no jet mass fluctuation has been done
        self.fit_WJetsNorm()
        ### fit data in the mlvj low sideband with two different models
        self.fit_mlvj_in_Mj_sideband("_WJets01_xww", "_sb_lo", self.MODEL_4_mlvj_alter, 1)
        self.fit_mlvj_in_Mj_sideband("_WJets0_xww", "_sb_lo", self.MODEL_4_mlvj, 1)

        ### Prepare the workspace and datacards
        self.prepare_limit("sideband_correction_method1", 1, 0, 0)
        ### finale plot and check of the workspace
        self.read_workspace(1)

    ##### Analysis with no shape uncertainty on alpha
    def analysis_sideband_correction_method1_without_shape_and_psmodel_systermatic(self):
        #### fit all the MC samples
        self.fit_AllSamples_Mj_and_Mlvj()
        #### take the real data
        self.get_data()
        #### fit WJets just with one shape parametrization
        self.fit_WJetsNormalization_in_Mj_signal_region("_WJets0_xww")
        #self.fit_WJetsNormalization_in_Mj_signal_region("_WJets0")
        #### fit sb lo with just one parametrization
        self.fit_mlvj_in_Mj_sideband("_WJets0_xww", "_sb_lo", self.MODEL_4_mlvj, 1)
        #### prepare limit
        self.prepare_limit("sideband_correction_method1", 1, 0, 0)
        #### read the workspace
        self.read_workspace(1)

   ###### Analysis fitting just signal lineshape in mlvj
    def fit_signal_only(self):
        #### fit signal MC samples
        self.fit_AllSamples_Mj_and_Mlvj()

def calculate_chi2(hist, rrv_x, mplot_orig, ndof, ismj):
    pulls = array('d',[])
    print "############### calculate chi2 (new) ########################"
    hpull = mplot_orig.pullHist()
    #bins = 0
    bins_ = 0
    x = ROOT.Double(0.)
    y = ROOT.Double(0)
    for ipoint in range(0, hpull.GetN()):
        hpull.GetPoint(ipoint, x, y)
    hist.get(bins_)
    hist.weightError(RooAbsData.SumW2)
    print x, y, bins_, hist.get(bins_).getRealValue(rrv_x.GetName()), hist.weight(), hist.weightError(RooAbsData.SumW2)
    if hist.weight() != 0:
        pulls.append(y)
    #print x, y, hist.GetBinCenter(bins_), hist.GetBinContent(bins_)
        #if not(ismj) and y != 0 and TMath.Abs(y) < 4: pulls.append(y)
    #elif ismj and y != 0 and (x < 65 or x > 135):
    # pulls.append(y)
        # bins+=1
    #else: print "Bin %f is empty!" %x
    bins_+=1
    chi2 = 0
    for p in pulls:
        chi2+=(p*p)

    #if ismj:
    # ndof = ndof - (hpull.GetN()-bins)
    # print hpull.GetN(), bins, ndof
    print "Chi2/ndof = %f/%f = %f" %(chi2, ndof, chi2/ndof)
    return chi2, ndof

### print the parameters of a given pdf --> only non constant ones
def ShowParam_Pdf(model_pdf, argset_notparameter):
    print "########### Show Parameters of a input model  ############"
    model_pdf.Print()
    parameters = model_pdf.getParameters(argset_notparameter)
    par = parameters.createIterator()
    par.Reset()
    param = par.Next()
    while (param):
        if not param.isConstant():
            param.Print()
            if (param.getVal()-param.getMin())< param.getError()*1 or (param.getMax()- param.getVal())< param.getError()*1:
                param.Print()
        param = par.Next()

### fix a pdf in a different way --> for RooAbsPdf
def fix_Pdf(model_pdf, argset_notparameter):
    print "########### Fixing a RooAbsPdf for mlvj or mj  ############"
    parameters = model_pdf.getParameters(argset_notparameter)
    par = parameters.createIterator()
    par.Reset()
    param = par.Next()
    while (param):
        param.setConstant(kTRUE)
        param.Print()
        param = par.Next()



def pre_limit_sb_correction_without_systermatic(channel, signal_sample,
        in_mlvj_signal_region_min = 500, in_mlvj_signal_region_max = 700,
        in_mj_min=30, in_mj_max=140, in_mlvj_min=700, in_mlvj_max=5000,
        fit_model = "ErfExp_v1", fit_model_alter = "ErfPow_v1"):
    """funtion to run the analysis with less systematics"""

    print "#################### pre_limit_sb_correction_without_systermatic: channel %s, signal %s, max and min signal region %f-%f, max and min mJ %f-%f, max and min mlvj %f-%f, fit model %s and alternate %s ######################"%(channel, signal_sample, in_mlvj_signal_region_min, in_mlvj_signal_region_max, in_mj_min, in_mj_max, in_mlvj_min, in_mlvj_max, fit_model, fit_model_alter)
    boostedW_fitter = DoFit(channel, signal_sample, in_mlvj_signal_region_min, in_mlvj_signal_region_max, in_mj_min, in_mj_max, in_mlvj_min, in_mlvj_max, fit_model, fit_model_alter)
    boostedW_fitter.analysis_sideband_correction_method1_without_shape_and_psmodel_systermatic()



def pre_limit_sb_correction(method, channel, signal_sample = "BulkG_c0p2_M1000",
        in_mlvj_signal_region_min=500, in_mlvj_signal_region_max=700,
        in_mj_min=30, in_mj_max=140, in_mlvj_min=700, in_mlvj_max=5000,
        fit_model = "ErfExp_v1", fit_model_alter = "ErfPow_v1"):
    """funtion to run the analysis with full systematics"""

    print "#################### pre_limit_sb_correction: channel %s, signal %s, max and min signal region %f-%f, max and min mJ %f-%f, max and min mlvj %f-%f, fit model %s and alternate %s ######################"%(channel, signal_sample, in_mlvj_signal_region_min, in_mlvj_signal_region_max, in_mj_min, in_mj_max, in_mlvj_min, in_mlvj_max, fit_model, fit_model_alter)

    boostedW_fitter = DoFit(channel, signal_sample, in_mlvj_signal_region_min, in_mlvj_signal_region_max, in_mj_min, in_mj_max, in_mlvj_min, in_mlvj_max, fit_model, fit_model_alter)
    getattr(boostedW_fitter, "analysis_sideband_correction_%s"%(method))()

### funtion to run the analysis without systematic
def pre_limit_simple(channel):
    """quickly do limit for one signal point"""
    print "######################### pre_limit_simple for %s sampel"%(channel)

    pre_limit_sb_correction_without_systermatic(channel, "BulkGravWW750", 650, 850, 40, 150, 600, 1500, "ExpN", "ExpTail")
    #pre_limit_sb_correction_without_systermatic(channel, "BulkGravWW750", 650, 850, 40, 150, 600, 1500, "ExpTail", "ExpN")
    #pre_limit_sb_correction_without_systermatic(channel, "BulkGravWW4000", 3900, 4100, 40, 150, 800, 5000, "ExpN", "Exp")
    #pre_limit_sb_correction_without_systermatic(channel, "BulkGravWW4500", 4400, 4600, 40, 150, 800, 5000, "ExpTail", "ExpN")
    #pre_limit_sb_correction_without_systermatic(channel, "BulkGravWW4000", 3900, 4100, 40, 150, 800, 5000, "ExpN", "ExpTail")
    #pre_limit_sb_correction_without_systermatic(channel, "WprimeWZ1000", 900, 1100, 40, 150, 600, 1500, "ExpN", "ExpTail")

def fit_signal(method, channel, signal_sample,  in_mlvj_min, in_mlvj_max): # the WJets M_lvj shape and normalization are from sb_correction
    """only fit Signal"""
    boostedW_fitter = DoFit(channel, signal_sample, in_mlvj_min, in_mlvj_max, 40, 150, in_mlvj_min, in_mlvj_max)
    boostedW_fitter.fit_Signal()

def control_single(channel):
    """quickly draw control plots"""
    print "control_single for %s sampel"%(channel)
    boostedW_fitter = DoFit(channel, "BulkGravWW750")
    boostedW_fitter.ControlPlots()


### function to check the workspace once it has already created
def check_workspace(channel, signal):
    """read workspace of el or mu channel, and draw the M_lvj plot"""
    boostedW_fitter = DoFit(channel, signal, 500, 700, 40, 140, 400, 1400, "ExpN", "ExpTail")#all the param is useless. only the pre-fit will be included in the picture name
    boostedW_fitter.read_workspace()
    #boostedW_fitter.read_workspace_postfit()

def combine_el_mu(channel):
    """read workspaces of el and mu channel, and draw the combined M_lvj plot"""
    #boostedW_fitter = DoFit("em", "BulkGravWW750", 650, 850, 40, 150, 600, 1500, "ExpN", "Pow")
    #boostedW_fitter = DoFit("em", "BulkGravWW4500", 4400, 4600, 40, 150, 800, 5000, "ExpN", "Pow")
    #boostedW_fitter = DoFit("em", "WprimeWZ-HVT-A800", 700, 900, 40, 150, 600, 1500, "ExpN", "Pow")
    boostedW_fitter = DoFit("em", "WprimeWZ-HVT-A4500", 4400, 4600, 40, 150, 800, 5000, "ExpN", "Pow")
    boostedW_fitter.combine_2workspaces(1)
    #boostedW_fitter.read_postfit_workspaces()

#### Main Code
if __name__ == '__main__':

    #if options.fitwtagger:
    #    print 'fitwtagger for %s sample'%(options.opt_channel)
    #    control_sample(options.opt_channel)#mu for muon sample; el for el sample

    #if options.fitwtaggersim:
    #    print 'fitwtagger for el+mu sample'
    #    control_sample_simultaneous()#mu for muon sample; el for el sample

    if options.control:
        print 'control for %s sample'%(options.opt_channel)
        control_single(options.opt_channel)

    if options.fitsignal:
        print "fitsignal"
        fit_signal("method1", options.opt_channel, "BulkGravWW600", 200, 1500)
        fit_signal("method1", options.opt_channel, "BulkGravWW700", 200, 1500)
        fit_signal("method1", options.opt_channel, "BulkGravWW750", 200, 1500)
        fit_signal("method1", options.opt_channel, "BulkGravWW800", 200, 1500)
        fit_signal("method1", options.opt_channel, "BulkGravWW900", 200, 1500)
        fit_signal("method1", options.opt_channel, "BulkGravWW1000", 200, 1500)
        fit_signal("method1", options.opt_channel, "BulkGravWW1200", 200, 2200)
        fit_signal("method1", options.opt_channel, "BulkGravWW1400", 400, 2400)
        fit_signal("method1", options.opt_channel, "BulkGravWW1800", 800, 2800)
        fit_signal("method1", options.opt_channel, "BulkGravWW2000", 1000, 3000)
        fit_signal("method1", options.opt_channel, "BulkGravWW2500", 1500, 3500)
        fit_signal("method1", options.opt_channel, "BulkGravWW3000", 2000, 4000)
        fit_signal("method1", options.opt_channel, "BulkGravWW3500", 2500, 4500)
        fit_signal("method1", options.opt_channel, "BulkGravWW4000", 3000, 5000)
        fit_signal("method1", options.opt_channel, "BulkGravWW4500", 3500, 5500)

    if options.check:
        print '################# check workspace for %s sample'%(options.opt_channel)
        check_workspace(options.opt_channel, "BulkGravWW750")
        #check_workspace(options.opt_channel, "BulkGravWW4500")

    if options.combine:
        print '################# check workspace for %s sample'%(options.opt_channel)
        combine_el_mu("em")

    if options.simple and (not options.fitwtagger) and (not options.fitwtaggersim) and (not options.multi) and (not options.control) and (not options.check) and (not options.combine) and (not options.fitsignal):
        print '################# simple mode for %s sample'%(options.opt_channel)
        pre_limit_simple(options.opt_channel)

    ### real function called by the command line parsing some arguments as argv
    if options.multi  and (not options.fitsignal):
        print '################# multi mode for %s sample'%(options.opt_channel)
        pre_limit_sb_correction("method1", sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]), sys.argv[9], sys.argv[10])

