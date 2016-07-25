#! /usr/bin/env python
import os
import ROOT

from optparse import OptionParser

parser = OptionParser()
parser.add_option('--vclean', help='clean all the so files', type=int, default=0.)

(options, args) = parser.parse_args()


if __name__ == "__main__":

  ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include");

  inPath = os.getenv("PWD")

  ##os.chdir(inPath+"/PlotStyle");

  ##if options.vclean : os.system("rm Util_cxx.so ; rm PlotUtils_cxx.so");
  ##
  ##ROOT.gROOT.ProcessLine(".L Util.cxx+");
  ##ROOT.gSystem.Load("Util_cxx.so");

  ##ROOT.gROOT.ProcessLine(".L PlotUtils.cxx+");
  ##ROOT.gSystem.Load("PlotUtils_cxx.so");
 
  os.chdir(inPath+"/PDFs");

  if options.vclean : os.system("rm PdfDiagonalizer_cc.so ; rm HWWLVJRooPdfs_cxx.so ; rm MakePdf_cxx.so; rm Util_cxx.so ");

  ROOT.gROOT.ProcessLine(".L PdfDiagonalizer.cc+");
  ROOT.gSystem.Load("PdfDiagonalizer_cc.so");

  ROOT.gROOT.ProcessLine(".L HWWLVJRooPdfs.cxx+");
  ROOT.gSystem.Load("HWWLVJRooPdfs_cxx.so");

  ROOT.gROOT.ProcessLine(".L Util.cxx+");
  ROOT.gSystem.Load("Util_cxx.so");

  ##os.chdir(inPath+"/BiasStudy");
  ##
  ##if options.vclean : os.system("rm BiasUtils_cxx.so");
  ##ROOT.gROOT.ProcessLine(".L BiasUtils.cxx+");
  ##ROOT.gSystem.Load("BiasUtils_cxx.so");

  ##os.chdir(inPath+"/FitUtils");

  ##if options.vclean : os.system("rm FitUtils_cxx.so");
  ##ROOT.gROOT.ProcessLine(".L FitUtils.cxx+");
  ##ROOT.gSystem.Load("FitUtils_cxx.so");
  
