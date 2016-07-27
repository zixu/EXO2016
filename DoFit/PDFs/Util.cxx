/*
 * =====================================================================================
 *
 *       Filename:  Util.cxx
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/26/2012 03:27:25 PM CST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   Zijun Xu (xuzijun123@gmail.com), 
 *        Company:  
 *
 * =====================================================================================
 */


#include <algorithm>
#include <vector>
#include <string>
#include "Math/Math.h"
#include "Math/QuantFuncMathCore.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooPlot.h"
#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2.h"
#include "TChain.h"
#include "TString.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TIterator.h"
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooCurve.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooMCStudy.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include <iostream>


#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH2.h"
#include "RooFitResult.h"
#include "TStyle.h"
#include "TDirectory.h"

using namespace std;
using namespace RooFit;

/// function used to draw an error band around a RooAbsPdf -> used to draw the band after each fit
void draw_error_band(RooAbsData &rdata, RooAbsPdf &rpdf, RooRealVar &rrv_number_events , RooFitResult *rfres, RooPlot *mplot, Int_t kcolor=6, std::string opt ="F", Int_t number_point=100, const Int_t number_errorband=2000){

	TRandom3 rand(1234);

	/// get the observables of the pdf --> mj or mlvj depends on which bands you are drawing
	RooArgSet* argset_obs=rpdf.getObservables(rdata);
	/// create an iterator on the observables
	TIterator *par=argset_obs->createIterator();
	par->Reset();
	/// extract the RooRealVar 
	RooRealVar *rrv_x=(RooRealVar*)par->Next();
	rrv_x->Print();
	rpdf.Print("v");
	/// Get the pdf pramters 
	rpdf.getParameters(RooArgSet(*rrv_x))->Print("v");
	rrv_number_events.Print();
	/// Define min and max x range
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	/// define the step to sample the function 
	Double_t delta_x=(x_max-x_min)/number_point;
	/// Original Bin width
	Double_t width_x = mplot->getFitRangeBinW();
	/// number of events  
	Double_t number_events_mean = rrv_number_events.getVal();
	Double_t number_events_sigma= rrv_number_events.getError();
	/// Create a Graph for the central background prediction 
	TGraph *bkgpred = new TGraph(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		rrv_x->setVal(x_min+delta_x*i); 
		bkgpred->SetPoint( i , x_min+delta_x*i , rrv_number_events.getVal()*rpdf.getVal(*rrv_x)*width_x );
	}
	bkgpred->SetLineWidth(2);
	bkgpred->SetLineColor(kcolor);

	/// Set of parameters
	RooArgSet* par_pdf  = rpdf.getParameters(RooArgSet(*rrv_x)) ;

	/// Build a envelope using number_errorband toys
	TGraph* syst[number_errorband];
	for(int j=0;j<number_errorband;j++){
		syst[j] = new TGraph(number_point+1);
		/// paramters value are randomized using rfres and this can be done also if they are not decorrelate
		RooArgList par_tmp = rfres->randomizePars();
		*par_pdf = par_tmp;
		Double_t number_events_tmp = rand.Gaus(number_events_mean,number_events_sigma); /// new poisson random number of events
		for(int i =0 ; i<=number_point ; i++){
			rrv_x->setVal(x_min+delta_x*i); 
			syst[j]->SetPoint( i , x_min+delta_x*i , number_events_tmp*rpdf.getVal(*rrv_x)*width_x);
		}
	}

	RooArgList par_tmp = rfres->floatParsFinal();
	*par_pdf = par_tmp;

	std::vector<double> val;
	val.resize(number_errorband);
	// Try to build and find max and minimum for each point --> not the curve but the value to do a real envelope -> take one 2sigma interval
	TGraph *ap=new TGraph(number_point+1);
	TGraph *am=new TGraph(number_point+1);
	TGraphAsymmErrors* errorband=new TGraphAsymmErrors(number_point+1);
	ap->SetName("error_up");
	am->SetName("error_dn");
	errorband->SetName("errorband");
	for(int i =0 ; i<= number_point ; i++){
		for(int j=0;j<number_errorband;j++){
			val[j]=(syst[j])->GetY()[i];
		}
		std::sort(val.begin(),val.end());
		ap->SetPoint(i, x_min+delta_x*i,val[Int_t(0.16*number_errorband)]);
		am->SetPoint(i, x_min+delta_x*i,val[Int_t(0.84*number_errorband)]);
		errorband->SetPoint(i, x_min+delta_x*i,bkgpred->GetY()[i] );
		errorband->SetPointError(i, 0.,0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)],val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]);
	}
	ap->SetLineWidth(2);
	ap->SetLineColor(kcolor);
	am->SetLineWidth(2);
	am->SetLineColor(kcolor);
	errorband->SetFillColor(kBlack);
	errorband->SetFillStyle(3013);

	if( TString(opt).Contains("F") ) mplot->addObject(errorband,"E3");             
	if( TString(opt).Contains("L") ){ am->SetMarkerStyle(1);  ap->SetMarkerStyle(1); mplot->addObject(am); mplot->addObject(ap); }
}

/// Variation of the previous method using a workspace -> used to draw the band for the final extrapolation -> take in input decorrelated parameters
void draw_error_band( RooAbsPdf &rpdf, std::string xaxis_name, RooRealVar &rrv_number_events , RooArgList &paras, RooWorkspace &ws, RooPlot *mplot, Int_t kcolor=6, std::string opt="F", Int_t number_point=100, const Int_t number_errorband=2000){


	std::cout << "I AM HERE : draw_error_band " <<  std::endl;
	TRandom3 rand(1234);
	/// get observables , pdf and number of events
	RooRealVar *rrv_x=ws.var(xaxis_name.c_str());
	rpdf.Print("v");
	rpdf.getParameters(RooArgSet(*rrv_x))->Print("v");
	rrv_number_events.Print();

	/// Define the sampling of the input pdf
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	Double_t delta_x=(x_max-x_min)/number_point;
	Double_t width_x=mplot->getFitRangeBinW();

	Double_t number_events_mean = rrv_number_events.getVal();
	Double_t number_events_sigma= rrv_number_events.getError();

	/// TGraph for the central bkg prediction 
	TGraph *bkgpred=new TGraph(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		rrv_x->setVal(x_min+delta_x*i); 
		bkgpred->SetPoint( i , x_min+delta_x*i , rrv_number_events.getVal()*rpdf.getVal(*rrv_x)*width_x );
	}
	bkgpred->SetLineWidth(2);
	bkgpred->SetLineColor(kcolor);

	/// Define the curve in each toy and fill them  not using the randomized par but vaying them by hand -> to be decorrelated
	TGraph* syst[number_errorband];
	for(int j=0;j<number_errorband;j++){
		for(Int_t ipara=0;ipara<paras.getSize();ipara++){
			ws.var(paras[ipara].GetName())->setConstant(0);
			ws.var(paras[ipara].GetName())->setVal( rand.Gaus(0.,ws.var(paras[ipara].GetName())->getError()) );
		}

		Double_t number_events_tmp = rand.Gaus(number_events_mean,number_events_sigma);
		syst[j]=new TGraph(number_point+1);
		for(int i =0 ; i<=number_point ; i++){
			rrv_x->setVal(x_min+delta_x*i); 
			syst[j]->SetPoint( i , x_min+delta_x*i , number_events_tmp*rpdf.getVal(*rrv_x)*width_x);
		}
	}

	/// Now look for the envelop at 2sigma CL
	std::vector<double> val;
	val.resize(number_errorband);
	TGraph *ap = new TGraph(number_point+1);
	TGraph *am = new TGraph(number_point+1);
	TGraphAsymmErrors* errorband=new TGraphAsymmErrors(number_point+1);
	ap->SetName("error_up");
	am->SetName("error_dn");
	errorband->SetName("errorband");
	for(int i =0 ; i<= number_point ; i++){
		for(int j=0;j<number_errorband;j++){
			val[j]=(syst[j])->GetY()[i];
		}
		std::sort(val.begin(),val.end());
		ap->SetPoint(i, x_min+delta_x*i,val[Int_t(0.16*number_errorband)]);
		am->SetPoint(i, x_min+delta_x*i,val[Int_t(0.84*number_errorband)]);
		errorband->SetPoint(i, x_min+delta_x*i,bkgpred->GetY()[i] );
		errorband->SetPointError(i, 0.,0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)],val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]);
	}
	ap->SetLineWidth(2);
	ap->SetLineColor(kcolor);
	am->SetLineWidth(2);
	am->SetLineColor(kcolor);
	errorband->SetFillColor(kBlack);
	errorband->SetFillStyle(3013);

	if( TString(opt).Contains("F") ) mplot->addObject(errorband,"E3");
	if( TString(opt).Contains("L") ){ am->SetMarkerStyle(1); ap->SetMarkerStyle(1); mplot->addObject(am); mplot->addObject(ap); }
}


//void draw_error_band( RooAbsPdf &rpdf, std::string xaxis_name, RooRealVar &rrv_number_events , RooArgList &paras, RooWorkspace &ws, RooPlot *mplot, RooPlot *mplotP, RooDataHist *datahist, Int_t kcolor=6, std::string opt="F", Int_t number_point=100, const Int_t number_errorband=2000){
void draw_error_band( RooAbsPdf &rpdf, std::string xaxis_name, RooRealVar &rrv_number_events , RooArgList &paras, RooWorkspace &ws, RooPlot *mplot, RooPlot *mplotP, TH1 *hdata, Int_t kcolor=6, std::string opt="F", Int_t number_point=100, const Int_t number_errorband=2000){
	//void draw_error_band( RooAbsPdf &rpdf, std::string xaxis_name, RooRealVar &rrv_number_events , RooArgList &paras, RooWorkspace &ws, RooPlot *mplot, RooPlot *mplotP, Int_t kcolor=6, std::string opt="F", Int_t number_point=100, const Int_t number_errorband=2000){

	//}
	//else{ 
	//std::cout << "I AM HERE " << paras[ipara].GetName() << std::endl;
	//ws.var(paras[ipara].GetName())->setConstant(1);
	//}

	std::cout << "I AM HERE : draw_error_band " <<  std::endl;
	//Int_t hpoints = datahist->numEntries();
	//std::cout << hpoints<<endl;
	//TH1 *hdata=datahist->createHistogram("rrv_mass_lvj",hpoints);
	Int_t hpoints = hdata->GetNbinsX();


	TRandom3 rand(1234);
	/// get observables , pdf and number of events
	RooRealVar *rrv_x=ws.var(xaxis_name.c_str());
	rpdf.Print("v");
	rpdf.getParameters(RooArgSet(*rrv_x))->Print("v");
	rrv_number_events.Print();

	/// Define the sampling of the input pdf
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	Double_t delta_x=(x_max-x_min)/number_point;
	Double_t width_x=mplot->getFitRangeBinW();

	Double_t number_events_mean = rrv_number_events.getVal();
	Double_t number_events_sigma= rrv_number_events.getError();

	/// TGraph for the central bkg prediction 
	TGraph *bkgpred=new TGraph(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		rrv_x->setVal(x_min+delta_x*i); 
		bkgpred->SetPoint( i , x_min+delta_x*i , rrv_number_events.getVal()*rpdf.getVal(*rrv_x)*width_x );
	}
	bkgpred->SetLineWidth(2);
	bkgpred->SetLineColor(kcolor);

	/// Define the curve in each toy and fill them  not using the randomized par but vaying them by hand -> to be decorrelated
	TGraph* syst[number_errorband];
	for(int j=0;j<number_errorband;j++){
		for(Int_t ipara=0;ipara<paras.getSize();ipara++){
			ws.var(paras[ipara].GetName())->setConstant(0);
			ws.var(paras[ipara].GetName())->setVal( rand.Gaus(0.,ws.var(paras[ipara].GetName())->getError()) );
		}

		Double_t number_events_tmp = rand.Gaus(number_events_mean,number_events_sigma);
		syst[j]=new TGraph(number_point+1);
		for(int i =0 ; i<=number_point ; i++){
			rrv_x->setVal(x_min+delta_x*i); 
			syst[j]->SetPoint( i , x_min+delta_x*i , number_events_tmp*rpdf.getVal(*rrv_x)*width_x);
		}
	}

	/// Now look for the envelop at 2sigma CL
	std::vector<double> val;
	val.resize(number_errorband);
	TGraph *ap = new TGraph(number_point+1);
	TGraph *am = new TGraph(number_point+1);
	TGraphAsymmErrors* errorband=new TGraphAsymmErrors(number_point+1);
	TGraphAsymmErrors* errorbandP=new TGraphAsymmErrors(number_point-1);
	ap->SetName("error_up");
	am->SetName("error_dn");
	errorband->SetName("errorband");
	Double_t x_tmp = 0.;
	for(int i =0 ; i<= number_point ; i++){
		for(int j=0;j<number_errorband;j++){
			val[j]=(syst[j])->GetY()[i];
		}
		std::sort(val.begin(),val.end());
		ap->SetPoint(i, x_min+delta_x*i,val[Int_t(0.16*number_errorband)]);
		am->SetPoint(i, x_min+delta_x*i,val[Int_t(0.84*number_errorband)]);
		errorband->SetPoint(i, x_min+delta_x*i,bkgpred->GetY()[i] );
		x_tmp = x_min+delta_x*i;
		errorband->SetPointError(i, 0.,0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)],val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]);
		double dataEYlow, dataEYhigh , hdata_err, alpha, N, L, U;
		for (int k=1;k<=hpoints;k++) {
			if(hdata->GetBinLowEdge(k)<x_tmp && x_tmp <  hdata->GetBinLowEdge(k)+hdata->GetBinWidth(k) ) {

				alpha = 1 - 0.6827;
				N = hdata->GetBinContent(k);
				if (N==0) {L = 0;}
				else { L = ROOT::Math::gamma_quantile(alpha/2,N,1.);}

				U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1);
				cout<<"N="<<N<<", L="<<L<<", U="<<U<<endl;
				dataEYlow = N-L;
				dataEYhigh= U-N;
				if ( hdata->GetBinContent(k) - bkgpred->GetY()[i] >0 ) {
					hdata_err = dataEYlow;
				} else {
					hdata_err = dataEYhigh;
				}
				cout<<"i="<<i<<" x_tmp="<<x_tmp<<" dataEYlow="<<dataEYlow<<" dataEYhigh="<<dataEYhigh<<endl;
				//cout<<"hdata_err="<<hdata_err<<endl;
				//std::cout<< "Hello World "<<i<<"   "<< x_tmp<<"\t"<<hdata->GetBinError(k)<<"  "<<hdata_err<<std::endl;
				errorbandP->SetPoint(i, x_tmp,0.); 
				if(hdata_err>0 ) {
					errorbandP->SetPointError(i, 0. , 0., (bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)])/ hdata_err ,(val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i])/ hdata_err); 
					//cout<<i<<" x_tmp= "<<x_tmp<<" error= "<<(bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)])/ hdata_err<<", "<<(val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i])/ hdata_err<<endl; 
					cout<<i<<" x_tmp= "<<x_tmp<<" error= "<<(bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)])<<", "<<(val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i])<<" hdata= "<< hdata_err<<endl; 
					//cout<<i<<" "<<bkgpred->GetY()[i]<<", "<<val[Int_t(0.84*number_errorband)]<<", "<<val[Int_t(0.16*number_errorband)]<<", "<<hdata_err<<endl; 
				}else{
					errorbandP->SetPointError(i, 0. , 0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)] ,val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]); 
				}                      
			}
		} 
	}
	ap->SetLineWidth(2);
	ap->SetLineColor(kcolor);
	am->SetLineWidth(2);
	am->SetLineColor(kcolor);
	errorband->SetFillColor(kBlack);
	errorband->SetFillStyle(3013);
	errorbandP->SetFillColor(kOrange-2);
	//        errorbandP->SetFillStyle(3013);


	if( TString(opt).Contains("F") ) mplot->addObject(errorband,"E3");
	if( TString(opt).Contains("F") ) mplotP->addObject(errorbandP,"E3");
	if( TString(opt).Contains("L") ){ am->SetMarkerStyle(1); ap->SetMarkerStyle(1);  mplot->addObject(am); mplot->addObject(ap); }
}

//void draw_error_band( RooAbsPdf &rpdf, std::string xaxis_name, RooRealVar &rrv_number_events , RooArgList &paras, RooWorkspace &ws, RooPlot *mplot, RooPlot *mplotP, RooDataHist *datahist, Int_t kcolor=6, std::string opt="F", Int_t number_point=100, const Int_t number_errorband=2000){
//
//	Int_t hpoints = datahist->numEntries();
//	//std::cout << hpoints<<endl;
//	TH1 *hdata=datahist->createHistogram("rrv_mass_lvj",hpoints);
//	
//	draw_error_band( rpdf, xaxis_name, rrv_number_events, paras, ws, mplot, mplotP, hdata, kcolor, opt, number_point, number_errorband);
//}

/// Draw error band giving directly the extended Pdf
void draw_error_band_extendPdf(RooAbsData &rdata, RooExtendPdf &rpdf, RooFitResult *rfres, RooPlot *mplot, Int_t kcolor=6, std::string opt="F", Int_t number_point=100, const Int_t number_errorband=2000){

	TRandom3 rand(1234);
	/// Take the observable for the pdf
	RooArgSet* argset_obs=rpdf.getObservables(rdata);
	TIterator *par=argset_obs->createIterator();
	par->Reset();
	RooRealVar *rrv_x=(RooRealVar*)par->Next();
	rrv_x->Print();
	/// Define the sampling
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	Double_t delta_x=(x_max-x_min)/number_point;
	Double_t width_x=mplot->getFitRangeBinW();

	/// Central value for the bkg prediction 
	TGraph *bkgpred=new TGraph(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		rrv_x->setVal(x_min+delta_x*i); 
		bkgpred->SetPoint( i , x_min+delta_x*i , rpdf.expectedEvents(*rrv_x)*rpdf.getVal(*rrv_x)*width_x );
	}
	bkgpred->SetLineWidth(2);
	bkgpred->SetLineColor(kcolor);

	/// Take the parameters
	RooArgSet* par_pdf  = rpdf.getParameters(RooArgSet(*rrv_x)) ;
	par_pdf->Print("v");

	/// Make the envelope
	TGraph* syst[number_errorband];
	for(int j=0;j<number_errorband;j++){
		syst[j]=new TGraph(number_point+1);
		RooArgList par_tmp = rfres->randomizePars();
		*par_pdf = par_tmp;
		for(int i =0 ; i<=number_point ; i++){
			rrv_x->setVal(x_min+delta_x*i); 
			syst[j]->SetPoint( i , x_min+delta_x*i , rpdf.expectedEvents(*rrv_x)*rpdf.getVal(*rrv_x)*width_x);
		}
	}

	RooArgList par_tmp = rfres->floatParsFinal();
	*par_pdf = par_tmp;

	/// now extract the error curve at 2sigma 
	std::vector<double> val;
	val.resize(number_errorband);
	TGraph *ap=new TGraph(number_point+1);
	TGraph *am=new TGraph(number_point+1);
	TGraphAsymmErrors* errorband=new TGraphAsymmErrors(number_point+1);
	ap->SetName("error_up");
	am->SetName("error_dn");
	errorband->SetName("errorband");
	for(int i =0 ; i<= number_point ; i++){
		for(int j=0;j<number_errorband;j++){
			val[j]=(syst[j])->GetY()[i];
		}
		std::sort(val.begin(),val.end());
		ap->SetPoint(i, x_min+delta_x*i,val[Int_t(0.16*number_errorband)]);
		am->SetPoint(i, x_min+delta_x*i,val[Int_t(0.84*number_errorband)]);
		errorband->SetPoint(i, x_min+delta_x*i,bkgpred->GetY()[i] );
		errorband->SetPointError(i, 0.,0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)],val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]);
	}
	ap->SetLineWidth(2);
	ap->SetLineColor(kcolor);
	am->SetLineWidth(2);
	am->SetLineColor(kcolor);
	errorband->SetFillColor(kBlack);
	errorband->SetFillStyle(3013);

	if( TString(opt).Contains("F") ) mplot->addObject(errorband,"E3"); 
	if( TString(opt).Contains("L") ){ am->SetMarkerStyle(1);ap->SetMarkerStyle(1);  mplot->addObject(am); mplot->addObject(ap); }
}

/// Draw error band giving directly the a general Pdf
void draw_error_band_extendPdf(RooAbsData &rdata, RooAbsPdf &rpdf, RooFitResult *rfres, RooPlot *mplot, Int_t kcolor=6, std::string opt="F", Int_t number_point=100, const Int_t number_errorband=2000){

	TRandom3 rand(1234);
	/// Take the observable for the pdf
	RooArgSet* argset_obs=rpdf.getObservables(rdata);
	TIterator *par=argset_obs->createIterator();
	par->Reset();
	RooRealVar *rrv_x=(RooRealVar*)par->Next();
	rrv_x->Print();
	/// Define the sampling
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	Double_t delta_x=(x_max-x_min)/number_point;
	Double_t width_x=mplot->getFitRangeBinW();

	/// Central value for the bkg prediction 
	TGraph *bkgpred=new TGraph(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		rrv_x->setVal(x_min+delta_x*i); 
		bkgpred->SetPoint( i , x_min+delta_x*i , rpdf.expectedEvents(*rrv_x)*rpdf.getVal(*rrv_x)*width_x );
	}
	bkgpred->SetLineWidth(2);
	bkgpred->SetLineColor(kcolor);

	/// Take the parameters
	RooArgSet* par_pdf  = rpdf.getParameters(RooArgSet(*rrv_x)) ;
	par_pdf->Print("v");

	/// Make the envelope
	TGraph* syst[number_errorband];
	for(int j=0;j<number_errorband;j++){
		syst[j]=new TGraph(number_point+1);
		RooArgList par_tmp = rfres->randomizePars();
		*par_pdf = par_tmp;
		for(int i =0 ; i<=number_point ; i++){
			rrv_x->setVal(x_min+delta_x*i); 
			syst[j]->SetPoint( i , x_min+delta_x*i , rpdf.expectedEvents(*rrv_x)*rpdf.getVal(*rrv_x)*width_x);
		}
	}

	RooArgList par_tmp = rfres->floatParsFinal();
	*par_pdf = par_tmp;

	/// now extract the error curve at 2sigma 
	std::vector<double> val;
	val.resize(number_errorband);
	TGraph *ap=new TGraph(number_point+1);
	TGraph *am=new TGraph(number_point+1);
	TGraphAsymmErrors* errorband=new TGraphAsymmErrors(number_point+1);
	ap->SetName("error_up");
	am->SetName("error_dn");
	errorband->SetName("errorband");
	for(int i =0 ; i<= number_point ; i++){
		for(int j=0;j<number_errorband;j++){
			val[j]=(syst[j])->GetY()[i];
		}
		std::sort(val.begin(),val.end());
		ap->SetPoint(i, x_min+delta_x*i,val[Int_t(0.16*number_errorband)]);
		am->SetPoint(i, x_min+delta_x*i,val[Int_t(0.84*number_errorband)]);
		errorband->SetPoint(i, x_min+delta_x*i,bkgpred->GetY()[i] );
		errorband->SetPointError(i, 0.,0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)],val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]);
	}
	ap->SetLineWidth(2);
	ap->SetLineColor(kcolor);
	am->SetLineWidth(2);
	am->SetLineColor(kcolor);
	errorband->SetFillColor(kBlack);
	errorband->SetFillStyle(3013);

	if( TString(opt).Contains("F") ) mplot->addObject(errorband,"E3"); 
	if( TString(opt).Contains("L") ){am->SetMarkerStyle(1); ap->SetMarkerStyle(1);  mplot->addObject(am); mplot->addObject(ap); }
}


/// Error band creator for a Decorellated Pdf starting from Workspace
void draw_error_band_Decor( std::string pdf_name, std::string xaxis_name, RooArgList &paras, RooWorkspace &ws,RooRealVar &rrv_shape_scale , RooPlot *mplot, Int_t kcolor=6, std::string opt="F", Int_t number_point=100, const Int_t number_errorband=2000){

	TRandom3 rand(1234);

	RooRealVar *rrv_x=ws.var(xaxis_name.c_str());
	rrv_x->Print();
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	Double_t delta_x=(x_max-x_min)/number_point;
	Double_t width_x=mplot->getFitRangeBinW();

	/// factor to scale the normalized shape
	Double_t shape_scale = rrv_shape_scale.getVal();
	Double_t shape_scale_error = rrv_shape_scale.getError();

	/// Central bkg prediction 
	TGraph *bkgpred=new TGraph(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		rrv_x->setVal(x_min+delta_x*i); 
		bkgpred->SetPoint( i , x_min+delta_x*i , shape_scale*ws.pdf(pdf_name.c_str())->getVal() );
	}
	bkgpred->SetLineWidth(2);
	bkgpred->SetLineColor(kcolor+3);

	/// error band -> each parameter can be randomly gaus generated
	TGraph* syst[number_errorband];
	for(int j=0;j<number_errorband;j++){
		for(Int_t ipara=0;ipara<paras.getSize();ipara++){
			ws.var(paras[ipara].GetName())->setVal( rand.Gaus(0.,1.) );
		}

		/// Change the scaling value
		Double_t shape_scale_tmp=rand.Gaus(shape_scale,shape_scale_error);
		syst[j]=new TGraph(number_point+1);
		for(int i =0 ; i<=number_point ; i++){
			rrv_x->setVal(x_min+delta_x*i); 
			syst[j]->SetPoint( i , x_min+delta_x*i ,shape_scale_tmp*ws.pdf(pdf_name.c_str())->getVal()*width_x);
		}
	}

	std::vector<double> val;
	val.resize(number_errorband);
	TGraph *ap=new TGraph(number_point+1);
	TGraph *am=new TGraph(number_point+1);
	TGraphAsymmErrors* errorband=new TGraphAsymmErrors(number_point+1);
	ap->SetName("error_up");
	am->SetName("error_dn");
	errorband->SetName("errorband");
	for(int i =0 ; i<= number_point ; i++){
		for(int j=0;j<number_errorband;j++){
			val[j]=(syst[j])->GetY()[i];
		}
		std::sort(val.begin(),val.end());
		ap->SetPoint(i, x_min+delta_x*i,val[Int_t(0.16*number_errorband)]);
		am->SetPoint(i, x_min+delta_x*i,val[Int_t(0.84*number_errorband)]);
		errorband->SetPoint(i, x_min+delta_x*i,bkgpred->GetY()[i] );
		errorband->SetPointError(i, 0.,0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)],val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]);
	}
	ap->SetLineWidth(2);
	ap->SetLineColor(kcolor);
	am->SetLineWidth(2);
	am->SetLineColor(kcolor);
	errorband->SetFillColor(kBlack);
	errorband->SetFillStyle(3013);
	errorband->SetName("Uncertainty");


	if( TString(opt).Contains("F") ){
		mplot->addObject(errorband,"E3"); 
		mplot->addObject(bkgpred);
	}
	if( TString(opt).Contains("L") ){
		ap->SetMarkerStyle(1); 
		am->SetMarkerStyle(1); 
		mplot->addObject(am); mplot->addObject(ap); 
		mplot->addObject(bkgpred);
	}

	for(Int_t ipara=0;ipara<paras.getSize();ipara++){
		ws.var(paras[ipara].GetName())->setVal(0.);
	}

} 

/// Just the shape and don't touch the normalization 
void draw_error_band_shape_Decor( std::string pdf_name, std::string xaxis_name, RooArgList &paras, RooWorkspace &ws,Double_t sigma , RooPlot *mplot, Int_t kcolor=6,std::string opt="F", Int_t fillstyle=3013,std::string uncertainty_title="", Int_t number_point=100, const Int_t number_errorband=2000){

	TRandom3 rand(1234);
	/// take the observable
	RooRealVar *rrv_x=ws.var(xaxis_name.c_str());
	rrv_x->Print();
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	Double_t delta_x=(x_max-x_min)/number_point;
	Double_t width_x=mplot->getFitRangeBinW();

	/// bkg prediction central value
	TGraph *bkgpred=new TGraph(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		rrv_x->setVal(x_min+delta_x*i); 
		bkgpred->SetPoint( i , x_min+delta_x*i , ws.pdf(pdf_name.c_str())->getVal(*rrv_x)*width_x );
	}
	bkgpred->SetLineWidth(2);
	bkgpred->SetLineColor(kcolor+3);

	// make the envelope
	TGraph* syst[number_errorband];
	for(int j=0;j<number_errorband;j++){
		for(Int_t ipara=0;ipara<paras.getSize();ipara++){
			ws.var(paras[ipara].GetName())->setVal( rand.Gaus(0.,sigma) ); // choose how many sigma on the parameters you wamt
		}
		syst[j]=new TGraph(number_point+1);
		for(int i =0 ; i<=number_point ; i++){
			rrv_x->setVal(x_min+delta_x*i); 
			syst[j]->SetPoint( i , x_min+delta_x*i ,ws.pdf(pdf_name.c_str())->getVal(*rrv_x)*width_x);
		}
	}

	std::vector<double> val;
	val.resize(number_errorband);
	TGraph *ap=new TGraph(number_point+1);
	TGraph *am=new TGraph(number_point+1);
	TGraphAsymmErrors* errorband=new TGraphAsymmErrors(number_point+1);
	ap->SetName("error_up");
	am->SetName("error_dn");
	errorband->SetName("errorband");
	for(int i =0 ; i<= number_point ; i++){
		for(int j=0;j<number_errorband;j++){
			val[j]=(syst[j])->GetY()[i];
		}
		std::sort(val.begin(),val.end());
		ap->SetPoint(i, x_min+delta_x*i,val[Int_t(0.16*number_errorband)]);
		am->SetPoint(i, x_min+delta_x*i,val[Int_t(0.84*number_errorband)]);
		errorband->SetPoint(i, x_min+delta_x*i,bkgpred->GetY()[i] );
		errorband->SetPointError(i, 0.,0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)],val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]);
		//if (i==55)cout<<val[Int_t(0.16*number_errorband)]<<"  "<<bkgpred->GetY()[i] <<"  "<<val[Int_t(0.84*number_errorband)] <<endl;
	}
	ap->SetLineWidth(2);
	ap->SetLineColor(kcolor);
	am->SetLineWidth(2);
	am->SetLineColor(kcolor);
	errorband->SetFillStyle(fillstyle);
	errorband->SetFillColor(kcolor);
	errorband->SetName(Form("%s %g#sigma",uncertainty_title.c_str(),sigma));

	if( TString(opt).Contains("F") ){ mplot->addObject(errorband,"E3"); }
	if( TString(opt).Contains("L") ){  am->SetMarkerStyle(1); ap->SetMarkerStyle(1);mplot->addObject(am); mplot->addObject(ap); }
	for(Int_t ipara=0;ipara<paras.getSize();ipara++){
		ws.var(paras[ipara].GetName())->setVal(0.);
	}
} 


/// Calculate the error when intgrating a pdf in a range -> take the error not as a single fit result but from toys randomizing the parameters
double Calc_error_extendPdf(RooAbsData &rdata, RooExtendPdf &rpdf, RooFitResult *rfres, std::string range, const Int_t calc_times=2000){

	TRandom3 rand(1234);
	/// Get the observable on the x-axis
	RooArgSet* argset_obs=rpdf.getObservables(rdata);
	TIterator *par=argset_obs->createIterator();
	par->Reset();
	RooRealVar *rrv_x=(RooRealVar*)par->Next();
	rrv_x->Print();

	/// Create integral of the whole function and in a given range
	RooAbsReal* fullInt = rpdf.createIntegral(*rrv_x,*rrv_x);
	RooAbsReal* signalInt = rpdf.createIntegral(*rrv_x,*rrv_x,range.c_str());
	double fullInt_var=fullInt->getVal();
	double signalInt_var=signalInt->getVal()/fullInt_var;

	double signal_number_media = signalInt_var*rpdf.expectedEvents(*rrv_x);

	/// Take the parameters and randomize them within uncertainty and do a lot of toys
	RooArgSet* par_pdf  = rpdf.getParameters(RooArgSet(*rrv_x)) ;
	par_pdf->Print("v");

	std::vector<double> val;
	val.resize(calc_times);
	for(int j=0;j<calc_times;j++){
		RooArgList par_tmp = rfres->randomizePars();
		*par_pdf = par_tmp;

		fullInt_var=fullInt->getVal();
		signalInt_var=signalInt->getVal()/fullInt_var;
		signal_number_media=signalInt_var*rpdf.expectedEvents(*rrv_x);
		val[j]=signal_number_media;
	}

	std::sort(val.begin(),val.end());
	return (val[Int_t(0.84*calc_times)]-val[Int_t(0.16*calc_times)])/2.; /// return a doble value 
}

/// useful for couting analysis
double Calc_error( std::string rpdfname, std::string xaxis_name , RooArgList &paras, RooWorkspace &ws,std::string range, const Int_t calc_times=2000){

	TRandom3 rand(1234);
	// Get the observable
	RooRealVar *rrv_x=ws.var(xaxis_name.c_str()); 
	rrv_x->Print();
	RooAbsPdf*rpdf=ws.pdf(rpdfname.c_str());

	RooAbsReal* fullInt = rpdf->createIntegral(*rrv_x,*rrv_x);
	RooAbsReal* signalInt = rpdf->createIntegral(*rrv_x,*rrv_x,range.c_str());
	double fullInt_var=fullInt->getVal();
	double signalInt_var=signalInt->getVal()/fullInt_var;

	double signal_number_media=signalInt_var;

	std::vector<double> val;
	val.resize(calc_times);
	for(int j=0;j<calc_times;j++){
		for(Int_t ipara=0;ipara<paras.getSize();ipara++){
			ws.var(paras[ipara].GetName())->setConstant(0);
			ws.var(paras[ipara].GetName())->setVal( rand.Gaus(0.,ws.var(paras[ipara].GetName())->getError()) );
		}

		fullInt_var=fullInt->getVal();
		signalInt_var=signalInt->getVal()/fullInt_var;

		signal_number_media=signalInt_var;
		val[j]=signal_number_media;
	}
	for(Int_t ipara=0;ipara<paras.getSize();ipara++){ ws.var(paras[ipara].GetName())->setVal(0.); }

	std::sort(val.begin(),val.end());
	double number_error=(val[Int_t(0.84*calc_times)]-val[Int_t(0.16*calc_times)])/2./signal_number_media;
	return number_error;
}


void MCStudy(RooAbsPdf* model, RooRealVar x, TString picname, Int_t nTop=400, Int_t nEvtPerSample=0){

	RooArgSet* parameters_General= model->getParameters( RooArgSet(x));
	parameters_General->Print("v");

	RooRealVar *para0;
	RooRealVar *para1;
	TIterator* par=parameters_General->createIterator();
	par->Reset();
	RooRealVar* param=(RooRealVar*)par->Next();
	for (Int_t i=0; param; i++ ){
		param->Print();
		if(i==0) para0=param;
		if(i==1) para1=param;
		param=(RooRealVar*)par->Next();
	}
	para0->Print("v");
	para1->Print("v");
	
	RooMCStudy* mcstudy = new RooMCStudy(*model,x,Binned(kTRUE),Silence(),Extended(), FitOptions(Save(kTRUE),PrintEvalErrors(0))) ;
	// G e n e r a t e   a n d   f i t   e v e n t s
	// ---------------------------------------------

	// Generate and fit 1000 samples of Poisson(nExpected) events
	mcstudy->generateAndFit(nTop,nEvtPerSample);


	// E x p l o r e   r e s u l t s   o f   s t u d y 
	// ------------------------------------------------

	// Make plots of the distributions of para0, the error on para0 and the pull of para0
	RooPlot* frame1 = mcstudy->plotParam(*para0,Bins(40)) ;
	RooPlot* frame2 = mcstudy->plotError(*para0,Bins(40)) ;
	RooPlot* frame3 = mcstudy->plotPull(*para0,Bins(40),FitGauss(kTRUE)) ;

	RooPlot* frame4 = mcstudy->plotParam(*para1,Bins(40)) ;
	RooPlot* frame5 = mcstudy->plotError(*para1,Bins(40)) ;
	RooPlot* frame6 = mcstudy->plotPull(*para1,Bins(40),FitGauss(kTRUE)) ;


	// Plot distribution of minimized likelihood
	RooPlot* frame7 = mcstudy->plotNLL(Bins(40)) ;

	// Access some of the saved fit results from individual toys
	TH2* corrHist000 = mcstudy->fitResult(0)->correlationHist("c000") ;
	TH2* corrHist127 = mcstudy->fitResult(127)->correlationHist("c127") ;

	// Draw all plots on a canvas
	TStyle tmpstyle;
	tmpstyle.SetPalette(1) ;
	tmpstyle.SetOptStat(0) ;
	TCanvas* c = new TCanvas("rf801_mcstudy","rf801_mcstudy",900,900) ;
	c->Divide(3,3) ;
	c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
	c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
	c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;
	c->cd(4) ; gPad->SetLeftMargin(0.15) ; frame4->GetYaxis()->SetTitleOffset(1.4) ; frame4->Draw() ;
	c->cd(5) ; gPad->SetLeftMargin(0.15) ; frame5->GetYaxis()->SetTitleOffset(1.4) ; frame5->Draw() ;
	c->cd(6) ; gPad->SetLeftMargin(0.15) ; frame6->GetYaxis()->SetTitleOffset(1.4) ; frame6->Draw() ;
	c->cd(7) ; gPad->SetLeftMargin(0.15) ; corrHist000->GetYaxis()->SetTitleOffset(1.4) ; corrHist000->Draw("colz") ;
	c->cd(8) ; gPad->SetLeftMargin(0.15) ; corrHist127->GetYaxis()->SetTitleOffset(1.4) ; corrHist127->Draw("colz") ;

	c->SaveAs(picname);
}
