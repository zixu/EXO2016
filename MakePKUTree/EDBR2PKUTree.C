#define EDBR2PKUTree_cxx
#include "EDBR2PKUTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//*****PU WEIGHT***************
vector<Double_t> generate_weights(TH1* data_npu_estimated, Int_t isForSynch){
	// see SimGeneral/MixingModule/python/mix_2015_25ns_Startup_PoissonOOTPU_cfi.pyy; copy and paste from there:
	const Double_t npu_probs[52] = {
		4.8551E-07,
		1.74806E-06,
		3.30868E-06,
		1.62972E-05,
		4.95667E-05,
		0.000606966,
		0.003307249,
		0.010340741,
		0.022852296,
		0.041948781,
		0.058609363,
		0.067475755,
		0.072817826,
		0.075931405,
		0.076782504,
		0.076202319,
		0.074502547,
		0.072355135,
		0.069642102,
		0.064920999,
		0.05725576,
		0.047289348,
		0.036528446,
		0.026376131,
		0.017806872,
		0.011249422,
		0.006643385,
		0.003662904,
		0.001899681,
		0.00095614,
		0.00050028,
		0.000297353,
		0.000208717,
		0.000165856,
		0.000139974,
		0.000120481,
		0.000103826,
		8.88868E-05,
		7.53323E-05,
		6.30863E-05,
		5.21356E-05,
		4.24754E-05,
		3.40876E-05,
		2.69282E-05,
		2.09267E-05,
		1.5989E-05,
		4.8551E-06,
		2.42755E-06,
		4.8551E-07,
		2.42755E-07,
		1.21378E-07,
		4.8551E-08
	};

	if (isForSynch==0) { //OFFICIAL RECIPE
		vector<Double_t> result(52);
		Double_t s = 0.0;
		for(Int_t npu=0; npu<52; ++npu){
			Double_t npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));              
			result[npu] = npu_estimated / npu_probs[npu];
			s += npu_estimated;
		}
		// normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
		for(Int_t npu=0; npu<52; ++npu){
			result[npu] /= s;
		}
		return result;
	}
	else { //THIS IS FOR THE SYNCH ONLY. THIS IS NOT THE OFFICIAL RECIPE!
		vector<Double_t> result(60);
		for(Int_t npu=0; npu<60; ++npu){
			if (data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu))==NULL)
			  result[npu] = 0.;
			else {
				Double_t npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));            
				result[npu] = npu_estimated;
			}
		}
		return result;
	}

}


void EDBR2PKUTree::Loop(TString channelname, Double_t XS, Double_t totaleventnumber, Int_t IsData) {

	std::vector<Double_t> weights_pu1; //these are made with our recipe
	std::vector<Double_t> weights_pu2; //these are made with the official recipe

	TFile* pileupFile1 = TFile::Open("pileupDataRun2015D_72mb.root");  
	TH1F* pileupHisto1 = (TH1F*)pileupFile1->Get("pileup");  
	weights_pu1 = generate_weights(pileupHisto1,0);
	pileupFile1->Close();

	//  TFile* pileupFile2 = TFile::Open("puweights.root");  
	TFile* pileupFile2 = TFile::Open("PUxSynch.root");  
	TH1F *pileupHisto2 = (TH1F*)pileupFile2->Get("puweights");
	weights_pu2 = generate_weights(pileupHisto2,1);
	pileupFile2->Close();

	//TFile * input1 = new TFile ("puweights.root");
	//TH1F* hR1= (TH1F*)input1->Get("puweights");
	//zixu
	TFile * input1 = new TFile ("puweight.root");	
	TH1F* hR1= (TH1F*)input1->Get("h2");
	//TFile * input1 = new TFile ("test_mu.root");
	//TH1F* hR1= (TH1F*)input1->Get("hRatio"); //"pileup");//hRatio");


	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();

	Double_t n_deltaRlepjet = 0; 
	Double_t n_delPhijetlep = 0; 
	Double_t ntau = 0;
	Double_t number_qq = 0; 
	Double_t nmassVhad = 0; 
	Double_t nptVlepJEC = 0;
	Double_t nID_e = 0;
	Double_t npt_e = 0;
	Double_t nmet_e = 0; 
	Double_t nnum_bJet_e = 0; 
	Double_t n_delPhijetmet = 0; 

	Double_t nID_mu = 0;
	Double_t npt_mu = 0;
	Double_t nmet_mu = 0; 
	Double_t nnum_bJet_mu = 0; 
	//Double_t nbtb_mu = 0; 

	Double_t nptVhad = 0;
	Double_t yields = 0;
	//TLorentzVector jetV, genjetV;
	//some constants inside this analysis
	Double_t pi_2=1.57079632679;
	Long64_t npp = fChain->GetEntries("theWeight>0.");
	Long64_t nmm = fChain->GetEntries("theWeight<0.");
	cout<<"npp="<<npp<<" nmm="<<nmm<<" totaleventnumber="<<totaleventnumber<<endl;

	Double_t nn;
	Double_t eff_and_pu_Weight;
	Double_t eff_and_pu_Weight1;
	Float_t Identical_lumiWeight = XS;//All the events inside a sample are same lumiweight
	//Float_t Identical_lumiWeight = XS/totaleventnumber;//All the events inside a sample are same lumiweight

	Long64_t nbytes = 0, nb = 0;
	//for (Long64_t jentry=0; jentry<10;jentry++)
	for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;

		nb = fChain->GetEntry(jentry); 
		nbytes += nb;
		pfMET             = Float_t(met);
		pfMETPhi          = Float_t(metPhi);
		l_pt              = Float_t(ptlep1);
		l_eta             = Float_t(etalep1);
		l_phi             = Float_t(philep1);
		ptVhad            = Float_t(ptVhad);
		jet_eta           = Float_t(yVhad);
		jet_phi           = Float_t(phiVhad);
		jet_mass_pruned   = Float_t(massVhadJEC);
		jetAK8_mass       = Float_t(jetAK8_mass);
		jet_mass_softdrop = Float_t(sdropJEC);
		jet_tau2tau1      = Float_t(tau21);
		W_pt              = Float_t(ptVlepJEC);
		W_eta             = Float_t(yVlep);
		W_phi             = Float_t(phiVlep);
		m_lvj             = Float_t(candMass);
		fjet2_pt          = Float_t(jet2_pt);
		fjet2_btag        = Float_t(jet2_btag);
		fjet3_pt          = Float_t(jet3_pt);
		fjet3_btag        = Float_t(jet3_btag);

		// GEN-RECO match
		//deltaRleplep = deltaR(etalep1,philep1,etalep2,philep2);
		//deltaRWlepGen = deltaR(etaGenVlep, phiGenVlep, yVlep, phiVlep);
		//jetV.SetPtEtaPhiM(ptVhad, yVhad, phiVhad, massVhad);
		//genjetV.SetPtEtaPhiM(ptGenVhad, etaGenVhad, phiGenVhad, massGenVhad);
		//deltaRWhadGen = deltaR(etaGenVhad, phiGenVhad, yVhad, phiVhad);
		Double_t deltaRWhadGen = sqrt(pow(etaGenVhad-yVhad,2) + pow(phiGenVhad-phiVhad,2));

		//Weight Calculation
		Int_t bin = hR1->FindBin(npT);
		pileupWeight = hR1->GetBinContent(bin);		

		eff_and_pu_Weight = 0;
		eff_and_pu_Weight1 = 0;
		if(IsData>0) {
			if(npT < weights_pu1.size()){
				eff_and_pu_Weight = weights_pu1[npT];
			}
			if(npT < weights_pu2.size()){
				eff_and_pu_Weight1 = weights_pu2[npT];
			}
		}
		//cout << "pileupWeight:"<<pileupWeight<< " eff_and_pu_Weight:" << eff_and_pu_Weight << " eff_and_pu_Weight1:" << eff_and_pu_Weight1 << endl;
		if(theWeight>0) nn=1;
		else nn= -1;
		if(npp>0) lumiWeight=Identical_lumiWeight/(npp-nmm);
		else lumiWeight=Identical_lumiWeight/nentries;
		weight=lumiWeight*triggerWeight*eff_and_pu_Weight*nn;
		//weight=lumiWeight*triggerWeight*pileupWeight*nn;
		if (IsData>1 ) weight = weight*1.21;

		//lumiWeight=Identical_lumiWeight;
		//if(npp>0) weight=lumiWeight*triggerWeight*pileupWeight/(npp-nmm)*nn*0.04024;//0.00559;
		//else weight=lumiWeight*triggerWeight*pileupWeight/nentries*0.04024;//0.00559;
		if ( IsData==0 ) weight=1;
		//Weight Calculation Done

		//number of bjet calculation
		num_bJet=0.;
		num_bJet_loose=0.;
		num_bJet_tight=0.;
		for(Int_t i=0; i<8; i++)  {
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.890 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8[i]>=0.8 ) {num_bJet=num_bJet+1;}
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.605 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8[i]>=0.8 ) {num_bJet_loose=num_bJet_loose+1;}
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.970 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8[i]>=0.8 ) {num_bJet_tight=num_bJet_tight+1;}
		}
		nbtag=num_bJet;
		//number of bjet calculation Done

		Int_t nLooseLep=nLooseEle+nLooseMu;//the tight Lep included

		Double_t isAnaHP=1.;
		Double_t isAnaLP=1.;
		Double_t isTTBarControl=1.;
		Int_t tmp_categoryID_channel=0;
		if( channelname=="el" ){
		tmp_categoryID_channel=-1;// -1 for el; 1 for mu

			//HP: 0<tau21<=0.5;
			if (isAnaHP>0 && lep==11 && HLT_Ele105>0 && isHEEP>0 && nLooseLep==1){ nID_e = nID_e+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptlep1>120 && fabs(etalep1)<2.5){ npt_e = npt_e+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && MET_et>80) { nmet_e = nmet_e+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptVlepJEC > 200.) { nptVlepJEC = nptVlepJEC +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0 ){ nptVhad = nptVhad+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && num_bJet<1){ nnum_bJet_e = nnum_bJet_e +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && deltaRlepjet>pi_2) {n_deltaRlepjet = n_deltaRlepjet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetmet) >2.0)  {n_delPhijetmet = n_delPhijetmet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetlep)>2.0) {n_delPhijetlep = n_delPhijetlep +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && tau21>0. && tau21<0.5) {ntau = ntau+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && massVhadJEC>30 && massVhadJEC <150)// && m_lvj>200 && m_lvj<5000)
			{
				nmassVhad = nmassVhad +1;
				yields = yields + weight;
				(*file_cutflow)<<"event:"<<event<<endl;
			} else{ isAnaHP=-1; }

			//LP: 0.5<tau21<=0.75;
			if ( lep==11 && HLT_Ele105>0 && isHEEP>0 && nLooseLep==1 && ptlep1>120 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200.  && ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0 && num_bJet<1 && deltaRlepjet>pi_2 && fabs(delPhijetmet) >2.0 && fabs(delPhijetlep)>2.0 && tau21>0.5 && tau21<=0.75 && (( massVhadJEC >30 &&  massVhadJEC< 150 )) )
			{ isAnaLP=1.; } 
			else{ isAnaLP=-1.; }

			//TTbar control
			if ( lep==11 && HLT_Ele105>0 && isHEEP>0 && nLooseLep==1 && ptlep1>120 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200. &&ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0  && num_bJet>0 &&  massVhadJEC>30 && massVhadJEC <150)
			{ isTTBarControl=1.; } 
			else{ isTTBarControl=-1.; }
		}
		else if( channelname=="mu" ){
		tmp_categoryID_channel=1;// -1 for el; 1 for mu
			//HP: 0<tau21<=0.5;
			if (isAnaHP>0 && lep==13 && HLT_Mu45_v1>0 && isHighPt>0 && trackIso/ptlep1<0.1 && fabs(etalep1)<2.1 && nLooseLep==1 ) { nID_mu = nID_mu+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptlep1>53){ npt_mu = npt_mu+1; }else{ isAnaHP=-1; }
			if (isAnaHP>0 && MET_et>40) { nmet_mu = nmet_mu+1; }else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptVlepJEC>200) { nptVlepJEC = nptVlepJEC +1;} else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0){ nptVhad = nptVhad+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && num_bJet<1){ nnum_bJet_mu = nnum_bJet_mu +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && deltaRlepjet>pi_2) {n_deltaRlepjet = n_deltaRlepjet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetmet) >2.0)  {n_delPhijetmet = n_delPhijetmet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetlep)>2.0) {n_delPhijetlep = n_delPhijetlep +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && tau21>0. && tau21<0.5) {ntau = ntau+1;} else{ isAnaHP=-1; }
			if (isAnaHP>0 && (( massVhadJEC >30&& massVhadJEC< 150 )))// && m_lvj>100 && m_lvj<5000 )
			{ 
				nmassVhad = nmassVhad +1; 
				(*file_cutflow)<<"event:"<<event<<endl;
			} else{ isAnaHP=-1; }

			//LP: 0.5<tau21<=0.75;
			if (lep==13 && HLT_Mu45_v1>0 && isHighPt>0 && trackIso/ptlep1<0.1 && fabs(etalep1)<2.1 && nLooseLep==1 && ptlep1>53 && MET_et>40 && ptVlepJEC>200 && ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0 && num_bJet<1 && deltaRlepjet>pi_2 && fabs(delPhijetmet) >2.0 && fabs(delPhijetlep)>2.0 && tau21>0.5 && tau21<0.75 && (( massVhadJEC >30&& massVhadJEC< 150 )))
			{ isAnaLP=1.; } 
			else{ isAnaLP=-1.; }

			//TTbar control
			if (lep==13 && HLT_Mu45_v1>0 && isHighPt>0 && trackIso/ptlep1<0.1 && fabs(etalep1)<2.1 && nLooseLep==1 && ptlep1>53 && MET_et>40 && ptVlepJEC>200 && ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0 && num_bJet>0 && massVhadJEC>30 && massVhadJEC <150)
			{ isTTBarControl=1.; } 
			else{ isTTBarControl=-1.; }

		}else{
			cout<<"We don't know channelname:"<<channelname<<endl;
		}

		Int_t tmp_categoryID_eventselection=0;
		if(isAnaHP>0)tmp_categoryID_eventselection=1;
		else if(isAnaLP>0)tmp_categoryID_eventselection=2;
		else if(isTTBarControl>0)tmp_categoryID_eventselection=3;
		else tmp_categoryID_eventselection=100;

		CategoryID=tmp_categoryID_channel* tmp_categoryID_eventselection;

		isMatch=1.;
		if(deltaRWhadGen >= 0.3) isMatch=-1;
		//cout << "massVhad" << massVhad << "jet_mass_pruned " << jet_mass_pruned << endl;
		if(tau21<=0){vTagID=2;}
		else if(tau21>0.45 && tau21<=0.60){vTagID=1;}
		else if(tau21>0.60 && tau21<=0.75){vTagID=0;}
		else if(tau21>0.75 && tau21<=1){vTagID=-1;}
		else {vTagID=-2;}

		if(TMath::Abs(CategoryID)<10) ExTree->Fill();
	}

	if(channelname=="el"){ 
		std::cout << "nID_e" << nID_e << "; npt_e" << npt_e << "; nmet_e" << nmet_e << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<<"; nnum_bJet_e" << nnum_bJet_e <<"; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad << "; number_qq" << number_qq << "; yields " << yields << std::endl;
		(*file_cutflow) << "nID_e" << nID_e << "; npt_e" << npt_e << "; nmet_e" << nmet_e << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<<"; nnum_bJet_e" << nnum_bJet_e <<"; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad << "; number_qq" << number_qq << std::endl;
	}
	if(channelname=="mu"){
		std::cout << "nID_mu" << nID_mu << "; npt_mu" << npt_mu << "; nmet_mu" << nmet_mu << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<< "; nnum_bJet_mu" << nnum_bJet_mu << "; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad<< "; number_qq" << number_qq << std::endl;
		(*file_cutflow) << "nID_mu" << nID_mu << "; npt_mu" << npt_mu << "; nmet_mu" << nmet_mu << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<< "; nnum_bJet_mu" << nnum_bJet_mu << "; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad<< "; number_qq" << number_qq << std::endl;
	}

}
