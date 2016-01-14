void MakePseudoData(TString channel) {

	Double_t lumi=1.;

	//Get old file, old tree and set top branch address
	TFile *oldfile = new TFile(channel+"_PKUTree_allBkg_xww.root");
	TTree *oldtree = (TTree*)oldfile->Get("PKUTree");
	Long64_t nentries = oldtree->GetEntries();

	Double_t weight=-100;
	oldtree->SetBranchAddress("weight",&weight);

	//Create a new file + a clone of old tree in new file
	TFile *newfile = new TFile(channel+"_PKUTree_pdata.root","recreate");
	TTree *newtree = oldtree->CloneTree(0);

	TRandom3 r3(0);

	for (Long64_t i=0;i<nentries; i++) {
		oldtree->GetEntry(i);
		//if (event->GetNtrack() > 605) newtree->Fill();
		//event->Clear();
		Double_t tmp=r3.Uniform(0.,1.);
		if (weight*lumi>tmp){
			weight=1.;
			newtree->Fill(); 
			cout<<"weight="<<weight<<" r3="<<tmp<<endl;
		}
		weight=-100;
	}
	newtree->Print();
	newtree->AutoSave();
	delete oldfile;
	delete newfile;
}
