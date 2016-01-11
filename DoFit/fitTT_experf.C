
void fitTT_experf(){

    gSystem->Load("libRooFit");
    using namespace RooFit;
    gROOT->ProcessLine(".L ./PDFs/HWWLVJRooPdfs_cxx.so");

    //TFile * file = new TFile("TTBarDataSet/el_PKUTree_pdata.root");
    TFile * file = new TFile("TTBarDataSet/el_out_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola.root");
    TTree * t = (TTree *)file->Get("PKUTree");

    float jet_tau2tau1,  jet_mass_softdrop;
    double massVhad;
    double isMatch;

    t->SetBranchAddress("massVhad", &massVhad);
    t->SetBranchAddress("jet_mass_softdrop", & jet_mass_softdrop);
    t->SetBranchAddress("jet_tau2tau1", &jet_tau2tau1);
    t->SetBranchAddress("isMatch", &isMatch);

    RooRealVar x("massVhad","massVhad", 40, 130);
    RooRealVar c_ErfExp("c_ErfExp", "c_ErfExp", -0.1, -10., 0.);
    RooRealVar off_ErfExp("off_ErfExp", "off_ErfExp", 50., 10., 100.);
    RooRealVar sigma_ErfExp("sigma_ErfExp", "sigma_ErfExp", 10., 0., 100.);
    RooErfExpPdf erfExp("erfExp", "erfExp", x, c_ErfExp, off_ErfExp, sigma_ErfExp);
     
    RooRealVar mean1("mean1", "mean1", 75, 60, 90);
    RooRealVar sigma1("sigma1", "sigma1", 16., 1., 40);
    RooGaussian gaus("gaus", "gaus", x, mean1, sigma1);

    RooRealVar  f("f", "f", 0.5, 0., 1.);
    //RooAddPdf model("model", "model", RooArgList(erfExp, gaus), f);
    RooAddPdf model("model", "model", RooArgList(erfExp, gaus), f);

    //RooDataSet data("data", "data", t, x);
    RooDataSet data("data", "data",  x);
    for(Int_t i=0; i< t->GetEntries(); i++){
	t->GetEntry(i);
	//if( jet_tau2tau1>0.5 && isMatch>0){
	//if(isMatch>0 && jet_tau2tau1<0.5){
	//if(isMatch>0 && jet_tau2tau1>0.5){
	//if(isMatch<0 && jet_tau2tau1<0.5){
	if(isMatch<0 && jet_tau2tau1>0.5){
	    x = massVhad;
            data.add(x);
	}
    }

//    model.fitTo(data);
    model.fitTo(data);

    RooPlot * frame = x.frame();
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components(erfExp), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Components(gaus), LineStyle(kDashed), LineColor(kBlue));

    TCanvas * c = new TCanvas("c","c");
    frame->Draw();
    //c->Print("TT_pass_nomatch.png");
    c->Print("TT_fail_nomatch.png");
}



