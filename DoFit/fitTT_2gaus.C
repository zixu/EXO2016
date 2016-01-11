
void fitTT_2gaus(){

    gSystem->Load("libRooFit");
    using namespace RooFit;

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
//    RooRealVar c_ErfExp("c_ErfExp", "c_ErfExp", -0.1, -10., 0.);
//    RooRealVar off_ErfExp("off_ErfExp", "off_ErfExp", 50., 10., 100.);
//    RooRealVar sigma_ErfExp("sigma_ErfExp", "sigma_ErfExp", 10., 0., 100.);
//    RooErfExpPdf erfExp("erfExp", "erfExp", x, c_ErfExp, off_ErfExp, sigma_ErfExp);
//    RooRealVar a0("a0","a0",0.5,0.,1.) ;
//    RooRealVar a1("a1","a1",0.2,0.,1.) ;
//    RooChebychev cheby("cheby", "cheby", x, RooArgSet(a0, a1));
     
    RooRealVar mean1("mean1", "mean1", 80, 70, 150);
    RooRealVar sigma1("sigma1", "sigma1", 5., 1., 100);
    RooGaussian gaus("gaus", "gaus", x, mean1, sigma1);

    RooRealVar mean2("mean2", "mean2", 80, 70, 150);
    RooRealVar sigma2("sigma2", "sigma2", 5., 1., 100);
    RooGaussian gaus2("gaus2", "gaus2", x, mean2, sigma2);

    RooRealVar  f("f", "f", 0.5, 0., 1.);
    //RooAddPdf model("model", "model", RooArgList(erfExp, gaus), f);
    RooAddPdf model("model", "model", RooArgList(gaus, gaus2), f);

    //RooDataSet data("data", "data", t, x);
    RooDataSet data("data", "data",  x);
    for(Int_t i=0; i< t->GetEntries(); i++){
	t->GetEntry(i);
	//if( jet_tau2tau1>0.5 && isMatch>0){
	if(isMatch>0 && jet_tau2tau1<0.5){
	//if(isMatch>0 && jet_tau2tau1>0.5){
	//if(isMatch<0 && jet_tau2tau1<0.5){
	//if(isMatch<0 && jet_tau2tau1>0.5){
	    x = massVhad;
            data.add(x);
	}
    }

//    model.fitTo(data);
    model.fitTo(data);

    RooPlot * frame = x.frame();
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components(gaus2), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Components(gaus), LineStyle(kDashed), LineColor(kBlue));

    TCanvas * c = new TCanvas("c","c");
    frame->Draw();
    c->Print("TT_pass_match.png");
}



