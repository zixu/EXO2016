void fitTT_experf(TString filename, double tau21min, double tau21max, double match, TString outname){

    gSystem->Load("libRooFit");
    using namespace RooFit;
    gROOT->ProcessLine(".L ./PDFs/HWWLVJRooPdfs_cxx.so");

    TFile * file = new TFile(filename);
    //TFile * file = new TFile("TTBarDataSet/el_out_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola.root");
    //TFile * file = new TFile("TTBarDataSet/mu_out_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
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
    RooRealVar off_ErfExp("off_ErfExp", "off_ErfExp", 50., 0., 140.);
    RooRealVar sigma_ErfExp("sigma_ErfExp", "sigma_ErfExp", 25., 0., 200.);
    RooErfExpPdf model("model", "model", x, c_ErfExp, off_ErfExp, sigma_ErfExp);

    //RooRealVar mean1("mean1", "mean1", 90, 75, 105);
    //RooRealVar sigma1("sigma1", "sigma1", 23., 1., 40);
    //RooGaussian gaus("gaus", "gaus", x, mean1, sigma1);

    //RooRealVar  f("f", "f", 0.5, 0., 1.);
    //RooAddPdf model("model", "model", RooArgList( gaus, erfExp), f);

    RooDataSet data("data", "data",  x);
    for(Int_t i=0; i< t->GetEntries(); i++){
        t->GetEntry(i);
        if(isMatch==match && jet_tau2tau1>tau21min && jet_tau2tau1< tau21max){
            x = massVhad;
            data.add(x);
        }
    }

    model.fitTo(data);

    RooPlot * frame = x.frame();
    data.plotOn(frame);
    model.plotOn(frame);
    //model.plotOn(frame, Components(erfExp), LineStyle(kDashed), LineColor(kRed));
    //model.plotOn(frame, Components(gaus), LineStyle(kDashed), LineColor(kBlue));

    TCanvas * c = new TCanvas("c","c");
    frame->Draw();
    c->Print(outname+"TT.png");
}

void fitTT_2gaus(TString filename, double tau21min, double tau21max, double match){

    gSystem->Load("libRooFit");
    using namespace RooFit;

    TFile * file = new TFile(filename);
    //TFile * file = new TFile("TTBarDataSet/el_out_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola.root");
    TTree * t = (TTree *)file->Get("PKUTree");

    float jet_tau2tau1,  jet_mass_softdrop;
    double massVhad;
    double isMatch;

    t->SetBranchAddress("massVhad", &massVhad);
    t->SetBranchAddress("jet_mass_softdrop", & jet_mass_softdrop);
    t->SetBranchAddress("jet_tau2tau1", &jet_tau2tau1);
    t->SetBranchAddress("isMatch", &isMatch);

    RooRealVar x("massVhad","massVhad", 40, 130);
    RooRealVar mean1("mean1", "mean1", 97, 50, 120);
    RooRealVar sigma1("sigma1", "sigma1", 28., 1., 150);
    RooGaussian gaus("gaus", "gaus", x, mean1, sigma1);

    RooRealVar mean2("mean2", "mean2", 79, 50, 100);
    RooRealVar sigma2("sigma2", "sigma2", 7., 1., 150);
    RooGaussian gaus2("gaus2", "gaus2", x, mean2, sigma2);

    RooRealVar  f("f", "f", 0.5, 0., 1.);
    RooAddPdf model("model", "model", RooArgList(gaus2, gaus), f);

    RooDataSet data("data", "data",  x);
    for(Int_t i=0; i< t->GetEntries(); i++){
        t->GetEntry(i);
        //if(isMatch>0 && jet_tau2tau1<0.5){
        if(isMatch==match && jet_tau2tau1>tau21min && jet_tau2tau1< tau21max){
            x = massVhad;
            data.add(x);
        }
    }

    model.fitTo(data);

    RooPlot * frame = x.frame();
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components(gaus2), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Components(gaus), LineStyle(kDashed), LineColor(kBlue));

    TCanvas * c = new TCanvas("c","c");
    frame->Draw();
    c->Print("pass_match_TT.png");
}

void fitTT_GausCheb(TString filename, double tau21min, double tau21max, double match){

    gSystem->Load("libRooFit");
    using namespace RooFit;

    TFile * file = new TFile(filename);
    //TFile * file = new TFile("TTBarDataSet/el_out_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola.root");
    TTree * t = (TTree *)file->Get("PKUTree");

    float jet_tau2tau1,  jet_mass_softdrop;
    double massVhad;
    double isMatch;

    t->SetBranchAddress("massVhad", &massVhad);
    t->SetBranchAddress("jet_mass_softdrop", & jet_mass_softdrop);
    t->SetBranchAddress("jet_tau2tau1", &jet_tau2tau1);
    t->SetBranchAddress("isMatch", &isMatch);

    RooRealVar x("massVhad","massVhad", 40, 130);
    RooRealVar a0("a0","a0",-0.3,-2.,2.) ;
    RooRealVar a1("a1","a1",0.2,-1.,1.) ;
    RooChebychev cheby("cheby", "cheby", x, RooArgSet(a0, a1));
     
    RooRealVar mean1("mean1", "mean1", 76, 50, 90);
    RooRealVar sigma1("sigma1", "sigma1", 11., 1., 40);
    RooGaussian gaus("gaus", "gaus", x, mean1, sigma1);

    RooRealVar  f("f", "f", 0.5, 0., 1.);
    RooAddPdf model("model", "model", RooArgList( gaus, cheby), f);

    RooDataSet data("data", "data",  x);
    for(Int_t i=0; i< t->GetEntries(); i++){
	t->GetEntry(i);
        if(isMatch==match && jet_tau2tau1>tau21min && jet_tau2tau1< tau21max){
	//if(isMatch>0 && jet_tau2tau1>0.5){
	//if(isMatch<0 && jet_tau2tau1<0.5){
	//if(isMatch<0 && jet_tau2tau1>0.5){
	    x = massVhad;
            data.add(x);
	}
    }

    model.fitTo(data);

    RooPlot * frame = x.frame();
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components(cheby), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Components(gaus), LineStyle(kDashed), LineColor(kBlue));

    TCanvas * c = new TCanvas("c","c");
    frame->Draw();
    c->Print("fail_match_TT.png");
}

void fitTT(){
    gSystem->Load("libRooFit");
    using namespace RooFit;
    gROOT->ProcessLine(".L ./PDFs/HWWLVJRooPdfs_cxx.so");
    TString openfilename;
    openfilename = "TTBarDataSet/mu_out_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root";
    //openfilename = "TTBarDataSet/mu_out_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola.root";
    fitTT_experf(openfilename, 0., 0.5, -1., "pass_nomatch_");
    fitTT_experf(openfilename, 0.5, 1.0, -1., "fail_nomatch_");
    fitTT_2gaus(openfilename, 0., 0.5, 1.);
    fitTT_GausCheb(openfilename, 0.5, 1.0, 1.);

}

