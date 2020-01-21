{
  gROOT->Reset();

  TFile* fin = new TFile("../inputs/TwoMuZH_nSelectedAODCaloJetTag_log.root","READ");
  TH1F* bkgtotal = (TH1F*)fin->Get("bkgtotal");
  TH1F* light = (TH1F*)fin->Get("light");
  TH1F* heavy = (TH1F*)fin->Get("heavy");
  TH1F* other = (TH1F*)fin->Get("other");
  TH1F* Sig_MS55ct100 = (TH1F*)fin->Get("Sig_MS55ct100");
  TH1F* Data = (TH1F*)fin->Get("Data");
  Data->Add(Sig_MS55ct100,5.0);
  
  TFile* fout = new TFile("fout.root", "RECREATE");
  bkgtotal->Write("bkgtotal");
  light->Write("light");
  heavy->Write("heavy");
  other->Write("other");
  Sig_MS55ct100->Write("Sig_MS55ct100");
  Data->Write("Data");
  fout->Close();
}
