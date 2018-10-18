void build_ws(){

  //Options
  bool blind = true;
  TString signal_string = "Sig_MS40ct100";

  //Some setup
  TString data_string = "Data";
  if(blind) data_string = "bkgtotal";

  // As usual, load the combine library to get access to the RooParametricHist
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");

  // Output file and workspace
  TFile *fOut = new TFile("param_ws.root","RECREATE");
  RooWorkspace wspace("wspace","wspace");

  //Search in ntags, define ntags as our variable
  RooRealVar ntags("ntags", "ntags", -.5, 2.5);
  RooArgList vars(ntags);


  //---------------------------------------------------------------------------------------------------------------
  //EleMu Control Region
  //---------------------------------------------------------------------------------------------------------------
  TFile* f_elemu = TFile::Open("../inputs/EleMuOSOF_nSelectedAODCaloJetTag_GH.root", "READ");
  TH1F* data_elemu_th1_file = (TH1F*)f_elemu->Get(data_string);
  
  TH1F data_elemu_th1("data_obs_elemu","Data observed in EleMu", 3, -0.5, 2.5);
  data_elemu_th1.SetBinContent(1, data_elemu_th1_file->GetBinContent(1));
  data_elemu_th1.SetBinContent(2, data_elemu_th1_file->GetBinContent(2));
  data_elemu_th1.SetBinContent(3, data_elemu_th1_file->Integral(3,6));//assumes bin 6 is overflow TODO

  RooDataHist data_elemu_hist("data_obs_elemu", "Data observed in EleMu", vars, &data_elemu_th1);
  wspace.import(data_elemu_hist);  

  //Heavy background is freely floating
  //Create one parameter per bin representing the yield 
  RooRealVar heavy_elemu_bin1("heavy_elemu_bin1", "Heavy background yield in EleMu, bin 1", 5126.56, 1e3, 5e4);
  RooRealVar heavy_elemu_bin2("heavy_elemu_bin2", "Heavy background yield in EleMu, bin 2", 98.25, 50, 200);
  RooRealVar heavy_elemu_bin3("heavy_elemu_bin3", "Heavy background yield in EleMu, bin 3", 1.29, 0, 5);
  RooArgList heavy_elemu_bins;
  heavy_elemu_bins.add(heavy_elemu_bin1);
  heavy_elemu_bins.add(heavy_elemu_bin2);
  heavy_elemu_bins.add(heavy_elemu_bin3);

  //Create a RooParametricHist which contains yields.  Last argument is just for binning
  // better to have ntags as vars?? 
  RooParametricHist p_heavy_elemu("heavy_elemu", "Heavy PDF in EleMu Control Region", ntags, heavy_elemu_bins, data_elemu_th1);
  //Always include _norm term which should be sum of yields
  RooAddition p_heavy_elemu_norm("heavy_elemu_norm", "Total number of heavy events in EleMu Control Region", heavy_elemu_bins);

  //For now at least assume this region is 100% heavy background TODO


  //---------------------------------------------------------------------------------------------------------------
  //OnePho Control Region
  //---------------------------------------------------------------------------------------------------------------
  TFile* f_onepho = TFile::Open("../inputs/OnePho_nSelectedAODCaloJetTag_GH.root", "READ");
  TH1F* data_onepho_th1_file = (TH1F*)f_onepho->Get(data_string);
  
  TH1F data_onepho_th1("data_obs_onepho","Data observed in OnePho", 3, -0.5, 2.5);
  data_onepho_th1.SetBinContent(1, data_onepho_th1_file->GetBinContent(1));
  data_onepho_th1.SetBinContent(2, data_onepho_th1_file->GetBinContent(2));
  data_onepho_th1.SetBinContent(3, data_onepho_th1_file->Integral(3,6));//assumes bin 6 is overflow

  RooDataHist data_onepho_hist("data_obs_onepho", "Data observed in OnePho", vars, &data_onepho_th1);
  wspace.import(data_onepho_hist);

  //Light background is freely floating
  //Create one parameter per bin representing the yield 
  RooRealVar light_onepho_bin1("light_onepho_bin1", "Light background yield in OnePho, bin 1", 288366, 200000, 400000);
  RooRealVar light_onepho_bin2("light_onepho_bin2", "Light background yield in OnePho, bin 2", 1377.44, 500, 2000);
  RooRealVar light_onepho_bin3("light_onepho_bin3", "Light background yield in OnePho, bin 3", 0.5, 0, 5);
  RooArgList light_onepho_bins;
  light_onepho_bins.add(light_onepho_bin1);
  light_onepho_bins.add(light_onepho_bin2);
  light_onepho_bins.add(light_onepho_bin3);
  RooParametricHist p_light_onepho("light_onepho", "Light PDF in OnePho Control Region", ntags, light_onepho_bins, data_onepho_th1);
  RooAddition p_light_onepho_norm("light_onepho_norm", "Total number of light events in OnePho Control Region", light_onepho_bins);
  
  //Now add heavy_onepho, function of heavy_elemu
  
  //Start with transfer factor -- assume all bins share one
  RooRealVar rrv_heavy_elemu_to_onepho("rrv_heavy_elemu_to_onepho", "heavy EleMu to OnePho",1);
  RooFormulaVar tf_heavy_elemu_to_onepho("tf_heavy_elemu_to_onepho", "heavy EleMu to OnePho transfer factor", " 0.138834*TMath::Power(1+0.5,@0)", RooArgList(rrv_heavy_elemu_to_onepho) );//Leave as gaussian for now TODO
  
  RooFormulaVar heavy_onepho_bin1("heavy_onepho_bin1", "Heavy background yield in OnePho, bin 1", "@0*@1", RooArgList(tf_heavy_elemu_to_onepho, heavy_elemu_bin1));
  RooFormulaVar heavy_onepho_bin2("heavy_onepho_bin2", "Heavy background yield in OnePho, bin 2", "@0*@1", RooArgList(tf_heavy_elemu_to_onepho, heavy_elemu_bin2));
  RooFormulaVar heavy_onepho_bin3("heavy_onepho_bin3", "Heavy background yield in OnePho, bin 3", "@0*@1", RooArgList(tf_heavy_elemu_to_onepho, heavy_elemu_bin3));
  RooArgList heavy_onepho_bins;
  heavy_onepho_bins.add(heavy_onepho_bin1);
  heavy_onepho_bins.add(heavy_onepho_bin2);
  heavy_onepho_bins.add(heavy_onepho_bin3);
  RooParametricHist p_heavy_onepho("heavy_onepho", "Heavy PDF in OnePho Control Region", ntags, heavy_onepho_bins, data_onepho_th1);
  RooAddition p_heavy_onepho_norm("heavy_onepho_norm", "Total number of heavy events in OnePho Control Region", heavy_onepho_bins);

  //---------------------------------------------------------------------------------------------------------------
  //TwoMuZH Signal Region
  //---------------------------------------------------------------------------------------------------------------

  //Data
  TFile* f_twomuzh = TFile::Open("../inputs/TwoMuZH_nSelectedAODCaloJetTag_GH.root", "READ");
  TH1F* data_twomuzh_th1_file = (TH1F*)f_twomuzh->Get(data_string);
  
  TH1F data_twomuzh_th1("data_obs_twomuzh","Data observed in TwoMuZH", 3, -0.5, 2.5);
  data_twomuzh_th1.SetBinContent(1, data_twomuzh_th1_file->GetBinContent(1));
  data_twomuzh_th1.SetBinContent(2, data_twomuzh_th1_file->GetBinContent(2));
  data_twomuzh_th1.SetBinContent(3, data_twomuzh_th1_file->Integral(3,6));//assumes bin 6 is overflow

  RooDataHist data_twomuzh_hist("data_obs_twomuzh", "Data observed in TwoMuZH", vars, &data_twomuzh_th1);
  wspace.import(data_twomuzh_hist);

  //Heavy as function of heavy_elemu
  RooRealVar rrv_heavy_elemu_to_twomuzh("rrv_heavy_elemu_to_twomuzh", "heavy EleMu to TwoMuZH",1);
  RooFormulaVar tf_heavy_elemu_to_twomuzh("tf_heavy_elemu_to_twomuzh", "heavy EleMu to TwoMuZH transfer factor", "0.733396*TMath::Power(1+0.5,@0)", RooArgList(rrv_heavy_elemu_to_twomuzh) );//Leave as gaussian for now TODO
  
  RooFormulaVar heavy_twomuzh_bin1("heavy_twomuzh_bin1", "Heavy background yield in TwoMuZH, bin 1", "@0*@1", RooArgList(tf_heavy_elemu_to_twomuzh, heavy_elemu_bin1));
  RooFormulaVar heavy_twomuzh_bin2("heavy_twomuzh_bin2", "Heavy background yield in TwoMuZH, bin 2", "@0*@1", RooArgList(tf_heavy_elemu_to_twomuzh, heavy_elemu_bin2));
  RooFormulaVar heavy_twomuzh_bin3("heavy_twomuzh_bin3", "Heavy background yield in TwoMuZH, bin 3", "@0*@1", RooArgList(tf_heavy_elemu_to_twomuzh, heavy_elemu_bin3));
  RooArgList heavy_twomuzh_bins;
  heavy_twomuzh_bins.add(heavy_twomuzh_bin1);
  heavy_twomuzh_bins.add(heavy_twomuzh_bin2);
  heavy_twomuzh_bins.add(heavy_twomuzh_bin3);
  RooParametricHist p_heavy_twomuzh("heavy_twomuzh", "Heavy PDF in TwoMuZH Control Region", ntags, heavy_twomuzh_bins, data_twomuzh_th1);
  RooAddition p_heavy_twomuzh_norm("heavy_twomuzh_norm", "Total number of heavy events in TwoMuZH Control Region", heavy_twomuzh_bins);

  //Light as function of light_onepho
  RooRealVar rrv_light_onepho_to_twomuzh("rrv_light_onepho_to_twomuzh", "light OnePhoo to TwoMuZH",1);
  RooFormulaVar tf_light_onepho_to_twomuzh("tf_light_onepho_to_twomuzh", "light OnePhoo to TwoMuZH transfer factor", "0.576097*TMath::Power(1+0.5,@0)", RooArgList(rrv_light_onepho_to_twomuzh) );//Leave as gaussian for now TODO
  
  RooFormulaVar light_twomuzh_bin1("light_twomuzh_bin1", "Light background yield in TwoMuZH, bin 1", "@0*@1", RooArgList(tf_light_onepho_to_twomuzh, light_onepho_bin1));
  RooFormulaVar light_twomuzh_bin2("light_twomuzh_bin2", "Light background yield in TwoMuZH, bin 2", "@0*@1", RooArgList(tf_light_onepho_to_twomuzh, light_onepho_bin2));
  RooFormulaVar light_twomuzh_bin3("light_twomuzh_bin3", "Light background yield in TwoMuZH, bin 3", "@0*@1", RooArgList(tf_light_onepho_to_twomuzh, light_onepho_bin3));
  RooArgList light_twomuzh_bins;
  light_twomuzh_bins.add(light_twomuzh_bin1);
  light_twomuzh_bins.add(light_twomuzh_bin2);
  light_twomuzh_bins.add(light_twomuzh_bin3);
  RooParametricHist p_light_twomuzh("light_twomuzh", "Light PDF in TwoMuZH Control Region", ntags, light_twomuzh_bins, data_twomuzh_th1);
  RooAddition p_light_twomuzh_norm("light_twomuzh_norm", "Total number of light events in TwoMuZH Control Region", light_twomuzh_bins);

  //Signal
  TH1F* signal_twomuzh_th1_file = (TH1F*)f_twomuzh->Get(signal_string);
  
  TH1F signal_twomuzh_th1("signal_twomuzh","Signal yield in TwoMuZH", 3, -0.5, 2.5);
  signal_twomuzh_th1.SetBinContent(1, signal_twomuzh_th1_file->GetBinContent(1));
  signal_twomuzh_th1.SetBinContent(2, signal_twomuzh_th1_file->GetBinContent(2));
  signal_twomuzh_th1.SetBinContent(3, signal_twomuzh_th1_file->Integral(3,6));//assumes bin 6 is overflow

  RooDataHist signal_twomuzh_hist("signal_twomuzh", "Signal yield in TwoMuZH", vars, &signal_twomuzh_th1);
  wspace.import(signal_twomuzh_hist);

  //---------------------------------------------------------------------------------------------------------------
  //TwoEeleZH Signal Region
  //---------------------------------------------------------------------------------------------------------------
  //Skip for now

 
  // Final stuff -- import RooParametricHists and _norms
  //EleMu
  wspace.import(p_heavy_elemu);
  wspace.import(p_heavy_elemu_norm, RooFit::RecycleConflictNodes());
  //OnePho
  wspace.import(p_light_onepho);
  wspace.import(p_light_onepho_norm, RooFit::RecycleConflictNodes());
  wspace.import(p_heavy_onepho);
  wspace.import(p_heavy_onepho_norm, RooFit::RecycleConflictNodes());
  //TwoMuZH
  wspace.import(p_light_twomuzh);
  wspace.import(p_light_twomuzh_norm, RooFit::RecycleConflictNodes());
  wspace.import(p_heavy_twomuzh);
  wspace.import(p_heavy_twomuzh_norm, RooFit::RecycleConflictNodes());

  fOut->cd();
  wspace.Write();

}
