#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TSystem.h"

#include "RooArgList.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAddition.h"
#include "RooDataHist.h"

#include "HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"


//Global options
TString signal_string = "Sig_MS40ct100";
TString data_string = "bkgtotal"; //"Data";


//---------------------------------------------------------------------------------------------------------------
//TwoMuZH Region
//---------------------------------------------------------------------------------------------------------------
void build_twomuzh(RooWorkspace* wspace, TString light_est = "OnePho"){

  RooRealVar* ntags = wspace->var("ntags");
  RooArgList vars(*ntags);


  //Data
  TFile* f_twomuzh = TFile::Open("../inputs/TwoMuZH_nSelectedAODCaloJetTag_GH.root", "READ");
  TH1F* data_twomuzh_th1_file = (TH1F*)f_twomuzh->Get(data_string);
  TH1F data_twomuzh_th1("data_obs_twomuzh","Data observed in TwoMuZH", 3, -0.5, 2.5);
  data_twomuzh_th1.SetBinContent(1, data_twomuzh_th1_file->GetBinContent(1));
  data_twomuzh_th1.SetBinContent(2, data_twomuzh_th1_file->GetBinContent(2));
  data_twomuzh_th1.SetBinContent(3, data_twomuzh_th1_file->Integral(3,6));//assumes bin 6 is overflow
  RooDataHist data_twomuzh_hist("data_obs_twomuzh", "Data observed in TwoMuZH", vars, &data_twomuzh_th1);
  wspace->import(data_twomuzh_hist);


  //Heavy as function of heavy_elemu
  RooRealVar* heavy_elemu_bin1 = wspace->var("heavy_elemu_bin1");
  RooRealVar* heavy_elemu_bin2 = wspace->var("heavy_elemu_bin2");
  RooRealVar* heavy_elemu_bin3 = wspace->var("heavy_elemu_bin3");

  RooRealVar rrv_heavy_elemu_to_twomuzh_bin1("rrv_heavy_elemu_to_twomuzh_bin1", "heavy EleMu to TwoMuZH bin 1",1);
  RooRealVar rrv_heavy_elemu_to_twomuzh_bin2("rrv_heavy_elemu_to_twomuzh_bin2", "heavy EleMu to TwoMuZH bin 2",1);
  RooRealVar rrv_heavy_elemu_to_twomuzh_bin3("rrv_heavy_elemu_to_twomuzh_bin3", "heavy EleMu to TwoMuZH bin 3",1);

  RooFormulaVar tf_heavy_elemu_to_twomuzh_bin1("tf_heavy_elemu_to_twomuzh_bin1", "heavy EleMu to TwoMuZH transfer factor bin 1", 
					  "0.733396*TMath::Power(1+0.5,@0)", RooArgList(rrv_heavy_elemu_to_twomuzh_bin1) );
  RooFormulaVar tf_heavy_elemu_to_twomuzh_bin2("tf_heavy_elemu_to_twomuzh_bin2", "heavy EleMu to TwoMuZH transfer factor bin 2", 
					  "0.733396*TMath::Power(1+0.5,@0)", RooArgList(rrv_heavy_elemu_to_twomuzh_bin2) );
  RooFormulaVar tf_heavy_elemu_to_twomuzh_bin3("tf_heavy_elemu_to_twomuzh_bin3", "heavy EleMu to TwoMuZH transfer factor bin 3", 
					  "0.733396*TMath::Power(1+0.5,@0)", RooArgList(rrv_heavy_elemu_to_twomuzh_bin3) );
  
  RooFormulaVar heavy_twomuzh_bin1("heavy_twomuzh_bin1", "Heavy background yield in TwoMuZH, bin 1", 
				   "@0*@1", RooArgList(tf_heavy_elemu_to_twomuzh_bin1, *heavy_elemu_bin1));
  RooFormulaVar heavy_twomuzh_bin2("heavy_twomuzh_bin2", "Heavy background yield in TwoMuZH, bin 2", 
				   "@0*@1", RooArgList(tf_heavy_elemu_to_twomuzh_bin2, *heavy_elemu_bin2));
  RooFormulaVar heavy_twomuzh_bin3("heavy_twomuzh_bin3", "Heavy background yield in TwoMuZH, bin 3", 
				   "@0*@1", RooArgList(tf_heavy_elemu_to_twomuzh_bin3, *heavy_elemu_bin3));
  RooArgList heavy_twomuzh_bins;
  heavy_twomuzh_bins.add(heavy_twomuzh_bin1);
  heavy_twomuzh_bins.add(heavy_twomuzh_bin2);
  heavy_twomuzh_bins.add(heavy_twomuzh_bin3);
  RooParametricHist p_heavy_twomuzh("heavy_twomuzh", "Heavy PDF in TwoMuZH Region", *ntags, heavy_twomuzh_bins, data_twomuzh_th1);
  RooAddition p_heavy_twomuzh_norm("heavy_twomuzh_norm", "Total number of heavy events in TwoMuZH Region", heavy_twomuzh_bins);

  wspace->import(p_heavy_twomuzh);
  wspace->import(p_heavy_twomuzh_norm, RooFit::RecycleConflictNodes());


  //Light as function of light_onepho
  if(light_est=="OnePho"){
    RooRealVar* light_onepho_bin1 = wspace->var("light_onepho_bin1");
    RooRealVar* light_onepho_bin2 = wspace->var("light_onepho_bin2");
    RooRealVar* light_onepho_bin3 = wspace->var("light_onepho_bin3");

    RooRealVar rrv_light_onepho_to_twomuzh_bin1("rrv_light_onepho_to_twomuzh_bin1", "light OnePhoo to TwoMuZH_bin1",1);
    RooRealVar rrv_light_onepho_to_twomuzh_bin2("rrv_light_onepho_to_twomuzh_bin2", "light OnePhoo to TwoMuZH_bin2",1);
    RooRealVar rrv_light_onepho_to_twomuzh_bin3("rrv_light_onepho_to_twomuzh_bin3", "light OnePhoo to TwoMuZH_bin3",1);

    RooFormulaVar tf_light_onepho_to_twomuzh_bin1("tf_light_onepho_to_twomuzh_bin1", "light OnePhoo to TwoMuZH transfer factor bin 1", 
					     "0.576097*TMath::Power(1+0.5,@0)", RooArgList(rrv_light_onepho_to_twomuzh_bin1) );
    RooFormulaVar tf_light_onepho_to_twomuzh_bin2("tf_light_onepho_to_twomuzh_bin2", "light OnePhoo to TwoMuZH transfer factor bin 2", 
					     "0.576097*TMath::Power(1+0.5,@0)", RooArgList(rrv_light_onepho_to_twomuzh_bin2) );
    RooFormulaVar tf_light_onepho_to_twomuzh_bin3("tf_light_onepho_to_twomuzh_bin3", "light OnePhoo to TwoMuZH transfer factor bin 3", 
					     "0.576097*TMath::Power(1+0.5,@0)", RooArgList(rrv_light_onepho_to_twomuzh_bin3) );
    
    RooFormulaVar light_twomuzh_bin1("light_twomuzh_bin1", "Light background yield in TwoMuZH, bin 1", 
				     "@0*@1", RooArgList(tf_light_onepho_to_twomuzh_bin1, *light_onepho_bin1));
    RooFormulaVar light_twomuzh_bin2("light_twomuzh_bin2", "Light background yield in TwoMuZH, bin 2", 
				     "@0*@1", RooArgList(tf_light_onepho_to_twomuzh_bin2, *light_onepho_bin2));
    RooFormulaVar light_twomuzh_bin3("light_twomuzh_bin3", "Light background yield in TwoMuZH, bin 3", 
				     "@0*@1", RooArgList(tf_light_onepho_to_twomuzh_bin3, *light_onepho_bin3));
    RooArgList light_twomuzh_bins;
    light_twomuzh_bins.add(light_twomuzh_bin1);
    light_twomuzh_bins.add(light_twomuzh_bin2);
    light_twomuzh_bins.add(light_twomuzh_bin3);
    RooParametricHist p_light_twomuzh("light_twomuzh", "Light PDF in TwoMuZH Region", *ntags, light_twomuzh_bins, data_twomuzh_th1);
    RooAddition p_light_twomuzh_norm("light_twomuzh_norm", "Total number of light events in TwoMuZH Region", light_twomuzh_bins);
    
    wspace->import(p_light_twomuzh);
    wspace->import(p_light_twomuzh_norm, RooFit::RecycleConflictNodes());
  }
  else if(light_est == "DY"){
    RooRealVar* light_twomudy_bin1 = wspace->var("light_twomudy_bin1");
    RooRealVar* light_twomudy_bin2 = wspace->var("light_twomudy_bin2");
    RooRealVar* light_twomudy_bin3 = wspace->var("light_twomudy_bin3");

    RooRealVar rrv_light_twomudy_to_twomuzh_bin1("rrv_light_twomudy_to_twomuzh_bin1", "light TwoMuDY to TwoMuZH bin 1",1);
    RooRealVar rrv_light_twomudy_to_twomuzh_bin2("rrv_light_twomudy_to_twomuzh_bin2", "light TwoMuDY to TwoMuZH bin 2",1);
    RooRealVar rrv_light_twomudy_to_twomuzh_bin3("rrv_light_twomudy_to_twomuzh_bin3", "light TwoMuDY to TwoMuZH bin 3",1);

    RooFormulaVar tf_light_twomudy_to_twomuzh_bin1("tf_light_twomudy_to_twomuzh_bin1", "light TwoMuDY to TwoMuZH transfer factor bin 1", 
					      "0.576097*TMath::Power(1+0.5,@0)", RooArgList(rrv_light_twomudy_to_twomuzh_bin1) );
    RooFormulaVar tf_light_twomudy_to_twomuzh_bin2("tf_light_twomudy_to_twomuzh_bin2", "light TwoMuDY to TwoMuZH transfer factor bin 2", 
					      "0.576097*TMath::Power(1+0.5,@0)", RooArgList(rrv_light_twomudy_to_twomuzh_bin2) );
    RooFormulaVar tf_light_twomudy_to_twomuzh_bin3("tf_light_twomudy_to_twomuzh_bin3", "light TwoMuDY to TwoMuZH transfer factor bin 3", 
					      "0.576097*TMath::Power(1+0.5,@0)", RooArgList(rrv_light_twomudy_to_twomuzh_bin3) );
    
    RooFormulaVar light_twomuzh_bin1("light_twomuzh_bin1", "Light background yield in TwoMuZH, bin 1", 
				     "@0*@1", RooArgList(tf_light_twomudy_to_twomuzh_bin1, *light_twomudy_bin1));
    RooFormulaVar light_twomuzh_bin2("light_twomuzh_bin2", "Light background yield in TwoMuZH, bin 2", 
				     "@0*@1", RooArgList(tf_light_twomudy_to_twomuzh_bin2, *light_twomudy_bin2));
    RooFormulaVar light_twomuzh_bin3("light_twomuzh_bin3", "Light background yield in TwoMuZH, bin 3", 
				     "@0*@1", RooArgList(tf_light_twomudy_to_twomuzh_bin3, *light_twomudy_bin3));
    RooArgList light_twomuzh_bins;
    light_twomuzh_bins.add(light_twomuzh_bin1);
    light_twomuzh_bins.add(light_twomuzh_bin2);
    light_twomuzh_bins.add(light_twomuzh_bin3);
    RooParametricHist p_light_twomuzh("light_twomuzh", "Light PDF in TwoMuZH Region", *ntags, light_twomuzh_bins, data_twomuzh_th1);
    RooAddition p_light_twomuzh_norm("light_twomuzh_norm", "Total number of light events in TwoMuZH Region", light_twomuzh_bins);
    
    wspace->import(p_light_twomuzh);
    wspace->import(p_light_twomuzh_norm, RooFit::RecycleConflictNodes());
  }
  else {
    cout << "Invalid option for light estimate" << endl;
    return; 
  }
  
  //Signal
  TH1F* signal_twomuzh_th1_file = (TH1F*)f_twomuzh->Get(signal_string);
  TH1F signal_twomuzh_th1("signal_twomuzh","Signal yield in TwoMuZH", 3, -0.5, 2.5);
  signal_twomuzh_th1.SetBinContent(1, signal_twomuzh_th1_file->GetBinContent(1));
  signal_twomuzh_th1.SetBinContent(2, signal_twomuzh_th1_file->GetBinContent(2));
  signal_twomuzh_th1.SetBinContent(3, signal_twomuzh_th1_file->Integral(3,6));//assumes bin 6 is overflow
  RooDataHist signal_twomuzh_hist("signal_twomuzh", "Signal yield in TwoMuZH", vars, &signal_twomuzh_th1);
  wspace->import(signal_twomuzh_hist, RooFit::RecycleConflictNodes());

}


//---------------------------------------------------------------------------------------------------------------
//EleMu Region
//---------------------------------------------------------------------------------------------------------------
void build_elemu(RooWorkspace *wspace){

  RooRealVar* ntags = wspace->var("ntags");
  RooArgList vars(*ntags);

  //For now assume this region is 100% heavy background (TODO with RooDataHist like signal)

  //Data
  TFile* f_elemu = TFile::Open("../inputs/EleMuOSOF_nSelectedAODCaloJetTag_GH.root", "READ");
  TH1F* data_elemu_th1_file = (TH1F*)f_elemu->Get(data_string);
  TH1F data_elemu_th1("data_obs_elemu","Data observed in EleMu", 3, -0.5, 2.5);
  data_elemu_th1.SetBinContent(1, data_elemu_th1_file->GetBinContent(1));
  data_elemu_th1.SetBinContent(2, data_elemu_th1_file->GetBinContent(2));
  data_elemu_th1.SetBinContent(3, data_elemu_th1_file->Integral(3,6));//assumes bin 6 is overflow TODO
  RooDataHist data_elemu_hist("data_obs_elemu", "Data observed in EleMu", vars, &data_elemu_th1);
  wspace->import(data_elemu_hist, RooFit::RecycleConflictNodes());  

  
  //Heavy background is freely floating
  //Create one parameter per bin representing the yield 
  RooRealVar heavy_elemu_bin1("heavy_elemu_bin1", "Heavy background yield in EleMu, bin 1", 5126.56, 1e3, 5e4);
  RooRealVar heavy_elemu_bin2("heavy_elemu_bin2", "Heavy background yield in EleMu, bin 2", 98.25, 50, 200);
  RooRealVar heavy_elemu_bin3("heavy_elemu_bin3", "Heavy background yield in EleMu, bin 3", 1.29, 0, 5);
  RooArgList heavy_elemu_bins;
  heavy_elemu_bins.add(heavy_elemu_bin1);
  heavy_elemu_bins.add(heavy_elemu_bin2);
  heavy_elemu_bins.add(heavy_elemu_bin3);
  RooParametricHist p_heavy_elemu("heavy_elemu", "Heavy PDF in EleMu Region", *ntags, heavy_elemu_bins, data_elemu_th1);
  RooAddition p_heavy_elemu_norm("heavy_elemu_norm", "Total number of heavy events in EleMu Region", heavy_elemu_bins);
  wspace->import(p_heavy_elemu, RooFit::RecycleConflictNodes());
  wspace->import(p_heavy_elemu_norm, RooFit::RecycleConflictNodes());
  
}


//---------------------------------------------------------------------------------------------------------------
//OnePho Region
//---------------------------------------------------------------------------------------------------------------
void build_onepho(RooWorkspace* wspace){
  
  RooRealVar* ntags = wspace->var("ntags");
  RooArgList vars(*ntags);


  //Data
  TFile* f_onepho = TFile::Open("../inputs/OnePho_nSelectedAODCaloJetTag_GH.root", "READ");
  TH1F* data_onepho_th1_file = (TH1F*)f_onepho->Get(data_string);
  TH1F data_onepho_th1("data_obs_onepho","Data observed in OnePho", 3, -0.5, 2.5);
  data_onepho_th1.SetBinContent(1, data_onepho_th1_file->GetBinContent(1));
  data_onepho_th1.SetBinContent(2, data_onepho_th1_file->GetBinContent(2));
  data_onepho_th1.SetBinContent(3, data_onepho_th1_file->Integral(3,6));//assumes bin 6 is overflow
  RooDataHist data_onepho_hist("data_obs_onepho", "Data observed in OnePho", vars, &data_onepho_th1);
  wspace->import(data_onepho_hist, RooFit::RecycleConflictNodes());


  //Light background is freely floating
  //Create one parameter per bin representing the yield 
  RooRealVar light_onepho_bin1("light_onepho_bin1", "Light background yield in OnePho, bin 1", 288366, 200000, 400000);
  RooRealVar light_onepho_bin2("light_onepho_bin2", "Light background yield in OnePho, bin 2", 1377.44, 500, 2000);
  RooRealVar light_onepho_bin3("light_onepho_bin3", "Light background yield in OnePho, bin 3", 0.5, 0, 5);
  RooArgList light_onepho_bins;
  light_onepho_bins.add(light_onepho_bin1);
  light_onepho_bins.add(light_onepho_bin2);
  light_onepho_bins.add(light_onepho_bin3);
  RooParametricHist p_light_onepho("light_onepho", "Light PDF in OnePho Region", *ntags, light_onepho_bins, data_onepho_th1);
  RooAddition p_light_onepho_norm("light_onepho_norm", "Total number of light events in OnePho Region", light_onepho_bins);
  wspace->import(p_light_onepho, RooFit::RecycleConflictNodes());
  wspace->import(p_light_onepho_norm, RooFit::RecycleConflictNodes());

  
  //Get heavy 
  RooRealVar* heavy_elemu_bin1 = wspace->var("heavy_elemu_bin1");
  RooRealVar* heavy_elemu_bin2 = wspace->var("heavy_elemu_bin2");
  RooRealVar* heavy_elemu_bin3 = wspace->var("heavy_elemu_bin3");

  RooRealVar rrv_heavy_elemu_to_onepho_bin1("rrv_heavy_elemu_to_onepho_bin1", "heavy EleMu to OnePho bin 1",1);
  RooRealVar rrv_heavy_elemu_to_onepho_bin2("rrv_heavy_elemu_to_onepho_bin2", "heavy EleMu to OnePho bin 2",1);
  RooRealVar rrv_heavy_elemu_to_onepho_bin3("rrv_heavy_elemu_to_onepho_bin3", "heavy EleMu to OnePho bin 3",1);
 
  RooFormulaVar tf_heavy_elemu_to_onepho_bin1("tf_heavy_elemu_to_onepho_bin1", "heavy EleMu to OnePho transfer factor bin 1", 
					 "0.138834*TMath::Power(1+0.5,@0)", RooArgList(rrv_heavy_elemu_to_onepho_bin1) );
  RooFormulaVar tf_heavy_elemu_to_onepho_bin2("tf_heavy_elemu_to_onepho_bin2", "heavy EleMu to OnePho transfer factor bin 2", 
					 "0.138834*TMath::Power(1+0.5,@0)", RooArgList(rrv_heavy_elemu_to_onepho_bin2) );
  RooFormulaVar tf_heavy_elemu_to_onepho_bin3("tf_heavy_elemu_to_onepho_bin3", "heavy EleMu to OnePho transfer factor bin 3", 
					 "0.138834*TMath::Power(1+0.5,@0)", RooArgList(rrv_heavy_elemu_to_onepho_bin3) );

  RooFormulaVar heavy_onepho_bin1("heavy_onepho_bin1", "Heavy background yield in OnePho, bin 1", 
				  "@0*@1", RooArgList(tf_heavy_elemu_to_onepho_bin1, *heavy_elemu_bin1));
  RooFormulaVar heavy_onepho_bin2("heavy_onepho_bin2", "Heavy background yield in OnePho, bin 2", 
				  "@0*@1", RooArgList(tf_heavy_elemu_to_onepho_bin2, *heavy_elemu_bin2));
  RooFormulaVar heavy_onepho_bin3("heavy_onepho_bin3", "Heavy background yield in OnePho, bin 3", 
				  "@0*@1", RooArgList(tf_heavy_elemu_to_onepho_bin3, *heavy_elemu_bin3));
  RooArgList heavy_onepho_bins;
  heavy_onepho_bins.add(heavy_onepho_bin1);
  heavy_onepho_bins.add(heavy_onepho_bin2);
  heavy_onepho_bins.add(heavy_onepho_bin3);
  RooParametricHist p_heavy_onepho("heavy_onepho", "Heavy PDF in OnePho Region", *ntags, heavy_onepho_bins, data_onepho_th1);
  RooAddition p_heavy_onepho_norm("heavy_onepho_norm", "Total number of heavy events in OnePho Region", heavy_onepho_bins);

  wspace->import(p_heavy_onepho, RooFit::RecycleConflictNodes());
  wspace->import(p_heavy_onepho_norm, RooFit::RecycleConflictNodes());

}



//---------------------------------------------------------------------------------------------------------------
//TwoMuDY Region
//---------------------------------------------------------------------------------------------------------------
void build_twomudy(RooWorkspace* wspace){
  
  RooRealVar* ntags = wspace->var("ntags");
  RooArgList vars(*ntags);


  //Data
  TFile* f_twomudy = TFile::Open("../inputs/TwoMuDY_nSelectedAODCaloJetTag_GH.root", "READ");
  TH1F* data_twomudy_th1_file = (TH1F*)f_twomudy->Get(data_string);
  TH1F data_twomudy_th1("data_obs_twomudy","Data observed in TwoMuDY", 3, -0.5, 2.5);
  data_twomudy_th1.SetBinContent(1, data_twomudy_th1_file->GetBinContent(1));
  data_twomudy_th1.SetBinContent(2, data_twomudy_th1_file->GetBinContent(2));
  data_twomudy_th1.SetBinContent(3, data_twomudy_th1_file->Integral(3,6));//assumes bin 6 is overflow
  RooDataHist data_twomudy_hist("data_obs_twomudy", "Data observed in TwoMuDY", vars, &data_twomudy_th1);
  wspace->import(data_twomudy_hist, RooFit::RecycleConflictNodes());


  //Light background is freely floating
  //Create one parameter per bin representing the yield 
  RooRealVar light_twomudy_bin1("light_twomudy_bin1", "Light background yield in TwoMuDY, bin 1", 288366, 200000, 400000);
  RooRealVar light_twomudy_bin2("light_twomudy_bin2", "Light background yield in TwoMuDY, bin 2", 1377.44, 500, 2000);
  RooRealVar light_twomudy_bin3("light_twomudy_bin3", "Light background yield in TwoMuDY, bin 3", 0.5, 0, 5);
  RooArgList light_twomudy_bins;
  light_twomudy_bins.add(light_twomudy_bin1);
  light_twomudy_bins.add(light_twomudy_bin2);
  light_twomudy_bins.add(light_twomudy_bin3);
  RooParametricHist p_light_twomudy("light_twomudy", "Light PDF in TwoMuDY Region", *ntags, light_twomudy_bins, data_twomudy_th1);
  RooAddition p_light_twomudy_norm("light_twomudy_norm", "Total number of light events in TwoMuDY Region", light_twomudy_bins);
  wspace->import(p_light_twomudy, RooFit::RecycleConflictNodes());
  wspace->import(p_light_twomudy_norm, RooFit::RecycleConflictNodes());

  
  //Get heavy 
  RooRealVar* heavy_elemu_bin1 = wspace->var("heavy_elemu_bin1");
  RooRealVar* heavy_elemu_bin2 = wspace->var("heavy_elemu_bin2");
  RooRealVar* heavy_elemu_bin3 = wspace->var("heavy_elemu_bin3");


  RooRealVar rrv_heavy_elemu_to_twomudy_bin1("rrv_heavy_elemu_to_twomudy_bin1", "heavy EleMu to TwoMuDY bin 1",1);
  RooRealVar rrv_heavy_elemu_to_twomudy_bin2("rrv_heavy_elemu_to_twomudy_bin2", "heavy EleMu to TwoMuDY bin 2",1);
  RooRealVar rrv_heavy_elemu_to_twomudy_bin3("rrv_heavy_elemu_to_twomudy_bin3", "heavy EleMu to TwoMuDY bin 3",1);

  RooFormulaVar tf_heavy_elemu_to_twomudy_bin1("tf_heavy_elemu_to_twomudy_bin1", "heavy EleMu to TwoMuDY transfer factor bin 1", 
					  "0.138834*TMath::Power(1+0.5,@0)", RooArgList(rrv_heavy_elemu_to_twomudy_bin1) );
  RooFormulaVar tf_heavy_elemu_to_twomudy_bin2("tf_heavy_elemu_to_twomudy_bin2", "heavy EleMu to TwoMuDY transfer factor bin 2", 
					  "0.138834*TMath::Power(1+0.5,@0)", RooArgList(rrv_heavy_elemu_to_twomudy_bin2) );
  RooFormulaVar tf_heavy_elemu_to_twomudy_bin3("tf_heavy_elemu_to_twomudy_bin3", "heavy EleMu to TwoMuDY transfer factor bin 3", 
					  "0.138834*TMath::Power(1+0.5,@0)", RooArgList(rrv_heavy_elemu_to_twomudy_bin3) );
  
  RooFormulaVar heavy_twomudy_bin1("heavy_twomudy_bin1", "Heavy background yield in TwoMuDY, bin 1", 
				   "@0*@1", RooArgList(tf_heavy_elemu_to_twomudy_bin1, *heavy_elemu_bin1));
  RooFormulaVar heavy_twomudy_bin2("heavy_twomudy_bin2", "Heavy background yield in TwoMuDY, bin 2", 
				   "@0*@1", RooArgList(tf_heavy_elemu_to_twomudy_bin2, *heavy_elemu_bin2));
  RooFormulaVar heavy_twomudy_bin3("heavy_twomudy_bin3", "Heavy background yield in TwoMuDY, bin 3", 
				   "@0*@1", RooArgList(tf_heavy_elemu_to_twomudy_bin3, *heavy_elemu_bin3));
  RooArgList heavy_twomudy_bins;
  heavy_twomudy_bins.add(heavy_twomudy_bin1);
  heavy_twomudy_bins.add(heavy_twomudy_bin2);
  heavy_twomudy_bins.add(heavy_twomudy_bin3);
  RooParametricHist p_heavy_twomudy("heavy_twomudy", "Heavy PDF in TwoMuDY Region", *ntags, heavy_twomudy_bins, data_twomudy_th1);
  RooAddition p_heavy_twomudy_norm("heavy_twomudy_norm", "Total number of heavy events in TwoMuDY Region", heavy_twomudy_bins);

  wspace->import(p_heavy_twomudy, RooFit::RecycleConflictNodes());
  wspace->import(p_heavy_twomudy_norm, RooFit::RecycleConflictNodes());

}


void build_ws(){

  // As usual, load the combine library to get access to the RooParametricHist
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");

  // Output file and workspace
  TFile *fOut = new TFile("param_ws.root","RECREATE");
  RooWorkspace* wspace = new RooWorkspace("wspace","wspace");

  //Search in ntags, define ntags as our variable
  RooRealVar ntags("ntags", "ntags", -.5, 2.5);
  wspace->import(ntags, RooFit::RecycleConflictNodes());

  //TwoMuZH+EleMu+OnePho top-nontop
  //build_elemu(wspace);
  //build_onepho(wspace);
  //build_twomuzh(wspace);
  
  //TwoMuZH+EleMu+TwoMuDY top-nontop
  build_elemu(wspace);
  build_twomudy(wspace);
  build_twomuzh(wspace, "DY");

  //TwoMuZH+EleMu+TwoMuZH flavor split
  
  wspace->Print("v");

  fOut->cd();
  wspace->Write();

}
