#include <iostream>
#include <math.h>

#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPad.h"

#include "RooArgList.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAddition.h"
#include "RooDataHist.h"

#include "HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"


#include <TSystem.h>
#include <TTree.h>
#include <TLatex.h>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>
#include <Math/GaussIntegrator.h>
#include <Math/IntegratorOptions.h>
//ROOFIT INCLUDES
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooMinimizer.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooExtendPdf.h>
#include <RooStats/SPlot.h>
#include <RooStats/ModelConfig.h>
#include <RooGenericPdf.h>
#include <RooFormulaVar.h>
#include <RooBernstein.h>
#include <RooMinuit.h>
#include <RooNLLVar.h>
#include <RooRandom.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooArgList.h>


//LOCAL
#include "CommandLineInput.hh"

using namespace std;


//Todo:
// - make constraint extraction integrate 2 and up
// - add sig eff and other systematics


//Global options
bool _additional_method_sys = true;
TString dummy_syst = "0.01";
TString xvar = "nSelectedAODCaloJetTag";
//TString xvar = "nSelectedAODCaloJetTagSB7";
//TString xvar_suffix = "_GH";
//TString xvar_suffix = "_log_eemumu_zpt65";
//TString xvar_suffix = "_log_eemumu";
TString xvar_suffix = "_log";
TString signal_string = "Sig_MS55ct100";
//TString data_string = "bkgtotal"; //"Data";
TString data_string = "Data";
vector<TString> sys_vec;


void plot_h(TH1F* h, TString name){
  TCanvas c(name, name, 640, 480);
  h->SetTitle(name);
  h->Draw("HIST E");
  c.SaveAs("../analysis/plots/"+name+".pdf");
}


void plot_syst(TH1F* h, TH1F* hu, TH1F* hd, TString name, bool doLog){

  double max = h->GetMaximum();
  if(hu->GetMaximum() > max) max = hu->GetMaximum();
  if(hd->GetMaximum() > max) max = hd->GetMaximum();
  h->SetMaximum(max*1.25);

  if(doLog){
    h->SetMinimum(1e-3);
  }

  h->SetLineWidth(3);
  hu->SetLineColor(kGreen+1);
  hd->SetLineColor(kRed+1);
  TCanvas c(name, name, 640, 480);
  gPad->SetLogy(doLog);
  h->SetTitle(name);
  h->Draw("HIST E");
  hu->Draw("HIST E SAME");
  hd->Draw("HIST E SAME");
  c.SaveAs("../analysis/plots/"+name+".pdf");

}


TString translate(TString in){
  TString out = "";
  //if(in=="elemu") out = "EleMuOSOF";
  if(in=="elemu") out = "EleMu_Combined";
  //else if(in == "elemul") out = "EleMuOSOFL";//combine emu region
  else if(in == "twomuzh") out = "TwoMuZH";
  else if(in == "twomudy") out = "TwoMuDY";
  else if(in == "twoelezh") out = "TwoEleZH";
  else if(in == "twoeledy") out = "TwoEleDY";
  else{
    cout << "Invalid translate!" << endl;
    out = "";
  }
  return out;
}


void build_tf(RooWorkspace* wspace, TString process, TString from_name, TString to_name, vector<TString> sys_vec){

  bool missing3 = false;

  TString full_name = process+"_"+from_name+"_to_"+to_name;

  TFile* f_from = TFile::Open("../inputs/"+translate(from_name)+"_"+xvar+xvar_suffix+".root", "READ");
  TH1F*  h_from = (TH1F*)f_from->Get(process);

  TFile* f_to = TFile::Open("../inputs/"+translate(to_name)+"_"+xvar+xvar_suffix+".root", "READ");
  TH1F*  h_to = (TH1F*)f_to->Get(process);

  TH1F* h_r = (TH1F*)h_to->Clone("h_r_"+full_name);
  h_r->Divide(h_from);

  plot_h(h_r,"h_r_"+full_name);

  TString s_bin1="", s_bin2="", s_bin3="";
  TString s_bin1_err="", s_bin2_err="", s_bin3_err="";
  s_bin1     += h_r->GetBinContent(1);
  s_bin1_err += h_r->GetBinError(1)/h_r->GetBinContent(1);//relative error
  s_bin2     += h_r->GetBinContent(2);
  s_bin2_err += h_r->GetBinError(2)/h_r->GetBinContent(2);//relative error
  if(h_r->GetBinContent(3)>0){
    s_bin3     += h_r->GetBinContent(3);
    s_bin3_err += h_r->GetBinError(3)/h_r->GetBinContent(3);//relative error
    if(s_bin3_err.Atof() > 1.0)
      {
        s_bin3 = s_bin2;
        s_bin3_err = "0.99";
      }
  }else{
    missing3 = true;
    cout << endl; cout << "*** WARNING *** : Using 1-tag ratio for " << process << " " << from_name << " to " << to_name << endl; cout << endl;
    s_bin3=s_bin2;
    s_bin3_err=s_bin2_err;//relative error
    //s_bin3="0";
    //s_bin3_err="0.5";//relative error
  }

  //TEST TF DUMMY UNCERTAINTY
  if( from_name=="twoeledy" )
  {
    //s_bin1_err = s_bin2_err = s_bin3_err = "0.01";
    //s_bin1 += "*1.06371";//2016 correction on pt
    //s_bin2 += "*1.06371";//2016 correction on pt
    //s_bin3 = "0.68945861";//add-hoc, mu and ele combination
    //s_bin3_err = "0.83582353";
    //s_bin3 += "*1.06371";//2016 correction on pt
  }
  else if( from_name=="twomudy" )
  {
    //s_bin1_err = s_bin2_err = s_bin3_err = "0.01";
    //s_bin1 += "*0.995138";//2016 correction on pt
    //s_bin2 += "*0.995138";//2016 correction on pt
    //s_bin3 = "0.68945861";//add-hoc, mu and ele combination
    //s_bin3_err = "0.83582353";
    //s_bin3 += "*0.995138";//2016 correction on pt
  }

  RooRealVar rrv_bin1("rrv_"+full_name+"_bin1", process+" "+from_name+" to "+to_name+" bin 1",0);
  RooRealVar rrv_bin2("rrv_"+full_name+"_bin2", process+" "+from_name+" to "+to_name+" bin 2",0);
  RooRealVar rrv_bin3("rrv_"+full_name+"_bin3", process+" "+from_name+" to "+to_name+" bin 3",0);

  RooArgList ral_bin1 = RooArgList(rrv_bin1);
  RooArgList ral_bin2 = RooArgList(rrv_bin2);
  RooArgList ral_bin3 = RooArgList(rrv_bin3);

  TString rfv_bin1 = s_bin1+"*TMath::Power("+s_bin1_err+"+1,@0)";
  TString rfv_bin2 = s_bin2+"*TMath::Power("+s_bin2_err+"+1,@0)";
  TString rfv_bin3 = s_bin3+"*TMath::Power("+s_bin3_err+"+1,@0)";

  int sys_cnt = 1;
  for(unsigned int i=0; i<sys_vec.size(); i++){

    if(from_name=="twoeledy" && sys_vec[i]=="MES"){
      cout << "SKIP SYS " << full_name << " " << sys_vec[i] << endl;
      continue;
    }
    if(from_name=="twomudy" && sys_vec[i]=="EGS"){
      cout << "SKIP SYS " << full_name << " " << sys_vec[i] << endl;
      continue;
    }

    TFile* f_from_up = TFile::Open("../inputs/"+translate(from_name)+"_"+xvar+xvar_suffix+"_"+sys_vec[i]+"Up.root", "READ");
    TH1F* h_from_up = (TH1F*)f_from_up->Get(process);

    TFile* f_to_up = TFile::Open("../inputs/"+translate(to_name)+"_"+xvar+xvar_suffix+"_"+sys_vec[i]+"Up.root", "READ");
    TH1F* h_to_up = (TH1F*)f_to_up->Get(process);

    TH1F* h_r_up = (TH1F*)h_to_up->Clone("h_r_"+full_name+"_"+sys_vec[i]+"_up");
    h_r_up->Divide(h_from_up);

    TFile* f_from_down = TFile::Open("../inputs/"+translate(from_name)+"_"+xvar+xvar_suffix+"_"+sys_vec[i]+"Down.root", "READ");
    TH1F* h_from_down = (TH1F*)f_from_down->Get(process);

    TFile* f_to_down = TFile::Open("../inputs/"+translate(to_name)+"_"+xvar+xvar_suffix+"_"+sys_vec[i]+"Down.root", "READ");
    TH1F* h_to_down = (TH1F*)f_to_down->Get(process);

    TH1F* h_r_down = (TH1F*)h_to_down->Clone("h_r_"+full_name+"_"+sys_vec[i]+"_down");
    h_r_down->Divide(h_from_down);

    plot_syst(h_from, h_from_up, h_from_down, "h_ntag_"+process+"_"+from_name+"_"+sys_vec[i], true);
    plot_syst(h_to, h_to_up, h_to_down, "h_ntag_"+process+"_"+to_name+"_"+sys_vec[i], true);
    plot_syst(h_r, h_r_up, h_r_down, "h_r_"+full_name+"_"+sys_vec[i], false);

    //Symmetrize -- (up-down)/2
    //TH1F* h_r_symm = (TH1F*)h_r_up->Clone("h_r_symm_"+full_name+"_"+sys_vec[i]);
    //h_r_symm->Add(h_r_down, -1);
    //h_r_symm->Scale(0.5);

    TH1F* h_r_symm = (TH1F*)h_r_up->Clone("h_r_symm_"+full_name+"_"+sys_vec[i]);
    bool sys_samedir[6]={0,0,0,0,0,0};
    for(int j=1; j<=h_r_up->GetNbinsX(); j++){

      double nom  = h_r->GetBinContent(j);
      double up   = h_r_up->GetBinContent(j);
      double down = h_r_down->GetBinContent(j);

      //Skip error because not used.
      //This is an error on an error.
      double sym = 0;
      double sym_err = 0;

      //both higher than nom
      if(up>nom && down>nom){
	sys_samedir[j-1]=1;
	if(up>down) sym = up-nom;//positive
	else sym = down-nom;//positive
      }
      //both lower than nom
      else if(up<nom && down<nom){
	sys_samedir[j-1]=1;
	if(up<down) sym = up-nom;//negative
	else sym = down-nom;//negative
      }
      //opposite directions w.r.t. nom
      else{
	sys_samedir[j-1]=0;
	sym = 0.5 * (up-down);//positive if up>nom and down, else negative
      }

      cout<< "RAW SYS " << full_name << sys_vec[i] << " " << sym << endl;
      cout<< "   up " << up << " nom " << nom << " down " << down << endl;

      h_r_symm->SetBinContent(j,sym);
      h_r_symm->SetBinError(j,sym_err);
    }

    //Make relative
    TH1F* h_r_symm_rel = (TH1F*)h_r_symm->Clone("h_r_symm_rel_"+full_name+"_"+sys_vec[i]);
    h_r_symm_rel->Divide(h_r);

    //Now the RooFit part
    //Might be cleaner to switch this to a loop over bins!
    TString s_bin1_syst = "", s_bin2_syst = "", s_bin3_syst = "";
    TString s_bin1_sign = "", s_bin2_sign = "", s_bin3_sign = "";
    TString s_bin1_abs = "", s_bin1_abs_end = "";
    TString s_bin2_abs = "", s_bin2_abs_end = "";
    TString s_bin3_abs = "", s_bin3_abs_end = "";

    if(sys_samedir[0]==1){
      s_bin1_abs = "TMath::Abs(";
      s_bin1_abs_end = ")";
    }
    if(sys_samedir[1]==1){
      s_bin2_abs = "TMath::Abs(";
      s_bin2_abs_end = ")";
    }
    if(sys_samedir[2]==1){
      s_bin3_abs = "TMath::Abs(";
      s_bin3_abs_end = ")";
    }

    if(h_r_symm_rel->GetBinContent(1)<0) s_bin1_sign="-";
    if(h_r_symm_rel->GetBinContent(2)<0) s_bin2_sign="-";

    s_bin1_syst += h_r_symm_rel->GetBinContent(1);
    s_bin2_syst += h_r_symm_rel->GetBinContent(2);
    //TEMP
    if(missing3 && h_r_symm_rel->GetBinContent(2)>0){
      s_bin3_syst = dummy_syst;
    }
    else if(missing3 && h_r_symm_rel->GetBinContent(2)<0){
      s_bin3_syst = "-"; s_bin3_syst = dummy_syst;
      s_bin3_sign="-";
    }
    else{
      s_bin3_syst += h_r_symm_rel->GetBinContent(3);
      if(h_r_symm_rel->GetBinContent(3)<0) s_bin3_sign="-";
    }

    rfv_bin1 += "*TMath::Power( TMath::Abs("; rfv_bin1+=s_bin1_syst; rfv_bin1+=")+1,"; rfv_bin1+=s_bin1_sign; rfv_bin1+=s_bin1_abs; rfv_bin1+="@"; rfv_bin1+=sys_cnt; rfv_bin1+=s_bin1_abs_end; rfv_bin1+=")";
    rfv_bin2 += "*TMath::Power( TMath::Abs("; rfv_bin2+=s_bin2_syst; rfv_bin2+=")+1,"; rfv_bin2+=s_bin2_sign; rfv_bin2+=s_bin2_abs; rfv_bin2+="@"; rfv_bin2+=sys_cnt; rfv_bin2+=s_bin2_abs_end; rfv_bin2+=")";
    rfv_bin3 += "*TMath::Power( TMath::Abs("; rfv_bin3+=s_bin3_syst; rfv_bin3+=")+1,"; rfv_bin3+=s_bin3_sign; rfv_bin3+=s_bin3_abs; rfv_bin3+="@"; rfv_bin3+=sys_cnt; rfv_bin3+=s_bin3_abs_end; rfv_bin3+=")";

    RooRealVar* rrv_syst = wspace->var("rrv_"+sys_vec[i]);
    ral_bin1.add(RooArgList(*rrv_syst));
    ral_bin2.add(RooArgList(*rrv_syst));
    ral_bin3.add(RooArgList(*rrv_syst));

    sys_cnt++;
  }//loop over systematics

  cout << "RFV1: " << rfv_bin1 << endl;
  ral_bin1.Print();
  cout << "RFV2: " << rfv_bin2 << endl;
  ral_bin2.Print();
  cout << "RFV3: " << rfv_bin3 << endl;
  ral_bin3.Print();

  RooFormulaVar tf_bin1("tf_"+full_name+"_bin1", process+" "+from_name+" to "+to_name+" transfer factor bin 1", rfv_bin1, ral_bin1);
  RooFormulaVar tf_bin2("tf_"+full_name+"_bin2", process+" "+from_name+" to "+to_name+" transfer factor bin 2", rfv_bin2, ral_bin2);
  RooFormulaVar tf_bin3("tf_"+full_name+"_bin3", process+" "+from_name+" to "+to_name+" transfer factor bin 3", rfv_bin3, ral_bin3);

  wspace->import(tf_bin1, RooFit::RecycleConflictNodes());
  wspace->import(tf_bin2, RooFit::RecycleConflictNodes());
  wspace->import(tf_bin3, RooFit::RecycleConflictNodes());

}


//---------------------------------------------------------------------------------------------------------------
//TwoMuZH Region
//---------------------------------------------------------------------------------------------------------------
void build_twomuzh(RooWorkspace* wspace, TString light_est = "DY"){

  RooRealVar* ntags = wspace->var("ntags");
  RooArgList vars(*ntags);


  //Data
  TFile* f_twomuzh = TFile::Open("../inputs/TwoMuZH_"+xvar+xvar_suffix+".root", "READ");
  TH1F* data_twomuzh_th1_file = (TH1F*)f_twomuzh->Get(data_string);
  TH1F data_twomuzh_th1("data_obs_twomuzh","Data observed in TwoMuZH", 3, -0.5, 2.5);
  data_twomuzh_th1.SetBinContent(1, data_twomuzh_th1_file->GetBinContent(1));
  data_twomuzh_th1.SetBinContent(2, data_twomuzh_th1_file->GetBinContent(2));
  data_twomuzh_th1.SetBinContent(3, data_twomuzh_th1_file->Integral(3,6));//assumes bin 6 is overflow
  RooDataHist data_twomuzh_hist("data_obs_twomuzh", "Data observed in TwoMuZH", vars, &data_twomuzh_th1);
  wspace->import(data_twomuzh_hist);

  //----------------------------------------------------
  //setting up to add additional systematic uncertainty
  //----------------------------------------------------
  //std::string add_logN_systematic = "TMath::Power(1+0.1,@0)";
  std::string add_logN_systematic[] = {"TMath::Power(1+0.15,@0)","TMath::Power(1+0.15,@0)","TMath::Power(1+0.15,@0)"};
  RooRealVar rrv_twomuzh_add_sys_bin1("rrv_twomuzh_add_sys_bin1","rrv_twomuzh_add_sys_bin1", 0.0);
  RooRealVar rrv_twomuzh_add_sys_bin2("rrv_twomuzh_add_sys_bin2","rrv_twomuzh_add_sys_bin2", 0.0);
  RooRealVar rrv_twomuzh_add_sys_bin3("rrv_twomuzh_add_sys_bin3","rrv_twomuzh_add_sys_bin3", 0.0);

  //Other contamination
  TH1F* other_twomuzh_th1_file = (TH1F*)f_twomuzh->Get("other");
  TH1F other_twomuzh_th1("other_twomuzh","Other yield in TwoMuZH", 3, -0.5, 2.5);
  other_twomuzh_th1.SetBinContent(1, other_twomuzh_th1_file->GetBinContent(1));
  other_twomuzh_th1.SetBinContent(2, other_twomuzh_th1_file->GetBinContent(2));
  other_twomuzh_th1.SetBinContent(3, other_twomuzh_th1_file->Integral(3,6));//assumes bin 6 is overflow TODO
  RooDataHist other_twomuzh_hist("other_twomuzh", "Other yield in TwoMuZH", vars, &other_twomuzh_th1);
  wspace->import(other_twomuzh_hist);


  //Heavy as function of heavy_elemu
  build_tf(wspace, "heavy", "elemu", "twomuzh", sys_vec);
  RooFormulaVar* tf_heavy_elemu_to_twomuzh_bin1 = (RooFormulaVar*)wspace->arg("tf_heavy_elemu_to_twomuzh_bin1");
  RooFormulaVar* tf_heavy_elemu_to_twomuzh_bin2 = (RooFormulaVar*)wspace->arg("tf_heavy_elemu_to_twomuzh_bin2");
  RooFormulaVar* tf_heavy_elemu_to_twomuzh_bin3 = (RooFormulaVar*)wspace->arg("tf_heavy_elemu_to_twomuzh_bin3");

  RooRealVar* heavy_elemu_bin1 = wspace->var("heavy_elemu_bin1");
  RooRealVar* heavy_elemu_bin2 = wspace->var("heavy_elemu_bin2");
  RooRealVar* heavy_elemu_bin3 = wspace->var("heavy_elemu_bin3");

  //---------------------------------------------------------
  //Include additional systematic uncertainty for the method
  //---------------------------------------------------------
  RooFormulaVar heavy_twomuzh_logN_add_sys_bin1("heavy_twomuzh_logN_add_sys_bin1", "additional log-normal for method, bin 1",
  add_logN_systematic[0].c_str(), RooArgList(rrv_twomuzh_add_sys_bin1));
  RooFormulaVar heavy_twomuzh_logN_add_sys_bin2("heavy_twomuzh_logN_add_sys_bin2", "additional log-normal for method, bin 2",
  add_logN_systematic[1].c_str(), RooArgList(rrv_twomuzh_add_sys_bin2));
  RooFormulaVar heavy_twomuzh_logN_add_sys_bin3("heavy_twomuzh_logN_add_sys_bin3", "additional log-normal for method, bin 3",
  add_logN_systematic[2].c_str(), RooArgList(rrv_twomuzh_add_sys_bin3));


  RooFormulaVar* heavy_twomuzh_bin1;
  RooFormulaVar* heavy_twomuzh_bin2;
  RooFormulaVar* heavy_twomuzh_bin3;
  if(_additional_method_sys)
  {
    heavy_twomuzh_bin1 = new RooFormulaVar("heavy_twomuzh_bin1", "Heavy background yield in twomuzh, bin 1",
             "@0*@1*@2", RooArgList(*tf_heavy_elemu_to_twomuzh_bin1, *heavy_elemu_bin1, heavy_twomuzh_logN_add_sys_bin1));
    heavy_twomuzh_bin2 = new RooFormulaVar("heavy_twomuzh_bin2", "Heavy background yield in twomuzh, bin 2",
             "@0*@1*@2", RooArgList(*tf_heavy_elemu_to_twomuzh_bin2, *heavy_elemu_bin2, heavy_twomuzh_logN_add_sys_bin2));
    heavy_twomuzh_bin3 = new RooFormulaVar("heavy_twomuzh_bin3", "Heavy background yield in twomuzh, bin 3",
             "@0*@1*@2", RooArgList(*tf_heavy_elemu_to_twomuzh_bin3, *heavy_elemu_bin3, heavy_twomuzh_logN_add_sys_bin3));
  }
  else
  {
    heavy_twomuzh_bin1 = new RooFormulaVar("heavy_twomuzh_bin1", "Heavy background yield in twomuzh, bin 1",
  				   "@0*@1", RooArgList(*tf_heavy_elemu_to_twomuzh_bin1, *heavy_elemu_bin1));
    heavy_twomuzh_bin2 = new RooFormulaVar("heavy_twomuzh_bin2", "Heavy background yield in twomuzh, bin 2",
  				   "@0*@1", RooArgList(*tf_heavy_elemu_to_twomuzh_bin2, *heavy_elemu_bin2));
    heavy_twomuzh_bin3 = new RooFormulaVar("heavy_twomuzh_bin3", "Heavy background yield in twomuzh, bin 3",
  				   "@0*@1", RooArgList(*tf_heavy_elemu_to_twomuzh_bin3, *heavy_elemu_bin3));
  }

  RooArgList heavy_twomuzh_bins;
  heavy_twomuzh_bins.add(*heavy_twomuzh_bin1);
  heavy_twomuzh_bins.add(*heavy_twomuzh_bin2);
  heavy_twomuzh_bins.add(*heavy_twomuzh_bin3);
  RooParametricHist p_heavy_twomuzh("heavy_twomuzh", "Heavy PDF in TwoMuZH Region", *ntags, heavy_twomuzh_bins, data_twomuzh_th1);
  RooAddition p_heavy_twomuzh_norm("heavy_twomuzh_norm", "Total number of heavy events in TwoMuZH Region", heavy_twomuzh_bins);

  wspace->import(p_heavy_twomuzh);
  wspace->import(p_heavy_twomuzh_norm, RooFit::RecycleConflictNodes());


  //Light as function of light_onepho
  if(light_est=="OnePho"){
    build_tf(wspace, "light", "onepho", "twomuzh", sys_vec);
    RooFormulaVar* tf_light_onepho_to_twomuzh_bin1 = (RooFormulaVar*)wspace->arg("tf_light_onepho_to_twomuzh_bin1");
    RooFormulaVar* tf_light_onepho_to_twomuzh_bin2 = (RooFormulaVar*)wspace->arg("tf_light_onepho_to_twomuzh_bin2");
    RooFormulaVar* tf_light_onepho_to_twomuzh_bin3 = (RooFormulaVar*)wspace->arg("tf_light_onepho_to_twomuzh_bin3");

    RooRealVar* light_onepho_bin1 = wspace->var("light_onepho_bin1");
    RooRealVar* light_onepho_bin2 = wspace->var("light_onepho_bin2");
    RooRealVar* light_onepho_bin3 = wspace->var("light_onepho_bin3");

    RooFormulaVar light_twomuzh_bin1("light_twomuzh_bin1", "Light background yield in TwoMuZH, bin 1",
				     "@0*@1", RooArgList(*tf_light_onepho_to_twomuzh_bin1, *light_onepho_bin1));
    RooFormulaVar light_twomuzh_bin2("light_twomuzh_bin2", "Light background yield in TwoMuZH, bin 2",
				     "@0*@1", RooArgList(*tf_light_onepho_to_twomuzh_bin2, *light_onepho_bin2));
    RooFormulaVar light_twomuzh_bin3("light_twomuzh_bin3", "Light background yield in TwoMuZH, bin 3",
				     "@0*@1", RooArgList(*tf_light_onepho_to_twomuzh_bin3, *light_onepho_bin3));
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
    build_tf(wspace, "light", "twomudy", "twomuzh", sys_vec);
    RooFormulaVar* tf_light_twomudy_to_twomuzh_bin1 = (RooFormulaVar*)wspace->arg("tf_light_twomudy_to_twomuzh_bin1");
    RooFormulaVar* tf_light_twomudy_to_twomuzh_bin2 = (RooFormulaVar*)wspace->arg("tf_light_twomudy_to_twomuzh_bin2");
    RooFormulaVar* tf_light_twomudy_to_twomuzh_bin3 = (RooFormulaVar*)wspace->arg("tf_light_twomudy_to_twomuzh_bin3");

    RooRealVar* light_twomudy_bin1 = wspace->var("light_twomudy_bin1");
    RooRealVar* light_twomudy_bin2 = wspace->var("light_twomudy_bin2");
    RooRealVar* light_twomudy_bin3 = wspace->var("light_twomudy_bin3");

    //---------------------------------------------------------
    //Include additional systematic uncertainty for the method
    //---------------------------------------------------------
    RooFormulaVar light_twomuzh_logN_add_sys_bin1("light_twomuzh_logN_add_sys_bin1", "additional log-normal for method, bin 1",
    add_logN_systematic[0].c_str(), RooArgList(rrv_twomuzh_add_sys_bin1));
    RooFormulaVar light_twomuzh_logN_add_sys_bin2("light_twomuzh_logN_add_sys_bin2", "additional log-normal for method, bin 2",
    add_logN_systematic[1].c_str(), RooArgList(rrv_twomuzh_add_sys_bin2));
    RooFormulaVar light_twomuzh_logN_add_sys_bin3("light_twomuzh_logN_add_sys_bin3", "additional log-normal for method, bin 3",
    add_logN_systematic[2].c_str(), RooArgList(rrv_twomuzh_add_sys_bin3));

    RooFormulaVar* light_twomuzh_bin1;
    RooFormulaVar* light_twomuzh_bin2;
    RooFormulaVar* light_twomuzh_bin3;
    if(_additional_method_sys)
    {
      light_twomuzh_bin1 = new RooFormulaVar("light_twomuzh_bin1", "Light background yield in twomuzh, bin 1",
               "@0*@1*@2", RooArgList(*tf_light_twomudy_to_twomuzh_bin1, *light_twomudy_bin1, light_twomuzh_logN_add_sys_bin1));
      light_twomuzh_bin2 = new RooFormulaVar("light_twomuzh_bin2", "Light background yield in twomuzh, bin 2",
               "@0*@1*@2", RooArgList(*tf_light_twomudy_to_twomuzh_bin2, *light_twomudy_bin2, light_twomuzh_logN_add_sys_bin2));
      light_twomuzh_bin3 = new RooFormulaVar("light_twomuzh_bin3", "Light background yield in twomuzh, bin 3",
               "@0*@1*@2", RooArgList(*tf_light_twomudy_to_twomuzh_bin3, *light_twomudy_bin3, light_twomuzh_logN_add_sys_bin3));
    }
    else
    {
      light_twomuzh_bin1 = new RooFormulaVar("light_twomuzh_bin1", "Light background yield in twomuzh, bin 1",
  				     "@0*@1", RooArgList(*tf_light_twomudy_to_twomuzh_bin1, *light_twomudy_bin1));
      light_twomuzh_bin2 = new RooFormulaVar("light_twomuzh_bin2", "Light background yield in twomuzh, bin 2",
  				     "@0*@1", RooArgList(*tf_light_twomudy_to_twomuzh_bin2, *light_twomudy_bin2));
      light_twomuzh_bin3 = new RooFormulaVar("light_twomuzh_bin3", "Light background yield in twomuzh, bin 3",
  				     "@0*@1", RooArgList(*tf_light_twomudy_to_twomuzh_bin3, *light_twomudy_bin3));
    }


    RooArgList light_twomuzh_bins;
    light_twomuzh_bins.add(*light_twomuzh_bin1);
    light_twomuzh_bins.add(*light_twomuzh_bin2);
    light_twomuzh_bins.add(*light_twomuzh_bin3);
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
  wspace->import(signal_twomuzh_hist);

}


//---------------------------------------------------------------------------------------------------------------
//TwoEleZH Region
//---------------------------------------------------------------------------------------------------------------
void build_twoelezh(RooWorkspace* wspace, TString light_est = "DY"){

  RooRealVar* ntags = wspace->var("ntags");
  RooArgList vars(*ntags);


  //Data
  TFile* f_twoelezh = TFile::Open("../inputs/TwoEleZH_"+xvar+xvar_suffix+".root", "READ");
  TH1F* data_twoelezh_th1_file = (TH1F*)f_twoelezh->Get(data_string);
  TH1F data_twoelezh_th1("data_obs_twoelezh","Data observed in TwoEleZH", 3, -0.5, 2.5);
  data_twoelezh_th1.SetBinContent(1, data_twoelezh_th1_file->GetBinContent(1));
  data_twoelezh_th1.SetBinContent(2, data_twoelezh_th1_file->GetBinContent(2));
  data_twoelezh_th1.SetBinContent(3, data_twoelezh_th1_file->Integral(3,6));//assumes bin 6 is overflow
  RooDataHist data_twoelezh_hist("data_obs_twoelezh", "Data observed in TwoEleZH", vars, &data_twoelezh_th1);
  wspace->import(data_twoelezh_hist);


  //----------------------------------------------------
  //setting up to add additional systematic uncertainty
  //----------------------------------------------------
  //std::string add_logN_systematic = "TMath::Power(1+0.1,@0)";
  std::string add_logN_systematic = "TMath::Power(1+0.3,@0)";
  RooRealVar rrv_twoelezh_add_sys_bin1("rrv_twoelezh_add_sys_bin1","rrv_twoelezh_add_sys_bin1", 1.0);
  RooRealVar rrv_twoelezh_add_sys_bin2("rrv_twoelezh_add_sys_bin2","rrv_twoelezh_add_sys_bin2", 1.0);
  RooRealVar rrv_twoelezh_add_sys_bin3("rrv_twoelezh_add_sys_bin3","rrv_twoelezh_add_sys_bin3", 1.0);

  //Other contamination
  TH1F* other_twoelezh_th1_file = (TH1F*)f_twoelezh->Get("other");
  TH1F other_twoelezh_th1("other_twoelezh","Other yield in TwoEleZH", 3, -0.5, 2.5);
  other_twoelezh_th1.SetBinContent(1, other_twoelezh_th1_file->GetBinContent(1));
  other_twoelezh_th1.SetBinContent(2, other_twoelezh_th1_file->GetBinContent(2));
  other_twoelezh_th1.SetBinContent(3, other_twoelezh_th1_file->Integral(3,6));//assumes bin 6 is overflow TODO
  RooDataHist other_twoelezh_hist("other_twoelezh", "Other yield in TwoEleZH", vars, &other_twoelezh_th1);
  wspace->import(other_twoelezh_hist);


  //Heavy as function of heavy_elemu
  build_tf(wspace, "heavy", "elemu", "twoelezh", sys_vec);
  RooFormulaVar* tf_heavy_elemu_to_twoelezh_bin1 = (RooFormulaVar*)wspace->arg("tf_heavy_elemu_to_twoelezh_bin1");
  RooFormulaVar* tf_heavy_elemu_to_twoelezh_bin2 = (RooFormulaVar*)wspace->arg("tf_heavy_elemu_to_twoelezh_bin2");
  RooFormulaVar* tf_heavy_elemu_to_twoelezh_bin3 = (RooFormulaVar*)wspace->arg("tf_heavy_elemu_to_twoelezh_bin3");

  RooRealVar* heavy_elemu_bin1 = wspace->var("heavy_elemu_bin1");
  RooRealVar* heavy_elemu_bin2 = wspace->var("heavy_elemu_bin2");
  RooRealVar* heavy_elemu_bin3 = wspace->var("heavy_elemu_bin3");

  //---------------------------------------------------------
  //Include additional systematic uncertainty for the method
  //---------------------------------------------------------
  RooFormulaVar heavy_twoelezh_logN_add_sys_bin1("heavy_twoelezh_logN_add_sys_bin1", "additional log-normal for method, bin 1",
  add_logN_systematic.c_str(), RooArgList(rrv_twoelezh_add_sys_bin1));
  RooFormulaVar heavy_twoelezh_logN_add_sys_bin2("heavy_twoelezh_logN_add_sys_bin2", "additional log-normal for method, bin 2",
  add_logN_systematic.c_str(), RooArgList(rrv_twoelezh_add_sys_bin2));
  RooFormulaVar heavy_twoelezh_logN_add_sys_bin3("heavy_twoelezh_logN_add_sys_bin3", "additional log-normal for method, bin 3",
  add_logN_systematic.c_str(), RooArgList(rrv_twoelezh_add_sys_bin3));


  RooFormulaVar* heavy_twoelezh_bin1;
  RooFormulaVar* heavy_twoelezh_bin2;
  RooFormulaVar* heavy_twoelezh_bin3;
  if(_additional_method_sys)
  {
    heavy_twoelezh_bin1 = new RooFormulaVar("heavy_twoelezh_bin1", "Heavy background yield in TwoEleZH, bin 1",
             "@0*@1*@2", RooArgList(*tf_heavy_elemu_to_twoelezh_bin1, *heavy_elemu_bin1, heavy_twoelezh_logN_add_sys_bin1));
    heavy_twoelezh_bin2 = new RooFormulaVar("heavy_twoelezh_bin2", "Heavy background yield in TwoEleZH, bin 2",
             "@0*@1*@2", RooArgList(*tf_heavy_elemu_to_twoelezh_bin2, *heavy_elemu_bin2, heavy_twoelezh_logN_add_sys_bin2));
    heavy_twoelezh_bin3 = new RooFormulaVar("heavy_twoelezh_bin3", "Heavy background yield in TwoEleZH, bin 3",
             "@0*@1*@2", RooArgList(*tf_heavy_elemu_to_twoelezh_bin3, *heavy_elemu_bin3, heavy_twoelezh_logN_add_sys_bin3));
  }
  else
  {
    heavy_twoelezh_bin1 = new RooFormulaVar("heavy_twoelezh_bin1", "Heavy background yield in TwoEleZH, bin 1",
  				   "@0*@1", RooArgList(*tf_heavy_elemu_to_twoelezh_bin1, *heavy_elemu_bin1));
    heavy_twoelezh_bin2 = new RooFormulaVar("heavy_twoelezh_bin2", "Heavy background yield in TwoEleZH, bin 2",
  				   "@0*@1", RooArgList(*tf_heavy_elemu_to_twoelezh_bin2, *heavy_elemu_bin2));
    heavy_twoelezh_bin3 = new RooFormulaVar("heavy_twoelezh_bin3", "Heavy background yield in TwoEleZH, bin 3",
  				   "@0*@1", RooArgList(*tf_heavy_elemu_to_twoelezh_bin3, *heavy_elemu_bin3));
  }
  //std::cout << "EleZH heavy bin1: " << heavy_twoelezh_bin1.Print() << std::endl;
  RooArgList heavy_twoelezh_bins;
  heavy_twoelezh_bins.add(*heavy_twoelezh_bin1);
  heavy_twoelezh_bins.add(*heavy_twoelezh_bin2);
  heavy_twoelezh_bins.add(*heavy_twoelezh_bin3);
  RooParametricHist p_heavy_twoelezh("heavy_twoelezh", "Heavy PDF in TwoEleZH Region", *ntags, heavy_twoelezh_bins, data_twoelezh_th1);
  RooAddition p_heavy_twoelezh_norm("heavy_twoelezh_norm", "Total number of heavy events in TwoEleZH Region", heavy_twoelezh_bins);

  wspace->import(p_heavy_twoelezh);
  wspace->import(p_heavy_twoelezh_norm, RooFit::RecycleConflictNodes());


  //Light as function of light_onepho
  if(light_est=="OnePho"){
    build_tf(wspace, "light", "onepho", "twoelezh", sys_vec);
    RooFormulaVar* tf_light_onepho_to_twoelezh_bin1 = (RooFormulaVar*)wspace->arg("tf_light_onepho_to_twoelezh_bin1");
    RooFormulaVar* tf_light_onepho_to_twoelezh_bin2 = (RooFormulaVar*)wspace->arg("tf_light_onepho_to_twoelezh_bin2");
    RooFormulaVar* tf_light_onepho_to_twoelezh_bin3 = (RooFormulaVar*)wspace->arg("tf_light_onepho_to_twoelezh_bin3");

    RooRealVar* light_onepho_bin1 = wspace->var("light_onepho_bin1");
    RooRealVar* light_onepho_bin2 = wspace->var("light_onepho_bin2");
    RooRealVar* light_onepho_bin3 = wspace->var("light_onepho_bin3");

    RooFormulaVar light_twoelezh_bin1("light_twoelezh_bin1", "Light background yield in TwoEleZH, bin 1",
				     "@0*@1", RooArgList(*tf_light_onepho_to_twoelezh_bin1, *light_onepho_bin1));
    RooFormulaVar light_twoelezh_bin2("light_twoelezh_bin2", "Light background yield in TwoEleZH, bin 2",
				     "@0*@1", RooArgList(*tf_light_onepho_to_twoelezh_bin2, *light_onepho_bin2));
    RooFormulaVar light_twoelezh_bin3("light_twoelezh_bin3", "Light background yield in TwoEleZH, bin 3",
				     "@0*@1", RooArgList(*tf_light_onepho_to_twoelezh_bin3, *light_onepho_bin3));
    RooArgList light_twoelezh_bins;
    light_twoelezh_bins.add(light_twoelezh_bin1);
    light_twoelezh_bins.add(light_twoelezh_bin2);
    light_twoelezh_bins.add(light_twoelezh_bin3);
    RooParametricHist p_light_twoelezh("light_twoelezh", "Light PDF in TwoEleZH Region", *ntags, light_twoelezh_bins, data_twoelezh_th1);
    RooAddition p_light_twoelezh_norm("light_twoelezh_norm", "Total number of light events in TwoEleZH Region", light_twoelezh_bins);

    wspace->import(p_light_twoelezh);
    wspace->import(p_light_twoelezh_norm, RooFit::RecycleConflictNodes());
  }
  else if(light_est == "DY"){
    build_tf(wspace, "light", "twoeledy", "twoelezh", sys_vec);
    RooFormulaVar* tf_light_twoeledy_to_twoelezh_bin1 = (RooFormulaVar*)wspace->arg("tf_light_twoeledy_to_twoelezh_bin1");
    RooFormulaVar* tf_light_twoeledy_to_twoelezh_bin2 = (RooFormulaVar*)wspace->arg("tf_light_twoeledy_to_twoelezh_bin2");
    RooFormulaVar* tf_light_twoeledy_to_twoelezh_bin3 = (RooFormulaVar*)wspace->arg("tf_light_twoeledy_to_twoelezh_bin3");

    RooRealVar* light_twoeledy_bin1 = wspace->var("light_twoeledy_bin1");
    RooRealVar* light_twoeledy_bin2 = wspace->var("light_twoeledy_bin2");
    RooRealVar* light_twoeledy_bin3 = wspace->var("light_twoeledy_bin3");

    //---------------------------------------------------------
    //Include additional systematic uncertainty for the method
    //---------------------------------------------------------
    RooFormulaVar light_twoelezh_logN_add_sys_bin1("light_twoelezh_logN_add_sys_bin1", "additional log-normal for method, bin 1",
    add_logN_systematic.c_str(), RooArgList(rrv_twoelezh_add_sys_bin1));
    RooFormulaVar light_twoelezh_logN_add_sys_bin2("light_twoelezh_logN_add_sys_bin2", "additional log-normal for method, bin 2",
    add_logN_systematic.c_str(), RooArgList(rrv_twoelezh_add_sys_bin2));
    RooFormulaVar light_twoelezh_logN_add_sys_bin3("light_twoelezh_logN_add_sys_bin3", "additional log-normal for method, bin 3",
    add_logN_systematic.c_str(), RooArgList(rrv_twoelezh_add_sys_bin3));


    RooFormulaVar* light_twoelezh_bin1;
    RooFormulaVar* light_twoelezh_bin2;
    RooFormulaVar* light_twoelezh_bin3;
    if(_additional_method_sys)
    {
      light_twoelezh_bin1 = new RooFormulaVar("light_twoelezh_bin1", "Light background yield in TwoEleZH, bin 1",
               "@0*@1*@2", RooArgList(*tf_light_twoeledy_to_twoelezh_bin1, *light_twoeledy_bin1, light_twoelezh_logN_add_sys_bin1));
      light_twoelezh_bin2 = new RooFormulaVar("light_twoelezh_bin2", "Light background yield in TwoEleZH, bin 2",
               "@0*@1*@2", RooArgList(*tf_light_twoeledy_to_twoelezh_bin2, *light_twoeledy_bin2, light_twoelezh_logN_add_sys_bin2));
      light_twoelezh_bin3 = new RooFormulaVar("light_twoelezh_bin3", "Light background yield in TwoEleZH, bin 3",
               "@0*@1*@2", RooArgList(*tf_light_twoeledy_to_twoelezh_bin3, *light_twoeledy_bin3, light_twoelezh_logN_add_sys_bin3));
    }
    else
    {
      light_twoelezh_bin1 = new RooFormulaVar("light_twoelezh_bin1", "Light background yield in TwoEleZH, bin 1",
  				     "@0*@1", RooArgList(*tf_light_twoeledy_to_twoelezh_bin1, *light_twoeledy_bin1));
      light_twoelezh_bin2 = new RooFormulaVar("light_twoelezh_bin2", "Light background yield in TwoEleZH, bin 2",
  				     "@0*@1", RooArgList(*tf_light_twoeledy_to_twoelezh_bin2, *light_twoeledy_bin2));
      light_twoelezh_bin3 = new RooFormulaVar("light_twoelezh_bin3", "Light background yield in TwoEleZH, bin 3",
  				     "@0*@1", RooArgList(*tf_light_twoeledy_to_twoelezh_bin3, *light_twoeledy_bin3));
    }

    RooArgList light_twoelezh_bins;
    light_twoelezh_bins.add(*light_twoelezh_bin1);
    light_twoelezh_bins.add(*light_twoelezh_bin2);
    light_twoelezh_bins.add(*light_twoelezh_bin3);
    RooParametricHist p_light_twoelezh("light_twoelezh", "Light PDF in TwoEleZH Region", *ntags, light_twoelezh_bins, data_twoelezh_th1);
    RooAddition p_light_twoelezh_norm("light_twoelezh_norm", "Total number of light events in TwoEleZH Region", light_twoelezh_bins);

    wspace->import(p_light_twoelezh);
    wspace->import(p_light_twoelezh_norm, RooFit::RecycleConflictNodes());
  }
  else {
    cout << "Invalid option for light estimate" << endl;
    return;
  }

  //Signal
  TH1F* signal_twoelezh_th1_file = (TH1F*)f_twoelezh->Get(signal_string);
  TH1F signal_twoelezh_th1("signal_twoelezh","Signal yield in TwoEleZH", 3, -0.5, 2.5);
  signal_twoelezh_th1.SetBinContent(1, signal_twoelezh_th1_file->GetBinContent(1));
  signal_twoelezh_th1.SetBinContent(2, signal_twoelezh_th1_file->GetBinContent(2));
  signal_twoelezh_th1.SetBinContent(3, signal_twoelezh_th1_file->Integral(3,6));//assumes bin 6 is overflow
  RooDataHist signal_twoelezh_hist("signal_twoelezh", "Signal yield in TwoEleZH", vars, &signal_twoelezh_th1);
  wspace->import(signal_twoelezh_hist);

}


//---------------------------------------------------------------------------------------------------------------
//EleMu Region
//---------------------------------------------------------------------------------------------------------------
void build_elemu(RooWorkspace *wspace){

  RooRealVar* ntags = wspace->var("ntags");
  RooArgList vars(*ntags);


  //Data
  TFile* f_elemu = TFile::Open("../inputs/EleMu_Combined_"+xvar+xvar_suffix+".root", "READ");
  TH1F* data_elemu_th1_file = (TH1F*)f_elemu->Get(data_string);
  TH1F data_elemu_th1("data_obs_elemu","Data observed in EleMu", 3, -0.5, 2.5);
  data_elemu_th1.SetBinContent(1, data_elemu_th1_file->GetBinContent(1));
  data_elemu_th1.SetBinContent(2, data_elemu_th1_file->GetBinContent(2));
  data_elemu_th1.SetBinContent(3, data_elemu_th1_file->Integral(3,6));//assumes bin 6 is overflow TODO
  RooDataHist data_elemu_hist("data_obs_elemu", "Data observed in EleMu", vars, &data_elemu_th1);
  wspace->import(data_elemu_hist);


  //Light contamination
  TH1F* light_elemu_th1_file = (TH1F*)f_elemu->Get("light");
  TH1F light_elemu_th1("light_elemu","Light yield in EleMu", 3, -0.5, 2.5);
  light_elemu_th1.SetBinContent(1, light_elemu_th1_file->GetBinContent(1));
  light_elemu_th1.SetBinContent(2, light_elemu_th1_file->GetBinContent(2));
  light_elemu_th1.SetBinContent(3, light_elemu_th1_file->Integral(3,6));//assumes bin 6 is overflow TODO
  RooDataHist light_elemu_hist("light_elemu", "Light yield in EleMu", vars, &light_elemu_th1);
  wspace->import(light_elemu_hist);


  //Other contamination
  TH1F* other_elemu_th1_file = (TH1F*)f_elemu->Get("other");
  TH1F other_elemu_th1("other_elemu","Other yield in EleMu", 3, -0.5, 2.5);
  other_elemu_th1.SetBinContent(1, other_elemu_th1_file->GetBinContent(1));
  other_elemu_th1.SetBinContent(2, other_elemu_th1_file->GetBinContent(2));
  other_elemu_th1.SetBinContent(3, other_elemu_th1_file->Integral(3,6));//assumes bin 6 is overflow TODO
  RooDataHist other_elemu_hist("other_elemu", "Other yield in EleMu", vars, &other_elemu_th1);
  wspace->import(other_elemu_hist);


  //Heavy background
  //Create one parameter per bin representing the yield
  TH1F* h_heavy_elemu = (TH1F*)f_elemu->Get("heavy");
  RooRealVar heavy_elemu_bin1("heavy_elemu_bin1", "Heavy background yield in EleMu, bin 1",
			      h_heavy_elemu->GetBinContent(1), h_heavy_elemu->GetBinContent(1)*0.5,h_heavy_elemu->GetBinContent(1)*1.5 );
  RooRealVar heavy_elemu_bin2("heavy_elemu_bin2", "Heavy background yield in EleMu, bin 2",
			      h_heavy_elemu->GetBinContent(2), h_heavy_elemu->GetBinContent(2)*0.5,h_heavy_elemu->GetBinContent(2)*1.5 );
  RooRealVar heavy_elemu_bin3("heavy_elemu_bin3", "Heavy background yield in EleMu, bin 3",
			      h_heavy_elemu->GetBinContent(3), 0 ,(h_heavy_elemu->GetBinContent(3)+10)*1.5 );
  RooArgList heavy_elemu_bins;
  heavy_elemu_bins.add(heavy_elemu_bin1);
  heavy_elemu_bins.add(heavy_elemu_bin2);
  heavy_elemu_bins.add(heavy_elemu_bin3);
  RooParametricHist p_heavy_elemu("heavy_elemu", "Heavy PDF in EleMu Region", *ntags, heavy_elemu_bins, data_elemu_th1);
  RooAddition p_heavy_elemu_norm("heavy_elemu_norm", "Total number of heavy events in EleMu Region", heavy_elemu_bins);
  wspace->import(p_heavy_elemu);
  wspace->import(p_heavy_elemu_norm, RooFit::RecycleConflictNodes());

}


//---------------------------------------------------------------------------------------------------------------
//EleMu-Low PT Region
//---------------------------------------------------------------------------------------------------------------
void build_elemul(RooWorkspace *wspace){

  RooRealVar* ntags = wspace->var("ntags");
  RooArgList vars(*ntags);


  //Data
  TFile* f_elemul = TFile::Open("../inputs/EleMuOSOFL_"+xvar+xvar_suffix+".root", "READ");
  TH1F* data_elemul_th1_file = (TH1F*)f_elemul->Get(data_string);
  TH1F data_elemul_th1("data_obs_elemul","Data observed in EleMuL", 3, -0.5, 2.5);
  data_elemul_th1.SetBinContent(1, data_elemul_th1_file->GetBinContent(1));
  data_elemul_th1.SetBinContent(2, data_elemul_th1_file->GetBinContent(2));
  data_elemul_th1.SetBinContent(3, data_elemul_th1_file->Integral(3,6));//assumes bin 6 is overflow TODO
  RooDataHist data_elemul_hist("data_obs_elemul", "Data observed in EleMuL", vars, &data_elemul_th1);
  wspace->import(data_elemul_hist);


  //Light contamination
  TH1F* light_elemul_th1_file = (TH1F*)f_elemul->Get("light");
  TH1F light_elemul_th1("light_elemul","Light yield in EleMuL", 3, -0.5, 2.5);
  light_elemul_th1.SetBinContent(1, light_elemul_th1_file->GetBinContent(1));
  light_elemul_th1.SetBinContent(2, light_elemul_th1_file->GetBinContent(2));
  light_elemul_th1.SetBinContent(3, light_elemul_th1_file->Integral(3,6));//assumes bin 6 is overflow TODO
  RooDataHist light_elemul_hist("light_elemul", "Light yield in EleMuL", vars, &light_elemul_th1);
  wspace->import(light_elemul_hist);


  //Other contamination
  TH1F* other_elemul_th1_file = (TH1F*)f_elemul->Get("other");
  TH1F other_elemul_th1("other_elemul","Other yield in EleMuL", 3, -0.5, 2.5);
  other_elemul_th1.SetBinContent(1, other_elemul_th1_file->GetBinContent(1));
  other_elemul_th1.SetBinContent(2, other_elemul_th1_file->GetBinContent(2));
  other_elemul_th1.SetBinContent(3, other_elemul_th1_file->Integral(3,6));//assumes bin 6 is overflow TODO
  RooDataHist other_elemul_hist("other_elemul", "Other yield in EleMuL", vars, &other_elemul_th1);
  wspace->import(other_elemul_hist);


  //Heavy background
  //Create one parameter per bin representing the yield
  TH1F* h_heavy_elemul = (TH1F*)f_elemul->Get("heavy");
  RooRealVar heavy_elemul_bin1("heavy_elemul_bin1", "Heavy background yield in EleMuL, bin 1",
			      h_heavy_elemul->GetBinContent(1), h_heavy_elemul->GetBinContent(1)*0.5,h_heavy_elemul->GetBinContent(1)*1.5 );
  RooRealVar heavy_elemul_bin2("heavy_elemul_bin2", "Heavy background yield in EleMuL, bin 2",
			      h_heavy_elemul->GetBinContent(2), h_heavy_elemul->GetBinContent(2)*0.5,h_heavy_elemul->GetBinContent(2)*1.5 );
  RooRealVar heavy_elemul_bin3("heavy_elemul_bin3", "Heavy background yield in EleMuL, bin 3",
			       h_heavy_elemul->GetBinContent(3), 0, (h_heavy_elemul->GetBinContent(3)+10)*1.5);
  RooArgList heavy_elemul_bins;
  heavy_elemul_bins.add(heavy_elemul_bin1);
  heavy_elemul_bins.add(heavy_elemul_bin2);
  heavy_elemul_bins.add(heavy_elemul_bin3);
  RooParametricHist p_heavy_elemul("heavy_elemul", "Heavy PDF in EleMuL Region", *ntags, heavy_elemul_bins, data_elemul_th1);
  RooAddition p_heavy_elemul_norm("heavy_elemul_norm", "Total number of heavy events in EleMuL Region", heavy_elemul_bins);
  wspace->import(p_heavy_elemul);
  wspace->import(p_heavy_elemul_norm, RooFit::RecycleConflictNodes());

}


//---------------------------------------------------------------------------------------------------------------
//OnePho Region
//---------------------------------------------------------------------------------------------------------------
void build_onepho(RooWorkspace* wspace){

  RooRealVar* ntags = wspace->var("ntags");
  RooArgList vars(*ntags);


  //Data
  TFile* f_onepho = TFile::Open("../inputs/OnePho_"+xvar+xvar_suffix+".root", "READ");
  TH1F* data_onepho_th1_file = (TH1F*)f_onepho->Get(data_string);
  TH1F data_onepho_th1("data_obs_onepho","Data observed in OnePho", 3, -0.5, 2.5);
  data_onepho_th1.SetBinContent(1, data_onepho_th1_file->GetBinContent(1));
  data_onepho_th1.SetBinContent(2, data_onepho_th1_file->GetBinContent(2));
  data_onepho_th1.SetBinContent(3, data_onepho_th1_file->Integral(3,6));//assumes bin 6 is overflow
  RooDataHist data_onepho_hist("data_obs_onepho", "Data observed in OnePho", vars, &data_onepho_th1);
  wspace->import(data_onepho_hist);


  //Other contamination
  TH1F* other_onepho_th1_file = (TH1F*)f_onepho->Get("other");
  TH1F other_onepho_th1("other_onepho","Other yield in OnePho", 3, -0.5, 2.5);
  other_onepho_th1.SetBinContent(1, other_onepho_th1_file->GetBinContent(1));
  other_onepho_th1.SetBinContent(2, other_onepho_th1_file->GetBinContent(2));
  other_onepho_th1.SetBinContent(3, other_onepho_th1_file->Integral(3,6));//assumes bin 6 is overflow TODO
  RooDataHist other_onepho_hist("other_onepho", "Other yield in OnePho", vars, &other_onepho_th1);
  wspace->import(other_onepho_hist);


  //Light background
  //Create one parameter per bin representing the yield
  TH1F* h_light_onepho = (TH1F*)f_onepho->Get("light");
  RooRealVar light_onepho_bin1("light_onepho_bin1", "Light background yield in OnePho, bin 1",
			       h_light_onepho->GetBinContent(1), h_light_onepho->GetBinContent(1)*0.5, h_light_onepho->GetBinContent(1)*1.5);
  RooRealVar light_onepho_bin2("light_onepho_bin2", "Light background yield in OnePho, bin 2",
			       h_light_onepho->GetBinContent(2), h_light_onepho->GetBinContent(2)*0.5, h_light_onepho->GetBinContent(2)*1.5);
  RooRealVar light_onepho_bin3("light_onepho_bin3", "Light background yield in OnePho, bin 3",
			       h_light_onepho->GetBinContent(3), 0, (h_light_onepho->GetBinContent(3)+10)*1.5);
  RooArgList light_onepho_bins;
  light_onepho_bins.add(light_onepho_bin1);
  light_onepho_bins.add(light_onepho_bin2);
  light_onepho_bins.add(light_onepho_bin3);
  RooParametricHist p_light_onepho("light_onepho", "Light PDF in OnePho Region", *ntags, light_onepho_bins, data_onepho_th1);
  RooAddition p_light_onepho_norm("light_onepho_norm", "Total number of light events in OnePho Region", light_onepho_bins);
  wspace->import(p_light_onepho);
  wspace->import(p_light_onepho_norm, RooFit::RecycleConflictNodes());


  //Get heavy
  build_tf(wspace, "heavy", "elemu", "onepho", sys_vec);
  RooFormulaVar* tf_heavy_elemu_to_onepho_bin1 = (RooFormulaVar*)wspace->arg("tf_heavy_elemu_to_onepho_bin1");
  RooFormulaVar* tf_heavy_elemu_to_onepho_bin2 = (RooFormulaVar*)wspace->arg("tf_heavy_elemu_to_onepho_bin2");
  RooFormulaVar* tf_heavy_elemu_to_onepho_bin3 = (RooFormulaVar*)wspace->arg("tf_heavy_elemu_to_onepho_bin3");

  RooRealVar* heavy_elemu_bin1 = wspace->var("heavy_elemu_bin1");
  RooRealVar* heavy_elemu_bin2 = wspace->var("heavy_elemu_bin2");
  RooRealVar* heavy_elemu_bin3 = wspace->var("heavy_elemu_bin3");

  RooFormulaVar heavy_onepho_bin1("heavy_onepho_bin1", "Heavy background yield in OnePho, bin 1",
				  "@0*@1", RooArgList(*tf_heavy_elemu_to_onepho_bin1, *heavy_elemu_bin1));
  RooFormulaVar heavy_onepho_bin2("heavy_onepho_bin2", "Heavy background yield in OnePho, bin 2",
				  "@0*@1", RooArgList(*tf_heavy_elemu_to_onepho_bin2, *heavy_elemu_bin2));
  RooFormulaVar heavy_onepho_bin3("heavy_onepho_bin3", "Heavy background yield in OnePho, bin 3",
				  "@0*@1", RooArgList(*tf_heavy_elemu_to_onepho_bin3, *heavy_elemu_bin3));

  RooArgList heavy_onepho_bins;
  heavy_onepho_bins.add(heavy_onepho_bin1);
  heavy_onepho_bins.add(heavy_onepho_bin2);
  heavy_onepho_bins.add(heavy_onepho_bin3);
  RooParametricHist p_heavy_onepho("heavy_onepho", "Heavy PDF in OnePho Region", *ntags, heavy_onepho_bins, data_onepho_th1);
  RooAddition p_heavy_onepho_norm("heavy_onepho_norm", "Total number of heavy events in OnePho Region", heavy_onepho_bins);

  wspace->import(p_heavy_onepho);
  wspace->import(p_heavy_onepho_norm, RooFit::RecycleConflictNodes());

}



//---------------------------------------------------------------------------------------------------------------
//TwoMuDY Region
//---------------------------------------------------------------------------------------------------------------
void build_twomudy(RooWorkspace* wspace){

  RooRealVar* ntags = wspace->var("ntags");
  RooArgList vars(*ntags);


  //Data
  TFile* f_twomudy = TFile::Open("../inputs/TwoMuDY_"+xvar+xvar_suffix+".root", "READ");
  TH1F* data_twomudy_th1_file = (TH1F*)f_twomudy->Get(data_string);
  TH1F data_twomudy_th1("data_obs_twomudy","Data observed in TwoMuDY", 3, -0.5, 2.5);
  data_twomudy_th1.SetBinContent(1, data_twomudy_th1_file->GetBinContent(1));
  data_twomudy_th1.SetBinContent(2, data_twomudy_th1_file->GetBinContent(2));
  data_twomudy_th1.SetBinContent(3, data_twomudy_th1_file->Integral(3,6));//assumes bin 6 is overflow
  RooDataHist data_twomudy_hist("data_obs_twomudy", "Data observed in TwoMuDY", vars, &data_twomudy_th1);
  wspace->import(data_twomudy_hist);


  //Other contamination
  TH1F* other_twomudy_th1_file = (TH1F*)f_twomudy->Get("other");
  TH1F other_twomudy_th1("other_twomudy","Other yield in TwoMuDY", 3, -0.5, 2.5);
  other_twomudy_th1.SetBinContent(1, other_twomudy_th1_file->GetBinContent(1));
  other_twomudy_th1.SetBinContent(2, other_twomudy_th1_file->GetBinContent(2));
  other_twomudy_th1.SetBinContent(3, other_twomudy_th1_file->Integral(3,6));//assumes bin 6 is overflow TODO
  RooDataHist other_twomudy_hist("other_twomudy", "Other yield in TwoMuDY", vars, &other_twomudy_th1);
  wspace->import(other_twomudy_hist);


  //Light background
  //Create one parameter per bin representing the yield
  TH1F* h_light_twomudy = (TH1F*)f_twomudy->Get("light");
  RooRealVar light_twomudy_bin1("light_twomudy_bin1", "Light background yield in TwoMuDY, bin 1",
				h_light_twomudy->GetBinContent(1), h_light_twomudy->GetBinContent(1)*0.5, h_light_twomudy->GetBinContent(1)*1.5);
  RooRealVar light_twomudy_bin2("light_twomudy_bin2", "Light background yield in TwoMuDY, bin 2",
				h_light_twomudy->GetBinContent(2), h_light_twomudy->GetBinContent(2)*0.5, h_light_twomudy->GetBinContent(2)*1.5);
  RooRealVar light_twomudy_bin3("light_twomudy_bin3", "Light background yield in TwoMuDY, bin 3",
				h_light_twomudy->GetBinContent(3), 0, (h_light_twomudy->GetBinContent(3)+10)*1.5);
  RooArgList light_twomudy_bins;
  light_twomudy_bins.add(light_twomudy_bin1);
  light_twomudy_bins.add(light_twomudy_bin2);
  light_twomudy_bins.add(light_twomudy_bin3);
  RooParametricHist p_light_twomudy("light_twomudy", "Light PDF in TwoMuDY Region", *ntags, light_twomudy_bins, data_twomudy_th1);
  RooAddition p_light_twomudy_norm("light_twomudy_norm", "Total number of light events in TwoMuDY Region", light_twomudy_bins);
  wspace->import(p_light_twomudy);
  wspace->import(p_light_twomudy_norm, RooFit::RecycleConflictNodes());


  //Get heavy
  build_tf(wspace, "heavy", "elemu", "twomudy", sys_vec);//ONE EMU CONTROL REGION -- NO Low EMU region anymore
  RooFormulaVar* tf_heavy_elemu_to_twomudy_bin1 = (RooFormulaVar*)wspace->arg("tf_heavy_elemu_to_twomudy_bin1");
  RooFormulaVar* tf_heavy_elemu_to_twomudy_bin2 = (RooFormulaVar*)wspace->arg("tf_heavy_elemu_to_twomudy_bin2");
  RooFormulaVar* tf_heavy_elemu_to_twomudy_bin3 = (RooFormulaVar*)wspace->arg("tf_heavy_elemu_to_twomudy_bin3");

  RooRealVar* heavy_elemu_bin1 = wspace->var("heavy_elemu_bin1");
  RooRealVar* heavy_elemu_bin2 = wspace->var("heavy_elemu_bin2");
  RooRealVar* heavy_elemu_bin3 = wspace->var("heavy_elemu_bin3");

  RooFormulaVar heavy_twomudy_bin1("heavy_twomudy_bin1", "Heavy background yield in TwoMuDY, bin 1",
				   "@0*@1", RooArgList(*tf_heavy_elemu_to_twomudy_bin1, *heavy_elemu_bin1));
  RooFormulaVar heavy_twomudy_bin2("heavy_twomudy_bin2", "Heavy background yield in TwoMuDY, bin 2",
				   "@0*@1", RooArgList(*tf_heavy_elemu_to_twomudy_bin2, *heavy_elemu_bin2));
  RooFormulaVar heavy_twomudy_bin3("heavy_twomudy_bin3", "Heavy background yield in TwoMuDY, bin 3",
				   "@0*@1", RooArgList(*tf_heavy_elemu_to_twomudy_bin3, *heavy_elemu_bin3));

  RooArgList heavy_twomudy_bins;
  heavy_twomudy_bins.add(heavy_twomudy_bin1);
  heavy_twomudy_bins.add(heavy_twomudy_bin2);
  heavy_twomudy_bins.add(heavy_twomudy_bin3);
  RooParametricHist p_heavy_twomudy("heavy_twomudy", "Heavy PDF in TwoMuDY Region", *ntags, heavy_twomudy_bins, data_twomudy_th1);
  RooAddition p_heavy_twomudy_norm("heavy_twomudy_norm", "Total number of heavy events in TwoMuDY Region", heavy_twomudy_bins);

  wspace->import(p_heavy_twomudy);
  wspace->import(p_heavy_twomudy_norm, RooFit::RecycleConflictNodes());


  //Signal
  TH1F* signal_twomudy_th1_file = (TH1F*)f_twomudy->Get(signal_string);
  TH1F signal_twomudy_th1("signal_twomudy","Signal yield in TwoMuDY", 3, -0.5, 2.5);
  signal_twomudy_th1.SetBinContent(1, signal_twomudy_th1_file->GetBinContent(1));
  signal_twomudy_th1.SetBinContent(2, signal_twomudy_th1_file->GetBinContent(2));
  signal_twomudy_th1.SetBinContent(3, signal_twomudy_th1_file->Integral(3,6));//assumes bin 6 is overflow
  cout << "TwoMuDY Sig Integral " <<  signal_twomudy_th1_file->Integral() << endl;
  RooDataHist signal_twomudy_hist("signal_twomudy", "Signal yield in TwoMuDY", vars, &signal_twomudy_th1);
  wspace->import(signal_twomudy_hist);

}


//---------------------------------------------------------------------------------------------------------------
//TwoEleDY Region
//---------------------------------------------------------------------------------------------------------------
void build_twoeledy(RooWorkspace* wspace){

  RooRealVar* ntags = wspace->var("ntags");
  RooArgList vars(*ntags);


  //Data
  TFile* f_twoeledy = TFile::Open("../inputs/TwoEleDY_"+xvar+xvar_suffix+".root", "READ");
  TH1F* data_twoeledy_th1_file = (TH1F*)f_twoeledy->Get(data_string);
  TH1F data_twoeledy_th1("data_obs_twoeledy","Data observed in TwoEleDY", 3, -0.5, 2.5);
  data_twoeledy_th1.SetBinContent(1, data_twoeledy_th1_file->GetBinContent(1));
  data_twoeledy_th1.SetBinContent(2, data_twoeledy_th1_file->GetBinContent(2));
  data_twoeledy_th1.SetBinContent(3, data_twoeledy_th1_file->Integral(3,6));//assumes bin 6 is overflow
  /*
  std::cout << "=============" << std::endl;
  std::cout << "fname: " << "../inputs/TwoEleDY_"+xvar+xvar_suffix+".root" << std::endl;
  std::cout << data_twoeledy_th1.GetBinContent(1) << " " << data_twoeledy_th1.GetBinContent(2)
  << " " << data_twoeledy_th1.GetBinContent(3) << std::endl;
  std::cout << "=============" << std::endl;
  */
  RooDataHist data_twoeledy_hist("data_obs_twoeledy", "Data observed in TwoEleDY", vars, &data_twoeledy_th1);
  wspace->import(data_twoeledy_hist);


  //Other contamination
  TH1F* other_twoeledy_th1_file = (TH1F*)f_twoeledy->Get("other");
  TH1F other_twoeledy_th1("other_twoeledy","Other yield in TwoEleDY", 3, -0.5, 2.5);
  other_twoeledy_th1.SetBinContent(1, other_twoeledy_th1_file->GetBinContent(1));
  other_twoeledy_th1.SetBinContent(2, other_twoeledy_th1_file->GetBinContent(2));
  other_twoeledy_th1.SetBinContent(3, other_twoeledy_th1_file->Integral(3,6));//assumes bin 6 is overflow TODO
  RooDataHist other_twoeledy_hist("other_twoeledy", "Other yield in TwoEleDY", vars, &other_twoeledy_th1);
  wspace->import(other_twoeledy_hist);


  //Light background
  //Create one parameter per bin representing the yield
  TH1F* h_light_twoeledy = (TH1F*)f_twoeledy->Get("light");
  RooRealVar light_twoeledy_bin1("light_twoeledy_bin1", "Light background yield in TwoEleDY, bin 1",
				h_light_twoeledy->GetBinContent(1), h_light_twoeledy->GetBinContent(1)*0.5, h_light_twoeledy->GetBinContent(1)*1.5);
  RooRealVar light_twoeledy_bin2("light_twoeledy_bin2", "Light background yield in TwoEleDY, bin 2",
				h_light_twoeledy->GetBinContent(2), h_light_twoeledy->GetBinContent(2)*0.5, h_light_twoeledy->GetBinContent(2)*1.5);
  RooRealVar light_twoeledy_bin3("light_twoeledy_bin3", "Light background yield in TwoEleDY, bin 3",
				 h_light_twoeledy->GetBinContent(3), 0, (h_light_twoeledy->GetBinContent(3)+10)*1.5);
  RooArgList light_twoeledy_bins;
  light_twoeledy_bins.add(light_twoeledy_bin1);
  light_twoeledy_bins.add(light_twoeledy_bin2);
  light_twoeledy_bins.add(light_twoeledy_bin3);
  RooParametricHist p_light_twoeledy("light_twoeledy", "Light PDF in TwoEleDY Region", *ntags, light_twoeledy_bins, data_twoeledy_th1);
  RooAddition p_light_twoeledy_norm("light_twoeledy_norm", "Total number of light events in TwoEleDY Region", light_twoeledy_bins);
  wspace->import(p_light_twoeledy);
  wspace->import(p_light_twoeledy_norm, RooFit::RecycleConflictNodes());


  //Get heavy
  RooRealVar* heavy_elemul_bin1 = wspace->var("heavy_elemul_bin1");
  RooRealVar* heavy_elemul_bin2 = wspace->var("heavy_elemul_bin2");
  RooRealVar* heavy_elemul_bin3 = wspace->var("heavy_elemul_bin3");

  build_tf(wspace, "heavy", "elemul", "twoeledy", sys_vec);
  RooFormulaVar* tf_heavy_elemul_to_twoeledy_bin1 = (RooFormulaVar*)wspace->arg("tf_heavy_elemul_to_twoeledy_bin1");
  RooFormulaVar* tf_heavy_elemul_to_twoeledy_bin2 = (RooFormulaVar*)wspace->arg("tf_heavy_elemul_to_twoeledy_bin2");
  RooFormulaVar* tf_heavy_elemul_to_twoeledy_bin3 = (RooFormulaVar*)wspace->arg("tf_heavy_elemul_to_twoeledy_bin3");

  RooFormulaVar heavy_twoeledy_bin1("heavy_twoeledy_bin1", "Heavy background yield in TwoEleDY, bin 1",
				   "@0*@1", RooArgList(*tf_heavy_elemul_to_twoeledy_bin1, *heavy_elemul_bin1));
  RooFormulaVar heavy_twoeledy_bin2("heavy_twoeledy_bin2", "Heavy background yield in TwoEleDY, bin 2",
				   "@0*@1", RooArgList(*tf_heavy_elemul_to_twoeledy_bin2, *heavy_elemul_bin2));
  RooFormulaVar heavy_twoeledy_bin3("heavy_twoeledy_bin3", "Heavy background yield in TwoEleDY, bin 3",
				   "@0*@1", RooArgList(*tf_heavy_elemul_to_twoeledy_bin3, *heavy_elemul_bin3));

  RooArgList heavy_twoeledy_bins;
  heavy_twoeledy_bins.add(heavy_twoeledy_bin1);
  heavy_twoeledy_bins.add(heavy_twoeledy_bin2);
  heavy_twoeledy_bins.add(heavy_twoeledy_bin3);
  RooParametricHist p_heavy_twoeledy("heavy_twoeledy", "Heavy PDF in TwoEleDY Region", *ntags, heavy_twoeledy_bins, data_twoeledy_th1);
  RooAddition p_heavy_twoeledy_norm("heavy_twoeledy_norm", "Total number of heavy events in TwoEleDY Region", heavy_twoeledy_bins);

  wspace->import(p_heavy_twoeledy);
  wspace->import(p_heavy_twoeledy_norm, RooFit::RecycleConflictNodes());


  //Signal
  std::cout << "[INFO]: Signal Model Name: " << signal_string << std::endl;
  TH1F* signal_twoeledy_th1_file = (TH1F*)f_twoeledy->Get(signal_string);
  TH1F signal_twoeledy_th1("signal_twoeledy","Signal yield in TwoEleDY", 3, -0.5, 2.5);
  signal_twoeledy_th1.SetBinContent(1, signal_twoeledy_th1_file->GetBinContent(1));
  signal_twoeledy_th1.SetBinContent(2, signal_twoeledy_th1_file->GetBinContent(2));
  signal_twoeledy_th1.SetBinContent(3, signal_twoeledy_th1_file->Integral(3,6));//assumes bin 6 is overflow
  RooDataHist signal_twoeledy_hist("signal_twoeledy", "Signal yield in TwoEleDY", vars, &signal_twoeledy_th1);
  wspace->import(signal_twoeledy_hist);

}



int main( int argc, char* argv[] ){

  // As usual, load the combine library to get access to the RooParametricHist
  //gSystem->Load("libHiggsAnalysisCombinedLimit.so");

  gROOT->Reset();
  //------------------------------
  // Getting Specific Signal Model
  //------------------------------
  std::string  signal_model = ParseCommandLine( argc, argv, "-signal_model=" );
  if (  signal_model == "" )
    {
      std::cerr << "[ERROR]: please provide a signal model --signal_model=<signal_model_name>" << std::endl;
      return -1;
    }

  signal_string = signal_model.c_str();
  //sys_vec.push_back("TagVars");
  //sys_vec.push_back("AMax");//use
  //sys_vec.push_back("IPSig");//use
  //sys_vec.push_back("TA");//use
  //sys_vec.push_back("EGS");
  //sys_vec.push_back("MES");
  //sys_vec.push_back("JES");

  // Output file and workspace
  TFile *fOut = new TFile("param_ws.root","RECREATE");
  RooWorkspace* wspace = new RooWorkspace("wspace","wspace");

  //Search in ntags, define ntags as our variable
  RooRealVar ntags("ntags", "ntags", -.5, 2.5);
  wspace->import(ntags, RooFit::RecycleConflictNodes());

  //Make RRVs for systematics
  for(unsigned int i=0; i<sys_vec.size(); i++){
    RooRealVar r("rrv_"+sys_vec[i], "rrv_"+sys_vec[i], 1, 0, 5);
    wspace->import(r, RooFit::RecycleConflictNodes());
  }

  //TwoMuZH+EleMu+OnePho top-nontop
  //build_elemu(wspace);
  //build_onepho(wspace);
  //build_twomuzh(wspace);

  //ZH+EleMu+EleMuL+DY top-nontop-other-signal
  std::cout << "===========MAKING EMU REGION===============" << std::endl;
  build_elemu(wspace);
  //build_elemul(wspace);
  std::cout << "===========MAKING Z-lowPT REGION===============" << std::endl;
  build_twomudy(wspace);
  //build_twoeledy(wspace);
  std::cout << "===========MAKING Z-highPT REGION===============" << std::endl;
  build_twomuzh(wspace, "DY");
  //build_twoelezh(wspace, "DY");

  wspace->Print("v");

  fOut->cd();
  wspace->Write();

  return 0;
}
