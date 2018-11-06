#include "TFile.h"
#include "TH1F.h"
#include "TDirectory.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include "TPad.h"

void make_one_plot(TString channel_name, TString fit_name){
  
  TFile* fin = TFile::Open("../model/fitDiagnostics.root");
  
  bool draw_signal = false;
  if( (channel_name.Contains("ZH")||channel_name.Contains("DY")) && fit_name!= "shapes_fit_b") draw_signal = true;

  TH1F* h_heavy=0;
  TH1F* h_light=0;
  TH1F* h_other=0;
  TH1F* h_signal=0;
  h_heavy = (TH1F*)fin->Get(fit_name+"/"+channel_name+"/heavy");
  h_light = (TH1F*)fin->Get(fit_name+"/"+channel_name+"/light");
  h_other = (TH1F*)fin->Get(fit_name+"/"+channel_name+"/other");
  if(draw_signal) h_signal = (TH1F*)fin->Get(fit_name+"/"+channel_name+"/signal");

  TGraphAsymmErrors* gr_data = (TGraphAsymmErrors*)fin->Get(fit_name+"/"+channel_name+"/data");
 
  //Style
  h_heavy->SetFillColor(kGreen+1);
  h_heavy->SetFillStyle(1001);
  h_light->SetFillColor(kAzure-3);
  h_light->SetFillStyle(1001);
  h_other->SetFillColor(kGray+1);
  h_other->SetFillStyle(1001);
  if(draw_signal){
    h_signal->SetFillColor(kRed);
    h_signal->SetFillStyle(1001);
  }
  gr_data->SetMarkerStyle(8);
  gr_data->SetMarkerSize(1);
  
  TCanvas c("c_"+channel_name, "c_"+channel_name, 640, 480);
  gStyle->SetOptStat(0);
  gPad->SetLogy(1);
  THStack *bgstack = new THStack("bgstack","");
  if(h_heavy->Integral()>h_light->Integral()){
    bgstack->Add(h_other);
    bgstack->Add(h_light);
    bgstack->Add(h_heavy);
  }
  else{
    bgstack->Add(h_other);
    bgstack->Add(h_heavy);
    bgstack->Add(h_light);
  }
  if(draw_signal) bgstack->Add(h_signal);

  bgstack->SetMinimum(0.1);
  bgstack->Draw("hist");
  gr_data->Draw("P");

  c.SaveAs("c_"+fit_name+"_"+channel_name+".pdf");
}

void make_plots(){
  //"shapes_fit_s"; //shapes_prefit, shapes_fit_b
  
  if(0){
    //Prefit
    make_one_plot("TwoMuZH", "shapes_prefit");
    make_one_plot("EleMu", "shapes_prefit");
    make_one_plot("OnePho", "shapes_prefit");
    
    //B only fit
    make_one_plot("TwoMuZH", "shapes_fit_b");
    make_one_plot("EleMu", "shapes_fit_b");
    make_one_plot("OnePho", "shapes_fit_b");
    
    //S+B fit
    make_one_plot("TwoMuZH", "shapes_fit_s");
    make_one_plot("EleMu", "shapes_fit_s");
    make_one_plot("OnePho", "shapes_fit_s");
  }
  else{
    //Prefit
    make_one_plot("TwoMuZH", "shapes_prefit");
    make_one_plot("EleMu", "shapes_prefit");
    make_one_plot("TwoMuDY", "shapes_prefit");
    
    //B only fit
    make_one_plot("TwoMuZH", "shapes_fit_b");
    make_one_plot("EleMu", "shapes_fit_b");
    make_one_plot("TwoMuDY", "shapes_fit_b");
    
    //S+B fit
    make_one_plot("TwoMuZH", "shapes_fit_s");
    make_one_plot("EleMu", "shapes_fit_s");
    make_one_plot("TwoMuDY", "shapes_fit_s");
  }
}

