#include <iostream>
#include <math.h>   

#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TMath.h"

using namespace std;


float r_err(float n, float dn, float d, float dd){
  return (n/d)*TMath::Sqrt(dn*dn/(n*n) + (dd*dd)/(d*d));
}

float r_err_n0(float n, float dn, float d, float dd){
  return dn/d;
}


void wavg(std::vector<TH1F*> hvec, int bin, float& a, float& da){
  float num = 0;
  float den = 0;
  for(unsigned int i=0; i<hvec.size(); i++){
    if(hvec[i]->GetBinContent(bin)<1e-9) continue; //skip zeros for now
    float w = 1.0 / ( hvec[i]->GetBinError(bin)*hvec[i]->GetBinError(bin) );
    num+=w*hvec[i]->GetBinContent(bin);
    den+=w;
  }
  a = num/den;
  da = 1.0 / TMath::Sqrt(den);
}


TString translate(TString in){
  TString out = "";
  if(in=="elemu") out = "EleMuOSOF";
  else if(in == "onepho") out = "OnePho";
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


void print_hist(TH1F* h, TString name){
  std::cout << name << " ";
  for(int i=1; i<=h->GetNbinsX(); i++){
    std::cout << h->GetBinContent(i) << " ";
  }
  std::cout << endl;
}


void print_hist_err(TH1F* h, TString name){
  std::cout << name << " ";
  for(int i=1; i<=h->GetNbinsX(); i++){
    std::cout << h->GetBinContent(i) << " +/- " << h->GetBinError(i) << " " ;
  }
  std::cout << endl;
}


TH1F* make_h_r(TString process, TString from_name, TString to_name, TString SB){

  TString full_name = process+"_"+from_name+"_to_"+to_name;  
  if(SB.Contains("SB")){ full_name += "_"; full_name += SB; }
  
  TFile* f_from = TFile::Open("../inputs/"+translate(from_name)+"_nSelectedAODCaloJetTag"+SB+"_GH.root", "READ");
  TH1F* h_from = (TH1F*)f_from->Get(process);
  
  TFile* f_to = TFile::Open("../inputs/"+translate(to_name)+"_nSelectedAODCaloJetTag"+SB+"_GH.root", "READ");
  TH1F* h_to = (TH1F*)f_to->Get(process);
  
  TH1F* h_r = (TH1F*)h_to->Clone("h_r_"+full_name);
  h_r->Divide(h_from);  
  
  h_r->SetTitle(SB);

  /*
  print_hist(h_to, full_name+" TO");
  print_hist(h_from, full_name+" FROM");
  print_hist(h_r, full_name);
  */

  /*
  //2 tag bin
  double n=h_to->GetBinContent(3);
  double dn=h_to->GetBinError(3);
  double d=h_from->GetBinContent(3);
  double dd=h_from->GetBinError(3);
  if(d>1e-9){//only print if d>0
    cout << "R " << n/d << " " << r_err(n, dn, d, dd) << endl;
    if(n<1e-9){//only try this if n=0
      cout << "Ro " << n/d << " " << r_err_n0(n, 1.8, d, dd) << endl;
    }
  }
  */

  return h_r;
}


void sb_data(TString from_name, TString to_name){

  TString SB = "SB"; //SB, SBL
  TString full_name = "DataMC_"+from_name+"_to_"+to_name;
  
  std::vector<TH1F*> rvec_mc;
  rvec_mc.push_back(make_h_r("bkgtotal", from_name, to_name, SB+"1"));
  rvec_mc.push_back(make_h_r("bkgtotal", from_name, to_name, SB+"2"));
  rvec_mc.push_back(make_h_r("bkgtotal", from_name, to_name, SB+"3"));
  rvec_mc.push_back(make_h_r("bkgtotal", from_name, to_name, SB+"4"));
  rvec_mc.push_back(make_h_r("bkgtotal", from_name, to_name, SB+"5"));
  rvec_mc.push_back(make_h_r("bkgtotal", from_name, to_name, SB+"6"));
  rvec_mc.push_back(make_h_r("bkgtotal", from_name, to_name, SB+"7"));

  std::vector<TH1F*> rvec_data;
  rvec_data.push_back(make_h_r("Data", from_name, to_name, SB+"1"));
  rvec_data.push_back(make_h_r("Data", from_name, to_name, SB+"2"));
  rvec_data.push_back(make_h_r("Data", from_name, to_name, SB+"3"));
  rvec_data.push_back(make_h_r("Data", from_name, to_name, SB+"4"));
  rvec_data.push_back(make_h_r("Data", from_name, to_name, SB+"5"));
  rvec_data.push_back(make_h_r("Data", from_name, to_name, SB+"6"));
  rvec_data.push_back(make_h_r("Data", from_name, to_name, SB+"7"));

  for(unsigned int i=0; i<rvec_data.size(); i++){
    
    rvec_mc[i]->SetTitle(full_name+"_"+rvec_mc[i]->GetTitle());
    rvec_mc[i]->SetLineWidth(2);
    rvec_mc[i]->SetLineColor(kRed);
    rvec_mc[i]->SetFillStyle(0);
    rvec_data[i]->SetLineWidth(2);
    rvec_data[i]->SetLineColor(kBlack);
    rvec_data[i]->SetFillStyle(0);

    TString name = rvec_data[i]->GetTitle(); name+= +"_data_mc";
    TH1F* h_data_mc = (TH1F*)rvec_data[i]->Clone(name);
    h_data_mc->Divide(rvec_mc[i]);
    h_data_mc->SetTitle("Data/MC");

    float max = rvec_mc[i]->GetMaximum();
    if(rvec_data[i]->GetMaximum()>max) max = rvec_data[i]->GetMaximum();
    if(max>10){
      max = 2;
      rvec_mc[i]->SetMinimum(-.2);
    }
    rvec_mc[i]->SetMaximum(1.25*max);

    TLegend *leg;
    leg = new TLegend(0.15,0.78,0.4,0.88);
    leg->SetBorderSize(0);
    leg->SetNColumns(2);
    leg->SetFillStyle(0);
    leg->AddEntry(rvec_mc[i], "MC", "l");
    leg->AddEntry(rvec_data[i], "Data", "l");
    
    TCanvas c("c", "c", 640, 2*480);
    c.Divide(1,2);
    c.cd(1);
    rvec_mc[i]->Draw("HIST E1");
    rvec_data[i]->Draw("HIST E1 SAME");
    leg->Draw();
    c.cd(2);
    h_data_mc->Draw("HIST E1");
    c.SaveAs(full_name+"_"+rvec_data[i]->GetTitle()+".pdf");
  }

}

void plot_tf(TString process, TString from_name, TString to_name){
  
  TString SB = "SB";//SB, SBL
  TString full_name = process+"_"+from_name+"_to_"+to_name;

  ////////////////
  //Prepare TFs
  /////////////////
  TH1F* h_r = make_h_r(process, from_name, to_name, "");//SIG
  h_r->SetTitle(full_name);
  print_hist_err(h_r, full_name);

  std::vector<TH1F*> rvec;
  rvec.push_back(make_h_r(process, from_name, to_name, SB+"1"));
  rvec.push_back(make_h_r(process, from_name, to_name, SB+"2"));
  rvec.push_back(make_h_r(process, from_name, to_name, SB+"3"));
  rvec.push_back(make_h_r(process, from_name, to_name, SB+"4"));
  rvec.push_back(make_h_r(process, from_name, to_name, SB+"5"));
  rvec.push_back(make_h_r(process, from_name, to_name, SB+"6"));
  rvec.push_back(make_h_r(process, from_name, to_name, SB+"7"));


  /////////////////////
  // Weighted average
  /////////////////////
  float a0, a1, a2;
  float da0, da1, da2;
  wavg(rvec,1,a0,da0);
  wavg(rvec,2,a1,da1);
  wavg(rvec,3,a2,da2);
  cout << "Weighted average of "+SB+"s, 0 tag: " << a0 << " +/- " << da0 << endl;
  cout << "Weighted average of "+SB+"s, 1 tag: " << a1 << " +/- " << da1 << endl;
  cout << "Weighted average of "+SB+"s, 2 tag: " << a2 << " +/- " << da2 << endl;


  ///////////////
  // Plot 
  ///////////////
  h_r->SetTitle(full_name);
  h_r->SetLineColor(kBlack);
  h_r->SetLineWidth(3);
  h_r->SetMarkerStyle(8);
  h_r->SetMarkerSize(1);
  for(unsigned int i=0; i<rvec.size(); i++){
    rvec[i]->SetLineColor(2+i);
    rvec[i]->SetLineWidth(2);
  }

  //Compute max
  double max = h_r->GetMaximum();
  for(unsigned int i=0; i<rvec.size(); i++){
    if(rvec[i]->GetMaximum()>max) max = rvec[i]->GetMaximum();
  }
  h_r->SetMaximum(max*1.5);

  TLegend *leg;
  leg = new TLegend(0.15,0.7,0.6,0.88);
  leg->SetBorderSize(0);
  leg->SetNColumns(2);
  //leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->AddEntry(h_r, "SIG", "lpe");
  for(unsigned int i=0; i<rvec.size(); i++){
    leg->AddEntry(rvec[i], rvec[i]->GetTitle(), "l");
  }

  TCanvas c("c", "c", 640, 480);
  h_r->Draw("P0 E1");
  for(unsigned int i=0; i<rvec.size(); i++){
    rvec[i]->Draw("HIST E1 SAME");
  }
  leg->Draw();
  c.SaveAs(SB+"_"+full_name+".pdf");

}


void tf(){

  //plot_tf("heavy", "elemu", "twomuzh");
  plot_tf("light", "twomudy", "twomuzh");
  //plot_tf("heavy", "elemu", "twoelezh");
  plot_tf("light", "twoeledy", "twoelezh");
  //plot_tf("heavy", "elemu", "twomudy");
  //plot_tf("heavy", "elemu", "twoeledy");

  //plot_tf("light", "onepho", "twomuzh");
  //plot_tf("light", "onepho", "twoelezh");

}
