#include <iostream>

#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"

using namespace std;


float r_err(float n, float dn, float d, float dd){
  return (n/d)*TMath::Sqrt(dn*dn/(n*n) + (dd*dd)/(d*d));
}


void print_hist(TH1F* h, TString name){
  std::cout << name << " ";
  for(int i=1; i<=h->GetNbinsX(); i++){
    std::cout << h->GetBinContent(i) << " ";
  }
  std::cout << endl;
}


void extract_ratios(){

  TFile* f_elemu   = TFile::Open("../inputs/EleMuOSOF_nSelectedAODCaloJetTag_GH.root", "READ");
  TH1F* h_heavy_elemu = (TH1F*)f_elemu->Get("heavy");

  TFile* f_onepho  = TFile::Open("../inputs/OnePho_nSelectedAODCaloJetTag_GH.root", "READ");
  TH1F* h_heavy_onepho = (TH1F*)f_onepho->Get("heavy");
  TH1F* h_light_onepho = (TH1F*)f_onepho->Get("light");

  TFile* f_twomudy = TFile::Open("../inputs/TwoMuDY_nSelectedAODCaloJetTag_GH.root", "READ");
  TH1F* h_heavy_twomudy = (TH1F*)f_twomudy->Get("heavy");
  TH1F* h_light_twomudy = (TH1F*)f_twomudy->Get("light");

  TFile* f_twomuzh = TFile::Open("../inputs/TwoMuZH_nSelectedAODCaloJetTag_GH.root", "READ");
  TH1F* h_heavy_twomuzh = (TH1F*)f_twomuzh->Get("heavy");
  TH1F* h_light_twomuzh = (TH1F*)f_twomuzh->Get("light");
  
  //print
  print_hist(h_heavy_elemu, "h_heavy_elemu");
  print_hist(h_heavy_onepho, "h_heavy_onepho");
  print_hist(h_light_onepho, "h_light_onepho");
  print_hist(h_heavy_twomudy, "h_heavy_twomudy");
  print_hist(h_light_twomudy, "h_light_twomudy");
  print_hist(h_heavy_twomuzh, "h_heavy_twomuzh");
  print_hist(h_light_twomuzh, "h_light_twomuzh");
  std::cout << std::endl;

  //heavy elemu-to-onepho
  Double_t n_heavy_onepho_err;
  float n_heavy_onepho = h_heavy_onepho->IntegralAndError(1,6,n_heavy_onepho_err);

  Double_t n_heavy_elemu_err;
  float n_heavy_elemu = h_heavy_elemu->IntegralAndError(1,6,n_heavy_elemu_err);

  float r_heavy_elemu_to_onepho = n_heavy_onepho/n_heavy_elemu;
  float r_heavy_elemu_to_onepho_err = r_err(n_heavy_onepho, n_heavy_onepho_err, n_heavy_elemu, n_heavy_elemu_err);

  std::cout << "heavy onepho " << n_heavy_onepho << " +/- " << n_heavy_onepho_err  << std::endl;
  std::cout << "heavy elemu "  << n_heavy_elemu  << " +/- " << n_heavy_elemu       << std::endl;
  std::cout << "heavy elemu to onepho " << r_heavy_elemu_to_onepho << " +/- " << r_heavy_elemu_to_onepho_err << std::endl;
  std::cout << std::endl;

  //heavy elemu-to-twomuzh
  Double_t n_heavy_twomuzh_err;
  float n_heavy_twomuzh = h_heavy_twomuzh->IntegralAndError(1,6,n_heavy_twomuzh_err);

  float r_heavy_elemu_to_twomuzh = n_heavy_twomuzh/n_heavy_elemu;
  float r_heavy_elemu_to_twomuzh_err = r_err(n_heavy_twomuzh, n_heavy_twomuzh_err, n_heavy_elemu, n_heavy_elemu_err);

  std::cout << "heavy twomuzh " << n_heavy_twomuzh << " +/- " << n_heavy_twomuzh_err  << std::endl;
  std::cout << "heavy elemu "  << n_heavy_elemu  << " +/- " << n_heavy_elemu       << std::endl;
  std::cout << "heavy elemu to twomuzh " << r_heavy_elemu_to_twomuzh << " +/- " << r_heavy_elemu_to_twomuzh_err << std::endl;
  std::cout << std::endl;

  //light onepho-to-twomuzh
  Double_t n_light_twomuzh_err;
  float n_light_twomuzh = h_light_twomuzh->IntegralAndError(1,6,n_light_twomuzh_err);

  Double_t n_light_onepho_err;
  float n_light_onepho = h_light_onepho->IntegralAndError(1,6,n_light_onepho_err);

  float r_light_onepho_to_twomuzh = n_light_twomuzh/n_light_onepho;
  float r_light_onepho_to_twomuzh_err = r_err(n_light_twomuzh, n_light_twomuzh_err, n_light_onepho, n_light_onepho_err);

  std::cout << "light twomuzh " << n_light_twomuzh << " +/- " << n_light_twomuzh_err  << std::endl;
  std::cout << "light onepho "  << n_light_onepho  << " +/- " << n_light_onepho       << std::endl;
  std::cout << "light onepho to twomuzh " << r_light_onepho_to_twomuzh << " +/- " << r_light_onepho_to_twomuzh_err << std::endl;
  std::cout << std::endl;

}
