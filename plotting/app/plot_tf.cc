//C++
#include <iostream>
#include <fstream>
#include <map>
#include <stdlib.h>
#include <utility>
//ROOT
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
//LOCAL INCLUDES
#include <CommandLineInput.hh>
#include <helper_functions.hh>

const bool _info  = true;
const bool _debug = false;

/*
const float lumi = 5;
//Axis
const float axisTitleSize = 0.06;
const float axisTitleOffset = .8;

const float axisTitleSizeRatioX   = 0.18;
const float axisLabelSizeRatioX   = 0.12;
const float axisTitleOffsetRatioX = 0.94;

const float axisTitleSizeRatioY   = 0.15;
const float axisLabelSizeRatioY   = 0.108;
const float axisTitleOffsetRatioY = 0.32;

//Margins
const float leftMargin   = 0.12;
const float rightMargin  = 0.05;
const float topMargin    = 0.07;
const float bottomMargin = 0.12;

//CMS STANDARD
TString CMSText = "CMS";
//TString extraText   = "";
TString extraText   = "Preliminary";
//TString lumiText = "2.32 fb^{-1} (13 TeV)";
//TString lumiText = "35.9 fb^{-1} (13 TeV)";
TString lumiText = "16.2 fb^{-1} (13 TeV)";
*/
bool AddCMS( TCanvas* C );

int main( int argc, char** argv )
{

  //-----------------
  //Input File List
  //-----------------
  std::string inputDir = ParseCommandLine( argc, argv, "-inputDir=" );
  if (  inputDir == "" )
  {
    std::cerr << "[ERROR]: please provide an inputDir. Use --inputDir=" << std::endl;
    return -1;
  }

  if( _info )
  {
    std::cout << "[INFO]: input directory: " << inputDir << std::endl;
  }

  //-----------------------
  TString dy_mumu_zh  = "/TwoMuZH/GH/TwoMuZH_nSelectedAODCaloJetTag_GH.root";
  TString dy_ee_zh    = "/TwoEleZH/GH/TwoEleZH_nSelectedAODCaloJetTag_GH.root";
  //-----------------------
  TString dy_mumu_dy  = "/TwoMuDY/GH/TwoMuDY_nSelectedAODCaloJetTag_GH.root";
  TString dy_ee_dy    = "/TwoEleDY/GH/TwoEleDY_nSelectedAODCaloJetTag_GH.root";
  //-----------------------
  TString top_emu_zh  = "/EleMuOSOF/GH/EleMuOSOF_nSelectedAODCaloJetTag_GH.root";
  TString topl_emu_dy = "/EleMuOSOFL/GH/EleMuOSOFL_nSelectedAODCaloJetTag_GH.root";



  TFile* ftmp = new TFile("ftmp.root", "recreate");
  //------------------
  //Getting ROOT files
  //------------------
  TString tmp_file = inputDir.c_str() + dy_mumu_zh;

  TFile* fin = new TFile( tmp_file, "READ");
  if( fin == NULL )
  {
    std::cerr << "could not open root file: " << tmp_file << std::endl;
    return -1;
  }


  //get light from mumuZH
  std::cout << tmp_file << std::endl;
  TH1F* h_dy_mumu_zh = (TH1F*)fin->Get("light");
  delete fin;
  //get light from eeZH
  tmp_file = inputDir.c_str() + dy_ee_zh;
  std::cout << tmp_file << std::endl;
  fin = new TFile( tmp_file, "READ");
  TH1F* h_dy_ee_zh   = (TH1F*)fin->Get("light");
  delete fin;

  //get light from mumuDY
  tmp_file = inputDir.c_str() + dy_mumu_dy;
  std::cout << tmp_file << std::endl;
  fin = new TFile( tmp_file, "READ");
  TH1F* h_dy_mumu_dy = (TH1F*)fin->Get("light");
  delete fin;
  //get light from eeDY
  tmp_file = inputDir.c_str() + dy_ee_dy;
  std::cout << tmp_file << std::endl;
  fin = new TFile( tmp_file, "READ");
  TH1F* h_dy_ee_dy   = (TH1F*)fin->Get("light");
  delete fin;

  create_tf_plot(h_dy_mumu_zh, h_dy_mumu_dy, "tf_z_mumu");


  return 0;
}