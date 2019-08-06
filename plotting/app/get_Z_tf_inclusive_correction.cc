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

struct Limit
{
  double obs;
  double exp0p025;
  double exp0p16;
  double exp0p5;
  double exp0p84;
  double exp0p975;
};


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
  std::string inputFile = ParseCommandLine( argc, argv, "-inputFile=" );
  if (  inputFile == "" )
  {
    std::cerr << "[ERROR]: please provide an inputList. Use --inputFile=" << std::endl;
    return -1;
  }

  if( _info )
  {
    std::cout << "[INFO]: input file: " << inputFile << std::endl;
  }

  //-----------------
  //Getting ROOT file
  //-----------------
  TFile* fin = new TFile( inputFile.c_str(), "READ");

  if( fin == NULL )
  {
    std::cerr << "could not open root file: " << inputFile << std::endl;
    return -1;
  }

  TH1F* data  = (TH1F*)fin->Get("Data");
  TH1F* light = (TH1F*)fin->Get("light");
  TH1F* heavy = (TH1F*)fin->Get("heavy");
  TH1F* other = (TH1F*)fin->Get("other");

  //-------------------------
  //get control region yields
  //-------------------------
  //(binning is 10 GeV from 0-500)
  double dataCR_yield  = data->Integral(1,10);//from 10 to 100 GeV
  double lightCR_yield = light->Integral(1,10);//from 10 to 100 GeV
  double heavyCR_yield = heavy->Integral(1,10);//from 10 to 100 GeV
  double otherCR_yield = other->Integral(1,10);//from 10 to 100 GeV

  //-------------------------
  //get signal region yields
  //-------------------------
  //(binning is 10 GeV from 0-500)
  double dataSR_yield  = data->Integral(11,50);//from 100 to 500 GeV
  double lightSR_yield = light->Integral(11,50);//from 100 to 500 GeV
  double heavySR_yield = heavy->Integral(11,50);//from 100 to 500 GeV
  double otherSR_yield = other->Integral(11,50);//from 100 to 500 GeV

  double low_pt_relative_correction  = (dataCR_yield-heavyCR_yield-otherCR_yield)/lightCR_yield-1.;
  double high_pt_relative_correction = (dataSR_yield-heavySR_yield-otherSR_yield)/lightSR_yield-1.;

  std::cout << "low_pt_relative_correction: " << low_pt_relative_correction << std::endl;
  std::cout << "high_pt_relative_correction: " << high_pt_relative_correction << std::endl;
  std::cout << "==================================="<< std::endl;
  std::cout << "TF correction -> " << (1+high_pt_relative_correction)/(1+low_pt_relative_correction) << std::endl;
  std::cout << "==================================="<< std::endl;


  return 0;
}
