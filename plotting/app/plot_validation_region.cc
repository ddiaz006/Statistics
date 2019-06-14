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
#include "CommandLineInput.hh"

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
  //std::ifstream ifs ( inputList.c_str(), std::ifstream::in );//input file list

  //----------------------
  //validation region name
  //----------------------
  std::string vrName = ParseCommandLine( argc, argv, "-vrName=" );
  if (  vrName == "" )
  {
    std::cerr << "[ERROR]: please provide validation region name. Use --vrName=" << std::endl;
    return -1;
  }


  if( _info )
  {
    std::cout << "[INFO]: input file: " << inputFile << std::endl;
    std::cout << "[INFO]: validation region name: " << vrName << std::endl;
  }

  TFile* fin = new TFile( inputFile.c_str(), "READ");
  TDirectory* dir_prefit_two_ee_dy = (TDirectory*)(((TDirectory*)fin->Get("shapes_prefit"))->Get("TwoEleDY"));

  TH1F* light_two_ee_dy = (TH1F*)dir_prefit_two_ee_dy->Get("light");
  TH1F* heavy_two_ee_dy = (TH1F*)dir_prefit_two_ee_dy->Get("heavy");
  TH1F* other_two_ee_dy = (TH1F*)dir_prefit_two_ee_dy->Get("other");
  TH1F* bkg_total_two_ee_dy = (TH1F*)dir_prefit_two_ee_dy->Get("total_background");
  TH1F* signal_two_ee_dy = (TH1F*)dir_prefit_two_ee_dy->Get("total_signal");
  TGraphAsymmErrors* data_two_ee_dy = (TGraphAsymmErrors*)dir_prefit_two_ee_dy->Get("data");



  TCanvas* c = new TCanvas( "c", "c", 2119, 33, 800, 700 );
  c->SetHighLightColor(2);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetLeftMargin( leftMargin );
  c->SetRightMargin( rightMargin );
  c->SetTopMargin( topMargin );
  c->SetBottomMargin( bottomMargin );
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);
  //c->SetLogy();
  //c->SetLogx();

  bkg_total_two_ee_dy->Draw("HISTO");
  data_two_ee_dy->Draw("P");

  gStyle->SetPaintTextFormat("4.3f");

  TLegend* leg = new TLegend( 0.51, 0.75-5*0.065, 0.85, 0.75-0.05, NULL, "brNDC" );
  leg->SetBorderSize(0);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetTextSize(0.04);

  //leg->AddEntry( gTheory, " NLO+NNL theory", "l");
  //leg->AddEntry( gObs, " Observed limit (95% CL)", "l" );
  //leg->AddEntry( gExp, " Median expected limit", "l" );
  //leg->AddEntry( gOneS, " 68% expected", "f" );
  //leg->AddEntry( gTwoS, " 95% expected", "f" );
  //leg->Draw("SAME");

  //95% CL label
  float cmsx = 0.81;
  float cmsy = 0.63-0.05;
  float cmsSize = 0.04;
  float cmsTextFont = 41;  // default is helvetic-bold
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);
  latex.SetTextAlign(31);
  latex.SetTextSize(cmsSize);
  latex.SetTextFont(cmsTextFont);
  //latex.DrawLatex(cmsx, cmsy, "95% CL upper limits");

  TLatex latex2;
  cmsx = 0.29;
  cmsy = 0.88;
  latex2.SetNDC();
  latex2.SetTextSize(0.038);
  latex2.SetTextFont(42);
  //HH
  //latex2.DrawLatex(cmsx, cmsy, "pp #rightarrow #tilde{#chi}^{0,#pm}_{i} #tilde{#chi}^{0,#pm}_{j} #rightarrow  #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} + X_{soft}; #tilde{#chi}^{0}_{1} #rightarrow H #tilde{G} (100%)");
  //latex2.DrawLatex(cmsx+0.263, cmsy-0.07, "m_{#tilde{#chi}^{0}_{2}} #approx m_{#tilde{#chi}^{#pm}_{1}} #approx m_{#tilde{#chi}^{0}_{1}};  m_{#tilde{G}} = 1 GeV");
  //HZ
  //latex2.DrawLatex(cmsx, cmsy, "pp #rightarrow #tilde{#chi}^{0,#pm}_{i} #tilde{#chi}^{0,#pm}_{j} #rightarrow  #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} + X_{soft}; #tilde{#chi}^{0}_{1} #rightarrow H #tilde{G} (50%)");
  latex2.DrawLatex(cmsx+0.15, cmsy, "pp #rightarrow ZH #rightarrow  Z #chi^{0}#chi^{0}; #chi}^{0} #rightarrow b#bar{b} (85%)");
  latex2.DrawLatex(cmsx+0.455, cmsy-0.07, "m_{#chi^{0}} = 40 GeV");
  //latex2.DrawLatex(cmsx+0.25, cmsy-0.13, "m_{#tilde{#chi}^{0}_{2}} #approx m_{#tilde{#chi}^{#pm}_{1}} #approx m_{#tilde{#chi}^{0}_{1}};  m_{#tilde{G}} = 1 GeV");
  //1D WH
  //latex2.DrawLatex(cmsx+0.15, cmsy, "pp #rightarrow #tilde{#chi}^{#pm}_{1} #tilde{#chi}^{0}_{2} ; #tilde{#chi}^{#pm}_{1} #rightarrow W^{#pm} #tilde{#chi}^{0}_{1}, #tilde{#chi}^{0}_{2} #rightarrow H #tilde{#chi}^{0}_{1}");
  //std::cout << "hola " << latex2.GetTextFont() << std::endl;
  TLatex latex3;
  latex3.SetNDC();
  latex3.SetTextSize(0.038);
  latex3.SetTextFont(42);
  //latex2.DrawLatex(0.2, 0.66, "#bf{EWP Analysis}");
  AddCMS(c);

  //c->SetLogx();
  c->SaveAs("NarrowResLimit_BIAS_fix.pdf");
  c->SaveAs("NarrowResLimit_BIAS_fix.C");

  //gObs->GetXaxis()->SetRangeUser(0, 30);
  //gObs->Write("gObs");
  //gExp->GetXaxis()->SetRangeUser(0, 30);
  //gExp->Write("gExp");
  //gOneS->Write("gOneS");
  //gTwoS->Write("gTwoS");

  //out->Close();
  return 0;
}


bool AddCMS( TCanvas* C )
{
  C->cd();
  float lumix = 0.955;
  float lumiy = 0.945;
  float lumifont = 42;

  float cmsx = 0.25;
  float cmsy = 0.94;
  float cmsTextFont   = 61;  // default is helvetic-bold
  float extrax = cmsx +0.20;
  float extray = cmsy;
  //float extrax = cmsx + 0.078;
  //float extray = cmsy - 0.04;
  float extraTextFont = 52;  // default is helvetica-italics
  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;
  float cmsSize = 0.06;
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);
  float extraTextSize = extraOverCmsTextSize*cmsSize;
  latex.SetTextFont(lumifont);
  latex.SetTextAlign(31);
  latex.SetTextSize(cmsSize);
  latex.DrawLatex(lumix, lumiy,lumiText);

  latex.SetTextFont(cmsTextFont);
  latex.SetTextAlign(31);
  latex.SetTextSize(cmsSize);
  latex.DrawLatex(cmsx, cmsy, CMSText);

  latex.SetTextFont(extraTextFont);
  latex.SetTextAlign(31);
  latex.SetTextSize(extraTextSize);
  latex.DrawLatex(extrax, extray, extraText);
  return true;
};
