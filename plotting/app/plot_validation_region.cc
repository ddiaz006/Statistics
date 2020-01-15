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

  //-----------------
  //Getting ROOT file
  //-----------------
  TFile* fin = new TFile( inputFile.c_str(), "READ");

  //----------------------------------------------------------
  //GETTING PRE-FIT RESULTS FROM COMBINE
  //----------------------------------------------------------
  //---------------------
  //TwoEleDY
  //---------------------
  TDirectory* dir_prefit_two_ee_dy = (TDirectory*)(((TDirectory*)fin->Get("shapes_prefit"))->Get("TwoEleDY"));
  TH1F* light_two_ee_dy = (TH1F*)dir_prefit_two_ee_dy->Get("light");
  TH1F* heavy_two_ee_dy = (TH1F*)dir_prefit_two_ee_dy->Get("heavy");
  TH1F* other_two_ee_dy = (TH1F*)dir_prefit_two_ee_dy->Get("other");
  TH1F* bkg_total_two_ee_dy = (TH1F*)dir_prefit_two_ee_dy->Get("total_background");
  TH1F* signal_two_ee_dy = (TH1F*)dir_prefit_two_ee_dy->Get("total_signal");
  TGraphAsymmErrors* data_two_ee_dy = (TGraphAsymmErrors*)dir_prefit_two_ee_dy->Get("data");

  THStack* stack_two_ee_dy = new THStack( "two_ee_dy_prefit" , "two_ee_dy combine-prefit" );
  create_stack( stack_two_ee_dy, light_two_ee_dy, heavy_two_ee_dy, other_two_ee_dy);
  create_ratio_plot(data_two_ee_dy, stack_two_ee_dy, bkg_total_two_ee_dy,
     Form("two_ee_dy_prefit_%s",vrName.c_str()), light_two_ee_dy, heavy_two_ee_dy, other_two_ee_dy);

  //---------------------
  //TwoMuDY
  //---------------------
  TDirectory* dir_prefit_two_mumu_dy = (TDirectory*)(((TDirectory*)fin->Get("shapes_prefit"))->Get("TwoMuDY"));
  TH1F* light_two_mumu_dy = (TH1F*)dir_prefit_two_mumu_dy->Get("light");
  TH1F* heavy_two_mumu_dy = (TH1F*)dir_prefit_two_mumu_dy->Get("heavy");
  TH1F* other_two_mumu_dy = (TH1F*)dir_prefit_two_mumu_dy->Get("other");
  TH1F* bkg_total_two_mumu_dy = (TH1F*)dir_prefit_two_mumu_dy->Get("total_background");
  TH1F* signal_two_mumu_dy = (TH1F*)dir_prefit_two_mumu_dy->Get("total_signal");
  TGraphAsymmErrors* data_two_mumu_dy = (TGraphAsymmErrors*)dir_prefit_two_mumu_dy->Get("data");

  THStack* stack_two_mumu_dy = new THStack( "two_mumu_dy_prefit" , "two_mumu_dy combine-prefit" );
  create_stack( stack_two_mumu_dy, light_two_mumu_dy, heavy_two_mumu_dy, other_two_mumu_dy);
  create_ratio_plot(data_two_mumu_dy, stack_two_mumu_dy, bkg_total_two_mumu_dy,
     Form("two_mumu_dy_prefit_%s",vrName.c_str()), light_two_mumu_dy, heavy_two_mumu_dy, other_two_mumu_dy);

  //---------------------
  //EleMu
  //---------------------
  TDirectory* dir_prefit_elemu_dy = (TDirectory*)(((TDirectory*)fin->Get("shapes_prefit"))->Get("EleMu"));
  TH1F* light_elemu_dy = (TH1F*)dir_prefit_elemu_dy->Get("light");
  TH1F* heavy_elemu_dy = (TH1F*)dir_prefit_elemu_dy->Get("heavy");
  TH1F* other_elemu_dy = (TH1F*)dir_prefit_elemu_dy->Get("other");
  TH1F* bkg_total_elemu_dy = (TH1F*)dir_prefit_elemu_dy->Get("total_background");
  TH1F* signal_elemu_dy = (TH1F*)dir_prefit_elemu_dy->Get("total_signal");
  TGraphAsymmErrors* data_elemu_dy = (TGraphAsymmErrors*)dir_prefit_elemu_dy->Get("data");

  THStack* stack_elemu_dy = new THStack( "elemu_dy_prefit" , "elemu_dy combine-prefit" );
  create_stack( stack_elemu_dy, light_elemu_dy, heavy_elemu_dy, other_elemu_dy);
  create_ratio_plot(data_elemu_dy, stack_elemu_dy, bkg_total_elemu_dy,
    Form("elemu_dy_prefit_%s",vrName.c_str()), light_elemu_dy, heavy_elemu_dy, other_elemu_dy);

  //---------------------
  //EleMuL
  //---------------------
  TDirectory* dir_prefit_elemul_dy = (TDirectory*)(((TDirectory*)fin->Get("shapes_prefit"))->Get("EleMuL"));
  TH1F* light_elemul_dy = (TH1F*)dir_prefit_elemul_dy->Get("light");
  TH1F* heavy_elemul_dy = (TH1F*)dir_prefit_elemul_dy->Get("heavy");
  TH1F* other_elemul_dy = (TH1F*)dir_prefit_elemul_dy->Get("other");
  TH1F* bkg_total_elemul_dy = (TH1F*)dir_prefit_elemul_dy->Get("total_background");
  TH1F* signal_elemul_dy = (TH1F*)dir_prefit_elemul_dy->Get("total_signal");
  TGraphAsymmErrors* data_elemul_dy = (TGraphAsymmErrors*)dir_prefit_elemul_dy->Get("data");

  THStack* stack_elemul_dy = new THStack( "elemul_dy_prefit" , "elemul_dy combine-prefit" );
  create_stack( stack_elemul_dy, light_elemul_dy, heavy_elemul_dy, other_elemul_dy);
  create_ratio_plot(data_elemul_dy, stack_elemul_dy, bkg_total_elemul_dy,
     Form("elemul_dy_prefit_%s",vrName.c_str()), light_elemul_dy, heavy_elemul_dy, other_elemul_dy);

  //---------------------
  //TwoEleZH
  //---------------------
  TDirectory* dir_prefit_two_ee_zh = (TDirectory*)(((TDirectory*)fin->Get("shapes_prefit"))->Get("TwoEleZH"));
  TH1F* light_two_ee_zh = (TH1F*)dir_prefit_two_ee_zh->Get("light");
  TH1F* heavy_two_ee_zh = (TH1F*)dir_prefit_two_ee_zh->Get("heavy");
  TH1F* other_two_ee_zh = (TH1F*)dir_prefit_two_ee_zh->Get("other");
  TH1F* bkg_total_two_ee_zh = (TH1F*)dir_prefit_two_ee_zh->Get("total_background");
  TH1F* signal_two_ee_zh = (TH1F*)dir_prefit_two_ee_zh->Get("total_signal");
  TGraphAsymmErrors* data_two_ee_zh = (TGraphAsymmErrors*)dir_prefit_two_ee_zh->Get("data");

  THStack* stack_two_ee_zh = new THStack( "two_ee_zh_prefit" , "two_ee_zh combine-prefit" );
  create_stack( stack_two_ee_zh, light_two_ee_zh, heavy_two_ee_zh, other_two_ee_zh);
  create_ratio_plot(data_two_ee_zh, stack_two_ee_zh, bkg_total_two_ee_zh,
     Form("two_ee_zh_prefit_%s",vrName.c_str()), light_two_ee_zh, heavy_two_ee_zh, other_two_ee_zh);

  //---------------------
  //TwoMuZH
  //---------------------
  TDirectory* dir_prefit_two_mumu_zh = (TDirectory*)(((TDirectory*)fin->Get("shapes_prefit"))->Get("TwoMuZH"));
  TH1F* light_two_mumu_zh = (TH1F*)dir_prefit_two_mumu_zh->Get("light");
  TH1F* heavy_two_mumu_zh = (TH1F*)dir_prefit_two_mumu_zh->Get("heavy");
  TH1F* other_two_mumu_zh = (TH1F*)dir_prefit_two_mumu_zh->Get("other");
  TH1F* bkg_total_two_mumu_zh = (TH1F*)dir_prefit_two_mumu_zh->Get("total_background");
  TH1F* signal_two_mumu_zh = (TH1F*)dir_prefit_two_mumu_zh->Get("total_signal");
  TGraphAsymmErrors* data_two_mumu_zh = (TGraphAsymmErrors*)dir_prefit_two_mumu_zh->Get("data");

  THStack* stack_two_mumu_zh = new THStack( "two_mumu_zh_prefit" , "two_mumu_zh combine-prefit" );
  create_stack( stack_two_mumu_zh, light_two_mumu_zh, heavy_two_mumu_zh, other_two_mumu_zh);
  create_ratio_plot(data_two_mumu_zh, stack_two_mumu_zh, bkg_total_two_mumu_zh,
     Form("two_mumu_zh_prefit_%s",vrName.c_str()), light_two_mumu_zh, heavy_two_mumu_zh, other_two_mumu_zh);


  //----------------------------------------------------------
  //GETTING POST-FIT RESULTS FROM COMBINE
  //----------------------------------------------------------
  //---------------------
  //TwoEleDY
  //---------------------
  TDirectory* dir_postfit_two_ee_dy = (TDirectory*)(((TDirectory*)fin->Get("shapes_fit_b"))->Get("TwoEleDY"));
  TH1F* light_two_ee_dy_pf = (TH1F*)dir_postfit_two_ee_dy->Get("light");
  TH1F* heavy_two_ee_dy_pf = (TH1F*)dir_postfit_two_ee_dy->Get("heavy");
  TH1F* other_two_ee_dy_pf = (TH1F*)dir_postfit_two_ee_dy->Get("other");
  TH1F* bkg_total_two_ee_dy_pf = (TH1F*)dir_postfit_two_ee_dy->Get("total_background");
  TH1F* signal_two_ee_dy_pf = (TH1F*)dir_postfit_two_ee_dy->Get("total_signal");
  TGraphAsymmErrors* data_two_ee_dy_pf = (TGraphAsymmErrors*)dir_postfit_two_ee_dy->Get("data");

  THStack* stack_two_ee_dy_pf = new THStack( "two_ee_dy_postfit" , "two_ee_dy combine-postfit" );
  create_stack( stack_two_ee_dy_pf, light_two_ee_dy_pf, heavy_two_ee_dy_pf, other_two_ee_dy_pf);
  create_ratio_plot(data_two_ee_dy_pf, stack_two_ee_dy_pf, bkg_total_two_ee_dy_pf,
    Form("two_ee_dy_postfit_%s",vrName.c_str()), light_two_ee_dy_pf, heavy_two_ee_dy_pf, other_two_ee_dy_pf);

  //---------------------
  //TwoMuDY
  //---------------------
  TDirectory* dir_postfit_two_mumu_dy = (TDirectory*)(((TDirectory*)fin->Get("shapes_fit_b"))->Get("TwoMuDY"));
  TH1F* light_two_mumu_dy_pf = (TH1F*)dir_postfit_two_mumu_dy->Get("light");
  TH1F* heavy_two_mumu_dy_pf = (TH1F*)dir_postfit_two_mumu_dy->Get("heavy");
  TH1F* other_two_mumu_dy_pf = (TH1F*)dir_postfit_two_mumu_dy->Get("other");
  TH1F* bkg_total_two_mumu_dy_pf = (TH1F*)dir_postfit_two_mumu_dy->Get("total_background");
  TH1F* signal_two_mumu_dy_pf = (TH1F*)dir_postfit_two_mumu_dy->Get("total_signal");
  TGraphAsymmErrors* data_two_mumu_dy_pf = (TGraphAsymmErrors*)dir_postfit_two_mumu_dy->Get("data");

  THStack* stack_two_mumu_dy_pf = new THStack( "two_mumu_dy_postfit" , "two_mumu_dy combine-postfit" );
  create_stack( stack_two_mumu_dy_pf, light_two_mumu_dy_pf, heavy_two_mumu_dy_pf, other_two_mumu_dy_pf);
  create_ratio_plot(data_two_mumu_dy_pf, stack_two_mumu_dy_pf, bkg_total_two_mumu_dy_pf,
     Form("two_mumu_dy_postfit_%s",vrName.c_str()), light_two_mumu_dy_pf, heavy_two_mumu_dy_pf, other_two_mumu_dy_pf);

  //---------------------
  //EleMu
  //---------------------
  TDirectory* dir_postfit_elemu_dy = (TDirectory*)(((TDirectory*)fin->Get("shapes_fit_b"))->Get("EleMu"));
  TH1F* light_elemu_dy_pf = (TH1F*)dir_postfit_elemu_dy->Get("light");
  TH1F* heavy_elemu_dy_pf = (TH1F*)dir_postfit_elemu_dy->Get("heavy");
  TH1F* other_elemu_dy_pf = (TH1F*)dir_postfit_elemu_dy->Get("other");
  TH1F* bkg_total_elemu_dy_pf = (TH1F*)dir_postfit_elemu_dy->Get("total_background");
  TH1F* signal_elemu_dy_pf = (TH1F*)dir_postfit_elemu_dy->Get("total_signal");
  TGraphAsymmErrors* data_elemu_dy_pf = (TGraphAsymmErrors*)dir_postfit_elemu_dy->Get("data");

  THStack* stack_elemu_dy_pf = new THStack( "elemu_dy_postfit" , "elemu_dy combine-postfit" );
  create_stack( stack_elemu_dy_pf, light_elemu_dy_pf, heavy_elemu_dy_pf, other_elemu_dy_pf);
  create_ratio_plot(data_elemu_dy_pf, stack_elemu_dy_pf, bkg_total_elemu_dy_pf,
     Form("elemu_dy_postfit_%s",vrName.c_str()), light_elemu_dy_pf, heavy_elemu_dy_pf, other_elemu_dy_pf);

  //---------------------
  //EleMuL
  //---------------------
  TDirectory* dir_postfit_elemul_dy = (TDirectory*)(((TDirectory*)fin->Get("shapes_fit_b"))->Get("EleMuL"));
  TH1F* light_elemul_dy_pf = (TH1F*)dir_postfit_elemul_dy->Get("light");
  TH1F* heavy_elemul_dy_pf = (TH1F*)dir_postfit_elemul_dy->Get("heavy");
  TH1F* other_elemul_dy_pf = (TH1F*)dir_postfit_elemul_dy->Get("other");
  TH1F* bkg_total_elemul_dy_pf = (TH1F*)dir_postfit_elemul_dy->Get("total_background");
  TH1F* signal_elemul_dy_pf = (TH1F*)dir_postfit_elemul_dy->Get("total_signal");
  TGraphAsymmErrors* data_elemul_dy_pf = (TGraphAsymmErrors*)dir_postfit_elemul_dy->Get("data");

  THStack* stack_elemul_dy_pf = new THStack( "elemul_dy_postfit" , "elemul_dy combine-postfit" );
  create_stack( stack_elemul_dy_pf, light_elemul_dy_pf, heavy_elemul_dy_pf, other_elemul_dy_pf);
  create_ratio_plot(data_elemul_dy_pf, stack_elemul_dy_pf, bkg_total_elemul_dy_pf,
     Form("elemul_dy_postfit_%s",vrName.c_str()), light_elemul_dy_pf, heavy_elemul_dy_pf, other_elemul_dy_pf);

  //---------------------
  //TwoEleZH
  //---------------------
  TDirectory* dir_postfit_two_ee_zh = (TDirectory*)(((TDirectory*)fin->Get("shapes_fit_b"))->Get("TwoEleZH"));
  TH1F* light_two_ee_zh_pf = (TH1F*)dir_postfit_two_ee_zh->Get("light");
  TH1F* heavy_two_ee_zh_pf = (TH1F*)dir_postfit_two_ee_zh->Get("heavy");
  TH1F* other_two_ee_zh_pf = (TH1F*)dir_postfit_two_ee_zh->Get("other");
  TH1F* bkg_total_two_ee_zh_pf = (TH1F*)dir_postfit_two_ee_zh->Get("total_background");
  TH1F* signal_two_ee_zh_pf = (TH1F*)dir_postfit_two_ee_zh->Get("total_signal");
  TGraphAsymmErrors* data_two_ee_zh_pf = (TGraphAsymmErrors*)dir_postfit_two_ee_zh->Get("data");

  THStack* stack_two_ee_zh_pf = new THStack( "two_ee_zh_postfit" , "two_ee_zh combine-postfit" );
  create_stack( stack_two_ee_zh_pf, light_two_ee_zh_pf, heavy_two_ee_zh_pf, other_two_ee_zh_pf);
  create_ratio_plot(data_two_ee_zh_pf, stack_two_ee_zh_pf, bkg_total_two_ee_zh_pf,
     Form("two_ee_zh_postfit_%s",vrName.c_str()), light_two_ee_zh_pf, heavy_two_ee_zh_pf, other_two_ee_zh_pf);

  //---------------------
  //TwoMuZH
  //---------------------
  TDirectory* dir_postfit_two_mumu_zh = (TDirectory*)(((TDirectory*)fin->Get("shapes_fit_b"))->Get("TwoMuZH"));
  TH1F* light_two_mumu_zh_pf = (TH1F*)dir_postfit_two_mumu_zh->Get("light");
  TH1F* heavy_two_mumu_zh_pf = (TH1F*)dir_postfit_two_mumu_zh->Get("heavy");
  TH1F* other_two_mumu_zh_pf = (TH1F*)dir_postfit_two_mumu_zh->Get("other");
  TH1F* bkg_total_two_mumu_zh_pf = (TH1F*)dir_postfit_two_mumu_zh->Get("total_background");
  TH1F* signal_two_mumu_zh_pf = (TH1F*)dir_postfit_two_mumu_zh->Get("total_signal");
  TGraphAsymmErrors* data_two_mumu_zh_pf = (TGraphAsymmErrors*)dir_postfit_two_mumu_zh->Get("data");

  THStack* stack_two_mumu_zh_pf = new THStack( "two_mumu_zh_postfit" , "two_mumu_zh combine-postfit" );
  create_stack( stack_two_mumu_zh_pf, light_two_mumu_zh_pf, heavy_two_mumu_zh_pf, other_two_mumu_zh_pf);
  create_ratio_plot(data_two_mumu_zh_pf, stack_two_mumu_zh_pf, bkg_total_two_mumu_zh_pf,
     Form("two_mumu_zh_postfit_%s",vrName.c_str()), light_two_mumu_zh_pf, heavy_two_mumu_zh_pf, other_two_mumu_zh_pf);


  //---------------------------
  //pre-fit
  //---------------------------
  std::cout << "=================================================" << std::endl;
  std::cout << "Z pre-fit yiedls: " << light_two_mumu_zh->GetBinContent(1)
  << " " << light_two_mumu_zh->GetBinContent(2)
  << " " << light_two_mumu_zh->GetBinContent(3) << std::endl;

  std::cout << "tt pre-fit yiedls: " << heavy_two_mumu_zh->GetBinContent(1)
  << " " << heavy_two_mumu_zh->GetBinContent(2)
  << " " << heavy_two_mumu_zh->GetBinContent(3) << std::endl;


  std::cout << "=================================================" << std::endl; 
  //----------------------------
  //post-fit
  //----------------------------
  std::cout << "Z postfit yiedls: " << light_two_mumu_zh_pf->GetBinContent(1)
  << " " << light_two_mumu_zh_pf->GetBinContent(2)
  << " " << light_two_mumu_zh_pf->GetBinContent(3) << std::endl;

  std::cout << "tt postfit yiedls: " << heavy_two_mumu_zh_pf->GetBinContent(1)
  << " " << heavy_two_mumu_zh_pf->GetBinContent(2)
  << " " << heavy_two_mumu_zh_pf->GetBinContent(3) << std::endl;
  /*
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

  stack_two_ee_dy->Draw();
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
  */
  return 0;
}
