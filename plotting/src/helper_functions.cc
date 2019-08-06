#include <iostream>
#include <TLegend.h>
#include <helper_functions.hh>

//------------------------
//CMS COSMETICS DEFINITION
//------------------------
//const float lumi = 5;
//Axis
//const float axisTitleSize = 0.06;
//const float axisTitleOffset = .8;

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


bool create_stack(THStack* stack, TH1F* light, TH1F* heavy, TH1F* other)
{
  //cosmetics for light
  light->SetFillColor(kGreen-3);
  light->SetLineColor(kGreen-3);
  //cosmetics for heavy
  heavy->SetFillColor(kMagenta-4);
  heavy->SetLineColor(kMagenta-4);
  //cosmetics for other
  other->SetFillColor(kAzure-4);
  other->SetLineColor(kAzure-4);

  stack->Add( other, "histo" );
  stack->Add( heavy, "histo" );
  stack->Add( light, "histo" );

  //stack cosmetics
  stack->SetTitle("");
  return true;
};

bool create_ratio_plot(TGraphAsymmErrors* data, THStack* stack, TH1F* total_bkg,
   TString plot_name, TH1F* light, TH1F* heavy, TH1F* other)
{
  //cosmetics for data
  data->SetLineColor(kBlack);
  data->SetMarkerColor(kBlack);
  data->SetLineWidth(2);
  data->SetMarkerSize(1.5);
  data->SetMarkerStyle(20);


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
  c->cd();

  TPad *pad1 = new TPad("pad1","pad1", 0.0, 0.3, 1., 1.);
  pad1->SetBottomMargin(0);
  pad1->SetRightMargin( rightMargin );
  pad1->SetLeftMargin( leftMargin );
  pad1->Draw();
  pad1->cd();

  stack->Draw();
  stack->SetMinimum(1e-3);
  //Get maximum for data and compare to max in bkg histo
  double x_b0,y_b0;
  if(data->GetPoint(0,x_b0,y_b0) > total_bkg->GetBinContent(1))
  {
    stack->SetMaximum(1000.*y_b0);
  }
  else
  {
    stack->SetMaximum(1000.*total_bkg->GetBinContent(1));
  }

  data->Draw("P");
  pad1->SetLogy();

  //-----------------
  //Legend
  //------------------
  TLegend* leg = new TLegend( 0.73, 0.65, 0.93, 0.88, NULL, "brNDC" );
  leg->SetBorderSize(0);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetTextSize(0.04);
  leg->AddEntry( data, "  data", "lep" );
  leg->AddEntry( light, "  Z/#gamma^{*}", "f" );
  leg->AddEntry( heavy, "  t#bar{t}+t", "f" );
  leg->AddEntry( other, "  other", "f" );

  leg->Draw();


  //------------------
  //Ratio Plot
  //------------------
  c->cd();
  TPad *pad2 = new TPad("pad2","pad2", .0, 0.0, 1., 0.29);
  pad2->SetTopMargin(0.04);
  pad2->SetTopMargin(0.008);
  pad2->SetBottomMargin(0.4);
  pad2->SetRightMargin( rightMargin );
  pad2->SetLeftMargin( leftMargin );
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  TH1F* data_hist = new TH1F("data_hist", "data_hist", 3, -0.5,2.5);
  for(int i = 0; i < data_hist->GetNbinsX(); i++)
  {
    double x,y;
    data->GetPoint(i,x,y);
    data_hist->SetBinContent(i+1, y);
  }
  TH1F* ratio = new TH1F( *data_hist );
  ratio->Divide( total_bkg );
  //Set Uncertainty to only the poisson, remove the one from bkg
  for ( int i  = 1; i <= ratio->GetNbinsX(); i++ )
  {
    double rel_unc = data_hist->GetBinError(i)/data_hist->GetBinContent(i);
    ratio->SetBinError(i,rel_unc);
  }
  ratio->SetMarkerStyle( 20 );
  ratio->SetMarkerSize( 1.0 );
  ratio->GetXaxis()->SetTitleSize( axisTitleSizeRatioX );
  ratio->GetXaxis()->SetLabelSize( axisLabelSizeRatioX );
  ratio->GetXaxis()->SetTitleOffset( axisTitleOffsetRatioX );
  ratio->GetYaxis()->SetTitleSize( axisTitleSizeRatioY );
  ratio->GetYaxis()->SetLabelSize( axisLabelSizeRatioY );
  ratio->GetYaxis()->SetTitleOffset( axisTitleOffsetRatioY );
  ratio->SetMarkerColor( kBlack );
  ratio->SetLineColor( kBlack );
  ratio->SetLineWidth( 2 );
  ratio->GetYaxis()->SetRangeUser( 0.0, 3.0 );
  ratio->SetTitle("");
  ratio->GetYaxis()->SetTitle("data / mc");
  ratio->GetXaxis()->SetTitle("nJet_{tags}");
  ratio->GetYaxis()->CenterTitle( true );
  ratio->GetYaxis()->SetNdivisions( 10, false );
  ratio->SetStats( 0 );
  ratio->Draw("E1");



  TH1F* ratio2 = new TH1F( *data_hist );
  TH1F* ratio3 = new TH1F( *data_hist );
  //Set Uncertainty to only the poisson, remove the one from bkg
  for ( int i  = 1; i <= ratio2->GetNbinsX(); i++ )
  {
    ratio2->SetBinContent(i,1.0);
    double rel_unc = total_bkg->GetBinError(i)/total_bkg->GetBinContent(i);
    ratio2->SetBinError(i,rel_unc);

    ratio3->SetBinContent(i,1.0);
    //ratio3->SetBinError(i,0.0);
  }
  //ratio2->SetMarkerStyle( 20 );
  //ratio2->SetMarkerSize( 1.5 );
  ratio2->GetXaxis()->SetTitleSize( axisTitleSizeRatioX );
  ratio2->GetXaxis()->SetLabelSize( axisLabelSizeRatioX );
  ratio2->GetXaxis()->SetTitleOffset( axisTitleOffsetRatioX );
  ratio2->GetYaxis()->SetTitleSize( axisTitleSizeRatioY );
  ratio2->GetYaxis()->SetLabelSize( axisLabelSizeRatioY );
  ratio2->GetYaxis()->SetTitleOffset( axisTitleOffsetRatioY );
  ratio2->SetMarkerColor( kBlue );
  ratio2->SetLineColor( kBlue );
  ratio2->SetFillColorAlpha(kBlue-9, 0.35);
  //ratio2->SetFillColor( kBlue-9 );
  ratio2->SetLineWidth( 2 );
  ratio2->GetYaxis()->SetRangeUser( 0.0, 3.0 );
  ratio2->SetTitle("");
  ratio2->GetYaxis()->SetTitle("data / mc");
  ratio2->GetXaxis()->SetTitle("nJet_{tags}");
  ratio2->GetYaxis()->CenterTitle( true );
  ratio2->GetYaxis()->SetNdivisions( 10, false );
  ratio2->SetStats( 0 );

  //ratio3->SetMarkerStyle( 20 );
  //ratio3->SetMarkerSize( 1.5 );
  ratio3->GetXaxis()->SetTitleSize( axisTitleSizeRatioX );
  ratio3->GetXaxis()->SetLabelSize( axisLabelSizeRatioX );
  ratio3->GetXaxis()->SetTitleOffset( axisTitleOffsetRatioX );
  ratio3->GetYaxis()->SetTitleSize( axisTitleSizeRatioY );
  ratio3->GetYaxis()->SetLabelSize( axisLabelSizeRatioY );
  ratio3->GetYaxis()->SetTitleOffset( axisTitleOffsetRatioY );
  ratio3->SetMarkerColor( kBlue );
  ratio3->SetLineColor( kBlue );
  //ratio3->SetFillColorAlpha(kBlue-9, 0.35);
  //ratio2->SetFillColor( kBlue-9 );
  ratio3->SetLineWidth( 2 );
  ratio3->GetYaxis()->SetRangeUser( 0.0, 3.0 );
  ratio3->SetTitle("");
  ratio3->GetYaxis()->SetTitle("data / mc");
  ratio3->GetXaxis()->SetTitle("nJet_{tags}");
  ratio3->GetYaxis()->CenterTitle( true );
  ratio3->GetYaxis()->SetNdivisions( 10, false );
  ratio3->SetStats( 0 );


  ratio3->Draw("HIST");
  ratio2->Draw("E2+SAME");
  ratio->Draw("E1+SAME");

  c->SaveAs(plot_name+".pdf");
  c->SaveAs(plot_name+".C");

  for ( int i  = 1; i <= ratio2->GetNbinsX(); i++ )
  {
    std::cout << "ratio bin " << i << " :" << ratio->GetBinContent(i) << std::endl;
  }
  return true;
};

bool create_tf_plot(TH1F* h_to, TH1F* h_from, TString plot_name)
{
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
  c->cd();
  //-----------------
  //Legend
  //------------------
  TLegend* leg = new TLegend( 0.73, 0.65, 0.93, 0.88, NULL, "brNDC" );
  leg->SetBorderSize(0);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetTextSize(0.04);
  //leg->AddEntry( data, "  data", "lep" );
  //leg->AddEntry( light, "  Z/#gamma^{*}", "f" );
  //leg->AddEntry( heavy, "  t#bar{t}+t", "f" );
  //leg->AddEntry( other, "  other", "f" );

  leg->Draw();


  //------------------
  //Ratio Plot
  //------------------
  c->cd();


  TH1F* ratio = new TH1F( *h_to );
  ratio->Divide( h_from );
  //Set Uncertainty to only the poisson, remove the one from bkg
  ratio->SetMarkerStyle( 20 );
  ratio->SetMarkerSize( 1.0 );
  ratio->GetXaxis()->SetTitleSize( axisTitleSizeRatioX );
  ratio->GetXaxis()->SetLabelSize( axisLabelSizeRatioX );
  ratio->GetXaxis()->SetTitleOffset( axisTitleOffsetRatioX );
  ratio->GetYaxis()->SetTitleSize( axisTitleSizeRatioY );
  ratio->GetYaxis()->SetLabelSize( axisLabelSizeRatioY );
  ratio->GetYaxis()->SetTitleOffset( axisTitleOffsetRatioY );
  ratio->SetMarkerColor( kBlack );
  ratio->SetLineColor( kBlack );
  ratio->SetLineWidth( 2 );
  ratio->GetYaxis()->SetRangeUser( 0.0, 3.0 );
  ratio->SetTitle("");
  ratio->GetYaxis()->SetTitle("data / mc");
  ratio->GetXaxis()->SetTitle("nJet_{tags}");
  ratio->GetYaxis()->CenterTitle( true );
  ratio->GetYaxis()->SetNdivisions( 10, false );
  ratio->SetStats( 0 );
  ratio->Draw("E1");


  //ratio->Draw("E1+SAME");

  c->SaveAs(plot_name+".pdf");
  c->SaveAs(plot_name+".C");

  return true;
};

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
