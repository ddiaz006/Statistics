#ifndef HELPER_FUNCTIONS_HH
#define HELPER_FUNCTIONS_HH
#include <iostream>
#include <THStack.h>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLatex.h>

bool create_stack(THStack* stack, TH1F* light, TH1F* heavy, TH1F* other);
bool create_ratio_plot(TGraphAsymmErrors* data, THStack* stack, TH1F* total_bkg,
  TString plot_name, TH1F* light, TH1F* heavy, TH1F* other);
bool create_tf_plot(TH1F* h_to, TH1F* h_from, TString plot_name);
bool AddCMS( TCanvas* C );


#endif
