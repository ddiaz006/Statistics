#include <iostream>
#include <string>
//ROOT
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
//LOCAL
#include "CommandLineInput.hh"


double _tag_var_sys[]   = {0.01,-0.02,-0.04};
double _ipsig_sys[]     = {0.01,-0.02,-0.04};
double _ta_sys[]        = {0.01,-0.02,-0.04};
double _amax_sys[]      = {0.01,-0.02,-0.04};
double _egs_sys[]       = {0.007,0.007,0.007};
double _mes_sys[]       = {0.007,0.007,0.007};
double _jes_sys[]       = {0.02,0.02,0.1};

int main( int argc, char* argv[] )
{

  // As usual, load the combine library to get access to the RooParametricHist
  //gSystem->Load("libHiggsAnalysisCombinedLimit.so");

  gROOT->Reset();
  //------------------------------
  // Getting Input ROOT File
  //------------------------------
  std::string  input_file = ParseCommandLine( argc, argv, "-input_file=" );
  if (  input_file == "" )
  {
    std::cerr << "[ERROR]: please provide a signal model --input_file=<input_file>" << std::endl;
    return -1;
  }

  //------------------------------
  // Getting OUT ROOT File Name
  //------------------------------
  std::string  output_file = ParseCommandLine( argc, argv, "-output_file=" );
  if (  output_file == "" )
  {
    std::cerr << "[ERROR]: please provide a signal model --output_file=<output_file>" << std::endl;
    return -1;
  }

  //---------------------------
  //open input file in read mode
  //---------------------------
  TFile* fin  = new TFile(input_file.c_str(), "READ");
  //get histograms
  TH1F* light = (TH1F*)fin->Get("light");
  TH1F* heavy = (TH1F*)fin->Get("heavy");
  TH1F* other = (TH1F*)fin->Get("other");
  //TH1F* data  = (TH1F*)fin->Get("Data");
  //TH1F* signal = (TH1F*)fin->Get("Sig_MS55ct100");

  //-----------------
  //TagVar Systematic
  //-----------------
  //light
  TH1F* light_TagVarsUp   = new TH1F(*light);
  TH1F* light_TagVarsDown = new TH1F(*light);
  //heavy
  TH1F* heavy_TagVarsUp   = new TH1F(*heavy);
  TH1F* heavy_TagVarsDown = new TH1F(*heavy);
  //other
  TH1F* other_TagVarsUp   = new TH1F(*other);
  TH1F* other_TagVarsDown = new TH1F(*other);

  //-----------------
  //AMax Systematic
  //-----------------
  //light
  TH1F* light_AMaxUp   = new TH1F(*light);
  TH1F* light_AMaxDown = new TH1F(*light);
  //heavy
  TH1F* heavy_AMaxUp   = new TH1F(*heavy);
  TH1F* heavy_AMaxDown = new TH1F(*heavy);
  //other
  TH1F* other_AMaxUp   = new TH1F(*other);
  TH1F* other_AMaxDown = new TH1F(*other);

  //-----------------
  //IPSIG Systematic
  //-----------------
  //light
  TH1F* light_IPSigUp   = new TH1F(*light);
  TH1F* light_IPSigDown = new TH1F(*light);
  //heavy
  TH1F* heavy_IPSigUp   = new TH1F(*heavy);
  TH1F* heavy_IPSigDown = new TH1F(*heavy);
  //other
  TH1F* other_IPSigUp   = new TH1F(*other);
  TH1F* other_IPSigDown = new TH1F(*other);

  //-----------------
  //TA Systematic
  //-----------------
  //light
  TH1F* light_TAUp   = new TH1F(*light);
  TH1F* light_TADown = new TH1F(*light);
  //heavy
  TH1F* heavy_TAUp   = new TH1F(*heavy);
  TH1F* heavy_TADown = new TH1F(*heavy);
  //other
  TH1F* other_TAUp   = new TH1F(*other);
  TH1F* other_TADown = new TH1F(*other);

  //-----------------
  // EGS Systematic
  //-----------------
  //light
  TH1F* light_EGSUp   = new TH1F(*light);
  TH1F* light_EGSDown = new TH1F(*light);
  //heavy
  TH1F* heavy_EGSUp   = new TH1F(*heavy);
  TH1F* heavy_EGSDown = new TH1F(*heavy);
  //other
  TH1F* other_EGSUp   = new TH1F(*other);
  TH1F* other_EGSDown = new TH1F(*other);

  //-----------------
  //MES Systematic
  //-----------------
  //light
  TH1F* light_MESUp   = new TH1F(*light);
  TH1F* light_MESDown = new TH1F(*light);
  //heavy
  TH1F* heavy_MESUp   = new TH1F(*heavy);
  TH1F* heavy_MESDown = new TH1F(*heavy);
  //other
  TH1F* other_MESUp   = new TH1F(*other);
  TH1F* other_MESDown = new TH1F(*other);

  //-----------------
  //JES Systematic
  //-----------------
  //light
  TH1F* light_JESUp   = new TH1F(*light);
  TH1F* light_JESDown = new TH1F(*light);
  //heavy
  TH1F* heavy_JESUp   = new TH1F(*heavy);
  TH1F* heavy_JESDown = new TH1F(*heavy);
  //other
  TH1F* other_JESUp   = new TH1F(*other);
  TH1F* other_JESDown = new TH1F(*other);

  //set systematic uncertainty already symmetric
  for ( int i = 1; i <= 3; i++ )
  {
    //light
    light_TagVarsUp->SetBinContent(i, light_TagVarsUp->GetBinContent(i)*(1.0+_tag_var_sys[i-1]));
    light_TagVarsDown->SetBinContent(i, light_TagVarsDown->GetBinContent(i)*(1.0-_tag_var_sys[i-1]));
    //heavy
    heavy_TagVarsUp->SetBinContent(i, heavy_TagVarsUp->GetBinContent(i)*(1.0+_tag_var_sys[i-1]));
    heavy_TagVarsDown->SetBinContent(i, heavy_TagVarsDown->GetBinContent(i)*(1.0-_tag_var_sys[i-1]));
    //other
    other_TagVarsUp->SetBinContent(i, other_TagVarsUp->GetBinContent(i)*(1.0+_tag_var_sys[i-1]));
    other_TagVarsDown->SetBinContent(i, other_TagVarsDown->GetBinContent(i)*(1.0-_tag_var_sys[i-1]));

    //----------------------
    //AMAX
    //----------------------
    //light
    light_AMaxUp->SetBinContent(i, light_AMaxUp->GetBinContent(i)*(1.0+_amax_sys[i-1]));
    light_AMaxDown->SetBinContent(i, light_AMaxDown->GetBinContent(i)*(1.0-_amax_sys[i-1]));
    //heavy
    heavy_AMaxUp->SetBinContent(i, heavy_AMaxUp->GetBinContent(i)*(1.0+_amax_sys[i-1]));
    heavy_AMaxDown->SetBinContent(i, heavy_AMaxDown->GetBinContent(i)*(1.0-_amax_sys[i-1]));
    //other
    other_AMaxUp->SetBinContent(i, other_AMaxUp->GetBinContent(i)*(1.0+_amax_sys[i-1]));
    other_AMaxDown->SetBinContent(i, other_AMaxDown->GetBinContent(i)*(1.0-_amax_sys[i-1]));

    //----------------------
    //IPSIG
    //----------------------
    //light
    light_IPSigUp->SetBinContent(i, light_IPSigUp->GetBinContent(i)*(1.0+_ipsig_sys[i-1]));
    light_IPSigDown->SetBinContent(i, light_IPSigDown->GetBinContent(i)*(1.0-_ipsig_sys[i-1]));
    //heavy
    heavy_IPSigUp->SetBinContent(i, heavy_IPSigUp->GetBinContent(i)*(1.0+_ipsig_sys[i-1]));
    heavy_IPSigDown->SetBinContent(i, heavy_IPSigDown->GetBinContent(i)*(1.0-_ipsig_sys[i-1]));
    //other
    other_IPSigUp->SetBinContent(i, other_IPSigUp->GetBinContent(i)*(1.0+_ipsig_sys[i-1]));
    other_IPSigDown->SetBinContent(i, other_IPSigDown->GetBinContent(i)*(1.0-_ipsig_sys[i-1]));

    //----------------------
    //TA
    //----------------------
    //light
    light_TAUp->SetBinContent(i, light_TAUp->GetBinContent(i)*(1.0+_ta_sys[i-1]));
    light_TADown->SetBinContent(i, light_TADown->GetBinContent(i)*(1.0-_ta_sys[i-1]));
    //heavy
    heavy_TAUp->SetBinContent(i, heavy_TAUp->GetBinContent(i)*(1.0+_ta_sys[i-1]));
    heavy_TADown->SetBinContent(i, heavy_TADown->GetBinContent(i)*(1.0-_ta_sys[i-1]));
    //other
    other_TAUp->SetBinContent(i, other_TAUp->GetBinContent(i)*(1.0+_ta_sys[i-1]));
    other_TADown->SetBinContent(i, other_TADown->GetBinContent(i)*(1.0-_ta_sys[i-1]));

    //----------------------
    //EGS
    //----------------------
    //light
    light_EGSUp->SetBinContent(i, light_EGSUp->GetBinContent(i)*(1.0+_egs_sys[i-1]));
    light_EGSDown->SetBinContent(i, light_EGSDown->GetBinContent(i)*(1.0-_egs_sys[i-1]));
    //heavy
    heavy_EGSUp->SetBinContent(i, heavy_EGSUp->GetBinContent(i)*(1.0+_egs_sys[i-1]));
    heavy_EGSDown->SetBinContent(i, heavy_EGSDown->GetBinContent(i)*(1.0-_egs_sys[i-1]));
    //other
    other_EGSUp->SetBinContent(i, other_EGSUp->GetBinContent(i)*(1.0+_egs_sys[i-1]));
    other_EGSDown->SetBinContent(i, other_EGSDown->GetBinContent(i)*(1.0-_egs_sys[i-1]));

    //----------------------
    //MES
    //----------------------
    //light
    light_MESUp->SetBinContent(i, light_MESUp->GetBinContent(i)*(1.0+_mes_sys[i-1]));
    light_MESDown->SetBinContent(i, light_MESDown->GetBinContent(i)*(1.0-_mes_sys[i-1]));
    //heavy
    heavy_MESUp->SetBinContent(i, heavy_MESUp->GetBinContent(i)*(1.0+_mes_sys[i-1]));
    heavy_MESDown->SetBinContent(i, heavy_MESDown->GetBinContent(i)*(1.0-_mes_sys[i-1]));
    //other
    other_MESUp->SetBinContent(i, other_MESUp->GetBinContent(i)*(1.0+_mes_sys[i-1]));
    other_MESDown->SetBinContent(i, other_MESDown->GetBinContent(i)*(1.0-_mes_sys[i-1]));

    //----------------------
    //JES
    //----------------------
    //light
    light_JESUp->SetBinContent(i, light_JESUp->GetBinContent(i)*(1.0+_jes_sys[i-1]));
    light_JESDown->SetBinContent(i, light_JESDown->GetBinContent(i)*(1.0-_jes_sys[i-1]));
    //heavy
    heavy_JESUp->SetBinContent(i, heavy_JESUp->GetBinContent(i)*(1.0+_jes_sys[i-1]));
    heavy_JESDown->SetBinContent(i, heavy_JESDown->GetBinContent(i)*(1.0-_jes_sys[i-1]));
    //other
    other_JESUp->SetBinContent(i, other_JESUp->GetBinContent(i)*(1.0+_jes_sys[i-1]));
    other_JESDown->SetBinContent(i, other_JESDown->GetBinContent(i)*(1.0-_jes_sys[i-1]));
  }

  //signal_string = signal_model.c_str();
  TFile* fout  = new TFile(output_file.c_str(), "RECREATE");
  //light->Write();
  //heavy->Write();
  //other->Write();
  //MakeMockupSystematics
  light_TagVarsUp->Write("light_TagVarsUp");
  light_TagVarsDown->Write("light_TagVarsDown");
  heavy_TagVarsUp->Write("heavy_TagVarsUp");
  heavy_TagVarsDown->Write("heavy_TagVarsDown");
  other_TagVarsUp->Write("other_TagVarsUp");
  other_TagVarsDown->Write("other_TagVarsDown");
  //AMAX
  light_AMaxUp->Write("light_AMaxUp");
  light_AMaxDown->Write("light_AMaxDown");
  heavy_AMaxUp->Write("heavy_AMaxUp");
  heavy_AMaxDown->Write("heavy_AMaxDown");
  other_AMaxUp->Write("other_AMaxUp");
  other_AMaxDown->Write("other_AMaxDown");
  //IPSIG
  light_IPSigUp->Write("light_IPSigUp");
  light_IPSigDown->Write("light_IPSigDown");
  heavy_IPSigUp->Write("heavy_IPSigUp");
  heavy_IPSigDown->Write("heavy_IPSigDown");
  other_IPSigUp->Write("other_IPSigUp");
  other_IPSigDown->Write("other_IPSigDown");
  //TA
  light_TAUp->Write("light_TAUp");
  light_TADown->Write("light_TADown");
  heavy_TAUp->Write("heavy_TAUp");
  heavy_TADown->Write("heavy_TADown");
  other_TAUp->Write("other_TAUp");
  other_TADown->Write("other_TADown");
  //EGS
  light_EGSUp->Write("light_EGSUp");
  light_EGSDown->Write("light_EGSDown");
  heavy_EGSUp->Write("heavy_EGSUp");
  heavy_EGSDown->Write("heavy_EGSDown");
  other_EGSUp->Write("other_EGSUp");
  other_EGSDown->Write("other_EGSDown");
  //MES
  light_MESUp->Write("light_MESUp");
  light_MESDown->Write("light_MESDown");
  heavy_MESUp->Write("heavy_MESUp");
  heavy_MESDown->Write("heavy_MESDown");
  other_MESUp->Write("other_MESUp");
  other_MESDown->Write("other_MESDown");
  //JES
  light_JESUp->Write("light_JESUp");
  light_JESDown->Write("light_JESDown");
  heavy_JESUp->Write("heavy_JESUp");
  heavy_JESDown->Write("heavy_JESDown");
  other_JESUp->Write("other_JESUp");
  other_JESDown->Write("other_JESDown");

  //signal->Write();
  fout->Close();

  return 0;
}
