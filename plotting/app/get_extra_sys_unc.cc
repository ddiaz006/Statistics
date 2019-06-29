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

bool AddCMS( TCanvas* C );

int main( int argc, char** argv )
{

  //-----------------
  //Input File List
  //-----------------
  std::string inputList = ParseCommandLine( argc, argv, "-inputList=" );
  if (  inputList == "" )
  {
    std::cerr << "[ERROR]: please provide an inputList. Use --inputList=" << std::endl;
    return -1;
  }
  //std::ifstream ifs ( inputList.c_str(), std::ifstream::in );//input file list

  //----------------------
  //validation region name
  //----------------------
  //std::string vrName = ParseCommandLine( argc, argv, "-vrName=" );
  //if (  vrName == "" )
  //{
  //  std::cerr << "[ERROR]: please provide validation region name. Use --vrName=" << std::endl;
  //  return -1;
  //}


  if( _info )
  {
    std::cout << "[INFO]: input file: " << inputList << std::endl;
    //std::cout << "[INFO]: validation region name: " << vrName << std::endl;
  }

  //--------------------------------
  //Array with bin extra uncertainty
  //--------------------------------
  double err_ee[]   = {0.0, 0.0, 0.0};
  double err_mumu[] = {0.0, 0.0, 0.0};

  TFile* fin;
  std::ifstream ifs ( inputList.c_str(), std::ifstream::in);
  if (ifs.is_open())
  {
    while (ifs.good())
    {
      std::string current_fname;
      ifs >> current_fname;
      if (ifs.eof()) break;
      //-----------------
      //Getting ROOT file
      //-----------------
      fin = new TFile( current_fname.c_str(), "READ");
      //----------------------------------------------------------
      //GETTING POST-FIT RESULTS FROM COMBINE
      //----------------------------------------------------------
      //---------------------
      //TwoEleZH
      //---------------------
      TDirectory* dir_postfit_two_ee_zh = (TDirectory*)(((TDirectory*)fin->Get("shapes_fit_b"))->Get("TwoEleZH"));
      TH1F* bkg_total_two_ee_zh_pf = (TH1F*)dir_postfit_two_ee_zh->Get("total_background");
      TGraphAsymmErrors* data_two_ee_zh_pf = (TGraphAsymmErrors*)dir_postfit_two_ee_zh->Get("data");
      for( unsigned int i = 1; i <= bkg_total_two_ee_zh_pf->GetNbinsX(); i++ )
      {
        //data treatment
        double y_data, x_data;
        data_two_ee_zh_pf->GetPoint(i-1, x_data, y_data);
        //relative difference
        double delta_y     = y_data - bkg_total_two_ee_zh_pf->GetBinContent(i);
        double delta_y_rel = delta_y/bkg_total_two_ee_zh_pf->GetBinContent(i);
        double y_unc       = bkg_total_two_ee_zh_pf->GetBinError(i);

        if( fabs( delta_y ) > y_unc && y_data != 0.0 )
        {
          std::cout << "ee->"<< i-1 << " " << delta_y_rel<< std::endl;
          err_ee[i-1] += pow(delta_y_rel,2.0);
        }
      }


      //---------------------
      //TwoMuZH
      //---------------------
      TDirectory* dir_postfit_two_mumu_zh = (TDirectory*)(((TDirectory*)fin->Get("shapes_fit_b"))->Get("TwoMuZH"));
      TH1F* bkg_total_two_mumu_zh_pf = (TH1F*)dir_postfit_two_mumu_zh->Get("total_background");
      TGraphAsymmErrors* data_two_mumu_zh_pf = (TGraphAsymmErrors*)dir_postfit_two_mumu_zh->Get("data");
      for( unsigned int i = 1; i <= bkg_total_two_mumu_zh_pf->GetNbinsX(); i++ )
      {
        //data treatment
        double y_data, x_data;
        data_two_mumu_zh_pf->GetPoint(i-1, x_data, y_data);
        //relative difference
        double delta_y     = y_data - bkg_total_two_mumu_zh_pf->GetBinContent(i);
        double delta_y_rel = delta_y/bkg_total_two_mumu_zh_pf->GetBinContent(i);
        double y_unc       = bkg_total_two_mumu_zh_pf->GetBinError(i);
        //std::cout << i-1 << " " << delta_y_rel<< std::endl;
        if( fabs( delta_y ) > y_unc && y_data != 0.0 )
        {
          std::cout << "mumu->"<< i-1 << " " << delta_y_rel<< std::endl;
          err_mumu[i-1] += pow(delta_y_rel,2.0);
        }
      }


      //Close current root file
      fin->Close();
    }
  }
  else {
    // show message:
    std::cerr << "[ERROR]: Error opening file: " << inputList << std::endl;
    std::cout << "Exiting program!!!" << std::endl;
    return -1;
  }

  //--------------------
  //ee extra uncertainty
  //--------------------
  std::cout << "ee extra_unc bins (0,1,2) is: (" << sqrt(err_ee[0]) << ", "
  << sqrt(err_ee[1])  << ", " << sqrt(err_ee[2]) << ")"<< std::endl;
  //----------------------
  //mumu extra uncertainty
  //----------------------
  std::cout << "mumu extra_unc bins (0,1,2) is: (" << sqrt(err_mumu[0]) << ", "
  << sqrt(err_mumu[1])  << ", " << sqrt(err_mumu[2]) << ")"<< std::endl;

  return 0;
}
