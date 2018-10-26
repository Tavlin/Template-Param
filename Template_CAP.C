#include "CommonHeader.h"

/*******************************************************************************
Start of the main function

current_path is needed for the function to be calleable for variable datat sets

templatemethod == 1 uses the backgrund templates from the 3 to 8 method
               == 2 uses the  -||-                        Next Neighbours method

*******************************************************************************/
void Template_CAP(std::string current_path, int templatemethod){
  TString safePath = gDirectory->GetPath();            // retrieve neutral path

  /*****************************************************************************
  open ESD histo which is inside TLists inside a rootfile
  There get the Histo with Number of Events (minimum Bias)
  *****************************************************************************/
  TFile* ESDFile_MC               = NULL; // ESD File from the MC
  TFile* ESDFile_data             = NULL; // ESD File from the Data
  TFile* BkgFile                  = NULL; // Self makde lower stat. Bkg File
  TList* lGammaCalo_data          = NULL; // TLists inside the ESD File
  TList* lCutNumber_data          = NULL; // TLists inside the ESD File
  TList* lESD_data                = NULL; // Innerst TList for ESD histos from hInvMass_Data
  TList* lGammaCalo_MC            = NULL;
  TList* lCutNumber_MC            = NULL;
  TList* lMC_MC                   = NULL; // Innerst TList for MC histos from MC
  TList* lTrue_MC                 = NULL; // Innerst TList for True histos from MC
  TH1D* hNEvents                  = NULL; // histo containing number of Events
  TH1D* hMC_Pi0InAcc_Pt           = NULL; // acceptance histo
  TH2D* hTrueDoubleCounting_Pi0   = NULL; // 2D Histo including Doublecounting
  TH1D* hTrueDoubleCounting_Pi0_Pro;      // X-Projection of current
                                          // pT for Double Counting

  /*****************************************************************************
  Access the ESD File form the MC simulation for two Histograms:
  1st. the Doublecounting Histogram
  2nd. the MC Histogram of Pi0 in acceptance as a function of pT
  *****************************************************************************/
  ESDFile_MC    = SafelyOpenRootfile("./../Daten/" + current_path + ".root");
  if (ESDFile_MC->IsOpen() ) printf("ESDFile_MC opened successfully\n");

  lGammaCalo_MC          = (TList*) ESDFile_MC->Get("GammaCalo_503");
  lCutNumber_MC          = (TList*) lGammaCalo_MC->FindObject("Cut Number 00010113_1111112067032220000_01631031000000d0");
  lMC_MC                 = (TList*) lCutNumber_MC->FindObject("00010113_1111112067032220000_01631031000000d0 MC histograms");
  lTrue_MC               = (TList*) lCutNumber_MC->FindObject("00010113_1111112067032220000_01631031000000d0 True histograms");
  hMC_Pi0InAcc_Pt         = (TH1D*)  lMC_MC->FindObject("MC_Pi0InAcc_Pt");
  hTrueDoubleCounting_Pi0 = (TH2D*)  lTrue_MC->FindObject("ESD_TrueDoubleCountPi0_InvMass_Pt");
  ESDFile_MC->Close();

  /*****************************************************************************
  Access the ESD File form the data for one Histograms:
  1st. the NEvents Histogram to get the number of Events with MinBias Trigger
  *****************************************************************************/
  ESDFile_data    = SafelyOpenRootfile("./../Daten/GammaCalo-data_503.root");
  if (ESDFile_data->IsOpen() ) printf("ESDFile_data opened successfully\n");

  lGammaCalo_data        = (TList*) ESDFile_data->Get("GammaCalo_503");
  lCutNumber_data        = (TList*) lGammaCalo_data->FindObject("Cut Number 00010113_1111112067032220000_01631031000000d0");
  lESD_data              = (TList*) lCutNumber_data->FindObject("00010113_1111112067032220000_01631031000000d0 ESD histograms");
  hNEvents               = (TH1D*)  lESD_data->FindObject("NEvents");

  Double_t NEvents  = hNEvents->GetBinContent(1);   // retrieve NEents MinBias
  ESDFile_data->Close();

  gDirectory->Cd(safePath.Data());                     // for saftey resetting path
}
