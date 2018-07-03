#include "CommonHeader.h"


void MixedBGComp(void){

  const Int_t nbins = 45;
  const Int_t ndrawpoints = 1.e5;

  //////////////////////////////////////////////////////////////////////////////
  // setting up the canvas to draw on. Will later be changed for the chi2 pic
  TCanvas *c1 = new TCanvas("c1","",1200,1100);
  c1->cd();
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.09);
  c1->SetRightMargin(0.005);
  c1->SetLeftMargin(0.09);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);

  //////////////////////////////////////////////////////////////////////////////
  // setting up the 2 Huistograms to compare chi2 from the to fit methods as
  // well as peak factor comp. between pol 1 and double temp fit and the ratio
  // of BG. scaling factor and the Peak scaling factor

  TH1F* hBGComp = new TH1F("hBGComp", "", 26, fBinsPi013TeVEMCPt);
  SetHistoStandardSettings(hBGComp);


  //////////////////////////////////////////////////////////////////////////////
  // going over all pt bins despite first one, which is some framework bs.
  for (int k = 1; k < 26; k++) {

    ////////////////////////////////////////////////////////////////////////////
    // open MC histo path
    TFile* MCFile = new TFile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_MC_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root", "READ");
    if (MCFile->IsOpen() ) printf("MCFile opened successfully\n");

    ////////////////////////////////////////////////////////////////////////////
    // retrieve MC histograms
    MCBG = (TH1F*) MCFile->Get(Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02i",k));

    ////////////////////////////////////////////////////////////////////////////
    // open Data path
    TFile* DataFile = new TFile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root", "READ");
    if (DataFile->IsOpen() ) printf("DataFile opened successfully\n");

    ////////////////////////////////////////////////////////////////////////////
    // retrieve data histogram
    DataBG = (TH1F*) DataFile->Get(Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02i",k));

    ////////////////////////////////////////////////////////////////////////////
    // making things look good
    // change y title offset to fit and look nicely
    DataBG->GetYaxis()->SetTitleOffset(1.2);
    DataBG->GetYaxis()->SetLabelOffset(0.006);
    DataBG->SetTitleSize(0.03, "xy");
    DataBG->SetLabelSize(0.03, "xy");
    DataBG->SetYTitle("d#it{N}/d#it{M}_{#gamma#gamma} (#it{c}^{2}/GeV)");
    DataBG->SetMarkerStyle(20);
    DataBG->SetMarkerSize(1.5);
    DataBG->SetTitle("");
    MCBG->SetLineColor(kTeal-7);
    MCBG->SetMarkerColor(kTeal-7);
    MCBG->SetMarkerStyle(34);
    MCBG->SetMarkerSize(1.5);

    DataBG->Divide(MCBG);
    DataBG->GetYaxis()->SetRangeUser(0,20);
    c1->Update();
    DataBG->Draw("");


    c1->SaveAs(Form("MixedBGComp/DataDivMC%02i.png",k));
    c1->Clear();

    DataBG->Delete();
    MCBG->Delete();

  }
  delete(c1);
}
