#include "CommonHeader.h"
#include "TFractionFitter.h"


Double_t mc_full_func1(Double_t x){
  return mc_full_clone1->GetBinContent(mc_full->FindBin(x));
}

Double_t mc_full_func2(Double_t x){
  return mc_full_clone2->GetBinContent(mc_full->FindBin(x));
}

Double_t PeakAKorrBG(Double_t x){
  return korrBG_clone1->GetBinContent(korrBG->FindBin(x));
}

Double_t mc_full_func42(Double_t x){
  return mc_full_clone42->GetBinContent(mc_full->FindBin(x));
}

Double_t PeakAKorrBG42(Double_t x){
  return korrBG_clone42->GetBinContent(korrBG->FindBin(x));
}

void IterTemPlot(void){

  TString str;
  const Int_t nbins = 45;
  const Int_t ndrawpoints = 1.e5;
  const int n_iter = 4;
  TString doubletempstring = "Double template param.";
  TString pol1string = "Peak template + 1^{st} ord. pol. param.";

  //////////////////////////////////////////////////////////////////////////////
  // setting up the canvas to draw on. Will later be changed for the chi2 pic
  TCanvas *c1 = new TCanvas("c1","",1200,1000);
  c1->cd();
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.09);
  c1->SetRightMargin(0.02);
  c1->SetLeftMargin(0.09);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);

  TCanvas *canInvMass = new TCanvas("canInvMass","",1200,1200);
  TPad *pad1InvMass = new TPad("pad1InvMass","",0.0,0.33,1.0,1.0);
  pad1InvMass->SetTopMargin(0.05);
  pad1InvMass->SetLeftMargin(0.09);
  pad1InvMass->SetBottomMargin(0.0);
  pad1InvMass->SetRightMargin(0.02);
  TPad *pad2InvMass = new TPad("pad2InvMass","",0.0,0.0,1.0,0.33);
  pad2InvMass->SetTopMargin(0.0);
  pad2InvMass->SetLeftMargin(0.09);
  pad2InvMass->SetBottomMargin(0.3);
  pad2InvMass->SetRightMargin(0.02);

  TF1* fit_eq_double_temp = new TF1("fit_eq_double_temp", "PeakAKorrBG(x)*[1] + mc_full_func1(x)*[0]", 0.0,0.4);
  fit_eq_double_temp->SetNpx(ndrawpoints);
  fit_eq_double_temp->SetNumberFitPoints(nbins);
  fit_eq_double_temp->SetLineColor(kTeal-7);
  fit_eq_double_temp->SetLineWidth(4);


  TF1* fit_eq_1 = new TF1("fit_eq_1", "mc_full_func2(x)*[0]+[2]+x*[3]",0.0,0.4);
  fit_eq_1->SetNpx(ndrawpoints);
  fit_eq_1->SetNumberFitPoints(nbins);
  fit_eq_1->SetLineColor(kRed);
  fit_eq_1->SetLineWidth(4);

  //////////////////////////////////////////////////////////////////////////////
  // setting up the 2 Histograms to compare chi2 from the to fit methods as
  // well as peak factor comp. between pol 1 and double temp fit and the ratio
  // of BG. scaling factor and the Peak scaling factor

  TH1F* hChi2_dt
  TH1F* hChi2_pol1
  TH1F* hPeakRatio
  TH1F* hBGtoPeak
  TH1F* hRatioDoubleTemp;
  TH1F* hRatioPol1;
  TH1F* hData;
  TH1F* hDTPeak;
  TH1F* hDTBG;
  TH1F* hPol1Peak;
  TF1* fpol1;

  //////////////////////////////////////////////////////////////////////////////
  // going over all pt bins despite first one, which is some framework bs.
  for (int k = 1; k < 26; k++) {


    TFile* IterTemp = SafelyOpenRootfile("IterTemp.root");
    if (IterTemp->IsOpen() ) printf("IterTemp opened successfully\n");
    

  }
}
