#include "CommonHeader.h"

////////////////////////////////////////////////////////////////////////////////
// Function for Double Template Param
Double_t funcCorrBackFitting(Double_t *x,  Double_t *par){
  Double_t xx = x[0];
  return (Double_t) par[0]*hBackStackup->GetBinContent(hBackStackup->FindBin(xx));
}

const Int_t ndrawpoints      = 1.e5;


////////////////////////////////////////////////////////////////////////////////
// Function to add similar looking background histos together for increase in
// statistics. This is needed for better Chi2 fitting since errors are too big
// otherwise.
// i is the number of the bin on which you add the background from the bins j to
// (including) k, but excluding i if i is between j and k
void BackgroundAdding(int i = 1, int j = 2, int k = 9, TFile* file = NULL, TString PicFormat = "png"){
  TH1D* hBack = NULL;
  TH1D* (aBackStackup[k-j]);
  TH1D* aPeak2[k-j];
  TH1D* aBack[k-j];


  /////////////////////////////////////////////////////////////////////////////
  // setting up the canvas to draw on. Will later be changed for the chi2 pic
  TCanvas *c1 = new TCanvas("c1","",1200,1000);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.09);
  c1->SetRightMargin(0.02);
  c1->SetLeftMargin(0.09);
  c1->SetTicky();
  c1->SetTickx();
  c1->SetLogz(1);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);

  hBack               = (TH1D*) file->Get(Form("hCorrBack_bin%02i",i));
  aBack[0] = (TH1D*) hBack->Clone(Form("aBack%02d", 0));

  int backscale = 0;

  for(int m = j; m <= k; m++){
    if(m != i){

      hBackStackup = (TH1D*) file->Get(Form("hCorrBack_bin%02i",m));
      aBackStackup[m] = (TH1D*) file->Get(Form("hCorrBack_bin%02i",m));
      ////////////////////////////////////////////////////////////////////////////
      // fit function
      TF1* fit = new TF1("fit", &funcCorrBackFitting, 0.0 ,0.4, 1);
      fit->SetNpx(ndrawpoints);
      fit->SetNumberFitPoints(numberbins);
      fit->SetLineColor(kRed);
      fit->SetLineWidth(3);
      hBack->Fit("fit", "QM0P","", 0.1, 0.3);
      hBackStackup->Scale(fit->GetParameter(0));
      aBackStackup[m]->Scale(fit->GetParameter(0));

      if(m == j){
        aBack[m] = (TH1D*) aBack[0]->Clone(Form("aBack%02d", m));
      }
      else{
        aBack[m] = (TH1D*) aBack[m-1]->Clone(Form("aBack%02d", m));
      }
      aBack[m]->Add(aBackStackup[m], 1);
      backscale++;
      delete fit;
    }
    else{
      continue;
    }
  }
  hBack->Clear();
  aBack[k]->Scale(1./(Double_t) backscale);
  hBack = (TH1D*) aBack[k]->Clone("hBack");


  SetHistoStandardSettings(hBack);

  c1->cd();
  hBack->Draw("");


  c1->Update();
  c1->SaveAs(Form("BackGroundFitting/PiledUpBack.png"));
  c1->Clear();

  delete c1;
}
/******************************************************************************/



void BackgroundTestFit(int i, TFile* file, TString PicFormat){
  // setting j and k right depending on the binning
  int j = i;
  int k = i;
  if(i == 1){
    j = 2;
  }
  else if(i == 2){
    j = 1;
  }
  else{
    j = i-2;
  }
  if(i > numberbins-4){
    exit(42);
  }
  else if(i == numberbins-5){
    k = i;
  }
  else if(i == numberbins-6){
    k = i+1;
  }
  else{
    k = i+2;
  }
  /////////////////////////////////////////////////////////////////////////////
  // setting up the canvas to draw on. Will later be changed for the chi2 pic
  TCanvas *c1 = new TCanvas("c1","",1200,1000);
  c1->cd();
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.09);
  c1->SetRightMargin(0.02);
  c1->SetLeftMargin(0.09);
  c1->SetTicky();
  c1->SetTickx();
  c1->SetLogz(1);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);

  TH1D* hBack           = NULL;

  hBack                 = (TH1D*) file->Get(Form("hCorrBack_bin%02i",i));

  for(int m = k; m >= j; m--){
    if(m != i){

      hBackStackup = (TH1D*) file->Get(Form("hCorrBack_bin%02i",m));

      // rebinning of the mid histo and the two to the left from rebin 4 to
      // rebin 8
      if(GetNBinningFromHistogram(hBackStackup) < GetNBinningFromHistogram(hBack)){
        hBack->Rebin(2);
        if(m < i){
          hBackStackup->Rebin(2);
        }
      }

      // rebin of the one histogram to the right if needed
      if(m == i+1){
        if(GetNBinningFromHistogram(hBackStackup) > GetNBinningFromHistogram(hBack)){
          hBackStackup->Rebin(2);
        }
      }

      ////////////////////////////////////////////////////////////////////////////
      // fit function
      TF1* fit = new TF1("fit", &funcCorrBackFitting, 0.0 ,0.4, 1);
      fit->SetNpx(ndrawpoints);
      fit->SetNumberFitPoints(numberbins);
      fit->SetLineColor(kRed);
      fit->SetLineWidth(3);
      hBack->Fit("fit", "QM0P","", 0.1, 0.3);
      hBackStackup->Scale(fit->GetParameter(0));

      SetHistoStandardSettings(hBackStackup);
      hBackStackup->SetMarkerStyle(1);
      if(m == i-2){
        hBackStackup->SetLineColor(kRed);
        hBackStackup->SetMarkerColor(kRed);
      }

      if(m == i-1){
        hBackStackup->SetLineColor(kTeal-7);
        hBackStackup->SetMarkerColor(kTeal-7);
      }

      if(m == i+1){
        hBackStackup->SetLineColor(kMagenta+2);
        hBackStackup->SetMarkerColor(kMagenta+2);
      }

      if(m == i+2){
        hBackStackup->SetLineColor(kBlue+2);
        hBackStackup->SetMarkerColor(kBlue+2);
      }

      if(m == k){
        c1->cd();
        c1->Clear();
        SetHistoStandardSettings(hBack);
        hBack->SetMarkerStyle(1);
        hBack->Draw("AXIS");
        hBack->DrawClone("SAME");
        hBack->DrawClone("SAME HIST");
      }
      hBackStackup->DrawClone("SAME");
      hBackStackup->DrawClone("SAME HIST");


      c1->Update();

      delete fit;
    }
    else{
      continue;
    }
  }
  c1->SaveAs(Form("BackGroundFitting/Fit%02d." + PicFormat,i));
  c1->Clear();
  delete hBack;
  delete c1;

}

////////////////////////////////////////////////////////////////////////////////
// Start of the Main
void BackGroundFitting(TString PicFormat = "png"){
  TString sPath = gDirectory->GetPath();            // retrieve neutral path

  ////////////////////////////////////////////////////////////////////////////
  // open MC histo path
  TFile* IterTemp = SafelyOpenRootfile("IterTemp.root");
  if (IterTemp->IsOpen() ) printf("IterTemp opened successfully\n");
  gDirectory->Cd(sPath.Data());

  BackgroundAdding(1, 2, 9, IterTemp, PicFormat);
  for (int n = 1; n < numberbins; n++) {
    BackgroundTestFit(n, IterTemp, PicFormat);
  }


  IterTemp->Close();

}
