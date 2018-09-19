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
void BackgroundAdding(int i, TFile* file, TString PicFormat){
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
  if(i >= numberbins-4){
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
  int scalefactor = 0;
  if(i == 1){
    scalefactor = k-j+1;
  }
  else{
    scalefactor = k-j;
  }
  TH1D* (aBackStackup[k-j]);
  TFile* BackFile;
  TH1D* hBack           = NULL;
  TList *list = new TList;

  hBack                 = (TH1D*) file->Get(Form("hCorrBack_bin%02d",i));
  list->Add(hBack);

  for(int m = k; m >= j; m--){
    if(m != i){
      std::cout << "m = " << m << '\n';

      aBackStackup[k-m] = (TH1D*) file->Get(Form("hCorrBack_bin%02d",m));

      // rebinning of the mid histo and the two to the left from rebin 4 to
      // rebin 8
      std::cout << "Bins not mid = " << aBackStackup[k-m]->fNcells  << '\n';
      std::cout << "Bins mid = " << hBack->fNcells  << '\n';
      if(aBackStackup[k-m]->fNcells < hBack->fNcells){
        std::cout << "inside upper rebinning" << '\n';
        hBack->Rebin(2);
      }

      // rebin of the one histogram to the right if needed
      if(aBackStackup[k-m]->fNcells > hBack->fNcells){
        std::cout << "inside lower rebinning" << '\n';
        (aBackStackup[k-m])->Rebin(2);
      }


      ////////////////////////////////////////////////////////////////////////////
      // fit function
      if(m != k){
        hBackStackup->Reset();
      }
      hBackStackup = NULL;
      hBackStackup = (TH1D*) (aBackStackup[k-m])->Clone("hBackStackup");
      TF1* fit = new TF1("fit", &funcCorrBackFitting, 0.0 ,0.4, 1);
      fit->SetNpx(ndrawpoints);
      fit->SetNumberFitPoints(numberbins);
      fit->SetLineColor(kRed);
      fit->SetLineWidth(3);
      hBack->Fit("fit", "QM0P","", 0.1, 0.2);
      aBackStackup[k-m]->Scale(fit->GetParameter(0));

      SetHistoStandardSettings(aBackStackup[k-m]);
      aBackStackup[k-m]->SetMarkerStyle(1);
      if(m == i-2){
        aBackStackup[k-m]->SetLineColor(kRed);
        aBackStackup[k-m]->SetMarkerColor(kRed);
      }

      if(m == i-1){
        aBackStackup[k-m]->SetLineColor(kTeal-7);
        aBackStackup[k-m]->SetMarkerColor(kTeal-7);
      }

      if(m == i+1){
        aBackStackup[k-m]->SetLineColor(kMagenta+2);
        aBackStackup[k-m]->SetMarkerColor(kMagenta+2);
      }

      if(m == i+2){
        aBackStackup[k-m]->SetLineColor(kBlue+2);
        aBackStackup[k-m]->SetMarkerColor(kBlue+2);
      }

      list->Add(aBackStackup[k-m]);
      delete fit;
    }
    else{
      continue;
    }
  }

  TH1D* hPilledUpBack = (TH1D*) hBack->Clone("hPilledUpBack");
  hPilledUpBack->Reset();
  hPilledUpBack->Merge(list);
  hPilledUpBack->Scale(1./(Double_t)scalefactor);

  TString sPath = gDirectory->GetPath();

  if(i == 1){
    BackFile      = new TFile("BackFile.root", "RECREATE");
  }
  else{
    BackFile      = new TFile("BackFile.root", "UPDATE");
  }

  hPilledUpBack->Write(Form("hPilledUpBack_Bin%02d", i));

  BackFile->Close();

  gDirectory->Cd(sPath.Data());

  PicFormat = "eps";
  for(int m = k; m >= j; m--){
    aBackStackup[k-m]= NULL;
  }
  delete hPilledUpBack;
  delete list;
  delete hBack;
  std::cout << "" << '\n';

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
      TLegend* legBkg = new TLegend(0.6, 0.6, 0.9, 0.9);
      SetLegendSettigns(legBkg);
      hBackStackup = (TH1D*) file->Get(Form("hCorrBack_bin%02i",m));
      legBkg->AddEntry(hBack, "actual bin Bkg");
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
      hBack->Fit("fit", "QM0P","", 0.1, 0.2);
      hBackStackup->Scale(fit->GetParameter(0));

      SetHistoStandardSettings(hBackStackup);
      hBackStackup->SetMarkerStyle(1);
      if(m == i-2){
        hBackStackup->SetLineColor(kRed);
        hBackStackup->SetMarkerColor(kRed);
        legBkg->AddEntry(hBackStackup, "bin-2 Bkg");
      }

      if(m == i-1){
        hBackStackup->SetLineColor(kTeal-7);
        hBackStackup->SetMarkerColor(kTeal-7);
        legBkg->AddEntry(hBackStackup, "bin-1 Bkg");
      }

      if(m == i+1){
        hBackStackup->SetLineColor(kMagenta+2);
        hBackStackup->SetMarkerColor(kMagenta+2);
        legBkg->AddEntry(hBackStackup, "bin+1 Bkg");
      }

      if(m == i+2){
        hBackStackup->SetLineColor(kBlue+2);
        hBackStackup->SetMarkerColor(kBlue+2);
        legBkg->AddEntry(hBackStackup, "bin+2 Bkg");
      }

      if(m == k){
        c1->cd();
        c1->Clear();
        SetHistoStandardSettings(hBack);
        hBack->SetMarkerStyle(1);
        hBack->SetLineWidth(3);
        hBack->Draw("AXIS");
        // hBack->DrawClone("SAME");
        hBack->DrawClone("SAME HIST");
      }
      // hBackStackup->DrawClone("SAME");
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


  for (int n = 1; n < numberbins - 2; n++) {
    std::cout << "n = " << n << '\n';
    BackgroundTestFit(n, IterTemp, PicFormat);
    BackgroundAdding(n, IterTemp, PicFormat);
  }


  IterTemp->Close();

}
