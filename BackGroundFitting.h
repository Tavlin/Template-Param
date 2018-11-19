#ifndef BackGroundFitting_H
#define BackGroundFitting_H
#include "chi2test.h"

////////////////////////////////////////////////////////////////////////////////
// Function for Double Template Param
Double_t funcCorrBackFitting(Double_t *x,  Double_t *par){
  Double_t xx = x[0];
  Double_t f = par[0]*hBackStackup->GetBinContent(hBackStackup->FindBin(xx));
  return f;
}

Double_t funcOpeningAngleCut(Double_t *x, Double_t *par){
  Double_t xx = x[0];
  return xx/(sqrt(2.*(1.-cos(0.017/2.))));
}

const Int_t ndrawpoints      =  1.e5;
Double_t tempup              = -1.e10;
Double_t templow             =  1.e10;

/**
 * Function to create the needed corr. bkg. templates for the NN and the 3 to 8
 * method.
 * The templates will be saved ion a File called "CorrBkgFileNoRebin.root"
 * Name of the templates: hCorrBkgNoRebinBin%02d !
 */
void CorrBkgCreation(void){

  TString safePath = gDirectory->GetPath();            // retrieve neutral path
  TFile* MCFile         = NULL;
  TFile* DataFile       = NULL;
  TH1D* hSignalTemplate = NULL;
  TH1D* hInvMassData    = NULL;
  TH1D* hCorrBkg        = NULL;

  TFile* CorrBkgFile    = new TFile("CorrBkgFileNoRebin.root", "RECREATE");

  gDirectory->Cd(safePath.Data());                  // for saftey resetting path



  for (int k = 1; k < numberbins; k++) {

    /**
     * Open the file which contains the MC output of the framework's work so to
     * say.
     */
    MCFile              =
    SafelyOpenRootfile("/data4/mhemmer/Documents/BachelorArbeit/GammaCalo_503_JJ_2017_Rebin1/00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_MC_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");

    /**
     * Histogram from the MC simulation which contains the only the true Pi0s
     * coming from y y; y_conv y; and double y_conv
     */
    hSignalTemplate     = (TH1D*) MCFile->Get(Form("Mapping_TrueMeson_InvMass_in_Pt_Bin%02d",k));

    /**
     * Open the file which contains the data output of the framework's work so to
     * say.
     */
    // DataFile            =
    // SafelyOpenRootfile("/data4/mhemmer/Documents/BachelorArbeit/GammaCalo_503_data_2017_Rebin1/00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_MC_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");

    /**
     * Histogram from the data which was created through the normal analysis
     * method: same events - scaled mixed events.
     */
    hInvMassData        = (TH1D*) MCFile->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",k));

    hCorrBkg            = (TH1D*) hInvMassData->Clone("hCorrBkg");
    hCorrBkg->          Add(hSignalTemplate, -1);

    for(int i = 1; i < hCorrBkg->GetNbinsX(); i++){
      hCorrBkg->SetBinError(i, sqrt(pow(hInvMassData->GetBinError(i),2.) - pow(hSignalTemplate->GetBinError(i),2.) ));
    }

    gDirectory = CorrBkgFile;
    // hCorrBkg->          GetXaxis()->SetRangeUser(0.0, 0.3);
    hCorrBkg->          Write(Form("hCorrBkgNoRebinBin%02d", k));

    gDirectory->Cd(safePath.Data());                  // for saftey resetting path

    MCFile->            Close();
    // DataFile->          Close();
  }
  CorrBkgFile->Close();
}


void GetLowerBounds(void){

  TFile* OACFile        = SafelyOpenRootfile("./OAC_ToyMCMerged.root");
  TH2D* hMinv_pT_ratio  = (TH2D*)  OACFile->Get("hMinv_pT_ratio");
  for (int i = 1; i < numberbins; i++) {
    TH1D* hMinv_pT_ratio_projectionX = NULL;
    hMinv_pT_ratio_projectionX = (TH1D*) hMinv_pT_ratio->ProjectionX(Form("hMinv_pT_ratio_projectionX[%02d,%02d]", i, i+1), i, i+1);
    Int_t lowerbinnumber = hMinv_pT_ratio_projectionX->FindFirstBinAbove(8.) + 5;

    std::cout << "lower bin = " << hMinv_pT_ratio_projectionX->GetXaxis()->GetBinCenter(lowerbinnumber)  << '\n';
  }

  OACFile->Close();
}


/**
 * Function to add neighboring  background histos together for increase in
 * statistics. This is needed for better Chi2 fitting since errors are too big
 * otherwise.
 * @param i         current pT-Bin
 * @param b         number of Neighbours
 * @param file      File where the normal corr. bkg. histos are.
 * @param PicFormat Formattype of the pictures. (.eps or .png)
 */
TH1D* BackgroundAdding(int i){

  TString sPath = gDirectory->GetPath();

  TFile* CorrBkgFile    = SafelyOpenRootfile("./CorrBkgFileNoRebin.root");

  // setting j and k right depending on the binning
  const int b = 6;
  int j = i-(b/2);
  int k = i+(b/2);
  if(k >= numberbins){
    k = numberbins-1;
  }
  if(j <= 0){
    j = 1;
  }
  int scalefactor = k-j;
  if(i == 1 || i == numberbins-1){
    scalefactor = 3;
  }
  TH1D* (aBackStackup[k-j]);
  TFile* BackFile;
  TH1D* hBack           = NULL;
  TH1D* hRatio          = NULL;
  TList *list = new TList;

  TF1* fPol0            = NULL;


  hBack                 = (TH1D*) CorrBkgFile->Get(Form("hCorrBkgNoRebinBin%02d",i));
  hBack->Rebin(fBinsPi013TeVEMCPtRebin[i]);

  // comment *out* if you want fit without original bin!
  // list->Add(hBack);
  // comment *in* if you want fit without original bin!
  scalefactor -= 1;

  for(int m = k; m >= j; m--){
    if(m != i){
      aBackStackup[k-m] = (TH1D*) CorrBkgFile->Get(Form("hCorrBkgNoRebinBin%02d",m));
      aBackStackup[k-m]->Rebin(fBinsPi013TeVEMCPtRebin[i]);

      ////////////////////////////////////////////////////////////////////////////
      // fit function
      if(m != k){
        hBackStackup->Reset();
      }
      hBackStackup = NULL;
      hBackStackup = (TH1D*) (aBackStackup[k-m])->Clone("hBackStackup");
      Int_t NFitPoints = hBackStackup->FindBin(0.3) - hBackStackup->FindBin(0.1);
      std::cout << "Before Fit" << '\n';
      TF1* fit = new TF1("fit", &funcCorrBackFitting, 0.0 ,0.3, 1);
      fit->SetNpx(ndrawpoints);
      fit->SetParameter(0, 1.);
      fit->SetNumberFitPoints(NFitPoints);
      fit->SetLineColor(kRed);
      fit->SetLineWidth(3);
      hBack->Fit("fit", "QM0P","", 0.1, 0.3);
      std::cout << "After Fit" << '\n';
      aBackStackup[k-m]->Scale(fit->GetParameter(0));

      aBackStackup[k-m]->SetMarkerStyle(1);
      list->Add(aBackStackup[k-m]);
      delete fit;
    }
    else{
      continue;
    }
  }

  TH1D* hUpperBkg     = (TH1D*) hBack->Clone("hUpperBkg");
  TH1D* hLowerBkg     = (TH1D*) hBack->Clone("hLowerBkg");
  TH1D* hPilledUpBack = (TH1D*) hBack->Clone("hPilledUpBack");

  Int_t Nbins = hBack->GetEntries();
  for(int o = 0; o < Nbins; o++){
    for(int l = k; l >= j; l--){
      if(l != i){
        if(tempup < aBackStackup[k-l]->GetBinContent(o)){
          tempup = aBackStackup[k-l]->GetBinContent(o);
        }
        if(templow > aBackStackup[k-l]->GetBinContent(o)){
          templow = aBackStackup[k-l]->GetBinContent(o);
        }
      }
    }
    hUpperBkg->SetBinContent(o,tempup);
    hUpperBkg->SetBinError(o, 0.0);
    hLowerBkg->SetBinContent(o, templow);
    hLowerBkg->SetBinError(o, 0.0);
    tempup = -1.e10;
    templow = 1.e10;
  }

  hPilledUpBack->Reset();
  hPilledUpBack->Merge(list);
  hPilledUpBack->Scale(1./(Double_t)scalefactor);

  hRatio      = (TH1D*) hBack->Clone();
  hRatio->Divide(hBack, hPilledUpBack, 1, 1);

  fPol0        = new TF1(Form("fPol0_bin%02d",i), "[0]", 0.0, 0.3);
  fPol0->SetParLimits(0, 0.0, 5.0);
  fPol0->SetParLimits(0, 0., 100.);

  for (int t = 1; t < hRatio->GetNbinsX(); t++) {
    if(hRatio->GetBinError(t) > 40 || fabs(hRatio->GetBinContent(t)) > 40){
      hRatio->SetBinContent(t, 0.0);
      hRatio->SetBinError(t, 0.0);
    }
  }
  hRatio->Fit(fPol0,"QM0P", "",  0.1, 0.2);


  if(i == 1){
    BackFile      = new TFile("BackFileNN.root", "RECREATE");
  }
  else{
    BackFile      = new TFile("BackFileNN.root", "UPDATE");
  }

  // hPilledUpBack->Write(Form("hPilledUpBack_Bin%02d_with%02d_bins", i, b));
  // hUpperBkg->Write(Form("hUpperBkg_Bin%02d_with%02d_bins", i, b));
  // hLowerBkg->Write(Form("hLowerBkg_Bin%02d_with%02d_bins", i, b));
  hPilledUpBack->Write(Form("hPilledUpBack_Bin%02d", i));
  hRatio->Write(Form("hRatio_Bin%02d", i));
  fPol0->Write(Form("fPol0_Bin%02d", i));

  // aBackStackup[k-m]= NULL;
  hRatio = NULL;
  delete fPol0;


  hUpperBkg = NULL;
  hLowerBkg = NULL;
  BackFile->Close();

  gDirectory->Cd(sPath.Data());


  delete list;
  hBack = NULL;
  std::cout << "" << '\n';

  return hPilledUpBack;

}

/**
 * Function to add  background histos from 3rd to 8th pT-bin together for
 * increase in statistics. This is needed for better Chi2 fitting since errors
 * are too big otherwise.
 * @param i    current pT-bin
 */
TH1D* BackGround3to8(int i){

  TFile* OACFile        = SafelyOpenRootfile("./OAC_ToyMCMerged.root");
  TH2D* hMinv_pT_ratio  = (TH2D*)  OACFile->Get("hMinv_pT_ratio");

  TFile* CorrBkgFile    = SafelyOpenRootfile("./CorrBkgFileNoRebin.root");

  // setting j and k right depending on the binning
  int j = 3;
  int k = 8;
  int scalefactor;
  if(i > 2 && i <= 8){
    scalefactor = 5;
  }
  else{
    scalefactor = 6;
  }
  TH1D* (aBackStackup[k-j]);
  TFile* BackFile;
  TH1D* hBack           = NULL;
  TH1D* hRatio          = NULL;
  TF1*  fPol0           = NULL;
  TList *list = new TList;
  hBack                 = (TH1D*) CorrBkgFile->Get(Form("hCorrBkgNoRebinBin%02d",i));
  hBack->Rebin(fBinsPi013TeVEMCPtRebin[i]);
  hMinv_pT_ratio->RebinX(fBinsPi013TeVEMCPtRebin[i]);
  // needs rescaling since it is 8 histos merged and also for the rebinnig
  hMinv_pT_ratio->Scale(1./(8.*fBinsPi013TeVEMCPtRebin[i]));

  for(int m = k; m >= j; m--){
    if(m != i){

      aBackStackup[k-m] = (TH1D*) CorrBkgFile->Get(Form("hCorrBkgNoRebinBin%02d",m));
      aBackStackup[k-m]->Rebin(fBinsPi013TeVEMCPtRebin[i]);
      for (int w = 1; w < aBackStackup[k-m]->FindBin(0.3); w++) {
        if(fabs(hMinv_pT_ratio->GetBinContent(w, m+1)) >= 1.e-2){
          aBackStackup[k-m]->SetBinContent(w, aBackStackup[k-m]->GetBinContent(w) * 1./hMinv_pT_ratio->GetBinContent(w, m+1));
          aBackStackup[k-m]->SetBinError(w, aBackStackup[k-m]->GetBinError(w) * 1./hMinv_pT_ratio->GetBinContent(w, m+1));
        }
        else{
          aBackStackup[k-m]->SetBinContent(w, 0.0);
          aBackStackup[k-m]->SetBinError(w, 0.0);
        }
      }
      list->Add(aBackStackup[k-m]);

    }
    else{
      continue;
    }
  }

  TH1D* hPilledUpBack = (TH1D*) hBack->Clone("hPilledUpBack");
  hPilledUpBack->Reset();
  hPilledUpBack->Merge(list);
  hPilledUpBack->Scale(1./(Double_t)scalefactor);
  hPilledUpBack->GetXaxis()->SetRangeUser(0.0, 0.3);

  Double_t sum = 0;
  Double_t OAC_scaling = 0;
  for (int b = 1; b <= (800/fBinsPi013TeVEMCPtRebin[i])*3./8.+1; b++) {
    sum = 0;
    // for (int c = 3; c <= 8; c++) {
    //   if(i == c){
    //     continue;
    //   }
    //   else{
    //     sum += hMinv_pT_ratio->GetBinContent(b, c+1);
    //   }
    // }
    // std::cout << "BinContent = " << hMinv_pT_ratio->GetBinContent(b, i+1) << '\n';
    if(hMinv_pT_ratio->GetBinContent(b, i+1) != 0){
      OAC_scaling = 0;
      OAC_scaling = hMinv_pT_ratio->GetBinContent(b, i+1);  //sum;
      hPilledUpBack->SetBinContent(b,hPilledUpBack->GetBinContent(b)*OAC_scaling);
      hPilledUpBack->SetBinError(b,hPilledUpBack->GetBinError(b)*OAC_scaling);
    }
    else{
      hPilledUpBack->SetBinContent(b, 0.0);
      hPilledUpBack->SetBinError(b, 0.0);
    }
    // std::cout << "BinContent = " << hPilledUpBack->GetBinContent(b) << '\n';
  }

  hRatio      = (TH1D*) hBack->Clone();
  hRatio->Divide(hBack, hPilledUpBack, 1, 1);

  fPol0        = new TF1(Form("fPol0_bin%02d",i), "[0]", 0.0, 0.3);
  fPol0->SetParLimits(0, 0.0, 5.0);

  for (int t = 1; t < hRatio->GetNbinsX(); t++) {
    if(hRatio->GetBinError(t) > 40 || fabs(hRatio->GetBinContent(t)) > 40){
      hRatio->SetBinContent(t, 0.0);
      hRatio->SetBinError(t, 0.0);
    }
  }
  hRatio->Fit(fPol0,"QM0P", "",  0.1, 0.2);



  TString sPath = gDirectory->GetPath();

  if(i == 1){
    BackFile      = new TFile("BackFile3to8.root", "RECREATE");
  }
  else{
    BackFile      = new TFile("BackFile3to8.root", "UPDATE");
  }

  hRatio->Write(Form("hRatio_Bin%02d", i));
  fPol0->Write(Form("fPol0_Bin%02d", i));
  BackFile->Close();

  gDirectory->Cd(sPath.Data());


  delete list;
  delete hBack;
  delete fPol0;
  OACFile->Close();
  CorrBkgFile->Close();

  return hPilledUpBack;
}
#endif
