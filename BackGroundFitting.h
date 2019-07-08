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
void CorrBkgCreation(const char* strFrameworkOutputNoRebinning){

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
    SafelyOpenRootfile(strFrameworkOutputNoRebinning);

    /**
     * Histogram from the MC simulation which contains the only the true Pi0s
     * coming from y y; y_conv y; and double y_conv
     */
    hSignalTemplate     = (TH1D*) MCFile->Get(Form("Mapping_TrueFullMeson_InvMass_in_Pt_Bin%02d",k));

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
    Int_t lowerbinnumber = hMinv_pT_ratio_projectionX->FindFirstBinAbove(0.0);

    std::cout << "bin " << i << " minv_min = " << hMinv_pT_ratio_projectionX->GetXaxis()->GetBinCenter(lowerbinnumber)  << '\n';
  }

  OACFile->Close();
}


/**
 * Function to add neighboring  background histos together for increase in
 * statistics. This is needed for better Chi2 fitting since errors are too big
 * otherwise.
 * @param i         current pT-Bin
 */
TH1D* BackgroundAdding(int i){

  TString sPath = gDirectory->GetPath();

  TFile* CorrBkgFile    = SafelyOpenRootfile("./CorrBkgFileNoRebin.root");

  // setting j and k right depending on the binning
  const double b = 3.;
  int j = i;
  int k = i+1;
  while(fBinsPi013TeVEMCPt[k]-fBinsPi013TeVEMCPt[j] < b){
    if(j < 2){
      j = 1;
      k++;
    }
    else if(k >= numberbins-2){
      k = numberbins-2;
      j--;
    }
    else{
      k++;
      j--;
    }
  }
  if(i >= 38){
    k = 38;
    j = 37;
  }
  std::cout << fBinsPi013TeVEMCPt[k]-fBinsPi013TeVEMCPt[j] << '\n';
  std::cout << "k and j are " << k <<  " and " << j << '\n';
  int scalefactor = k-j;

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
      TF1* fit = new TF1(Form("fit%d", m), &funcCorrBackFitting, 0.0 ,0.3, 1);
      fit->SetNpx(ndrawpoints);
      fit->SetParameter(0, 1.);
      fit->SetNumberFitPoints(NFitPoints);
      fit->SetLineColor(kRed);
      fit->SetLineWidth(3);
      hBack->Fit(Form("fit%d", m), "QM0","", 0.1, 0.3);
      aBackStackup[k-m]->Scale(fit->GetParameter(0));

      aBackStackup[k-m]->SetMarkerStyle(1);
      list->Add(aBackStackup[k-m]);
      delete fit;
    }
    else{
      continue;
    }
  }
  // std::cout << "outside the loop" << '\n';

  // std::cout << "before upper and lower ranges" << '\n';
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

  fPol0        = new TF1(Form("fPol0_bin%02d",i), "pol0", 0.0, 0.3);
  fPol0->SetParameter(0, 1.0);
  // fPol0->SetParLimits(0, -1., 100.);

  for (int t = 1; t < hRatio->GetNbinsX(); t++) {
    if(hRatio->GetBinError(t) > 40 || fabs(hRatio->GetBinContent(t)) > 40){
      hRatio->SetBinContent(t, 0.0);
      hRatio->SetBinError(t, 0.0);
    }
  }
  hRatio->Fit(Form("fPol0_bin%02d",i), "QM0", "",  0.08, 0.2);

  // std::cout << "before making outputfile" << '\n';
  if(i == 1){
    BackFile      = new TFile("BackFileNN.root", "RECREATE");
  }
  else{
    BackFile      = new TFile("BackFileNN.root", "UPDATE");
  }
  // std::cout << "after making outputfile" << '\n';

  // hPilledUpBack->Write(Form("hPilledUpBack_Bin%02d_with%02d_bins", i, b));
  // hUpperBkg->Write(Form("hUpperBkg_Bin%02d_with%02d_bins", i, b));
  // hLowerBkg->Write(Form("hLowerBkg_Bin%02d_with%02d_bins", i, b));
  hBack->Write(Form("hBack_Bin%02d", i));
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
TH1D* BackGround3to8(int i, const int TEMPLATEMETHOD){

  TFile* OACFile        = SafelyOpenRootfile("./OAC_ToyMCMerged.root");
  std::cout << "OACFile safetly opend" << '\n';
  TH2D* hMinv_pT_ratio  = (TH2D*)  OACFile->Get("hMinv_pT_ratio");

  TFile* CorrBkgFile    = SafelyOpenRootfile("./CorrBkgFileNoRebin.root");
  std::cout << "CorrBkgFile safetly opend" << '\n';

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
  switch (TEMPLATEMETHOD) {
    case 3:
      hBack->Rebin(fBinsPi013TeVEMCPtSmallerRebin[i]);
      hMinv_pT_ratio->RebinX(fBinsPi013TeVEMCPtSmallerRebin[i]);
      // needs rescaling since it is 8 histos merged and also for the rebinnig
      hMinv_pT_ratio->Scale(1./(8.0*fBinsPi013TeVEMCPtSmallerRebin[i]));
      break;
    case 4:
      hBack->Rebin(fBinsPi013TeVEMCPtHigherRebin[i]);
      hMinv_pT_ratio->RebinX(fBinsPi013TeVEMCPtHigherRebin[i]);
      // needs rescaling since it is 8 histos merged and also for the rebinnig
      hMinv_pT_ratio->Scale(1./(8.0*fBinsPi013TeVEMCPtHigherRebin[i]));
      break;
    default:
      hBack->Rebin(fBinsPi013TeVEMCPtRebin[i]);
      hMinv_pT_ratio->RebinX(fBinsPi013TeVEMCPtRebin[i]);
      // needs rescaling since it is 8 histos merged and also for the rebinnig
      hMinv_pT_ratio->Scale(1./(8.0*fBinsPi013TeVEMCPtRebin[i]));
      break;
  }

  for(int m = k; m >= j; m--){
    if(m != i){

      aBackStackup[k-m] = (TH1D*) CorrBkgFile->Get(Form("hCorrBkgNoRebinBin%02d",m));
      switch (TEMPLATEMETHOD) {
        case 3:
          aBackStackup[k-m]->Rebin(fBinsPi013TeVEMCPtSmallerRebin[i]);
          break;
        case 4:
          aBackStackup[k-m]->Rebin(fBinsPi013TeVEMCPtHigherRebin[i]);
          break;
        default:
          aBackStackup[k-m]->Rebin(fBinsPi013TeVEMCPtRebin[i]);
          break;
      }
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
  switch (TEMPLATEMETHOD) {
    case 3:
      for (int b = 1; b <= (150/fBinsPi013TeVEMCPtSmallerRebin[i]); b++) { //NO IDEA ABOUT THAT +1!!!
        sum = 0;
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
      }
      break;
    case 4:
      for (int b = 1; b <= (150/fBinsPi013TeVEMCPtHigherRebin[i]); b++) { //NO IDEA ABOUT THAT +1!!!
        sum = 0;
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
      }
      break;
    default:
      for (int b = 1; b <= (150/fBinsPi013TeVEMCPtRebin[i]); b++) { //NO IDEA ABOUT THAT +1!!!
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
      break;
  }
  std::cout << "number bins = " << hMinv_pT_ratio->GetNbinsX() << '\n';
  hRatio      = (TH1D*) hBack->Clone();
  hRatio->Divide(hBack, hPilledUpBack, 1, 1);

  fPol0        = new TF1(Form("fPol0_bin%02d",i), "[0]", 0.0, 0.3);
  fPol0->SetParLimits(0, 0.0, 5.0);

  std::cout << "Smashing where? 1" << '\n';
  for (int t = 1; t < hRatio->GetNbinsX(); t++) {
    if(hRatio->GetBinError(t) > 40 || fabs(hRatio->GetBinContent(t)) > 40){
      hRatio->SetBinContent(t, 0.0);
      hRatio->SetBinError(t, 0.0);
    }
  }
  hRatio->Fit(fPol0,"QM0", "",  0.1, 0.2);



  TString sPath = gDirectory->GetPath();

  switch (TEMPLATEMETHOD) {
    case 3:
      if(i == 1){
        BackFile      = new TFile("BackFile3to8_SmallerRebin.root", "RECREATE");
      }
      else{
        BackFile      = new TFile("BackFile3to8_SmallerRebin.root", "UPDATE");
      }
      break;
    case 4:
      if(i == 1){
        BackFile      = new TFile("BackFile3to8_HigherRebin.root", "RECREATE");
      }
      else{
        BackFile      = new TFile("BackFile3to8_HigherRebin.root", "UPDATE");
      }
      break;
    default:
      if(i == 1){
        BackFile      = new TFile("BackFile3to8.root", "RECREATE");
      }
      else{
        BackFile      = new TFile("BackFile3to8.root", "UPDATE");
      }
      break;
  }
  std::cout << "Smashing where? 2" << '\n';

  hRatio->Write(Form("hRatio_Bin%02d", i));
  fPol0->Write(Form("fPol0_Bin%02d", i));
  BackFile->Close();
  std::cout << "Smashing where? 3" << '\n';

  gDirectory->Cd(sPath.Data());


  delete list;
  delete hBack;
  delete fPol0;
  OACFile->Close();
  CorrBkgFile->Close();

  std::cout << "Smashing where? 4" << '\n';

  return hPilledUpBack;
}
#endif
