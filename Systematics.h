#ifndef Systematics
#define Systematics

#include "CommonHeader.h"

void systematics(int templatemethod, TFile* OutputFile){

  TFile* InputFile      = NULL;
  TFile* NNMethodFile   = NULL;
  TFile* SingleBkgFile  = NULL;

  if(templatemethod == 2){
    InputFile      = SafelyOpenRootfile("OutputFileBetterBkgNN.root");
  }
  else if(templatemethod == 1){
    InputFile      = SafelyOpenRootfile("OutputFileBetterBkg3to8.root");
    NNMethodFile   = SafelyOpenRootfile("OutputFileBetterBkgNN.root");
    SingleBkgFile  = SafelyOpenRootfile("OutputFileNormal.root");
  }
  else if(templatemethod == 3){
    InputFile      = SafelyOpenRootfile("OutputFileBetterBkgPulse.root");
  }
  else if(templatemethod == 4){
    InputFile      = SafelyOpenRootfile("OutputFileNormal.root");
  }
  else if(templatemethod == 5){
    InputFile      = SafelyOpenRootfile("OutputFileOneTemplate.root");
  }
  else{
    std::cerr << "templatemethod not found!" << '\n';
    exit(1);
  }

  TH1D* hCorrYieldME                      = NULL;
  TH1D* hCorrYieldME_StatError            = NULL;
  TH1D* hCorrectedYieldNormEff_StatError  = NULL;
  TH1D* hCorrectedYieldNormEff            = NULL;
  TH1D* hCorrectedYieldNormEff_Ratio      = NULL;
  TH1D* hCorrYieldME_Ratio                = NULL;
  TH1D* hPi0_gen                          = NULL;
  TH1D* hCorrYield_HigherInt              = NULL;
  TH1D* hCorrYield_SmallInt               = NULL;
  TH1D* hCorrYield_IntSysError            = NULL;
  TH1D* hCorrYield_HigherFit              = NULL;
  TH1D* hCorrYield_SmallFit               = NULL;
  TH1D* hCorrYield_FitSysErrro            = NULL;
  TH1D* hCorrYield_LowerRebinning         = NULL;
  TH1D* hCorrYield_HigherRebinning        = NULL;
  TH1D* hCorrYield_RebinningSysError      = NULL;
  TH1D* hCorrYield_SysError               = NULL;
  TH1D* hCorrYield_RelativSyserror        = NULL;
  TH1D* hCorrYield_SyserrorRatio          = NULL;
  TH1D* hCorrYield_NNMethod               = NULL;
  TH1D* hCorrYield_SingleBkg              = NULL;
  TH1D* hCorrYield_CorrBkgSysError        = NULL;
  TFile* FStatUnc                         = NULL;

  TF1* fitBylikin13TeV = new TF1("fitBylikin13TeV", "[0]*exp(-(sqrt(x^(2)+0.135^(2))-0.135)/[1])+[2]/((1+x*x/[3])^([4]))", 1.4, 12.);
  fitBylikin13TeV->SetParameters(13, 0.1, 2, 0.7, 2.9);
  fitBylikin13TeV->SetLineWidth(3);
  fitBylikin13TeV->SetLineColor(kBlack);

  TF1* fitBylikin13TeV_3to8 = new TF1("fitBylikin13TeV_3to8", "[0]*exp(-(sqrt(x^(2)+0.135^(2))-0.135)/[1])+[2]/((1+x*x/[3])^([4]))", 1.4, 12.);
  fitBylikin13TeV_3to8->SetParameters(13, 0.1, 2, 0.7, 2.9);
  fitBylikin13TeV_3to8->SetLineWidth(3);
  fitBylikin13TeV_3to8->SetLineColor(kBlue+2);

  TF1 *ftsallis13TeV = new TF1("ftsallis13TeV", "[0]/(2*3.1415) * ( ([1]-1) * ([1]-2) )/( [1]*[2] *( [1]*[2] + [3] *([1]-2)) )* pow(( 1 + ( ( pow(([3]*[3]+x*x),0.5) -[3]) /( [1]*[2] ) ) ), -[1])", 1.4, 12.);
  ftsallis13TeV->SetParameter(0, 9.4);           // A
  ftsallis13TeV->SetParameter(1, 7.169);         // n
  ftsallis13TeV->SetParameter(2, 0.159);         // T
  ftsallis13TeV->FixParameter(3, 0.1349766);     // M
  ftsallis13TeV->SetLineWidth(3);
  ftsallis13TeV->SetLineColor(kBlack);


  std::vector<Double_t> vCountSys;
  std::vector<Double_t> vParamSys;
  std::vector<Double_t> vBGFitSys;
  std::vector<Double_t> vFinalSys;
  std::vector<Double_t> vCorrBkgSys;

  // TFile* JoshuasFile      = SafelyOpenRootfile("Pi0_data_GammaConvV1Correction_00010113_1111112067032220000_01631031000000d0.root");

  hCorrYieldME                  = (TH1D*) InputFile->Get("hYield_dt_chi2map_corrected");
  hPi0_gen                      = (TH1D*) InputFile->Get("MC_Meson_genPt");
  hCorrectedYieldNormEff        = (TH1D*) InputFile->Get("hCorrectedYieldNormEff");

  if(templatemethod ==1){

  hCorrYield_NNMethod           = (TH1D*) NNMethodFile->  Get("hYield_dt_chi2map_corrected");
  hCorrYield_SingleBkg          = (TH1D*) SingleBkgFile-> Get("hYield_dt_chi2map_corrected");
  }

  /**
   * CHANGED
   * for MC we want ratio to generated, else we want ratio to FW Yield
   */
  hCorrYieldME_Ratio            = (TH1D*) hCorrYieldME->Clone("hCorrYieldME_Ratio");
  hCorrYieldME_Ratio->          Divide(hCorrYieldME_Ratio,      hCorrectedYieldNormEff, 1, 1, "B");

  hCorrectedYieldNormEff_Ratio  = (TH1D*) hCorrectedYieldNormEff->Clone("hCorrectedYieldNormEff_Ratio");
  hCorrectedYieldNormEff_Ratio->Divide(hCorrectedYieldNormEff,  hPi0_gen, 1, 1, "B");


  hCorrYieldME_StatError            = (TH1D*) hCorrYieldME->Clone("hCorrYieldME_StatError");
  hCorrectedYieldNormEff_StatError  = (TH1D*) hCorrectedYieldNormEff->Clone("hCorrectedYieldNormEff_StatError");
  hCorrYieldME_StatError->SetXTitle(pt_str);
  hCorrYieldME_StatError->SetYTitle("relative stat. Unsicherheit (%)");
  hCorrectedYieldNormEff_StatError->SetXTitle(pt_str);
  hCorrectedYieldNormEff_StatError->SetYTitle("relative stat. Unsicherheit (%)");

  for(int k = 2; k < numberbins; k++){
    hCorrYieldME_StatError->SetBinContent(k, hCorrYieldME_StatError->GetBinError(k)/(Double_t)hCorrYieldME_StatError->GetBinContent(k)*100.);
    hCorrYieldME_StatError->SetBinError(k,0.);

    hCorrectedYieldNormEff_StatError->SetBinContent(k, hCorrectedYieldNormEff_StatError->GetBinError(k)/(Double_t)hCorrectedYieldNormEff_StatError->GetBinContent(k)*100.);
    hCorrectedYieldNormEff_StatError->SetBinError(k,0.);
  }
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // hCorrectedYieldNormEff->Rebin(39, "", fBinsPi013TeVEMCPt);
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  hCorrectedYieldNormEff->Fit(fitBylikin13TeV, "QM0P", "", 1.4, 12.);
  hCorrYieldME->Fit(fitBylikin13TeV_3to8, "QM0P", "", 1.4, 12.);
  hCorrectedYieldNormEff->Fit(ftsallis13TeV, "QM0P", "", 1.4, 12.);

  TF1 *fBylikinRatio = new TF1("fBylikinRatio", "fitBylikin13TeV_3to8/fitBylikin13TeV", 1.4, 12.);
  fBylikinRatio->SetLineWidth(3);
  fBylikinRatio->SetLineColor(kBlue+2);


  hCorrYield_SysError           = (TH1D*) hCorrYieldME->Clone("hCorrYield_SysError");
  hCorrYield_CorrBkgSysError    = (TH1D*) hCorrYieldME->Clone("hCorrYield_CorrBkgSysError");
  hCorrYield_IntSysError        = (TH1D*) hCorrYieldME->Clone("hCorrYield_IntSysError");
  hCorrYield_FitSysErrro        = (TH1D*) hCorrYieldME->Clone("hCorrYield_FitSysErrro");
  hCorrYield_RebinningSysError  = (TH1D*) hCorrYieldME->Clone("hCorrYield_RebinningSysError");

  //////////////////////////////////////////////////////////////////////////////
  // First the counting range variation:
  TFile* OutputFileBetterBkg3to8_HigherInt  = SafelyOpenRootfile("OutputFileBetterBkg3to8_HigherInt.root");
  if (OutputFileBetterBkg3to8_HigherInt->    IsOpen() ) printf("OutputFileBetterBkg3to8_HigherInt opened successfully\n");

  hCorrYield_HigherInt                      = (TH1D*) OutputFileBetterBkg3to8_HigherInt->Get("hYield_dt_chi2map_corrected");


  TFile* OutputFileBetterBkg3to8_SmallInt   = SafelyOpenRootfile("OutputFileBetterBkg3to8_SmallInt.root");
  if (OutputFileBetterBkg3to8_SmallInt->    IsOpen() ) printf("OutputFileBetterBkg3to8_SmallInt opened successfully\n");

  hCorrYield_SmallInt                       = (TH1D*) OutputFileBetterBkg3to8_SmallInt->Get("hYield_dt_chi2map_corrected");

  //////////////////////////////////////////////////////////////////////////////
  // 2nd the param range variation:
  TFile* OutputFileBetterBkg3to8_HigherFit  = SafelyOpenRootfile("OutputFileBetterBkg3to8_HigherFit.root");
  if (OutputFileBetterBkg3to8_HigherFit->   IsOpen() ) printf("OutputFileBetterBkg3to8_HigherFit successfully\n");

  hCorrYield_HigherFit                      = (TH1D*) OutputFileBetterBkg3to8_HigherFit->Get("hYield_dt_chi2map_corrected");

  TFile* OutputFileBetterBkg3to8_SmallFit   = SafelyOpenRootfile("OutputFileBetterBkg3to8_SmallFit.root");
  if (OutputFileBetterBkg3to8_SmallFit->      IsOpen() ) printf("OutputFileBetterBkg3to8_SmallFit successfully\n");

  hCorrYield_SmallFit                       = (TH1D*) OutputFileBetterBkg3to8_SmallFit->Get("hYield_dt_chi2map_corrected");



  //////////////////////////////////////////////////////////////////////////////
  // 3rd the Rebinning  variation
  //
  TFile* OutputFileBetterBkg3to8_LowerRebinning   = SafelyOpenRootfile("OutputFileBetterBkg3to8_LowerRebinning.root");
  if (OutputFileBetterBkg3to8_LowerRebinning->    IsOpen() ) printf("OutputFileBetterBkg3to8_LowerRebinning successfully\n");

  hCorrYield_LowerRebinning                       = (TH1D*) OutputFileBetterBkg3to8_LowerRebinning->Get("hYield_dt_chi2map_corrected");

  TFile* OutputFileBetterBkg3to8_HigherRebinning  = SafelyOpenRootfile("OutputFileBetterBkg3to8_HigherRebinning.root");
  if (OutputFileBetterBkg3to8_HigherRebinning->   IsOpen() ) printf("OutputFileBetterBkg3to8_HigherRebinning opened successfully\n");

  hCorrYield_HigherRebinning                      = (TH1D*) OutputFileBetterBkg3to8_HigherRebinning->Get("hYield_dt_chi2map_corrected");


  if(templatemethod == 1){
    Double_t temp = 0;
    for (int i = 2; i < numberbins; i++) {

      temp = 0;
      // check for biggest diff. in param vari
      // if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_HigherFit->GetBinContent(i)) > temp){
      temp += fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_HigherFit->GetBinContent(i));
      // }
      // if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_SmallFit->GetBinContent(i)) > temp){
      temp += fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_SmallFit->GetBinContent(i));
      // }
      // pushing biggest difference back
      temp = temp/2.;
      vParamSys.push_back(temp);
    }
    for (int i = 2; i < numberbins; i++) {
      // resetting temp
      temp = 0;

      // chech for biggest diff. in count vari
      // if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_HigherInt->GetBinContent(i)) > temp){
      temp += fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_HigherInt->GetBinContent(i));
      // }
      // if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_SmallInt->GetBinContent(i)) > temp){
      temp += fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_SmallInt->GetBinContent(i));
      // }
      temp = temp/2.;
      // pushing biggest difference back
      vCountSys.push_back(temp);
    }
    for (int i = 2; i < numberbins; i++) {
      temp = 0;

      // chech for biggest diff. in BGGit vari
      // if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_LowerRebinning->GetBinContent(i)) > temp){
      temp += fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_LowerRebinning->GetBinContent(i));
      // }
      // if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_HigherRebinning->GetBinContent(i)) > temp){
      temp += fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_HigherRebinning->GetBinContent(i));
      // }

      temp = temp/2.;
      // pushing biggest difference back
      vBGFitSys.push_back(temp);
    }

    //////////////////////////////////////////////////////////////////////////////
    // different corr. bkg templates
    for (int i = 2; i < numberbins; i++) {
      temp = 0;

      // chech for biggest diff. in BGGit vari
      // if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_NNMethod->GetBinContent(i)) > temp){
      temp += fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_NNMethod->GetBinContent(i));
      // }
      // if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_SingleBkg->GetBinContent(i)) > temp){
      temp += fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_SingleBkg->GetBinContent(i));
      // }

      temp = temp/2.;
      // pushing biggest difference back
      vCorrBkgSys.push_back(temp);
    }


    for (int i = 2; i < numberbins; i++) {
      temp = sqrt(1./4.*(pow(vParamSys[i-2], 2.) + pow(vCountSys[i-2], 2.) + pow(vBGFitSys[i-2], 2.)+ pow(vCorrBkgSys[i-2], 2.)));
      vFinalSys.push_back(temp);
      hCorrYield_SysError->         SetBinError(i, vFinalSys[i-2]);
      hCorrYield_IntSysError->      SetBinError(i, vCountSys[i-2]);
      hCorrYield_FitSysErrro->      SetBinError(i, vParamSys[i-2]);
      hCorrYield_RebinningSysError->SetBinError(i, vBGFitSys[i-2]);
      hCorrYield_CorrBkgSysError->  SetBinError(i, vCorrBkgSys[i-2]);
      temp = 0;
    }

    hCorrYield_RelativSyserror            = (TH1D*) hCorrYield_SysError->Clone("hCorrYield_RelativSyserror");
    hCorrYield_RelativSyserror->SetXTitle(pt_str);
    hCorrYield_RelativSyserror->SetYTitle("relative sys. Unsicherheit (%)");

    hCorrYield_IntSysError->SetXTitle(pt_str);
    hCorrYield_IntSysError->SetYTitle("relative sys. Unsicherheit (%)");

    hCorrYield_FitSysErrro->SetXTitle(pt_str);
    hCorrYield_FitSysErrro->SetYTitle("relative sys. Unsicherheit (%)");

    hCorrYield_RebinningSysError->SetXTitle(pt_str);
    hCorrYield_RebinningSysError->SetYTitle("relative sys. Unsicherheit (%)");

    hCorrYield_CorrBkgSysError->SetXTitle(pt_str);
    hCorrYield_CorrBkgSysError->SetYTitle("relative sys. Unsicherheit (%)");


    for(int k = 2; k < numberbins; k++){
      hCorrYield_RelativSyserror->SetBinContent(k, hCorrYield_RelativSyserror->GetBinError(k)/(Double_t)hCorrYield_RelativSyserror->GetBinContent(k)*100.);
      hCorrYield_RelativSyserror->SetBinError(k,0.);

      hCorrYield_IntSysError->SetBinContent(k, hCorrYield_IntSysError->GetBinError(k)/(Double_t)hCorrYield_IntSysError->GetBinContent(k)*100.);
      hCorrYield_IntSysError->SetBinError(k,0.);

      hCorrYield_FitSysErrro->SetBinContent(k, hCorrYield_FitSysErrro->GetBinError(k)/(Double_t)hCorrYield_FitSysErrro->GetBinContent(k)*100.);
      hCorrYield_FitSysErrro->SetBinError(k,0.);

      hCorrYield_RebinningSysError->SetBinContent(k, hCorrYield_RebinningSysError->GetBinError(k)/(Double_t)hCorrYield_RebinningSysError->GetBinContent(k)*100.);
      hCorrYield_RebinningSysError->SetBinError(k,0.);

      hCorrYield_CorrBkgSysError->SetBinContent(k, hCorrYield_CorrBkgSysError->GetBinError(k)/(Double_t)hCorrYield_CorrBkgSysError->GetBinContent(k)*100.);
      hCorrYield_CorrBkgSysError->SetBinError(k,0.);

    }

    hCorrYield_SyserrorRatio            = (TH1D*) hCorrYield_SysError->Clone("hCorrYield_SyserrorRatio");
    hCorrYield_SyserrorRatio->          Divide(hCorrYield_SyserrorRatio, hCorrectedYieldNormEff, 1, 1, "B");
  }


  gDirectory = OutputFile;          // changing directory to the output file
  hCorrYieldME_Ratio->              Write("hCorrYieldME_Ratio");
  hCorrYieldME_StatError->          Write("hCorrYieldME_StatError");
  hCorrectedYieldNormEff_StatError->Write("hCorrectedYieldNormEff_StatError");
  hCorrectedYieldNormEff_Ratio->    Write("hCorrectedYieldNormEff_Ratio");
  if(templatemethod == 1){
    hCorrYield_SysError->             Write("hCorrYield_SysError");
    hCorrYield_IntSysError->          Write("hCorrYield_IntSysError");
    hCorrYield_FitSysErrro->          Write("hCorrYield_FitSysErrro");
    hCorrYield_RebinningSysError->    Write("hCorrYield_RebinningSysError");
    hCorrYield_RelativSyserror->      Write("hCorrYield_RelativSyserror");
    hCorrYield_SyserrorRatio->        Write("hCorrYield_SyserrorRatio");
    hCorrYield_CorrBkgSysError->      Write("hCorrYield_CorrBkgSysError");
  }



delete fitBylikin13TeV;
delete fitBylikin13TeV_3to8;
delete ftsallis13TeV;
delete fBylikinRatio;

InputFile->Close();
OutputFileBetterBkg3to8_HigherInt->Close();
OutputFileBetterBkg3to8_SmallInt->Close();
OutputFileBetterBkg3to8_HigherFit->Close();
OutputFileBetterBkg3to8_SmallFit->Close();
OutputFileBetterBkg3to8_LowerRebinning->Close();
OutputFileBetterBkg3to8_HigherRebinning->Close();
if(templatemethod == 1){
  NNMethodFile->Close();
  SingleBkgFile->Close();
}


}

#endif
