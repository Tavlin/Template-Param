#ifndef Systematics
#define Systematics

#include "CommonHeader.h"

void systematics(int templatemethod, TFile* OutputFile){

  TFile* InputFile = NULL;

  if(templatemethod == 2){
    InputFile      = SafelyOpenRootfile("OutputFileBetterBkgNN.root");
  }
  else if(templatemethod == 1){
    InputFile      = SafelyOpenRootfile("OutputFileBetterBkg3to8.root");
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
  TH1D* hCorrYieldME_Ratio                = NULL;
  // TH1D* hCorrYield_count0d06to0d225      = NULL;
  // TH1D* hCorrYield_count0d1to0d225       = NULL;
  // TH1D* hCorrYield_count0d02to0d185      = NULL;
  // TH1D* hCorrYield_count0d02to0d285      = NULL;
  // TH1D* hCorrYield_countsyserror         = NULL;
  // TH1D* hCorrYield_param0d06to0d225      = NULL;
  // TH1D* hCorrYield_param0d1to0d225       = NULL;
  // TH1D* hCorrYield_param0d02to0d185      = NULL;
  // TH1D* hCorrYield_param0d02to0d285      = NULL;
  // TH1D* hCorrYield_paramsyserror         = NULL;
  // TH1D* hCorrYield_BGFitRange0d25        = NULL;
  // TH1D* hCorrYield_BGFitRange0d29        = NULL;
  // TH1D* hCorrYield_BGLeft                = NULL;
  // TH1D* hCorrYield_BGFitsyserror         = NULL;
  // TH1D* hCorrYield_syserror              = NULL;
  // TH1D* hCorrYield_RelativSyserror       = NULL;
  // TFile* FStatUnc                        = NULL;

  TF1* fitBylikin13TeV = new TF1("fitBylikin13TeV", "[0]*exp(-(sqrt(x^(2)+0.135^(2))-0.135)/[1])+[2]/((1+x*x/[3])^([4]))", 1.4, 12.);
  fitBylikin13TeV->SetParameters(13, 0.1, 2, 0.7, 2.9);
  fitBylikin13TeV->SetLineWidth(3);
  fitBylikin13TeV->SetLineColor(kBlack);

  TF1* fitBylikin13TeV_3to8 = new TF1("fitBylikin13TeV_3to8", "[0]*exp(-(sqrt(x^(2)+0.135^(2))-0.135)/[1])+[2]/((1+x*x/[3])^([4]))", 1.4, 12.);
  fitBylikin13TeV_3to8->SetParameters(13, 0.1, 2, 0.7, 2.9);
  fitBylikin13TeV_3to8->SetLineWidth(3);
  fitBylikin13TeV_3to8->SetLineColor(kBlue+2);

  TF1 *ftsallis13TeV = new TF1("ftsallis13TeV", "[0]/(2*3.1415) * ( ([1]-1) * ([1]-2) )/( [1]*[2] *( [1]*[2] + [3] *([1]-2)) )* pow(( 1 + ( ( pow(([3]*[3]+x*x),0.5) -[3]) /( [1]*[2] ) ) ), -[1])", 1.4, 12.);
  ftsallis13TeV->SetParameter(0, 9.4); // A
  ftsallis13TeV->SetParameter(1, 7.169); // n
  ftsallis13TeV->SetParameter(2, 0.159); // T
  ftsallis13TeV->FixParameter(3, 0.1349766); // M
  ftsallis13TeV->SetLineWidth(3);
  ftsallis13TeV->SetLineColor(kBlack);


  std::vector<Double_t> vCountSys;
  std::vector<Double_t> vParamSys;
  std::vector<Double_t> vBGFitSys;
  std::vector<Double_t> vFinalSys;

  // TFile* JoshuasFile      = SafelyOpenRootfile("Pi0_data_GammaConvV1Correction_00010113_1111112067032220000_01631031000000d0.root");

  hCorrYieldME            = (TH1D*) InputFile->Get("hYield_dt_chi2map_corrected");

  hCorrectedYieldNormEff  = (TH1D*) InputFile->Get("hCorrectedYieldNormEff");

  hCorrYieldME_Ratio = (TH1D*) hCorrYieldME->Clone("hCorrYieldME_Ratio");
  hCorrYieldME_Ratio->Divide(hCorrYieldME_Ratio, hCorrectedYieldNormEff, 1, 1, "B");


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

  gDirectory = OutputFile;          // changing directory to the output file
  hCorrYieldME_Ratio->              Write("hCorrYieldME_Ratio");
  hCorrYieldME_StatError->          Write("hCorrYieldME_StatError");
  hCorrectedYieldNormEff_StatError->Write("hCorrectedYieldNormEff_StatError");

  // hCorrYield_syserror           = (TH1D*) hCorrYieldME->Clone("hCorrYield_syserror");
  // hCorrYield_countsyserror      = (TH1D*) hCorrYieldME->Clone("hCorrYield_countsyserror");
  // hCorrYield_paramsyserror      = (TH1D*) hCorrYieldME->Clone("hCorrYield_paramsyserror");
  // hCorrYield_BGFitsyserror      = (TH1D*) hCorrYieldME->Clone("hCorrYield_BGFitsyserror");
  //
  // //////////////////////////////////////////////////////////////////////////////
  // // First the counting range variation:
  // TFile* IterTemp_count0d06to0d225 = SafelyOpenRootfile("IterTemp_count0d06to0d225.root");
  // if (IterTemp_count0d06to0d225->IsOpen() ) printf("IterTemp_count0d06to0d225 opened successfully\n");
  //
  // hCorrYield_count0d06to0d225 = (TH1D*) IterTemp_count0d06to0d225->Get("hYield_dt_chi2map_corrected");
  //
  // TFile* IterTemp_count0d1to0d225 = SafelyOpenRootfile("IterTemp_count0d1to0d225.root");
  // if (IterTemp_count0d1to0d225->IsOpen() ) printf("IterTemp_count0d1to0d225 opened successfully\n");
  //
  // hCorrYield_count0d1to0d225 = (TH1D*) IterTemp_count0d1to0d225->Get("hYield_dt_chi2map_corrected");
  //
  // TFile* IterTemp_count0d02to0d185 = SafelyOpenRootfile("IterTemp_count0d02to0d185.root");
  // if (IterTemp_count0d02to0d185->IsOpen() ) printf("IterTemp_count0d02to0d185 opened successfully\n");
  //
  // hCorrYield_count0d02to0d185 = (TH1D*) IterTemp_count0d02to0d185->Get("hYield_dt_chi2map_corrected");
  //
  // TFile* IterTemp_count0d02to0d285 = SafelyOpenRootfile("IterTemp_count0d02to0d285.root");
  // if (IterTemp_count0d02to0d285->IsOpen() ) printf("IterTemp_count0d02to0d285 opened successfully\n");
  //
  // hCorrYield_count0d02to0d285 = (TH1D*) IterTemp_count0d02to0d285->Get("hYield_dt_chi2map_corrected");
  //
  //
  //
  //
  // //////////////////////////////////////////////////////////////////////////////
  // // 2nd the param range variation:
  // TFile* IterTemp_param0d06to0d225 = SafelyOpenRootfile("IterTemp_param0d06to0d225.root");
  // if (IterTemp_param0d06to0d225->IsOpen() ) printf("IterTemp_param0d06to0d225 opened successfully\n");
  //
  // hCorrYield_param0d06to0d225 = (TH1D*) IterTemp_param0d06to0d225->Get("hYield_dt_chi2map_corrected");
  //
  // TFile* IterTemp_param0d1to0d225 = SafelyOpenRootfile("IterTemp_param0d1to0d225.root");
  // if (IterTemp_param0d1to0d225->IsOpen() ) printf("IterTemp_param0d1to0d225 opened successfully\n");
  //
  // hCorrYield_param0d1to0d225 = (TH1D*) IterTemp_param0d1to0d225->Get("hYield_dt_chi2map_corrected");
  //
  // TFile* IterTemp_param0d02to0d185 = SafelyOpenRootfile("IterTemp_param0d02to0d185.root");
  // if (IterTemp_param0d02to0d185->IsOpen() ) printf("IterTemp_param0d02to0d185 opened successfully\n");
  //
  // hCorrYield_param0d02to0d185 = (TH1D*) IterTemp_param0d02to0d185->Get("hYield_dt_chi2map_corrected");
  //
  // TFile* IterTemp_param0d02to0d285 = SafelyOpenRootfile("IterTemp_param0d02to0d285.root");
  // if (IterTemp_param0d02to0d285->IsOpen() ) printf("IterTemp_param0d02to0d285 opened successfully\n");
  //
  // hCorrYield_param0d02to0d285 = (TH1D*) IterTemp_param0d02to0d285->Get("hYield_dt_chi2map_corrected");
  //
  //
  //
  // //////////////////////////////////////////////////////////////////////////////
  // // 3rd the bg fitting variation
  // TFile* IterTemp_BGFitRange0d25 = SafelyOpenRootfile("IterTemp_BGFitRange0d25.root");
  // if (IterTemp_BGFitRange0d25->IsOpen() ) printf("IterTemp_BGFitRange0d25 opened successfully\n");
  //
  // hCorrYield_BGFitRange0d25 = (TH1D*) IterTemp_BGFitRange0d25->Get("hYield_dt_chi2map_corrected");
  //
  // TFile* IterTemp_BGLeft = SafelyOpenRootfile("IterTemp_BGLeft.root");
  // if (IterTemp_BGLeft->IsOpen() ) printf("IterTemp_BGLeft opened successfully\n");
  //
  // hCorrYield_BGLeft = (TH1D*) IterTemp_BGFitRange0d25->Get("hYield_dt_chi2map_corrected");
  //
  // TFile* IterTemp_BGFitRange0d29 = SafelyOpenRootfile("IterTemp_BGFitRange0d29.root");
  // if (IterTemp_BGFitRange0d29->IsOpen() ) printf("IterTemp_BGFitRange0d29 opened successfully\n");
  //
  // hCorrYield_BGFitRange0d29 = (TH1D*) IterTemp_BGFitRange0d29->Get("hYield_dt_chi2map_corrected");
  //
  //
  // Double_t temp = 0;
  // for (int i = 2; i < numberbins; i++) {
  //
  //   temp = 0;
  //   // check for biggest diff. in param vari
  //   if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_param0d02to0d285->GetBinContent(i)) > temp){
  //     temp = fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_param0d02to0d285->GetBinContent(i));
  //   }
  //   if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_param0d06to0d225->GetBinContent(i)) > temp){
  //     temp = fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_param0d06to0d225->GetBinContent(i));
  //   }
  //   if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_param0d1to0d225->GetBinContent(i)) > temp){
  //     temp = fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_param0d1to0d225->GetBinContent(i));
  //   }
  //   if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_param0d02to0d185->GetBinContent(i)) > temp){
  //     temp = fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_param0d02to0d185->GetBinContent(i));
  //   }
  //   // pushing biggest difference back
  //   vParamSys.push_back(temp);
  //
  //   // resetting temp
  //   temp = 0;
  //
  //   // chech for biggest diff. in count vari
  //   if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_count0d06to0d225->GetBinContent(i)) > temp){
  //     temp = fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_count0d06to0d225->GetBinContent(i));
  //   }
  //   if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_count0d1to0d225->GetBinContent(i)) > temp){
  //     temp = fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_count0d1to0d225->GetBinContent(i));
  //   }
  //   if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_count0d02to0d185->GetBinContent(i)) > temp){
  //     temp = fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_count0d02to0d185->GetBinContent(i));
  //   }
  //   if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_count0d02to0d285->GetBinContent(i)) > temp){
  //     temp = fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_count0d02to0d285->GetBinContent(i));
  //   }
  //   // pushing biggest difference back
  //   vCountSys.push_back(temp);
  //
  //   // chech for biggest diff. in BGGit vari
  //
  //   if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_BGFitRange0d25->GetBinContent(i)) > temp){
  //     temp = fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_BGFitRange0d25->GetBinContent(i));
  //   }
  //   if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_BGLeft->GetBinContent(i)) > temp){
  //     temp = fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_BGLeft->GetBinContent(i));
  //   }
  //   // if(fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_BGFitRange0d29->GetBinContent(i)) > temp){
  //   //   temp = fabs(hCorrYieldME->GetBinContent(i)-hCorrYield_BGFitRange0d29->GetBinContent(i));
  //   // }
  //
  //   // pushing biggest difference back
  //   vBGFitSys.push_back(temp);
  //
  //   temp = sqrt(pow(vParamSys[i-2], 2.) + pow(vCountSys[i-2], 2.) + pow(vBGFitSys[i-2], 2.));
  //   vFinalSys.push_back(temp);
  //   hCorrYield_syserror->SetBinError(i, vFinalSys[i-2]);
  //   hCorrYield_countsyserror->SetBinError(i, vCountSys[i-2]);
  //   hCorrYield_paramsyserror->SetBinError(i, vParamSys[i-2]);
  //   hCorrYield_BGFitsyserror->SetBinError(i, vBGFitSys[i-2]);
  //   temp = 0;
  // }
  //


delete fitBylikin13TeV;
delete fitBylikin13TeV_3to8;
delete ftsallis13TeV;
delete fBylikinRatio;

InputFile->Close();



}

#endif
