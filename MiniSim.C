#include "CommonHeader.h"

Double_t marvinsbinning[40]                      =   {0.0,  1.4,   1.6,   1.8,   2.0,   2.2,
                                                      2.4,   2.6,   2.8,   3.0,   3.2,
                                                      3.4,   3.6,   3.8,   4.0,   4.2,
                                                      4.4,   4.6,   4.8,   5.0,   5.2,
                                                      5.4,   5.6,   5.8,   6.0,   6.2,
                                                      6.4,   6.6,   6.8,   7.0,   7.5,
                                                      8.0,   8.5,   9.0,   9.5,  10.0,
                                                     12.0,  14.0,  16.0,  20.0};

Double_t joshuasbinning[106]                      = {      // size: 106
                                                            0.00, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7,
                                                            1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2,
                                                            2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7,
                                                            2.75, 2.8, 2.85, 2.9, 2.95, 3.0, 3.10, 3.2, 3.30, 3.4,
                                                            3.50, 3.6, 3.70, 3.8, 3.90, 4.0, 4.10, 4.2, 4.30, 4.4,
                                                            4.50, 4.6, 4.70, 4.8, 4.90, 5.0, 5.10, 5.2, 5.30, 5.4,
                                                            5.50, 5.6, 5.70, 5.8, 5.90, 6.0, 6.10, 6.2, 6.30, 6.4,
                                                            6.50, 6.6, 6.70, 6.8, 6.90, 7.0, 7.10, 7.2, 7.30, 7.4,
                                                            7.50, 7.6, 7.70, 7.8, 7.90, 8.0, 8.20, 8.4, 8.60, 8.8,
                                                            9.00, 9.2, 9.40, 9.6, 9.80, 10.0,10.5,11.0, 11.5,12.0,
                                                            13.0, 14.0, 15.0, 16.0, 18.0, 20.0 };

// Mini Sim to check if different binning results in need for bin shifting in
// the corrected yield
void MiniSim(){
  gRandom->SetSeed(time(NULL));

  // Bylikin Function
  TF1* fitBylikin13TeV = new TF1("fitBylikin13TeV", "[0]*exp(-(sqrt(x^(2)+0.135^(2))-0.135)/[1])+[2]/((1+x*x/[3])^([4]))", 1.4, 12.);
  fitBylikin13TeV->SetParameters(13, 0.1, 2, 0.7, 2.9);
  fitBylikin13TeV->SetLineWidth(3);
  fitBylikin13TeV->SetLineColor(kBlack);

  // Getting the corrected Yield for the Fit from Joshuas Analysis
  TFile* JoshuasFile     = SafelyOpenRootfile("Pi0_data_GammaConvV1Correction_00010113_1111112067032220000_01631031000000d0.root");
  if (JoshuasFile->IsOpen() ) printf("JoshuasFile opened successfully\n");
  TH1D* hCorrectedYieldNormEff = (TH1D*) JoshuasFile->Get("CorrectedYieldNormEff");

  hCorrectedYieldNormEff->Fit(fitBylikin13TeV, "QM0P", "", 1.4, 12.);

  TH1D* hCorrYield_MarvinBinning = new TH1D("hCorrYield_MarvinBinning", "", 39, marvinsbinning);
  TH1D* hCorrYield_JoshuaBinning = new TH1D("hCorrYield_JoshuaBinning", "", 105, joshuasbinning);

  SetHistoStandardSettings(hCorrYield_MarvinBinning);
  SetHistoStandardSettings(hCorrYield_JoshuaBinning);

  for (int i = 0; i < 1.e8; i++) {
    Double_t fillvalue = fitBylikin13TeV->GetRandom(1.4, 12.);
    hCorrYield_MarvinBinning->Fill(fillvalue /*, fitBylikin13TeV->Eval(fillvalue)*/);
    hCorrYield_JoshuaBinning->Fill(fillvalue /*, fitBylikin13TeV->Eval(fillvalue)*/);
  }
  hCorrYield_MarvinBinning->Scale(1./1.e8, "width");
  hCorrYield_JoshuaBinning->Scale(1./1.e8, "width");

  hCorrYield_MarvinBinning->SetLineColor(kBlue+2);
  hCorrYield_MarvinBinning->SetMarkerColor(kBlue+2);

  hCorrYield_MarvinBinning->SaveAs("hCorrYield_MarvinBinning.root");
  hCorrYield_JoshuaBinning->SaveAs("hCorrYield_JoshuaBinning.root");

}
