#include "CommonHeader.h"


}
Double_t chi2test(TH1F* hData, TH1F* hSignal, TH1F* hCorrback){
  Double_t chi2 = 0;
  Double_t A_c = 0;               // Area of hData
  Double_t A_b = 0;               // Area of corr. back.
  Double_t A_a = 0;               // Area of signal

  TH2F* hChi2map = new TH2F("hChi2map", "", 1000, 0., 1., 1000, 0., 1.);
  TH1F* hData_clone = (TH1F*) hData->Clone("hData");
  TH1F* hSignal_clone = (TH1F*) hSignal->Clone("hSignal");
  TH1F* hCorrback_clone = (TH1F*) hCorrback->Clone("hCorrback");

  //////////////////////////////////////////////////////////////////////////////
  // Setting all the bins with pT > 0.3 GeV/c to 0
  //////////////////////////////////////////////////////////////////////////////
  for (int i = 75; i < 200; i++) {
    hData_clone->SetBinContent(i,0.);
    hData_clone->SetBinError(i,0.);
    hSignal_clone->SetBinContent(i,0.);
    hSignal_clone->SetBinError(i,0.);
    hCorrback_clone->SetBinContent(i,0.);
    hCorrback_clone->SetBinError(i,0.);
  }
  A_c = hData_clone->Integral();
  A_b = hCorrback_clone->Integral();
  A_a = hSignal_clone->Integral();

  hCorrback_clone->Scale(A_c/A_b);
  hSignal_clone->Scale(A_c/A_a);

  for (int ix = 0; ix < 1000; ix++) {
    for (int iy = 0; iy < 1000; iy++) {
      for (int i = 0; i < 75; i++) {
        chi2 += pow(,2.)/(pow(,2.)+pow(,2.));
      }
    }
  }



  delete hChi2map;
  delete hData_clone;
  delete hSignal_clone;
  delete hCorrback_clone;

  return chi2;
}
