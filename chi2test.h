#include "CommonHeader.h"

TH2D* chi2test(TH1F* hData, TH1F* hSignal, TH1F* hCorrback){
  Double_t chi2 = 0;
  Double_t A_c = 0;               // Area of hData
  Double_t A_b = 0;               // Area of corr. back.
  Double_t A_a = 0;               // Area of signal
  Double_t dx = 1.e-3;
  Double_t dy = 1.e-3;
  Double_t temp_error = 0;
  const int bin0dot3 = 75;

  TH2D* hChi2map = new TH2D("hChi2map", "", 1000, 0., 1., 1000, 0., 1.);
  SetHistoStandardSettings2(hChi2map);

  hChi2map->SetXTitle("normalized signal scaling factor");
  hChi2map->SetYTitle("normalized corr. back. scaling factor");
  hChi2map->SetZTitle("#chi^{2}/dof");

  TH1F* hData_clone = (TH1F*) hData->Clone("hData");
  TH1F* hSignal_clone = (TH1F*) hSignal->Clone("hSignal");
  TH1F* hCorrback_clone = (TH1F*) hCorrback->Clone("hCorrback");

  //////////////////////////////////////////////////////////////////////////////
  // Setting all the bins with pT > 0.3 GeV/c to 0
  //////////////////////////////////////////////////////////////////////////////
  for (int i = bin0dot3; i < 200; i++) {
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
      chi2 = 0;
      for (int i = 0; i < bin0dot3; i++) {

        temp_error = sqrt(pow(dx*ix*hSignal_clone->GetBinError(i), 2.)
        +pow(dy*iy*hCorrback_clone->GetBinError(i), 2.));

        chi2 += pow(dx*ix*hSignal_clone->GetBinContent(i)
        +dy*iy*hCorrback_clone->GetBinContent(i)
        -hData_clone->GetBinContent(i),2.)
        /(pow(temp_error,2.)+pow(hData_clone->GetBinError(i),2.));
      }
      chi2 /= (Double_t)(bin0dot3 -2.);
      hChi2map->Fill(dx*ix, dy*iy, chi2);
    }
  }


  return hChi2map;

  delete hChi2map;
  delete hData_clone;
  delete hSignal_clone;
  delete hCorrback_clone;

}

////////////////////////////////////////////////////////////////////////////////
// Gets the 2D Histo and Searches for Chi2_min
////////////////////////////////////////////////////////////////////////////////
Double_t findChi2_min(TH2D *h1){
  Double_t val = 100;
  for (int ix = 0; ix < 1000; ix++) {
    for (int iy = 0; iy < 1000; iy++) {
      if(h1->GetBinContent(ix,iy) < val)
      {
        val = h1->GetBinContent(ix,iy);
      }
    }
  }
  return val;
}

////////////////////////////////////////////////////////////////////////////////
// Gets the 2D Histo and Searches for the 1simga Error around Chi2_min+1
////////////////////////////////////////////////////////////////////////////////
TH2D *getErrorHist(TString name, TH2D *h1, double val)
{
  TH2D *res = (TH2D *) h1->Clone(name);
  res->Reset();
  int nBinsX = h1->GetNbinsX();
  int nBinsY = h1->GetNbinsY();
  for (int ix = 1; ix <= nBinsX; ix++) {
    for (int iy = 1; iy <= nBinsY; iy++) {
      if (h1->GetBinContent(ix,iy) <= val){res->SetBinContent(ix,iy,2.);}
    }
  }
  return res;
}

////////////////////////////////////////////////////////////////////////////////
// Gets the 2D Histo with 2. values and returns the simga values?
////////////////////////////////////////////////////////////////////////////////
std::vector<double> getErrors(TH2D* h1, double xStart, double yStart)
{
  //returns vector with errors as
  //vec[0] = xLow
  //vec[1] = xhigh
  //vec[2] = yLow
  //vec[3] = yhigh
  std::vector<double> vecOut;
  double errXlow_tmp = 0.;
  double errXhigh_tmp = 0.;
  double errYlow_tmp = 0.;
  double errYhigh_tmp = 0.;
  double errXlow =  xStart;
  double errXhigh = xStart;
  double errYlow =  yStart;
  double errYhigh = yStart;

  for (int ix = 10; ix <= h1->GetNbinsX(); ix++) {
    for (int iy = 10; iy <= h1->GetNbinsY(); iy++) {
      if (h1->GetBinContent(ix,iy)){

        errXlow_tmp = h1->GetXaxis()->GetBinCenter(ix);
        errXhigh_tmp = h1->GetXaxis()->GetBinCenter(ix);
        errYlow_tmp = h1->GetYaxis()->GetBinCenter(iy);
        errYhigh_tmp = h1->GetYaxis()->GetBinCenter(iy);

        if (errXlow_tmp < errXlow){errXlow = errXlow_tmp;}
        if (errXhigh_tmp > errXhigh){errXhigh = errXhigh_tmp;}
        if (errYlow_tmp < errYlow){errYlow = errYlow_tmp;}
        if (errYhigh_tmp > errYhigh){errYhigh = errYhigh_tmp;}
      }
    }
  }
  vecOut.push_back(errXlow);
  vecOut.push_back(errXhigh);
  vecOut.push_back(errYlow);
  vecOut.push_back(errYhigh);
  return vecOut;
}
