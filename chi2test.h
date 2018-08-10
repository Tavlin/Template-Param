#include "CommonHeader.h"


Double_t chi2_selfmade(TH1D* h1, TH1D* h2, TH1D* h3, Double_t &ndf, Double_t a,
                       Double_t b){
  Double_t chi2 = 0;
  Double_t temp_error = 0;
  Int_t lowerfitrange = h3->FindBin(0.085);
  Int_t upperfitrange = h3->FindBin(0.225);

  for (int j = lowerfitrange; j <= upperfitrange; j++) {

    if(h1->GetBinContent(j) != 0 && h1->GetBinError(j) != 0
      && h2->GetBinError(j) != 0 && h2->GetBinContent(j) != 0){

        temp_error = sqrt(pow(a*h1->GetBinError(j), 2.)
        +pow(b*h2->GetBinError(j), 2.));

        chi2 += pow(a*h1->GetBinContent(j) + b*h2->GetBinContent(j)
        -h3->GetBinContent(j),2.)
        /(pow(temp_error,2.)+pow(h3->GetBinError(j),2.));
      }

    else{
      ndf -= 1;
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  // constraint for parameter b. b should not be too big!
  chi2 += pow(b-a-0.8, 2.)/pow(0.8, 2.);
  if(chi2 == pow(b-a-0.8, 2.)/pow(0.8, 2.)){
    return 1000;
  }
  else{
    return chi2;
  }
}


TH2D* chi2test(TH1D* hData, TH1D* hSignal, TH1D* hCorrback, Double_t &chi2_min,
  Double_t &signalAreaScaling, Double_t &corrbackAreaScaling, Double_t &x_min,
  Double_t &y_min, Double_t &ndf){
  Double_t chi2_min_temp = 10.e10;
  Double_t A_c = 0;                         // Area of hData
  Double_t A_b = 0;                         // Area of corr. back.
  Double_t A_a = 0;                         // Area of signal
  const Double_t dx = 0.005;                // Stepsize in x
  const Double_t dy = 0.01;                // Stepsize in y
  Double_t temp_error = 0;                  // Fehlervariable fuer die Templates
  int binnumber2D = 500;                    // Binzahl ~ Feinheit der Suche
  const int bin0dot3 = 75;                  // Binzahl wo 0.3 GeV/c liegt

  TH2D* hChi2map;

  hChi2map = new TH2D("hChi2map", "", binnumber2D, 0.0, 2.5, binnumber2D, -1.0, 4.0);
  SetHistoStandardSettings2(hChi2map);

  hChi2map->SetXTitle("signal scaling factor");
  hChi2map->SetYTitle("corr. back. scaling factor");
  hChi2map->SetZTitle("#chi^{2}");

  TH1D* hData_clone = (TH1D*) hData->Clone("hData");
  TH1D* hSignal_clone = (TH1D*) hSignal->Clone("hSignal");
  TH1D* hCorrback_clone = (TH1D*) hCorrback->Clone("hCorrback");

  Int_t lowerfitrange = hData_clone->FindBin(0.085);
  Int_t upperfitrange = hData_clone->FindBin(0.225);

  //////////////////////////////////////////////////////////////////////////////
  // Setting all the bins with pT > 0.3 GeV/c to 0
  //////////////////////////////////////////////////////////////////////////////
  for (int i = 0; i < 200; i++) {
    if(i < lowerfitrange || i > upperfitrange){
      hData_clone->SetBinContent(i,0.);
      hData_clone->SetBinError(i,0.);
      hSignal_clone->SetBinContent(i,0.);
      hSignal_clone->SetBinError(i,0.);
      hCorrback_clone->SetBinContent(i,0.);
      hCorrback_clone->SetBinError(i,0.);
    }
  }
  // A_c = hData_clone->Integral();
  // A_b = hCorrback_clone->Integral();
  // A_a = hSignal_clone->Integral();

  // if(A_b < 0){
  //   hCorrback_clone->Scale(-1.*(A_c/A_b));
  //   corrbackAreaScaling = -1.*(A_c/A_b);
  // }
  // else{
  //   hCorrback_clone->Scale(A_c/A_b);
  //   corrbackAreaScaling = A_c/A_b;
  // }
  // hSignal_clone->Scale(A_c/A_a);
  // signalAreaScaling = A_c/A_a;

  corrbackAreaScaling = 1.;
  signalAreaScaling = 1.;

  for (int ix = 0; ix < binnumber2D; ix++) {
    for (int iy = 0; iy < binnumber2D; iy++) {
      ndf = 61-16;
      Double_t chi2 = 0;
      chi2 = chi2_selfmade(hSignal_clone, hCorrback_clone, hData_clone, ndf,
                           dx*(Double_t)ix, dy*(Double_t)(iy-100));

      if(chi2 < chi2_min_temp){
        chi2_min_temp = chi2;
        x_min = (Double_t)ix*dx;
        y_min = (Double_t)(iy-100)*dy;
      }
      hChi2map->SetBinContent(ix+1, iy+1, chi2);
    }
  }
  chi2_min = chi2_min_temp;
  std::cout << "chi^2_min = " << chi2_min << std::endl;
  std::cout << "ndf = " << ndf << std::endl;


  return hChi2map;

}

////////////////////////////////////////////////////////////////////////////////
// Gets the 2D Histo and Searches for Chi2_min
////////////////////////////////////////////////////////////////////////////////
std::vector<double> GetMinmumTH2 (TH2D *h){
  std::vector<double> vec;
  int nBinsX = h->GetNbinsX();
  int nBinsY = h->GetNbinsY();
  double posX = 0.;
  double posY = 0.;
  double val = h->GetBinContent(1,1);
  double val_tmp = 0.;
  for (int ix = 1; ix <= nBinsX; ++ix){
    for (int iy = 1; iy <= nBinsY; ++iy){
      val_tmp = h->GetBinContent(ix,iy);
      if(val_tmp < val){
        val = val_tmp;
        posX = h->GetXaxis()->GetBinCenter(ix);
        posY = h->GetYaxis()->GetBinCenter(iy);
        // cout << "Pos X: " << posX << "\tPos Y: " << posY << "\t val: " << val << endl;
      }
    }
  }
  vec.push_back(posX);
  vec.push_back(posY);
  vec.push_back(val);
  // cout << "Min X: " << gr->GetX()[1] << "\tMin Y: " << gr->GetY()[1] << endl;
  return vec;
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
