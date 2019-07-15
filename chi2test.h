#ifndef CHI2TEST
#define CHI2TEST
#include "CommonHeader.h"

/**
 * Function called by Chi2MapFunction to actually calculated Chi^2 and Ndf
 * @param  h1             data histogram containing (same - scaled mixed) events
 * @param  h2             MC histogram containing MC "true" pi^0 peak
 * @param  h3             correlated background histogram
 * @param  ndf            temp. variable to obtain ndf
 * @param  a              variable for the signal scaling
 * @param  b              variable for the corr. bkg scaling
 * @param  TEMPLATEMETHOD telling the function which method is currently
 * used
 * @param  binnumber      current p_T bin
 * corr. bkg. scaling factor from the 3 to 8 method after one iteration.
 * function which is used as 1 sigma uncertainty as a constraint
 * @return                Chi^2
 */
Double_t Chi2Calc(TH1D* h1, TH1D* h2, TH1D* h3, Double_t &ndf, Double_t a,
                       Double_t b, const int TEMPLATEMETHOD, int binnumber){
  Double_t chi2 = 0;
  Double_t temp_error = 0;
  Int_t lowerfitrange = 0;
  Int_t upperfitrange = 0;

  switch (TEMPLATEMETHOD) {
    case  7:
      lowerfitrange = h3->FindBin(lowerparamrange_narrow);
      upperfitrange = h3->FindBin(upperparamrange_narrow);
      break;
    case  8:
      lowerfitrange = h3->FindBin(lowerparamrange_wide);
      upperfitrange = h3->FindBin(upperparamrange_wide);
      break;
    default:
      lowerfitrange = h3->FindBin(lowerparamrange);
      upperfitrange = h3->FindBin(upperparamrange);
      break;
  }


  /*
  loop over all bins inside the fitting range
  check if the MC histograms have meaninungful entries in these bins

  if so, calculate the chi squared value for this bin and added it's vaule to
  the sum of chi squared values for all bins

  if not, decrease number of degrees of freedom by one

  at the end return the sum of chi squared over all bins
   */
  for (int j = lowerfitrange; j <= upperfitrange; j++) {

    if(h1->GetBinContent(j) != 0 && h1->GetBinError(j) != 0
      && h2->GetBinError(j) != 0 && h2->GetBinContent(j) != 0){

        temp_error = sqrt(
          pow(a*h1->GetBinError(j), 2.)
          +pow(b*h2->GetBinError(j), 2.)
        );

        chi2 += pow(a*h1->GetBinContent(j) + b*h2->GetBinContent(j)
        -h3->GetBinContent(j),2.)
        /(pow(temp_error,2.)+pow(h3->GetBinError(j),2.));
      }

    else{
      ndf -= 1;
    }
  }
  if(a >= 0 && b>= 0)
  {
    return chi2;
  }
  else
  {
    return 0.0;
  }
}

/**
 * Self written function which creates the so called Chi2Map.
 * @param hInvMass_Data       data histogram containing same event - scaled
 * mixed event
 * @param hPeak_MC            MC histogram containing MC truth Pi0 Peak
 * @param hCorrBkg            correlated background histogram
 * @param temp_chi2_dt        temp. variable to obtain min. Chi^2
 * @param signalAreaScaling   temp. variable for the signalAreaScaling
 * @param corrbackAreaScaling temp. variable for the corrbackAreaScaling
 * @param x_min               temp. variable for the signal scaling
 * @param y_min               temp. variable for the corr. bkg scaling
 * @param ndf                 temp. variable to obtain ndf
 * @param TEMPLATEMETHOD      telling the function which method is currently
 * used
 * @param fBinsPi013TeVEMCPt  pT binning
 * @param k                   current PT bin
 */
 TH2D* Chi2MapFunction(TH1D* hData, TH1D* hSignal, TH1D* hCorrback, Double_t &chi2_min,
                       Double_t &signalAreaScaling, Double_t &corrbackAreaScaling,
                       Double_t &x_min, Double_t &y_min, Double_t &ndf,
                       int TEMPLATEMETHOD, Double_t pT, int binnumber,
                       Double_t NEvents_data, Double_t NEvents_MC)
 {

   Double_t chi2_min_temp = 10.e10;
   Double_t chi2 = 0;
   Double_t A_c        = 0;                   // Area of the same - scaled mixed event
   Double_t A_b        = 0;                   // Area of corr. back. template
   Double_t A_a        = 0;                   // Area of the signal template
   Int_t numbersteps   = 101;                 // == the number of bins in the TH2
   Double_t temp_error = 0;                   // Fehlervariable fuer die Templates
   Double_t x_min_temp = 0;
   Double_t y_min_temp = 0;
   TString safePath = gDirectory->GetPath();  // retrieve neutral path


   TH1D* hData_clone = (TH1D*) hData->Clone("hData");
   TH1D* hSignal_clone = (TH1D*) hSignal->Clone("hSignal");
   TH1D* hCorrback_clone = (TH1D*) hCorrback->Clone("hCorrback");

   Int_t lowerfitrange = 0;
   Int_t upperfitrange = 0;
   switch (TEMPLATEMETHOD) {
     case  7:
       lowerfitrange = hData->FindBin(lowerparamrange_narrow);
       upperfitrange = hData->FindBin(upperparamrange_narrow);
       break;
     case  8:
       lowerfitrange = hData->FindBin(lowerparamrange_wide);
       upperfitrange = hData->FindBin(upperparamrange_wide);
       break;
     default:
       lowerfitrange = hData->FindBin(lowerparamrange);
       upperfitrange = hData->FindBin(upperparamrange);
       break;
   }

   /**
    * Setting all the bins with pT > 0.3 GeV/c to 0
    * @param i loop variable indicating
    */
   for (int i = 1; i < hData_clone->FindBin(0.3); i++)
   {
     if(i < lowerfitrange || i > upperfitrange){
       hData_clone->SetBinContent(i,0.);
       hData_clone->SetBinError(i,0.);
       hSignal_clone->SetBinContent(i,0.);
       hSignal_clone->SetBinError(i,0.);
       hCorrback_clone->SetBinContent(i,0.);
       hCorrback_clone->SetBinError(i,0.);
     }
   }

   /**
    * Areascaling method as general idea to compensate for difference in number
    * of Events in MC and data. Can theortically still be used, but is currently
    * not used as standard.
    */
   // A_c = hData_clone->Integral();
   // A_b = hCorrback_clone->Integral();
   // A_a = hSignal_clone->Integral();
   //
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


   /**
   * Chi2Map Creation:
   * First Setting the Stepsize in which Chi2 will be calculated
   * Then create the Chi2Maps
   */

   TH2D* hChi2map      = NULL;
   const int MAXnIterationsChi2Fit  = 2;

   Double_t PMValue = 5.0*pow(10.0, -2);

   hChi2map    = new TH2D("hChi2map", "",
   numbersteps, -1.0*PMValue , 10.0+PMValue,
   numbersteps, -1.0*PMValue , 10.0+PMValue);
   SetHistoStandardSettings2(hChi2map);

   /*****************************************************************************
   **************************First CHI2MAP CALCULATION***************************
   *****************************************************************************/

   Double_t x_point = 0.0;
   Double_t y_point = 0.0;
   for (int ix = 0; ix < numbersteps; ix++) {
     for (int iy = 0; iy < numbersteps; iy++) {
       x_point = hChi2map->GetXaxis()->GetBinCenter(ix);
       y_point = hChi2map->GetYaxis()->GetBinCenter(iy);
       ndf = upperfitrange-lowerfitrange-2;
       chi2 = 0;

       chi2 = Chi2Calc(hSignal_clone, hCorrback_clone, hData_clone, ndf,
                       x_point, y_point, TEMPLATEMETHOD, binnumber);

       if(chi2 < chi2_min_temp && chi2 > 0.0){
         chi2_min_temp = chi2;
         x_min_temp = max(x_point, 0.0);
         y_min_temp = max(y_point, 0.0);
       }

          hChi2map->SetBinContent(ix+1, iy+1, chi2); // since bin (0,0) is underflow!
      }
    }

    chi2_min = chi2_min_temp;
    x_min = x_min_temp;
    y_min = y_min_temp;
    std::cout << "*********************************************************************" << '\n';
    std::cout << "| \u03C7\u00B2 _min = " << chi2_min << std::endl;
    std::cout << "| number of degrees of freedom = " << ndf << std::endl;
    std::cout << "| x_min = " << x_min << '\n';
    std::cout << "| y_min = " << y_min << '\n';
    std::cout << "*********************************************************************" << '\n';

    /****************************************************************************
    **************************LOOP CHI2MAP CALCULATION***************************
    ****************************************************************************/
   for(int nIterationsChi2Fit = 1; nIterationsChi2Fit <= MAXnIterationsChi2Fit; nIterationsChi2Fit++){
     PMValue = 5.0*pow(10.0, -(2+nIterationsChi2Fit));

     // hChi2map      = NULL;
     if(corrbackAreaScaling == 1. && signalAreaScaling == 1.){
       hChi2map    = new TH2D("hChi2map", "",
       numbersteps,
       hChi2map->GetXaxis()->GetBinCenter(hChi2map->GetXaxis()->FindBin(x_min)-5) - PMValue,
       hChi2map->GetXaxis()->GetBinCenter(hChi2map->GetXaxis()->FindBin(x_min)+5) + PMValue,
       numbersteps,
       hChi2map->GetYaxis()->GetBinCenter(hChi2map->GetYaxis()->FindBin(y_min)-5) - PMValue,
       hChi2map->GetYaxis()->GetBinCenter(hChi2map->GetYaxis()->FindBin(y_min)+5) + PMValue);
       SetHistoStandardSettings2(hChi2map);
     }
     else if(corrbackAreaScaling != 1. || signalAreaScaling != 1.){
       hChi2map    = new TH2D("hChi2map", "",
       numbersteps,
       hChi2map->GetXaxis()->GetBinCenter(hChi2map->GetXaxis()->FindBin(x_min)-5) - PMValue,
       hChi2map->GetXaxis()->GetBinCenter(hChi2map->GetXaxis()->FindBin(x_min)+5) + PMValue,
       numbersteps,
       hChi2map->GetYaxis()->GetBinCenter(hChi2map->GetYaxis()->FindBin(y_min)-5) - PMValue,
       hChi2map->GetYaxis()->GetBinCenter(hChi2map->GetYaxis()->FindBin(y_min)+5) + PMValue);
       SetHistoStandardSettings2(hChi2map);
     }
     else{
       std::cerr << "Neither Areascaling nor no Areascaling oO? This doesn't make any sense. ABORT!" << '\n';
       exit(42);
     }

     hChi2map->SetXTitle("SF_{Signal}");
     hChi2map->SetYTitle("SF_{korr. Untergrund}");
     hChi2map->SetZTitle("#chi^{2}");

     for (int ix = 0; ix < numbersteps; ix++) {
       for (int iy = 0; iy < numbersteps; iy++) {
         x_point = hChi2map->GetXaxis()->GetBinCenter(ix);
         y_point = hChi2map->GetYaxis()->GetBinCenter(iy);
         ndf = upperfitrange-lowerfitrange-2;
         chi2 = 0;

         chi2 = Chi2Calc(hSignal_clone, hCorrback_clone, hData_clone, ndf,
                         x_point, y_point, TEMPLATEMETHOD, binnumber);

         if(chi2 < chi2_min_temp && chi2 > 0.0){
           chi2_min_temp = chi2;
           x_min_temp = max(x_point, 0.0);
           y_min_temp = max(y_point, 0.0);
         }

            hChi2map->SetBinContent(ix+1, iy+1, chi2); // since bin (0,0) is underflow!
        }
      }

      // if(nIterationsChi2Fit == MAXnIterationsChi2Fit)
      // {
        chi2_min = chi2_min_temp;
        x_min = x_min_temp;
        y_min = y_min_temp;
        std::cout << "*********************************************************************" << '\n';
        std::cout << "| \u03C7\u00B2 _min = " << chi2_min << std::endl;
        std::cout << "| number of degrees of freedom = " << ndf << std::endl;
        std::cout << "| x_min = " << x_min << '\n';
        std::cout << "| y_min = " << y_min << '\n';
        std::cout << "| Lower x Binedge = " << hChi2map->GetXaxis()->GetBinLowEdge(1) << '\n';
        std::cout << "| Lower x Binedge +1 = " << hChi2map->GetXaxis()->GetBinLowEdge(2) << '\n';
        std::cout << "| Upper x Binedge = " << hChi2map->GetXaxis()->GetBinLowEdge(numbersteps+1) << '\n';
        std::cout << "| Lower y Binedge = " << hChi2map->GetYaxis()->GetBinLowEdge(1) << '\n';
        std::cout << "| Lower y Binedge +1 = " << hChi2map->GetYaxis()->GetBinLowEdge(2) << '\n';
        std::cout << "| Upper y Binedge = " << hChi2map->GetYaxis()->GetBinLowEdge(numbersteps+1) << '\n';
        std::cout << "*********************************************************************" << '\n';
      // }

    }
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
      }
    }
  }
  vec.push_back(posX);
  vec.push_back(posY);
  vec.push_back(val);
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

  for (int ix = 1; ix <= h1->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= h1->GetNbinsY(); iy++) {
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
#endif
