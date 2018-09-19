#include "chi2test.h"                         // self made chi2 calc and mapping
////////////////////////////////////////////////////////////////////////////////
// Function for Signal Template + Pol 1 Param
Double_t mc_full_func2(Double_t *x,  Double_t *par){
  Double_t xx = x[0];
  return (Double_t) par[0]*mc_full_clone2->GetBinContent(hPeak_MC->FindBin(xx)) +
  par[1]+par[2]*xx;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Start of the Main
void IterTempCreation2(std::string current_path, int cutmode){

  TString sPath = gDirectory->GetPath();            // retrieve neutral path

  //////////////////////////////////////////////////////////////////////////////
  // open ESD histo which is inside TLists inside a rootfile
  // There get the Histo with Number of Events (minimum Bias)
  TFile* ESDFile_MC    = SafelyOpenRootfile("./../Daten/" + current_path + ".root");
  if (ESDFile_MC->IsOpen() ) printf("ESDFile_MC opened successfully\n");

  TFile* ESDFile_data             = NULL; // ESD File from the Data
  TFile* BkgFile                  = NULL; // Self makde lower stat. Bkg File
  TList* lGammaCalo_data          = NULL; // TLists inside the ESD File
  TList* lCutNumber_data          = NULL; // TLists inside the ESD File
  TList* lESD_data                = NULL; // Innerst TList for ESD histos from hInvMass_Data
  TList* lGammaCalo_MC            = NULL;
  TList* lCutNumber_MC            = NULL;
  TList* lESD_MC                  = NULL; // Innerst TList for ESD histos from MC
  TList* lMC_MC                   = NULL; // Innerst TList for MC histos from MC
  TList* lTrue_MC                 = NULL; // Innerst TList for True histos from MC
  TH1D* hNEvents                  = NULL; // histo containing number of Events
  TH1D* hMC_Pi0InAcc_Pt           = NULL; // acceptance histo
  TH2D* hTrueDoubleCounting_Pi0   = NULL; // 2D Histo including Doublecounting
  TH1D* hTrueDoubleCounting_Pi0_Pro;      // X-Projection of current
                                          // pT for Double Counting


  ESDFile_data    = SafelyOpenRootfile("./../Daten/GammaCalo-data_503.root");
  if (ESDFile_data->IsOpen() ) printf("ESDFile_data opened successfully\n");

  lGammaCalo_data        = (TList*) ESDFile_data->Get("GammaCalo_503");
  lCutNumber_data        = (TList*) lGammaCalo_data->FindObject("Cut Number 00010113_1111112067032220000_01631031000000d0");
  lESD_data              = (TList*) lCutNumber_data->FindObject("00010113_1111112067032220000_01631031000000d0 ESD histograms");

  lGammaCalo_MC          = (TList*) ESDFile_MC->Get("GammaCalo_503");
  lCutNumber_MC          = (TList*) lGammaCalo_MC->FindObject("Cut Number 00010113_1111112067032220000_01631031000000d0");
  lESD_MC                = (TList*) lCutNumber_MC->FindObject("00010113_1111112067032220000_01631031000000d0 ESD histograms");
  lMC_MC                 = (TList*) lCutNumber_MC->FindObject("00010113_1111112067032220000_01631031000000d0 MC histograms");
  lTrue_MC               = (TList*) lCutNumber_MC->FindObject("00010113_1111112067032220000_01631031000000d0 True histograms");
  hNEvents                = (TH1D*)  lESD_data->FindObject("NEvents");
  hMC_Pi0InAcc_Pt         = (TH1D*)  lMC_MC->FindObject("MC_Pi0InAcc_Pt");
  hTrueDoubleCounting_Pi0 = (TH2D*)  lTrue_MC->FindObject("ESD_TrueDoubleCountPi0_InvMass_Pt");

  //////////////////////////////////////////////////////////////////////////////


  Double_t NEvents  = hNEvents->GetBinContent(1);   // retrieve NEents MinBias
  delete hNEvents;
  ESDFile_MC->Close();
  ESDFile_data->Close();
  //////////////////////////////////////////////////////////////////////////////


  gDirectory->Cd(sPath.Data());                     // for saftey resetting path

  //////////////////////////////////////////////////////////////////////////////
  // declaring global variables
  TString str;                                      // contains the pT range
  Double_t int_error           = 0;                 // contains the Yiled errors
  Double_t temp_ndf            = 0;                 // temp variable containig
                                                    // ndf for Iter Method
  Double_t int_value           = 0;                 // contains -||- values
  // const Int_t numberbins    = numberbins;
  const Int_t ndrawpoints      = 1.e5;              // # points for TF drawing
  // const int n_iter          = 4;
  const int epsilon            = 1.e-6;             // presicion in Iter Method
  int nIterStep                = 0;                 // number of Iterationsteps
  int iterMethodBool           = 1;                 // bool value for IterMethod
                                                    // while loop

  std::vector<Double_t> vInIntRangePercent;         // vector containig needed
                                                    // correction for the
                                                    // efficiency cuz of int
                                                    // boundaries

  std::vector<Double_t> vChi2_DT_Chi2Map;           // vecotr containig Chi2
                                                    // from Chi2MapMethod

  std::vector<Double_t> vNDF_DT_Chi2Map;            // ector containing Ndf
                                                    // from Chi2MapMethod

  Double_t temp_chi2_dt        = 0;                 // temp variable which will
                                                    // hold the current Chi2 and
                                                    // be pushed back into
                                                    // vChi2_DT_Chi2Map

  std::vector<Double_t> vSignalAreaScaling;         // vecotr containig Area-
                                                    // scaling factor for the
                                                    // signal
  Double_t signalAreaScaling   = 0;                 // temp variable which will
                                                    // hold the current Area-
                                                    // scaling factor for the
                                                    // Signal and be pushed
                                                    // back into
                                                    // vSignalAreaScaling

  std::vector<Double_t> vCorrbackAreaScaling;       // vecotr containig Area-
                                                    // scaling factor for the
                                                    // corr. Background
  Double_t corrbackAreaScaling = 0;                 // temp variable which will
                                                    // hold the current Area-
                                                    // scaling factor for the
                                                    // corr. Background and be
                                                    // pushed back into
                                                    // vCorrbackAreaScaling

  std::vector<Double_t> v_x_min;                    // vecotr containig the
                                                    // scaling factor for the
                                                    // Signal
  Double_t x_min               = 0;                 // temp variable which will
                                                    // hold the current scaling
                                                    // factor for the Signal
                                                    // and be pushed back into
                                                    // v_x_min

  std::vector<Double_t> v_y_min;                    // vecotr containig the
                                                    // scaling factor for the
                                                    // corr. Background
  Double_t y_min               = 0;                 // temp variable which will
                                                    // hold the current scaling
                                                    // factor for the corr.
                                                    // Background and be pushed
                                                    // back into v_y_min


  Double_t ndf                 = 0;                 // contains the current ndf

  std::vector<std::vector<Double_t>> vsigma_dt;     // matrix containig the
                                                    // upper and lower errors
                                                    // for the scaling parameter
                                                    // from the Chi2Map Method

  std::vector<Double_t> vChi2_Pol1_Iter;            // vector containig Chi/ndf
                                                    // from the Pol1 Iter method

  TFile *IterTemp;                                  // FilePointer where things
                                                    // will be safed to

  TH2D* hChi2_2D[numberbins];                       // Array of pointers to 2D
                                                    // histos, which contain the
                                                    // Chi2 map

  TH2D* hChi2_2D_sigma[numberbins];                 // Array of TH2D* which
                                                    // contain the 1sigma range


  TH1D* hChi2_pol1_iter[numberbins];                // Array of TH1D* which
                                                    // contain Chi2/ndf for the
                                                    // Pol1 + temp IterMethod

  //////////////////////////////////////////////////////////////////////////////
  // Histos conating the uncorrected yields from the different methods
  TH1D* hYield_dt_chi2map_uncorr  = new TH1D("hYield_dt_chi2map_uncorr",
                                            "", numberbins, fBinsPi013TeVEMCPt);

  TH1D* hYield_pol1_uncorr        = new TH1D("hYield_pol1_uncorr",
                                            "", numberbins, fBinsPi013TeVEMCPt);
  //////////////////////////////////////////////////////////////////////////////
  // setting up the 2 Histograms to compare chi2 from the to fit methods as
  // well as peak factor comp. between pol 1 and double temp fit and the ratio
  // of BG. scaling factor and the Peak scaling factor

  TH1D* hchi2_pol1                = new TH1D("hchi2_pol1",
                                            "", numberbins, fBinsPi013TeVEMCPt);

  TH1D* ha_pol1                   = new TH1D("ha_pol1",
                                            "", numberbins, fBinsPi013TeVEMCPt);
  //////////////////////////////////////////////////////////////////////////////
  // Final Pol1 Iter
  TH1D* hPol1;

  //////////////////////////////////////////////////////////////////////////////
  // Data Clone Histos which will be Integrated to get the uncorrected Yield
  TH1D* data_clone_for_int_dt_chi2map;
  TH1D* data_clone_for_int_pol1;

  //////////////////////////////////////////////////////////////////////////////
  // Standard Settings for all The TH1 histos
  SetHistoStandardSettings(hchi2_pol1);
  SetHistoStandardSettings(ha_pol1);
  SetHistoStandardSettings(hYield_dt_chi2map_uncorr);
  SetHistoStandardSettings(hYield_pol1_uncorr);

  //////////////////////////////////////////////////////////////////////////////
  // Giving the Histos their unique colors depending on the method
  hYield_dt_chi2map_uncorr->SetLineColor(kMagenta+2);
  hYield_dt_chi2map_uncorr->SetMarkerColor(kMagenta+2);
  hYield_pol1_uncorr->SetLineColor(kRed);
  hYield_pol1_uncorr->SetMarkerColor(kRed);


  //////////////////////////////////////////////////////////////////////////////
  // open True Yield Path
  TH1D* CorrectedYieldTrueEff;
  TFile* FData_corrected = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1Correction_00010113_1111112067032220000_01631031000000d0.root");
  if (FData_corrected->IsOpen() ) printf("FData_corrected opened successfully\n");

  CorrectedYieldTrueEff = (TH1D*) FData_corrected->Get(Form("CorrectedYieldNormEff"));
  //////////////////////////////////////////////////////////////////////////////
  // going over all pt bins despite first one, which is some framework bs.
  for (int k = 1; k < numberbins; k++) {
    std::cout << "starte bin " << k << " reading and wrinting!" << std::endl << std::endl;
    TFile* MCFile = NULL;
    TFile* DataFile = NULL;

    ////////////////////////////////////////////////////////////////////////////
    // open MC histo path
    MCFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_MC_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");
    if (MCFile->IsOpen() ) printf("MCFile opened successfully\n");


    ////////////////////////////////////////////////////////////////////////////
    // retrieve MC histograms
    hInvMass_MC = (TH1D*) MCFile->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",k)); //fHistoMappingSignalInvMass_in_Pt_Bin
    hPeak_MC = (TH1D*) MCFile->Get(Form("Mapping_TrueMeson_InvMass_in_Pt_Bin%02d",k));


    DataFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");
    if (DataFile->IsOpen() ) printf("DataFile opened successfully\n");
    hInvMass_Data = (TH1D*) DataFile->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",k));

    ////////////////////////////////////////////////////////////////////////////
    // Getting the purposed corr Background

    BkgFile = SafelyOpenRootfile("./BackFile.root");
    if (BkgFile->IsOpen() ) printf("BkgFile opened successfully\n");
    hCorrBkg = (TH1D*) BkgFile->Get(Form("hPilledUpBack_Bin%02d",k));

    // hCorrBkg = (TH1D*) hInvMass_MC->Clone("hCorrBkg");
    // hCorrBkg->Add(hPeak_MC,-1);


    hPeak_MC->GetXaxis()->SetRangeUser(0.,0.4);
    hInvMass_Data->GetXaxis()->SetRangeUser(0.,0.4);
    hCorrBkg->GetXaxis()->SetRangeUser(0.,0.4);


    hTrueDoubleCounting_Pi0_Pro = hTrueDoubleCounting_Pi0->ProjectionX(Form("hTrueDoubleCounting_Pi0_Pro_bin%02d",k),hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k]),
    hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k+1]));
    Double_t InIntRangeDoubleCounting = (Double_t)hTrueDoubleCounting_Pi0_Pro->Integral(hTrueDoubleCounting_Pi0_Pro->FindBin(lowercountrange),
                                                   hTrueDoubleCounting_Pi0_Pro->FindBin(uppercountrange))/
                                                   (Double_t)hPeak_MC->Integral(
                                                    hPeak_MC->FindBin(lowercountrange),
                                                    hPeak_MC->FindBin(uppercountrange));

    Double_t InIntRangePercent = (Double_t)hPeak_MC->Integral(hPeak_MC->FindBin(lowercountrange),
                                                   hPeak_MC->FindBin(uppercountrange))/
                                                   (Double_t)hMC_Pi0InAcc_Pt->Integral(
                                                     hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k]),
                                                     hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k+1])-1);

    vInIntRangePercent.push_back(InIntRangePercent);

    ////////////////////////////////////////////////////////////////////////////
    // Function for ndf
    TF1* f_ChiOverNdf = new TF1("f_ChiOverNdf", "[0]", 0.0 ,100.);
    f_ChiOverNdf->SetNpx(ndrawpoints);
    f_ChiOverNdf->SetNumberFitPoints(numberbins);
    f_ChiOverNdf->SetLineColor(kBlue+2);
    f_ChiOverNdf->SetLineWidth(4);

    //////////////////////////////////////////////////////////////////////////
    // Fix! Changes < in TLatex to #leq
    str = hInvMass_MC->GetTitle();
    TString str_copy = str.Copy();
    str_copy.ReplaceAll("<","#leq");
    str.Replace(0,20,str_copy,23);

    SetHistoStandardSettings(hInvMass_MC);
    SetHistoStandardSettings(hPeak_MC);
    SetHistoStandardSettings(hCorrBkg);

    hInvMass_Data->SetTitle(str);

    ////////////////////////////////////////////////////////////////////////////
    // normal lame pol 1 fit with template
    TF1* fit_eq_1 = new TF1("fit_eq_1", &mc_full_func2, 0.0, 0.4, 3);
    fit_eq_1->SetNpx(ndrawpoints);
    fit_eq_1->SetNumberFitPoints(numberbins);
    fit_eq_1->SetLineColor(kRed);
    fit_eq_1->SetLineWidth(4);

    //////////////////////////////////////////////////////////////////////////
    // making things look good
    // change y title offset to fit and look nicely
    hInvMass_Data->GetYaxis()->SetTitleOffset(1.2);
    hInvMass_Data->GetYaxis()->SetLabelOffset(0.006);
    hInvMass_Data->SetTitleSize(0.03, "xy");
    hInvMass_Data->SetLabelSize(0.03, "xy");
    hInvMass_Data->SetYTitle("d#it{N}/d#it{m}_{inv} ((GeV/#it{c}^{2})^{-1})");
    hInvMass_Data->SetXTitle("d#it{m}_{inv}");
    hInvMass_Data->SetMarkerStyle(24);
    hInvMass_Data->SetMarkerSize(1.5);
    hInvMass_Data->SetTitle("");
    hInvMass_Data->SetLineWidth(3);
    hCorrBkg->SetLineColor(kCyan+3);
    hCorrBkg->SetMarkerColor(kCyan+3);
    hCorrBkg->SetMarkerStyle(21);
    hCorrBkg->SetMarkerSize(1.5);
    hPeak_MC->SetLineColor(kGreen+3);
    hPeak_MC->SetMarkerColor(kGreen+3);
    hPeak_MC->SetMarkerStyle(33);
    hPeak_MC->SetMarkerSize(2);

    //////////////////////////////////////////////////////////////////////////
    //clone for pol 1 fit
    TH1D* data_clone2 = (TH1D*) hInvMass_Data->Clone("data_clone2");
    TH1D* data_clone3 = (TH1D*) hInvMass_Data->Clone("data_clone3");


    // clearing the vectors
    vChi2_Pol1_Iter.clear();
    vChi2_Pol1_Iter.resize(0);
    iterMethodBool = 1;

    // resetting # of IterSteps to 0
    nIterStep = 0;

    ////////////////////////////////////////////////////////////////////////////
    // while block where the Iter Method is done
    while(iterMethodBool){

      //////////////////////////////////////////////////////////////////////////
      //clone for pol 1 fit
      mc_full_clone2  = (TH1D*) hPeak_MC->Clone("mc_full_clone2");

      ////////////////////////////////////////////////////////////////////////////
      // creating the new root file to safe all the related histograms and fits
      // in it.
      if(k == 1 && nIterStep == 0){
        IterTemp      = new TFile("IterTempBetterBkg.root", "RECREATE");
      }

      else{
        IterTemp      = new TFile("IterTempBetterBkg.root", "UPDATE");
      }
      //////////////////////////////////////////////////////////////////////////
      //fit pol 1 + temp
      TFitResultPtr r_pol1_temp1 = data_clone2->Fit("fit_eq_1", "QM0PS","", lowerparamrange, upperparamrange);
      vChi2_Pol1_Iter.push_back(r_pol1_temp1->Chi2() / r_pol1_temp1->Ndf());

      //////////////////////////////////////////////////////////////////////////
      // scale for pol 1 + temp
      mc_full_clone2->Scale(r_pol1_temp1->Parameter(0));

      //////////////////////////////////////////////////////////////////////////
      // reset data_clone histos and then calculate their new errors
      gDirectory->Cd(sPath.Data());
      data_clone2 = (TH1D*) hInvMass_Data->Clone("data_clone2");
      for(int j = 13; j < 63; j++){
        data_clone2->SetBinError(j,sqrt(pow(data_clone2->GetBinError(j),2.)
        + pow((r_pol1_temp1->Parameter(0)*hPeak_MC->GetBinError(j)),2.)));
      }

      IterTemp->Close();
      if(nIterStep >=1){
        if(fabs(vChi2_Pol1_Iter[nIterStep-1]-vChi2_Pol1_Iter[nIterStep] <= epsilon)){

          //////////////////////////////////////////////////////////////////////
          // for the last iteration step don't reset the clones, instead calc
          // errors for the hInvMass_Data histos that will be used in the final Fit
          // afterwards

          for(int j = 0; j < 75; j++){
          data_clone3->SetBinError(j,sqrt(data_clone3->GetBinError(j) *
          data_clone3->GetBinError(j) + mc_full_clone2->GetBinError(j) *
          mc_full_clone2->GetBinError(j)));

          iterMethodBool = 0;
        }
      }
    }
      nIterStep++;
    }

    ////////////////////////////////////////////////////////////////////////////
    // making the CHi2 monitoring histos!
    hChi2_pol1_iter[k-1]        = new TH1D(Form("1hChi2_pol1_iter_bin%02d",k-1),"",
                                      nIterStep-1, 0.5, (Double_t) nIterStep-0.5);
    SetHistoStandardSettings(hChi2_pol1_iter[k-1]);


    for(int i = 1; i < nIterStep; i++){
      hChi2_pol1_iter[k-1]->SetBinContent(i, vChi2_Pol1_Iter[i-1]);
      hChi2_pol1_iter[k-1]->SetYTitle("#chi^{2}/ndf");
      hChi2_pol1_iter[k-1]->SetXTitle("Iterationstep");
      hChi2_pol1_iter[k-1]->SetLineColor(kRed);
      hChi2_pol1_iter[k-1]->SetMarkerColor(kRed);
    }

    gDirectory->Cd(sPath.Data());

    ///////////////////////////////////////////////////////////////////////////
    // final pol 1 + temp fit
    mc_full_clone2 = (TH1D*) hPeak_MC->Clone("mc_full_clone2");
    TFitResultPtr r_pol1_temp = data_clone3->Fit("fit_eq_1", "M0S","", lowerparamrange, upperparamrange);
    TH1D* mc_full_clone4 = (TH1D*) hPeak_MC->Clone("mc_full_clone4");
    mc_full_clone4->Scale(r_pol1_temp->Parameter(0));
    mc_full_clone4->SetLineColor(kRed);
    mc_full_clone4->SetMarkerColor(kRed);

    TF1* fpol1 = new TF1("fpol1", "[0]+x*[1]", 0.0, 0.4);
    fpol1->SetParameter(0, fit_eq_1->GetParameter(1));
    fpol1->SetParameter(1, fit_eq_1->GetParameter(2));
    fpol1->SetLineColor(kTeal-7);
    fpol1->SetLineWidth(3);

    ////////////////////////////////////////////////////////////////////////////
    //  making full histogram of the Pol 1 Param.
    hPol1 = (TH1D*) mc_full_clone4->Clone("hPol1");
    hPol1->Add(fpol1);
    hPol1->SetMarkerStyle(20);
    hPol1->SetMarkerSize(1.5);
    hPol1->SetMarkerColor(kRed+1);
    hPol1->SetLineColor(kRed+1);

    ////////////////////////////////////////////////////////////////////////////
    // Double Template with Chi2 Map Part
    ////////////////////////////////////////////////////////////////////////////
    hChi2_2D[k-1] = chi2test(hInvMass_Data, hPeak_MC, hCorrBkg, temp_chi2_dt, signalAreaScaling, corrbackAreaScaling, x_min, y_min, ndf);

    hChi2_2D_sigma[k-1] = getErrorHist(Form("hChi2_2D_sigma_bin%02d",k), hChi2_2D[k-1],temp_chi2_dt+1);
    vChi2_DT_Chi2Map.push_back(temp_chi2_dt);
    vNDF_DT_Chi2Map.push_back(ndf);
    vSignalAreaScaling.push_back(signalAreaScaling);
    vCorrbackAreaScaling.push_back(corrbackAreaScaling);
    v_x_min.push_back(x_min);
    v_y_min.push_back(y_min);
    vsigma_dt.push_back(getErrors(hChi2_2D_sigma[k-1], x_min, y_min));
    f_ChiOverNdf->SetParameter(0,temp_chi2_dt/ndf);

    temp_chi2_dt = 0;

    hInvMass_Data->SetTitle(str);
    IterTemp = new TFile("IterTempBetterBkg.root","UPDATE");
    gDirectory = IterTemp;
    hInvMass_Data->Write(Form("data_bin%02d",k));
    data_clone3->Write(Form("data_addedErrosPol1_bin%02d",k));
    mc_full_clone4->Write(Form("mc_peak_pol1_bin%02d",k));
    fpol1->Write(Form("fpol1_bin%02d",k));
    hPol1->Write(Form("hPol1_bin%02d",k));
    hChi2_pol1_iter[k-1]->Write(Form("hChi2_pol1_iter_bin%02d",k));
    hChi2_2D[k-1]->Write(Form("hChi2_2Dbin%02d",k));
    hChi2_2D_sigma[k-1]->Write(Form("hChi2_2D_sigma_bin%02d",k));
    hPeak_MC->Write(Form("hSignal_bin%02d",k));
    hCorrBkg->Write(Form("hCorrBack_bin%02d",k));
    f_ChiOverNdf->Write(Form("f_ChiOverNdf%02d",k));

    gDirectory->Cd(sPath.Data());

    ////////////////////////////////////////////////////////////////////////////
    // getting the chi2 of the current pT bin
    hchi2_pol1->SetBinContent(k+1,r_pol1_temp->Chi2()/temp_ndf);
    hchi2_pol1->SetBinError(k+1, sqrt(2./(temp_ndf+3)));

    data_clone_for_int_dt_chi2map = (TH1D*) hInvMass_Data->Clone("data_clone_for_int_dt_chi2map");
    hCorrBkg->Scale(y_min*corrbackAreaScaling);
    data_clone_for_int_dt_chi2map->Add(hCorrBkg, -1);
    int_value = data_clone_for_int_dt_chi2map->IntegralAndError(data_clone_for_int_dt_chi2map->GetXaxis()->FindBin(lowercountrange), data_clone_for_int_dt_chi2map->GetXaxis()->FindBin(uppercountrange), int_error);
    hYield_dt_chi2map_uncorr->SetBinContent(k+1, int_value);
    hYield_dt_chi2map_uncorr->SetBinError(k+1, int_error);


    data_clone_for_int_pol1 = (TH1D*) hInvMass_Data->Clone("hYield_pol1_uncorr");
    data_clone_for_int_pol1->Add(fpol1, -1);
    int_value = data_clone_for_int_pol1->IntegralAndError(data_clone_for_int_pol1->GetXaxis()->FindBin(lowercountrange), data_clone_for_int_pol1->GetXaxis()->FindBin(uppercountrange), int_error);
    hYield_pol1_uncorr->SetBinContent(k+1, int_value);
    hYield_pol1_uncorr->SetBinError(k+1, int_error);

    ////////////////////////////////////////////////////////////////////////////
    // getting the peakratio of the current pT bin
    Double_t peakscale_pol1, peakerr_pol1;

    peakscale_pol1 = r_pol1_temp->Parameter(0);
    peakerr_pol1 = r_pol1_temp->Error(0);

    ha_pol1->SetBinContent(k+1, peakscale_pol1);
    ha_pol1->SetBinError(k+1, peakerr_pol1);

    ////////////////////////////////////////////////////////////////////////////
    // garbage collection part 1
    delete mc_full_clone2;
    delete mc_full_clone4;
    delete data_clone2;
    delete data_clone3;
    delete hPeak_MC;
    delete hCorrBkg;
    delete hInvMass_Data;
    delete fpol1;
    delete f_ChiOverNdf;
    delete fit_eq_1;

    if(k < numberbins-1){
      // IterTemp->Close();
      MCFile->Close();
      DataFile->Close();
    }
    IterTemp->Close();
    std::cout << "bin number" << k << "Ende" << std::endl << std::endl;
    delete hChi2_pol1_iter[k-1];
    delete hChi2_2D[k-1];
  }
  //////////////////////////////////////////////////////////////////////////////
  // end of the for loop over all 1 <= k < nbins
  //////////////////////////////////////////////////////////////////////////////
  // Chi2Map Histos
  TH1D* hChi2_DT_Chi2map = new TH1D("hChi2_DT_Chi2map", "", numberbins, fBinsPi013TeVEMCPt);
  hChi2_DT_Chi2map->SetYTitle("#chi^{2}/ndf");
  hChi2_DT_Chi2map->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hChi2_DT_Chi2map->SetLineWidth(3);

  TH1D* hSignalAreaScaling = new TH1D("hSignalAreaScaling", "", numberbins, fBinsPi013TeVEMCPt);
  hSignalAreaScaling->SetYTitle("signal areascaling factor");
  hSignalAreaScaling->SetXTitle("#it{p}_{T} (GeV/#it{c})");

  TH1D* hCorrbackAreaScaling = new TH1D("hCorrbackAreaScaling", "", numberbins, fBinsPi013TeVEMCPt);
  hCorrbackAreaScaling->SetYTitle("corr. bkg. areascaling factor");
  hCorrbackAreaScaling->SetXTitle("#it{p}_{T} (GeV/#it{c})");

  TH1D* h_x_min = new TH1D("h_x_min", "", numberbins, fBinsPi013TeVEMCPt);
  h_x_min->SetYTitle("signal scaling factor");
  h_x_min->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  TH1D* h_y_min = new TH1D("h_y_min", "", numberbins, fBinsPi013TeVEMCPt);
  h_y_min->SetYTitle("corr. back. scaling factor");
  h_y_min->SetXTitle("#it{p}_{T} (GeV/#it{c})");


  TH1D* hErrXlow = new TH1D("hErrXlow", "", numberbins, fBinsPi013TeVEMCPt);
  hErrXlow->SetYTitle("lower signal scaling factor uncertainty");
  hErrXlow->SetXTitle("#it{p}_{T} (GeV/#it{c})");

  TH1D* hErrXhigh = new TH1D("hErrXhigh", "", numberbins, fBinsPi013TeVEMCPt);
  hErrXhigh->SetYTitle("upper signal scaling factor uncertainty");
  hErrXhigh->SetXTitle("#it{p}_{T} (GeV/#it{c})");;

  TH1D* hErrYlow = new TH1D("hErrYlow", "", numberbins, fBinsPi013TeVEMCPt);
  hErrYlow->SetYTitle("lower corr. back. scaling factor uncertainty");
  hErrYlow->SetXTitle("#it{p}_{T} (GeV/#it{c})");

  TH1D* hErrYhigh = new TH1D("hErrYhigh", "", numberbins, fBinsPi013TeVEMCPt);
  hErrYhigh->SetYTitle("upper corr. back. scaling factor uncertainty");
  hErrYhigh->SetXTitle("#it{p}_{T} (GeV/#it{c})");


  for (int k = 1; k < numberbins; k++) {
    hChi2_DT_Chi2map->SetBinContent(k+1, vChi2_DT_Chi2Map[k-1]/vNDF_DT_Chi2Map[k-1]);
    hChi2_DT_Chi2map->SetBinError(k+1, sqrt(2./vNDF_DT_Chi2Map[k-1]));
    hSignalAreaScaling->SetBinContent(k, vSignalAreaScaling[k-1]);
    hCorrbackAreaScaling->SetBinContent(k, vCorrbackAreaScaling[k-1]);
    h_x_min->SetBinContent(k+1, v_x_min[k-1]);
    h_y_min->SetBinContent(k+1, v_y_min[k-1]);
    hErrXlow->SetBinContent(k+1, vsigma_dt[k-1][0]);
    hErrXhigh->SetBinContent(k+1, vsigma_dt[k-1][1]);
    hErrYlow->SetBinContent(k+1, vsigma_dt[k-1][2]);
    hErrYhigh->SetBinContent(k+1, vsigma_dt[k-1][3]);
    h_x_min->SetBinError(k+1,
    max(hErrXhigh->GetBinContent(k+1)-h_x_min->GetBinContent(k+1),
    h_x_min->GetBinContent(k+1) - hErrXlow->GetBinContent(k+1)));
    h_y_min->SetBinError(k+1,
    max(hErrYhigh->GetBinContent(k+1) - h_y_min->GetBinContent(k+1),
    h_y_min->GetBinContent(k+1) - hErrYlow->GetBinContent(k+1)));
  }
  //////////////////////////////////////////////////////////////////////////////
  // Writing of Chi2(pT)
  IterTemp = new TFile("IterTempBetterBkg.root","UPDATE");

  hchi2_pol1->SetLineColor(kRed);
  hchi2_pol1->SetLineWidth(3);
  hchi2_pol1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hchi2_pol1->GetYaxis()->SetTitle("#chi^{2}/ndf");


  hchi2_pol1->Write(Form("hchi2_pol1"));
  ha_pol1->Write("Pol1PeakFactor");

  ////////////////////////////////////////////////////////////////////////////
  // open MC histo path for the correction histos
  TFile* CorrectionFile = NULL;
  CorrectionFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_MC_GammaConv_OnlyCorrectionFactor_00010113_1111112067032220000_01631031000000d0.root");
  if (CorrectionFile->IsOpen() ) printf("CorrectionFile opened successfully\n");

  TH1D* hAcc    = (TH1D*) CorrectionFile->Get(Form("fMCMesonAccepPt"));
  TH1D* hEffi   = (TH1D*) CorrectionFile->Get(Form("TrueMesonEffiPt"));

  // correction for 2pi, BR, NEvents, Y, Binwidth
  hYield_pol1_uncorr->Scale(1./(NEvents*2*M_PI*1.6*0.98798),"width");
  hYield_pol1_uncorr->SetYTitle(rawyield);
  hYield_pol1_uncorr->SetXTitle(pt_str);

  hYield_dt_chi2map_uncorr->Scale(1./(NEvents*2*M_PI*1.6*0.98798),"width");
  hYield_dt_chi2map_uncorr->SetYTitle(rawyield);
  hYield_dt_chi2map_uncorr->SetXTitle(pt_str);

  // open Data File for the Yield coming from the Framework and for the Chi2
  // from the framework method
  TFile* DataFile = NULL;
  DataFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");
  if (DataFile->IsOpen() ) printf("DataFile opened successfully\n");

  TH1D* hYield_framework = (TH1D*) DataFile->Get(Form("histoYieldMeson"));
  TH1D* histoChi2_0 = (TH1D*) DataFile->Get(Form("histoChi2_0"));
  SetHistoStandardSettings(histoChi2_0);

  hYield_framework->Divide(hEffi);

  // correcting with bin center (1/pT)
  for (int i = 2; i <= numberbins; i++) {
    hYield_pol1_uncorr->SetBinContent(i,hYield_pol1_uncorr->GetBinContent(i)/hYield_pol1_uncorr->GetBinCenter(i));
    hYield_dt_chi2map_uncorr->SetBinContent(i,hYield_dt_chi2map_uncorr->GetBinContent(i)/hYield_dt_chi2map_uncorr->GetBinCenter(i));
    hYield_pol1_uncorr->SetBinError(i,hYield_pol1_uncorr->GetBinError(i)/hYield_pol1_uncorr->GetBinCenter(i));
    hYield_dt_chi2map_uncorr->SetBinError(i,hYield_dt_chi2map_uncorr->GetBinError(i)/hYield_dt_chi2map_uncorr->GetBinCenter(i));
    hYield_framework->SetBinContent(i,hYield_framework->GetBinContent(i)/hYield_framework->GetBinCenter(i));
    hYield_framework->SetBinError(i,hYield_framework->GetBinError(i)/hYield_framework->GetBinCenter(i));
  }

  //////////////////////////////////////////////////////////////////////////////
  // correcting yield for the acceptance
  TH1D * hYield_pol1_acceptance_corrected       =  (TH1D*) hYield_pol1_uncorr->Clone("hYield_pol1_acceptance_corrected");
  hYield_pol1_acceptance_corrected->Divide(hAcc);

  TH1D * hYield_dt_chi2map_acceptance_corrected =  (TH1D*) hYield_dt_chi2map_uncorr->Clone("hYield_dt_chi2map_acceptance_corrected");
  hYield_dt_chi2map_acceptance_corrected->Divide(hAcc);

  //////////////////////////////////////////////////////////////////////////////
  // correcting yield for the efficiency
  TH1D * hYield_pol1_corrected       =  (TH1D*) hYield_pol1_acceptance_corrected->Clone("hYield_pol1_corrected");

  TH1D * hYield_dt_chi2map_corrected =  (TH1D*) hYield_dt_chi2map_acceptance_corrected->Clone("hYield_dt_chi2map_corrected");

  for (int i = 2; i <= numberbins; i++) {
    hYield_pol1_corrected->SetBinContent(i,hYield_pol1_acceptance_corrected->GetBinContent(i)/vInIntRangePercent[i-2]);
    hYield_dt_chi2map_corrected->SetBinContent(i,hYield_dt_chi2map_acceptance_corrected->GetBinContent(i)/vInIntRangePercent[i-2]);
    hYield_pol1_corrected->SetBinError(i,hYield_pol1_acceptance_corrected->GetBinError(i)/vInIntRangePercent[i-2]);
    hYield_dt_chi2map_corrected->SetBinError(i,hYield_dt_chi2map_acceptance_corrected->GetBinError(i)/vInIntRangePercent[i-2]);
  }

  hYield_pol1_corrected->SetYTitle(strCorrectedYield);
  hYield_dt_chi2map_corrected->SetYTitle(strCorrectedYield);



  hYield_pol1_uncorr->Write("hYield_pol1_uncorr");
  hYield_dt_chi2map_uncorr->Write("hYield_dt_chi2map_uncorr");
  hYield_pol1_acceptance_corrected->Write("hYield_pol1_acceptance_corrected");
  hYield_dt_chi2map_acceptance_corrected->Write("hYield_dt_chi2map_acceptance_corrected");
  hYield_pol1_corrected->Write("hYield_pol1_corrected");
  hYield_dt_chi2map_corrected->Write("hYield_dt_chi2map_corrected");
  hYield_framework->Write("hYield_framework");
  CorrectedYieldTrueEff->Write("hCorrectedYieldTrueEff");
  hChi2_DT_Chi2map->Write("hChi2_DT_Chi2map");
  hSignalAreaScaling->Write("hSignalAreaScaling");
  hCorrbackAreaScaling->Write("hCorrbackAreaScaling");
  h_x_min->Write("h_x_min");
  h_y_min->Write("h_y_min");
  hErrXlow->Write("hErrXlow");
  hErrXhigh->Write("hErrXhigh");
  hErrYlow->Write("hErrYlow");
  hErrYhigh->Write("hErrYhigh");
  histoChi2_0->Write("histoChi2_0");

  //////////////////////////////////////////////////////////////////////////////
  // Garbage collection part 2
  delete hchi2_pol1;
  delete hYield_dt_chi2map_uncorr;
  delete hYield_pol1_uncorr;
  delete hYield_pol1_acceptance_corrected;
  delete hYield_dt_chi2map_acceptance_corrected;
  delete hPol1;
  delete hChi2_DT_Chi2map;
  delete hSignalAreaScaling;
  delete hCorrbackAreaScaling;
  delete h_x_min;
  delete h_y_min;
  delete hErrXlow;
  delete hErrXhigh;
  delete hErrYlow;
  delete hErrYhigh;
  delete hTrueDoubleCounting_Pi0_Pro;
  delete hTrueDoubleCounting_Pi0;
  delete hYield_framework;
  delete histoChi2_0;

  vInIntRangePercent.clear();
  vChi2_DT_Chi2Map.clear();
  vNDF_DT_Chi2Map.clear();
  vSignalAreaScaling.clear();
  vCorrbackAreaScaling.clear();
  v_x_min.clear();
  v_y_min.clear();
  vsigma_dt.clear();
  vChi2_Pol1_Iter.clear();
  IterTemp->Close();
  DataFile->Close();
  BkgFile->Close();
}
