#include "CommonHeader.h"

/**
 * [the main function]
 * @param  current_path   [needed for the function to be calleable for variable datat sets]
 * @param  templatemethod [== 1 uses the backgrund templates from the 3 to 8 method,
 *                         == 2 uses the backgrund templates from the Next Neighbours method]
 * @return                [void]
 */
void Template_CAP(std::string current_path, int templatemethod){


  TString safePath = gDirectory->GetPath();            // retrieve neutral path

  /*****************************************************************************
  Declartion of global variables which are used in this programm
  *****************************************************************************/

  TString str;                                      // contains the pT range

  Double_t int_error           = 0;                 // contains the Yiled errors

  Double_t int_value           = 0;                 // contains -||- values

  const Int_t ndrawpoints      = 1.e5;              // # points for TF1 drawing

  std::vector<Double_t> vInIntRangePercent;         // vector containig needed
                                                    // correction for the
                                                    // efficiency cuz of integral
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

  // matter of change
  Double_t x_min               = 0;                 // temp variable which will
                                                    // hold the current scaling
                                                    // factor for the Signal
                                                    // and be pushed back into
                                                    // v_x_min

  std::vector<Double_t> v_y_min;                    // vecotr containig the
                                                    // scaling factor for the
                                                    // corr. Background

  // matter of change
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

  // matter of change
  TFile *IterTemp;                                  // FilePointer where things
                                                    // will be safed to

  TH2D* hChi2_2D[numberbins];                       // Array of pointers to 2D
                                                    // histos, which contain the
                                                    // Chi2 map

  TH2D* hChi2_2D_sigma[numberbins];                 // Array of TH2D* which
                                                    // contain the 1sigma range


  /**
   * Histogram containig the uncorrected Yield obtained via the selfmade Chi2Map
   * Method
   */
  TH1D* hYield_dt_chi2map_uncorr  = new TH1D("hYield_dt_chi2map_uncorr",
                                            "", numberbins, fBinsPi013TeVEMCPt);

  /**
   * Data histogram clone which will be integrated to get the uncorrected Yield
   * for the Chi2Map
   */
  TH1D* data_clone_for_int_dt_chi2map;

  /**
   * Making the histogram look nice
   */
  SetHistoStandardSettings(hYield_dt_chi2map_uncorr);
  hYield_dt_chi2map_uncorr->SetLineColor(kMagenta+2);
  hYield_dt_chi2map_uncorr->SetMarkerColor(kMagenta+2);

  /**
   * open ESD histo which is inside TLists inside a rootfile. There get the
   * Histo with Number of Events (minimum Bias)
   */
  TFile* ESDFile_MC               = NULL; // ESD File from the MC
  TFile* ESDFile_data             = NULL; // ESD File from the Data
  TFile* BkgFile                  = NULL; // Self makde lower stat. Bkg File
  TList* lGammaCalo_data          = NULL; // TLists inside the ESD File
  TList* lCutNumber_data          = NULL; // TLists inside the ESD File
  TList* lESD_data                = NULL; // Innerst TList for ESD histos from hInvMass_Data
  TList* lGammaCalo_MC            = NULL;
  TList* lCutNumber_MC            = NULL;
  TList* lMC_MC                   = NULL; // Innerst TList for MC histos from MC
  TList* lTrue_MC                 = NULL; // Innerst TList for True histos from MC
  TH1D* hNEvents                  = NULL; // histo containing number of Events
  TH1D* hMC_Pi0InAcc_Pt           = NULL; // acceptance histo
  TH2D* hTrueDoubleCounting_Pi0   = NULL; // 2D Histo including Doublecounting
  TH1D* hTrueDoubleCounting_Pi0_Pro;      // X-Projection of current
                                          // pT for Double Counting

  /**
   * Access the ESD File form the MC simulation for two Histograms:
   * 1st. the Doublecounting Histogram
   * 2nd. the MC Histogram of Pi0 in acceptance as a function of pT
   */
  ESDFile_MC    = SafelyOpenRootfile("./../Daten/" + current_path + ".root");
  if (ESDFile_MC->IsOpen() ) printf("ESDFile_MC opened successfully\n");

  lGammaCalo_MC          = (TList*) ESDFile_MC->Get("GammaCalo_503");
  lCutNumber_MC          = (TList*) lGammaCalo_MC->FindObject("Cut Number 00010113_1111112067032220000_01631031000000d0");
  lMC_MC                 = (TList*) lCutNumber_MC->FindObject("00010113_1111112067032220000_01631031000000d0 MC histograms");
  lTrue_MC               = (TList*) lCutNumber_MC->FindObject("00010113_1111112067032220000_01631031000000d0 True histograms");
  hMC_Pi0InAcc_Pt         = (TH1D*)  lMC_MC->FindObject("MC_Pi0InAcc_Pt");
  hTrueDoubleCounting_Pi0 = (TH2D*)  lTrue_MC->FindObject("ESD_TrueDoubleCountPi0_InvMass_Pt");
  ESDFile_MC->Close();

  /**
   * Access the ESD File form the data for one Histograms:
   * 1st. the NEvents Histogram to get the number of Events with MinBias Trigger
   */
  ESDFile_data    = SafelyOpenRootfile("./../Daten/GammaCalo-data_503.root");
  if (ESDFile_data->IsOpen() ) printf("ESDFile_data opened successfully\n");

  lGammaCalo_data        = (TList*) ESDFile_data->Get("GammaCalo_503");
  lCutNumber_data        = (TList*) lGammaCalo_data->FindObject("Cut Number 00010113_1111112067032220000_01631031000000d0");
  lESD_data              = (TList*) lCutNumber_data->FindObject("00010113_1111112067032220000_01631031000000d0 ESD histograms");
  hNEvents               = (TH1D*)  lESD_data->FindObject("NEvents");

  Double_t NEvents  = hNEvents->GetBinContent(1);   // retrieve NEents MinBias
  ESDFile_data->Close();

  gDirectory->Cd(safePath.Data());                     // for saftey resetting path

  //////////////////////////////////////////////////////////////////////////////
  // open True Yield Path
  TH1D* CorrectedYieldNormEff;
  TFile* FData_corrected = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1Correction_00010113_1111112067032220000_01631031000000d0.root");
  if (FData_corrected->IsOpen() ) printf("FData_corrected opened successfully\n");

  CorrectedYieldNormEff = (TH1D*) FData_corrected->Get(Form("CorrectedYieldNormEff"));
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

    // Bkg made up with 6 NN Bins
    if(mario == 2){
      BkgFile = SafelyOpenRootfile("./BackFile.root");
      if (BkgFile->IsOpen() ) printf("BkgFile opened successfully\n");
      hCorrBkg = (TH1D*) BkgFile->Get(Form("hPilledUpBack_Bin%02d_with%02d_bins",k, numberneighbours));
    }

    // Bkg made up with bins 3 to 8
    if(mario == 1){
      BkgFile = SafelyOpenRootfile("./BackGround3to8.root");
      if (BkgFile->IsOpen() ) printf("BkgFile opened successfully\n");
      hCorrBkg = (TH1D*) BkgFile->Get(Form("hPilledUpBack_Bin%02d_enhanced",k));
    }

    // normal Bkg
    if(mario == 0){
      hCorrBkg = (TH1D*) hInvMass_MC->Clone("hCorrBkg");
      hCorrBkg->Add(hPeak_MC,-1);
    }


    hPeak_MC->GetXaxis()->SetRangeUser(0.,0.3);
    hInvMass_Data->GetXaxis()->SetRangeUser(0.,0.3);
    hCorrBkg->GetXaxis()->SetRangeUser(0.,0.3);


    hTrueDoubleCounting_Pi0_Pro = hTrueDoubleCounting_Pi0->ProjectionX(Form("hTrueDoubleCounting_Pi0_Pro_bin%02d",k),hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k]),
    hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k+1]));
    Double_t InIntRangeDoubleCounting = (Double_t)hTrueDoubleCounting_Pi0_Pro->Integral(hTrueDoubleCounting_Pi0_Pro->FindBin(lowercountrange[k]),
                                                   hTrueDoubleCounting_Pi0_Pro->FindBin(uppercountrange))/
                                                   (Double_t)hPeak_MC->Integral(
                                                    hPeak_MC->FindBin(lowercountrange[k]),
                                                    hPeak_MC->FindBin(uppercountrange));

    Double_t InIntRangePercent = (Double_t)hPeak_MC->Integral(hPeak_MC->FindBin(lowercountrange[k]),
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
        if(mario == 2){
          IterTemp      = new TFile("IterTempBetterBkgNN.root", "RECREATE");
        }
        if(mario == 1){
          IterTemp      = new TFile("IterTempBetterBkg3to8.root", "RECREATE");
        }
        if(mario == 0){
          IterTemp      = new TFile("IterTemp.root", "RECREATE");
        }
      }

      else{
        if(mario == 2){
          IterTemp      = new TFile("IterTempBetterBkgNN.root", "UPDATE");
        }
        if(mario == 1){
          IterTemp      = new TFile("IterTempBetterBkg3to8.root", "UPDATE");
        }
        if(mario == 0){
          IterTemp      = new TFile("IterTemp.root", "UPDATE");
        }
      }
      //////////////////////////////////////////////////////////////////////////
      //fit pol 1 + temp
      TFitResultPtr r_pol1_temp1 = data_clone2->Fit("fit_eq_1", "QM0PS","", lowerparamrange[k], upperparamrange);
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
    TFitResultPtr r_pol1_temp = data_clone3->Fit("fit_eq_1", "M0S","", lowerparamrange[k], upperparamrange);
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
    hChi2_2D[k-1] = chi2test(hInvMass_Data, hPeak_MC, hCorrBkg, temp_chi2_dt,
      signalAreaScaling, corrbackAreaScaling, x_min, y_min, ndf, mario,
      fBinsPi013TeVEMCPt[k], k);

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
    if(mario == 2){
      IterTemp      = new TFile("IterTempBetterBkgNN.root", "UPDATE");
    }
    if(mario == 1){
      IterTemp      = new TFile("IterTempBetterBkg3to8.root", "UPDATE");
    }
    if(mario == 0){
      IterTemp      = new TFile("IterTemp.root", "UPDATE");
    }
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
    temp_ndf = ndf;
    hchi2_pol1->SetBinContent(k+1,r_pol1_temp->Chi2()/temp_ndf);
    hchi2_pol1->SetBinError(k+1, sqrt(2./(temp_ndf+3)));

    data_clone_for_int_dt_chi2map = (TH1D*) hInvMass_Data->Clone("data_clone_for_int_dt_chi2map");
    hCorrBkg->Scale(y_min*corrbackAreaScaling);
    data_clone_for_int_dt_chi2map->Add(hCorrBkg, -1);
    int_value = data_clone_for_int_dt_chi2map->
    IntegralAndError(data_clone_for_int_dt_chi2map->GetXaxis()->
    FindBin(lowercountrange[k]),
    data_clone_for_int_dt_chi2map->GetXaxis()->FindBin(uppercountrange),
    int_error);
    hYield_dt_chi2map_uncorr->SetBinContent(k+1, int_value);
    hYield_dt_chi2map_uncorr->SetBinError(k+1, int_error);


    data_clone_for_int_pol1 = (TH1D*) hInvMass_Data->Clone("hYield_pol1_uncorr");
    data_clone_for_int_pol1->Add(fpol1, -1);
    int_value = data_clone_for_int_pol1->IntegralAndError(data_clone_for_int_pol1->GetXaxis()->FindBin(lowercountrange[k]), data_clone_for_int_pol1->GetXaxis()->FindBin(uppercountrange), int_error);
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



}
