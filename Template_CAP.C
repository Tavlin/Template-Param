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
  TFile* MCFile                   = NULL; // File containig the MC outcome
  TFile* DataFile                 = NULL; // File containig the data outcome

  TList* lGammaCalo_data          = NULL; // TLists inside the ESD File (data)
  TList* lCutNumber_data          = NULL; // TLists inside the ESD File (data)
  TList* lESD_data                = NULL; // Innerst TList for ESD histos from hInvMass_Data

  TList* lGammaCalo_MC            = NULL; // TLists inside the ESD File (MC)
  TList* lCutNumber_MC            = NULL; // TLists inside the ESD File (MC)
  TList* lMC_MC                   = NULL; // Innerst TList for MC histos from MC
  TList* lTrue_MC                 = NULL; // Innerst TList for True histos from MC

  TH1D* hNEvents                  = NULL; // histo containing number of Events
  TH1D* hMC_Pi0InAcc_Pt           = NULL; // acceptance histo
  TH2D* hTrueDoubleCounting_Pi0   = NULL; // 2D Histo including Doublecounting
  TH1D* CorrectedYieldNormEff;    = NULL; // Corrected Yield from the Framwork

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

  gDirectory->Cd(safePath.Data());                  // for saftey resetting path

  /**
   * Open the file which contains the corrected true norm efficiency Yield from
   * the framework.
   */
  TFile* FData_corrected = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1Correction_00010113_1111112067032220000_01631031000000d0.root");
  if (FData_corrected->IsOpen() ) printf("FData_corrected opened successfully\n");

  CorrectedYieldNormEff = (TH1D*) FData_corrected->Get(Form("CorrectedYieldNormEff"));

  FData_corrected->Close();

  /**
   * For loop to loop over all pT bins definded by fBinsPi013TeVEMCPt. Bin 0 is
   * excluded since it only contains 0 <= pT (GeV/c) < 1.4 where no data should
   * be present, since we have an energy cut at 0.7 GeV per Cluster.
   */
  for (int k = 1; k < numberbins; k++) {

    std::cout << "Start bin  " << k << " reading and wrinting!" << "\n\n";

    /**
     * resetting the Filepointers just for "saftey" reasons
     */
    MCFile = NULL;
    DataFile = NULL;

    /**
     * Open the file which contains the MC output of the framework's work so to
     * say.
     */
    MCFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_MC_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");
    if (MCFile->IsOpen() ) printf("MCFile opened successfully\n");


    /**
     * Histogram from the MC simulation which was created through the normal
     * analysis method: same events - scaled mixed events.
     */
    hInvMass_MC = (TH1D*) MCFile->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",k));

    /**
     * Histogram from the MC simulation which contains the only the true Pi0s
     * coming from y y; y_conv y; and double y_conv
     */
    hPeak_MC = (TH1D*) MCFile->Get(Form("Mapping_TrueMeson_InvMass_in_Pt_Bin%02d",k));

    /**
    * Open the file which contains the data output of the framework's work so to
    * say.
     */
    DataFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");
    if (DataFile->IsOpen() ) printf("DataFile opened successfully\n");

    /**
    * Histogram from the data which was created through the normal analysis
    * method: same events - scaled mixed events.
     */
    hInvMass_Data = (TH1D*) DataFile->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",k));

    /**
     * Obtaining the correlated background histograms for the corresponding
     * method.
     */
    if(templatemethod == 1){
      BkgFile = SafelyOpenRootfile("./BackGround3to8.root");
      if (BkgFile->IsOpen() ) printf("BkgFile opened successfully\n");
      hCorrBkg = (TH1D*) BkgFile->Get(Form("hPilledUpBack_Bin%02d_enhanced",k));
    }

    if(templatemethod == 2){
      BkgFile = SafelyOpenRootfile("./BackFile.root");
      if (BkgFile->IsOpen() ) printf("BkgFile opened successfully\n");
      hCorrBkg = (TH1D*) BkgFile->Get(Form("hPilledUpBack_Bin%02d_with%02d_bins",k, numberneighbours));
    }

    /**
     * Calculating the % of how many Pi0s are in the True Histo from the MC in
     * the current pT range, divided by the amount of produced Pi0s which should
     * hit the Detector (here: EMCal)
     * This is the value of our efficiency in the corresponding pT interval
     */
    Double_t InIntRangePercent = (Double_t)hPeak_MC->Integral(hPeak_MC->FindBin(lowercountrange[k]),
                                                   hPeak_MC->FindBin(uppercountrange))/
                                                   (Double_t)hMC_Pi0InAcc_Pt->Integral(
                                                     hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k]),
                                                     hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k+1])-1);

    vInIntRangePercent.push_back(InIntRangePercent);

    /**
     * little fix for the string which contains the pT intervall comming from
     * the framework.
     * Like: xy < pT (GeV/c) < yz ---> xy <= pT (GeV/c)
     */
    str = hInvMass_MC->GetTitle();
    TString str_copy = str.Copy();
    str_copy.ReplaceAll("<","#leq");
    str.Replace(0,20,str_copy,23);

    hInvMass_Data->SetTitle(str);

    /**
     * creating the new root file(s) to safe all the related histograms and fits
     * in it.
     */
    if(k == 1){
      if(templatemethod == 2){
        IterTemp      = new TFile("IterTempBetterBkgNN.root", "RECREATE");
      }
      if(templatemethod == 1){
        IterTemp      = new TFile("IterTempBetterBkg3to8.root", "RECREATE");
    }

    else{
      if(templatemethod == 2){
        IterTemp      = new TFile("IterTempBetterBkgNN.root", "UPDATE");
      }
      if(templatemethod == 1){
        IterTemp      = new TFile("IterTempBetterBkg3to8.root", "UPDATE");
    }

    gDirectory->Cd(sPath.Data()); // reset path so no nwe generated pointer is
                                  // connected to the Files above.

    ////////////////////////////////////////////////////////////////////////////
    // Double Template with Chi2 Map Part
    ////////////////////////////////////////////////////////////////////////////
    /**
     * Self written function which creates the so called Chi2Map.
     * @param hInvMass_Data       [description]
     * @param hPeak_MC            [description]
     * @param hCorrBkg            [description]
     * @param temp_chi2_dt        [description]
     * @param signalAreaScaling   [description]
     * @param corrbackAreaScaling [description]
     * @param x_min               [description]
     * @param y_min               [description]
     * @param ndf                 [description]
     * @param mario               [description]
     * @param fBinsPi013TeVEMCPt  [description]
     * @param k                   [description]
     */
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
    if(templatemethod == 2){
      IterTemp      = new TFile("IterTempBetterBkgNN.root", "UPDATE");
    }
    if(templatemethod == 1){
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
