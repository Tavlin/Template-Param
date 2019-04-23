#include "BackGroundFitting.h"
#include "chi2test.h"
#include "Systematics.h"

/**
 * [the main function]
 * @param  TEMPLATEMETHOD [==  0 uses the background templates from the 3 to 8 method,
 *                         ==  1 uses the background templates from the Next Neighbours method
 *                         ==  2 uses the background templates from the same pT intervals as the signal templates
 *                         ==  3 uses lower rebinning output              with TEMPLATEMETHOD 1
 *                         ==  4 uses higher rebinning output             with TEMPLATEMETHOD 1
 *                         ==  5 uses narrow uncorr. bck fit range output with TEMPLATEMETHOD 1
 *                         ==  6 uses wide uncorr. bck   fit range output with TEMPLATEMETHOD 1
 *                         ==  7 uses narrow             fit range output with TEMPLATEMETHOD 1
 *                         ==  8 uses wide               fit range output with TEMPLATEMETHOD 1
 *                         ==  9 uses narrow           count range output with TEMPLATEMETHOD 1
 *                         == 10 uses wide             count range output with TEMPLATEMETHOD 1]
 *
 * @param  ESD_MC_PATH                        == path to the ESD MC
 * @param  ESD_DATA_PATH                      == path to the ESD data
 * @param  FrameworkOutput_Normal             == path to the normal framework output
 * @param  FrameworkOutput_NoRebinning        == path to the framework output without rebinning
 * @param  FrameworkOutput_HigherRebinning    == path to the framework output with higher rebinning
 * @param  FrameworkOutput_LowerRebinning     == path to the framework output with lower rebinning
 * @param  FrameworkOutput_NarrowUncorrBckFit == path to the framework output with narrow uncorr back fit range
 * @param  FrameworkOutput_WideUncorrBckFit   == path to the framework output with wide uncorr back fit range
 * @param  CutString                          == cut string (since the cut string log can contain more then one and this is made only for one cut string!)
 */
void Template_CAP(const int TEMPLATEMETHOD,
                  const char* ESD_MC_PATH,                        const char* ESD_DATA_PATH,
                  const char* FrameworkOutput_Normal,             const char* FrameworkOutput_NoRebinning,
                  const char* FrameworkOutput_HigherRebinning,    const char* FrameworkOutput_LowerRebinning,
                  const char* FrameworkOutput_NarrowUncorrBckFit, const char* FrameworkOutput_WideUncorrBckFit,
                  const char* CutString){


  TString safePath = gDirectory->GetPath();         // retrieve neutral path
  TString strESD_MC_PATH                = ESD_MC_PATH;
  TString strESD_DATA_PATH              = ESD_DATA_PATH;
  TString strCutString                  = CutString;
  TString strFrameworkOutput            = "";
  TString strFrameworkOutputNoRebinning = "";
  TString strMCWithoutCorrection        = "";
  TString strDataWithoutCorrection      = "";
  TString strDataCorrection             = "";
  TString strMCOnlyCorrectionFactor     = "";
  TString strMCCorrectionHistos         = "";

  strFrameworkOutputNoRebinning = FrameworkOutput_NoRebinning;
  strFrameworkOutputNoRebinning = strFrameworkOutputNoRebinning + strCutString + "/13TeV/Pi0_MC_GammaConvV1WithoutCorrection_" + strCutString + ".root";

  switch (TEMPLATEMETHOD){
    case 1:
      strFrameworkOutput = FrameworkOutput_Normal;
      break;
    case 2:
      strFrameworkOutput = FrameworkOutput_Normal;
      break;
    case 3:
      strFrameworkOutput = FrameworkOutput_LowerRebinning;
      break;
    case 4:
      strFrameworkOutput = FrameworkOutput_HigherRebinning;
      break;
    case 5:
      strFrameworkOutput = FrameworkOutput_NarrowUncorrBckFit;
      break;
    case 6:
      strFrameworkOutput = FrameworkOutput_WideUncorrBckFit;
      break;
    case 7:
      strFrameworkOutput = FrameworkOutput_Normal;
      break;
    case 8:
      strFrameworkOutput = FrameworkOutput_Normal;
      break;
    case 9:
      strFrameworkOutput = FrameworkOutput_Normal;
      break;
    case 10:
      strFrameworkOutput = FrameworkOutput_Normal;
      break;
    default:
      strFrameworkOutput = FrameworkOutput_Normal;
      break;
  }

  strMCWithoutCorrection    = strFrameworkOutput + strCutString + "/13TeV/Pi0_MC_GammaConvV1WithoutCorrection_" + strCutString + ".root";
  strDataWithoutCorrection  = strFrameworkOutput + strCutString + "/13TeV/Pi0_data_GammaConvV1WithoutCorrection_" + strCutString + ".root";
  strDataCorrection         = strFrameworkOutput + strCutString + "/13TeV/Pi0_data_GammaConvV1Correction_" + strCutString + ".root";
  strMCOnlyCorrectionFactor = strFrameworkOutput + strCutString + "/13TeV/Pi0_MC_GammaConv_OnlyCorrectionFactor_" + strCutString + ".root";
  strMCCorrectionHistos     = strFrameworkOutput + strCutString + "/13TeV/Pi0_MC_GammaConvV1WithoutCorrection_" + strCutString + ".root";

  /*****************************************************************************
  Declartion of global variables which are used in this programm
  *****************************************************************************/

  TString str;                                      // contains the pT range

  Double_t int_error           = 0;                 // contains the pi0 yield errors

  Double_t int_value           = 0;                 // contains -||- values

  const Int_t ndrawpoints      = 1.e5;              // # points for TF1 drawing

  std::vector<Double_t> vMyEffi;                    // vector containig needed
                                                    // correction for the
                                                    // efficiency cuz of integral
                                                    // boundaries

  std::vector<Double_t> vMyEffiUncer;               // vector containig needed
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
  TFile *OutputFile;                                // FilePointer where things
                                                    // will be safed to

  TH2D* hChi2Map[numberbins];                       // Array of pointers to 2D
                                                    // histos, which contain the
                                                    // Chi2 map

  TH2D* hChi2Map_sigma[numberbins];                 // Array of TH2D* which
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
   * open ESD histo which is inside TLists inside a rootfile. There get the
   * Histo with Number of Events (minimum Bias)
   */
  TFile* ESDFile_MC               = NULL; // ESD File from the MC
  TFile* ESDFile_data             = NULL; // ESD File from the Data
  TFile* BkgFile                  = NULL; // Self makde lower stat. Bkg File
  TFile* MCFile                   = NULL; // File containig the MC outcome
  TFile* DataFile                 = NULL; // File containig the data outcome
  TFile* CorrBkgFile              = NULL; // File that contains the corr. bkg
  TFile* CorrectionFile           = NULL; // File which contains the correction-
                                          // histograms for effi and acceptance

  TList* lGammaCalo_data          = NULL; // TLists inside the ESD File (data)
  TList* lCutNumber_data          = NULL; // TLists inside the ESD File (data)
  TList* lESD_data                = NULL; // Innerst TList for ESD histos from hInvMass_Data

  TList* lGammaCalo_MC            = NULL; // TLists inside the ESD File (MC)
  TList* lCutNumber_MC            = NULL; // TLists inside the ESD File (MC)
  TList* lMC_MC                   = NULL; // Innerst TList for MC histos from MC
  TList* lTrue_MC                 = NULL; // Innerst TList for True histos from MC
  TList* lESD_MC                  = NULL; // Innerst TList for ESD histos for NEvents_data

  TH1D* hNEvents_data             = NULL; // histo containing number of Events in Data
  TH1D* hNEvents_MC               = NULL; // histo containing number of Events in MC
  TH1D* hMC_Pi0InAcc_Pt           = NULL; // acceptance histo
  TH2D* hTrueDoubleCounting_Pi0   = NULL; // 2D Histo including Doublecounting
  TH1D* CorrectedYieldNormEff     = NULL; // Corrected Yield from the Framwork

  TH1D* hYieldMC                  = new TH1D("hYieldMC",
                                            "", numberbins, fBinsPi013TeVEMCPt);
                                          // Yield if we only count the signal
                                          // template not data-corr_back
                                          //
  TH1D* hYieldMC_acc_corr        = NULL;
  TH1D* hYieldMC_effi_corr       = NULL;
  TH1D* hPi0_gen                 = NULL;

  /**
   * Open the MC file which contains the correction histograms for efficiency and
   * acceptance. Obtaining those two directly afterwards.
   */
  CorrectionFile = SafelyOpenRootfile(strMCOnlyCorrectionFactor.Data());
  if (CorrectionFile->IsOpen() ) printf("CorrectionFile opened successfully\n");

  TH1D* hAcc    = (TH1D*) CorrectionFile->Get(Form("fMCMesonAccepPt"));
  TH1D* hEffi   = (TH1D*) CorrectionFile->Get(Form("TrueMesonEffiPt"));

  /**
   * Access the ESD File form the MC simulation for two Histograms:
   * 1st. the Doublecounting Histogram
   * 2nd. the MC Histogram of Pi0 in acceptance as a function of pT
   */
  ESDFile_MC    = SafelyOpenRootfile(strESD_MC_PATH.Data());
  if (ESDFile_MC->IsOpen() ) printf("ESDFile_MC opened successfully\n");

  lGammaCalo_MC          = (TList*) ESDFile_MC->Get("GammaCalo_2020");
  lCutNumber_MC          = (TList*) lGammaCalo_MC->FindObject("Cut Number 00010113_4117911067032230000_01631031000000d0");
  lMC_MC                 = (TList*) lCutNumber_MC->FindObject("00010113_4117911067032230000_01631031000000d0 MC histograms");
  lTrue_MC               = (TList*) lCutNumber_MC->FindObject("00010113_4117911067032230000_01631031000000d0 True histograms");
  lESD_MC                 = (TList*) lCutNumber_MC->FindObject("00010113_4117911067032230000_01631031000000d0 ESD histograms");
  hMC_Pi0InAcc_Pt         = (TH1D*)  lMC_MC->FindObject("MC_Pi0InAcc_Pt");
  hTrueDoubleCounting_Pi0 = (TH2D*)  lTrue_MC->FindObject("ESD_TrueDoubleCountPi0_InvMass_Pt");
  hNEvents_MC             = (TH1D*)  lESD_MC->FindObject("NEvents");
  Double_t NEvents_MC  =
  hNEvents_MC->GetBinContent(1) +(
    hNEvents_MC->GetBinContent(1)/(
      hNEvents_MC->GetBinContent(1)+hNEvents_MC->GetBinContent(5)
    )
  )*hNEvents_MC->GetBinContent(6);
  ESDFile_MC->Close();

  /**
   * Access the ESD File form the data for one Histograms:
   * 1st. the NEvents_data Histogram to get the number of Events with MinBias Trigger
   */
  ESDFile_data    = SafelyOpenRootfile(strESD_DATA_PATH.Data());
  if (ESDFile_data->IsOpen() ) printf("ESDFile_data opened successfully\n");

  lGammaCalo_data        = (TList*) ESDFile_data->Get("GammaCalo_2020");
  lCutNumber_data        = (TList*) lGammaCalo_data->FindObject("Cut Number 00010113_4117911067032230000_01631031000000d0");
  lESD_data              = (TList*) lCutNumber_data->FindObject("00010113_4117911067032230000_01631031000000d0 ESD histograms");
  hNEvents_data          = (TH1D*)  lESD_data->FindObject("NEvents");

  Double_t NEvents_data  = //hNEvents_data->GetBinContent(1);   // retrieve NEents MinBias
  hNEvents_data->GetBinContent(1) +(
    hNEvents_data->GetBinContent(1)/(
      hNEvents_data->GetBinContent(1)+hNEvents_data->GetBinContent(5)
    )
  )*hNEvents_data->GetBinContent(6);
  ESDFile_data->Close();

  gDirectory->Cd(safePath.Data());                  // for saftey resetting path

  /**
   * Calls the CorrBkgCreation function from BackGroundFitting.h This function
   * should make the needed normal corr. bkg templates for the other two
   * methods.
   */
  CorrBkgCreation(strFrameworkOutputNoRebinning.Data());

  /**
   * For loop to loop over all pT bins definded by fBinsPi013TeVEMCPt. Bin 0 is
   * excluded since it only contains 0 <= pT (GeV/c) < 1.4 where no data should
   * be present, since we have an energy cut at 0.7 GeV per Cluster.
   */
  for (int k = 1; k <= numberbins; ++k) {

    if(k >= 39){
      std::cout << "k ist zu gross!" << '\n';
      continue;
    }
    std::cout << "Start bin  " << k << " reading and wrinting!" << "\n\n";

    /**
     * resetting the Filepointers and TH1 pointers just for "saftey" reasons
     */
    MCFile        = NULL;
    DataFile      = NULL;
    hInvMass_MC   = NULL;
    hPeak_MC      = NULL;
    hInvMass_Data = NULL;
    CorrBkgFile   = NULL;
    hCorrBkg      = NULL;

    /**
     * Open the file which contains the MC output of the framework's work so to
     * say.
     */
    MCFile = SafelyOpenRootfile(strMCWithoutCorrection.Data());
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
    hPeak_MC = (TH1D*) MCFile->Get(Form("Mapping_TrueFullMeson_InvMass_in_Pt_Bin%02d",k));

    /**
    * Open the file which contains the data output of the framework's work so to
    * say.
     */
    DataFile = SafelyOpenRootfile(strDataWithoutCorrection.Data());
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
    switch (TEMPLATEMETHOD) {
      case 1:
        hCorrBkg = NULL;
        hCorrBkg = BackgroundAdding(k);                // from BackGroundFitting.h
        break;

      case 2:
        CorrBkgFile = SafelyOpenRootfile("CorrBkgFileNoRebin.root");
        hCorrBkg    = NULL;
        hCorrBkg    = (TH1D*) CorrBkgFile->Get(Form("hCorrBkgNoRebinBin%02d", k));
        hCorrBkg->Rebin(fBinsPi013TeVEMCPtRebin[k]);
        break;

      default:
        // CorrBkgFile = SafelyOpenRootfile("CorrBkgFileNoRebin.root");
        // hCorrBkg    = NULL;
        // hCorrBkg    = (TH1D*) CorrBkgFile->Get(Form("hCorrBkgNoRebinBin%02d", k));
        // hCorrBkg->Rebin(fBinsPi013TeVEMCPtRebin[k]);
        hCorrBkg = NULL;
        hCorrBkg = BackGround3to8(k, TEMPLATEMETHOD);   	             // from BackGroundFitting.h
        std::cout << "After BackGround3to8 1" << '\n';
        break;
    }

    std::cout << "After BackGround3to8 2" << '\n';

    /**
     * little fix for the string which contains the pT intervall comming from
     * the framework.
     * Like: xy < pT (GeV/c) < yz ---> xy <= pT (GeV/c)
     */
    str = hInvMass_MC->GetTitle();
    TString str_copy = str.Copy();
    str_copy.ReplaceAll("<","#leq");
    str.Replace(0,20,str_copy,23);

    std::cout << "After BackGround3to8 3" << '\n';

    hInvMass_Data->SetTitle(str);

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
    hChi2Map[k-1] = Chi2MapFunction(hInvMass_Data, hPeak_MC, hCorrBkg, temp_chi2_dt,
      signalAreaScaling, corrbackAreaScaling, x_min, y_min, ndf, TEMPLATEMETHOD,
      fBinsPi013TeVEMCPt[k], k, NEvents_data, NEvents_MC);

    /**
     * Function from Sebastian to calculate the 1 sigma region around the min
     * Chi^2. The return value is again a TH2D!
     * @param k                 current PT bin
     * @param hChi2Map[k-1]     Chi2Map from above
     * @param temp_chi2_dt+1    temp. variable containig min. Chi^2
     */
    hChi2Map_sigma[k-1] = getErrorHist(Form("hChi2_2D_sigma_bin%02d",k),
    hChi2Map[k-1] ,temp_chi2_dt+1);

    /**
     * Adding all the information we want to monitor in the corresponding
     * vectors
     */
    vChi2_DT_Chi2Map.push_back(temp_chi2_dt);
    vNDF_DT_Chi2Map.push_back(ndf);
    vSignalAreaScaling.push_back(signalAreaScaling);
    vCorrbackAreaScaling.push_back(corrbackAreaScaling);
    v_x_min.push_back(x_min);
    v_y_min.push_back(y_min);
    vsigma_dt.push_back(getErrors(hChi2Map_sigma[k-1], x_min, y_min));

    temp_chi2_dt = 0;           // resetting the temp. variable for min. Chi^2

    /**
     * creating the new root file(s) to safe all the related histograms and fits
     * in it.
     */
    if(k == 1){
      switch (TEMPLATEMETHOD){
        case 1:
          OutputFile      = new TFile("OutputFileBetterBkgNN.root", "RECREATE");
          break;
        case 2:
          OutputFile      = new TFile("OutputFileNormal.root", "RECREATE");
          break;
        case 3:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_LowerRebinning.root", "RECREATE");
          break;
        case 4:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_HigherRebinning.root", "RECREATE");
          break;
        case 5:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_UncorrBckSmallFitRange.root", "RECREATE");
          break;
        case 6:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_UncorrBckWideFitRange.root", "RECREATE");
          break;
        case 7:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_SmallFit.root", "RECREATE");
          break;
        case 8:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_HigherFit.root", "RECREATE");
          break;
        case 9:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_SmallInt.root", "RECREATE");
          break;
        case 10:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_HigherInt.root", "RECREATE");
          break;
        default:
          OutputFile      = new TFile("OutputFileBetterBkg3to8.root", "RECREATE");
          break;
      }
    }

    else{
      switch (TEMPLATEMETHOD){
        case 1:
          OutputFile      = new TFile("OutputFileBetterBkgNN.root", "UPDATE");
          break;
        case 2:
          OutputFile      = new TFile("OutputFileNormal.root", "UPDATE");
          break;
        case 3:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_LowerRebinning.root", "UPDATE");
          break;
        case 4:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_HigherRebinning.root", "UPDATE");
          break;
        case 5:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_UncorrBckSmallFitRange.root", "UPDATE");
          break;
        case 6:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_UncorrBckWideFitRange.root", "UPDATE");
          break;
        case 7:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_SmallFit.root", "UPDATE");
          break;
        case 8:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_HigherFit.root", "UPDATE");
          break;
        case 9:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_SmallInt.root", "UPDATE");
          break;
        case 10:
          OutputFile      = new TFile("OutputFileBetterBkg3to8_HigherInt.root", "UPDATE");
          break;
        default:
          OutputFile      = new TFile("OutputFileBetterBkg3to8.root", "UPDATE");
          break;
      }
    }

    gDirectory->Cd(safePath.Data());  // reset path so no nwe generated pointer is
                                      // connected to the Files above.

    gDirectory = OutputFile;          // changing directory to the output file

    /**
     * wrinting all the wanted histograms for plotting purposes in the output
     * file. part 1
     */
    hInvMass_Data->       Write(Form("data_bin%02d",k));
    hChi2Map[k-1]->       Write(Form("hChi2MapBin%02d",k));
    hChi2Map_sigma[k-1]-> Write(Form("hChi2_2D_sigma_bin%02d",k));
    hPeak_MC->            Write(Form("hSignal_bin%02d",k));
    hCorrBkg->            Write(Form("hCorrBack_bin%02d",k));
    OutputFile->Close();

    gDirectory->cd(safePath.Data()); // resetting directory again.



    for(int bin = 1; bin <= hCorrBkg->GetNbinsX(); bin++){
      hCorrBkg->SetBinContent(bin, hCorrBkg->GetBinContent(bin)*y_min*corrbackAreaScaling);
      hCorrBkg->SetBinError(bin, sqrt(pow(hCorrBkg->GetBinError(bin)*y_min*corrbackAreaScaling,2.)+
                                      pow(hCorrBkg->GetBinContent(bin)*max(vsigma_dt[k-1][3] - y_min,
                                      y_min - vsigma_dt[k-1][2])*corrbackAreaScaling,2.)));
    }

    /**
     * Calculating the % of how many Pi0s are in the True Histo from the MC in
     * the current pT range, divided by the amount of produced Pi0s which should
     * hit the Detector (here: EMCal)
     * This is the value of our efficiency in the corresponding pT interval
     */
    Int_t lowerEffirange = 0;
    Int_t upperEffirange = 0;
    switch (TEMPLATEMETHOD) {
      case  9:
        lowerEffirange = hPeak_MC->FindBin(lowercountrange_narrow);
        upperEffirange = hPeak_MC->FindBin(uppercountrange_narrow);
        break;
      case 10:
        lowerEffirange = hPeak_MC->FindBin(lowercountrange_wide);
        upperEffirange = hPeak_MC->FindBin(uppercountrange_wide);
        break;
      default:
        lowerEffirange = hPeak_MC->FindBin(lowercountrange);
        upperEffirange = hPeak_MC->FindBin(uppercountrange);
        break;
    }
    Double_t IntInRangeError  = 0;
    Double_t IntInRange       = hPeak_MC->IntegralAndError(lowerEffirange, upperEffirange, IntInRangeError);
    Double_t IntAcceptedError = 0;
    Double_t IntAccepted      = hMC_Pi0InAcc_Pt->IntegralAndError(
      hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k]),
      hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k+1]) - 1,
      IntAcceptedError
    );
    Double_t myEffi = (IntInRange)/(IntAccepted);
    Double_t myEffi_Uncer = sqrt(
      pow(IntInRangeError/IntAccepted, 2.) +
      pow((IntInRange*IntAcceptedError)/pow(IntAccepted, 2.) , 2.)
    );

    vMyEffi.push_back(myEffi);
    vMyEffiUncer.push_back(myEffi_Uncer);



    /**
     * Data histogram clone which will be integrated to get the uncorrected Yield
     * for the Chi2Map
     */
    data_clone_for_int_dt_chi2map = (TH1D*) hInvMass_Data->Clone("data_clone_for_int_dt_chi2map");

    // hCorrBkg->Scale(y_min*corrbackAreaScaling);       //sclaing of the corr. bkg.

    /**
     * Subtracting the scaled corr. bkg. from the same - scaled mixed event data
     * histogram (clone).
     * @param hCorrBkg (now scaled) corr. bkg histogram
     */
    data_clone_for_int_dt_chi2map->Add(hCorrBkg, -1);

    /**
     * Obtaining the integral value and error for the uncorrected yield in the
     * current pT bin.
     * Then giving theses values into the uncorrected Yield histogram.
     */
    std::cout << "lowerEffirange = " << lowerEffirange << '\n';
    int_value = data_clone_for_int_dt_chi2map->IntegralAndError(
      lowerEffirange, upperEffirange, int_error);

    hYield_dt_chi2map_uncorr->SetBinContent(k+1, int_value);
    hYield_dt_chi2map_uncorr->SetBinError(k+1, int_error);

    hYieldMC->SetBinContent(
      k+1, hPeak_MC->Integral(lowerEffirange, upperEffirange )
    );

    /**
     * Garbage collection part 1.
     */
    hPeak_MC      = NULL;
    hCorrBkg      = NULL;
    hInvMass_Data = NULL;
    hChi2Map[k-1] = NULL;

    MCFile->Close();
    DataFile->Close();
    if(TEMPLATEMETHOD == 4){
      CorrBkgFile->Close();
    }

    std::cout << "bin number " << k << " reading and writing... DONE!" << "\n\n";
  }
  /**
   * end of the for loop over all 1 <= k < nbins
   */


  /**
   * Creating all the monitoring histograms:
   * @hChi2Map_Chi2_pT        Chi^2(pT)
   * @hSignalAreaScaling      Signal area scaling (pT)
   * @hCorrbackAreaScaling    Corr. bkg. area scaling (pT)
   * @h_x_min                 signal scaling factor (pT)
   * @h_y_min                 corr. bkg. scaling dactor(pT)
   * @hErrXlow                lower signal scaling factor uncertainty (pT)
   * @hErrXhigh               upper signal scaling factor uncertainty (pT)
   * @hErrYlow                lower corr. bkg. scaling factor uncertainty (pT)
   * @hErrYhigh               upper corr. bkg. scaling factor uncertainty (pT)
   */
  TH1D* hChi2Map_Chi2_pT = new TH1D("hChi2Map_Chi2_pT", "", numberbins, fBinsPi013TeVEMCPt);
  hChi2Map_Chi2_pT->SetYTitle("#chi^{2}/ndf");
  hChi2Map_Chi2_pT->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hChi2Map_Chi2_pT->SetLineWidth(3);

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

  TH1D* hEfficiency = new TH1D("hEfficiency", "", numberbins, fBinsPi013TeVEMCPt);
  hEfficiency->SetYTitle("#epsilon_{rek}");
  hEfficiency->SetXTitle("#it{p}_{T} (GeV/#it{c})");

  for (int k = 1; k < numberbins; ++k) {
    hChi2Map_Chi2_pT->SetBinContent(k+1, vChi2_DT_Chi2Map[k-1]/vNDF_DT_Chi2Map[k-1]);
    hChi2Map_Chi2_pT->SetBinError(k+1, sqrt(2./vNDF_DT_Chi2Map[k-1]));
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
    hEfficiency->SetBinContent(k+1, vMyEffi[k-1]);
    hEfficiency->SetBinError(k+1, vMyEffiUncer[k-1]);
  }

  /**
   * reopening the new root file(s) to safe all the related histograms in it.
   */
   switch (TEMPLATEMETHOD){
     case 1:
       OutputFile      = new TFile("OutputFileBetterBkgNN.root", "UPDATE");
       break;
     case 2:
       OutputFile      = new TFile("OutputFileNormal.root", "UPDATE");
       break;
     case 3:
       OutputFile      = new TFile("OutputFileBetterBkg3to8_LowerRebinning.root", "UPDATE");
       break;
     case 4:
       OutputFile      = new TFile("OutputFileBetterBkg3to8_HigherRebinning.root", "UPDATE");
       break;
     case 5:
       OutputFile      = new TFile("OutputFileBetterBkg3to8_UncorrBckSmallFitRange.root", "UPDATE");
       break;
     case 6:
       OutputFile      = new TFile("OutputFileBetterBkg3to8_UncorrBckWideFitRange.root", "UPDATE");
       break;
     case 7:
       OutputFile      = new TFile("OutputFileBetterBkg3to8_SmallFit.root", "UPDATE");
       break;
     case 8:
       OutputFile      = new TFile("OutputFileBetterBkg3to8_HigherFit.root", "UPDATE");
       break;
     case 9:
       OutputFile      = new TFile("OutputFileBetterBkg3to8_SmallInt.root", "UPDATE");
       break;
     case 10:
       OutputFile      = new TFile("OutputFileBetterBkg3to8_HigherInt.root", "UPDATE");
       break;
     default:
       OutputFile      = new TFile("OutputFileBetterBkg3to8.root", "UPDATE");
       break;
   }


  /**
   * Open the file which contains the corrected true norm efficiency Yield from
   * the framework.
   */
  TFile* FData_corrected = SafelyOpenRootfile(strDataCorrection.Data());
  if (FData_corrected->IsOpen() ) printf("FData_corrected opened successfully\n");

  TFile* FCorrectedMC = SafelyOpenRootfile(strMCCorrectionHistos.Data());
  hPi0_gen = (TH1D*) FCorrectedMC->Get(Form("MC_Meson_genPt"));

  CorrectedYieldNormEff = (TH1D*) FData_corrected->Get(Form("CorrectedYieldNormEff"));

  /**
  * Open the data file which contains the uncorrected Yield aswell as the Chi^2
  * from the function parametrisation method from the framework.
   */
  DataFile = NULL;
  DataFile = SafelyOpenRootfile(strDataWithoutCorrection.Data());
  if (DataFile->IsOpen() ) printf("DataFile opened successfully\n");

  TH1D* hYield_framework = (TH1D*) DataFile->Get(Form("histoYieldMeson"));
  TH1D* histoChi2_0 = (TH1D*) DataFile->Get(Form("histoChi2_0"));
  SetHistoStandardSettings(histoChi2_0);

  /**
   * correcting the yields with division by pT (bincenter)
   * then correcting Yields with the acceptance
   * then with the efficiency.
   * i starts at 2 since 1 would be 0 <= pT (GeV/c) < 1.4!!!
   */
  for (int i = 2; i <= numberbins; i++) {
    hYield_dt_chi2map_uncorr->SetBinContent(i,hYield_dt_chi2map_uncorr->GetBinContent(i)/hYield_dt_chi2map_uncorr->GetBinCenter(i));
    hYield_dt_chi2map_uncorr->SetBinError(i,hYield_dt_chi2map_uncorr->GetBinError(i)/hYield_dt_chi2map_uncorr->GetBinCenter(i));
    hYield_framework->SetBinContent(i,hYield_framework->GetBinContent(i)/hYield_framework->GetBinCenter(i));
    hYield_framework->SetBinError(i,hYield_framework->GetBinError(i)/hYield_framework->GetBinCenter(i));
    hYieldMC->SetBinContent(i,hYieldMC->GetBinContent(i)/hYieldMC->GetBinCenter(i));
    hYieldMC->SetBinError(i,hYieldMC->GetBinError(i)/hYieldMC->GetBinCenter(i));
    hPi0_gen->SetBinContent(i,hPi0_gen->GetBinContent(i)/hPi0_gen->GetBinCenter(i));
    hPi0_gen->SetBinError(i,hPi0_gen->GetBinError(i)/hPi0_gen->GetBinCenter(i));
  }

  // correction for 2pi, BR, NEvents_data, Y, Binwidth
  /**
   * strMCCorrectionHistos.Data() for the number of Events
   * 2*Pi
   * the rapidity
   * the branching ratio for pi0 to decay into two photons
   * the bin width
   */
  hYield_dt_chi2map_uncorr->Scale(1./(NEvents_data*2*M_PI*1.6*0.98798),"width");
  hYield_dt_chi2map_uncorr->SetYTitle(rawyield);
  hYield_dt_chi2map_uncorr->SetXTitle(pt_str);

  /**
   * normalize the uncorrected framework yield.
   */
  hYield_framework->Scale(1./(NEvents_data*2*M_PI*1.6*0.98798));
  hPi0_gen->Scale(1./(NEvents_MC*2*M_PI*1.6*0.98798));
  hYieldMC->Scale(1./(NEvents_data*2*M_PI*1.6*0.98798), "width");
  hYieldMC->SetYTitle(rawyield);
  hYieldMC->SetXTitle(pt_str);


  /**
   * correcting Yields with the acceptance.
   */
  TH1D * hYield_dt_chi2map_acceptance_corrected =  (TH1D*) hYield_dt_chi2map_uncorr->Clone("hYield_dt_chi2map_acceptance_corrected");
  hYield_dt_chi2map_acceptance_corrected->Divide(hYield_dt_chi2map_uncorr, hAcc, 1, 1);
  hYieldMC_acc_corr = (TH1D*) hYieldMC->Clone("hYieldMC_acc_corr");
  hYieldMC_acc_corr->Divide(hYieldMC, hAcc, 1, 1);

  /**
   * correcting Yields with the efficiency.
   */
  TH1D * hYield_dt_chi2map_corrected =  (TH1D*) hYield_dt_chi2map_acceptance_corrected->Clone("hYield_dt_chi2map_corrected");
  hYieldMC_effi_corr = (TH1D*) hYieldMC_acc_corr->Clone("hYieldMC_acc_corr");

  hYield_dt_chi2map_corrected->Divide(hYield_dt_chi2map_acceptance_corrected, hEfficiency, 1, 1);
  hYieldMC_effi_corr->Divide(hYield_dt_chi2map_acceptance_corrected, hEfficiency, 1 , 1);
  // for (int i = 2; i <= numberbins; i++) {
  //   hYield_dt_chi2map_corrected->SetBinContent(
  //     i,hYield_dt_chi2map_acceptance_corrected->GetBinContent(i)/hEfficiency->GetBinContent(i)
  //   );
  //   hYield_dt_chi2map_corrected->
  //   SetBinError(i,
  //     sqrt(
  //       pow(hYield_dt_chi2map_acceptance_corrected->GetBinError(i)/hEfficiency->GetBinContent(i), 2.) +
  //       pow(hYield_dt_chi2map_acceptance_corrected->GetBinContent(i)*hEfficiency->GetBinError(i)/
  //       pow(hEfficiency->GetBinContent(i),  2.), 2.)));
  //   hYieldMC_effi_corr->SetBinContent(i,hYieldMC_effi_corr->GetBinContent(i)/hEfficiency->GetBinContent(i));
  // }

  hYield_dt_chi2map_corrected->SetYTitle(strCorrectedYield);

  std::cout << "2nd bin ratio:" <<  hYield_dt_chi2map_corrected->GetBinContent(3)/CorrectedYieldNormEff->GetBinContent(3) << '\n';


  /**
   * wrinting all the wanted histograms for plotting purposes in the output
   * file. part 2
   */
  hYield_dt_chi2map_uncorr->              Write("hYield_dt_chi2map_uncorr");
  hYield_dt_chi2map_acceptance_corrected->Write("hYield_dt_chi2map_acceptance_corrected");
  hYield_dt_chi2map_corrected->           Write("hYield_dt_chi2map_corrected");
  hYield_framework->                      Write("hYield_framework");
  CorrectedYieldNormEff->                 Write("hCorrectedYieldNormEff");
  hYieldMC->                              Write("hYieldMC");
  hYieldMC_acc_corr->                     Write("hYieldMC_acc_corr");
  hYieldMC_effi_corr->                    Write("hYieldMC_effi_corr");
  hChi2Map_Chi2_pT->                      Write("hChi2Map_Chi2_pT");
  hSignalAreaScaling->                    Write("hSignalAreaScaling");
  hCorrbackAreaScaling->                  Write("hCorrbackAreaScaling");
  h_x_min->                               Write("h_x_min");
  h_y_min->                               Write("h_y_min");
  hErrXlow->                              Write("hErrXlow");
  hErrXhigh->                             Write("hErrXhigh");
  hErrYlow->                              Write("hErrYlow");
  hErrYhigh->                             Write("hErrYhigh");
  histoChi2_0->                           Write("histoChi2_0");
  hEfficiency->                           Write("hEfficiency");
  hEffi->                                 Write("TrueMesonEffiPt");
  hAcc->                                  Write("hAcc");
  hPi0_gen->                              Write("MC_Meson_genPt");

  /**
   * Garbage collection part 2.
   */
  delete hYield_dt_chi2map_uncorr;
  delete hYield_dt_chi2map_acceptance_corrected;
  delete hChi2Map_Chi2_pT;
  delete hSignalAreaScaling;
  delete hCorrbackAreaScaling;
  delete h_x_min;
  delete h_y_min;
  delete hErrXlow;
  delete hErrXhigh;
  delete hErrYlow;
  delete hErrYhigh;
  delete hYield_framework;
  delete histoChi2_0;
  delete hEfficiency;

  /**
   * clearing all the vectors and freeing memory. Maybe not needed. I dunno
   */
  vMyEffi.clear();
  vMyEffiUncer.clear();
  vChi2_DT_Chi2Map.clear();
  vNDF_DT_Chi2Map.clear();
  vSignalAreaScaling.clear();
  vCorrbackAreaScaling.clear();
  v_x_min.clear();
  v_y_min.clear();
  vsigma_dt.clear();

  if(TEMPLATEMETHOD == 0){
    std::cout << "\nSTART SYSTEMATIC ANALYSIS\n\n" << '\n';
    systematics(OutputFile);
    std::cout << "\n\nFINISHED SYSTEMATIC ANALYSIS" << '\n';
  }

  /**
   * Closing all the files which were opend for the Yields.
   */
  // CorrectionFile->Close();
  // OutputFile->Close();
  // DataFile->Close();
  // FData_corrected->Close();

}
