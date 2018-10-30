#include "BackGroundFitting.h"
#include "chi2test.h"

/**
 * [the main function]
 * @param  current_path   [needed for the function to be calleable for variable datat sets]
 * @param  templatemethod [== 1 uses the backgrund templates from the 3 to 8 method,
 *                         == 2 uses the backgrund templates from the Next Neighbours method]
 *                         == 3 uses the fPulse function on the 3 to 8 Output to refine it.
 *                         only useable AFTER 3 to 8 root file was created!
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
  TFile* CorrectionFile           = NULL; // File which contains the correction-
                                          // histograms for effi and acceptance

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
  TH1D* CorrectedYieldNormEff     = NULL; // Corrected Yield from the Framwork

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
   * Calls the CorrBkgCreation function from BackGroundFitting.h This function
   * should make the needed normal corr. bkg templates for the other two
   * methods.
   */
  CorrBkgCreation();

  /**
   * For loop to loop over all pT bins definded by fBinsPi013TeVEMCPt. Bin 0 is
   * excluded since it only contains 0 <= pT (GeV/c) < 1.4 where no data should
   * be present, since we have an energy cut at 0.7 GeV per Cluster.
   */
  for (int k = 1; k < numberbins; ++k) {

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
    hCorrBkg      = NULL;

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
    if(templatemethod == 1 || templatemethod == 3){
      hCorrBkg = NULL;
      hCorrBkg = BackGround3to8(k);
    }

    else if(templatemethod == 2){
      hCorrBkg = BackgroundAdding(k);
    }
    else{
      std::cerr << "templatemethod not found!" << '\n';
      exit(2);
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
     * @param templatemethod      telling the function which method is currently
     * used
     * @param fBinsPi013TeVEMCPt  pT binning
     * @param k                   current PT bin
     */
    hChi2Map[k-1] = Chi2MapFunction(hInvMass_Data, hPeak_MC, hCorrBkg, temp_chi2_dt,
      signalAreaScaling, corrbackAreaScaling, x_min, y_min, ndf, templatemethod,
      fBinsPi013TeVEMCPt[k], k);

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
      if(templatemethod == 2){
        OutputFile      = new TFile("OutputFileBetterBkgNNforAdrian.root", "RECREATE");
      }
      else if(templatemethod == 1){
        OutputFile      = new TFile("OutputFileBetterBkg3to8.root", "RECREATE");
      }
      else if(templatemethod == 3){
        OutputFile      = new TFile("OutputFileBetterBkgPulse.root", "RECREATE");
      }
      else{
        std::cerr << "templatemethod not found!" << '\n';
        exit(1);
      }
    }

    else{
      if(templatemethod == 2){
        OutputFile      = new TFile("OutputFileBetterBkgNNforAdrian.root", "UPDATE");
      }
      else if(templatemethod == 1){
        OutputFile      = new TFile("OutputFileBetterBkg3to8.root", "UPDATE");
      }
      else if(templatemethod == 3){
        OutputFile      = new TFile("OutputFileBetterBkgPulse.root", "UPDATE");
      }
      else{
        std::cerr << "templatemethod not found!" << '\n';
        exit(1);
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
    hChi2Map[k-1]->       Write(Form("hChi2_2Dbin%02d",k));
    hChi2Map_sigma[k-1]-> Write(Form("hChi2_2D_sigma_bin%02d",k));
    hPeak_MC->            Write(Form("hSignal_bin%02d",k));
    hCorrBkg->            Write(Form("hCorrBack_bin%02d",k));
    OutputFile->Close();

    gDirectory->Cd(safePath.Data()); // resetting directory again.

    /**
     * Data histogram clone which will be integrated to get the uncorrected Yield
     * for the Chi2Map
     */
    data_clone_for_int_dt_chi2map = (TH1D*) hInvMass_Data->Clone("data_clone_for_int_dt_chi2map");

    hCorrBkg->Scale(y_min*corrbackAreaScaling);       //sclaing of the corr. bkg.

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
    int_value = data_clone_for_int_dt_chi2map->
    IntegralAndError(data_clone_for_int_dt_chi2map->GetXaxis()->
    FindBin(lowercountrange[k]),
    data_clone_for_int_dt_chi2map->GetXaxis()->FindBin(uppercountrange),
    int_error);

    hYield_dt_chi2map_uncorr->SetBinContent(k+1, int_value);
    hYield_dt_chi2map_uncorr->SetBinError(k+1, int_error);

    /**
     * Garbage collection part 1.
     */
    hPeak_MC      = NULL;
    hCorrBkg      = NULL;
    hInvMass_Data = NULL;
    hChi2Map[k-1] = NULL;

    MCFile->Close();
    DataFile->Close();

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
  }

  /**
   * reopening the new root file(s) to safe all the related histograms in it.
   */
  if(templatemethod == 2){
    OutputFile      = new TFile("OutputFileBetterBkgNNforAdrian.root", "UPDATE");
  }
  else if(templatemethod == 1){
    OutputFile      = new TFile("OutputFileBetterBkg3to8.root", "UPDATE");
  }
  else if(templatemethod == 3){
    OutputFile      = new TFile("OutputFileBetterBkgPulse.root", "UPDATE");
  }
  else{
    std::cerr << "templatemethod not found!" << '\n';
    exit(2);
  }


  /**
   * Open the file which contains the corrected true norm efficiency Yield from
   * the framework.
   */
  TFile* FData_corrected = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1Correction_00010113_1111112067032220000_01631031000000d0.root");
  if (FData_corrected->IsOpen() ) printf("FData_corrected opened successfully\n");

  CorrectedYieldNormEff = (TH1D*) FData_corrected->Get(Form("CorrectedYieldNormEff"));


  /**
  * Open the MC file which contains the correction histograms for efficiency and
  * acceptance. Obtaining those two directly afterwards.
   */
  CorrectionFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_MC_GammaConv_OnlyCorrectionFactor_00010113_1111112067032220000_01631031000000d0.root");
  if (CorrectionFile->IsOpen() ) printf("CorrectionFile opened successfully\n");



  TH1D* hAcc    = (TH1D*) CorrectionFile->Get(Form("fMCMesonAccepPt"));
  TH1D* hEffi   = (TH1D*) CorrectionFile->Get(Form("TrueMesonEffiPt"));

  // correction for 2pi, BR, NEvents, Y, Binwidth
  /**
   * Correction for the number of Events
   * 2*Pi
   * the rapidity
   * the branching ratio for pi0 to decay into two photons
   * the bin width
   */
  hYield_dt_chi2map_uncorr->Scale(1./(NEvents*2*M_PI*1.6*0.98798),"width");
  hYield_dt_chi2map_uncorr->SetYTitle(rawyield);
  hYield_dt_chi2map_uncorr->SetXTitle(pt_str);


  /**
  * Open the data file which contains the uncorrected Yield aswell as the Chi^2
  * from the function parametrisation method from the framework.
   */
  DataFile = NULL;
  DataFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");
  if (DataFile->IsOpen() ) printf("DataFile opened successfully\n");

  TH1D* hYield_framework = (TH1D*) DataFile->Get(Form("histoYieldMeson"));
  TH1D* histoChi2_0 = (TH1D*) DataFile->Get(Form("histoChi2_0"));
  SetHistoStandardSettings(histoChi2_0);

  /**
   * correcting the uncorrected framework yield with the efficiency.
   */
  hYield_framework->Divide(hEffi);

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
  }

  /**
   * correcting Yields with the acceptance.
   */
  TH1D * hYield_dt_chi2map_acceptance_corrected =  (TH1D*) hYield_dt_chi2map_uncorr->Clone("hYield_dt_chi2map_acceptance_corrected");
  hYield_dt_chi2map_acceptance_corrected->Divide(hAcc);

  /**
   * correcting Yields with the efficiency.
   */
  TH1D * hYield_dt_chi2map_corrected =  (TH1D*) hYield_dt_chi2map_acceptance_corrected->Clone("hYield_dt_chi2map_corrected");

  for (int i = 2; i <= numberbins; i++) {
    hYield_dt_chi2map_corrected->SetBinContent(i,hYield_dt_chi2map_acceptance_corrected->GetBinContent(i)/vInIntRangePercent[i-2]);
    hYield_dt_chi2map_corrected->SetBinError(i,hYield_dt_chi2map_acceptance_corrected->GetBinError(i)/vInIntRangePercent[i-2]);
  }

  hYield_dt_chi2map_corrected->SetYTitle(strCorrectedYield);


  /**
   * wrinting all the wanted histograms for plotting purposes in the output
   * file. part 2
   */
  hYield_dt_chi2map_uncorr->              Write("hYield_dt_chi2map_uncorr");
  hYield_dt_chi2map_acceptance_corrected->Write("hYield_dt_chi2map_acceptance_corrected");
  hYield_dt_chi2map_corrected->           Write("hYield_dt_chi2map_corrected");
  hYield_framework->                      Write("hYield_framework");
  CorrectedYieldNormEff->                 Write("hCorrectedYieldNormEff");
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

  /**
   * clearing all the vectors and freeing memory. Maybe not needed. I dunno
   */
  vInIntRangePercent.clear();
  vChi2_DT_Chi2Map.clear();
  vNDF_DT_Chi2Map.clear();
  vSignalAreaScaling.clear();
  vCorrbackAreaScaling.clear();
  v_x_min.clear();
  v_y_min.clear();
  vsigma_dt.clear();

  /**
   * Closing all the files which were opend for the Yields.
   */
  CorrectionFile->Close();
  OutputFile->Close();
  DataFile->Close();
  FData_corrected->Close();

}
