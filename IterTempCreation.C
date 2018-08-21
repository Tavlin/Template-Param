// #include "CommonHeader.h"
#include "chi2test.h"                         // self made chi2 calc and mapping
#include "TFractionFitter.h"

////////////////////////////////////////////////////////////////////////////////
// Function for Double Template Param
Double_t mc_full_func1(Double_t *x,  Double_t *par){
  Double_t xx = x[0];
  return (Double_t) par[0]*mc_full_clone1->GetBinContent(mc_full->FindBin(xx)) +
  par[1]*korrBG_clone1->GetBinContent(korrBG->FindBin(xx));
}
////////////////////////////////////////////////////////////////////////////////
// Function for Double Template Param for drawing purposes since the histos from
// above get scaled
Double_t mc_full_func42(Double_t *x,  Double_t *par){
  Double_t xx = x[0];
  return (Double_t) par[0]*mc_full_clone42->GetBinContent(mc_full->FindBin(xx)) +
  par[1]*korrBG_clone42->GetBinContent(korrBG->FindBin(xx));
}

////////////////////////////////////////////////////////////////////////////////
// Function for Signal Template + Pol 1 Param
Double_t mc_full_func2(Double_t *x,  Double_t *par){
  Double_t xx = x[0];
  return (Double_t) par[0]*mc_full_clone2->GetBinContent(mc_full->FindBin(xx)) +
  par[1]+par[2]*xx;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Start of the Main
void IterTempCreation(std::string current_path){

  TString sPath = gDirectory->GetPath();            // retrieve neutral path

  //////////////////////////////////////////////////////////////////////////////
  // open ESD histo which is inside TLists inside a rootfile
  // There get the Histo with Number of Events (minimum Bias)
  TFile* ESDFile_MC    = SafelyOpenRootfile("./../Daten/" + current_path + ".root");
  if (ESDFile_MC->IsOpen() ) printf("ESDFile_MC opened successfully\n");


  TFile* ESDFile_data    = SafelyOpenRootfile("./../Daten/GammaCalo-data_503.root");
  if (ESDFile_data->IsOpen() ) printf("ESDFile_data opened successfully\n");

  TList* lGammaCalo_data        = (TList*) ESDFile_data->Get("GammaCalo_503");
  TList* lCutNumber_data        = (TList*) lGammaCalo_data->FindObject("Cut Number 00010113_1111112067032220000_01631031000000d0");
  TList* lESD_data              = (TList*) lCutNumber_data->FindObject("00010113_1111112067032220000_01631031000000d0 ESD histograms");

  TList* lGammaCalo_MC          = (TList*) ESDFile_MC->Get("GammaCalo_503");
  TList* lCutNumber_MC          = (TList*) lGammaCalo_MC->FindObject("Cut Number 00010113_1111112067032220000_01631031000000d0");
  TList* lESD_MC                = (TList*) lCutNumber_MC->FindObject("00010113_1111112067032220000_01631031000000d0 ESD histograms");
  TList* lMC_MC                 = (TList*) lCutNumber_MC->FindObject("00010113_1111112067032220000_01631031000000d0 MC histograms");
  TList* lTrue_MC               = (TList*) lCutNumber_MC->FindObject("00010113_1111112067032220000_01631031000000d0 True histograms");
  TH1D* hNEvents                = (TH1D*)  lESD_data->FindObject("NEvents");
  TH1D* hMC_Pi0InAcc_Pt         = (TH1D*)  lMC_MC->FindObject("MC_Pi0InAcc_Pt");
  TH2D* hTrueDoubleCounting_Pi0 = (TH2D*)  lTrue_MC->FindObject("ESD_TrueDoubleCountPi0_InvMass_Pt");
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

  std::vector<Double_t> vChi2_DT_Iter;              // vector containig Chi/ndf
                                                    // from the Double Temp Iter
                                                    // method

  std::vector<Double_t> a_dt_iter;                  // sclaing paramter for the
                                                    // signal for the double
                                                    // template Iter Method

  std::vector<Double_t> b_dt_iter;                  // sclaing paramter for the
                                                    // corr. Background for the
                                                    // double template Iter
                                                    // Method

  std::vector<Double_t> vChi2_DT_Iter_Selfcalc;      // vector containig Chi/ndf
                                                    // from the Double Temp Iter
                                                    // method but self
                                                    // calculated

  std::vector<Double_t> vChi2_Pol1_Iter;            // vector containig Chi/ndf
                                                    // from the Pol1 Iter method

  std::vector<Double_t> vChi2_DT_Iter_Test;          // vector containig Chi/ndf
                                                    // from Double Temp Iter
                                                    // via the Chi2Test Function
  // const int numberbins = 26;
  TFile *IterTemp;                                  // FilePointer where things
                                                    // will be safed to

  TH2D* hChi2_2D[numberbins];                       // Array of pointers to 2D
                                                    // histos, which contain the
                                                    // Chi2 map

  TH2D* hChi2_2D_sigma[numberbins];                 // Array of TH2D* which
                                                    // contain the 1sigma range

  TH1D* hChi2_dt_iter[numberbins];                  // Array of TH1D* which
                                                    // contain Chi2/ndf for the
                                                    // double temp IterMethod

  TH1D* ha_dt_iter[numberbins];                     // Array of TH1D* which
                                                    // which contain the signal
                                                    // scaling factor for the
                                                    // double temp IterMethod
                                                    // depending on the IterStep

  TH1D* hb_dt_iter[numberbins];                     // Array of TH1D* which
                                                    // which contain the corr
                                                    // background scaling factor
                                                    // for the double temp
                                                    // IterMethod depending on
                                                    // the IterStep

  TH1D* hChi2_dt_iter_selfcalc[numberbins];         // Array of TH1D* which
                                                    // contain Chi2/ndf for the
                                                    // double temp IterMethod
                                                    // but self calculated with
                                                    // with own chi2 function

  TH1D* hChi2_pol1_iter[numberbins];                // Array of TH1D* which
                                                    // contain Chi2/ndf for the
                                                    // Pol1 + temp IterMethod

  TH1D* hChi2_dt_iter_test[numberbins];             // Array of TH1D* which
                                                    // contain Chi2/ndf for the
                                                    // double temp with the
                                                    // Chi2Test Function

  //////////////////////////////////////////////////////////////////////////////
  // Histos conating the uncorrected yields from the different methods
  TH1D* hYield_dt_uncorr          = new TH1D("hYield_dt_uncorr",
                                            "", numberbins, fBinsPi013TeVEMCPt);

  TH1D* hYield_dt_chi2map_uncorr  = new TH1D("hYield_dt_chi2map_uncorr",
                                            "", numberbins, fBinsPi013TeVEMCPt);

  TH1D* hYield_pol1_uncorr        = new TH1D("hYield_pol1_uncorr",
                                            "", numberbins, fBinsPi013TeVEMCPt);

  //////////////////////////////////////////////////////////////////////////////
  // setting up the 2 Histograms to compare chi2 from the to fit methods as
  // well as peak factor comp. between pol 1 and double temp fit and the ratio
  // of BG. scaling factor and the Peak scaling factor

  TH1D* hChi2_DT_Iter             = new TH1D("hChi2_DT_Iter",
                                            "", numberbins, fBinsPi013TeVEMCPt);

  TH1D* hchi2_pol1                = new TH1D("hchi2_pol1",
                                            "", numberbins, fBinsPi013TeVEMCPt);

  TH1D* ha_DT                     = new TH1D("ha_DT",
                                            "", numberbins, fBinsPi013TeVEMCPt);

  TH1D* hb_DT                     = new TH1D("hb_DT",
                                            "", numberbins, fBinsPi013TeVEMCPt);

  TH1D* ha_pol1                   = new TH1D("ha_pol1",
                                            "", numberbins, fBinsPi013TeVEMCPt);

  TH1D* hpeakratio                = new TH1D("hpeakratio",
                                            "", numberbins, fBinsPi013TeVEMCPt);

  TH1D* hpeakcomp                 = new TH1D("hpeakcomp",
                                            "", numberbins, fBinsPi013TeVEMCPt);


  //////////////////////////////////////////////////////////////////////////////
  // Histos for the normalized pull plots
  TH1D* hRatioDoubleTemp;
  TH1D* hRatioDoubleTemp_chi2map;
  TH1D* hRatioPol1;


  TH1D* hDataCloneforCHi2;                          // Data Clone Histo for the
                                                    // Chi2Test Function Method

  TH1D* hDoubleTemp;                                // Final Double Temp Iter
                                                    // Method for the pT bins

  TH1D* hPol1;                                      // Final Pol1 + Temp Iter
                                                    // Method for the pT bins

  TH1D* hTrueDoubleCounting_Pi0_Pro;                // X-Projection of current
                                                    // pT for Double Counting

  //////////////////////////////////////////////////////////////////////////////
  // Data Clone Histos which will be Integrated to get the uncorrected Yield
  TH1D* data_clone_for_int_dt;
  TH1D* data_clone_for_int_dt_chi2map;
  TH1D* data_clone_for_int_pol1;

  //////////////////////////////////////////////////////////////////////////////
  // Standard Settings for all The TH1 histos
  SetHistoStandardSettings(hChi2_DT_Iter);
  SetHistoStandardSettings(hchi2_pol1);
  SetHistoStandardSettings(hpeakratio);
  SetHistoStandardSettings(hpeakcomp);
  SetHistoStandardSettings(ha_DT);
  SetHistoStandardSettings(hb_DT);
  SetHistoStandardSettings(ha_pol1);
  SetHistoStandardSettings(hYield_dt_uncorr);
  SetHistoStandardSettings(hYield_dt_chi2map_uncorr);
  SetHistoStandardSettings(hYield_pol1_uncorr);

  //////////////////////////////////////////////////////////////////////////////
  // Giving the Histos their unique colors depending on the method
  hYield_dt_uncorr->SetLineColor(kTeal-7);
  hYield_dt_uncorr->SetMarkerColor(kTeal-7);
  hYield_dt_chi2map_uncorr->SetLineColor(kMagenta+2);
  hYield_dt_chi2map_uncorr->SetMarkerColor(kMagenta+2);
  hYield_pol1_uncorr->SetLineColor(kRed);
  hYield_pol1_uncorr->SetMarkerColor(kRed);


  //////////////////////////////////////////////////////////////////////////////
  // open True Yield Path
  TFile* FData_corrected = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1Correction_00010113_1111112067032220000_01631031000000d0.root");
  if (FData_corrected->IsOpen() ) printf("FData_corrected opened successfully\n");

  TH1D* CorrectedYieldTrueEff = (TH1D*) FData_corrected->Get(Form("CorrectedYieldNormEff"));

  //////////////////////////////////////////////////////////////////////////////
  // going over all pt bins despite first one, which is some framework bs.
  for (int k = 1; k < numberbins; k++) {
    std::cout << "starte bin " << k << " reading and wrinting!" << std::endl << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    // open MC histo path
    TFile* MCFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_MC_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");
    if (MCFile->IsOpen() ) printf("MCFile opened successfully\n");


    ////////////////////////////////////////////////////////////////////////////
    // retrieve MC histograms
    data_MC = (TH1D*) MCFile->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02i",k));
    mc_full = (TH1D*) MCFile->Get(Form("Mapping_TrueMeson_InvMass_in_Pt_Bin%02i",k));


    TFile* DataFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");
    if (DataFile->IsOpen() ) printf("DataFile opened successfully\n");
    data = (TH1D*) DataFile->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02i",k));


    mc_full->GetXaxis()->SetRangeUser(0.,0.4);
    data->GetXaxis()->SetRangeUser(0.,0.4);

    hTrueDoubleCounting_Pi0_Pro = hTrueDoubleCounting_Pi0->ProjectionX(Form("hTrueDoubleCounting_Pi0_Pro_bin%02d",k),hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k]),
    hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k+1]));
    Double_t InIntRangeDoubleCounting = (Double_t)hTrueDoubleCounting_Pi0_Pro->Integral(hTrueDoubleCounting_Pi0_Pro->FindBin(lowerparamrange),
                                                   hTrueDoubleCounting_Pi0_Pro->FindBin(upperparamrange))/
                                                   (Double_t)mc_full->Integral(
                                                    mc_full->FindBin(lowerparamrange),
                                                    mc_full->FindBin(upperparamrange));

    Double_t InIntRangePercent = (Double_t)mc_full->Integral(mc_full->FindBin(lowerparamrange),
                                                   mc_full->FindBin(upperparamrange))/
                                                   (Double_t)hMC_Pi0InAcc_Pt->Integral(
                                                     hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k]),
                                                     hMC_Pi0InAcc_Pt->FindBin(fBinsPi013TeVEMCPt[k+1]));

    std::cout << "InIntRangePercent = " << InIntRangePercent << '\n';
    std::cout << "InIntRangeDoubleCounting = " << InIntRangeDoubleCounting << '\n';
    vInIntRangePercent.push_back(InIntRangePercent);

    ////////////////////////////////////////////////////////////////////////////
    // Getting the purposed corr Background
    korrBG = (TH1D*) data_MC->Clone("korrBG");
    korrBG->Add(mc_full,-1);

    ////////////////////////////////////////////////////////////////////////////
    // Function for ndf
    TF1* f_ChiOverNdf = new TF1("f_ChiOverNdf", "[0]", 0.0 ,100.);
    f_ChiOverNdf->SetNpx(ndrawpoints);
    f_ChiOverNdf->SetNumberFitPoints(numberbins);
    f_ChiOverNdf->SetLineColor(kBlue+2);
    f_ChiOverNdf->SetLineWidth(4);

    ////////////////////////////////////////////////////////////////////////////
    // using the correlated BG from MC as Template
    TF1* fit_eq_double_temp = new TF1("fit_eq_double_temp", &mc_full_func1, 0.0,0.4, 2);
    fit_eq_double_temp->SetNpx(ndrawpoints);
    fit_eq_double_temp->SetNumberFitPoints(numberbins);
    fit_eq_double_temp->SetLineColor(kTeal-7);
    fit_eq_double_temp->SetLineWidth(4);

    ////////////////////////////////////////////////////////////////////////////
    // second TF1 for drawing only since the other function gets scaled two
    // times since the external scaling of the templates and the internal
    // scaling
    TF1* fit_eq_double_temp42 = new TF1("fit_eq_double_temp42", &mc_full_func42, 0.0,0.4, 2);
    fit_eq_double_temp42->SetNpx(ndrawpoints);
    fit_eq_double_temp42->SetNumberFitPoints(numberbins);
    fit_eq_double_temp42->SetLineColor(kTeal-7);
    fit_eq_double_temp42->SetLineWidth(4);

    //////////////////////////////////////////////////////////////////////////
    // Fix! Changes < in TLatex to #leq
    str = data_MC->GetTitle();
    TString str_copy = str.Copy();
    str_copy.ReplaceAll("<","#leq");
    str.Replace(0,20,str_copy,23);

    SetHistoStandardSettings(data_MC);
    SetHistoStandardSettings(mc_full);
    SetHistoStandardSettings(korrBG);
    SetHistoStandardSettings(korrBG);

    data->SetTitle(str);

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
    data->GetYaxis()->SetTitleOffset(1.2);
    data->GetYaxis()->SetLabelOffset(0.006);
    data->SetTitleSize(0.03, "xy");
    data->SetLabelSize(0.03, "xy");
    data->SetYTitle("d#it{N}/d#it{M}_{#gamma#gamma} (#it{c}^{2}/GeV)");
    data->SetMarkerStyle(24);
    data->SetMarkerSize(1.5);
    data->SetTitle("");
    data->SetLineWidth(3);
    korrBG->SetLineColor(kCyan+3);
    korrBG->SetMarkerColor(kCyan+3);
    korrBG->SetMarkerStyle(21);
    korrBG->SetMarkerSize(1.5);
    mc_full->SetLineColor(kGreen+3);
    mc_full->SetMarkerColor(kGreen+3);
    mc_full->SetMarkerStyle(33);
    mc_full->SetMarkerSize(2);

    //////////////////////////////////////////////////////////////////////////
    // clone for 2 temp fit
    TH1D* data_clone1 = (TH1D*) data->Clone("data_clone1");
    TH1D* data_clone4 = (TH1D*) data->Clone("data_clone4");

    //////////////////////////////////////////////////////////////////////////
    //clone for pol 1 fit
    TH1D* data_clone2 = (TH1D*) data->Clone("data_clone2");
    TH1D* data_clone3 = (TH1D*) data->Clone("data_clone3");

    // clearing the vectors
    vChi2_DT_Iter.clear();
    vChi2_Pol1_Iter.clear();
    vChi2_DT_Iter.resize(0);
    a_dt_iter.clear();
    a_dt_iter.resize(0);
    b_dt_iter.clear();
    b_dt_iter.resize(0);
    vChi2_Pol1_Iter.resize(0);
    vChi2_DT_Iter_Test.clear();
    vChi2_DT_Iter_Test.resize(0);
    vChi2_DT_Iter_Selfcalc.clear();
    vChi2_DT_Iter_Selfcalc.resize(0);
    iterMethodBool = 1;

    // resetting # of IterSteps to 0
    nIterStep = 0;


    while(iterMethodBool){

      //////////////////////////////////////////////////////////////////////////
      //clone for 2 temp fit
      mc_full_clone1  = (TH1D*) mc_full->Clone("mc_full_clone1");
      mc_full_clone42 = (TH1D*) mc_full->Clone("mc_full_clone42");
      korrBG_clone1   = (TH1D*) korrBG->Clone("korrBG_clone1");
      korrBG_clone42  = (TH1D*) korrBG->Clone("korrBG_clone42");

      //////////////////////////////////////////////////////////////////////////
      //clone for pol 1 fit
      mc_full_clone2  = (TH1D*) mc_full->Clone("mc_full_clone2");

      ////////////////////////////////////////////////////////////////////////////
      // creating the new root file to safe all the related histograms and fits
      // in it.
      if(k == 1 && nIterStep == 0){
        IterTemp      = new TFile("IterTemp.root", "RECREATE");
      }

      else{
        IterTemp      = new TFile("IterTemp.root", "UPDATE");
      }

      //////////////////////////////////////////////////////////////////////////
      // fit 2 temp
      TFitResultPtr r_double_temp1 = data_clone1->Fit("fit_eq_double_temp", "QM0PS","", lowerparamrange, upperparamrange);
      data_clone1->Fit("fit_eq_double_temp42", "QM0PS","", lowerparamrange, upperparamrange);
      vChi2_DT_Iter.push_back(r_double_temp1->Chi2() / r_double_temp1->Ndf());
      a_dt_iter.push_back(r_double_temp1->Parameter(0));
      b_dt_iter.push_back(r_double_temp1->Parameter(1));


      //////////////////////////////////////////////////////////////////////////
      //fit pol 1 + temp
      TFitResultPtr r_pol1_temp1 = data_clone2->Fit("fit_eq_1", "QM0PS","", lowerparamrange, upperparamrange);
      vChi2_Pol1_Iter.push_back(r_pol1_temp1->Chi2() / r_pol1_temp1->Ndf());

      //////////////////////////////////////////////////////////////////////////
      // scale 2 temp
      mc_full_clone1->Scale(r_double_temp1->Parameter(0));
      korrBG_clone1->Scale(r_double_temp1->Parameter(1));

      temp_ndf = r_double_temp1->Ndf();

      vChi2_DT_Iter_Selfcalc.push_back(chi2_selfmade(mc_full, korrBG, data,
                                                    temp_ndf,
                                                    r_double_temp1->Parameter(0),
                                                    r_double_temp1->Parameter(1))
                                                  /(Double_t)temp_ndf);


      //////////////////////////////////////////////////////////////////////////
      // Wrinting after the scaling:
      // mc_full_clone1->Write(Form("mc_full_clone_afterScaling_bin%02d_iter%d", k, nIterStep));
      // korrBG_clone1->Write(Form("korrBG_clone_afterScaling_bin%02d_iter%d", k, nIterStep));



      //////////////////////////////////////////////////////////////////////////
      // scale for pol 1 + temp
      mc_full_clone2->Scale(r_pol1_temp1->Parameter(0));

      //////////////////////////////////////////////////////////////////////////
      // test chi2 for monitoring
      gDirectory->Cd(sPath.Data());
      hDoubleTemp = (TH1D*) mc_full_clone1->Clone("hDoubleTemp");
      hDoubleTemp->Add(korrBG_clone1);
      hDataCloneforCHi2 = (TH1D*) data->Clone("hDataCloneforCHi2");

      for(int i = 0; i < 200; i++){
        if( i < 13 || i > 63){
          hDataCloneforCHi2->SetBinContent(i,0);
          hDataCloneforCHi2->SetBinError(i,0);
          hDoubleTemp->SetBinContent(i,0);
          hDoubleTemp->SetBinError(i,0);
        }
        ////////////////////////////////////////////////////////////////////////
        // calc the correct error for the Chi2Test Histo
        else{
          hDoubleTemp->SetBinError(i,sqrt(pow((r_double_temp1->Parameter(0)*mc_full->GetBinError(i)),2.)
          + pow((r_double_temp1->Parameter(1)*korrBG->GetBinError(i)),2.)));
        }
      }
      vChi2_DT_Iter_Test.push_back(hDataCloneforCHi2->Chi2Test(hDoubleTemp, "WW CHI2/NDF", 0));

      //////////////////////////////////////////////////////////////////////////
      // reset data_clone histos and then calculate their new errors

      gDirectory->Cd(sPath.Data());
      data_clone1 = (TH1D*) data->Clone("data_clone1");
      data_clone2 = (TH1D*) data->Clone("data_clone2");
      hDoubleTemp =( TH1D*) mc_full_clone1->Clone("hDoubleTemp");
      hDoubleTemp->Add(korrBG_clone1);
      for(int j = 13; j < 63; j++){
        data_clone1->SetBinError(j,sqrt(pow(data_clone1->GetBinError(j),2.)
        + pow((r_double_temp1->Parameter(0)*mc_full->GetBinError(j)),2.)
        + pow((r_double_temp1->Parameter(1)*korrBG->GetBinError(j)),2.)));

        data_clone2->SetBinError(j,sqrt(pow(data_clone2->GetBinError(j),2.)
        + pow((r_double_temp1->Parameter(0)*mc_full->GetBinError(j)),2.)));
      }

      IterTemp->Close();
      if(nIterStep >=1){
        if(fabs(vChi2_DT_Iter[nIterStep-1]-vChi2_DT_Iter[nIterStep] <= epsilon) &&
           fabs(vChi2_Pol1_Iter[nIterStep-1]-vChi2_Pol1_Iter[nIterStep] <= epsilon)){

          //////////////////////////////////////////////////////////////////////
          // for the last iteration step don't reset the clones, instead calc
          // errors for the data histos that will be used in the final Fit
          // afterwards

          for(int j = 0; j < 75; j++){
          data_clone4->SetBinError(j,sqrt(data_clone4->GetBinError(j) *
          data_clone4->GetBinError(j) + mc_full_clone1->GetBinError(j) *
          mc_full_clone1->GetBinError(j) + korrBG_clone1->GetBinError(j) *
          korrBG_clone1->GetBinError(j)));

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

    hChi2_dt_iter[k-1]          = new TH1D(Form("hChi2_dt_iter_bin%02d",k-1),"",
                                      nIterStep-1, 0.5, (Double_t) nIterStep-0.5);
    SetHistoStandardSettings(hChi2_dt_iter[k-1]);

    ha_dt_iter[k-1]             = new TH1D(Form("ha_dt_iter_bin%02d",k-1),"",
                                      nIterStep-1 ,0.5, (Double_t )nIterStep-0.5);
    SetHistoStandardSettings(ha_dt_iter[k-1]);

    hb_dt_iter[k-1]             = new TH1D(Form("hb_dt_iter_bin%02d",k-1),"",
                                      nIterStep-1, 0.5, (Double_t) nIterStep-0.5);
    SetHistoStandardSettings(hb_dt_iter[k-1]);


    hChi2_dt_iter_selfcalc[k-1] = new TH1D(Form("hChi2_dt_iter_selfcalc_bin%02d",k-1),"",
                                      nIterStep-2, 1.5, (Double_t) nIterStep-0.5);
    SetHistoStandardSettings(hChi2_dt_iter_selfcalc[k-1]);

    hChi2_pol1_iter[k-1]        = new TH1D(Form("1hChi2_pol1_iter_bin%02d",k-1),"",
                                      nIterStep-1, 0.5, (Double_t) nIterStep-0.5);
    SetHistoStandardSettings(hChi2_pol1_iter[k-1]);

    hChi2_dt_iter_test[k-1]     = new TH1D(Form("hChi2_dt_iter_test_bin%02d",k-1),"",
                                      nIterStep-1, 0.5, (Double_t) nIterStep-0.5);
    SetHistoStandardSettings(hChi2_dt_iter_test[k-1]);

    for(int i = 1; i < nIterStep; i++){
      hChi2_dt_iter[k-1]->SetBinContent(i, vChi2_DT_Iter[i-1]);
      ha_dt_iter[k-1]->SetBinContent(i, a_dt_iter[i-1]);
      hb_dt_iter[k-1]->SetBinContent(i, b_dt_iter[i-1]);
      hChi2_dt_iter_selfcalc[k-1]->SetBinContent(i, vChi2_DT_Iter_Selfcalc[i-1]);
      hChi2_pol1_iter[k-1]->SetBinContent(i, vChi2_Pol1_Iter[i-1]);
      hChi2_dt_iter_test[k-1]->SetBinContent(i, vChi2_DT_Iter_Test[i-1]);

      hChi2_dt_iter[k-1]->SetYTitle("#chi^{2}/ndf");
      hChi2_dt_iter[k-1]->SetXTitle("Iterationstep");
      hChi2_dt_iter[k-1]->SetLineColor(kTeal-7);
      hChi2_dt_iter[k-1]->SetMarkerColor(kTeal-7);

      hChi2_dt_iter_selfcalc[k-1]->SetYTitle("#chi^{2}/ndf");
      hChi2_dt_iter_selfcalc[k-1]->SetXTitle("Iterationstep");
      hChi2_dt_iter_selfcalc[k-1]->SetLineColor(kViolet+3);
      hChi2_dt_iter_selfcalc[k-1]->SetMarkerColor(kViolet+3);

      hChi2_pol1_iter[k-1]->SetYTitle("#chi^{2}/ndf");
      hChi2_pol1_iter[k-1]->SetXTitle("Iterationstep");
      hChi2_pol1_iter[k-1]->SetLineColor(kRed);
      hChi2_pol1_iter[k-1]->SetMarkerColor(kRed);

      hChi2_dt_iter_test[k-1]->SetYTitle("#chi^{2}/ndf");
      hChi2_dt_iter_test[k-1]->SetXTitle("Iterationstep");
      hChi2_dt_iter_test[k-1]->SetLineColor(kBlue+2);
      hChi2_dt_iter_test[k-1]->SetMarkerColor(kBlue+2);
    }

    gDirectory->Cd(sPath.Data());

    ///////////////////////////////////////////////////////////////////////////
    // final double temp fit
    mc_full_clone1 = (TH1D*) mc_full->Clone("mc_full_clone1");
    korrBG_clone1 = (TH1D*) korrBG->Clone("korrBG_clone1");



    TFitResultPtr r_double_temp = data_clone4->Fit("fit_eq_double_temp", "M0S","",lowerparamrange , upperparamrange);
    TH1D* mc_full_clone3 = (TH1D*) mc_full->Clone("mc_full_clone3");
    TH1D* korrBG_clone3 = (TH1D*) korrBG->Clone("korrBG_clone3");
    Double_t corrbackerror[101];
    for(int i = 0; i <= 100; i++){
      corrbackerror[i] = sqrt(pow(korrBG_clone3->GetBinError(i)
      *r_double_temp->Parameter(1),2.)+ pow(korrBG_clone3->GetBinContent(i)
      *r_double_temp->Error(1),2.));
    }
    mc_full_clone3->Scale(r_double_temp->Parameter(0));
    korrBG_clone3->Scale(r_double_temp->Parameter(1));

    ////////////////////////////////////////////////////////////////////////////
    //  making full histogram of the DoubleTemplate Param.
    hDoubleTemp = (TH1D*) mc_full_clone3->Clone("hDoubleTemp");
    hDoubleTemp->Add(korrBG_clone3);
    hDoubleTemp->SetMarkerStyle(20);
    hDoubleTemp->SetMarkerSize(1.5);
    hDoubleTemp->SetMarkerColor(kTeal-7);
    hDoubleTemp->SetLineColor(kTeal-7);
    ////////////////////////////////////////////////////////////////////////////

    for(int i = 0; i <= 100; i++){
      korrBG_clone3->SetBinError(i,corrbackerror[i]);
    }

    ///////////////////////////////////////////////////////////////////////////
    // final pol 1 + temp fit
    mc_full_clone2 = (TH1D*) mc_full->Clone("mc_full_clone2");
    TFitResultPtr r_pol1_temp = data_clone3->Fit("fit_eq_1", "M0S","", lowerparamrange, upperparamrange);
    TH1D* mc_full_clone4 = (TH1D*) mc_full->Clone("mc_full_clone4");
    mc_full_clone4->Scale(r_pol1_temp->Parameter(0));
    mc_full_clone4->SetLineColor(kRed);
    mc_full_clone4->SetMarkerColor(kRed);

    TF1* fpol1 = new TF1("fpol1", "[0]+x*[1]", 0.0, 0.4);
    fpol1->SetParameter(0, fit_eq_1->GetParameter(1));
    fpol1->SetParameter(1, fit_eq_1->GetParameter(2));
    fpol1->SetLineColor(kTeal-7);
    fpol1->SetLineWidth(3);

    ////////////////////////////////////////////////////////////////////////////
    //  making full histogram of the DoubleTemplate Param.
    hPol1 = (TH1D*) mc_full_clone4->Clone("hPol1");
    hPol1->Add(fpol1);
    hPol1->SetMarkerStyle(20);
    hPol1->SetMarkerSize(1.5);
    hPol1->SetMarkerColor(kRed+1);
    hPol1->SetLineColor(kRed+1);


    ///////////////////////////////////////////////////////////////////////////
    // Picture of double template fit with chi2 and factors as well as ratio of
    // the factors for the template.
    data_clone4->SetLineWidth(3);
    data_clone3->SetLineWidth(3);
    data_clone3->SetLineColor(kGray+2);
    data_clone3->SetMarkerColor(kGray+2);
    data_clone4->SetTitle("");
    fpol1->SetLineColor(kRed);

    hRatioDoubleTemp = (TH1D*) data_clone4->Clone("RatioDoubleTemp");
    hRatioPol1 = (TH1D*) data_clone3->Clone("RatioPol1");
    hRatioDoubleTemp->Add(fit_eq_double_temp,-1.);
    hRatioPol1->Add(fit_eq_1, -1.);
    for (int i = 0; i < 101; i++) {
      hRatioDoubleTemp->SetBinContent(i,hRatioDoubleTemp->GetBinContent(i)/hRatioDoubleTemp->GetBinError(i));
      hRatioPol1->SetBinContent(i,hRatioPol1->GetBinContent(i)/hRatioPol1->GetBinError(i));
      hRatioDoubleTemp->SetBinError(i,0);
      hRatioPol1->SetBinError(i,0);

    }
    hRatioDoubleTemp->GetYaxis()->SetRangeUser(-5.,5.);
    hRatioDoubleTemp->GetYaxis()->SetNdivisions(509);
    hRatioPol1->GetYaxis()->SetRangeUser(-4.5,4.5);
    SetHistoStandardSettings(hRatioDoubleTemp,0.,0.,0.09);
    SetHistoStandardSettings(hRatioPol1,0.,0.,0.09);
    hRatioDoubleTemp->SetYTitle("(data-param)/#sigma");
    hRatioPol1->SetYTitle("(data-param)/#sigma");
    hRatioDoubleTemp->GetYaxis()->SetTitleOffset(0.4);
    hRatioPol1->GetYaxis()->SetTitleOffset(0.4);
    hRatioDoubleTemp->SetMarkerStyle(24);
    hRatioDoubleTemp->SetMarkerSize(1.5);
    hRatioPol1->SetMarkerStyle(24);
    hRatioPol1->SetMarkerSize(1.5);

    hRatioDoubleTemp->SetLineColor(kTeal-7);
    hRatioDoubleTemp->SetMarkerColor(kTeal-7);
    hRatioPol1->SetLineColor(kRed);
    hRatioPol1->SetMarkerColor(kRed);

    ////////////////////////////////////////////////////////////////////////////
    // testitesti
    ////////////////////////////////////////////////////////////////////////////
    hChi2_2D[k-1] = chi2test(data, mc_full, korrBG, temp_chi2_dt, signalAreaScaling, corrbackAreaScaling, x_min, y_min, ndf);

    hChi2_2D_sigma[k-1] = getErrorHist(Form("hChi2_2D_sigma_bin%02",k), hChi2_2D[k-1],temp_chi2_dt+1);
    vChi2_DT_Chi2Map.push_back(temp_chi2_dt);
    vNDF_DT_Chi2Map.push_back(ndf);
    vSignalAreaScaling.push_back(signalAreaScaling);
    vCorrbackAreaScaling.push_back(corrbackAreaScaling);
    v_x_min.push_back(x_min);
    v_y_min.push_back(y_min);
    vsigma_dt.push_back(getErrors(hChi2_2D_sigma[k-1], x_min, y_min));

    f_ChiOverNdf->SetParameter(0,temp_chi2_dt/ndf);

    hRatioDoubleTemp_chi2map = (TH1D*) data_clone4->Clone("hRatioDoubleTemp_chi2map");
    for (int j = 0; j < 75; j++) {
        Double_t temp_error = sqrt(pow(x_min*mc_full->GetBinError(j), 2.)
        +pow(y_min*korrBG->GetBinError(j), 2.));

        Double_t chi2 = ((Double_t)x_min*mc_full->GetBinContent(j) + (Double_t)y_min*korrBG->GetBinContent(j)
        -data->GetBinContent(j))
        /sqrt((pow(temp_error,2.)+pow(data->GetBinError(j),2.)));
      hRatioDoubleTemp_chi2map->SetBinContent(j, chi2/(Double_t)ndf);
      }


    temp_chi2_dt = 0;





    data->SetTitle(str);
    IterTemp = new TFile("IterTemp.root","UPDATE");
    gDirectory = IterTemp;
    data->Write(Form("data_bin%02d",k));
    data_clone3->Write(Form("data_addedErrosPol1_bin%02d",k));
    data_clone4->Write(Form("data_addedErrosDT_bin%02d",k));
    mc_full_clone4->Write(Form("mc_peak_pol1_bin%02d",k));
    mc_full_clone3->Write(Form("mc_full_DT_bin%02d", k));
    korrBG_clone3->Write(Form("korrBG_bin%02d", k));
    fpol1->Write(Form("fpol1_bin%02d",k));
    hRatioDoubleTemp->Write(Form("hRatioDoubleTemp_bin%02d", k));
    hRatioDoubleTemp_chi2map->Write(Form("hRatioDoubleTemp_chi2map_bin%02d", k));
    hRatioPol1->Write(Form("hRatioPol1_bin%02d", k));
    hDoubleTemp->Write(Form("hDoubleTemp_bin%02d",k));
    hPol1->Write(Form("hPol1_bin%02d",k));
    hChi2_dt_iter[k-1]->Write(Form("hChi2_dt_iter_bin%02d",k));
    ha_dt_iter[k-1]->Write(Form("ha_dt_iter_bin%02d",k));
    hb_dt_iter[k-1]->Write(Form("hb_dt_iter_bin%02d",k));
    hChi2_dt_iter_selfcalc[k-1]->Write(Form("hChi2_dt_iter_selfcalc_bin%02d",k));
    hChi2_pol1_iter[k-1]->Write(Form("hChi2_pol1_iter_bin%02d",k));
    hChi2_dt_iter_test[k-1]->Write(Form("hChi2_dt_iter_test_bin%02d",k));
    hChi2_2D[k-1]->Write(Form("hChi2_2Dbin%02d",k));
    hChi2_2D_sigma[k-1]->Write(Form("hChi2_2D_sigma_bin%02d",k));
    mc_full->Write(Form("hSignal_bin%02d",k));
    korrBG->Write(Form("hCorrBack_bin%02d",k));
    f_ChiOverNdf->Write(Form("f_ChiOverNdf%02d",k));

    gDirectory->Cd(sPath.Data());
    ////////////////////////////////////////////////////////////////////////////
    // getting the chi2 of the current pT bin
    hchi2_pol1->SetBinContent(k+1,r_pol1_temp->Chi2()/temp_ndf);
    hChi2_DT_Iter->SetBinContent(k+1,r_double_temp->Chi2()/temp_ndf);



    data_clone_for_int_dt = (TH1D*) data->Clone("hYield_dt_uncorr");
    data_clone_for_int_dt->Add(korrBG_clone3, -1);
    int_value = data_clone_for_int_dt->IntegralAndError(data_clone_for_int_dt->GetXaxis()->FindBin(0.085), data_clone_for_int_dt->GetXaxis()->FindBin(0.225), int_error);
    hYield_dt_uncorr->SetBinContent(k+1, int_value);
    hYield_dt_uncorr->SetBinError(k+1, int_error);

    data_clone_for_int_dt_chi2map = (TH1D*) data->Clone("data_clone_for_int_dt_chi2map");
    korrBG->Scale(y_min);
    data_clone_for_int_dt_chi2map->Add(korrBG, -1);
    int_value = data_clone_for_int_dt_chi2map->IntegralAndError(data_clone_for_int_dt_chi2map->GetXaxis()->FindBin(0.085), data_clone_for_int_dt_chi2map->GetXaxis()->FindBin(0.225), int_error);
    hYield_dt_chi2map_uncorr->SetBinContent(k+1, int_value);
    hYield_dt_chi2map_uncorr->SetBinError(k+1, int_error);


    data_clone_for_int_pol1 = (TH1D*) data->Clone("hYield_pol1_uncorr");
    data_clone_for_int_pol1->Add(fpol1, -1);
    int_value = data_clone_for_int_pol1->IntegralAndError(data_clone_for_int_pol1->GetXaxis()->FindBin(0.085), data_clone_for_int_pol1->GetXaxis()->FindBin(0.225), int_error);
    hYield_pol1_uncorr->SetBinContent(k+1, int_value);
    hYield_pol1_uncorr->SetBinError(k+1, int_error);

    ////////////////////////////////////////////////////////////////////////////
    // getting the peakratio of the current pT bin
    Double_t peakratio, peakratioerr, peakscale, bgscale, peakerr, bgerr;
    Double_t peakscale_pol1, peakerr_pol1, peakratio_to_pol1, peakratioerr_to_pol1;
    peakscale = r_double_temp->Parameter(0);
    bgscale = r_double_temp->Parameter(1);
    peakerr = r_double_temp->Error(0);
    bgerr = r_double_temp->Error(1);

    peakscale_pol1 = r_pol1_temp->Parameter(0);
    peakerr_pol1 = r_pol1_temp->Error(0);

    peakratio = bgscale/peakscale;
    peakratio_to_pol1 = peakscale/peakscale_pol1;

    peakratioerr = sqrt(pow(bgerr/peakscale, 2.)+pow(bgscale*peakerr/pow(peakscale, 2.), 2.));
    peakratioerr_to_pol1 = sqrt(pow(peakerr/peakscale_pol1,2.0)+pow(peakscale*peakerr_pol1/pow(peakscale_pol1,2.0),2.0));

    hpeakratio->SetBinContent(k+1, peakratio);
    hpeakratio->SetBinError(k+1, peakratioerr);

    hpeakcomp->SetBinContent(k+1, peakratio_to_pol1);
    hpeakcomp->SetBinError(k+1, peakratioerr_to_pol1);

    ha_DT->SetBinContent(k+1, peakscale);
    ha_DT->SetBinError(k+1, peakerr);
    hb_DT->SetBinContent(k+1, bgscale);
    hb_DT->SetBinError(k+1, bgerr);
    ha_pol1->SetBinContent(k+1, peakscale_pol1);
    ha_pol1->SetBinError(k+1, peakerr_pol1);



    ////////////////////////////////////////////////////////////////////////////
    // garbage collection part 1
    delete mc_full_clone1;
    delete mc_full_clone2;
    delete mc_full_clone3;
    delete mc_full_clone4;
    delete korrBG_clone1;
    delete korrBG_clone3;
    delete data_clone1;
    delete data_clone2;
    delete data_clone3;
    delete data_clone4;
    delete mc_full;
    delete korrBG;
    delete data;
    delete fpol1;
    delete fit_eq_double_temp;
    delete f_ChiOverNdf;
    delete fit_eq_1;

    if(k < numberbins-1){
      // IterTemp->Close();
      MCFile->Close();
      DataFile->Close();
    }
    IterTemp->Close();
    std::cout << "bin number" << k << "Ende" << std::endl << std::endl;
    delete hChi2_dt_iter[k-1];
    delete ha_dt_iter[k-1];
    delete hb_dt_iter[k-1];
    delete hChi2_pol1_iter[k-1];
    delete hChi2_dt_iter_test[k-1];
    delete hChi2_2D[k-1];
  }

  TH1D* hChi2_DT_Chi2map = new TH1D("hChi2_DT_Chi2map", "", numberbins, fBinsPi013TeVEMCPt);
  hChi2_DT_Chi2map->SetYTitle("#chi^{2}/ndf");
  hChi2_DT_Chi2map->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hChi2_DT_Chi2map->SetLineWidth(3);

  TH1D* hSignalAreaScaling = new TH1D("hSignalAreaScaling", "", numberbins, fBinsPi013TeVEMCPt);
  hSignalAreaScaling->SetYTitle("signal areascaling factor");
  hSignalAreaScaling->SetXTitle("#it{p}_{T} (GeV/#it{c})");

  TH1D* hCorrbackAreaScaling = new TH1D("hCorrbackAreaScaling", "", numberbins, fBinsPi013TeVEMCPt);
  hCorrbackAreaScaling->SetYTitle("corr. back. areascaling factor");
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


  /////////////////////////////////////////////////////////////////////////////
  // Writing of Chi2(pT)
  IterTemp = new TFile("IterTemp.root","UPDATE");

  hchi2_pol1->SetLineColor(kRed);
  hchi2_pol1->SetLineWidth(3);
  hChi2_DT_Iter->SetLineWidth(3);
  hpeakratio->SetLineWidth(3);
  hpeakcomp->SetLineWidth(3);
  hchi2_pol1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hChi2_DT_Iter->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hchi2_pol1->GetYaxis()->SetTitle("#chi^{2}/ndf");
  hChi2_DT_Iter->GetYaxis()->SetTitle("#chi^{2}/ndf");
  hpeakratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hpeakratio->GetYaxis()->SetTitle("b_{double}/a_{double}");
  hpeakcomp->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hpeakcomp->GetYaxis()->SetTitle("a_{pol1}/ a_{double}");
  hChi2_DT_Iter->GetYaxis()->SetRangeUser(0.0,5.0);


  ////////////////////////////////////////////////////////////////////////////
  // open MC histo path
  TFile* CorrectionFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_MC_GammaConv_OnlyCorrectionFactor_00010113_1111112067032220000_01631031000000d0.root");
  if (CorrectionFile->IsOpen() ) printf("CorrectionFile opened successfully\n");

  TH1D* hAcc    = (TH1D*) CorrectionFile->Get(Form("fMCMesonAccepPt"));
  TH1D* hEffi   = (TH1D*) CorrectionFile->Get(Form("TrueMesonEffiPt"));


  //
  // //////////////////////////////////////////////////////////////////////////////
  // // sclaing/ correcting the efficienc
  // for (int i = 0; i < numberbins; i++) {
  //   hEffi->SetBinContent(i+2, hEffi->GetBinContent(i+2)*vInIntRangePercent[i]);
  //   hEffi->SetBinError(i+2, hEffi->GetBinError(i+2)*vInIntRangePercent[i]);
  // }

  hChi2_DT_Iter->Write(Form("hChi2_DT_Iter"));
  hchi2_pol1->Write(Form("hchi2_pol1"));
  hpeakratio->Write("hpeakratio");
  hpeakcomp->Write("hpeakcomp");
  ha_DT->Write("DoubleTemplatePeakFactor");
  hb_DT->Write("DoubleTemplatecorrBGFactor");
  ha_pol1->Write("Pol1PeakFactor");

  hYield_dt_uncorr->Scale(1./(NEvents*2*M_PI*1.6*0.98798),"width");
  // hYield_dt_uncorr->Scale(1, "width");
  // hYield_dt_uncorr->Scale(1./0.98798);      //branching ratio scaling
  hYield_dt_uncorr->SetYTitle(rawyield);
  hYield_dt_uncorr->SetXTitle(pt_str);

  hYield_pol1_uncorr->Scale(1./(NEvents*2*M_PI*1.6*0.98798),"width");
  // hYield_pol1_uncorr->Scale(1, "width");
  // hYield_pol1_uncorr->Scale(1./0.98798);      //branching ratio scaling
  hYield_pol1_uncorr->SetYTitle(rawyield);
  hYield_pol1_uncorr->SetXTitle(pt_str);

  hYield_dt_chi2map_uncorr->Scale(1./(NEvents*2*M_PI*1.6*0.98798),"width");
  // hYield_dt_chi2map_uncorr->Scale(1, "width");
  // hYield_dt_chi2map_uncorr->Scale(1./0.98798);      //branching ratio scaling
  hYield_dt_chi2map_uncorr->SetYTitle(rawyield);
  hYield_dt_chi2map_uncorr->SetXTitle(pt_str);



  TFile* DataFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");
  if (DataFile->IsOpen() ) printf("DataFile opened successfully\n");
  TH1D* hYield_framework = (TH1D*) DataFile->Get(Form("histoYieldMeson"));
  TH1D* histoChi2_0 = (TH1D*) DataFile->Get(Form("histoChi2_0"));
  SetHistoStandardSettings(histoChi2_0);
  if(hYield_framework) printf("histoYieldMeson ist da\n");
  else printf("histoYieldMeson ist nicht da\n");
  hYield_framework->Divide(hEffi);

  for (int i = 2; i <= numberbins; i++) {
    hYield_dt_uncorr->SetBinContent(i,hYield_dt_uncorr->GetBinContent(i)/hYield_dt_uncorr->GetBinCenter(i));
    hYield_pol1_uncorr->SetBinContent(i,hYield_pol1_uncorr->GetBinContent(i)/hYield_pol1_uncorr->GetBinCenter(i));
    hYield_dt_chi2map_uncorr->SetBinContent(i,hYield_dt_chi2map_uncorr->GetBinContent(i)/hYield_dt_chi2map_uncorr->GetBinCenter(i));
    hYield_dt_uncorr->SetBinError(i,hYield_dt_uncorr->GetBinError(i)/hYield_dt_uncorr->GetBinCenter(i));
    hYield_pol1_uncorr->SetBinError(i,hYield_pol1_uncorr->GetBinError(i)/hYield_pol1_uncorr->GetBinCenter(i));
    hYield_dt_chi2map_uncorr->SetBinError(i,hYield_dt_chi2map_uncorr->GetBinError(i)/hYield_dt_chi2map_uncorr->GetBinCenter(i));
    hYield_framework->SetBinContent(i,hYield_framework->GetBinContent(i)/hYield_framework->GetBinCenter(i));
    hYield_framework->SetBinError(i,hYield_framework->GetBinError(i)/hYield_framework->GetBinCenter(i));
  }

  //////////////////////////////////////////////////////////////////////////////
  // correcting yield for the acceptance
  TH1D * hYield_dt_acceptance_corrected         =  (TH1D*) hYield_dt_uncorr->Clone("hYield_dt_acceptance_corrected");
  hYield_dt_acceptance_corrected->Divide(hAcc);

  TH1D * hYield_pol1_acceptance_corrected       =  (TH1D*) hYield_pol1_uncorr->Clone("hYield_pol1_acceptance_corrected");
  hYield_pol1_acceptance_corrected->Divide(hAcc);

  TH1D * hYield_dt_chi2map_acceptance_corrected =  (TH1D*) hYield_dt_chi2map_uncorr->Clone("hYield_dt_chi2map_acceptance_corrected");
  hYield_dt_chi2map_acceptance_corrected->Divide(hAcc);

  //////////////////////////////////////////////////////////////////////////////
  // correcting yield for the efficiency
  TH1D * hYield_dt_corrected         =  (TH1D*) hYield_dt_acceptance_corrected->Clone("hYield_dt_corrected");
  // hYield_dt_corrected->Divide(hEffi);

  TH1D * hYield_pol1_corrected       =  (TH1D*) hYield_pol1_acceptance_corrected->Clone("hYield_pol1_corrected");
  // hYield_pol1_corrected->Divide(hEffi);

  TH1D * hYield_dt_chi2map_corrected =  (TH1D*) hYield_dt_chi2map_acceptance_corrected->Clone("hYield_dt_chi2map_corrected");
  // hYield_dt_chi2map_corrected->Divide(hEffi);

  for (int i = 2; i <= numberbins; i++) {
    hYield_dt_corrected->SetBinContent(i,hYield_dt_acceptance_corrected->GetBinContent(i)/vInIntRangePercent[i-2]);
    hYield_pol1_corrected->SetBinContent(i,hYield_pol1_acceptance_corrected->GetBinContent(i)/vInIntRangePercent[i-2]);
    hYield_dt_chi2map_corrected->SetBinContent(i,hYield_dt_chi2map_acceptance_corrected->GetBinContent(i)/vInIntRangePercent[i-2]);
    hYield_dt_corrected->SetBinError(i,hYield_dt_acceptance_corrected->GetBinError(i)/vInIntRangePercent[i-2]);
    hYield_pol1_corrected->SetBinError(i,hYield_pol1_acceptance_corrected->GetBinError(i)/vInIntRangePercent[i-2]);
    hYield_dt_chi2map_corrected->SetBinError(i,hYield_dt_chi2map_acceptance_corrected->GetBinError(i)/vInIntRangePercent[i-2]);
  }

  hYield_dt_corrected->SetYTitle(strCorrectedYield);
  hYield_pol1_corrected->SetYTitle(strCorrectedYield);
  hYield_dt_chi2map_corrected->SetYTitle(strCorrectedYield);



  hYield_dt_uncorr->Write("hYield_dt_uncorr");
  hYield_pol1_uncorr->Write("hYield_pol1_uncorr");
  hYield_dt_chi2map_uncorr->Write("hYield_dt_chi2map_uncorr");
  hYield_dt_acceptance_corrected->Write("hYield_dt_acceptance_corrected");
  hYield_pol1_acceptance_corrected->Write("hYield_pol1_acceptance_corrected");
  hYield_dt_chi2map_acceptance_corrected->Write("hYield_dt_chi2map_acceptance_corrected");
  hYield_dt_corrected->Write("hYield_dt_corrected");
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
  delete hpeakratio;
  delete hpeakcomp;
  delete hChi2_DT_Iter;
  delete hchi2_pol1;
  delete hRatioDoubleTemp;
  delete hRatioPol1;
  delete hYield_dt_uncorr;
  delete hYield_dt_chi2map_uncorr;
  delete hYield_pol1_uncorr;
  delete hYield_dt_acceptance_corrected;
  delete hYield_pol1_acceptance_corrected;
  delete hYield_dt_chi2map_acceptance_corrected;
  delete hDoubleTemp;
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
  vChi2_DT_Iter.clear();
  vChi2_Pol1_Iter.clear();
  vChi2_DT_Iter_Test.clear();
  IterTemp->Close();
  DataFile->Close();

}
