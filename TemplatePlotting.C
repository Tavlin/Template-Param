#include "Plotting_Patrick.h"
#include "CommonHeader.h"
//MCTemplatesAnData

// wpsid = which picture should I draw
void TemplatePlotting(TString wpsid = "all", TString PicFormat = "png",
                      const char* FrameworkOutput = "", const char*CutString = "" ,
                      TString SaveAppendix = ""){

  TGaxis::SetMaxDigits(3);
 /**
  * loop over templatemethod
  * templatemethod == 0 for 3 to 8 Method
  * templatemethod == 1 for NN method
  * templatemethod == 2 for normal
  */
  for (int templatemethod = 0; templatemethod <= 2; templatemethod++) {

    TString str;
    const Int_t nbins = numberbins;
    const Int_t ndrawpoints = 1.e5;
    Double_t somelist[2] = {1., 2.};

    //////////////////////////////////////////////////////////////////////////////
    // setting up the 2 Histograms to compare chi2 from the to fit methods as
    // well as peak factor comp. between pol 1 and double temp fit and the ratio
    // of BG. scaling factor and the Peak scaling factor
    auto OAhists = new TObjArray();
    auto OAratios = new TObjArray();
    TCanvas* c1                       = NULL;
    TCanvas* c2                       = NULL;
    TCanvas* cChi2Map                 = NULL;
    TH1D* hData                       = NULL;
    TH1D* hPol1PeakFactor             = NULL;
    TH1D* hYield_dt_chi2map_uncorr    = NULL;
    TH2D* hChi2_2D                    = NULL;
    TH2D* hChi2_2D_sigma              = NULL;
    TH1D* hChi2Map_Chi2_pT            = NULL;
    TH1D* histoChi2_0                 = NULL;
    TH1D* hSignalAreaScaling          = NULL;
    TH1D* hCorrbackAreaScaling        = NULL;
    TH1D* hYield_dt_chi2map_corrected = NULL;
    TH1D* hYield_pol1_corrected       = NULL;
    TH1D* hCorrectedYieldNormEff      = NULL;
    TH1D* h_x_min                     = NULL;
    TH1D* h_y_min                     = NULL;
    TH1D* hErrXlow                    = NULL;
    TH1D* hErrXhigh                   = NULL;
    TH1D* hErrYlow                    = NULL;
    TH1D* hErrYhigh                   = NULL;
    TH1D* hSignal                     = NULL;
    TH1D* hSignal_scaled              = NULL;
    TH1D* hAdded                      = NULL;
    TH1D* hCorrBack                   = NULL;
    TH1D* hCorrBack_scaled            = NULL;
    TH1D* hGG                         = NULL;
    TH1D* hGC                         = NULL;
    TH1D* hCC                         = NULL;
    TH1D* hInvMass_MC                 = NULL;
    TH1D* hEfficiency                 = NULL;
    TH1D* hAcc                        = NULL;
    TH1D* hCorrYieldME_StatError            = NULL;
    TH1D* hCorrectedYieldNormEff_StatError  = NULL;
    TH1D* hCorrYieldME_Ratio                = NULL;

    TH1D* hRatio_Bkg                  = NULL;
    TH1D* hCorrBackNoRebin            = NULL;
    TF1*  fPol0                       = NULL;
    TF1* f_ChiOverNdf                 = NULL;
    TLine* fitrange2                  = NULL;
    TLine* line_0                     = NULL;
    TLine* line_one                   = NULL;
    TLine* lChi2MinX                  = NULL;
    TLine* lChi2MinY                  = NULL;

    TLegend* legpT                    = NULL;
    TLegend* legTemplat               = NULL;
    TLegend* legTemplatChi2           = NULL;
    TLegend* legTemplateYield         = NULL;
    TLegend* legTemplateRatio         = NULL;
    TLegend* legpTRatio               = NULL;

    Double_t line_y = 0;

    TFile* MCWithOutFile = NULL;

    TString strCutString              = CutString;
    TString strFrameworkOutput        = FrameworkOutput;
    TString strMCWithoutCorrection    = strFrameworkOutput + strCutString + "/13TeV/Pi0_MC_GammaConvV1WithoutCorrection_" + strCutString + ".root";

    MCWithOutFile = SafelyOpenRootfile(strMCWithoutCorrection.Data());
    TFile* OutputFile   = NULL;
    TFile* BackFileNN   = NULL;
    TFile* BackFile3to8 = NULL;
    TFile* CorrBkgFile  = NULL;

    /*
      TLegend which displays most important info as header!
      There are so many different ones, since many pictues have different places
      where they can go.
     */
    TLegend* legSystem = new TLegend(0.1, 0.94, 0.7, 0.98);
    legSystem->AddEntry((TObject*) 0, "ALICE, pp #sqrt{#it{s}} = 13 TeV, #pi^{0} #rightarrow #gamma#gamma with EMCal", "");

    TLegend* legSystemChi2Map = new TLegend(0.03, 0.94, 0.6, 0.98);
    legSystemChi2Map->AddEntry((TObject*) 0, "ALICE, pp #sqrt{#it{s}} = 13 TeV, #pi^{0} #rightarrow #gamma#gamma with EMCal", "");

    legTemplat = new TLegend(0.09, 0.65, 0.38, 0.85);
    legTemplat->AddEntry((TObject*) 0, "templates:", "");
    legTemplat->AddEntry((TObject*) 0, "PYTHIA 8", "");
    legTemplat->AddEntry((TObject*) 0, "Monash 2013", "");
    legTemplat->AddEntry((TObject*) 0, "GEANT 3", "");

    legTemplatChi2 = new TLegend(0.2, 0.65, 0.6, 0.85);
    legTemplatChi2->AddEntry((TObject*) 0, "templates:", "");
    legTemplatChi2->AddEntry((TObject*) 0, "PYTHIA 8", "");
    legTemplatChi2->AddEntry((TObject*) 0, "Monash 2013", "");
    legTemplatChi2->AddEntry((TObject*) 0, "GEANT 3", "");

    legTemplateYield = new TLegend(0.5, 0.65, 0.8, 0.85);
    legTemplateYield->AddEntry((TObject*) 0, "templates:", "");
    legTemplateYield->AddEntry((TObject*) 0, "PYTHIA 8", "");
    legTemplateYield->AddEntry((TObject*) 0, "Monash 2013", "");
    legTemplateYield->AddEntry((TObject*) 0, "GEANT 3", "");

    CorrBkgFile       = SafelyOpenRootfile("CorrBkgFileNoRebin.root");
    if (CorrBkgFile->IsOpen() ) printf("CorrBkgFile opened successfully\n");

    BackFileNN        = SafelyOpenRootfile("BackFileNN.root");
    if (BackFileNN->IsOpen() ) printf("BackFileNN opened successfully\n");

    BackFile3to8      = SafelyOpenRootfile("BackFile3to8.root");
    if (BackFile3to8->IsOpen() ) printf("BackFile3to8 opened successfully\n");

    /*
    Open the correct file corresponding to the current templatemethod
     */
    switch (templatemethod) {
      case 1:
        OutputFile      = SafelyOpenRootfile("OutputFileBetterBkgNN.root");
        std::cout << "_____________________________________________________________________" << '\n';
        std::cout << "| Start Plotting Better Bkg NN Method" << '\n';
        break;
      case 2:
        OutputFile      = SafelyOpenRootfile("OutputFileNormal.root");
        std::cout << "_____________________________________________________________________" << '\n';
        std::cout << "| Start Plotting 'Normal' Method" << '\n';
        break;
      default:
        OutputFile      = SafelyOpenRootfile("OutputFileBetterBkg3to8.root");
        std::cout << "_____________________________________________________________________" << '\n';
        std::cout << "| Start Plotting Better Bkg 3 to 8 Method" << '\n';
        break;
    }
    if (OutputFile->IsOpen() ) printf("OutputFile opened successfully\n");

    /**
     * uncorrected Yield obtained via my (template) method
     */
    hYield_dt_chi2map_uncorr    = (TH1D*) OutputFile->Get(Form("hYield_dt_chi2map_uncorr"));
    SetHistogramProperties(hYield_dt_chi2map_uncorr, "pt", rawyield, 0, 1.4, 14.);

    /**
     * Chi^2/ndf (pT) from my method
     */
    hChi2Map_Chi2_pT            = (TH1D*)OutputFile->Get("hChi2Map_Chi2_pT");
    SetHistogramProperties(hChi2Map_Chi2_pT, "pt", "#chi^{2}/ndf", 0, 1.4, 14.);

    /**
     * Chi^2/ndf from the framework method with function parametrization
     */
    histoChi2_0                 = (TH1D*) OutputFile->Get("histoChi2_0");
    SetHistogramProperties(histoChi2_0, "pt", "#chi^{2}/ndf", 5, 1.4, 14.);

    /**
     * Signal Area Scaling factor. Currently disabled so NO USE!
     */
    hSignalAreaScaling          = (TH1D*)OutputFile->Get("hSignalAreaScaling");
    SetHistogramProperties(hSignalAreaScaling, "pt", "FS_{Signal}", 0, 1.4, 14.);

    /**
     * Correlated Background Area Scaling factor. Currently disabled so NO USE!
     */
    hCorrbackAreaScaling        = (TH1D*)OutputFile->Get("hCorrbackAreaScaling");
    SetHistogramProperties(hCorrbackAreaScaling, "pt", "FS_{korr. Untergrund}", 0, 1.4, 14.);

    /**
     * corrected Yield obtained via my (template) method
     */
    hYield_dt_chi2map_corrected = (TH1D*)OutputFile->Get("hYield_dt_chi2map_corrected");
    SetHistogramProperties(hYield_dt_chi2map_corrected, "pt", strCorrectedYield, 0, 1.4, 14.);

    /**
     * corrected Yield with the framework method (function parametrization)
     */
    hCorrectedYieldNormEff      = (TH1D*) OutputFile->Get("hCorrectedYieldNormEff");
    SetHistogramProperties(hCorrectedYieldNormEff, "pt", strCorrectedYield, 5, 1.4, 14.);

    /**
     * Signal Template Scaling factor
     */
    h_x_min = (TH1D*)OutputFile->Get("h_x_min");
    SetHistogramProperties(h_x_min, "pt", "SF_{signal}", 2, 1.4, 14.);

    /**
     * corr. bkg. template sclaing factor
     */
    h_y_min                           = (TH1D*)OutputFile->Get("h_y_min");
    SetHistogramProperties(h_y_min, "pt", "SF_{corr. background}", 8, 1.4, 14.);


    /*
    My own efficiency; NOT from the after burner!
     */
    hEfficiency                       = (TH1D*) OutputFile->Get("hEfficiency");
    SetHistogramProperties(hEfficiency, "pt", "correction factors", 0, 1.4, 14.);

    /*
    The normal acceptance (from the afterburner)
     */
    hAcc                              = (TH1D*) OutputFile->Get("hAcc");
    SetHistogramProperties(hAcc, "pt", "correction factors", 1, 1.4, 14.);

    /*
    Ratio: my Yield / AfterBurner Yield
     */
    hCorrYieldME_Ratio                = (TH1D*) OutputFile->Get("hCorrYieldME_Ratio");
    SetHistogramProperties(hCorrYieldME_Ratio, "pt", "ratio", 4, 1.4, 14.);

    /*
    Statistical uncertainties of my yield
     */
    hCorrYieldME_StatError            = (TH1D*) OutputFile->Get("hCorrYieldME_StatError");
    SetHistogramProperties(hCorrYieldME_StatError, "pt", StatUn_Str, 0, 1.4, 14.);

    /*
    Statistical uncertainties of the yield from the after burner
     */
    hCorrectedYieldNormEff_StatError  = (TH1D*) OutputFile->Get("hCorrectedYieldNormEff_StatError");
    SetHistogramProperties(hCorrectedYieldNormEff_StatError, "pt", StatUn_Str, 5, 1.4, 14.);

    /*
    one TLine at 0.0 for the minv dependent pictures
     */
    line_0 = new TLine(0.0, 0.0, 0.3, 0.0);
    line_0->SetLineWidth(3);
    line_0->SetLineStyle(1);
    line_0->SetLineColor(kBlack);

    /*
    one TLine at 0.0 for the pT dependent pictures
     */
    line_one = new TLine(1.4, 1.0, 14.0, 1.0);
    line_one->SetLineWidth(2);
    line_one->SetLineStyle(3);
    line_one->SetLineColor(kBlack);




    //////////////////////////////////////////////////////////////////////////////
    // going over all pt bins despite first one, which is some framework bs.
    for (int k = 1; k < numberbins-1; k++) {


        /***********only possible if NOT lighweight train output!***********/
        // hGG                      = (TH1D*) MCWithOutFile->Get(Form("Mapping_TrueMesonCaloPhoton_InvMass_in_Pt_Bin%02d", k));
        // SetHistogramProperties(hGG, "minv", count_str, 2, 0.0, 0.3);
        //
        // hGC                      = (TH1D*) MCWithOutFile->Get(Form("Mapping_TrueMesonMixedCaloConvPhoton_InvMass_in_Pt_Bin%02d", k));
        // SetHistogramProperties(hGC, "minv", count_str, 4, 0.0, 0.3);
        //
        // hCC                      = (TH1D*) MCWithOutFile->Get(Form("Mapping_TrueMesonCaloConvPhoton_InvMass_in_Pt_Bin%02d", k));
        // SetHistogramProperties(hCC, "minv", count_str, 7, 0.0, 0.3);


        /**
         * Same - scaled mixed event from MC
         */
        hInvMass_MC              = (TH1D*) MCWithOutFile->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d", k));
        SetHistogramProperties(hInvMass_MC, "minv", count_str, 5, 0.0, 0.3);


        /**
         * Same - scaled mixed event from data
         */
        hData                    = (TH1D*) OutputFile->Get(Form("data_bin%02i",k));
        SetHistogramProperties(hData, "minv", count_str, 5, 0.0, 0.3);
        str                      = hData->GetTitle(); // pT range string
        hData->SetTitle("");

        /**
         * Chi2Map from my method
         */
        hChi2_2D                 = (TH2D*) OutputFile->Get(Form("hChi2MapBin%02d",k));

        /**
         * Chi2Map that only contains the 1 sigma range
         */
        hChi2_2D_sigma           = (TH2D*) OutputFile->Get(Form("hChi2_2D_sigma_bin%02d",k));


        /**
         * Signal Template UNSCALED!
         */
        hSignal                  = (TH1D*) OutputFile->Get(Form("hSignal_bin%02d",k));
        SetHistogramProperties(hSignal, "minv", count_str, 1, 0.0, 0.3);

        /**
         * corr. Bkg template UNSCALED!
         */
        hCorrBack                = (TH1D*) OutputFile->Get(Form("hCorrBack_bin%02d",k));
        SetHistogramProperties(hCorrBack, "minv", count_str, 3, 0.0, 0.3);


        /*
        scaling of the signal and correalted background templates
         */
        hSignal_scaled = (TH1D*) hSignal->Clone("hSignal_scaled");
        hCorrBack_scaled = (TH1D*) hCorrBack->Clone("hCorrBack_scaled");


        hSignal_scaled->       Scale(hSignalAreaScaling->  GetBinContent(k+1)*h_x_min->GetBinContent(k+1));
        hCorrBack_scaled->     Scale(hCorrbackAreaScaling->GetBinContent(k+1)*h_y_min->GetBinContent(k+1));

        /*
        Histogram containing the combination of scaled (signal + correlated
        background)
         */
        hAdded                = (TH1D*) hSignal_scaled->   Clone("hAdded");
        hAdded->              Add(hCorrBack_scaled);
        SetHistogramProperties(hAdded, "minv", count_str, 8, 0.0, 0.3);

        /*
        Setting up the pT Legend!
         */
        str = Form("%.1lf #leq #it{p}_{T} /(GeV/#it{c}) < %.1lf", fBinsPi013TeVEMCPt[k], fBinsPi013TeVEMCPt[k+1]);
        legpT = new TLegend(0.09, 0.85, 0.38, 0.9);
        legpT->AddEntry((TObject*) 0, str, "");

        /**
         * Getting the two histograms which show the Ratio of the corr. bkg Template
         * of the normal method and the NN/3to8 method
         */
        if(templatemethod == 0){

          /*
          get the correalted background without rebinning and rebin them with
          the proper value of the current pT bin.
           */

          hCorrBackNoRebin    = NULL;
          hCorrBackNoRebin    = (TH1D*) CorrBkgFile->Get(Form("hCorrBkgNoRebinBin%02d",k));
          hCorrBackNoRebin->Rebin(fBinsPi013TeVEMCPtRebin[k-1]);

          hRatio_Bkg = NULL;
          fPol0      = NULL;

          /*
          Get the ratio plot for the correlated background (combined/single)
           */
          hRatio_Bkg = (TH1D*) BackFileNN->Get(Form("hRatio_Bin%02d", k));
          SetHistogramProperties(hCorrBackNoRebin, "minv", count_str, 1, 0.0, 0.3);
          SetHistogramProperties(hRatio_Bkg, "minv", "ratio", 5, 0.0, 0.3);

          /*
          a constant fit to the ratio above
           */
          fPol0      = (TF1*)  BackFileNN->Get(Form("fPol0_Bin%02d", k));

          hRatio_Bkg->GetYaxis()->SetRangeUser(fPol0->GetParameter(0)-1.54, fPol0->GetParameter(0)+1.54);

          OAhists->Clear();
          OAratios->Clear();

          /*
          TLegends for the correlated background templates in comparison.
           */
          TLegend* legCorrBkgComp = new TLegend(0.15, 0.75, 0.5, 0.85);
          legCorrBkgComp->AddEntry(hCorrBack, "corr. background (combined)", "p");
          legCorrBkgComp->AddEntry(hCorrBackNoRebin, "corr. background (single)", "p");

          TLegend* legTemplate = new TLegend(0.6, 0.5, 0.8, 0.75);
          legTemplate->AddEntry((TObject*) 0, "templates:", "");
          legTemplate->AddEntry((TObject*) 0, "PYTHIA 8", "");
          legTemplate->AddEntry((TObject*) 0, "Monash 2013", "");
          legTemplate->AddEntry((TObject*) 0, "GEANT 3", "");

          TLegend* legpTBckRatio = new TLegend(0.6, 0.8, 0.8, 0.85);
          legpTBckRatio->AddEntry((TObject*) 0, str, "");

          /*
          plotting part using Plotting_Patrick.h
           */
          hCorrBack->GetYaxis()->SetRangeUser(-1.5e3, 10.5e3);
          OAhists->Add(hCorrBack);
          OAhists->Add(hCorrBackNoRebin);

          OAhists->Add(legCorrBkgComp);
          OAhists->Add(legSystem);
          OAhists->Add(legpTBckRatio);
          OAhists->Add(legTemplate);
          OAratios->Add(hRatio_Bkg);
          OAratios->Add(fPol0);

          c2 = makeCanvas(OAhists, OAratios, "notimeThickSquare", 0, 0);

          c2->Update();
          c2->SaveAs(Form("BetterBkg3to8/BackgroundWithRatio%02d" + SaveAppendix + "." + PicFormat,k));
          c2->Clear();

          delete legCorrBkgComp;
        }
        else if(templatemethod == 1){
          hCorrBackNoRebin    = NULL;
          hCorrBackNoRebin    = (TH1D*) CorrBkgFile->Get(Form("hCorrBkgNoRebinBin%02d", k));
          hCorrBackNoRebin->Rebin(fBinsPi013TeVEMCPtRebin[k-1]);

          /*
          TLegend which displays most important info as header!
           */
          TLegend* legSystemBWR = new TLegend(0.03, 0.85, 0.5, 0.9);
          legSystemBWR->AddEntry((TObject*) 0, "ALICE, pp #sqrt{#it{s}} = 13 TeV, #pi^{0} #rightarrow #gamma#gamma with EMCal", "");


          hRatio_Bkg = NULL;
          fPol0      = NULL;

          hRatio_Bkg = (TH1D*) BackFile3to8->Get(Form("hRatio_Bin%02d", k));
          fPol0      = (TF1*)  BackFile3to8->Get(Form("fPol0_Bin%02d", k));

          SetHistogramProperties(hCorrBackNoRebin, "minv", count_str, 3, 0.0, 0.3);
          SetHistogramProperties(hCorrBack, "minv", count_str, 8, 0.0, 0.3);
          SetHistogramProperties(hRatio_Bkg, "minv", "#frac{single}{combined}", 5, 0.0, 0.3);
          hRatio_Bkg->GetYaxis()->SetRangeUser(fPol0->GetParameter(0)-1.5, fPol0->GetParameter(0)+1.5);
          fPol0->SetLineStyle(1);

          OAhists->Clear();
          OAratios->Clear();

          TLegend* legCorrBkgComp = new TLegend(0.15, 0.75, 0.5, 0.85);
          legCorrBkgComp->AddEntry(hCorrBack, "corr. background (combined)", "p");
          legCorrBkgComp->AddEntry(hCorrBackNoRebin, "corr. background (single)", "p");

          TLegend* legTemplate = new TLegend(0.6, 0.5, 0.8, 0.75);
          legTemplate->AddEntry((TObject*) 0, "templates:", "");
          legTemplate->AddEntry((TObject*) 0, "PYTHIA 8", "");
          legTemplate->AddEntry((TObject*) 0, "Monash 2013", "");
          legTemplate->AddEntry((TObject*) 0, "GEANT 3", "");

          TLegend* legpTBckRatio = new TLegend(0.6, 0.8, 0.8, 0.85);
          legpTBckRatio->AddEntry((TObject*) 0, str, "");

          hCorrBack->GetYaxis()->SetRangeUser(-1.5e3, 10.5e3);
          OAhists->Add(hCorrBack);
          OAhists->Add(hCorrBackNoRebin);

          OAhists->Add(legCorrBkgComp);
          OAhists->Add(legTemplate);
          OAhists->Add(legpTBckRatio);

          OAratios->Add(hRatio_Bkg);
          OAratios->Add(fPol0);
          OAratios->Add(legSystemBWR);

          c2 = makeCanvas(OAhists, OAratios, "notimeThickSquare", 0, 0);


          c2->Update();
          c2->SaveAs(Form("BetterBkgNN/BackgroundWithRatio%02d" + SaveAppendix + "." + PicFormat,k));
          c2->Clear();

          delete legCorrBkgComp;
          delete legSystemBWR;

        }
        else{
          hRatio_Bkg = NULL;
          fPol0      = NULL;
        }


        if(wpsid == "all" || wpsid.Contains("chi2map")){
          ////////////////////////////////////////////////////////////////////////
          // Drawing Chi2 maps

          TGaxis::SetMaxDigits(1); //setting maxnumber of digits to 1 which may
                                   //need to be changed back to 3

          /*
          Setting the axis titles for the Chi2Map
           */
          hChi2_2D->SetXTitle("#lambda_{S}");
          hChi2_2D->SetYTitle("#lambda_{CB}");
          hChi2_2D->SetZTitle("#chi^{2}");

          /*
          Set Up and Drawing of the Chi2Map which is needed first to be able
          to draw the white Chi²_min + 1 line which is saved in a different
          histogram.
           */
          OAhists->Clear();
          OAhists->Add(hChi2_2D);
          OAhists->Add(legSystemChi2Map);
          legpT->SetTextColor(kWhite);
          OAhists->Add(legpT);
          legTemplat->SetTextColor(kWhite);
          OAhists->Add(legTemplat);


          cChi2Map = NULL;
          lChi2MinY = new TLine(h_x_min->GetBinContent(k+1),
                                h_y_min->GetBinContent(k+1)-0.001,
                                h_x_min->GetBinContent(k+1),
                                h_y_min->GetBinContent(k+1)+0.001);

          lChi2MinY->SetLineWidth(5);
          lChi2MinY->SetLineColor(kBlack);

          lChi2MinX = new TLine(h_x_min->GetBinContent(k+1)-0.001,
                                h_y_min->GetBinContent(k+1),
                                h_x_min->GetBinContent(k+1)+0.001,
                                h_y_min->GetBinContent(k+1));
          lChi2MinX->SetLineWidth(5);
          lChi2MinX->SetLineColor(kBlack);
          cChi2Map = makeCanvas(OAhists, 0, "colznotimesquare", 0, 0);
          cChi2Map->cd();
          cChi2Map->Update();
          /*
          After drawing the Chi2Map once we can now draw the 1 sigma range which
          is not an option in Plotting_Patrick.h
           */
          hChi2_2D_sigma->SetLineColor(kWhite);
          hChi2_2D_sigma->SetLineWidth(2);
          hChi2_2D_sigma->SetContour(2, somelist);
          hChi2_2D_sigma->Draw("SAME CONT3");
          lChi2MinX->Draw("SAME");
          lChi2MinY->Draw("SAME");

          TGaxis::SetMaxDigits(1);
          cChi2Map->Update();

          cChi2Map->Update();
          if(templatemethod == 1){
            cChi2Map->SaveAs(Form("BetterBkgNN/Chi2Map%02d" + SaveAppendix + "." + PicFormat,k));
          }
          if(templatemethod == 0){
            cChi2Map->SaveAs(Form("BetterBkg3to8/Chi2Map%02d" + SaveAppendix + "." + PicFormat,k));
          }
          if(templatemethod == 2){
            cChi2Map->SaveAs(Form("Normal/Chi2Map%02d" + SaveAppendix + "." + PicFormat,k));
          }
          cChi2Map->Clear();
          TGaxis::SetMaxDigits(3); // Resetting MaxDigits to 3
          delete lChi2MinX;
          delete lChi2MinY;
        }


        /**
         * Drawing of the parametrization result using templates
         * one time only the components with the data, one time the complete
         * parametrization comapred with the data.
         */
        if(wpsid == "all" || wpsid.Contains("param")){

          OAhists->Clear();
          OAhists->Add(hData);
          OAhists->Add(hSignal_scaled);
          OAhists->Add(hCorrBack_scaled);
          OAhists->Add(legSystem);
          legpT->SetTextColor(kBlack);
          OAhists->Add(legpT);
          legTemplat->SetTextColor(kBlack);
          OAhists->Add(legTemplat);

          c1 = makeCanvas(OAhists, 0, "notimethickhorizontal", 0, 0);
          c1->cd();
          c1->Update();

          TLegend* lParamResultParts  = new TLegend(0.59, 0.8, 0.85, 0.9);
          TLegend* lParamResultParts2 = new TLegend(0.53, 0.75, 0.82, 0.8);
          TLegend* lParamResultParts3 = new TLegend(0.59, 0.65, 0.85, 0.75);

          lParamResultParts->   AddEntry(hData, "signal", "p");
          lParamResultParts->   AddEntry((TObject*) 0, "+ corr. background:", "");
          lParamResultParts2->  AddEntry((TObject*) 0, "template:", "");
          lParamResultParts3->  AddEntry(hSignal_scaled,   "signal",      "p");
          lParamResultParts3->  AddEntry(hCorrBack_scaled, "corr. background",  "p");


          double line_y = gPad->GetUymax()*0.999;
          TLine* paramrange = new TLine(lowerparamrange/*[k-1]*/, line_y, upperparamrange, line_y);
          paramrange->SetLineColor(kAzure+10);
          paramrange->SetLineWidth(7);

          lParamResultParts3->   AddEntry(paramrange, "parametrisation range",  "l");

          OAhists->Clear();
          OAhists->Add(hData);
          OAhists->Add(hSignal_scaled);
          OAhists->Add(hCorrBack_scaled);
          OAhists->Add(lParamResultParts);
          OAhists->Add(lParamResultParts2);
          OAhists->Add(lParamResultParts3);
          OAhists->Add(paramrange);
          OAhists->Add(line_0);
          OAhists->Add(legSystem);
          OAhists->Add(legpT);
          OAhists->Add(legTemplat);

          c1 = makeCanvas(OAhists, 0, "notimethickhorizontal", 0, 0);

          // DrawLabelALICE(0.6, 0.9, 0.02, 40, str);

          c1->Update();

          c1->Update();
          if(templatemethod == 1){
            c1->SaveAs(Form("BetterBkgNN/ParamResultParts_Bin%02d" + SaveAppendix + "." + PicFormat,k));
          }
          if(templatemethod == 0){
            c1->SaveAs(Form("BetterBkg3to8/ParamResultParts_Bin%02d" + SaveAppendix + "." + PicFormat,k));
          }
          if(templatemethod == 2){
            c1->SaveAs(Form("Normal/ParamResultParts_Bin%02d" + SaveAppendix + "." + PicFormat,k));
          }
          c1->Clear();

          TLegend* lParamResult  = new TLegend(0.59, 0.85, 0.83, 0.9);
          TLegend* lParamResult4  = new TLegend(0.555, 0.8, 0.78, 0.85);
          TLegend* lParamResult2 = new TLegend(0.49, 0.75, 0.83, 0.8);
          TLegend* lParamResult3 = new TLegend(0.59, 0.6, 0.83, 0.75);

          lParamResult->   AddEntry(hData, "signal", "p");
          lParamResult4->   AddEntry((TObject*) 0, "+ corr. background", "");
          lParamResult2->  AddEntry((TObject*) 0, "parametrisation:", "");
          lParamResult3->  AddEntry(hAdded,   "T_{comb}",      "p");
          lParamResult3->  AddEntry(hCorrBack_scaled, "#lambda_{CB}T_{CB}",  "p");
          lParamResult3->  AddEntry(paramrange, "range",  "l");

          hCorrBack_scaled->SetMarkerStyle(33);
          hCorrBack_scaled->SetMarkerSize(1.7);

          OAhists->Clear();
          OAhists->Add(hData);
          OAhists->Add(hAdded);
          OAhists->Add(hCorrBack_scaled);
          OAhists->Add(paramrange);
          OAhists->Add(lParamResult);
          OAhists->Add(lParamResult4);
          OAhists->Add(lParamResult2);
          OAhists->Add(lParamResult3);
          OAhists->Add(paramrange);
          OAhists->Add(line_0);
          OAhists->Add(legSystem);
          OAhists->Add(legpT);
          OAhists->Add(legTemplat);

          c1 = makeCanvas(OAhists, 0, "notimethickSquare", 0, 0);

          c1->Update();

          if(templatemethod == 1){
            c1->SaveAs(Form("BetterBkgNN/ParamResult_Bin%02d" + SaveAppendix + "." + PicFormat,k));
          }
          if(templatemethod == 0){
            c1->SaveAs(Form("BetterBkg3to8/ParamResult_Bin%02d" + SaveAppendix + "." + PicFormat,k));
          }
          if(templatemethod == 2){
            c1->SaveAs(Form("Normal/ParamResult_Bin%02d" + SaveAppendix + "." + PicFormat,k));
          }

          c1->Clear();

          delete lParamResultParts;
          delete lParamResult;
          delete paramrange;

          /***********only possible if NOT lighweight train output!***********/
          // if(templatemethod == 0){
          //   TLegend* lGammas = new TLegend(0.63, 0.7, 0.8, 0.9);
          //   lGammas->AddEntry(hSignal, "MC true #pi^{0} signal", "p");
          //   lGammas->AddEntry(hGG, "#gamma#gamma", "p");
          //   lGammas->AddEntry(hGC, "#gamma#gamma_{conv}", "p");
          //   lGammas->AddEntry(hCC, "#gamma_{conv}#gamma_{conv}",  "p");
          //
          //
          //   OAhists->Clear();
          //   OAhists->Add(hSignal);
          //   OAhists->Add(hGG);
          //   OAhists->Add(hGC);
          //   OAhists->Add(hCC);
          //   OAhists->Add(lGammas);
          //   OAhists->Add(legSystem);
          //   OAhists->Add(line_0);
          //   OAhists->Add(legpT);
          //   OAhists->Add(legTemplat);
          //
          //   c1 = makeCanvas(OAhists, 0, "notimethicksquare", 0, 0);
          //
          //   c1->Update();
          //   // DrawLabelALICE(0.6, 0.9, 0.02, 40, str);
          //   c1->Update();
          //
          //   c1->SaveAs(Form("BetterBkg3to8/PeakTemplateMotivation%02d" + SaveAppendix + "." + PicFormat,k));
          //
          // }

          if(templatemethod == 2){
            c1->Clear();

            TLegend* lTemplates  = new TLegend(0.56, 0.85, 0.8, 0.9);
            TLegend* lTemplates2 = new TLegend(0.46, 0.8, 0.75, 0.85);
            TLegend* lTemplates3 = new TLegend(0.56, 0.7, 0.8, 0.8);

            lTemplates->AddEntry(hInvMass_MC, "MC_{S+CB}", "p");
            lTemplates2->AddEntry((TObject*) 0, "template:" , "");
            lTemplates3->AddEntry(hSignal, "T_{S}", "p");
            lTemplates3->AddEntry(hCorrBack, "T_{CB}", "p");

            OAhists->Clear();
            OAhists->Add(hInvMass_MC);
            OAhists->Add(hSignal);
            OAhists->Add(hCorrBack);
            OAhists->Add(lTemplates);
            OAhists->Add(lTemplates2);
            OAhists->Add(lTemplates3);
            OAhists->Add(legSystem);
            OAhists->Add(line_0);
            OAhists->Add(legpT);
            OAhists->Add(legTemplat);

            c1 = makeCanvas(OAhists, 0, "notimetSquare", 0, 0);

            c1->Update();
            c1->SaveAs(Form("Normal/EntstehungUntergrund%02d" + SaveAppendix + "." + PicFormat,k));
          }
          c1->Clear();


        }
      delete legpT;
    }

    /**
     * Drawing of Chi^2/ndf (pT) used as comparison between my method and the
     * framework function parametrization
     */
    if(wpsid == "all" || wpsid.Contains("chi2")){
      TLegend* leg2 = new TLegend(0.3,0.75,0.6,0.9);
      SetLegendSettigns(leg2, 40);
      leg2->SetHeader("parametrization with:");
      leg2->AddEntry(hChi2Map_Chi2_pT, "template", "l");
      leg2->AddEntry(histoChi2_0, "function", "l");

      // hChi2Map_Chi2_pT->GetYaxis()->SetRangeUser(0.0, 4.0);

      OAhists->Clear();
      if(hChi2Map_Chi2_pT->GetMaximum()>histoChi2_0->GetMaximum()){
        OAhists->Add(hChi2Map_Chi2_pT);
        OAhists->Add(histoChi2_0);
      }

      else{
        OAhists->Add(histoChi2_0);
        OAhists->Add(hChi2Map_Chi2_pT);
      }
      OAhists->Add(leg2);
      OAhists->Add(legSystem);
      OAhists->Add(legTemplat);

      c1 = makeCanvas(OAhists, 0, "notimethickhorizontalLines", 0, 0);
      // DrawLabelALICE(0.6, 0.9, 0.03, 40);


      c1->Update();
      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkgNN/Chi2Comparison" + SaveAppendix + "." + PicFormat));
      }
      if(templatemethod == 0){
        c1->SaveAs(Form("BetterBkg3to8/Chi2Comparison" + SaveAppendix + "." + PicFormat));
      }
      if(templatemethod == 2){
        c1->SaveAs(Form("Normal/Chi2Comparison" + SaveAppendix + "." + PicFormat));
      }
      c1->Clear();
      c1->SetLogx(0);

      OAhists->Clear();
      OAhists->Add(hChi2Map_Chi2_pT);
      OAhists->Add(legSystem);
      OAhists->Add(legTemplatChi2);

      c1 = makeCanvas(OAhists, 0, "notimethickhorizontalLines", 0, 0);

      c1->Update();
      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkgNN/Chi2NoComp" + SaveAppendix + "." + PicFormat));
      }
      if(templatemethod == 0){
        c1->SaveAs(Form("BetterBkg3to8/Chi2NoComp" + SaveAppendix + "." + PicFormat));
      }
      if(templatemethod == 2){
        c1->SaveAs(Form("Normal/Chi2NoComp" + SaveAppendix + "." + PicFormat));
      }
      c1->Clear();

      delete leg2;
    }

    /**
     * Drawing of the uncorrected Yield from my method
     */
    if(wpsid == "all" || wpsid.Contains("uncorryield")){

      TLegend* leg = new TLegend(0.2,0.2,0.4,0.4);
      hYield_dt_chi2map_uncorr->GetXaxis()->SetRangeUser(1.4, 14.);

      OAhists->Clear();
      OAhists->Add(hYield_dt_chi2map_uncorr);
      OAhists->Add(legSystem);
      OAhists->Add(legTemplateYield);

      c1 = makeCanvas(OAhists, 0, "notimeThickHorizontalogy", 0, 0);

      // DrawLabelALICE(0.55, 0.9, 0.018, 40);
      c1->Update();
      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkgNN/UncorrYields" + SaveAppendix + "." + PicFormat));
      }
      if(templatemethod == 0){
        c1->SaveAs(Form("BetterBkg3to8/UncorrYields" + SaveAppendix + "." + PicFormat));
      }
      if(templatemethod == 2){
        c1->SaveAs(Form("Normal/UncorrYields" + SaveAppendix + "." + PicFormat));
      }
      delete leg;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Factor comp between Chi2 map and Iterative method
    /**
     * Drawing of the scaling factors
     * as well as the ratio of corr. bkg/ signal scaling
     * @param wpsid [description]
     */
    if(wpsid == "all" || wpsid.Contains("factorcomp")){

      TFile* PulseFile = NULL;
      TH1D* hConvInter = NULL;
      TF1* fPulse      = NULL;


      OAhists->Clear();
      OAhists->Add(h_y_min);
      OAhists->Add(legSystem);
      OAhists->Add(legTemplat);

      c1 = makeCanvas(OAhists, 0, "notimethickhorizontal", 0, 0);

      c1->Update();
      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkgNN/BGFactorComp" + SaveAppendix + "." + PicFormat));
      }
      if(templatemethod == 0){
        c1->SaveAs(Form("BetterBkg3to8/BGFactorComp" + SaveAppendix + "." + PicFormat));
      }
      if(templatemethod == 2){
        c1->SaveAs(Form("Normal/BGFactorComp" + SaveAppendix + "." + PicFormat));
      }
      c1->Clear();

      TH1D* SFratio = (TH1D*) h_y_min->Clone("SFratio");
      SFratio->Divide(h_x_min);
      SFratio->SetYTitle("ratio");
      SetHistogramProperties(h_x_min, "", "#lambda_{S}", 2, 1.4, 14.);
      SetHistogramProperties(h_y_min, "", "#lambda_{CB}", 4, 1.4, 14.);

      TLegend* SFLegend = new TLegend(0.6,0.75,0.85,0.9);
      SFLegend->SetHeader("Skalierungsfaktor:");
      SFLegend->AddEntry(h_y_min, "korr. Untergrund", "p");
      SFLegend->AddEntry(h_x_min, "Signal", "p");
      SFratio->SetMarkerColor(kGreen-3);
      SFratio->SetLineColor(kGreen-3);

      OAhists->Clear();
      OAhists->Add(h_y_min);
      OAhists->Add(h_x_min);
      OAhists->Add(SFLegend);
      OAhists->Add(legSystem);
      OAhists->Add(legTemplat);
      OAratios->Clear();
      OAratios->Add(SFratio);

      c1 = makeCanvas(OAhists, OAratios, "notimethick", 0, 0);

      c1->Update();
      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkgNN/b_to_a_ratio" + SaveAppendix + "." + PicFormat));
      }
      if(templatemethod == 0){
        c1->SaveAs(Form("BetterBkg3to8/b_to_a_ratio" + SaveAppendix + "." + PicFormat));
      }
      if(templatemethod == 2){
        c1->SaveAs(Form("Normal/b_to_a_ratio" + SaveAppendix + "." + PicFormat));
      }
      c1->Clear();

      h_y_min->SetYTitle("scaling factors");
      h_y_min->SetXTitle("#it{p}_{T} (GeV/#it{c}^{2})");

      OAhists->Clear();
      OAhists->Add(h_y_min);
      OAhists->Add(h_x_min);
      OAhists->Add(SFLegend);
      OAhists->Add(legSystem);
      OAhists->Add(legTemplatChi2);
      OAratios->Clear();

      c1 = makeCanvas(OAhists, 0, "notimeThickHorizontal", 0, 0);

      c1->Update();
      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkgNN/SF" + SaveAppendix + "." + PicFormat));
      }
      if(templatemethod == 0){
        c1->SaveAs(Form("BetterBkg3to8/SF" + SaveAppendix + "." + PicFormat));
      }
      if(templatemethod == 2){
        c1->SaveAs(Form("Normal/SF" + SaveAppendix + "." + PicFormat));
      }
      c1->Clear();

      SetHistogramProperties(h_x_min, "pt", "SF_{signal}", 2, 1.4, 14.);
      SetHistogramProperties(h_y_min, "pt", "SF_{corr. background}", 4, 1.4, 14.);

      TLegend* lCorrection = new TLegend(0.15,0.7,0.6,0.9);
      lCorrection->AddEntry(hEfficiency, "efficiency" , "l");
      lCorrection->AddEntry(hAcc,  "acceptance", "l");

      OAhists->Clear();
      OAhists->Add(hEfficiency);
      OAhists->Add(hAcc);
      OAhists->Add(lCorrection);
      OAhists->Add(legSystem);

      c1 = makeCanvas(OAhists, 0, "notimethickhorizontal", 0, 0);
      c1->Update();

      if(templatemethod == 0){
        c1->SaveAs(Form("BetterBkg3to8/Korrekturfaktoren" + SaveAppendix + "." + PicFormat));
      }
      c1->Clear();


    }

    delete hData;
    delete hPol1PeakFactor;
    delete line_0;
    delete line_one;
    delete hChi2_2D;
    delete hChi2_2D_sigma;
    delete hChi2Map_Chi2_pT;
    delete hSignalAreaScaling;
    delete hCorrbackAreaScaling;
    delete h_x_min;
    delete h_y_min;
    delete hErrXlow;
    delete hErrXhigh;
    delete hErrYlow;
    delete hErrYhigh;
    delete hSignal;
    delete hCorrBack;
    delete f_ChiOverNdf;
    delete hYield_dt_chi2map_uncorr;
    delete hYield_dt_chi2map_corrected;
    delete hYield_pol1_corrected;
    delete hCorrectedYieldNormEff;
    delete histoChi2_0;

    OutputFile->Close();
  }
}
