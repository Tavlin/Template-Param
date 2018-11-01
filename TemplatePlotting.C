#include "CommonHeader.h"
//MCTemplatesAnData

// wpsid = which picture should I draw
void TemplatePlotting(TString wpsid = "all", TString PicFormat = "png"){
  int binnumber = -1;
 /**
  * loop over templatemethod
  * templatemethod == 0 for Next Neighbor method
  * templatemethod == 1 for 3 to 8 Method
  * templatemethod == 2 for BetterBkg3to8Pulse method
  */
  for (int templatemethod = 1; templatemethod <= 4; templatemethod++) {

    TString str;
    const Int_t nbins = numberbins;
    const Int_t ndrawpoints = 1.e5;
    Double_t somelist[2] = {1., 2.};
    TString TempStr = "template";

    //////////////////////////////////////////////////////////////////////////////
    // setting up the canvas to draw on. Will later be changed for the chi2 pic
    TCanvas *c1 = new TCanvas("c1","",1540,1417);
    c1->cd();
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.09);
    c1->SetRightMargin(0.02);
    c1->SetLeftMargin(0.09);
    c1->SetTicky();
    c1->SetTickx();
    c1->SetLogz(1);
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);

    TCanvas *c2 = new TCanvas("c2","",1540,1417);
    c2->cd();
    c2->SetTopMargin(0.05);
    c2->SetBottomMargin(0.09);
    c2->SetRightMargin(0.15);
    c2->SetLeftMargin(0.09);
    c2->SetTicky();
    c2->SetTickx();
    c2->SetLogz(1);

    TCanvas *canInvMass = new TCanvas("canInvMass","",1200,1600);
    TPad *pad1InvMass = new TPad("pad1InvMass","",0.0,0.33,1.0,1.0);
    pad1InvMass->SetTopMargin(0.05);
    pad1InvMass->SetLeftMargin(0.15);
    pad1InvMass->SetBottomMargin(0.0);
    pad1InvMass->SetRightMargin(0.02);
    TPad *pad2InvMass = new TPad("pad2InvMass","",0.0,0.0,1.0,0.33);
    pad2InvMass->SetTopMargin(0.0);
    pad2InvMass->SetLeftMargin(0.15);
    pad2InvMass->SetBottomMargin(0.3);
    pad2InvMass->SetRightMargin(0.02);
    pad2InvMass->SetTicky();

    //////////////////////////////////////////////////////////////////////////////
    // setting up the 2 Histograms to compare chi2 from the to fit methods as
    // well as peak factor comp. between pol 1 and double temp fit and the ratio
    // of BG. scaling factor and the Peak scaling factor

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
    TH1D* hCorrBack                   = NULL;
    TH1D* hGG                         = NULL;
    TH1D* hGC                         = NULL;
    TH1D* hCC                         = NULL;
    TH1D* hInvMass_MC                 = NULL;
    TH1D* hEfficiency                 = NULL;
    TH1D* hAcc                        = NULL;
    TF1* f_ChiOverNdf                 = NULL;
    TLine* fitrange2                  = NULL;
    TLine* line_0                     = NULL;
    TLine* line_p1                    = NULL;
    TLine* line_m1                    = NULL;
    TLine* line_p3                    = NULL;
    TLine* line_m3                    = NULL;
    TLine* line_one                   = NULL;

    Double_t line_y = 0;

    TFile* MCWithOutFile = NULL;
    MCWithOutFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_MC_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");
    TFile* OutputFile = NULL;
    if(templatemethod == 2){
      OutputFile      = SafelyOpenRootfile("OutputFileBetterBkgNNforAdrian.root");
    }
    else if(templatemethod == 1){
      OutputFile      = SafelyOpenRootfile("OutputFileBetterBkg3to8.root");
    }
    else if(templatemethod == 3){
      OutputFile      = SafelyOpenRootfile("OutputFileBetterBkgPulse.root");
    }
    else if(templatemethod == 4){
      OutputFile      = SafelyOpenRootfile("OutputFileNormal.root");
    }
    else{
      std::cerr << "No Outputfile loaded!" << std::endl;
    }
    if (OutputFile->IsOpen() ) printf("OutputFile opened successfully\n");

    /**
     * uncorrected Yield obtained via my (template) method
     */
    hYield_dt_chi2map_uncorr    = (TH1D*) OutputFile->Get(Form("hYield_dt_chi2map_uncorr"));

    /**
     * Chi^2/ndf (pT) from my method
     */
    hChi2Map_Chi2_pT            = (TH1D*)OutputFile->Get("hChi2Map_Chi2_pT");

    /**
     * Chi^2/ndf from the framework method with function parametrization
     */
    histoChi2_0                 = (TH1D*) OutputFile->Get("histoChi2_0");

    /**
     * Signal Area Scaling factor. Currently disabled so NO USE!
     */
    hSignalAreaScaling          = (TH1D*)OutputFile->Get("hSignalAreaScaling");

    /**
     * Correlated Background Area Scaling factor. Currently disabled so NO USE!
     */
    hCorrbackAreaScaling        = (TH1D*)OutputFile->Get("hCorrbackAreaScaling");

    /**
     * corrected Yield obtained via my (template) method
     */
    hYield_dt_chi2map_corrected = (TH1D*)OutputFile->Get("hYield_dt_chi2map_corrected");

    /**
     * corrected Yield with the framework method (function parametrization)
     */
    hCorrectedYieldNormEff      = (TH1D*) OutputFile->Get("hCorrectedYieldNormEff");

    /**
     * Signal Template Scaling factor
     */
    h_x_min = (TH1D*)OutputFile->Get("h_x_min");

    /**
     * corr. bkg. template sclaing factor
     * @param [name] [description]
     */
    h_y_min = (TH1D*)OutputFile->Get("h_y_min");
    hEfficiency = (TH1D*) OutputFile->Get("hEfficiency");
    hAcc = (TH1D*) OutputFile->Get("hAcc");



    line_0 = new TLine(0.0, 0.0, 0.4, 0.0);
    line_0->SetLineWidth(3);
    line_0->SetLineStyle(1);
    line_0->SetLineColor(kBlack);
    line_p1 = new TLine(0.0, 1.0, 0.4, 1.0);
    line_p1->SetLineWidth(2);
    line_p1->SetLineStyle(9);
    line_p1->SetLineColor(kGray+2);
    line_m1 = new TLine(0.0, -1.0, 0.4, -1.0);
    line_m1->SetLineWidth(2);
    line_m1->SetLineStyle(9);
    line_m1->SetLineColor(kGray+2);
    line_p3 = new TLine(0.0, 3.0, 0.4, 3.0);
    line_p3->SetLineWidth(2);
    line_p3->SetLineStyle(2);
    line_p3->SetLineColor(kGray+2);
    line_m3 = new TLine(0.0, -3.0, 0.4, -3.0);
    line_m3->SetLineWidth(2);
    line_m3->SetLineStyle(2);
    line_m3->SetLineColor(kGray+2);
    line_one = new TLine(1.4, 1.0, 12.0, 1.0);
    line_one->SetLineWidth(2);
    line_one->SetLineStyle(3);
    line_one->SetLineColor(kBlack);




    //////////////////////////////////////////////////////////////////////////////
    // going over all pt bins despite first one, which is some framework bs.
    for (int k = 1; k < numberbins-1; k++) {

      if(binnumber <=  0 || binnumber > numberbins){

        hGG                      = (TH1D*) MCWithOutFile->Get(Form("Mapping_TrueMesonCaloPhoton_InvMass_in_Pt_Bin%02d", k));
        hGC                      = (TH1D*) MCWithOutFile->Get(Form("Mapping_TrueMesonMixedCaloConvPhoton_InvMass_in_Pt_Bin%02d", k));
        hCC                      = (TH1D*) MCWithOutFile->Get(Form("Mapping_TrueMesonCaloConvPhoton_InvMass_in_Pt_Bin%02d", k));
        hInvMass_MC              = (TH1D*) MCWithOutFile->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d", k));

        /**
         * Same - scaled mixed event from data
         */
        hData                    = (TH1D*) OutputFile->Get(Form("data_bin%02i",k));
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

        /**
         * corr. Bkg template UNSCALED!
         */
        hCorrBack                = (TH1D*) OutputFile->Get(Form("hCorrBack_bin%02d",k));

        SetHistoStandardSettings(hData,      1.2, 1., 40, black);
        SetHistoStandardSettings(hInvMass_MC,1.2, 1., 40, black);
        SetHistoStandardSettings(hSignal,    1.2, 1., 40, teal-7);
        SetHistoStandardSettings(hCorrBack,  1.2, 1., 40, pink-2);

        hData->SetYTitle("d#it{N}/d#it{m}_{inv} (GeV/#it{c}^2)^{-1}");
        hSignal->SetYTitle("d#it{N}/d#it{m}_{inv} (GeV/#it{c}^2)^{-1}");
        hCorrBack->SetYTitle("d#it{N}/d#it{m}_{inv} (GeV/#it{c}^2)^{-1}");
        hData->SetXTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
        hSignal->SetXTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
        hCorrBack->SetXTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
        hInvMass_MC->SetYTitle("d#it{N}/d#it{m}_{inv} (GeV/#it{c}^2)^{-1}");
        hInvMass_MC->SetXTitle("#it{m}_{inv} (GeV/#it{c}^{2})");


        if(wpsid == "all" || wpsid.Contains("chi2map")){
          ////////////////////////////////////////////////////////////////////////
          // Drawing Chi2 maps
          c2->cd();

          TLegend* lChi2Map = new TLegend(0.4,0.85,0.7,0.9);
          SetLegendSettigns(lChi2Map, 40);
          lChi2Map->AddEntry((TObject*) 0x0, str, ""),

          SetHistoStandardSettings2(hChi2_2D);
          hChi2_2D->GetZaxis()->SetTitleOffset(0.5);
          hChi2_2D->SetXTitle("SF_{Signal}");
          hChi2_2D->SetYTitle("SF_{korr. Untergrund}");
          hChi2_2D->Draw("COLZ");
          hChi2_2D_sigma->SetLineColor(kWhite);
          hChi2_2D_sigma->SetLineWidth(2);
          hChi2_2D_sigma->SetContour(2, somelist);
          hChi2_2D_sigma->Draw("SAME CONT3");
          lChi2Map->Draw("SAME");
          c2->Update();

          c2->Update();
          if(templatemethod == 2){
            c2->SaveAs(Form("BetterBkgNN/Chi2Map%02d." + PicFormat,k));
          }
          if(templatemethod == 1){
            c2->SaveAs(Form("BetterBkg3to8/Chi2Map%02d." + PicFormat,k));
          }
          if(templatemethod == 3){
            c2->SaveAs(Form("BetterBkg3to8Pulse/Chi2Map%02d." + PicFormat,k));
          }
          if(templatemethod == 4){
            c2->SaveAs(Form("Normal/Chi2Map%02d." + PicFormat,k));
          }
          c2->Clear();
        }

        /**
         * Drawing of the parametrization result using templates
         * one time only the components with the data, one time the complete
         * parametrization comapred with the data.
         */
        if(wpsid == "all" || wpsid.Contains("param")){
          c1->cd();
          TLegend* lParamResultParts = new TLegend(0.2,0.5,0.4,0.63);
          SetLegendSettigns(lParamResultParts, 40);
          lParamResultParts->   AddEntry(hData, "Daten", "l");
          TH1D* hAdded          = NULL;
          TH1D* hSignal_Clone   = NULL;
          TH1D* hCorrBack_Clone = NULL;
          hSignal_Clone         = (TH1D*) hSignal->         Clone("hSignal_Clone");
          hCorrBack_Clone       = (TH1D*) hCorrBack->       Clone("hCorrBack_Clone");
          hSignal_Clone->       Scale(hSignalAreaScaling->  GetBinContent(k)*h_x_min->GetBinContent(k+1));
          // for(int i = 1; i < hCorrBack_Clone->fNcells -2; i++){
          //   hCorrBack_Clone->SetBinContent(i, hCorrBack_Clone->GetBinContent(i)*hSignalAreaScaling->GetBinContent(k)*h_y_min->GetBinContent(k+1));
          //   hCorrBack_Clone->SetBinError(i, sqrt(pow(hCorrBack_Clone->GetBinContent(i)*hSignalAreaScaling->GetBinContent(k)*h_y_min->GetBinError(k+1),2.)+
          //   pow(hCorrBack_Clone->GetBinError(i)*hSignalAreaScaling->GetBinContent(k)*h_y_min->GetBinContent(k+1),2));
          // }
          hCorrBack_Clone->     Scale(hCorrbackAreaScaling->GetBinContent(k)*h_y_min->GetBinContent(k+1));
          hAdded                = (TH1D*) hSignal_Clone->   Clone("hAdded");
          hAdded->              Add(hCorrBack_Clone);
          lParamResultParts->   AddEntry((TObject*) 0x0, "Template:", "");
          lParamResultParts->   AddEntry(hSignal_Clone,   "Signal",      "l");
          lParamResultParts->   AddEntry(hCorrBack_Clone, "korr. Untergrund.",  "l");

          hData->             Draw("AXIS");
          c1->Update();

          double line_y = gPad->GetUymax()*0.995;
          TLine* paramrange = new TLine(lowerparamrange[k-1], line_y, upperparamrange, line_y);
          paramrange->SetLineColor(kAzure+10);
          paramrange->SetLineWidth(7);
          paramrange->Draw("SAME");

          lParamResultParts->   AddEntry(paramrange, "Parametrisierungsbereich",  "l");

          hData->             Draw("SAME");
          hSignal_Clone->     Draw("SAME EP");
          hSignal_Clone->     Draw("SAME HIST");
          hCorrBack_Clone->   Draw("SAME EP");
          hCorrBack_Clone->   Draw("SAME HIST");
          paramrange->        Draw("SAME");
          lParamResultParts-> Draw("SAME");
          DrawLabelALICE(0.6, 0.9, 0.02, 40, str);

          c1->Update();

          c1->Update();
          if(templatemethod == 2){
            c1->SaveAs(Form("BetterBkgNN/ParamResultParts_Bin%02d." + PicFormat,k));
          }
          if(templatemethod == 1){
            c1->SaveAs(Form("BetterBkg3to8/ParamResultParts_Bin%02d." + PicFormat,k));
          }
          if(templatemethod == 3){
            c1->SaveAs(Form("BetterBkg3to8Pulse/ParamResultParts_Bin%02d." + PicFormat,k));
          }
          if(templatemethod == 4){
            c1->SaveAs(Form("Normal/ParamResultParts_Bin%02d." + PicFormat,k));
          }
          c1->Clear();
          TLegend* lParamResult = new TLegend(0.6,0.5,0.9,0.63);
          SetLegendSettigns(lParamResult, 40);
          lParamResult->AddEntry(hData, "Daten", "l");
          lParamResult->AddEntry((TObject*) 0x0, "Parametrisierung:", "");
          lParamResult->AddEntry(hAdded, "Template", "l");
          lParamResult->AddEntry(paramrange, "Bereich",  "l");


          c1->Update();
          hData->             Draw("AXIS");
          hData->             Draw("SAME");
          hAdded->            Draw("SAME");
          hAdded->            Draw("SAME HIST");
          lParamResult->      Draw("SAME");
          paramrange->        Draw("SAME");
          DrawLabelALICE(0.6, 0.9, 0.02, 40, str);
          c1->Update();

          if(templatemethod == 2){
            c1->SaveAs(Form("BetterBkgNN/ParamResult_Bin%02d." + PicFormat,k));
          }
          if(templatemethod == 1){
            c1->SaveAs(Form("BetterBkg3to8/ParamResult_Bin%02d." + PicFormat,k));
          }
          if(templatemethod == 3){
            c1->SaveAs(Form("BetterBkg3to8Pulse/ParamResult_Bin%02d." + PicFormat,k));
          }
          if(templatemethod == 4){
            c1->SaveAs(Form("Normal/ParamResult_Bin%02d." + PicFormat,k));
          }

          c1->Clear();

          delete lParamResultParts;
          delete lParamResult;
          delete paramrange;

          hSignal->SetMarkerSize(2);
          hSignal->SetYTitle("d#it{N}/d#it{m}_{inv} (GeV/#it{c}^2)^{-1}");
          hSignal->SetXTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
          hSignal->SetLineColor(kBlack);
          hGG->SetMarkerSize(2);
          hGC->SetMarkerSize(2);
          hCC->SetMarkerSize(2);
          hGC->SetMarkerStyle(34);
          hCC->SetMarkerStyle(33);

          TLegend* lGammas = new TLegend(0.2,0.7,0.4,0.9);
          SetLegendSettigns(lGammas, 40);
          lGammas->AddEntry(hSignal, "Summe", "p");
          lGammas->AddEntry(hGG, "#gamma#gamma", "p");
          lGammas->AddEntry(hGC, "#gamma#gamma_{conv}", "p");
          lGammas->AddEntry(hCC, "#gamma_{conv}#gamma_{conv}",  "p");

          c1->Update();
          hSignal->Draw("AXIS");
          hSignal->Draw("SAME");
          hGG->Draw("SAME");
          hGC->Draw("SAME");
          hCC->Draw("SAME");
          lGammas->Draw("SAME");
          DrawLabelALICE(0.6, 0.9, 0.02, 40, str);
          c1->Update();

          if(templatemethod == 4){
            c1->SaveAs(Form("BetterBkg3to8/PeakTemplateMotivation%02d." + PicFormat,k));
          }
          c1->Clear();

          TLegend* lTemplates = new TLegend(0.2,0.6,0.4,0.8);
          SetLegendSettigns(lTemplates, 40);
          lTemplates->AddEntry(hInvMass_MC, "Signal", "l");
          lTemplates->AddEntry((TObject*) 0x0, "+ korr. Untergrund", "");
          lTemplates->AddEntry((TObject*) 0x0, "+ unkorr. Untergrund", "");
          lTemplates->AddEntry(hSignal, "Signal Template", "l");
          lTemplates->AddEntry(hCorrBack, "korr. Untergrund Template", "l");

          c1->Update();
          hInvMass_MC-> Draw("AXIS");
          hInvMass_MC-> Draw("SAME");
          hSignal->     Draw("SAME");
          hCorrBack->   Draw("SAME");
          lTemplates->  Draw("SAME");
          DrawLabelALICE(0.6, 0.9, 0.02, 40, str);
          c1->Update();

          if(templatemethod == 4){
            c1->SaveAs(Form("BetterBkg3to8/EntstehungUntergrund%02d." + PicFormat,k));
          }
          c1->Clear();


        }
      }
    }

    /**
     * Drawing of Chi^2/ndf (pT) used as comparison between my method and the
     * framework function parametrization
     */
    if(wpsid == "all" || wpsid.Contains("chi2")){
      TLegend* leg2 = new TLegend(0.6,0.75,0.9,0.9);
      SetLegendSettigns(leg2, 40);
      leg2->SetHeader("Parametrisierung mit:");
      leg2->AddEntry(hChi2Map_Chi2_pT, "Template", "l");
      leg2->AddEntry(histoChi2_0, "Funktionen", "l");


      SetHistoStandardSettings(hChi2Map_Chi2_pT, 1.2, 1., 40, magenta-2);
      SetHistoStandardSettings(histoChi2_0     , 1.2, 1., 40, black);
      hChi2Map_Chi2_pT->GetXaxis()->SetRangeUser(1.4, 12.);
      histoChi2_0->     GetXaxis()->SetRangeUser(1.4, 12.);

      c1->cd();
      c1->Clear();
      c1->SetLogx(1);
      histoChi2_0->     Draw("AXIS");
      line_one->        Draw("same");
      histoChi2_0->     Draw("SAME HIST");
      hChi2Map_Chi2_pT->Draw("SAME HIST");
      hChi2Map_Chi2_pT->Draw("SAME P");
      leg2->            Draw("same");
      DrawLabelALICE(0.2, 0.9, 0.018, 40);


      c1->Update();
      if(templatemethod == 2){
        c1->SaveAs(Form("BetterBkgNN/Chi2Comparison." + PicFormat));
      }
      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkg3to8/Chi2Comparison." + PicFormat));
      }
      if(templatemethod == 3){
        c1->SaveAs(Form("BetterBkg3to8Pulse/Chi2Comparison." + PicFormat));
      }
      if(templatemethod == 4){
        c1->SaveAs(Form("Normal/Chi2Comparison." + PicFormat));
      }
      c1->Clear();
      c1->SetLogx(0);

      delete leg2;
    }

    /**
     * Drawing of the uncorrected Yield from my method
     */
    if(wpsid == "all" || wpsid.Contains("uncorryield")){

      TLegend* leg = new TLegend(0.2,0.2,0.4,0.4);
      SetLegendSettigns(leg);
      leg->AddEntry(hYield_dt_chi2map_uncorr, TempStr + " with chi2map" , "lp");
      hYield_dt_chi2map_uncorr->GetXaxis()->SetRangeUser(1.4, 12.);

      c1->cd();
      c1->SetLeftMargin(0.11);
      c1->Clear();
      c1->SetLogy(1);
      hYield_dt_chi2map_uncorr->Draw("AXIS");
      hYield_dt_chi2map_uncorr->Draw("samelp");
      leg->Draw("same");

      DrawLabelALICE(0.55, 0.9, 0.018, 40);
      c1->Update();
      if(templatemethod == 2){
        c1->SaveAs(Form("BetterBkgNN/UncorrYields." + PicFormat));
      }
      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkg3to8/UncorrYields." + PicFormat));
      }
      if(templatemethod == 3){
        c1->SaveAs(Form("BetterBkg3to8Pulse/UncorrYields." + PicFormat));
      }
      if(templatemethod == 4){
        c1->SaveAs(Form("Normal/UncorrYields." + PicFormat));
      }
      c1->Clear();
      c1->SetLogy(0);
      c1->SetLeftMargin(0.09);
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


      c1->cd();

      h_y_min->SetLineColor(kMagenta+2);
      h_y_min->SetMarkerColor(kMagenta+2);
      h_y_min->SetMarkerStyle(25);
      h_y_min->SetMarkerSize(1.5);
      h_y_min->GetXaxis()->SetRangeUser(-1.4, 12.);

      h_y_min->Draw();

      c1->Update();
      if(templatemethod == 2){
        c1->SaveAs(Form("BetterBkgNN/BGFactorComp." + PicFormat));
      }
      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkg3to8/BGFactorComp." + PicFormat));
      }
      if(templatemethod == 3){
        c1->SaveAs(Form("BetterBkg3to8Pulse/BGFactorComp." + PicFormat));
      }
      if(templatemethod == 4){
        c1->SaveAs(Form("Normal/BGFactorComp." + PicFormat));
      }
      c1->Clear();

      if(templatemethod == 3){
        PulseFile = SafelyOpenRootfile("Pulse.root");
        hConvInter = (TH1D*) PulseFile->Get("hConvInter");
        SetHistoStandardSettings(hConvInter);
        hConvInter->SetFillColor(2);
        fPulse      = (TF1*)  PulseFile->Get("fPulse");

        TLegend* leg = new TLegend(0.15,0.7,0.6,0.9);
        SetLegendSettigns(leg, 40);
        leg->AddEntry(h_y_min, "korr. Untergrund. Skalierungsfaktor" , "lp");
        leg->AddEntry(fPulse,  "Pulsefunktion", "l");
        leg->AddEntry(hConvInter, "Konfidenzintervall", "e");

        c1->Update();
        h_y_min->Draw("AXIS");
        fPulse->Draw("SAME");
        hConvInter->Draw("SAME E3");
        c1->Update();
        c1->SaveAs(Form("BetterBkg3to8Pulse/BkgConfidenceIntervall." + PicFormat));

        delete leg;

      }

      h_y_min->Divide(h_x_min);

      c1->cd();

      h_y_min->SetLineColor(kMagenta+2);
      h_y_min->SetMarkerColor(kMagenta+2);
      h_y_min->SetMarkerStyle(25);
      h_y_min->SetMarkerSize(1.5);
      h_y_min->GetYaxis()->SetRangeUser(-0.2, 5.);



      h_y_min->Draw();
      h_y_min->SetYTitle("back. scaling/signal scaling");

      c1->Update();
      if(templatemethod == 2){
        c1->SaveAs(Form("BetterBkgNN/b_to_a_ratio." + PicFormat));
      }
      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkg3to8/b_to_a_ratio." + PicFormat));
      }
      if(templatemethod == 3){
        c1->SaveAs(Form("BetterBkg3to8Pulse/b_to_a_ratio." + PicFormat));
      }
      if(templatemethod == 4){
        c1->SaveAs(Form("Normal/b_to_a_ratio." + PicFormat));
      }
      c1->Clear();

      h_y_min->Add(h_x_min, -1);

      TLegend* lCorrection = new TLegend(0.15,0.7,0.6,0.9);
      SetLegendSettigns(lCorrection, 40);
      lCorrection->AddEntry(hEfficiency, "Effizienz" , "l");
      lCorrection->AddEntry(hAcc,  "Akzeptanz", "l");

      hEfficiency->SetXTitle(pt_str);
      hEfficiency->SetYTitle("Korrekturfaktor");

      hEfficiency->SetLineColor(kBlack);
      hEfficiency->SetLineWidth(5);
      hAcc->SetMarkerColor(kRed+2);
      hAcc->SetLineWidth(5);
      hAcc->SetLineColor(kRed+2);


      hEfficiency->Draw("AXIS");
      hEfficiency->GetYaxis()->SetTitleOffset(0.09);
      hEfficiency->Draw("SAME");
      hAcc->Draw("SAME");
      lCorrection->Draw("SAME");
      c1->Update();
      c1->SetLeftMargin(1.2);

      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkg3to8/Korrekturfaktoren." + PicFormat));
      }
      c1->Clear();


    }

    delete pad1InvMass;
    delete pad2InvMass;
    delete canInvMass;
    delete hData;
    delete hPol1PeakFactor;
    delete line_0;
    delete line_p1;
    delete line_m1;
    delete line_p3;
    delete line_m3;
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

    delete c1;
    delete c2;
    OutputFile->Close();
  }
}
