#include "CommonHeader.h"
//MCTemplatesAnData

// wpsid = which picture should I draw
void IterTempPlot(int binnumber = 3, TString wpsid = "all", TString PicFormat = "png", TString SaveFile = "MCTemplatesAnData"){

  TString str;
  const Int_t nbins = numberbins;
  const Int_t ndrawpoints = 1.e5;
  const int n_iter = 4;
  Double_t somelist[2] = {1., 2.};
  TString doubletempstring = "double temp.";
  TString pol1string = "signal temp. + 1^{st} ord. pol.";

  //////////////////////////////////////////////////////////////////////////////
  // setting up the canvas to draw on. Will later be changed for the chi2 pic
  TCanvas *c1 = new TCanvas("c1","",1200,1000);
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

  TCanvas *c2 = new TCanvas("c2","",1200,1000);
  c2->cd();
  c2->SetTopMargin(0.05);
  c2->SetBottomMargin(0.09);
  c2->SetRightMargin(0.15);
  c2->SetLeftMargin(0.09);
  c2->SetTicky();
  c2->SetTickx();
  c2->SetLogz(1);

  TCanvas *canInvMass = new TCanvas("canInvMass","",1200,1200);
  TPad *pad1InvMass = new TPad("pad1InvMass","",0.0,0.33,1.0,1.0);
  pad1InvMass->SetTopMargin(0.05);
  pad1InvMass->SetLeftMargin(0.09);
  pad1InvMass->SetBottomMargin(0.0);
  pad1InvMass->SetRightMargin(0.02);
  TPad *pad2InvMass = new TPad("pad2InvMass","",0.0,0.0,1.0,0.33);
  pad2InvMass->SetTopMargin(0.0);
  pad2InvMass->SetLeftMargin(0.09);
  pad2InvMass->SetBottomMargin(0.3);
  pad2InvMass->SetRightMargin(0.02);
  pad2InvMass->SetTicky();

  // TF1* fit_eq_double_temp = new TF1("fit_eq_double_temp", "PeakAKorrBG(x)*[1] + mc_full_func1(x)*[0]", 0.0,0.4);
  // fit_eq_double_temp->SetNpx(ndrawpoints);
  // fit_eq_double_temp->SetNumberFitPoints(nbins);
  // fit_eq_double_temp->SetLineColor(kTeal-7);
  // fit_eq_double_temp->SetLineWidth(4);


  // TF1* fit_eq_1 = new TF1("fit_eq_1", "mc_full_func1(x)*[0]+[2]+x*[3]",0.0,0.4);
  // fit_eq_1->SetNpx(ndrawpoints);
  // fit_eq_1->SetNumberFitPoints(nbins);
  // fit_eq_1->SetLineColor(kRed);
  // fit_eq_1->SetLineWidth(4);

  //////////////////////////////////////////////////////////////////////////////
  // setting up the 2 Histograms to compare chi2 from the to fit methods as
  // well as peak factor comp. between pol 1 and double temp fit and the ratio
  // of BG. scaling factor and the Peak scaling factor

  TH1D* hChi2_DT_Iter = NULL;
  TH1D* hChi2_pol1 = NULL;
  TH1D* hPeakRatio = NULL;
  TH1D* hPeakComp = NULL;
  TH1D* hRatioDoubleTemp = NULL;
  TH1D* hRatioDoubleTemp_chi2map = NULL;
  TH1D* hRatioPol1 = NULL;
  TH1D* hData = NULL;
  TH1D* hData_Pol1Error = NULL;
  TH1D* hData_DTError = NULL;
  TH1D* hDTPeak = NULL;
  TH1D* hDTBG = NULL;
  TH1D* hPol1Peak = NULL;
  TH1D* hDoubleTemplatePeakFactor = NULL;
  TH1D* hDoubleTemplatecorrBGFactor = NULL;
  TH1D* hPol1PeakFactor = NULL;
  TH1D* hYield_dt_uncorr = NULL;
  TH1D* hYield_pol1_uncorr = NULL;
  TH1D* hYield_dt_chi2map_uncorr = NULL;
  TH1D* hDoubleTemp = NULL;
  TH1D* hPol1 = NULL;
  TH1D* hChi2_dt_iter = NULL;
  TH1D* hChi2_pol1_iter = NULL;
  TH1D* hChi2_dt_iter_test = NULL;
  TH1D* hChi2_dt_iter_selfcalc = NULL;
  TH2D* hChi2_2D = NULL;
  TH2D* hChi2_2D_sigma = NULL;
  TH1D* hChi2_DT_Chi2map = NULL;
  TH1D* histoChi2_0 = NULL;
  TH1D* hSignalAreaScaling = NULL;
  TH1D* hCorrbackAreaScaling = NULL;
  TH1D* hYield_dt_corrected = NULL;
  TH1D* hYield_dt_chi2map_corrected = NULL;
  TH1D* hYield_pol1_corrected = NULL;
  TH1D* hCorrectedYieldTrueEff = NULL;
  TH1D* h_x_min = NULL;
  TH1D* h_y_min = NULL;
  TH1D* hErrXlow = NULL;
  TH1D* hErrXhigh = NULL;
  TH1D* hErrYlow = NULL;
  TH1D* hErrYhigh = NULL;
  TH1D* hSignal = NULL;
  TH1D* hCorrBack = NULL;
  TF1* fpol1 = NULL;
  TF1* f_ChiOverNdf = NULL;
  TLine* fitrange2 = NULL;
  TLine* line_0 = NULL;
  TLine* line_p1 = NULL;
  TLine* line_m1 = NULL;
  TLine* line_p3 = NULL;
  TLine* line_m3 = NULL;
  TLine* line_one = NULL;

  Double_t line_y = 0;


  TFile* IterTemp = SafelyOpenRootfile("IterTemp.root");
  if (IterTemp->IsOpen() ) printf("IterTemp opened successfully\n");

  hChi2_DT_Iter = (TH1D*) IterTemp->Get(Form("hChi2_DT_Iter"));
  hChi2_pol1 = (TH1D*) IterTemp->Get(Form("hchi2_pol1"));
  hPeakRatio = (TH1D*) IterTemp->Get(Form("hpeakratio"));
  hPeakComp = (TH1D*) IterTemp->Get(Form("hpeakcomp"));
  hDoubleTemplatePeakFactor = (TH1D*) IterTemp->Get(Form("DoubleTemplatePeakFactor"));
  hDoubleTemplatecorrBGFactor = (TH1D*) IterTemp->Get(Form("DoubleTemplatecorrBGFactor"));
  hPol1PeakFactor = (TH1D*) IterTemp->Get(Form("Pol1PeakFactor"));
  hYield_dt_uncorr = (TH1D*) IterTemp->Get(Form("hYield_dt_uncorr"));
  hYield_pol1_uncorr = (TH1D*) IterTemp->Get(Form("hYield_pol1_uncorr"));
  hYield_dt_chi2map_uncorr = (TH1D*) IterTemp->Get(Form("hYield_dt_chi2map_uncorr"));
  hChi2_DT_Chi2map = (TH1D*)IterTemp->Get("hChi2_DT_Chi2map");
  histoChi2_0 = (TH1D*) IterTemp->Get("histoChi2_0");
  hSignalAreaScaling = (TH1D*)IterTemp->Get("hSignalAreaScaling");
  hCorrbackAreaScaling = (TH1D*)IterTemp->Get("hCorrbackAreaScaling");
  hYield_dt_corrected = (TH1D*)IterTemp->Get("hYield_dt_corrected");
  hYield_pol1_corrected = (TH1D*)IterTemp->Get("hYield_pol1_corrected");
  hYield_dt_chi2map_corrected = (TH1D*)IterTemp->Get("hYield_dt_chi2map_corrected");
  hCorrectedYieldTrueEff = (TH1D*) IterTemp->Get("hCorrectedYieldTrueEff");
  h_x_min = (TH1D*)IterTemp->Get("h_x_min");
  h_y_min = (TH1D*)IterTemp->Get("h_y_min");


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
  line_one = new TLine(0.0, 1.0, 20.0, 1.0);
  line_one->SetLineWidth(2);
  line_one->SetLineStyle(3);
  line_one->SetLineColor(kBlack);




  //////////////////////////////////////////////////////////////////////////////
  // going over all pt bins despite first one, which is some framework bs.
  for (int k = 1; k < numberbins; k++) {

    if(binnumber <=  0 || binnumber > numberbins){
      hData                    = (TH1D*) IterTemp->Get(Form("data_bin%02i",k));
      str                      = hData->GetTitle();
      hData->SetTitle("");
      hData_Pol1Error          = (TH1D*) IterTemp->Get(Form("data_addedErrosPol1_bin%02i",k));
      hData_DTError            = (TH1D*) IterTemp->Get(Form("data_addedErrosDT_bin%02i",k));
      hPol1Peak                = (TH1D*) IterTemp->Get(Form("mc_peak_pol1_bin%02i",k));
      hDTPeak                  = (TH1D*) IterTemp->Get(Form("mc_full_DT_bin%02i",k));
      hDTBG                    = (TH1D*) IterTemp->Get(Form("korrBG_bin%02i",k));
      fpol1                    = (TF1*) IterTemp->Get(Form("fpol1_bin%02i",k));
      hRatioDoubleTemp         = (TH1D*) IterTemp->Get(Form("hRatioDoubleTemp_bin%02i",k));
      hRatioDoubleTemp_chi2map = (TH1D*) IterTemp->Get(Form("hRatioDoubleTemp_chi2map_bin%02i",k));
      hRatioPol1               = (TH1D*) IterTemp->Get(Form("hRatioPol1_bin%02i",k));
      // mc_full_clone1        = (TH1D*) IterTemp->Get(Form("mc_full_clone_beforeIterFit_bin%02d",k));
      // korrBG_clone1         = (TH1D*) IterTemp->Get(Form("korrBG_clone_beforeIterFit_bin%02d",k));
      hDoubleTemp              = (TH1D*) IterTemp->Get(Form("hDoubleTemp_bin%02d",k));
      hPol1                    = (TH1D*) IterTemp->Get(Form("hPol1_bin%02d",k));
      hChi2_dt_iter            = (TH1D*) IterTemp->Get(Form("hChi2_dt_iter_bin%02d",k));
      hChi2_dt_iter_test       = (TH1D*) IterTemp->Get(Form("hChi2_dt_iter_test_bin%02d",k));
      hChi2_pol1_iter          = (TH1D*) IterTemp->Get(Form("hChi2_pol1_iter_bin%02d",k));
      hChi2_2D                 = (TH2D*) IterTemp->Get(Form("hChi2_2Dbin%02d",k));
      hChi2_2D_sigma           = (TH2D*) IterTemp->Get(Form("hChi2_2D_sigma_bin%02d",k));
      hSignal                  = (TH1D*) IterTemp->Get(Form("hSignal_bin%02d",k));
      hCorrBack                = (TH1D*) IterTemp->Get(Form("hCorrBack_bin%02d",k));
      f_ChiOverNdf             = (TF1*) IterTemp->Get(Form("f_ChiOverNdf%02d",k));
      hChi2_dt_iter_selfcalc   = (TH1D*) IterTemp->Get(Form("hChi2_dt_iter_selfcalc_bin%02d",k));

      hChi2_dt_iter->SetLineWidth(3);
      hChi2_pol1_iter->SetLineWidth(3);
      hChi2_dt_iter_test->SetLineColor(kRed);

      // mc_full_clone1->SetName("mc_full_clone1");
      // korrBG_clone1->SetName("korrBG_clone1");

      // fit_eq_double_temp->SetParameter(0,hDoubleTemplatePeakFactor->GetBinContent(k+1));
      // fit_eq_double_temp->SetParameter(1,hDoubleTemplatecorrBGFactor->GetBinContent(k+1));
      // fit_eq_1->SetParameter(0,hPol1PeakFactor->GetBinContent(k+1));
      // fit_eq_1->SetParameter(2,fpol1->GetParameter(0));
      // fit_eq_1->SetParameter(3,fpol1->GetParameter(1));



      if(wpsid == "all" || wpsid.Contains("paramcomp")){
        canInvMass->cd();
        pad1InvMass->Draw();
        pad2InvMass->Draw("same");
        pad1InvMass->cd();

        TLegend* leg = new TLegend(0.5,0.35,0.9,0.55);
        SetLegendSettigns(leg, 0.03*3./2.);
        leg->AddEntry(hData, strData, "p");
        leg->AddEntry((TObject*)0x0, "parametrization:", "");
        leg->AddEntry(hPol1, pol1string, "p");
        leg->AddEntry(hDoubleTemp, doubletempstring, "p");

        SetHistoStandardSettings(hData, 0., 0., 0.03*3./2.);
        hData->GetYaxis()->SetRangeUser(
          1.5*hData->GetBinContent(hData->GetMinimumBin()),
          1.1*hData->GetBinContent(hData->GetMaximumBin()));

        hData->SetTitle("");
        hData->GetYaxis()->SetTitleOffset(0.9);
        hData->Draw("p");
        hPol1->Draw("same");
        hDoubleTemp->Draw("same");
        canInvMass->Update();
        line_y = gPad->GetUymax()*0.995;
        fitrange2 = new TLine(lowerparamrange, line_y, upperparamrange, line_y);
        fitrange2->SetLineColor(kAzure+10);
        fitrange2->SetLineWidth(7);
        fitrange2->Draw("same");
        leg->AddEntry(fitrange2, "range", "l");
        leg->Draw("same");
        DrawLabelALICE(0.5, 0.9, 0.035, 0.03*3./2., str);
        pad1InvMass->Update();


        pad2InvMass->cd();
        // hRatioDoubleTemp->DrawCopy("P");
        line_0->Draw("same");
        line_p1->Draw("same");
        line_m1->Draw("same");
        line_p3->Draw("same");
        line_m3->Draw("same");
        hRatioDoubleTemp->DrawCopy("P");
        hRatioPol1->DrawCopy("SAME P");
        pad2InvMass->Update();

        canInvMass->Update();
        canInvMass->SaveAs(Form(SaveFile + "/DataFitWithMCCompIter%02i." + PicFormat,k));
        canInvMass->Clear("D");

        delete leg;

      }

      if(wpsid == "all" || wpsid.Contains("dtcomp")){
        ////////////////////////////////////////////////////////////////////////
        // Comparison between the two double temp methods with data
        c1->cd();

        TH1D* hSignal_Clone = NULL;
        TH1D* hCorrBack_Clone = NULL;
        hSignal_Clone = (TH1D*) hSignal->Clone("hSignal_Clone");
        hCorrBack_Clone = (TH1D*) hCorrBack->Clone("hCorrBack_Clone");
        hSignal_Clone->Scale(hSignalAreaScaling->GetBinContent(k)*h_x_min->GetBinContent(k+1));
        hCorrBack_Clone->Scale(hCorrbackAreaScaling->GetBinContent(k)*h_y_min->GetBinContent(k+1));
        hSignal_Clone->Add(hCorrBack_Clone);
        hSignal_Clone->SetMarkerStyle(20);
        hSignal_Clone->SetMarkerSize(1.5);
        hSignal_Clone->SetMarkerColor(kMagenta+2);
        hSignal_Clone->SetLineColor(kMagenta+2);

        TLegend* leg = new TLegend(0.5,0.35,0.9,0.55);
        SetLegendSettigns(leg, 0.03);
        leg->AddEntry(hData, strData, "p");
        leg->AddEntry((TObject*)0x0, "parametrization:", "");
        leg->AddEntry(hSignal_Clone, doubletempstring + " with chi2map", "p");
        leg->AddEntry(hDoubleTemp, doubletempstring, "p");

        SetHistoStandardSettings(hData, 0., 0., 0.03);
        hData->SetTitle("");
        hData->GetYaxis()->SetTitleOffset(0.9);
        hData->Draw("p");
        hDoubleTemp->Draw("same");
        hSignal_Clone->Draw("same");
        c1->Update();
        line_y = gPad->GetUymax()*0.995;
        fitrange2 = new TLine(lowerparamrange, line_y, upperparamrange, line_y);
        fitrange2->SetLineColor(kAzure+10);
        fitrange2->SetLineWidth(7);
        fitrange2->Draw("same");
        leg->AddEntry(fitrange2, "range", "l");
        leg->Draw("same");
        DrawLabelALICE(0.5, 0.9, 0.03, 0.03, str);
        c1->Update();
        c1->SaveAs(Form(SaveFile + "/DataFitDTComp%02i." + PicFormat,k));
        c1->Clear("D");


        delete leg;
        // delete hSignal_Clone;
        // delete hCorrBack_Clone;
      }

      if(wpsid == "all" || wpsid.Contains("bgcomp")){
        //////////////////////////////////////////////////////////////////////
        // Drawing both corr. BG versions to Data with normal errors

        TH1D* hCorrBack_Chi2Map = (TH1D*) hCorrBack->Clone("hCorrBack_Chi2Map");
        hCorrBack_Chi2Map->Scale(h_y_min->GetBinContent(k+1));
        hCorrBack_Chi2Map->SetLineColor(kMagenta+2);
        hCorrBack_Chi2Map->SetMarkerColor(kMagenta+2);

        c1->cd();
        TLegend* leg = new TLegend(0.5,0.5,0.9,0.63);
        SetLegendSettigns(leg, 0.03);
        leg->AddEntry(fpol1, "1^{st} ord. pol.", "l");
        leg->AddEntry(hDTBG, "scaled corr. back. temp.", "p");
        leg->AddEntry(hCorrBack_Chi2Map, "scaled corr. back. temp with chi2map");

        hData->GetXaxis()->SetTitleSize(0.03);
        hData->GetYaxis()->SetTitleSize(0.03);
        hData->GetXaxis()->SetLabelSize(0.03);
        hData->GetYaxis()->SetLabelSize(0.03);

        hCorrBack_Chi2Map->Draw("");
        fpol1->Draw("same");
        hDTBG->Draw("same");
        c1->Update();
        line_y = gPad->GetUymax()*0.995;
        TLine* fitrange = new TLine(lowerparamrange, line_y, upperparamrange, line_y);

        fitrange->SetLineColor(kAzure+10);
        fitrange->SetLineWidth(7);
        fitrange->Draw("same");
        leg->AddEntry(fitrange, "param. range", "l");
        leg->Draw("same");
        DrawLabelALICE(0.5, 0.9, 0.02, 0.03, str);

        c1->Update();
        c1->SaveAs(Form(SaveFile + "/CorrBGComp%02i." + PicFormat,k));
        c1->Clear();

        delete leg;
        delete fitrange;
      }

      if(wpsid == "all" || wpsid.Contains("monitoring")){
        ////////////////////////////////////////////////////////////////////////
        // Drawing monitoring plots for Chi^2 for both fits
        c1->cd();

        TLegend* leg = new TLegend(0.5,0.5,0.9,0.63);
        SetLegendSettigns(leg, 0.03);
        leg->AddEntry(hChi2_dt_iter_selfcalc, doubletempstring + " self calculated", "l");
        leg->AddEntry(f_ChiOverNdf, doubletempstring + " Chi2Map");

        hChi2_dt_iter_selfcalc->Draw("HIST");
        f_ChiOverNdf->Draw("SAME");
        leg->Draw("same");
        DrawLabelALICE(0.5, 0.45, 0.02, 0.03, str);


        c1->Update();
        c1->SaveAs(Form(SaveFile + "/Monitoring%02i." + PicFormat,k));
        c1->Clear();

        delete leg;

        ////////////////////////////////////////////////////////////////////////
        // Drawing monitoring plots for Chi^2 for both fits
        c1->cd();

        TLegend* leg3 = new TLegend(0.5,0.5,0.9,0.63);
        SetLegendSettigns(leg3, 0.03);
        leg3->AddEntry(hChi2_dt_iter, doubletempstring, "l");
        leg3->AddEntry(hChi2_dt_iter_selfcalc, doubletempstring + " self calculated", "l");
        leg3->AddEntry(hChi2_dt_iter_test, "Chi2Test Function" , "l");
        hChi2_dt_iter->GetYaxis()->SetRangeUser(0.,1.5*hChi2_dt_iter->GetMaximum());
        hChi2_dt_iter->GetXaxis()->SetRangeUser(0.5, 4.5);

        hChi2_dt_iter->Draw("HIST");
        hChi2_dt_iter_selfcalc->Draw("SAME HIST");
        hChi2_dt_iter_test->Draw("SAME HIST");
        c1->Update();
        leg3->Draw("same");
        DrawLabelALICE(0.5, 0.9, 0.02, 0.03, str);

        c1->Update();
        c1->SaveAs(Form(SaveFile + "/MonitoringLesserXRangep%02i." + PicFormat,k));
        c1->Clear();

        delete leg3;
        hChi2_dt_iter->GetXaxis()->SetRangeUser(0,0);


        ////////////////////////////////////////////////////////////////////////
        // Drawing monitoring plots for Chi^2 for both fits
        c1->cd();

        TLegend* leg2 = new TLegend(0.5,0.5,0.9,0.63);
        SetLegendSettigns(leg2, 0.03);
        leg2->AddEntry(hChi2_dt_iter, doubletempstring, "l");
        leg2->AddEntry(hChi2_dt_iter_selfcalc, doubletempstring + " self calculated", "l");
        leg2->AddEntry(hChi2_dt_iter_test, "Chi2Test Function" , "l");
        hChi2_dt_iter->GetYaxis()->SetRangeUser(0.,1.5*hChi2_dt_iter->GetMaximum());

        hChi2_dt_iter->Draw("HIST");
        hChi2_dt_iter_selfcalc->Draw("SAME HIST");
        hChi2_dt_iter_test->Draw("SAME HIST");
        c1->Update();
        leg2->Draw("same");
        DrawLabelALICE(0.5, 0.9, 0.02, 0.03, str);

        c1->Update();
        c1->SaveAs(Form(SaveFile + "/MonitoringChi2DTComp%02i." + PicFormat,k));
        c1->Clear();

        delete leg2;

      }

      if(wpsid == "all" || wpsid.Contains("chi2map")){
        ////////////////////////////////////////////////////////////////////////
        // Drawing Chi2 maps
        c2->cd();

        hChi2_2D->Draw("colz");
        hChi2_2D_sigma->SetLineColor(kWhite);
        hChi2_2D_sigma->SetLineWidth(2);
        hChi2_2D_sigma->SetContour(2, somelist);
        hChi2_2D_sigma->Draw("same cont3");
        c2->Update();

        c2->Update();
        c2->SaveAs(Form(SaveFile + "/Chi2Map%02i." + PicFormat,k));
        c2->Clear();
      }

        if(wpsid == "all" || wpsid.Contains("chi2map")){
          ////////////////////////////////////////////////////////////////////////
          // Drawing the Plot coming from the chi2map data
          c1->cd();
          TH1D* hSignal_Clone = NULL;
          TH1D* hCorrBack_Clone = NULL;
          hSignal_Clone = (TH1D*) hSignal->Clone("hSignal_Clone");
          hCorrBack_Clone = (TH1D*) hCorrBack->Clone("hCorrBack_Clone");
          hSignal_Clone->Scale(hSignalAreaScaling->GetBinContent(k)*h_x_min->GetBinContent(k+1));
          hCorrBack_Clone->Scale(hCorrbackAreaScaling->GetBinContent(k)*h_y_min->GetBinContent(k+1));
          hSignal_Clone->Add(hCorrBack_Clone);

          hData->Draw("");
          hSignal_Clone->Draw("same");

          c1->Update();

          c1->Update();
          c1->SaveAs(Form(SaveFile + "/ParamResultWithChi2Map%02i." + PicFormat,k));
          c1->Clear();

          // delete hSignal_Clone;
          // delete hCorrBack_Clone;
        }
    }

      else{
        if(binnumber == k){
          hData                    = (TH1D*) IterTemp->Get(Form("data_bin%02i",k));
          str                      = hData->GetTitle();
          hData->SetTitle("");
          hData_Pol1Error          = (TH1D*) IterTemp->Get(Form("data_addedErrosPol1_bin%02i",k));
          hData_DTError            = (TH1D*) IterTemp->Get(Form("data_addedErrosDT_bin%02i",k));
          hPol1Peak                = (TH1D*) IterTemp->Get(Form("mc_peak_pol1_bin%02i",k));
          hDTPeak                  = (TH1D*) IterTemp->Get(Form("mc_full_DT_bin%02i",k));
          hDTBG                    = (TH1D*) IterTemp->Get(Form("korrBG_bin%02i",k));
          fpol1                    = (TF1*) IterTemp->Get(Form("fpol1_bin%02i",k));
          hRatioDoubleTemp         = (TH1D*) IterTemp->Get(Form("hRatioDoubleTemp_bin%02i",k));
          hRatioDoubleTemp_chi2map = (TH1D*) IterTemp->Get(Form("hRatioDoubleTemp_chi2map_bin%02i",k));
          hRatioPol1               = (TH1D*) IterTemp->Get(Form("hRatioPol1_bin%02i",k));
          // mc_full_clone1        = (TH1D*) IterTemp->Get(Form("mc_full_clone_beforeIterFit_bin%02d",k));
          // korrBG_clone1         = (TH1D*) IterTemp->Get(Form("korrBG_clone_beforeIterFit_bin%02d",k));
          hDoubleTemp              = (TH1D*) IterTemp->Get(Form("hDoubleTemp_bin%02d",k));
          hPol1                    = (TH1D*) IterTemp->Get(Form("hPol1_bin%02d",k));
          hChi2_dt_iter            = (TH1D*) IterTemp->Get(Form("hChi2_dt_iter_bin%02d",k));
          hChi2_dt_iter_test       = (TH1D*) IterTemp->Get(Form("hChi2_dt_iter_test_bin%02d",k));
          hChi2_pol1_iter          = (TH1D*) IterTemp->Get(Form("hChi2_pol1_iter_bin%02d",k));
          hChi2_2D                 = (TH2D*) IterTemp->Get(Form("hChi2_2Dbin%02d",k));
          hChi2_2D_sigma           = (TH2D*) IterTemp->Get(Form("hChi2_2D_sigma_bin%02d",k));
          hSignal                  = (TH1D*) IterTemp->Get(Form("hSignal_bin%02d",k));
          hCorrBack                = (TH1D*) IterTemp->Get(Form("hCorrBack_bin%02d",k));
          f_ChiOverNdf             = (TF1*) IterTemp->Get(Form("f_ChiOverNdf%02d",k));
          hChi2_dt_iter_selfcalc   = (TH1D*) IterTemp->Get(Form("hChi2_dt_iter_selfcalc_bin%02d",k));
          // mc_full_clone1->SetName("mc_full_clone1");
          // korrBG_clone1->SetName("korrBG_clone1");

          // fit_eq_double_temp->SetParameter(0,hDoubleTemplatePeakFactor->GetBinContent(k+1));
          // fit_eq_double_temp->SetParameter(1,hDoubleTemplatecorrBGFactor->GetBinContent(k+1));
          // fit_eq_1->SetParameter(0,hPol1PeakFactor->GetBinContent(k+1));
          // fit_eq_1->SetParameter(2,fpol1->GetParameter(0));
          // fit_eq_1->SetParameter(3,fpol1->GetParameter(1));


          hChi2_dt_iter->SetLineWidth(3);
          hChi2_pol1_iter->SetLineWidth(3);
          hChi2_dt_iter_test->SetLineColor(kRed);


          if(wpsid == "all" || wpsid.Contains("paramcomp")){


            canInvMass->cd();
            pad1InvMass->Draw();
            pad2InvMass->Draw("same");
            pad1InvMass->cd();

            TLegend* leg = new TLegend(0.5,0.35,0.9,0.55);
            SetLegendSettigns(leg, 0.03*3./2.);
            leg->AddEntry(hData, strData, "p");
            leg->AddEntry((TObject*)0x0, "parametrization:", "");
            leg->AddEntry(hPol1, pol1string, "p");
            leg->AddEntry(hDoubleTemp, doubletempstring, "p");

            SetHistoStandardSettings(hData, 0., 0., 0.03*3./2.);
            hData->GetYaxis()->SetRangeUser(
              1.5*hData->GetBinContent(hData->GetMinimumBin()),
              1.1*hData->GetBinContent(hData->GetMaximumBin()));

            hData->SetTitle("");
            hData->GetYaxis()->SetTitleOffset(0.9);
            hData->Draw("p");
            hPol1->Draw("same");
            hDoubleTemp->Draw("same");
            canInvMass->Update();
            line_y = gPad->GetUymax()*0.995;
            fitrange2 = new TLine(lowerparamrange, line_y, upperparamrange, line_y);
            fitrange2->SetLineColor(kAzure+10);
            fitrange2->SetLineWidth(7);
            fitrange2->Draw("same");
            leg->AddEntry(fitrange2, "range", "l");
            leg->Draw("same");
            DrawLabelALICE(0.5, 0.9, 0.035, 0.03*3./2., str);
            pad1InvMass->Update();


            pad2InvMass->cd();
            hRatioDoubleTemp->DrawCopy("P");
            line_0->Draw("same");
            line_p1->Draw("same");
            line_m1->Draw("same");
            line_p3->Draw("same");
            line_m3->Draw("same");
            hRatioDoubleTemp->DrawCopy("SAME P");
            hRatioPol1->DrawCopy("SAME P");
            pad2InvMass->Update();

            canInvMass->Update();
            canInvMass->SaveAs(Form(SaveFile + "/DataFitWithMCCompIter%02i." + PicFormat,k));
            canInvMass->Clear("D");

            delete leg;

          }

          if(wpsid == "all" || wpsid.Contains("dtcomp")){
            ////////////////////////////////////////////////////////////////////////
            // Comparison between the two double temp methods with data
            c1->cd();

            TH1D* hSignal_Clone = NULL;
            TH1D* hCorrBack_Clone = NULL;
            hSignal_Clone = (TH1D*) hSignal->Clone("hSignal_Clone");
            hCorrBack_Clone = (TH1D*) hCorrBack->Clone("hCorrBack_Clone");
            hSignal_Clone->Scale(hSignalAreaScaling->GetBinContent(k)*h_x_min->GetBinContent(k+1));
            hCorrBack_Clone->Scale(hCorrbackAreaScaling->GetBinContent(k)*h_y_min->GetBinContent(k+1));
            hSignal_Clone->Add(hCorrBack_Clone);
            hSignal_Clone->SetMarkerStyle(20);
            hSignal_Clone->SetMarkerSize(1.5);
            hSignal_Clone->SetMarkerColor(kMagenta+2);
            hSignal_Clone->SetLineColor(kMagenta+2);

            TLegend* leg = new TLegend(0.5,0.35,0.9,0.55);
            SetLegendSettigns(leg, 0.03);
            leg->AddEntry(hData, strData, "p");
            leg->AddEntry((TObject*)0x0, "parametrization:", "");
            leg->AddEntry(hSignal_Clone, doubletempstring + " with chi2map", "p");
            leg->AddEntry(hDoubleTemp, doubletempstring, "p");

            SetHistoStandardSettings(hData, 0., 0., 0.03);
            hData->SetTitle("");
            hData->GetYaxis()->SetTitleOffset(0.9);
            hData->Draw("p");
            hDoubleTemp->Draw("same");
            hSignal_Clone->Draw("same");
            c1->Update();
            line_y = gPad->GetUymax()*0.995;
            fitrange2 = new TLine(lowerparamrange, line_y, upperparamrange, line_y);
            fitrange2->SetLineColor(kAzure+10);
            fitrange2->SetLineWidth(7);
            fitrange2->Draw("same");
            leg->AddEntry(fitrange2, "range", "l");
            leg->Draw("same");
            DrawLabelALICE(0.5, 0.9, 0.03, 0.03, str);
            c1->Update();
            c1->SaveAs(Form(SaveFile + "/DataFitDTComp%02i." + PicFormat,k));
            c1->Clear("D");


            delete leg;
            // delete hSignal_Clone;
            // delete hCorrBack_Clone;
          }


          if(wpsid == "all" || wpsid.Contains("bgcomp")){
            //////////////////////////////////////////////////////////////////////
            // Drawing both corr. BG versions to Data with normal errors

            TH1D* hCorrBack_Chi2Map = (TH1D*) hCorrBack->Clone("hCorrBack_Chi2Map");
            hCorrBack_Chi2Map->Scale(h_y_min->GetBinContent(k+1));
            hCorrBack_Chi2Map->SetLineColor(kMagenta+2);
            hCorrBack_Chi2Map->SetMarkerColor(kMagenta+2);

            c1->cd();
            TLegend* leg = new TLegend(0.5,0.5,0.9,0.63);
            SetLegendSettigns(leg, 0.03);
            leg->AddEntry(fpol1, "1^{st} ord. pol.", "l");
            leg->AddEntry(hDTBG, "scaled corr. back. temp.", "p");
            leg->AddEntry(hCorrBack_Chi2Map, "scaled corr. back. temp with chi2map");

            hData->GetXaxis()->SetTitleSize(0.03);
            hData->GetYaxis()->SetTitleSize(0.03);
            hData->GetXaxis()->SetLabelSize(0.03);
            hData->GetYaxis()->SetLabelSize(0.03);

            hCorrBack_Chi2Map->Draw("");
            fpol1->Draw("same");
            hDTBG->Draw("same");
            c1->Update();
            line_y = gPad->GetUymax()*0.995;
            TLine* fitrange = new TLine(lowerparamrange, line_y, upperparamrange, line_y);

            fitrange->SetLineColor(kAzure+10);
            fitrange->SetLineWidth(7);
            fitrange->Draw("same");
            leg->AddEntry(fitrange, "param. range", "l");
            leg->Draw("same");
            DrawLabelALICE(0.5, 0.9, 0.02, 0.03, str);

            c1->Update();
            c1->SaveAs(Form(SaveFile + "/CorrBGComp%02i." + PicFormat,k));
            c1->Clear();

            delete leg;
            delete fitrange;
          }

          if(wpsid == "all" || wpsid.Contains("monitoring")){
            ////////////////////////////////////////////////////////////////////////
            // Drawing monitoring plots for Chi^2 for both fits
            c1->cd();

            TLegend* leg = new TLegend(0.5,0.5,0.9,0.63);
            SetLegendSettigns(leg, 0.03);
            leg->AddEntry(hChi2_dt_iter_selfcalc, doubletempstring + " self calculated", "l");
            leg->AddEntry(f_ChiOverNdf, doubletempstring + " Chi2Map");

            hChi2_dt_iter_selfcalc->Draw("HIST");
            f_ChiOverNdf->Draw("SAME");
            leg->Draw("same");
            DrawLabelALICE(0.5, 0.45, 0.02, 0.03, str);


            c1->Update();
            c1->SaveAs(Form(SaveFile + "/Monitoring%02i." + PicFormat,k));
            c1->Clear();

            delete leg;

            ////////////////////////////////////////////////////////////////////////
            // Drawing monitoring plots for Chi^2 for both fits
            c1->cd();

            TLegend* leg3 = new TLegend(0.5,0.5,0.9,0.63);
            SetLegendSettigns(leg3, 0.03);
            leg3->AddEntry(hChi2_dt_iter, doubletempstring, "l");
            leg3->AddEntry(hChi2_dt_iter_selfcalc, doubletempstring + " self calculated", "l");
            leg3->AddEntry(hChi2_dt_iter_test, "Chi2Test Function" , "l");
            hChi2_dt_iter->GetYaxis()->SetRangeUser(0.,1.5*hChi2_dt_iter->GetMaximum());
            hChi2_dt_iter->GetXaxis()->SetRangeUser(0.5, 4.5);

            hChi2_dt_iter->Draw("HIST");
            hChi2_dt_iter_selfcalc->Draw("SAME HIST");
            hChi2_dt_iter_test->Draw("SAME HIST");
            c1->Update();
            leg3->Draw("same");
            DrawLabelALICE(0.5, 0.9, 0.02, 0.03, str);

            c1->Update();
            c1->SaveAs(Form(SaveFile + "/MonitoringChi2DTCompLesserXRange%02i." + PicFormat,k));
            c1->Clear();

            delete leg3;
            hChi2_dt_iter->GetXaxis()->SetRangeUser(0,0);


            ////////////////////////////////////////////////////////////////////////
            // Drawing monitoring plots for Chi^2 for both fits
            c1->cd();

            TLegend* leg2 = new TLegend(0.5,0.5,0.9,0.63);
            SetLegendSettigns(leg2, 0.03);
            leg2->AddEntry(hChi2_dt_iter, doubletempstring, "l");
            leg2->AddEntry(hChi2_dt_iter_selfcalc, doubletempstring + " self calculated", "l");
            leg2->AddEntry(hChi2_dt_iter_test, "Chi2Test Function" , "l");
            hChi2_dt_iter->GetYaxis()->SetRangeUser(0.,1.5*hChi2_dt_iter->GetMaximum());

            hChi2_dt_iter->Draw("HIST");
            hChi2_dt_iter_selfcalc->Draw("SAME HIST");
            hChi2_dt_iter_test->Draw("SAME HIST");
            c1->Update();
            leg2->Draw("same");
            DrawLabelALICE(0.5, 0.9, 0.02, 0.03, str);

            c1->Update();
            c1->SaveAs(Form(SaveFile + "/MonitoringChi2DTComp%02i." + PicFormat,k));
            c1->Clear();

            delete leg2;

          }

          if(wpsid == "all" || wpsid.Contains("chi2map")){
            ////////////////////////////////////////////////////////////////////////
            // Drawing the Plot coming from the chi2map data
            c1->cd();
            TH1D* hSignal_Clone = NULL;
            TH1D* hCorrBack_Clone = NULL;
            hSignal_Clone = (TH1D*) hSignal->Clone("hSignal_Clone");
            hCorrBack_Clone = (TH1D*) hCorrBack->Clone("hCorrBack_Clone");
            hSignal_Clone->Scale(hSignalAreaScaling->GetBinContent(k)*h_x_min->GetBinContent(k+1));
            hCorrBack_Clone->Scale(hCorrbackAreaScaling->GetBinContent(k)*h_y_min->GetBinContent(k+1));
            hSignal_Clone->Add(hCorrBack_Clone);

            hData->Draw("");
            hSignal_Clone->Draw("same");

            c1->Update();

            c1->Update();
            c1->SaveAs(Form(SaveFile + "/ParamResultWithChi2Map%02i." + PicFormat,k));
            c1->Clear();

            // delete hSignal_Clone;
            // delete hCorrBack_Clone;
          }
        }
      }
  }

  // Drawing of Chi^2 comparison bwtween the two fits
  if(wpsid == "all" || wpsid.Contains("chi2")){

    TLegend* leg = new TLegend(0.6,0.75,0.9,0.9);
    SetLegendSettigns(leg);
    leg->AddEntry(hChi2_DT_Iter, doubletempstring, "l");
    leg->AddEntry(hChi2_pol1, pol1string, "l");
    Double_t chi2_dt_mean = 0;
    Double_t chi2_pol1_mean = 0;
    for (int i = 0; i < numberbins; i++) {
      chi2_dt_mean += hChi2_DT_Iter->GetBinContent(i);
      chi2_pol1_mean += hChi2_pol1->GetBinContent(i);
    }
    chi2_dt_mean /= (Double_t)numberbins;
    chi2_pol1_mean /= (Double_t)numberbins;
    hChi2_DT_Iter->SetMarkerStyle(1);
    hChi2_DT_Iter->SetMarkerSize(1);
    hChi2_pol1->SetMarkerStyle(1);
    hChi2_pol1->SetMarkerSize(1);

    c1->cd();
    c1->Clear();
    hChi2_DT_Iter->Draw("AXIS");
    line_one->Draw("same");
    hChi2_DT_Iter->Draw("SAME HIST");
    hChi2_DT_Iter->Draw("SAME P");
    hChi2_pol1->Draw("SAME HIST");
    hChi2_pol1->Draw("SAME P");
    leg->Draw("same");
    DrawLabelALICE(0.2, 0.9, 0.018, 0.03);

    TLatex* tex = new TLatex();
    SetLatexSettings(tex);
    tex->DrawLatexNDC(0.6,0.7,Form("#LT#chi^{2}_{dt}/ndf#GT = %1.2lf",chi2_dt_mean));
    tex->DrawLatexNDC(0.6,0.65,Form("#LT#chi^{2}_{pol1}/ndf#GT = %1.2lf",chi2_pol1_mean));

    c1->Update();
    c1->SaveAs(Form(SaveFile + "/Chi2." + PicFormat));
    c1->Clear();

    TLegend* leg2 = new TLegend(0.6,0.75,0.9,0.9);
    SetLegendSettigns(leg2);
    leg2->AddEntry(hChi2_pol1, pol1string, "l");
    leg2->AddEntry(hChi2_DT_Chi2map, doubletempstring, "l");
    leg2->AddEntry(histoChi2_0, "parametrization with function", "l");


    hChi2_pol1->SetLineColor(kRed);
    hChi2_DT_Chi2map->SetLineColor(kMagenta+2);
    histoChi2_0->SetLineWidth(3);

    c1->cd();
    c1->Clear();
    c1->SetLogx(1);
    histoChi2_0->Draw("AXIS");
    line_one->Draw("same");
    histoChi2_0->Draw("SAME HIST");
    hChi2_pol1->Draw("SAME HIST");
    hChi2_pol1->Draw("SAME P");
    hChi2_DT_Chi2map->Draw("SAME HIST");
    hChi2_DT_Chi2map->Draw("SAME P");
    leg2->Draw("same");
    DrawLabelALICE(0.2, 0.9, 0.018, 0.03);


    c1->Update();
    c1->SaveAs(Form(SaveFile + "/Chi2Iter_vs_Map." + PicFormat));
    c1->Clear();
    c1->SetLogx(0);






    delete leg;
    delete leg2;
    delete tex;
  }
  // draing of b_double/a_double temp
  if(wpsid == "all" || wpsid.Contains("bgtopeak")){
    c1->cd();
    c1->Clear();

    hPeakRatio->Draw("");
    line_one->Draw("same");
    DrawLabelALICE(0.34, 0.9, 0.018, 0.03);
    c1->Update();
    c1->SaveAs(Form(SaveFile + "/corr_BG_to_peak." + PicFormat));
    c1->Clear();
  }

  // draing of a_pol1/a_double temp
  if(wpsid == "all" || wpsid.Contains("peakcomp")){
    c1->cd();
    c1->Clear();

    hPeakComp->Draw("");
    line_one->Draw("same");
    DrawLabelALICE(0.15, 0.9, 0.018, 0.03);
    c1->Update();
    c1->SaveAs(Form(SaveFile + "/Peakcomp." + PicFormat));
    c1->Clear();
  }
  // drawing uncorrected yields
  if(wpsid == "all" || wpsid.Contains("uncorryield")){

    TLegend* leg = new TLegend(0.2,0.2,0.4,0.4);
    SetLegendSettigns(leg);
    leg->AddEntry(hYield_dt_uncorr, doubletempstring, "lp");
    leg->AddEntry(hYield_dt_chi2map_uncorr, doubletempstring + " with chi2map" , "lp");
    leg->AddEntry(hYield_pol1_uncorr, pol1string, "lp");

    c1->cd();
    c1->SetLeftMargin(0.11);
    c1->Clear();
    c1->SetLogy(1);
    hYield_pol1_uncorr->GetYaxis()->SetRangeUser(1.e-8,1.e-2);
    hYield_pol1_uncorr->GetYaxis()->SetTitleOffset(1.5);
    hYield_pol1_uncorr->Draw("lp");
    hYield_dt_uncorr->Draw("samelp");
    hYield_dt_chi2map_uncorr->Draw("samelp");
    leg->Draw("same");

    DrawLabelALICE(0.55, 0.9, 0.018, 0.03);
    c1->Update();
    c1->SaveAs(Form(SaveFile + "/UncorrYields." + PicFormat));
    c1->Clear();
    c1->SetLogy(0);
    c1->SetLeftMargin(0.09);
    delete leg;
  }

  // drawing uncorrected yields
  if(wpsid == "all" || wpsid.Contains("corryield")){

    hCorrectedYieldTrueEff->SetMarkerSize(1.5);

    // TLegend* leg = new TLegend(0.2,0.2,0.4,0.4);
    // SetLegendSettigns(leg);
    // leg->AddEntry(hYield_dt_corrected, doubletempstring, "lp");
    // leg->AddEntry(hYield_dt_chi2map_corrected, doubletempstring + " with chi2map" , "lp");
    // leg->AddEntry(hYield_pol1_corrected, pol1string, "lp");
    //
    // c1->cd();
    // c1->SetLeftMargin(0.11);
    // c1->Clear();
    // c1->SetLogy(1);
    // hYield_pol1_corrected->GetYaxis()->SetTitleOffset(1.5);
    // hYield_pol1_corrected->Draw("lp");
    // hYield_dt_corrected->Draw("samelp");
    // hYield_dt_chi2map_corrected->Draw("samelp");
    // leg->Draw("same");
    //
    // DrawLabelALICE(0.55, 0.9, 0.018, 0.03);
    // c1->Update();
    // c1->SaveAs(Form(SaveFile + "/CorrYields." + PicFormat));
    // c1->Clear();
    // c1->SetLogy(0);
    // c1->SetLeftMargin(0.09);
    // delete leg;

    ////////////////////////////////////////////////////////////////////////////
    // drwaing yields + ratios
    canInvMass->cd();
    pad1InvMass->Draw();
    pad2InvMass->Draw("same");
    pad1InvMass->cd();

    pad1InvMass->SetLogy(1);

    TLegend* leg = new TLegend(0.5,0.35,0.9,0.55);
    SetLegendSettigns(leg, 0.03*3./2.);
    leg->AddEntry(hYield_dt_corrected, doubletempstring, "lp");
    leg->AddEntry(hYield_dt_chi2map_corrected, doubletempstring + " with chi2map" , "lp");
    leg->AddEntry(hYield_pol1_corrected, pol1string, "lp");
    leg->AddEntry(hCorrectedYieldTrueEff, "standard method from framework", "lp");

    hYield_dt_corrected->GetYaxis()->SetTitleOffset(1.2);
    hYield_dt_corrected->GetYaxis()->SetRangeUser(1.e-7-5.e-8,1.e-1);
    hYield_dt_corrected->GetXaxis()->SetRangeUser(0.0, 16.0);
    hYield_dt_corrected->Draw("p");
    hYield_dt_chi2map_corrected->Draw("same");
    hYield_pol1_corrected->Draw("same");
    hCorrectedYieldTrueEff->Draw("same");
    leg->Draw("same");
    canInvMass->Update();
    DrawLabelALICE(0.5, 0.9, 0.035, 0.03*3./2., "");
    pad1InvMass->Update();

    TH1D* hYield_dt_corrected_ratio = (TH1D*) hCorrectedYieldTrueEff->Clone("hYield_dt_corrected_ratio");
    hYield_dt_corrected_ratio->Divide(hYield_dt_corrected);
    hYield_dt_corrected_ratio->SetLineColor(kTeal-7);
    hYield_dt_corrected_ratio->SetMarkerColor(kTeal-7);
    hYield_dt_corrected_ratio->SetYTitle("Ratio");

    TH1D* hYield_dt_chi2map_corrected_ratio = (TH1D*) hCorrectedYieldTrueEff->Clone("hYield_dt_chi2map_corrected_ratio");
    hYield_dt_chi2map_corrected_ratio->Divide(hYield_dt_chi2map_corrected);
    hYield_dt_chi2map_corrected_ratio->SetLineColor(kMagenta+2);
    hYield_dt_chi2map_corrected_ratio->SetMarkerColor(kMagenta+2);
    hYield_dt_chi2map_corrected_ratio->SetYTitle("Ratio");

    TH1D* hYield_pol1_corrected_ratio = (TH1D*) hCorrectedYieldTrueEff->Clone("hYield_pol1_corrected_ratio");
    hYield_pol1_corrected_ratio->Divide(hYield_pol1_corrected);
    hYield_pol1_corrected_ratio->SetLineColor(kRed);
    hYield_pol1_corrected_ratio->SetMarkerColor(kRed);
    hYield_pol1_corrected_ratio->SetYTitle("Ratio");

    pad2InvMass->cd();
    TLine* line_ratio1 = new TLine(0.0, 1.0, 16.0, 1.0);
    line_ratio1->SetLineWidth(2);
    line_ratio1->SetLineStyle(3);

    hYield_dt_chi2map_corrected_ratio->GetXaxis()->SetRangeUser(0.0, 16.0);
    hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetRangeUser(0.69, 1.25);

    hYield_dt_chi2map_corrected_ratio->DrawCopy("P");
    line_ratio1->Draw("SAME");
    hYield_dt_corrected_ratio->DrawCopy("SAME P");
    hYield_dt_chi2map_corrected_ratio->DrawCopy("SAME P");
    hYield_pol1_corrected_ratio->DrawCopy("SAME P");
    pad2InvMass->Update();

    canInvMass->Update();
    canInvMass->SaveAs(Form(SaveFile + "/CorrectedYieldComp." + PicFormat));
    canInvMass->Clear("D");

    delete leg;
    delete hYield_dt_corrected_ratio;
    delete hYield_dt_chi2map_corrected_ratio;
    delete hYield_pol1_corrected_ratio;
  }

  //////////////////////////////////////////////////////////////////////////////
  // constraint plots
  // if(wpsid == "all" || wpsid.Contains("constraint")){
  //   TH1D* h_y_min_clone = (TH1D*) h_y_min->Clone("h_y_min_clone");
  //   h_y_min_clone->Add(h_x_min,-1);
  //   h_y_min_clone->Draw();
  //   c1->Update();
  //   c1->SaveAs(Form(SaveFile + "/ConstraintPloti." + PicFormat));
  //   c1->Clear();
  //   delete h_y_min_clone;
  // }


  //////////////////////////////////////////////////////////////////////////////
  // Factor comp between Chi2 map and Iterative method
  if(wpsid == "all" || wpsid.Contains("factorcomp")){
    c1->cd();

    h_x_min->SetLineColor(kMagenta+2);
    h_x_min->SetMarkerColor(kMagenta+2);
    h_x_min->SetMarkerStyle(25);
    h_x_min->SetMarkerSize(1.5);
    hDoubleTemplatePeakFactor->SetMarkerStyle(25);

    TLegend* leg = new TLegend(0.15,0.7,0.6,0.9);
    SetLegendSettigns(leg);
    leg->AddEntry(hDoubleTemplatePeakFactor, "signal scaling factor" + doubletempstring, "lp");
    leg->AddEntry(h_x_min, "signal scaling factor" +  doubletempstring + " with chi2map" , "lp");


    h_x_min->Draw();
    hDoubleTemplatePeakFactor->Draw("same");
    leg->Draw("same");

    c1->Update();
    c1->SaveAs(Form(SaveFile + "/SignalFactorComp." + PicFormat));
    c1->Clear();

    delete leg;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Factor comp between Chi2 map and Iterative method
  if(wpsid == "all" || wpsid.Contains("factorcomp")){
    c1->cd();

    h_y_min->SetLineColor(kMagenta+2);
    h_y_min->SetMarkerColor(kMagenta+2);
    h_y_min->SetMarkerStyle(25);
    h_y_min->SetMarkerSize(1.5);
    hDoubleTemplatecorrBGFactor->SetMarkerStyle(25);
    h_y_min->GetYaxis()->SetRangeUser(-0.2, 5.);

    TLegend* leg = new TLegend(0.15,0.7,0.6,0.9);
    SetLegendSettigns(leg);
    leg->AddEntry(hDoubleTemplatecorrBGFactor, "corr. BG. scaling factor" + doubletempstring, "lp");
    leg->AddEntry(h_y_min, "corr. BG. scaling factor" +  doubletempstring + " with chi2map" , "lp");


    h_y_min->Draw();
    hDoubleTemplatecorrBGFactor->Draw("same");
    leg->Draw("same");

    c1->Update();
    c1->SaveAs(Form(SaveFile + "/BGFactorComp." + PicFormat));
    c1->Clear();

    delete leg;

    hDoubleTemplatecorrBGFactor->Add(hDoubleTemplatePeakFactor,-1);
    double mean = 0;
    for (int i = 1; i < 18; i++) {
      mean += hDoubleTemplatecorrBGFactor->GetBinContent(i);
    }
    mean /= 17.;

    double variance = 0;
    for (int i = 1; i < 18; i++) {
      variance += pow(hDoubleTemplatecorrBGFactor->GetBinContent(i)-mean, 2.);
    }
    variance = variance/17.;
    variance = sqrt(variance);

    std::cout << "IterMethod:" << '\n';
    std::cout << "mean = " << mean << '\n';
    std::cout << "variance = " << variance << '\n';

    h_y_min->Add(h_x_min, -1);
    mean = 0;
    variance = 0;
    for (int i = 1; i < 17; i++) {
      mean += h_y_min->GetBinContent(i);
    }
    mean /= 16.;

    for (int i = 1; i < 17; i++) {
      variance += pow(h_y_min->GetBinContent(i)-mean, 2.);
    }
    variance = variance/16.;
    variance = sqrt(variance);

    std::cout << "Chi2MapMethod:" << '\n';
    std::cout << "mean = " << mean << '\n';
    std::cout << "variance = " << variance << '\n';

  }




  delete pad1InvMass;
  delete pad2InvMass;
  delete canInvMass;
  // delete fit_eq_double_temp;
  // delete fit_eq_1;
  delete hChi2_DT_Iter;
  delete hChi2_pol1;
  delete hPeakRatio;
  delete hPeakComp;
  delete hRatioDoubleTemp;
  delete hRatioPol1;
  delete hData;
  delete hData_Pol1Error;
  delete hData_DTError;
  delete hDTPeak;
  delete hDTBG;
  delete hPol1Peak;
  delete hDoubleTemplatePeakFactor;
  delete hDoubleTemplatecorrBGFactor;
  delete hPol1PeakFactor;
  delete fpol1;
  // delete fitrange2;
  delete line_0;
  delete line_p1;
  delete line_m1;
  delete line_p3;
  delete line_m3;
  delete line_one;
  delete hDoubleTemp;
  delete hPol1;
  delete hChi2_dt_iter;
  delete hChi2_dt_iter_test;
  delete hChi2_dt_iter_selfcalc;
  delete hChi2_pol1_iter;
  delete hChi2_2D;
  delete hChi2_2D_sigma;
  delete hChi2_DT_Chi2map;
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
  delete hYield_dt_uncorr;
  delete hYield_pol1_uncorr;
  delete hYield_dt_chi2map_uncorr;
  delete hYield_dt_corrected;
  delete hYield_dt_chi2map_corrected;
  delete hYield_pol1_corrected;
  delete hCorrectedYieldTrueEff;
  delete histoChi2_0;


  delete c1;
  delete c2;
  IterTemp->Close();
}
