#include "CommonHeader.h"

// wpsid = which picture should I draw
void IterTempPlot(int binnumber = 3, TString wpsid = "all"){

  TString str;
  const Int_t nbins = 45;
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

  TH1D* hChi2ndf_dt = NULL;
  TH1D* hChi2_pol1 = NULL;
  TH1D* hPeakRatio = NULL;
  TH1D* hPeakComp = NULL;
  TH1D* hRatioDoubleTemp = NULL;
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
  TH1D* hDoubleTemp = NULL;
  TH1D* hPol1 = NULL;
  TH1D* hChi2_dt_iter = NULL;
  TH1D* hChi2_pol1_iter = NULL;
  TH2D* hChi2_2D = NULL;
  TH2D* hChi2_2D_sigma = NULL;
  TH1D* hChi2_dt = NULL;
  TH1D* hSignalAreaScaling = NULL;
  TH1D* hCorrbackAreaScaling = NULL;
  TH1D* h_x_min = NULL;
  TH1D* h_y_min = NULL;
  TH1D* hErrXlow = NULL;
  TH1D* hErrXhigh = NULL;
  TH1D* hErrYlow = NULL;
  TH1D* hErrYhigh = NULL;
  TH1D* hSignal = NULL;
  TH1D* hCorrBack = NULL;
  TF1* fpol1 = NULL;
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

  hChi2ndf_dt = (TH1D*) IterTemp->Get(Form("hchi2_dt"));
  hChi2_pol1 = (TH1D*) IterTemp->Get(Form("hchi2_pol1"));
  hPeakRatio = (TH1D*) IterTemp->Get(Form("hpeakratio"));
  hPeakComp = (TH1D*) IterTemp->Get(Form("hpeakcomp"));
  hDoubleTemplatePeakFactor = (TH1D*) IterTemp->Get(Form("DoubleTemplatePeakFactor"));
  hDoubleTemplatecorrBGFactor = (TH1D*) IterTemp->Get(Form("DoubleTemplatecorrBGFactor"));
  hPol1PeakFactor = (TH1D*) IterTemp->Get(Form("Pol1PeakFactor"));
  hYield_dt_uncorr = (TH1D*) IterTemp->Get(Form("hYield_dt_uncorr"));
  hYield_pol1_uncorr = (TH1D*) IterTemp->Get(Form("hYield_pol1_uncorr"));
  hChi2_dt = (TH1D*)IterTemp->Get("hChi2_dt");
  hSignalAreaScaling = (TH1D*)IterTemp->Get("hSignalAreaScaling");
  hCorrbackAreaScaling = (TH1D*)IterTemp->Get("hCorrbackAreaScaling");
  h_x_min = (TH1D*)IterTemp->Get("h_x_min");
  h_y_min = (TH1D*)IterTemp->Get("h_y_min");
  hErrXlow = (TH1D*)IterTemp->Get("hErrXlow");
  hErrXhigh = (TH1D*)IterTemp->Get("hErrXhigh");
  hErrYlow = (TH1D*)IterTemp->Get("hErrYlow");
  hErrYhigh = (TH1D*)IterTemp->Get("hErrYhigh");


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
  line_one = new TLine(0.0, 1.0, 21.0, 1.0);
  line_one->SetLineWidth(3);
  line_one->SetLineStyle(1);
  line_one->SetLineColor(kBlack);




  //////////////////////////////////////////////////////////////////////////////
  // going over all pt bins despite first one, which is some framework bs.
  for (int k = 1; k < numberbins; k++) {

    if(binnumber <=  0 || binnumber > numberbins){
      hData = (TH1D*) IterTemp->Get(Form("data_bin%02i",k));
      str = hData->GetTitle();
      hData->SetTitle("");
      hData_Pol1Error = (TH1D*) IterTemp->Get(Form("data_addedErrosPol1_bin%02i",k));
      hData_DTError = (TH1D*) IterTemp->Get(Form("data_addedErrosDT_bin%02i",k));
      hPol1Peak = (TH1D*) IterTemp->Get(Form("mc_peak_pol1_bin%02i",k));
      hDTPeak = (TH1D*) IterTemp->Get(Form("mc_full_DT_bin%02i",k));
      hDTBG = (TH1D*) IterTemp->Get(Form("korrBG_bin%02i",k));
      fpol1 = (TF1*) IterTemp->Get(Form("fpol1_bin%02i",k));
      hRatioDoubleTemp = (TH1D*) IterTemp->Get(Form("hRatioDoubleTemp_bin%02i",k));
      hRatioPol1 = (TH1D*) IterTemp->Get(Form("hRatioPol1_bin%02i",k));
      // mc_full_clone1 = (TH1D*) IterTemp->Get(Form("mc_full_clone_beforeIterFit_bin%02d",k));
      // korrBG_clone1 = (TH1D*) IterTemp->Get(Form("korrBG_clone_beforeIterFit_bin%02d",k));
      hDoubleTemp = (TH1D*) IterTemp->Get(Form("hDoubleTemp_bin%02d",k));
      hPol1 = (TH1D*) IterTemp->Get(Form("hPol1_bin%02d",k));
      hChi2_dt_iter = (TH1D*) IterTemp->Get(Form("hChi2_dt_iter_bin%02d",k));
      hChi2_pol1_iter = (TH1D*) IterTemp->Get(Form("hChi2_pol1_iter_bin%02d",k));
      hChi2_2D = (TH2D*) IterTemp->Get(Form("hChi2_2Dbin%02d",k));
      hChi2_2D_sigma = (TH2D*) IterTemp->Get(Form("hChi2_2D_sigma_bin%02d",k));
      hSignal = (TH1D*) IterTemp->Get(Form("hSignal_bin%02d",k));
      hCorrBack = (TH1D*) IterTemp->Get(Form("hCorrBack_bin%02d",k));
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
        canInvMass->SaveAs(Form("MCTemplatesAnData/DataFitWithMCCompIter%02i.png",k));
        canInvMass->Clear("D");

        delete leg;

      }
      if(wpsid == "all" || wpsid.Contains("bgcomp")){
        //////////////////////////////////////////////////////////////////////
        // Drawing both corr. BG versions to Data with normal errors
        c1->cd();
        TLegend* leg = new TLegend(0.5,0.5,0.9,0.63);
        SetLegendSettigns(leg, 0.03);
        leg->AddEntry(fpol1, "1^{st} ord. pol.", "l");
        leg->AddEntry(hDTBG, "scaled corr. back. temp.", "p");
        hData->GetXaxis()->SetTitleSize(0.03);
        hData->GetYaxis()->SetTitleSize(0.03);
        hData->GetXaxis()->SetLabelSize(0.03);
        hData->GetYaxis()->SetLabelSize(0.03);

        hDTBG->Draw("");
        fpol1->Draw("same");
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
        c1->SaveAs(Form("MCTemplatesAnData/CorrBGComp%02i.png",k));
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
        leg->AddEntry(hChi2_dt_iter, doubletempstring, "p");
        leg->AddEntry(hChi2_pol1_iter, pol1string, "p");
        hChi2_dt_iter->GetYaxis()->SetRangeUser(0.,1.5*hChi2_dt_iter->GetMaximum());

        hChi2_dt_iter->Draw("HIST");
        hChi2_pol1_iter->Draw("SAME HIST");
        c1->Update();
        leg->Draw("same");
        DrawLabelALICE(0.5, 0.9, 0.02, 0.03, str);

        c1->Update();
        c1->SaveAs(Form("MCTemplatesAnData/Monitoring%02i.png",k));
        c1->Clear();

        delete leg;
      }

      if(wpsid == "all" || wpsid.Contains("chi2map")){
        ////////////////////////////////////////////////////////////////////////
        // Drawing Chi2 maps
        c1->cd();

        hChi2_2D->Draw("colz");
        hChi2_2D_sigma->SetLineColor(kWhite);
        hChi2_2D_sigma->SetLineWidth(2);
        hChi2_2D_sigma->SetContour(2, somelist);
        hChi2_2D_sigma->Draw("same cont3");
        c1->Update();

        c1->Update();
        c1->SaveAs(Form("MCTemplatesAnData/Chi2Map%02i.png",k));
        c1->Clear();
      }

        if(wpsid == "all" || wpsid.Contains("chi2map")){
          ////////////////////////////////////////////////////////////////////////
          // Drawing the Plot coming from the chi2map data
          c1->cd();
          hSignal->Scale(hSignalAreaScaling->GetBinContent(k)*h_x_min->GetBinContent(k));
          hCorrBack->Scale(hCorrbackAreaScaling->GetBinContent(k)*h_y_min->GetBinContent(k));
          hSignal->Add(hCorrBack);

          hData->Draw("");
          hSignal->Draw("same");

          c1->Update();

          c1->Update();
          c1->SaveAs(Form("MCTemplatesAnData/ParamResultWithChi2Map%02i.png",k));
          c1->Clear();
        }
    }

      else{
        if(binnumber == k){
          hData = (TH1D*) IterTemp->Get(Form("data_bin%02i",k));
          str = hData->GetTitle();
          hData->SetTitle("");
          hData_Pol1Error = (TH1D*) IterTemp->Get(Form("data_addedErrosPol1_bin%02i",k));
          hData_DTError = (TH1D*) IterTemp->Get(Form("data_addedErrosDT_bin%02i",k));
          hPol1Peak = (TH1D*) IterTemp->Get(Form("mc_peak_pol1_bin%02i",k));
          hDTPeak = (TH1D*) IterTemp->Get(Form("mc_full_DT_bin%02i",k));
          hDTBG = (TH1D*) IterTemp->Get(Form("korrBG_bin%02i",k));
          fpol1 = (TF1*) IterTemp->Get(Form("fpol1_bin%02i",k));
          hRatioDoubleTemp = (TH1D*) IterTemp->Get(Form("hRatioDoubleTemp_bin%02i",k));
          hRatioPol1 = (TH1D*) IterTemp->Get(Form("hRatioPol1_bin%02i",k));
          // mc_full_clone1 = (TH1D*) IterTemp->Get(Form("mc_full_clone_beforeIterFit_bin%02d",k));
          // korrBG_clone1 = (TH1D*) IterTemp->Get(Form("korrBG_clone_beforeIterFit_bin%02d",k));
          hDoubleTemp = (TH1D*) IterTemp->Get(Form("hDoubleTemp_bin%02d",k));
          hPol1 = (TH1D*) IterTemp->Get(Form("hPol1_bin%02d",k));
          hChi2_dt_iter = (TH1D*) IterTemp->Get(Form("hChi2_dt_iter_bin%02d",k));
          hChi2_pol1_iter = (TH1D*) IterTemp->Get(Form("hChi2_pol1_iter_bin%02d",k));
          hChi2_2D = (TH2D*) IterTemp->Get(Form("hChi2_2Dbin%02d",k));
          hChi2_2D_sigma = (TH2D*) IterTemp->Get(Form("hChi2_2D_sigma_bin%02d",k));
          hSignal = (TH1D*) IterTemp->Get(Form("hSignal_bin%02d",k));
          hCorrBack = (TH1D*) IterTemp->Get(Form("hCorrBack_bin%02d",k));
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
            canInvMass->SaveAs(Form("MCTemplatesAnData/DataFitWithMCCompIter%02i.png",k));
            canInvMass->Clear("D");

            delete leg;

          }
          if(wpsid == "all" || wpsid.Contains("bgcomp")){
            //////////////////////////////////////////////////////////////////////
            // Drawing both corr. BG versions to Data with normal errors
            c1->cd();
            TLegend* leg = new TLegend(0.5,0.5,0.9,0.64);
            SetLegendSettigns(leg, 0.03);
            leg->AddEntry(hPol1, "1^{st} ord. pol.", "l");
            leg->AddEntry(hDoubleTemp, "scaled corr. back. temp.", "p");
            hData->GetXaxis()->SetTitleSize(0.03);
            hData->GetYaxis()->SetTitleSize(0.03);
            hData->GetXaxis()->SetLabelSize(0.03);
            hData->GetYaxis()->SetLabelSize(0.03);

            hDTBG->Draw("");
            fpol1->Draw("same");
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
            c1->SaveAs(Form("MCTemplatesAnData/CorrBGComp%02i.png",k));
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
            leg->AddEntry(hChi2_dt_iter, doubletempstring, "p");
            leg->AddEntry(hChi2_pol1_iter, pol1string, "p");
            hChi2_dt_iter->GetYaxis()->SetRangeUser(0.,1.5*hChi2_dt_iter->GetMaximum());

            hChi2_dt_iter->Draw("HIST");
            hChi2_pol1_iter->Draw("SAME HIST");
            c1->Update();
            leg->Draw("same");
            DrawLabelALICE(0.5, 0.9, 0.02, 0.03, str);

            c1->Update();
            c1->SaveAs(Form("MCTemplatesAnData/Monitoring%02i.png",k));
            c1->Clear();

            delete leg;
          }

          if(wpsid == "all" || wpsid.Contains("chi2map")){
            ////////////////////////////////////////////////////////////////////////
            // Drawing chi2maps
            c1->cd();

            hChi2_2D->Draw("colz");
            hChi2_2D_sigma->SetLineColor(kWhite);
            hChi2_2D_sigma->SetLineWidth(2);
            hChi2_2D_sigma->SetContour(2, somelist);
            hChi2_2D_sigma->Draw("same cont3");
            c1->Update();

            c1->Update();
            c1->SaveAs(Form("MCTemplatesAnData/Chi2Map%02i.png",k));
            c1->Clear();
          }
        }
      }
  }

  // Drawing of Chi^2 comparison bwtween the two fits
  if(wpsid == "all" || wpsid.Contains("chi2")){

    TLegend* leg = new TLegend(0.6,0.75,0.9,0.9);
    SetLegendSettigns(leg);
    leg->AddEntry(hChi2ndf_dt, doubletempstring, "l");
    leg->AddEntry(hChi2_pol1, pol1string, "l");
    Double_t chi2_dt_mean = 0;
    Double_t chi2_pol1_mean = 0;
    for (int i = 0; i < numberbins; i++) {
      chi2_dt_mean += hChi2ndf_dt->GetBinContent(i);
      chi2_pol1_mean += hChi2_pol1->GetBinContent(i);
    }
    chi2_dt_mean /= (Double_t)numberbins;
    chi2_pol1_mean /= (Double_t)numberbins;

    c1->cd();
    c1->Clear();
    hChi2ndf_dt->Draw("");
    hChi2_pol1->Draw("same");
    line_one->Draw("same");
    leg->Draw("same");
    DrawLabelALICE(0.2, 0.9, 0.018, 0.03);

    TLatex* tex = new TLatex();
    SetLatexSettings(tex);
    tex->DrawLatexNDC(0.6,0.7,Form("#LT#chi^{2}_{dt}/ndf#GT = %1.2lf",chi2_dt_mean));
    tex->DrawLatexNDC(0.6,0.65,Form("#LT#chi^{2}_{pol1}/ndf#GT = %1.2lf",chi2_pol1_mean));

    c1->Update();
    c1->SaveAs(Form("MCTemplatesAnData/Chi2.png"));
    c1->Clear();
    delete leg;
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
    c1->SaveAs(Form("MCTemplatesAnData/corr_BG_to_peak.png"));
    c1->Clear();
  }

  // draing of a_pol1/a_double temp
  if(wpsid == "all" || wpsid.Contains("peakcomp")){
    c1->cd();
    c1->Clear();

    hPeakComp->Draw("");
    line_one->Draw("same");
    DrawLabelALICE(0.13, 0.9, 0.018, 0.03);
    c1->Update();
    c1->SaveAs(Form("MCTemplatesAnData/Peakcomp.png"));
    c1->Clear();
  }
  // drawing uncorrected yields
  if(wpsid == "all" || wpsid.Contains("uncorryield")){

    TLegend* leg = new TLegend(0.6,0.75,0.9,0.9);
    SetLegendSettigns(leg);
    leg->AddEntry(hYield_dt_uncorr, doubletempstring, "lp");
    leg->AddEntry(hYield_pol1_uncorr, pol1string, "lp");

    c1->cd();
    c1->Clear();
    c1->SetLogy(1);
    hYield_pol1_uncorr->GetYaxis()->SetRangeUser(1.e0,1.e6);
    hYield_pol1_uncorr->Draw("lp");
    hYield_dt_uncorr->Draw("samelp");
    leg->Draw("same");

    DrawLabelALICE(0.3, 0.9, 0.018, 0.03);
    c1->Update();
    c1->SaveAs(Form("MCTemplatesAnData/UncorrYields.png"));
    c1->Clear();
    c1->SetLogy(0);
    delete leg;
  }


  delete pad1InvMass;
  delete pad2InvMass;
  delete canInvMass;
  // delete fit_eq_double_temp;
  // delete fit_eq_1;
  delete hChi2ndf_dt;
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
  delete hChi2_pol1_iter;
  delete hChi2_2D;
  delete hChi2_2D_sigma;
  delete hChi2_dt;
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
  IterTemp->Close();
}
