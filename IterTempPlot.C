#include "CommonHeader.h"
//MCTemplatesAnData

// wpsid = which picture should I draw
void IterTempPlot(int binnumber = 3, TString wpsid = "all", TString PicFormat = "png", TString SaveFile = "MCTemplatesAnData"){

  /*
  loop over mario
  mario == 0 for Next Neighbor method
  mario == 1 for 3 to 8 Method
  mario == 2 for Normal method
  */
  for (int mario = 0; mario < 3; mario++) {

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

    TH1D* hChi2_pol1 = NULL;
    TH1D* hData = NULL;
    TH1D* hData_Pol1Error = NULL;
    TH1D* hData_DTError = NULL;
    TH1D* hDTPeak = NULL;
    TH1D* hPol1Peak = NULL;
    TH1D* hPol1PeakFactor = NULL;
    TH1D* hYield_pol1_uncorr = NULL;
    TH1D* hYield_dt_chi2map_uncorr = NULL;
    TH1D* hPol1 = NULL;
    TH1D* hChi2_pol1_iter = NULL;
    TH2D* hChi2_2D = NULL;
    TH2D* hChi2_2D_sigma = NULL;
    TH1D* hChi2_DT_Chi2map = NULL;
    TH1D* histoChi2_0 = NULL;
    TH1D* hSignalAreaScaling = NULL;
    TH1D* hCorrbackAreaScaling = NULL;
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


    TFile* IterTemp = NULL;
    if(mario == 0){
      IterTemp      = SafelyOpenRootfile("IterTempBetterBkgNN.root");
    }
    if(mario == 1){
      IterTemp      = SafelyOpenRootfile("IterTempBetterBkg3to8.root");
    }
    if(mario == 2){
      IterTemp      = SafelyOpenRootfile("IterTemp.root");
    }
    if (IterTemp->IsOpen() ) printf("IterTemp opened successfully\n");

    hChi2_pol1 = (TH1D*) IterTemp->Get(Form("hchi2_pol1"));
    hPol1PeakFactor = (TH1D*) IterTemp->Get(Form("Pol1PeakFactor"));
    hYield_pol1_uncorr = (TH1D*) IterTemp->Get(Form("hYield_pol1_uncorr"));
    hYield_dt_chi2map_uncorr = (TH1D*) IterTemp->Get(Form("hYield_dt_chi2map_uncorr"));
    hChi2_DT_Chi2map = (TH1D*)IterTemp->Get("hChi2_DT_Chi2map");
    histoChi2_0 = (TH1D*) IterTemp->Get("histoChi2_0");
    hSignalAreaScaling = (TH1D*)IterTemp->Get("hSignalAreaScaling");
    hCorrbackAreaScaling = (TH1D*)IterTemp->Get("hCorrbackAreaScaling");
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
        fpol1                    = (TF1*)  IterTemp->Get(Form("fpol1_bin%02i",k));
        hPol1                    = (TH1D*) IterTemp->Get(Form("hPol1_bin%02d",k));
        hChi2_pol1_iter          = (TH1D*) IterTemp->Get(Form("hChi2_pol1_iter_bin%02d",k));
        hChi2_2D                 = (TH2D*) IterTemp->Get(Form("hChi2_2Dbin%02d",k));
        hChi2_2D_sigma           = (TH2D*) IterTemp->Get(Form("hChi2_2D_sigma_bin%02d",k));
        hSignal                  = (TH1D*) IterTemp->Get(Form("hSignal_bin%02d",k));
        hCorrBack                = (TH1D*) IterTemp->Get(Form("hCorrBack_bin%02d",k));
        f_ChiOverNdf             = (TF1*)  IterTemp->Get(Form("f_ChiOverNdf%02d",k));

        hChi2_pol1_iter->SetLineWidth(3);


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
          leg->AddEntry(hCorrBack_Chi2Map, "scaled corr. back. temp with chi2map");

          hData->GetXaxis()->SetTitleSize(0.03);
          hData->GetYaxis()->SetTitleSize(0.03);
          hData->GetXaxis()->SetLabelSize(0.03);
          hData->GetYaxis()->SetLabelSize(0.03);

          hCorrBack_Chi2Map->Draw("");
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
          if(mario == 0){
            c1->SaveAs(Form("BetterBkgNN/CorrBGComp_Bin%02d." + PicFormat, k));
          }
          if(mario == 1){
            c1->SaveAs(Form("BetterBkg3to8/CorrBGComp_Bin%02d." + PicFormat, k));
          }
          if(mario == 2){
            c1->SaveAs(Form("Normal/CorrBGComp_Bin%02d." + PicFormat, k));
          }
          c1->Clear();

          delete leg;
          delete fitrange;
        }

        if(wpsid == "all" || wpsid.Contains("chi2map")){
          ////////////////////////////////////////////////////////////////////////
          // Drawing Chi2 maps
          c2->cd();

          hChi2_2D->GetXaxis()->SetTitleSize(0.035);
          hChi2_2D->GetYaxis()->SetTitleSize(0.035);
          hChi2_2D->GetXaxis()->SetLabelSize(0.035);
          hChi2_2D->GetYaxis()->SetLabelSize(0.035);

          hChi2_2D->Draw("colz");
          hChi2_2D_sigma->SetLineColor(kWhite);
          hChi2_2D_sigma->SetLineWidth(2);
          hChi2_2D_sigma->SetContour(2, somelist);
          hChi2_2D_sigma->Draw("same cont3");
          c2->Update();

          c2->Update();
          if(mario == 0){
            c2->SaveAs(Form("BetterBkgNN/Chi2Map%02d." + PicFormat,k));
          }
          if(mario == 1){
            c2->SaveAs(Form("BetterBkg3to8/Chi2Map%02d." + PicFormat,k));
          }
          if(mario == 2){
            c2->SaveAs(Form("Normal/Chi2Map%02d." + PicFormat,k));
          }
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
            if(mario == 0){
              c1->SaveAs(Form("BetterBkgNN/ParamResultWithChi2Map_Bin%02d." + PicFormat,k));
            }
            if(mario == 1){
              c1->SaveAs(Form("BetterBkg3to8/ParamResultWithChi2Map_Bin%02d." + PicFormat,k));
            }
            if(mario == 2){
              c1->SaveAs(Form("Normal/ParamResultWithChi2Map_Bin%02d." + PicFormat,k));
            }
            c1->Clear();

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
            fpol1                    = (TF1*) IterTemp->Get(Form("fpol1_bin%02i",k));
            hPol1                    = (TH1D*) IterTemp->Get(Form("hPol1_bin%02d",k));
            hChi2_pol1_iter          = (TH1D*) IterTemp->Get(Form("hChi2_pol1_iter_bin%02d",k));
            hChi2_2D                 = (TH2D*) IterTemp->Get(Form("hChi2_2Dbin%02d",k));
            hChi2_2D_sigma           = (TH2D*) IterTemp->Get(Form("hChi2_2D_sigma_bin%02d",k));
            hSignal                  = (TH1D*) IterTemp->Get(Form("hSignal_bin%02d",k));
            hCorrBack                = (TH1D*) IterTemp->Get(Form("hCorrBack_bin%02d",k));
            f_ChiOverNdf             = (TF1*) IterTemp->Get(Form("f_ChiOverNdf%02d",k));

            hChi2_pol1_iter->SetLineWidth(3);

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
              leg->AddEntry(hCorrBack_Chi2Map, "scaled corr. back. temp with chi2map");

              hData->GetXaxis()->SetTitleSize(0.03);
              hData->GetYaxis()->SetTitleSize(0.03);
              hData->GetXaxis()->SetLabelSize(0.03);
              hData->GetYaxis()->SetLabelSize(0.03);

              hCorrBack_Chi2Map->Draw("");
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
              if(mario == 0){
                c1->SaveAs(Form("BetterBkgNN/CorrBGComp%02i." + PicFormat,k));
              }
              if(mario == 1){
                c1->SaveAs(Form("BetterBkg3to8/CorrBGComp%02i." + PicFormat,k));
              }
              if(mario == 2){
                c1->SaveAs(Form("Normal/CorrBGComp%02i." + PicFormat,k));
              }
              c1->Clear();

              delete leg;
              delete fitrange;
            }
          }
        }
    }

    // Drawing of Chi^2 comparison bwtween the two fits
    if(wpsid == "all" || wpsid.Contains("chi2")){

      c1->Update();
      if(mario == 0){
        c1->SaveAs(Form("BetterBkgNN/Chi2." + PicFormat));
      }
      if(mario == 1){
        c1->SaveAs(Form("BetterBkg3to8/Chi2." + PicFormat));
      }
      if(mario == 2){
        c1->SaveAs(Form("Normal/Chi2." + PicFormat));
      }
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
      if(mario == 0){
        c1->SaveAs(Form("BetterBkgNN/Chi2Iter_vs_Map." + PicFormat));
      }
      if(mario == 1){
        c1->SaveAs(Form("BetterBkg3to8/Chi2Iter_vs_Map." + PicFormat));
      }
      if(mario == 2){
        c1->SaveAs(Form("Normal/Chi2Iter_vs_Map." + PicFormat));
      }
      c1->Clear();
      c1->SetLogx(0);

      delete leg2;
    }

    // drawing uncorrected yields
    if(wpsid == "all" || wpsid.Contains("uncorryield")){

      TLegend* leg = new TLegend(0.2,0.2,0.4,0.4);
      SetLegendSettigns(leg);
      leg->AddEntry(hYield_dt_chi2map_uncorr, doubletempstring + " with chi2map" , "lp");
      leg->AddEntry(hYield_pol1_uncorr, pol1string, "lp");

      c1->cd();
      c1->SetLeftMargin(0.11);
      c1->Clear();
      c1->SetLogy(1);
      hYield_pol1_uncorr->GetYaxis()->SetRangeUser(1.e-8,1.e-2);
      hYield_pol1_uncorr->GetYaxis()->SetTitleOffset(1.6);
      hYield_pol1_uncorr->Draw("lp");
      hYield_dt_chi2map_uncorr->Draw("samelp");
      leg->Draw("same");

      DrawLabelALICE(0.55, 0.9, 0.018, 0.03);
      c1->Update();
      if(mario == 0){
        c1->SaveAs(Form("BetterBkgNN/UncorrYields." + PicFormat));
      }
      if(mario == 1){
        c1->SaveAs(Form("BetterBkg3to8/UncorrYields." + PicFormat));
      }
      if(mario == 2){
        c1->SaveAs(Form("Normal/UncorrYields." + PicFormat));
      }
      c1->Clear();
      c1->SetLogy(0);
      c1->SetLeftMargin(0.09);
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
      h_y_min->GetYaxis()->SetRangeUser(-0.2, 5.);

      TLegend* leg = new TLegend(0.15,0.7,0.6,0.9);
      SetLegendSettigns(leg);
      leg->AddEntry(h_y_min, "corr. BG. scaling factor" +  doubletempstring + " with chi2map" , "lp");


      h_y_min->Draw();
      leg->Draw("same");

      c1->Update();
      if(mario == 0){
        c1->SaveAs(Form("BetterBkgNN/BGFactorComp." + PicFormat));
      }
      if(mario == 1){
        c1->SaveAs(Form("BetterBkg3to8/BGFactorComp." + PicFormat));
      }
      if(mario == 2){
        c1->SaveAs(Form("Normal/BGFactorComp." + PicFormat));
      }
      c1->Clear();

      delete leg;

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
      if(mario == 0){
        c1->SaveAs(Form("BetterBkgNN/b_to_a_ratio." + PicFormat));
      }
      if(mario == 1){
        c1->SaveAs(Form("BetterBkg3to8/b_to_a_ratio." + PicFormat));
      }
      if(mario == 2){
        c1->SaveAs(Form("Normal/b_to_a_ratio." + PicFormat));
      }
      c1->Clear();

      h_y_min->Add(h_x_min, -1);
      Double_t mean = 0;
      Double_t variance = 0;
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
    delete hChi2_pol1;
    delete hData;
    delete hData_Pol1Error;
    delete hData_DTError;
    delete hDTPeak;
    delete hPol1Peak;
    delete hPol1PeakFactor;
    delete fpol1;
    delete line_0;
    delete line_p1;
    delete line_m1;
    delete line_p3;
    delete line_m3;
    delete line_one;
    delete hPol1;
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
    delete hYield_pol1_uncorr;
    delete hYield_dt_chi2map_uncorr;
    delete hYield_dt_chi2map_corrected;
    delete hYield_pol1_corrected;
    delete hCorrectedYieldTrueEff;
    delete histoChi2_0;

    delete c1;
    delete c2;
    IterTemp->Close();
  }
}
