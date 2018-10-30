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
  for (int templatemethod = 1; templatemethod <= 3; templatemethod++) {

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

    TH1D* hData = NULL;
    TH1D* hPol1PeakFactor = NULL;
    TH1D* hYield_dt_chi2map_uncorr = NULL;
    TH2D* hChi2_2D = NULL;
    TH2D* hChi2_2D_sigma = NULL;
    TH1D* hChi2Map_Chi2_pT = NULL;
    TH1D* histoChi2_0 = NULL;
    TH1D* hSignalAreaScaling = NULL;
    TH1D* hCorrbackAreaScaling = NULL;
    TH1D* hYield_dt_chi2map_corrected = NULL;
    TH1D* hYield_pol1_corrected = NULL;
    TH1D* hCorrectedYieldNormEff = NULL;
    TH1D* h_x_min = NULL;
    TH1D* h_y_min = NULL;
    TH1D* hErrXlow = NULL;
    TH1D* hErrXhigh = NULL;
    TH1D* hErrYlow = NULL;
    TH1D* hErrYhigh = NULL;
    TH1D* hSignal = NULL;
    TH1D* hCorrBack = NULL;
    TF1* f_ChiOverNdf = NULL;
    TLine* fitrange2 = NULL;
    TLine* line_0 = NULL;
    TLine* line_p1 = NULL;
    TLine* line_m1 = NULL;
    TLine* line_p3 = NULL;
    TLine* line_m3 = NULL;
    TLine* line_one = NULL;

    Double_t line_y = 0;


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
    else{
      std::cerr << "No Outputfile loaded!" << std::endl;
    }
    if (OutputFile->IsOpen() ) printf("OutputFile opened successfully\n");

    hYield_dt_chi2map_uncorr = (TH1D*) OutputFile->Get(Form("hYield_dt_chi2map_uncorr"));
    hChi2Map_Chi2_pT = (TH1D*)OutputFile->Get("hChi2Map_Chi2_pT");
    histoChi2_0 = (TH1D*) OutputFile->Get("histoChi2_0");
    hSignalAreaScaling = (TH1D*)OutputFile->Get("hSignalAreaScaling");
    hCorrbackAreaScaling = (TH1D*)OutputFile->Get("hCorrbackAreaScaling");
    hYield_dt_chi2map_corrected = (TH1D*)OutputFile->Get("hYield_dt_chi2map_corrected");
    hCorrectedYieldNormEff = (TH1D*) OutputFile->Get("hCorrectedYieldNormEff");
    h_x_min = (TH1D*)OutputFile->Get("h_x_min");
    h_y_min = (TH1D*)OutputFile->Get("h_y_min");


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
    for (int k = 1; k < numberbins-1; k++) {

      if(binnumber <=  0 || binnumber > numberbins){
        hData                    = (TH1D*) OutputFile->Get(Form("data_bin%02i",k));
        str                      = hData->GetTitle();
        hData->SetTitle("");
        hChi2_2D                 = (TH2D*) OutputFile->Get(Form("hChi2MapBin%02d",k));
        hChi2_2D_sigma           = (TH2D*) OutputFile->Get(Form("hChi2_2D_sigma_bin%02d",k));
        hSignal                  = (TH1D*) OutputFile->Get(Form("hSignal_bin%02d",k));
        hCorrBack                = (TH1D*) OutputFile->Get(Form("hCorrBack_bin%02d",k));

        SetHistoStandardSettings(hData,     1.2, 1., 35, kBlack);
        SetHistoStandardSettings(hSignal,   1.2, 1., 35, kTeal-7);
        SetHistoStandardSettings(hCorrBack, 1.2, 1., 35, kPink-2);


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
          if(templatemethod == 2){
            c2->SaveAs(Form("BetterBkgNN/Chi2Map%02d." + PicFormat,k));
          }
          if(templatemethod == 1){
            c2->SaveAs(Form("BetterBkg3to8/Chi2Map%02d." + PicFormat,k));
          }
          if(templatemethod == 3){
            c2->SaveAs(Form("BetterBkg3to8Pulse/Chi2Map%02d." + PicFormat,k));
          }
          c2->Clear();
        }

        if(wpsid == "all" || wpsid.Contains("chi2map")){
          ////////////////////////////////////////////////////////////////////////
          // Drawing the Plot coming from the chi2map data
          c1->cd();
          TLegend* lParamResultParts = new TLegend(0.6,0.5,0.9,0.63);
          SetLegendSettigns(lParamResultParts, 35);
          lParamResultParts->AddEntry(hData, "Data", "l");
          TH1D* hAdded = NULL;
          TH1D* hSignal_Clone = NULL;
          TH1D* hCorrBack_Clone = NULL;
          hSignal_Clone   = (TH1D*) hSignal->Clone("hSignal_Clone");
          hCorrBack_Clone = (TH1D*) hCorrBack->Clone("hCorrBack_Clone");
          hSignal_Clone->Scale(hSignalAreaScaling->GetBinContent(k)*h_x_min->GetBinContent(k+1));
          hCorrBack_Clone->Scale(hCorrbackAreaScaling->GetBinContent(k)*h_y_min->GetBinContent(k+1));
          hAdded          = (TH1D*) hSignal_Clone->Clone("hAdded");
          hAdded->Add(hCorrBack_Clone);
          // hSignal_Clone->Add(hCorrBack_Clone);
          lParamResultParts->AddEntry(hSignal_Clone, "Signal template", "l");
          lParamResultParts->AddEntry(hCorrBack_Clone, "Bkg. template", "l");

          hData->Draw("");
          hSignal_Clone->Draw("SAME")
          hSignal_Clone->Draw("SAME HIST");
          hCorrBack_Clone->Draw("SAME");
          hCorrBack_Clone->Draw("SAME HIST");
          lParamResultParts->Draw("SAME");
          DrawLabelALICE(0.6, 0.9, 0.02, 35, str);

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
          c1->Clear();
          TLegend* lParamResult = new TLegend(0.6,0.5,0.9,0.63);
          SetLegendSettigns(lParamResulParts, 35);
          lParamResulParts->AddEntry(hData, "Data", "l");
          lParamResulParts->AddEntry(hAdded, "Templateparametrization", "l");

          c1->Update();
          hData->Draw("AXIS");
          hData->Draw("SAME");
          hAdded->Draw("SAME");
          DrawLabelALICE(0.6, 0.9, 0.02, 35, str);
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

          c1->Clear();


          delete lParamResultParts;
          delete lParamResult;

        }
      }
    }

    // Drawing of Chi^2 comparison bwtween the two fits
    if(wpsid == "all" || wpsid.Contains("chi2")){
      TLegend* leg2 = new TLegend(0.6,0.75,0.9,0.9);
      SetLegendSettigns(leg2);
      leg2->AddEntry(hChi2Map_Chi2_pT, doubletempstring, "l");
      leg2->AddEntry(histoChi2_0, "parametrization with function", "l");


      hChi2Map_Chi2_pT->SetLineColor(kMagenta+2);
      histoChi2_0->SetLineWidth(3);

      c1->cd();
      c1->Clear();
      c1->SetLogx(1);
      histoChi2_0->Draw("AXIS");
      line_one->Draw("same");
      histoChi2_0->Draw("SAME HIST");
      hChi2Map_Chi2_pT->Draw("SAME HIST");
      hChi2Map_Chi2_pT->Draw("SAME P");
      leg2->Draw("same");
      DrawLabelALICE(0.2, 0.9, 0.018, 0.03);


      c1->Update();
      if(templatemethod == 0){
        c1->SaveAs(Form("BetterBkgNN/Chi2Iter_vs_Map." + PicFormat));
      }
      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkg3to8/Chi2Iter_vs_Map." + PicFormat));
      }
      if(templatemethod == 2){
        c1->SaveAs(Form("BetterBkg3to8Pulse/Chi2Iter_vs_Map." + PicFormat));
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
      hYield_dt_chi2map_uncorr->GetXaxis()->SetRangeUser(1.4, 12.);

      c1->cd();
      c1->SetLeftMargin(0.11);
      c1->Clear();
      c1->SetLogy(1);
      hYield_dt_chi2map_uncorr->Draw("AXIS");
      hYield_dt_chi2map_uncorr->Draw("samelp");
      leg->Draw("same");

      DrawLabelALICE(0.55, 0.9, 0.018, 0.03);
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
      if(templatemethod == 2){
        c1->SaveAs(Form("BetterBkgNN/BGFactorComp." + PicFormat));
      }
      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkg3to8/BGFactorComp." + PicFormat));
      }
      if(templatemethod == 3){
        c1->SaveAs(Form("BetterBkg3to8Pulse/BGFactorComp." + PicFormat));
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
      if(templatemethod == 2){
        c1->SaveAs(Form("BetterBkgNN/b_to_a_ratio." + PicFormat));
      }
      if(templatemethod == 1){
        c1->SaveAs(Form("BetterBkg3to8/b_to_a_ratio." + PicFormat));
      }
      if(templatemethod == 3){
        c1->SaveAs(Form("BetterBkg3to8Pulse/b_to_a_ratio." + PicFormat));
      }
      c1->Clear();

      h_y_min->Add(h_x_min, -1);
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
