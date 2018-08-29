#include "CommonHeader.h"

// wpsid = which picture should I draw
void dpg(int binnumber = 3, TString wpsid = "all", TString PicFormat = "png", TString SaveFile = "DPG"){

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
  c2->SetBottomMargin(0.1);
  c2->SetRightMargin(0.15);
  c2->SetLeftMargin(0.09);
  c2->SetTicky();
  c2->SetTickx();
  c2->SetLogz(1);

  TCanvas *canInvMass = new TCanvas("canInvMass","",1200,1600);
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

  hPol1PeakFactor = (TH1D*) IterTemp->Get(Form("Pol1PeakFactor"));
  hSignalAreaScaling = (TH1D*)IterTemp->Get("hSignalAreaScaling");
  hCorrbackAreaScaling = (TH1D*)IterTemp->Get("hCorrbackAreaScaling");
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
  line_one->SetLineWidth(3);
  line_one->SetLineStyle(1);
  line_one->SetLineColor(kBlack);

  //////////////////////////////////////////////////////////////////////////////
  // going over all pt bins despite first one, which is some framework bs.
  for (int k = 1; k < numberbins; k++) {

    if(binnumber <=  0 || binnumber > numberbins){
      hData                    = (TH1D*) IterTemp->Get(Form("data_bin%02i",k));
      str                      = hData->GetTitle();
      hData->SetTitle("");
      hPol1Peak                = (TH1D*) IterTemp->Get(Form("mc_peak_pol1_bin%02i",k));
      hDTPeak                  = (TH1D*) IterTemp->Get(Form("mc_full_DT_bin%02i",k));
      hDTBG                    = (TH1D*) IterTemp->Get(Form("korrBG_bin%02i",k));
      fpol1                    =  (TF1*) IterTemp->Get(Form("fpol1_bin%02i",k));
      hRatioPol1               = (TH1D*) IterTemp->Get(Form("hRatioPol1_bin%02i",k));
      hPol1                    = (TH1D*) IterTemp->Get(Form("hPol1_bin%02d",k));
      hSignal                  = (TH1D*) IterTemp->Get(Form("hSignal_bin%02d",k));
      hCorrBack                = (TH1D*) IterTemp->Get(Form("hCorrBack_bin%02d",k));

    if(wpsid == "all" || wpsid.Contains("Chi2mapParam")){
      ////////////////////////////////////////////////////////////////////////
      // Comparison between the two double temp methods with data
      c1->cd();

      TH1D* hChi2Map_Param = NULL;
      TH1D* hSignal_Clone = NULL;
      TH1D* hCorrBack_Clone = NULL;
      hChi2Map_Param = (TH1D*) hSignal->Clone("hSignal_Clone");
      hSignal_Clone = (TH1D*) hSignal->Clone("hSignal_Clone");
      hCorrBack_Clone = (TH1D*) hCorrBack->Clone("hCorrBack_Clone");
      hChi2Map_Param->Scale(hSignalAreaScaling->GetBinContent(k)*h_x_min->GetBinContent(k+1));
      hSignal_Clone->Scale(hSignalAreaScaling->GetBinContent(k)*h_x_min->GetBinContent(k+1));
      hCorrBack_Clone->Scale(hCorrbackAreaScaling->GetBinContent(k)*h_y_min->GetBinContent(k+1));
      hChi2Map_Param->Add(hCorrBack_Clone);

      hSignal_Clone->SetMarkerStyle(20);
      hSignal_Clone->SetMarkerSize(1.5);


      hData->SetMarkerColor(kBlack);
      hData->SetLineColor(kBlack);
      hData->SetMarkerStyle(21);
      hChi2Map_Param->SetMarkerColor(kRed);
      hChi2Map_Param->SetLineColor(kRed);
      hChi2Map_Param->SetMarkerStyle(1);
      hChi2Map_Param->SetMarkerSize(1);
      hCorrBack_Clone->SetMarkerColor(kMagenta+2);
      hCorrBack_Clone->SetLineColor(kMagenta+2);
      hCorrBack_Clone->SetMarkerStyle(1);
      hCorrBack_Clone->SetMarkerSize(1);

      TLegend* leg = new TLegend(0.6,0.9,0.9,0.65);
      SetLegendSettigns(leg, 0.035);
      leg->AddEntry(hData, strData, "p");
      leg->AddEntry(hChi2Map_Param, "signal + back. temp.", "l");
      leg->AddEntry((TObject*) 0x0, "parametrization", "");
      leg->AddEntry(hCorrBack_Clone, "scaled back. template", "l");

      SetHistoStandardSettings(hData, 0., 0., 0.035);

      hData->GetXaxis()->SetRangeUser(0.0, 0.3);
      hData->SetTitle("");
      hData->GetYaxis()->SetTitleOffset(1.1);
      hData->GetXaxis()->SetTitleOffset(1.);
      hData->Draw("AXIS");
      hChi2Map_Param->Draw("same HIST");
      hChi2Map_Param->Draw("SAME P");
      hCorrBack_Clone->Draw("same HIST");
      hCorrBack_Clone->Draw("SAME P");
      hData->Draw("SAME P");
      c1->Update();
      line_y = gPad->GetUymax()*0.995;
      fitrange2 = new TLine(lowerparamrange, line_y, upperparamrange, line_y);
      fitrange2->SetLineColor(kAzure+10);
      fitrange2->SetLineWidth(7);
      fitrange2->Draw("same");
      leg->AddEntry(fitrange2, "parametrization range", "l");
      leg->Draw("same");
      DrawLabelALICE(0.12, 0.9, 0.03, 0.035, str);
      c1->Update();
      c1->SaveAs(Form(SaveFile + "/Chi2mapParam%02i." + PicFormat,k));
      c1->Clear("D");


      delete leg;
    }

    if(wpsid == "all" || wpsid.Contains("Pol1Param")){
      ////////////////////////////////////////////////////////////////////////
      // Comparison between the two double temp methods with data
      c1->cd();

      TLegend* leg = new TLegend(0.6,0.9,0.9,0.65);
      SetLegendSettigns(leg, 0.035);
      leg->AddEntry(hData, strData, "p");
      leg->AddEntry(hPol1, "signal temp. + 1^{st} ord. pol.", "l");
      leg->AddEntry((TObject*) 0x0, "parametrization", "");
      leg->AddEntry(fpol1, " 1^{st} order polynomial", "l");

      hData->GetXaxis()->SetRangeUser(0.0, 0.3);
      hData->SetMarkerColor(kBlack);
      hData->SetLineColor(kBlack);
      hData->SetMarkerStyle(21);
      hPol1->SetMarkerColor(kRed);
      hPol1->SetLineColor(kRed);
      hPol1->SetMarkerStyle(1);
      hPol1->SetMarkerSize(1);
      fpol1->SetLineColor(kMagenta+2);


      SetHistoStandardSettings(hData, 0., 0., 0.035);
      hData->GetXaxis()->SetTitleOffset(1.);
      hData->SetTitle("");
      hData->GetYaxis()->SetTitleOffset(1.1);
      hData->Draw("AXIS");
      hPol1->Draw("SAME HIST");
      hPol1->Draw("SAME P");
      hData->Draw("SAME P");
      fpol1->Draw("same");
      c1->Update();
      line_y = gPad->GetUymax()*0.995;
      fitrange2 = new TLine(lowerparamrange, line_y, upperparamrange, line_y);
      fitrange2->SetLineColor(kAzure+10);
      fitrange2->SetLineWidth(7);
      fitrange2->Draw("same");
      leg->AddEntry(fitrange2, "parametrization range", "l");
      leg->Draw("same");
      DrawLabelALICE(0.12, 0.9, 0.03, 0.035, str);
      c1->Update();
      c1->SaveAs(Form(SaveFile + "/Pol1Param%02i." + PicFormat,k));
      c1->Clear("D");


      delete leg;
    }
  }
    else{
      break;
    }
  }

  delete hChi2_DT_Iter;
  delete hChi2_pol1;
  delete hPeakRatio;
  delete hPeakComp;
  delete hRatioDoubleTemp;
  delete hRatioDoubleTemp_chi2map;
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
  delete hYield_dt_uncorr;
  delete hYield_pol1_uncorr;
  delete hYield_dt_chi2map_uncorr;
  delete hDoubleTemp;
  delete hPol1;
  delete hChi2_dt_iter;
  delete hChi2_pol1_iter;
  delete hChi2_dt_iter_test;
  delete hChi2_dt_iter_selfcalc;
  delete hChi2_2D;
  delete hChi2_2D_sigma;
  delete hChi2_DT_Chi2map;
  delete histoChi2_0;
  delete hSignalAreaScaling;
  delete hCorrbackAreaScaling;
  delete hYield_dt_corrected;
  delete hYield_dt_chi2map_corrected;
  delete hYield_pol1_corrected;
  delete hCorrectedYieldTrueEff;
  delete h_x_min;
  delete h_y_min;
  delete hErrXlow;
  delete hErrXhigh;
  delete hErrYlow;
  delete hErrYhigh;
  delete hSignal;
  delete hCorrBack;
  delete fpol1;
  delete f_ChiOverNdf;
  delete fitrange2;
  delete line_0;
  delete line_p1;
  delete line_m1;
  delete line_p3;
  delete line_m3;
  delete line_one;

}
