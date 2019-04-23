#include "CommonHeader.h"

void CorrBackFitPlot(TString PicFormat = "eps"){

  TH1D* hCorrBkgNormal = NULL;
  TH1D* hCorrBkg3to8   = NULL;
  TH1D* hCorrBkgRatio  = NULL;

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


  TCanvas *canInvMass = new TCanvas("canInvMass","",1540,1417);
  TPad *pad1InvMass = new TPad("pad1InvMass","",0.0,0.33,1.0,1.0);
  pad1InvMass->SetTopMargin(0.05);
  pad1InvMass->SetLeftMargin(0.13);
  pad1InvMass->SetBottomMargin(0.0);
  pad1InvMass->SetRightMargin(0.02);
  TPad *pad2InvMass = new TPad("pad2InvMass","",0.0,0.0,1.0,0.33);
  pad2InvMass->SetTopMargin(0.0);
  pad2InvMass->SetLeftMargin(0.13);
  pad2InvMass->SetBottomMargin(0.2);
  pad2InvMass->SetRightMargin(0.02);
  pad2InvMass->SetTicky();


  TFile* OutputFile3to8   = SafelyOpenRootfile("OutputFileBetterBkg3to8.root");
  TFile* OutputFilePulse  = SafelyOpenRootfile("Pulse.root");
  TFile* OutputFileNormal = SafelyOpenRootfile("OutputFileNormal.root");

  TH1D* hCorrBkgScaling   = (TH1D*) OutputFile3to8->Get("h_y_min");
  TH1D* hConvInter        = (TH1D*) OutputFilePulse->Get("hConvInter");
  TF1*  fPulse            = (TF1*) OutputFilePulse->Get("fPulse");

  hCorrBkgScaling->SetMarkerSize(3);
  hCorrBkgScaling->SetMarkerStyle(20);
  hCorrBkgScaling->SetLineWidth(5);
  hCorrBkgScaling->SetMarkerColor(kBlue);
  hCorrBkgScaling->SetLineColor(kBlue);

  c1->cd();
  SetHistoStandardSettings(hConvInter);
  hConvInter->SetFillColor(2);

  TLegend* leg = new TLegend(0.3,0.7,0.7,0.9);
  SetLegendSettigns(leg, 40);
  hConvInter->SetMarkerStyle(1);
  hConvInter->SetMarkerSize(1);
  hCorrBkgScaling->SetYTitle("SF_{korr. Untergrund}");

  hCorrBkgScaling->GetXaxis()->SetRangeUser(1.4, 12.);
  fPulse->GetXaxis()->SetRangeUser(1.4, 12.);
  hConvInter->GetXaxis()->SetRangeUser(1.4, 12.);
  leg->AddEntry(hCorrBkgScaling, "SF_{korr. Untergrund}" , "lp");
  leg->AddEntry(fPulse,  "Konfidenzintervall", "l");

  c1->Update();
  hCorrBkgScaling->Draw("AXIS");
  hConvInter->Draw("SAME E3");
  hCorrBkgScaling->Draw("SAME PE");
  leg->Draw("SAME");
  c1->Update();
  c1->SaveAs(Form("BetterBkg3to8Pulse/BkgConfidenceIntervall2." + PicFormat));
  c1->Clear();

  for(int k = 1; k < numberbins - 1; k++){
    TString str = Form("%.1lf #leq #it{p}_{T} < %.1lf",fBinsPi013TeVEMCPt[k], fBinsPi013TeVEMCPt[k+1]);
    TString pTBig = Form("%.1lf #leq #it{p}_{T} < %.1lf",fBinsPi013TeVEMCPt[3], fBinsPi013TeVEMCPt[8]);
    hCorrBkgNormal = NULL;
    hCorrBkg3to8   = NULL;
    hCorrBkgRatio  = NULL;


    hCorrBkgNormal = (TH1D*) OutputFileNormal->Get(Form("hCorrBack_bin%02d",k));
    hCorrBkg3to8   = (TH1D*) OutputFile3to8->Get(Form("hCorrBack_bin%02d",k));

    hCorrBkgNormal->SetLineColor(kBlue);
    hCorrBkgNormal->SetMarkerColor(kBlue);
    hCorrBkgNormal->SetLineWidth(5);
    hCorrBkg3to8->SetLineColor(kRed);
    hCorrBkg3to8->SetMarkerColor(kRed);
    hCorrBkg3to8->SetLineWidth(5);

    TLegend* leg2 = new TLegend(0.25,0.05,0.3,0.2);
    SetLegendSettigns(leg2, 40);
    leg2->SetHeader("korr. Untergrund Template");
    leg2->AddEntry(hCorrBkgNormal, str , "l");
    leg2->AddEntry(hCorrBkg3to8,  pTBig, "l");

    canInvMass->cd();
    pad1InvMass->Draw();
    pad2InvMass->Draw("same");
    pad1InvMass->cd();

    pad1InvMass->Update();
    hCorrBkg3to8->Draw("AXIS");
    hCorrBkgNormal->DrawCopy("SAME HiST");
    hCorrBkg3to8->DrawCopy("SAME HIST");
    hCorrBkgNormal->DrawCopy("SAME");
    hCorrBkg3to8->DrawCopy("SAME");
    leg2->Draw("SAME");
    DrawLabelALICE(0.6, 0.9, 0.03, 40, "");

    pad1InvMass->Update();
    pad2InvMass->cd();
    hCorrBkgRatio = (TH1D*) hCorrBkg3to8->Clone("hCorrBkgRatio");
    SetHistoStandardSettings(hCorrBkgRatio);
    hCorrBkgRatio->Divide(hCorrBkg3to8, hCorrBkgNormal, 1 , 1, "B");
    hCorrBkgRatio->SetYTitle("Ratio");
    hCorrBkgRatio->Draw("AXIS");
    hCorrBkgRatio->Draw("SAME");
    hCorrBkgRatio->GetYaxis()->SetRangeUser(0., hCorrBkgRatio->GetBinContent(hCorrBkgRatio->FindBin(0.135))*3.);

    pad2InvMass->Update();
    canInvMass->Update();

    canInvMass->SaveAs(Form("Untergrund/BkgRatio%02d." + PicFormat, k));
    canInvMass->Clear("D");




    pad2InvMass->cd();


  }

  delete leg;
  OutputFile3to8->Close();
  OutputFilePulse->Close();
  OutputFileNormal->Close();
}
