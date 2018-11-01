#include "CommonHeader.h"

void CorrBackFitPlot(TString PicFormat = "eps"){

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

  TFile* OutputFile3to8   = SafelyOpenRootfile("OutputFileBetterBkg3to8.root");
  TFile* OutputFilePulse  = SafelyOpenRootfile("Pulse.root");

  TH1D* hCorrBkgScaling   = (TH1D*) OutputFile3to8->Get("h_y_min");
  TH1D* hConvInter        = (TH1D*) OutputFilePulse->Get("hConvInter");
  TF1*  fPulse            = (TF1*) OutputFilePulse->Get("fPulse");

  hCorrBkgScaling->SetMarkerSize(3);
  hCorrBkgScaling->SetMarkerStyle(20);
  hCorrBkgScaling->SetLineWidth(5);
  hCorrBkgScaling->SetMarkerColor(kBlue);
  hCorrBkgScaling->SetLineColor(kBlue);


  SetHistoStandardSettings(hConvInter);
  hConvInter->SetFillColor(2);

  TLegend* leg = new TLegend(0.15,0.7,0.6,0.9);
  SetLegendSettigns(leg, 40);
  leg->AddEntry(hCorrBkgScaling, "korr. Untergrund. Skalierungsfaktor" , "lp");
  leg->AddEntry(fPulse,  "Pulsefunktion", "l");
  leg->AddEntry(hConvInter, "Konfidenzintervall", "e");

  c1->Update();
  hCorrBkgScaling->Draw("AXIS");
  // fPulse->Draw("SAME");
  hConvInter->Draw("SAME E3");
  hCorrBkgScaling->Draw("SAME PE");
  c1->Update();
  c1->SaveAs(Form("BetterBkg3to8Pulse/BkgConfidenceIntervall2." + PicFormat));
  c1->Clear();

  delete leg;
}
