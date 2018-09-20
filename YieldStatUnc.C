#include "CommonHeader.h"

void YieldStatUnc(TString PicFormat = "png"){
  TH1D* hCorrYield_RelativStaterror      = NULL;
  TFile* FStatUnc                        = NULL;


  TH2D* hYieldStatUnc_Map = new TH2D("hYieldStatUnc_Map", "",numberbins, fBinsPi013TeVEMCPt, ((numberbins-1)/2), 1., numberbins);
  SetHistoStandardSettings2(hYieldStatUnc_Map);
  hYieldStatUnc_Map->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hYieldStatUnc_Map->GetXaxis()->SetRangeUser(1.4, 12.0);
  hYieldStatUnc_Map->SetYTitle("number of surrounding bins");
  hYieldStatUnc_Map->SetZTitle("relative statistical uncertainty (%)");

  FStatUnc           = SafelyOpenRootfile("FStatUnc.root");
  if (FStatUnc->IsOpen() ) printf("FStatUnc opened successfully\n");

  for (int i = 2; i < numberbins; i+=2) {
    hCorrYield_RelativStaterror = (TH1D*) FStatUnc->Get(Form("hCorrYield_RelativStaterror_with%02d_bins", i));
    for (int j = 0; j < hCorrYield_RelativStaterror->fNcells; j++) {
      hYieldStatUnc_Map->SetBinContent(j, i/2, hCorrYield_RelativStaterror->GetBinContent(j));
    }
  }

  TCanvas* cStatUnc = new TCanvas("cStatUnc", "", 1600,1600);
  TLine* l5 = new TLine(1.4, 5.0, 12.0, 5.0);
  l5->SetLineWidth(2);
  l5->SetLineColor(kBlack);
  TLine* l7 = new TLine(1.4, 7.0, 12.0, 7.0);
  l7->SetLineWidth(2);
  l7->SetLineColor(kBlack);
  TLine* l9 = new TLine(1.4, 9.0, 12.0, 9.0);
  l9->SetLineWidth(2);
  l9->SetLineColor(kBlack);
  SetCanvasStandardSettings(cStatUnc);
  cStatUnc->SetLeftMargin(0.15);
  cStatUnc->SetRightMargin(0.12);
  cStatUnc->SetBottomMargin(0.12);
  cStatUnc->SetLogz(1);
  cStatUnc->cd();

  hYieldStatUnc_Map->Draw("AXIS COLZ");
  hYieldStatUnc_Map->Draw("SAME COLZ");
  hYieldStatUnc_Map->Draw("SAME AXIS COLZ");
  l5->Draw("SAME");
  l7->Draw("SAME");
  l9->Draw("SAME");

  cStatUnc->Update();
  cStatUnc->SaveAs(Form("Systematics/StatUncertaintyMap." + PicFormat));
  cStatUnc->Clear();


  delete l5;
  delete l7;
  delete l9;
  delete hYieldStatUnc_Map;
  delete cStatUnc;


}
