#include "CommonHeader.h"

void Systematics(TString PicFormat = "png"){

  TH1D* hCorrYield                   = NULL;
  TH1D* hCorrYield_count0d06to0d225  = NULL;
  TH1D* hCorrYield_count0d1to0d225   = NULL;
  TH1D* hCorrYield_count0d02to0d185  = NULL;
  TH1D* hCorrYield_count0d02to0d285  = NULL;
  TH1D* hCorrYield_countsyserror     = NULL;
  TH1D* hCorrYield_param0d06to0d225  = NULL;
  TH1D* hCorrYield_param0d1to0d225   = NULL;
  TH1D* hCorrYield_param0d02to0d185  = NULL;
  TH1D* hCorrYield_param0d02to0d285  = NULL;
  TH1D* hCorrYield_paramsyserror     = NULL;
  TH1D* hCorrYield_BGFitRange0d25    = NULL;
  TH1D* hCorrYield_BGFitRange0d29    = NULL;
  TH1D* hCorrYield_BGLeft            = NULL;
  TH1D* hCorrYield_BGFitsyserror     = NULL;
  TH1D* hCorrYield_syserror          = NULL;
  TH1D* hCorrYield_RelativSyserror   = NULL;
  TH1D* hCorrectedYieldTrueEff       = NULL;

  std::vector<Double_t> vCountSys;
  std::vector<Double_t> vParamSys;
  std::vector<Double_t> vBGFitSys;
  std::vector<Double_t> vFinalSys;

  TFile* IterTemp           = SafelyOpenRootfile("IterTemp.root");
  if (IterTemp->IsOpen() ) printf("IterTemp opened successfully\n");

  hCorrYield                = (TH1D*) IterTemp->Get("hYield_dt_chi2map_corrected");

  hCorrectedYieldTrueEff    = (TH1D*) IterTemp->Get("hCorrectedYieldTrueEff");

  hCorrYield_syserror       = (TH1D*) hCorrYield->Clone("hCorrYield_syserror");
  hCorrYield_countsyserror  = (TH1D*) hCorrYield->Clone("hCorrYield_countsyserror");
  hCorrYield_paramsyserror  = (TH1D*) hCorrYield->Clone("hCorrYield_paramsyserror");
  hCorrYield_BGFitsyserror  = (TH1D*) hCorrYield->Clone("hCorrYield_BGFitsyserror");

  //////////////////////////////////////////////////////////////////////////////
  // First the counting range variation:
  TFile* IterTemp_count0d06to0d225 = SafelyOpenRootfile("IterTemp_count0d06to0d225.root");
  if (IterTemp_count0d06to0d225->IsOpen() ) printf("IterTemp_count0d06to0d225 opened successfully\n");

  hCorrYield_count0d06to0d225 = (TH1D*) IterTemp_count0d06to0d225->Get("hYield_dt_chi2map_corrected");

  TFile* IterTemp_count0d1to0d225 = SafelyOpenRootfile("IterTemp_count0d1to0d225.root");
  if (IterTemp_count0d1to0d225->IsOpen() ) printf("IterTemp_count0d1to0d225 opened successfully\n");

  hCorrYield_count0d1to0d225 = (TH1D*) IterTemp_count0d1to0d225->Get("hYield_dt_chi2map_corrected");

  TFile* IterTemp_count0d02to0d185 = SafelyOpenRootfile("IterTemp_count0d02to0d185.root");
  if (IterTemp_count0d02to0d185->IsOpen() ) printf("IterTemp_count0d02to0d185 opened successfully\n");

  hCorrYield_count0d02to0d185 = (TH1D*) IterTemp_count0d02to0d185->Get("hYield_dt_chi2map_corrected");

  TFile* IterTemp_count0d02to0d285 = SafelyOpenRootfile("IterTemp_count0d02to0d285.root");
  if (IterTemp_count0d02to0d285->IsOpen() ) printf("IterTemp_count0d02to0d285 opened successfully\n");

  hCorrYield_count0d02to0d285 = (TH1D*) IterTemp_count0d02to0d285->Get("hYield_dt_chi2map_corrected");




  //////////////////////////////////////////////////////////////////////////////
  // 2nd the param range variation:
  TFile* IterTemp_param0d06to0d225 = SafelyOpenRootfile("IterTemp_param0d06to0d225.root");
  if (IterTemp_param0d06to0d225->IsOpen() ) printf("IterTemp_param0d06to0d225 opened successfully\n");

  hCorrYield_param0d06to0d225 = (TH1D*) IterTemp_param0d06to0d225->Get("hYield_dt_chi2map_corrected");

  TFile* IterTemp_param0d1to0d225 = SafelyOpenRootfile("IterTemp_param0d1to0d225.root");
  if (IterTemp_param0d1to0d225->IsOpen() ) printf("IterTemp_param0d1to0d225 opened successfully\n");

  hCorrYield_param0d1to0d225 = (TH1D*) IterTemp_param0d1to0d225->Get("hYield_dt_chi2map_corrected");

  TFile* IterTemp_param0d02to0d185 = SafelyOpenRootfile("IterTemp_param0d02to0d185.root");
  if (IterTemp_param0d02to0d185->IsOpen() ) printf("IterTemp_param0d02to0d185 opened successfully\n");

  hCorrYield_param0d02to0d185 = (TH1D*) IterTemp_param0d02to0d185->Get("hYield_dt_chi2map_corrected");

  TFile* IterTemp_param0d02to0d285 = SafelyOpenRootfile("IterTemp_param0d02to0d285.root");
  if (IterTemp_param0d02to0d285->IsOpen() ) printf("IterTemp_param0d02to0d285 opened successfully\n");

  hCorrYield_param0d02to0d285 = (TH1D*) IterTemp_param0d02to0d285->Get("hYield_dt_chi2map_corrected");



  //////////////////////////////////////////////////////////////////////////////
  // 3rd the bg fitting variation
  TFile* IterTemp_BGFitRange0d25 = SafelyOpenRootfile("IterTemp_BGFitRange0d25.root");
  if (IterTemp_BGFitRange0d25->IsOpen() ) printf("IterTemp_BGFitRange0d25 opened successfully\n");

  hCorrYield_BGFitRange0d25 = (TH1D*) IterTemp_BGFitRange0d25->Get("hYield_dt_chi2map_corrected");

  TFile* IterTemp_BGLeft = SafelyOpenRootfile("IterTemp_BGLeft.root");
  if (IterTemp_BGLeft->IsOpen() ) printf("IterTemp_BGLeft opened successfully\n");

  hCorrYield_BGLeft = (TH1D*) IterTemp_BGFitRange0d25->Get("hYield_dt_chi2map_corrected");

  TFile* IterTemp_BGFitRange0d29 = SafelyOpenRootfile("IterTemp_BGFitRange0d29.root");
  if (IterTemp_BGFitRange0d29->IsOpen() ) printf("IterTemp_BGFitRange0d29 opened successfully\n");

  hCorrYield_BGFitRange0d29 = (TH1D*) IterTemp_BGFitRange0d29->Get("hYield_dt_chi2map_corrected");


  Double_t temp = 0;
  for (int i = 2; i < numberbins; i++) {

    temp = 0;
    // check for biggest diff. in param vari
    if(fabs(hCorrYield->GetBinContent(i)-hCorrYield_param0d02to0d285->GetBinContent(i)) > temp){
      temp = fabs(hCorrYield->GetBinContent(i)-hCorrYield_param0d02to0d285->GetBinContent(i));
    }
    if(fabs(hCorrYield->GetBinContent(i)-hCorrYield_param0d06to0d225->GetBinContent(i)) > temp){
      temp = fabs(hCorrYield->GetBinContent(i)-hCorrYield_param0d06to0d225->GetBinContent(i));
    }
    if(fabs(hCorrYield->GetBinContent(i)-hCorrYield_param0d1to0d225->GetBinContent(i)) > temp){
      temp = fabs(hCorrYield->GetBinContent(i)-hCorrYield_param0d1to0d225->GetBinContent(i));
    }
    if(fabs(hCorrYield->GetBinContent(i)-hCorrYield_param0d02to0d185->GetBinContent(i)) > temp){
      temp = fabs(hCorrYield->GetBinContent(i)-hCorrYield_param0d02to0d185->GetBinContent(i));
    }
    // pushing biggest difference back
    vParamSys.push_back(temp);

    // resetting temp
    temp = 0;

    // chech for biggest diff. in count vari
    if(fabs(hCorrYield->GetBinContent(i)-hCorrYield_count0d06to0d225->GetBinContent(i)) > temp){
      temp = fabs(hCorrYield->GetBinContent(i)-hCorrYield_count0d06to0d225->GetBinContent(i));
    }
    if(fabs(hCorrYield->GetBinContent(i)-hCorrYield_count0d1to0d225->GetBinContent(i)) > temp){
      temp = fabs(hCorrYield->GetBinContent(i)-hCorrYield_count0d1to0d225->GetBinContent(i));
    }
    if(fabs(hCorrYield->GetBinContent(i)-hCorrYield_count0d02to0d185->GetBinContent(i)) > temp){
      temp = fabs(hCorrYield->GetBinContent(i)-hCorrYield_count0d02to0d185->GetBinContent(i));
    }
    if(fabs(hCorrYield->GetBinContent(i)-hCorrYield_count0d02to0d285->GetBinContent(i)) > temp){
      temp = fabs(hCorrYield->GetBinContent(i)-hCorrYield_count0d02to0d285->GetBinContent(i));
    }
    // pushing biggest difference back
    vCountSys.push_back(temp);

    // chech for biggest diff. in BGGit vari

    if(fabs(hCorrYield->GetBinContent(i)-hCorrYield_BGFitRange0d25->GetBinContent(i)) > temp){
      temp = fabs(hCorrYield->GetBinContent(i)-hCorrYield_BGFitRange0d25->GetBinContent(i));
    }
    if(fabs(hCorrYield->GetBinContent(i)-hCorrYield_BGLeft->GetBinContent(i)) > temp){
      temp = fabs(hCorrYield->GetBinContent(i)-hCorrYield_BGLeft->GetBinContent(i));
    }
    // if(fabs(hCorrYield->GetBinContent(i)-hCorrYield_BGFitRange0d29->GetBinContent(i)) > temp){
    //   temp = fabs(hCorrYield->GetBinContent(i)-hCorrYield_BGFitRange0d29->GetBinContent(i));
    // }

    // pushing biggest difference back
    vBGFitSys.push_back(temp);

    temp = sqrt(pow(vParamSys[i-2], 2.) + pow(vCountSys[i-2], 2.) + pow(vBGFitSys[i-2], 2.));
    vFinalSys.push_back(temp);
    hCorrYield_syserror->SetBinError(i, vFinalSys[i-2]);
    hCorrYield_countsyserror->SetBinError(i, vCountSys[i-2]);
    hCorrYield_paramsyserror->SetBinError(i, vParamSys[i-2]);
    hCorrYield_BGFitsyserror->SetBinError(i, vBGFitSys[i-2]);
    temp = 0;
  }

  hCorrYield_RelativSyserror = (TH1D*) hCorrYield_syserror->Clone("hCorrYield_RelativSyserror");
  for(int k = 2; k < numberbins; k++){
    hCorrYield_RelativSyserror->SetBinContent(k, hCorrYield_syserror->GetBinError(k)/(Double_t)hCorrYield->GetBinContent(k)*100.);
    std::cout << "rel Error = " <<  hCorrYield_RelativSyserror->GetBinContent(k) << '\n';
    hCorrYield_RelativSyserror->SetBinError(k,0.);
    hCorrYield_countsyserror->SetBinContent(k, hCorrYield_countsyserror->GetBinError(k)/(Double_t)hCorrYield->GetBinContent(k)*100.);
    hCorrYield_paramsyserror->SetBinContent(k, hCorrYield_paramsyserror->GetBinError(k)/(Double_t)hCorrYield->GetBinContent(k)*100.);
    hCorrYield_BGFitsyserror->SetBinContent(k, hCorrYield_BGFitsyserror->GetBinError(k)/(Double_t)hCorrYield->GetBinContent(k)*100.);
    hCorrYield_countsyserror->SetBinError(k,0.);
    hCorrYield_paramsyserror->SetBinError(k,0.);
    hCorrYield_BGFitsyserror->SetBinError(k,0.);
  }

  //////////////////////////////////////////////////////////////////////////////
  // setting up canvas to draw Yield plus relative Systematic error
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

  hCorrYield_syserror->SetMarkerSize(1.5);
  hCorrYield->SetMarkerSize(1.5);

  ////////////////////////////////////////////////////////////////////////////
  // drwaing yields + ratios
  canInvMass->cd();
  pad1InvMass->Draw();
  pad2InvMass->Draw("same");
  pad1InvMass->cd();
  pad1InvMass->SetTickx();
  pad1InvMass->SetTicky();

  pad1InvMass->SetLogy(1);

  TLegend* leg_yield = new TLegend(0.2,0.07,0.35,0.3);
  SetLegendSettigns(leg_yield, 0.025*3./2.);
  leg_yield->AddEntry(hCorrYield, "signal + back. temp." , "lp");
  leg_yield->AddEntry(hCorrectedYieldTrueEff, "standard method", "lp");
  hCorrectedYieldTrueEff->SetMarkerSize(1.5);


  hCorrYield->GetYaxis()->SetTitleOffset(1.7);
  hCorrYield->GetYaxis()->SetRangeUser(1.e-7-5.e-8,1.e-1);
  hCorrYield->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYield->GetXaxis()->SetTitleSize(0.025*3./2.);
  hCorrYield->GetYaxis()->SetTitleSize(0.025*3./2.);
  hCorrYield->GetXaxis()->SetLabelSize(0.025*3./2.);
  hCorrYield->GetYaxis()->SetLabelSize(0.025*3./2.);
  hCorrYield->SetMarkerColor(kRed+3);
  hCorrYield->SetLineColor(kRed+3);
  hCorrYield_syserror->SetMarkerColor(kRed+3);
  hCorrYield_syserror->SetLineColor(kRed+3);
  hCorrYield_syserror->SetFillColor(kGray+2);
  hCorrYield_syserror->SetFillStyle(1001);

  hCorrYield->Draw("AXIS");
  hCorrYield_syserror->Draw("SAME E2");
  hCorrYield->Draw("SAME");
  leg_yield->Draw("SAME");
  canInvMass->Update();
  DrawLabelALICE(0.64, 0.85, 0.035, 0.025*3./2., "");
  pad1InvMass->Update();


  TH1D* hYield_dt_chi2map_corrected_ratio = (TH1D*) hCorrectedYieldTrueEff->Clone("hYield_dt_chi2map_corrected_ratio");
  hYield_dt_chi2map_corrected_ratio->Divide(hCorrYield);
  hYield_dt_chi2map_corrected_ratio->SetLineColor(kRed);
  hYield_dt_chi2map_corrected_ratio->SetMarkerColor(kRed);
  hYield_dt_chi2map_corrected_ratio->SetYTitle("Ratio");


  pad2InvMass->cd();
  TLine* line_ratio1 = new TLine(1.4, 1.0, 12.0, 1.0);
  line_ratio1->SetLineWidth(2);
  line_ratio1->SetLineStyle(3);

  hYield_dt_chi2map_corrected_ratio->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hYield_dt_chi2map_corrected_ratio->GetXaxis()->SetTitleOffset(1.0);
  hYield_dt_chi2map_corrected_ratio->GetXaxis()->SetLabelOffset(0.008);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetTitleOffset(0.8);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetLabelOffset(0.008);

  hYield_dt_chi2map_corrected_ratio->GetXaxis()->SetRangeUser(1.4, 12.0);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetRangeUser(0.69, 1.25);
  hYield_dt_chi2map_corrected_ratio->GetXaxis()->SetTitleSize(0.025*3./1.);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetTitleSize(0.025*3./1.);
  hYield_dt_chi2map_corrected_ratio->GetXaxis()->SetLabelSize(0.025*3./1.);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetLabelSize(0.025*3./1.);

  hYield_dt_chi2map_corrected_ratio->GetXaxis()->SetTitleFont(42);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetTitleFont(42);
  hYield_dt_chi2map_corrected_ratio->GetXaxis()->SetLabelFont(42);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetLabelFont(42);

  hYield_dt_chi2map_corrected_ratio->DrawCopy("AXIS");
  line_ratio1->Draw("SAME");
  hYield_dt_chi2map_corrected_ratio->DrawCopy("SAME P");
  pad2InvMass->Update();

  pad2InvMass->SetTickx();
  pad2InvMass->SetTicky();

  pad2InvMass->Update();

  canInvMass->Update();
  canInvMass->SaveAs(Form("Systematics/CorrectedYieldComp." + PicFormat));
  canInvMass->Clear("D");

  delete line_ratio1;
  delete leg_yield;

  //////////////////////////////////////////////////////////////////////////////
  // setting up the canvas to draw on. Will later be changed for the chi2 pic
  TCanvas *c1 = new TCanvas("c1","",1200,600);
  c1->cd();
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.18);
  c1->SetRightMargin(0.02);
  c1->SetLeftMargin(0.10);
  c1->SetTicky();
  c1->SetTickx();
  c1->SetLogz(1);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);


  // pad2InvMass->SetLogy(1);
  Double_t meanSys = 0;
  for(int i = 1; i < numberbins-3; i++){
    meanSys += hCorrYield_RelativSyserror->GetBinContent(i);
  }
  meanSys /= (Double_t)(numberbins-4);

  TLine* line_MeanSys = new TLine(1.4, meanSys, 12.0, meanSys);
  line_MeanSys->SetLineWidth(2);
  line_MeanSys->SetLineStyle(3);

  TLegend* leg = new TLegend(0.15,0.7,0.5,0.9);
  SetLegendSettigns(leg, 0.079);
  // leg->SetHeader("systematic uncertainties");
  leg->AddEntry(line_MeanSys, Form("mean value: %1.2lf %%", meanSys) , "l");
  leg->AddEntry(hCorrYield_RelativSyserror, "sum" , "l");

  TLegend* leg2 = new TLegend(0.55,0.6,0.9,0.9);
  SetLegendSettigns(leg2, 0.079);
  leg2->AddEntry(hCorrYield_paramsyserror, "param. range" , "l");
  leg2->AddEntry(hCorrYield_countsyserror, "integration range" , "l");
  leg2->AddEntry(hCorrYield_BGFitsyserror, "uncorr. back. variation" , "l");

  hCorrYield_RelativSyserror->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hCorrYield_RelativSyserror->SetYTitle("rel. syst. uncertainty (%)");
  hCorrYield_countsyserror->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hCorrYield_countsyserror->SetYTitle("rel. syst. uncertainty (%)");
  hCorrYield_paramsyserror->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hCorrYield_paramsyserror->SetYTitle("rel. syst. uncertainty (%)");
  hCorrYield_BGFitsyserror->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hCorrYield_BGFitsyserror->SetYTitle("rel. syst. uncertainty (%)");

  hCorrYield_RelativSyserror->GetXaxis()->SetTitleOffset(0.9);
  hCorrYield_RelativSyserror->GetXaxis()->SetLabelOffset(0.008);
  hCorrYield_RelativSyserror->GetYaxis()->SetTitleOffset(0.5);
  hCorrYield_RelativSyserror->GetYaxis()->SetLabelOffset(0.008);

  hCorrYield_RelativSyserror->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYield_paramsyserror->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYield_countsyserror->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYield_BGFitsyserror->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYield_RelativSyserror->GetYaxis()->SetRangeUser(0, 9.9);
  hCorrYield_RelativSyserror->GetXaxis()->SetTitleSize(0.085);
  hCorrYield_RelativSyserror->GetYaxis()->SetTitleSize(0.085);
  hCorrYield_RelativSyserror->GetXaxis()->SetLabelSize(0.085);
  hCorrYield_RelativSyserror->GetYaxis()->SetLabelSize(0.085);

  hCorrYield_RelativSyserror->GetXaxis()->SetTitleFont(42);
  hCorrYield_RelativSyserror->GetYaxis()->SetTitleFont(42);
  hCorrYield_RelativSyserror->GetXaxis()->SetLabelFont(42);
  hCorrYield_RelativSyserror->GetYaxis()->SetLabelFont(42);

  hCorrYield_countsyserror->SetLineColor(kMagenta);
  hCorrYield_paramsyserror->SetLineColor(kBlue+2);
  hCorrYield_BGFitsyserror->SetLineColor(kTeal-7);
  hCorrYield_RelativSyserror->SetLineColor(kRed);

  hCorrYield_RelativSyserror->DrawCopy("AXIS");
  line_MeanSys->Draw("SAME");
  hCorrYield_countsyserror->DrawCopy("SAME HIST");
  hCorrYield_paramsyserror->DrawCopy("SAME HIST");
  hCorrYield_BGFitsyserror->DrawCopy("SAME HIST");
  hCorrYield_RelativSyserror->DrawCopy("SAME HIST");
  leg->Draw("SAME");
  leg2->Draw("SAME");

  c1->Update();
  c1->SaveAs("Systematics/RelativeSystematics." + PicFormat);
  c1->Clear();

  delete leg;
  delete leg2;
  delete hCorrYield;
  // delete hCorrYield_count0d06to0d225;
  // delete hCorrYield_count0d1to0d225;
  // delete hCorrYield_count0d02to0d185;
  // delete hCorrYield_count0d02to0d285;
  // delete hCorrYield_param0d06to0d225;
  // delete hCorrYield_param0d1to0d225;
  // delete hCorrYield_param0d02to0d185;
  // delete hCorrYield_param0d02to0d285;
  // delete hCorrYield_BGFitRange0d25;
  // delete hCorrYield_BGLeft;
  // delete hCorrYield_syserror;
  delete hCorrectedYieldTrueEff;

}
