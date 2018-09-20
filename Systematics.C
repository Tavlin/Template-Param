#include "CommonHeader.h"

void Systematics(TString PicFormat = "png", int numberneighbours = 4){

  TH1D* hCorrYield                       = NULL;
  TH1D* hCorrYield_count0d06to0d225      = NULL;
  TH1D* hCorrYield_count0d1to0d225       = NULL;
  TH1D* hCorrYield_count0d02to0d185      = NULL;
  TH1D* hCorrYield_count0d02to0d285      = NULL;
  TH1D* hCorrYield_countsyserror         = NULL;
  TH1D* hCorrYield_param0d06to0d225      = NULL;
  TH1D* hCorrYield_param0d1to0d225       = NULL;
  TH1D* hCorrYield_param0d02to0d185      = NULL;
  TH1D* hCorrYield_param0d02to0d285      = NULL;
  TH1D* hCorrYield_paramsyserror         = NULL;
  TH1D* hCorrYield_BGFitRange0d25        = NULL;
  TH1D* hCorrYield_BGFitRange0d29        = NULL;
  TH1D* hCorrYield_BGLeft                = NULL;
  TH1D* hCorrYield_BGFitsyserror         = NULL;
  TH1D* hCorrYield_syserror              = NULL;
  TH1D* hCorrYield_RelativSyserror       = NULL;
  TH1D* hCorrYield_RelativStaterror      = NULL;
  TH1D* hCorrectedYieldTrueEff           = NULL;
  TH1D* hCorrectedYieldTrueEff_StatError = NULL;
  TFile* FStatUnc                        = NULL;

  std::vector<Double_t> vCountSys;
  std::vector<Double_t> vParamSys;
  std::vector<Double_t> vBGFitSys;
  std::vector<Double_t> vFinalSys;

  TFile* IterTemp           = SafelyOpenRootfile("IterTempBetterBkg.root");
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
  hCorrYield_RelativStaterror = (TH1D*) hCorrYield->Clone("hCorrYield_RelativStaterror");
  hCorrectedYieldTrueEff_StatError = (TH1D*) hCorrectedYieldTrueEff->Clone("hCorrectedYieldTrueEff_StatError");

  for(int k = 2; k < numberbins; k++){
    hCorrYield_RelativSyserror->SetBinContent(k, hCorrYield_syserror->GetBinError(k)/(Double_t)hCorrYield->GetBinContent(k)*100.);
    std::cout << "rel Error = " <<  hCorrYield_RelativSyserror->GetBinContent(k) << '\n';
    hCorrYield_RelativStaterror->SetBinContent(k, hCorrYield_RelativStaterror->GetBinError(k)/(Double_t)hCorrYield->GetBinContent(k)*100.);
    hCorrYield_RelativStaterror->SetBinError(k,0.);
    hCorrectedYieldTrueEff_StatError->SetBinContent(k, hCorrectedYieldTrueEff->GetBinError(k)/(Double_t)hCorrectedYieldTrueEff->GetBinContent(k)*100.);
    hCorrectedYieldTrueEff_StatError->SetBinError(k,0.);
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
  TCanvas *canError = new TCanvas("canError","",2000,1000);
  TPad *pad1Error = new TPad("pad1Error","",0.0,0.53,1.0,1.0);
  pad1Error->SetTopMargin(0.05);
  pad1Error->SetLeftMargin(0.15);
  pad1Error->SetBottomMargin(0.0);
  pad1Error->SetRightMargin(0.02);
  pad1Error->SetTicky();
  pad1Error->SetTickx();
  TPad *pad2Error = new TPad("pad2Error","",0.0,0.0,1.0,0.53);
  pad2Error->SetTopMargin(0.0);
  pad2Error->SetLeftMargin(0.15);
  pad2Error->SetBottomMargin(0.18);
  pad2Error->SetRightMargin(0.02);
  pad2Error->SetTicky();
  pad2Error->SetTickx();

  ///////////////////////////////////////////////////////////////////////////
  // drwaing yields + ratios
  canError->cd();
  pad1Error->Draw();
  pad2Error->Draw("same");
  pad1Error->cd();



  // pad2InvMass->SetLogy(1);
  Double_t meanSys = 0;
  for(int i = 1; i < numberbins-3; i++){
    meanSys += hCorrYield_RelativSyserror->GetBinContent(i);
  }
  meanSys /= (Double_t)(numberbins-4);

  TLine* line_MeanSys = new TLine(1.4, meanSys, 12.0, meanSys);
  line_MeanSys->SetLineWidth(2);
  line_MeanSys->SetLineStyle(3);

  TLegend* leg = new TLegend(0.18,0.7,0.5,0.9);
  SetLegendSettigns(leg, 0.079);
  // leg->SetHeader("systematic uncertainties");
  leg->AddEntry(line_MeanSys, Form("mean value: %1.2lf %%", meanSys) , "l");
  leg->AddEntry(hCorrYield_RelativSyserror, "sum" , "l");

  TLegend* leg2 = new TLegend(0.55,0.6,0.9,0.9);
  SetLegendSettigns(leg2, 0.079);
  leg2->AddEntry(hCorrYield_paramsyserror, "param. range" , "l");
  leg2->AddEntry(hCorrYield_countsyserror, "integration range" , "l");
  leg2->AddEntry(hCorrYield_BGFitsyserror, "uncorr. bkg. variation" , "l");

  hCorrYield_RelativSyserror->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hCorrYield_RelativSyserror->SetYTitle("rel. syst. uncertainty (%)");
  hCorrYield_countsyserror->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hCorrYield_countsyserror->SetYTitle("rel. syst. uncertainty (%)");
  hCorrYield_paramsyserror->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hCorrYield_paramsyserror->SetYTitle("rel. syst. uncertainty (%)");
  hCorrYield_BGFitsyserror->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hCorrYield_BGFitsyserror->SetYTitle("rel. syst. uncertainty (%)");

  hCorrYield_RelativSyserror->SetLineWidth(3);
  hCorrYield_countsyserror->SetLineWidth(3);
  hCorrYield_paramsyserror->SetLineWidth(3);
  hCorrYield_BGFitsyserror->SetLineWidth(3);

  hCorrYield_RelativSyserror->GetXaxis()->SetTitleOffset(1.2);
  hCorrYield_RelativSyserror->GetXaxis()->SetLabelOffset(0.008);
  hCorrYield_RelativSyserror->GetYaxis()->SetTitleOffset(0.9);
  hCorrYield_RelativSyserror->GetYaxis()->SetLabelOffset(0.008);

  hCorrYield_RelativSyserror->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYield_paramsyserror->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYield_countsyserror->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYield_BGFitsyserror->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYield_RelativSyserror->GetYaxis()->SetRangeUser(-0.5, 9.9);
  hCorrYield_RelativSyserror->GetXaxis()->SetTitleSize(42);
  hCorrYield_RelativSyserror->GetYaxis()->SetTitleSize(42);
  hCorrYield_RelativSyserror->GetXaxis()->SetLabelSize(42);
  hCorrYield_RelativSyserror->GetYaxis()->SetLabelSize(42);

  hCorrYield_RelativSyserror->GetXaxis()->SetTitleFont(43);
  hCorrYield_RelativSyserror->GetYaxis()->SetTitleFont(43);
  hCorrYield_RelativSyserror->GetXaxis()->SetLabelFont(43);
  hCorrYield_RelativSyserror->GetYaxis()->SetLabelFont(43);

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

  pad1Error->Update();

  pad2Error->cd();

  hCorrYield_RelativStaterror->SetLineColor(kRed);
  hCorrYield_RelativStaterror->SetMarkerColor(kRed);
  hCorrectedYieldTrueEff_StatError->SetLineColor(kBlack);
  hCorrectedYieldTrueEff_StatError->SetMarkerColor(kBlack);

  TLegend* leg_stat = new TLegend(0.18,0.7,0.5,0.9);
  SetLegendSettigns(leg_stat, 0.079);
  leg_stat->SetTextFont(43);
  leg_stat->SetTextSize(42);
  leg_stat->AddEntry(hCorrectedYieldTrueEff_StatError, "standard method", "l");
  leg_stat->AddEntry(hCorrYield_RelativStaterror, "this method" , "l");
  hCorrectedYieldTrueEff_StatError->SetYTitle("rel. stat. uncertainty (%)");
  hCorrectedYieldTrueEff_StatError->SetXTitle("#it{p}_{T} (GeV/#it{c})");

  hCorrectedYieldTrueEff_StatError->SetLineWidth(3);
  hCorrYield_RelativStaterror->SetLineWidth(3);

  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetTitleOffset(1.8);
  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetLabelOffset(0.008);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetTitleOffset(0.9);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetLabelOffset(0.008);

  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYield_RelativStaterror->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetRangeUser(-0.09, 9.9);
  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetTitleSize(42);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetTitleSize(42);
  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetLabelSize(42);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetLabelSize(42);

  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetTitleFont(43);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetTitleFont(43);
  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetLabelFont(43);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetLabelFont(43);

  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetNdivisions(309);

  hCorrectedYieldTrueEff_StatError->DrawCopy("AXIS");
  hCorrYield_RelativStaterror->DrawCopy("SAME HIST");
  hCorrectedYieldTrueEff_StatError->DrawCopy("SAME HIST");
  leg_stat->Draw("");
  pad2Error->Update();


  canError->Update();
  canError->SaveAs("Systematics/ErrorPlot." + PicFormat);
  canError->Clear("D");

  delete leg;
  delete leg2;



  //////////////////////////////////////////////////////////////////////////////
  // setting up canvas to draw Yield plus relative Systematic error
  TCanvas *canInvMass = new TCanvas("canInvMass","",1600,1600);
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
  leg_yield->SetHeader("method:");
  leg_yield->AddEntry(hCorrYield, "templates parametrization" , "lp");
  leg_yield->AddEntry(hCorrectedYieldTrueEff, "function parametrization", "lp");
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
  hCorrectedYieldTrueEff->GetYaxis()->SetTitleOffset(1.7);
  hCorrectedYieldTrueEff->GetYaxis()->SetRangeUser(1.e-7-5.e-8,1.e-1);
  hCorrectedYieldTrueEff->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrectedYieldTrueEff->GetXaxis()->SetTitleSize(0.025*3./2.);
  hCorrectedYieldTrueEff->GetYaxis()->SetTitleSize(0.025*3./2.);
  hCorrectedYieldTrueEff->GetXaxis()->SetLabelSize(0.025*3./2.);
  hCorrectedYieldTrueEff->GetYaxis()->SetLabelSize(0.025*3./2.);
  hCorrectedYieldTrueEff->SetMarkerColor(kBlack);
  hCorrectedYieldTrueEff->SetLineColor(kBlack);
  hCorrYield_syserror->SetMarkerColor(kRed+3);
  hCorrYield_syserror->SetLineColor(kRed+3);
  hCorrYield_syserror->SetFillColor(kGray+2);
  hCorrYield_syserror->SetFillStyle(1001);

  hCorrYield->DrawCopy("AXIS");
  hCorrectedYieldTrueEff->Draw("SAME");
  // hCorrYield_syserror->DrawCopy("SAME E2");          // sys Error draw
  hCorrYield->DrawCopy("SAME");
  leg_yield->Draw("SAME");
  canInvMass->Update();
  DrawLabelALICE(0.6, 0.85, 0.035, 0.025*3./2., "");
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

  TLegend* leg4 = new TLegend(0.15,0.7,0.5,0.9);
  SetLegendSettigns(leg4, 0.079);
  // leg->SetHeader("systematic uncertainties");
  leg4->AddEntry(line_MeanSys, Form("mean value: %1.2lf %%", meanSys) , "l");
  leg4->AddEntry(hCorrYield_RelativSyserror, "sum" , "l");

  TLegend* leg5 = new TLegend(0.55,0.6,0.9,0.9);
  SetLegendSettigns(leg5, 0.079);
  leg5->AddEntry(hCorrYield_paramsyserror, "param. range" , "l");
  leg5->AddEntry(hCorrYield_countsyserror, "integration range" , "l");
  leg5->AddEntry(hCorrYield_BGFitsyserror, "uncorr. bkg. variation" , "l");

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
  leg4->Draw("SAME");
  leg5->Draw("SAME");

  c1->Update();
  c1->SaveAs("Systematics/RelativeSystematics." + PicFormat);
  c1->Clear();

  delete leg4;
  delete leg5;

  //////////////////////////////////////////////////////////////////////////////
  // setting up canvas to draw Yield plus relative Systematic error
  TCanvas *canYield = new TCanvas("canYield","",1000,1600);
  TPad *pad1Yield = new TPad("pad1Yield","",0.0,0.50,1.0,1.0);
  pad1Yield->SetTopMargin(0.05);
  pad1Yield->SetLeftMargin(0.21);
  pad1Yield->SetBottomMargin(0.0);
  pad1Yield->SetRightMargin(0.02);
  pad1Yield->SetTicky();
  pad1Yield->SetTickx();
  pad1Yield->SetLogy(1);
  TPad *pad2Yield = new TPad("pad2Yield","",0.0,0.28,1.0,0.50);
  pad2Yield->SetTopMargin(0.0);
  pad2Yield->SetLeftMargin(0.21);
  pad2Yield->SetBottomMargin(0.0);
  pad2Yield->SetRightMargin(0.02);
  pad2Yield->SetTicky();
  pad2Yield->SetTickx();
  TPad *pad3Yield = new TPad("pad3Yield","",0.0,0.0,1.0,0.28);
  pad3Yield->SetTopMargin(0.0);
  pad3Yield->SetLeftMargin(0.21);
  pad3Yield->SetBottomMargin(0.3);
  pad3Yield->SetRightMargin(0.02);
  pad3Yield->SetTicky();
  pad3Yield->SetTickx();

  canYield->cd();
  pad1Yield->Draw();
  pad2Yield->Draw("same");
  pad3Yield->Draw("same");

  pad1Yield->cd();

  TLegend* leg_yield2 = new TLegend(0.24,0.07,0.35,0.3);
  SetLegendSettigns(leg_yield2, 0.025*2.);
  leg_yield2->SetTextFont(43);
  leg_yield2->SetTextSize(42);
  leg_yield2->AddEntry(hCorrYield, "signal + corr. bkg. temp." , "lp");
  leg_yield2->AddEntry(hCorrectedYieldTrueEff, "standard method", "lp");


  hCorrYield->GetXaxis()->SetTitleOffset(1.4);
  hCorrYield->GetXaxis()->SetLabelOffset(0.008);
  hCorrYield->GetYaxis()->SetTitleOffset(3.5);
  hCorrYield->GetYaxis()->SetLabelOffset(0.008);
  hCorrYield->GetXaxis()->SetTitleSize(42);
  hCorrYield->GetYaxis()->SetTitleSize(42);
  hCorrYield->GetXaxis()->SetLabelSize(42);
  hCorrYield->GetYaxis()->SetLabelSize(42);

  hCorrYield->GetXaxis()->SetTitleFont(43);
  hCorrYield->GetYaxis()->SetTitleFont(43);
  hCorrYield->GetXaxis()->SetLabelFont(43);
  hCorrYield->GetYaxis()->SetLabelFont(43);


  hCorrYield->DrawCopy("AXIS");
  hCorrYield_syserror->DrawCopy("SAME E2");
  hCorrYield->DrawCopy("SAME");
  leg_yield2->Draw("SAME");
  canYield->Update();
  DrawLabelALICE(0.55, 0.85, 0.035, 0.025*2., "");
  pad1Yield->Update();

  pad2Yield->cd();

  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetRangeUser(0.76,1.26);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetNdivisions(505);
  hYield_dt_chi2map_corrected_ratio->SetXTitle("");
  hYield_dt_chi2map_corrected_ratio->GetXaxis()->SetLabelOffset(0.008);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetTitleOffset(3.5);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetLabelOffset(0.008);
  hYield_dt_chi2map_corrected_ratio->GetXaxis()->SetTitleSize(42);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetTitleSize(42);
  hYield_dt_chi2map_corrected_ratio->GetXaxis()->SetLabelSize(42);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetLabelSize(42);

  hYield_dt_chi2map_corrected_ratio->GetXaxis()->SetTitleFont(43);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetTitleFont(43);
  hYield_dt_chi2map_corrected_ratio->GetXaxis()->SetLabelFont(43);
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetLabelFont(43);

  hYield_dt_chi2map_corrected_ratio->DrawCopy("AXIS");
  line_ratio1->Draw("SAME");
  hYield_dt_chi2map_corrected_ratio->DrawCopy("SAME P");
  pad2Yield->Update();

  pad3Yield->cd();

  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetNdivisions(505);

  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetTitleOffset(3.5);
  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetLabelOffset(0.008);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetTitleOffset(3.5);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetLabelOffset(0.008);
  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetTitleSize(42);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetTitleSize(42);
  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetLabelSize(42);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetLabelSize(42);

  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetTitleFont(43);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetTitleFont(43);
  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetLabelFont(43);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetLabelFont(43);

  TLegend* leg_stat_yield = new TLegend(0.24,0.7,0.5,0.9);
  SetLegendSettigns(leg_stat_yield, 0.079);
  leg_stat_yield->SetTextFont(43);
  leg_stat_yield->SetTextSize(42);
  leg_stat_yield->AddEntry(hCorrectedYieldTrueEff_StatError, "standard method", "l");
  leg_stat_yield->AddEntry(hCorrYield_RelativStaterror, "this method" , "l");

  hCorrectedYieldTrueEff_StatError->DrawCopy("AXIS");
  hCorrYield_RelativStaterror->DrawCopy("SAME HIST");
  hCorrectedYieldTrueEff_StatError->DrawCopy("SAME HIST");
  leg_stat_yield->Draw("");
  pad3Yield->Update();

  canYield->Update();
  canYield->SaveAs(Form("Systematics/CorrectedYieldCompWithStat." + PicFormat));
  canYield->Clear("D");


  delete leg_yield2;
  delete leg_stat_yield;
  delete leg_stat;
  delete line_ratio1;

  TString sPath = gDirectory->GetPath();
  if(numberneighbours == 2){
    FStatUnc      = new TFile("FStatUnc.root", "RECREATE");
  }

  else{
    FStatUnc      = new TFile("FStatUnc.root", "UPDATE");
  }

  hCorrYield_RelativStaterror->Write(Form("hCorrYield_RelativStaterror_with%02d_bins", numberneighbours));
  gDirectory->Cd(sPath.Data());

  TCanvas* cStatUnc = new TCanvas("cStatUnc", "", 2500,1000);
  SetCanvasStandardSettings(cStatUnc);
  cStatUnc->cd();

  TLegend* leg_stat_yield2 = new TLegend(0.20,0.75,0.5,0.9);
  SetLegendSettigns(leg_stat_yield2, 0.079);
  leg_stat_yield2->SetTextFont(43);
  leg_stat_yield2->SetTextSize(42);
  leg_stat_yield2->SetHeader("method:");
  leg_stat_yield2->AddEntry(hCorrectedYieldTrueEff_StatError, "function parametrization", "l");
  leg_stat_yield2->AddEntry(hCorrYield_RelativStaterror, "templates parametrization" , "l");

  hCorrectedYieldTrueEff_StatError->GetXaxis()->SetTitleOffset(1.2);
  hCorrectedYieldTrueEff_StatError->GetYaxis()->SetTitleOffset(1.2);

  hCorrectedYieldTrueEff_StatError->DrawCopy("AXIS");
  hCorrYield_RelativStaterror->DrawCopy("SAME HIST");
  hCorrectedYieldTrueEff_StatError->DrawCopy("SAME HIST");
  leg_stat_yield2->Draw("");

  cStatUnc->Update();
  cStatUnc->SaveAs(Form("Systematics/StatUncertainty." + PicFormat));
  cStatUnc->Clear();

  delete leg_stat_yield2;

  delete hCorrYield;
  delete cStatUnc;
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
