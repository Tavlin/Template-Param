#include "CommonHeader.h"

void Systematics(TString PicFormat = "png", int numberneighbours = 4){

  TH1D* hCorrYieldNormal                 = NULL;
  TH1D* hCorrYieldBetterBkgNN            = NULL;
  TH1D* hCorrYieldBetterBkg3to8          = NULL;
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
  TH1D* hCorrectedYieldNormEff           = NULL;
  TH1D* hCorrectedYieldNormEff_StatError = NULL;
  TFile* FStatUnc                        = NULL;

  TF1* fitBylikin13TeV = new TF1("fitBylikin13TeV", "[0]*exp(-(sqrt(x^(2)+0.135^(2))-0.135)/[1])+[2]/((1+x*x/[3])^([4]))", 1.4, 12.);
  fitBylikin13TeV->SetParameters(13, 0.1, 2, 0.7, 2.9);
  fitBylikin13TeV->SetLineWidth(3);
  fitBylikin13TeV->SetLineColor(kBlack);

  TF1* fitBylikin13TeV_3to8 = new TF1("fitBylikin13TeV_3to8", "[0]*exp(-(sqrt(x^(2)+0.135^(2))-0.135)/[1])+[2]/((1+x*x/[3])^([4]))", 1.4, 12.);
  fitBylikin13TeV_3to8->SetParameters(13, 0.1, 2, 0.7, 2.9);
  fitBylikin13TeV_3to8->SetLineWidth(3);
  fitBylikin13TeV_3to8->SetLineColor(kBlue+2);

  TF1 *ftsallis13TeV = new TF1("ftsallis13TeV", "[0]/(2*3.1415) * ( ([1]-1) * ([1]-2) )/( [1]*[2] *( [1]*[2] + [3] *([1]-2)) )* pow(( 1 + ( ( pow(([3]*[3]+x*x),0.5) -[3]) /( [1]*[2] ) ) ), -[1])", 1.4, 12.);
  ftsallis13TeV->SetParameter(0, 9.4); // A
  ftsallis13TeV->SetParameter(1, 7.169); // n
  ftsallis13TeV->SetParameter(2, 0.159); // T
  ftsallis13TeV->FixParameter(3, 0.1349766); // M
  ftsallis13TeV->SetLineWidth(3);
  ftsallis13TeV->SetLineColor(kBlack);


  std::vector<Double_t> vCountSys;
  std::vector<Double_t> vParamSys;
  std::vector<Double_t> vBGFitSys;
  std::vector<Double_t> vFinalSys;

  TFile* IterTempNormal         = SafelyOpenRootfile("IterTemp.root");
  if (IterTempNormal->IsOpen() ) printf("IterTempNormal opened successfully\n");

  TFile* IterTempBetterBkgNN    = SafelyOpenRootfile("IterTempBetterBkgNN.root");
  if (IterTempBetterBkgNN->IsOpen() ) printf("IterTempBetterBkgNN opened successfully\n");

  TFile* IterTempBetterBkg3to8  = SafelyOpenRootfile("IterTempBetterBkg3to8_WithFit.root");
  if (IterTempBetterBkg3to8->IsOpen() ) printf("IterTempBetterBkg3to8 opened successfully\n");

  TFile* JoshuasFile  = SafelyOpenRootfile("Pi0_data_GammaConvV1Correction_00010113_1111112067032220000_01631031000000d0.root");
  if (JoshuasFile->IsOpen() ) printf("JoshuasFile opened successfully\n");

  hCorrYieldNormal              = (TH1D*) IterTempNormal->Get("hYield_dt_chi2map_corrected");
  hCorrYieldBetterBkgNN         = (TH1D*) IterTempBetterBkgNN->Get("hYield_dt_chi2map_corrected");
  hCorrYieldBetterBkg3to8       = (TH1D*) IterTempBetterBkg3to8->Get("hYield_dt_chi2map_corrected");

  hCorrectedYieldNormEff        = (TH1D*) IterTempNormal->Get("hCorrectedYieldNormEff");
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // hCorrectedYieldNormEff->Rebin(39, "", fBinsPi013TeVEMCPt);
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  hCorrectedYieldNormEff->Fit(fitBylikin13TeV, "QM0P", "", 1.4, 12.);
  hCorrYieldBetterBkg3to8->Fit(fitBylikin13TeV_3to8, "QM0P", "", 1.4, 12.);
  hCorrectedYieldNormEff->Fit(ftsallis13TeV, "QM0P", "", 1.4, 12.);

  TF1 *fBylikinRatio = new TF1("fBylikinRatio", "fitBylikin13TeV_3to8/fitBylikin13TeV", 1.4, 12.);
  fBylikinRatio->SetLineWidth(3);
  fBylikinRatio->SetLineColor(kBlue+2);

  hCorrYield_syserror           = (TH1D*) hCorrYieldNormal->Clone("hCorrYield_syserror");
  hCorrYield_countsyserror      = (TH1D*) hCorrYieldNormal->Clone("hCorrYield_countsyserror");
  hCorrYield_paramsyserror      = (TH1D*) hCorrYieldNormal->Clone("hCorrYield_paramsyserror");
  hCorrYield_BGFitsyserror      = (TH1D*) hCorrYieldNormal->Clone("hCorrYield_BGFitsyserror");

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
    if(fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_param0d02to0d285->GetBinContent(i)) > temp){
      temp = fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_param0d02to0d285->GetBinContent(i));
    }
    if(fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_param0d06to0d225->GetBinContent(i)) > temp){
      temp = fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_param0d06to0d225->GetBinContent(i));
    }
    if(fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_param0d1to0d225->GetBinContent(i)) > temp){
      temp = fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_param0d1to0d225->GetBinContent(i));
    }
    if(fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_param0d02to0d185->GetBinContent(i)) > temp){
      temp = fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_param0d02to0d185->GetBinContent(i));
    }
    // pushing biggest difference back
    vParamSys.push_back(temp);

    // resetting temp
    temp = 0;

    // chech for biggest diff. in count vari
    if(fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_count0d06to0d225->GetBinContent(i)) > temp){
      temp = fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_count0d06to0d225->GetBinContent(i));
    }
    if(fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_count0d1to0d225->GetBinContent(i)) > temp){
      temp = fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_count0d1to0d225->GetBinContent(i));
    }
    if(fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_count0d02to0d185->GetBinContent(i)) > temp){
      temp = fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_count0d02to0d185->GetBinContent(i));
    }
    if(fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_count0d02to0d285->GetBinContent(i)) > temp){
      temp = fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_count0d02to0d285->GetBinContent(i));
    }
    // pushing biggest difference back
    vCountSys.push_back(temp);

    // chech for biggest diff. in BGGit vari

    if(fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_BGFitRange0d25->GetBinContent(i)) > temp){
      temp = fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_BGFitRange0d25->GetBinContent(i));
    }
    if(fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_BGLeft->GetBinContent(i)) > temp){
      temp = fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_BGLeft->GetBinContent(i));
    }
    // if(fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_BGFitRange0d29->GetBinContent(i)) > temp){
    //   temp = fabs(hCorrYieldNormal->GetBinContent(i)-hCorrYield_BGFitRange0d29->GetBinContent(i));
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
  hCorrYield_RelativStaterror = (TH1D*) hCorrYieldNormal->Clone("hCorrYield_RelativStaterror");
  hCorrectedYieldNormEff_StatError = (TH1D*) hCorrectedYieldNormEff->Clone("hCorrectedYieldNormEff_StatError");

  for(int k = 2; k < numberbins; k++){
    hCorrYield_RelativSyserror->SetBinContent(k, hCorrYield_syserror->GetBinError(k)/(Double_t)hCorrYieldNormal->GetBinContent(k)*100.);
    std::cout << "rel Error = " <<  hCorrYield_RelativSyserror->GetBinContent(k) << '\n';
    hCorrYield_RelativStaterror->SetBinContent(k, hCorrYield_RelativStaterror->GetBinError(k)/(Double_t)hCorrYieldNormal->GetBinContent(k)*100.);
    hCorrYield_RelativStaterror->SetBinError(k,0.);
    hCorrectedYieldNormEff_StatError->SetBinContent(k, hCorrectedYieldNormEff->GetBinError(k)/(Double_t)hCorrectedYieldNormEff->GetBinContent(k)*100.);
    hCorrectedYieldNormEff_StatError->SetBinError(k,0.);
    hCorrYield_RelativSyserror->SetBinError(k,0.);
    hCorrYield_countsyserror->SetBinContent(k, hCorrYield_countsyserror->GetBinError(k)/(Double_t)hCorrYieldNormal->GetBinContent(k)*100.);
    hCorrYield_paramsyserror->SetBinContent(k, hCorrYield_paramsyserror->GetBinError(k)/(Double_t)hCorrYieldNormal->GetBinContent(k)*100.);
    hCorrYield_BGFitsyserror->SetBinContent(k, hCorrYield_BGFitsyserror->GetBinError(k)/(Double_t)hCorrYieldNormal->GetBinContent(k)*100.);
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
  hCorrectedYieldNormEff_StatError->SetLineColor(kBlack);
  hCorrectedYieldNormEff_StatError->SetMarkerColor(kBlack);

  TLegend* leg_stat = new TLegend(0.18,0.7,0.5,0.9);
  SetLegendSettigns(leg_stat, 0.079);
  leg_stat->SetTextFont(43);
  leg_stat->SetTextSize(42);
  leg_stat->AddEntry(hCorrectedYieldNormEff_StatError, "standard method", "l");
  leg_stat->AddEntry(hCorrYield_RelativStaterror, "this method" , "l");
  hCorrectedYieldNormEff_StatError->SetYTitle("rel. stat. uncertainty (%)");
  hCorrectedYieldNormEff_StatError->SetXTitle("#it{p}_{T} (GeV/#it{c})");

  hCorrectedYieldNormEff_StatError->SetLineWidth(3);
  hCorrYield_RelativStaterror->SetLineWidth(3);

  hCorrectedYieldNormEff_StatError->GetXaxis()->SetTitleOffset(1.8);
  hCorrectedYieldNormEff_StatError->GetXaxis()->SetLabelOffset(0.008);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetTitleOffset(0.9);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetLabelOffset(0.008);

  hCorrectedYieldNormEff_StatError->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYield_RelativStaterror->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetRangeUser(-0.09, 9.9);
  hCorrectedYieldNormEff_StatError->GetXaxis()->SetTitleSize(42);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetTitleSize(42);
  hCorrectedYieldNormEff_StatError->GetXaxis()->SetLabelSize(42);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetLabelSize(42);

  hCorrectedYieldNormEff_StatError->GetXaxis()->SetTitleFont(43);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetTitleFont(43);
  hCorrectedYieldNormEff_StatError->GetXaxis()->SetLabelFont(43);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetLabelFont(43);

  hCorrectedYieldNormEff_StatError->GetXaxis()->SetNdivisions(309);

  hCorrectedYieldNormEff_StatError->DrawCopy("AXIS");
  hCorrYield_RelativStaterror->DrawCopy("SAME HIST");
  hCorrectedYieldNormEff_StatError->DrawCopy("SAME HIST");
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
  hCorrYieldNormal->SetMarkerSize(1.5);

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
  leg_yield->SetHeader("parametrization method:");
  leg_yield->AddEntry(hCorrYieldNormal, "templates (normal)" , "lp");
  leg_yield->AddEntry(hCorrYieldBetterBkg3to8, "templates (3 to 8)" , "lp");
  leg_yield->AddEntry(hCorrYieldBetterBkgNN, "templates (next neighbours)" , "lp");
  leg_yield->AddEntry(hCorrectedYieldNormEff, "function", "lp");
  hCorrectedYieldNormEff->SetMarkerSize(1.5);


  hCorrYieldBetterBkg3to8->GetYaxis()->SetTitleOffset(1.7);
  hCorrYieldBetterBkg3to8->GetYaxis()->SetRangeUser(1.e-7-5.e-8,1.e-1);
  hCorrYieldBetterBkg3to8->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYieldBetterBkg3to8->GetXaxis()->SetTitleSize(0.025*3./2.);
  hCorrYieldBetterBkg3to8->GetYaxis()->SetTitleSize(0.025*3./2.);
  hCorrYieldBetterBkg3to8->GetXaxis()->SetLabelSize(0.025*3./2.);
  hCorrYieldBetterBkg3to8->GetYaxis()->SetLabelSize(0.025*3./2.);
  hCorrYieldBetterBkg3to8->SetMarkerColor(kBlue+2);
  hCorrYieldBetterBkg3to8->SetLineColor(kBlue+2);
  hCorrYieldNormal->GetYaxis()->SetTitleOffset(1.7);
  hCorrYieldNormal->GetYaxis()->SetRangeUser(1.e-7-5.e-8,1.e-1);
  hCorrYieldNormal->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYieldNormal->GetXaxis()->SetTitleSize(0.025*3./2.);
  hCorrYieldNormal->GetYaxis()->SetTitleSize(0.025*3./2.);
  hCorrYieldNormal->GetXaxis()->SetLabelSize(0.025*3./2.);
  hCorrYieldNormal->GetYaxis()->SetLabelSize(0.025*3./2.);
  hCorrYieldNormal->SetMarkerColor(kRed+3);
  hCorrYieldNormal->SetLineColor(kRed+3);
  hCorrYieldBetterBkgNN->GetYaxis()->SetTitleOffset(1.7);
  hCorrYieldBetterBkgNN->GetYaxis()->SetRangeUser(1.e-7-5.e-8,1.e-1);
  hCorrYieldBetterBkgNN->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrYieldBetterBkgNN->GetXaxis()->SetTitleSize(0.025*3./2.);
  hCorrYieldBetterBkgNN->GetYaxis()->SetTitleSize(0.025*3./2.);
  hCorrYieldBetterBkgNN->GetXaxis()->SetLabelSize(0.025*3./2.);
  hCorrYieldBetterBkgNN->GetYaxis()->SetLabelSize(0.025*3./2.);
  hCorrYieldBetterBkgNN->SetMarkerColor(kGreen+3);
  hCorrYieldBetterBkgNN->SetLineColor(kGreen+3);
  hCorrectedYieldNormEff->GetYaxis()->SetTitleOffset(1.7);
  hCorrectedYieldNormEff->GetYaxis()->SetRangeUser(1.e-7-5.e-8,1.e-1);
  hCorrectedYieldNormEff->GetXaxis()->SetRangeUser(1.4, 12.0);
  hCorrectedYieldNormEff->GetXaxis()->SetTitleSize(0.025*3./2.);
  hCorrectedYieldNormEff->GetYaxis()->SetTitleSize(0.025*3./2.);
  hCorrectedYieldNormEff->GetXaxis()->SetLabelSize(0.025*3./2.);
  hCorrectedYieldNormEff->GetYaxis()->SetLabelSize(0.025*3./2.);
  hCorrectedYieldNormEff->SetMarkerColor(kBlack);
  hCorrectedYieldNormEff->SetLineColor(kBlack);
  hCorrYield_syserror->SetMarkerColor(kRed+3);
  hCorrYield_syserror->SetLineColor(kRed+3);
  hCorrYield_syserror->SetFillColor(kGray+2);
  hCorrYield_syserror->SetFillStyle(1001);

  hCorrYieldNormal->DrawCopy("AXIS");
  hCorrectedYieldNormEff->Draw("SAME");
  // hCorrYield_syserror->DrawCopy("SAME E2");          // sys Error draw
  hCorrYieldNormal->DrawCopy("SAME");
  hCorrYieldBetterBkg3to8->DrawCopy("SAME");
  hCorrYieldBetterBkgNN->DrawCopy("SAME");
  leg_yield->Draw("SAME");
  canInvMass->Update();
  DrawLabelALICE(0.6, 0.85, 0.035, 0.025*3./2., "");
  pad1InvMass->Update();


  TH1D* hYield_dt_chi2map_corrected_ratio = (TH1D*) hCorrectedYieldNormEff->Clone("hYield_dt_chi2map_corrected_ratio");
  hYield_dt_chi2map_corrected_ratio->Divide(hCorrYieldNormal);
  hYield_dt_chi2map_corrected_ratio->SetLineColor(kRed);
  hYield_dt_chi2map_corrected_ratio->SetMarkerColor(kRed);
  hYield_dt_chi2map_corrected_ratio->SetYTitle("Ratio");

  TH1D* hYield_ratio_NN_to_norm = (TH1D*) hCorrYieldBetterBkgNN->Clone("hYield_dt_chi2map_corrected_ratio");
  hYield_ratio_NN_to_norm->Divide(hCorrYieldNormal);
  hYield_ratio_NN_to_norm->SetLineColor(kGreen+3);
  hYield_ratio_NN_to_norm->SetMarkerColor(kGreen+3);
  hYield_ratio_NN_to_norm->SetYTitle("Ratio");

  TH1D* hYield_ratio_3to8_to_norm = (TH1D*) hCorrYieldBetterBkg3to8->Clone("hYield_dt_chi2map_corrected_ratio");
  hYield_ratio_3to8_to_norm->Divide(hCorrYieldNormal);
  hYield_ratio_3to8_to_norm->SetLineColor(kBlue+2);
  hYield_ratio_3to8_to_norm->SetMarkerColor(kBlue+2);
  hYield_ratio_3to8_to_norm->SetYTitle("Ratio");

  TH1D* hYield_ratio_3to8_to_NN = (TH1D*) hCorrYieldBetterBkg3to8->Clone("hYield_dt_chi2map_corrected_ratio");
  hYield_ratio_3to8_to_NN->Divide(hCorrYieldBetterBkgNN);
  hYield_ratio_3to8_to_NN->SetLineColor(kBlack);
  hYield_ratio_3to8_to_NN->SetMarkerColor(kBlack);
  hYield_ratio_3to8_to_NN->SetYTitle("Ratio");



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
  hYield_dt_chi2map_corrected_ratio->GetYaxis()->SetRangeUser(0.84, 1.17);
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
  hYield_ratio_NN_to_norm->DrawCopy("SAME P");
  hYield_ratio_3to8_to_norm->DrawCopy("SAME P");
  hYield_ratio_3to8_to_NN->DrawCopy("SAME P");
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
  leg_yield2->AddEntry(hCorrYieldNormal, "signal + corr. bkg. temp." , "lp");
  leg_yield2->AddEntry(hCorrectedYieldNormEff, "standard method", "lp");


  hCorrYieldNormal->GetXaxis()->SetTitleOffset(1.4);
  hCorrYieldNormal->GetXaxis()->SetLabelOffset(0.008);
  hCorrYieldNormal->GetYaxis()->SetTitleOffset(3.5);
  hCorrYieldNormal->GetYaxis()->SetLabelOffset(0.008);
  hCorrYieldNormal->GetXaxis()->SetTitleSize(42);
  hCorrYieldNormal->GetYaxis()->SetTitleSize(42);
  hCorrYieldNormal->GetXaxis()->SetLabelSize(42);
  hCorrYieldNormal->GetYaxis()->SetLabelSize(42);

  hCorrYieldNormal->GetXaxis()->SetTitleFont(43);
  hCorrYieldNormal->GetYaxis()->SetTitleFont(43);
  hCorrYieldNormal->GetXaxis()->SetLabelFont(43);
  hCorrYieldNormal->GetYaxis()->SetLabelFont(43);


  hCorrYieldNormal->DrawCopy("AXIS");
  hCorrYield_syserror->DrawCopy("SAME E2");
  hCorrYieldNormal->DrawCopy("SAME");
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

  hCorrectedYieldNormEff_StatError->GetYaxis()->SetNdivisions(505);

  hCorrectedYieldNormEff_StatError->GetXaxis()->SetTitleOffset(3.5);
  hCorrectedYieldNormEff_StatError->GetXaxis()->SetLabelOffset(0.008);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetTitleOffset(3.5);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetLabelOffset(0.008);
  hCorrectedYieldNormEff_StatError->GetXaxis()->SetTitleSize(42);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetTitleSize(42);
  hCorrectedYieldNormEff_StatError->GetXaxis()->SetLabelSize(42);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetLabelSize(42);

  hCorrectedYieldNormEff_StatError->GetXaxis()->SetTitleFont(43);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetTitleFont(43);
  hCorrectedYieldNormEff_StatError->GetXaxis()->SetLabelFont(43);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetLabelFont(43);

  TLegend* leg_stat_yield = new TLegend(0.24,0.7,0.5,0.9);
  SetLegendSettigns(leg_stat_yield, 0.079);
  leg_stat_yield->SetTextFont(43);
  leg_stat_yield->SetTextSize(42);
  leg_stat_yield->AddEntry(hCorrectedYieldNormEff_StatError, "standard method", "l");
  leg_stat_yield->AddEntry(hCorrYield_RelativStaterror, "this method" , "l");

  hCorrectedYieldNormEff_StatError->DrawCopy("AXIS");
  hCorrYield_RelativStaterror->DrawCopy("SAME HIST");
  hCorrectedYieldNormEff_StatError->DrawCopy("SAME HIST");
  leg_stat_yield->Draw("");
  pad3Yield->Update();

  canYield->Update();
  canYield->SaveAs(Form("Systematics/CorrectedYieldCompWithStat." + PicFormat));
  canYield->Clear("D");


  delete leg_yield2;
  delete leg_stat_yield;
  delete leg_stat;

  TString sPath = gDirectory->GetPath();
  if(numberneighbours == 2){
    FStatUnc      = new TFile("FStatUnc.root", "RECREATE");
  }

  else{
    FStatUnc      = new TFile("FStatUnc.root", "UPDATE");
  }

  hCorrYield_RelativStaterror->Write(Form("hCorrYield_RelativStaterror_with%02d_bins", numberneighbours));
  gDirectory->Cd(sPath.Data());

  TCanvas* cStatUnc = new TCanvas("cStatUnc", "", 2*1058,2*319);
  SetCanvasStandardSettings(cStatUnc);
  cStatUnc->cd();
  cStatUnc->SetLeftMargin(0.12);

  TLegend* leg_stat_yield2 = new TLegend(0.20,0.7,0.5,0.9);
  SetLegendSettigns(leg_stat_yield2, 0.079);
  leg_stat_yield2->SetTextFont(43);
  leg_stat_yield2->SetTextSize(42);
  leg_stat_yield2->SetHeader("method:");
  leg_stat_yield2->AddEntry(hCorrectedYieldNormEff_StatError, "function parametrization", "l");
  leg_stat_yield2->AddEntry(hCorrYield_RelativStaterror, "templates parametrization" , "l");

  hCorrectedYieldNormEff_StatError->GetXaxis()->SetTitleOffset(1.0);
  hCorrectedYieldNormEff_StatError->GetYaxis()->SetTitleOffset(0.5);

  hCorrectedYieldNormEff_StatError->DrawCopy("AXIS");
  hCorrYield_RelativStaterror->DrawCopy("SAME HIST");
  hCorrectedYieldNormEff_StatError->DrawCopy("SAME HIST");
  leg_stat_yield2->Draw("");

  cStatUnc->Update();
  cStatUnc->SaveAs(Form("Systematics/StatUncertainty." + PicFormat));
  cStatUnc->Clear();

  delete leg_stat_yield2;





  ////////////////////////////////////////////////////////////////////////////
  // drwaing yields + ratio to Bylikin
  canInvMass->cd();
  pad1InvMass->Draw();
  pad2InvMass->Draw("same");
  pad1InvMass->cd();
  pad1InvMass->SetTickx();
  pad1InvMass->SetTicky();

  pad1InvMass->SetLogy(1);

  TLegend* leg_yield_bylikin = new TLegend(0.2,0.07,0.35,0.3);
  SetLegendSettigns(leg_yield_bylikin, 0.025*3./2.);
  leg_yield_bylikin->SetHeader("parametrization method:");
  leg_yield_bylikin->AddEntry(hCorrYieldNormal, "templates (normal)" , "lp");
  leg_yield_bylikin->AddEntry(hCorrYieldBetterBkg3to8, "templates (3 to 8)" , "lp");
  leg_yield_bylikin->AddEntry(hCorrYieldBetterBkgNN, "templates (next neighbours)" , "lp");
  leg_yield_bylikin->AddEntry(fitBylikin13TeV, "Bylikin", "lp");
  hCorrectedYieldNormEff->SetMarkerSize(1.5);


  hCorrYieldNormal->DrawCopy("AXIS");
  hCorrectedYieldNormEff->DrawCopy("SAME");
  fitBylikin13TeV->Draw("SAME");
  fitBylikin13TeV_3to8->Draw("SAME");
  hCorrYieldNormal->DrawCopy("SAME");
  hCorrYieldBetterBkg3to8->DrawCopy("SAME");
  hCorrYieldBetterBkgNN->DrawCopy("SAME");
  leg_yield_bylikin->Draw("SAME");
  canInvMass->Update();
  DrawLabelALICE(0.6, 0.85, 0.035, 0.025*3./2., "");
  pad1InvMass->Update();


  TH1D* hYield_ratio_func_to_bylikin = (TH1D*) hCorrectedYieldNormEff->Clone("hYield_dt_chi2map_corrected_ratio");
  hYield_ratio_func_to_bylikin->Divide(fitBylikin13TeV);
  hYield_ratio_func_to_bylikin->SetLineColor(kBlack);
  hYield_ratio_func_to_bylikin->SetMarkerColor(kBlack);
  hYield_ratio_func_to_bylikin->SetMarkerStyle(20);
  hYield_ratio_func_to_bylikin->SetMarkerSize(1.5);
  hYield_ratio_func_to_bylikin->SetYTitle("Ratio");

  TH1D* hYield_ratio_norm_to_bylikin = (TH1D*) hCorrYieldNormal->Clone("hYield_dt_chi2map_corrected_ratio");
  hYield_ratio_norm_to_bylikin->Divide(hCorrectedYieldNormEff);
  hYield_ratio_norm_to_bylikin->SetLineColor(kRed);
  hYield_ratio_norm_to_bylikin->SetMarkerColor(kRed);
  hYield_ratio_norm_to_bylikin->SetYTitle("Ratio");

  TH1D* hYield_ratio_NN_to_bylikin = (TH1D*) hCorrYieldBetterBkgNN->Clone("hYield_dt_chi2map_corrected_ratio");
  hYield_ratio_NN_to_bylikin->Divide(hCorrectedYieldNormEff);
  hYield_ratio_NN_to_bylikin->SetLineColor(kGreen+3);
  hYield_ratio_NN_to_bylikin->SetMarkerColor(kGreen+3);
  hYield_ratio_NN_to_bylikin->SetYTitle("Ratio");

  TH1D* hYield_ratio_3to8_to_bylikin = (TH1D*) hCorrYieldBetterBkg3to8->Clone("hYield_dt_chi2map_corrected_ratio");
  hYield_ratio_3to8_to_bylikin->Divide(hCorrectedYieldNormEff);
  hYield_ratio_3to8_to_bylikin->SetLineColor(kBlue+2);
  hYield_ratio_3to8_to_bylikin->SetMarkerColor(kBlue+2);
  hYield_ratio_3to8_to_bylikin->SetYTitle("Ratio");

  Int_t nb = 36;
  Int_t inside = 0;
  for(int y = 1; y <= 36; y++){
    Double_t bc = hYield_ratio_3to8_to_bylikin->GetBinContent(y);
    if(bc < 1.0){
      if(bc + hYield_ratio_3to8_to_bylikin->GetBinError(y) >= 1.0){
        inside++;
      }
    }
    else{
      if(bc - hYield_ratio_3to8_to_bylikin->GetBinError(y) <= 1.0){
        inside++;
      }
    }
  }
std::cout << "percantage of Bins inside of 1sigma range: " << (Double_t)inside/(Double_t)nb << '\n';
std::cout << "number of Bins inside of 1sigma range: " << inside << '\n';

  hYield_ratio_norm_to_bylikin->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hYield_ratio_norm_to_bylikin->GetXaxis()->SetTitleOffset(1.0);
  hYield_ratio_norm_to_bylikin->GetXaxis()->SetLabelOffset(0.008);
  hYield_ratio_norm_to_bylikin->GetYaxis()->SetTitleOffset(0.8);
  hYield_ratio_norm_to_bylikin->GetYaxis()->SetLabelOffset(0.008);

  //changed
  hYield_ratio_norm_to_bylikin->GetXaxis()->SetRangeUser(1.4, 12.0);
  hYield_ratio_norm_to_bylikin->GetYaxis()->SetRangeUser(0.84, 1.17); //(0.84, 1.17)
  hYield_ratio_norm_to_bylikin->GetXaxis()->SetTitleSize(0.025*3./1.);
  hYield_ratio_norm_to_bylikin->GetYaxis()->SetTitleSize(0.025*3./1.);
  hYield_ratio_norm_to_bylikin->GetXaxis()->SetLabelSize(0.025*3./1.);
  hYield_ratio_norm_to_bylikin->GetYaxis()->SetLabelSize(0.025*3./1.);

  hYield_ratio_norm_to_bylikin->GetXaxis()->SetTitleFont(42);
  hYield_ratio_norm_to_bylikin->GetYaxis()->SetTitleFont(42);
  hYield_ratio_norm_to_bylikin->GetXaxis()->SetLabelFont(42);
  hYield_ratio_norm_to_bylikin->GetYaxis()->SetLabelFont(42);

  pad2InvMass->cd();

  hYield_ratio_norm_to_bylikin->DrawCopy("AXIS");
  line_ratio1->Draw("SAME");
  // fBylikinRatio->Draw("SAME");
  hYield_ratio_norm_to_bylikin->DrawCopy("SAME P");
  hYield_ratio_NN_to_bylikin->DrawCopy("SAME P");
  hYield_ratio_3to8_to_bylikin->DrawCopy("SAME P");
  // hYield_ratio_func_to_bylikin->DrawCopy("SAME P");
  pad2InvMass->Update();

  pad2InvMass->SetTickx();
  pad2InvMass->SetTicky();

  pad2InvMass->Update();

  canInvMass->Update();
  canInvMass->SaveAs(Form("Systematics/CorrectedYieldBylikin." + PicFormat));
  canInvMass->Clear("D");

  delete leg_yield_bylikin;

  ////////////////////////////////////////////////////////////////////////////
  // drwaing yields + ratio to Tsallis
  canInvMass->cd();
  pad1InvMass->Draw();
  pad2InvMass->Draw("same");
  pad1InvMass->cd();
  pad1InvMass->SetTickx();
  pad1InvMass->SetTicky();

  pad1InvMass->SetLogy(1);

  TLegend* leg_yield_tsallis = new TLegend(0.2,0.07,0.35,0.3);
  SetLegendSettigns(leg_yield_tsallis, 0.025*3./2.);
  leg_yield_tsallis->SetHeader("parametrization method:");
  leg_yield_tsallis->AddEntry(hCorrYieldNormal, "templates (normal)" , "lp");
  leg_yield_tsallis->AddEntry(hCorrYieldBetterBkg3to8, "templates (3 to 8)" , "lp");
  leg_yield_tsallis->AddEntry(hCorrYieldBetterBkgNN, "templates (next neighbours)" , "lp");
  leg_yield_tsallis->AddEntry(ftsallis13TeV, "Tsallis", "lp");
  hCorrectedYieldNormEff->SetMarkerSize(1.5);


  hCorrYieldNormal->DrawCopy("AXIS");
  ftsallis13TeV->Draw("SAME");
  hCorrYieldNormal->DrawCopy("SAME");
  hCorrYieldBetterBkgNN->DrawCopy("SAME");
  leg_yield_tsallis->Draw("SAME");
  canInvMass->Update();
  DrawLabelALICE(0.6, 0.85, 0.035, 0.025*3./2., "");
  pad1InvMass->Update();


  TH1D* hYield_ratio_func_to_tsallis = (TH1D*) hCorrectedYieldNormEff->Clone("hYield_dt_chi2map_corrected_ratio");
  hYield_ratio_func_to_tsallis->Divide(fitBylikin13TeV);
  hYield_ratio_func_to_tsallis->SetLineColor(kBlack);
  hYield_ratio_func_to_tsallis->SetMarkerColor(kBlack);
  hYield_ratio_func_to_tsallis->SetYTitle("Ratio");

  TH1D* hYield_ratio_norm_to_tsallis = (TH1D*) hCorrYieldNormal->Clone("hYield_dt_chi2map_corrected_ratio");
  hYield_ratio_norm_to_tsallis->Divide(hCorrectedYieldNormEff);
  hYield_ratio_norm_to_tsallis->SetLineColor(kRed);
  hYield_ratio_norm_to_tsallis->SetMarkerColor(kRed);
  hYield_ratio_norm_to_tsallis->SetYTitle("Ratio");

  TH1D* hYield_ratio_NN_to_tsallis = (TH1D*) hCorrYieldBetterBkgNN->Clone("hYield_dt_chi2map_corrected_ratio");
  hYield_ratio_NN_to_tsallis->Divide(hCorrectedYieldNormEff);
  hYield_ratio_NN_to_tsallis->SetLineColor(kGreen+3);
  hYield_ratio_NN_to_tsallis->SetMarkerColor(kGreen+3);
  hYield_ratio_NN_to_tsallis->SetYTitle("Ratio");

  TH1D* hYield_ratio_3to8_to_tsallis = (TH1D*) hCorrYieldBetterBkg3to8->Clone("hYield_dt_chi2map_corrected_ratio");
  hYield_ratio_3to8_to_tsallis->Divide(hCorrectedYieldNormEff);
  hYield_ratio_3to8_to_tsallis->SetLineColor(kBlue+2);
  hYield_ratio_3to8_to_tsallis->SetMarkerColor(kBlue+2);
  hYield_ratio_3to8_to_tsallis->SetYTitle("Ratio");


  hYield_ratio_norm_to_tsallis->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hYield_ratio_norm_to_tsallis->GetXaxis()->SetTitleOffset(1.0);
  hYield_ratio_norm_to_tsallis->GetXaxis()->SetLabelOffset(0.008);
  hYield_ratio_norm_to_tsallis->GetYaxis()->SetTitleOffset(0.8);
  hYield_ratio_norm_to_tsallis->GetYaxis()->SetLabelOffset(0.008);

  hYield_ratio_norm_to_tsallis->GetXaxis()->SetRangeUser(1.4, 12.0);
  //changed
  hYield_ratio_norm_to_tsallis->GetYaxis()->SetRangeUser(0.84, 1.15); //(0.84, 1.15)
  hYield_ratio_norm_to_tsallis->GetXaxis()->SetTitleSize(0.025*3./1.);
  hYield_ratio_norm_to_tsallis->GetYaxis()->SetTitleSize(0.025*3./1.);
  hYield_ratio_norm_to_tsallis->GetXaxis()->SetLabelSize(0.025*3./1.);
  hYield_ratio_norm_to_tsallis->GetYaxis()->SetLabelSize(0.025*3./1.);

  hYield_ratio_norm_to_tsallis->GetXaxis()->SetTitleFont(42);
  hYield_ratio_norm_to_tsallis->GetYaxis()->SetTitleFont(42);
  hYield_ratio_norm_to_tsallis->GetXaxis()->SetLabelFont(42);
  hYield_ratio_norm_to_tsallis->GetYaxis()->SetLabelFont(42);

  pad2InvMass->cd();

  hYield_ratio_norm_to_tsallis->DrawCopy("AXIS");
  line_ratio1->Draw("SAME");
  hYield_ratio_norm_to_tsallis->DrawCopy("SAME P");
  hYield_ratio_NN_to_tsallis->DrawCopy("SAME P");
  hYield_ratio_3to8_to_tsallis->DrawCopy("SAME P");
  // hYield_ratio_func_to_tsallis->DrawCopy("SAME P");
  pad2InvMass->Update();

  pad2InvMass->SetTickx();
  pad2InvMass->SetTicky();

  pad2InvMass->Update();

  canInvMass->Update();
  canInvMass->SaveAs(Form("Systematics/CorrectedYieldTsallis." + PicFormat));
  canInvMass->Clear("D");

  delete ftsallis13TeV;

  delete line_ratio1;
  delete hCorrYieldNormal;
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
  delete hCorrectedYieldNormEff;

}
