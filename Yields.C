#include "Plotting_Patrick.h"
#include "CommonHeader.h"

void Yields(TString PicFormat = "png", TString SaveAppendix = ""){
    TGaxis::SetMaxDigits(3);
    auto OAhists = new TObjArray();
    auto OAratios = new TObjArray();

    TFile* BetterBkgNN    = SafelyOpenRootfile("OutputFileBetterBkgNN.root");

    TFile* BetterBkg3to8  = SafelyOpenRootfile("OutputFileBetterBkg3to8.root");

    TFile* BetterBkgPulse = SafelyOpenRootfile("OutputFileBetterBkgPulse.root");

    TFile* Normal         = SafelyOpenRootfile("OutputFileNormal.root");

    // TFile* NormalWC       = SafelyOpenRootfile("OutputFileNormalWithConstraint.root");

    TFile* OneTemplate    = SafelyOpenRootfile("OutputFileOneTemplate.root");

    TH1D* hCorrYieldME_BetterBkgNN              = NULL;
    TH1D* hCorrYieldME_StatError_BetterBkgNN    = NULL;
    TH1D* hCorrYieldME_BetterBkg3to8            = NULL;
    TH1D* hCorrYieldME_StatError_BetterBkg3to8  = NULL;
    TH1D* hCorrYieldME_BetterBkgPulse           = NULL;
    TH1D* hCorrYieldME_StatError_BetterBkgPulse = NULL;
    TH1D* hCorrYieldME_Normal                   = NULL;
    TH1D* hCorrYieldME_StatError_Normal         = NULL;
    TH1D* hCorrYieldME_Ratio_BetterBkgNN        = NULL;
    TH1D* hCorrYieldME_Ratio_BetterBkg3to8      = NULL;
    TH1D* hCorrYieldME_Ratio_BetterBkgPulse     = NULL;
    TH1D* hCorrYieldME_Ratio_Normal             = NULL;
    TH1D* hCorrectedYieldNormEff_StatError      = NULL;
    TH1D* hCorrectedYieldNormEff                = NULL;
    TH1D* hCorrectedYieldNormEff_Ratio          = NULL;
    TH1D* hPi0_gen                              = NULL;
    // TH1D* hCorrYieldME_NormalWC                 = NULL;
    // TH1D* hCorrYieldME_StatError_NormalWC       = NULL;
    // TH1D* hCorrYieldME_Ratio_NormalWC           = NULL;
    TH1D* Chi2_pT_BetterBkgNN                   = NULL;
    TH1D* Chi2_pT_BetterBkg3to8                 = NULL;
    TH1D* Chi2_pT_BetterBkgPulse                = NULL;
    TH1D* Chi2_pT_Normal                        = NULL;
    // TH1D* Chi2_pT_NormalWC                      = NULL;
    TH1D* Chi2_pT_Framework                     = NULL;

    TH1D* hCorrYieldME_OneTemplate              = NULL;
    TH1D* hCorrYieldME_StatError_OneTemplate    = NULL;
    TH1D* hCorrYieldME_Ratio_OneTemplate        = NULL;
    TH1D* Chi2_pT_OneTemplate                   = NULL;


    TH1D* hYield_dt_chi2map_uncorr              = NULL;
    TH1D* hYield_framework                      = NULL;
    TH1D* hEfficiency                           = NULL;
    TH1D* TrueMesonEffiPt                       = NULL;

    TH1D*hCorrYield_SysError                    = NULL;
    TH1D*hCorrYield_RelativSyserror             = NULL;
    TH1D*hCorrYield_SyserrorRatio               = NULL;

    TCanvas* c1                                 = NULL;

    /**
     * Ratio of uncorrected yields and efficiency
     */
    hYield_dt_chi2map_uncorr              = (TH1D*) OneTemplate->Get("hYield_dt_chi2map_uncorr");
    SetHistogramProperties(hYield_dt_chi2map_uncorr, "pt", rawyield, 1, 1.4, 12.);

    hYield_framework                      = (TH1D*) OneTemplate->Get("hYield_framework");
    SetHistogramProperties(hYield_framework, "pt", rawyield, 1, 1.4, 12.);

    hEfficiency                           = (TH1D*) OneTemplate->Get("hEfficiency");
    SetHistogramProperties(hEfficiency, "pt", "Effizienz", 4, 1.4, 12.);

    TrueMesonEffiPt                       = (TH1D*) OneTemplate->Get("TrueMesonEffiPt");
    SetHistogramProperties(TrueMesonEffiPt, "pt", "Effizienz", 4, 1.4, 12.);

    TH1D* UncorrYieldRatio = (TH1D*) hYield_dt_chi2map_uncorr->Clone("UncorrYieldRatio");
    UncorrYieldRatio->Divide(hYield_dt_chi2map_uncorr, hYield_framework, 1, 1);

    TH1D* EffiRatio = (TH1D*) hEfficiency->Clone("EffiRatio");
    EffiRatio->Divide(hEfficiency, TrueMesonEffiPt, 1, 1);

    SetHistogramProperties(UncorrYieldRatio, "pt", "Ratio", 4, 1.4, 12.);
    SetHistogramProperties(EffiRatio, "pt", "Ratio", 1, 1.4, 12.);

    TLegend* legUncorrYieldAndEffiRatio = new TLegend(0.6, 0.7, 0.9, 0.9);
    legUncorrYieldAndEffiRatio->AddEntry(UncorrYieldRatio, "Unkorr. Yield Ratio", "l");
    legUncorrYieldAndEffiRatio->AddEntry(EffiRatio, "Effizienzratio", "l");


    hPi0_gen                              = (TH1D*) Normal->Get("MC_Meson_genPt");
    SetHistogramProperties(hCorrYieldME_BetterBkgNN, "pt", strCorrectedYield, 5, 1.4, 12.);

    hCorrYieldME_BetterBkgNN              = (TH1D*) BetterBkgNN->Get("hYield_dt_chi2map_corrected");
    SetHistogramProperties(hCorrYieldME_BetterBkgNN, "pt", strCorrectedYield, 1, 1.4, 12.);

    hCorrYieldME_StatError_BetterBkgNN    = (TH1D*) BetterBkgNN->Get("hCorrYieldME_StatError");
    SetHistogramProperties(hCorrYieldME_StatError_BetterBkgNN, "pt", StatUn_Str, 1, 1.4, 12.);

    hCorrYieldME_BetterBkg3to8            = (TH1D*) BetterBkg3to8->Get("hYield_dt_chi2map_corrected");
    SetHistogramProperties(hCorrYieldME_BetterBkg3to8, "pt", strCorrectedYield, 2, 1.4, 12.);

    hCorrYield_SysError                   = (TH1D*) BetterBkg3to8->Get("hCorrYield_SysError");
    SetHistogramProperties(hCorrYield_SysError, "pt", strCorrectedYield, 2, 1.4, 12.);

    hCorrYield_RelativSyserror            = (TH1D*) BetterBkg3to8->Get("hCorrYield_RelativSyserror");
    SetHistogramProperties(hCorrYield_RelativSyserror, "pt", "relative sys. Unsicherheit (%)", 2, 1.4, 12.);

    hCorrYieldME_StatError_BetterBkg3to8  = (TH1D*) BetterBkg3to8->Get("hCorrYieldME_StatError");
    SetHistogramProperties(hCorrYieldME_StatError_BetterBkg3to8, "pt", StatUn_Str, 2, 1.4, 12.);

    hCorrYieldME_BetterBkgPulse           = (TH1D*) BetterBkgPulse->Get("hYield_dt_chi2map_corrected");
    SetHistogramProperties(hCorrYieldME_BetterBkgPulse, "pt", strCorrectedYield, 0, 1.4, 12.);

    hCorrYieldME_StatError_BetterBkgPulse = (TH1D*) BetterBkgPulse->Get("hCorrYieldME_StatError");
    SetHistogramProperties(hCorrYieldME_StatError_BetterBkgPulse, "pt", StatUn_Str, 0, 1.4, 12.);

    hCorrYieldME_Normal                   = (TH1D*) Normal->Get("hYield_dt_chi2map_corrected");
    SetHistogramProperties(hCorrYieldME_Normal, "pt", strCorrectedYield, 4, 1.4, 12.);

    hCorrYieldME_StatError_Normal         = (TH1D*) Normal->Get("hCorrYieldME_StatError");
    SetHistogramProperties(hCorrYieldME_StatError_Normal, "pt", StatUn_Str, 4, 1.4, 12.);

    hCorrYieldME_Ratio_BetterBkgNN        = (TH1D*) BetterBkgNN->Get("hCorrYieldME_Ratio");
    SetHistogramProperties(hCorrYieldME_Ratio_BetterBkgNN, "pt", "Ratio", 1, 1.4, 12.);

    hCorrYieldME_Ratio_BetterBkg3to8      = (TH1D*) BetterBkg3to8->Get("hCorrYieldME_Ratio");
    SetHistogramProperties(hCorrYieldME_Ratio_BetterBkg3to8, "pt", "Ratio", 2, 1.4, 12.);

    hCorrYield_SyserrorRatio              = (TH1D*) BetterBkg3to8->Get("hCorrYield_SyserrorRatio");
    SetHistogramProperties(hCorrYield_SyserrorRatio, "pt", "Ratio", 2, 1.4, 12.);

    hCorrYieldME_Ratio_BetterBkgPulse     = (TH1D*) BetterBkgPulse->Get("hCorrYieldME_Ratio");
    SetHistogramProperties(hCorrYieldME_Ratio_BetterBkgPulse, "pt", "Ratio", 0, 1.4, 12.);

    hCorrYieldME_Ratio_Normal             = (TH1D*) Normal->Get("hCorrYieldME_Ratio");
    SetHistogramProperties(hCorrYieldME_Ratio_Normal, "pt", "Ratio", 4, 1.4, 12.);

    hCorrectedYieldNormEff_StatError      = (TH1D*) Normal->Get("hCorrectedYieldNormEff_StatError");
    SetHistogramProperties(hCorrectedYieldNormEff_StatError, "pt", StatUn_Str, 3, 1.4, 12.);

    hCorrectedYieldNormEff                = (TH1D*) Normal->Get("hCorrectedYieldNormEff");
    SetHistogramProperties(hCorrectedYieldNormEff, "pt", strCorrectedYield, 3, 1.4, 12.);

    hCorrectedYieldNormEff_Ratio          = (TH1D*) Normal->Get("hCorrectedYieldNormEff_Ratio");
    SetHistogramProperties(hCorrectedYieldNormEff_Ratio, "pt", "Ratio", 3, 1.4, 12.);

    // hCorrYieldME_NormalWC                 = (TH1D*) NormalWC->Get("hYield_dt_chi2map_corrected");
    // SetHistogramProperties(hCorrYieldME_NormalWC, "pt", strCorrectedYield, 3, 1.4, 12.);

    // hCorrYieldME_StatError_NormalWC       = (TH1D*) NormalWC->Get("hCorrYieldME_StatError");
    // SetHistogramProperties(hCorrYieldME_StatError_NormalWC, "pt", "Ratio", 3, 1.4, 12.);

    // hCorrYieldME_Ratio_NormalWC           = (TH1D*) NormalWC->Get("hCorrYieldME_Ratio");
    // SetHistogramProperties(hCorrYieldME_Ratio_NormalWC, "pt", "Ratio", 3, 1.4, 12.);


    Chi2_pT_BetterBkgNN    = (TH1D*) BetterBkgNN->Get("hChi2Map_Chi2_pT");
    SetHistogramProperties(Chi2_pT_BetterBkgNN, "pt", "#chi^{2}/ndf", 1, 1.4, 12.);

    Chi2_pT_BetterBkg3to8  = (TH1D*) BetterBkg3to8->Get("hChi2Map_Chi2_pT");
    SetHistogramProperties(Chi2_pT_BetterBkg3to8, "pt", "#chi^{2}/ndf", 2, 1.4, 12.);

    Chi2_pT_BetterBkgPulse = (TH1D*) BetterBkgPulse->Get("hChi2Map_Chi2_pT");
    SetHistogramProperties(Chi2_pT_BetterBkgPulse, "pt", "#chi^{2}/ndf", 0, 1.4, 12.);

    Chi2_pT_Normal         = (TH1D*) Normal->Get("hChi2Map_Chi2_pT");
    SetHistogramProperties(Chi2_pT_Normal, "pt", "#chi^{2}/ndf", 4, 1.4, 12.);

    // Chi2_pT_NormalWC       = (TH1D*) NormalWC->Get("hChi2Map_Chi2_pT");
    // SetHistogramProperties(Chi2_pT_NormalWC, "pt", "#chi^{2}/ndf", 3, 1.4, 12.);

    Chi2_pT_Framework      = (TH1D*) Normal->Get("histoChi2_0");
    SetHistogramProperties(Chi2_pT_Framework, "pt", "#chi^{2}/ndf", 5, 1.4, 12.);

    hCorrYieldME_OneTemplate = (TH1D*) OneTemplate->Get("hYield_dt_chi2map_corrected");
    SetHistogramProperties(hCorrYieldME_OneTemplate, "pt", strCorrectedYield, 6, 1.4, 12.);

    hCorrYieldME_StatError_OneTemplate = (TH1D*) OneTemplate->Get("hCorrYieldME_StatError");
    SetHistogramProperties(hCorrYieldME_StatError_OneTemplate, "pt", StatUn_Str, 6, 1.4, 12.);

    hCorrYieldME_Ratio_OneTemplate = (TH1D*) OneTemplate->Get("hCorrYieldME_Ratio");
    SetHistogramProperties(hCorrYieldME_Ratio_OneTemplate, "pt", "Ratio", 6, 1.4, 12.);

    Chi2_pT_OneTemplate = (TH1D*) OneTemplate->Get("hChi2Map_Chi2_pT");
    SetHistogramProperties(Chi2_pT_OneTemplate, "pt", "#chi^{2}/ndf", 6, 1.4, 12.);





    TLegend* legYields = new TLegend(0.4, 0.5, 0.8, 0.9);
    // legYields->SetHeader("Methode:");
    // legYields->AddEntry(hPi0_gen, "normiertes gen. Spektrum", "p");
    legYields->AddEntry(hCorrectedYieldNormEff, "Funktionsparametrisierung", "p");
    legYields->AddEntry(hCorrYieldME_Normal, "Templates (normal)", "p");
    legYields->AddEntry(hCorrYieldME_OneTemplate, "Template", "p");
    legYields->AddEntry(hCorrYieldME_BetterBkg3to8, "Templates (3 bis 8)", "p");
    legYields->AddEntry(hCorrYieldME_BetterBkgPulse, "Templates (Pulsefunktion)", "p");
    legYields->AddEntry(hCorrYieldME_BetterBkgNN, "Templates (naechste Nachbarn)", "p");

    TLegend* legYieldNN = new TLegend(0.22, 0.01, 0.47, 0.25);
    // legYieldNN->SetHeader("Methode:");
    // legYieldNN->AddEntry(hPi0_gen, "normiertes gen. Spektrum", "p");
    legYieldNN->AddEntry(hCorrectedYieldNormEff, "Funktionsparametrisierung", "p");
    legYieldNN->AddEntry(hCorrYieldME_BetterBkgNN, "Templates (naechste Nachbarn)", "p");

    TLegend* legYield3to8 = new TLegend(0.22, 0.01, 0.47, 0.25);
    // legYield3to8->SetHeader("Methode:");
    // legYield3to8->AddEntry(hPi0_gen, "normiertes gen. Spektrum", "p");
    legYield3to8->AddEntry(hCorrectedYieldNormEff, "Funktionsparametrisierung", "p");
    legYield3to8->AddEntry(hCorrYieldME_BetterBkg3to8, "Templates (3 bis 8)", "p");

    TLegend* legYieldPulse = new TLegend(0.22, 0.01, 0.47, 0.25);
    // legYieldPulse->SetHeader("Methode:");
    // legYieldPulse->AddEntry(hPi0_gen, "normiertes gen. Spektrum", "p");
    legYieldPulse->AddEntry(hCorrectedYieldNormEff, "Funktionsparametrisierung", "p");
    legYieldPulse->AddEntry(hCorrYieldME_BetterBkgPulse, "Templates (Pulsefunktion)", "p");

    TLegend* legYieldNormal = new TLegend(0.22, 0.01, 0.47, 0.25);
    // legYieldNormal->SetHeader("Methode:");
    // legYieldNormal->AddEntry(hPi0_gen, "normiertes gen. Spektrum", "p");
    legYieldNormal->AddEntry(hCorrectedYieldNormEff, "Funktionsparametrisierung", "p");
    legYieldNormal->AddEntry(hCorrYieldME_Normal, "Templates (normal)", "p");

    TLegend* legYieldOneTemplate = new TLegend(0.22, 0.01, 0.47, 0.25);
    // legYieldOneTemplate->SetHeader("Methode:");
    // legYieldOneTemplate->AddEntry(hPi0_gen, "normiertes gen. Spektrum", "p");
    legYieldOneTemplate->AddEntry(hCorrectedYieldNormEff, "Funktionsparametrisierung", "p");
    legYieldOneTemplate->AddEntry(hCorrYieldME_OneTemplate, " einzelnes Template", "p");

    TLegend* legYieldStatUncer= new TLegend(0.5, 0.2, 0.8, 0.4);
    // legYieldStatUncer->SetHeader("Methode:");
    legYieldStatUncer->AddEntry(hCorrectedYieldNormEff_StatError, "Funktionsparametrisierung", "l");
    legYieldStatUncer->AddEntry(hCorrYieldME_StatError_Normal, "Templates (normal)", "l");
    legYieldStatUncer->AddEntry(hCorrYieldME_StatError_OneTemplate, "Template", "l");
    legYieldStatUncer->AddEntry(hCorrYieldME_StatError_BetterBkg3to8, "Templates (3 bis 8)", "l");
    legYieldStatUncer->AddEntry(hCorrYieldME_StatError_BetterBkgPulse, "Templates (Pulsefunktion)", "l");
    legYieldStatUncer->AddEntry(hCorrYieldME_StatError_BetterBkgNN, "Templates (naechste Nachbarn)", "l");

    TLine* lOne = new TLine(1.4, 1.0, 12., 1.0);

    OAhists->Clear();
    OAratios->Clear();

    OAhists->Add(hCorrYieldME_Normal);
    OAhists->Add(hCorrYieldME_OneTemplate);
    OAhists->Add(hCorrYieldME_BetterBkg3to8);
    OAhists->Add(hCorrYieldME_BetterBkgPulse);
    OAhists->Add(hCorrYieldME_BetterBkgNN);
    OAhists->Add(hCorrectedYieldNormEff);
    // OAhists->Add(hPi0_gen);
    OAhists->Add(legYields);

    OAratios->Add(hCorrYieldME_Ratio_Normal);
    OAratios->Add(hCorrYieldME_Ratio_OneTemplate);
    OAratios->Add(hCorrYieldME_Ratio_BetterBkg3to8);
    OAratios->Add(hCorrYieldME_Ratio_BetterBkgPulse);
    OAratios->Add(hCorrYieldME_Ratio_BetterBkgNN);
    // OAratios->Add(hCorrectedYieldNormEff_Ratio);
    OAratios->Add(lOne);


    c1 = NULL;
    c1 = makeCanvas(OAhists, OAratios, "notimethicklogY", 0, 0);
    c1->SaveAs("Yields/Yields" + SaveAppendix + "." + PicFormat);
    c1->SaveAs("Yields/Yields" + SaveAppendix + ".root");

    /**
     * Plotting NN method Yield in comp. with framework yield
     */
    OAhists->Clear();
    OAratios->Clear();

    OAhists->Add(hCorrYieldME_BetterBkgNN);
    OAhists->Add(hCorrectedYieldNormEff);
    OAhists->Add(legYieldNN);

    OAratios->Add(hCorrYieldME_Ratio_BetterBkgNN);
    OAratios->Add(lOne);

    c1->Clear();
    c1 = NULL;
    c1 = makeCanvas(OAhists, OAratios, "notimethicklogY", 0, 0);
    c1->SaveAs("Yields/YieldNN" + SaveAppendix + "." + PicFormat);
    c1->SaveAs("Yields/YieldNN" + SaveAppendix + ".root");

    /**
     * Plotting of 3 to 8 method Yield in comp. with framework yield
     */
    OAhists->Clear();
    OAratios->Clear();

    OAhists->Add(hCorrYieldME_BetterBkg3to8);
    OAhists->Add(hCorrectedYieldNormEff);
    OAhists->Add(legYield3to8);

    OAratios->Add(hCorrYieldME_Ratio_BetterBkg3to8);
    OAratios->Add(lOne);

    c1->Clear();
    c1 = NULL;
    c1 = makeCanvas(OAhists, OAratios, "notimethicklogY", 0, 0);
    c1->SaveAs("Yields/Yield3to8" + SaveAppendix + "." + PicFormat);
    c1->SaveAs("Yields/Yield3to8" + SaveAppendix + ".root");

    /**
     * Plotting of Pulse method Yield in comp. with framework yield
     */
    OAhists->Clear();
    OAratios->Clear();

    OAhists->Add(hCorrYieldME_BetterBkgPulse);
    OAhists->Add(hCorrectedYieldNormEff);
    OAhists->Add(legYieldPulse);

    OAratios->Add(hCorrYieldME_Ratio_BetterBkgPulse);
    OAratios->Add(lOne);

    c1->Clear();
    c1 = NULL;
    c1 = makeCanvas(OAhists, OAratios, "notimethicklogY", 0, 0);
    c1->SaveAs("Yields/YieldPulse" + SaveAppendix + "." + PicFormat);
    c1->SaveAs("Yields/YieldPulse" + SaveAppendix + ".root");


    /**
     * Plotting of Normal method Yield in comp. with framework yield
     */
    OAhists->Clear();
    OAratios->Clear();

    OAhists->Add(hCorrYieldME_Normal);
    OAhists->Add(hCorrectedYieldNormEff);
    OAhists->Add(legYieldNormal);

    OAratios->Add(hCorrYieldME_Ratio_Normal);
    OAratios->Add(lOne);

    c1->Clear();
    c1 = NULL;
    c1 = makeCanvas(OAhists, OAratios, "notimethicklogY", 0, 0);
    c1->SaveAs("Yields/YieldNormal" + SaveAppendix + "." + PicFormat);
    c1->SaveAs("Yields/YieldNormal" + SaveAppendix + ".root");

    /**
     * Plotting of OneTemplate method Yield in comp. with framework yield
     */
    OAhists->Clear();
    OAratios->Clear();

    OAhists->Add(hCorrYieldME_OneTemplate);
    OAhists->Add(hCorrectedYieldNormEff);
    OAhists->Add(legYieldOneTemplate);

    OAratios->Add(hCorrYieldME_Ratio_OneTemplate);
    OAratios->Add(lOne);

    c1->Clear();
    c1 = NULL;
    c1 = makeCanvas(OAhists, OAratios, "notimethicklogY", 0, 0);
    c1->SaveAs("Yields/YieldOneTemplate" + SaveAppendix + "." + PicFormat);
    c1->SaveAs("Yields/YieldOneTemplate" + SaveAppendix + ".root");

    /**
     * Plotting relative stat. uncertainty
     */
    OAhists->Clear();
    OAratios->Clear();

    // hCorrectedYieldNormEff_StatError->GetYaxis()->SetRangeUser(0.0, 6.05);

    OAhists->Add(hCorrectedYieldNormEff_StatError);
    OAhists->Add(hCorrYieldME_StatError_BetterBkg3to8);
    OAhists->Add(hCorrYieldME_StatError_BetterBkgPulse);
    OAhists->Add(hCorrYieldME_StatError_Normal);
    OAhists->Add(hCorrYieldME_StatError_OneTemplate);
    OAhists->Add(hCorrYieldME_StatError_BetterBkgNN);
    OAhists->Add(legYieldStatUncer);

    c1->Clear();
    c1 = NULL;
    c1 = makeCanvas(OAhists, 0, "notimethickHorizontal", 0, 0);
    c1->SaveAs("Yields/YieldsStatUncer" + SaveAppendix + "." + PicFormat);
    c1->SaveAs("Yields/YieldsStatUncer" + SaveAppendix + ".root");

    OAhists->Clear();
    OAratios->Clear();

    TLegend* legChi2 = new TLegend(0.15, 0.7, 0.4, 0.95);
    legChi2->SetHeader("Methode:");
    legChi2->AddEntry(Chi2_pT_Framework, "Funktionsparametrisierung", "l");
    legChi2->AddEntry(Chi2_pT_Normal, "Templates (normal)", "l");
    legChi2->AddEntry(Chi2_pT_OneTemplate, "Template", "l");
    legChi2->AddEntry(Chi2_pT_BetterBkg3to8, "Templates (3 bis 8)", "l");
    legChi2->AddEntry(Chi2_pT_BetterBkgPulse, "Templates (Pulsefunktion)", "l");
    legChi2->AddEntry(Chi2_pT_BetterBkgNN, "Templates (naechste Nachbarn)", "l");

    OAhists->Add(Chi2_pT_Framework);
    OAhists->Add(Chi2_pT_Normal);
    OAhists->Add(Chi2_pT_OneTemplate);
    OAhists->Add(Chi2_pT_BetterBkg3to8);
    OAhists->Add(Chi2_pT_BetterBkgPulse);
    OAhists->Add(Chi2_pT_BetterBkgNN);
    OAhists->Add(legChi2);

    c1->Clear();
    c1 = NULL;
    c1 = makeCanvas(OAhists, 0, "notimethickHorizontalLines", 0, 0);
    c1->SaveAs("Yields/Chi2Comp" + SaveAppendix + "." + PicFormat);
    c1->SaveAs("Yields/Chi2Comp" + SaveAppendix + ".root");

    OAhists->Clear();
    OAratios->Clear();

    OAhists->Add(UncorrYieldRatio);
    OAhists->Add(EffiRatio);
    OAhists->Add(legUncorrYieldAndEffiRatio);

    c1->Clear();
    c1 = NULL;
    c1 = makeCanvas(OAhists, 0, "notimethickHorizontalLines", 0, 0);
    c1->SaveAs("Yields/YieldAndEffiRatio" + SaveAppendix + "." + PicFormat);
    c1->SaveAs("Yields/YieldAndEffiRatio" + SaveAppendix + ".root");

    OAhists->Clear();
    OAratios->Clear();

    OAhists->Add(hCorrYield_SysError);
    OAhists->Add(hCorrYieldME_BetterBkg3to8);
    OAhists->Add(hCorrectedYieldNormEff);
    OAhists->Add(legYield3to8);

    OAratios->Add(hCorrYield_SyserrorRatio);
    OAratios->Add(hCorrYieldME_Ratio_BetterBkg3to8);
    OAratios->Add(lOne);

    c1->Clear();
    c1 = NULL;
    c1 = makeCanvas(OAhists, OAratios, "notimelogYSystematics", 0, 0);
    c1->SaveAs("Yields/Yield3to8WithSys" + SaveAppendix + "." + PicFormat);
    c1->SaveAs("Yields/Yield3to8WithSys" + SaveAppendix + ".root");

    OAhists->Clear();
    OAratios->Clear();


    OAhists->Add(hCorrYield_RelativSyserror);

    c1->Clear();
    c1 = NULL;
    c1 = makeCanvas(OAhists, 0, "notimethickHorizontal", 0, 0);
    c1->SaveAs("Yields/YieldsSysUncer" + SaveAppendix + "." + PicFormat);
    c1->SaveAs("Yields/YieldsSysUncer" + SaveAppendix + ".root");

    OAhists->Clear();
    OAratios->Clear();





    delete OAhists;
    delete OAratios;

}
