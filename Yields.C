#include "Plotting_Patrick.h"
#include "CommonHeader.h"

void Yields(TString PicFormat = "png", TString SaveAppendix = ""){
    TGaxis::SetMaxDigits(3);
    TObjArray* OAhists = new TObjArray();
    TObjArray* OAratios = new TObjArray();

    TFile* BetterBkgNN    = SafelyOpenRootfile("OutputFileBetterBkgNN.root");

    TFile* BetterBkg3to8  = SafelyOpenRootfile("OutputFileBetterBkg3to8.root");

    TFile* Normal         = SafelyOpenRootfile("OutputFileNormal.root");


    TH1D* hCorrYieldME_BetterBkgNN              = NULL;
    TH1D* hCorrYieldME_StatError_BetterBkgNN    = NULL;
    TH1D* hCorrYieldME_BetterBkg3to8            = NULL;
    TH1D* hCorrYieldME_StatError_BetterBkg3to8  = NULL;
    TH1D* hCorrYieldME_BetterBkgPulse           = NULL;
    TH1D* hCorrYieldME_Normal                   = NULL;
    TH1D* hCorrYieldME_StatError_Normal         = NULL;
    TH1D* hCorrYieldME_Ratio_BetterBkgNN        = NULL;
    TH1D* hCorrYieldME_Ratio_BetterBkg3to8      = NULL;
    TH1D* hCorrYieldME_Ratio_BetterBkgPulse     = NULL;
    TH1D* hCorrYieldME_Ratio_Normal             = NULL;
    TH1D* hCorrYield_LeftBkgSmallSysUncertainty = NULL;
    TH1D* hCorrYield_LeftBkgWideSysUncertainty  = NULL;
    TH1D* hCorrectedYieldNormEff_StatError      = NULL;
    TH1D* hCorrectedYieldNormEff                = NULL;
    TH1D* hCorrectedYieldNormEff_Ratio          = NULL;
    TH1D* hPi0_gen                              = NULL;
    TH1D* Chi2_pT_BetterBkgNN                   = NULL;
    TH1D* Chi2_pT_BetterBkg3to8                 = NULL;
    TH1D* Chi2_pT_BetterBkgPulse                = NULL;
    TH1D* Chi2_pT_Normal                        = NULL;

    TH1D* Chi2_pT_Framework                     = NULL;



    TH1D* hYield_dt_chi2map_uncorr              = NULL;
    TH1D* hYield_framework                      = NULL;
    TH1D* hEfficiency                           = NULL;
    TH1D* TrueMesonEffiPt                       = NULL;

    TH1D*hCorrYield_SysError                    = NULL;
    TH1D*hCorrYield_RelativSyserror             = NULL;
    TH1D*hCorrYield_SyserrorRatio               = NULL;

    TH1D* hCorrYield_HigherFitSysUncertainty        = NULL;
    TH1D* hCorrYield_SmallFitSysUncertainty         = NULL;
    TH1D* hCorrYield_HigherIntSysUncertainty        = NULL;
    TH1D* hCorrYield_SmallIntSysUncertainty         = NULL;
    TH1D* hCorrYield_LowerRebinningSysUncertainty   = NULL;
    TH1D* hCorrYield_HigherRebinningSysUncertainty  = NULL;
    TH1D* hCorrYield_NNMethodSysUncertainty         = NULL;
    TH1D* hCorrYield_SingleBkgSysUncertainty        = NULL;

    TLegend* legSystem = new TLegend(0.1, 0.94, 0.7, 0.98);
    legSystem->AddEntry((TObject*) 0, "ALICE, pp bei #sqrt{#it{s}} = 13 TeV, #pi^{0} #rightarrow #gamma#gamma mit EMCal", "");


    hCorrYield_HigherFitSysUncertainty        = (TH1D*) BetterBkg3to8->Get("hCorrYield_HigherFitSysUncertainty");
    SetHistogramProperties(hCorrYield_HigherFitSysUncertainty, "pt", "relative sys. Abweichung (%)", 2, 1.4, 12.);

    hCorrYield_SmallFitSysUncertainty         = (TH1D*) BetterBkg3to8->Get("hCorrYield_SmallFitSysUncertainty");
    SetHistogramProperties(hCorrYield_SmallFitSysUncertainty, "pt", "relative sys. Abweichung (%)", 1, 1.4, 12.);

    hCorrYield_HigherIntSysUncertainty        = (TH1D*) BetterBkg3to8->Get("hCorrYield_HigherIntSysUncertainty");
    SetHistogramProperties(hCorrYield_HigherIntSysUncertainty, "pt", "relative sys. Abweichung (%)", 2, 1.4, 12.);

    hCorrYield_SmallIntSysUncertainty         = (TH1D*) BetterBkg3to8->Get("hCorrYield_SmallIntSysUncertainty");
    SetHistogramProperties(hCorrYield_SmallIntSysUncertainty, "pt", "relative sys. Abweichung (%)", 1, 1.4, 12.);

    hCorrYield_LowerRebinningSysUncertainty   = (TH1D*) BetterBkg3to8->Get("hCorrYield_LowerRebinningSysUncertainty");
    SetHistogramProperties(hCorrYield_LowerRebinningSysUncertainty, "pt", "relative sys. Abweichung (%)", 1, 1.4, 12.);

    hCorrYield_HigherRebinningSysUncertainty  = (TH1D*) BetterBkg3to8->Get("hCorrYield_HigherRebinningSysUncertainty");
    SetHistogramProperties(hCorrYield_HigherRebinningSysUncertainty, "pt", "relative sys. Abweichung (%)", 2, 1.4, 12.);

    hCorrYield_NNMethodSysUncertainty         = (TH1D*) BetterBkg3to8->Get("hCorrYield_NNMethodSysUncertainty");
    SetHistogramProperties(hCorrYield_NNMethodSysUncertainty, "pt", "relative sys. Abweichung (%)", 2, 1.4, 12.);

    hCorrYield_SingleBkgSysUncertainty        = (TH1D*) BetterBkg3to8->Get("hCorrYield_SingleBkgSysUncertainty");
    SetHistogramProperties(hCorrYield_SingleBkgSysUncertainty, "pt", "relative sys. Abweichung (%)", 1, 1.4, 12.);

    hCorrYield_LeftBkgSmallSysUncertainty     = (TH1D*) BetterBkg3to8->Get("hCorrYield_LeftBkgSmallSysUncertainty");
    SetHistogramProperties(hCorrYield_LeftBkgSmallSysUncertainty, "pt", "relative sys. Abweichung (%)", 1, 1.4, 12.);

    hCorrYield_LeftBkgWideSysUncertainty      = (TH1D*) BetterBkg3to8->Get("hCorrYield_LeftBkgWideSysUncertainty");
    SetHistogramProperties(hCorrYield_LeftBkgWideSysUncertainty, "pt", "relative sys. Abweichung (%)", 2, 1.4, 12.);


    TCanvas* c1                                 = NULL;

    /*
    Draw the relative systematic uncertainty for the different Variations
     */
     OAhists->Clear();

     TLine* ZeroLine = new TLine(1.4, 0.0, 12.0, 0.0);
     ZeroLine->SetLineStyle(2);
     ZeroLine->SetLineWidth(3);

     TLegend* legFitSysUncertainty = new TLegend(0.13, 0.7, 0.5, 0.9);
     legFitSysUncertainty->AddEntry((TObject*) 0, "Parametrisierungsbereich:", "");
     legFitSysUncertainty->AddEntry(hCorrYield_HigherFitSysUncertainty, "breite Grenzen", "p");
     legFitSysUncertainty->AddEntry(hCorrYield_SmallFitSysUncertainty, "schmale Grenzen", "p");

     hCorrYield_SmallFitSysUncertainty->GetYaxis()->SetRangeUser(-25.7, 25.7);

     OAhists->Add(hCorrYield_SmallFitSysUncertainty);
     OAhists->Add(hCorrYield_HigherFitSysUncertainty);
     OAhists->Add(ZeroLine);
     OAhists->Add(legSystem);
     OAhists->Add(legFitSysUncertainty);

    //  c1->Clear();
     c1 = NULL;
     c1 = makeCanvas(OAhists, 0, "notimethickHorizontalPoints", 0, 0);
     c1->SaveAs("Yields/YieldsSysUncerFitRange" + SaveAppendix + "." + PicFormat);
     c1->SaveAs("Yields/YieldsSysUncerFitRange" + SaveAppendix + ".root");

     OAhists->Clear();

     /*************************************************************************/

     TLegend* legIntSysUncertainty = new TLegend(0.13, 0.7, 0.5, 0.9);
     legIntSysUncertainty->AddEntry((TObject*) 0, "Z#ddot{a}hlbereich:", "");
     legIntSysUncertainty->AddEntry(hCorrYield_HigherIntSysUncertainty, "breite Grenzen", "p");
     legIntSysUncertainty->AddEntry(hCorrYield_SmallIntSysUncertainty, "schmale Grenzen", "p");

     hCorrYield_SmallIntSysUncertainty->GetYaxis()->SetRangeUser(-25.7, 25.7);

     OAhists->Add(hCorrYield_SmallIntSysUncertainty);
     OAhists->Add(hCorrYield_HigherIntSysUncertainty);
     OAhists->Add(ZeroLine);
     OAhists->Add(legSystem);
     OAhists->Add(legIntSysUncertainty);


     c1->Clear();
     c1 = NULL;
     c1 = makeCanvas(OAhists, 0, "notimethickHorizontalPoints", 0, 0);
     c1->SaveAs("Yields/YieldsSysUncerIntRange" + SaveAppendix + "." + PicFormat);
     c1->SaveAs("Yields/YieldsSysUncerIntRange" + SaveAppendix + ".root");

     OAhists->Clear();

     /*************************************************************************/

     TLegend* legRebinningSysUncertainty = new TLegend(0.13, 0.7, 0.5, 0.9);
     legRebinningSysUncertainty->AddEntry((TObject*) 0, "Intervallbreite #Delta #it{m}_{inv}:", "");
     legRebinningSysUncertainty->AddEntry(hCorrYield_HigherRebinningSysUncertainty, "vergr#ddot{o}#betaert", "p");
     legRebinningSysUncertainty->AddEntry(hCorrYield_LowerRebinningSysUncertainty, "verkleinert", "p");

     hCorrYield_LowerRebinningSysUncertainty->GetYaxis()->SetRangeUser(-25.7, 25.7);

     OAhists->Add(hCorrYield_LowerRebinningSysUncertainty);
     OAhists->Add(hCorrYield_HigherRebinningSysUncertainty);
     OAhists->Add(ZeroLine);
     OAhists->Add(legSystem);
     OAhists->Add(legRebinningSysUncertainty);

     c1->Clear();
     c1 = NULL;
     c1 = makeCanvas(OAhists, 0, "notimethickHorizontalPoints", 0, 0);
     c1->SaveAs("Yields/YieldsSysUncerRebinning" + SaveAppendix + "." + PicFormat);
     c1->SaveAs("Yields/YieldsSysUncerRebinning" + SaveAppendix + ".root");

     OAhists->Clear();

     /*************************************************************************/

     TLegend* legBkgVariationSysUncertainty = new TLegend(0.13, 0.7, 0.5, 0.9);
     legBkgVariationSysUncertainty->AddEntry((TObject*) 0, "templates des korrelierten Untergrund:", "");
     legBkgVariationSysUncertainty->AddEntry(hCorrYield_NNMethodSysUncertainty, "n#ddot{a}chste Nachbarn", "p");
     legBkgVariationSysUncertainty->AddEntry(hCorrYield_SingleBkgSysUncertainty, "einzelne", "p");

     hCorrYield_SingleBkgSysUncertainty->GetYaxis()->SetRangeUser(-25.7, 25.7);

     OAhists->Add(hCorrYield_SingleBkgSysUncertainty);
     OAhists->Add(hCorrYield_NNMethodSysUncertainty);
     OAhists->Add(ZeroLine);
     OAhists->Add(legSystem);
     OAhists->Add(legBkgVariationSysUncertainty);


     c1->Clear();
     c1 = NULL;
     c1 = makeCanvas(OAhists, 0, "notimethickHorizontalPoints", 0, 0);
     c1->SaveAs("Yields/YieldsSysUncerBkgVariation" + SaveAppendix + "." + PicFormat);
     c1->SaveAs("Yields/YieldsSysUncerBkgVariation" + SaveAppendix + ".root");

     OAhists->Clear();

     /*************************************************************************/

     TLegend* legLeftBkgVariationSysUncertainty = new TLegend(0.13, 0.7, 0.5, 0.9);
     legLeftBkgVariationSysUncertainty->AddEntry((TObject*) 0, "Skalierungsbereich f#ddot{u}r:", "");
     legLeftBkgVariationSysUncertainty->AddEntry((TObject*) 0, "den unkorrelierten Untergrund:", "");
     legLeftBkgVariationSysUncertainty->AddEntry(hCorrYield_LeftBkgWideSysUncertainty, "0.225 #leq #it{m}_{inv} /(GeV/#it{c}^{2}) < 0.3", "p");
     legLeftBkgVariationSysUncertainty->AddEntry(hCorrYield_LeftBkgSmallSysUncertainty, "0.19 #leq #it{m}_{inv} /(GeV/#it{c}^{2}) < 0.265", "p");

     hCorrYield_LeftBkgSmallSysUncertainty->GetYaxis()->SetRangeUser(-25.7, 25.7);

     OAhists->Add(hCorrYield_LeftBkgSmallSysUncertainty);
     OAhists->Add(hCorrYield_LeftBkgWideSysUncertainty);
     OAhists->Add(ZeroLine);
     OAhists->Add(legSystem);
     OAhists->Add(legLeftBkgVariationSysUncertainty);


     c1->Clear();
     c1 = NULL;
     c1 = makeCanvas(OAhists, 0, "notimethickHorizontalPoints", 0, 0);
     c1->SaveAs("Yields/YieldsSysUncerUncorrBkgVariation" + SaveAppendix + "." + PicFormat);
     c1->SaveAs("Yields/YieldsSysUncerUncorrBkgVariation" + SaveAppendix + ".root");

     OAhists->Clear();

     /**************************************************************************
     ***************************************************************************
     **************************************************************************/

     delete legFitSysUncertainty;
     delete legIntSysUncertainty;
     delete legRebinningSysUncertainty;
     delete legBkgVariationSysUncertainty;
     delete legSystem;


     delete OAhists;
     delete OAratios;

}
