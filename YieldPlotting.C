#include "Plotting_Patrick.h"
#include "CommonHeader.h"
//MCTemplatesAnData

// wpsid = which picture should I draw
void YieldPlotting(TString wpsid = "all", TString PicFormat = "pdf", std::string Data = "", TString SaveAppendix = ""){

  TGaxis::SetMaxDigits(3);
  TObjArray* OAhists = new TObjArray();
  TObjArray* OAratios = new TObjArray();

  TFile* BetterBkg3to8File  = SafelyOpenRootfile("OutputFileBetterBkg3to8.root");

  TCanvas* c1                                 = NULL; // Canvas to draw
  TH1D* hCorrYieldME_BetterBkg3to8            = NULL; // Yield from 3to8 Method
  TH1D* hCorrYield_Framework                  = NULL; // Yield from Framework
  TH1D* hCorrYieldME_Ratio_BetterBkg3to8      = NULL; // Ratio me/Framework
  TH1D* hCorrYield_SysError                   = NULL; // Yield with systematics.
  TH1D* hCorrYield_RelativSyserror            = NULL; // rel. systematics
  TH1D* hCorrYieldME_StatError                = NULL; // relative stat. errors
  TH1D* hCorrectedYieldNormEff_StatError      = NULL; // -||- for framework

  hCorrYieldME_BetterBkg3to8        = (TH1D*) BetterBkg3to8File->Get("hYield_dt_chi2map_corrected");
  hCorrYield_Framework              = (TH1D*) BetterBkg3to8File->Get("hCorrectedYieldNormEff");
  hCorrYieldME_Ratio_BetterBkg3to8  = (TH1D*) BetterBkg3to8File->Get("hCorrYieldME_Ratio");
  hCorrYield_SysError               = (TH1D*) BetterBkg3to8File->Get("hCorrYield_SysError");
  hCorrYield_RelativSyserror        = (TH1D*) BetterBkg3to8File->Get("hCorrYield_RelativSyserror");
  hCorrYieldME_StatError            = (TH1D*) BetterBkg3to8File->Get("hCorrYieldME_StatError");
  hCorrectedYieldNormEff_StatError  = (TH1D*) BetterBkg3to8File->Get("hCorrectedYieldNormEff_StatError");



  SetHistogramProperties(hCorrYieldME_BetterBkg3to8,        "pt", strCorrectedYield, 8, 1.4, 14.0);
  SetHistogramProperties(hCorrYield_Framework,              "pt", strCorrectedYield, 5, 1.4, 14.0);
  SetHistogramProperties(hCorrYieldME_Ratio_BetterBkg3to8,  "pt", "#frac{template}{traditional}", 0, 1.4, 14.0);
  SetHistogramProperties(hCorrYield_SysError,               "pt", strCorrectedYield, 8, 1.4, 14.0);
  SetHistogramProperties(hCorrYield_RelativSyserror,        "pt", "rel. systematic uncertainties (%)", 8, 1.4, 14.0);
  SetHistogramProperties(hCorrYieldME_StatError,            "pt", "rel. statistical uncertainties (%)", 8, 1.4, 14.0);
  SetHistogramProperties(hCorrectedYieldNormEff_StatError,  "pt", "rel. statistical uncertainties (%)", 5, 1.4, 14.0);


  /*
    TLegend which displays most important info as header!
   */
  TLegend* legSystem = new TLegend(0.01, 0.94, 0.63, 0.98);
  legSystem->AddEntry((TObject*) 0, "ALICE, pp #sqrt{#it{s}} = 13 TeV, #pi^{0} #rightarrow #gamma#gamma with EMCal", "");

  TLegend* legTemplat = new TLegend(0.2, 0.05, 0.5, 0.3);
  legTemplat->AddEntry((TObject*) 0, "templates:", "");
  legTemplat->AddEntry((TObject*) 0, "PYTHIA 8", "");
  legTemplat->AddEntry((TObject*) 0, "Monash 2013", "");
  legTemplat->AddEntry((TObject*) 0, "GEANT 3", "");

  TLine* OneLine = new TLine(1.4, 1.0, 14.0, 1.0);
  OneLine->SetLineColor(kBlack);
  OneLine->SetLineWidth(3);
  OneLine->SetLineStyle(2);

  /*
    Plotting of the corrected Yield with both uncertainties and Ratio
   */

  TLegend* legYields = new TLegend(0.6, 0.7, 0.8, 0.9);
  legYields->AddEntry((TObject*) 0, "method:", "");
  legYields->AddEntry(hCorrYieldME_BetterBkg3to8, "template", "l");
  legYields->AddEntry(hCorrYield_Framework, "traditional", "l");

  hCorrYield_SysError->GetYaxis()->SetRangeUser(1.5E-7, 4.0E-2);

  OAhists->Add(hCorrYield_SysError); //seems to be buggy
  OAhists->Add(hCorrYieldME_BetterBkg3to8);
  OAhists->Add(hCorrYield_Framework);
  OAhists->Add(legSystem);
  OAhists->Add(legTemplat);
  OAhists->Add(legYields);

  hCorrYieldME_Ratio_BetterBkg3to8->GetYaxis()->SetRangeUser(0.87, 1.148);

  OAratios->Add(hCorrYieldME_Ratio_BetterBkg3to8);
  OAratios->Add(OneLine);

  c1 = NULL;
  c1 = makeCanvas(OAhists, OAratios, "notimeSystematicsLogYThick", 0, 0);

  c1->Update();
  c1->SaveAs(Form("Yields/KorrigierteYields" + SaveAppendix + "." + PicFormat));
  c1->Clear();
  c1 = NULL;
  OAhists->Clear();
  OAratios->Clear();


  /****************************************************************************/


  TLegend* legSystemStat = new TLegend(0.1, 0.94, 0.7, 0.98);
  legSystemStat->AddEntry((TObject*) 0, "ALICE, pp #sqrt{#it{s}} = 13 TeV, #pi^{0} #rightarrow #gamma#gamma with EMCal", "");

  TLegend* legTemplatStat = new TLegend(0.1, 0.68, 0.3, 0.88);
  legTemplatStat->AddEntry((TObject*) 0, "templates:", "");
  legTemplatStat->AddEntry((TObject*) 0, "PYTHIA 8", "");
  legTemplatStat->AddEntry((TObject*) 0, "Monash 2013", "");
  legTemplatStat->AddEntry((TObject*) 0, "GEANT 3", "");


  TLegend* legStat = new TLegend(0.8, 0.22, 0.9, 0.4);
  legStat->AddEntry((TObject*) 0, "method:", "");
  legStat->AddEntry(hCorrYieldME_StatError, "template", "l");
  legStat->AddEntry(hCorrectedYieldNormEff_StatError, "traditional", "l");

  /*
    Plotting of the statistical uncertainties of the Yields
   */

  hCorrectedYieldNormEff_StatError->GetYaxis()->SetRangeUser(0.0, 5.03);

  OAhists->Add(hCorrectedYieldNormEff_StatError);
  OAhists->Add(hCorrYieldME_StatError);
  OAhists->Add(legSystemStat);
  OAhists->Add(legStat);
  OAhists->Add(legTemplatStat);

  c1 = makeCanvas(OAhists, 0, "notimeHorizontalThick", 0, 0);

  c1->Update();
  c1->SaveAs(Form("Yields/StatistischeUnsicherheitVergleich" + SaveAppendix + "." + PicFormat));
  c1->Clear();
  c1 = NULL;
  OAhists->Clear();
  OAratios->Clear();
  /****************************************************************************/

  TLegend* legTemplatSys = new TLegend(0.6, 0.66, 0.9, 0.9);
  legTemplatSys->AddEntry((TObject*) 0, "templates:", "");
  legTemplatSys->AddEntry((TObject*) 0, "PYTHIA 8", "");
  legTemplatSys->AddEntry((TObject*) 0, "Monash 2013", "");
  legTemplatSys->AddEntry((TObject*) 0, "GEANT 3", "");

  /*
    Plotting of the systematic uncertainties of the Yield
   */

  OAhists->Add(hCorrYield_RelativSyserror);
  OAhists->Add(legTemplatSys);
  OAhists->Add(legSystemStat);

  c1 = makeCanvas(OAhists, 0, "notimethickHorizontalLines", 0, 0);

  c1->Update();
  c1->SaveAs(Form("Yields/SystematischeUnsicherheit" + SaveAppendix + "." + PicFormat));
  c1->Clear();
  c1 = NULL;
  OAhists->Clear();
  OAratios->Clear();
  /****************************************************************************/

  /*
    Plotting of the corrected Yield with both uncertainties
   */

  TLegend* legSystemYield = new TLegend(0.1, 0.94, 0.7, 0.98);
  legSystemYield->AddEntry((TObject*) 0, "ALICE, pp #sqrt{#it{s}} = 13 TeV, #pi^{0} #rightarrow #gamma#gamma with EMCal", "");
  hCorrYield_SysError->GetYaxis()->SetRangeUser(1.5E-7, 4.0E-2);
  hCorrYieldME_BetterBkg3to8->GetYaxis()->SetRangeUser(1.5E-7, 4.0E-2);

  OAhists->Add(hCorrYield_SysError);
  OAhists->Add(hCorrYieldME_BetterBkg3to8);
  OAhists->Add(legSystemYield);
  OAhists->Add(legTemplatSys);

  c1 = makeCanvas(OAhists, 0, "notimeSystematicsLogYHorizontal", 0, 0);

  c1->Update();
  c1->SaveAs(Form("Yields/KorrigierterYield" + SaveAppendix + "." + PicFormat));
  c1->Clear();
  c1 = NULL;
  OAhists->Clear();
  OAratios->Clear();
  /****************************************************************************/

  /*
    Plotting of the corrected Yield with only statistical uncertainties
   */

  OAhists->Add(hCorrYieldME_BetterBkg3to8);
  OAhists->Add(legSystemYield);
  OAhists->Add(legTemplatSys);

  c1 = makeCanvas(OAhists, 0, "notimeLogYHorizontal", 0, 0);

  c1->Update();
  c1->SaveAs(Form("Yields/KorrigierterYieldNurStat" + SaveAppendix + "." + PicFormat));
  c1->Clear();
  c1 = NULL;
  OAhists->Clear();
  OAratios->Clear();
  /****************************************************************************/

  delete legSystem;
  delete legTemplat;
  delete OneLine;

}
