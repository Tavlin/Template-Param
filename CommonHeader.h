#include "TLatex.h"
#include "stddef.h"
#include "TLegend.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TList.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TTree.h"
#include "TGraph.h"
#include "TRandom.h"
#include "TFitResult.h"
#include "TROOT.h"
#include "TArrayD.h"
#include "TObjArray.h"
#include <TSystem.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>

Double_t fBinsPi013TeVEMCPt[40]                  =   {0.0,  1.4,   1.6,   1.8,   2.0,   2.2,
                                                      2.4,   2.6,   2.8,   3.0,   3.2,
                                                      3.4,   3.6,   3.8,   4.0,   4.2,
                                                      4.4,   4.6,   4.8,   5.0,   5.2,
                                                      5.4,   5.6,   5.8,   6.0,   6.2,
                                                      6.4,   6.6,   6.8,   7.0,   7.5,
                                                      8.0,   8.5,   9.0,   9.5,  10.0,
                                                     12.0,  14.0,  16.0,  20.0};

// Double_t fBinsPi013TeVEMCPt[37]               = { 0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4,   1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
//                                                     3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10.0,  12.0, 16.0, 20.0, 25.0, 30.0, 35.0, 40.0};


// Double_t fBinsPi013TeVEMCPt[46]                  =   {0.0,  1.4,   1.6,   1.8,   2.0,   2.2,
//                                                       2.4,   2.6,   2.8,   3.0,   3.2,
//                                                       3.4,   3.6,   3.8,   4.0,   4.2,
//                                                       4.4,   4.6,   4.8,   5.0,   5.2,
//                                                       5.4,   5.6,   5.8,   6.0,   6.2,
//                                                       6.4,   6.6,   6.8,   7.0,   7.2,
//                                                       7.4,   7.6,   7.8,   8.0,   8.4,
//                                                       8.8,   9.2,   9.6,  10.0,  11.0,
//                                                      12.0,  13.0,  14.0,  16.0,  20.0};


const Int_t kMaxHit             = 2000;
const Double_t lowerparamrange  = 0.05;         //0.085;
const Double_t upperparamrange  = 0.225;
const Double_t lowercountrange  = 0.05;         //0.085;
const Double_t uppercountrange  = 0.225;
const int numberbins            = 39;
TH1D* hInvMass_Data             = NULL;         //data histogram
TH1D* hInvMass_MC               = NULL;
TH1D* mc_photon                 = NULL;         // gamma gamma MC histogram
TH1D* mc_MixedConvPhoton        = NULL;         // gamma gamma_conv MC histogram
TH1D* mc_ConvPhoton             = NULL;         // gamma_conv gamma_conv MC histogram
TH1D* testhisto2                = NULL;
TH1D* hPeak_MC                  = NULL;
TH1D* hCorrBkg                  = NULL;
TH1D* mc_full_clone1            = NULL;
TH1D* mc_full_clone2            = NULL;
TH1D* korrBG_clone1             = NULL;
TH1D* MCBG                      = NULL;
TH1D* DataBG                    = NULL;
TLegend *legiter                = NULL;
TLatex* chi_and_param42         = NULL;
TH1D* mc_full_clone42           = NULL;
TH1D* korrBG_clone42            = NULL;

////////////////////////////////////////////////////////////////////////////////
// Backgrund fitting stuff
TH1D* hBackStackup  = NULL;         // Histo of back. for adding on main back
                                    // histo
// TH1D* hPeak1         = NULL;        // peak histo, needed to subtract from
                                    // data to retrieve the background
TH1D* hPeak2         = NULL;        // peak histo, needed to subtract from
                                    // data to retrieve the background
// TH1D* hBack         = NULL;         // main background histo
/******************************************************************************/


void drawchi_and_param42(TLatex* tex,TFitResultPtr r ){
  tex->DrawLatexNDC(0.15,0.85,
  Form("#frac{#chi^{2}}{ndf} = %.2lf ",r->Chi2() / r->Ndf()));

  tex->DrawLatexNDC(0.17,0.70,
  Form("a = %.2lf ",r->Parameter(0)));

  tex->DrawLatexNDC(0.17,0.65,
  Form("b = %.2lf ",r->Parameter(1)));

  tex->DrawLatexNDC(0.17,0.60,
  Form("#frac{b}{a} = %.2lf ",r->Parameter(1) / r->Parameter(0)));
}


TString strData = "data (signal + corr. back.)";
TString MCInfo = "#splitline{pp, #sqrt{#it{s}} = 13TeV}{MC Monasch 13}";
TString massaxis = "#it{M}_{#gamma #gamma} (GeV/#it{c}^{2})";
TString rawyield = "#frac{1}{2#pi N_{evt}} #frac{d#it{N}}{d#it{p}_{T}}";
TString strCorrectedYield = "#frac{1}{2#pi #it{N}_{evt} #it{p}_{T}} #frac{d^{2}#it{N}}{d#it{y}d#it{p}_{T}}";
TString DoubleTempStr = "signal + corr. back.";
TString Pol1Str = "Peak temp. + 1^{st} ord. pol param.";
TLatex* texMCInfo = new TLatex();

// Double_t myFunc(Double_t x){
//   return (mc_photon->GetBinContent(mc_photon->FindBin(x))
//   +mc_MixedConvPhoton->GetBinContent(mc_MixedConvPhoton->FindBin(x))
//   +mc_ConvPhoton->GetBinContent(mc_ConvPhoton->FindBin(x)));
// }

Double_t* GetBinningFromHistogram(TH1D* hist){
  if(!hist) return 0;
  TArrayD* dArray = (TArrayD*)hist->GetXaxis()->GetXbins();
  return dArray->GetArray();
}

Int_t GetNBinningFromHistogram(TH1D* hist){
  if(!hist) return 0;
  TArrayD* dArray = (TArrayD*)hist->GetXaxis()->GetXbins();
  return dArray->GetSize();
}

void RotateToLabSystem(const float& theta, const float& phi,
		       const float& p1, const float& p2, const float& p3,
		       float& p1rot, float& p2rot, float& p3rot) {

  Float_t st = sin(theta);
  Float_t ct = cos(theta);
  Float_t sp = sin(phi);
  Float_t cp = cos(phi);

  p1rot = cp*ct*p1 - p2*sp + cp*p3*st;
  p2rot = cp*p2 + ct*p1*sp + p3*sp*st;
  p3rot = ct*p3 - p1*st;

}

// open rootfile safetly from Patrick Reichelt
TFile* /*LmHelper::*/SafelyOpenRootfile(const std::string filename)
{
  /// Opens a rootfile without affecting the active path, which otherwise would point into the file, often causing trouble.
  //save current path before opening rootfile.
  TString sPath = gDirectory->GetPath();

  TFile* ffile = 0x0;
  // check if file is already open.
  if ( gROOT->GetListOfFiles()->FindObject(filename.data()) ) {
    ffile = gROOT->GetFile(filename.data()); // avoid to open same file twice
  }
  if (ffile && ffile->IsOpen()) return ffile;

  ffile = TFile::Open(filename.data()); // gives root error and returns 0x0 on fail.
  //ffile = new TFile(filename.data()); // gives root error on fail.
  //ffile->OpenFile(filename.data());   // gives no error on fail.
  if (!ffile) printf(Form("SafelyOpenRootfile(): file '%s' not found.", filename.data()));

  // change to previous path again, so that it will be possible to close the file later without crash.
  // otherwise heap based objects will be created in memory that will be freed when the file is closed.
  gDirectory->Cd(sPath.Data());

  // alternatively one can do hist->SetDirectory(0); // according to Dario (Analysis Tutorial 26.06.2015)
  // but this seems not to work if the class object (which owns the hist) was created in the file path.

  return ffile;
}

class DataTree{
  private:
    Float_t pxdata[kMaxHit];
    Float_t pydata[kMaxHit];
    Float_t pzdata[kMaxHit];
    Int_t iNCluster;
    Int_t NEvents;
    TTree* tClusters;

  public:
    DataTree(TFile* fDaten){
      tClusters = (TTree*) fDaten->Get("ntu");
      NEvents = tClusters->GetEntries();
      tClusters->SetBranchAddress("nhit", &iNCluster);
      tClusters->SetBranchAddress("px", pxdata);
      tClusters->SetBranchAddress("py", pydata);
      tClusters->SetBranchAddress("pz", pzdata);
    }

    void GetEvent(Int_t iEvt){
      tClusters->GetEntry(iEvt);
    }

    Float_t GetPX(Int_t iEvt, Int_t iHit){
      GetEvent(iEvt);
      return pxdata[iHit];
    }

    Float_t GetPY(Int_t iEvt, Int_t iHit){
      GetEvent(iEvt);
      return pydata[iHit];
    }

    Float_t GetPZ(Int_t iEvt, Int_t iHit){
      GetEvent(iEvt);
      return pzdata[iHit];
    }
    Float_t GetClusterID(Int_t iEvt){
      GetEvent(iEvt);
      return iNCluster;
    }

    Int_t GetNEvents(){
      return NEvents;
    }
};




void SetCanvasStandardSettings(TCanvas *cCanv){
  gStyle->SetOptStat(kFALSE); // <- das hier macht dies box rechts oben weg
  cCanv->SetTopMargin(0.025);
  cCanv->SetBottomMargin(0.15);
  cCanv->SetRightMargin(0.05);
  cCanv->SetLeftMargin(0.15);
  cCanv->SetTickx();
  cCanv->SetTicky();
  cCanv->SetLogy(0);
  cCanv->SetLogx(0);
}


void SetHistoStandardSettings(TH1* histo, Double_t XOffset = 1.2, Double_t YOffset = 1., Double_t textSize = 0.03){
  histo->SetStats(0);
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetXaxis()->SetTitleSize(textSize);
  histo->GetYaxis()->SetTitleSize(textSize);
  histo->GetXaxis()->SetLabelSize(textSize);
  histo->GetYaxis()->SetLabelSize(textSize);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetTitleFont(42);
  histo->GetXaxis()->SetTitleFont(42);

  histo->SetTitle("");
  // histo->SetXTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
  histo->GetXaxis()->SetTitleOffset(1.4);
  // histo->SetYTitle("#frac{d#it{N}_{#gamma #gamma}}{d#it{m}_{inv}} (GeV/#it{c}^{2})^{-1}");
  histo->GetYaxis()->SetTitleOffset(1.4);
  histo->Sumw2();
  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(1.5);
  histo->SetLineWidth(2);
  histo->SetLineColor(kBlack);
  histo->SetMarkerColor(kBlack);
}



void SetHistoStandardSettings2(TH2* histo, Double_t XOffset = 1.2, Double_t YOffset = 1., Double_t textSize = 0.03){
  histo->SetStats(0);
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetXaxis()->SetTitleSize(textSize);
  histo->GetYaxis()->SetTitleSize(textSize);
  histo->GetXaxis()->SetLabelSize(textSize);
  histo->GetYaxis()->SetLabelSize(textSize);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetTitleFont(42);
  histo->GetXaxis()->SetTitleFont(42);

  histo->GetZaxis()->SetTitleOffset(YOffset);
  histo->GetZaxis()->SetLabelSize(textSize);
  histo->GetZaxis()->SetLabelFont(42);
  histo->GetZaxis()->SetTitleFont(42);
  // histo->Sumw2();


  // histo->SetTitle("");
  // histo->SetXTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
  // histo->GetXaxis()->SetTitleOffset(1.4);
  // histo->SetYTitle("#it{p}_{T} (GeV/#it{c})");
  // histo->GetYaxis()->SetTitleOffset(1.4);
  // histo->SetZTitle("#it{counts}");
  // histo->GetZaxis()->SetTitleOffset(1.4);
  // histo->GetZaxis()->SetRangeUser(1.e-10,100.);

}


void SetLegendSettigns(TLegend* leg, Double_t textSize = 0.03){
  leg->SetTextFont(42);
  leg->SetTextSize(textSize);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->SetLineColor(0);
  leg->SetMargin(0.15);
  leg->SetBorderSize(0);
}

void SetLatexSettings(TLatex* tex, Double_t textSize = 0.03){
  tex->SetTextSize(textSize);
  tex->SetTextFont(42);
  }

// gStyle->SetCanvasColor(0);
// gStyle->SetPadColor(0);
// gStyle->SetCanvasBorderMode(0);
// gStyle->SetPadBorderMode(0);
//
// gStyle->SetTitleXOffset(1.4);
// gStyle->SetTitleYOffset(1.8);
//
// gStyle->SetPadLeftMargin(0.17);
// gStyle->SetPadRightMargin(0.1);      // 0.1 = root default
// gStyle->SetPadTopMargin(0.1);
// gStyle->SetPadBottomMargin(0.14);




void printProgress (Double_t progress)
{
  int barWidth = 50;
  std::cout.flush();
  std::cout << "["<< int(progress * 100.0) << "%]" << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "|";
      else std::cout << " ";
  }
  std::cout << "]\r";
}


TString sigma_minv_str = TString("#it{#sigma} (GeV/#it{c}^{2})^{-1}");
TString minv_str = TString("#it{m}_{inv} (GeV/#it{c}^{2})");
TString pt_str = TString("#it{p}_{T} (GeV/#it{c})");
TString dNdmin_str = TString("#frac{d#it{N}_{#gamma #gamma}}{d#it{m}_{inv}} (GeV/#it{c}^{2})^{-1}");
TString poweek_str = TString("Powerweek Daten");
TString pi0togamma_str = TString("#pi^{0} #rightarrow #gamma #gamma");







//// Advanced

void DrawLabelALICE(Double_t startTextX = 0.13, Double_t startTextY = 0.9, Double_t textHeight = 0.04, Double_t textSize = 0.03, TString str = " "){
  TString textAlice       = "ALICE work in progress";
  TString textEvents      = "Data";
  TLatex *alice           = NULL;
  TLatex *pt              = NULL;
  TLatex *energy          = NULL;
  TLatex *detprocess      = NULL;
  TLatex *detprocess2     = NULL;
  TLatex *Template        = NULL;
  TLatex *Template2       = NULL;

  Double_t differenceText  = textHeight*1.7;
  if(str == "" || str == " "){
    alice            = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));

    pt               = new TLatex(startTextX, (startTextY-1.5*differenceText), str);

    energy           = new TLatex(startTextX, (startTextY-1.5*differenceText), "pp, #sqrt{#it{s}} = 13 TeV");

    detprocess       = new TLatex(startTextX, (startTextY-2.5*differenceText), "#pi^{0}#rightarrow#gamma#gamma");

    detprocess2      = new TLatex(startTextX, (startTextY-3.5*differenceText), "#gamma's rec. with EMCal ");

    Template         = new TLatex(startTextX, (startTextY-4.5*differenceText), "Templates:");

    Template2        = new TLatex(startTextX, (startTextY-5.5*differenceText), "Pythia 8 Monash 2013");
  }

  else{
    alice            = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));

    pt               = new TLatex(startTextX, (startTextY-1.5*differenceText), str);

    energy           = new TLatex(startTextX, (startTextY-2.5*differenceText), "pp, #sqrt{#it{s}} = 13 TeV");

    detprocess       = new TLatex(startTextX, (startTextY-3.5*differenceText), "#pi^{0}#rightarrow#gamma#gamma");

    detprocess2      = new TLatex(startTextX, (startTextY-4.5*differenceText), "#gamma's rec. with EMCal ");

    Template         = new TLatex(startTextX, (startTextY-5.5*differenceText), "Templates:");

    Template2        = new TLatex(startTextX, (startTextY-6.5*differenceText), "Pythia 8 Monash 2013");
  }

  alice->SetNDC();
  alice->SetTextColor(1);
  alice->SetTextFont(42);
  alice->SetTextSize(textSize);
  alice->DrawClone();

  energy->SetNDC();
  energy->SetTextColor(1);
  energy->SetTextSize(textSize);
  energy->SetTextFont(42);
  energy->DrawClone();

  Template->SetNDC();
  Template->SetTextColor(1);
  Template->SetTextSize(textSize);
  Template->SetTextFont(42);
  Template->DrawClone();

  Template2->SetNDC();
  Template2->SetTextColor(1);
  Template2->SetTextSize(textSize);
  Template2->SetTextFont(42);
  Template2->DrawClone();


  detprocess->SetNDC();
  detprocess->SetTextColor(1);
  detprocess->SetTextSize(textSize);
  detprocess->SetTextFont(42);
  detprocess->DrawClone();

  detprocess2->SetNDC();
  detprocess2->SetTextColor(1);
  detprocess2->SetTextSize(textSize);
  detprocess2->SetTextFont(42);
  detprocess2->DrawClone();

  pt->SetNDC();
  pt->SetTextColor(1);
  pt->SetTextSize(textSize);
  pt->SetTextFont(42);
  if(!(str == "" || str == " ")){
    pt->DrawClone();
  }


  delete alice;
  delete energy;
  delete Template;
  delete detprocess;
  delete pt;

}
