// printer.h file
#ifndef CommonHeader_H
#define CommonHeader_H

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
#include "TVirtualFitter.h"
#include <TSystem.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <time.h>


/**
 * pT Binning comming from the framework. If its changed there, you need to
 * change it here!
 */
Double_t fBinsPi013TeVEMCPt[40]                  =   {0.0,  1.4,   1.6,   1.8,   2.0,   2.2,
                                                      2.4,   2.6,   2.8,   3.0,   3.2,
                                                      3.4,   3.6,   3.8,   4.0,   4.2,
                                                      4.4,   4.6,   4.8,   5.0,   5.2,
                                                      5.4,   5.6,   5.8,   6.0,   6.2,
                                                      6.4,   6.6,   6.8,   7.0,   7.5,
                                                      8.0,   8.5,   9.0,   9.5,  10.0,
                                                     12.0,  14.0,  16.0,  20.0};

Int_t fBinsPi013TeVEMCPtRebin[39]               =   {4, 4, 4, 4, 4,    4, 4, 4, 4, 4,
                                                     4, 4, 4, 4, 4,    4, 4, 4, 4, 4,
                                                     4, 4, 4, 4, 4,    4, 4, 4, 4, 4,
                                                     8, 8, 8, 8, 8,   10,10,10,10};
/**
 * parametrisation- and counting intervall boarders.
 */
// Double_t lowerparamrange[38] = {0.0395, 0.0395, 0.0395, 0.0395, 0.0395,
//                                 0.0395, 0.0395, 0.0395, 0.0395, 0.0395,
//                                 0.0395, 0.0395, 0.0395, 0.0395, 0.0415,
//                                 0.0425, 0.0445, 0.0465, 0.0475, 0.0495,
//                                 0.0515, 0.0525, 0.0545, 0.0555, 0.0575,
//                                 0.0595, 0.0605, 0.0625, 0.0645, 0.0675,
//                                 0.0715, 0.0755, 0.0795, 0.0835, 0.0895,
//                                 0.1015, 0.1175, 0.1365};
Double_t lowerparamrange[38] = {0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1};
                                // {0.0195, 0.0205, 0.0225, 0.0235, 0.0255,
                                //  0.0275, 0.0285, 0.0305, 0.0315, 0.0335,
                                //  0.0355, 0.0365, 0.0385, 0.0395, 0.0415,
                                //  0.0425, 0.0445, 0.0465, 0.0475, 0.0495,
                                //  0.0515, 0.0525, 0.0545, 0.0555, 0.0575,
                                //  0.0595, 0.0605, 0.0625, 0.0645, 0.0675,
                                //  0.0715, 0.0755, 0.0795, 0.0835, 0.0895,
                                //  0.1015, 0.1175, 0.1365};
const Double_t upperparamrange  = 0.2;


// Double_t lowercountrange[38] = {0.0395, 0.0395, 0.0395, 0.0395, 0.0395,
//                                 0.0395, 0.0395, 0.0395, 0.0395, 0.0395,
//                                 0.0395, 0.0395, 0.0395, 0.0395, 0.0415,
//                                 0.0425, 0.0445, 0.0465, 0.0475, 0.0495,
//                                 0.0515, 0.0525, 0.0545, 0.0555, 0.0575,
//                                 0.0595, 0.0605, 0.0625, 0.0645, 0.0675,
//                                 0.0715, 0.0755, 0.0795, 0.0835, 0.0895,
//                                 0.1015, 0.1175, 0.1365};

Double_t lowercountrange[38] = {0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1};
const Double_t uppercountrange  = 0.2;

const Int_t kMaxHit             = 2000;
const int numberbins            = 39;           // number of actual used pT Bins
TH1D* hInvMass_Data             = NULL;         // data histogram
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

/**
 * definition of the color coding. This is needed to give the
 * SetHistoStandardSettings Function the option to change the color of the
 * histogram, since EColor, the type of the kBlack and so furth won't work.
 */
const int black   = 1;
const int teal    = 840;
const int red     = 632;
const int pink    = 900;
const int magenta = 616;
const int blue    = 600;
const int gray    = 920;

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


/**
 * functions for the mini toy MC for the opening angle cut.

 */
Double_t fCalcP(Double_t p1x, Double_t p1y, Double_t p1z){
  return sqrt(pow(p1x,2.)+pow(p1y,2.)+pow(p1z,2.));
}

Double_t fCalcTheta(Double_t p1x, Double_t p1y, Double_t p1z, Double_t p2x, Double_t p2y, Double_t p2z){
  return acos((p1x*p2x+ p1y*p2y + p1z*p2z)/(fCalcP(p1x, p1y, p1z)*fCalcP(p2x, p2y, p2z)));
}

Double_t fCalcInvMass(Double_t p1, Double_t p2, Double_t theta){
  return sqrt(2.*p1*p2*(1-cos(theta)));
}

Double_t fCalcPT(Double_t p1x, Double_t p1y, Double_t p2x, Double_t p2y){
  return sqrt(pow((p1x+p2x),2.)+pow((p1y+p2y),2.));
}

void fSmear(Double_t &px, Double_t &py, Double_t &pz){
  // one cell == 6 cm
  // distance to EMCal == 400 cm
  // -> 6/400 = 3sigma -> 1sigma = 1/200
  gRandom->Gaus(px, px/200.);
  gRandom->Gaus(py, py/200.);
  gRandom->Gaus(pz, pz/200.);
}

/**
 * function to draw Chi^2 and sacling factors within a canvas/pad.
 */
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

/**
 * String collection:
 * This collection is there, so that you do not need to write all these nasty
 * long strings when you give some axis a title.
 */
TString strData = "data (signal + corr. back.)";
TString MCInfo = "#splitline{pp, #sqrt{#it{s}} = 13TeV}{MC Monasch 13}";
TString rawyield = "#frac{1}{2#pi N_{evt}} #frac{d#it{N}}{d#it{p}_{T}}";
TString strCorrectedYield = "#frac{1}{2#pi #it{N}_{evt} #it{p}_{T}} #frac{d^{2}#it{N}}{d#it{y}d#it{p}_{T}}";
TString DoubleTempStr = "signal + corr. back.";
TString Pol1Str = "Peak temp. + 1^{st} ord. pol param.";
TString count_str = TString("d#it{N}/d#it{m}_{inv} ((GeV/#it{c}^{2})^{-1})");
TString StatUn_Str = TString("relative stat. Unsicherheit (%)");
TString sigma_minv_str = TString("#it{#sigma} (GeV/#it{c}^{2})^{-1}");
TString minv_str = TString("#it{m}_{inv} (GeV/#it{c}^{2})");
TString pt_str = TString("#it{p}_{T} (GeV/#it{c})");
TString dNdmin_str = TString("#frac{d#it{N}_{#gamma #gamma}}{d#it{m}_{inv}} (GeV/#it{c}^{2})^{-1}");
TString poweek_str = TString("Powerweek Daten");
TString pi0togamma_str = TString("#pi^{0} #rightarrow #gamma #gamma");

TLatex* texMCInfo = new TLatex();

/**
 * functions written by Patrick Huhn
 * used to get the Binning Array or the number of bins of an given TH1D.
 */
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

/**
 * Toy MC Stuff. System Rotation from CMS to Lab.
 */
void RotateToLabSystem(const double& theta, const double& phi,
		       const double& p1, const double& p2, const double& p3,
		       double& p1rot, double& p2rot, double& p3rot) {

  Double_t st = sin(theta);
  Double_t ct = cos(theta);
  Double_t sp = sin(phi);
  Double_t cp = cos(phi);

  p1rot = cp*ct*p1 - p2*sp + cp*p3*st;
  p2rot = cp*p2 + ct*p1*sp + p3*sp*st;
  p3rot = ct*p3 - p1*st;

}

// open rootfile safetly from Patrick Reichelt
/**
 * function to open root files safely which means, without changing the
 * gDirectory, which influences pointers. So if you create a pointer to a TH1D
 * after opening a root file normally, and then close the file, the pointer
 * would be broken. This is not good. We do not want this.
 * @filename        name of the file whoch you want to open
 * @return          return TFile Pointer to the file with corresponding filename
 */
TFile* SafelyOpenRootfile(const std::string filename)
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
  // if (!ffile) printf(Form("SafelyOpenRootfile(): file '%s' not found.", filename.data()));

  // change to previous path again, so that it will be possible to close the file later without crash.
  // otherwise heap based objects will be created in memory that will be freed when the file is closed.
  gDirectory->Cd(sPath.Data());

  // alternatively one can do hist->SetDirectory(0); // according to Dario (Analysis Tutorial 26.06.2015)
  // but this seems not to work if the class object (which owns the hist) was created in the file path.

  return ffile;
}


/**
 * DataTree Class used for the Toy MC aswell.
 */
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



/**
 * Functions to set Standard settings on the histograms
 * They probably all need some work done on them.
 */

/**
 * Sets the Standard Canvas Settings up.
 * @param cCanv pointer to the TCanvas you want set up.
 */
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

/**
 * Sets the Standard TH1 Settings up.
 * @param histo    pointer to TH1
 * @param XOffset  X-Title Offset
 * @param YOffset  Y-Title Offset
 * @param textSize Text size
 */
void SetHistoStandardSettings(TH1* histo, Double_t XOffset = 1.2, Double_t YOffset = 1., Double_t textSize = 35, int color = 1){
  histo->SetStats(0);
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetXaxis()->SetTitleSize(textSize);
  histo->GetYaxis()->SetTitleSize(textSize);
  histo->GetXaxis()->SetLabelSize(textSize);
  histo->GetYaxis()->SetLabelSize(textSize);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetTitleFont(43);

  histo->SetTitle("");
  histo->Sumw2();
  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(1.5);
  histo->SetLineWidth(3);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
}



void SetHistoStandardSettings2(TH2* histo, Double_t XOffset = 1.2, Double_t YOffset = 1., Double_t textSize = 35){
  histo->SetStats(0);
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetXaxis()->SetTitleSize(textSize);
  histo->GetYaxis()->SetTitleSize(textSize);
  histo->GetXaxis()->SetLabelSize(textSize);
  histo->GetYaxis()->SetLabelSize(textSize);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetTitleFont(43);

  histo->GetZaxis()->SetTitleOffset(YOffset);
  histo->GetZaxis()->SetLabelSize(textSize);
  histo->GetZaxis()->SetLabelFont(43);
  histo->GetZaxis()->SetTitleFont(43);

  histo->Sumw2();
}


void SetLegendSettigns(TLegend* leg, Double_t textSize = 35){
  leg->SetTextFont(43);
  leg->SetTextSize(textSize);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->SetLineColor(0);
  leg->SetMargin(0.15);
  leg->SetBorderSize(0);
}

void SetLatexSettings(TLatex* tex, Double_t textSize = 35){
  tex->SetTextSize(textSize);
  tex->SetTextFont(43);
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



/**
 * function from the Toy MC which prints out the progress in a fancy little way
 * into your command prompt
 * made by Adrian
 */
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


/**
 * Function to draw the ALICE label onto your canvas/pad
 * @param startTextX X-Startpostion
 * @param startTextY Y-Startpostion
 * @param textHeight influences the gap between the lines
 * @param textSize   Size of the Text Font
 * @param str        TString containing the pT range
 */
void DrawLabelALICE(Double_t startTextX = 0.13, Double_t startTextY = 0.9, Double_t textHeight = 0.04, Double_t textSize = 35, TString str = " "){
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
  alice->SetTextFont(43);
  alice->SetTextSize(textSize);
  alice->DrawClone();

  energy->SetNDC();
  energy->SetTextColor(1);
  energy->SetTextSize(textSize);
  energy->SetTextFont(43);
  energy->DrawClone();

  Template->SetNDC();
  Template->SetTextColor(1);
  Template->SetTextSize(textSize);
  Template->SetTextFont(43);
  Template->DrawClone();

  Template2->SetNDC();
  Template2->SetTextColor(1);
  Template2->SetTextSize(textSize);
  Template2->SetTextFont(43);
  Template2->DrawClone();


  detprocess->SetNDC();
  detprocess->SetTextColor(1);
  detprocess->SetTextSize(textSize);
  detprocess->SetTextFont(43);
  detprocess->DrawClone();

  detprocess2->SetNDC();
  detprocess2->SetTextColor(1);
  detprocess2->SetTextSize(textSize);
  detprocess2->SetTextFont(43);
  detprocess2->DrawClone();

  pt->SetNDC();
  pt->SetTextColor(1);
  pt->SetTextSize(textSize);
  pt->SetTextFont(43);
  if(!(str == "" || str == " ")){
    pt->DrawClone();
  }


  delete alice;
  delete energy;
  delete Template;
  delete detprocess;
  delete pt;

}
#endif
