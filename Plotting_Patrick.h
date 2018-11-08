#ifndef Plotting_Patrick
#define Plotting_Patrick

#include "TH1D.h"
#include "TLatex.h"
#include "TObjArray.h"
#include "TPad.h"
#include "TCanvas.h"

#include "TLegend.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom.h"
#include "TImage.h"
#include "TTimeStamp.h"
#include "TMath.h"

#include <iostream>

Double_t minXValue = 0.15;
Double_t maxXValue = 99.9;
Double_t maxYValue = -1;
Double_t minYValue = 1e20;
Double_t fPtPlotRange = 49.9;

TCanvas *makeCanvas(TObjArray *histArray, TObjArray *ratioArray = 0,const char *controlstring="", Short_t *colorArray = 0, Short_t *markerArray = 0);
void SetHistogramProperties(TH1D* &h, TString XTitle, TString YTitle, Int_t ColorIndex);

Int_t GetNiceColor(Int_t i)
{//
  Int_t niceColors[] = {kRed+1, kGreen-3, kBlue+1, kViolet-4, kOrange-3, kBlack, kCyan-2, kGray+2, kOrange+2, kBlack, kYellow+3};
  return niceColors[i%11];
}

void PlotArray(TObjArray *arraytoplot, const char *controlstring ,Short_t *colorArray, Short_t *markerArray, Double_t relativeTextSize,Float_t lableSize,Float_t titleSize,  Int_t lableFont,Int_t titleFont,Float_t titleOffsetY,Float_t titleOffsetX,Float_t margin,Double_t xMin,Double_t xMax,Double_t yMin = 0,Double_t yMax = 0){
  TString control(controlstring);
  Int_t nbrEntries = arraytoplot->GetEntries();
  for(Int_t hh=0;hh<=nbrEntries;hh++){
    if(!arraytoplot->At(hh)){break;}

    if (arraytoplot->At(hh)->InheritsFrom("TH1") && !arraytoplot->At(hh)->InheritsFrom("TH2")){
      cout<<"| Drawing a TH1 "<<hh<<endl;
      TH1D *hist = (TH1D*) arraytoplot->At(hh);
      hist->SetStats(0);
      hist->GetXaxis()->SetLabelSize(lableSize*relativeTextSize);
      hist->GetXaxis()->SetTitleSize(titleSize*relativeTextSize);
      hist->GetXaxis()->SetLabelFont(lableFont);
      hist->GetXaxis()->SetTitleFont(titleFont);
      hist->GetXaxis()->SetTitleOffset(titleOffsetX);

      hist->GetYaxis()->SetLabelSize(lableSize*relativeTextSize);
      hist->GetYaxis()->SetTitleSize(titleSize*relativeTextSize);
      hist->GetYaxis()->SetLabelFont(lableFont);
      hist->GetYaxis()->SetTitleFont(titleFont);
      hist->GetYaxis()->SetTitleOffset(titleOffsetY);
      hist->SetTitle("");

      hist->SetAxisRange(xMin,xMax,"X");



      TString xTitle(hist->GetXaxis()->GetTitle());
      if (xTitle.Contains("pt")  || xTitle.Contains("p_{T} (GeV/c)")){
        hist->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      }else if(xTitle.Contains("mcPt")){
        hist->GetXaxis()->SetTitle("#it{p}_{T}^{MC} (GeV/#it{c})");
      }else if (xTitle.Contains("minv")  || xTitle.Contains("m_{inv} (GeV/c^2)")) {
        hist->GetXaxis()->SetTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
      }



      // hist->SetMarkerStyle(markerArray[hh]);
      // hist->SetMarkerColor(colorArray[hh]);
      // hist->SetLineColor(colorArray[hh]);


      Double_t mrkSize = hist->GetMarkerSize();
      if(hist->GetMarkerStyle()== 33 || hist->GetMarkerStyle()== 27) hist->SetMarkerSize(1.4*mrkSize);
      if(hist->GetMarkerStyle()== 29 || hist->GetMarkerStyle()== 30) hist->SetMarkerSize(1.4*mrkSize);
      if(hist->GetMarkerStyle()== 5  || hist->GetMarkerStyle()== 5 ) hist->SetMarkerSize(1.4*mrkSize);
      if(hist->GetMarkerStyle()== 34 || hist->GetMarkerStyle()== 28) hist->SetMarkerSize(1.2*mrkSize);


      /// Set tick length independent of pad size
      Double_t pad_width  = gPad->XtoPixel(gPad->GetX2());
      Double_t pad_height = gPad->YtoPixel(gPad->GetY1());
      hist->SetTickLength(15/pad_height,"X");
      hist->SetTickLength(15/pad_width,"Y");


      if(control.Contains("Thick")||control.Contains("thick")){
        hist->SetLineWidth(4);
      }else if(control.Contains("Thin")||control.Contains("thin")){
        hist->SetLineWidth(1);
      }else{
        hist->SetLineWidth(2);
      }

      if(markerArray[hh]==0){
        hist->SetLineWidth(0);
      }
      if(control.Contains("Lines")||control.Contains("lines")){
        cout<<"|         - Lines No errors"<<endl;
        if(markerArray){ hist->SetLineStyle(markerArray[hh]);}
        if(hh==0) { hist->Draw("hist");}
        else{hist->Draw("Same hist"); }
      }else if(markerArray[hh]==999){
        hist->SetFillColor(kGray);
        hist->SetMarkerSize(0);
        if(hh==0) {hist->Draw("E2");}
        else{hist->Draw("Same E2");}
      }else{
        if(hh==0) {hist->Draw("");}
        else{hist->Draw("Same");}
      }
    }else if (arraytoplot->At(hh)->InheritsFrom("TH2")){
      cout<<"| Drawing a TH2 "<<hh<<endl;
      TH2D *hist2 = (TH2D*) arraytoplot->At(hh);
      hist2->SetStats(0);
      hist2->GetXaxis()->SetLabelSize(lableSize*relativeTextSize);
      hist2->GetXaxis()->SetTitleSize(titleSize*relativeTextSize);
      hist2->GetXaxis()->SetLabelFont(lableFont);
      hist2->GetXaxis()->SetTitleFont(titleFont);
      hist2->GetXaxis()->SetTitleOffset(titleOffsetX);

      hist2->GetYaxis()->SetLabelSize(lableSize*relativeTextSize);
      hist2->GetYaxis()->SetTitleSize(titleSize*relativeTextSize);
      hist2->GetYaxis()->SetLabelFont(lableFont);
      hist2->GetYaxis()->SetTitleFont(titleFont);
      hist2->GetYaxis()->SetTitleOffset(titleOffsetY);
      hist2->SetTitle("");

      hist2->SetAxisRange(xMin,xMax,"X");
      hist2->SetAxisRange(yMin,yMax,"Y");


      TString xTitle(hist2->GetXaxis()->GetTitle());
      if (xTitle.Contains("pt")  || xTitle.Contains("p_{T} (GeV/c)")){
        hist2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      }else if(xTitle.Contains("mcPt")){
        hist2->GetXaxis()->SetTitle("#it{p}_{T}^{MC} (GeV/#it{c})");
      }


      /// Set tick length independent of pad size
      Double_t pad_width  = gPad->XtoPixel(gPad->GetX2());
      Double_t pad_height = gPad->YtoPixel(gPad->GetY1());
      hist2->SetTickLength(15/pad_height,"X");
      hist2->SetTickLength(15/pad_width,"Y");

      hist2->Draw("COLZ");

    }else if (arraytoplot->At(hh)->InheritsFrom("TLine")){
      cout<<"| Drawing line: "<<arraytoplot->At(hh)->GetName()<<endl;
      TLine *line = (TLine*) arraytoplot->At(hh);
      if (markerArray[hh]<11){line->SetLineStyle(markerArray[hh]);}
      line->Draw("Same");
    }else if (arraytoplot->At(hh)->InheritsFrom("TPaveText")){
      cout<<"| Drawing a TPaveText "<<hh<<endl;
      TPaveText *pt = (TPaveText*) arraytoplot->At(hh);
      if (! pt) cout<<"|Error with TPaveText"<<endl;
      cout<<"| Rel. text size: "<<relativeTextSize<<endl;
      pt->SetFillStyle(0);
      pt->SetFillColor(0);
      pt->SetBorderSize(0);
      pt->SetTextFont(titleFont);
      pt->SetTextSize(relativeTextSize);
      pt->Draw();
    }else if (arraytoplot->At(hh)->InheritsFrom("TLegend")){
      cout<<"| Drawing a Legend "<<hh<<endl;
      cout<<"| Rel. text size: "<<relativeTextSize<<endl;
      TLegend *leghh = (TLegend*) arraytoplot->At(hh);
      leghh->SetTextFont(titleFont);
      leghh->SetTextSize(relativeTextSize);
      leghh->SetFillColor(kWhite);
      leghh->SetFillStyle(0);
      leghh->SetBorderSize(0);
      leghh->Draw("SAME");
    }else if (arraytoplot->At(hh)->InheritsFrom("TBox")){
      cout<<"| Drawing a TBox"<<endl;
      TBox *box = (TBox*) arraytoplot->At(hh);
      box->Draw();
    }else if (arraytoplot->At(hh)->InheritsFrom("TGraphAsymmErrors")){
      cout<<"| Drawing a TGraphAsymmErrors "<<hh<<endl;
      TGraphAsymmErrors *asym = (TGraphAsymmErrors*) arraytoplot->At(hh);

      asym->SetLineColor(colorArray[hh]);
      asym->SetLineStyle(markerArray[hh]);
      asym->SetLineWidth(4);
      asym->Draw("SAME");
    }else if (arraytoplot->At(hh)->InheritsFrom("TGraphErrors")){
      cout<<"| Drawing a TGraph "<<hh<<endl;
      TGraph *grph = (TGraphAsymmErrors*) arraytoplot->At(hh);
        cout<<"| - As sytematic error"<<endl;
        grph->Draw("SAME []2");
    }else if (arraytoplot->At(hh)->InheritsFrom("TGraph")){
      cout<<"| Drawing a TGraph "<<hh<<endl;
      TGraph *grph = (TGraphAsymmErrors*) arraytoplot->At(hh);

      if(markerArray[hh] == 999){
        cout<<"| - As sytematic error"<<endl;
        grph->SetFillColorAlpha(colorArray[hh],0.6);
        //         grph->Draw("SAME a3");
        grph->Draw("SAME []");
        //         grph->Draw("SAME 5");
      }else{
        grph->SetMarkerColor(colorArray[hh]);
        grph->SetMarkerStyle(markerArray[hh]);
        Double_t mrkSize = grph->GetMarkerSize();
        if(grph->GetMarkerStyle()== 33 || grph->GetMarkerStyle()== 27)    grph->SetMarkerSize(1.4*mrkSize);
        if(grph->GetMarkerStyle()== 29 || grph->GetMarkerStyle()== 30)    grph->SetMarkerSize(1.4*mrkSize);
        if(grph->GetMarkerStyle()== 34 || grph->GetMarkerStyle()== 28)    grph->SetMarkerSize(1.2*mrkSize);
        grph->Draw("SAME P");
      }
    }else if (arraytoplot->At(hh)->InheritsFrom("TF1")){
      cout<<"| Drawing a TF1 "<<hh<<endl;
      TF1 *fun = (TF1*) arraytoplot->At(hh);
      fun->SetLineColor(colorArray[hh]);
      if (markerArray[hh]< 10){
        fun->SetLineStyle(markerArray[hh]);
      }else{
        fun->SetLineStyle(3);
      }
      if(control.Contains("Thick")||control.Contains("thick")){
        fun->SetLineWidth(4);
      }else if(control.Contains("Thin")||control.Contains("thin")){
        fun->SetLineWidth(1);
      }else{
        fun->SetLineWidth(2);
      }
      fun->Draw("SAME");
    }
    gPad->RedrawAxis();
  }
}

TCanvas *makeCanvas(TObjArray *histArray, TObjArray *ratioArray,const char *controlstring, Short_t *colorArray, Short_t *markerArray){


  /// This fuction draws a canvas, which contains two panels. The upper one takes up the
  /// most space and is designed to show spectra.
  /// The lower panel is designed to show the ratio of the histograms shown in the upper
  /// panel.The histograms and the ratios as well as the ledgend are created befor and
  /// are part of the parameters.
  ///
  /// The function expects a <b>TObjArray</b> for the Historgrams and a <b>Short_t</b> Array for the
  /// colors
  ///
  /// \param TObjArray Array of histograms
  /// \param TObjArray Array of Ratios
  /// \param Char String ("LogY"&"LogX" to set Axis)
  /// \param Short_t Array of color integers
  /// \param Short_t Array of marker integers
  ///
  /// \return TCanvas with the drawn histograms
  ///
  /// <b>String posibilities:</b>
  /// \arg LogX     Set X-axis logarythmic
  /// \arg LogY     Set Y-axis logarythmic
  /// \arg Square   Produce a squared canvas
  ///
  ///
  ///
  cout<<" ____________________________________________________________________________"<<endl;
  cout<<"|"<<endl;
  cout<<"| Plots generated via makeCanvas"<<endl;

  // Define common colors for plots
  // Particle composition
  Int_t colKaon = kBlue-7;
  Int_t mkrKaon = 20;
  Int_t mkrKaonO = 24;
  Int_t colPion = kBlue+1;
  Int_t mkrPion = 21;
  Int_t mkrPionO = 25;
  Int_t colProt = kRed+1;
  Int_t mkrProt = 33;
  Int_t mkrProtO = 27;
  Int_t colSigp = kCyan+3;
  Int_t mkrSigm = 34;
  Int_t mkrSigmO = 28;
  Int_t colSigm = kGreen+3;
  Int_t mkrSigp = 29;
  Int_t mkrSigpO = 30;
  Int_t colRest = kGray +1;
  Int_t mkrRest = 5;
  Int_t mkrRestO = 5;

  Int_t colPrim = kCyan+3;
  Int_t mrkPrim = 33;
  Int_t mrkPrimO = 27;

  Int_t colSec  = kGreen+3;
  Int_t mrkSec = 34;
  Int_t mrkSecO = 28;

  // Centrality
  Int_t col05   = kRed+1;
  Int_t col510  = kRed-1;
  Int_t col1020 = kOrange-3;
  Int_t col2030 = kYellow-2;
  Int_t col3040 = kGreen-3;
  Int_t col4050 = kGreen+3;
  Int_t col5060 = kCyan+3;
  Int_t col6070 = kBlue-7;
  Int_t col7080 = kBlue+1;

  // Energies
  Int_t col276TeV = kBlue+1;
  Int_t mrk276TeV = 33;
  Int_t mrk276TeVO = 27;
  Int_t col502TeV = kGreen-2;
  Int_t mrk502TeV = 21;
  Int_t mrk502TeVO = 25;
  Int_t col544TeV = kOrange-3;
  Int_t mrk544TeV = 34;
  Int_t mrk544TeVO = 28;

  // Data MC comparisons
  Int_t colData = kBlue-7;
  Int_t mrkData = 20;
  Int_t mrkDataO= 24;
  Int_t colMC   = kRed+1;
  Int_t mrkMC   = 21;
  Int_t mrkMCO  = 25;

  gStyle->SetTextFont(43);
  TString control(controlstring);

  //   Short_t colorArray[]={kRed,kOrange-3,kGreen-3,kGreen+3,kCyan+3,kBlue};
    Short_t defaultColorArray[14] = {0};
  if(!colorArray) {
      Short_t copyArray[14]={kBlue-7,kRed+1,kGreen-3,6,8,9,11,12,13,14,15,16,17,18};
      for (int i =0; i<14; ++i) {
          defaultColorArray[i] = copyArray[i];
      }

  } /// \remark If no colors are specified, use default: blue, red, green, magenta, wierd green, wierd blue and various grey
    Short_t defaultMarkerArray[14] = {0};
  if(!markerArray) {
      Short_t copyArray[14]={20,21,22,23,24,25,26,27,28,29,30,20,20,20};
      for (int i =0; i<14; ++i) {
          defaultMarkerArray[i] = copyArray[i];
      }
  } /// \remark If no markers are specified, use default.
  // if(control.Contains("OwnColors")){//DefineOwnColors();
  //     Short_t copyArray[14]={4001,4002,4003,4004,4005,4006,4007,4008,4009,4010,4011,4012,17,18};
  //     for (int i =0; i<14; ++i) {
  //         defaultColorArray[i] = copyArray[i];
  //     }
  // }

  Int_t nHist = histArray->GetEntries();
  if(nHist == 0){cout<<"| ERROR: Arrays empty"<<endl; return nullptr;}

  Float_t lableSize = 1.40; // Factor the label will be larger than the relative textsize
  Float_t titleSize = 1.40; // Factor the title will be larger than the relative textsize
  Float_t textSizeFactor = 12000;
  Int_t lableFont= 43;
  Int_t titleFont= 43;

  Float_t canvasWidth;
  Float_t canvasHeight;
  Float_t padFraction;
  Float_t titleOffsetY;
  Float_t titleOffsetX;
  Float_t sideMargin;
  Float_t topMargin;
  Float_t lowMargin;

  if(control.Contains("square")||control.Contains("Square")||control.Contains("SQUARE")){
    cout<<"| - Square Canvas"<<endl;
    canvasWidth=700;
    canvasHeight=700;
    padFraction = 0.5;
    titleOffsetY=1.4;
    titleOffsetX=1.5;
    sideMargin=0.12;
    topMargin=sideMargin-0.05;
    lowMargin=sideMargin+0.05;
    if(ratioArray){
      textSizeFactor =0.5 * 12000;
      sideMargin=0.13;
      topMargin = sideMargin/(padFraction)-0.05;
      lowMargin = sideMargin/(padFraction)+0.05;
      titleOffsetY=1.4;
      titleOffsetX=2.8;
    }
  }else if (control.Contains("horizontal")||control.Contains("Horizontal")||control.Contains("HORIZONTAL")){
    cout<<"| - Horizontal Canvas"<<endl;
    canvasWidth=600*1.414213562;
    canvasHeight=600;
    padFraction = 0.5;
    titleOffsetY=1.05;
    titleOffsetX=1.02;
    sideMargin = 0.09*1.414213562;
    topMargin=0.06;
    lowMargin=0.12;
    if (control.Contains("BinLabels")){
      titleOffsetY=1.05;
      titleOffsetX=2.1;
      lowMargin=0.19;
    }
    if(ratioArray){
    }
  }else if (control.Contains("A4")||control.Contains("a4")){
    cout<<"| - DIN-A4 shaped Canvas"<<endl;
    canvasWidth=600;
    canvasHeight=600*1.414213562;
    padFraction = 0.25;
    titleOffsetY=1.4;
    titleOffsetX=1.4;
    sideMargin = 0.11;
    topMargin=sideMargin * 1.414213562 - 0.05;
    lowMargin=sideMargin * 1.414213562 + 0.05;
    if(ratioArray){
    }
  }else if (control.Contains("OldAspect")||control.Contains("oldaspect")){
    cout<<"| - Improved aspect ratio Canvas"<<endl;
    canvasWidth=700;
    canvasHeight=800;
    padFraction = 0.28;
    titleOffsetY=1.7;
    titleOffsetX=3.5;
    sideMargin=0.12;
    topMargin=0.11;
    lowMargin=0.13;
    if(ratioArray){
      titleOffsetY=1.7;
      titleOffsetX=4.1;
      topMargin=0.08;
      lowMargin=0.30;
    }
  }else{
    cout<<"| - Improved aspect ratio Canvas"<<endl;
    textSizeFactor =1.8 * 12000;
    canvasWidth=800;
    canvasHeight=800*1.41421356;
    padFraction = 0.25;
    titleOffsetY=2;
    titleOffsetX=4.5;
    sideMargin = 0.14;
    topMargin = 0.12 * 1.414213562/2;
    lowMargin = 0.12 * 1.414213562/2;
    if(ratioArray){
      titleOffsetY=2;
      titleOffsetX=4.5;
      topMargin = (0.12 * 1.414213562/2)/(1-padFraction);
      lowMargin = (0.12 * 1.414213562/2)/(padFraction);
    }
  }



  if (control.Contains("CMYK")||control.Contains("cmyk")){gStyle->SetColorModelPS(1);}

  TRandom *random = new TRandom(histArray->Hash());

    TString title("default");
    title += Form("%f",random->Rndm());
  /// Check if first element of th histogram array is indeed an Histogram.
    if(!histArray->At(0)) return nullptr;
    if (!histArray->At(0)->InheritsFrom("TH1")){cout<<"|ERROR: First entry has to be an Histogram"<<endl; return nullptr;}
    TH1D *hist0 = (TH1D*) histArray->At(0);
    if(!hist0){cout<<"| ERROR: Histogram empty"<<endl; return nullptr;}

    Double_t xMin = hist0->GetXaxis()->GetBinLowEdge(hist0->GetXaxis()->GetFirst())+0.0000001;
    Double_t xMax = hist0->GetXaxis()->GetBinUpEdge(hist0->GetXaxis()->GetLast())-0.0000001;
    cout<<"| - Axis range: "<<xMin<<"-"<<xMax<<endl;
    hist0->GetYaxis()->CenterTitle(0);
    if(control.Contains("CenterTitle")||control.Contains("centertitle")) hist0->GetYaxis()->CenterTitle(1);

  /// Check if the first Element of the ratio array is an histogram.
  if(ratioArray){
    if(ratioArray->At(0)){
      if (!ratioArray->At(0)->InheritsFrom("TH1")){cout<<"|ERROR: First entry has to be an Histogram"<<endl;}
      TH1D *ratio0 = (TH1D*) ratioArray->At(0);
      if(!ratio0){cout<<"| ERROR: Histogram empty"<<endl;}
      ratio0->GetXaxis()->SetTickLength(0.06);
      ratio0->GetYaxis()->SetNdivisions(506);
      ratio0->SetLineWidth(2);
      ratio0->GetYaxis()->CenterTitle(0);

      /// Automatic range + zoom
      if(control.Contains("Zoom")){
        Double_t ratio0Maximum = ratio0->GetBinContent(ratio0->GetMaximumBin());
        Double_t ratio0Minimum = 1;
        for(Int_t ii = 0; ii < ratio0->GetNbinsX(); ii++){
          Double_t val = ratio0->GetBinContent(ii);
          if (val > 0){if (val < ratio0Minimum){ratio0Minimum = val;}}
        }
        Double_t zoom = 0.05;
        if(control.Contains("Zoom"))     zoom = 0.15;
        if(control.Contains("Zoom+"))    zoom = 0.30;
        if(control.Contains("Zoom++"))   zoom = 0.45;
        if(control.Contains("Zoom+++"))  zoom = 0.60;
        if(control.Contains("Zoom++++")) zoom = 0.75;
        if(control.Contains("Zoom+++++"))zoom = 0.90;
        if ( ( 1 - ratio0Minimum) > ( ratio0Maximum - 1) ){zoom = zoom*(1- ratio0Minimum); ratio0->SetAxisRange(ratio0Minimum+zoom,2-ratio0Minimum-zoom,"Y");}
        else{zoom = zoom*(ratio0Maximum -1); if ( ratio0Maximum > 2){ratio0->SetAxisRange(zoom,ratio0Maximum - zoom,"Y");}else{ ratio0->SetAxisRange(1-(ratio0Maximum-1)+zoom,ratio0Maximum - zoom,"Y");}}
      }
    }
  }

  /// Create Canvas in a given size
  TCanvas *newCanvas = new TCanvas(title,title,10,10,canvasWidth+10,canvasHeight+10);
  newCanvas->SetFillStyle(4000);
  newCanvas->cd();

    TPad *upperPad = nullptr;

  /// If ratioArray exist, split the Canvas in two and create a ratio panel
  if(ratioArray){
    /// Create upper pad for the histograms \remark upper pad has 72%
    upperPad = new TPad("upperPad","Distribution" ,0 ,padFraction ,1 ,1);
    TPad *lowerPad = new TPad("lowerPad","Ratio"        ,0 ,0           ,1 ,padFraction);

    upperPad->SetFillStyle(4000);
    upperPad->SetTopMargin(topMargin);
    upperPad->SetLeftMargin(sideMargin);
    upperPad->SetRightMargin(sideMargin);
    upperPad->SetBottomMargin(0.0);
    upperPad->SetTickx(1);
    upperPad->SetTicky(1);
    upperPad->Draw();

    lowerPad->SetFillStyle(4000);
    lowerPad->SetTopMargin(0.0);
    lowerPad->SetLeftMargin(sideMargin);
    lowerPad->SetRightMargin(sideMargin);
    lowerPad->SetBottomMargin(lowMargin);
    lowerPad->SetTickx(1);
    lowerPad->SetTicky(1);
    lowerPad->Draw();


    if(control.Contains("logX")||control.Contains("logx")||control.Contains("LogX")||control.Contains("LOGX")){
      upperPad->SetLogx(1); lowerPad->SetLogx(1); cout<<"| - LogX"<<endl;
    }
    if(control.Contains("logY")||control.Contains("logy")||control.Contains("LogY")||control.Contains("LOGY")){
      upperPad->SetLogy(1);  cout<<"| - LogY"<<endl;
    }
    if(control.Contains("RatioGridY")){lowerPad->SetGridy(1); cout<<"| - Ratio: GridY"<<endl;}

    newCanvas->cd();
    upperPad->cd();

    Double_t relativeTextSize;
    Double_t pad_width;
    Double_t pad_height;

    /// Calculate the relative text size, independent of the size of the pad
    pad_width  = gPad->XtoPixel(gPad->GetX2());
    pad_height = gPad->YtoPixel(gPad->GetY1());
    if (pad_width < pad_height){relativeTextSize = 1/pad_width;}
    else{relativeTextSize = 1/pad_height;}
    relativeTextSize = textSizeFactor * relativeTextSize;

    /// Use the standardtised PlotArray function to loop over the histogram an to draw the elements.
    if(!markerArray && !colorArray){PlotArray(histArray,controlstring,defaultColorArray,defaultMarkerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, sideMargin, xMin,xMax);}
    else if(!markerArray){PlotArray(histArray,controlstring,colorArray,defaultMarkerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, sideMargin, xMin,xMax);}
    else if(!colorArray){ PlotArray(histArray,controlstring,defaultColorArray,markerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, sideMargin, xMin,xMax);}
    else{PlotArray(histArray,controlstring,colorArray,markerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, sideMargin, xMin,xMax);}

    /// Go to the ratio pad and repeat the procedure
    lowerPad->cd();
    if(!markerArray && !colorArray){PlotArray(ratioArray,controlstring,defaultColorArray,defaultMarkerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, sideMargin, xMin,xMax);}
    else if(!markerArray){PlotArray(ratioArray,controlstring,colorArray,defaultMarkerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, sideMargin, xMin,xMax);}
    else if(!colorArray){PlotArray(ratioArray,controlstring,defaultColorArray,markerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, sideMargin, xMin,xMax);}
    else{PlotArray(ratioArray,controlstring,colorArray,markerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, sideMargin, xMin,xMax);}
    upperPad->cd();

  }else{// If no ratios are given:

    newCanvas->SetLeftMargin(sideMargin);
    newCanvas->SetRightMargin(sideMargin);
    newCanvas->SetTopMargin(topMargin);
    newCanvas->SetBottomMargin(lowMargin);


    newCanvas->SetTickx(1);
    newCanvas->SetTicky(1);


    if(control.Contains("logX")||control.Contains("logx")||control.Contains("LogX")||control.Contains("LOGX")){
      newCanvas->SetLogx(1); cout<<"| - LogX"<<endl;}
    if(control.Contains("logY")||control.Contains("logy")||control.Contains("LogY")||control.Contains("LOGY")){
      newCanvas->SetLogy(1); cout<<"| - LogY"<<endl;}
    if(control.Contains("gridY")||control.Contains("gridy")||control.Contains("GridY")||control.Contains("GRIDY")){
      if(!control.Contains("atio") ) newCanvas->SetGridy(1);cout<<"| - Grid Y"<<endl;
    }

    Double_t relativeTextSize;
    Double_t pad_width  = gPad->XtoPixel(gPad->GetX2());
    Double_t pad_height = gPad->YtoPixel(gPad->GetY1());
    if (pad_width < pad_height){relativeTextSize = 1/pad_width;}
    else{relativeTextSize = 1/pad_height;}
    relativeTextSize = textSizeFactor * relativeTextSize;

    if(!markerArray && !colorArray){PlotArray(histArray,controlstring,defaultColorArray,defaultMarkerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, sideMargin, xMin,xMax);}
    else if(!markerArray){PlotArray(histArray,controlstring,colorArray,defaultMarkerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, sideMargin, xMin,xMax);}
    else if(!colorArray){PlotArray(histArray,controlstring,defaultColorArray,markerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, sideMargin, xMin,xMax);}
    else{PlotArray(histArray,controlstring,colorArray,markerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, sideMargin, xMin,xMax);}
  }

  /// Draw extras
  /// Timestamp
  if(!control.Contains("Notime")&&!control.Contains("NoTime")&&!control.Contains("notime")&&!control.Contains("NOTIME")){
    TTimeStamp time;
    TString stTime;
    stTime=Form("#color[16]{%i}",time.GetDate());
    TLatex *timeStamp = new TLatex();
    timeStamp->SetNDC();
    timeStamp->SetTextAlign(22);
    timeStamp->SetTextFont(63);
    timeStamp->SetTextSizePixels(10);
    timeStamp->DrawLatex(0.96,0.98,stTime);
  }else{
    cout<<"|         - No Timestamp"<<endl;
  }
  /// Text "This work"
  if(control.Contains("ThisWork")||control.Contains("Thiswork")||control.Contains("thiswork")||control.Contains("THISWORK")){
    TLatex *thsWrk = new TLatex();
    thsWrk->SetNDC();
    thsWrk->SetTextAlign(22);
    thsWrk->SetTextFont(43);
    thsWrk->SetTextSizePixels(20);
    thsWrk->DrawLatex(0.83,0.94,"#color[16]{this work}");
  }
  /// Unified or work logo
  if(control.Contains("Logo") ||control.Contains("logo") ||control.Contains("LOGO")){
    Double_t scale = 0.1;
    Double_t x1 = 0.7; Double_t y1 = 0.7;
    if(control.Contains("square")||control.Contains("Square")||control.Contains("SQUARE")){ x1 = 0.7;  y1 = 0.7;}
    else if (control.Contains("horizontal")||control.Contains("Horizontal")||control.Contains("HORIZONTAL")){ x1 = 0.75;  y1 = 0.8;}
    else{ x1 = 0.7;  y1 = 0.8;}
    Double_t x2 = x1+scale;
    Double_t y2 = y1+scale;

    Char_t* logopath = nullptr;
    if(control.Contains("Unified") ||control.Contains("unified") ||control.Contains("UNIFIED")){logopath = "./graphics/einhorn.pdf";}
    if(control.Contains("Work") ||control.Contains("work") ||control.Contains("WORK")){logopath = "./graphics/wrkprgs.png";}
    TImage* logo = TImage::Open(logopath);
    TPad *pad_logo = new TPad("pad_logo","pad_logo",x1,y1,x2,y2);
    pad_logo->SetMargin(0,0,0,0);
    newCanvas->cd();
    pad_logo->Draw();
    pad_logo->cd();
    logo->Draw();
  }

  if(ratioArray && upperPad){upperPad->cd();}
  cout<<"|_____________________________________________________________________________"<<endl;
  cout<<" "<<endl;
  return newCanvas;
}

void SetHistogramProperties(TH1D* &h, TString XTitle, TString YTitle,
                            Int_t ColorIndex, Double_t XMin, Double_t XMax)
{
  if(!h) {Printf("Can't read Histogram!"); return;}
  h->SetTitle("");
  h->SetLineWidth(3);
  h->SetTitleFont(43);
  h->SetTitleSize(0);
  h->SetXTitle(XTitle);
  h->SetYTitle(YTitle);
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetTitleOffset(1.5);
  h->SetTitleSize(24, "X");
  h->SetTitleSize(24, "Y");
  h->SetTitleFont(43, "X");
  h->SetTitleFont(43, "Y");
  h->SetLabelFont(43, "X");
  h->SetLabelFont(43, "Y");
  h->SetLabelSize(24, "X");
  h->SetLabelSize(24, "Y");
  h->SetMarkerSize(1.5);
  h->SetStats(kFALSE);
  h->GetXaxis()->SetRangeUser(XMin, XMax);
  h->SetMarkerColor(GetNiceColor(ColorIndex));
  h->SetMarkerStyle(kFullCircle);
  h->SetLineColor(GetNiceColor(ColorIndex));
  Double_t minValue = 1e10;
  Double_t maxValue = -1;

  for(Int_t bin=1; bin<h->GetNbinsX()+1; bin++)
  {
    Double_t binContent = h->GetBinContent(bin);
    if(binContent<=0) continue; /// reject empty bins
    if(h->GetBinLowEdge(bin)<minXValue) continue; /// only in visable area
    if(h->GetBinLowEdge(bin)+h->GetBinWidth(bin)>maxXValue) continue;
    if(binContent>maxValue) maxValue = binContent;
    if(binContent<minValue) minValue = binContent;
  }
  h->GetYaxis()->SetRangeUser(0.95*minValue, 1.05*maxValue);
  if(minValue<minYValue) minYValue = minValue;
  if(maxValue>maxYValue) maxYValue = maxValue;
}

#endif
