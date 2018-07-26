#include "CommonHeader.h"
#include "TFractionFitter.h"


Double_t mc_full_func1(Double_t x){
  return mc_full_clone1->GetBinContent(mc_full->FindBin(x));
}

Double_t mc_full_func2(Double_t x){
  return mc_full_clone2->GetBinContent(mc_full->FindBin(x));
}

Double_t PeakAKorrBG(Double_t x){
  return korrBG_clone1->GetBinContent(korrBG->FindBin(x));
}

Double_t mc_full_func42(Double_t x){
  return mc_full_clone42->GetBinContent(mc_full->FindBin(x));
}

Double_t PeakAKorrBG42(Double_t x){
  return korrBG_clone42->GetBinContent(korrBG->FindBin(x));
}

void Iter_Temp_Fit(void){

  TString str;
  const Int_t nbins = 45;
  const Int_t ndrawpoints = 1.e5;
  const int n_iter = 4;
  TString doubletempstring = "Double template param.";
  TString pol1string = "Peak template + 1^{st} ord. pol. param.";

  //////////////////////////////////////////////////////////////////////////////
  // setting up the canvas to draw on. Will later be changed for the chi2 pic
  TCanvas *c1 = new TCanvas("c1","",1200,1000);
  c1->cd();
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.09);
  c1->SetRightMargin(0.02);
  c1->SetLeftMargin(0.09);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);

  TCanvas *canInvMass = new TCanvas("canInvMass","",1200,1200);
  TPad *pad1InvMass = new TPad("pad1InvMass","",0.0,0.33,1.0,1.0);
  pad1InvMass->SetTopMargin(0.05);
  pad1InvMass->SetLeftMargin(0.09);
  pad1InvMass->SetBottomMargin(0.0);
  pad1InvMass->SetRightMargin(0.02);
  TPad *pad2InvMass = new TPad("pad2InvMass","",0.0,0.0,1.0,0.33);
  pad2InvMass->SetTopMargin(0.0);
  pad2InvMass->SetLeftMargin(0.09);
  pad2InvMass->SetBottomMargin(0.3);
  pad2InvMass->SetRightMargin(0.02);


  //////////////////////////////////////////////////////////////////////////////
  // setting up the 2 Huistograms to compare chi2 from the to fit methods as
  // well as peak factor comp. between pol 1 and double temp fit and the ratio
  // of BG. scaling factor and the Peak scaling factor

  TH1F* hchi2_dt = new TH1F("hchi2_dt", "", 26, fBinsPi013TeVEMCPt);
  TH1F* hchi2_pol1 = new TH1F("hchi2_pol1", "", 26, fBinsPi013TeVEMCPt);
  TH1F* hpeakratio = new TH1F("hpeakratio", "", 26, fBinsPi013TeVEMCPt);
  TH1F* hpeakcomp = new TH1F("hpeakcomp", "", 26, fBinsPi013TeVEMCPt);
  TH1F* hRatioDoubleTemp;
  TH1F* hRatioPol1;
  TH1F* hDoubleTemp;
  SetHistoStandardSettings(hchi2_dt);
  SetHistoStandardSettings(hchi2_pol1);
  SetHistoStandardSettings(hpeakratio);
  SetHistoStandardSettings(hpeakcomp);

  chi_and_param42 = new TLatex();
  SetLatexSettings(chi_and_param42);

  //////////////////////////////////////////////////////////////////////////////
  // going over all pt bins despite first one, which is some framework bs.
  for (int k = 1; k < 26; k++) {

    ////////////////////////////////////////////////////////////////////////////
    // open MC histo path
    TFile* MCFile = new TFile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_MC_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root", "READ");
    if (MCFile->IsOpen() ) printf("MCFile opened successfully\n");

    ////////////////////////////////////////////////////////////////////////////
    // retrieve MC histograms
    data_MC = (TH1F*) MCFile->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02i",k));
    mc_full = (TH1F*) MCFile->Get(Form("Mapping_TrueFullMeson_InvMass_in_Pt_Bin%02i",k));


    TFile* DataFile = new TFile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root", "READ");
    if (DataFile->IsOpen() ) printf("DataFile opened successfully\n");
    data = (TH1F*) DataFile->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02i",k));


    mc_full->GetXaxis()->SetRangeUser(0.,0.4);
    data->GetXaxis()->SetRangeUser(0.,0.4);
    ////////////////////////////////////////////////////////////////////////////
    // Getinng the purposed corr Background
    korrBG = (TH1F*) data_MC->Clone("korrBG");
    korrBG->Add(mc_full,-1);


    ////////////////////////////////////////////////////////////////////////////
    // using the correlated BG from MC as Template
    TF1* fit_eq_double_temp = new TF1("fit_eq_double_temp", "PeakAKorrBG(x)*[1] + mc_full_func1(x)*[0]", 0.0,0.4);
    fit_eq_double_temp->SetNpx(ndrawpoints);
    fit_eq_double_temp->SetNumberFitPoints(nbins);
    fit_eq_double_temp->SetLineColor(kTeal-7);
    fit_eq_double_temp->SetLineWidth(4);

    ////////////////////////////////////////////////////////////////////////////
    // second TF1 for drawing only since the other function gets scaled two
    // times since the external scaling of the templates and the internal
    // scaling
    TF1* fit_eq_double_temp42 = new TF1("fit_eq_double_temp42", "PeakAKorrBG42(x)*[1] + mc_full_func42(x)*[0]", 0.0,0.4);
    fit_eq_double_temp42->SetNpx(ndrawpoints);
    fit_eq_double_temp42->SetNumberFitPoints(nbins);
    fit_eq_double_temp42->SetLineColor(kTeal-7);
    fit_eq_double_temp42->SetLineWidth(4);

    ///////////////////////////////////////////////////////////////////////////
    // Fix! Changes < in TLatex to #leq
    str = data_MC->GetTitle();
    TString str_copy = str.Copy();
    str_copy.ReplaceAll("<","#leq");
    str.Replace(0,20,str_copy,23);
    cout << str << endl;

    SetHistoStandardSettings(data_MC);
    SetHistoStandardSettings(mc_full);
    SetHistoStandardSettings(korrBG);
    SetHistoStandardSettings(korrBG);

    ////////////////////////////////////////////////////////////////////////////
    // normal lame pol 1 fit with template
    TF1* fit_eq_1 = new TF1("fit_eq_1", "mc_full_func2(x)*[0]+[2]+x*[3]",0.0,0.4);
    fit_eq_1->SetNpx(ndrawpoints);
    fit_eq_1->SetNumberFitPoints(nbins);
    fit_eq_1->SetLineColor(kRed);
    fit_eq_1->SetLineWidth(4);

    //////////////////////////////////////////////////////////////////////////
    // making things look good
    // change y title offset to fit and look nicely
    data->GetYaxis()->SetTitleOffset(1.2);
    data->GetYaxis()->SetLabelOffset(0.006);
    data->SetTitleSize(0.03, "xy");
    data->SetLabelSize(0.03, "xy");
    data->GetYaxis()->SetRangeUser(
      1.5*data->GetBinContent(data->GetMinimumBin()),
      1.1*data->GetBinContent(data->GetMaximumBin()));
    data->SetYTitle("d#it{N}/d#it{M}_{#gamma#gamma} (#it{c}^{2}/GeV)");
    data->SetMarkerStyle(20);
    data->SetMarkerSize(1.5);
    data->SetTitle("");
    data->SetLineWidth(3);
    korrBG->SetLineColor(kCyan+3);
    korrBG->SetMarkerColor(kCyan+3);
    korrBG->SetMarkerStyle(21);
    korrBG->SetMarkerSize(1.5);
    mc_full->SetLineColor(kGreen+3);
    mc_full->SetMarkerColor(kGreen+3);
    mc_full->SetMarkerStyle(33);
    mc_full->SetMarkerSize(2);

    //////////////////////////////////////////////////////////////////////////
    // clone for 2 temp fit
    TH1F* data_clone1 = (TH1F*) data->Clone("data_clone1");
    TH1F* data_clone4 = (TH1F*) data->Clone("data_clone4");

    //////////////////////////////////////////////////////////////////////////
    //clone for pol 1 fit
    TH1F* data_clone2 = (TH1F*) data->Clone("data_clone2");
    TH1F* data_clone3 = (TH1F*) data->Clone("data_clone3");



    for (int iter = 0; iter <= n_iter; iter++) {

      //////////////////////////////////////////////////////////////////////////
      //clone for 2 temp fit
      mc_full_clone1 = (TH1F*) mc_full->Clone("mc_full_clone1");
      mc_full_clone42 = (TH1F*) mc_full->Clone("mc_full_clone42");
      korrBG_clone1 = (TH1F*) korrBG->Clone("korrBG_clone1");
      korrBG_clone42 = (TH1F*) korrBG->Clone("korrBG_clone42");

      //////////////////////////////////////////////////////////////////////////
      //clone for pol 1 fit
      mc_full_clone2 = (TH1F*) mc_full->Clone("mc_full_clone2");

      //////////////////////////////////////////////////////////////////////////
      // StartDrawing before something happend

      legiter = new TLegend(0.5,0.43,0.9,0.63); // was 0.55
      SetLegendSettigns(legiter);
      legiter->AddEntry(data_clone1, "data: (signal + corr. back.)", "lp");
      legiter->AddEntry((TObject*) 0x0, "template:", "");
      legiter->AddEntry(mc_full_clone1, "signal", "p");
      legiter->AddEntry(korrBG_clone1, "corr. back. ", "p");



      c1->cd();
      c1->Clear();
      c1->Update();

      data_clone4->Draw();
      mc_full_clone1->Draw("same");
      korrBG_clone1->Draw("same");
      DrawLabelALICE(0.5, 0.9, 0.02, 0.03, str);
      c1->Update();

      ////////////////////////////////////////////////////////////////////////////
      // creating TLine which represents the range in which the fit will be made
      double line_y = gPad->GetUymax()*0.995;
      cout << line_y << endl;
      TLine* fitrange_in_iter = new TLine(0.054,line_y,0.252,line_y);

      fitrange_in_iter->SetLineColor(kAzure+10);
      fitrange_in_iter->SetLineWidth(7);
      fitrange_in_iter->Draw("same");
      legiter->AddEntry(fitrange_in_iter, "Param. range", "l");
      legiter->Draw();

      c1->Update();
      c1->SaveAs(Form("IterationProgress/bin%02i_%02iBeforeFit.png",k,01));
      c1->Clear();

      //////////////////////////////////////////////////////////////////////////
      // fit 2 temp
      TFitResultPtr r_double_temp1 = data_clone1->Fit("fit_eq_double_temp", "QM0PS","", 0.054, 0.252);
      data_clone1->Fit("fit_eq_double_temp42", "QM0PS","", 0.054, 0.252);

      //////////////////////////////////////////////////////////////////////////
      //fit pol 1 + temp
      TFitResultPtr r_pol1_temp1 = data_clone2->Fit("fit_eq_1", "QM0PS","", 0.054, 0.252);

      c1->Clear();
      c1->Update();
      data_clone4->Draw();
      //////////////////////////////////////////////////////////////////////////
      // scale 2 temp
      mc_full_clone1->Scale(r_double_temp1->Parameter(0));
      korrBG_clone1->Scale(r_double_temp1->Parameter(1));


      //////////////////////////////////////////////////////////////////////////
      // Drawing after the scaling:
      hDoubleTemp = (TH1F*) mc_full_clone1->Clone("hDoubleTemp");
      hDoubleTemp->Add(korrBG_clone1);
      hDoubleTemp->SetLineColor(kTeal-7);
      hDoubleTemp->SetMarkerColor(kTeal-7);
      hDoubleTemp->SetMarkerStyle(20);
      hDoubleTemp->SetMarkerSize(1.5);

      TLegend* legiter2 = new TLegend(0.5,0.43,0.9,0.63);
      SetLegendSettigns(legiter2);
      legiter2->AddEntry(data_clone1, "data: (signal + corr. back.)", "lp");
      legiter2->AddEntry((TObject*) 0x0, "template:", "");
      legiter2->AddEntry(mc_full_clone1, "signal", "p");
      legiter2->AddEntry(korrBG_clone1, "corr. back. ", "p");
      legiter2->AddEntry(hDoubleTemp, "signal + corr. back.", "p");
      legiter2->AddEntry(fitrange_in_iter, "Param. range", "l");

      c1->Update();
      hDoubleTemp->Draw("same");
      mc_full_clone1->Draw("same");
      korrBG_clone1->Draw("same");
      fitrange_in_iter->Draw("same");
      drawchi_and_param42(chi_and_param42, r_double_temp1);
      legiter2->Draw();
      DrawLabelALICE(0.5, 0.9, 0.02, 0.03, str);

      c1->Update();
      c1->SaveAs(Form("IterationProgress/bin%02i_%02iAfterFitWithScaling_iter%02i.png",k,02, iter));
      c1->Clear();

      //////////////////////////////////////////////////////////////////////////
      // scale for pol 1 + temp
      mc_full_clone2->Scale(r_pol1_temp1->Parameter(0));

      //////////////////////////////////////////////////////////////////////////
      // reset data_clone histos and then calculate their new errors
      if(iter < n_iter){
        data_clone1 = (TH1F*) data->Clone("data_clone1");
        data_clone2 = (TH1F*) data->Clone("data_clone2");
        for(int j = 0; j < 75; j++){
          data_clone1->SetBinError(j,sqrt(data_clone1->GetBinError(j) *
          data_clone1->GetBinError(j) + mc_full_clone1->GetBinError(j) *
          mc_full_clone1->GetBinError(j) + korrBG_clone1->GetBinError(j) *
          korrBG_clone1->GetBinError(j)));

          data_clone2->SetBinError(j,sqrt(data_clone2->GetBinError(j) *
          data_clone2->GetBinError(j) + mc_full_clone2->GetBinError(j) *
          mc_full_clone2->GetBinError(j)));
        }
        //////////////////////////////////////////////////////////////////////////
        // Drawing after the scaling and error calculation:

        data_clone1->Draw();
        mc_full_clone1->Draw("same");
        korrBG_clone1->Draw("same");
        hDoubleTemp->Draw("same");
        fitrange_in_iter->Draw("same");
        drawchi_and_param42(chi_and_param42, r_double_temp1);
        legiter2->Draw();
        DrawLabelALICE(0.5, 0.9, 0.02, 0.03, str);

        c1->Update();
        c1->SaveAs(Form("IterationProgress/bin%02i_%02iAfterFitWithScalingAndErros_iter%02i.png",k,03, iter));
        c1->Clear();
      }

      //////////////////////////////////////////////////////////////////////////
      // Drawing after the scaling and error calculation:

      //////////////////////////////////////////////////////////////////////////
      // for the last iteration step don't reset the clones, instead calc errors
      // for the data histos that will be used in the final Fit afterwards
      if(iter == n_iter){
        for(int j = 0; j < 75; j++){
          data_clone4->SetBinError(j,sqrt(data_clone4->GetBinError(j) *
          data_clone4->GetBinError(j) + mc_full_clone1->GetBinError(j) *
          mc_full_clone1->GetBinError(j) + korrBG_clone1->GetBinError(j) *
          korrBG_clone1->GetBinError(j)));

          data_clone3->SetBinError(j,sqrt(data_clone3->GetBinError(j) *
          data_clone3->GetBinError(j) + mc_full_clone2->GetBinError(j) *
          mc_full_clone2->GetBinError(j)));
        }

        //////////////////////////////////////////////////////////////////////////
        // Drawing after the scaling and error calculation:

        data_clone4->Draw();
        mc_full_clone1->Draw("same");
        korrBG_clone1->Draw("same");
        hDoubleTemp->Draw("same");
        fitrange_in_iter->Draw("same");
        drawchi_and_param42(chi_and_param42, r_double_temp1);
        legiter2->Draw();
        DrawLabelALICE(0.5, 0.9, 0.02, 0.03, str);

        c1->Update();
        c1->SaveAs(Form("IterationProgress/bin%02i_%02iAfterFitWithScalingAndErros_iter%02i.png",k,03, iter));
        c1->Clear();
      }
      delete fitrange_in_iter;
      delete legiter2;
    }

    ///////////////////////////////////////////////////////////////////////////
    // final  2 temp fit
    mc_full_clone1 = (TH1F*) mc_full->Clone("mc_full_clone1");
    korrBG_clone1 = (TH1F*) korrBG->Clone("korrBG_clone1");
    TFitResultPtr r_double_temp = data_clone4->Fit("fit_eq_double_temp", "M0PS","", 0.054, 0.252);
    TH1F* mc_full_clone3 = (TH1F*) mc_full->Clone("mc_full_clone3");
    TH1F* korrBG_clone3 = (TH1F*) korrBG->Clone("korrBG_clone3");
    mc_full_clone3->Scale(r_double_temp->Parameter(0));
    korrBG_clone3->Scale(r_double_temp->Parameter(1));

    ///////////////////////////////////////////////////////////////////////////
    // final pol 1 + temp fit
    mc_full_clone2 = (TH1F*) mc_full->Clone("mc_full_clone2");
    TFitResultPtr r_pol1_temp = data_clone3->Fit("fit_eq_1", "M0PS","", 0.054, 0.252);
    TH1F* mc_full_clone4 = (TH1F*) mc_full->Clone("mc_full_clone4");
    mc_full_clone4->Scale(r_pol1_temp->Parameter(0));
    mc_full_clone4->SetLineColor(kRed);
    mc_full_clone4->SetMarkerColor(kRed);

    TF1* fpol1 = new TF1("fpol1", "[0]+x*[1]", 0.0, 0.3);
    fpol1->SetParameter(0, fit_eq_1->GetParameter(2));
    fpol1->SetParameter(1, fit_eq_1->GetParameter(3));
    fpol1->SetLineColor(kTeal-7);
    fpol1->SetLineWidth(3);


    // ///////////////////////////////////////////////////////////////////////////
    // // create TLatex for extra information  (pT range)
    // TLatex* tex = new TLatex();
    // data_MC->SetTitle("");
    // tex->SetTextSize(0.03);
    // tex->SetTextFont(42);
    //
    // TLegend *leg = new TLegend(0.52,0.65,0.8,0.85); // was 0.55
    // SetLegendSettigns(leg);
    // leg->AddEntry(data, "Data same - scaled mixed evt.", "lp");
    // leg->AddEntry(fit_eq_double_temp, doubletempstring, "l");
    // c1->Clear();
    //
    // ///////////////////////////////////////////////////////////////////////////
    // // Picture of double template fit with chi2 and factors as well as ratio of
    // // the factors for the template.
    // data_clone4->SetLineWidth(3);
    // data_clone3->SetLineWidth(3);
    // data_clone3->SetLineColor(kGray+2);
    // data_clone3->SetMarkerColor(kGray+2);
    // data_clone4->SetTitle("");
    // data_clone4->Draw("p");
    //
    // fit_eq_double_temp->Draw("same");
    //
    // c1->Update();
    //
    // ////////////////////////////////////////////////////////////////////////////
    // // creating TLine which represents the range in which the fit will be made
    // double line_y = gPad->GetUymax()*0.995;
    // cout << line_y << endl;
    // TLine* fitrange = new TLine(0.054,line_y,0.252,line_y);
    //
    // fitrange->SetLineColor(kAzure+10);
    // fitrange->SetLineWidth(7);
    // fitrange->Draw("same");
    //
    // leg->AddEntry(fitrange, "Param. range", "l");
    // leg->Draw("same");
    // // tex->DrawLatexNDC(0.35,0.9, str);
    //
    // TLatex* chi_and_param = new TLatex();
    // SetLatexSettings(chi_and_param);
    // chi_and_param->DrawLatexNDC(0.18,0.85,
    // Form("#frac{#chi^{2}_{double temp}}{ndf} = %.2lf ",r_double_temp->Chi2() /
    // r_double_temp->Ndf()));
    //
    // // chi_and_param->DrawLatexNDC(0.18,0.75,
    // // Form("a_{double} = %.2lf ",r_double_temp->Parameter(0)));
    // //
    // // chi_and_param->DrawLatexNDC(0.18,0.70,
    // // Form("b_{double} = %.2lf ",r_double_temp->Parameter(1)));
    // //
    // // chi_and_param->DrawLatexNDC(0.18,0.65,
    // // Form("b_{double}/a_{double} = %.2lf ",r_double_temp->Parameter(1) /
    // // r_double_temp->Parameter(0)));
    // DrawLabelALICE(0.55, 0.5, 0.02, 0.03, str);
    //
    //
    // c1->SaveAs(Form("MCTemplatesAnData/DataFitWithMCIter%02i.png",k));
    // c1->Clear();
    //
    // fpol1->SetLineColor(kRed);
    // TLegend *leg2 = new TLegend(0.52,0.52,0.8,0.87); // was 0.55
    // SetLegendSettigns(leg2);
    // leg2->SetTextSize(0.06);
    // leg2->AddEntry(data_clone4, "Data same - scaled mixed evt.", "lp");
    // leg2->AddEntry(fit_eq_double_temp, doubletempstring, "l");
    // leg2->AddEntry(fit_eq_1, pol1string);
    //
    // c1->Clear();
    //
    // ///////////////////////////////////////////////////////////////////////////
    // // Comparison Picture between 1st Order polynomial fit and double template
    // // fit to data.
    // // Also drawn are the chi2 of the 2 fits and the factors coming from the
    // // double template fit and their ratio
    //
    // canInvMass->cd();
    // pad1InvMass->Draw();
    // pad2InvMass->Draw("same");
    // pad1InvMass->cd();
    //
    // SetHistoStandardSettings(data_clone4, 0., 0., 0.03*3./2.);
    // data_clone4->GetYaxis()->SetRangeUser(
    //   1.4*data_clone4->GetBinContent(data_clone4->GetMinimumBin()),
    //   1.4*data_clone4->GetBinContent(data_clone4->GetMaximumBin()));
    //
    // data_clone4->SetTitle("");
    // data_clone4->SetLineColor(kBlue+2);
    // data_clone4->SetMarkerColor(kBlue+2);
    // data_clone4->GetYaxis()->SetTitleOffset(0.9);
    // data_clone4->Draw("p");
    // data_clone3->Draw("samep");
    // fit_eq_double_temp->Draw("same");
    // fit_eq_1->Draw("same");
    //
    // canInvMass->Update();
    // line_y = gPad->GetUymax()*0.995;
    // cout << line_y << endl;
    // TLine* fitrange2 = new TLine(0.054,line_y,0.252,line_y);
    // fitrange2->SetLineColor(kAzure+10);
    // fitrange2->SetLineWidth(7);
    // fitrange2->Draw("same");
    // leg2->AddEntry(fitrange2, "Param. range", "l");
    // SetLegendSettigns(leg2, 0.03*3./2.);
    // leg2->Draw("same");
    // // tex->DrawLatexNDC(0.35,0.9, str);
    //
    // TLatex* chi_and_param2 = new TLatex();
    // SetLatexSettings(chi_and_param2, 0.03);
    // chi_and_param2->DrawLatexNDC(0.18,0.67,
    // Form("#frac{#chi^{2}_{pol1 + temp}}{ndf} = %.2lf ",r_pol1_temp->Chi2() /
    // r_pol1_temp->Ndf()));
    //
    // chi_and_param2->DrawLatexNDC(0.18,0.83,
    // Form("#frac{#chi^{2}_{double temp}}{ndf} = %.2lf ",r_double_temp->Chi2() /
    // r_double_temp->Ndf()));
    //
    // // chi_and_param2->DrawLatexNDC(0.19,0.55,
    // // Form("a_{double} = %.2lf",r_double_temp->Parameter(0)));
    // //
    // // chi_and_param2->DrawLatexNDC(0.19,0.45,
    // // Form("a_{pol1} = %.2lf",r_pol1_temp->Parameter(0)));
    // DrawLabelALICE(0.55, 0.48, 0.04, 0.03*3./2., str);
    // pad1InvMass->Update();
    //
    //
    // pad2InvMass->cd();
    // hRatioDoubleTemp = (TH1F*) data_clone4->Clone("RatioDoubleTemp");
    // hRatioPol1 = (TH1F*) data_clone3->Clone("RatioPol1");
    // hRatioDoubleTemp->Add(fit_eq_double_temp,-1.);
    // hRatioPol1->Add(fit_eq_1, -1.);
    // for (int i = 0; i < 75; i++) {
    //   hRatioDoubleTemp->SetBinContent(i,hRatioDoubleTemp->GetBinContent(i)/hRatioDoubleTemp->GetBinError(i));
    //   hRatioPol1->SetBinContent(i,hRatioPol1->GetBinContent(i)/hRatioPol1->GetBinError(i));
    //   hRatioDoubleTemp->SetBinError(i,0);
    //   hRatioPol1->SetBinError(i,0);
    //
    // }
    // hRatioDoubleTemp->GetYaxis()->SetRangeUser(-5.,5.);
    // hRatioDoubleTemp->GetYaxis()->SetNdivisions(509);
    // hRatioPol1->GetYaxis()->SetRangeUser(-4.5,4.5);
    // SetHistoStandardSettings(hRatioDoubleTemp,0.,0.,0.09);
    // SetHistoStandardSettings(hRatioPol1,0.,0.,0.09);
    // hRatioDoubleTemp->SetYTitle("(data-param)/#sigma(data)");
    // hRatioPol1->SetYTitle("(data-param)/#sigma(data)");
    // hRatioDoubleTemp->GetYaxis()->SetTitleOffset(0.4);
    // hRatioPol1->GetYaxis()->SetTitleOffset(0.4);
    //
    // hRatioDoubleTemp->SetLineColor(kTeal-7);
    // hRatioDoubleTemp->SetMarkerColor(kTeal-7);
    // hRatioPol1->SetLineColor(kRed);
    // hRatioPol1->SetMarkerColor(kRed);
    //
    // TLine* line_0 = new TLine(0.0, 0.0, 0.4, 0.0);
    // line_0->SetLineWidth(3);
    // line_0->SetLineStyle(1);
    // line_0->SetLineColor(kBlack);
    // TLine* line_p1 = new TLine(0.0, 1.0, 0.4, 1.0);
    // line_p1->SetLineWidth(2);
    // line_p1->SetLineStyle(3);
    // line_p1->SetLineColor(kGray+2);
    // TLine* line_m1 = new TLine(0.0, -1.0, 0.4, -1.0);
    // line_m1->SetLineWidth(2);
    // line_m1->SetLineStyle(3);
    // line_m1->SetLineColor(kGray+2);
    // TLine* line_p3 = new TLine(0.0, 3.0, 0.4, 3.0);
    // line_p3->SetLineWidth(2);
    // line_p3->SetLineStyle(2);
    // line_p3->SetLineColor(kGray+2);
    // TLine* line_m3 = new TLine(0.0, -3.0, 0.4, -3.0);
    // line_m3->SetLineWidth(2);
    // line_m3->SetLineStyle(2);
    // line_m3->SetLineColor(kGray+2);
    //
    //
    //
    // hRatioDoubleTemp->DrawCopy("");
    // line_0->Draw("same");
    // line_p1->Draw("same");
    // line_m1->Draw("same");
    // line_p3->Draw("same");
    // line_m3->Draw("same");
    // hRatioDoubleTemp->DrawCopy("same");
    // hRatioPol1->DrawCopy("same");
    // pad2InvMass->Update();
    //
    // canInvMass->Update();
    // canInvMass->SaveAs(Form("MCTemplatesAnData/DataFitWithMCCompIter%02i.png",k));
    // canInvMass->Clear("D");
    //
    //
    //
    //
    //
    //
    // c1->cd();
    // data_clone4->SetLineColor(kGray+3);
    // data_clone4->SetMarkerColor(kGray+3);
    //
    // ///////////////////////////////////////////////////////////////////////////
    // // drawing Double Template fit to data with extra Error from iterartion
    // c1->Update();
    // TLatex* tex_double_temp = new TLatex();
    // SetLatexSettings(tex_double_temp);
    //
    // TLegend *leg3 = new TLegend(0.52,0.65,0.8,0.85); // was 0.55
    // SetLegendSettigns(leg3);
    // leg3->AddEntry(data_clone4, "Data same - scaled mixed evt.", "lp");
    // leg3->AddEntry(mc_full_clone1, "scaled peak template", "p");
    // leg3->AddEntry(korrBG_clone1, "scaled corr. BG template ", "p");
    //
    // SetHistoStandardSettings(data_clone4);
    // data_clone4->Draw("p");
    // mc_full_clone3->Draw("same");
    // korrBG_clone3->Draw("same");
    // leg3->Draw("same");
    // fitrange2->Draw("same");
    // // tex->DrawLatexNDC(0.35,0.9, str);
    // tex_double_temp->DrawLatexNDC(0.15, 0.85, "Double template parametrization");
    // DrawLabelALICE(0.55, 0.5, 0.02, 0.03, str);
    //
    // c1->Update();
    // c1->SaveAs(Form("MCTemplatesAnData/DataDoubleTempComp%02i.png",k));
    // c1->Clear();
    //
    // ///////////////////////////////////////////////////////////////////////////
    // // drawing 1st order polynomial plus peak scaled to data with extra error
    // // from iterration
    // c1->Update();
    // TLatex* tex_pol1 = new TLatex();
    // SetLatexSettings(tex_pol1);
    //
    // TLegend *leg4 = new TLegend(0.52,0.65,0.8,0.85); // was 0.55
    // SetLegendSettigns(leg4);
    // leg4->AddEntry(data_clone4, "Data same - scaled mixed evt.", "lp");
    // leg4->AddEntry(mc_full_clone4, "scaled peak template", "p");
    // leg4->AddEntry(fpol1, "1^{st} ord. pol. ", "l");
    // leg4->AddEntry(fit_eq_1, pol1string);
    // leg4->AddEntry(fitrange2, "Param. range", "l");
    //
    // data_clone3->SetTitle("");
    // data_clone3->Draw("p");
    // mc_full_clone4->Draw("same");
    // fpol1->Draw("same");
    // fit_eq_1->Draw("same");
    // leg4->Draw("same");
    // fitrange2->Draw("same");
    // // tex->DrawLatexNDC(0.35,0.9, str);
    // tex_pol1->DrawLatexNDC(0.15, 0.85, "#splitline{Peak template +}{1^{st} ord. pol. param}");
    // DrawLabelALICE(0.55, 0.5, 0.02, 0.03, str);
    //
    // c1->Update();
    // c1->SaveAs(Form("MCTemplatesAnData/DataPol1Comp%02i.png",k));
    // c1->Clear();
    //
    // ////////////////////////////////////////////////////////////////////////////
    // // Drawing double temp fit + data with normal errors:
    // c1->Update();
    //
    //
    // TLegend *leg5 = new TLegend(0.52,0.65,0.8,0.85); // was 0.55
    // SetLegendSettigns(leg5);
    // leg5->AddEntry(data, "Data same - scaled mixed evt.", "lp");
    // leg5->AddEntry(mc_full_clone3, "scaled peak template", "p");
    // leg5->AddEntry(korrBG_clone3, "scaled corr. BG. template", "p");
    // leg5->AddEntry(fit_eq_double_temp, doubletempstring, "l");
    // leg5->AddEntry(fitrange2, "Param. range", "l");
    //
    // data->Draw("");
    // mc_full_clone3->Draw("same");
    // korrBG_clone3->Draw("same");
    // fit_eq_double_temp->Draw("same");
    // leg5->Draw("same");
    // fitrange2->Draw("same");
    // // tex->DrawLatexNDC(0.35,0.9, str);
    // chi_and_param2->SetTextSize(0.03);
    // chi_and_param2->DrawLatexNDC(0.15,0.85,
    // Form("#frac{#chi^{2}_{double temp}}{ndf} = %.2lf ",r_double_temp->Chi2() /
    // r_double_temp->Ndf()));
    //
    // // chi_and_param2->DrawLatexNDC(0.17,0.70,
    // // Form("a_{double} = %.2lf ",r_double_temp->Parameter(0)));
    // //
    // // chi_and_param2->DrawLatexNDC(0.17,0.65,
    // // Form("b_{double} = %.2lf ",r_double_temp->Parameter(1)));
    // //
    // // chi_and_param2->DrawLatexNDC(0.17,0.60,
    // // Form("b_{double}/a_{double} = %.2lf ",r_double_temp->Parameter(1) /
    // // r_double_temp->Parameter(0)));
    // DrawLabelALICE(0.55, 0.5, 0.02, 0.03, str);
    //
    // c1->Update();
    // c1->SaveAs(Form("MCTemplatesAnData/DataDoubleTempFit_NormalErrors%02i.png",k));
    // c1->Clear();
    //
    // ////////////////////////////////////////////////////////////////////////////
    // // Drawing double temp fit + data with iter. errors:
    // c1->Update();
    //
    //
    // TLegend *leg6 = new TLegend(0.52,0.65,0.8,0.85); // was 0.55
    // SetLegendSettigns(leg6);
    // leg6->AddEntry(data_clone4, "Data same - scaled mixed evt.", "lp");
    // leg6->AddEntry(mc_full_clone3, "scaled peak template", "p");
    // leg6->AddEntry(korrBG_clone3, "scaled corr. BG. template ", "p");
    // leg6->AddEntry(fit_eq_double_temp, doubletempstring, "l");
    // leg6->AddEntry(fitrange2, "Param. range", "l");
    //
    // data_clone4->Draw("");
    // data_clone4->GetYaxis()->SetTitleSize(0.03);
    // data_clone4->GetXaxis()->SetTitleSize(0.03);
    // mc_full_clone3->Draw("same");
    // korrBG_clone3->Draw("same");
    // fit_eq_double_temp->Draw("same");
    // leg6->Draw("same");
    // fitrange2->Draw("same");
    // // tex->DrawLatexNDC(0.35, 0.9, str);
    // tex->DrawLatexNDC(0.15,0.85, "Added MC-template");
    // tex->DrawLatexNDC(0.15,0.8, "uncertainty to Data");
    // tex->DrawLatexNDC(0.15, 0.75, Form("#frac{#chi^{2}}{ndf} = %.2lf", r_double_temp->Chi2()/r_double_temp->Ndf()));
    // DrawLabelALICE(0.55, 0.5, 0.02, 0.03, str);
    //
    // c1->Update();
    // c1->SaveAs(Form("MCTemplatesAnData/DataDoubleTempFit_AddedErrors%02i.png",k));
    // c1->Clear();
    //
    //
    // ////////////////////////////////////////////////////////////////////////////
    // // Drawing first ord. pol. fit + data with normal errors:
    // c1->Update();
    //
    // mc_full_clone4->SetLineColor(kRed);
    // mc_full_clone4->SetMarkerColor(kRed);
    // TLegend *leg7 = new TLegend(0.52,0.65,0.8,0.85); // was 0.55
    // SetLegendSettigns(leg7);
    // leg7->AddEntry(data, "Data same - scaled mixed evt.", "lp");
    // leg7->AddEntry(mc_full_clone4, "scaled peak template", "p");
    // leg7->AddEntry(fpol1, "1^{st} ord. pol.", "l");
    // leg7->AddEntry(fit_eq_1, pol1string);
    // leg7->AddEntry(fitrange2, "Param. range", "l");
    //
    // data->Draw("");
    // mc_full_clone4->Draw("same");
    // fpol1->Draw("same");
    // fit_eq_1->Draw("same");
    // leg7->Draw("same");
    // fitrange2->Draw("same");
    // // tex->DrawLatexNDC(0.35,0.9, str);
    //
    // chi_and_param2->DrawLatexNDC(0.15,0.85,
    // Form("#frac{#chi^{2}_{double temp}}{ndf} = %.2lf ",r_pol1_temp->Chi2() /
    // r_pol1_temp->Ndf()));
    //
    // // chi_and_param2->DrawLatexNDC(0.60,0.55,
    // // Form("a_{pol1} = %.2lf ",r_pol1_temp->Parameter(0)));
    // DrawLabelALICE(0.55, 0.5, 0.02, 0.03, str);
    //
    //
    // c1->Update();
    // c1->SaveAs(Form("MCTemplatesAnData/DataPol1Fit_NormalErrors%02i.png",k));
    // c1->Clear();
    //
    // ////////////////////////////////////////////////////////////////////////////
    // // Drawing first ord. pol. fit + data with iter. errors:
    // c1->Update();
    //
    //
    // TLegend *leg8 = new TLegend(0.52,0.65,0.8,0.85); // was 0.55
    // SetLegendSettigns(leg8);
    // leg8->AddEntry(data_clone4, "Data same - scaled mixed evt.", "lp");
    // leg8->AddEntry(mc_full_clone4, "scaled peak template", "p");
    // leg8->AddEntry(fpol1, "1^{st} ord. pol.", "l");
    // leg8->AddEntry(fit_eq_1, pol1string);
    // leg8->AddEntry(fitrange2, "Param. range", "l");
    //
    // data_clone4->Draw("");
    // mc_full_clone4->Draw("same");
    // fpol1->Draw("same");
    // fit_eq_1->Draw("same");
    // leg8->Draw("same");
    // fitrange2->Draw("same");
    // // tex->DrawLatexNDC(0.35,0.9, str);
    // tex->DrawLatexNDC(0.15,0.85, "Added MC-template");
    // tex->DrawLatexNDC(0.15,0.8, "uncertainty to Data");
    // tex->DrawLatexNDC(0.15, 0.75, Form("#frac{#chi^{2}}{ndf} = %.2lf", r_pol1_temp->Chi2()/r_pol1_temp->Ndf()));
    // DrawLabelALICE(0.55, 0.5, 0.02, 0.03, str);
    //
    // c1->Update();
    // c1->SaveAs(Form("MCTemplatesAnData/DataPol1Fit_AddedErrors%02i.png",k));
    // c1->Clear();
    //
    // ////////////////////////////////////////////////////////////////////////////
    // // Drawing both corr. BG versions to Data with normal errors
    // c1->Update();
    //
    //
    // TLegend *leg9 = new TLegend(0.52,0.65,0.8,0.85); // was 0.55
    // SetLegendSettigns(leg9);
    // leg9->AddEntry(data, "Data same - scaled mixed evt.", "lp");
    // leg9->AddEntry(korrBG_clone3, "scaled corr. BG. template", "p");
    // leg9->AddEntry(fpol1, "1^{st} ord. pol.", "l");
    // leg9->AddEntry(fitrange2, "Param. range", "l");
    //
    // data->Draw("");
    // data->GetYaxis()->SetTitleOffset(1.5);
    // korrBG_clone3->Draw("same");
    // fpol1->SetLineColor(kRed);
    // fpol1->Draw("same");
    // leg9->Draw("same");
    // fitrange2->Draw("same");
    // // tex->DrawLatexNDC(0.35,0.9, str);
    // tex->DrawLatexNDC(0.15,0.8, "Corr. BG. comparison");
    // DrawLabelALICE(0.55, 0.5, 0.02, 0.03, str);
    //
    //
    // c1->Update();
    // c1->SaveAs(Form("MCTemplatesAnData/CorrBGComp%02i.png",k));
    // fpol1->SetLineColor(kTeal-7);
    // c1->Clear();
    //
    // ////////////////////////////////////////////////////////////////////////////
    // // getting the chi2 of the current pT bin
    // hchi2_pol1->SetBinContent(k+1,r_pol1_temp->Chi2()/r_pol1_temp->Ndf());
    // hchi2_dt->SetBinContent(k+1,r_double_temp->Chi2()/r_double_temp->Ndf());
    //
    // ////////////////////////////////////////////////////////////////////////////
    // // getting the peakratio of the current pT bin
    // Double_t peakratio, peakratioerr, peakscale, bgscale, peakerr, bgerr;
    // Double_t peakscale_pol1, peakerr_pol1, peakratio_to_pol1, peakratioerr_to_pol1;
    // peakscale = r_double_temp->Parameter(0);
    // bgscale = r_double_temp->Parameter(1);
    // peakerr = r_double_temp->Error(0);
    // bgerr = r_double_temp->Error(1);
    //
    // peakscale_pol1 = r_pol1_temp->Parameter(0);
    // peakerr_pol1 = r_pol1_temp->Error(0);
    //
    // peakratio = bgscale/peakscale;
    // peakratio_to_pol1 = peakscale/peakscale_pol1;
    //
    // peakratioerr = sqrt(pow(bgerr/peakscale, 2.)+pow(bgscale*peakerr/pow(peakscale, 2.), 2.));
    // peakratioerr_to_pol1 = sqrt(pow(peakerr/peakscale_pol1,2.0)+pow(peakscale*peakerr_pol1/pow(peakscale_pol1,2.0),2.0));
    //
    // hpeakratio->SetBinContent(k+1, peakratio);
    // hpeakratio->SetBinError(k+1, peakratioerr);
    //
    // hpeakcomp->SetBinContent(k+1, peakratio_to_pol1);
    // hpeakcomp->SetBinError(k+1, peakratioerr_to_pol1);



    ////////////////////////////////////////////////////////////////////////////
    // garbage collection since cint does not have one?
    delete mc_full_clone1;
    delete mc_full_clone2;
    delete mc_full_clone3;
    delete mc_full_clone4;
    delete korrBG_clone1;
    delete korrBG_clone3;
    delete data_clone1;
    delete data_clone2;
    delete data_clone3;
    delete data_clone4;
    delete mc_full;
    delete korrBG;
    delete data;
    delete hDoubleTemp;
    // delete leg;
    // delete leg2;
    // delete leg3;
    // delete leg4;
    // delete leg5;
    // delete leg6;
    // delete leg7;
    // delete leg8;
    // delete leg9;
    // delete chi_and_param;
    // delete chi_and_param2;
    // delete tex;
    // delete fitrange;
    delete fpol1;
    delete fit_eq_double_temp;
    delete fit_eq_1;
    // delete fitrange2;
    // delete tex_double_temp;
    // delete tex_pol1;

    if(k < 25){
      MCFile->Close();
      DataFile->Close();
    }



  }

  // /////////////////////////////////////////////////////////////////////////////
  // // Drawing of Chi2(pT)
  // TLine* oneline = new TLine(0.0, 1., 21.0, 1.0);
  // oneline->SetLineWidth(3);
  // oneline->SetLineStyle(2);
  // oneline->SetLineColor(kGray+3);
  //
  // c1->Update();
  //
  //
  // hchi2_pol1->SetLineColor(kRed);
  // hchi2_pol1->SetLineWidth(3);
  // hchi2_dt->SetLineWidth(3);
  // hpeakratio->SetLineWidth(3);
  // hpeakcomp->SetLineWidth(3);
  // hchi2_pol1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  // hchi2_dt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  // hchi2_pol1->GetYaxis()->SetTitle("#chi^{2}/ndf");
  // hchi2_dt->GetYaxis()->SetTitle("#chi^{2}/ndf");
  // hpeakratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  // hpeakratio->GetYaxis()->SetTitle("b_{double}/a_{double}");
  // hpeakcomp->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  // hpeakcomp->GetYaxis()->SetTitle("a_{pol1}/ a_{double}");
  // hchi2_dt->GetYaxis()->SetRangeUser(0.0,5.0);
  //
  // TLegend* legchi = new TLegend(0.35,0.55,0.9,0.65);
  // SetLegendSettigns(legchi);
  // legchi->AddEntry(hchi2_dt, doubletempstring, "l");
  // legchi->AddEntry(hchi2_pol1, pol1string, "l");
  //
  // // c1->SetTopMargin(0.05);
  // // c1->SetBottomMargin(0.14);
  // // c1->SetRightMargin(0.005);
  // // c1->SetLeftMargin(0.14);
  //
  //
  // hchi2_dt->Draw("");
  // hchi2_pol1->Draw("same");
  // legchi->Draw("same");
  // oneline->Draw("same");
  // DrawLabelALICE(0.2, 0.9, 0.018, 0.03);
  //
  // c1->SaveAs(Form("MCTemplatesAnData/Chi2.png"));
  // c1->Clear();
  //
  // c1->Update();
  // hpeakratio->Draw("");
  // oneline->Draw("same");
  // DrawLabelALICE(0.34, 0.9, 0.018, 0.03);
  // c1->SaveAs(Form("MCTemplatesAnData/Peakratio.png"));
  // c1->Clear();
  //
  //
  // c1->Update();
  // hpeakcomp->Draw("");
  // oneline->Draw("same");
  // DrawLabelALICE(0.13, 0.9, 0.018, 0.03);
  // c1->SaveAs(Form("MCTemplatesAnData/Peakcomp.png"));
  // c1->Clear();


  delete hpeakratio;
  delete hpeakcomp;
  delete hchi2_dt;
  delete hchi2_pol1;
  delete legiter;
  delete chi_and_param42;
  delete hRatioDoubleTemp;
  delete hRatioPol1;
  delete pad1InvMass;
  delete pad2InvMass;
  delete canInvMass;
  delete c1;
  // delete oneline;
  // delete legchi;

}
