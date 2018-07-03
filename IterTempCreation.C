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

void IterTempCreation(void){

  TString sPath = gDirectory->GetPath();

  TString str;
  const Int_t nbins = 45;
  const Int_t ndrawpoints = 1.e5;
  const int n_iter = 4;
  const int numberbins = 26;
  TFile *IterTemp;

  //////////////////////////////////////////////////////////////////////////////
  // setting up the 2 Huistograms to compare chi2 from the to fit methods as
  // well as peak factor comp. between pol 1 and double temp fit and the ratio
  // of BG. scaling factor and the Peak scaling factor

  TH1F* hchi2_dt = new TH1F("hchi2_dt", "", numberbins, fBinsPi013TeVEMCPt);
  TH1F* hchi2_pol1 = new TH1F("hchi2_pol1", "", numberbins, fBinsPi013TeVEMCPt);

  TH1F* ha_DT = new TH1F("ha_DT", "", numberbins, fBinsPi013TeVEMCPt);
  TH1F* hb_DT = new TH1F("hb_DT", "", numberbins, fBinsPi013TeVEMCPt);
  TH1F* ha_pol1 = new TH1F("ha_pol1", "", numberbins, fBinsPi013TeVEMCPt);


  TH1F* hpeakratio = new TH1F("hpeakratio", "", numberbins, fBinsPi013TeVEMCPt);
  TH1F* hpeakcomp = new TH1F("hpeakcomp", "", numberbins, fBinsPi013TeVEMCPt);
  TH1F* hRatioDoubleTemp;
  TH1F* hRatioPol1;
  SetHistoStandardSettings(hchi2_dt);
  SetHistoStandardSettings(hchi2_pol1);
  SetHistoStandardSettings(hpeakratio);
  SetHistoStandardSettings(hpeakcomp);
  SetHistoStandardSettings(ha_DT);
  SetHistoStandardSettings(hb_DT);
  SetHistoStandardSettings(ha_pol1);


//////////////////////////////////////////////////////////////////////////////
  // going over all pt bins despite first one, which is some framework bs.
  for (int k = 1; k < numberbins; k++) {
    cout << "starte bin " << k << " reading and wrinting!" << endl << endl;

    ////////////////////////////////////////////////////////////////////////////
    // open MC histo path
    TFile* MCFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_MC_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");
    if (MCFile->IsOpen() ) printf("MCFile opened successfully\n");

    ////////////////////////////////////////////////////////////////////////////
    // retrieve MC histograms
    data_MC = (TH1F*) MCFile->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02i",k));
    mc_full = (TH1F*) MCFile->Get(Form("Mapping_TrueFullMeson_InvMass_in_Pt_Bin%02i",k));


    TFile* DataFile = SafelyOpenRootfile("./00010113_1111112067032220000_01631031000000d0/13TeV/Pi0_data_GammaConvV1WithoutCorrection_00010113_1111112067032220000_01631031000000d0.root");
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

    //////////////////////////////////////////////////////////////////////////
    // Fix! Changes < in TLatex to #leq
    str = data_MC->GetTitle();
    TString str_copy = str.Copy();
    str_copy.ReplaceAll("<","#leq");
    str.Replace(0,20,str_copy,23);

    SetHistoStandardSettings(data_MC);
    SetHistoStandardSettings(mc_full);
    SetHistoStandardSettings(korrBG);
    SetHistoStandardSettings(korrBG);

    data_MC->SetTitle(str);

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
    data->SetYTitle("d#it{N}/d#it{M}_{#gamma#gamma} (#it{c}^{2}/GeV)");
    data->SetMarkerStyle(20);
    data->SetMarkerSize(1.5);
    data->SetTitle("");
    korrBG->SetLineColor(kTeal-7);
    korrBG->SetMarkerColor(kTeal-7);
    korrBG->SetMarkerStyle(34);
    korrBG->SetMarkerSize(2);
    mc_full->SetLineColor(kTeal-7);
    mc_full->SetMarkerColor(kTeal-7);
    mc_full->SetMarkerStyle(33);
    mc_full->SetMarkerSize(2.5);

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

      ////////////////////////////////////////////////////////////////////////////
      // creating the new root file to safe all the related histograms and fits
      // in it.
      if(k == 1 && iter == 0){
        IterTemp = new TFile("IterTemp.root", "RECREATE");
          // if (IterTemp->IsOpen() ) printf("IterTemp opened successfully\n");
      }

      else{
        IterTemp = new TFile("IterTemp.root", "UPDATE");
          // if (IterTemp->IsOpen() ) printf("IterTemp opened successfully\n");
      }

      if(iter == 0){
        // mc_full_clone1->Write(Form("mc_full_clone_beforeIterFit_bin%02d_iter%d",k, iter));
        // korrBG_clone1->Write(Form("korrBG_clone_beforeIterFit_bin%02d_iter%d",k, iter));
      }

      //////////////////////////////////////////////////////////////////////////
      // fit 2 temp
      TFitResultPtr r_double_temp1 = data_clone1->Fit("fit_eq_double_temp", "QM0PS","", 0.054, 0.252);
      data_clone1->Fit("fit_eq_double_temp42", "QM0PS","", 0.054, 0.252);


      //////////////////////////////////////////////////////////////////////////
      //fit pol 1 + temp
      TFitResultPtr r_pol1_temp1 = data_clone2->Fit("fit_eq_1", "QM0PS","", 0.054, 0.252);

      //////////////////////////////////////////////////////////////////////////
      // scale 2 temp
      mc_full_clone1->Scale(r_double_temp1->Parameter(0));
      korrBG_clone1->Scale(r_double_temp1->Parameter(1));


      //////////////////////////////////////////////////////////////////////////
      // Wrinting after the scaling:
      // mc_full_clone1->Write(Form("mc_full_clone_afterScaling_bin%02d_iter%d", k, iter));
      // korrBG_clone1->Write(Form("korrBG_clone_afterScaling_bin%02d_iter%d", k, iter));



      //////////////////////////////////////////////////////////////////////////
      // scale for pol 1 + temp
      mc_full_clone2->Scale(r_pol1_temp1->Parameter(0));

      //////////////////////////////////////////////////////////////////////////
      // reset data_clone histos and then calculate their new errors
      if(iter < n_iter){
        gDirectory->Cd(sPath.Data());
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
        // Writing after the scaling and error calculation:
        gDirectory = IterTemp;
        // data_clone1->Write(Form("data_clone_afterFitWithScalingAndErros_bin%02d_iter%d", k, iter));

      }

      //////////////////////////////////////////////////////////////////////////
      // Writing after the scaling and error calculation:

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
        // Writing after the scaling and error calculation:
        // data_clone4->Write(Form("data_clone_afterFitWithScalingAndErros_bin%02d_iter%d", k, iter));

      }
      IterTemp->Close();
    }
    gDirectory->Cd(sPath.Data());
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


    ///////////////////////////////////////////////////////////////////////////
    // Picture of double template fit with chi2 and factors as well as ratio of
    // the factors for the template.
    data_clone4->SetLineWidth(3);
    data_clone3->SetLineWidth(3);
    data_clone3->SetLineColor(kGray+2);
    data_clone3->SetMarkerColor(kGray+2);
    data_clone4->SetTitle("");
    fpol1->SetLineColor(kRed);

    hRatioDoubleTemp = (TH1F*) data_clone4->Clone("RatioDoubleTemp");
    hRatioPol1 = (TH1F*) data_clone3->Clone("RatioPol1");
    hRatioDoubleTemp->Add(fit_eq_double_temp,-1.);
    hRatioPol1->Add(fit_eq_1, -1.);
    for (int i = 0; i < 75; i++) {
      hRatioDoubleTemp->SetBinContent(i,hRatioDoubleTemp->GetBinContent(i)/hRatioDoubleTemp->GetBinError(i));
      hRatioPol1->SetBinContent(i,hRatioPol1->GetBinContent(i)/hRatioPol1->GetBinError(i));
      hRatioDoubleTemp->SetBinError(i,0);
      hRatioPol1->SetBinError(i,0);

    }
    hRatioDoubleTemp->GetYaxis()->SetRangeUser(-5.,5.);
    hRatioDoubleTemp->GetYaxis()->SetNdivisions(509);
    hRatioPol1->GetYaxis()->SetRangeUser(-4.5,4.5);
    SetHistoStandardSettings(hRatioDoubleTemp,0.,0.,0.09);
    SetHistoStandardSettings(hRatioPol1,0.,0.,0.09);
    hRatioDoubleTemp->SetYTitle("(data-param)/#sigma(data)");
    hRatioPol1->SetYTitle("(data-param)/#sigma(data)");
    hRatioDoubleTemp->GetYaxis()->SetTitleOffset(0.4);
    hRatioPol1->GetYaxis()->SetTitleOffset(0.4);

    hRatioDoubleTemp->SetLineColor(kTeal-7);
    hRatioDoubleTemp->SetMarkerColor(kTeal-7);
    hRatioPol1->SetLineColor(kRed);
    hRatioPol1->SetMarkerColor(kRed);

    
    IterTemp = new TFile("IterTemp.root","UPDATE");
    gDirectory = IterTemp;
    data->Write(Form("data_bin%02d",k));
    data_clone4->Write(Form("data_addedErrosPol1_bin%02d",k));
    data_clone3->Write(Form("data_addedErrosDT_bin%02d",k));
    mc_full_clone4->Write(Form("mc_peak_pol1_bin%02d",k));
    mc_full_clone3->Write(Form("mc_full_DT_bin%02d", k));
    korrBG_clone3->Write(Form("korrBG_bin%02d", k));
    fpol1->Write(Form("fpol1_bin%02d",k));
    hRatioDoubleTemp->Write(Form("hRatioDoubleTemp_bin%02d", k));
    hRatioPol1->Write(Form("hRatioPol1_bin%02d", k));


    gDirectory->Cd(sPath.Data());
    ////////////////////////////////////////////////////////////////////////////
    // getting the chi2 of the current pT bin
    hchi2_pol1->SetBinContent(k+1,r_pol1_temp->Chi2()/r_pol1_temp->Ndf());
    hchi2_dt->SetBinContent(k+1,r_double_temp->Chi2()/r_double_temp->Ndf());

    ////////////////////////////////////////////////////////////////////////////
    // getting the peakratio of the current pT bin
    Double_t peakratio, peakratioerr, peakscale, bgscale, peakerr, bgerr;
    Double_t peakscale_pol1, peakerr_pol1, peakratio_to_pol1, peakratioerr_to_pol1;
    peakscale = r_double_temp->Parameter(0);
    bgscale = r_double_temp->Parameter(1);
    peakerr = r_double_temp->Error(0);
    bgerr = r_double_temp->Error(1);

    peakscale_pol1 = r_pol1_temp->Parameter(0);
    peakerr_pol1 = r_pol1_temp->Error(0);

    peakratio = bgscale/peakscale;
    peakratio_to_pol1 = peakscale/peakscale_pol1;

    peakratioerr = sqrt(pow(bgerr/peakscale, 2.)+pow(bgscale*peakerr/pow(peakscale, 2.), 2.));
    peakratioerr_to_pol1 = sqrt(pow(peakerr/peakscale_pol1,2.0)+pow(peakscale*peakerr_pol1/pow(peakscale_pol1,2.0),2.0));

    hpeakratio->SetBinContent(k+1, peakratio);
    hpeakratio->SetBinError(k+1, peakratioerr);

    hpeakcomp->SetBinContent(k+1, peakratio_to_pol1);
    hpeakcomp->SetBinError(k+1, peakratioerr_to_pol1);

    ha_DT->SetBinContent(k+1, peakscale);
    ha_DT->SetBinError(k+1, peakerr);
    hb_DT->SetBinContent(k+1, bgscale);
    hb_DT->SetBinError(k+1, bgerr);
    ha_pol1->SetBinContent(k+1, peakscale_pol1);
    ha_pol1->SetBinError(k+1, peakerr_pol1);



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
    delete fpol1;
    delete fit_eq_double_temp;
    delete fit_eq_1;

    if(k < numberbins-1){
      // IterTemp->Close();
      MCFile->Close();
      DataFile->Close();
    }
    IterTemp->Close();
    cout << "bin number" << k << "Ende" << endl << endl;
  }

  /////////////////////////////////////////////////////////////////////////////
  // Writing of Chi2(pT)
  IterTemp = new TFile("IterTemp.root","UPDATE");

  hchi2_pol1->SetLineColor(kRed);
  hchi2_pol1->SetLineWidth(3);
  hchi2_dt->SetLineWidth(3);
  hpeakratio->SetLineWidth(3);
  hpeakcomp->SetLineWidth(3);
  hchi2_pol1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hchi2_dt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hchi2_pol1->GetYaxis()->SetTitle("#chi^{2}/ndf");
  hchi2_dt->GetYaxis()->SetTitle("#chi^{2}/ndf");
  hpeakratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hpeakratio->GetYaxis()->SetTitle("b_{double}/a_{double}");
  hpeakcomp->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hpeakcomp->GetYaxis()->SetTitle("a_{pol1}/ a_{double}");
  hchi2_dt->GetYaxis()->SetRangeUser(0.0,5.0);


  hchi2_dt->Write(Form("hchi2_dt"));
  hchi2_pol1->Write(Form("hchi2_pol1"));
  hpeakratio->Write("hpeakratio");
  hpeakcomp->Write("hpeakcomp");
  ha_DT->Write("DoubleTemplatePeakFactor");
  hb_DT->Write("DoubleTemplatecorrBGFactor");
  ha_pol1->Write("Pol1PeakFactor");

  delete hpeakratio;
  delete hpeakcomp;
  delete hchi2_dt;
  delete hchi2_pol1;
  delete hRatioDoubleTemp;
  delete hRatioPol1;

  IterTemp->Close();

}
