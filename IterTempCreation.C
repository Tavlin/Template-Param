// #include "CommonHeader.h"
#include "chi2test.h"
#include "TFractionFitter.h"


Double_t mc_full_func1(Double_t *x,  Double_t *par){
  Double_t xx = x[0];
  return (Double_t) par[0]*mc_full_clone1->GetBinContent(mc_full->FindBin(xx)) +
  par[1]*korrBG_clone1->GetBinContent(korrBG->FindBin(xx));
}

Double_t mc_full_func42(Double_t *x,  Double_t *par){
  Double_t xx = x[0];
  return (Double_t) par[0]*mc_full_clone42->GetBinContent(mc_full->FindBin(xx)) +
  par[1]*korrBG_clone42->GetBinContent(korrBG->FindBin(xx));
}

Double_t mc_full_func2(Double_t *x,  Double_t *par){
  Double_t xx = x[0];
  return (Double_t) par[0]*mc_full_clone2->GetBinContent(mc_full->FindBin(xx)) +
  par[1]+par[2]*xx;
}

// Double_t PeakAKorrBG(Double_t *x,  Double_t *par){
//   Double_t xx = x[0];
//   return (Double_t) par[0]*korrBG_clone1->GetBinContent(korrBG->FindBin(xx));
// }
//
// Double_t mc_full_func42(Double_t *x,  Double_t *par){
//   Double_t xx = x[0];
//   return (Double_t) par[]*mc_full_clone42->GetBinContent(mc_full->FindBin(xx));
// }
//
// Double_t PeakAKorrBG42(Double_t *x,  Double_t *par){
//   Double_t xx = x[0];
//   return (Double_t) par[]*korrBG_clone42->GetBinContent(korrBG->FindBin(xx));
// }

void IterTempCreation(void){

  TString sPath = gDirectory->GetPath();

  TString str;
  Double_t int_error = 0;
  Double_t int_value = 0;
  const Int_t nbins = 45;
  const Int_t ndrawpoints = 1.e5;
  const int n_iter = 4;
  const int epsilon = 1.e-6;
  int iter = 0;
  int chi2_test_vari = 1;
  std::vector<Double_t> chi2_dt_iter;
  std::vector<Double_t> chi2_pol1_iter;
  std::vector<Double_t> chi2_dt_iter_test;
  // const int numberbins = 26;
  TFile *IterTemp;

  TH2D* testtest[numberbins];

  TH1F* hChi2_dt_iter[numberbins];
  TH1F* hChi2_pol1_iter[numberbins];
  // TH1F* hChi2_dt_iter_test[numberbins];
  TH1F* hYield_dt_uncorr = new TH1F("hYield_dt_uncorr","",numberbins, fBinsPi013TeVEMCPt);
  TH1F* hYield_pol1_uncorr = new TH1F("hYield_pol1_uncorr","",numberbins, fBinsPi013TeVEMCPt);
  // TF1* f1 = new TF1("f1", "1", 0.0, 0.3);
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
  TH1F* hDataCloneforCHi2;
  TH1F* hDoubleTemp;
  TH1F* hPol1;
  TH1F* data_clone_for_int_dt;
  TH1F* data_clone_for_int_pol1;

  SetHistoStandardSettings(hchi2_dt);
  SetHistoStandardSettings(hchi2_pol1);
  SetHistoStandardSettings(hpeakratio);
  SetHistoStandardSettings(hpeakcomp);
  SetHistoStandardSettings(ha_DT);
  SetHistoStandardSettings(hb_DT);
  SetHistoStandardSettings(ha_pol1);
  SetHistoStandardSettings(hYield_dt_uncorr);
  SetHistoStandardSettings(hYield_pol1_uncorr);
  hYield_dt_uncorr->SetLineColor(kTeal-7);
  hYield_dt_uncorr->SetMarkerColor(kTeal-7);
  hYield_pol1_uncorr->SetLineColor(kRed);
  hYield_pol1_uncorr->SetMarkerColor(kRed);


//////////////////////////////////////////////////////////////////////////////
  // going over all pt bins despite first one, which is some framework bs.
  for (int k = 1; k < numberbins; k++) {
    std::cout << "starte bin " << k << " reading and wrinting!" << std::endl << std::endl;

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
    // Getting the purposed corr Background
    korrBG = (TH1F*) data_MC->Clone("korrBG");
    korrBG->Add(mc_full,-1);


    ////////////////////////////////////////////////////////////////////////////
    // using the correlated BG from MC as Template
    TF1* fit_eq_double_temp = new TF1("fit_eq_double_temp", &mc_full_func1, 0.0,0.4, 2);
    fit_eq_double_temp->SetNpx(ndrawpoints);
    fit_eq_double_temp->SetNumberFitPoints(nbins);
    fit_eq_double_temp->SetLineColor(kTeal-7);
    fit_eq_double_temp->SetLineWidth(4);

    ////////////////////////////////////////////////////////////////////////////
    // second TF1 for drawing only since the other function gets scaled two
    // times since the external scaling of the templates and the internal
    // scaling
    TF1* fit_eq_double_temp42 = new TF1("fit_eq_double_temp42", &mc_full_func42, 0.0,0.4, 2);
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

    data->SetTitle(str);

    ////////////////////////////////////////////////////////////////////////////
    // normal lame pol 1 fit with template
    TF1* fit_eq_1 = new TF1("fit_eq_1", &mc_full_func2, 0.0, 0.4, 3);
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
    data->SetMarkerStyle(24);
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

    // clearing the vectors
    chi2_dt_iter.clear();
    chi2_pol1_iter.clear();
    chi2_dt_iter.resize(0);
    chi2_pol1_iter.resize(0);
    chi2_dt_iter_test.clear();
    chi2_dt_iter_test.resize(0);
    chi2_test_vari = 1;
    iter = 0;


    while(chi2_test_vari){
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

      // if(iter == 0){
      //   mc_full_clone1->Write(Form("mc_full_clone_beforeIterFit_bin%02d",k));
      //   korrBG_clone1->Write(Form("korrBG_clone_beforeIterFit_bin%02d",k));
      // }

      //////////////////////////////////////////////////////////////////////////
      // fit 2 temp
      TFitResultPtr r_double_temp1 = data_clone1->Fit("fit_eq_double_temp", "QM0PS","", lowerparamrange, upperparamrange);
      data_clone1->Fit("fit_eq_double_temp42", "QM0PS","", lowerparamrange, upperparamrange);
      chi2_dt_iter.push_back(r_double_temp1->Chi2() / r_double_temp1->Ndf());


      //////////////////////////////////////////////////////////////////////////
      //fit pol 1 + temp
      TFitResultPtr r_pol1_temp1 = data_clone2->Fit("fit_eq_1", "QM0PS","", lowerparamrange, upperparamrange);
      chi2_pol1_iter.push_back(r_pol1_temp1->Chi2() / r_pol1_temp1->Ndf());

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
      // test chi2 for monitoring
      // gDirectory->Cd(sPath.Data());
      // hDoubleTemp = (TH1F*) mc_full_clone1->Clone("hDoubleTemp");
      // hDoubleTemp->Add(korrBG_clone1);
      // hDataCloneforCHi2 = (TH1F*) data->Clone("hDataCloneforCHi2");
      // gDirectory = IterTemp;
      // hDoubleTemp->Write(Form("hDoubleTemp_bin%02d_iter%02d",k,iter));
      // gDirectory->Cd(sPath.Data());
      // for(int i = 0; i < 200; i++){
      //   if( i < 13 || i > 63){
      //     hDataCloneforCHi2->SetBinContent(i,0);
      //     hDataCloneforCHi2->SetBinError(i,0);
      //     hDoubleTemp->SetBinContent(i,0);
      //     hDoubleTemp->SetBinError(i,0);
      //   }
      //   ////////////////////////////////////////////////////////////////////////
      //   // desperate try to calc the correct correlated error
      //   else{
      //     hDoubleTemp->SetBinError(i,sqrt(pow((r_double_temp1->Parameter(0)*mc_full->GetBinError(i)),2.)
      //     + pow((r_double_temp1->Parameter(1)*korrBG->GetBinError(i)),2.)));
      //   // std::cout << "hDoubleTemp Error:" << hDoubleTemp->GetBinError(i) << std::endl << std::endl;
      //   }
      // }
      // chi2_dt_iter_test.push_back(hDataCloneforCHi2->Chi2Test(hDoubleTemp, "WW CHI2/NDF", 0));

      //////////////////////////////////////////////////////////////////////////
      // reset data_clone histos and then calculate their new errors

      gDirectory->Cd(sPath.Data());
      data_clone1 = (TH1F*) data->Clone("data_clone1");
      data_clone2 = (TH1F*) data->Clone("data_clone2");
      hDoubleTemp =( TH1F*) mc_full_clone1->Clone("hDoubleTemp");
      hDoubleTemp->Add(korrBG_clone1);
      for(int j = 13; j < 63; j++){
        data_clone1->SetBinError(j,sqrt(pow(data_clone1->GetBinError(j),2.)
        + pow((r_double_temp1->Parameter(0)*mc_full->GetBinError(j)),2.)
        + pow((r_double_temp1->Parameter(1)*korrBG->GetBinError(j)),2.)));

        data_clone2->SetBinError(j,sqrt(pow(data_clone2->GetBinError(j),2.)
        + pow((r_double_temp1->Parameter(0)*mc_full->GetBinError(j)),2.)));
      }
      //////////////////////////////////////////////////////////////////////////
      // Writing after the scaling and error calculation:

      // gDirectory = IterTemp;

      // data_clone1->Write(Form("data_clone_afterFitWithScalingAndErros_bin%02d_iter%d", k, iter));


        //////////////////////////////////////////////////////////////////////////
        // Writing after the scaling and error calculation:
        // data_clone4->Write(Form("data_clone_afterFitWithScalingAndErros_bin%02d_iter%d", k, iter));

      IterTemp->Close();
      if(iter >=1){
        if(fabs(chi2_dt_iter[iter-1]-chi2_dt_iter[iter] <= epsilon) &&
           fabs(chi2_pol1_iter[iter-1]-chi2_pol1_iter[iter] <= epsilon)){
            //  std::cout << "iter = " << iter << "chi2 = " << chi2_dt_iter[iter] << std::endl;
          //////////////////////////////////////////////////////////////////////
          // Writing after the scaling and error calculation:

          //////////////////////////////////////////////////////////////////////
          // for the last iteration step don't reset the clones, instead calc
          // errors for the data histos that will be used in the final Fit
          // afterwards

          for(int j = 0; j < 75; j++){
          data_clone4->SetBinError(j,sqrt(data_clone4->GetBinError(j) *
          data_clone4->GetBinError(j) + mc_full_clone1->GetBinError(j) *
          mc_full_clone1->GetBinError(j) + korrBG_clone1->GetBinError(j) *
          korrBG_clone1->GetBinError(j)));

          data_clone3->SetBinError(j,sqrt(data_clone3->GetBinError(j) *
          data_clone3->GetBinError(j) + mc_full_clone2->GetBinError(j) *
          mc_full_clone2->GetBinError(j)));

          chi2_test_vari = 0;
        }
      }
    }
      iter++;
    }
    ////////////////////////////////////////////////////////////////////////////
    // making the CHi2 monitoring histos!

    hChi2_dt_iter[k-1] = new TH1F(Form("hChi2_dt_iter_bin%02d",k-1),"",iter,0.5,(Double_t)iter+0.5);
    SetHistoStandardSettings(hChi2_dt_iter[k-1]);
    hChi2_pol1_iter[k-1] = new TH1F(Form("1hChi2_pol1_iter_bin%02d",k-1),"",iter,0.5,(Double_t)iter+0.5);
    SetHistoStandardSettings(hChi2_pol1_iter[k-1]);
    // hChi2_dt_iter_test[k-1] = new TH1F(Form("hChi2_dt_iter_test_bin%02d",k-1),"",iter,0.5,(Double_t)iter+0.5);
    // SetHistoStandardSettings(hChi2_dt_iter_test[k-1]);

    for(int i = 0; i < iter; i++){
      hChi2_dt_iter[k-1]->SetBinContent(i+1,chi2_dt_iter[i]);
      hChi2_pol1_iter[k-1]->SetBinContent(i+1,chi2_pol1_iter[i]);
      // hChi2_dt_iter_test[k-1]->SetBinContent(i+1,chi2_dt_iter_test[i]);
      hChi2_dt_iter[k-1]->SetYTitle("#chi^{2}/ndf");
      hChi2_dt_iter[k-1]->SetXTitle("Iterationstep");
      hChi2_dt_iter[k-1]->SetLineColor(kTeal-7);
      hChi2_dt_iter[k-1]->SetMarkerColor(kTeal-7);
      hChi2_pol1_iter[k-1]->SetYTitle("#chi^{2}/ndf");
      hChi2_pol1_iter[k-1]->SetXTitle("Iterationstep");
      hChi2_pol1_iter[k-1]->SetLineColor(kRed);
      hChi2_pol1_iter[k-1]->SetMarkerColor(kRed);
      // hChi2_dt_iter_test[k-1]->SetYTitle("#chi^{2}/ndf");
      // hChi2_dt_iter_test[k-1]->SetXTitle("Iterationstep");
      // hChi2_dt_iter_test[k-1]->SetLineColor(kBlue+2);
      // hChi2_dt_iter_test[k-1]->SetMarkerColor(kBlue+2);
    }
    gDirectory->Cd(sPath.Data());
    ///////////////////////////////////////////////////////////////////////////
    // final  2 temp fit
    mc_full_clone1 = (TH1F*) mc_full->Clone("mc_full_clone1");
    korrBG_clone1 = (TH1F*) korrBG->Clone("korrBG_clone1");



    TFitResultPtr r_double_temp = data_clone4->Fit("fit_eq_double_temp", "M0S","",lowerparamrange , upperparamrange);
    TH1F* mc_full_clone3 = (TH1F*) mc_full->Clone("mc_full_clone3");
    TH1F* korrBG_clone3 = (TH1F*) korrBG->Clone("korrBG_clone3");
    Double_t corrbackerror[101];
    for(int i = 0; i <= 100; i++){
      corrbackerror[i] = sqrt(pow(korrBG_clone3->GetBinError(i)
      *r_double_temp->Parameter(1),2.)+ pow(korrBG_clone3->GetBinContent(i)
      *r_double_temp->Error(1),2.));
    }
    mc_full_clone3->Scale(r_double_temp->Parameter(0));
    korrBG_clone3->Scale(r_double_temp->Parameter(1));

    ////////////////////////////////////////////////////////////////////////////
    //  making full histogram of the DoubleTemplate Param.
    hDoubleTemp = (TH1F*) mc_full_clone3->Clone("hDoubleTemp");
    hDoubleTemp->Add(korrBG_clone3);
    hDoubleTemp->SetMarkerStyle(20);
    hDoubleTemp->SetMarkerSize(1.5);
    hDoubleTemp->SetMarkerColor(kTeal-7);
    hDoubleTemp->SetLineColor(kTeal-7);
    ////////////////////////////////////////////////////////////////////////////

    for(int i = 0; i <= 100; i++){
      korrBG_clone3->SetBinError(i,corrbackerror[i]);
    }

    ///////////////////////////////////////////////////////////////////////////
    // final pol 1 + temp fit
    mc_full_clone2 = (TH1F*) mc_full->Clone("mc_full_clone2");
    TFitResultPtr r_pol1_temp = data_clone3->Fit("fit_eq_1", "M0S","", lowerparamrange, upperparamrange);
    TH1F* mc_full_clone4 = (TH1F*) mc_full->Clone("mc_full_clone4");
    mc_full_clone4->Scale(r_pol1_temp->Parameter(0));
    mc_full_clone4->SetLineColor(kRed);
    mc_full_clone4->SetMarkerColor(kRed);

    TF1* fpol1 = new TF1("fpol1", "[0]+x*[1]", 0.0, 0.4);
    fpol1->SetParameter(0, fit_eq_1->GetParameter(1));
    fpol1->SetParameter(1, fit_eq_1->GetParameter(2));
    fpol1->SetLineColor(kTeal-7);
    fpol1->SetLineWidth(3);

    ////////////////////////////////////////////////////////////////////////////
    //  making full histogram of the DoubleTemplate Param.
    hPol1 = (TH1F*) mc_full_clone4->Clone("hDoubleTemp");
    hPol1->Add(fpol1);
    hPol1->SetMarkerStyle(20);
    hPol1->SetMarkerSize(1.5);
    hPol1->SetMarkerColor(kRed+1);
    hPol1->SetLineColor(kRed+1);


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
    for (int i = 0; i < 101; i++) {
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
    hRatioDoubleTemp->SetYTitle("(data-param)/#sigma");
    hRatioPol1->SetYTitle("(data-param)/#sigma");
    hRatioDoubleTemp->GetYaxis()->SetTitleOffset(0.4);
    hRatioPol1->GetYaxis()->SetTitleOffset(0.4);
    hRatioDoubleTemp->SetMarkerStyle(24);
    hRatioDoubleTemp->SetMarkerSize(1.5);
    hRatioPol1->SetMarkerStyle(24);
    hRatioPol1->SetMarkerSize(1.5);

    hRatioDoubleTemp->SetLineColor(kTeal-7);
    hRatioDoubleTemp->SetMarkerColor(kTeal-7);
    hRatioPol1->SetLineColor(kRed);
    hRatioPol1->SetMarkerColor(kRed);

    ////////////////////////////////////////////////////////////////////////////
    // testitesti
    ////////////////////////////////////////////////////////////////////////////
    testtest[k-1] = chi2test(data, mc_full, korrBG);

    data->SetTitle(str);
    IterTemp = new TFile("IterTemp.root","UPDATE");
    gDirectory = IterTemp;
    data->Write(Form("data_bin%02d",k));
    data_clone3->Write(Form("data_addedErrosPol1_bin%02d",k));
    data_clone4->Write(Form("data_addedErrosDT_bin%02d",k));
    mc_full_clone4->Write(Form("mc_peak_pol1_bin%02d",k));
    mc_full_clone3->Write(Form("mc_full_DT_bin%02d", k));
    korrBG_clone3->Write(Form("korrBG_bin%02d", k));
    fpol1->Write(Form("fpol1_bin%02d",k));
    hRatioDoubleTemp->Write(Form("hRatioDoubleTemp_bin%02d", k));
    hRatioPol1->Write(Form("hRatioPol1_bin%02d", k));
    hDoubleTemp->Write(Form("hDoubleTemp_bin%02d",k));
    hPol1->Write(Form("hPol1_bin%02d",k));
    hChi2_dt_iter[k-1]->Write(Form("hChi2_dt_iter_bin%02d",k));
    hChi2_pol1_iter[k-1]->Write(Form("hChi2_pol1_iter_bin%02d",k));
    // hChi2_dt_iter_test[k-1]->Write(Form("hChi2_dt_iter_test_bin%02d",k));
    testtest[k-1]->Write(Form("hTestTest_bin%02d",k));

    gDirectory->Cd(sPath.Data());
    ////////////////////////////////////////////////////////////////////////////
    // getting the chi2 of the current pT bin
    hchi2_pol1->SetBinContent(k+1,r_pol1_temp->Chi2()/r_pol1_temp->Ndf());
    hchi2_dt->SetBinContent(k+1,r_double_temp->Chi2()/r_double_temp->Ndf());



    data_clone_for_int_dt = (TH1F*) data->Clone("hYield_dt_uncorr");
    data_clone_for_int_dt->Add(korrBG_clone3, -1);
    int_value = data_clone_for_int_dt->IntegralAndError(0, data_clone_for_int_dt->GetXaxis()->FindBin(0.3), int_error);
    hYield_dt_uncorr->SetBinContent(k+1, int_value);
    hYield_dt_uncorr->SetBinError(k+1, int_error);

    data_clone_for_int_pol1 = (TH1F*) data->Clone("hYield_pol1_uncorr");
    data_clone_for_int_pol1->Add(fpol1, -1);
    int_value = data_clone_for_int_pol1->IntegralAndError(0, data_clone_for_int_pol1->GetXaxis()->FindBin(0.3), int_error);
    hYield_pol1_uncorr->SetBinContent(k+1, int_value);
    hYield_pol1_uncorr->SetBinError(k+1, int_error);

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
    std::cout << "bin number" << k << "Ende" << std::endl << std::endl;
    delete hChi2_dt_iter[k-1];
    delete hChi2_pol1_iter[k-1];
    // delete hChi2_dt_iter_test[k-1];
    delete testtest[k-1];
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
  hYield_dt_uncorr->Scale(1,"width");
  hYield_pol1_uncorr->Scale(1,"width");
  hYield_dt_uncorr->SetYTitle(rawyield);
  hYield_dt_uncorr->SetXTitle(massaxis);
  hYield_pol1_uncorr->SetYTitle(rawyield);
  hYield_pol1_uncorr->SetXTitle(massaxis);
  hYield_dt_uncorr->Write("hYield_dt_uncorr");
  hYield_pol1_uncorr->Write("hYield_pol1_uncorr");

  delete hpeakratio;
  delete hpeakcomp;
  delete hchi2_dt;
  delete hchi2_pol1;
  delete hRatioDoubleTemp;
  delete hRatioPol1;
  delete hYield_dt_uncorr;
  delete hYield_pol1_uncorr;
  delete hDoubleTemp;
  delete hPol1;
  IterTemp->Close();

}
