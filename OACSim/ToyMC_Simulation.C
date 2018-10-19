#include "CommonHeader.h"

/*
ToyMC simulating particle decay with and without Opening Angle Cut (OAC). Also
builds ratio from sim with OAC divided without OAC. This will be used to
calculate the scaling factor for the correlated background templates to scale
their lower m_inv bins.
*/
void ToyMC_Simulation(int PID){
  gRandom->SetSeed(PID);
  //////////////////////////////////////////////////////////////////////////////
  //                  Declaration of global variables                         //
  //                                                                          //
  // 2D-Histograms for minv vs pT
  TH2D* hMinv_pT_wo_OAC = new TH2D("hMinv_pT_wo_OAC","", 300, 0.0, 0.3, 39, fBinsPi013TeVEMCPt);
  SetHistoStandardSettings2(hMinv_pT_wo_OAC);
  TH2D* hMinv_pT_w_OAC  = new TH2D("hMinv_pT_w_OAC", "", 300, 0.0, 0.3, 39, fBinsPi013TeVEMCPt);
  SetHistoStandardSettings2(hMinv_pT_w_OAC);
  TH2D* hMinv_pT_ratio  = new TH2D("hMinv_pT_ratio", "", 300, 0.0, 0.3, 39, fBinsPi013TeVEMCPt);
  SetHistoStandardSettings2(hMinv_pT_ratio);
  //                                                                          //
  //                                                                          //
  // TFile where the simulations result will be saved
  TFile* SaveFile = new TFile(Form("OAC_ToyMC_PID%d.root",PID), "RECREATE");
  //                                                                          //
  //                                                                          //
  // Bool variable to check wether the sim was already made without the OAC or
  // not
  Bool_t bWO_OAC_done = false;
  //                                                                          //
  //                                                                          //
  // Int variable to count how many times the histograms were filled to have
  // evenly filled. The max value is const. and stored in maxfillcount.
  Int_t fillcount = 0;
  Int_t loopcount = 0;
  Int_t count_WO_OAC =0;
  Int_t count_OAC =0;
  const Int_t maxfillcount = 2.e7;
  //                                                                          //
  //                                                                          //
  // OAC min angle in rad
  const Double_t OA = 0.018;
  //                                                                          //
  //                                                                          //
  // Simulation loop where two photons with random momenta are generated
  while(fillcount < maxfillcount){

    loopcount++;
    printProgress(((Double_t)fillcount)/((Double_t)maxfillcount));

    Double_t m = gRandom->Uniform(1.e-9, 0.3);
    Double_t pi = TMath::Pi();

    // set pT, rapidity, ...
    Double_t pt_lab = gRandom->Uniform(1.e-9 ,20.);
    Double_t phi_lab = gRandom->Uniform(2.*pi);
    Double_t y_lab = gRandom->Uniform(-0.8 ,0.8);

    Double_t mt_lab = sqrt(m*m + pt_lab*pt_lab);
    Double_t e_lab = mt_lab * cosh(y_lab);
    Double_t px_lab = pt_lab * cos(phi_lab);
    Double_t py_lab = pt_lab * sin(phi_lab);
    Double_t pz_lab = mt_lab * sinh(y_lab);
    Double_t p_lab  = sqrt(e_lab*e_lab - m*m);
    Double_t theta_lab = atan2(pt_lab,pz_lab);

    // draw cosine of the decay angle in the CMS from uniform distribution
    Double_t cos_theta_star = gRandom->Uniform(0.,1.);
    Double_t phi_star = gRandom->Uniform(2.*pi);
    Double_t sin_theta_star = sqrt(1. - pow(cos_theta_star,2.0));

    // rest frame of the pi0
    Double_t p1x_star = m/2. * sin_theta_star * cos(phi_star);
    Double_t p1y_star = m/2. * sin_theta_star * sin(phi_star);
    Double_t p1z_star = m/2. * cos_theta_star;

    Double_t p2x_star = - p1x_star;
    Double_t p2y_star = - p1y_star;
    Double_t p2z_star = - p1z_star;

    Double_t beta = 1./sqrt(1+(m*m)/(p_lab*p_lab));
    Double_t gamma = 1./sqrt(1.-beta*beta);

    // Lorentz transform of the momentum vectors of the two decay photons
    Double_t p1x = p1x_star;
    Double_t p1y = p1y_star;
    Double_t p1z = gamma*(p1z_star + beta * m/2.);

    Double_t p2x = p2x_star;
    Double_t p2y = p2y_star;
    Double_t p2z = gamma*(p2z_star + beta * m/2.);

    // 3D rotation to lab system
    Double_t p1xrot, p1yrot, p1zrot;
    Double_t p2xrot, p2yrot, p2zrot;

    RotateToLabSystem(theta_lab,phi_lab,p1x,p1y,p1z,p1xrot,p1yrot,p1zrot);
    RotateToLabSystem(theta_lab,phi_lab,p2x,p2y,p2z,p2xrot,p2yrot,p2zrot);

    p1x = p1xrot;
    p1y = p1yrot;
    p1z = p1zrot;

    p2x = p2xrot;
    p2y = p2yrot;
    p2z = p2zrot;

    // momentum/position uncertainty
    fSmear(p1x, p1y, p1z);
    fSmear(p2x, p2y, p2z);

    Double_t p1  = fCalcP(p1x, p1y, p1z);
    Double_t p2  = fCalcP(p2x, p2y, p2z);

    // only accept photons with 0.7 < E (GeV/c) < 20.
    if(fabs(p1) > 0.7 && fabs(p2) > 0.7 && fabs(p1) < 20. && fabs(p2) < 20.){
      Double_t theta = fCalcTheta(p1x, p1y, p1z, p2x, p2y, p2z);
      m = fCalcInvMass(p1, p2, theta);
      if(m > 0.0 && m <= 0.3){
        if(bWO_OAC_done == false){
          Double_t pT    = fCalcPT(p1x, p1y, p2x, p2y);
          hMinv_pT_wo_OAC->Fill(m, pT);
          fillcount++;
          count_WO_OAC++;

          // change to filling the histogram with OAC
          if(fillcount >= maxfillcount/2.){
            bWO_OAC_done = true;
          }
        }
        // only if we are about to fill the histo with OAC
        else{
          Double_t pT    = fCalcPT(p1x, p1y, p2x, p2y);
          fillcount++;
          count_OAC++;
          if(theta > 0.017){
            hMinv_pT_w_OAC->Fill(m, pT);
          }
        }
      }
    }
  }
  std::cout << '\n' << "loopcount - fillcount = " << loopcount - fillcount << '\n';
  std::cout << "count_WO_OAC = " << count_WO_OAC << '\n';
  std::cout << "count_OAC = " << count_OAC << '\n';
  // hMinv_pT_ratio = (TH2D*) hMinv_pT_w_OAC->Clone("hMinv_pT_ratio");
  // hMinv_pT_ratio->Divide(hMinv_pT_wo_OAC);
  hMinv_pT_ratio->Divide(hMinv_pT_w_OAC, hMinv_pT_wo_OAC, 1, 1, "");
  hMinv_pT_ratio->Write("hMinv_pT_ratio");
  delete hMinv_pT_wo_OAC;
  delete hMinv_pT_w_OAC;
  delete hMinv_pT_ratio;
  SaveFile->Close();


}
