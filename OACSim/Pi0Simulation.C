#include <CommenHeader.h>

void Pi0Simulation(TString AddName = "") {

  // Wir definieren ein Canvas auf das wir malen können
  TCanvas *cExtractSignal = new TCanvas("cExtractSignal","",800,800);
  // Wir stellen ein paar grundlegende Settings ein
  SetCanvasStandardSettings(cExtractSignal);// (diese Funktion ist in ExtractSignal.h definiert)

  Float_t m = 0.135; // pi0 mass

  // generate a certain number of pi0
  const Int_t Npi0 = 100000;

  // pT distribution
  TF1* fpt = new TF1("fpt","x*exp(-x/0.2)",0.,10.);

  // rapidity distribution
  TF1* fy = new TF1("fy", "gaus", -0.5, 0.5);
  fy->SetParameters(1., 0., 4.);

  // histograms for generated and accepted pi0's
  TH1F* hNPi0_gen_pt = new TH1F("hNPi0_gen_pt","generated pi0 pT spectrum",20,0.,10.);
  SetHistoStandardSettings(hNPi0_gen_pt);
  TH1F* hNPi0_gen_minv = new TH1F("hNPi0_gen_minv","generated pi0 minv spectrum",100,0.,0.3);
  SetHistoStandardSettings(hNPi0_gen_minv);
  TH2F* hNPi0_gen_minv_pt = new TH2F("hNPi0_gen_minv_pt","generated pi0: minv vs. pT",100,0.,0.3,20,0.,10.);
  SetHistoStandardSettings2(hNPi0_gen_minv_pt);

  // pi0 accepted by VCal (eta coverage |eta| < 0.5)
  TH1F* hNPi0_acc = new TH1F("hNPi0_acc","accepted pi0 pT spectrum",20,0.,10.);
  SetHistoStandardSettings(hNPi0_acc);

  for (int ip=0; ip < Npi0; ip++) {
    printProgress( ((Double_t)ip) / ((Double_t)Npi0) );

    Float_t pi = TMath::Pi();

    // set pT, rapidity, ...
    Float_t pt_lab = fpt->GetRandom();
    Float_t phi_lab = gRandom->Uniform(2.*pi);
    Float_t y_lab = fy->GetRandom();

    Float_t mt_lab = sqrt(m*m + pt_lab*pt_lab);
    Float_t e_lab = mt_lab * cosh(y_lab);
    Float_t px_lab = pt_lab * cos(phi_lab);
    Float_t py_lab = pt_lab * sin(phi_lab);
    Float_t pz_lab = mt_lab * sinh(y_lab);
    Float_t p_lab  = sqrt(e_lab*e_lab - m*m);
    Float_t theta_lab = atan2(pt_lab,pz_lab);

    // draw cosine of the decay angle in the CMS from uniform distribution
    Float_t cos_theta_star = gRandom->Uniform(0.,1.);
    Float_t phi_star = gRandom->Uniform(2.*pi);
    Float_t sin_theta_star = sqrt(1. - pow(cos_theta_star,2.0));

    // rest frame of the pi0
    Float_t p1x_star = m/2. * sin_theta_star * cos(phi_star);
    Float_t p1y_star = m/2. * sin_theta_star * sin(phi_star);
    Float_t p1z_star = m/2. * cos_theta_star;

    Float_t p2x_star = - p1x_star;
    Float_t p2y_star = - p1y_star;
    Float_t p2z_star = - p1z_star;

    Float_t beta = 1./sqrt(1+(m*m)/(p_lab*p_lab));
    Float_t gamma = 1./sqrt(1.-beta*beta);

    // Lorentz transform of the momentum vectors of the two decay photons
    Float_t p1x = p1x_star;
    Float_t p1y = p1y_star;
    Float_t p1z = gamma*(p1z_star + beta * m/2.);

    Float_t p2x = p2x_star;
    Float_t p2y = p2y_star;
    Float_t p2z = gamma*(p2z_star + beta * m/2.);

    // 3D rotation to lab system
    Float_t p1xrot, p1yrot, p1zrot;
    Float_t p2xrot, p2yrot, p2zrot;

    RotateToLabSystem(theta_lab,phi_lab,p1x,p1y,p1z,p1xrot,p1yrot,p1zrot);
    RotateToLabSystem(theta_lab,phi_lab,p2x,p2y,p2z,p2xrot,p2yrot,p2zrot);

    p1x = p1xrot;
    p1y = p1yrot;
    p1z = p1zrot;

    p2x = p2xrot;
    p2y = p2yrot;
    p2z = p2zrot;

    // calculate phi and eta for the two decay photons
    Float_t pt1 = sqrt(p1x*p1x + p1y*p1y);
    Float_t theta1 = atan2(pt1,p1z);
    Float_t eta1 = -log(tan(theta1/2.));

    Float_t pt2 = sqrt(p2x*p2x + p2y*p2y);
    Float_t theta2 = atan2(pt2,p2z);
    Float_t eta2 = -log(tan(theta2/2.));

///////////////////////////////////////////////

    // Berechne die invariante Masse
    // fCalcInvMass im header
    // Float_t minv = fCalcInvMass(p1x,p1y,p1z,p2x,p2y,p2z);

///////////////////////////////////////////////

    // Implementiere eine Energieverschmierung
    // und berechne den Transversalimpuls neu
    // pt_lab = fCalcPT(p1x,p1y,p2x,p2y);

///////////////////////////////////////////////

    // Implementiere eine endliche Detektorakzeptanz in phi

///////////////////////////////////////////////

    hNPi0_gen_pt->Fill(pt_lab);
    if (fabs(eta1) < 0.5 && fabs(eta2) < 0.5) hNPi0_acc->Fill(pt_lab);

  }


  // hNPi0_gen_pt->Draw("");
  // hNPi0_acc->Draw("same");
  // cExtractSignal->SaveAs(Form("Simulation/InvarianteMasseSignal%s.pdf", AddName.Data()));

}
