#include "CommonHeader.h"

TH1D* normpull(TH1D* a, TH1D* b, TString name){
  TH1D* ret = (TH1D*) a->Clone(name);
  Double_t uncertainty[a->fNcells-2];

  for (int j = 0; j < a->fNcells-2; j++) {
    uncertainty[j] = sqrt(pow(a->GetBinError(j),2.)+pow(b->GetBinError(j),2.));
  }
  for (int i = 0; i < a->fNcells-2; i++) {
    if(a->GetBinContent(i) != 0 && b->GetBinContent(i) != 0){
      ret->SetBinContent(i,fabs((a->GetBinContent(i)-b->GetBinContent(i))/uncertainty[i]));
    }
    else{
      ret->SetBinContent(i, 0.0);
    }
    ret->SetBinError(i, 0.0);
  }

  return ret;
}
