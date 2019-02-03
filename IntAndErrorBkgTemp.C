#include "CommonHeader.h"

void IntAndErrorBkgTemp(){

  TFile* InputFile = SafelyOpenRootfile("OutputFileNormal.root");
  TH1D* hCorrBkg   = NULL;

  Double_t int_value = 0;
  Double_t error_value = 0;


  ofstream OutputFile;
  OutputFile.open ("IntAndError.txt");
  OutputFile << std::setw(10) << "Binnumber\t\t";
  OutputFile << std::setw(10) << "Int Value\t\t";
  OutputFile << std::setw(10) << "Error Value.\n";

  for (int k = 1; k < numberbins; ++k) {

    if(k >= 39){
      std::cout << "k ist zu gross!" << '\n';
      continue;
    }
    std::cout << "Start bin  " << k << " reading and wrinting!" << "\n";

    hCorrBkg = (TH1D*) InputFile->Get(Form("hCorrBack_bin%02d",k));

    Int_t lowerEffirange = hCorrBkg->FindBin(0.0);
    Int_t upperEffirange = hCorrBkg->FindBin(0.3);
    int_value            = hCorrBkg->IntegralAndError(lowerEffirange, upperEffirange, error_value);

    // if(fabs(error_value)>fabs(int_value)){
      OutputFile << std::setw(10) <<  k << "\t\t";
      OutputFile << std::setw(10) <<  int_value << "\t\t";
      OutputFile << std::setw(10) <<  error_value << "\n";
    // }

    std::cout << "End bin  " << k << " reading and wrinting!" << "\n\n";
  }
  OutputFile.close();
}
