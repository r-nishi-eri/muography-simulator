#include "TH2.h"
#include "sddata.h"


// You have to load det**.root file before running this script.
void tubanalysis_flux(int nDivPhi, double phiMin, double phiMax, int nDivCosz, double coszMin, double coszMax){

  // ask some calculation configurations
  double detR;
  cout << "Detector Radius (m)?: " ;
  cin >> detR;
  double detH;
  cout << "Detector Height (m)?: ";
  cin >> detH;    
  double T;
  cout << "Elapsed Time (sec)?: ";
  cin >> T; cin.ignore();

  TH2F *h2 = new TH2F ("h2", "# particles", nDivPhi, phiMin, phiMax, nDivCosz, coszMin, coszMax);

  string str_contents = "cz:phi";
  
  string str_cut;
  cout << "cut?: ";
  getline(cin, str_cut);
  if (str_cut == "x") str_cut = "";
  
  string str_option;
  cout << "option?: ";
  getline(cin, str_option);
  if (str_cut == "x") str_cut = "";
  
  ofstream fout;
  char ans;
  cout << "want to output ascii file? (y/n): " << endl;
  cin >> ans;
  if (ans == 'y'){
    char fnameout[256];
    cout << "filename ?:" << endl;
    cin >> fnameout;
    fout.open(fnameout); 
  }
  
  TCanvas *canvas = new TCanvas ("canvas", "canvas", 1000, 1000);
  canvas->Divide(2,2);
  
  // scatter plot
  canvas->cd(1);
  char cut_scat[1024];
  sprintf(cut_scat, "%s&&%f<phi&&phi<%f&&%f<cz&&cz<%f", str_cut.c_str(), phiMin, phiMax, coszMin, coszMax);
  sd_data->Draw(str_contents.c_str(), cut_scat, "scat");
  
  // raw plot
  canvas->cd(2)->SetLogz();
  sd_data->Draw((str_contents+">>h2").c_str(), str_cut.c_str(), "COLZ,TEXT");
  
  // converted flux
  TH2F *h3 = (TH2F*) h2->Clone();
  h3->SetName("h3");  
  h3->SetTitle("flux (m-2 sec-1 sr-1)");
  
  int numPhi =   h2->GetXaxis()->GetNbins();
  int numTheta = h2->GetYaxis()->GetNbins();
  
  for (int iTheta = 1; iTheta <= numTheta; iTheta++){
    for (int iPhi = 1; iPhi <= numPhi; iPhi++){
      
      double nObs = h2->GetBinContent(iPhi,iTheta);
      double phimin = phiMin + (phiMax - phiMin) / double(numPhi) * double(iPhi - 1);
      double phimax = phiMin + (phiMax - phiMin) / double(numPhi) * double(iPhi);
      double cos_phimin = cos(phimin);
      double cos_phimax = cos(phimax);
      double costhetamin = coszMin + (coszMax - coszMin) / double(numTheta) * double(iTheta - 1);
      double costhetamax = coszMin + (coszMax - coszMin) / double(numTheta) * double(iTheta);
      double thetamin = acos(costhetamax);
      double thetamax = acos(costhetamin);
      
      
      double s0 = 2.0 * detR * TMath::Pi() * detH;
      double s_factor1 = (thetamax - thetamin) / 2.0 - (sin(2.0*thetamax) - sin(2.0*thetamin)) / 4.0;
      double s_factor2;
      if      (phimax <= -TMath::Pi() / 2.0                               ) s_factor2 = sin(phimin) - sin(phimax);
      else if (phimin <= -TMath::Pi() / 2.0 && -TMath::Pi() / 2.0 < phimax) s_factor2 = sin(phimin) + sin(phimax) + 2.0;
      else if (-TMath::Pi() / 2.0 < phimin && phimax <= +TMath::Pi() / 2.0) s_factor2 = sin(phimax) - sin(phimin);
      else if (phimin <= +TMath::Pi() / 2.0 && +TMath::Pi() / 2.0 < phimax) s_factor2 = -sin(phimin) - sin(phimax) + 2.0;
      else                                                                  s_factor2 = sin(phimin) - sin(phimax); 
      
      double flux = nObs / T / (s0 * s_factor1 * s_factor2);
      double fluxerr = sqrt(nObs) / T / (s0 * s_factor1 * s_factor2);
      h3->SetBinContent(iPhi,iTheta,flux);
      if (ans == 'y'){
	fout << phimin << " " <<  phimax << " " << costhetamin << " " << costhetamax << " " << nObs << " " <<  flux  << " " << fluxerr << endl;
      }
    }
  }
  canvas->cd(3);
  h3->Draw(str_option.c_str());
  canvas->cd(3)->SetLogz();
  canvas->Draw();
  
}
