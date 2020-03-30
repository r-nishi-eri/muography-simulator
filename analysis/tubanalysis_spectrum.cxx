
#include "TH1.h"
#include "TH2.h"
#include "sddata.h"

#include "TMath.h"

using namespace TMath;




// Rebinning to log scale
void BinLogX(TH1 *h){
  TAxis *axis = h->GetXaxis();
  Int_t bins = axis->GetNbins();
  
  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t *new_bins = new Axis_t[bins+1];
  
  for(Int_t i = 0 ; i <= bins; i++){
    new_bins[i] = from * Exp((Log(to)-Log(from))*(Float_t)i/(Float_t)bins);
  }
  axis->Set(bins, new_bins);
  delete new_bins;
}





// You have to load det**.root file before running this script.
void tubanalysis_spectrum(char* filenameRoot, double phiMin, double phiMax, double czMin, double czMax, int nDiv, double kEmin, double kEmax, double detRadius, double detHeight, double detTime, char* cut, char* filenameOut){

  TFile *file = new TFile(filenameRoot, "READ");

  // get configurations from the arguments
  double detR = detRadius;
  double detH = detHeight;
  double T = detTime;

  TH1F *h1 = new TH1F ("h1", "energy spectrum", nDiv, kEmin, kEmax);
  BinLogX(h1);

  string str_cut = cut;
  
  ofstream fout;
  fout.open(filenameOut); 
  fout << "# phiMin: " << phiMin << endl
       << "# phiMax: " << phiMax << endl
       << "# czMin:  " << czMin << endl
       << "# czMax:  " << czMax << endl
       << "# cut: " << str_cut << endl
       << "# nDiv: " << nDiv << endl
       << "# kEmin: " << kEmin << endl
       << "# kEmax: " << kEmax << endl;

  

  TCanvas *canvas1 = new TCanvas ("canvas", "number hit spectrum", 600, 600);
  canvas1->cd();
  canvas1->SetLogx();
  char cutrange[1024];
  sprintf(cutrange, "%s && %f<=phi && phi<=%f && %f<=cz && cz<=%f", str_cut.c_str(), phiMin, phiMax, czMin, czMax);
  sd_data->Draw("kE>>h1", cutrange);


  double cos_phimin = cos(phiMin);
  double cos_phimax = cos(phiMax);
  double thetamin = acos(czMax);
  double thetamax = acos(czMin);

  double s0 = 2.0 * detR * Pi() * detH;
  double s_factor1 = (thetamax - thetamin) / 2.0 - (sin(2.0*thetamax) - sin(2.0*thetamin)) / 4.0;
  double s_factor2;
  if      (phiMax <= -Pi() / 2.0                        ) s_factor2 =  sin(phiMin) - sin(phiMax);
  else if (phiMin <= -Pi() / 2.0 && -Pi() / 2.0 < phiMax) s_factor2 =  sin(phiMin) + sin(phiMax) + 2.0;
  else if (-Pi() / 2.0 < phiMin && phiMax <= +Pi() / 2.0) s_factor2 =  sin(phiMax) - sin(phiMin);
  else if (phiMin <= +Pi() / 2.0 && +Pi() / 2.0 < phiMax) s_factor2 = -sin(phiMin) - sin(phiMax) + 2.0;
  else                                                    s_factor2 =  sin(phiMin) - sin(phiMax); 
    
  TH1F *hflux = h1->Clone("hflux");

  
  TAxis *axis = h1->GetXaxis();
  for (int iBin = 1; iBin <= axis->GetNbins() ; iBin++){
    
    double npars = h1->GetBinContent(iBin);
    double nparsbelow = (iBin >= 2 ? h1->Integral(1, iBin-1) : 0);
    double nparsabove = h1->Integral(iBin, axis->GetNbins());
    double de = axis->GetBinWidth(iBin);
    double ecenter = axis->GetBinCenter(iBin);
    double emin = axis->GetBinLowEdge(iBin);
    double emax = axis->GetBinUpEdge (iBin);
        
    double flux = npars / T / (s0 * s_factor1 * s_factor2);
    double fluxerr = sqrt(npars) / T / (s0 * s_factor1 * s_factor2);
    double fluxbelow = nparsbelow / T / (s0 * s_factor1 * s_factor2);
    double fluxbelowerr = sqrt(nparsbelow) / T / (s0 * s_factor1 * s_factor2);
    double fluxabove = nparsabove / T / (s0 * s_factor1 * s_factor2);
    double fluxaboveerr = sqrt(nparsabove) / T / (s0 * s_factor1 * s_factor2);

    hflux->SetBinContent(iBin, flux);
    hflux->SetBinError(iBin, fluxerr);

    fout << emin << " " 
	 << emax << " " 
	 << npars << " " 
	 << flux << " " 
	 << fluxerr << " " 
	 << flux / de << " " 
	 << fluxerr / de << " " 
	 << fluxbelow << " " 
	 << fluxbelowerr << " " 
	 << fluxabove << " " 
	 << fluxaboveerr << endl;
  }


  TCanvas *canvas2 = new TCanvas ("canvas2", "energy spectrum", 600, 600);
  canvas2->cd();
  canvas2->SetLogx();
  canvas2->SetLogy();
  hflux->Draw();

  fout.close();

  
}
