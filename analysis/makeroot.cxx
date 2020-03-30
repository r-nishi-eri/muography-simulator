
#include "TH2.h"

void makeroot (char *dirname, char *filename, int idstart, int idend, char *filenameout){

  //TNtuple *sd_data = new TNtuple("sd_data", "sd_data", "code:x:y:z:vx:vy:vz:tE:pricode:prix:priy:priz:privx:privy:privz:pritE");
  TTree* tree = new TTree("sd_data","Sensitive Detector Data");

  Int_t charge;
  Int_t code;
  Double_t x;
  Double_t y;
  Double_t z;
  Double_t vx;
  Double_t vy;
  Double_t vz;
  Double_t kE;
  Double_t mom;
  Int_t pricharge;
  Int_t pricode;
  Double_t prix;
  Double_t priy;
  Double_t priz;
  Double_t privx;
  Double_t privy;
  Double_t privz;
  Double_t prikE;
  Double_t scat;
  Int_t eID;
  Int_t bID;
  Double_t cz;
  Double_t phi;

  tree->Branch("charge",&charge,"charge/I");
  tree->Branch("code",&code,"code/I");
  tree->Branch("x",&x,"x/D");
  tree->Branch("y",&y,"y/D");
  tree->Branch("z",&z,"z/D");
  tree->Branch("vx",&vx,"vx/D");
  tree->Branch("vy",&vy,"vy/D");
  tree->Branch("vz",&vz,"vz/D");
  tree->Branch("kE",&kE,"kE/D");
  tree->Branch("mom",&mom,"mom/D");
  tree->Branch("pricharge",&pricharge,"pricharge/I");
  tree->Branch("pricode",&pricode,"pricode/I");
  tree->Branch("prix",&prix,"prix/D");
  tree->Branch("priy",&priy,"priy/D");
  tree->Branch("priz",&priz,"priz/D");
  tree->Branch("privx",&privx,"privx/D");
  tree->Branch("privy",&privy,"privy/D");
  tree->Branch("privz",&privz,"privz/D");
  tree->Branch("prikE",&prikE,"prikE/D");
  tree->Branch("scat",&scat,"scat/D");
  tree->Branch("eID",&eID,"eID/I");
  tree->Branch("bID",&bID,"bID/I");
  tree->Branch("cz",&cz,"cz/D");
  tree->Branch("phi",&phi,"phi/D");


  for(int id = idstart; id <= idend; id++){
    char fnamae [1024];
    sprintf(fnamae, "%s/%s-%d.dat", dirname, filename, id);
    ifstream fin(fnamae);
    while(fin >> charge >> code >> x >> y >> z >> vx >> vy >> vz >> kE >> 
	  pricharge >> pricode >> prix >> priy >> priz >> privx >> privy >> privz >> prikE >> 
	  eID >> bID){

      /*
      cout << charge << " "
	   << code << " " 
	   << x << " " 
	   << y << " " 
	   << z << " " 
	   << vx << " " 
	   << vy << " " 
	   << vz << " " 
	   << tE << " " 
	   << pricharge << " " 
	   << pricode << " " 
	   << prix << " " 
	   << priy << " " 
	   << priz << " " 
	   << privx << " " 
	   << privy << " " 
	   << privz << " " 
	   << pritE << " " 
	   << eID << " " 
	   << bID << endl;
      */

      // convert kE -> mom
      double mass = 0; // GeV
      if (code == 11 || code == -11) mass = 0.000510988910; 
      if (code == 13 || code == -13) mass = 0.1056583715;
      if (code == 15 || code == -15) mass = 1.77682;
      if (code == 2212 || code == -2212) mass = 0.938272013;
      if (code == 2112) mass = 0.939565346;
      if (code == 111) mass = 0.1349766;
      if (code == 211 || code == -211) mass = 0.13957018;
      if (code == 311) mass = 0.497648;
      if (code == 321 || code == -321) mass = 0.493667;
      mom = sqrt(kE * kE + 2.0 * kE * mass);

      cz = -vz;
      double phidet = atan2(y,x);
      phi = (cz > 0 ? atan2(vy, vx) : atan2(-vy,-vx)) - phidet;
      if (phi > +TMath::Pi()) phi -= TMath::Pi() * 2;
      if (phi < -TMath::Pi()) phi += TMath::Pi() * 2;
      scat = TMath::ACos(vx * privx + vy * privy + vz * privz);

      tree->Fill();
    }
  }

  TFile *f = new TFile(filenameout, "recreate");
  tree->Write();
  f->Close();  
  
}
