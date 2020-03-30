//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul  7 12:11:54 2014 by ROOT version 5.34/14
// from TTree sd_data/Sensitive Detector Data
// found on file: /work/nisiyama/result/msm/det81.root
//////////////////////////////////////////////////////////

#ifndef sddata_h
#define sddata_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class sddata {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           charge;
   Int_t           code;
   Double_t        x;
   Double_t        y;
   Double_t        z;
   Double_t        vx;
   Double_t        vy;
   Double_t        vz;
   Double_t        kE;
   Double_t        mom;
   Int_t           pricharge;
   Int_t           pricode;
   Double_t        prix;
   Double_t        priy;
   Double_t        priz;
   Double_t        privx;
   Double_t        privy;
   Double_t        privz;
   Double_t        prikE;
   Double_t        scat;
   Int_t           eID;
   Int_t           bID;
   Double_t        cz;
   Double_t        phi;

   // List of branches
   TBranch        *b_charge;   //!
   TBranch        *b_code;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_kE;   //!
   TBranch        *b_mom;   //!
   TBranch        *b_pricharge;   //!
   TBranch        *b_pricode;   //!
   TBranch        *b_prix;   //!
   TBranch        *b_priy;   //!
   TBranch        *b_priz;   //!
   TBranch        *b_privx;   //!
   TBranch        *b_privy;   //!
   TBranch        *b_privz;   //!
   TBranch        *b_prikE;   //!
   TBranch        *b_scat;   //!
   TBranch        *b_eID;   //!
   TBranch        *b_bID;   //!
   TBranch        *b_cz;   //!
   TBranch        *b_phi;   //!

   sddata(TTree *tree=0);
   virtual ~sddata();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef sddata_cxx
sddata::sddata(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/work/nisiyama/result/msm/det81.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/work/nisiyama/result/msm/det81.root");
      }
      f->GetObject("sd_data",tree);

   }
   Init(tree);
}

sddata::~sddata()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t sddata::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t sddata::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void sddata::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("code", &code, &b_code);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("kE", &kE, &b_kE);
   fChain->SetBranchAddress("mom", &mom, &b_mom);
   fChain->SetBranchAddress("pricharge", &pricharge, &b_pricharge);
   fChain->SetBranchAddress("pricode", &pricode, &b_pricode);
   fChain->SetBranchAddress("prix", &prix, &b_prix);
   fChain->SetBranchAddress("priy", &priy, &b_priy);
   fChain->SetBranchAddress("priz", &priz, &b_priz);
   fChain->SetBranchAddress("privx", &privx, &b_privx);
   fChain->SetBranchAddress("privy", &privy, &b_privy);
   fChain->SetBranchAddress("privz", &privz, &b_privz);
   fChain->SetBranchAddress("prikE", &prikE, &b_prikE);
   fChain->SetBranchAddress("scat", &scat, &b_scat);
   fChain->SetBranchAddress("eID", &eID, &b_eID);
   fChain->SetBranchAddress("bID", &bID, &b_bID);
   fChain->SetBranchAddress("cz", &cz, &b_cz);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   Notify();
}

Bool_t sddata::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void sddata::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t sddata::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef sddata_cxx
