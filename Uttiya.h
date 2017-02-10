//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 10 14:18:47 2017 by ROOT version 5.34/26
// from TTree miniTree/myTree
// found on file: DY_Sample7m.root
//////////////////////////////////////////////////////////

#ifndef Uttiya_h
#define Uttiya_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class Uttiya {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<double>  *ElePx;
   vector<double>  *ElePy;
   vector<double>  *ElePz;
   vector<double>  *ElePt;
   vector<double>  *EleEta;
   vector<double>  *ElePhi;
   vector<double>  *EleTheta;
   vector<double>  *EleE;
   vector<double>  *EleCharge;
   vector<double>  *EleIdLoose;
   vector<double>  *EleIdTight;
   vector<double>  *EleIdRobustLoose;
   vector<double>  *EleIdRobustTight;
   vector<double>  *EleTrackIso;
   vector<double>  *EleEcalIso;
   vector<double>  *EleCaloIso;
   vector<double>  *EleNHits;
   vector<double>  *EleNChi2;
   vector<double>  *EleRecoDxy;
   vector<double>  *EleRecoDz;
   vector<double>  *JetPt;
   vector<double>  *JetEta;
   vector<double>  *JetPhi;
   vector<double>  *JetPx;
   vector<double>  *JetPy;
   vector<double>  *JetPz;
   vector<double>  *JetE;
   vector<double>  *jptJetPt;
   vector<double>  *jptJetEta;
   vector<double>  *jptJetPhi;
   vector<double>  *jptJetPx;
   vector<double>  *jptJetPy;
   vector<double>  *jptJetPz;
   vector<double>  *jptJetE;
   vector<double>  *pfJetPt;
   vector<double>  *pfJetEta;
   vector<double>  *pfJetPhi;
   vector<double>  *pfJetPx;
   vector<double>  *pfJetPy;
   vector<double>  *pfJetPz;
   vector<double>  *pfJetE;

   // List of branches
   TBranch        *b_ElePx;   //!
   TBranch        *b_ElePy;   //!
   TBranch        *b_ElePz;   //!
   TBranch        *b_ElePt;   //!
   TBranch        *b_EleEta;   //!
   TBranch        *b_ElePhi;   //!
   TBranch        *b_EleTheta;   //!
   TBranch        *b_EleE;   //!
   TBranch        *b_EleCharge;   //!
   TBranch        *b_EleIdLoose;   //!
   TBranch        *b_EleIdTight;   //!
   TBranch        *b_EleIdRobustLoose;   //!
   TBranch        *b_EleIdRobustTight;   //!
   TBranch        *b_EleTrackIso;   //!
   TBranch        *b_EleEcalIso;   //!
   TBranch        *b_EleCaloIso;   //!
   TBranch        *b_EleNHits;   //!
   TBranch        *b_EleNChi2;   //!
   TBranch        *b_EleRecoDxy;   //!
   TBranch        *b_EleRecoDz;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetPx;   //!
   TBranch        *b_JetPy;   //!
   TBranch        *b_JetPz;   //!
   TBranch        *b_JetE;   //!
   TBranch        *b_jptJetPt;   //!
   TBranch        *b_jptJetEta;   //!
   TBranch        *b_jptJetPhi;   //!
   TBranch        *b_jptJetPx;   //!
   TBranch        *b_jptJetPy;   //!
   TBranch        *b_jptJetPz;   //!
   TBranch        *b_jptJetE;   //!
   TBranch        *b_pfJetPt;   //!
   TBranch        *b_pfJetEta;   //!
   TBranch        *b_pfJetPhi;   //!
   TBranch        *b_pfJetPx;   //!
   TBranch        *b_pfJetPy;   //!
   TBranch        *b_pfJetPz;   //!
   TBranch        *b_pfJetE;   //!

   Uttiya(TTree *tree=0);
   virtual ~Uttiya();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Uttiya_cxx
Uttiya::Uttiya(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DY_Sample7m.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("DY_Sample7m.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("DY_Sample7m.root:/analyzeBasicPat");
      dir->GetObject("miniTree",tree);

   }
   Init(tree);
}

Uttiya::~Uttiya()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Uttiya::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Uttiya::LoadTree(Long64_t entry)
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

void Uttiya::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ElePx = 0;
   ElePy = 0;
   ElePz = 0;
   ElePt = 0;
   EleEta = 0;
   ElePhi = 0;
   EleTheta = 0;
   EleE = 0;
   EleCharge = 0;
   EleIdLoose = 0;
   EleIdTight = 0;
   EleIdRobustLoose = 0;
   EleIdRobustTight = 0;
   EleTrackIso = 0;
   EleEcalIso = 0;
   EleCaloIso = 0;
   EleNHits = 0;
   EleNChi2 = 0;
   EleRecoDxy = 0;
   EleRecoDz = 0;
   JetPt = 0;
   JetEta = 0;
   JetPhi = 0;
   JetPx = 0;
   JetPy = 0;
   JetPz = 0;
   JetE = 0;
   jptJetPt = 0;
   jptJetEta = 0;
   jptJetPhi = 0;
   jptJetPx = 0;
   jptJetPy = 0;
   jptJetPz = 0;
   jptJetE = 0;
   pfJetPt = 0;
   pfJetEta = 0;
   pfJetPhi = 0;
   pfJetPx = 0;
   pfJetPy = 0;
   pfJetPz = 0;
   pfJetE = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ElePx", &ElePx, &b_ElePx);
   fChain->SetBranchAddress("ElePy", &ElePy, &b_ElePy);
   fChain->SetBranchAddress("ElePz", &ElePz, &b_ElePz);
   fChain->SetBranchAddress("ElePt", &ElePt, &b_ElePt);
   fChain->SetBranchAddress("EleEta", &EleEta, &b_EleEta);
   fChain->SetBranchAddress("ElePhi", &ElePhi, &b_ElePhi);
   fChain->SetBranchAddress("EleTheta", &EleTheta, &b_EleTheta);
   fChain->SetBranchAddress("EleE", &EleE, &b_EleE);
   fChain->SetBranchAddress("EleCharge", &EleCharge, &b_EleCharge);
   fChain->SetBranchAddress("EleIdLoose", &EleIdLoose, &b_EleIdLoose);
   fChain->SetBranchAddress("EleIdTight", &EleIdTight, &b_EleIdTight);
   fChain->SetBranchAddress("EleIdRobustLoose", &EleIdRobustLoose, &b_EleIdRobustLoose);
   fChain->SetBranchAddress("EleIdRobustTight", &EleIdRobustTight, &b_EleIdRobustTight);
   fChain->SetBranchAddress("EleTrackIso", &EleTrackIso, &b_EleTrackIso);
   fChain->SetBranchAddress("EleEcalIso", &EleEcalIso, &b_EleEcalIso);
   fChain->SetBranchAddress("EleCaloIso", &EleCaloIso, &b_EleCaloIso);
   fChain->SetBranchAddress("EleNHits", &EleNHits, &b_EleNHits);
   fChain->SetBranchAddress("EleNChi2", &EleNChi2, &b_EleNChi2);
   fChain->SetBranchAddress("EleRecoDxy", &EleRecoDxy, &b_EleRecoDxy);
   fChain->SetBranchAddress("EleRecoDz", &EleRecoDz, &b_EleRecoDz);
   fChain->SetBranchAddress("JetPt", &JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", &JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetPx", &JetPx, &b_JetPx);
   fChain->SetBranchAddress("JetPy", &JetPy, &b_JetPy);
   fChain->SetBranchAddress("JetPz", &JetPz, &b_JetPz);
   fChain->SetBranchAddress("JetE", &JetE, &b_JetE);
   fChain->SetBranchAddress("jptJetPt", &jptJetPt, &b_jptJetPt);
   fChain->SetBranchAddress("jptJetEta", &jptJetEta, &b_jptJetEta);
   fChain->SetBranchAddress("jptJetPhi", &jptJetPhi, &b_jptJetPhi);
   fChain->SetBranchAddress("jptJetPx", &jptJetPx, &b_jptJetPx);
   fChain->SetBranchAddress("jptJetPy", &jptJetPy, &b_jptJetPy);
   fChain->SetBranchAddress("jptJetPz", &jptJetPz, &b_jptJetPz);
   fChain->SetBranchAddress("jptJetE", &jptJetE, &b_jptJetE);
   fChain->SetBranchAddress("pfJetPt", &pfJetPt, &b_pfJetPt);
   fChain->SetBranchAddress("pfJetEta", &pfJetEta, &b_pfJetEta);
   fChain->SetBranchAddress("pfJetPhi", &pfJetPhi, &b_pfJetPhi);
   fChain->SetBranchAddress("pfJetPx", &pfJetPx, &b_pfJetPx);
   fChain->SetBranchAddress("pfJetPy", &pfJetPy, &b_pfJetPy);
   fChain->SetBranchAddress("pfJetPz", &pfJetPz, &b_pfJetPz);
   fChain->SetBranchAddress("pfJetE", &pfJetE, &b_pfJetE);
   Notify();
}

Bool_t Uttiya::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Uttiya::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Uttiya::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Uttiya_cxx
