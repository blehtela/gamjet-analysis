//////////////////////////////////////////////////////////
//							// 
// Same structure as GamHistosFill.h but simplified.	//
// Created a MakeClass Summer24RunsClass on 18.02.25	//
// 							//
//////////////////////////////////////////////////////////

#ifndef CalcGenWeight_h
#define CalcGenWeight_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
#include <cstdio>
#include <map>
#include <string>
using namespace std;

// Header file for the classes stored in the TTree if any.
class CalcGenWeight {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   string          pthtbin;
   string          version;
   string          _filename; // file name for debugging purposes
   static const bool debugFiles = true;

  
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   Long64_t        genEventCount;
   Double_t        genEventSumw;
   Double_t        genEventSumw2;
   Int_t           nLHEScaleSumw;
   Double_t        LHEScaleSumw[9];   //[nLHEScaleSumw]
   Int_t           nLHEPdfSumw;
   Double_t        LHEPdfSumw[103];   //[nLHEPdfSumw]
   Int_t           nPSSumw;
   Double_t        PSSumw[4];   //[nPSSumw]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_genEventCount;   //!
   TBranch        *b_genEventSumw;   //!
   TBranch        *b_genEventSumw2;   //!
   TBranch        *b_nLHEScaleSumw;   //!
   TBranch        *b_LHEScaleSumw;   //!
   TBranch        *b_nLHEPdfSumw;   //!
   TBranch        *b_LHEPdfSumw;   //!
   TBranch        *b_nPSSumw;   //!
   TBranch        *b_PSSumw;   //!


   CalcGenWeight(TTree *tree=0, string pthtbinname="X");
   virtual ~CalcGenWeight();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef CalcGenWeight_cxx
CalcGenWeight::CalcGenWeight(TTree *tree, string pthtbinname) : fChain(0), pthtbin(pthtbinname)
{

  // Use data set to decide on active branches
  string& bin = pthtbinname;
  //isSummer24 = (bin=="PTG10to100-HT40to100" || bin=="PTG10to100-HT40to100"); 
  
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.


// That's how it was coming out of MakeClass, but i adapted it to match with  GamHistosFill version
/*
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://hip-cms-se.csc.fi/store/user/rverma/cms-jerc-run3/Skim/GamJet/2024/MCSummer24/GJets/date-14Feb2025_time-111503_commit-40144c3/MCSummer24_2024_GamJet_GJetsHT40to100PTG10to100_Skim_5of75.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://hip-cms-se.csc.fi/store/user/rverma/cms-jerc-run3/Skim/GamJet/2024/MCSummer24/GJets/date-14Feb2025_time-111503_commit-40144c3/MCSummer24_2024_GamJet_GJetsHT40to100PTG10to100_Skim_5of75.root");
      }
      f->GetObject("Runs",tree);

   }
   Init(tree);
*/

   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("Runs",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain of trees.
      TChain * chain = new TChain("Runs","");
      chain->Add("root://hip-cms-se.csc.fi/store/user/rverma/cms-jerc-run3/Skim/GamJet/2024/MCSummer24/GJets/date-14Feb2025_time-111503_commit-40144c3/MCSummer24_2024_GamJet_GJetsHT40to100PTG10to100_Skim_5of75.root/Runs"); //not so sure about this? (analoguous to GamHistosFill, but not from MakeClass)
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
} //end of constructor


CalcGenWeight::~CalcGenWeight()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CalcGenWeight::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CalcGenWeight::LoadTree(Long64_t entry)
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

void CalcGenWeight::Init(TTree *tree)
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

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("genEventCount", &genEventCount, &b_genEventCount);
   fChain->SetBranchAddress("genEventSumw", &genEventSumw, &b_genEventSumw);
   fChain->SetBranchAddress("genEventSumw2", &genEventSumw2, &b_genEventSumw2);
   fChain->SetBranchAddress("nLHEScaleSumw", &nLHEScaleSumw, &b_nLHEScaleSumw);
   fChain->SetBranchAddress("LHEScaleSumw", LHEScaleSumw, &b_LHEScaleSumw);
   fChain->SetBranchAddress("nLHEPdfSumw", &nLHEPdfSumw, &b_nLHEPdfSumw);
   fChain->SetBranchAddress("LHEPdfSumw", LHEPdfSumw, &b_LHEPdfSumw);
   fChain->SetBranchAddress("nPSSumw", &nPSSumw, &b_nPSSumw);
   fChain->SetBranchAddress("PSSumw", PSSumw, &b_PSSumw);
   Notify();

   //in GamHistosFill.h made some safety reset of all branches, not here...
}

Bool_t CalcGenWeight::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

  //kept the followin from GamHistosFill code
  if (fChain && fChain->GetCurrentFile()) {
    _filename = fChain->GetCurrentFile()->GetName();
    if (debugFiles)
      cout << endl << "Opened file: " << _filename << endl << flush;
  }

   return kTRUE;
}

void CalcGenWeight::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CalcGenWeight::Cut(Long64_t entry)
{
  if (entry) {}; // suppress compiler error //not sure what this was about, kept it...
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef CalcGenWeight_cxx
