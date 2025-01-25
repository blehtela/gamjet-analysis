//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 28 16:43:19 2021 by ROOT version 6.24/00
// from TChain Events/
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
// This class has been modified to remove extra branches
// and make a hybrid of 2016, 2017 and 2018 data and MC
// in GamHistosFill2016.h and GamHistosFill2018.h
//////////////////////////////////////////////////////////

#ifndef GamHistosFill_h
#define GamHistosFill_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1D.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include <iostream>
#include <cstdio>
#include <map>
#include <string>
using namespace std;

//typedef for lumimap
typedef std::map<int,double> LumiMap;


// Header file for the classes stored in the TTree if any.

class GamHistosFill {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   int             isMC;     // data=0, PythiaPtFlat=1, MadGraphHT=2
   bool            is16apv, is16fgh, is16;
   bool            is17;
   bool            is18;
   bool            is22;
   bool            is22v10;
   bool            is23;
   bool            is24;
   bool            isRun2, isRun3;
   bool            isQCD;
   bool            isMG;
   string          dataset;
   string	   puera;  //data era used for pu reweighting
   string          version;
   string          _filename; // file name for debugging purposes
   static const bool debugFiles = true;
  
// Fixed size dimensions of array or collections stored in the TTree if any.
   // nCorrT1MetJet, nFatJet, nJet, nSoftActivityJet, nSubJet
   static const int nJetMax = 200;
   // nPhoton, nFsrPhoton, nGenIsolatedPhoton
   static const int nPhotonMax = 300; //200;//100;//20; //needed to go to higher number, as 24F had an event with more than 200 photons...
   // nElectron, nGenDressedLepton
   static const int nElectronMax = 10;
   static const int nTauMax = 10;
   static const int nMuonMax = 20;

   // nGenJet, nGenJetAK8, nSubGenJetAK8
   static const int nGenJetMax = 100;
   static const int nGenPartMax = 200;
   static const int nIsoTrackMax = 10;
   static const int nLHEPartMax = 10;

   static const int nTrigObjMax = 50;
   static const int nSVMax = 100;

   //static const int nPSWeightMax = 4; //THIS SHOULD BE USED PER DEFAULT - FOUR WEIGHT VARIATIONS
   static const int nPSWeightMax = 44; //THIS IS JUST A TEST for winter2024 madgraph samples


// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          bunchCrossing; // NanoV12
  
   Float_t         ChsMET_phi; // Run2
   Float_t         ChsMET_pt;  // Run2
   Float_t         RawPuppiMET_phi; // Run3
   Float_t         RawPuppiMET_pt;  // Run3


   //UInt_t          nJet; // NanoV11,10
   Int_t           nJet; // NanoV12
   Float_t         Jet_pt[nJetMax];
   Float_t         Jet_eta[nJetMax];
   Float_t         Jet_phi[nJetMax];
   Float_t         Jet_mass[nJetMax];

   Float_t         Jet_rawFactor[nJetMax];
   Float_t         Jet_area[nJetMax];
   //Int_t           Jet_jetId[nJetMax]; // NanoV10,11
   UChar_t         Jet_jetId[nJetMax];   // NanoV12

   Float_t         Jet_btagDeepB[nJetMax]; // Run2
   Float_t         Jet_btagDeepC[nJetMax]; // Run2
   Float_t         Jet_qgl[nJetMax];       // Run2

   Float_t         Jet_btagDeepFlavB[nJetMax];   // Run3
   Float_t         Jet_btagDeepFlavCvB[nJetMax]; // Run3
   Float_t         Jet_btagDeepFlavCvL[nJetMax]; // Run3
   Float_t         Jet_btagDeepFlavQG[nJetMax];  // Run3

   Float_t         Jet_chHEF[nJetMax];
   Float_t         Jet_neHEF[nJetMax];
   Float_t         Jet_neEmEF[nJetMax];
   Float_t         Jet_chEmEF[nJetMax];
   Float_t         Jet_muEF[nJetMax];
   //Float_t         Jet_chFPV0EF[nJetMax];   //[nJet]

   Float_t         Jet_hfEmEF[nJetMax];   //[nJet]
   Float_t         Jet_hfHEF[nJetMax];   //[nJet]

  
   //Int_t           Jet_nConstituents[nJetMax];   //[nJet]
   //Int_t           Jet_nElectrons[nJetMax];   //[nJet]
   //Int_t           Jet_nMuons[nJetMax];   //[nJet]

   //Int_t           Jet_genJetIdx[nJetMax]; // NanoV11,10
   Short_t           Jet_genJetIdx[nJetMax]; // NanoV12
  
   //UInt_t          nPhoton; // NanoV11,10
   Int_t           nPhoton; // NanoV12
   Float_t         Photon_pt[nPhotonMax];
   Float_t         Photon_eta[nPhotonMax];
   Float_t         Photon_phi[nPhotonMax];
   Float_t         Photon_mass[nPhotonMax]; // Run2

   Float_t         Photon_hoe[nPhotonMax];
   //Int_t           Photon_cutBased[nPhotonMax]; // NanoV11,10
   UChar_t         Photon_cutBased[nPhotonMax];   // NanoV12
   Int_t           Photon_cutBasedBitmap[nPhotonMax];   // EXTRA
   //Int_t           Photon_jetIdx[nPhotonMax];   // NanoV10,V11
   Short_t         Photon_jetIdx[nPhotonMax];     // NanoV12

   UChar_t         Photon_seedGain[nPhotonMax];
   Float_t         Photon_r9[nPhotonMax];
   Float_t         Photon_eCorr[nPhotonMax];
   Float_t         Photon_energyErr[nPhotonMax];
  
   Float_t         fixedGridRhoFastjetAll; // Run2,Run3
   Float_t         Rho_fixedGridRhoFastjetAll; // Run2,Run3
   Float_t         Rho_fixedGridRhoFastjetCentral; // Run2,Run3
   Float_t         Rho_fixedGridRhoFastjetCentralChargedPileUp; // Run2,Run3

   //UInt_t          nTrigObj; // NanoV11,10
   //Int_t          nTrigObj; // NanoV12
   //Float_t         TrigObj_pt[nTrigObjMax];   //[nTrigObj]
   //Float_t         TrigObj_eta[nTrigObjMax];   //[nTrigObj]
   //Float_t         TrigObj_phi[nTrigObjMax];   //[nTrigObj]
   //Float_t         TrigObj_l1pt[nTrigObjMax];   //[nTrigObj]
   //Float_t         TrigObj_l1pt_2[nTrigObjMax];   //[nTrigObj]
   //Float_t         TrigObj_l2pt[nTrigObjMax];   //[nTrigObj]
   //Int_t           TrigObj_id[nTrigObjMax];   //[nTrigObj]
   //Int_t           TrigObj_l1iso[nTrigObjMax];   //[nTrigObj]
   //Int_t           TrigObj_l1charge[nTrigObjMax];   //[nTrigObj]
   //Int_t           TrigObj_filterBits[nTrigObjMax];   //[nTrigObj]
   //UInt_t          nOtherPV;
   //Float_t         OtherPV_z[nSVMax];   //[nOtherPV]
   //Float_t         PV_ndof;
   //Float_t         PV_x;
   //Float_t         PV_y;
   //Float_t         PV_z;
   //Float_t         PV_chi2;
   //Float_t         PV_score;

   //Int_t           PV_npvs; // NanoV11,10
   UChar_t           PV_npvs; // NanoV12
   //Int_t           PV_npvsGood; // NanoV11,10
   UChar_t           PV_npvsGood; // NanoV12
 
   Bool_t          Flag_HBHENoiseFilter;
   //Bool_t          Flag_HBHENoiseIsoFilter;	//this and some others commented out since w44 (see git 25.01.2025)
   //Bool_t          Flag_CSCTightHaloFilter;
   //Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter;
   //Bool_t          Flag_CSCTightHalo2015Filter;
   //Bool_t          Flag_globalTightHalo2016Filter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   //Bool_t          Flag_HcalStripHaloFilter;
   //Bool_t          Flag_hcalLaserEventFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   //Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_eeBadScFilter;
   //Bool_t          Flag_ecalLaserCorrFilter;
   //Bool_t          Flag_trkPOGFilters;
   //Bool_t          Flag_chargedHadronTrackResolutionFilter;
   //Bool_t          Flag_muonBadTrackFilter;
   //Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadPFMuonDzFilter;
   Bool_t          Flag_hfNoisyHitsFilter;
   //Bool_t          Flag_BadChargedCandidateSummer16Filter;
   //Bool_t          Flag_BadPFMuonSummer16Filter;
   //Bool_t          Flag_trkPOG_manystripclus53X;
   //Bool_t          Flag_trkPOG_toomanystripclus53X;
   //Bool_t          Flag_trkPOG_logErrorTooManyClusters;
   //Bool_t          Flag_METFilters;

   // Triggers from 2016
   Bool_t          HLT_Photon250_NoHE;
   Bool_t          HLT_Photon300_NoHE;

   // Triggers from 2016
   Bool_t          HLT_Photon22;
   Bool_t          HLT_Photon30;
   Bool_t          HLT_Photon36;
   Bool_t          HLT_Photon50;
   Bool_t          HLT_Photon75;
   Bool_t          HLT_Photon90;
   Bool_t          HLT_Photon120;
   Bool_t          HLT_Photon175;
   Bool_t          HLT_Photon165_HE10;
   Bool_t          HLT_Photon22_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon30_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon36_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon50_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon90_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon120_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon165_R9Id90_HE10_IsoM;

   //Bool_t          HLT_HISinglePhoton10;
   //Bool_t          HLT_HISinglePhoton15;
   //Bool_t          HLT_HISinglePhoton20;
   //Bool_t          HLT_HISinglePhoton40;
   //Bool_t          HLT_HISinglePhoton60;
   //Bool_t          HLT_EcalCalibration;
   //Bool_t          HLT_HcalCalibration;

   //Bool_t          HLT_Photon500;
   //Bool_t          HLT_Photon600;

   // Triggers and branches from 24
   //Bool_t          HLT_Photon50EB_TightID_TightIso_PFJet30;
   Bool_t          HLT_Photon50EB_TightID_TightIso;
   //Bool_t          HLT_Photon55EB_TightID_TightIso; //commented out since w44
   Bool_t          HLT_Photon75EB_TightID_TightIso;
   Bool_t          HLT_Photon90EB_TightID_TightIso;
   //Bool_t          HLT_Photon110EB_TightID_TightIso;
   //Bool_t          HLT_Photon110EB_TightID_TightIso_PFJet30;

   // Triggers and branches from 22-23
   Bool_t          HLT_Photon100EBHE10;
   Bool_t          HLT_Photon30EB_TightID_TightIso;
  
   // Triggers and branches from 2018
   Bool_t          HLT_Photon100EB_TightID_TightIso;
   Bool_t          HLT_Photon110EB_TightID_TightIso;
   Bool_t          HLT_Photon120EB_TightID_TightIso;
   Bool_t          HLT_Photon200;
   Bool_t          HLT_Photon20_HoverELoose;
   Bool_t          HLT_Photon30_HoverELoose;
   Bool_t          HLT_Photon150;
   Bool_t          HLT_Photon33;
   Bool_t          HLT_Photon20;

   // Triggers from 2017
   Bool_t          HLT_Photon40_HoverELoose;
   Bool_t          HLT_Photon50_HoverELoose;
   Bool_t          HLT_Photon60_HoverELoose;

   // Extras for MC
   //UInt_t          nGenJet; // NanoV11,10
   Int_t          nGenJet; // NanoV12
   Float_t         GenJet_eta[nGenJetMax];   //[nGenJet]
   Float_t         GenJet_mass[nGenJetMax];   //[nGenJet]
   Float_t         GenJet_phi[nGenJetMax];   //[nGenJet]
   Float_t         GenJet_pt[nGenJetMax];   //[nGenJet]
   //Int_t           GenJet_partonFlavour[nGenJetMax]; //[nGenJet] // NanoV10,11
   Short_t         GenJet_partonFlavour[nGenJetMax]; //[nGenJet] // NanoV12

   UInt_t          nGenIsolatedPhoton;
   Float_t         GenIsolatedPhoton_pt[nPhotonMax];
   Float_t         GenIsolatedPhoton_eta[nPhotonMax];
   Float_t         GenIsolatedPhoton_phi[nPhotonMax];
   Float_t         GenIsolatedPhoton_mass[nPhotonMax];

   Float_t         LHE_HT;
  
   Float_t         genWeight;
   Float_t         Pileup_nTrueInt;
  //UInt_t          nPSWeight; // NanoV11,10
   Int_t           nPSWeight; // NanoV12
   Float_t         PSWeight[nPSWeightMax];   //[nPSWeight]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!

   TBranch        *b_ChsMET_phi;   //!
   TBranch        *b_ChsMET_pt;   //!
  //TBranch        *b_ChsMET_sumEt;   //!

   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_btagDeepB;   //!
   TBranch        *b_Jet_btagDeepC;   //!
   //TBranch        *b_Jet_btagDeepCvB;   //!
   //TBranch        *b_Jet_btagDeepCvL;   //!
   TBranch        *b_Jet_btagDeepFlavB;   //!
  //TBranch        *b_Jet_btagDeepFlavC;   //!
   TBranch        *b_Jet_btagDeepFlavCvB;   //!
   TBranch        *b_Jet_btagDeepFlavCvL;   //!
   TBranch        *b_Jet_btagDeepFlavQG;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_chFPV0EF;   //!
   //TBranch        *b_Jet_chFPV1EF;   //!
   //TBranch        *b_Jet_chFPV2EF;   //!
   //TBranch        *b_Jet_chFPV3EF;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_muEF;   //!
   TBranch        *b_Jet_muonSubtrFactor;   //!
   TBranch        *b_Jet_neEmEF;   //!
   TBranch        *b_Jet_neHEF;   //!
   TBranch        *b_Jet_hfEmEF;   //! new
   TBranch        *b_Jet_hfHEF;   //! new
 
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   //TBranch        *b_Jet_puIdDisc;   //!
   TBranch        *b_Jet_qgl;   //!
   TBranch        *b_Jet_rawFactor;   //!
   TBranch        *b_Jet_jetId;   //!
   TBranch        *b_Jet_nConstituents;   //!
   TBranch        *b_Jet_nElectrons;   //!
   TBranch        *b_Jet_nMuons;   //!
  //TBranch        *b_Jet_puId;   //!

   TBranch        *b_Jet_genJetIdx;   //!

   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_pt;   //!
   TBranch        *b_MET_significance;   //!
   TBranch        *b_MET_sumEt;   //!
   //TBranch        *b_MET_sumPtUnclustered;   //!

   TBranch        *b_nPhoton;   //!
   TBranch        *b_Photon_energyErr;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_hoe;   //!
   TBranch        *b_Photon_mass;   //!
   TBranch        *b_Photon_mvaID;   //!
   TBranch        *b_Photon_mvaID_Fall17V1p1;   //!
  //TBranch        *b_Photon_pfRelIso03_all;   //!
  //TBranch        *b_Photon_pfRelIso03_chg;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_r9;   //!
   TBranch        *b_Photon_sieie;   //!
  //TBranch        *b_Photon_charge;   //!
   TBranch        *b_Photon_cutBased;   //!
   TBranch        *b_Photon_cutBasedBitmap;   // EXTRA
  //TBranch        *b_Photon_cutBased_Fall17V1Bitmap;   //!
   TBranch        *b_Photon_electronIdx;   //!
   TBranch        *b_Photon_jetIdx;   //!
  //TBranch        *b_Photon_pdgId;   //!
   TBranch        *b_Photon_vidNestedWPBitmap;   //!
   TBranch        *b_Photon_electronVeto;   //!
   TBranch        *b_Photon_isScEtaEB;   //!
   TBranch        *b_Photon_isScEtaEE;   //!
   TBranch        *b_Photon_mvaID_WP80;   //!
   TBranch        *b_Photon_mvaID_WP90;   //!
   TBranch        *b_Photon_pixelSeed;   //!
   TBranch        *b_Photon_seedGain;   //!

   TBranch        *b_RawMET_phi;   //!
   TBranch        *b_RawMET_pt;   //!
   TBranch        *b_RawMET_sumEt;   //!

   TBranch        *b_RawPuppiMET_phi;   //!
   TBranch        *b_RawPuppiMET_pt;   //!
   TBranch        *b_RawPuppiMET_sumEt;   //!

   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_fixedGridRhoFastjetCentral;   //!
  //TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralChargedPileUp;   //!
  //TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!

   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_l1pt;   //!
   TBranch        *b_TrigObj_l1pt_2;   //!
   TBranch        *b_TrigObj_l2pt;   //!
   TBranch        *b_TrigObj_id;   //!
   TBranch        *b_TrigObj_l1iso;   //!
   TBranch        *b_TrigObj_l1charge;   //!
   TBranch        *b_TrigObj_filterBits;   //!
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_OtherPV_z;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!

   TBranch        *b_Flag_HBHENoiseFilter;   //!
   //TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   //TBranch        *b_Flag_CSCTightHaloFilter;   //!
   //TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter;   //!
   //TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   //TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   //TBranch        *b_Flag_HcalStripHaloFilter;   //!
   //TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   //TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   //TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   //TBranch        *b_Flag_trkPOGFilters;   //!
   //TBranch        *b_Flag_chargedHadronTrackResolutionFilter;   //!
   //TBranch        *b_Flag_muonBadTrackFilter;   //!
   //TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_BadPFMuonDzFilter;   //!
   TBranch        *b_Flag_hfNoisyHitsFilter;   //!
   //TBranch        *b_Flag_BadChargedCandidateSummer16Filter;   //!
   //TBranch        *b_Flag_BadPFMuonSummer16Filter;   //!
   //TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   //TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   //TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   //TBranch        *b_Flag_METFilters;   //!

   // Extras for MC
   TBranch        *b_genWeight;
   TBranch        *b_Pileup_nTrueInt;
   TBranch        *b_nPSWeight;   //!
   TBranch        *b_PSWeight;   //!
   TBranch        *b_nGenJet;   //!
   TBranch        *b_GenJet_eta;   //!
   TBranch        *b_GenJet_mass;   //!
   TBranch        *b_GenJet_phi;   //!
   TBranch        *b_GenJet_pt;   //!
   TBranch        *b_GenJet_partonFlavour;   //!
   TBranch        *b_nGenIsolatedPhoton;   //!
   TBranch        *b_GenIsolatedPhoton_eta;   //!
   TBranch        *b_GenIsolatedPhoton_mass;   //!
   TBranch        *b_GenIsolatedPhoton_phi;   //!
   TBranch        *b_GenIsolatedPhoton_pt;   //!

   TBranch        *b_LHE_HT;   //!

   // Triggers from 2016
   TBranch        *b_HLT_Photon250_NoHE;   //!
   TBranch        *b_HLT_Photon300_NoHE;   //!

   TBranch        *b_HLT_Photon22;   //!
   TBranch        *b_HLT_Photon30;   //!
   TBranch        *b_HLT_Photon36;   //!
   TBranch        *b_HLT_Photon50;   //!
   TBranch        *b_HLT_Photon75;   //!
   TBranch        *b_HLT_Photon90;   //!
   TBranch        *b_HLT_Photon120;   //!
   TBranch        *b_HLT_Photon175;   //!
   TBranch        *b_HLT_Photon165_HE10;   //!
   TBranch        *b_HLT_Photon22_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon30_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon36_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon50_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon75_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon90_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon120_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon165_R9Id90_HE10_IsoM;   //!

   //TBranch        *b_HLT_HISinglePhoton10;   //!
   //TBranch        *b_HLT_HISinglePhoton15;   //!
   //TBranch        *b_HLT_HISinglePhoton20;   //!
   //TBranch        *b_HLT_HISinglePhoton40;   //!
   //TBranch        *b_HLT_HISinglePhoton60;   //!
   //TBranch        *b_HLT_EcalCalibration;   //!
   //TBranch        *b_HLT_HcalCalibration;   //!

   //TBranch        *b_HLT_Photon500;   //!
   //TBranch        *b_HLT_Photon600;   //!


   // Triggers (and branches) from 24
   TBranch        *b_HLT_Photon50EB_TightID_TightIso_PFJet30;   //!
   TBranch        *b_HLT_Photon50EB_TightID_TightIso;   //!
   //TBranch        *b_HLT_Photon55EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon75EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon90EB_TightID_TightIso;   //!
   //TBranch        *b_HLT_Photon110EB_TightID_TightIso;   //!
   //TBranch        *b_HLT_Photon110EB_TightID_TightIso_PFJet30;   //!


   // Triggers (and branches) from 22-23
   TBranch        *b_HLT_Photon100EBHE10;
   TBranch        *b_HLT_Photon30EB_TightID_TightIso;
  
   // Triggers (and branches) from 2018
   TBranch        *b_Photon_eCorr;
   TBranch        *b_HLT_Photon100EB_TightID_TightIso;
   TBranch        *b_HLT_Photon110EB_TightID_TightIso;
   TBranch        *b_HLT_Photon120EB_TightID_TightIso;
   TBranch        *b_HLT_Photon200;
   TBranch        *b_HLT_Photon150;
   TBranch        *b_HLT_Photon33;
   TBranch        *b_HLT_Photon20;
   TBranch        *b_HLT_Photon20_HoverELoose;
   TBranch        *b_HLT_Photon30_HoverELoose;
   
   // Triggers from 2017
   TBranch        *b_HLT_Photon40_HoverELoose;
   TBranch        *b_HLT_Photon50_HoverELoose;
   TBranch        *b_HLT_Photon60_HoverELoose;

  GamHistosFill(TTree *tree=0, int itype=1, string datasetname="X", string pueraname="",
		string versionname="vX");
   virtual ~GamHistosFill();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   // Code originally from jetphys/HistosFill.C
   void PrintInfo(string info, bool printcout);
   bool LoadJSON(string json);
   void LoadPU();
   LumiMap LoadLumi(string filename); //added function for writing lumimap based on brilcalc input file
   std::map<int, std::map<int, int> > _json;
   map<string, map<int, TH1D*> > _pu;
   map<string, map<int, double> >  _lumi;
};

#endif

#ifdef GamHistosFill_cxx
GamHistosFill::GamHistosFill(TTree *tree, int itype, string datasetname, string pueraname, string versionname)
  : fChain(0), isMC(itype), dataset(datasetname), puera(pueraname), version(versionname)
{

  // Use data set to decide on active branches
  string& ds = datasetname;
  is16apv = (ds=="2016B" || ds=="2016C" || ds=="2016D" || ds=="2016BCD" ||
	     ds=="2016E" || ds=="2016F" || ds=="2016EF" || ds=="2016BCDEF" ||
	     ds=="2016APVP8" || ds=="2016APVQCD");
  is16fgh = (ds=="2016FG" || ds=="2016H" || ds=="2016FGH" ||
	     ds=="2016P8" || ds=="2016QCD");
  is16 = (is16apv || is16fgh);
  is17 = (ds=="2017B" || ds=="2017C" || ds=="2017D" || ds=="2017E" ||
	  ds=="2017F" || ds=="2017BCDEF" || ds=="2017P8" || ds=="2017QCD");
  is18 = (ds=="2018A" || ds=="2018B" || ds=="2018C" || ds=="2018D" || 
	  ds=="2018A1" || ds=="2018A2" ||
	  ds=="2018D1" || ds=="2018D2" || ds=="2018D3" || ds=="2018D4" ||
	  ds=="2018ABCD" || ds=="2018P8" || ds=="2018QCD");
  is22 = (ds=="2022C" || ds=="2022D" || ds=="2022E" || ds=="2022F" || ds=="2022G" || 
	  ds=="2022Cnib1" || ds=="2022Dnib1" || ds=="2022Enib1" || ds=="2022Fnib1" || ds=="2022Gnib1" || //nibs and fibs
	  ds=="2022P8" || ds=="2022QCD" || ds=="2022EEP8" ||
	  ds=="2022EEQCD");
  //is22v10 = (ds=="2022F" || ds=="2022G" ||
  is22v10 = false;//(ds=="2022G");
  is23 = (ds=="2023B" || ds=="2023Cv123" || ds=="2023Cv4" || ds=="2023D" ||
	  ds=="2023Bnib1" || ds=="2023Cv123nib1" || ds=="2023Cv4nib1" || ds=="2023Cv4nib2" || ds=="2023Dnib1" || //nibs and fibs
          ds=="2023Cv123X" || ds=="2023Cv4X" || ds=="2023DX" || //for testing wX23 and wX22
	  ds=="2023P8X" || ds=="2023QCDX" || ds=="2023P8-BPixX"|| ds=="2023QCDX" || ds=="2023QCD-BPixX" || //for testing wX23 and wX22
	  ds=="2023P8" || ds=="2023QCD" || ds=="2023P8-BPix"|| ds=="2023QCD" || ds=="2023QCD-BPix"); //added 2023P8_BPix
  is24 = (ds=="2024B-PromptReco-v1" || ds=="2024B" || ds=="2024C" || ds=="2024D" || ds=="2024Ev1" || ds=="2024Ev2" || 
					ds=="2024F" || ds=="2024Fa" || ds=="2024Fb" || ds=="2024Fc" || ds=="2024Fd" || ds=="2024F-ECALCC-HCALDI-skim" ||
					ds=="2024G" || ds=="2024Ga" || ds=="2024Gb" || ds=="2024Gc" || ds=="2024Gd" || ds=="2024Ge" || ds=="2024Gtest" ||
					ds=="2024H" || ds=="2024Hskim" || ds=="2024Ha" || ds=="2024Hb" || ds=="2024Hc" || ds=="2024Hd" || ds=="2024He" || ds=="2024Htest" ||
					ds=="2024I" || ds=="2024Iskim" || ds=="2024Iv1" || ds=="2024Iv2" || ds=="2024Itest" ||
  		ds=="2024Bnib1" || ds=="2024Cnib1" || ds=="2024Dnib1" || ds=="2024Ev1nib1" || ds=="2024Ev2nib1" || 
					ds=="2024Fnib1" || ds=="2024Fnib2" || ds=="2024Fnib3" || ds=="2024F-ECALCC-HCALDI-nib1" || ds=="2024F-ECALCC-HCALDI-nib2" || ds=="2024F-ECALCC-HCALDI-nib3" ||
					ds=="2024Gnib1" || ds=="2024Gnib2" || ds=="2024Hnib1" || ds=="2024Inib1" ||
					ds=="2024B-ECALRATIO" || ds=="2024C-ECALRATIO" || ds=="2024C-ECALR-HCALDI" || ds=="2024C-ECALCC-HCALDI" ||
					ds=="winter2024P8" ||  //in use since 07.08.2024
					ds=="winter2024P8a" ||
					ds=="winter2024P8b" ||
					ds=="winter2024P8c" ||
					ds=="winter2024P8-test" ||
					ds=="winter2024P8-v14" ||  //in use since 21.08.2024
          ds=="2024P8" || ds=="2024QCD" || 
					ds=="2024QCDa" || ds=="2024QCDb" || ds=="2024QCDc" || ds=="2024QCDd" || ds=="2024QCDe" || ds=="2024QCDf" || 
					ds=="2024QCD-v14"); //added these already, even though no summer MC for 2024 yet, use winter2024P8 right now.
  isQCD = (ds=="2016QCD" || ds=="2016QCDAPV" || ds=="2017QCD" ||
	   ds=="2018QCD" || ds=="2022QCD" || ds=="2022EEQCD" ||
           ds=="2023QCD" || ds=="2023QCD-BPix" || 
           ds=="2024QCD" || ds=="2024QCD-v14" ||
					 ds=="2024QCDa" || ds=="2024QCDb" || ds=="2024QCDc" || ds=="2024QCDd" || ds=="2024QCDe" || ds=="2024QCDf"); //added 2024QCD already here
  isMG = (ds=="winter2024P8" || ds=="winter2024P8a" || ds=="winter2024P8b" || ds=="winter2024P8c" || 
					 ds=="2022P8" || ds=="2022EEP8" || 
					 ds=="2022QCD" || ds=="2022EEQCD" || //added winter24 here, too.. checking
           ds=="2023P8" || ds=="2023QCD" || ds=="2023P8-BPix" || ds=="2023QCD-BPix" || //); //should 2023P8 and 2023P8_BPix be added here, too? (for correct weight in HT bins)
					 ds=="winter2024P8" || ds=="winter2024P8-test" || ds=="winter2024P8-v14" || ds=="2024QCD" || 
					 ds=="2024QCD-v14" || ds=="2024QCDa" || ds=="2024QCDb" || ds=="2024QCDc" || ds=="2024QCDd" || ds=="2024QCDe" || ds=="2024QCDf" ||
           ds=="2024P8");  //added already here
  isRun3 = (is22 || is23 || is24);
  isRun2 = (is16  || is17 || is18);
  assert(is16 || is17 || is18 || is22 || is23 || is24);
  
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("Events",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("Events","");
      chain->Add("../nanoV8EGammaLocal/dasgo_Plus_xrdcp/SinglePhoton_Run2016H-UL2016_MiniAODv1_NanoAODv2-v1_NANOAOD/FB28C338-9F05-274B-812E-B1F77B49C09D.root/Events");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

GamHistosFill::~GamHistosFill()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GamHistosFill::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GamHistosFill::LoadTree(Long64_t entry)
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

void GamHistosFill::Init(TTree *tree)
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
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);

   if (isRun2) {
     fChain->SetBranchAddress("ChsMET_phi", &ChsMET_phi, &b_ChsMET_phi);
     fChain->SetBranchAddress("ChsMET_pt", &ChsMET_pt, &b_ChsMET_pt);
     //fChain->SetBranchAddress("ChsMET_sumEt", &ChsMET_sumEt, &b_ChsMET_sumEt);
   }

   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   if (!(is22 || is23 || is24)) fChain->SetBranchAddress("Jet_btagDeepB", Jet_btagDeepB, &b_Jet_btagDeepB);
   if (!(is22 || is23 || is24)) fChain->SetBranchAddress("Jet_btagDeepC", Jet_btagDeepC, &b_Jet_btagDeepC);
   //fChain->SetBranchAddress("Jet_btagDeepCvB", Jet_btagDeepCvB, &b_Jet_btagDeepCvB);
   //fChain->SetBranchAddress("Jet_btagDeepCvL", Jet_btagDeepCvL, &b_Jet_btagDeepCvL);
   if (is22 || is23 || is24) {
     fChain->SetBranchAddress("Jet_btagDeepFlavB", Jet_btagDeepFlavB, &b_Jet_btagDeepFlavB);
     //fChain->SetBranchAddress("Jet_btagDeepFlavC", Jet_btagDeepFlavC, &b_Jet_btagDeepFlavC);
     fChain->SetBranchAddress("Jet_btagDeepFlavCvB", Jet_btagDeepFlavCvB, &b_Jet_btagDeepFlavCvB);
     fChain->SetBranchAddress("Jet_btagDeepFlavCvL", Jet_btagDeepFlavCvL, &b_Jet_btagDeepFlavCvL);
     fChain->SetBranchAddress("Jet_btagDeepFlavQG", Jet_btagDeepFlavQG, &b_Jet_btagDeepFlavQG);
   }
   fChain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   //if (!(is22 || is23 || is24)) fChain->SetBranchAddress("Jet_chFPV0EF", Jet_chFPV0EF, &b_Jet_chFPV0EF);
   //fChain->SetBranchAddress("Jet_chFPV1EF", Jet_chFPV1EF, &b_Jet_chFPV1EF);
   //fChain->SetBranchAddress("Jet_chFPV2EF", Jet_chFPV2EF, &b_Jet_chFPV2EF);
   //fChain->SetBranchAddress("Jet_chFPV3EF", Jet_chFPV3EF, &b_Jet_chFPV3EF);
   fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   //fChain->SetBranchAddress("Jet_muonSubtrFactor", Jet_muonSubtrFactor, &b_Jet_muonSubtrFactor);
   fChain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);

   if (is24) fChain->SetBranchAddress("Jet_hfEmEF", Jet_hfEmEF, &b_Jet_hfEmEF);   //new in 2024
   if (is24) fChain->SetBranchAddress("Jet_hfHEF", Jet_hfHEF, &b_Jet_hfHEF);      //new in 2024

   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   //fChain->SetBranchAddress("Jet_puIdDisc", Jet_puIdDisc, &b_Jet_puIdDisc);
   if (!(is22 || is23 || is24)) fChain->SetBranchAddress("Jet_qgl", Jet_qgl, &b_Jet_qgl);
   fChain->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);
   fChain->SetBranchAddress("Jet_jetId", Jet_jetId, &b_Jet_jetId);
   //fChain->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   //fChain->SetBranchAddress("Jet_nElectrons", Jet_nElectrons, &b_Jet_nElectrons);
   //fChain->SetBranchAddress("Jet_nMuons", Jet_nMuons, &b_Jet_nMuons);
   //fChain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);

   if (isMC)
     fChain->SetBranchAddress("Jet_genJetIdx", Jet_genJetIdx, &b_Jet_genJetIdx);

   //fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   //fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   //fChain->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
   //fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   //fChain->SetBranchAddress("MET_sumPtUnclustered", &MET_sumPtUnclustered, &b_MET_sumPtUnclustered);

   fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   //if (!is16)
   if (!(is16 || is22 || is23 || is24)) 
     fChain->SetBranchAddress("Photon_eCorr", Photon_eCorr, &b_Photon_eCorr);
   else
     b_Photon_eCorr = 0;
   fChain->SetBranchAddress("Photon_energyErr", Photon_energyErr, &b_Photon_energyErr);
   fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_hoe", Photon_hoe, &b_Photon_hoe);
   if (!(is22 || is23 || is24)) fChain->SetBranchAddress("Photon_mass", Photon_mass, &b_Photon_mass);
   //fChain->SetBranchAddress("Photon_mvaID", Photon_mvaID, &b_Photon_mvaID);
   //fChain->SetBranchAddress("Photon_mvaID_Fall17V1p1", Photon_mvaID_Fall17V1p1, &b_Photon_mvaID_Fall17V1p1);
   //fChain->SetBranchAddress("Photon_pfRelIso03_all", Photon_pfRelIso03_all, &b_Photon_pfRelIso03_all);
   //fChain->SetBranchAddress("Photon_pfRelIso03_chg", Photon_pfRelIso03_chg, &b_Photon_pfRelIso03_chg);
   fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_r9", Photon_r9, &b_Photon_r9);
   //fChain->SetBranchAddress("Photon_sieie", Photon_sieie, &b_Photon_sieie);
   //fChain->SetBranchAddress("Photon_charge", Photon_charge, &b_Photon_charge);
   if (is17 && isMC && isQCD)
     fChain->SetBranchAddress("Photon_cutBasedBitmap", Photon_cutBased, &b_Photon_cutBasedBitmap);
   else 
     fChain->SetBranchAddress("Photon_cutBased", Photon_cutBased, &b_Photon_cutBased);
   //fChain->SetBranchAddress("Photon_cutBasedBitmap", Photon_cutBased, &b_Photon_cutBasedBitmap);

   //fChain->SetBranchAddress("Photon_cutBased_Fall17V1Bitmap", Photon_cutBased_Fall17V1Bitmap, &b_Photon_cutBased_Fall17V1Bitmap);
   //fChain->SetBranchAddress("Photon_electronIdx", Photon_electronIdx, &b_Photon_electronIdx);
   fChain->SetBranchAddress("Photon_jetIdx", Photon_jetIdx, &b_Photon_jetIdx);
   //fChain->SetBranchAddress("Photon_pdgId", Photon_pdgId, &b_Photon_pdgId);
   //fChain->SetBranchAddress("Photon_vidNestedWPBitmap", Photon_vidNestedWPBitmap, &b_Photon_vidNestedWPBitmap);
   //fChain->SetBranchAddress("Photon_electronVeto", Photon_electronVeto, &b_Photon_electronVeto);
   //fChain->SetBranchAddress("Photon_isScEtaEB", Photon_isScEtaEB, &b_Photon_isScEtaEB);
   //fChain->SetBranchAddress("Photon_isScEtaEE", Photon_isScEtaEE, &b_Photon_isScEtaEE);
   //fChain->SetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80, &b_Photon_mvaID_WP80);
   //fChain->SetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90, &b_Photon_mvaID_WP90);
   //fChain->SetBranchAddress("Photon_pixelSeed", Photon_pixelSeed, &b_Photon_pixelSeed);
   fChain->SetBranchAddress("Photon_seedGain", Photon_seedGain, &b_Photon_seedGain);

   //fChain->SetBranchAddress("RawMET_phi", &RawMET_phi, &b_RawMET_phi);
   //fChain->SetBranchAddress("RawMET_pt", &RawMET_pt, &b_RawMET_pt);
   //fChain->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt, &b_RawMET_sumEt);

   if (isRun3) {
     fChain->SetBranchAddress("RawPuppiMET_phi", &RawPuppiMET_phi, &b_RawPuppiMET_phi);
     fChain->SetBranchAddress("RawPuppiMET_pt", &RawPuppiMET_pt, &b_RawPuppiMET_pt);
     //fChain->SetBranchAddress("RawPuppiMET_sumEt", &RawPuppiMET_sumEt, &b_RawPuppiMET_sumEt);
   }

   if (!isRun3)
     fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   if (isRun3)
     fChain->SetBranchAddress("Rho_fixedGridRhoFastjetAll", &Rho_fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
     //fChain->SetBranchAddress("Rho_fixedGridRhoAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   		fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentral", &Rho_fixedGridRhoFastjetCentral, &b_fixedGridRhoFastjetCentral);
   //fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   		fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentralChargedPileUp", &Rho_fixedGridRhoFastjetCentralChargedPileUp, &b_fixedGridRhoFastjetCentralChargedPileUp);
   //fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);

   //fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   //fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   ////fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   //fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   //fChain->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
   //fChain->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
   //fChain->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);
   //fChain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   //fChain->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
   //fChain->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
   //fChain->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
   //fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   //fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   //fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   //fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   //fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   //fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   //fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   //fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);

   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   //fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   //fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   //fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   //fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   //fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   //fChain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   //fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   //fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   //fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   //fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   //fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   //fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   //fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   if (isRun3) fChain->SetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter, &b_Flag_BadPFMuonDzFilter);
   if (isRun3) fChain->SetBranchAddress("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter, &b_Flag_hfNoisyHitsFilter);
   //fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter, &b_Flag_BadChargedCandidateSummer16Filter);
   //fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter, &b_Flag_BadPFMuonSummer16Filter);
   //fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   //fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   //fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   //fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);

   if (isMC) {
     fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
     fChain->SetBranchAddress("nPSWeight", &nPSWeight, &b_nPSWeight);
     fChain->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight);
     fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);

     fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
     fChain->SetBranchAddress("GenJet_eta", GenJet_eta, &b_GenJet_eta);
     fChain->SetBranchAddress("GenJet_mass", GenJet_mass, &b_GenJet_mass);
     fChain->SetBranchAddress("GenJet_phi", GenJet_phi, &b_GenJet_phi);
     fChain->SetBranchAddress("GenJet_pt", GenJet_pt, &b_GenJet_pt);
     fChain->SetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour, &b_GenJet_partonFlavour);

     if (is22 || is23 || is24)
       fChain->SetBranchAddress("LHE_HT", &LHE_HT, &b_LHE_HT);

     if (!isQCD) {
       fChain->SetBranchAddress("nGenIsolatedPhoton", &nGenIsolatedPhoton, &b_nGenIsolatedPhoton);
       fChain->SetBranchAddress("GenIsolatedPhoton_eta", GenIsolatedPhoton_eta, &b_GenIsolatedPhoton_eta);
       fChain->SetBranchAddress("GenIsolatedPhoton_mass", GenIsolatedPhoton_mass, &b_GenIsolatedPhoton_mass);
       fChain->SetBranchAddress("GenIsolatedPhoton_phi", GenIsolatedPhoton_phi, &b_GenIsolatedPhoton_phi);
       fChain->SetBranchAddress("GenIsolatedPhoton_pt", GenIsolatedPhoton_pt, &b_GenIsolatedPhoton_pt);
     } // isQCD
   } // isMC


   // Safety resets of all branches
   b_HLT_Photon250_NoHE = 0;
   b_HLT_Photon300_NoHE = 0;

   b_HLT_Photon22 = 0;
   b_HLT_Photon30 = 0;
   b_HLT_Photon36 = 0;
   b_HLT_Photon50 = 0;
   b_HLT_Photon75 = 0;
   b_HLT_Photon90 = 0;
   b_HLT_Photon120 = 0;
   b_HLT_Photon150 = 0;
   b_HLT_Photon175 = 0;
   b_HLT_Photon200 = 0;
   //b_HLT_Photon500 = 0;
   //b_HLT_Photon600 = 0;

   b_HLT_Photon20_HoverELoose = 0;
   b_HLT_Photon30_HoverELoose = 0;
   b_HLT_Photon40_HoverELoose = 0;
   b_HLT_Photon50_HoverELoose = 0;
   b_HLT_Photon60_HoverELoose = 0;

   b_HLT_Photon22_R9Id90_HE10_IsoM = 0;
   b_HLT_Photon30_R9Id90_HE10_IsoM = 0;
   b_HLT_Photon36_R9Id90_HE10_IsoM = 0;
   b_HLT_Photon50_R9Id90_HE10_IsoM = 0;
   b_HLT_Photon75_R9Id90_HE10_IsoM = 0;
   b_HLT_Photon90_R9Id90_HE10_IsoM = 0;
   b_HLT_Photon120_R9Id90_HE10_IsoM = 0;
   b_HLT_Photon165_R9Id90_HE10_IsoM = 0;

   b_HLT_Photon165_HE10 = 0;

   b_HLT_Photon50EB_TightID_TightIso = 0;
   //b_HLT_Photon55EB_TightID_TightIso = 0;
   b_HLT_Photon75EB_TightID_TightIso = 0;
   b_HLT_Photon90EB_TightID_TightIso = 0;

   b_HLT_Photon100EB_TightID_TightIso = 0;
   b_HLT_Photon110EB_TightID_TightIso = 0;
   b_HLT_Photon120EB_TightID_TightIso = 0;

   b_HLT_Photon100EBHE10 = 0;
   b_HLT_Photon30EB_TightID_TightIso = 0;

   if (is16) {

   // 250 only in 2016 (2017-2018 has 300_NoHE)
   if (!isMC) {
   fChain->SetBranchAddress("HLT_Photon250_NoHE", &HLT_Photon250_NoHE, &b_HLT_Photon250_NoHE);
   } // !isMC
   fChain->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE, &b_HLT_Photon300_NoHE);

   // 2017-2018 thresholds (vs 2016) 20 (22), 33 (30), x (36)
   if (!isMC) {
   fChain->SetBranchAddress("HLT_Photon22", &HLT_Photon22, &b_HLT_Photon22);
   fChain->SetBranchAddress("HLT_Photon30", &HLT_Photon30, &b_HLT_Photon30);
   fChain->SetBranchAddress("HLT_Photon36", &HLT_Photon36, &b_HLT_Photon36);
   } // !isMC
   fChain->SetBranchAddress("HLT_Photon50", &HLT_Photon50, &b_HLT_Photon50);
   fChain->SetBranchAddress("HLT_Photon75", &HLT_Photon75, &b_HLT_Photon75);
   fChain->SetBranchAddress("HLT_Photon90", &HLT_Photon90, &b_HLT_Photon90);
   fChain->SetBranchAddress("HLT_Photon120", &HLT_Photon120, &b_HLT_Photon120);
   fChain->SetBranchAddress("HLT_Photon175", &HLT_Photon175, &b_HLT_Photon175);

   // only in 2016
   if (is16 && !isMC) {
   fChain->SetBranchAddress("HLT_Photon165_HE10", &HLT_Photon165_HE10, &b_HLT_Photon165_HE10);
   } // !isMC

   // only in 2016 (medium 22,30,36 replaced by loose 20,30 in 2017-2018)
   if (is16 && !isMC) {
     fChain->SetBranchAddress("HLT_Photon22_R9Id90_HE10_IsoM", &HLT_Photon22_R9Id90_HE10_IsoM, &b_HLT_Photon22_R9Id90_HE10_IsoM);
     fChain->SetBranchAddress("HLT_Photon30_R9Id90_HE10_IsoM", &HLT_Photon30_R9Id90_HE10_IsoM, &b_HLT_Photon30_R9Id90_HE10_IsoM);
     fChain->SetBranchAddress("HLT_Photon36_R9Id90_HE10_IsoM", &HLT_Photon36_R9Id90_HE10_IsoM, &b_HLT_Photon36_R9Id90_HE10_IsoM);
   }

   // also in 2017-2018
   fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM, &b_HLT_Photon50_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM, &b_HLT_Photon75_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM, &b_HLT_Photon90_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM, &b_HLT_Photon120_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM, &b_HLT_Photon165_R9Id90_HE10_IsoM);
   
   // only in 2016 (HI pp reference runs?)
   //fChain->SetBranchAddress("HLT_HISinglePhoton10", &HLT_HISinglePhoton10, &b_HLT_HISinglePhoton10);
   //fChain->SetBranchAddress("HLT_HISinglePhoton15", &HLT_HISinglePhoton15, &b_HLT_HISinglePhoton15);
   //fChain->SetBranchAddress("HLT_HISinglePhoton20", &HLT_HISinglePhoton20, &b_HLT_HISinglePhoton20);
   //fChain->SetBranchAddress("HLT_HISinglePhoton40", &HLT_HISinglePhoton40, &b_HLT_HISinglePhoton40);
   //fChain->SetBranchAddress("HLT_HISinglePhoton60", &HLT_HISinglePhoton60, &b_HLT_HISinglePhoton60);

   // also in 2017-2018?
   //fChain->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration, &b_HLT_EcalCalibration);
   //fChain->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration, &b_HLT_HcalCalibration);

   // only in 2016 (backup?)
   //fChain->SetBranchAddress("HLT_Photon500", &HLT_Photon500, &b_HLT_Photon500);
   //fChain->SetBranchAddress("HLT_Photon600", &HLT_Photon600, &b_HLT_Photon600);

   } // is16

   if (is17 || is18 || is22 || is23 || is24) {

   // also in 2016
   fChain->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE, &b_HLT_Photon300_NoHE);

   // 2017-2018 thresholds (vs 2016) 20 (22), 33 (30), x (36)
   fChain->SetBranchAddress("HLT_Photon33", &HLT_Photon33, &b_HLT_Photon33);

   // also in 2016, except 150
   fChain->SetBranchAddress("HLT_Photon50", &HLT_Photon50, &b_HLT_Photon50);
   fChain->SetBranchAddress("HLT_Photon75", &HLT_Photon75, &b_HLT_Photon75);
   fChain->SetBranchAddress("HLT_Photon90", &HLT_Photon90, &b_HLT_Photon90);
   fChain->SetBranchAddress("HLT_Photon120", &HLT_Photon120, &b_HLT_Photon120);
   fChain->SetBranchAddress("HLT_Photon150", &HLT_Photon150, &b_HLT_Photon150);
   fChain->SetBranchAddress("HLT_Photon175", &HLT_Photon175, &b_HLT_Photon175);
   fChain->SetBranchAddress("HLT_Photon200", &HLT_Photon200, &b_HLT_Photon200);

   // also in 2016
   fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM, &b_HLT_Photon50_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM, &b_HLT_Photon75_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM, &b_HLT_Photon90_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM, &b_HLT_Photon120_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM, &b_HLT_Photon165_R9Id90_HE10_IsoM);

   // only in 2017-2018 (loose 20,30 replace medium 22,30,36 in 2016)
   fChain->SetBranchAddress("HLT_Photon20_HoverELoose", &HLT_Photon20_HoverELoose, &b_HLT_Photon20_HoverELoose);
   fChain->SetBranchAddress("HLT_Photon30_HoverELoose", &HLT_Photon30_HoverELoose, &b_HLT_Photon30_HoverELoose);

   } // is17 || is18 || is22 || is23 || is24
   if (is17 && !isMC) {

     fChain->SetBranchAddress("HLT_Photon40_HoverELoose", &HLT_Photon40_HoverELoose, &b_HLT_Photon40_HoverELoose);
     fChain->SetBranchAddress("HLT_Photon50_HoverELoose", &HLT_Photon50_HoverELoose, &b_HLT_Photon50_HoverELoose);
     fChain->SetBranchAddress("HLT_Photon60_HoverELoose", &HLT_Photon60_HoverELoose, &b_HLT_Photon60_HoverELoose);

   } // is17
   if (is18) {

   // only in 2018
   fChain->SetBranchAddress("HLT_Photon20", &HLT_Photon20, &b_HLT_Photon20);

   // only in 2018
   fChain->SetBranchAddress("HLT_Photon100EB_TightID_TightIso", &HLT_Photon100EB_TightID_TightIso, &b_HLT_Photon100EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon110EB_TightID_TightIso", &HLT_Photon110EB_TightID_TightIso, &b_HLT_Photon110EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon120EB_TightID_TightIso", &HLT_Photon120EB_TightID_TightIso, &b_HLT_Photon120EB_TightID_TightIso);

   } // is18
   if (is22 || is23 || is24) {

       //fChain->SetBranchAddress("HLT_Photon20", &HLT_Photon20, &b_HLT_Photon20);
       fChain->SetBranchAddress("HLT_Photon30EB_TightID_TightIso", &HLT_Photon30EB_TightID_TightIso, &b_HLT_Photon30EB_TightID_TightIso);
       fChain->SetBranchAddress("HLT_Photon100EBHE10", &HLT_Photon100EBHE10, &b_HLT_Photon100EBHE10);
       fChain->SetBranchAddress("HLT_Photon110EB_TightID_TightIso", &HLT_Photon110EB_TightID_TightIso, &b_HLT_Photon110EB_TightID_TightIso);

   } // is22 is23 is24
   if (is24){
       //fChain->SetBranchAddress("HLT_Photon50EB_TightID_TightIso_PFJet30", &HLT_Photon50EB_TightID_TightIso_PFJet30, &b_HLT_Photon50EB_TightID_TightIso_PFJet30);
       fChain->SetBranchAddress("HLT_Photon50EB_TightID_TightIso", &HLT_Photon50EB_TightID_TightIso, &b_HLT_Photon50EB_TightID_TightIso);
       //fChain->SetBranchAddress("HLT_Photon55EB_TightID_TightIso", &HLT_Photon55EB_TightID_TightIso, &b_HLT_Photon55EB_TightID_TightIso);
       fChain->SetBranchAddress("HLT_Photon75EB_TightID_TightIso", &HLT_Photon75EB_TightID_TightIso, &b_HLT_Photon75EB_TightID_TightIso);
       fChain->SetBranchAddress("HLT_Photon90EB_TightID_TightIso", &HLT_Photon90EB_TightID_TightIso, &b_HLT_Photon90EB_TightID_TightIso);

			 if(!isMC){//somehow not found for winter2024P8, so only use for data
       	//fChain->SetBranchAddress("HLT_Photon55EB_TightID_TightIso", &HLT_Photon55EB_TightID_TightIso, &b_HLT_Photon55EB_TightID_TightIso);
			 }

       //fChain->SetBranchAddress("HLT_Photon110EB_TightID_TightIso", &HLT_Photon110EB_TightID_TightIso, &b_HLT_Photon110EB_TightID_TightIso);
       //fChain->SetBranchAddress("HLT_Photon110EB_TightID_TightIso_PFJet30", &HLT_Photon110EB_TightID_TightIso_PFJet30, &b_HLT_Photon110EB_TightID_TightIso_PFJet30);
   } // is24

   Notify();
}

Bool_t GamHistosFill::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

  if (fChain && fChain->GetCurrentFile()) {
    _filename = fChain->GetCurrentFile()->GetName();
    if (debugFiles)
      cout << endl << "Opened file: " << _filename << endl << flush;
  }

   return kTRUE;
}

void GamHistosFill::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GamHistosFill::Cut(Long64_t entry)
{
  if (entry) {}; // suppress compiler error
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GamHistosFill_cxx
