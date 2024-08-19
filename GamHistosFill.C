#define GamHistosFill_cxx
#include "GamHistosFill.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"

#include <iostream>
//used for lumimap
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <map>

//typedef for lumimap
typedef std::map<int,double> LumiMap;



using namespace std;
//using namespace GamHistosFill;

#include "parsePileUpJSON.C"

bool _gh_debug = false;
bool _gh_debug100 = false;

bool doGamjet = true;
bool doGamjet2 = true;
bool doJetveto = true; //like in dijet: eta-phi maps
bool smearJets = false;

// Error counters
int cntErrDR(0);

// Classes to structure sets of histograms and profiles
struct BasicHistos {
  TH1D *hn;
  TH1D *hxsec;
  TProfile *prpt;
  TProfile *prbal;
  TProfile *prdb;
  TProfile *prmpf;
  TProfile *prmpf1;
  TProfile *prmpfn;
  TProfile *prmpfu;
  TProfile *prho;
  TProfile *pdjes;
  TProfile *pjes;
  TProfile *pres;
};

class gamjetHistos {
 public:
  
  TH1D *hpt13, *hpt13a, *hpt13j;
  TProfile *pptg, *pptj;
  TProfile *pres, *pjsf, *pm0, *pm2, *pmn, *pmu;

  // Resolution
  TProfile *pm0x, *pm2x;

  // Extra FRS studies
  TProfile *pmnu, *pmnx, *pmux, *pmnux;
  
  // Composition
  TProfile *prho, *pchf, *pnef, *pnhf;

  // Alternative pT bins
  TProfile *presa, *pm0a, *pm2a, *pmna, *pmua;
  TProfile *presj, *pm0j, *pm2j, *pmnj, *pmuj;
};

class gamjetHistos2 {
public:

  // Basic information about the trigger
  //string trg;
  //int trgpt;
  //double ptmin, ptmax, absetamin, absetamax;

  TH2D *h2pteta;
  TProfile2D *p2res, *p2m0, *p2m2, *p2mn, *p2mu;
  TProfile2D *p2m0x, *p2m2x;

  // Extra for FSR studies
  TProfile2D *p2mnu, *p2mnx, *p2mux, *p2mnux;
  //TH2D *h2ptetatc, *h2ptetapf;
  //TProfile2D *p2restc, *p2m0tc, *p2m2tc, *p2mntc, *p2mutc; // pT,tag (central)
  //TProfile2D *p2respf, *p2m0pf, *p2m2pf, *p2mnpf, *p2mupf; // pT,probe (forward)

  // Smearing controls
  TProfile2D *p2jsf;//, *p2jsftc, *p2jsfpf;
};

// Helper function to retrieve FactorizedJetCorrector 
FactorizedJetCorrector *getFJC(string l1="", string l2="", string res="",
			       string path="") {

  // Set default jet algo                                                       
  if (l1!="" && !(TString(l1.c_str()).Contains("_AK")))
    l1 += "_AK4PFchs";
  if (l2!="" && !(TString(l2.c_str()).Contains("_AK")))
    l2 += "_AK4PFchs";
  if (res!="" && !(TString(res.c_str()).Contains("_AK")))
    res += "_AK4PFchs";

  // Set default path
  if (path=="") path = "CondFormats/JetMETObjects/data";
  const char *cd = path.c_str();
  const char *cl1 = l1.c_str();
  const char *cl2 = l2.c_str();
  const char *cres = res.c_str();
  string s("");

  vector<JetCorrectorParameters> v;
  if (l1!=""){
    s = Form("%s/%s.txt",cd,cl1);
    cout << s << endl << flush;
    JetCorrectorParameters *pl1 = new JetCorrectorParameters(s);
    v.push_back(*pl1);
  }
  if (l2!="") {
    s = Form("%s/%s.txt",cd,cl2);
    cout << s << endl << flush;
    JetCorrectorParameters *pl2 = new JetCorrectorParameters(s);
    v.push_back(*pl2);
  }
  if (res!="") {
    s = Form("%s/%s.txt",cd,cres);
    cout << s << endl << flush;
    JetCorrectorParameters *pres = new JetCorrectorParameters(s);
    v.push_back(*pres);
  }
  FactorizedJetCorrector *jec = new FactorizedJetCorrector(v);

  return jec;
} // getJFC      

void GamHistosFill::Loop()
{
//   In a ROOT session, you can do:
//      root> .L GamHistosFill.C
//      root> GamHistosFill t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

  //Not working in v6.26.06 as is
  //ROOT.EnableImplicitMT(); // From Nico on Skype, to parallelize processing

  TStopwatch fulltime, laptime;
  fulltime.Start();
  TDatime bgn;
  int nlap(0);

  int _ntot(0), _nevents(0), _nbadevents_json(0), _nbadevents_trigger(0);
  int _nbadevents_veto(0);
  
  if (true) { // ProcessFast
    fChain->SetBranchStatus("*",0);  // disable all branches
    
    // Baseline triggers with very high prescale except Photon200
    if (is18) {// || is22 || is23) {
      fChain->SetBranchStatus("HLT_Photon20",1);
    }
    if (is16) {
      if (!isMC) fChain->SetBranchStatus("HLT_Photon22",1);
      if (!isMC) fChain->SetBranchStatus("HLT_Photon30",1);
      if (!isMC) fChain->SetBranchStatus("HLT_Photon36",1);
    }
    if (is17 || is18 || is22 || is23 || is24) {
      fChain->SetBranchStatus("HLT_Photon33",1);
      fChain->SetBranchStatus("HLT_Photon200",1);
    //fChain->SetBranchStatus("HLT_Photon500",1);
    //fChain->SetBranchStatus("HLT_Photon600",1);
    }
    fChain->SetBranchStatus("HLT_Photon50",1);
    fChain->SetBranchStatus("HLT_Photon75",1);
    fChain->SetBranchStatus("HLT_Photon90",1);
    fChain->SetBranchStatus("HLT_Photon120",1);
    fChain->SetBranchStatus("HLT_Photon175",1);

    // Presumably backup triggers for unprescaled Photon175 or Photon200
    if (is16 && !isMC) fChain->SetBranchStatus("HLT_Photon165_HE10",1);
    if (is16 && !isMC) fChain->SetBranchStatus("HLT_Photon250_NoHE",1);
    fChain->SetBranchStatus("HLT_Photon300_NoHE",1);
    
    // Triggers to recover pT=105-230 GeV range. Efficient, low prescale
    // Not active in some of 2018 data (~10%)?
    if (is18) {
      fChain->SetBranchStatus("HLT_Photon100EB_TightID_TightIso",1); // ok
      fChain->SetBranchStatus("HLT_Photon110EB_TightID_TightIso",1); // best
      fChain->SetBranchStatus("HLT_Photon120EB_TightID_TightIso",1); // backup
    }
    if (isRun3) {
      fChain->SetBranchStatus("HLT_Photon30EB_TightID_TightIso",1);
      fChain->SetBranchStatus("HLT_Photon110EB_TightID_TightIso",1);
      fChain->SetBranchStatus("HLT_Photon100EBHE10",1);
    }
    if (is24) { //new 2024 trigger paths
      fChain->SetBranchStatus("HLT_Photon50EB_TightID_TightIso",1);
      //fChain->SetBranchStatus("HLT_Photon55EB_TightID_TightIso",1);
      fChain->SetBranchStatus("HLT_Photon75EB_TightID_TightIso",1);
      fChain->SetBranchStatus("HLT_Photon90EB_TightID_TightIso",1);
			if(!isMC){ //did not work for winter2024P8, so do it only for data
	      fChain->SetBranchStatus("HLT_Photon55EB_TightID_TightIso",1);
			}
    }

    // Triggers to recover 60-105 GeV range. However, inefficient up to high pT
    // Possibly medium HLT cuts not fully consistent with tight ID?
    if (is16 && !isMC) {
      fChain->SetBranchStatus("HLT_Photon22_R9Id90_HE10_IsoM",1);
      fChain->SetBranchStatus("HLT_Photon30_R9Id90_HE10_IsoM",1);
      fChain->SetBranchStatus("HLT_Photon36_R9Id90_HE10_IsoM",1);
    }
    fChain->SetBranchStatus("HLT_Photon50_R9Id90_HE10_IsoM",1);
    fChain->SetBranchStatus("HLT_Photon75_R9Id90_HE10_IsoM",1);
    fChain->SetBranchStatus("HLT_Photon90_R9Id90_HE10_IsoM",1);
    fChain->SetBranchStatus("HLT_Photon120_R9Id90_HE10_IsoM",1);
    fChain->SetBranchStatus("HLT_Photon165_R9Id90_HE10_IsoM",1);
    
    // Triggers to recover 30-60 GeV range. Efficient above 30 GeV
    if (is17 || is18 || is22 || is23 || is24) {
      fChain->SetBranchStatus("HLT_Photon20_HoverELoose",1);
      fChain->SetBranchStatus("HLT_Photon30_HoverELoose",1);
    }
    if (is17 && !isMC) {
      fChain->SetBranchStatus("HLT_Photon40_HoverELoose",1);
      fChain->SetBranchStatus("HLT_Photon50_HoverELoose",1);
      fChain->SetBranchStatus("HLT_Photon60_HoverELoose",1);
    }
    
    // JSON filtering, PU reweighing
    if (isMC) fChain->SetBranchStatus("Pileup_nTrueInt");
    fChain->SetBranchStatus("run",1);
    if (!isMC) fChain->SetBranchStatus("luminosityBlock",1);
    if (!isMC && debugFiles) fChain->SetBranchStatus("event",1);

    if (isRun3) {
      // Same cut as for dijet sample
      //fChain->SetBranchStatus("Flag_METFilters",1);
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Run_3_recommendations
      fChain->SetBranchStatus("Flag_goodVertices");
      fChain->SetBranchStatus("Flag_globalSuperTightHalo2016Filter");
      fChain->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter");
      fChain->SetBranchStatus("Flag_BadPFMuonFilter");
      fChain->SetBranchStatus("Flag_BadPFMuonDzFilter");
      fChain->SetBranchStatus("Flag_hfNoisyHitsFilter");
      fChain->SetBranchStatus("Flag_eeBadScFilter");
      fChain->SetBranchStatus("Flag_ecalBadCalibFilter");
    }
    else {    
    // Event filters listing from Sami:
    // twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
    // #2018_2017_data_and_MC_UL
    fChain->SetBranchStatus("Flag_goodVertices",1);
    fChain->SetBranchStatus("Flag_globalSuperTightHalo2016Filter",1);
    fChain->SetBranchStatus("Flag_HBHENoiseFilter",1);
    fChain->SetBranchStatus("Flag_HBHENoiseIsoFilter",1);
    fChain->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter",1);
    fChain->SetBranchStatus("Flag_BadPFMuonFilter",1);
    //fChain->SetBranchStatus("Flag_BadPFMuonDzFilter",1); // not in nAOD?
    //inputList.("Flag_BadPFMuonDzFilter",1); // new in UL, but part of above?
    //inputList.("Flag_BadChargedCandidateFilter",0); // not recommended
    //inputList.("Flag_globalTightHalo2016Filter",0); // obsolete?
    //inputList.("Flag_CSCTightHaloFilter",0); // obsolete?
    if (!is16) {
      fChain->SetBranchStatus("Flag_ecalBadCalibFilter",1);//new (add for Sami)
    }
    //if dataset.isData:
    //if (!isMC) {
    fChain->SetBranchStatus("Flag_eeBadScFilter",1); // MC added 7 July 2021
    //}
    }    

    // MC weights
    if (isMC)  fChain->SetBranchStatus("genWeight",1);
    if (isMC)  fChain->SetBranchStatus("nPSWeight",1);
    if (isMC)  fChain->SetBranchStatus("PSWeight",1);
    //if (isMC && !isRun3)  fChain->SetBranchStatus("nPSWeight",1);
    //if (isMC && !isRun3)  fChain->SetBranchStatus("PSWeight",1);

    if (!isRun3)
      fChain->SetBranchStatus("fixedGridRhoFastjetAll",1);
    if (isRun3)
      fChain->SetBranchStatus("Rho_fixedGridRhoFastjetAll",1);
      //fChain->SetBranchStatus("Rho_fixedGridRhoAll",1);
    fChain->SetBranchStatus("PV_npvs",1);
    fChain->SetBranchStatus("PV_npvsGood",1);

    if (isRun2) {
      fChain->SetBranchStatus("ChsMET_pt",1);
      fChain->SetBranchStatus("ChsMET_phi",1);
    }
    if (isRun3) {
      fChain->SetBranchStatus("RawPuppiMET_pt",1);
      fChain->SetBranchStatus("RawPuppiMET_phi",1);
    }
    
    fChain->SetBranchStatus("nPhoton",1);
    fChain->SetBranchStatus("Photon_pt",1);
    fChain->SetBranchStatus("Photon_eta",1);
    fChain->SetBranchStatus("Photon_phi",1);
    if (!(is22||is23||is24)) fChain->SetBranchStatus("Photon_mass",1);
    fChain->SetBranchStatus("Photon_hoe",1);
    if (is17 && isMC && isQCD)
      fChain->SetBranchStatus("Photon_cutBasedBitmap",1);
    else
      fChain->SetBranchStatus("Photon_cutBased",1);
    fChain->SetBranchStatus("Photon_jetIdx",1);

    fChain->SetBranchStatus("Photon_seedGain",1);
    if (b_Photon_eCorr)
      fChain->SetBranchStatus("Photon_eCorr",1); // not in 2016
    fChain->SetBranchStatus("Photon_energyErr",1);
    fChain->SetBranchStatus("Photon_r9",1);
    
    fChain->SetBranchStatus("nJet",1);
    fChain->SetBranchStatus("Jet_pt",1);
    fChain->SetBranchStatus("Jet_eta",1);
    fChain->SetBranchStatus("Jet_phi",1);
    fChain->SetBranchStatus("Jet_mass",1);
    fChain->SetBranchStatus("Jet_rawFactor",1);
    fChain->SetBranchStatus("Jet_area",1);
    fChain->SetBranchStatus("Jet_jetId",1);

    // PF composition
    fChain->SetBranchStatus("Jet_chHEF",1);
    fChain->SetBranchStatus("Jet_neHEF",1);
    fChain->SetBranchStatus("Jet_neEmEF",1);
    fChain->SetBranchStatus("Jet_chEmEF",1);
    fChain->SetBranchStatus("Jet_muEF",1);
    //if (!isRun3) fChain->SetBranchStatus("Jet_chFPV0EF",1);

    //new: HF stuff (branches only available in 2024 onwards)
    if (is24) fChain->SetBranchStatus("Jet_hfEmEF",1);    //electromagnetic Energy Fraction in HF
    if (is24) fChain->SetBranchStatus("Jet_hfHEF",1);     //hadronic Energy Fraction in HF

    if (isMC) fChain->SetBranchStatus("Jet_genJetIdx",1);
    
    if (!isRun3) fChain->SetBranchStatus("Jet_btagDeepB",1);
    if (!isRun3) fChain->SetBranchStatus("Jet_btagDeepC",1);
    if (!isRun3) fChain->SetBranchStatus("Jet_qgl",1);
    //
    if (isRun3) fChain->SetBranchStatus("Jet_btagDeepFlavB",1);
    if (isRun3) fChain->SetBranchStatus("Jet_btagDeepFlavCvB",1);
    if (isRun3) fChain->SetBranchStatus("Jet_btagDeepFlavCvL",1);
    if (isRun3) fChain->SetBranchStatus("Jet_btagDeepFlavQG",1);

    //if (isMC && isRun3) {
    if (isMG) {
      fChain->SetBranchStatus("LHE_HT",1);
      cout << "Running over MadGraph HT bins" << endl;
    }
    /*
      fChain->SetBranchStatus("nLHEPart",1);
      fChain->SetBranchStatus("LHEPart_pt",1);
      fChain->SetBranchStatus("LHEPart_eta",1);
      fChain->SetBranchStatus("LHEPart_phi",1);
      fChain->SetBranchStatus("Jet_mass",1);
      fChain->SetBranchStatus("LHEPart_status",1);
      fChain->SetBranchStatus("LHEPart_pdgId",1);
    */
    
    // Photon pT correlates with Generator_binvar (pThat) in P8 PtFlat,
    // and with LHEPart_pt[2] (parton photon) in MadGraph HT samples, but
    // MadGraph has plenty of uncorrelated soft and hard photons as well
    if (isMC) {
      if (!isQCD) {
	fChain->SetBranchStatus("nGenIsolatedPhoton",1);
	fChain->SetBranchStatus("GenIsolatedPhoton_pt",1);
	fChain->SetBranchStatus("GenIsolatedPhoton_eta",1);
	fChain->SetBranchStatus("GenIsolatedPhoton_phi",1);
	fChain->SetBranchStatus("GenIsolatedPhoton_mass",1);
      }
      else
	nGenIsolatedPhoton = 0;
    }
    if (isMC) {
      fChain->SetBranchStatus("nGenJet",1);
      fChain->SetBranchStatus("GenJet_pt",1);
      fChain->SetBranchStatus("GenJet_eta",1);
      fChain->SetBranchStatus("GenJet_phi",1);
      fChain->SetBranchStatus("GenJet_mass",1);
      fChain->SetBranchStatus("GenJet_partonFlavour",1);
    }
  }
  
  // Select appropriate L1RC for type-I MET L1L2L3-RC calculation
  FactorizedJetCorrector *jecl1rc(0), *jec(0);
  string& ds = dataset;
  if (ds=="2016B")    jecl1rc = getFJC("Summer19UL16APV_RunBCD_V7_DATA_L1RC");
  if (ds=="2016C")    jecl1rc = getFJC("Summer19UL16APV_RunBCD_V7_DATA_L1RC");
  if (ds=="2016D")    jecl1rc = getFJC("Summer19UL16APV_RunBCD_V7_DATA_L1RC");
  if (ds=="2016BCD")  jecl1rc = getFJC("Summer19UL16APV_RunBCD_V7_DATA_L1RC");
  if (ds=="2016E")    jecl1rc = getFJC("Summer19UL16APV_RunEF_V7_DATA_L1RC");
  if (ds=="2016F")    jecl1rc = getFJC("Summer19UL16APV_RunEF_V7_DATA_L1RC");
  if (ds=="2016EF")   jecl1rc = getFJC("Summer19UL16APV_RunEF_V7_DATA_L1RC");
  //if (ds=="2016BCDEF")jecl1rc=getFJC("Summer19UL16APV_RunBCDEF_V7_DATA_L1RC");

  if (ds=="2016FG")   jecl1rc = getFJC("Summer19UL16_RunFGH_V7_DATA_L1RC");
  if (ds=="2016H")    jecl1rc = getFJC("Summer19UL16_RunFGH_V7_DATA_L1RC");
  if (ds=="2016FGH")  jecl1rc = getFJC("Summer19UL16_RunFGH_V7_DATA_L1RC");

  if (ds=="2017B")  jecl1rc = getFJC("Summer19UL17_RunB_V5_DATA_L1RC");
  if (ds=="2017C")  jecl1rc = getFJC("Summer19UL17_RunC_V5_DATA_L1RC");
  if (ds=="2017D")  jecl1rc = getFJC("Summer19UL17_RunD_V5_DATA_L1RC");
  if (ds=="2017E")  jecl1rc = getFJC("Summer19UL17_RunE_V5_DATA_L1RC");
  if (ds=="2017F")  jecl1rc = getFJC("Summer19UL17_RunF_V5_DATA_L1RC");
  //if (ds=="2017BCDEF")jecl1rc=getFJC("Summer19UL17_RunBCDEF_V5_DATA_L1RC");

  if (ds=="2018A")  jecl1rc = getFJC("Summer19UL18_RunA_V5_DATA_L1RC");
  if (ds=="2018B")  jecl1rc = getFJC("Summer19UL18_RunB_V5_DATA_L1RC");
  if (ds=="2018C")  jecl1rc = getFJC("Summer19UL18_RunC_V5_DATA_L1RC");
  if (ds=="2018D")  jecl1rc = getFJC("Summer19UL18_RunD_V5_DATA_L1RC");

  if (ds=="2018A1") jecl1rc = getFJC("Summer19UL18_RunA_V5_DATA_L1RC");
  if (ds=="2018A2") jecl1rc = getFJC("Summer19UL18_RunA_V5_DATA_L1RC");
  if (ds=="2018D1") jecl1rc = getFJC("Summer19UL18_RunD_V5_DATA_L1RC");
  if (ds=="2018D2") jecl1rc = getFJC("Summer19UL18_RunD_V5_DATA_L1RC");
  if (ds=="2018D3") jecl1rc = getFJC("Summer19UL18_RunD_V5_DATA_L1RC");
  if (ds=="2018D4") jecl1rc = getFJC("Summer19UL18_RunD_V5_DATA_L1RC");
  //if (ds=="2018ABCD")jecl1rc=getFJC("Summer19UL18_RunABCD_V5_DATA_L1RC");

  if (isRun3) {
    //if (isMC)
    //jecl1rc = getFJC("Winter23Prompt23_V2_MC_L1RC_AK4PFchs");
    //else
    //jecl1rc = getFJC("Winter23Prompt23_RunC_V2_DATA_L1RC_AK4PFchs");
    jecl1rc = 0;
  }
    
  if (ds=="2016P8")    jecl1rc = getFJC("Summer19UL16_V7_MC_L1RC");
  if (ds=="2016APVP8") jecl1rc = getFJC("Summer19UL16APV_V7_MC_L1RC");
  if (ds=="2017P8")    jecl1rc = getFJC("Summer19UL17_V5_MC_L1RC");
  if (ds=="2018P8")    jecl1rc = getFJC("Summer19UL18_V5_MC_L1RC");
  //
  if (ds=="2016QCD")    jecl1rc = getFJC("Summer19UL16_V7_MC_L1RC");
  if (ds=="2016APVQCD") jecl1rc = getFJC("Summer19UL16APV_V7_MC_L1RC");
  if (ds=="2017QCD")    jecl1rc = getFJC("Summer19UL17_V5_MC_L1RC");
  if (ds=="2018QCD")    jecl1rc = getFJC("Summer19UL18_V5_MC_L1RC");
  assert(jecl1rc || isRun3);

  // FactorizedJetCorrector for redoing JEC on the fly.
  //2022
  if (ds=="2022C") {
    jec = getFJC("", "Summer22Run3_V1_MC_L2Relative_AK4PUPPI",
		 //"Run22CD-22Sep2023_DATA_L2L3Residual_AK4PFPuppi");
		 "Summer22-22Sep2023_Run2022CD_V3_DATA_L2L3Residual_AK4PFPuppi");
  }
  if (ds=="2022D") {
    jec = getFJC("", "Summer22Run3_V1_MC_L2Relative_AK4PUPPI",
		 //"Run22CD-22Sep2023_DATA_L2L3Residual_AK4PFPuppi");
		 "Summer22-22Sep2023_Run2022CD_V3_DATA_L2L3Residual_AK4PFPuppi");
  }
  if (ds=="2022E") {
    jec = getFJC("", "Summer22EEVetoRun3_V1_MC_L2Relative_AK4PUPPI",
		 //"Run22E-22Sep2023_DATA_L2L3Residual_AK4PFPuppi");
		 "Summer22EE-22Sep2023_Run2022E_V3_DATA_L2L3Residual_AK4PFPuppi");		 
  }
  if (ds=="2022F") {
    jec = getFJC("", "Summer22EEVetoRun3_V1_MC_L2Relative_AK4PUPPI",
		 //"Run22F-Prompt_DATA_L2L3Residual_AK4PFPuppi");
		 "Summer22EEPrompt22_Run2022F_V3_DATA_L2L3Residual_AK4PFPuppi");
  }
  if (ds=="2022G") {
    jec = getFJC("", "Summer22EEVetoRun3_V1_MC_L2Relative_AK4PUPPI",
		 //"Run22G-Prompt_DATA_L2L3Residual_AK4PFPuppi");
		 "Summer22EEPrompt22_Run2022G_V3_DATA_L2L3Residual_AK4PFPuppi");
  }
  //22/23 MC
  if (ds=="2022P8" || ds=="2022QCD") {
    jec = getFJC("", "Summer22Run3_V1_MC_L2Relative_AK4PUPPI","");
  }
  if (ds=="2022EEP8" || ds=="2022EEQCD") {
    jec = getFJC("", "Summer22EEVetoRun3_V1_MC_L2Relative_AK4PUPPI", "");
  }
 /*
  if (dataset=="Summer23") { //includes thus both BPix corrected ones and the ones without it
    //got corrections from here: https://github.com/cms-jet/JECDatabase/tree/master/textFiles/Winter23Prompt23_V2_MC
    jec = getFJC("", "Winter23Prompt23_V2_MC_L2Relative_AK4PFPuppi", "");// use this for Summer23 MC (get this from github)
    //assert(false); // not yet available --> use the Winter23Prompt23_V2_MC_L2Relative_AK4PFPuppi, which is available.
  }
 */
  //MC2023 --> added for running 2023MC with and without BPix stuff
  if (ds=="2023P8" || ds=="2023QCD" || ds=="2023P8X") { //earlier called: Summer2023
    //jec = getFJC("", "Summer22Run3_V1_MC_L2Relative_AK4PUPPI", ""); //23rd of Feb2024 - investigating how plots look with 2022 corrections
    jec = getFJC("", "Summer23Run3_V1_MC_L2Relative_AK4PUPPI", ""); //16th of Feb2024, w4, w5 and onwards
		    //"Summer23Run3_V1_MC_L2Relative_AK4PUPPI", ""); //only used MC L2Rel for w2
		    //"Winter23Prompt23_V2_MC_L2Relative_AK4PFPuppi","");
    //assert(false); // not yet available --> use the Winter23Prompt23_V2_MC_L2Relative_AK4PFPuppi, which is available.
  }
  if (ds=="2023P8-BPix" || ds=="2023QCD-BPix" || ds=="2023P8-BPixX") { //earlier called: Summer2023, BPix separately!!
    //jec = getFJC("", "Summer22Run3_V1_MC_L2Relative_AK4PUPPI", ""); //23rd of Feb2024 - investigating how plots look with 2022 corrections
    jec = getFJC("", "Summer23BPixRun3_V3_MC_L2Relative_AK4PUPPI", "" ); //16th of Feb2024, w4, w5 and onwards 
		    //"Summer23BPixRun3_V3_MC_L2Relative_AK4PUPPI", ""); //only used MC L2Rel for w2
		    //"Winter23Prompt23_V2_MC_L2Relative_AK4PFPuppi",""); //old
  }
  //data2023
  if (ds=="2023B" || ds=="2023Cv123" || ds=="2023Cv123X") {//2023C --> no bpix issue
	//got Winter23 corrections from here: https://github.com/cms-jet/JECDatabase/tree/master/textFiles/Winter23Prompt23_V2_MC
    //jec = getFJC("", "Summer22Run3_V1_MC_L2Relative_AK4PUPPI", "Summer22Prompt23_Run2023Cv123_V3_DATA_L2L3Residual_AK4PFPUPPI"); //23rd of Feb2024 - investigating how plots look with 2022 corrections
    jec = getFJC("", "Summer23Run3_V1_MC_L2Relative_AK4PUPPI", 
			"Summer23Prompt23_Run2023Cv123_V2_DATA_L2L3Residual_AK4PFPuppi"); //29th of Feb2024, w7 and onwards
		  	//"Summer23Prompt23_Run2023Cv123_V1_DATA_L2L3Residual_AK4PFPuppi"); //18th of Feb2024, w5
		    	//"Summer23Prompt23_Run2023Cv123_V1_DATA_L2Residual_AK4PFPuppi"); //16th of Feb2024, w4; and 27th of Feb2024, w6
		 //"Summer23Run3_V1_MC_L2Relative_AK4PUPPI", ""); //only used MC L2Rel for w2
		 //"Winter23Prompt23_V2_MC_L2Relative_AK4PFPuppi", ""); //old
		 //"Summer22Prompt23_Run2023Cv123_V3_DATA_L2L3Residual_AK4PFPUPPI"); //even older
  }
  if (ds=="2023Cv4" || ds=="2023Cv4X") {//2023C --> no bpix issue
    //jec = getFJC("", "Summer22Run3_V1_MC_L2Relative_AK4PUPPI", "Summer22Prompt23_Run2023Cv4_V3_DATA_L2L3Residual_AK4PFPUPPI"); //23rd of Feb2024 - investigating how plots look with 2022 corrections
    jec = getFJC("","Summer23Run3_V1_MC_L2Relative_AK4PUPPI",  
			"Summer23Prompt23_Run2023Cv4_V2_DATA_L2L3Residual_AK4PFPuppi"); //29th of Feb2024, w7 and onwards
		    	//"Summer23Prompt23_Run2023Cv4_V1_DATA_L2L3Residual_AK4PFPuppi"); //18th of Feb2024, w5
		    	//"Summer23Prompt23_Run2023Cv4_V1_DATA_L2Residual_AK4PFPuppi"); //16th of Feb2024, w4; and 27th of Feb2024, w6
		 //"Summer23Run3_V1_MC_L2Relative_AK4PUPPI", ""); //only used MC L2Rel for w2
		 //"Winter23Prompt23_V2_MC_L2Relative_AK4PFPuppi", ""); //old
		 //"Summer22Prompt23_Run2023Cv4_V3_DATA_L2L3Residual_AK4PFPUPPI"); //even older
  }
  if (ds=="2023D" || ds=="2023DX") { //2023D needs BPix stuff!
    //jec = getFJC("", "Summer22Run3_V1_MC_L2Relative_AK4PUPPI", "Summer22Prompt23_Run2023D_V3_DATA_L2L3Residual_AK4PFPUPPI"); //23rd of Feb2024 - investigating how plots look with 2022 corrections
    jec = getFJC("", "Summer23BPixRun3_V3_MC_L2Relative_AK4PUPPI",  
                        "Summer23Prompt23_Run2023D_V2_DATA_L2L3Residual_AK4PFPuppi"); //9th of Mar2024, w8 (fixed this...)
			//"Summer23Prompt23_Run2023Cv4_V2_DATA_L2L3Residual_AK4PFPuppi"); //29th of Feb2024, w7
		    	//"Summer23Prompt23_Run2023D_V1_DATA_L2L3Residual_AK4PFPuppi"); //18th of Feb2024, w5
		    	//"Summer23Prompt23_Run2023D_V1_DATA_L2Residual_AK4PFPuppi"); //16th of Feb2024, w4; and 27th of Feb2024, w6
		 //"Summer23BPixRun3_V3_MC_L2Relative_AK4PUPPI", ""); //only used MC L2Rel for w2
		 //"Winter23Prompt23_V2_MC_L2Relative_AK4PFPuppi", ""); //old 
		 //"Summer22Prompt23_Run2023D_V3_DATA_L2L3Residual_AK4PFPUPPI"); //even older
  }
	//mc2024
  if (ds=="winter2024P8" || ds=="2024QCD") { //7th of Aug2024, w32 onwards; 14.8. for QCD w33
    jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "" ); //use this?
  }
  //data2024
  if (ds=="2024B-PromptReco-v1" || ds=="2024B" || ds=="2024C" || ds=="2024D") { //2023D needs BPix stuff, use this also for 2024B prompt data (12.4.24)
			jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Prompt24_Run2024BCD_V5M_DATA_L2L3Residual_AK4PFPuppi"); //w34 onwards (16.08.2024), V5M
    //jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Prompt24_Run2024BCD_V4M_DATA_L2L3Residual_AK4PFPuppi"); //w30 onwards (14.06.2024), V4M
    //jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Prompt24_Run2024BCD_V3M_DATA_L2L3Residual_AK4PFPuppi"); //w27 and w28 (starting 03.06.24) and onwards
    //jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Prompt24_Run2024BC_V2M_DATA_L2L3Residual_AK4PFPuppi"); //w17 and w18 (starting 10.05.24) and onwards
    //jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Prompt24_Run2024BC_V1M_DATA_L2L3Residual_AK4PFPuppi"); //w15 and w16 (starting 06.05.24)
    //jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Summer23BPixPrompt23_RunD_V1_DATA_L2L3Residual_AK4PFPuppi"); //Winter2024 L2Rel, and 2023D-L2L3Res (w12, 30.04.2024)
    //jec = getFJC("", "Summer23BPixPrompt23_V1_MC_L2Relative_AK4PFPuppi", "Summer23BPixPrompt23_RunD_V1_DATA_L2L3Residual_AK4PFPuppi"); //took the official ones from: (the one with V2 was an internal one from Mikko) --> should update also for 2023 stuff above (TO DO).
    //jec = getFJC("", "Summer23BPixRun3_V3_MC_L2Relative_AK4PUPPI", "Summer23Prompt23_Run2023D_V2_DATA_L2L3Residual_AK4PFPuppi"); //9th of Mar2024, w8 (fixed this...)
  }
  if (ds=="2024Ev1" || ds=="2024Ev2"){ //separated 24E corrections starting from V4M jecs
			jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Prompt24_Run2024E_V5M_DATA_L2L3Residual_AK4PFPuppi"); //w34
    //jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Prompt24_Run2024E_V4M_DATA_L2L3Residual_AK4PFPuppi"); //w30
	}
	if (ds=="2024B-ECALRATIO" || ds=="2024C-ECALRATIO") { //|| ds=="2024C-ECALR-HCALDI" ) { //for 2024 re-reco data
    jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Prompt24_Run2024CR_V4M_DATA_L2L3Residual_AK4PFPuppi"); //w30 (starting 14.06.24) and onwards, V4M
    //jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Prompt24_Run2024CR_V3M_DATA_L2L3Residual_AK4PFPuppi"); //w27 and w28 (starting 03.06.24) and onwards
	}
	if (ds=="2024C-ECALR-HCALDI" || ds=="2024F" || ds=="2024G") { //for 2024 re-reco data, but also for 2024F and onwards (fixed in w33)
			jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Prompt24_Run2024F_V5M_DATA_L2L3Residual_AK4PFPuppi"); //w34 (starting 16.08.24) and onwards, V5M
    //jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Prompt24_Run2024CS_V4M_DATA_L2L3Residual_AK4PFPuppi"); //w30 (starting 14.06.24) and onwards, V4M
		//jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Prompt24_Run2024CR_V3M_DATA_L2L3Residual_AK4PFPuppi"); 
	}
  assert(jec);

  
  string sera("");
  if (ds=="2016APVP8" || ds=="2016APVQCD") sera = "2016APV";
  if (ds=="2016P8" || ds=="2016QCD") sera = "2016FGH";
  if (ds=="2017P8" || ds=="2017QCD") sera = "2017";
  if (ds=="2018P8" || ds=="2018QCD") sera = "2018";
  if (ds=="2022P8" || ds=="2022QCD") sera = "2022";
  if (ds=="2022EEP8" || ds=="2022EEQCD") sera = "2022EE";
  if (ds=="2023P8" || ds=="2023QCD" || ds=="2023P8-BPix" || ds=="2023QCD-BPix") sera = "2023"; //added 2023P8-BPix
  //
  if (ds=="2016B"||ds=="2016C"||ds=="2016D"||ds=="2016BCD"||
      ds=="2016E"||ds=="2016F"||ds=="2016EF"||ds=="2016BCDEF") sera = "2016APV";
  if (ds=="2016FG"||ds=="2016H"||ds=="2016FGH") sera = "2016FGH";
  if (ds=="2017B"||ds=="2017C"||ds=="2017D"||ds=="2017E"||ds=="2017F"||
      ds=="2017BCDEF") sera = "2017";
  if (ds=="2018A"||ds=="2018A1"||ds=="2018A2"||ds=="2018B"||ds=="2018C"||
      ds=="2018D"||ds=="2018D1"||ds=="2018D2"||ds=="2018D3"||ds=="2018D4")
    sera = "2018";
  //
  if (ds=="2022P8" || ds=="2022QCD") sera = "2022";
  if (ds=="2022EEP8" || ds=="2022EEQCD") sera = "2022EE";
  if (ds=="2023P8" || ds=="2023QCD" || ds=="2023P8-BPix" || ds=="2023QCD-BPix") sera = "2023"; //added 2023P8-BPix
  if (ds=="2023P8X" || ds=="2023QCDX" || ds=="2023P8-BPixX" || ds=="2023QCD-BPixX") sera = "2023"; //added for w23X and w22X
  if (ds=="2022C" || ds=="2022D") sera ="2022";
  if (ds=="2022E" || ds=="2022F" || ds=="2022G") sera = "2022EE";
  if (ds=="2023B" || ds=="2023Cv123" || ds=="2023Cv4" || ds=="2023D") sera = "2023";
  if (ds=="2023Cv123X" || ds=="2023Cv4X" || ds=="2023DX") sera = "2023";
  if (ds=="2024B-PromptReco-v1" || ds=="2024B" || ds=="2024C" || ds=="2024D" || ds=="2024Ev1" || ds=="2024Ev2" || ds=="2024F" || ds=="2024G" ||
			ds=="2024B-ECALRATIO" || ds=="2024C-ECALRATIO" || ds=="2024C-ECALR-HCALDI") sera = "2024";
  if (ds=="winter2024P8" || ds=="2024QCD" || ds=="2024P8") sera = "2024"; //currently only winter2024P8 in use (w32), now also QCD (w33)
  assert(sera!="");

  // Load JSON files
  if (TString(ds.c_str()).Contains("2016"))
    LoadJSON("files/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt");
  if (TString(ds.c_str()).Contains("2017"))
    LoadJSON("files/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt");
  if (TString(ds.c_str()).Contains("2018"))
    LoadJSON("files/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt");
  if (TString(ds.c_str()).Contains("2022"))
    LoadJSON("files/Cert_Collisions2022_355100_362760_Golden.json");
  if (TString(ds.c_str()).Contains("2023"))
    LoadJSON("files/Cert_Collisions2023_366442_370790_Golden.json");
//for prompt data 2024B - UPDATE THIS REGULARLY
  if (TString(ds.c_str()).Contains("2024"))
    //LoadJSON("files/Collisions24_13p6TeV_378981_379355_DCSOnly_TkPx.json");
    //LoadJSON("files/Collisions24_13p6TeV_378981_379774_DCSOnly_TkPx.json"); //daily json from 22.4. -> used with w12
    //LoadJSON("files/Collisions24_13p6TeV_378981_380074_DCSOnly_TkPx.json"); //daily json from 28.4. (and 29.4. same) -> used also with w12
    //LoadJSON("files/Collisions24_13p6TeV_378981_380403_DCSOnly_TkPx.json"); //daily json from 06.05. around noon --> w15
    //LoadJSON("files/Collisions24_13p6TeV_378981_380627_DCSOnly_TkPx.json"); //daily json from 12.05. --> w17
    //LoadJSON("files/Collisions24_13p6TeV_378981_380649_DCSOnly_TkPx.json"); //daily json from 15.05. --> w19
    //LoadJSON("files/Collisions24_13p6TeV_378981_380883_DCSOnly_TkPx.json"); //daily json from 19.05. --> w21
    //LoadJSON("files/Collisions24_13p6TeV_378981_381053_DCSOnly_TkPx.json"); //daily json from 22.05. --> w23
    //LoadJSON("files/Collisions24_13p6TeV_378981_381199_DCSOnly_TkPx.json"); //daily json from 26.05. --> not used (went with newer one)
    //LoadJSON("files/Collisions24_13p6TeV_378981_381212_DCSOnly_TkPx.json"); //daily json from 27.05. --> w25
		//LoadJSON("files/Collisions24_13p6TeV_378981_381478_DCSOnly_TkPx.json"); //daily json from 03.06. --> w27
		//LoadJSON("files/Collisions24_13p6TeV_378981_381516_DCSOnly_TkPx.json"); //daily json from 04.06.?? --> used for Ev2 with w27


    //LoadJSON("files/Cert_Collisions2024_378981_379075_Golden.json"); //preliminary golden json (only until B?)
    //LoadJSON("files/Cert_Collisions2024_378981_379470_Golden.json"); //golden json from 30.4.
    //LoadJSON("files/Cert_Collisions2024_eraB_Golden.json");             //golden json for 2024B only from 03.05.
    //LoadJSON("files/Cert_Collisions2024_378981_379866_Golden.json");    //golden json from 03.05.
    //LoadJSON("files/Cert_Collisions2024_378981_379866_Golden.json");    //golden json from 06.05. --> w16, w18
    //LoadJSON("files/Cert_Collisions2024_378981_380115_Golden.json");    //golden json from 15.05. --> w20
    //LoadJSON("files/Cert_Collisions2024_378981_380470_Golden.json");  //golden json from 16.05. --> w22, w24
    //LoadJSON("files/Cert_Collisions2024_378981_380649_Golden.json");  //golden json from 27.05. --> w26, w28


    //bool golden=0; // --> not used anymore, but was IMPORTANT SWITCH, could also try to check that last run# in json name and in lumi name match, as long as using my naming.
	
		//hybrid JSON! from w29 onwards
		//LoadJSON("files/CombinedJSONS_GoldenRuns_378985to381152_DCSRuns_381153to381594_.json"); //hybrid json --> w29
		//LoadJSON("files/CombinedJSONS_GoldenRuns_378981to382329_DCSRuns_382343to378981_382330to382749_.json"); //hybrid json --> w31 (05.07.2024)
		//LoadJSON("files/CombinedJSONS_GoldenRuns_378985to383448_DCSRuns_378981to378985_383449to384128_.json"); //hybrid json --> w32 (07.08.2024) and w33 (11.08.2024)
		LoadJSON("files/CombinedJSONS_GoldenRuns_378981to383724_DCSRuns_383740to378981_383725to384491_.json"); //hybrid json --> w34 (16.08.2024)




  //Cert_Collisions2023_370354_370790_Golden.json");

  // Load pileup JSON
  parsePileUpJSON("files/pileup_ASCII_UL16-UL18.txt");

  // Load pileup profiles
  LoadPU();


  //Get recorded luminosity for different triggers, pb=in picobarn:
	//w32 also used for w33
  LumiMap lumi30, lumi50, lumi110, lumi200;
	lumi30 = LoadLumi("files/lumi2024_hybrid_photon30eb_pb_w34.csv"); //w34
	lumi50 = LoadLumi("files/lumi2024_hybrid_photon50eb_pb_w34.csv");
	lumi110 = LoadLumi("files/lumi2024_hybrid_photon110eb_pb_w34.csv");
	lumi200 = LoadLumi("files/lumi2024_hybrid_photon200_pb_w34.csv");



	/*
  LumiMap lumi30, lumi50, lumi110, lumi200;
  if(golden){
      lumi30 = LoadLumi("files/lumi2024_378981_380649_golden_photon30eb_pb.csv"); //could maybe add an assertion checking whather name matches with json used
      lumi50 = LoadLumi("files/lumi2024_378981_380649_golden_photon50eb_pb.csvv");
      lumi110 = LoadLumi("files/lumi2024_378981_380649_golden_photon110eb_pb.csvv");
      lumi200 = LoadLumi("files/lumi2024_378981_380649_golden_photon200_pb.csv");
  }
  else{
      lumi30 = LoadLumi("files/lumi2024_378981_381478_daily_photon30eb_pb.csv");
      lumi50 = LoadLumi("files/lumi2024_378981_381478_daily_photon50eb_pb.csv");
      lumi110 = LoadLumi("files/lumi2024_378981_381478_daily_photon110eb_pb.csv");
      lumi200 = LoadLumi("files/lumi2024_378981_381478_daily_photon200_pb.csv");
  }
	*/


  //cout << "Use lumi files: " << endl << flush;
  //PrintInfo(string("Use lumi files") + json + ":",true);


  
  // Load veto maps
  // JECDatabase/jet_veto_maps/Summer19UL16_V0/hotjets-UL16.root
  // JECDatabase/jet_veto_maps/Summer19UL17_V2/hotjets-UL17_v2.root
  // JECDatabase/jet_veto_maps/Summer19UL18_V1/hotjets-UL18.root
  TFile *fjv(0);
  if (TString(ds.c_str()).Contains("2016"))
    fjv = new TFile("files/hotjets-UL16.root","READ");
  if (TString(ds.c_str()).Contains("2017"))
    fjv = new TFile("files/hotjets-UL17_v2.root","READ");
  if (TString(ds.c_str()).Contains("2018"))
    fjv = new TFile("files/hotjets-UL18.root","READ");
  if (TString(ds.c_str()).Contains("2022")) {
    if (TString(ds.c_str()).Contains("2022C") || 
	TString(ds.c_str()).Contains("2022D") ||
	TString(ds.c_str()).Contains("2022P8") ||
      	TString(ds.c_str()).Contains("2022QCD"))
      fjv = new TFile("files/jetveto2022CD.root","READ");
    if (TString(ds.c_str()).Contains("2022E") || // incl. EEP8 
	TString(ds.c_str()).Contains("2022F") || 
	TString(ds.c_str()).Contains("2022G") ||
      	TString(ds.c_str()).Contains("2022EEP8") ||
	TString(ds.c_str()).Contains("2022EEQCD") )
      fjv = new TFile("files/jetveto2022EFG.root","READ");
  }
  if (TString(ds.c_str()).Contains("2023")) {
    if (TString(ds.c_str()).Contains("2023B") || 
	TString(ds.c_str()).Contains("2023C") ||
	(TString(ds.c_str()).Contains("2023P8") && TString(ds.c_str()).Contains("BPix")==false) ||  //need to make sure that it does not use this for BPix stuff
	(TString(ds.c_str()).Contains("2023QCD") && TString(ds.c_str()).Contains("BPix")==false) )
      fjv = new TFile("files/jetveto2023BC.root","READ");
 // ADD MONTE CARLO SETS HERE - ALSO NEED JETVETOMAPS	(below: BPix, above: no bpix)
    if (TString(ds.c_str()).Contains("2023D") ||
	TString(ds.c_str()).Contains("2023P8-BPix") || //overwrites the previous choice of fjv (in case of BPix it first sets the wrong one, as both strings contain 22023P8
	TString(ds.c_str()).Contains("2023QCD-BPix"))
      fjv = new TFile("files/jetveto2023D.root","READ");
    }
//for now also use jetvetomap 2023D for the new 2024B prompt reco data:
  if (TString(ds.c_str()).Contains("2024")) {
    if (TString(ds.c_str()).Contains("2024B-PromptReco-v1") ||
        TString(ds.c_str()).Contains("2024B") ||
        TString(ds.c_str()).Contains("2024C") ||
        TString(ds.c_str()).Contains("2024D") ||
        TString(ds.c_str()).Contains("2024Ev1") ||
        TString(ds.c_str()).Contains("2024Ev2") ||
        TString(ds.c_str()).Contains("2024B-ECALRATIO") ||
        TString(ds.c_str()).Contains("2024C-ECALRATIO") ||
        TString(ds.c_str()).Contains("2024C-ECALR-HCALDI") ||
        //TString(ds.c_str()).Contains("winter2024P8") || //also for MC now 2024.
        //TString(ds.c_str()).Contains("2024QCD")) //also for MC now 2024.
      //fjv = new TFile("files/jetveto2024BC_V1M.root","READ"); //updated this last on 06.05.
      //fjv = new TFile("files/jetveto2024BC_V2M.root","READ"); //updated this last on 10.05. (for w17, w18 and onwards)
      //fjv = new TFile("files/jetveto2024BCD_V3M.root","READ"); //updated this last on 03.06. (for w27, w28 and onwards)
				fjv = new TFile("files/jetveto2024BCDE.root","READ"); //V5M: updated this last on 16.08. (for w34 and onwards)
    if (TString(ds.c_str()).Contains("2024F") ||
        TString(ds.c_str()).Contains("2024G"))
        TString(ds.c_str()).Contains("winter2024P8") || //also for MC now 2024.
        TString(ds.c_str()).Contains("2024QCD")) //also for MC now 2024.
        fjv = new TFile("files/jetveto2024F.root","READ"); //V5M: updated this last on 16.08. (for w34 and onwards)
  }
  if (!fjv) cout << "Jetvetomap file not found for " << ds << endl << flush;
  assert(fjv);
  
  // Veto lists for different years (NB: extra MC map for UL16):
  // h2hot_ul16_plus_hbm2_hbp12_qie11 + h2hot_mc (for UL16)
  // h2hot_ul17_plus_hep17_plus_hbpw89 (UL17)
  // h2hot_ul18_plus_hem1516_and_hbp2m1 (UL18)
  // bpix added beginning 2023D, now also applying to 2024 data
  TH2D *h2jv = 0;
  TH2D *bpixjv = 0;
  if (TString(ds.c_str()).Contains("2016")) {
    h2jv = (TH2D*)fjv->Get("h2hot_ul16_plus_hbm2_hbp12_qie11");
    assert(h2jv);
    TH2D *h2mc = (TH2D*)fjv->Get("h2hot_mc");
    assert(h2mc);
    h2jv->Add(h2mc);
  }
  if (TString(ds.c_str()).Contains("2017"))
    h2jv = (TH2D*)fjv->Get("h2hot_ul17_plus_hep17_plus_hbpw89");
  if (TString(ds.c_str()).Contains("2018"))
    h2jv = (TH2D*)fjv->Get("h2hot_ul18_plus_hem1516_and_hbp2m1");
  if (TString(ds.c_str()).Contains("2022") ||
      TString(ds.c_str()).Contains("2023") ||
      TString(ds.c_str()).Contains("2024"))
    h2jv = (TH2D*)fjv->Get("jetvetomap");
  if (TString(ds.c_str()).Contains("2024"))
    bpixjv = (TH2D*)fjv->Get("jetvetomap_bpix"); //loading the bpix vetomap for all '24 stuff
  if (!h2jv) cout << "Jetvetomap histo not found for " << ds << endl << flush;
  assert(h2jv);
  if (!bpixjv && TString(ds.c_str()).Contains("2024")){ //need extra bpix veto map only for 2024 data (for 2023 handled differntly via 2023D, see above)
		cout << "Jetvetomap for bpix not found for " << ds << endl << flush; //
  	assert(bpixjv);
	}

  // Setup B and C tagging thresholds according to Z+jet settings (Sami)
  double bthr(0.7527), cthr(0.3985), frac(0.5);
  if (TString(ds.c_str()).Contains("2016")) {
    //btagDeepB->set("2016",0.2217,0.6321,0.8953);
    //btagDeepC->set("2016",-0.48,-0.1,-0.1+frac*(0.69+0.1));
    //bthr = 0.8953; // tight
    bthr = 0.6321; // medium
    cthr = -0.1+frac*(0.69+0.1);
  }
  if (TString(ds.c_str()).Contains("2017")) {
    //btagDeepB->set("2017",0.1522,0.4941,0.8001);
    //btagDeepC->set("2017",0.05,0.15,0.15+frac*(0.8-0.15));
    //bthr = 0.8001; // tight
    bthr = 0.4941; // medium
    cthr = 0.15+frac*(0.8-0.15);
  }
  if (TString(ds.c_str()).Contains("2018")) {
    //btagDeepB->set("2018",0.1241,0.4184,0.7527);
    //btagDeepC->set("2018",0.04,0.137,0.137+frac*(0.66-0.137));
    //bthr = 0.7527; // tight
    bthr = 0.4184; // medium
    cthr = 0.137+frac*(0.66-0.137);
  }
  if (TString(ds.c_str()).Contains("22")  ||
      TString(ds.c_str()).Contains("23")) {
    // Copy of 2018
    bthr = 0.4184;
    cthr = 0.137+frac*(0.66-0.137);
  }

  // Create histograms. Copy format from existing files from Lyon
  // Keep only histograms actually used by global fit (reprocess.C)
  TDirectory *curdir = gDirectory;
  TFile *fout = new TFile(Form("rootfiles/GamHistosFill_%s_%s_%s.root",
			       isMC ? "mc" : "data",
			       dataset.c_str(), version.c_str()),
			  "RECREATE");
  assert(fout && !fout->IsZombie());
  
  // Original gamma+jet binning
  //     old bin trigger edges  (20,30,60,85,*95*,105,130,230)
  //double vx[] = {15, 20, 25, 30, 35, 40, 50, 60, 70, 85, 105, 130, 175, 230,
  //		 300, 400, 500, 600, 700, 850, 1000, 1200, 1450, 1750};
  // 22-23 binning
  double vx[] = {15, 20, 25, 30, 35, 40, 50, 60, 75, 90, 110, 130, 175, 230,
  		 300, 400, 500, 600, 700, 850, 1000, 1200, 1450, 1750,
		 2100, 2500, 3000};
  /*
  //alternative binning used in w10 to combine upper bins for more statistics --> only for these histograms. --> starting version w10 (removed in version w12)
  double vx[] = {15, 20, 25, 30, 35, 40, 50, 60, 75, 90, 110, 130, 175, 230,
  		 300, 400, 500, 600, 700, 850, 1000, 1200, 1450, 3000}; //no overflow bin?
  */

  const int nx = sizeof(vx)/sizeof(vx[0])-1;
  
  // L2L3Res eta binning
  double vy[] = {-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.650,
		 -2.500, -2.322, -2.172, -1.930, -1.653, -1.479, -1.305,
		 -1.044, -0.783, -0.522, -0.261, 0.000, 0.261, 0.522, 0.783,
		 1.044, 1.305, 1.479, 1.653, 1.930, 2.172, 2.322, 2.500,
		 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};


	//TO DO: etabinning (barrel only) also to neg. values, used for 2Dprofiles checking gain 1/6/12
  double veta[] = {-1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.879, 
					-0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 
					0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 
					0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653};

  const int ny = sizeof(vy)/sizeof(vy[0])-1;
  const int nveta = sizeof(veta)/sizeof(veta[0])-1;


  string dir = (isMC ? "MC" : "DATA");
  
  vector<pair<double,double> > etas;
  etas.push_back(make_pair<double,double>(0,1.305));
  etas.push_back(make_pair<double,double>(0,2.500));
  
  vector<double> alphas;
  alphas.push_back(1.00); 
  alphas.push_back(0.30);
  alphas.push_back(0.20); 
  alphas.push_back(0.15);
  alphas.push_back(0.10);
  alphas.push_back(0.50); 

  // Setup HT bin weighting and monitoring
  TH1D *hxsec(0), *hnevt(0), *hsumw(0), *hLHE_HT(0), *hHT(0);
  //double vht_gam[] = {0, 25, 50, 100, 200, 300, 500, 700, 1000, 1500, 2000, 6500};
  double vht_gam[] = {0, 40, 70, 100, 200, 400, 600, 6500}; // G-4Jets
  const int nht_gam = sizeof(vht_gam)/sizeof(vht_gam[0])-1;
  int nMG_gam(0);
  double wMG_gam(0);
  if (isMG && !isQCD) {

    hxsec = new TH1D("hxsec",";H_{T} (GeV);pb",nht_gam,vht_gam);
    hnevt = new TH1D("hnevt",";H_{T} (GeV);N_{evt}",nht_gam,vht_gam);
    hsumw = new TH1D("hsumw",";H_{T} (GeV);Sum(weights)",nht_gam,vht_gam);
    hLHE_HT = new TH1D("hLHE_HT",";H_{T} (GeV);N_{evt} (unweighted)",
		       nht_gam,vht_gam);
    hHT = new TH1D("hHT",";H_{T} (GeV);N_{evt} (weighted)",2485,15,2500);

    // Reference number of events, retrieved manuallay with
    // TChain c("Events"); c.AddFile("<path to files>/*.root"); c.GetEntries();
    // Also re-calculated this code before event loop when needed
    //int vnevt[nht] = {0, 0, 11197186, 23002929, 17512439, 16405924, 14359110,
    //		      13473185, 4365993, 2944561, 1836165};
    int vnevt[7] = {1, 6862, 7213, 5825, 6575, 3185, 2815}; // 2022EEP8 (local)
    double vsumw[7] = {0, 5.277e+08, 3.126e+08, 2.698e+08, 7.937e+07, 4.976e+06, 1.596e+06}; // 2022EEP8 (local)
    for (int i = 0; i != nht_gam; ++i) {
      hnevt->SetBinContent(i+1, vnevt[i]);
      nMG_gam += vnevt[i];
      hsumw->SetBinContent(i+1, vsumw[i]);
      wMG_gam += vsumw[i];
    }
    cout << "Loaded (local) MadGraph event numbers ("
	 << nMG_gam << ", sumw=" << wMG_gam << ")" << endl << flush;
    
    // xsec from jetphys/settings.h_template
    //double vxsec[nht] = {0, 0, 246300000.*23700000./28060000., 23700000,
    //			 1547000, 322600, 29980, 6334, 1088, 99.11, 20.23};
    // pT^-4
    //double vxsec[nht] = {0, 1, 1.07e-1, 2.56e-2, 1.60e-3, 1.00e-4, 1.98e-5};
    // pT^-3
    //double vxsec[nht] = {0, 1, 1.86e-1, 6.40e-2, 8.00e-3, 1.00e-3, 2.97e-4};
    // pT^-2
    //double vxsec[nht] = {0, 1, 3.27e-1, 1.60e-1, 4.00e-2, 1.00e-2, 4.44e-3};
    // pT^-2 + fine tuning on hHT histogram (a*pT^b fits at boundary)
    //double k70 = 5056./3181.;
    //double k100 = 1475./831.3 * k70;
    //double k200 = 106.8/124.3 * k100;
    //double k400 = 9.339/21.11 * k200;
    //double k600 = 3.622/5.122 * k400;
    //double vxsec[nht] = {0, 1, k70*3.27e-1, k100*1.60e-1, k200*4.00e-2,
    //			 k400*1.00e-2, k600*4.44e-3};
  //double vxsec[7] = {0, 2.727, 1.417, 1.23, 0.2643, 0.02923, 0.009177}; // 2022EEP8 hand-adjusted, including vxsec[i] *= 4./1.467; for above

    // Values from Laurent Thomas, 20 Sep 2023, based on
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HowToGenXSecAnalyzer#Running_the_GenXSecAnalyzer_on_a
    double vxsec[nht_gam] = {0, 1.524e+04, 8.111e+03, 7.327e+03, 1.541e+03,
			     1.676e+02, 5.439e+01}; // xsec in pb
    cout << Form("double vxsec[%d] = {",nht_gam);
    for (int i = 0; i != nht_gam; ++i) {
      //vxsec[i] *= 1.7893701e-4; // match 4to70 to 2.727 for 2022C test file
      hxsec->SetBinContent(i+1, vxsec[i]);
      cout << Form("%s%1.4g",i==0 ? "" : ", ", vxsec[i]);
    }
    cout << "}; // 2022EEP8 hand-adjusted" << endl << flush;
  }

  // Setup HT bin weighting and monitoring for QCD
  double vht_qcd2[] =
    {0, 25, 50, 100, 200, 300, 500, 700, 1000, 1500, 2000, 6500};
  const int nht_qcd2 = sizeof(vht_qcd2)/sizeof(vht_qcd2[0])-1;
  double vht_qcd3[] =
    {0, 40, 70, 100, 200, 400, 600, 800, 1000, 1200, 1500, 2000, 6800};
  const int nht_qcd3 = sizeof(vht_qcd3)/sizeof(vht_qcd3[0])-1;
  const double *vht_qcd = (isRun3 ? &vht_qcd3[0] : &vht_qcd2[0]);
  const int nht_qcd = (isRun3 ? nht_qcd3 : nht_qcd2);
  int nMG_qcd(0);
  double wMG_qcd(0);
  if (isMG && isQCD) {
     
     hxsec = new TH1D("hxsec",";H_{T} (GeV);pb",nht_qcd,vht_qcd);
     hnevt = new TH1D("hnevt",";H_{T} (GeV);N_{evt}",nht_qcd,vht_qcd);
     hsumw = new TH1D("hsumw",";H_{T} (GeV);Sum of weights",nht_qcd,vht_qcd);
     hLHE_HT = new TH1D("hLHE_HT",";H_{T} (GeV);N_{evt} (unweighted)",
			nht_qcd,vht_qcd);
     //hLHE_HTw = new TH1D("hLHE_HTw",";H_{T} (GeV);N_{evt} (weighted)",
     //			 nht_qcd,vht_qcd);
     hHT = new TH1D("hHT",";H_{T} (GeV);N_{evt} (weighted)",2485,15,2500);

     // Reference number of events, retrieved manuallay with
     // TChain c("Events"); c.AddFile("<path to files>/*.root"); c.GetEntries();
     // Also re-calculated this code before event loop when needed
     int vnevt2[nht_qcd2] =
       {0, 0, 11197186, 23002929, 17512439, 16405924, 14359110,
	13473185, 4365993, 2944561, 1836165};
     //int vnevt3[nht3] = {0, 9929, 26573, 16411, 10495, 8260, 7929, 10082, 14390, 6548, 6250}; // Summer22MG, local files
     int vnevt3[nht_qcd3] =
       {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // Summer22MG, local files
     int vnwgt3[nht_qcd3] =
       {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // Summer22MG, local files
     const int *vnevt = (isRun3 ? &vnevt3[0] : &vnevt2[0]);
     const int *vnwgt = (isRun3 ? &vnwgt3[0] : &vnevt2[0]);
     for (int i = 0; i != nht_qcd; ++i) {
       hnevt->SetBinContent(i+1, vnevt[i]);
       hsumw->SetBinContent(i+1, vnwgt[i]);
       nMG_qcd += vnevt[i];
       wMG_qcd += vnwgt[i];
     }
     cout << "Loaded Hefaistos MadGraph event numbers ("
	  << nMG_qcd << ")" << endl << flush;
     
     // xsec from jetphys/settings.h_template
     double vxsec2[nht_qcd2] =
       {0, 0, 246300000.*23700000./28060000., 23700000,
	1547000, 322600, 29980, 6334, 1088, 99.11, 20.23};
     // Run3 xsec from Mikel Mendizabal, MatterMost 16 Oct 2023
     double vxsec3[nht_qcd3] =
       {0, 3.136e+08,
	5.883e+07, //5.807e+07 * 3478241.4 / 5305581.5, // HT 70to100
	2.520e+07, 1.936e+06, 9.728e+04,
	1.323e+04, //3.044e+04, //HT 60to800
	3.027e+03, 8.883e+02, 3.834e+02, 1.253e+02, 2.629e+01};
     const double *vxsec = (isRun3 ? &vxsec3[0] : &vxsec2[0]);
     for (int i = 0; i != nht_qcd; ++i) {
       hxsec->SetBinContent(i+1, vxsec[i]);
     }
   } // isMG && isQCD

  const double *vht = (isMG ? (isQCD ? &vht_qcd[0] : &vht_gam[0]) : 0);
  const int nht = (isMG ? (isQCD ? nht_qcd : nht_gam) : 0);
  int nMG = (isMG ? (isQCD ? nMG_qcd : nMG_gam) : 0);
  const int wMG = (isMG ? (isQCD ? wMG_qcd : wMG_gam) : 0);
  
  // PF composition plots stored in a separate directory
  fout->mkdir("pf");
  fout->cd("pf");  
  
  TH2D *h2pteta = new TH2D("h2pteta","",nx,vx,ny,vy);
  TProfile *pabseta = new TProfile("pabseta","",nx,vx);

  // 1D composition and response
  TProfile *pdb = new TProfile("pdb","",nx,vx);
  TProfile *pmpf = new TProfile("pmpf","",nx,vx);
  TProfile *pchf = new TProfile("pchf","",nx,vx);
  TProfile *pnhf = new TProfile("pnhf","",nx,vx);
  TProfile *pnef = new TProfile("pnef","",nx,vx);
  TProfile *pcef = new TProfile("pcef","",nx,vx);
  TProfile *pmuf = new TProfile("pmuf","",nx,vx);
  TProfile *ppuf = new TProfile("ppuf","",nx,vx);

  // for (int i = 0; i != 72; ++i) cout << Form("%1.3f, ",-TMath::Pi()+i*TMath::TwoPi()/72.); cout << endl;
  const int nphi = 72;
  const double vphi[nphi+1] = 
    {-3.142, -3.054, -2.967, -2.880, -2.793, -2.705, -2.618, -2.531, -2.443,
     -2.356, -2.269, -2.182, -2.094, -2.007, -1.920, -1.833, -1.745, -1.658,
     -1.571, -1.484, -1.396, -1.309, -1.222, -1.134, -1.047, -0.960, -0.873,
     -0.785, -0.698, -0.611, -0.524, -0.436, -0.349, -0.262, -0.175, -0.087,
     0.000, 0.087, 0.175, 0.262, 0.349, 0.436, 0.524, 0.611, 0.698, 0.785,
     0.873, 0.960, 1.047, 1.134, 1.222, 1.309, 1.396, 1.484, 1.571, 1.658,
     1.745, 1.833, 1.920, 2.007, 2.094, 2.182, 2.269, 2.356, 2.443, 2.531,
     2.618, 2.705, 2.793, 2.880, 2.967, 3.054,3.142};

  // 2D composition and response
  TProfile2D *p2db = new TProfile2D("p2db","",ny,vy,nphi,vphi);
  TProfile2D *p2mpf = new TProfile2D("p2mpf","",ny,vy,nphi,vphi);
  TProfile2D *p2chf = new TProfile2D("p2chf","",ny,vy,nphi,vphi);
  TProfile2D *p2nhf = new TProfile2D("p2nhf","",ny,vy,nphi,vphi);
  TProfile2D *p2nef = new TProfile2D("p2nef","",ny,vy,nphi,vphi);
  TProfile2D *p2cef = new TProfile2D("p2cef","",ny,vy,nphi,vphi);
  TProfile2D *p2muf = new TProfile2D("p2muf","",ny,vy,nphi,vphi);
  TProfile2D *p2puf = new TProfile2D("p2puf","",ny,vy,nphi,vphi);

  // Run control plots stored in a separate directory
  fout->mkdir("runs");
  fout->cd("runs");

  // Time stability of xsec
  //TH1D *pr30n = new TH1D("pr30n",";Run;N_{events};",16000,355000.5,371000.5);         //can remove this (keep updated version)
  //TH1D *pr110n = new TH1D("pr110n",";Run;N_{events};",16000,355000.5,371000.5);
  //TH1D *pr230n = new TH1D("pr230n",";Run;N_{events};",16000,355000.5,371000.5);
  //TH1D *prg1n = new TH1D("prg1n",";Run;N_{events};",16000,355000.5,371000.5);

	//double xmax = 383000.5; //need to update this regularly
	double xmax = 386000.5;
	double xmin = 355000.5;
	double histnx = xmax-xmin; //should be int of course
  //TH1D *pr30n = new TH1D("pr30n",";Run;N_{events};",26000,355000.5,383000.5); //updated all to xmin and xmax and number of bins
	TH1D *pr30n = new TH1D("pr30n",";Run;N_{events};",histnx,xmin,xmax);
  TH1D *pr50n = new TH1D("pr50n",";Run;N_{events};",histnx,xmin,xmax);
  TH1D *pr110n = new TH1D("pr110n",";Run;N_{events};",histnx,xmin,xmax);
  TH1D *pr230n = new TH1D("pr230n",";Run;N_{events};",histnx,xmin,xmax);
  TH1D *prg1n = new TH1D("prg1n",";Run;N_{events};",histnx,xmin,xmax);

  //time stability of xs, with actual N/lumi (added 17.5.24 for monitoring of new '24 data)
  TH1D *pr30xs = new TH1D("pr30xs",";Run;xs (pb);",histnx,xmin,xmax);
  TH1D *pr50xs = new TH1D("pr50xs",";Run;xs (pb);",histnx,xmin,xmax);
  TH1D *pr110xs = new TH1D("pr110xs",";Run;xs (pb);",histnx,xmin,xmax);
  TH1D *pr230xs = new TH1D("pr230xs",";Run;xs (pb);",histnx,xmin,xmax);
 
  
  // Time stability of JEC
  TProfile *pr30b = new TProfile("pr30b",";Run;BAL;",histnx,xmin,xmax);
  TProfile *pr50b = new TProfile("pr50b",";Run;BAL;",histnx,xmin,xmax);
  TProfile *pr110b = new TProfile("pr110b",";Run;BAL;",histnx,xmin,xmax);
  TProfile *pr230b = new TProfile("pr230b",";Run;BAL;",histnx,xmin,xmax);
  TProfile *prg1b = new TProfile("prg1b",";Run;BAL;",histnx,xmin,xmax);
  TProfile *pr30m = new TProfile("pr30m",";Run;MPF;",histnx,xmin,xmax);
  TProfile *pr50m = new TProfile("pr50m",";Run;MPF;",histnx,xmin,xmax);
  TProfile *pr110m = new TProfile("pr110m",";Run;MPF;",histnx,xmin,xmax);
  TProfile *pr230m = new TProfile("pr230m",";Run;MPF;",histnx,xmin,xmax);
  TProfile *prg1m = new TProfile("prg1m",";Run;MPF;",histnx,xmin,xmax);

  // Time stability of PF composition
  TProfile *pr30chf = new TProfile("pr30chf",";Run;CHF;",histnx,xmin,xmax);
  TProfile *pr50chf = new TProfile("pr50chf",";Run;CHF;",histnx,xmin,xmax);
  TProfile *pr110chf = new TProfile("pr110chf",";Run;CHF;",histnx,xmin,xmax);
  TProfile *pr230chf = new TProfile("pr230chf",";Run;CHF;",histnx,xmin,xmax);
  TProfile *prg1chf = new TProfile("prg1chf",";Run;CHF;",histnx,xmin,xmax);
  //
  TProfile *pr30nhf = new TProfile("pr30nhf",";Run;NHF;",histnx,xmin,xmax);
  TProfile *pr50nhf = new TProfile("pr50nhf",";Run;NHF;",histnx,xmin,xmax);
  TProfile *pr110nhf = new TProfile("pr110nhf",";Run;NHF;",histnx,xmin,xmax);
  TProfile *pr230nhf = new TProfile("pr230nhf",";Run;NHF;",histnx,xmin,xmax);
  TProfile *prg1nhf = new TProfile("prg1nhf",";Run;NHF;",histnx,xmin,xmax);
  //
  TProfile *pr30nef = new TProfile("pr30nef",";Run;NHF;",histnx,xmin,xmax);
  TProfile *pr50nef = new TProfile("pr50nef",";Run;NHF;",histnx,xmin,xmax);
  TProfile *pr110nef = new TProfile("pr110nef",";Run;NHF;",histnx,xmin,xmax);
  TProfile *pr230nef = new TProfile("pr230nef",";Run;NHF;",histnx,xmin,xmax);
  TProfile *prg1nef = new TProfile("prg1nef",";Run;NHF;",histnx,xmin,xmax);

  // - - - - - - - high jet eta investigations - - - - - - - - //
  //histograms for investigation of high jet eta (couldn't these be added as a loop???)
  //like: TH1D *pr30n, *pr50n, *pr110n
  //for(histi=0; histi<histos.size(); histi++){<create histos>}
  fout->mkdir("runs_high-jeteta");
  fout->cd("runs_high-jeteta");

  // Time stability of xsec
  TH1D *pr30n_eta3to4 = new TH1D("pr30n_eta3to4","3 <= #eta_jet{} < 4;Run;N_{events};",histnx,xmin,xmax);
  TH1D *pr50n_eta3to4 = new TH1D("pr50n_eta3to4",";Run;N_{events};",histnx,xmin,xmax);
  TH1D *pr110n_eta3to4 = new TH1D("pr110n_eta3to4",";Run;N_{events};",histnx,xmin,xmax);

  TH1D *pr30n_eta4to5 = new TH1D("pr30n_eta4to5","4 <= jet_eta < 5;Run;N_{events};",histnx,xmin,xmax);
  TH1D *pr50n_eta4to5 = new TH1D("pr50n_eta4to5",";Run;N_{events};",histnx,xmin,xmax);
  TH1D *pr110n_eta4to5 = new TH1D("pr110n_eta4to5",";Run;N_{events};",histnx,xmin,xmax);


  //time stability of xs, with actual N/lumi (added 17.5.24 for monitoring of new '24 data)
  TH1D *pr30xs_eta3to4 = new TH1D("pr30xs_eta3to4",";Run;xs (pb);",histnx,xmin,xmax);
  TH1D *pr50xs_eta3to4 = new TH1D("pr50xs_eta3to4",";Run;xs (pb);",histnx,xmin,xmax);
  TH1D *pr110xs_eta3to4 = new TH1D("pr110xs_eta3to4",";Run;xs (pb);",histnx,xmin,xmax);

  TH1D *pr30xs_eta4to5 = new TH1D("pr30xs_eta4to5",";Run;xs (pb);",histnx,xmin,xmax);
  TH1D *pr50xs_eta4to5 = new TH1D("pr50xs_eta4to5",";Run;xs (pb);",histnx,xmin,xmax);
  TH1D *pr110xs_eta4to5 = new TH1D("pr110xs_eta4to5",";Run;xs (pb);",histnx,xmin,xmax);
 
  
  // Time stability of JEC
  TProfile *pr30b_eta3to4 = new TProfile("pr30b_eta3to4",";Run;BAL;",histnx,xmin,xmax);
  TProfile *pr50b_eta3to4 = new TProfile("pr50b_eta3to4",";Run;BAL;",histnx,xmin,xmax);
  TProfile *pr110b_eta3to4 = new TProfile("pr110b_eta3to4",";Run;BAL;",histnx,xmin,xmax);
  TProfile *pr30m_eta3to4 = new TProfile("pr30m_eta3to4",";Run;MPF;",histnx,xmin,xmax);
  TProfile *pr50m_eta3to4 = new TProfile("pr50m_eta3to4",";Run;MPF;",histnx,xmin,xmax);
  TProfile *pr110m_eta3to4 = new TProfile("pr110m_eta3to4",";Run;MPF;",histnx,xmin,xmax);

  TProfile *pr30b_eta4to5 = new TProfile("pr30b_eta4to5",";Run;BAL;",histnx,xmin,xmax);
  TProfile *pr50b_eta4to5 = new TProfile("pr50b_eta4to5",";Run;BAL;",histnx,xmin,xmax);
  TProfile *pr110b_eta4to5 = new TProfile("pr110b_eta4to5",";Run;BAL;",histnx,xmin,xmax);
  TProfile *pr30m_eta4to5 = new TProfile("pr30m_eta4to5",";Run;MPF;",histnx,xmin,xmax);
  TProfile *pr50m_eta4to5 = new TProfile("pr50m_eta4to5",";Run;MPF;",histnx,xmin,xmax);
  TProfile *pr110m_eta4to5 = new TProfile("pr110m_eta4to5",";Run;MPF;",histnx,xmin,xmax);

  //Time stability of h vs em in HF
  TProfile *pr30hfEmEF_eta3to4 = new TProfile("pr30hfEmEF_eta3to4","HF em energy fraction (30EB);Run;N_{events};",histnx,xmin,xmax);
  TProfile *pr50hfEmEF_eta3to4 = new TProfile("pr50hfEmEF_eta3to4","HF em energy fraction (50EB);Run;N_{events};",histnx,xmin,xmax);
  TProfile *pr110hfEmEF_eta3to4 = new TProfile("pr110hfEmEF_eta3to4","HF em energy fraction (110EB);Run;N_{events};",histnx,xmin,xmax);
  TProfile *pr30hfHEF_eta3to4 = new TProfile("pr30hfHEF_eta3to4","HF had energy fraction (30EB);Run;N_{events};",histnx,xmin,xmax);
  TProfile *pr50hfHEF_eta3to4 = new TProfile("pr50hfHEF_eta3to4","HF had energy fraction (50EB);Run;N_{events};",histnx,xmin,xmax);
  TProfile *pr110hfHEF_eta3to4 = new TProfile("pr110hfHEF_eta3to4","HF had energy fraction (110EB);Run;N_{events};",histnx,xmin,xmax);

  TProfile *pr30hfEmEF_eta4to5 = new TProfile("pr30hfEmEF_eta4to5","HF em energy fraction (30EB);Run;N_{events};",histnx,xmin,xmax);
  TProfile *pr50hfEmEF_eta4to5 = new TProfile("pr50hfEmEF_eta4to5","HF em energy fraction (50EB);Run;N_{events};",histnx,xmin,xmax);
  TProfile *pr110hfEmEF_eta4to5 = new TProfile("pr110hfEmEF_eta4to5","HF em energy fraction (110EB);Run;N_{events};",histnx,xmin,xmax);
  TProfile *pr30hfHEF_eta4to5 = new TProfile("pr30hfHEF_eta4to5","HF had energy fraction (30EB);Run;N_{events};",histnx,xmin,xmax);
  TProfile *pr50hfHEF_eta4to5 = new TProfile("pr50hfHEF_eta4to5","HF had energy fraction (50EB);Run;N_{events};",histnx,xmin,xmax);
  TProfile *pr110hfHEF_eta4to5 = new TProfile("pr110hfHEF_eta4to5","HF had energy fraction (110EB);Run;N_{events};",histnx,xmin,xmax);




  // Time stability of PF composition -- not yet for high jet eta... do this later if needed
  /*
  TProfile *pr30chf_eta3to4 = new TProfile("pr30chf",";Run;CHF;",26000,355000.5,383000.5);
  TProfile *pr50chf_eta3to4 = new TProfile("pr50chf",";Run;CHF;",26000,355000.5,383000.5);
  TProfile *pr110chf_eta3to4 = new TProfile("pr110chf",";Run;CHF;",26000,355000.5,383000.5);
  //
  TProfile *pr30nhf_eta3to4 = new TProfile("pr30nhf",";Run;NHF;",26000,355000.5,383000.5);
  TProfile *pr50nhf_eta3to4 = new TProfile("pr50nhf",";Run;NHF;",26000,355000.5,383000.5);
  TProfile *pr110nhf_eta3to4 = new TProfile("pr110nhf",";Run;NHF;",26000,355000.5,383000.5);
  //
  TProfile *pr30nef_eta3to4 = new TProfile("pr30nef",";Run;NHF;",26000,355000.5,383000.5);
  TProfile *pr50nef_eta3to4 = new TProfile("pr50nef",";Run;NHF;",26000,355000.5,383000.5);
  TProfile *pr110nef_eta3to4 = new TProfile("pr110nef",";Run;NHF;",26000,355000.5,383000.5);
  */







  // Control plots stored in a separate directory
  fout->mkdir("control");
  fout->cd("control");
  
  // Cut flow controls
  TProfile *pcutflow1 = new TProfile("pcutflow1",";Cut;Pass",20,-0.5,19.5);
  TProfile *pcutflow2 = new TProfile("pcutflow2",";Cut;Pass",20,-0.5,19.5);
  TProfile *pcutflow3 = new TProfile("pcutflow3",";Cut;Pass",20,-0.5,19.5);

  const int nc = 18;
  string cuts[nc] = {"pass_trig", "pass_ngam", "pass_njet", "pass_gameta",
		     //"pass_dphi", "pass_jetid", "pass_veto", "pass_leak",
                     "pass_dphi", "pass_jetid", "pass_jetveto", "pass_gamveto", "pass_leak",
		     "pass_basic", "pass_bal", "pass_mpf", "pass_basic_ext",
		     "pass_jeteta", "pass_alpha100", "pass_all", "pass_gen"}; 
  for (int i = 0; i != nc; ++i) {
    pcutflow1->GetXaxis()->SetBinLabel(i+1, cuts[i].c_str());
    pcutflow2->GetXaxis()->SetBinLabel(i+1, cuts[i].c_str());
    pcutflow3->GetXaxis()->SetBinLabel(i+1, cuts[i].c_str());
  }

  // Follow up on problematic cuts
  TH1D *hdphi = new TH1D("hdphi",";#Delta#phi;N_{events}",
			 126,-TMath::TwoPi(),TMath::TwoPi());
  TH1D *hdr = new TH1D("hdr",";#DeltaR;N_{events}",100,0,10);
		       

  

  // 1D plots for mu per trigger
  //hmusmc,20,30,50,75,90,100,110,200
  TH1D *hmus  = new TH1D("hmus","",100,0,100);
  TH1D *hmusmc = new TH1D("hmusmc","",100,0,100);
  TH1D *hmus20 = new TH1D("hmus20","",100,0,100);
  TH1D *hmus30 = new TH1D("hmus30","",100,0,100);
  TH1D *hmus50 = new TH1D("hmus50","",100,0,100);
  TH1D *hmus75 = new TH1D("hmus75","",100,0,100);
  TH1D *hmus90 = new TH1D("hmus90","",100,0,100);
  TH1D *hmus100 = new TH1D("hmus100","",100,0,100);
  TH1D *hmus110 = new TH1D("hmus110","",100,0,100);
  TH1D *hmus200 = new TH1D("hmus200","",100,0,100);

  // 2D plots for mu vs photon pT
  TH2D *h2mus = new TH2D("h2mus","",nx,vx,100,0,100);

  // Plots of npvgood, npvall vs mu
  TProfile *pmuvsmu = new TProfile("pmuvsmu","",100,0,100);
  TProfile *prhovsmu = new TProfile("prhovsmu","",100,0,100);
  TProfile *pnpvgoodvsmu = new TProfile("pnpvgoodvsmu","",100,0,100);
  TProfile *pnpvallvsmu = new TProfile("pnpvallvsmu","",100,0,100);
  // Plots of photon corr, err, hoe, r9, vs mu 
  TProfile *pgainvsmu = new TProfile("pgainvsmu","",100,0,100);
  TProfile *pcorrvsmu = new TProfile("pcorrvsmu","",100,0,100);
  TProfile *perrvsmu = new TProfile("perrvsmu","",100,0,100);
  TProfile *phoevsmu = new TProfile("phoevsmu","",100,0,100);
  TProfile *pr9vsmu = new TProfile("pr9vsmu","",100,0,100);
  // ...and vs pT
  TProfile *pmuvspt = new TProfile("pmuvspt","",nx,vx);
  TProfile *prhovspt = new TProfile("prhovspt","",nx,vx);
  TProfile *pnpvgoodvspt = new TProfile("pnpvgoodvspt","",nx,vx);
  TProfile *pnpvallvspt = new TProfile("pnpvallvspt","",nx,vx);
  // ..and vs pT
  TProfile *pgain1vspt = new TProfile("pgain1vspt","",nx,vx);
  TProfile *pgain6vspt = new TProfile("pgain6vspt","",nx,vx);
  TProfile *pgain12vspt = new TProfile("pgain12vspt","",nx,vx);
  TProfile *pgainvspt = new TProfile("pgainvspt","",nx,vx);
  TProfile *pcorrvspt = new TProfile("pcorrvspt","",nx,vx);
  TProfile *perrvspt = new TProfile("perrvspt","",nx,vx);
  TH2D *h2hoevspt = new TH2D("h2hoevspt","",nx,vx,125,0,0.025);
  TProfile *phoevspt = new TProfile("phoevspt","",nx,vx);
  TH2D *h2r9vspt = new TH2D("h2r9vspt","",nx,vx,150,0.90,1.05);
  TProfile *pr9vspt = new TProfile("pr9vspt","",nx,vx);

	//new (w27+w28): 2D plots for gain vs pt and eta (nx = #xbins, vx = pt-xbins, ny=#ybins, vy=eta-ybins)
	//changed to narrower eta-bins called veta, #bins=nveta
  TProfile2D *pgain1vsptvseta = new TProfile2D("pgain1vsptvseta","",nx,vx,nveta,veta);
  TProfile2D *pgain6vsptvseta = new TProfile2D("pgain6vsptvseta","",nx,vx,nveta,veta);
  TProfile2D *pgain12vsptvseta = new TProfile2D("pgain12vsptvseta","",nx,vx,nveta,veta);
  TProfile2D *pgainvsptvseta = new TProfile2D("pgainvsptvseta","",nx,vx,nveta,veta);

	//new (w34): 2D plots (vs eta vs phi) for jet response (for 30GeV, 50GeV and 110GeV photon trigger)
	//TH2D *h2balvsetavsphi
	TProfile2D *p2bal50_jetetaphi = new TProfile2D("p2bal50_jetetaphi", ";#eta_{j1};#phi_{j1};BAL50", nveta, veta, 72, -TMath::Pi(), TMath::Pi()); //same #phibins as h2gametaphi (72); etabins as above (veta)
	//TH2D *h2mpf50_jetetaphi = new TH2D("h2mpf50_jetetaphi", ";#eta;#phi;MPF50", nveta, veta, 72, -TMath::Pi(), TMath::Pi());
 	TProfile2D *p2bal110_jetetaphi = new TProfile2D("p2bal110_jetetaphi", ";#eta_{j1};#phi_{j1};BAL110", nveta, veta, 72, -TMath::Pi(), TMath::Pi());
	//TH2D *h2mpf110_jetetaphi = new TH2D("h2mpf110_jetetaphi", ";#eta;#phi;MPF110", nveta, veta, 72, -TMath::Pi(), TMath::Pi());
 	TProfile2D *p2bal200_jetetaphi = new TProfile2D("p2bal200_jetetaphi", ";#eta_{j1};#phi_{j1};BAL200", nveta, veta, 72, -TMath::Pi(), TMath::Pi());
	//TH2D *h2mpf200_jetetaphi = new TH2D("h2mpf110_jetetaphi", ";#eta;#phi;MPF200", nveta, veta, 72, -TMath::Pi(), TMath::Pi());

	//new (w34): 2D plots (vs jet eta vs jet phi) for event rate (for 30GeV, 50GeV and 110GeV photon trigger)
	TH2D *h2n50_jetetaphi = new TH2D("h2n50_jetetaphi", "Rate for 50GeV #gamma-trigger;#eta_{j1};#phi_{j1};N_{evt}", nveta, veta, 72, -TMath::Pi(), TMath::Pi());
 	TH2D *h2n110_jetetaphi = new TH2D("h2n110_jetetaphi", "Rate for 110GeV #gamma-trigger;#eta_{j1};#phi_{j1};N_{evt}", nveta, veta, 72, -TMath::Pi(), TMath::Pi());
 	TH2D *h2n200_jetetaphi = new TH2D("h2n200_jetetaphi", "Rate for 200GeV #gamma-trigger;#eta_{j1};#phi_{j1};N_{evt}", nveta, veta, 72, -TMath::Pi(), TMath::Pi());

	//new (w34): 2D plots (vs PHOTON eta vs jet phi) for event rate (for 30GeV, 50GeV and 110GeV photon trigger)
	TH2D *h2n50_gametaphi = new TH2D("h2bal50_gagetaphi", "Rate for 50GeV #gamma-trigger;#eta_{#gamma};#phi_{#gamma};N_{evt}", nveta, veta, 72, -TMath::Pi(), TMath::Pi());
 	TH2D *h2n110_gametaphi = new TH2D("h2bal110_gametaphi", "Rate for 110GeV #gamma-trigger;#eta_{#gamma};#phi_{#gamma};N_{evt}", nveta, veta, 72, -TMath::Pi(), TMath::Pi());
 	TH2D *h2n200_gametaphi = new TH2D("h2bal200_gametaphi", "Rate for 200GeV #gamma-trigger;#eta_{#gamma};#phi_{#gamma};N_{evt}", nveta, veta, 72, -TMath::Pi(), TMath::Pi());



	//Narrow eta-binning (from dijet), but also to negative values
  double etabins[] = {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, 
											-2.853, -2.65, -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, 
											-1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 
											0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
       								1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 
       								3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};
  const int netabins = sizeof(etabins)/sizeof(etabins[0])-1;


	//new (w35): N_events over leading jet's eta for each photon trigger and for three different jet-pt-thresholds
	TH1D *h50n_j1eta_j1pt30 = new TH1D("h50n_j1eta_j1pt30","p_{T,jet1} > 30GeV;#eta_{j1};N_{events};",netabins,etabins); //j1pt > 30GeV
	TH1D *h50n_j1eta_j1pt40 = new TH1D("h50n_j1eta_j1pt40",";#eta_{j1};N_{events};",netabins,etabins); //j1pt > 40GeV
	TH1D *h50n_j1eta_j1pt50 = new TH1D("h50n_j1eta_j1pt50",";#eta_{j1};N_{events};",netabins,etabins); //j1pt > 50GeV

	TH1D *h110n_j1eta_j1pt30 = new TH1D("h110n_j1eta_j1pt30","p_{T,jet1} > 30GeV;#eta_{j1};N_{events};",netabins,etabins); //j1pt > 30GeV
	TH1D *h110n_j1eta_j1pt40 = new TH1D("h110n_j1eta_j1pt40",";#eta_{j1};N_{events};",netabins,etabins); //j1pt > 40GeV
	TH1D *h110n_j1eta_j1pt50 = new TH1D("h110n_j1eta_j1pt50",";#eta_{j1};N_{events};",netabins,etabins); //j1pt > 50GeV

	TH1D *h200n_j1eta_j1pt30 = new TH1D("h200n_j1eta_j1pt30","p_{T,jet1} > 30GeV;#eta_{j1};N_{events};",netabins,etabins); //j1pt > 30GeV
	TH1D *h200n_j1eta_j1pt40 = new TH1D("h200n_j1eta_j1pt40",";#eta_{j1};N_{events};",netabins,etabins); //j1pt > 40GeV
	TH1D *h200n_j1eta_j1pt50 = new TH1D("h200n_j1eta_j1pt50",";#eta_{j1};N_{events};",netabins,etabins); //j1pt > 50GeV






  // 2D plots for jet response
  TH2D *h2bal = new TH2D("h2bal","",nx,vx,200,0,4);
  TH2D *h2mpf = new TH2D("h2mpf","",nx,vx,300,-2,4);
  TH2D *h2balc = new TH2D("h2balc","",nx,vx,200,0,4);
  TH2D *h2mpfc = new TH2D("h2mpfc","",nx,vx,300,-2,4);
  TH2D *h2balc2 = new TH2D("h2balc2","",nx,vx,200,0,4);
  TH2D *h2mpfc2 = new TH2D("h2mpfc2","",nx,vx,300,-2,4);
  
  // 2D and profile for photon-jet pT
  TH2D *h2phoj = new TH2D("h2phoj","",nx,vx,240,-0.1,0.5);
  TProfile *pphoj = new TProfile("pphoj","",nx,vx);

  // Extras without zero suppression and for gain paths
  TH2D *h2phoj0 = new TH2D("h2phoj0","",nx,vx,140,-0.2,0.5);
  TProfile *pphoj0 = new TProfile("pphoj0","",nx,vx);
  TH2D *h2phoj1 = new TH2D("h2phoj1","",nx,vx,140,-0.2,0.5);
  TProfile *pphoj1 = new TProfile("pphoj1","",nx,vx);
  TH2D *h2phoj6 = new TH2D("h2phoj6","",nx,vx,140,-0.2,0.5);
  TProfile *pphoj6 = new TProfile("pphoj6","",nx,vx);
  TH2D *h2phoj12 = new TH2D("h2phoj12","",nx,vx,140,-0.2,0.5);
  TProfile *pphoj12 = new TProfile("pphoj12","",nx,vx);

  // Plots for photon properties (more in MC)
  TH2D *h2gametaphi = new TH2D("h2gametaphi","",30,-1.305,+1.305,
			       72,-TMath::Pi(),TMath::Pi());
  //TH2D *h2gametaphi2 = new TH2D("h2gametaphi2","",150,-1.305,+1.305,
  // Match barrel edge to 1.305 with 3.132, even though EC edge should be 3.139
  TH2D *h2gametaphi2 = new TH2D("h2gametaphi2","",360,-3.132,+3.132,
				360,-TMath::Pi(),TMath::Pi());
  //TH2D *h2gametaphi3 = new TH2D("h2gametaphi3","",150,-1.305,+1.305,
  //				720,-TMath::Pi(),TMath::Pi());
  //TH2D *h2gametaphi4 = new TH2D("h2gametaphi4","",150,-1.305,+1.305,
  //				1440,-TMath::Pi(),TMath::Pi());
  TH2D *h2ngam = new TH2D("h2ngam","",nx,vx,5,0,5);
  TH1D *hgen = new TH1D("hgen","",nx,vx);
  TH1D *hgam = new TH1D("hgam","",nx,vx);
  TH1D *hgamtrg = new TH1D("hgamtrg","",nx,vx);
  TProfile *peffgr = new TProfile("peffgr","",nx,vx);
  TProfile *peffid = new TProfile("peffid","",nx,vx);
  TProfile *pfake = new TProfile("pfake","",nx,vx);
  TProfile *pfakeqcd = new TProfile("pfakeqcd","",nx,vx); // for QCD bkg
  TProfile *pfakeqcd2 = new TProfile("pfakeqcd2","",nx,vx); // for QCD bkg
  TH2D *h2rgam = new TH2D("h2rgam","",nx,vx,350,0.80,1.15);
  TH2D *h2rgamqcd = new TH2D("h2rgamqcd","",nx,vx,350,0.80,1.15); // for QCD bkg
  TProfile *prgam = new TProfile("prgam","",nx,vx);
  TProfile *prgamqcd = new TProfile("prgamqcd","",nx,vx); // for QCD bkg
  TProfile *prgamqcd2 = new TProfile("prgamqcd2","",nx,vx); // for QCD bkg
  TProfile *prgamorigqcd = new TProfile("prgamorigqcd","",nx,vx); // fpr QCD bkg
  TProfile *prgamorigqcd2 = new TProfile("prgamorigqcd2","",nx,vx); // QCD
  TH2D *h2cgam = new TH2D("h2cgam","",nx,vx,100,0.90,1.10);
  TProfile *pcgam = new TProfile("pcgam","",nx,vx);

  // Plots for jet properties
  TH2D *h2gjet = new TH2D("h2gjet","",nx,vx,100,0.90,1.10);
  TProfile *pgjet = new TProfile("pgjet","",nx,vx);
  TH2D *h2rjet = new TH2D("h2rjet","",nx,vx,100,0.90,1.10);
  TProfile *prjet = new TProfile("prjet","",nx,vx);
  TH2D *h2rgen = new TH2D("h2rgen","",nx,vx,100,0.90,1.10);
  TProfile *prgen = new TProfile("prgen","",nx,vx);

  TH2D *h2gjet2 = new TH2D("h2gjet2","",nx,vx,100,0.90,1.10);
  TProfile *pgjet2 = new TProfile("pgjet2","",nx,vx);
  TH2D *h2rjet2 = new TH2D("h2rjet2","",nx,vx,100,0.90,1.10);
  TProfile *prjet2 = new TProfile("prjet2","",nx,vx);
  TH2D *h2rgen2 = new TH2D("h2rgen2","",nx,vx,100,0.90,1.10);
  TProfile *prgen2 = new TProfile("prgen2","",nx,vx);

  // Plots for jet flavor
  /*
  TProfile *prgenud = new TProfile("prgenud","",nx,vx);
  TProfile *prgens = new TProfile("prgens","",nx,vx);
  TProfile *prgenc = new TProfile("prgenc","",nx,vx);
  TProfile *prgenb = new TProfile("prgenb","",nx,vx);
  TProfile *prgeng = new TProfile("prgeng","",nx,vx);
  TProfile *prgeno = new TProfile("prgeno","",nx,vx);
  TProfile *pgjetud = new TProfile("pgjetud","",nx,vx);
  TProfile *pgjets = new TProfile("pgjets","",nx,vx);
  TProfile *pgjetc = new TProfile("pgjetc","",nx,vx);
  TProfile *pgjetb = new TProfile("pgjetb","",nx,vx);
  TProfile *pgjetg = new TProfile("pgjetg","",nx,vx);
  TProfile *pgjeto = new TProfile("pgjeto","",nx,vx);
  TProfile *pfud = new TProfile("pfud","",nx,vx);
  TProfile *pfs = new TProfile("pfs","",nx,vx);
  TProfile *pfc = new TProfile("pfc","",nx,vx);
  TProfile *pfb = new TProfile("pfb","",nx,vx);
  TProfile *pfg = new TProfile("pfg","",nx,vx);
  TProfile *pfo = new TProfile("pfo","",nx,vx);
  */

  // Plots for photon response in data
  TProfile *prbal = new TProfile("prbal","",nx,vx);
  TProfile *prmpf = new TProfile("prmpf","",nx,vx);
  TProfile *prbal0 = new TProfile("prbal0","",980,20,1000);
  TProfile *prmpf0 = new TProfile("prmpf0","",980,20,1000);
  TProfile *prbal1 = new TProfile("prbal1","",nx,vx);
  TProfile *prmpf1 = new TProfile("prmpf1","",nx,vx);
  TProfile *prbal6 = new TProfile("prbal6","",nx,vx);
  TProfile *prmpf6 = new TProfile("prmpf6","",nx,vx);
  TProfile *prbal12 = new TProfile("prbal12","",nx,vx);
  TProfile *prmpf12 = new TProfile("prmpf12","",nx,vx);

	//new (w27+w28): 2D plots for gain vs pt and eta (nx = #xbins, vx = pt-xbins, ny=#ybins, vy=eta-ybins)
	//changed to narrower eta-bins called veta, #bins=nveta
  TProfile2D *pr2bal = new TProfile2D("pr2bal","",nx,vx,nveta,veta);
  TProfile2D *pr2mpf = new TProfile2D("pr2mpf","",nx,vx,nveta,veta);
  TProfile2D *pr2bal1 = new TProfile2D("pr2bal1","",nx,vx,nveta,veta);
  TProfile2D *pr2mpf1 = new TProfile2D("pr2mpf1","",nx,vx,nveta,veta);
  TProfile2D *pr2bal6 = new TProfile2D("pr2bal6","",nx,vx,nveta,veta);
  TProfile2D *pr2mpf6 = new TProfile2D("pr2mpf6","",nx,vx,nveta,veta);
  TProfile2D *pr2bal12 = new TProfile2D("pr2bal12","",nx,vx,nveta,veta);
  TProfile2D *pr2mpf12 = new TProfile2D("pr2mpf12","",nx,vx,nveta,veta);


  
  // Plots for photon trigger efficiencies
  // TBD: need to create these more systematically with a loop (yep, agree... will look into this, - Bettina)
  TH1D *hgam0_data = new TH1D("hgam0_data","",197,15,1000);
  TH1D *hgam0_mc = new TH1D("hgam0_mc","",197,15,1000);
  TH1D *hgam0 = new TH1D("hgam0","",197,15,1000);
  TH1D *hgam20 = new TH1D("hgam20","",197,15,1000);
  TH1D *hgam22 = new TH1D("hgam22","",197,15,1000);
  TH1D *hgam30 = new TH1D("hgam30","",197,15,1000);
  TH1D *hgam33 = new TH1D("hgam33","",197,15,1000);
  TH1D *hgam36 = new TH1D("hgam36","",197,15,1000);
  TH1D *hgam50 = new TH1D("hgam50","",197,15,1000);
  TH1D *hgam75 = new TH1D("hgam75","",197,15,1000);
  TH1D *hgam90 = new TH1D("hgam90","",197,15,1000);
  TH1D *hgam120 = new TH1D("hgam120","",197,15,1000);
  TH1D *hgam150 = new TH1D("hgam150","",197,15,1000);
  TH1D *hgam175 = new TH1D("hgam175","",197,15,1000);
  TH1D *hgam200 = new TH1D("hgam200","",197,15,1000);
  TH1D *hgam300 = new TH1D("hgam300","",197,15,1000);
  TH1D *hgam30t = new TH1D("hgam30t","",197,15,1000);
  TH1D *hgam50t = new TH1D("hgam50t","",197,15,1000); //new (27.05.2024)
  TH1D *hgam100t = new TH1D("hgam100t","",197,15,1000);
  TH1D *hgam110t = new TH1D("hgam110t","",197,15,1000);
  TH1D *hgam120t = new TH1D("hgam120t","",197,15,1000);
  TH1D *hgam22m = new TH1D("hgam22m","",197,15,1000);
  TH1D *hgam30m = new TH1D("hgam30m","",197,15,1000);
  TH1D *hgam36m = new TH1D("hgam36m","",197,15,1000);
  TH1D *hgam50m = new TH1D("hgam50m","",197,15,1000);
  TH1D *hgam75m = new TH1D("hgam75m","",197,15,1000);
  TH1D *hgam90m = new TH1D("hgam90m","",197,15,1000);
  TH1D *hgam120m = new TH1D("hgam120m","",197,15,1000);
  TH1D *hgam165m = new TH1D("hgam165m","",197,15,1000);
  TH1D *hgam165h = new TH1D("hgam165h","",197,15,1000);
  TH1D *hgam100h = new TH1D("hgam100h","",197,15,1000);
  TH1D *hgam20l = new TH1D("hgam20l","",197,15,1000);
  TH1D *hgam30l = new TH1D("hgam30l","",197,15,1000);
  TH1D *hgam40l = new TH1D("hgam40l","",197,15,1000);
  TH1D *hgam50l = new TH1D("hgam50l","",197,15,1000);
  TH1D *hgam60l = new TH1D("hgam60l","",197,15,1000);
  TH1D *hgamtrig = new TH1D("hgamtrig","",197,15,1000);
  TH1D *hgamtrig_data = new TH1D("hgamtrig_data","",197,15,1000);
  TH1D *hgamtrig_mc = new TH1D("hgamtrig_mc","",197,15,1000);

  // Flavor plots stored in a separate directory
  fout->mkdir("flavor");
  fout->cd("flavor");

  map<string, double> mvar;
  map<string, map<string, map<string, TH1*> > > mp;
  string avar[] = {"counts","mpfchs1","ptchs","mpf1","mpfn","mpfu","rho",
		   "rjet","gjet","rgen"};
  string atag[] = {"i","b","c","q","g","n"};
  //string aflv[] = {"i","b","c","q","g","n"};
  string aflv[] = {"i","b","c","q","s","ud","g","n"};
  const int nvar = sizeof(avar)/sizeof(avar[0]);
  const int ntag = sizeof(atag)/sizeof(atag[0]);
  const int nflv = sizeof(aflv)/sizeof(aflv[0]);
  
  cout << "Creating flavor histograms/profiles in folder flavor/..." << endl;
  for (int ivar = 0; ivar != nvar; ++ivar) {
    for (int itag = 0; itag != ntag; ++itag) {
      for (int iflv = 0; iflv != nflv; ++iflv) {
	string& var = avar[ivar]; const char *cv = var.c_str();
	string& tag = atag[itag]; const char *ct = tag.c_str();
	string& flv = aflv[iflv]; const char *cf = flv.c_str();
	if (var=="counts")
	  mp[var][tag][flv] = new TH1D(Form("%s_g%s%s",cv,ct,cf),"",nx,vx);
	else
	  mp[var][tag][flv] = new TProfile(Form("%s_g%s%s",cv,ct,cf),"",nx,vx);
      } // for iflv
    } // for itag
  } // for ivar

  // Results similar to Multijet in dijet package
  gamjetHistos *hg(0);
  if (doGamjet) {

    // L3Res (gamma+jet) pT binning adapted and extended
    //const double vpt[] = {15, 20, 25, 30, 35,
    //			  40, 50, 60, 70, 85, 100, 125, 155, 180, 210, 250, 300,
    //			  350, 400, 500, 600, 800, 1000, 1200, 1500,
    //			  1800, 2100, 2400, 2700, 3000};
    // 22-23 binning
    const double vpt[] = {15, 20, 25, 30, 35,
			  40, 50, 60, 75, 90, 110, 130, 175, 230, 300,
			  400, 500, 600, 700, 850, 1000, 1200, 1450,
			  1750, 2100, 2500, 3000};

    double npt = sizeof(vpt)/sizeof(vpt[0])-1;

    if (debug) cout << "Setup doGamjet " << endl << flush;
      
    fout->mkdir("Gamjet");
    fout->cd("Gamjet");

    gamjetHistos *h = new gamjetHistos();
    hg = h;

    // Counting of events, and JEC L2L3Res+JERSF for undoing
    h->hpt13 = new TH1D("hpt13",";p_{T,#gamma};N_{events}",npt,vpt);
    h->hpt13a = new TH1D("hpt13a",";p_{T,avg};N_{events}",npt,vpt);
    h->hpt13j = new TH1D("hpt13j",";p_{T,jet};N_{events}",npt,vpt);
    
    h->pptg = new TProfile("pptg",";p_{T,#gamma};p_{T,#gamma}",npt,vpt);
    h->pptj = new TProfile("pptj",";p_{T,#gamma};p_{T,jet}",npt,vpt);

   // MPF decomposition for HDM method
    h->pres = new TProfile("pres",";p_{T,#gamma} (GeV);JES",npt,vpt);
    h->pjsf = new TProfile("pjsf",";p_{T,#gamma} (GeV);JERSF",npt,vpt);
    h->pm0 = new TProfile("pm0",";p_{T,#gamma} (GeV);MPF0",npt,vpt);
    h->pm2 = new TProfile("pm2",";p_{T,#gamma} (GeV);MPF2",npt,vpt);
    h->pmn = new TProfile("pmn",";p_{T,#gamma} (GeV);MPFn",npt,vpt);
    h->pmu = new TProfile("pmu",";p_{T,#gamma} (GeV);MPFu",npt,vpt);
    
    h->pm0x = new TProfile("pm0x",";p_{T,#gamma} (GeV);MPF0X (MPFX)", npt, vpt, "S");
    h->pm2x = new TProfile("pm2x",";;p_{T,#gamma} (GeV);MPF2X (DBX)", npt, vpt, "S");

    // Extra for FSR studies
    h->pmnu = new TProfile("pmnu",";p_{T,#gamma} (GeV);MPFnu", npt, vpt);
    h->pmnx = new TProfile("pmnx",";p_{T,#gamma} (GeV);MPFNX", npt, vpt, "S");
    h->pmux = new TProfile("pmux",";p_{T,#gamma} (GeV);MPFUX", npt, vpt, "S");
    h->pmnux = new TProfile("pmnux",";p_{T,#gamma} (GeV;MPFNUX", npt, vpt, "S");
    
    // Composition
    h->prho = new TProfile("prho",";p_{T,#gamma};#rho (GeV)",npt,vpt);
    h->pchf = new TProfile("pchf",";p_{T,#gamma};CHF",npt,vpt);
    h->pnhf = new TProfile("pnhf",";p_{T,#gamma};NHF",npt,vpt);
    h->pnef = new TProfile("pnef",";p_{T,#gamma};NEF",npt,vpt);

    // Alternative pT bins
    h->presa = new TProfile("presa",";p_{T,avp} (GeV);JES",npt,vpt);
    h->pm0a = new TProfile("pm0a",";p_{T,avp} (GeV);MPF0",npt,vpt);
    h->pm2a = new TProfile("pm2a",";p_{T,avp} (GeV);MPF2",npt,vpt);
    h->pmna = new TProfile("pmna",";p_{T,avp} (GeV);MPFn",npt,vpt);
    h->pmua = new TProfile("pmua",";p_{T,avp} (GeV);MPFu",npt,vpt);
    //
    h->presj = new TProfile("presj",";p_{T,jet} (GeV);JES",npt,vpt);
    h->pm0j = new TProfile("pm0j",";p_{T,jet} (GeV);MPF0",npt,vpt);
    h->pm2j = new TProfile("pm2j",";p_{T,jet} (GeV);MPF2",npt,vpt);
    h->pmnj = new TProfile("pmnj",";p_{T,jet} (GeV);MPFn",npt,vpt);
    h->pmuj = new TProfile("pmuj",";p_{T,jet} (GeV);MPFu",npt,vpt);
    
  } // doGamjet
  
  // Results similar to Dijet2 directory in dijet package
  gamjetHistos2 *hg2(0);
  if (doGamjet2) {

    // L2Res pT binning (central+forward hybrid)
    double vptd[] = 
    //{59.,85.,104.,170.,236., 302., 370., 460., 575.}; // central
    //{86., 110., 132., 204., 279., 373.} // forward
      {15, 21, 28, 37, 49,
       59, 86, 110, 132, 170, 204, 236, 279, 302, 373, 460, 575,
       638, 737, 846, 967, 1101, 1248,
       1410, 1588, 1784, 2000, 2238, 2500, 2787, 3103};
    // ^where did this original binning come from? How is that motivated? (I think it comes from dijet)

    double nptd = sizeof(vptd)/sizeof(vptd[0])-1;


    // L3Res (gamma+jet) pT binning adapted and extended
    const double vpt[] = {15, 20, 25, 30, 35,
			  40, 50, 60, 70, 85, 100, 125, 155, 180, 210, 250, 300,
			  350, 400, 500, 600, 800, 1000, 1200, 1500,
			  1800, 2100, 2400, 2700, 3000};

    double npt = sizeof(vpt)/sizeof(vpt[0])-1;

    // Current L2Res |eta| binning from Jindrich
    // https://indico.cern.ch/event/1263476/contributions/5311425/attachments/2612023/4513129/L2Res+HDM-March15.pdf
    double vxd[] =
    //{0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322,
    // 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
    // Newer L2Res |eta| binning from Mikel
    // https://indico.cern.ch/event/1335203/#7-update-on-l2res-for-2022-rer
    //  {0., 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.839, 4.013, 4.583, 5.191};
    //binning from dijet analysis (added on 22.04.2024), version w12:
    {0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
       1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 
       3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

    const int nxd = sizeof(vxd)/sizeof(vxd[0])-1;

    if (debug) cout << "Setup doGamjet2 " << endl << flush;
      
    fout->mkdir("Gamjet2");
    fout->cd("Gamjet2");

    gamjetHistos2 *h = new gamjetHistos2();
    hg2 = h;

    // Counting of events, and JEC L2L3Res+JERSF for undoing
    h->h2pteta = new TH2D("h2pteta",";#eta;p_{T,avp} (GeV);"
			  "N_{events}",nxd,vxd, nptd, vptd);
    h->p2res = new TProfile2D("p2res",";#eta;p_{T,avp} (GeV);JES",
			      nxd,vxd, nptd, vptd);
    h->p2jsf = new TProfile2D("p2jsf",";#eta;p_{T,avp} (GeV);JERSF",
			      nxd,vxd, nptd, vptd);
       
    // MPF decomposition for HDM method
    h->p2m0 = new TProfile2D("p2m0",";#eta;p_{T,avp} (GeV);MPF0",
			     nxd,vxd, nptd, vptd);
    h->p2m2 = new TProfile2D("p2m2",";#eta;p_{T,avp} (GeV);MPF2",
			     nxd,vxd, nptd, vptd);
    h->p2mn = new TProfile2D("p2mn",";#eta;p_{T,avp} (GeV);MPFn",
			     nxd,vxd, nptd, vptd);
    h->p2mu = new TProfile2D("p2mu",";#eta;p_{T,avp} (GeV);MPFu",
			     nxd,vxd, nptd, vptd);
    
    h->p2m0x = new TProfile2D("p2m0x",";#eta;p_{T,avp} (GeV);"
			      "MPF0X (MPFX)",nxd,vxd, nptd, vptd, "S");
    h->p2m2x = new TProfile2D("p2m2x",";#eta;p_{T,avp} (GeV);"
			      "MPF2X (DBX)",nxd,vxd, nptd, vptd, "S");

    // Extra for FSR studies
    h->p2mnu = new TProfile2D("p2mnu",";#eta;p_{T,avp} (GeV);MPFnu",
			      nxd,vxd, nptd, vptd);
    h->p2mnx = new TProfile2D("p2mnx",";#eta;p_{T,avp} (GeV);"
			      "MPFNX",nxd,vxd, nptd, vptd, "S");
    h->p2mux = new TProfile2D("p2mux",";#eta;p_{T,avp} (GeV);"
			      "MPFUX",nxd,vxd, nptd, vptd, "S");
    h->p2mnux = new TProfile2D("p2mnux",";#eta;p_{T,avp} (GeV);"
			       "MPFNUX",nxd,vxd, nptd, vptd, "S");

    /*
    h->h2ptetatc = new TH2D("h2ptetatc",";#eta;p_{T,tag} (GeV);"
			    "N_{events}",nxd,vxd, nptd, vptd);
    h->p2restc = new TProfile2D("p2restc",";#eta;p_{T,tag} (GeV);"
				"JES(probe)/JES(tag)",
				nxd,vxd, nptd, vptd);
    h->p2jsftc = new TProfile2D("p2jsftc ",";#eta;p_{T,tag} (GeV);"
				"JERSF(probe)/JERSF(tag)",
				nxd,vxd, nptd, vptd);
    h->p2m0tc = new TProfile2D("p2m0tc",";#eta;p_{T,tag} (GeV);MPF0",
			       nxd,vxd, nptd, vptd);
    h->p2m2tc = new TProfile2D("p2m2tc",";#eta;p_{T,tag} (GeV);MPF2",
			       nxd,vxd, nptd, vptd);
    h->p2mntc = new TProfile2D("p2mntc",";#eta;p_{T,tag} (GeV);MPFn",
			       nxd,vxd, nptd, vptd);
    h->p2mutc = new TProfile2D("p2mutc",";#eta;p_{T,tag} (GeV);MPFu",
			       nxd,vxd, nptd, vptd);
    
    h->h2ptetapf = new TH2D("h2ptetapf",";#eta;p_{T,probe} (GeV);"
			    "N_{events}",nxd,vxd, nptd, vptd);
    h->p2respf = new TProfile2D("p2respf",";#eta;p_{T,probe} (GeV);"
				"JES(probe)/JES(tag)",
				nxd,vxd, nptd, vptd);
    h->p2jsfpf = new TProfile2D("p2jsfpf",";#eta;p_{T,probe} (GeV);"
				"JERSF(probe)/JERSF(tag)",
				nxd,vxd, nptd, vptd);
    h->p2m0pf = new TProfile2D("p2m0pf",";#eta;p_{T,probe} (GeV);MPF0",
			       nxd,vxd, nptd, vptd);
    h->p2m2pf = new TProfile2D("p2m2pf",";#eta;p_{T,probe} (GeV);MPF2",
			       nxd,vxd, nptd, vptd);
    h->p2mnpf = new TProfile2D("p2mnpf",";#eta;p_{T,probe} (GeV);MPFn",
			       nxd,vxd, nptd, vptd);
    h->p2mupf = new TProfile2D("p2mupf",";#eta;p_{T,probe} (GeV);MPFu",
			       nxd,vxd, nptd, vptd);
    */    
  } // doGamjet2
  

  fout->cd();
  
  // Loop to create histograms and profiles
  // Match ordering to Lyon files (alpha->eta->data/MC) when creating
  // Although otherwise ordering is data/MC->eta->alpha
  // Add PS weight variations
  unsigned int nps = (isMC ? nPSWeightMax+1 : 1);
  //unsigned int nps = ((isMC && !isRun3) ? nPSWeightMax+1 : 1);
  map<int, map<int, map<int, BasicHistos*> > > mBasicHistos;
  for (unsigned int ialpha = 0; ialpha != alphas.size(); ++ialpha) {
    for (unsigned int ieta = 0; ieta != etas.size(); ++ieta) { 
    for (unsigned int ips = 0; ips != nps; ++ips) {

      // Select data/MC, alpha and eta bin
      const char *cd = dir.c_str();
      int ia = int(100*alphas[ialpha]);
      int iy1 = int(10*etas[ieta].first);
      int iy2 = int(10*etas[ieta].second);
      int iy = 100*int(iy1) + int(iy2);
      const char *cps = (ips==0 ? "" : Form("_ps%d",ips-1));
      
      // Counts
      TH1D *hn = new TH1D(Form("resp_MPFchs_%s_a%d_eta%02d_%02d_RawNEvents"
			       "_data_vs_pt%s",cd,ia,iy1,iy2,cps),"",nx,vx);
      TH1D *hxsec = new TH1D(Form("resp_MPFchs_%s_a%d_eta%02d_%02d_Xsec"
				  "_data_vs_pt%s",cd,ia,iy1,iy2,cps),"",nx,vx);
      
      // Response profiles
      string name = Form("resp_%%s_%s_a%d_eta%02d_%02d%s",cd,ia,iy1,iy2,cps);
      const char *cname = name.c_str();
      TProfile *prpt  =  new TProfile(Form(cname,"PtGam"),"",nx,vx);
      TProfile *prbal =  new TProfile(Form(cname,"PtBalchs"),"",nx,vx);
      TProfile *prdb =  new TProfile(Form(cname,"DBchs"),"",nx,vx);
      TProfile *prmpf =  new TProfile(Form(cname,"MPFchs"),"",nx,vx);
      TProfile *prmpf1 = new TProfile(Form(cname,"MPFR1chs"),"",nx,vx);
      TProfile *prmpfn = new TProfile(Form(cname,"MPFRnchs"),"",nx,vx);
      TProfile *prmpfu = new TProfile(Form(cname,"MpfRuchs"),"",nx,vx);
      TProfile *prho = new TProfile(Form(cname,"Rho_CHS"),"",nx,vx);
      TProfile *pdjes = new TProfile(Form(cname,"DeltaJES"),"",nx,vx);
      TProfile *pjes = new TProfile(Form(cname,"JES"),"",nx,vx);
      TProfile *pres = new TProfile(Form(cname,"RES"),"",nx,vx);

      // Store links to histograms and profiles into maps
      BasicHistos *pmh = new BasicHistos();
      BasicHistos& mh = (*pmh);
      mh.hn = hn;
      mh.hxsec = hxsec;
      mh.prpt = prpt;
      mh.prbal = prbal;
      mh.prdb = prdb;
      mh.prmpf = prmpf;
      mh.prmpf1 = prmpf1;
      mh.prmpfn = prmpfn;
      mh.prmpfu = prmpfu;
      mh.prho = prho;
      mh.pdjes = pdjes;
      mh.pjes = pjes;
      mh.pres = pres;
      mBasicHistos[iy][ia][ips] = pmh;
    } // for ips in PSWeight
    } // for ieta in etas
  } // for ialpha in alphas
  
  curdir->cd();
  
  TLorentzVector gam, gami, lhe, gengam, phoj, phoj0, phoj0off, jet, jet2, jetn;
  TLorentzVector gamorig; // for QCD bkg
  TLorentzVector met, met1, metn, metu, metnu, rawmet, corrmet, rawgam;
  TLorentzVector jeti, corrjets, rawjet, rawjets, rcjet, rcjets, rcoffsets;
  TLorentzVector geni, genjet, genjet2;
  TLorentzVector fox; // for isQCD
  TLorentzVector gamx; // for MPFX
  
  //Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nentries = fChain->GetEntries(); // Long startup time
  cout << "\nStarting loop over " << dataset << " with "
       << nentries << " entries" << endl;

  if (isMG && nentries!=nMG) {
    cout << "Nentries = "<<nentries<<", expected nMG = "<<nMG<<endl << flush;
     //assert(false);
    cout << "Recalculate HT bin counts prior to starting."
	 << " This will take a few minutes" << endl;
    hnevt->Reset();
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      b_LHE_HT->GetEntry(ientry); //read only this branch
      hnevt->Fill(LHE_HT);
      b_genWeight->GetEntry(ientry); //and this, too
      hsumw->Fill(LHE_HT, genWeight);
      if (jentry%1000000==0) cout << "." << flush;
      if (jentry%50000000==0 && jentry!=0) cout << "\nn="<<jentry<<endl<<flush;
    } // for jentry
    nMG = nentries;
    cout << "\nProcessed " << nMG << " entries" << endl << flush;
    //
    cout << Form("int vnevt[%d] = ",hnevt->GetNbinsX());
    for (int i = 1; i != hnevt->GetNbinsX()+1; ++i) {
      cout<<Form("%s%d",(i==1 ? "{" : ", "),int(hnevt->GetBinContent(i)+0.5));
    }
    cout << "}; // " << dataset << endl << flush;
    //
    cout << Form("double vsumw[%d] = ",hsumw->GetNbinsX());
    for (int i = 1; i != hsumw->GetNbinsX()+1; ++i) {
      cout<<Form("%s%1.4g",(i==1 ? "{" : ", "),hsumw->GetBinContent(i));
    }
    cout << "}; // " << dataset << endl << flush;
  } // isMC && nentries!=nMG
  
  //int skip = 21700000; // 2018A first events without 110EB
  //int skip = 55342793; // 2018A first events with 92 photon
  //int skip = 265126992; // 2018A first events with 191 photons, 23 jets
  //int skip = 14648973; // 2017C bad HDM

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) { //is this the start of the analysis loop? (trying to find event 0..)

    // Skip events, typically for debugging purposes
    //if (jentry<skip) continue;
    //if (_gh_debug && jentry%10000==0) cout << "," << endl << flush;
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    ++_ntot;
    if (jentry==100000 || jentry==1000000 || jentry==1000000 ||
	(jentry%1000000==0 && jentry<10000000) ||
	(jentry%10000000==0 && jentry!=0) ||
	jentry==nentries-1) {
      if (jentry==0) { laptime.Start(); }
      if (nentries!=0) {
	cout << Form("\nProcessed %lld events (%1.1f%%) in %1.0f sec. "
		     "(%1.0f sec. for last %d)",
		     jentry, 100.*jentry/nentries, fulltime.RealTime(),
		     laptime.RealTime(), nlap);
      }
      if (jentry!=0 && nlap!=0) {
	cout << Form("\nEstimated runtime:  %1.0f sec. "
		     " (%1.0f sec. for last %d)",
		     1.*nentries/jentry*fulltime.RealTime(),
		     1.*nentries/nlap*laptime.RealTime(),nlap) << flush;
	laptime.Reset();
	nlap = 0;
      }
      if (jentry==0) fulltime.Reset(); // Leave out initialization time
      fulltime.Continue();
      laptime.Continue();
    }
    if (jentry%10000==0) cout << "." << flush;
    ++nlap;

    //Safety resets for triggers only in 2024
    HLT_Photon50EB_TightID_TightIso = HLT_Photon55EB_TightID_TightIso =
      HLT_Photon75EB_TightID_TightIso = HLT_Photon90EB_TightID_TightIso =
      kFALSE;
    // Safety resets for triggers only in 2022-2023
    HLT_Photon30EB_TightID_TightIso = kFALSE;
    // Safety resets for tight triggers that are partly missing in 2018A, 2016
    HLT_Photon100EB_TightID_TightIso =  HLT_Photon110EB_TightID_TightIso =
      HLT_Photon120EB_TightID_TightIso = kFALSE;
    // Safety resets for triggers that are partly missing in 2016, or 2017-18
    HLT_Photon200 = HLT_Photon175 = HLT_Photon150 = HLT_Photon90 =
      HLT_Photon75 = HLT_Photon50 = HLT_Photon36 = HLT_Photon33 =
      HLT_Photon30 = HLT_Photon22 = HLT_Photon20 = kFALSE;
    // Safety resets for medium triggers
    HLT_Photon165_R9Id90_HE10_IsoM = HLT_Photon120_R9Id90_HE10_IsoM = 
      HLT_Photon90_R9Id90_HE10_IsoM = HLT_Photon75_R9Id90_HE10_IsoM = 
      HLT_Photon50_R9Id90_HE10_IsoM = HLT_Photon36_R9Id90_HE10_IsoM = 
      HLT_Photon30_R9Id90_HE10_IsoM = HLT_Photon22_R9Id90_HE10_IsoM = kFALSE;
    // Safety resets for loose triggers (only 2017-18)
    HLT_Photon60_HoverELoose = HLT_Photon50_HoverELoose =
      HLT_Photon40_HoverELoose =  HLT_Photon30_HoverELoose =
      HLT_Photon20_HoverELoose = kFALSE;
    HLT_Photon165_HE10 = kFALSE;

    if (!isMC) { // Fast trigger filtering (useful for data)

      if (b_HLT_Photon200 && !is16) // not in 2016
	b_HLT_Photon200->GetEntry(ientry);
      if (b_HLT_Photon175)
	b_HLT_Photon175->GetEntry(ientry);
      if (b_HLT_Photon150 && !is16) // not in 2016
	b_HLT_Photon150->GetEntry(ientry);
      if (b_HLT_Photon120)
	b_HLT_Photon120->GetEntry(ientry);
      if (b_HLT_Photon90)
	b_HLT_Photon90->GetEntry(ientry);
      if (b_HLT_Photon75)
	b_HLT_Photon75->GetEntry(ientry);
      if (b_HLT_Photon50)
	b_HLT_Photon50->GetEntry(ientry);
      if (b_HLT_Photon36 && is16) // 2016
	b_HLT_Photon36->GetEntry(ientry);
      if (b_HLT_Photon30 && is16) // 2016
	b_HLT_Photon30->GetEntry(ientry);
      if (b_HLT_Photon33 && !is16) // 2016
	b_HLT_Photon33->GetEntry(ientry);
      if (b_HLT_Photon22 && is16) // 2016
	b_HLT_Photon22->GetEntry(ientry);
      //if (b_HLT_Photon20 && (is18 || is22 || is23)) // 2018
      if (b_HLT_Photon20 && is18) // 2018
	b_HLT_Photon20->GetEntry(ientry);

      // Only in (most of) 2018, and now some in 22-23
      if (b_HLT_Photon120EB_TightID_TightIso && is18) // not in 2016-17, 2018A
	b_HLT_Photon120EB_TightID_TightIso->GetEntry(ientry);
      if (b_HLT_Photon110EB_TightID_TightIso && (is18 || is22 || is23 || is24)) // not in 2016-17, 2018A
	b_HLT_Photon110EB_TightID_TightIso->GetEntry(ientry);
      if (b_HLT_Photon100EB_TightID_TightIso && is18) // not in 2016-17, 2018A
	b_HLT_Photon100EB_TightID_TightIso->GetEntry(ientry);

      // Only 24
      if (b_HLT_Photon50EB_TightID_TightIso && is24)
	b_HLT_Photon50EB_TightID_TightIso->GetEntry(ientry);
      if (b_HLT_Photon55EB_TightID_TightIso && is24)
	b_HLT_Photon55EB_TightID_TightIso->GetEntry(ientry);
      if (b_HLT_Photon75EB_TightID_TightIso && is24)
	b_HLT_Photon75EB_TightID_TightIso->GetEntry(ientry);
      if (b_HLT_Photon90EB_TightID_TightIso && is24)
	b_HLT_Photon90EB_TightID_TightIso->GetEntry(ientry);

      // Only 22-23-24
      if (b_HLT_Photon100EBHE10 && isRun3)
	b_HLT_Photon100EBHE10->GetEntry(ientry);
      if (b_HLT_Photon30EB_TightID_TightIso && isRun3)
	b_HLT_Photon30EB_TightID_TightIso->GetEntry(ientry);

      // Only 2016
      if (b_HLT_Photon22_R9Id90_HE10_IsoM && is16) // only in 2016
	b_HLT_Photon22_R9Id90_HE10_IsoM->GetEntry(ientry);
      if (b_HLT_Photon30_R9Id90_HE10_IsoM && is16) // only in 2016
	b_HLT_Photon30_R9Id90_HE10_IsoM->GetEntry(ientry);
      if (b_HLT_Photon36_R9Id90_HE10_IsoM && is16) // only in 2016
	b_HLT_Photon36_R9Id90_HE10_IsoM->GetEntry(ientry);
      if (b_HLT_Photon165_HE10 && is16)
	b_HLT_Photon165_HE10->GetEntry(ientry);

      // In all years
      if (b_HLT_Photon165_R9Id90_HE10_IsoM)
	b_HLT_Photon165_R9Id90_HE10_IsoM->GetEntry(ientry);
      if (b_HLT_Photon120_R9Id90_HE10_IsoM)
	b_HLT_Photon120_R9Id90_HE10_IsoM->GetEntry(ientry);
      if (b_HLT_Photon90_R9Id90_HE10_IsoM)
	b_HLT_Photon90_R9Id90_HE10_IsoM->GetEntry(ientry);
      if (b_HLT_Photon75_R9Id90_HE10_IsoM)
	b_HLT_Photon75_R9Id90_HE10_IsoM->GetEntry(ientry);
      if (b_HLT_Photon50_R9Id90_HE10_IsoM)
	b_HLT_Photon50_R9Id90_HE10_IsoM->GetEntry(ientry);
							
      // Only in 2017
      if (b_HLT_Photon60_HoverELoose && is17) // only in 2017
	b_HLT_Photon60_HoverELoose->GetEntry(ientry);
      if (b_HLT_Photon50_HoverELoose && is17) // only in 2017
	b_HLT_Photon50_HoverELoose->GetEntry(ientry);
      if (b_HLT_Photon40_HoverELoose && is17) // only in 2017
	b_HLT_Photon40_HoverELoose->GetEntry(ientry);

      // 2017 and 2018
      if (b_HLT_Photon30_HoverELoose && !is16) // not in 2016
	b_HLT_Photon30_HoverELoose->GetEntry(ientry);
      if (b_HLT_Photon20_HoverELoose && !is16) // not in 2016
	b_HLT_Photon20_HoverELoose->GetEntry(ientry);
      
      if ((isRun3 &&
	   !(HLT_Photon200 ||
	     HLT_Photon175 || 
	     HLT_Photon150 || 
	     HLT_Photon120 || 
	     HLT_Photon90 || 
	     HLT_Photon75 || 
	     HLT_Photon50 || 
	     HLT_Photon33 || 
	     //HLT_Photon20 ||
	     //HLT_Photon120EB_TightID_TightIso ||
	     HLT_Photon110EB_TightID_TightIso ||
	     //HLT_Photon100EB_TightID_TightIso ||
	     HLT_Photon100EBHE10 ||
             HLT_Photon90EB_TightID_TightIso ||
	     HLT_Photon75EB_TightID_TightIso ||
	     HLT_Photon55EB_TightID_TightIso ||
	     HLT_Photon50EB_TightID_TightIso ||
	     HLT_Photon30EB_TightID_TightIso ||
	     HLT_Photon90_R9Id90_HE10_IsoM ||
	     HLT_Photon75_R9Id90_HE10_IsoM ||
	     HLT_Photon50_R9Id90_HE10_IsoM ||
	     HLT_Photon30_HoverELoose ||
	     HLT_Photon20_HoverELoose)) ||
	  (is18 &&
	   !(HLT_Photon200 ||
	     HLT_Photon175 || 
	     HLT_Photon150 || 
	     HLT_Photon120 || 
	     HLT_Photon90 || 
	     HLT_Photon75 || 
	     HLT_Photon50 || 
	     HLT_Photon33 || 
	     HLT_Photon20 ||
	     HLT_Photon120EB_TightID_TightIso ||
	     HLT_Photon110EB_TightID_TightIso ||
	     HLT_Photon100EB_TightID_TightIso ||
	     HLT_Photon90_R9Id90_HE10_IsoM ||
	     HLT_Photon75_R9Id90_HE10_IsoM ||
	     HLT_Photon50_R9Id90_HE10_IsoM ||
	     HLT_Photon30_HoverELoose ||
	     HLT_Photon20_HoverELoose)) ||
	   (is17 &&
	   !(HLT_Photon200 ||
	     HLT_Photon175 || 
	     HLT_Photon150 || 
	     HLT_Photon120 || 
	     HLT_Photon90 || 
	     HLT_Photon75 || 
	     HLT_Photon50 || 
	     HLT_Photon33 || 
	     HLT_Photon165_R9Id90_HE10_IsoM ||
	     HLT_Photon120_R9Id90_HE10_IsoM ||
	     HLT_Photon90_R9Id90_HE10_IsoM ||
	     HLT_Photon75_R9Id90_HE10_IsoM ||
	     HLT_Photon50_R9Id90_HE10_IsoM ||
	     HLT_Photon60_HoverELoose ||
	     HLT_Photon50_HoverELoose ||
	     HLT_Photon40_HoverELoose ||
	     HLT_Photon30_HoverELoose ||
	     HLT_Photon20_HoverELoose)) ||
	  (is16 &&
	   !(HLT_Photon175 || 
	     HLT_Photon120 || 
	     HLT_Photon90 || 
	     HLT_Photon75 || 
	     HLT_Photon50 || 
	     HLT_Photon36 || 
	     HLT_Photon30 || 
	     HLT_Photon22 || 
	     //HLT_Photon165_HE10 ||
	     HLT_Photon165_R9Id90_HE10_IsoM ||
	     HLT_Photon120_R9Id90_HE10_IsoM ||
	     HLT_Photon90_R9Id90_HE10_IsoM ||
	     HLT_Photon75_R9Id90_HE10_IsoM ||
	     HLT_Photon50_R9Id90_HE10_IsoM ||
	     HLT_Photon36_R9Id90_HE10_IsoM ||
	     HLT_Photon30_R9Id90_HE10_IsoM ||
	     HLT_Photon22_R9Id90_HE10_IsoM))) {
	++_nbadevents_trigger;
	continue;
      }
    } // !isMC

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    // Safety check for rho being NaN
    if (!(fixedGridRhoFastjetAll>=0 && fixedGridRhoFastjetAll<150))
      fixedGridRhoFastjetAll = 34; // average conditions
    
    // Sanity check PS weights
    if (!isMC) { nPSWeight = 0; }
    //if (!isMC || is22 || is23) { nPSWeight = 0; }
		//cout << "\n" << "nPSWeight: " << nPSWeight << endl; //test
		//cout << "nPSWeightMax: " << nPSWeightMax << endl; //test

		//check first event and whether it has 4 weight variations or 44
		/*
		if(jentry==0){
			cout << "\n" << "nPSWeight: " << nPSWeight << endl;
			if(nPSWeight==4){
				nPSWeightMax = 4;
			}
			else if(nPSWeight==44){
				nPSWeightMax = 44;
			}
		}
		*/
    assert(nPSWeight<=nPSWeightMax);

    // Does the run/LS pass the latest JSON selection?
    if (!isMC && _json[run][luminosityBlock]==0) {
      //_badjson.insert(pair<int, int>(run, lbn));
      ++_nbadevents_json;
      continue;
    }
    else 
      ++_nevents;


    // Safety checks for array sizes (segfault in 2018A)
    if (nJet > nJetMax || nPhoton > nPhotonMax) {
      cout << "Error: nJet="<<nJet<<" or nPhoton="<<nPhoton
	   << " exceeding maximum. Skip event " << jentry << endl << flush;
      exit(0);
    }
    
    // Re-JEC + re-MET, if needed
    // Skipped for now      
    
    // Select leading photon. Use tight cut-based ID and PF relative isolation
    // Temporary: select photon based on LHE photon match
    int iGamGen(-1), iGam(-1), nGam(0);
    int iGamOrig(-1); // for QCD bkg
    gengam.SetPtEtaPhiM(0,0,0,0);
    gam.SetPtEtaPhiM(0,0,0,0);
    rawgam.SetPtEtaPhiM(0,0,0,0);
    phoj.SetPtEtaPhiM(0,0,0,0);
    phoj0.SetPtEtaPhiM(0,0,0,0);

    // Gen-photon
    if (isMC && nGenIsolatedPhoton>0) {
      gengam.SetPtEtaPhiM(GenIsolatedPhoton_pt[0],GenIsolatedPhoton_eta[0],
			  GenIsolatedPhoton_phi[0],GenIsolatedPhoton_mass[0]);
    }

    // Select tight photons and photon matching gen photon
    for (int i = 0; i != nPhoton; ++i) {

      if (isRun3) Photon_mass[i] = 0;
      gami.SetPtEtaPhiM(Photon_pt[i],  Photon_eta[i],
			Photon_phi[i], Photon_mass[i]);
      
      // Photon matching gen photon
      if (gengam.Pt()>0 && gengam.DeltaR(gami)<0.2 && iGamGen==-1) {
	iGamGen = i;
      } 
      
      // Leading tight photon(s)
      // R9>0.94 to avoid bias wrt R9Id90 triggers and from photon conversions
      if (Photon_pt[i]>15 && Photon_cutBased[i]==3 && Photon_hoe[i]<0.02148 && // id cut should be >=3, in case there's tighter ones
	  Photon_r9[i]>0.94 && Photon_r9[i]<1.00) { //added photon cut R9<1.00
	++nGam;
	if (iGam==-1) {
	  iGam = i;
	  gam = gami;
	}
      } // tight photon
    } // for i in nPhoton

    // Correct photon for gain1 and MPF for "footprint" (photon vs PFgamma)
		// TO DO: in the same way as for gain1, we should also correct for jump in the end of 24C
    rawgam = gam;
    if (iGam!=-1 && Photon_seedGain[iGam]==1 && !isMC) {
      //gam *= 1./1.01;
      // minitools/drawGainVsPt.C (add R_6/12+R_1/6, take MPF+statTowardsDB)
      if (!isRun3) gam *= 1./1.011; // MPF=1.13+/-0.04%, DB=1.05+/-0.08%
      if ( isRun3) gam *= 1./1.017; // MPF=1.74+/-0.07%, DB=1.41+/-0.16%
    }
    if (iGam!=-1 && !isRun3) {
      // [0]+log(x)*([1]+log(x)*[2]) in range [15,1750] to MC pphoj0
      //1  p0           4.57516e-02   3.91871e-04   1.09043e-07   4.17033e-05
      //2  p1          -1.27462e-02   1.50968e-04   2.08432e-08   3.92715e-03
      //3  p2           1.07760e-03   1.45293e-05   3.93020e-09   1.65460e-01
      double x = max(60.,rawgam.Pt());
      double f = 4.57516e-02 + log(x) * ( -1.27462e-02 + log(x) * 1.07760e-03);
      rawgam *= (1+f);
    }

    //TO DO
    /*
    //create new script to calculate the factor needed for scale correction
    if (is24 && run=< NUMBER && iGam!=-1 && !isMC) {
        gam *= 1./1.01 //from comparison Z vs gamma, to account for the "step"/"jump", to permille precision
    }
    */
     
   
    // Photon-jet: uncorrected jet minus (uncorr.) photon minus L1RC
    if (iGam!=-1 && Photon_jetIdx[iGam]!=-1) {
      //assert(Photon_jetIdx[iGam]>=-1);
      //assert(Photon_jetIdx[iGam]<nJet);
      int idx = Photon_jetIdx[iGam];
      if (!(idx<nJet && idx>=0)) {
	cout << endl << "idx = " << idx << " nJet = " << nJet << endl << flush;
	cout << "Skip event " << event << " in LS " << luminosityBlock
	     << " in  run " << run << " in File: " << _filename << endl;
	continue;
	//assert(idx<nJet);
      }
      phoj.SetPtEtaPhiM(Jet_pt[idx], Jet_eta[idx], Jet_phi[idx], Jet_mass[idx]);
      phoj *= (1-Jet_rawFactor[idx]);
      if (rawgam.DeltaR(phoj)<0.4) { // does not always hold in Run3
	//phoj.Pt() >= rawgam.Pt()) { // not always true in Run3 (esp. 2022C)
	// 2022 data is missing proper Puppi photon protection for jets
	// (but MET ok?)
	//double r22 = max(0.15,min(1.0,(rawgam.Pt()-20.)/180.));
	//phoj -= (is22v10 ? r22*rawgam : rawgam);
	phoj -= rawgam; // NanoV12
      }
      else {
	if (cntErrDR<10) {
	  cout << endl << "entry " << jentry << ", rawgam.DeltaR(phoj) = "
	       << rawgam.DeltaR(phoj) << endl << flush;
	  cout << "Skip event " << event << " in LS " << luminosityBlock
	       << " in  run " << run << " in File: " << _filename << endl;
	  ++cntErrDR;
	  if (cntErrDR==10) {
	    cout << "Stop reporting rawgam.DeltaR, silently skip." << endl;
	  }
	}
	phoj.SetPtEtaPhiM(0,0,0,0);
	continue;
      }
      phoj0 = phoj;

      // Calculate L1RC correction
      double corrl1rc(1.); // isRun3
      if (isRun2) {
	jecl1rc->setJetPt(phoj.Pt());
	jecl1rc->setJetEta(phoj.Eta());
	jecl1rc->setJetA(Jet_area[idx]);
	jecl1rc->setRho(fixedGridRhoFastjetAll);
	corrl1rc = jecl1rc->getCorrection();
      }
      phoj *= corrl1rc;

      // Calculate L1RC correction without "zero suppression"
      double refpt = 30; // phoj.Pt~0 leads to negative offset cutoff
      double corrl1rc0(1.); // isRun3
      if (isRun2) {
	jecl1rc->setJetPt(refpt);
	jecl1rc->setJetEta(phoj0.Eta());
	jecl1rc->setJetA(Jet_area[idx]);
	jecl1rc->setRho(fixedGridRhoFastjetAll);
	corrl1rc0 = jecl1rc->getCorrection();
      }
      double off0 = (corrl1rc0 - 1) * refpt; // corr*ptref = (ptref-off)
      phoj0off.SetPtEtaPhiM(off0,phoj0.Eta(),phoj0.Phi(),0.);
      phoj0 -= phoj0off;
    }

    // For QCD background, emulate a photon+jet event by replacing
    // one of the leading jets with genjet taking the place of the photon
    int iFox(-1);
    fox.SetPtEtaPhiM(0,0,0,0);
    //if (isQCD && iGam==-1 && nJet>=2) {
    if (isQCD && nJet>=2) {
      // Save original good photon, if one was found
      iGamOrig = iGam;
      gamorig = gam;

      iFox = (jentry%2); // "random" selection from two leading jets
      // Jet_genJetIdx would be great, but only there for UL18 nAOD? Maybe there
      int k = Jet_genJetIdx[iFox];
      if (k>=0 && k<nGenJet) {
	gam.SetPtEtaPhiM(GenJet_pt[k], GenJet_eta[k], GenJet_phi[k],
			 GenJet_mass[k]);
	// NB: should remove UE clustered into gam. In Minsuk's rho_ZB_new.pdf
	// QCD_CP5 has about 3.5 GeV/A of UE offset at generator level
	double area = Jet_area[iFox];
	gam *= (gam.Pt()>0 ? 1 - 3.5*area/gam.Pt() : 1.);
	rawgam = gam;

	fox.SetPtEtaPhiM(Jet_pt[iFox], Jet_eta[iFox], Jet_phi[iFox],
			 Jet_mass[iFox]);
	fox *= (1-Jet_rawFactor[iFox]);
	// Calculate L1RC correction
	double corrl1rc(1.); // isRun3
	if (isRun2) {
	  jecl1rc->setJetPt(fox.Pt());
	  jecl1rc->setJetEta(fox.Eta());
	  jecl1rc->setJetA(Jet_area[iFox]);
	  jecl1rc->setRho(fixedGridRhoFastjetAll);
	  corrl1rc = jecl1rc->getCorrection();
	}
	fox *= corrl1rc;
	// NB2: should also remove UE clustered into fox. In Minsuk's plot
	// QCD_CP5 has about 2.5 GeV/A of UE offset at FullSim level
	fox *= (fox.Pt()>0 ? 1 - 2.5*area/fox.Pt() : 1.);

	// NB3: For consistency with gamma+jet, phoj should still have reco UE
	phoj.SetPtEtaPhiM(2.5*area,0,0,0);
	phoj0 = phoj;
	nGam = 1;
      }
    } // isQCD
  
    // Event weights (1 for MadGraph)
    //bool isMC = (run==1);
    assert((isMC && run==1) || (!isMC && run!=1));
    double w = (isMC ? genWeight : 1);
    if (isMG) {
      int iht = hxsec->FindBin(LHE_HT);
      double xsec = hxsec->GetBinContent(iht);
      //double nevt = hnevt->GetBinContent(iht);
      //double wht = (nevt ? xsec / nevt : 1);
      double sumw = hsumw->GetBinContent(iht);
      double wht = (sumw ? xsec / sumw : 1);
      w *= wht;
      hLHE_HT->Fill(LHE_HT); // cross-check hnevt afterwards
      hHT->Fill(LHE_HT, w); // cross-check HT spectrum smoothness afterwards
    }

    //bool doPtHatFilter = true;
    //if (doPtHatFilter && isMC) {
    //if ( isMG && 2.*Pileup_pthatmax>LHE_HT) continue;
    //if (!isMG && Pileup_pthatmax>Generator_binvar) continue;
    //}
    
    // Pileup
    double TruePUrms(0);
    if (!isMC) Pileup_nTrueInt = getTruePU(run,luminosityBlock,&TruePUrms);
    double ptgam = gam.Pt();

    // Trigger selection. Take care to match pT bin edges
    // {15, 20, 25, 30, 35, 40, 50, 60, 70, 85, 105, 130, 175, 230,
    //  300, 400, 500, 600, 700, 850, 1000, 1200, 1450, 1750};
    // NB: Photon90 threshold could be 95, Photon175 coud be 185, if bins split?
    double pt = ptgam; // shorthand to make trigger selection more readable
    int itrg(0); // choose trigger for PU reweighing as side effect (hack...)
    bool pass_trig = 
      ((is16 && 
	((HLT_Photon175                  && pt>=230            && (itrg=175)) ||
	 //(HLT_Photon165_HE10             && pt>=175 && pt<230) || // not in MC
	 (HLT_Photon165_R9Id90_HE10_IsoM && pt>=175 && pt<230  && (itrg=165)) ||
	 (HLT_Photon120_R9Id90_HE10_IsoM && pt>=130 && pt<175  && (itrg=120)) ||
	 (HLT_Photon90_R9Id90_HE10_IsoM  && pt>=105 && pt<130  && (itrg=90)) ||
	 (HLT_Photon75_R9Id90_HE10_IsoM  && pt>=85  && pt<105  && (itrg=75)) ||
	 (HLT_Photon50_R9Id90_HE10_IsoM  && pt>=60  && pt<85   && (itrg=50)) ||
	 (HLT_Photon36_R9Id90_HE10_IsoM  && pt>=40  && pt<60   && (itrg=36)) ||
	 (HLT_Photon30_R9Id90_HE10_IsoM  && pt>=35  && pt<40   && (itrg=30)) ||
	 (HLT_Photon22_R9Id90_HE10_IsoM  && pt>=20  && pt<35   && (itrg=22)) ||
	 (isMC                           && pt>=40  && pt<60   && (itrg=36)) ||
	 (isMC                           && pt>=35  && pt<40   && (itrg=30)) ||
	 (isMC                           && pt>=20  && pt<35   && (itrg=22))
	 )) ||
       (is17 &&
	((HLT_Photon200                  && pt>=230            && (itrg=200)) ||
	 (HLT_Photon165_R9Id90_HE10_IsoM && pt>=175 && pt<230  && (itrg=165)) ||
	 (HLT_Photon120_R9Id90_HE10_IsoM && pt>=130 && pt<175  && (itrg=120)) ||
	 (HLT_Photon90_R9Id90_HE10_IsoM  && pt>=105 && pt<130  && (itrg=90)) ||
	 (HLT_Photon75_R9Id90_HE10_IsoM  && pt>=85  && pt<105  && (itrg=75)) ||
	 (HLT_Photon50_R9Id90_HE10_IsoM  && pt>=60  && pt<85   && (itrg=50)) ||
	 //(HLT_Photon33                   && pt>=35  && pt<60 ) ||
	 (HLT_Photon30_HoverELoose       && pt>=35  && pt<60   && (itrg=30)) ||
	 (HLT_Photon20_HoverELoose       && pt>=20  && pt<35   && (itrg=20)) ||
	 //(HLT_Photon20                     && pt>=20  && pt<60 )
	 (isMC                           && pt>=35  && pt<60   && (itrg=30)) ||
	 (isMC                           && pt>=20  && pt<35   && (itrg=20))
	 )) ||
       (is18 &&
	((HLT_Photon200                    && pt>=230           && (itrg=200))||
	 (HLT_Photon110EB_TightID_TightIso && pt>=130 && pt<230 && (itrg=110))||
	 (HLT_Photon100EB_TightID_TightIso && pt>=105 && pt<130 && (itrg=100))||
	 (HLT_Photon90_R9Id90_HE10_IsoM    && pt>=95  && pt<105 && (itrg=90)) ||
	 (HLT_Photon75_R9Id90_HE10_IsoM    && pt>=85  && pt<95  && (itrg=75)) ||
	 (HLT_Photon50_R9Id90_HE10_IsoM    && pt>=60  && pt<85  && (itrg=50)) ||
	 //(HLT_Photon33                     && pt>=35  && pt<60 ) ||
	 (HLT_Photon30_HoverELoose         && pt>=35  && pt<60  && (itrg=30)) ||
	 (HLT_Photon20_HoverELoose         && pt>=20  && pt<35  && (itrg=20)) ||
	 //(HLT_Photon20                   && pt>=20  && pt<35 )
	 (isMC                           && pt>=35  && pt<60  && (itrg=30)) ||
	 (isMC                           && pt>=20  && pt<35  && (itrg=20))
	 )) ||
       // Push triggers to the limit for 22-23 (2022C bad 75,90)
       (isRun3 && (is22 || is23) &&
	((HLT_Photon200                 && pt>=230            && (itrg=200)) ||
	 (HLT_Photon110EB_TightID_TightIso && pt>=110&&pt<230 && (itrg=110)) ||
	 (HLT_Photon90_R9Id90_HE10_IsoM && pt>=90  && pt<110  && (itrg=90))  ||
	 (HLT_Photon75_R9Id90_HE10_IsoM && pt>=75  && pt<90   && (itrg=75))  ||
	 (HLT_Photon50_R9Id90_HE10_IsoM && pt>=50  && pt<75   && (itrg=50))  ||
       //(HLT_Photon50_R9Id90_HE10_IsoM && pt>=50  && pt<130  && (itrg=50))  ||
	 (HLT_Photon30EB_TightID_TightIso  && pt>=30&&pt<50   && (itrg=30))  ||
	 (HLT_Photon20_HoverELoose      && pt>=20  && pt<30   && (itrg=20))
	 //|| (true && (itrg=1))// trigger bypass for EGamma on photonTrigs.C
	 )) ||
       // Updated menu in 2024 with high rate isolated triggers
       (isRun3 && (is24) &&
	((HLT_Photon200                    && pt>=230         && (itrg=200)) ||
	 (HLT_Photon110EB_TightID_TightIso && pt>=110&&pt<230 && (itrg=110)) ||
	 (HLT_Photon50EB_TightID_TightIso  && pt>=50 &&pt<230 && (itrg=50))  ||
	 (HLT_Photon30EB_TightID_TightIso  && pt>=30 &&pt<50  && (itrg=30))  ||
	 (HLT_Photon20_HoverELoose         && pt>=20 && pt<30 && (itrg=20))
	 || (true && (itrg=1))// trigger bypass for EGamma on photonTrigs.C
	 ))
       );
  
    // Select trigger pT bins by hand for QCD. Error prone...
    if (isQCD && !pass_trig) {
      pass_trig = 
	((is16 && 
	  ((pt>=230            && (itrg=175)) ||
	   (pt>=175 && pt<230  && (itrg=165)) ||
	   (pt>=130 && pt<175  && (itrg=120)) ||
	   (pt>=105 && pt<130  && (itrg=90)) ||
	   (pt>=85  && pt<105  && (itrg=75)) ||
	   (pt>=60  && pt<85   && (itrg=50)) ||
	   (pt>=40  && pt<60   && (itrg=36)) ||
	   (pt>=35  && pt<40   && (itrg=30)) ||
	   (pt>=20  && pt<35   && (itrg=22))
	   )) ||
	 (is17 &&
	  ((pt>=230            && (itrg=200)) ||
	   (pt>=175 && pt<230  && (itrg=165)) ||
	   (pt>=130 && pt<175  && (itrg=120)) ||
	   (pt>=105 && pt<130  && (itrg=90)) ||
	   (pt>=85  && pt<105  && (itrg=75)) ||
	   (pt>=60  && pt<85   && (itrg=50)) ||
	   (pt>=35  && pt<60   && (itrg=30)) ||
	   (pt>=20  && pt<35   && (itrg=20))
	   )) ||
	 (is18 &&
	  ((pt>=230           && (itrg=200))||
	   (pt>=130 && pt<230 && (itrg=110))||
	   (pt>=105 && pt<130 && (itrg=100))||
	   (pt>=95  && pt<105 && (itrg=90)) ||
	   (pt>=85  && pt<95  && (itrg=75)) ||
	   (pt>=60  && pt<85  && (itrg=50)) ||
	   (pt>=35  && pt<60  && (itrg=30)) ||
	   (pt>=20  && pt<35  && (itrg=20))
	   )) ||
	 (isRun3 &&
	  ((pt>=230            && (itrg=200)) ||
	   (pt>=110&&pt<230 && (itrg=110)) ||
	   (pt>=90  && pt<110  && (itrg=90))  ||
	   (pt>=75  && pt<90   && (itrg=75))  ||
	   (pt>=50  && pt<75   && (itrg=50))  ||
	   (pt>=30&&pt<50   && (itrg=30))  ||
	   (pt>=20  && pt<30   && (itrg=20))
	   ))
	 );
    } // isQCD

    assert(itrg>0 || !pass_trig);

    // Reweight MC pileup (except for 22-23)
    if (isMC && pass_trig && !isRun3) {
      TH1D *hm = _pu[dataset][1]; assert(hm);
      TH1D *hd = _pu[sera][itrg];
      if (!hd) cout << "Missing _pu[sera="<<sera<<"][itrg="<<itrg<<"]"
		    << endl << flush;
      assert(hd);
      assert(hm->GetNbinsX()==hd->GetNbinsX());
      int k = hm->FindBin(Pileup_nTrueInt);
      assert(hm->GetBinLowEdge(k)==hd->GetBinLowEdge(k));
      double nm  = hm->GetBinContent(k);
      assert(nm>0); // should never get here if hm made from fullMC
      double nd  = hd->GetBinContent(k);
      double wt = (nm>0 ? nd / nm : 0);
      w *= wt;
    }
    // Normalize data luminosity (except for 22-23)
    if (!isMC && pass_trig && !isRun3) {
      double lumi = _lumi[sera][itrg];
      assert(lumi>0);
      w *= 1./lumi;
    }

    // Select leading jets. Just exclude photon, don't apply JetID yet
    Float_t         Jet_resFactor[nJetMax]; // Custom addition
    Float_t         Jet_deltaJES[nJetMax]; // Custom addition
    Float_t         Jet_CF[nJetMax]; // Custom addition
    int iJet(-1), iJet2(-1), nJets(0);
    double djes(1), jes(1), res(1);
    jet.SetPtEtaPhiM(0,0,0,0);
    jet2.SetPtEtaPhiM(0,0,0,0);
    jetn.SetPtEtaPhiM(0,0,0,0);
    // Also calculate corrected type-I chsMET and HDM inputs
    corrjets.SetPtEtaPhiM(0,0,0,0);
    rawjets.SetPtEtaPhiM(0,0,0,0);
    rcjets.SetPtEtaPhiM(0,0,0,0);
    rcoffsets.SetPtEtaPhiM(0,0,0,0);
    for (int i = 0; i != nJet; ++i) {
      
      // Redo JEC on the fly (should be no previous use of corrected jets)
      if (jec!=0) {
	
	double rawJetPt = Jet_pt[i] * (1.0 - Jet_rawFactor[i]);
	double rawJetMass = Jet_mass[i] * (1.0 - Jet_rawFactor[i]);
	jec->setJetPt(rawJetPt);
	jec->setJetEta(Jet_eta[i]);
	jec->setJetPhi(Jet_phi[i]); //added this line to make BPix work (should be marked as w3)
	if (isRun2) {
	  jec->setJetA(Jet_area[i]);
	  jec->setRho(fixedGridRhoFastjetAll);
	}
	//double corr = jec->getCorrection();
	vector<float> v = jec->getSubCorrections();
	double corr = v.back();
	double res = (v.size()>1 ? v[v.size()-1]/v[v.size()-2] : 1.);
	//Jet_RES[i] = 1./res;
	Jet_deltaJES[i] = (1./corr) / (1.0 - Jet_rawFactor[i]);
	Jet_pt[i] = corr * rawJetPt;
	Jet_mass[i] = corr * rawJetMass;
	Jet_rawFactor[i] = (1.0 - 1.0/corr);
	Jet_resFactor[i] = (1.0 - 1.0/res);
      }

      // Smear jets
      if (smearJets) {
	assert(false);
      }
      
      // Check that jet is not photon and pTcorr>15 GeV
      if (Jet_pt[i]>15 && (iGam==-1 || i != Photon_jetIdx[iGam]) &&
	  (!isQCD || i != iFox)) {
	
	//++nJets;
	jeti.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
 	if (gam.DeltaR(jeti)<0.2) continue; // should not happen, but does?
	++nJets;

	if (iJet==-1) { // Leading jet for balance
	  iJet = i;
	  jet = jeti;
	  djes = Jet_deltaJES[i];
	  jes = (1.-Jet_rawFactor[i]);
	  res = (1.-Jet_resFactor[i]);
	}
	else { // Subleading jets 
	  jetn += jeti;
	  if (iJet2==-1) { // First subleading jet for alpha
	    iJet2 = i;
	    jet2 = jeti;
	  }
	}
	
	// Calculate L1RC correction
	rawjet = (1-Jet_rawFactor[i]) * jeti;
	double corrl1rc(1.); // isRun3
	if (isRun2) {
	  jecl1rc->setJetPt(rawjet.Pt());
	  jecl1rc->setJetEta(rawjet.Eta());
	  jecl1rc->setJetA(Jet_area[i]);
	  jecl1rc->setRho(fixedGridRhoFastjetAll);
	  corrl1rc = jecl1rc->getCorrection();
	}
	rcjet = corrl1rc * rawjet;
	
	// Corrected type-I chsMET calculation
	corrjets += jeti;
	rawjets += rawjet;
	rcjets += rcjet;
	rcoffsets += (rawjet - rcjet);
      } // non-photon jet
    } // for i in nJet
    
    // Select genjet matching leading and subleading reco jet
    int iGenJet(-1), iGenJet2(-1);
    genjet.SetPtEtaPhiM(0,0,0,0);
    genjet2.SetPtEtaPhiM(0,0,0,0);
    if (isMC) {
      for (Int_t i = 0; i != nGenJet; ++i) {
	geni.SetPtEtaPhiM(GenJet_pt[i],GenJet_eta[i],GenJet_phi[i],
			  GenJet_mass[i]);
	if (iJet!=-1 && geni.DeltaR(jet)<0.4 && iGenJet==-1) {
	  iGenJet = i;
	  genjet = geni;
	}
	else if (iJet2!=-1 && geni.DeltaR(jet2)<0.4 && iGenJet2==-1) {
	  iGenJet2 = i;
	  genjet2 = geni;
	}
      } // for i in nGenJet
    } // isMC

    // Set MET vectors
    if (isRun3) {
      rawmet.SetPtEtaPhiM(RawPuppiMET_pt, 0, RawPuppiMET_phi, 0);
    }
    else {
      rawmet.SetPtEtaPhiM(ChsMET_pt, 0, ChsMET_phi, 0);
    }
    if (isQCD && iFox!=-1) rawmet += fox - gam; // fox=rawjet-PU, gam=genjet
    else rawmet += rawgam - gam; // replace PF photon with Reco photon
    met1 = -jet -gam;
    metn = -jetn;
    //corrmet = rawmet +rawjets -corrjets -rcoffsets;
    corrmet = rawmet +rcjets -corrjets; // same as above
    // Unclustered MET from rawMET by taking out all the hard stuff
    // metu = rawmet +gam +rawjets -rcoffsets;
    // metu = rawmet +gam +rcjets;
    // Or equally well, from corrMET (modulo rounding errors)
    metu = corrmet +gam +corrjets;
    metnu = metn + 1.1*metu;
    met = corrmet;
    
    // Make MET transverse
    corrmet.SetPz(0);
    met.SetPz(0);
    metn.SetPz(0);
    met1.SetPz(0);
    metu.SetPz(0);
    
    // Calculate basic variables
    double ptjet = jet.Pt();
    double abseta = fabs(jet.Eta());
    double pt2 = jet2.Pt();
    double pt2min = 30;
    double bal(0), mpf(0), mpf1(0), mpfn(0), mpfu(0), mpfnu(0);
    double mpfx(0), mpf1x(0), mpfnx(0), mpfux(0), mpfnux(0);
    if (ptgam>0) {
      bal = ptjet / ptgam;
      mpf = 1 + met.Vect().Dot(gam.Vect()) / (ptgam*ptgam);
      mpf1 = 1 + met1.Vect().Dot(gam.Vect()) / (ptgam*ptgam);
      mpfn = metn.Vect().Dot(gam.Vect()) / (ptgam*ptgam);
      mpfu = metu.Vect().Dot(gam.Vect()) / (ptgam*ptgam);
      mpfnu = metnu.Vect().Dot(gam.Vect()) / (ptgam*ptgam);
      //
      gamx.SetPtEtaPhiM(gam.Pt(),gam.Eta(),gam.Phi()+0.5*TMath::Pi(),0.);
      mpfx = 1 + met.Vect().Dot(gamx.Vect()) / (ptgam*ptgam);
      mpf1x = 1 + met1.Vect().Dot(gamx.Vect()) / (ptgam*ptgam);
      mpfnx = metn.Vect().Dot(gamx.Vect()) / (ptgam*ptgam);
      mpfux = metu.Vect().Dot(gamx.Vect()) / (ptgam*ptgam);
      mpfnux = metnu.Vect().Dot(gamx.Vect()) / (ptgam*ptgam);
    }
    
    // Sanity checks for HDM inputs
    if (!(fabs(mpf1+mpfn+mpfu-mpf)<1e-4)) {
      cout << "\nHDM input error: mpf=" << mpf << " mpf1=" << mpf1
	   << " mpfn=" << mpfn << " mpfu=" << mpfu << endl;
      cout << "Difference = " << mpf1+mpfn+mpfu-mpf << endl << flush;
      //assert(false);
      cout << "Skip entry " << jentry
	   << " ("<<run<<","<<luminosityBlock<<","<<event<<")"
	   << " in file " << _filename << endl << flush;
      continue;
    }

    // Event filters for 2016 and 2017+2018 data and MC
    // UL lists are separate, but all filter recommendations looked the same
    // Run3: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Run_3_recommendations
    bool pass_filt = 
      (//(isRun3 && Flag_METFilters>0) ||
       (isRun3 &&
	Flag_goodVertices &&
	Flag_globalSuperTightHalo2016Filter &&
	Flag_EcalDeadCellTriggerPrimitiveFilter &&
	Flag_BadPFMuonFilter &&
	Flag_BadPFMuonDzFilter &&
	Flag_hfNoisyHitsFilter &&
	Flag_eeBadScFilter &&
	Flag_ecalBadCalibFilter) ||
       (!isRun3 &&
	Flag_goodVertices &&
	Flag_globalSuperTightHalo2016Filter &&
	Flag_HBHENoiseFilter &&
	Flag_HBHENoiseIsoFilter &&
	Flag_EcalDeadCellTriggerPrimitiveFilter &&
	Flag_BadPFMuonFilter &&
	//Flag_BadPFMuonDzFilter && // new in UL, but not in nAOD?
	//Flag_BadChargedCandidateFilter && // not recommended
	//Flag_globalTightHalo2016Filter && // obsolete?
	//Flag_CSCTightHaloFilter // obsolete?
	(is16 || Flag_ecalBadCalibFilter) && //new in UL, not for UL16
	//(isMC || Flag_eeBadScFilter) // data only
	Flag_eeBadScFilter // MC added 7 July 2021
	));
    //) || isRun3; // pass_filt
    
    // Photon control plots
    h2ngam->Fill(ptgam, nGam, w);
    if (gengam.Pt()>0 && fabs(gengam.Eta()) < 1.3) {
      hgen->Fill(gengam.Pt(), w);
      peffgr->Fill(gengam.Pt(), iGamGen!=-1 ? 1 : 0, w);
      peffid->Fill(gengam.Pt(), iGam==iGamGen ? 1 : 0, w);
    }
    if (ptgam>0 && fabs(gam.Eta()) < 1.3 && pass_filt) {
      hgam->Fill(ptgam, w);
      if (isMC) pfake->Fill(ptgam, iGam!=iGamGen ? 1 : 0, w);
      if (isQCD) {
	bool hasorig = (iGamOrig!=-1 && gam.DeltaR(gamorig)<0.2);
	bool inwindow = (fabs(gamorig.Pt() / ptgam - 0.9) < 0.2); // [0.8,1.1]
	pfakeqcd->Fill(ptgam, hasorig ? 1 : 0, w);
	pfakeqcd2->Fill(ptgam, hasorig && inwindow ? 1 : 0, w);
	if (hasorig) {
	  h2rgamqcd->Fill(ptgam, gamorig.Pt() / ptgam, w);
	  prgamqcd->Fill(ptgam, gamorig.Pt() / ptgam, w);
	  if (inwindow) prgamqcd2->Fill(ptgam, gamorig.Pt() / ptgam, w);
	  prgamorigqcd->Fill(gamorig.Pt(), ptgam / gamorig.Pt(), w);
	  if (inwindow) prgamorigqcd2->Fill(gamorig.Pt(), ptgam / gamorig.Pt(), w);
	} // hasorig
      }
      if (iGam==iGamGen && gengam.Pt()>0) {
	h2rgam->Fill(gengam.Pt(), ptgam / gengam.Pt(), w);
	prgam->Fill(gengam.Pt(), ptgam / gengam.Pt(), w);
	h2cgam->Fill(ptgam, gengam.Pt() / ptgam, w);
	pcgam->Fill(ptgam, gengam.Pt() / ptgam, w);
      }

      // Plots for photon trigger efficiencies
      if (isMC)  hgam0_mc->Fill(ptgam, w);
      if (!isMC) hgam0_data->Fill(ptgam, w);
      
      hgam0 ->Fill(ptgam, w);
      // Backup high pT
      if (HLT_Photon300_NoHE)  hgam300->Fill(ptgam, w);
      if (HLT_Photon165_HE10)  hgam165h->Fill(ptgam, w);
      if (HLT_Photon100EBHE10) hgam100h->Fill(ptgam, w);
      // Main unprescaled trigger in 2018
      if (HLT_Photon200) hgam200->Fill(ptgam, w);
      // Run 1 style prescaled triggers
      if (HLT_Photon175) hgam175->Fill(ptgam, w);
      if (HLT_Photon150) hgam150->Fill(ptgam, w);
      if (HLT_Photon120) hgam120->Fill(ptgam, w);
      if (HLT_Photon90)  hgam90 ->Fill(ptgam, w);
      if (HLT_Photon75)  hgam75 ->Fill(ptgam, w);
      if (HLT_Photon50)  hgam50 ->Fill(ptgam, w);
      if (HLT_Photon36)  hgam36 ->Fill(ptgam, w);
      if (HLT_Photon33)  hgam33 ->Fill(ptgam, w);
      if (HLT_Photon30)  hgam30 ->Fill(ptgam, w);
      if (HLT_Photon22)  hgam22 ->Fill(ptgam, w);
      if (HLT_Photon20)  hgam20 ->Fill(ptgam, w);
      // 105-230 GeV intermediate range with tight triggers
      if (HLT_Photon120EB_TightID_TightIso) hgam120t->Fill(ptgam, w);
      if (HLT_Photon110EB_TightID_TightIso) hgam110t->Fill(ptgam, w);
      if (HLT_Photon100EB_TightID_TightIso) hgam100t->Fill(ptgam, w);
      if (HLT_Photon30EB_TightID_TightIso) hgam30t->Fill(ptgam, w);
      if (HLT_Photon50EB_TightID_TightIso) hgam50t->Fill(ptgam, w); //new, added hgam50t (27.05.2024)
      // 60-105 GeV with medium triggers. NB: conflicting ID?
      if (HLT_Photon165_R9Id90_HE10_IsoM) hgam165m->Fill(ptgam, w);
      if (HLT_Photon120_R9Id90_HE10_IsoM) hgam120m->Fill(ptgam, w);
      if (HLT_Photon90_R9Id90_HE10_IsoM)  hgam90m ->Fill(ptgam, w);
      if (HLT_Photon75_R9Id90_HE10_IsoM)  hgam75m ->Fill(ptgam, w);
      if (HLT_Photon50_R9Id90_HE10_IsoM)  hgam50m ->Fill(ptgam, w);
      if (HLT_Photon36_R9Id90_HE10_IsoM)  hgam36m ->Fill(ptgam, w);
      if (HLT_Photon30_R9Id90_HE10_IsoM)  hgam30m ->Fill(ptgam, w);
      if (HLT_Photon22_R9Id90_HE10_IsoM)  hgam22m ->Fill(ptgam, w);
      // 20-60 GeV with loose triggers
      if (HLT_Photon60_HoverELoose) hgam60l->Fill(ptgam, w);
      if (HLT_Photon50_HoverELoose) hgam50l->Fill(ptgam, w);
      if (HLT_Photon40_HoverELoose) hgam40l->Fill(ptgam, w);
      if (HLT_Photon30_HoverELoose) hgam30l->Fill(ptgam, w);
      if (HLT_Photon20_HoverELoose) hgam20l->Fill(ptgam, w);
    } // barrel photon

      // Summary of combined trigger efficiencies
      if (ptgam>0 && fabs(gam.Eta())<1.3 && pass_trig && pass_filt) {
	if (isMC)  hgamtrig_mc->Fill(ptgam, w);
	if (!isMC) hgamtrig_data->Fill(ptgam, w);
	hgamtrig->Fill(ptgam, w); // 5 GeV bins to match hgam[trgX]
	hgamtrg->Fill(ptgam, w); // wider binning to higher pT (=hgam)
      }
      //if (ptgam>=105 && fabs(gam.Eta())<1.3 && pass_trig) {
      if (ptgam>=110 && pass_trig && pass_filt) { //why this pt cut??
	//HLT_Photon110EB_TightID_TightIso) {
	h2gametaphi->Fill(gam.Eta(), gam.Phi(), w);
	h2gametaphi2->Fill(gam.Eta(), gam.Phi(), w);
	//h2gametaphi3->Fill(gam.Eta(), gam.Phi(), w);
	//h2gametaphi4->Fill(gam.Eta(), gam.Phi(), w);
      }
      
      bool pass_ngam = (nGam>=1);
      bool pass_njet = (nJets>=1);
      bool pass_gameta = (fabs(gam.Eta()) < 1.3);
      bool pass_dphi = (fabs(gam.DeltaPhi(jet)) > 2.7); // pi-0.44 as in KIT Z+j
      bool pass_jetid = (iJet!=-1 && Jet_jetId[iJet]>=4); // tightLepVeto
      //bool pass_veto = true; //for now: replaced by pass_jetveto and pass_gamveto, since applying jetvetomap also to photons
      bool pass_jetveto = true;
      bool pass_gamveto = true;
      if (true) { // jet veto
        int i1 = h2jv->GetXaxis()->FindBin(jet.Eta());
        int j1 = h2jv->GetYaxis()->FindBin(jet.Phi());
        //if (h2jv->GetBinContent(i1,j1)>0) { //use this if no bpix veto applied (in 2023: handle bpix veto differently)
			if(is24){
        if(h2jv->GetBinContent(i1,j1)>0 or bpixjv->GetBinContent(i1,j1)>0) { //use this when also vetoing bpix (2024 onwards)
          //++_nbadevents_veto;
	  //pass_veto = false;
          pass_jetveto = false;
				}
			} else{ //if not 2024
					if (h2jv->GetBinContent(i1,j1)>0) {
						pass_jetveto = false;
					}
				}
      } // jet veto
      if (true) { // photon veto
        int i1 = h2jv->GetXaxis()->FindBin(gam.Eta());
        int j1 = h2jv->GetYaxis()->FindBin(gam.Phi());
        //if (h2jv->GetBinContent(i1,j1)>0) { //use this if no bpix veto applied
				if(is24){
        	if(h2jv->GetBinContent(i1,j1)>0 or bpixjv->GetBinContent(i1,j1)>0) { //use this when also vetoing bpix
          	//++_nbadevents_veto;
	  				//pass_veto = false;
          	pass_gamveto = false; //reject also photons that end up in bad regions
					}
				}else{//if not 2024
					if (h2jv->GetBinContent(i1,j1)>0) {
						pass_gamveto = false;
					}
				}
      } // photon veto
      /* NOTE: REMOVE ALSO PHOTONS ACCORDING TO JETVETOMPA pass_veto_gam = false --> add it also to pass_basic*/
      if(pass_jetveto==false || pass_gamveto==false){//increase counter if event gets discarded either due to jet-veto or photon-veto
        ++_nbadevents_veto;
      }
      bool pass_leak = (phoj.Pt()<0.06*ptgam);// || isRun3);
      bool pass_basic = (pass_trig && pass_filt && pass_ngam && pass_njet &&
			 pass_gameta && pass_dphi && pass_jetid && 
                         pass_jetveto && pass_gamveto && //pass_veto &&
			 pass_leak); // add pass_gameta v19 / 202111122 !
      
      // Control plots for jet response
      bool pass_bal = (fabs(1-bal)<0.7);
      bool pass_mpf = (fabs(1-mpf)<0.7);
      bool pass_jeteta = (abseta < 1.3);
      bool pass_alpha100 = (pt2 < ptgam || pt2 < pt2min);      
      bool pass_basic_ext = (pass_basic && pass_bal && pass_mpf);
      bool pass_all = (pass_basic_ext && pass_jeteta && pass_alpha100);
      bool pass_gen = (iGenJet!=-1);

      //const int nc = 18;
      bool cut[nc] = {pass_trig, pass_ngam, pass_njet, pass_gameta,
		      //pass_dphi, pass_jetid, pass_veto, pass_leak,
		      pass_dphi, pass_jetid, pass_jetveto, pass_gamveto, pass_leak,
		      pass_basic, pass_bal, pass_mpf, pass_basic_ext,
		      pass_jeteta, pass_alpha100, pass_all, pass_gen};
      bool passcuts(true);
      for (int i = 0; i != nc; ++i) {
	pcutflow1->Fill(i, cut[i] ? 1 : 0, w);
	if (passcuts) {
	  pcutflow2->Fill(i, cut[i] ? 1 : 0, w);
	  passcuts = cut[i];
	}
	pcutflow3->Fill(i, passcuts ? 1 : 0, w);
      }
      if (pass_trig && pass_ngam && pass_njet && pass_gameta) {
	hdphi->Fill(gam.DeltaPhi(jet), w);
	hdr->Fill(gam.DeltaR(jet), w);
      }

      // 1) Time controls for JES and PF composition; 2) etaphi maps (for jet veto) with response over eta and phi; 3) etaphi maps for rate (n jets)
      if (pass_all) {
	if (itrg==30 && ptgam>32) { //check turn-on curves, changed offline cut to ptgam>32 (used to be at 30gev)
	  pr30n->Fill(run, w); 
	  // TO DO: make this pr30xs plot (below) without jet requirements (photon passes trigger + has required pT)
          pr30xs->Fill(run, lumi30[run] ? 1./lumi30[run] : 1.); //new, if lumi calculated for that run number, normalise, if not then just use weight=1.0
	  pr30b->Fill(run, bal, w); 
	  pr30m->Fill(run, mpf, w);
	  pr30chf->Fill(run, Jet_chHEF[iJet], w);
	  pr30nhf->Fill(run, Jet_neHEF[iJet], w);
	  pr30nef->Fill(run, Jet_neEmEF[iJet], w);
	}
        if (itrg==50 && ptgam>53) { //ptgam>53 (to avoid trouble with hlt scale) (used to be ptgam>50)
	  pr50n->Fill(run, w); 
          pr50xs->Fill(run, lumi50[run] ? 1./lumi50[run] : 1.); //new (can remove this when pr50n one above shows xs)
	  pr50b->Fill(run, bal, w); 
	  pr50m->Fill(run, mpf, w);
	  pr50chf->Fill(run, Jet_chHEF[iJet], w);
	  pr50nhf->Fill(run, Jet_neHEF[iJet], w);
	  pr50nef->Fill(run, Jet_neEmEF[iJet], w);

		//etaphi maps:
		p2bal50_jetetaphi->Fill(jet.Eta(), jet.Phi(), bal, w);
		//h2mpf50_jetetaphi->Fill(gam.Eta(), gam.Phi(), mpf, w);
		h2n50_jetetaphi->Fill(jet.Eta(), jet.Phi(), w); //event rate 
		h2n50_gametaphi->Fill(gam.Eta(), gam.Phi(), w); //event rate (photon eta, photon phi)

		//jet1 eta distribution for different jet-pt
		if(jet.Pt() > 30){h50n_j1eta_j1pt30->Fill(jet.Eta())};
		if(jet.Pt() > 40){h50n_j1eta_j1pt40->Fill(jet.Eta())};
		if(jet.Pt() > 50){h50n_j1eta_j1pt50->Fill(jet.Eta())};
		
	}
	if (itrg==110 && ptgam>120) { //offline cut ptgam > 120 (used to be 110)
	  pr110n->Fill(run, w);
          pr110xs->Fill(run, lumi110[run] ? 1./lumi110[run] : 1.); //new
	  pr110b->Fill(run, bal, w); 
	  pr110m->Fill(run, mpf, w);
	  pr110chf->Fill(run, Jet_chHEF[iJet], w);
	  pr110nhf->Fill(run, Jet_neHEF[iJet], w);
	  pr110nef->Fill(run, Jet_neEmEF[iJet], w);

		//etaphi maps:
		p2bal110_jetetaphi->Fill(jet.Eta(), jet.Phi(), bal, w);
		h2n110_jetetaphi->Fill(jet.Eta(), jet.Phi(), w); //event rate 
		h2n110_gametaphi->Fill(gam.Eta(), gam.Phi(), w); //event rate (photon eta, photon phi)

		//jet1 eta distribution for different jet-pt
		if(jet.Pt() > 30){h110n_j1eta_j1pt30->Fill(jet.Eta())};
		if(jet.Pt() > 40){h110n_j1eta_j1pt40->Fill(jet.Eta())};
		if(jet.Pt() > 50){h110n_j1eta_j1pt50->Fill(jet.Eta())};

	}
	if (itrg==200 && ptgam>230) {
	  pr230n->Fill(run, w);
          pr230xs->Fill(run, lumi200[run] ? 1./lumi200[run] : 1.); //new
	  pr230b->Fill(run, bal, w); 
	  pr230m->Fill(run, mpf, w);
	  pr230chf->Fill(run, Jet_chHEF[iJet], w);
	  pr230nhf->Fill(run, Jet_neHEF[iJet], w);
	  pr230nef->Fill(run, Jet_neEmEF[iJet], w);

		//etaphi maps:
		p2bal200_jetetaphi->Fill(jet.Eta(), jet.Phi(), bal, w);
		h2n200_jetetaphi->Fill(jet.Eta(), jet.Phi(), w); //event rate 
		h2n200_gametaphi->Fill(gam.Eta(), gam.Phi(), w); //event rate (photon eta, photon phi)

		//jet1 eta distribution for different jet-pt
		if(jet.Pt() > 30){h200n_j1eta_j1pt30->Fill(jet.Eta())};
		if(jet.Pt() > 40){h200n_j1eta_j1pt40->Fill(jet.Eta())};
		if(jet.Pt() > 50){h200n_j1eta_j1pt50->Fill(jet.Eta())};

	}
	if (iGam!=-1 && Photon_seedGain[iGam]==1) {
	  prg1n->Fill(run, w);
	  prg1b->Fill(run, bal, w); 
	  prg1m->Fill(run, mpf, w);
	  prg1chf->Fill(run, Jet_chHEF[iJet], w);
	  prg1nhf->Fill(run, Jet_neHEF[iJet], w);
	  prg1nef->Fill(run, Jet_neEmEF[iJet], w);
	}
      }

      //for high eta
      // Time controls for JES and PF composition in HIGH ETA range
      //NEW: (jet eta) 3<=eta<4 and TO DO:  4<=eta<5
      //do this eta investigation for all histograms, might not be so useful for composition plots, but check
      if (pass_basic_ext and pass_alpha100 and fabs(Jet_eta[iJet])>=3.0 and fabs(Jet_eta[iJet])<4.0) {
	if (itrg==30 && ptgam>32) {
	  pr30n_eta3to4->Fill(run, w); 
          pr30xs_eta3to4->Fill(run, lumi30[run] ? 1./lumi30[run] : 1.); //new, if lumi calculated for that run number, normalise, if not then just use weight=1.0
	  pr30b_eta3to4->Fill(run, bal, w); 
	  pr30m_eta3to4->Fill(run, mpf, w);
	  //pr30chf->Fill(run, Jet_chHEF[iJet], w);
	  //pr30nhf->Fill(run, Jet_neHEF[iJet], w);
	  //pr30nef->Fill(run, Jet_neEmEF[iJet], w);
					if(is24){ //HF branches only since 2024
          		pr30hfEmEF_eta3to4->Fill(run, Jet_hfEmEF[iJet], w);
          		pr30hfHEF_eta3to4->Fill(run, Jet_hfHEF[iJet], w);
					}
	}
        if (itrg==50 && ptgam>53) {
          pr50n_eta3to4->Fill(run, w);
          pr50xs_eta3to4->Fill(run, lumi50[run] ? 1./lumi50[run] : 1.);
          pr50b_eta3to4->Fill(run, bal, w);
          pr50m_eta3to4->Fill(run, mpf, w);
	  //pr50chf->Fill(run, Jet_chHEF[iJet], w);
	  //pr50nhf->Fill(run, Jet_neHEF[iJet], w);
	  //pr50nef->Fill(run, Jet_neEmEF[iJet], w);
					if(is24){ //HF branches only since 2024
							pr50hfEmEF_eta3to4->Fill(run, Jet_hfEmEF[iJet], w);
							pr50hfHEF_eta3to4->Fill(run, Jet_hfHEF[iJet], w);
					}
	}
	if (itrg==110 && ptgam>120) {
	  pr110n_eta3to4->Fill(run, w);
          pr110xs_eta3to4->Fill(run, lumi110[run] ? 1./lumi110[run] : 1.); //new
	  pr110b_eta3to4->Fill(run, bal, w); 
	  pr110m_eta3to4->Fill(run, mpf, w);
	  //pr110chf->Fill(run, Jet_chHEF[iJet], w);
	  //pr110nhf->Fill(run, Jet_neHEF[iJet], w);
	  //pr110nef->Fill(run, Jet_neEmEF[iJet], w);
					if(is24){ //HF branches only since 2024
						pr110hfEmEF_eta3to4->Fill(run, Jet_hfEmEF[iJet], w);
						pr110hfHEF_eta3to4->Fill(run, Jet_hfHEF[iJet], w);
				}
	}
        /*
	if (itrg==200 && ptgam>230) {
	  pr230n->Fill(run, w);
          pr230xs->Fill(run, lumi200[run] ? 1./lumi200[run] : 1.); //new
	  pr230b->Fill(run, bal, w); 
	  pr230m->Fill(run, mpf, w);
	  pr230chf->Fill(run, Jet_chHEF[iJet], w);
	  pr230nhf->Fill(run, Jet_neHEF[iJet], w);
	  pr230nef->Fill(run, Jet_neEmEF[iJet], w);
	}
	if (iGam!=-1 && Photon_seedGain[iGam]==1) {
	  prg1n->Fill(run, w);
	  prg1b->Fill(run, bal, w); 
	  prg1m->Fill(run, mpf, w);
	  prg1chf->Fill(run, Jet_chHEF[iJet], w);
	  prg1nhf->Fill(run, Jet_neHEF[iJet], w);
	  prg1nef->Fill(run, Jet_neEmEF[iJet], w);
	}
        */
      }

      //JET ETA between 4.0 and 5.0
      if (pass_basic_ext and pass_alpha100 and fabs(Jet_eta[iJet])>=4.0 and fabs(Jet_eta[iJet])<5.0) {
	if (itrg==30 && ptgam>32) {
	  pr30n_eta4to5->Fill(run, w); 
          pr30xs_eta4to5->Fill(run, lumi30[run] ? 1./lumi30[run] : 1.); //new, if lumi calculated for that run number, normalise, if not then just use weight=1.0
	  pr30b_eta4to5->Fill(run, bal, w); 
	  pr30m_eta4to5->Fill(run, mpf, w);
	  //pr30chf->Fill(run, Jet_chHEF[iJet], w);
	  //pr30nhf->Fill(run, Jet_neHEF[iJet], w);
	  //pr30nef->Fill(run, Jet_neEmEF[iJet], w);
					if(is24){ //HF branches only since 2024
          		pr30hfEmEF_eta4to5->Fill(run, Jet_hfEmEF[iJet], w);
          		pr30hfHEF_eta4to5->Fill(run, Jet_hfHEF[iJet], w);
					}
	}
        if (itrg==50 && ptgam>53) {
          pr50n_eta4to5->Fill(run, w);
          pr50xs_eta4to5->Fill(run, lumi50[run] ? 1./lumi50[run] : 1.);
          pr50b_eta4to5->Fill(run, bal, w);
          pr50m_eta4to5->Fill(run, mpf, w);
	  //pr50chf->Fill(run, Jet_chHEF[iJet], w);
	  //pr50nhf->Fill(run, Jet_neHEF[iJet], w);
	  //pr50nef->Fill(run, Jet_neEmEF[iJet], w);
					if(is24){ //HF branches only since 2024
          		pr50hfEmEF_eta4to5->Fill(run, Jet_hfEmEF[iJet], w);
          		pr50hfHEF_eta4to5->Fill(run, Jet_hfHEF[iJet], w);
					}
	}
	if (itrg==110 && ptgam>120) {
	  pr110n_eta4to5->Fill(run, w);
          pr110xs_eta4to5->Fill(run, lumi110[run] ? 1./lumi110[run] : 1.); //new
	  pr110b_eta4to5->Fill(run, bal, w); 
	  pr110m_eta4to5->Fill(run, mpf, w);
	  //pr110chf->Fill(run, Jet_chHEF[iJet], w);
	  //pr110nhf->Fill(run, Jet_neHEF[iJet], w);
	  //pr110nef->Fill(run, Jet_neEmEF[iJet], w);
					if(is24){ //HF branches only since 2024
          		pr110hfEmEF_eta4to5->Fill(run, Jet_hfEmEF[iJet], w);
          		pr110hfHEF_eta4to5->Fill(run, Jet_hfHEF[iJet], w);
					}
	}
        /*
	if (itrg==200 && ptgam>230) {
	  pr230n->Fill(run, w);
          pr230xs->Fill(run, lumi200[run] ? 1./lumi200[run] : 1.); //new
	  pr230b->Fill(run, bal, w); 
	  pr230m->Fill(run, mpf, w);
	  pr230chf->Fill(run, Jet_chHEF[iJet], w);
	  pr230nhf->Fill(run, Jet_neHEF[iJet], w);
	  pr230nef->Fill(run, Jet_neEmEF[iJet], w);
	}
	if (iGam!=-1 && Photon_seedGain[iGam]==1) {
	  prg1n->Fill(run, w);
	  prg1b->Fill(run, bal, w); 
	  prg1m->Fill(run, mpf, w);
	  prg1chf->Fill(run, Jet_chHEF[iJet], w);
	  prg1nhf->Fill(run, Jet_neHEF[iJet], w);
	  prg1nef->Fill(run, Jet_neEmEF[iJet], w);
	}
        */
      }







      if (_gh_debug100 && jentry<100) {
	cout << "Event " << jentry << " decisions" << endl;
	cout << Form("pass_ngam = %d, pass_njet = %d, pass_gameta = %d "
		     "pass_dphi = %d, pass_jetid = %d\n"
		     //"pass_veto = %d, pass_leak = %d, pass_basic = %d "
                     "pass_jetveto = %d, pass_gamveto = %d, pass_leak = %d, pass_basic = %d "
		     "pass_bal = %d, pass_mpf = %d, \n"
		     "pass_jeteta = %d, pass_alpha100 = %d, "
		     "pass_basic_ext = %d, "
		     "pass_gen = %d,\n"
		     "pass_trig = %d, pass_filt = %d",
		     pass_ngam, pass_njet, pass_gameta,
		     pass_dphi, pass_jetid,
		     //pass_veto, pass_leak, pass_basic,
                     pass_jetveto, pass_gamveto, pass_leak, pass_basic,
		     pass_bal, pass_mpf, pass_jeteta,
		     pass_alpha100, pass_basic_ext,
		     pass_gen, pass_trig, pass_filt) << endl << flush;
	cout << Form("Flags/Filters: goodVertices = %d, "
		     "globalSuperTightHalo2016 = %d,\n"
		     "HBHENoise = %d, "
		     "HBHENoiseIso = %d, "
		     "EcalDeadCellTriggerPrimitive = %d,\n"
		     "BadPFMuon = %d, "
		     "ecalBadCalib = %d, "
		     "eeBadSc = %d.\n",  
		     Flag_goodVertices,
		     Flag_globalSuperTightHalo2016Filter,
		     Flag_HBHENoiseFilter,
		     Flag_HBHENoiseIsoFilter,
		     Flag_EcalDeadCellTriggerPrimitiveFilter,
		     Flag_BadPFMuonFilter,
		     Flag_ecalBadCalibFilter,
		     Flag_eeBadScFilter) << endl << flush;
      }

      if (pass_basic && pass_jeteta && pass_alpha100) {
	h2bal->Fill(ptgam, bal, w);
	h2mpf->Fill(ptgam, mpf, w);
	if (pass_mpf) h2balc->Fill(ptgam, bal, w);
	if (pass_bal) h2mpfc->Fill(ptgam, mpf, w);
	if (pass_mpf && pass_bal) h2balc2->Fill(ptgam, bal, w);
	if (pass_mpf && pass_bal) h2mpfc2->Fill(ptgam, mpf, w);
	if (pass_basic_ext) {

	  if (pass_gen || !isMC) {
	    h2rjet->Fill(ptgam, jet.Pt() / ptgam, w);
	    prjet->Fill(ptgam, jet.Pt() / ptgam, w);
	  }
	  if (pass_gen) {
	    h2gjet->Fill(ptgam, genjet.Pt() / ptgam, w);
	    pgjet->Fill(ptgam, genjet.Pt() / ptgam, w);
	    h2rgen->Fill(genjet.Pt(), jet.Pt() / genjet.Pt(), w);
	    prgen->Fill(genjet.Pt(), jet.Pt() / genjet.Pt(), w);
	  }

	  int flv = (isMC ? GenJet_partonFlavour[iGenJet] : 99);
	  mvar["counts"] = 1;
	  mvar["mpfchs1"] = mpf;
	  mvar["ptchs"] = bal;
	  mvar["mpf1"] = mpf1;
	  mvar["mpfn"] = mpfn;
	  mvar["mpfu"] = mpfu;
	  mvar["rho"] = fixedGridRhoFastjetAll;
	  mvar["rjet"] = (ptgam!=0 ? jet.Pt() / ptgam : 0);
	  mvar["gjet"] = (ptgam!=0 ? genjet.Pt() / ptgam : 0);
	  mvar["rgen"] = (genjet.Pt()!=0 ? jet.Pt() / genjet.Pt() : 0);

	  if (isRun3) { // temporary patch
	    Jet_btagDeepB[iJet] = Jet_btagDeepFlavB[iJet];
	    Jet_btagDeepC[iJet] = 0.5*(Jet_btagDeepFlavCvB[iJet] +
				       Jet_btagDeepFlavCvL[iJet]);
	    Jet_qgl[iJet] = Jet_btagDeepFlavQG[iJet];
	  }
	  bool isb = (Jet_btagDeepB[iJet] > bthr);
	  bool isc = (Jet_btagDeepC[iJet] > cthr && !isb);
	  bool isq = (Jet_qgl[iJet]>=0.5 && Jet_qgl[iJet] && !isb && !isc);
	  bool isg = (Jet_qgl[iJet]>=0 && Jet_qgl[iJet]<0.5 && !isb && !isc);
	  bool isn = (!isb && !isc && !isq && !isg);
	  
	  for (int ivar = 0; ivar != nvar; ++ivar) {
	    for (int itag = 0; itag != ntag; ++itag) {
	      for (int iflv = 0; iflv != nflv; ++iflv) {

		string& svr = avar[ivar]; const char *cv = svr.c_str();
		string& stg = atag[itag]; const char *ct = stg.c_str();
		string& sfl = aflv[iflv]; const char *cf = sfl.c_str();

		if (((sfl=="i") ||
		     (sfl=="b" && abs(flv)==5) ||
		     (sfl=="c" && abs(flv)==4) ||
		     (sfl=="q" && abs(flv)<=3 && flv!=0) ||
		     (sfl=="s" && abs(flv)==3) ||
		     (sfl=="ud" && abs(flv)<=2 && flv!=0) ||
		     (sfl=="g" && flv==21) ||
		     (sfl=="n" && flv==0)) &&
		    ((stg=="i") ||
		     (stg=="b" && isb) ||
		     (stg=="c" && isc) ||
		     (stg=="q" && isq) ||
		     (stg=="g" && isg) ||
		     (stg=="n" && isn))) {

		  double x = ptgam;
		  if (svr=="rgen") x = genjet.Pt();
		  if ((svr=="rjet" || svr=="gjet") && iGenJet==-1) x = 0;
		  double var = mvar[svr];
		  TH1* h = mp[svr][stg][sfl];
		  if (!h) {
		    cout << "Missing "<<svr<<stg<<sfl<<endl<<flush;
		    assert(h);
		  }

		  if (svr=="counts")
		    ((TH1D*)h)->Fill(x, w);
		  else
		    ((TProfile*)h)->Fill(x, var, w);
		}
	      } // for iflv
	    } // for itag
	  } // for ivar

	  /*
	  double genpt = genjet.Pt();
	  double rgen = jet.Pt() / genjet.Pt();
	  double ggam = genjet.Pt() / gam.Pt();
	  if (abs(flv)==1) prgenud->Fill(genpt, rgen, w);
	  if (abs(flv)==2) prgenud->Fill(genpt, rgen, w);
	  if (abs(flv)==3) prgens->Fill(genpt, rgen, w);
	  if (abs(flv)==4) prgenc->Fill(genpt, rgen, w);
	  if (abs(flv)==5) prgenb->Fill(genpt, rgen, w);
	  if (flv==21)     prgeng->Fill(genpt, rgen, w);
	  if (flv==0)      prgeno->Fill(genpt, rgen, w);
	  if (abs(flv)==1) pgjetud->Fill(ptgam, ggam, w);
	  if (abs(flv)==2) pgjetud->Fill(ptgam, ggam, w);
	  if (abs(flv)==3) pgjets->Fill(ptgam, ggam, w);
	  if (abs(flv)==4) pgjetc->Fill(ptgam, ggam, w);
	  if (abs(flv)==5) pgjetb->Fill(ptgam, ggam, w);
	  if (flv==21)     pgjetg->Fill(ptgam, ggam, w);
	  if (flv==0)      pgjeto->Fill(ptgam, ggam, w);
	  pfud->Fill(ptgam, abs(flv)==1||abs(flv)==2 ? 1 : 0, w);
	  pfs->Fill(ptgam, abs(flv)==3 ? 1 : 0, w);
	  pfc->Fill(ptgam, abs(flv)==4 ? 1 : 0, w);
	  pfb->Fill(ptgam, abs(flv)==5 ? 1 : 0, w);
	  pfg->Fill(ptgam, flv==21 ? 1 : 0, w);
	  pfo->Fill(ptgam, flv==0 ? 1 : 0, w);
          */
	}
	if (pass_basic_ext && jet2.Pt()>0) {
	  if (iGenJet2!=-1 || !isMC) {
	    h2rjet2->Fill(ptgam, jet2.Pt() / ptgam, w);
	    prjet2->Fill(ptgam, jet2.Pt() / ptgam, w);
	  }
	  if (iGenJet2!=-1) {
	    h2gjet2->Fill(ptgam, genjet2.Pt() / ptgam, w);
	    pgjet2->Fill(ptgam, genjet2.Pt() / ptgam, w);
	    h2rgen2->Fill(genjet2.Pt(), jet2.Pt() / genjet2.Pt(), w);
	    prgen2->Fill(genjet2.Pt(), jet2.Pt() / genjet2.Pt(), w);
	  }
	}
      }

      // Control plots for trigger 
      if (ptgam>0 && pass_basic_ext && pass_jeteta && pass_alpha100) {
	//     optimal trigger edges: (20,30,(35),55,80,95,105,115,210)
	//     old bin trigger edges  (20,30,60,85,*95*,105,130,230)
	double pt = ptgam;
	double mu = Pileup_nTrueInt;
	if (isMC                             && pt>210)  hmusmc->Fill(mu, w);
	int nmax = (isMC ? 1 : 100);
	for (int i=0; i!=nmax; ++i) {
	  mu = gRandom->Gaus(Pileup_nTrueInt,TruePUrms);
	  double w1 = 0.01*w;
	  if (HLT_Photon20_HoverELoose         && pt>=20)  hmus20->Fill(mu,w1);
	  if (HLT_Photon30_HoverELoose         && pt>=30)  hmus30->Fill(mu,w1);
	  if (HLT_Photon50_R9Id90_HE10_IsoM    && pt>=55)  hmus50->Fill(mu,w1);
          if (HLT_Photon50EB_TightID_TightIso  && pt>=55)  hmus50->Fill(mu,w1);
	  if (HLT_Photon75_R9Id90_HE10_IsoM    && pt>=80)  hmus75->Fill(mu,w1);
          if (HLT_Photon75EB_TightID_TightIso  && pt>=55)  hmus75->Fill(mu,w1);
	  if (HLT_Photon90_R9Id90_HE10_IsoM    && pt>=95)  hmus90->Fill(mu,w1);
          if (HLT_Photon90EB_TightID_TightIso  && pt>=55)  hmus90->Fill(mu,w1);
	  if (HLT_Photon100EB_TightID_TightIso && pt>=105) hmus100->Fill(mu,w1);
	  if (HLT_Photon110EB_TightID_TightIso && pt>=115) hmus110->Fill(mu,w1);
	  if (HLT_Photon200                    && pt>=210) hmus200->Fill(mu,w1);
	} // for i in nmax
      } // control plots for triggers

      // Control plots for photon-jet 
      if (ptgam>0 && pass_basic_ext) {

	h2phoj->Fill(ptgam, phoj.Pt()/ptgam, w);
	pphoj->Fill(ptgam, phoj.Pt()/ptgam, w);

	double footprint =
	  (fabs(gam.DeltaPhi(phoj0))<0.4 ? 1 : -1)*phoj0.Pt()/ptgam;
	h2phoj0->Fill(ptgam, footprint, w);
	pphoj0->Fill(ptgam, footprint, w);

	if (iGam!=-1 && Photon_seedGain[iGam]==1) {
	  h2phoj1->Fill(ptgam, footprint, w);
	  pphoj1->Fill(ptgam, footprint, w);
	  if (pass_jeteta && pass_alpha100) {
	    prbal1->Fill(ptgam, bal, w);
	    prmpf1->Fill(ptgam, mpf, w);
	    pr2bal1->Fill(ptgam, gam.Eta(), bal, w); //new in w27+w28
	    pr2mpf1->Fill(ptgam, gam.Eta(), mpf, w);
	  }
	}
	if (iGam!=-1 && Photon_seedGain[iGam]==6) {
	  h2phoj6->Fill(ptgam, footprint, w);
	  pphoj6->Fill(ptgam, footprint, w);
	  if (pass_jeteta && pass_alpha100) {
	    prbal6->Fill(ptgam, bal, w);
	    prmpf6->Fill(ptgam, mpf, w);
	    pr2bal6->Fill(ptgam, gam.Eta(), bal, w);
	    pr2mpf6->Fill(ptgam, gam.Eta(), mpf, w);
	  }
	}
	if (iGam!=-1 && Photon_seedGain[iGam]==12) {
	  h2phoj12->Fill(ptgam, footprint, w);
	  pphoj12->Fill(ptgam, footprint, w);
	  if (pass_jeteta && pass_alpha100) {
	    prbal12->Fill(ptgam, bal, w);
	    prmpf12->Fill(ptgam, mpf, w);
	    pr2bal12->Fill(ptgam, gam.Eta(), bal, w);
	    pr2mpf12->Fill(ptgam, gam.Eta(), mpf, w);
	  }
	}

	if (pass_jeteta && pass_alpha100) {
	  prbal->Fill(ptgam, bal, w);
	  prmpf->Fill(ptgam, mpf, w);
	  pr2bal->Fill(ptgam, gam.Eta(), bal, w);
	  pr2mpf->Fill(ptgam, gam.Eta(), mpf, w);
	  prbal0->Fill(ptgam, bal, w);
	  prmpf0->Fill(ptgam, mpf, w);
	  
	  // PF composition plots
	  h2pteta->Fill(ptgam, Jet_eta[iJet], w);
	  pabseta->Fill(ptgam, fabs(Jet_eta[iJet]), w);

	  // 1D composition and response
	  pdb->Fill(ptgam, bal, w);
	  pmpf->Fill(ptgam, mpf, w);
	  pchf->Fill(ptgam, Jet_chHEF[iJet], w);
	  pnhf->Fill(ptgam, Jet_neHEF[iJet], w);
	  pnef->Fill(ptgam, Jet_neEmEF[iJet], w);
	  pcef->Fill(ptgam, Jet_chEmEF[iJet], w);
	  pmuf->Fill(ptgam, Jet_muEF[iJet], w);
	  //if (isRun3) Jet_chFPV0EF[iJet] = 0;
	  //ppuf->Fill(ptgam, Jet_chFPV0EF[iJet], w);
	  
	  // 2D composition and response
	  if (ptgam>230) {
	    double eta = Jet_eta[iJet];
	    double phi = Jet_phi[iJet];
	    p2db->Fill(eta, phi, bal, w);
	    p2mpf->Fill(eta, phi, mpf, w);
	    p2chf->Fill(eta, phi, Jet_chHEF[iJet], w);
	    p2nhf->Fill(eta, phi, Jet_neHEF[iJet], w);
	    p2nef->Fill(eta, phi, Jet_neEmEF[iJet], w);
	    p2cef->Fill(eta, phi, Jet_chEmEF[iJet], w);
	    p2muf->Fill(eta, phi, Jet_muEF[iJet], w);
	    //p2puf->Fill(eta, phi, Jet_chFPV0EF[iJet], w);
	  }

	  if (isMC) {
	    if (ptgam>=105 && ptgam<230)
	      hmus->Fill(Pileup_nTrueInt, w);
	    h2mus->Fill(ptgam, Pileup_nTrueInt, w);
	  }
	  else {
	    for (int i=0; i!=100; ++i) {
	      double mu = gRandom->Gaus(Pileup_nTrueInt,TruePUrms);
	      if (ptgam>=105 && ptgam<230)
		hmus->Fill(mu, 0.01*w);
	      h2mus->Fill(ptgam, mu, 0.01*w);
	    } // for i in 100
	  } // is MC
	  //if (ptgam>=130 && ptgam<175) {
	  //if ((is16 && ptgam>175) ||
	  //  (is17 && ptgam>230) ||
	  //  (is18 && ptgam>130)) {
	  if (ptgam>230 && iGam!=-1) {
	    pgainvsmu->Fill(Pileup_nTrueInt, Photon_seedGain[iGam], w);
	    if (b_Photon_eCorr) // safety for 2016
	      pcorrvsmu->Fill(Pileup_nTrueInt, Photon_eCorr[iGam], w);
	    perrvsmu->Fill(Pileup_nTrueInt, Photon_energyErr[iGam], w);
	    phoevsmu->Fill(Pileup_nTrueInt, Photon_hoe[iGam], w);
	    pr9vsmu->Fill(Pileup_nTrueInt, Photon_r9[iGam], w);
	  }
	  if (ptgam>230) {
	    pmuvsmu->Fill(Pileup_nTrueInt, Pileup_nTrueInt, w);
	    prhovsmu->Fill(Pileup_nTrueInt, fixedGridRhoFastjetAll, w);
	    pnpvgoodvsmu->Fill(Pileup_nTrueInt, PV_npvsGood, w);
	    pnpvallvsmu->Fill(Pileup_nTrueInt, PV_npvs, w);
	  } // high pT range
	  
	  if (iGam!=-1) {
	    pgain1vspt->Fill(ptgam, Photon_seedGain[iGam]==1 ? 1 : 0, w);
	    pgain6vspt->Fill(ptgam, Photon_seedGain[iGam]==6 ? 1 : 0, w);
	    pgain12vspt->Fill(ptgam, Photon_seedGain[iGam]==12 ? 1 : 0, w);
	    pgainvspt->Fill(ptgam, Photon_seedGain[iGam], w);

			//2D profiles (vs pt, vs eta)
	    pgain1vsptvseta->Fill(ptgam, gam.Eta(), Photon_seedGain[iGam]==1 ? 1 : 0, w);
	    pgain6vsptvseta->Fill(ptgam, gam.Eta(), Photon_seedGain[iGam]==6 ? 1 : 0, w);
	    pgain12vsptvseta->Fill(ptgam, gam.Eta(), Photon_seedGain[iGam]==12 ? 1 : 0, w);
	    pgainvsptvseta->Fill(ptgam, gam.Eta(), Photon_seedGain[iGam], w);

	    if (b_Photon_eCorr) // safety for 2016
	      pcorrvspt->Fill(ptgam, Photon_eCorr[iGam], w);
	    perrvspt->Fill(ptgam, Photon_energyErr[iGam], w);
	    h2hoevspt->Fill(ptgam, Photon_hoe[iGam], w);
	    phoevspt->Fill(ptgam, Photon_hoe[iGam], w);
	    h2r9vspt->Fill(ptgam, Photon_r9[iGam], w);
	    pr9vspt->Fill(ptgam, Photon_r9[iGam], w);
	  }
	  pmuvspt->Fill(ptgam, Pileup_nTrueInt, w);
	  prhovspt->Fill(ptgam, fixedGridRhoFastjetAll, w);
	  pnpvgoodvspt->Fill(ptgam, PV_npvsGood, w);
	  pnpvallvspt->Fill(ptgam, PV_npvs, w);
	} // barrel
      } // basic_ext cuts
      
      // Specific event selection for alpha and eta bins
      for (unsigned int ialpha = 0; ialpha != alphas.size(); ++ialpha) {
	for (unsigned int ieta = 0; ieta != etas.size(); ++ieta) { 
	for (int ips = 0; ips != nPSWeight+1; ++ips) { 

	  double alpha = alphas[ialpha];
	  double ymin = etas[ieta].first;
	  double ymax = etas[ieta].second;
	  double wps = (ips==0 ? 1. : PSWeight[ips-1]);
	  
	  // Specific event selection
	  bool pass_alpha = (pt2 < alpha*ptgam || pt2 < pt2min);
	  bool pass_eta = (abseta >= ymin && abseta < ymax);
	  bool pass = (pass_basic_ext && pass_alpha && pass_eta);
	  
	  if (pass) {
	    
	    // Retrieve pointers to correct set of histograms/profiles
	    int ia = int(100*alpha);
	    int iy = 100*int(10*ymin) + int(10*ymax);
	    
	    // Get reference instead of pointer so can use . and not ->
	    BasicHistos *pmh = mBasicHistos[iy][ia][ips]; assert(pmh);
	    BasicHistos& mh = (*pmh); assert(mh.hn);
	    
	    // Fill histograms (h*) and profiles (p*)
	    //assert(fabs(mpf1+mpfn+mpfu-mpf)<1e-4);
	    mh.hn->Fill(ptgam);
	    mh.hxsec->Fill(ptgam, w*wps);
	    mh.prpt->Fill(ptgam, ptgam, w*wps);
	    mh.prbal->Fill(ptgam, bal, w*wps);
	    mh.prdb->Fill(ptgam, mpf1, w*wps);
	    mh.prmpf->Fill(ptgam, mpf, w*wps);
	    mh.prmpf1->Fill(ptgam, mpf1, w*wps);
	    mh.prmpfn->Fill(ptgam, mpfn, w*wps);
	    mh.prmpfu->Fill(ptgam, mpfu, w*wps);
	    mh.prho->Fill(ptgam, fixedGridRhoFastjetAll, w);
	    mh.pdjes->Fill(ptgam, djes, w);
	    mh.pjes->Fill(ptgam, jes, w);
	    mh.pres->Fill(ptgam, res, w);
	  } // pass
	} // for ips in PSWeight
	} // for ieta in etas
      } // for ialpha in alphas

      if (doGamjet && hg) {

	gamjetHistos *h = hg;

	// Specific event selection
	bool pass_alpha = (pt2 < 1.00*ptgam || pt2 < pt2min);
	bool pass_eta = (fabs(Jet_eta[iJet]) < 1.3);
	bool pass = (pass_basic_ext && pass_alpha && pass_eta);

	if (pass) {
	
	  double jsf = (smearJets ? Jet_CF[iJet] : 1);
	  double ptjet = Jet_pt[iJet];
	  double ptavp = 0.5*(ptgam + ptjet); // Do better later, now pT,ave (not pT,avp)
	  h->hpt13->Fill(ptgam, w);
	  h->hpt13a->Fill(ptavp, w);
	  h->hpt13j->Fill(ptjet, w);
	  
	  h->pptg->Fill(ptgam, ptgam, w);
	  h->pptj->Fill(ptgam, ptjet, w);
	  
	  h->pres->Fill(ptgam, res, w);
	  h->pjsf->Fill(ptgam, jsf, w);
	  h->pm0->Fill(ptgam, mpf, w);
	  h->pm2->Fill(ptgam, mpf1, w);
	  h->pmn->Fill(ptgam, mpfn, w);
	  h->pmu->Fill(ptgam, mpfu, w);
	  
	  h->pm0x->Fill(ptgam, mpfx, w);
	  h->pm2x->Fill(ptgam, mpf1x, w);
	  
	  // Extras for FSR studies
	  h->pmnu->Fill(ptgam, mpfnu, w);
	  h->pmnx->Fill(ptgam, mpfnx, w);
	  h->pmux->Fill(ptgam, mpfux, w);
	  h->pmnux->Fill(ptgam, mpfnux, w);
	  
	  // Composition
	  h->prho->Fill(ptgam, fixedGridRhoFastjetAll, w);
	  h->pchf->Fill(ptgam, Jet_chHEF[iJet], w);
	  h->pnhf->Fill(ptgam, Jet_neHEF[iJet], w);
	  h->pnef->Fill(ptgam, Jet_neEmEF[iJet], w);
	  
	  // Alternative pT bins
	  h->presa->Fill(ptavp, res, w);
	  h->pm0a->Fill(ptavp, mpf, w);
	  h->pm2a->Fill(ptavp, mpf1, w);
	  h->pmna->Fill(ptavp, mpfn, w);
	  h->pmua->Fill(ptavp, mpfu, w);
	  //
	  h->presj->Fill(ptjet, res, w);
	  h->pm0j->Fill(ptjet, mpf, w);
	  h->pm2j->Fill(ptjet, mpf1, w);
	  h->pmnj->Fill(ptjet, mpfn, w);
	  h->pmuj->Fill(ptjet, mpfu, w);
	}
      } // doGamjet
      
      if (doGamjet2 && hg2) {

	gamjetHistos2 *h = hg2;

	// Specific event selection
	bool pass_alpha = (pt2 < 1.00*ptgam || pt2 < pt2min);
	//bool pass_eta = (abseta >= ymin && abseta < ymax);
	bool pass = (pass_basic_ext && pass_alpha);// && pass_eta);

	if (pass) {
	
	  //double res = Jet_RES[iprobe] / Jet_RES[itag];
	  double jsf = (smearJets ? Jet_CF[iJet] : 1);
	  
	  double abseta = fabs(Jet_eta[iJet]);
	  h->h2pteta->Fill(abseta, ptgam, w);
	  
	  h->p2res->Fill(abseta, ptgam, res, w);
	  h->p2jsf->Fill(abseta, ptgam, jsf, w);
	  h->p2m0->Fill(abseta, ptgam, mpf, w);
	  h->p2m2->Fill(abseta, ptgam, mpf1, w);
	  h->p2mn->Fill(abseta, ptgam, mpfn, w);
	  h->p2mu->Fill(abseta, ptgam, mpfu, w);
	  
	  h->p2m0x->Fill(abseta, ptgam, mpfx, w);
	  h->p2m2x->Fill(abseta, ptgam, mpf1x, w);
	  
	  // Extras for FSR studies
	  h->p2mnu->Fill(abseta, ptgam, mpfnu, w);
	  h->p2mnx->Fill(abseta, ptgam, mpfnx, w);
	  h->p2mux->Fill(abseta, ptgam, mpfux, w);
	  h->p2mnux->Fill(abseta, ptgam, mpfnux, w);
	}
      } // doGamjet2
      
    } // for jentry in nentries
    cout << endl << "Finished loop, writing file." << endl << flush;
    cout << "Processed " << _nevents << " events\n";
    cout << "Skipped " << _nbadevents_json << " events due to JSON ("
	 << (100.*_nbadevents_json/_nevents) << "%) \n"; //why not _ntot?
    cout << "Skipped " << _nbadevents_trigger << " events due to trigger ("
      	 << (100.*_nbadevents_trigger/_ntot) << "%) \n";
    cout << "Skipped " << _nbadevents_veto << " events due to veto ("
      	 << (100.*_nbadevents_veto/_nevents) << "%) \n";

    // Add extra plot for jet response vs photon pT
    if (isMC) {
      fout->cd("control");
      TH1D *hrgenvgen = prgen->ProjectionX("hrgenvgen");
      TH1D *hrgenvgam = prjet->ProjectionX("hrgenvgam");
      hrgenvgam->Divide(pgjet);
      TH1D *hrgen2vgen = prgen2->ProjectionX("hrgen2vgen");
      TH1D *hrgen2vgam = prjet2->ProjectionX("hrgen2vgam");
      hrgen2vgam->Divide(pgjet2);
      curdir->cd();
    }

    fout->Write();
    cout << "File written." << endl << flush;
    
    fout->Close();
    cout << "File closed." << endl << flush;

} // GamHistosFill::Loop()


// Code originally from jetphys/HistosFill.C
void GamHistosFill::PrintInfo(string info, bool printcout)
{
  //*ferr << info << endl << flush;
  if (printcout) cout << info << endl << flush;
}

// Code originally from jetphys/HistosFill.C
bool GamHistosFill::LoadJSON(string json)
{
  PrintInfo(string("Processing LoadJSON() with ") + json + " ...",true);
  ifstream file(json, ios::in);
  if (!file.is_open()) { assert(false); return false; }
  char c;
  string s, s2, s3, s4, s5;
  char s1[256];
  int rn(0), ls1(0), ls2(0), nrun(0), nls(0);
  file.get(c);
  if (c!='{') return false;
  while (file >> s and sscanf(s.c_str(),"\"%d\":",&rn)==1) {
    if (_gh_debug) PrintInfo(Form("\"%d\": ",rn),true);

    while (file.get(c) and c==' ') {};
    if (_gh_debug) { PrintInfo(Form("%c",c),true); assert(c=='['); }
    ++nrun;

    bool endrun = false;
    while (!endrun and file >> s >> s2 and (sscanf((s+s2).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3 or (file >> s3 and sscanf((s+s2+s3).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3) or 
            (file >> s4 and sscanf((s+s2+s3+s4).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3) or (file >> s5 and sscanf((s+s2+s3+s4+s5).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3))) {

      s2 = s1;
      if (s2=="]") { file >> s3; s2 += s3; }

      if (_gh_debug) PrintInfo(Form("[%d,%d,'%s']",ls1,ls2,s1),true);

      for (int ls = ls1; ls != ls2+1; ++ls) {
        _json[rn][ls] = 1;
        ++nls;
      }

      endrun = (s2=="]," || s2=="]}");
      if (_gh_debug and !endrun and s2!=",") { PrintInfo(string("s1: ")+s2,true); assert(s2==","); }
    } // while ls
    if (_gh_debug) PrintInfo("",true);

    if (s2=="]}") continue;
    else if (_gh_debug and s2!="],") PrintInfo(string("s2: ")+s2,true); assert(s2=="],");
  } // while run
  if (s2!="]}") { PrintInfo(string("s3: ")+s2,true); return false; }

  PrintInfo(string("Called LoadJSON() with ") + json + ":",true);
  PrintInfo(Form("Loaded %d good runs and %d good lumi sections",nrun,nls),true);
  return true;
} // LoadJSON

void GamHistosFill::LoadPU() {

  cout << endl << "GamHistosFill::LoadPU" << endl << flush;
  TDirectory *curdir = gDirectory;

  string eras[] =
    {"2016P8",/*"2016APVP8",*/"2016P8APV","2017P8", "2018P8",
     "2016QCD",/*"2016APVQCD",*/"2016QCDAPV","2017QCD", "2018QCD",
     "2016APV","2016FGH","2017","2018",
     "2022P8", "2022EEP8","2022QCD", "2022EEQCD"};
  //"2016BCD","2016EF","2016FGH",
  //"2017B","2017C","2017D","2017E","2017F",
  //"2018A","2018B","2018C","2018D"};
  const int neras = sizeof(eras)/sizeof(eras[0]);
  map<string, vector<string> > trigs;
  trigs["2016P8"].push_back("mc");
  trigs["2016APVP8"] = trigs["2017P8"] = trigs["2018P8"] = 
    trigs["2016QCD"] =  trigs["2016APVQCD"] = trigs["2017QCD"] =
    trigs["2018QCD"] = trigs["2018P8"] =
    trigs["2022P8"] = trigs["2022EEP8"] =
    trigs["2022QCD"] = trigs["2022EEQCD"] =
    trigs["2016P8"];

  trigs["2016APV"].push_back("HLT_Photon22_R9Id90_HE10_IsoM");
  trigs["2016APV"].push_back("HLT_Photon30_R9Id90_HE10_IsoM");
  trigs["2016APV"].push_back("HLT_Photon36_R9Id90_HE10_IsoM");
  trigs["2016APV"].push_back("HLT_Photon50_R9Id90_HE10_IsoM");
  trigs["2016APV"].push_back("HLT_Photon75_R9Id90_HE10_IsoM");
  trigs["2016APV"].push_back("HLT_Photon90_R9Id90_HE10_IsoM");
  trigs["2016APV"].push_back("HLT_Photon120_R9Id90_HE10_IsoM");
  trigs["2016APV"].push_back("HLT_Photon165_R9Id90_HE10_IsoM");
  //trigs["2016APV"].push_back("HLT_Photon165_HE10");
  trigs["2016APV"].push_back("HLT_Photon175");
  trigs["2016FGH"] = trigs["2016APV"];

  trigs["2017"].push_back("HLT_Photon20_HoverELoose");
  trigs["2017"].push_back("HLT_Photon30_HoverELoose");
  trigs["2017"].push_back("HLT_Photon50_R9Id90_HE10_IsoM");
  trigs["2017"].push_back("HLT_Photon75_R9Id90_HE10_IsoM");
  trigs["2017"].push_back("HLT_Photon90_R9Id90_HE10_IsoM");
  trigs["2017"].push_back("HLT_Photon120_R9Id90_HE10_IsoM");
  trigs["2017"].push_back("HLT_Photon165_R9Id90_HE10_IsoM");
  trigs["2017"].push_back("HLT_Photon200");

  trigs["2018"].push_back("HLT_Photon20_HoverELoose");
  trigs["2018"].push_back("HLT_Photon30_HoverELoose");
  trigs["2018"].push_back("HLT_Photon50_R9Id90_HE10_IsoM");
  trigs["2018"].push_back("HLT_Photon75_R9Id90_HE10_IsoM");
  trigs["2018"].push_back("HLT_Photon90_R9Id90_HE10_IsoM");
  trigs["2018"].push_back("HLT_Photon100EB_TightID_TightIso");
  trigs["2018"].push_back("HLT_Photon110EB_TightID_TightIso");
  trigs["2018"].push_back("HLT_Photon200");

  trigs["2022"].push_back("HLT_Photon20_HoverELoose");
  trigs["2022"].push_back("HLT_Photon30_HoverELoose");
  trigs["2022"].push_back("HLT_Photon30EB_TightID_TightIso");
  trigs["2022"].push_back("HLT_Photon50_R9Id90_HE10_IsoM");
  trigs["2022"].push_back("HLT_Photon75_R9Id90_HE10_IsoM");
  trigs["2022"].push_back("HLT_Photon90_R9Id90_HE10_IsoM");
  trigs["2022"].push_back("HLT_Photon100EBHE10");
  trigs["2022"].push_back("HLT_Photon110EB_TightID_TightIso");
  trigs["2022"].push_back("HLT_Photon200");

  // files/pileup.root updated with tchain.C on Hefaistos
  TFile *fmc = new TFile("files/pileup.root","READ");
  assert(fmc && !fmc->IsZombie());

  for (int i = 0; i != neras; ++i) {
    string se = eras[i]; const char *ce = se.c_str();
    for (unsigned int j = 0; j != trigs[se].size(); ++j) {
      string st = trigs[se][j]; const char *ct = st.c_str();

      // Read trigger threshold from trigger name
      int itrg(0);
      if (st=="mc") itrg = 1;
      else sscanf(ct,"HLT_Photon%d*",&itrg);
      
      TFile *fdt(0);
      TH1D *h(0);
      if (st=="mc") {
	h = (TH1D*)fmc->Get(Form("pileup_%s",ce));
	if (!h) cout << "Failed to find pileup_"<<ce<<endl<<flush;
	assert(h);
      }
      else {
	// data files from Laura (on CERNbox)
	fdt = new TFile(Form("pileup/%s/pu_%s.root",ce,ct),"READ");
	assert(fdt && !fdt->IsZombie());
	h = (TH1D*)fdt->Get("pileup");
	assert(h);
      }
      assert(h);

      curdir->cd();
      h = (TH1D*)h->Clone(Form("pileup_%s_%s",ce,ct));
      double lumi = h->Integral();
      h->Scale(1./lumi);
      _pu[se][itrg] = h;
      _lumi[se][itrg] = lumi;

      cout << Form("%s/%s (%d): %1.3g %s\n",ce,ct,itrg,
		   lumi,st=="mc" ? "events" : "fb-1");

      if (fdt) fdt->Close();
    } // for j in trigs
  } // for i in eras
  fmc->Close();
  cout << endl << flush;

  return;
} // LoadPU




//void lumimap() {
LumiMap GamHistosFill::LoadLumi(string filename){
    //the file(s) where lumi is stored for trigger
    fstream lumifile; //pointer to file

    //read in the .csv file for, e.g. 50EB, photon trigger
    //TFile* lumifile =;
    //lumifile.open("files/lumi2024_378981_380115_Golden_photon50eb.csv", ios::in);
    lumifile.open(filename, ios::in);
    LumiMap lumi;

    //reading (use vector of strings)
    vector<string> row;
    string currline, colentry, temp;

    //read the file
    int nrows=0;
    while(getline(lumifile, currline)){
        row.clear();

        //break words
        stringstream str(currline);

        //delimiter
        char delim = ',';
        //go through columns of this row
        while(getline(str, colentry, delim)){
            //add current column's entry to row vector
            row.push_back(colentry);
        }

        //show FIRST and LAST entry of row
        if(row.at(0)[0]!='#'){ //check if first entry isn't starting with # (as these would be comments)
            nrows+=1; //increase valid row counter by one
            //cout << "row (first element = run number): " << row.at(0) << endl; // first element
            //cout << "row (last element = recorded lumi): " << row.back() << endl;
            
            stringstream runstring;
            int runnum;
            runstring << row.at(0);
            runstring >> runnum;

            stringstream lumistring;
            double currlum; 
            lumistring << row.back();
            lumistring >> currlum;

            lumi.emplace(runnum, currlum); //fill value to lumimap
        }
    }//stop reading file

    //print the lumi map
    /*
    for(const auto& elem : lumi50){ //iterator through map's elements
        cout << "run#: " << elem.first << " lumi: " << elem.second << " " << endl;
    }
    */

    return lumi; //return the lumimap
}//end of lumimap
