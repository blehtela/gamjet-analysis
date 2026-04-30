// Run this script to compile CondFormats libraries. After this can easily run 
// root -l -b -q mk_GamHistosFill.C
// using R__LOAD_LIBRARY to load *.so
{

  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  //for JER smearing (since w80, 26.04.2026)
  //gROOT->ProcessLine(".L JetMETCorrections/Modules/src/JetResolution.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetResolutionObject.cc+"); //getting error here...
  gROOT->ProcessLine(".L JetMETCorrections/Modules/src/JetResolution.cc+");

  // For Gamjet code (v6.26.06)
  gROOT->ProcessLine(".L GamHistosFill.C+g");
  //gROOT->ProcessLine(".L DiGamHistosFill.C+g"); //commented out for now... do not need to compile every time

}
