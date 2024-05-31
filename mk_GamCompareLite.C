{
  //gROOT->ProcessLine(".exception");

  gROOT->ProcessLine(".L tools.C+g");

  // Link JEC libraries that are modified to work stand-alone
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  gROOT->ProcessLine(".L compareLite.C+g");
  gROOT->ProcessLine(".L drawCompareLite.C+g");

  //compareLite("TEST");
  //compareLite("2024B");

  drawCompareLite("2024B");
  
  // hadd compareLite_Run3.root compareLite_2022CD.root compareLite_2022E.root compareLite_2022FG.root compareLite_2023Cv123.root compareLite_2023Cv4.root compareLite_2023D.root

  /*
  drawCompareLiteIOVs();
  */
}
