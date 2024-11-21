// Purpose: Fill gamma+jet analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: June 6, 2021
//#include "CondFormats/JetMETObjects/src/Utilities.cc"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "GamHistosFill.h"

#include "TSystem.h"

#include <fstream>
#include <string>

#define GPU
//#define LOCAL

#ifdef LOCAL
// Compile these libraries into *.so first with root -l -b -q mk_CondFormats.C
// (works for 6.18.04?)
/*
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+)

R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+)
*/
//R__LOAD_LIBRARY(GamHistosFill.C+g)
// As in jetphys/mk2_histosFill.C:
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector_cc)

R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectionUncertainty_cc)

R__LOAD_LIBRARY(GamHistosFill_C)
#else
// (works for 6.26/10)
R__LOAD_LIBRARY(GamHistosFill_C.so)
#endif


void mk_GamHistosFill(string dataset = "X", string puera = "", string version = "w42") { //using w-version names for my code (Bettina).
//void mk_GamHistosFill(string inputlist = "X", string version = "w36") { //using w-version names for my code (Bettina). TO DO


//void mk_GamHistosFill(string dataset = "X", string version = "wX23") { //using version wX23 for summer23 corrections TEST (with one single file per iov)
//void mk_GamHistosFill(string dataset = "X", string version = "wX22full") { //using version wX22full for summer23 corrections TEST (with all files, but old jec)

	//map the different parts of the dataset to the actual main dataset, less confusion in GamHistosFill - TEST IT
/*
	string dataset;
	if(inputlist.Contains("2024G")) dataset = "2024G"; //covers 2024Ga, 2024Gb
*/



  // Settings
  bool addData = (dataset=="2016B"  || dataset=="2016C" || dataset=="2016D" || 
		  dataset=="2016E"  || dataset=="2016F" || 
		  dataset=="2016FG" || dataset=="2016H" || 
		  dataset=="2016BCD"|| dataset=="2016EF"|| dataset=="2016FGH" ||
		  //dataset=="2016BCDEF" ||
		  dataset=="2017B" || dataset=="2017C" || dataset=="2017D" || 
		  dataset=="2017E" || dataset=="2017F" ||
		  //dataset=="2017BCDEF" || 
		  dataset=="2018A"  || dataset=="2018B" ||
		  dataset=="2018C"  || dataset=="2018D" ||
		  dataset=="2018A1" || dataset=="2018A2" ||
		  dataset=="2018D1" || dataset=="2018D2" ||
		  dataset=="2018D3" || dataset=="2018D4" ||
                  //dataset=="2018ABCD");
		  dataset=="2022C"  || dataset=="2022D" || dataset=="2022E" ||
		  dataset=="2022F"  || dataset=="2022G" ||
		  dataset=="2023Cv123X" || dataset=="2023Cv4X" || dataset=="2023DX" || //for my test wX23
		  dataset=="2023B" || dataset=="2023Cv123" ||
		  dataset=="2023Cv4" || dataset=="2023D" ||
      dataset=="2024B" || dataset=="2024C" || dataset=="2024D" || dataset=="2024Ev1" || dataset=="2024Ev2" || 
			dataset=="2024F" || dataset=="2024Fa" || dataset=="2024Fb" || dataset=="2024Fc" || dataset=="2024Fd" || dataset=="2024F-ECALCC-HCALDI-skim" ||
			dataset=="2024G" || dataset=="2024Ga" || dataset=="2024Gb" || dataset=="2024Gc" || dataset=="2024Gd" || dataset=="2024Ge" || dataset=="2024Gtest" ||
			dataset=="2024H" || dataset=="2024Hskim" || dataset=="2024Ha" || dataset=="2024Hb" || dataset=="2024Hc" || dataset=="2024Hd" || dataset=="2024He" || dataset=="2024Htest" ||
			dataset=="2024I" || dataset=="2024Iskim" || dataset=="2024Iv1" || dataset=="2024Iv2" || dataset=="2024Itest" ||
									dataset=="2024B-ECALRATIO" || dataset=="2024C-ECALRATIO" || dataset=="2024C-ECALR-HCALDI" || dataset=="2024C-ECALCC-HCALDI" ||
                  dataset=="2024B-PromptReco-v1");

  bool addMC = (dataset=="2016P8" || dataset=="2017P8" || dataset=="2018P8" ||
		dataset=="2016APVP8" ||
		dataset=="2022P8" || //dataset=="2022QCD" ||
		dataset=="2022EEP8" || //dataset=="2022EEQCD" ||
                dataset=="2023P8X" || dataset=="2023P8-BPixX"|| //for my test wX23
		dataset=="2023P8" || //);// || dataset=="2023QCD");
  		dataset=="2023P8-BPix" || //added the BPix MC files
  		dataset=="winter2024P8" || dataset=="winter2024P8-test" || 
			dataset=="winter2024P8-v14" || //winter 2024 madgraph p8 (added on 07.08.2024)
			dataset=="winter2024P8a" || dataset=="winter2024P8b" || dataset=="winter2024P8c"); //winter 2024 madgraph p8 (added on 07.08.2024)


  bool addQCD = (dataset=="2016QCD" || dataset=="2016APVQCD" || 
		 dataset=="2017QCD" || dataset=="2018QCD" ||
		 dataset=="2022QCD" || dataset=="2022EEQCD" ||
                 dataset=="2023QCDX" || dataset=="2023QCD-BPixX" || //for my test wX23
		 dataset=="2023QCD" || dataset=="2023QCD-BPix" ||  //added BPix QCD MC
		 dataset=="2024QCD" || dataset=="2024QCD-v14" || //); //added BPix QCD MC
		dataset=="2024QCDa" || dataset=="2024QCDb" || dataset=="2024QCDc" || dataset=="2024QCDd" || 
		dataset=="2024QCDe" || dataset=="2024QCDf"); //these do use v14 puppi despite not in name here


  //cout << "Clean old shared objects and link files" << endl << flush;
  //gSystem->Exec("rm *.d");
  //gSystem->Exec("rm *.so");
  //gSystem->Exec("rm *.pcm");	

  string path = gSystem->pwd();
  cout << "Current path: " << path << endl << flush;
  /*
  gSystem->AddIncludePath(Form("-I%s",path.c_str()));
  gSystem->AddIncludePath(Form("-I%s/CondFormats/JetMETObjects/interface",path.c_str()));
  */
#ifdef GPU
  // Compile these libraries into *.so first with root -l -b -q mk_CondFormats.C
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  gROOT->ProcessLine(".L GamHistosFill.C+g");
#endif


  TChain *c = new TChain("Events");
  //c->AddFile("files/mc-ht-200to400.root");
  //c->AddFile("files/mc-ht-600toInf.root");
  //c->AddFile("files/mc-pt-15to6000.root");
  //c->AddFile("files/data-2018a.root");
  
  // Automatically figure out where we are running the job
  bool runGPU = (path=="/media/storage/gamjet" ||
		  path=="/home/bschilli/Cern/gamjet-analysis" || //); //on Vulcan, I have this home directory (should fix this code, as i have same directory also locally...)
		  path=="/media/Duo2/bschilli/gamjet-analysis" || //new directory on vulcan
	  	  path=="/media/storage/bschilli/gamjet-analysis"); //hefaistos
  bool runLocal = (path=="/Users/voutila/Dropbox/Cern/gamjet" ||
		   path=="/Users/voutila/Library/CloudStorage/Dropbox/Cern/gamjet" ||
		   path=="/Users/manvouti/Dropbox/Cern/gamjet" ||
		   path=="/home/bschilli/Cern/gamjet");
  if (!runLocal) assert(runGPU);

  if (runGPU) cout << "Running on Hefaistos or Vulcan (runGPU)" << endl;
  if (runLocal) cout << "Running on locally (runLocal)" << endl;

  if (puera==""){
    cout << "No pileup reweighting applied." << endl << flush;
  }
  else {
    cout << "Using pileup for era " << puera << " for PU reweighting." << endl << flush;
  }

  
  if (addData) {
    ifstream fin(runLocal ? Form("input_files/dataFiles_local_Run%s.txt",dataset.c_str()) : 
		Form("input_files/dataFiles_Run%s.txt",dataset.c_str()), ios::in);
		//Form("input_files/dataFiles_Run%s.txt",inputlist.c_str()), ios::in);
    string filename;
    cout << "Chaining data files:" << endl << flush;
    int nFiles(0), nFilesMax(827);//9999);
    while (fin >> filename && nFiles<nFilesMax) {
      ++nFiles;
      c->AddFile(filename.c_str());
    }
    cout << "Chained " << nFiles <<  " files" << endl << flush;
    
    GamHistosFill filler(c,0,dataset,puera,version);
    filler.Loop();
  }
  
  if (addMC) {
    ifstream fin(runLocal ? Form("input_files/mcFiles_local_%s.txt",dataset.c_str()) :
		 Form("input_files/mcFiles_%s.txt",dataset.c_str()), ios::in);
    string filename;
    cout << "Chaining MC files:" << endl << flush;
    int nFiles(0), nFilesMax(1437);//100);
    while (fin >> filename && nFiles<nFilesMax) {
      ++nFiles;
      c->AddFile(filename.c_str());
    }
    cout << "Chained " << nFiles <<  " files" << endl << flush;
  
    GamHistosFill filler(c,1,dataset,puera,version);
    filler.Loop();
  }

  if (addQCD) {
    ifstream fin(runLocal ? Form("input_files/mcFiles_local_%s.txt",dataset.c_str()) :
		 Form("input_files/mcFiles_%s.txt",dataset.c_str()), ios::in);
    string filename;
    cout << "Chaining QCD MC files:" << endl << flush;
    int nFiles(0), nFilesMax(5097);
    while (fin >> filename && nFiles<nFilesMax) {
      ++nFiles;
      c->AddFile(filename.c_str());
    }
    cout << "Chained " << nFiles <<  " files" << endl << flush;
  
    GamHistosFill filler(c,1,dataset,puera,version);
    filler.Loop();
  }

}
