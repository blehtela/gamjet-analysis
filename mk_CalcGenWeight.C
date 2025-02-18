// makefile for using CalcGenWeight code, based on mk_GamHistosFill.C
// author:  blehtela
// created: 18th of February 2025

#include "CalcGenWeight.h"

#include "TSystem.h"

#include <fstream>
#include <string>

#define GPU
//#define LOCAL

void mk_CalcGenWeight(string pthtbin = "X") {

  // Settings
  bool addBin = (pthtbin=="summer2024_PTG10to100-HT40to100" || pthtbin=="summer2024_PTG10to100-HT100to200" ||
		pthtbin=="summer2024_PTG10to100-HT200to400" || pthtbin=="summer2024_PTG10to100-HT400to600" ||
		pthtbin=="summer2024_PTG10to100-HT600to1000" || pthtbin=="summer2024_PTG10to100-HT1000toInf" ||
		pthtbin=="summer2024_PTG100to200-HT40to200" || pthtbin=="summer2024_PTG100to200-HT200to400" || 
		pthtbin=="summer2024_PTG100to200-HT400to600" || pthtbin=="summer2024_PTG100to200-HT600to1000" || 
		pthtbin=="summer2024_PTG100to200-HT1000toInf" || pthtbin=="summer2024_PTG200toInf-HT40to400" || 
		pthtbin=="summer2024_PTG200toInf-HT400to600" || pthtbin=="summer2024_PTG200toInf-HT600to1000" || 
		pthtbin=="summer2024_PTG200toInf-HT1000toInf")

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
  gROOT->ProcessLine(".L CalcGenWeight.C+g");
#endif


  TChain *chainruns = new TChain("Runs");
  //c->AddFile("files/mc-ht-200to400.root"); //for test... (never used with this)
  
  // Automatically figure out where we are running the job
  bool runGPU = (path=="/media/storage/gamjet" ||
		  path=="/home/bschilli/Cern/gamjet-analysis" || //); //on Vulcan, I have this home directory (should fix this code, as i have same directory also locally...)
		  path=="/media/Duo2/bschilli/gamjet-analysis" || //new directory on vulcan
	  	  path=="/media/storage/bschilli/gamjet-analysis"); //hefaistos
  bool runLocal = (path=="/home/bschilli/Cern/gamjet");
  bool runLxplus = (path=="/afs/cern.ch/user/b/blehtela/private/gamjet-analysis");
  if (!runLocal && !runGPU) assert(runLxplus);

  if (runGPU) cout << "Running on Hefaistos or Vulcan (runGPU)" << endl;
  if (runLocal) cout << "Running locally (runLocal)" << endl;
  if (runLxplus) cout << "Running on lxplus (runLxplus)" << endl;

  
  //had this originally in CalcGenWeight.C, but moved here
  if (addBin) {
    //ifstream fin(runLocal ? Form("input_files/dataFiles_local_Run%s.txt",dataset.c_str()) : Form("input_files/dataFiles_Run%s.txt",dataset.c_str()), ios::in);
    ifstream fin(Form("input_files/mcFiles_%s.txt",pthtbin.c_str()) 
    string filename;
    cout << "Chaining MC files for:" << pthtbin << endl << flush;
    int nfile(0), nFilesMax(1000);
    while (fin >> filename && nfile<nFilesMax) {
      ++nfile;
      chainruns->AddFile(filename.c_str());
    }
    cout << "Chained " << nFiles <<  " files\n" << endl << flush;
    
    CalcGenWeight filler(c,0,pthtbin); //why the 0 in second argument? what is it? first event index?
    filler.Loop();
  }
  
}
