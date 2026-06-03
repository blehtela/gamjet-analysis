// Calculating the overall GenWeight for a given sample
// author: blehtela
// created: 17th of February, 2025
#define CalcGenWeight_cxx
#include "CalcGenWeight.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std; 

//def
//void CalcGenWeight(string pthtbin);


//calling the fct several times
//void CalcGenWeight(string pthtbin = "X"){
//
//instead of the following, use calcAllGenWeights.py
/*
void CalcGenWeight(){
	CalcGenWeight("summer2024P8_PTG10to100-HT40to100");
	CalcGenWeight("summer2024P8_PTG10to100-HT100to200");
	CalcGenWeight("summer2024P8_PTG10to100-HT200to400");
	CalcGenWeight("summer2024P8_PTG10to100-HT400to600");
	CalcGenWeight("summer2024P8_PTG10to100-HT600to1000");
	CalcGenWeight("summer2024P8_PTG10to100-HT1000toInf");

	CalcGenWeight("summer2024P8_PTG100to200-HT40to200");
	CalcGenWeight("summer2024P8_PTG100to200-HT200to400");
	CalcGenWeight("summer2024P8_PTG100to200-HT400to600");
	CalcGenWeight("summer2024P8_PTG100to200-HT600to1000");
	CalcGenWeight("summer2024P8_PTG100to200-HT1000toInf");

	CalcGenWeight("summer2024P8_PTG200toInf-HT40to400");
	CalcGenWeight("summer2024P8_PTG200toInf-HT400to600");
	CalcGenWeight("summer2024P8_PTG200toInf-HT600to1000");
	CalcGenWeight("summer2024P8_PTG200toInf-HT1000toInf");
}
*/

//implementation
//void CalcGenWeight(string pthtbin="X"){
void CalcGenWeight::Loop() {
	//get files for current sample, e.g. MCSummer24_2024_GamJet_GJetsHT40to100PTG10to100
	//i call them input_files/mcFiles_summer2024P8_PTG10to100-HT40to100.txt
	//there are overall lists for the entire sample (processed with GamHistosFill, called: mcFiles_summer2024P8.txt

	// The following is now handled in mk_CalcGenWeight.C
	/*
	// TChain based on Runs Tree including all files for ONE pt-ht-bin
	TChain *chainrun = new TChain("Runs");
	ifstream fin(Form("input_files/mcFiles_%s.txt", pthtbin.c_str()));
	string filename;
	cout << "Chaining MC files for " << pthtbin << endl << flush;
	int nfile(0), nFilesMax(1000);
	while(fin >> filename && nfile<nFilesMax){
		++nfile;
		chainrun->AddFile(filename.c_str()); //add current file to the chain
	}

	//show file count to user
	cout << "Chained " << nFiles << " files\n" << endl << flush;
	*/

	cout << "now in loop!" << endl << flush;

	if(fChain == 0) return;

	//switch on the genEventSumw branch
	fChain->SetBranchStatus("*",0); //switch off all branches
	cout << "after set branchstatus * for all!" << endl << flush;
	fChain->SetBranchStatus("genEventSumw",1); //switch on branch with genEventSumw, the gen event weight
	cout << "after set branchstatus 1 for genEventSumw!" << endl << flush;
	fChain->SetBranchStatus("genEventCount",1); //switch on branch with genEventCount, the gen event count
	cout << "after set branchstatus 1 for genEventCount!" << endl << flush;


	//Loop through events to add the get genEventSumw for each event
	Long64_t nentries = fChain->GetEntries(); //this still gives the number of original files, even though running over skim?!
	cout << "\nStarting loop over " << pthtbin << " with "<< nentries << " entries" << endl;

	//double summedGenWeight(0.0);
	Double_t summedGenWeight(0);
	Double_t summedGenEventCount(0);

	for(Long64_t jentry=0; jentry<nentries; jentry++){
		Long64_t ientry = LoadTree(jentry);
		//cout << "ientry = " << ientry << endl << flush; //somehow does not go higher than 3...
		if(ientry<0) break;
		//b_genEventSumw->GetEntry(ientry); //read from branch (why from branch end not chain?)
		//b_genEventCount->GetEntry(ientry); //read from branch (why from branch end not chain?)
		fChain->GetEntry(ientry);
		summedGenWeight += genEventSumw;
		summedGenEventCount += genEventCount;
		//cout << "genEventCount now: " << genEventCount << endl << flush; 
		//cout << "genEventSumw now: " << genEventSumw << endl << flush; 
		//cout << "current total sum: " << summedGenWeight << endl << flush;
    	
		if (jentry%1000000==0) cout << "." << flush;
      		if (jentry%50000000==0 && jentry!=0) cout << "\nn="<<jentry<<endl<<flush;
	}

	cout << "\nProcessed " << nentries << " entries." << endl << flush;
	double prelimBinWeight = summedGenWeight/summedGenEventCount;
	const char* prelimBinWeightString = Form("%f", prelimBinWeight);
	/*
	cout << "double(prelimBinWeight): " << double(prelimBinWeight) << endl << flush;
	cout << "float(prelimBinWeight): " << float(prelimBinWeight) << endl << flush;
	cout << "summedGenWeight/summedGenEventCount: " << summedGenWeight/summedGenEventCount << endl << flush;
	cout << "float(summedGenWeight/summedGenEventCount): " << float(summedGenWeight/summedGenEventCount) << endl << flush;
	//cout << "setprecision(10)" << cout.setprecision(10) << prelimBinWeight << endl << flush;
	*/



	//summedGenWeight / summedEventCount
	//summedEventCount is not the same as nevents

	//write to file: sample name, summedGenWeight (or write it to console for now)
	cout << "Total sum of gen event weights in sample " << pthtbin << " is: " << summedGenWeight << endl << flush;
	cout << "Total sum of gen event counts in sample " << pthtbin << " is: " << summedGenEventCount << endl << flush;
	cout << "Normalised (by gen evt count) summedGenWeight in sample:  " << pthtbin << " is: " << prelimBinWeightString << endl << flush;
	cout << "----------------------------------------\n" << endl << flush;

	ofstream outputfile;
	//const char* filename = Form("genweight_%s.txt", pthtbin.c_str());
	const char* filename = Form("genweight_%s_newMay2026.txt", pthtbin.c_str());
	//outputfile.open(filename, app); //operations at end of file only (app = append)
	outputfile.open(filename, fstream::app); //operations at end of file only (app = append)

	//outputfile << pthtbin  << " " << nentries << " " << summedGenWeight << "\n"; //write bin, nevts, sumw to file
	const char* outputmessage = Form("File will contain: pthtbin (%s), nentries (#files in chain: %lld), summedGenWeight (%f), summedGenEventCount (%f), summedGenWeight/summedGenEventCount (pure sample bin weight: %f)", pthtbin.c_str(), nentries, summedGenWeight, summedGenEventCount, prelimBinWeight);
	cout << outputmessage << endl << flush;
	outputfile << "\npthtbin " << "nentries " << "summedGenWeight " << "summedGenEventCount " << "summedGenWeight/summedGenEventCount "<< "\n"; //write bin, nevts, sumw to file
	outputfile << pthtbin  << " " << nentries << " " << summedGenWeight << " " << summedGenEventCount << " " << prelimBinWeightString << "\n"; //write bin, nevts, sumw to file

	outputfile.close();
} //end of CalcGenWeight::Loop()
