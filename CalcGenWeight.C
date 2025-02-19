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

	if(fChain == 0) return;

	//switch on the genEventSumw branch
	fChain->SetBranchStatus("*",0); //switch off all branches
	fChain->SetBranchStatus("genEventSumw",1); //switch on branch with genEventSumw, the gen event weight

	//Loop through events to add the get genEventSumw for each event
	Long64_t nentries = fChain->GetEntries();
	cout << "\nStarting loop over " << pthtbin << " with "<< nentries << " entries" << endl;

	//double summedGenWeight(0.0);
	Double_t summedGenWeight(0);
	for(Long64_t jentry=0; jentry<nentries; jentry++){
		Long64_t ientry = LoadTree(jentry);
		if(ientry<0) break;
		b_genEventSumw->GetEntry(ientry); //read from branch
		summedGenWeight += genEventSumw;
		//cout << "genEventSumw now: " << genEventSumw << endl << flush; 
		//cout << "current total sum: " << summedGenWeight << endl << flush;
    	
		if (jentry%1000000==0) cout << "." << flush;
      		if (jentry%50000000==0 && jentry!=0) cout << "\nn="<<jentry<<endl<<flush;
	}

	cout << "\nProcessed " << nentries << " entries." << endl << flush;



	//write to file: sample name, summedGenWeight (or write it to console for now)
	cout << "Total sum of gen event weights in sample " << pthtbin << " is: " << summedGenWeight << endl << flush;
	cout << "----------------------------------------\n" << endl << flush;

	ofstream outputfile;
	const char* filename = Form("genweight_%s.txt", pthtbin.c_str());
	//outputfile.open(filename, app); //operations at end of file only (app = append)
	outputfile.open(filename, fstream::app); //operations at end of file only (app = append)


	outputfile << pthtbin << " " << summedGenWeight << "\n";

	outputfile.close();
} //end of CalcGenWeight::Loop()
