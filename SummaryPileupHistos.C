//Reading and saving all the PU histogram for gamjet into one single summary .root file 
#include <TFile.h>
#include <TChain.h>
#include <TH1D.h>
#include <iostream>
#include <ctime>



void SummaryPileupHistos(string ver = "w41") {
    TChain chain("Events");

		//list of data eras to be saved to summary file:
		string data_eras[] = {"2024B", "2024C", "2024D", "2024Ev1", "2024Ev2", "2024F", "2024G", "2024H", "2024I", "2024GH"};
		const int neras = sizeof(data_eras)/sizeof(data_eras[0]); //.size();


		//list of triggers to be considered (for now only 50EB, but have others too)
		string trigs[] = {"HLT_Photon50EB_TightID_TightIso"};
		const int ntrigs = sizeof(trigs)/sizeof(trigs[0]); //.size();
		string trg[] = {"photon50eb"}; //short names



		//For SUMMARY file (will contain: MC pre-reweighting, MC post-reweighting, data pu
		TFile *summaryFile = new TFile(Form("pileup/pu_summary_%s.root", ver.c_str()), "RECREATE");


		// MONTE CARLO PU DISTRIBUTIONS  --> read from files and save to summary file.
		TFile *mfp8 = new TFile(Form("pileup/2024/pileup_winter2024P8.root")); //winter2024P8
		TFile *mfqcd = new TFile(Form("pileup/2024/pileup_2024QCD.root")); //2024QCD

		//MC PU pre-reweighting
		TH1D *mcp8_pu = (TH1D*)mfp8->Get(Form("pileup_winter2024P8"));
		TH1D *mcqcd_pu = (TH1D*)mfqcd->Get(Form("pileup_2024QCD"));


		//SAVING original distribution, and scaled distribution
		//get number of MC events in pu profile
		Float_t mcp8_nevt = mcp8_pu->Integral();
		Float_t mcqcd_nevt = mcqcd_pu->Integral();

		summaryFile->cd();
		mcp8_pu->Write(); 	//write pu distribution for P8 (not normalised)
		mcqcd_pu->Write();	//write pu distribution for QCD (not normalised)

		mcp8_pu->Scale(1./mcp8_nevt); //normalise to 1
		mcp8_pu->Write(Form("pileup_winter2024P8_scaled")); //write normalised pu dist for MC (BEFORE reweighting)

		mcqcd_pu->Scale(1./mcqcd_nevt); //normalise to 1
		mcqcd_pu->Write(Form("pileup_2024QCD_scaled")); //write normalised pu dist for MC (BEFORE reweighting)



		// DATA PU DISTRIBUTION
		//read it, save it to summary file, normalise it (save again) and use it for reweighting
		//for the data distribution (just read the histo and save it into the new file)
		//TFile *df = new TFile(Form("pileup/pileup-histos_2024_%s.root",version);
		//histograms are called like this: MyDataPileupHistogram_2024G_photon50eb_w41.root

		for(int i=0; i<neras; ++i){
				for(int j=0; j<ntrigs; ++j){
						TFile *df = new TFile(Form("/media/Duo2/bschilli/pileup2024_%s/MyDataPileupHistogram_%s_%s_%s.root", 
																	ver.c_str(), data_eras[i].c_str(), trg[j].c_str(), ver.c_str()),"READ"); //test it with w41

						string era_hlt_name = Form("%s_%s",data_eras[i].c_str(), trigs[j].c_str()); //e.g. "2024F_Photon50EB_TightID_TightIso"
						//TH1D *data_pu = (TH1D*)df->Get("pileup_2024F_Photon50EB_TightID_TightIso");
						//TH1D *data_pu = (TH1D*)df->Get(Form("pileup_%s", era_hlt_name.c_str()));
						TH1D *data_pu = (TH1D*)df->Get(Form("pileup"));
						summaryFile->cd();
						data_pu->Write(Form("pileup_%s", era_hlt_name.c_str()));	//write existing pu histo for data to summary file

						//get the integral to know how many events we talk about
						Float_t data_nevt = data_pu->Integral(); //it will be an integer, but use floating point number to not mess up division later
						data_pu->Scale(1./data_nevt); //normalise histogram to one --> is already normalised??? why does it have only 120 entries?
						data_pu->Write(Form("pileup_%s_scaled", era_hlt_name.c_str()));

						/*
						//reweight both MC pu distributions with this era+trigger combo and store to file
						//REWEIGHTED MC PU DISTRIBUTIONS (PU reweighting to get reweighted pu profile)
						TH1D *puweight = (TH1D*)data_pu->Clone();
						//divide data pu dist by mc pu dist: (should be the normalised ones?)
						puweight->Divide(mcp8_pu); //--> for each pu-bin get a pu weight, which is #(evts in data with this number of pu)/ #(evt in mc with this number of pu)
						TH1D *mcp8_new = (TH1D*) mcp8_pu->Clone(Form("pileup_%s_pu-reweighted_scaled",dataset.c_str()));

						//for each bin in the pu distribution histogram, reweight it with the pu weight
						for(Long64_t i = 0; i < (mcpu_new->GetNbinsX()); i++){ //should go from 1 to 120
								cout << "Looping through " << mcpu_new->GetNbinsX() << " bins to do the reweighting of the MC distribution." << endl << flush;
								//the pu weight for the current (pu)bin:
								puweight->GetBinContent(i); //shouldn't this give exactly the original data histo again?
								//--> leave this out... we get the reweighted stuff from GamHistosFill (mu, rho, NPV etc).
						}
						*/
				}//end of loop through trigs
		}//end of loop through eras



    // Clean up
    //delete MCpileup;
	
		summaryFile->Close();

    //std::cout << "Pileup histogram created and saved to pileup/GJ-4Jets_Winter24_madgraph-p8-MC_V14_PUProfile.root" << std::endl;
    std::cout << Form("Pileup histograms saved to pileup/pu_summary_%s.root", ver.c_str()) << std::endl;
}
