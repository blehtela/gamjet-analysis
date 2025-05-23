//Calculating the PU histogram for gamjet
//Code based on Nestor's PileupHistogram() module
#include <TFile.h>
#include <TChain.h>
#include <TH1D.h>
#include <iostream>
#include <ctime>



void GamHistosPileup(string dataset = "summer2024P8", string version = "w48") { //moved default from winter2024P8 to summer2024P8
    TChain chain("Events");


    // 2024 MC (v14) files for Gamma+Jet
    		if(dataset=="summer2024P8"){	//updated on 17.03.2025 to do summer24 reweighting
			chain.Add("root://hip-cms-se.csc.fi/store/user/rverma/cms-jerc-run3/Skim/GamJet/2024/MCSummer24/GJets/date-14Feb2025_time-111503_commit-40144c3/MCSummer24_2024_GamJet_GJets*PTG*.root"); //make sure it's the pt-binned one (ptg and ht, actually)
		}
		else if(dataset=="summer2024QCD"){
			chain.Add("root://hip-cms-se.csc.fi/store/user/rverma/cms-jerc-run3/Skim/GamJet/2024/MCSummer24/QCD/date-14Feb2025_time-111503_commit-40144c3/MCSummer24_2024_GamJet_QCD4JetsHT*.root");
		}
		else if(dataset=="winter2024P8"){
			chain.Add("/eos/cms/store/group/phys_jetmet/rverma/cms-jerc-run3/Skim/15Nov2024/2024/GamJet/MC_2024_GamJet_G4Jets*V14*.root");
		}
		else if(dataset=="2024QCD"){
			chain.Add("/eos/cms/store/group/phys_jetmet/rverma/cms-jerc-run3/Skim/15Nov2024/2024/GamJet/MC_2024_GamJet_QCD*V14*.root");
		}
		else if(dataset=="2023P8"){
			chain.Add("/eos/cms/store/group/phys_jetmet/rverma/cms-jerc-run3/Skim/15Nov2024/2023/GamJet/MC_2023PreBPix_GamJet_G4JetsHT*.root");
		}
		else if(dataset=="2023P8-BPix"){
			chain.Add("/eos/cms/store/group/phys_jetmet/rverma/cms-jerc-run3/Skim/15Nov2024/2023/GamJet/MC_2023PostBPix_GamJet_G4JetsHT*.root");
		}
		else if(dataset=="2023QCD"){
			chain.Add("/eos/cms/store/group/phys_jetmet/rverma/cms-jerc-run3/Skim/15Nov2024/2023/GamJet/MC_2023PreBPix_GamJet_QCD4JetsHT*.root");
		}
		else if(dataset=="2023QCD-BPix"){
			chain.Add("/eos/cms/store/group/phys_jetmet/rverma/cms-jerc-run3/Skim/15Nov2024/2023/GamJet/MC_2023PostBPix_GamJet_QCD4JetsHT*.root");
		}
		else if(dataset=="2022P8"){
			chain.Add("/eos/cms/store/group/phys_jetmet/rverma/cms-jerc-run3/Skim/15Nov2024/2022/GamJet/MC_2022PreEE_GamJet_G4JetsHT*.root");
		}
		else if(dataset=="2022EEP8"){
			chain.Add("/eos/cms/store/group/phys_jetmet/rverma/cms-jerc-run3/Skim/15Nov2024/2022/GamJet/MC_2022PostEE_GamJet_G4JetsHT*.root");
		}
		else if(dataset=="2022QCD"){
			chain.Add("/eos/cms/store/group/phys_jetmet/rverma/cms-jerc-run3/Skim/15Nov2024/2022/GamJet/MC_2022PreEE_GamJet_QCD4JetsHT*.root");
		}
		else if(dataset=="2022EEQCD"){
			chain.Add("/eos/cms/store/group/phys_jetmet/rverma/cms-jerc-run3/Skim/15Nov2024/2022/GamJet/MC_2022PostEE_GamJet_QCD4JetsHT*.root");
		}



   //on Vulcan, use:
   // FOR WINTER2024P8, files at "/media/Duo2/JME_NANO_MC/2024/GJ-4Jets_HT-*" (folder with several files).
   // FOR QCD, files at "/media/Duo2/JME_NANO_MC/2024/QCD-4Jets_HT-*" (folder with several files).


    // Disable all branches by default
    chain.SetBranchStatus("*", 0);

    // Enable the "Pileup_nTrueInt" and "genWeight" branches -- do they only exist in MC?
    chain.SetBranchStatus("Pileup_nTrueInt", 1);
    chain.SetBranchStatus("genWeight", 1);

    //New histogram for the MC pileup distribution
    //TH1D *MCpileup = new TH1D("pileup", "PUProfile", 130, 0, 130); //should switch to this, but my data histos were with 120bins
		TH1D *MCpileup = new TH1D(Form("pileup_%s", dataset.c_str()), "pileup", 120, 0, 120);


/*
		//For SUMMARY file (will contain: MC pre-reweighting, MC post-reweighting, data pu
		TFile *summaryFile = new TFile(Form("pu_summary_%s_2024F",dataset.c_str()), "RECREATE"); //change the hardcoded era to flexible one


		// DATA PU DISTRIBUTION (currently only 50EB trigger)  --> do this in extra script!
		//read it, save it to summary file, normalise it (save again) and use it for reweighting
		//for the data distribution (just read the histo and save it into the new file)
		//TFile *df = new TFile(Form("pileup/pileup-histos_2024_%s.root",version);
		TFile *df = new TFile(Form("pileup-histos_2024_w36.root"),"READ"); //test it with w36

		//list of data eras to be saved to summary file:
		string data_eras[] = {"2024B", "2024C", "2024D", "2024Ev1", "2024Ev2", "2024F", "2024G", "2024H", "2024I", "2024GH"};
		int neras = data_eras.size();

		for(int i=0; i<neras; ++i){
			string era_hlt_name = Form("%s_Photon50EB_TightID_TightIso", data_eras[i].c_str());
			//TH1D *data_pu = (TH1D*)df->Get("pileup_2024F_Photon50EB_TightID_TightIso");
			TH1D *data_pu = (TH1D*)df->Get(Form("pileup_%s", era_hlt_name.c_str()));
			summaryFile->cd();
			data_pu->Write();	//write existing pu histo for data to summary file

			//get the integral to know how many events we talk about
			Float_t data_nevt = data_pu->Integral(); //it will be an integer, but use floating point number to not mess up division later
			data_pu->Scale(1./data_nevt); //normalise histogram to one --> is already normalised??? why does it have only 120 entries?
			data_pu->Write(Form("pileup_%s_scaled", era_hlt_name.c_str()));
		} 
*/




		//for the reweighted distribution (see lower)


    // Define variables to hold the values of Pileup_nTrueInt and genWeight
    Float_t Pileup_nTrueInt;
    Float_t genWeight;

    // Set the branch address to connect the variables with the branches in the tree
    chain.SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt);
    chain.SetBranchAddress("genWeight", &genWeight);

    // Loop over all events and fill the histogram with the weight applied
    Long64_t nEntries = chain.GetEntries();

    std::time_t startTime = std::time(nullptr);
    for (Long64_t i = 0; i < nEntries; i++) {
        chain.GetEntry(i);
        MCpileup->Fill(Pileup_nTrueInt, genWeight);  // Use weighted fill

        // Print progress every 1000000 events
        if (i % 1000000 == 0) {
            // Calculate elapsed time
            std::time_t currentTime = std::time(nullptr);
            double elapsedTime = std::difftime(currentTime, startTime);

            // Estimate remaining time
            double averageTimePerEvent = elapsedTime / (i + 1);
            double remainingEvents = nEntries - i;
            double estimatedTimeRemaining = averageTimePerEvent * remainingEvents;

            // Print progress and estimated time
            std::cout << "Processing event " << i << " of " << nEntries
                      << ". Elapsed time: " << elapsedTime << " seconds"
                      << ". Estimated time remaining: " << estimatedTimeRemaining << " seconds." << std::endl;
        }
    }

		//get number of MC events in pu profile
		Float_t mc_nevt = MCpileup->Integral();

		//divide number of events in data by num. of events in mc, to get normalisation factor for (or wait, maybe better normalise to one?)
		//yeah, normalise both histograms to one, by dividing by their integral

    // Specify where to save the output ROOT file
    //TFile outFile("luminosityscripts/PUWeights/Winter24MGV14_PUProfile_v2.root", "RECREATE");
		//TFile outFile("pileup/GJ-4Jets_Winter24_madgraph-p8-MC_V14_PUProfile.root", "RECREATE");
		TFile outFile(Form("pileup/pileup_%s_%s.root", dataset.c_str(), version.c_str()), "RECREATE"); //file that only contains ONE histogram --> original MC pileup distr.
    MCpileup->Write(); //write pu distribution (not normalised)


/*
		summaryFile->cd();
		MCpileup->Write();
		//MCpileup->Scale(data_nevt/mc_nevt);
		MCpileup->Scale(1./data_nevt); 
		MCpileup->Write(Form("pileup_%s_scaled",dataset.c_str())); //write normalised pu dist for MC (BEFORE reweighting)
*/



/*
		//PU reweighting to get reweighted pu profile
		//Float_t 
		TH1D *puweight = (TH1D*)data_pu->Clone();
		//divide data pu dist by mc pu dist: (should be the normalised ones?)
		puweight->Divide(MCpileup); //--> for each pu-bin get a pu weight, which is #(evts in data with this number of pu)/ #(evt in mc with this number of pu)
		TH1D *mcpu_new = (TH1D*) MCpileup->Clone(Form("pileup_%s_pu-reweighted_scaled",dataset.c_str()));

		//for each bin in the pu distribution histogram, reweight it with the pu weight
		for(Long64_t i = 0; i < (mcpu_new->GetNbinsX()); i++){ //should go from 1 to 120
			cout << "Looping through " << mcpu_new->GetNbinsX() << " bins to do the reweighting of the MC distribution." << endl << flush;

			//the pu weight for the current (pu)bin:
			puweight->GetBinContent(i); //shouldn't this give exactly the original data histo again?
		}
*/


    // Clean up
    outFile.Close();
    delete MCpileup;
	
		//summaryFile->Close();

    //std::cout << "Pileup histogram created and saved to pileup/GJ-4Jets_Winter24_madgraph-p8-MC_V14_PUProfile.root" << std::endl;
    std::cout << Form("Pileup histogram created and saved to pileup/pileup_%s_%s.root", dataset.c_str(), version.c_str()) << std::endl;
}
