// Purpose: Mix photon+jet and QCD using purity from drawPurityEstimates.C
// Output:  files/GamHistosMix_mc_[x]P8QCD_v[Y].root
//          This is used as input in GamHistosRatio.C
// Model:   Inspired by miquork/jecsys2020/recombine.C
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TF1.h"
#include "TKey.h"

#include <string>
#include <vector>
#include <iostream>

bool _debug = false;

// Scale EM jets (pi0->gamma+gamma) relative to prompt photons
// Possible causes: EM shower containment, pi>gamma+gamma is a bit wider
// also have primary energy coming from two lower energy photons
// Temporary: move this eventually to GamHistosFill.C for isQCD
bool scaleEMjet = true;
//double emJetScale = 0.99; // QCD [0.8,1.1] only
double emJetScale = 0.98; // QCD wider window and/or lower scale in data
// lower scale in data could connect to PF composition NEF in global fit?

// Take into account low gluon response in data, affecting QCD bkg wrt Z+jet
bool scaleGluonBkg = true;
//#include "../jecsys2020/Flavor.C"

void GamHistosMixEra(string epoch, string ver);
void recurseGamHistosFile(TDirectory *gamdir, TDirectory *qcddir,
			  TDirectory *outdir, int lvl = 0);

void GamHistosMix() {

  /*
  GamHistosMixEra("2016APV","v20");
  GamHistosMixEra("2016","v20");
  GamHistosMixEra("2017","v20");
  GamHistosMixEra("2018","v20");
  */

  // v27->v29->v30->v31
  /*
  GamHistosMixEra("2022","v32");
  GamHistosMixEra("2022EE","v32");
  GamHistosMixEra("Run3","v32");
  */

  //w1
  //GamHistosMixEra("2023","w14"); //make this contain MC without BPix issue
  //GamHistosMixEra("2023-BPix","w31"); //MC accounted for BPix issue
  //GamHistosMixEra("2023-BPix","w33"); //MC accounted for BPix issue
  //GamHistosMixEra("2024","w44");

  //2025
  GamHistosMixEra("2025", "w54"); //added on 02 Jun 2025


}

TF1 *_f1p(0);
TF1 *_f1qcd(0);
void GamHistosMixEra(string sepoch, string sver) {

  //const char *epoch = "2018";
  //const char *ver = "v19";
  const char *epoch = sepoch.c_str();
  string sqepoch = (sepoch=="2022EE"||sepoch=="Run3" ? "2022" : sepoch);//patch
  const char *qepoch = sqepoch.c_str();
  const char *ver = sver.c_str();
  //TFile *fout = new TFile(Form("rootfiles/GamHistosMix_mc_%sP8QCD_%s.root",
	//		       epoch,ver),"RECREATE");
	TFile *fout(0); //only for now with winter24p8 and 23-bpix QCD

  TFile *fgam(0);
  TFile *fqcd(0);

    //loading input files, needed if-condition due to different naming of BPix files. (TO DO: adjust this in GamHistosFill, so the if gets obsolete)
    if(TString(sepoch.c_str()).Contains("2023-BPix")){
        //fgam = new TFile(Form("rootfiles/GamHistosFill_mc_2023P8-BPix_%s.root", //needed to change this to different naming due to BPix samples.
        //                        ver),"READ");			                //e.g. GamHistosFill_mc_2023P8-BPix_w4.root
        fgam = new TFile(Form("rootfiles/GamHistosFill_mc_winter2024P8_%s.root", //for now the madgraph sample from 2024.. just interim result
                                ver),"READ");			 
	fout = new TFile(Form("rootfiles/GamHistosMix_mc_winter24P8_%sQCD_%s.root", epoch,ver),"RECREATE"); //only for now...



        assert(fgam && !fgam->IsZombie());
        //fqcd = new TFile(Form("rootfiles/GamHistosFill_mc_2023QCD-BPix_%s.root", //GamHistosFill_mc_2023QCD-BPix_w4.root
        //                        ver),"READ");
	fqcd = new TFile(Form("rootfiles/GamHistosFill_mc_2023QCD-BPix_w31.root"),"READ");

        assert(fqcd && !fqcd->IsZombie());
    }
    else{ //the usual way of loading the files (i.e. no BPix era)
        //TFile *fgam = new TFile(Form("rootfiles/GamHistosFill_mc_%sP8_%s.root", //this is how it usually is (GamHistosFill_mc_2023P8_w4.root)
        //fgam = new TFile(Form("rootfiles/GamHistosFill_mc_%sP8_%s.root",
				//       epoch,ver),"READ");
        //fgam = new TFile(Form("rootfiles/GamHistosFill_mc_winter%sP8_%s.root",epoch,ver),"READ"); //workaround as our mc is currently winter2024P8 instead of 2024P8
        //fgam = new TFile(Form("rootfiles/GamHistosFill_mc_winter%sP8-v14_%s.root",epoch,ver),"READ"); //for MC V14 stuff
        //fgam = new TFile(Form("rootfiles/GamHistosFill_mc_winter%sP8_pu-_%s.root",epoch,ver),"READ");
        fgam = new TFile(Form("rootfiles/GamHistosFill_mc_winter%sP8_no-pu_%s.root",epoch,ver),"READ");

	assert(fgam && !fgam->IsZombie());

	//fout = new TFile(Form("rootfiles/GamHistosMix_mc_winter24P8_%sQCD_%s.root", epoch,ver),"RECREATE"); //only for now... (delete when MC is called 2024P8 again)
	//fout = new TFile(Form("rootfiles/GamHistosMix_mc_winter24P8_%sQCD_V14_%s.root", epoch,ver),"RECREATE"); //for MC V14 stuff
	//fout = new TFile(Form("rootfiles/GamHistosMix_mc_winter2024P8_%sQCD_no-pu_%s.root", epoch,ver),"RECREATE"); //for MC V14 stuff
	fout = new TFile(Form("rootfiles/GamHistosMix_mc_winter%sP8_%sQCD_no-pu_%s.root",epoch,epoch,ver),"RECREATE"); //for MC V14 stuff



	//TFile *fqcd = new TFile(Form("rootfiles/GamHistosFill_mc_%sQCD_%s.root", //GamHistosFill_mc_2023QCD_w4.root
        //fqcd = new TFile(Form("rootfiles/GamHistosFill_mc_%sQCD-v14_%s.root",qepoch,ver),"READ"); //for the MC V14 stuff
        //fqcd = new TFile(Form("rootfiles/GamHistosFill_mc_%sQCD_%s.root",qepoch,ver),"READ"); //default
        //fqcd = new TFile(Form("rootfiles/GamHistosFill_mc_%sQCD_pu-_%s.root",qepoch,ver),"READ"); //for the MC V14 stuff
        fqcd = new TFile(Form("rootfiles/GamHistosFill_mc_%sQCD_no-pu_%s.root",qepoch,ver),"READ"); //for the MC V14 stuff

	assert(fqcd && !fqcd->IsZombie());
    }

  cout << "Calling GamHistosMix("<<epoch<<","<<ver<<");" << endl;
  cout << "Input files ";
  cout << fgam->GetName() << ",\n" << fqcd->GetName() << endl;
  cout << "Output file " << fout->GetName() << endl;
  cout << "Starting recursive loop. This may take a minute" << endl << flush;

  // Setup purity fit
  //assert(_f1p==0);
  // Hand-crafted "fit" to 2016+2018 control/pfakeqcd2x
  // processed through drawPurityEstimates.C
  //if (!_f1p) _f1p = new TF1("f1p","[0]+[1]*exp(-[2]*x)+[3]*x",15,3500);
  ////_f1p->SetParameters(7e-5,0.008,0.1,4e-8);
  //_f1p->SetParameters(0.1737, 0.9318, 0.01506, 0); // chi2/NDF=3.4/19

  // 2022 variant of manual "fit" of control/pfakeqcd2 through
  // drawPurityEstimates
  if (!_f1p) _f1p = new TF1("f1p","min([p0]+[p1]/x+[p2]/(x*x)+[p3]/(x*x*x),1.)",15,3500);
  //_f1p->SetParameters(0.02505, 9.455, -214.3, 3.88e+04); // chi2/NDF=8.0/19
  // Redo fit for pT>50 GeV that seems reliable for photon+jet xsec
  _f1p->SetParameters(0.01305, 4.573, -59.35, 2.299e+04); // chi2/NDF=2.2/14

  // From drawPurityEstimates.C:drawGluonResponse()
  if (!_f1qcd) _f1qcd = new TF1("f1qcd","[p0]+[p1]*pow(x,[p2])",15,1650);
  _f1qcd->SetParameters(1.0015, -0.11621, -0.56594);

  // Loop over all the directories recursively                                 
  recurseGamHistosFile(fgam,fqcd,fout);

  cout << endl;
  cout << "Recursive loop done." << endl;
  cout << "Output written in " << fout->GetName() << endl;
  fout->Close();
  cout << "Output file closed" << endl;
  fout->Delete();
  cout << "Output file pointer deleted" << endl << flush;

  fgam->Close();
  fgam->Delete();
  fqcd->Close();
  fqcd->Delete();
} // void GamHistosMix

void recurseGamHistosFile(TDirectory *gamdir, TDirectory *qcddir, 
			  TDirectory *outdir, int lvl) {
  
  // Automatically go through the list of keys (directories)                   
  TList *keys = gamdir->GetListOfKeys();
  TListIter itkey(keys);
  TObject *obj;
  TKey *key;

  while ( (key = dynamic_cast<TKey*>(itkey.Next())) ) {
    const char *c = Form("%%%ds%%s/%%s",lvl+1);
    if (_debug) cout << Form(c," ",gamdir->GetName(),key->GetName())
		     << endl << flush;
    obj = key->ReadObj(); assert(obj);
		if (TString(key->GetName()).Contains("_ps")) continue; //added this for w33, might need to remove it again


    // Found a subdirectory: copy it to output and go deeper                   
    if (obj->InheritsFrom("TDirectory")) {
      int loclvl = lvl;
      outdir->mkdir(obj->GetName());
      bool enteroutdir = outdir->cd(obj->GetName());
      assert(enteroutdir);
      TDirectory *outdir2 = outdir->GetDirectory(obj->GetName());
      assert(outdir2);
      outdir2->cd();

      TDirectory *gamdir2 = gamdir->GetDirectory(obj->GetName());
      TDirectory *qcddir2 = qcddir->GetDirectory(obj->GetName());

      if (loclvl>=0) loclvl++;

      //if (loclvl==1) 
      if ((string(gamdir2->GetName())!="runs" &&
	   string(gamdir2->GetName())!="Gamjet")) { // patch v30
	cout << endl << "Entering: " << gamdir2->GetName() << endl;
	recurseGamHistosFile(gamdir2, qcddir2, outdir2, loclvl);
      }
      else
	cout << endl << "Skip 'runs' and 'Gamjet' for v30" << endl;
    } // inherits from TDirectory  
    else if (obj->InheritsFrom("TH1")) { // Combine histograms

      outdir->cd();
      //string name = obj->GetName();
      //string trgname = indir->GetName();

      // Copy photon+jet histogram to output
      TObject *obj2 = obj->Clone(key->GetName()); // Copy input histo to output
      TH1D *hout(0);
      TProfile *pout(0);
      if (obj->InheritsFrom("TProfile")) {
	pout = dynamic_cast<TProfile*>(obj2);
	pout->SetName(Form("tmp_%s",key->GetName()));
	hout = pout->ProjectionX(key->GetName());
	//pout->Delete();
	pout->Reset();
      }
      // TProfile also inherits from TH1D, so need else if here!
      else if (obj->InheritsFrom("TH1D")) hout = dynamic_cast<TH1D*>(obj2);
      //else if (obj->InheritsFrom("TH1F")) hout = dynamic_cast<TH1F*>(obj2);
      //else if (obj->InheritsFrom("TH2D")) hout = dynamic_cast<TH2D*>(obj2);
      if (!hout) { // TH1F, TH2D, TH3D etc.
	obj2->Delete();
	continue;
      }

      // Start adding photon+jet signal and QCD backgrounds
      hout->Reset();
      TH1D *hingam = (TH1D*)gamdir->Get(key->GetName());
			cout << "key->GetName(): " << key->GetName() << endl;
      TH1D *hinqcd = (TH1D*)qcddir->Get(key->GetName()); assert(hinqcd);
      
      // Project TProfiles to TH1D for easier handling
      if (obj->InheritsFrom("TProfile")) {
	assert(hingam->InheritsFrom("TProfile"));
	TProfile *pingam = dynamic_cast<TProfile*>(hingam);
	hingam = pingam->ProjectionX(Form("p%s_gam",key->GetName()));

	assert(hinqcd->InheritsFrom("TProfile"));
	TProfile *pinqcd = dynamic_cast<TProfile*>(hinqcd);
	hinqcd = pinqcd->ProjectionX(Form("p%s_qcd",key->GetName()));
      }

      if (!hingam) {
	cout << outdir->GetName()<<"/"<<key->GetName()
	     << ": missing photon+jet from input"<<endl;
	continue;
      }
      if (!hinqcd) {
	cout << outdir->GetName()<<"/"<<key->GetName()
	     << ": missing QCD dijet from input"<<endl;
	continue;
      }
      
      if (hinqcd->GetNbinsX()!=hout->GetNbinsX()) {
	cout << outdir->GetName()<<"/"<<key->GetName()
	     << ": bin # mismatch ("
	     << hinqcd->GetNbinsX() << " vs " << hout->GetNbinsX() <<")"
	     << endl;
	continue;
      } // bin # mistmatch
	
      // Scale QCD "photon" i.e. EM jet based on MC best estimate
      if (scaleEMjet && obj->InheritsFrom("TProfile")) {
	TString s(key->GetName());
	if (s.Contains("resp_PtGam")    || s.Contains("resp_PtBalchs") ||
	    s.Contains("resp_DBchs")    || s.Contains("resp_MPFchs") ||
	    s.Contains("resp_MPFR1chs") || s.Contains("resp_MPFRnchs") ||
	    s.Contains("resp_MpfRuchs"))
	  hinqcd->Scale(1./emJetScale);
      }

      // Scale gluon jets in EM+jet wrt photon+jet for data
      if (scaleGluonBkg && obj->InheritsFrom("TProfile")) {
	TString s(key->GetName());
	if (s.Contains("resp_PtBalchs") || s.Contains("resp_DBchs") ||
	    s.Contains("resp_MPFchs") || s.Contains("resp_MPFR1chs"))
	  //hinqcd->Scale(_f1qcd);
	  for (int i = 1; i != hinqcd->GetNbinsX()+1; ++i) {
	    double k = _f1qcd->Eval(hinqcd->GetBinCenter(i));
	    hinqcd->SetBinContent(i, hinqcd->GetBinContent(i) * k);
	    hinqcd->SetBinError(i, hinqcd->GetBinError(i) * k);
	  }
      }

      for (int j = 1; j != hout->GetNbinsX()+1; ++j) {

	// Check that binning is ok
	if (hinqcd->GetBinLowEdge(j)!=hout->GetBinLowEdge(j)) {
	  cout << outdir->GetName()<<"/"<<key->GetName()
	       << Form(": bin mismatch ([%1.0f,%1.0f] vs [%1.0f,%1.0f]), ",
		       hinqcd->GetBinLowEdge(j), hinqcd->GetBinLowEdge(j+1),
		       hout->GetBinLowEdge(j), hout->GetBinLowEdge(j+1))
	       << endl;
	  continue;
	} // check binning

	// Combine results
	double pt = hout->GetBinCenter(j);
	double P = 1 - min(1.,max(0.,_f1p->Eval(pt)));
	hout->SetBinContent(j, (P * hingam->GetBinContent(j) +
				(1-P) * hinqcd->GetBinContent(j)));
	hout->SetBinError(j, sqrt(pow(P * hingam->GetBinError(j),2) +
				  pow((1-P) * hinqcd->GetBinError(j),2)));
      } // for j in hout
      
      // Save the stuff into an identical directory                            
      outdir->cd();
      assert(hout);
      // Hack TProfile back close to original
      if (obj->InheritsFrom("TProfile")) {
	assert(pout);
	hout->SetName(Form("tmp2_%s",key->GetName()));
	pout = (TProfile*)pout->Clone(key->GetName());
	pout->Reset();
	for (int i = 1; i != hout->GetNbinsX()+1; ++i) {
	  double pt = hout->GetBinCenter(i);
	  double y = hout->GetBinContent(i);
	  double ey = hout->GetBinError(i);
	  pout->Fill(pt, y+ey, 0.5);
	  pout->Fill(pt, y-ey, 0.5);
	} // for  i in pout
	pout->Write();
	pout->Delete();
      }
      else {
	hout->Write();
      }
      hout->Delete();
      gamdir->cd();
      // inherits from TH1
    }
    else
      assert(false);
    
    obj->Delete();
  } // while key      

} // recurseGamHistosFile




