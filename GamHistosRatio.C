// Purpose: produce data/MC ratios for JEC global fit
//          copy data and MC into the combination file
//          leave controls still in separate data and MC files
#include "TSystem.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TH2D.h"

#include <iostream>
using namespace std;

// Patches to input raw gamma+jet
// v15: moved these hard-coded into GamHistosFill, don't redo
bool addGain1 = false;//true;
//double gain1 = 0.99; // 1.02
bool addFootprint = false;//true;
bool patchGx = true; // patch Gxi to missing Gx

// Remove zeroes
void cleanGraph(TGraphErrors *g) {
  for (int i = g->GetN()-1; i != -1; --i) {
    if (g->GetY()[i]==0 && g->GetEY()[i]==0) g->RemovePoint(i);
  }
} // cleanGraph

// Replace pT in graph
void replacePt(TGraphErrors *g, TH1 *h) {
  for (int i = 0; i != g->GetN(); ++i) {
    double x = g->GetX()[i];
    int j = h->FindBin(x);
    double xnew = h->GetBinContent(j);
    double exnew = h->GetBinError(j);
    g->SetPoint(i, xnew, g->GetY()[i]);
    g->SetPointError(i, exnew, g->GetEY()[i]);
  }

  return;
}

void GamHistosRatios(string ver, string iov);

// v27->v29->v30->v31(L2L3Res_V2)->v32(L2L3Res_V3)
void GamHistosRatio(string ver = "w16") { 
  //GamHistosRatios(ver,"2016BCDEF");
  //GamHistosRatios(ver,"2016FGH");
  //GamHistosRatios(ver,"2017BCDEF");
  //GamHistosRatios(ver,"2018ABCD");
  //GamHistosRatios(ver,"Run2");

/*
  GamHistosRatios(ver,"2022C");
  GamHistosRatios(ver,"2022D");
  GamHistosRatios(ver,"2022E");
  //GamHistosRatios(ver,"2022CDE");
  GamHistosRatios(ver,"2022F");
  GamHistosRatios(ver,"2022G");
*/
  //GamHistosRatios(ver,"2022FG");
  //GamHistosRatios(ver,"2023Cv123");
  //GamHistosRatios(ver,"2023Cv4");

  //GamHistosRatios(ver,"2022CD");
  //GamHistosRatios(ver,"2022E");
  //GamHistosRatios(ver,"2022FG");

  //GamHistosRatios(ver,"2023Cv123");
  //GamHistosRatios(ver,"2023Cv4");
  //GamHistosRatios(ver,"2023D");
  GamHistosRatios(ver,"2024B");
  GamHistosRatios(ver,"2024C");
  GamHistosRatios(ver,"2024D");
  //GamHistosRatios(ver, "2024BCD");


  //GamHistosRatios(ver,"2023Cv4D");
  //GamHistosRatios(ver,"Run3");
  
}

void GamHistosRatios(string ver, string iov) {

  const char *cv = ver.c_str();
  const char *ci = iov.c_str();
  const char *cvo = cv; // patch to use older version of data, if needed

  TFile *fd(0), *fm(0), *fr(0);
  if (iov=="2016BCDEF") {
    // Merge files, if not already done (delete combination file to redo)
    //if (ps=="")
    gSystem->Exec(Form("hadd files/GamHistosFill_data_2016BCDEF_%s.root "
		       "files/GamHistosFill_data_2016BCD_%s.root "
		       "files/GamHistosFill_data_2016EF_%s.root ",
		       cv,cvo,cvo));

    fd = new TFile(Form("files/GamHistosFill_data_%s_%s.root",ci,cv),"READ");
    fm = new TFile(Form("files/GamHistosMix_mc_2016APVP8QCD_%s.root",cv),"READ");
    fr = new TFile(Form("files/GamHistosRatio_%s_P8QCDAPV_%s.root",ci,cv),
		   "RECREATE");//ps=="" ? "RECREATE" : "UPDATE");
  }
  if (iov=="2016FGH") {

    fd = new TFile(Form("files/GamHistosFill_data_%s_%s.root",ci,cvo),"READ");
    fm = new TFile(Form("files/GamHistosMix_mc_2016P8QCD_%s.root",cv),"READ");
    fr = new TFile(Form("files/GamHistosRatio_%s_P8QCD_%s.root",ci,cv),
		   "RECREATE");//ps=="" ? "RECREATE" : "UPDATE");
  }
  if (iov=="2017BCDEF") {
    // Merge files, if not already done (delete combination file to redo)
    //if (ps=="")
    gSystem->Exec(Form("hadd files/GamHistosFill_data_2017BCDEF_%s.root "
		       "files/GamHistosFill_data_2017B_%s.root "
		       "files/GamHistosFill_data_2017C_%s.root "
		       "files/GamHistosFill_data_2017D_%s.root "
		       "files/GamHistosFill_data_2017E_%s.root "
		       "files/GamHistosFill_data_2017F_%s.root ",
		       cv, cvo,cvo,cvo,cvo,cvo));

    fd = new TFile(Form("files/GamHistosFill_data_%s_%s.root",ci,cv),"READ");
    fm = new TFile(Form("files/GamHistosMix_mc_2017P8QCD_%s.root",cv),"READ");
    fr = new TFile(Form("files/GamHistosRatio_%s_P8QCD_%s.root",ci,cv),
		   "RECREATE");//ps=="" ? "RECREATE" : "UPDATE");
  
  }
  if (iov=="2018ABCD") {

    // Merge files, if not already done (delete combination file to redo)
    //if (ps=="")
    gSystem->Exec(Form("hadd files/GamHistosFill_data_2018ABCD_%s.root "
		       "files/GamHistosFill_data_2018A1_%s.root "
		       "files/GamHistosFill_data_2018A2_%s.root "
		       "files/GamHistosFill_data_2018B_%s.root "
		       "files/GamHistosFill_data_2018C_%s.root "
		       "files/GamHistosFill_data_2018D1_%s.root "
		       "files/GamHistosFill_data_2018D2_%s.root "
		       "files/GamHistosFill_data_2018D3_%s.root "
		       "files/GamHistosFill_data_2018D4_%s.root",
		       cv, cvo,cvo,cvo,cvo, cvo,cvo,cvo,cvo));
    
    fd = new TFile(Form("files/GamHistosFill_data_%s_%s.root",ci,cv),"READ");
    fm = new TFile(Form("files/GamHistosMix_mc_2018P8QCD_%s.root",cv),"READ");
    fr = new TFile(Form("files/GamHistosRatio_%s_P8QCD_%s.root",ci,cv),
		   "RECREATE");//ps=="" ? "RECREATE" : "UPDATE");
  }
  //
  if (iov=="2022C" || iov=="2022D" || iov=="2022CD" || iov=="2022CDE") {

    // Merge files, if not already done (delete combination file to redo)
    // Now updated into addAllIOVs.py
    /*
    if (iov=="2022CD")    
      gSystem->Exec(Form("hadd files/GamHistosFill_data_2022CD_%s.root "
			 "files/GamHistosFill_data_2022C_%s.root "
			 "files/GamHistosFill_data_2022D_%s.root",
			 cv, cv,cv));
    if (iov=="2022CDE")    
      gSystem->Exec(Form("hadd files/GamHistosFill_data_2022CDE_%s.root "
			 "files/GamHistosFill_data_2022C_%s.root "
			 "files/GamHistosFill_data_2022D_%s.root "
			 "files/GamHistosFill_data_2022E_%s.root",
			 cv, cv,cv,cv));
    */

    fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",ci,cv),"READ");
    //fm = new TFile(Form("files/GamHistosMix_mc_2018P8QCD_%s.root",cv),"READ");
    fm = new TFile(Form("rootfiles/GamHistosFill_mc_2022P8_%s.root",cv),"READ");
    //fr = new TFile(Form("files/GamHistosRatio_%s_P8QCD_%s.root",ci,cv),
    fr = new TFile(Form("rootfiles/GamHistosRatio_%s_P8_%s.root",ci,cv),"RECREATE");

    //fm = new TFile(Form("rootfiles/GamHistosMix_mc_2022P8QCD_%s.root",cv),"READ");
    //fr = new TFile(Form("rootfiles/GamHistosRatio_%s_P8QCD_%s.root",ci,cv),"RECREATE");
  }
  if (iov=="2022E" || iov=="2022F" || iov=="2022G" || iov=="2022FG") {

    // Merge files, if not already done (delete combination file to redo)
    /*
    if (iov=="2022FG")    
      gSystem->Exec(Form("hadd files/GamHistosFill_data_2022FG_%s.root "
			 "files/GamHistosFill_data_2022F_%s.root "
			 "files/GamHistosFill_data_2022G_%s.root",
			 cv, cv,cv));
    */
    fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",ci,cv),"READ");
    fm = new TFile(Form("rootfiles/GamHistosFill_mc_2022EEP8_%s.root",cv),"READ"); // no QCD???
    fr = new TFile(Form("rootfiles/GamHistosRatio_%s_P8_%s.root",ci,cv),"RECREATE"); // no QCD???
    //fm = new TFile(Form("rootfiles/GamHistosMix_mc_2022EEP8QCD_%s.root",cv),"READ");
    //fr = new TFile(Form("rootfiles/GamHistosRatio_%s_P8QCD_%s.root",ci,cv),"RECREATE");
  }
  /*
  if (iov=="2023Cv123" || iov=="2023Cv4" || iov=="2023Cv4D" ||
      iov=="2023D") {
    fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",ci,cv),"READ");
    //fm = new TFile(Form("rootfiles/GamHistosFill_mc_2022P8_%s.root",cv),"READ");
    //fr = new TFile(Form("rootfiles/GamHistosRatio_%s_P8_%s.root",ci,cv),"RECREATE");
    
    //fm = new TFile(Form("rootfiles/GamHistosMix_mc_2023P8QCD_%s.root",cv),"READ"); //gam+qcd MC
    //fm = new TFile(Form("rootfiles/GamHistosMix_mc_2023P8-BPixQCD_%s.root",cv),"READ"); //gam+qcd MC (BPix)
    //fm = new TFile(Form("rootfiles/GamHistosMix_mc_2022P8-BPix_%s.root",cv),"READ"); //gam MC (BPix)
    fm = new TFile(Form("rootfiles/GamHistosMix_mc_2023P8_%s.root",cv),"READ"); //gam MC: have no gamma+qcd mix yet for these, so need to do this now without qcd --> use one of the two lines above once we have it. (needed to hardcode version in this case... )
    //fr = new TFile(Form("rootfiles/GamHistosRatio_%s_P8QCD_%s.root",ci,cv),"RECREATE"); //switch back to this when we have QCD
    fr = new TFile(Form("rootfiles/GamHistosRatio_%s_P8-noQCD_%s.root",ci,cv),"RECREATE");
  }
  */
  //SHOULD RETURN TO THE ABOVE HANDLING AFTER THIS ONE, WHEN EVERYTHING HAS GONE TO VERSION W4 -- cannot because of BPIX issue
  //THE TWO IF CONDITIONS BELOW ARE ONLY FOR 16.02.2024 WHEN RUNNING WITHOUT QCD AND WITHOUT L2L3RES - TODO: update this.
      if (iov=="2023Cv123" || iov=="2023Cv4") {
        fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",ci,cv),"READ");
        //fm = new TFile(Form("rootfiles/GamHistosMix_mc_2023P8QCD_%s.root",cv),"READ"); 
        //fr = new TFile(Form("rootfiles/GamHistosRatio_%s_P8QCD_%s.root",ci,cv),"RECREATE");
        fm = new TFile(Form("rootfiles/GamHistosFill_mc_2023P8_%s.root",cv),"READ"); //NO QCD (remove later)
        fr = new TFile(Form("rootfiles/GamHistosRatio_%s_P8-noQCD_%s.root",ci,cv),"RECREATE"); //when no QCD in the mix (can remove it later)
      }
      //if (iov=="2023D") { //BPix issue (GamHistosMix_mc_2023-BPixP8QCD_w4.root) //should adjust the naming in the other programs!
      if (iov=="2023D" || iov=="2024B" || iov=="2024C" || iov=="2024D" || iov=="2024BCD") { //BPix issue, for now use this also for early 2024 data
        fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",ci,cv),"READ");
        fm = new TFile(Form("rootfiles/GamHistosMix_mc_2023-BPixP8QCD_%s.root",cv),"READ"); 
        fr = new TFile(Form("rootfiles/GamHistosRatio_%s_P8BPixQCD_%s.root",ci,cv),"RECREATE"); //changed name!!
        //fm = new TFile(Form("rootfiles/GamHistosFill_mc_2023P8-BPix_%s.root",cv),"READ"); //NO QCD (remove later)
        //fr = new TFile(Form("rootfiles/GamHistosRatio_%s_P8BPix-noQCD_%s.root",ci,cv),"RECREATE"); //when no QCD in the mix (can remove it later)
      }
  if (iov=="Run3") {

    // Merge files, if not already done (delete combination file to redo)
    /*
    gSystem->Exec(Form("hadd files/GamHistosFill_data_Run3_%s.root "
		       "files/GamHistosFill_data_2022C_%s.root "
		       "files/GamHistosFill_data_2022D_%s.root "
		       "files/GamHistosFill_data_2022E_%s.root "
		       "files/GamHistosFill_data_2022F_%s.root "
		       "files/GamHistosFill_data_2022G_%s.root "
		       "files/GamHistosFill_data_2023Cv123_%s.root "
		       "files/GamHistosFill_data_2023Cv4_%s.root",
		       cv, cv,cv,cv, cv,cv, cv,cv));
    gSystem->Exec(Form("hadd files/GamHistosFill_mc_2022XP8_%s.root "
		       "files/GamHistosFill_mc_2022EEP8_%s.root "
		       "files/GamHistosFill_mc_2022P8_%s.root",
		       cv, cv,cv));
    */

    fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",ci,cv),"READ");
    //fm = new TFile(Form("files/GamHistosFill_mc_2022XP8_%s.root",cv),"READ");
    //fm = new TFile(Form("rootfiles/GamHistosFill_mc_Run3P8_%s.root",cv),"READ");
    //fr = new TFile(Form("rootfiles/GamHistosRatio_%s_P8_%s.root",ci,cv),"RECREATE");
    fm = new TFile(Form("rootfiles/GamHistosMix_mc_Run3P8QCD_%s.root",cv),"READ");
    fr = new TFile(Form("rootfiles/GamHistosRatio_%s_P8QCD_%s.root",ci,cv),"RECREATE");
  }
  
  assert(fd && !fd->IsZombie());
  assert(fm && !fm->IsZombie());
  assert(fr && !fr->IsZombie());
  fr->mkdir("orig");
  fr->mkdir("flavor");
  fr->mkdir("pf");
  
  cout << "Merging files "<<fd->GetName()<<" and "<<fm->GetName() << endl;  
  cout << "Output file " << fr->GetName() << endl << flush;
  
  // Automatically go through the list of keys (profiles, histograms)
  TList *keys = fm->GetListOfKeys();
  TListIter itkey(keys);
  TObject *obj;
  TKey *key;
  
  bool debug = false;//true;
  while ( (key = dynamic_cast<TKey*>(itkey.Next())) ) {
    if (debug) cout << key->GetName() << endl << flush;
    obj = key->ReadObj(); assert(obj);

    fm->cd();
    if (obj->InheritsFrom("TProfile")) {

      TProfile *pm = (TProfile*)obj;
      TString mcname = obj->GetName();
      TString dataname = mcname;
      dataname.ReplaceAll("MC","DATA");
      dataname.ReplaceAll("_ps0",""); // v20
      dataname.ReplaceAll("_ps1",""); // v20
      dataname.ReplaceAll("_ps2",""); // v20
      dataname.ReplaceAll("_ps3",""); // v20
      TProfile *pd = (TProfile*)fd->Get(dataname.Data());
      if (!pd && dataname.Contains("Rho")) continue; // v19
      if (!pd && dataname.Contains("_ps")) continue; // v19
      assert(pd);

      fr->cd("orig");
      TH1D *hm = pm->ProjectionX(Form("h_%s",pm->GetName()));
      TH1D *hd = pd->ProjectionX(Form("h_%s",pd->GetName()));

      // Patch presumed Gain1 miscalibration in data
      if ((mcname.Contains("_MPFchs_") || mcname.Contains("_MPFR1chs_") ||
	   mcname.Contains("_PtBalchs_") || mcname.Contains("_DBchs_")) &&
	  addGain1) {
	const char *co = "control/pgain1vspt";
	TProfile *pd = (TProfile*)fd->Get(Form(co,"DATA")); assert(pd);
	//TProfile *pm = (TProfile*)fm->Get(Form(co,"MC")); assert(pm);
	hd->Add(pd,0.02);
      }

      // Patch presumed Photon vs PFphoton difference + footprint
      if ((mcname.Contains("_MPFchs_") || mcname.Contains("_MpfRuchs_")) &&
	  addFootprint) {
	const char *co = "control/pphoj0";
	TProfile *pd = (TProfile*)fd->Get(Form(co,"DATA")); assert(pd);
	TProfile *pm = (TProfile*)fm->Get(Form(co,"MC")); assert(pm);
	hd->Add(pd);
	hm->Add(pm);
      }

      TString rationame = mcname;
      rationame.ReplaceAll("_MC","");
      TH1D *hr = (TH1D*)hd->Clone(Form("h_%s",rationame.Data()));
      if (TString(key->GetName()).Contains("MPFR1chs") || 
	  TString(key->GetName()).Contains("MPFRnchs") || 
	  TString(key->GetName()).Contains("MpfRuchs"))
	hr->Add(hm,-1);
      else
	hr->Divide(hm);
      fm->cd();

      // Retrieve pT mapping
      char c[256];
      int a, y1, y2;
      if (debug) cout << "mcname: " << mcname << endl << flush;
      //sscanf(mcname.Data(),"resp_%[PtGamBalDBMPFpfR1nuchsRho]_MC_a%d_eta%d_%d",
      sscanf(mcname.Data(),"resp_%[PtGamBalDBMPFpfR1nuchsRhoDeltaJES]_MC_a%d_eta%d_%d",
	     c,&a,&y1,&y2);
      if (debug) cout << "c: " << c << endl << flush;
      TString mcptname = mcname;
      if (string(c)=="Rho") mcptname.ReplaceAll("Rho_CHS","PtGam");
      else mcptname.ReplaceAll(c,"PtGam");
      if (debug) cout << "PtGam(MC):" << mcptname << endl << flush;
      TString dataptname = dataname;
      if (string(c)=="Rho") dataptname.ReplaceAll("Rho_CHS","PtGam");
      else dataptname.ReplaceAll(c,"PtGam");
      if (debug) cout << "PtGam(data):" << dataptname << endl << flush;
      TProfile *pmpt = (TProfile*)fm->Get(mcptname.Data());   assert(pmpt);
      TProfile *pdpt = (TProfile*)fd->Get(dataptname.Data()); assert(pdpt);

      // Detailed calculation for data/MC average
      TH1D *hrpt = pmpt->ProjectionX("tmp");
      for (int i = 1; i != hrpt->GetNbinsX()+1; ++i) {

	double pt = 0.5*(pmpt->GetBinContent(i)+pdpt->GetBinContent(i));
	double ept = sqrt(pow(0.5*(pmpt->GetBinContent(i)-
				   pdpt->GetBinContent(i)),2) +
			  pow(pmpt->GetBinError(i),2) +
			  pow(pdpt->GetBinError(i),2));
	// Safeguards for empty bins in data or MC (e.g. due to trigger)
	if (pmpt->GetBinContent(i)==0) {
	  pt  = pdpt->GetBinContent(i);
	  ept = pdpt->GetBinError(i);
	}
	if (pdpt->GetBinContent(i)==0) {
	  pt  = pmpt->GetBinContent(i);
	  ept = pmpt->GetBinError(i);
	}
	hrpt->SetBinContent(i, pt);
	hrpt->SetBinError(i, ept);
      } // for i in hrpt

      fr->cd();
      TGraphErrors *gd = new TGraphErrors(hd);
      gd->SetNameTitle(dataname.Data(), dataname.Data());
      TGraphErrors *gm = new TGraphErrors(hm);
      gm->SetNameTitle(mcname.Data(), mcname.Data());
      TGraphErrors *gr = new TGraphErrors(hr);
      gr->SetNameTitle(rationame.Data(), rationame.Data());
      
      cleanGraph(gm);
      cleanGraph(gd);
      cleanGraph(gr);
      
      // Replace pT with more accurate estimate
      replacePt(gm,pmpt);
      replacePt(gd,pdpt);
      replacePt(gr,hrpt);
      delete hrpt;
      
      gm->SetLineColor(kRed);
      gd->SetLineColor(kBlue);
      gr->SetLineColor(kBlack);
      
      gm->Write(gm->GetName(),TObject::kOverwrite); 
      gd->Write(gd->GetName(),TObject::kOverwrite); 
      gr->Write(gr->GetName(),TObject::kOverwrite); 
    } // TProfile
    else if (obj->InheritsFrom("TH1D")) {

      TH1D *hm = (TH1D*)obj;
      TString mcname = obj->GetName();
      if (mcname.Contains("RawNEvents")) hm->Scale(750e3);

      TString dataname = mcname;
      dataname.ReplaceAll("MC","DATA");
      dataname.ReplaceAll("_ps0",""); // v20
      dataname.ReplaceAll("_ps1",""); // v20
      dataname.ReplaceAll("_ps2",""); // v20
      dataname.ReplaceAll("_ps3",""); // v20
      TH1D *hd = (TProfile*)fd->Get(dataname.Data());
      if (!hd && dataname.Contains("Xsec")) continue; // v19
      if (!hd && dataname.Contains("_ps")) continue; // v19
      if (!hd) {
	cout << "Missing " << dataname << endl << flush;
	if (dataname=="hxsec") continue; // run3 MC extra
	if (dataname=="hnevt") continue; // run3 MC extra
	if (dataname=="hsumw") continue; // run3 MC extra
	if (dataname=="hLHE_HT") continue; // run3 MC extra
	if (dataname=="hHT") continue; // run3 MC extra
      }
      assert(hd);

      TString rationame = mcname;
      rationame.ReplaceAll("_MC","");
      TH1D *hr = (TH1D*)hd->Clone(rationame.Data());
      hr->Divide(hm);

      fr->cd();
      //hm->SetDirectory(gDirectory);
      //hd->SetDirectory(gDirectory);
      //hr->SetDirectory(gDirectory);
      hm->Write(hm->GetName(),TObject::kOverwrite); 
      hd->Write(hd->GetName(),TObject::kOverwrite); 
      hr->Write(hr->GetName(),TObject::kOverwrite); 
    } // TH1D
  } // while key in itkey


  // Flavor responses in their own folder
  fm->cd("flavor");
  TDirectory *d1m = gDirectory;
  fd->cd("flavor");
  TDirectory *d1d = gDirectory;

  TList *fkeys = d1m->GetListOfKeys();
  TListIter itfkey(fkeys);
  
  while ( (key = dynamic_cast<TKey*>(itfkey.Next())) ) {
    if (debug) cout << key->GetName() << endl << flush;
    obj = key->ReadObj(); assert(obj);

    d1m->cd();
    if (obj->InheritsFrom("TProfile")) {
      TProfile *pm = (TProfile*)obj;
      TProfile *pd = (TProfile*)d1d->Get(key->GetName());
      if (!pd) continue; // v20 strange jets missing from data
      assert(pd);

      fr->cd("orig");

      TH1D *hm = pm->ProjectionX(Form("%s_mc",pm->GetName()));
      TH1D *hd = pd->ProjectionX(Form("%s_data",pd->GetName()));
      TH1D *hr = (TH1D*)hd->Clone(Form("%s_ratio",pd->GetName()));
      if (TString(key->GetName()).Contains("mpf1") || 
	  TString(key->GetName()).Contains("mpfn") || 
	  TString(key->GetName()).Contains("mpfu"))
	hr->Add(hm,-1);
      else
	hr->Divide(hm);

      fr->cd("flavor");

      TGraphErrors *gd = new TGraphErrors(hd);
      gd->SetNameTitle(hd->GetName(),hd->GetName());
      TGraphErrors *gm = new TGraphErrors(hm);
      gm->SetNameTitle(hm->GetName(),hm->GetName());
      TGraphErrors *gr = new TGraphErrors(hr);
      gr->SetNameTitle(hr->GetName(),hr->GetName());
      
      cleanGraph(gm);
      cleanGraph(gd);
      cleanGraph(gr);
      
      gm->SetLineColor(kRed);
      gd->SetLineColor(kBlue);
      gr->SetLineColor(kBlack);
      
      gm->Write(); 
      gd->Write();
      gr->Write();

      string s = key->GetName();
      if (patchGx && s[s.length()-1]=='i') { // ends with 'i' (e.g. counts_gii)
	string s2 = s.substr(0,s.length()-1); // remove final 'i' (counts_gi)
	gm->Write((s2+"_mc").c_str());
	gd->Write((s2+"_data").c_str());
	gr->Write((s2+"_ratio").c_str());
      }

      d1m->cd();
    } // TProfile
    else if (obj->InheritsFrom("TH1D")) {

      TH1D *hm = (TH1D*)obj;
      TH1D *hd = (TProfile*)d1d->Get(key->GetName());
      if (!hd) continue; // v20 strange jets missing from data
      assert(hd);
      TH1D *hr = (TH1D*)hd->Clone(Form("%s_ratio",key->GetName()));
      hr->Divide(hm);

      hm->SetName(Form("%s_mc",key->GetName()));
      hd->SetName(Form("%s_data",key->GetName()));

      fr->cd("flavor");
      hm->Write();
      hd->Write();
      hr->Write();

      string s = key->GetName();
      if (patchGx && s[s.length()-1]=='i') { // ends with 'i' (e.g. counts_gii)
	string s2 = s.substr(0,s.length()-1); // remove final 'i' (counts_gi)
	hm->Write((s2+"_mc").c_str());
	hd->Write((s2+"_data").c_str());
	hr->Write((s2+"_ratio").c_str());
      }
    } // TH1D
  } // while flavor


  // PF composition in their own folder
  fm->cd("pf");
  TDirectory *d2m = gDirectory;
  fd->cd("pf");
  TDirectory *d2d = gDirectory;

  TList *pkeys = d2m->GetListOfKeys();
  TListIter itpkey(pkeys);
  
  while ( (key = dynamic_cast<TKey*>(itpkey.Next())) ) {
    if (debug) cout << key->GetName() << endl << flush;
    obj = key->ReadObj(); assert(obj);

    d2m->cd();
    if (obj->InheritsFrom("TProfile2D")) {
      TProfile2D *p2m = (TProfile2D*)obj;
      TProfile2D *p2d = (TProfile2D*)d2d->Get(key->GetName());
      assert(p2d);

      fr->cd("orig");

      TH2D *h2m = p2m->ProjectionXY(Form("%s_mc",p2m->GetName()));
      TH2D *h2d = p2d->ProjectionXY(Form("%s_data",p2d->GetName()));
      TH2D *h2r = (TH2D*)h2d->Clone(Form("%s_ratio",p2d->GetName()));
      if (TString(key->GetName()).Contains("chf") || 
	  TString(key->GetName()).Contains("nhf") || 
	  TString(key->GetName()).Contains("nef") || 
	  TString(key->GetName()).Contains("cef") || 
	  TString(key->GetName()).Contains("muf") || 
	  TString(key->GetName()).Contains("puf"))
	h2r->Add(h2m,-1);
      else
	h2r->Divide(h2m);

      fr->cd("pf");

      h2m->Write(); 
      h2d->Write();
      h2r->Write();

      d2m->cd();
    }
    else if (obj->InheritsFrom("TProfile")) {
      TProfile *pm = (TProfile*)obj;
      TProfile *pd = (TProfile*)d2d->Get(key->GetName());
      assert(pd);

      fr->cd("orig");

      TH1D *hm = pm->ProjectionX(Form("%s_mc",pm->GetName()));
      TH1D *hd = pd->ProjectionX(Form("%s_data",pd->GetName()));
      TH1D *hr = (TH1D*)hd->Clone(Form("%s_ratio",pd->GetName()));
      if (TString(key->GetName()).Contains("chf") || 
	  TString(key->GetName()).Contains("nhf") || 
	  TString(key->GetName()).Contains("nef") || 
	  TString(key->GetName()).Contains("cef") || 
	  TString(key->GetName()).Contains("muf") || 
	  TString(key->GetName()).Contains("puf"))
	hr->Add(hm,-1);
      else
	hr->Divide(hm);

      fr->cd("pf");

      TGraphErrors *gd = new TGraphErrors(hd);
      gd->SetNameTitle(hd->GetName(),hd->GetName());
      TGraphErrors *gm = new TGraphErrors(hm);
      gm->SetNameTitle(hm->GetName(),hm->GetName());
      TGraphErrors *gr = new TGraphErrors(hr);
      gr->SetNameTitle(hr->GetName(),hr->GetName());
      
      cleanGraph(gm);
      cleanGraph(gd);
      cleanGraph(gr);
      
      gm->SetLineColor(kRed);
      gd->SetLineColor(kBlue);
      gr->SetLineColor(kBlack);
      
      gm->Write(); 
      gd->Write();
      gr->Write();

      d2m->cd();
    } // TProfile
    else if (obj->InheritsFrom("TH2D")) {

      TH2D *h2m = (TH2D*)obj;
      TH2D *h2d = (TH2D*)d2d->Get(key->GetName());
      assert(h2d);
      TH2D *h2r = (TH2D*)h2d->Clone(Form("%s_ratio",key->GetName()));
      h2r->Divide(h2m);

      h2m->SetName(Form("%s_mc",key->GetName()));
      h2d->SetName(Form("%s_data",key->GetName()));

      fr->cd("pf");
      h2m->Write();
      h2d->Write();
      h2r->Write();
    } // TH2D
  } // while pf composition


  // Rho plots in the control folder
  fm->cd("control");
  TDirectory *d3m = gDirectory;
  fd->cd("control");
  TDirectory *d3d = gDirectory;

  TList *rkeys = d3m->GetListOfKeys();
  TListIter itrkey(rkeys);
  
  while ( (key = dynamic_cast<TKey*>(itrkey.Next())) ) {
    if (debug) cout << key->GetName() << endl << flush;
    obj = key->ReadObj(); assert(obj);

    // Select the plots to keep further
    if (!(string(obj->GetName())=="prhovspt" ||
	  string(obj->GetName())=="prhovsmu" ||
	  string(obj->GetName())=="pnpvgoodvsmu" ||
	  string(obj->GetName())=="pnpvallvsmu"))
	continue;

    d3m->cd();
    if (obj->InheritsFrom("TProfile")) {
      TProfile *pm = (TProfile*)obj;
      TProfile *pd = (TProfile*)d3d->Get(key->GetName());
      assert(pd);

      fr->cd("orig");

      TH1D *hm = pm->ProjectionX(Form("%s_mc",pm->GetName()));
      TH1D *hd = pd->ProjectionX(Form("%s_data",pd->GetName()));
      TH1D *hr = (TH1D*)hd->Clone(Form("%s_ratio",pd->GetName()));
      hr->Divide(hm);

      fr->cd();

      TGraphErrors *gd = new TGraphErrors(hd);
      gd->SetNameTitle(hd->GetName(),hd->GetName());
      TGraphErrors *gm = new TGraphErrors(hm);
      gm->SetNameTitle(hm->GetName(),hm->GetName());
      TGraphErrors *gr = new TGraphErrors(hr);
      gr->SetNameTitle(hr->GetName(),hr->GetName());
      
      cleanGraph(gm);
      cleanGraph(gd);
      cleanGraph(gr);
      
      gm->SetLineColor(kRed);
      gd->SetLineColor(kBlue);
      gr->SetLineColor(kBlack);
      
      gm->Write(); 
      gd->Write();
      gr->Write();

      d3m->cd();
    } // TProfile
  } // while rho (control)

  cout << "Writing file " << fr->GetName() << endl << flush;
  fr->Write();
  fr->Close();
  cout << "File closed." << endl << flush;
} // GamHistosRatio

// J. Mnich, DG's talk on June 22, 2021
// Primary target of LHC: di-Higgs, ATLAS HH->gamma gamma b bar (4b?)
