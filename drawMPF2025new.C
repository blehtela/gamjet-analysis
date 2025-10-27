// My script created for plotting MPF for different minbiasXS used in pileup reweighting.
// Created by Bettina Lehtelä, March 2025.
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TLine.h"

#include "tdrstyle_mod22.C"
#include <iostream>

bool multiera = true; //draw all eras given in iov list
//string id = "w60"; //version of code used to produce the input files for this
string id = "w63"; 

// Forward declaration of function call
//void drawMPF2025new(string sobj, string var, string name, double x1, double x2); //, double z1, double z2); //start with only upper plot, will add lower one later
void drawMPF2025new(string sobj, string var, string name, double x1, double x2, double z1, double z2); //start with only upper plot, will add lower one later
//upper plot: distributions
//lower plot: reweighted MC / data (TO DO)
//double x1, double x2, double z1, double z2); //x1, x2 upper plot, z1, z2 lower plot

// Multiple calls to draw function
void drawMPF2025new() {
	/*
	drawMPF2025new("pf/pmpf","p_{T} in GeV","MPF", 15.,4000.,0.85,1.15); //20.,60.
	drawMPF2025new("pf/pdb","p_{T} in GeV","DB", 15.,4000.,0.7,1.3); //20.,60.
	*/
	drawMPF2025new("pf/pmpf","p_{T} in GeV","MPF", 15.,4000.,0.98,1.02); //20.,60.
	drawMPF2025new("pf/pdb","p_{T} in GeV","DB", 15.,4000.,0.97,1.03); //20.,60.

	//drawMPF2025new("pileup/h_npvgood","NPV_{good}","NPVgood",10.0,80.0);
	//drawMPF2025new("pileup/h_npvall","NPV_{all}","NPVall",10.0,80.0);
	//drawMPF2025new("pileup/h_mu","#mu","Mu",0.0,120.0); //only mc? could read it from data though, different file!
} // drawMPF2025new


//function to draw rho, mu, NPV distributions
//TO DO: modify this, to do it for several data eras at once..
//void drawMPF2025new(string sobj, string var, string name, double x1, double x2){
void drawMPF2025new(string sobj, string var, string name, double x1, double x2, double z1, double z2){
			     //double x1, double y2, double z1, double z2) {
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  //string iovs[] = {"2025Cv1", "2025Cv2", "2025D", "2025E", "2025F", "winter2025P8"};  //no qcd, no 25B
  string iovs[] = {"2025Cv1", "2025Cv2", "2025C-TrkRadDamage"};  //no qcd, no 25B



  const int niov = sizeof(iovs)/sizeof(iovs[0]); //but will be +1 for pu-off
  //string *iovs = &iovs[0];

  map<string,int> mcolor;
  map<string,int> mmarker;

  //MARKERS
  //2025
  mmarker["winter2025P8"] = kOpenCircle; //NOT reweighted
  mmarker["2025Cv1"] = kFullTriangleUp;
  mmarker["2025Cv2"] = kFullTriangleDown;
  mmarker["2025C-TrkRadDamage"] = kFullSquare;

  mmarker["2025D"] = kFullCircle;
  mmarker["2025E"] = kFullSquare;
  mmarker["2025F"] = kFullStar;
 

  //COLOURS
  //2025
  mcolor["winter2025P8"] = kGreen+2; //NOT reweighted
  mcolor["2025Cv1"] = kOrange+1;
  mcolor["2025Cv2"] = kPink+1;
  mcolor["2025C-TrkRadDamage"] = kTeal+5; //kGreen+1;
  mcolor["2025D"] = kBlue+1;
  mcolor["2025E"] = kGreen+1;
  mcolor["2025F"] = kCyan+1;
 
  
  const char *cvar = var.c_str(); //x-variable
  const char *cname = name.c_str();


  //TH1D *h = tdrHist("h",xlabel=cvar, ylabel="N_{events}",x1,x2); //axis range setting does not work here properly...
  TH1D *h = new TH1D("h", Form(";%s;%s",cvar,cname), x2-x1, x1, x2);  // N_{events}/N_{events-total}

  h->GetXaxis()->SetRangeUser(x1,x2);
  if(name=="DB"){
    h->GetYaxis()->SetRangeUser(0.85,1.15); //set this depending on the type of plot
  }
  else if(name=="MPF"){
    //h->GetYaxis()->SetRangeUser(0.9,1.1); //set this depending on the type of plot
    h->GetYaxis()->SetRangeUser(0.93,1.1); //set this depending on the type of plot
  }
  //TH1D *h2 = tdrHist("h2","data/MC",z1,z2); //TO DO
  TH1D *h2 = tdrHist("h2","data/data Cv1",z1,z2); //TO DO: JUST FOR NOW COMPARE WITH Cv1


  lumi_13TeV = "Run3"; // 4=13 TeV
  if (id!="") lumi_136TeV = Form("Run3 %s",id.c_str()); // 8=13.6 TeV
  TCanvas *c1 = tdrDiCanvas(Form("c1_%s",cname),h,h2,8,11); //when adding 2nd hist
  //TCanvas *c1 = tdrCanvas(Form("c1_%s",cname),h,8,11); //when only one

	//custom place for "Private" extraText (usually handled via tdrCanvas)
	TLatex *extratext = new TLatex();
  extratext->SetNDC();
	extratext->SetTextSize(0.039);
	extratext->SetTextFont(52);
	extratext->DrawLatex(0.3, 0.9, "Preliminary");
	double textpos = 0.00;


	//Add information to warn that the average is that of the histogram (still includes the hard primary vertex itself)
	TLatex *avginfo = new TLatex();
	avginfo->SetNDC(); //switch relativ coord. on
	avginfo->SetTextSize(0.031);
	avginfo->SetTextFont(42);
	avginfo->SetTextColor(kPink+7);
	//avginfo->DrawLatex(0.46, 0.87, "(mean #pm mean_err) still includes the primary vertex!"); //not for MPF plot



  c1->cd(1);
  TLine *l = new TLine(); //only for ratio
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray);
  gPad->SetLogx(); //not here

  //TLegend *leg = tdrLeg(0.46,0.850-niov*0.06,0.80,0.850); //after adding also average-info (x1,y1,x2,y2)
  double legx1 = 0.46;
  double legx2 = 0.80;
  //double legy1 = 0.850-niov*0.05;
  double legy2 = 0.85;

  if(name=="MPF"){
    legx1=0.3;
    //legx2=;
    //legy1=;
    legy2=0.35;
  }
  TLegend *leg = tdrLeg(legx1,legy2-niov*0.05,legx2,legy2);
  //TLegend *leg = tdrLeg(0.46,0.850-niov*0.05,0.80,0.850); //after adding also average-info (x1,y1,x2,y2)
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);


  //Ratio plot (data over MC) - TO DO
  c1->cd(2);
  gPad->SetLogx();
  l->DrawLine(h->GetXaxis()->GetXmin(),1,h->GetXaxis()->GetXmax(),1); //adjust!
  //l->DrawLine(h->GetXaxis()->GetXmin(),0,h->GetXaxis()->GetXmax(),0);

  for (int i = 0; i != niov; ++i) { //go through iovs, this does NOT include mc prior to reweighting. (note: should change this... only look at one era at a time)

    string iov = iovs[i];
    const char *ciov = iov.c_str(); //era
    //const char *cmc = mcs[i].c_str();
    const char *cid = id.c_str(); //code version, e.g. w38
    //const char *cdataera1 = dataera1.c_str();
    //const char *cdataera2 = dataera2.c_str();
    //const char *cmcera = mcera.c_str();


    TFile *fd(0), *fm(0);

    //draw this iov's distributions (four histos, for data mu from extra file?)
    //if(iov=="summer2024P8"){// Monte Carlo
    if(iov=="winter2025P8"){// Monte Carlo
    	fd = new TFile(Form("rootfiles/GamHistosFill_mc_%s_no-pu_w56.root",ciov)); //rootfiles/GamHistosFill_mc_winter2025P8_no-pu_w56.root
    	assert(fd && !fd->IsZombie());
    }
    else{ //data
    	fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",ciov,cid));
    }

    //always same MC for now
    //fm = new TFile(Form("rootfiles/GamHistosFill_mc_winter2025P8_no-pu_w56.root")); //rootfiles/GamHistosFill_mc_winter2025P8_no-pu_w56.root
    fm = new TFile(Form("rootfiles/GamHistosFill_data_2025Cv1_%s.root",cid)); // TO DO: JUST FOR NOW Cv1
 
 
    assert(fd && !fd->IsZombie());
    assert(fm && !fm->IsZombie());
    curdir->cd();
    
    cout << "Trying to get sobj.c_str() = " << sobj.c_str() << " for era " << ciov << " version " << cid  << endl << flush;
    TObject *ohisto = fd->Get(Form("%s",sobj.c_str())); assert(ohisto); //get histo object
    TObject *omc = fm->Get(Form("%s",sobj.c_str())); assert(omc); //get histo object
 
    
    TH1D *hd(0), *hm(0);
    hd = (TH1D*)ohisto;
    hm = (TH1D*)omc;
    assert(hd);
    assert(hm);
    cout << "\nGot histo from file " << fd << Form("(%s)",sobj.c_str()) << endl << flush;
    
    //normalise with its integral
    //hd->Scale(1./hd->Integral());
    
    double ymax = hd->GetMaximum(); //for setting correct y-range (not an optimal solution, but works for now)
    //h->GetYaxis()->SetRangeUser(0., ymax+0.2*ymax); //use this when investigating which ymax to set, otherwise: fixed in mymax.
    //if(ymax>textpos){textpos=ymax+0.1*ymax;}
    //extratext->DrawLatex((x2-x1)*0.2, ymax, "Preliminary"); //reposition text
    
    double xmean = hd->GetMean(1); //get the average along x-axis
    double xmeanerr = hd->GetMeanError(1); //get mean error
    //cout << Form("Mean for %s is ",sobj.c_str()) << xmean << endl << flush;

    //new: smaller marker
    //hd->SetMarkerSize(2);

    //RATIO
    TH1D *hr = (TH1D*)hd->Clone(Form("hr_%s_%s", cvar, ciov));
    hm->Sumw2(1);
    hr->Sumw2(1); //preserve err?
    hr->Divide(hm); //do the actual ratio, 1.0 for MC

    //NORMAL PLOT (no ratio)
    c1->cd(1);
    tdrDraw(hd,"Pz",(mmarker[iov] ? mmarker[iov] : kFullCircle), (mcolor[iov] ? mcolor[iov] : kBlack));
    if(name!="Mu" && iov!="winter2025P8"){ //for mu, do not add the data here, as this histo is empty, i.e. not drawn, will come from separate file (see lower)
	    //leg->AddEntry(hd,Form("%s (%.2f #pm %.2f)",ciov,xmean, xmeanerr),"PLE");
	    leg->AddEntry(hd,Form("%s data",ciov),"PLE");
    }
    else if(name=="Mu" && iov=="summer2024P8"){ //need to add this manually
	    //leg->AddEntry(hd,Form("%s (%.2f #pm %.2f)",ciov,xmean, xmeanerr),"PLE");
	    leg->AddEntry(hd,Form("%s (75.3mb)",ciov),"PLE");
    }
    else if(name!="Mu" && iov=="winter2025P8"){ //need to add this manually
	    //leg->AddEntry(hd,Form("%s (%.2f #pm %.2f)",ciov,xmean, xmeanerr),"PLE");
	    leg->AddEntry(hd,Form("%s MC (no pu-reweighting)",ciov),"PLE");
    }


    //additionally if mc
    /*
    if(iov=="summer2024P8"){
	//add for other minbiasXS manually (remove this part again when case is closed)
	fmalt = new TFile(Form("rootfiles/w48_xs69200/GamHistosFill_mc_%s_pu-%s-xs69200_%s.root",ciov,cdataera,cid)); //adjust minbiasXS if needed!!
	
	//hmpualt->Scale(1./hmpualt->Integral()); //normalise with its integral
	double xmean_pualt = hmpualt->GetMean(1); //get the average along x-axis
	double xmeanerr_pualt = hmpualt->GetMeanError(1); //get mean error
	

    	c1->cd(1);
    	//tdrDraw(hmpualt,"Pz", kOpenCircle, kOrange+2);
    	tdrDraw(hmpualt,"Pz", kOpenCircle, kBlue+2);
    	leg->AddEntry(hmpualt,Form("%s (69.2mb)",iov.c_str()),"PLE");


	//the not-reweighted MC (keep this always)
	TObject *ohistopuoff = fmpuoff->Get(Form("%s",sobj.c_str())); assert(ohistopuoff);
	hmpuoff = (TH1D*)ohistopuoff;
	assert(hmpuoff);
	cout << "Got histo from file " << fmpuoff << endl << flush;
	
	//hmpuoff->Scale(1./hmpuoff->Integral()); //normalise with its integral
	double xmean_puoff = hmpuoff->GetMean(1); //get the average along x-axis
	double xmeanerr_puoff = hmpuoff->GetMeanError(1); //get mean error
	

    	c1->cd(1);
    	tdrDraw(hmpuoff,"Pz", kOpenCircle, kRed+2);
    	leg->AddEntry(hmpuoff,Form("%s (not reweighted)",iov.c_str()),"PLE");


    }
    */ 
 

    c1->cd(2); //for ratio plot
    gPad->SetLogx();
    //tdrDraw(hr,"Pz",(mmarker[iov] ? mmarker[iov] : kFullCircle), (mcolor[iov] ? mcolor[iov] : kBlack));//ERRORS
    tdrDraw(hr,"HIST P",(mmarker[iov] ? mmarker[iov] : kFullCircle), (mcolor[iov] ? mcolor[iov] : kBlack));//ADD CORRECT ERRORS

  } // for each iov

	//extratext->DrawLatex((x2-x1)*0.2, textpos, "Preliminary"); //reposition text




    //Now still adding the mu distribution for data (different file)
    //these were created with brilcalc
    /*
    if(name=="Mu"){
	c1->cd(1);
    	TFile *fdatapu(0);
	//open the manually (well, with my own script) created file containing all currently used pu distributions (i.a. MyDataPileupHisto...)
   	//fdatapu = new TFile(Form("pileup/pileup_mc_data%s_%s.root",cdataera,id.c_str()), "READ"); //last without reweighting so far, w36
   	////fdatapu = new TFile(Form("pileup/2024/pu_summary_xs69200_%s.root",id.c_str()), "READ"); //REMEMBER TO UPDATE MINBIASXS
   	fdatapu = new TFile(Form("pileup/2024/pu_summary_xs75300_%s.root",id.c_str()), "READ"); //REMEMBER TO UPDATE MINBIASXS


	assert(fdatapu && !fdatapu->IsZombie());
	TH1D *hdatapu = (TH1D*)fdatapu->Get(Form("pileup_%s_HLT_Photon50EB_TightID_TightIso", cdataera));  //only trigger i look at right now
	//pileup_2024Ev2nib1_HLT_Photon50EB_TightID_TightIso_scaled --> should i use scaled version instead?, currently scaling here (see below)
	//TH1D *hdatapu = (TH1D*)fdatapu->Get(Form("pileup_%s_%s",ce,ct));  //with automatic reading of era (ce) and trigger (ct)
	//should I do fdatapu->Clone() instead?
	assert(hdatapu); 
	//hdatapu->Scale(1./hdatapu->Integral()); //normalise to unity
	double xmean_data = hdatapu->GetMean(1); //get the average along x-axis
	double xmeanerr_data = hdatapu->GetMeanError(1); //get mean error


	//tdrDraw(hdatapu,"Pz",kFullTriangleDown,kBlue+2,kSolid,-1,kNone);
	tdrDraw(hdatapu,"Pz",kFullTriangleDown,kBlack,kSolid,-1,kNone);
	leg->AddEntry(hdatapu,Form("%s (from brilcalc); %.2f #pm %.2f",cdataera, xmean_data, xmeanerr_data),"PLE");
    }
    */

    gPad->SetLogx();


    //finish everything.
    if (id!="")
	//c1->SaveAs(Form("pdf/drawMPF2025new_%s_%s_%s.pdf",cdataera,name.c_str(),id.c_str()));
	//c1->SaveAs(Form("pdf/drawMPF2025new_TEST_%s_%s.pdf",name.c_str(),id.c_str()));
	c1->SaveAs(Form("pdf/drawMPF2025new_compare-Cv1-Cv2-CTrkRadDamage_%s_%s.pdf",name.c_str(),id.c_str()));
    else 
	c1->SaveAs(Form("pdf/drawMPF2025new_compare-Cv1-Cv2-CTrkRadDamage_%s.pdf",name.c_str()));
} // void drawMPF2025new


