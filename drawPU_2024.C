//My script created for pileup investigations during autumn 2024.
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TLine.h"

#include "tdrstyle_mod22.C"

bool multiera = true; //draw all eras given in iov list
string id = "w38"; //version of code used to produce the input files for this


// Forward declaration of function call
void drawPU_2024(string sobj, string var, string name, double x1, double x2); //, double z1, double z2); //start with only upper plot, will add lower one later
//upper plot: distributions
//lower plot: reweighted MC / data (TO DO)
//double x1, double x2, double z1, double z2); //x1, x2 upper plot, z1, z2 lower plot

// Multiple calls to draw function
void drawPU_2024() {
	drawPU_2024("pileup/h_rho","#rho in GeV","Rho", 0.,80.); //20.,60.
	drawPU_2024("pileup/h_npvgood","NPV_{good}","NPVgood",10.0,80.0);
	drawPU_2024("pileup/h_npvall","NPV_{all}","NPVall",10.0,80.0);
	drawPU_2024("pileup/h_mu","#mu","Mu",0.0,120.0); //only mc? could read it from data though, different file!
} // drawPU_2024


//function to draw rho, mu, NPV distributions
//TO DO: modify this, to do it for several data eras at once..
void drawPU_2024(string sobj, string var, string name, double x1, double x2){
			     //double x1, double y2, double z1, double z2) {
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //string iovs[] = {"2024E", "2024F", "winter2024P8"};  // soon!! 
	//string iovs[] = {"2024F", "winter2024P8"}; 
	//string iovs[] = {"2024E", "winter2024P8"}; 
	string iovs[] = {"2024D", "winter2024P8"}; 
	string dataera = "2024D";
	string mcera = "winter2024P8";

	//string reweighted_mc = "w37";

  const int niov = sizeof(iovs)/sizeof(iovs[0]); //but will be +1 for pu-off
  //string *iovs = &iovs[0];

  map<string,int> mcolor;
  map<string,int> mmarker;
  //2024
	mmarker["2024E"] = kFullSquare;
	mmarker["2024F"] = kFullTriangleUp;
  mmarker["winter2024P8"] = kOpenCircle;

 	mcolor["2024E"] = kCyan+2;
  mcolor["2024F"] = kBlue+2;

 	mcolor["2024B-ECALRATIO"] = kBlue+3;
  mcolor["2024C-ECALRATIO"] = kGreen+3;
  mcolor["2024C-ECALR-HCALDI"] = kGreen+4;
	mcolor["2024C-ECALCC-HCALDI"] = kGreen+4;

  mcolor["winter2024P8"] = kGreen+2; //reweighted

	//for non-reweighted MC, will use kRed+2, but manually as has different code version (w38puoff)


	//different ymax used for plotting different variables, depending on name (but fixed to easily compare eras)
	map<string,double> mymax;
	mymax["Rho"] = 0.06;
	mymax["NPVgood"] = 0.06;
	mymax["NPVall"] = 0.06;
	mymax["Mu"] = 0.07;

 
  
  const char *cvar = var.c_str(); //x-variable
  const char *cname = name.c_str();
  const char *cdataera = dataera.c_str();
  const char *cmcera= mcera.c_str();


  //TH1D *h = tdrHist("h",xlabel=cvar, ylabel="N_{events}",x1,x2); //axis range setting does not work here properly...
	TH1D *h = new TH1D("h", Form(";%s;n_{events}",cvar), x2-x1, x1, x2);  // N_{events}/N_{events-total}

	h->GetXaxis()->SetRangeUser(x1,x2);
	//h->GetYaxis()->SetRangeUser(0.,1.); //set this later by getting the maximum
	h->GetYaxis()->SetRangeUser(0.,mymax[name]); //set this depending on the type of plot
  //TH1D *h2 = tdrHist("h2","MC (rew.) / data ",z1,z2); //TO DO

  lumi_13TeV = "Run3"; // 4=13 TeV
  if (id!="") lumi_136TeV = Form("Run3 %s",id.c_str()); // 8=13.6 TeV
  //TCanvas *c1 = tdrDiCanvas(Form("c1_%s",cname),h,h2,8,11); //when adding 2nd hist
  TCanvas *c1 = tdrCanvas(Form("c1_%s",cname),h,8,11);

	//custom place for "Private" extraText (usually handled via tdrCanvas)
	TLatex *extratext = new TLatex();
	extratext->SetTextSize(0.039);
	extratext->SetTextFont(52);
	//extratext->DrawLatex(60, 56, "Preliminary");
	double textpos = 0.00;


	//Add information to warn that the average is that of the histogram (still includes the hard primary vertex itself)
	TLatex *avginfo = new TLatex();
	avginfo->SetNDC(); //switch relativ coord. on
	avginfo->SetTextSize(0.031);
	avginfo->SetTextFont(42);
	avginfo->SetTextColor(kPink+7);
	//avginfo->DrawLatex(40, 0.053, "Avg still includes the primary vertex!");
	avginfo->DrawLatex(0.46, 0.87, "(mean #pm mean_err) still includes the primary vertex!");




  c1->cd(1);
  //TLine *l = new TLine(); //maybe later, not now, only for ratio
  //gPad->SetLogx(); //not here

  //TLegend *leg = tdrLeg(0.60,0.90-niov*0.08,0.80,0.90); //this is good for fewer graphs
	//TLegend *leg = tdrLeg(0.50,0.90-niov*0.03,0.60,0.90); //moved a bit to the left and smaller font
  TLegend *leg = tdrLeg(0.46,0.865-niov*0.06,0.80,0.865); //after adding also average-info (x1,y1,x2,y2)

  leg->SetTextSize(0.03);
  leg->SetTextFont(42);


	//Ratio plot (reweighted MC over data) - TO DO
/*
  c1->cd(2);
  gPad->SetLogx();
  l->DrawLine(h->GetXaxis()->GetXmin(),1,h->GetXaxis()->GetXmax(),1);
  l->DrawLine(h->GetXaxis()->GetXmin(),0,h->GetXaxis()->GetXmax(),0);
*/

  for (int i = 0; i != niov; ++i) { //go through iovs, this does NOT include mc prior to reweighting. (note: should change this... only look at one era at a time)

    string iov = iovs[i];
    const char *ciov = iov.c_str(); //era
    //const char *cmc = mcs[i].c_str();
    const char *cid = id.c_str(); //code version, e.g. w38
		const char *cdataera = dataera.c_str();
		const char *cmcera = mcera.c_str();


    TFile *fd(0), *fmpuoff(0);

		//draw this iov's distributions (four histos, for data mu from extra file?)
		if(iov=="winter2024P8"){// Monte Carlo
    	//fd = new TFile(Form("rootfiles/GamHistosFill_mc_%s_%s.root",ciov,ciov,cid));
			//get file like this: GamHistosFill_mc_winter2024P8_2024E-pu_w38.root
			fd = new TFile(Form("rootfiles/winter2024p8_reweighting-with-%s/GamHistosFill_mc_%s_%s-pu_%s.root",cdataera,ciov,cdataera,cid)); //different directory for each data-era used
			fmpuoff = new TFile(Form("rootfiles/GamHistosFill_mc_%s_%spuoff.root",ciov,cid)); //same era, but without pu reweighting
			assert(fd && !fmpuoff->IsZombie());
		}
		else{ //data
    	fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",ciov,cid));
		}
 
    assert(fd && !fd->IsZombie());
    curdir->cd();
    
    TObject *ohisto = fd->Get(Form("%s",sobj.c_str())); assert(ohisto); //get histo object
    
    TH1D *hd(0), *hmpuoff(0);
    hd = (TH1D*)ohisto;
    assert(hd);
		cout << "\nGot histo from file " << fd << Form("(%s)",sobj.c_str()) << endl << flush;

		//normalise with its integral
		hd->Scale(1./hd->Integral());

		double ymax = hd->GetMaximum(); //for setting correct y-range (not an optimal solution, but works for now)
		//h->GetYaxis()->SetRangeUser(0., ymax+0.2*ymax); //use this when investigating which ymax to set, otherwise: fixed in mymax.
		//if(ymax>textpos){textpos=ymax+0.1*ymax;}
		//extratext->DrawLatex((x2-x1)*0.2, ymax, "Preliminary"); //reposition text

		double xmean = hd->GetMean(1); //get the average along x-axis
		double xmeanerr = hd->GetMeanError(1); //get mean error
		//cout << Form("Mean for %s is ",sobj.c_str()) << xmean << endl << flush;


    c1->cd(1);
    tdrDraw(hd,"Pz",(mmarker[iov] ? mmarker[iov] : kFullCircle), (mcolor[iov] ? mcolor[iov] : kBlack));
		if(name!="Mu"){ //for mu, do not add the data here, as this histo is empty, i.e. not drawn, will come from separate file (see lower)
			leg->AddEntry(hd,Form("%s (%.2f #pm %.2f)",ciov,xmean, xmeanerr),"PLE");
		}
		else if(name=="Mu" && iov=="winter2024P8"){ //need to add this manually
			leg->AddEntry(hd,Form("%s (%.2f #pm %.2f)",ciov,xmean, xmeanerr),"PLE");
		}


		//additionally if mc
		if(iov=="winter2024P8"){
    	TObject *ohistopuoff = fmpuoff->Get(Form("%s",sobj.c_str())); assert(ohistopuoff);
			hmpuoff = (TH1D*)ohistopuoff;
			assert(hmpuoff);
			cout << "Got histo from file " << fmpuoff << endl << flush;

			hmpuoff->Scale(1./hmpuoff->Integral()); //normalise with its integral
			double xmean_puoff = hmpuoff->GetMean(1); //get the average along x-axis
			double xmeanerr_puoff = hmpuoff->GetMeanError(1); //get mean error
	

    	c1->cd(1);
    	tdrDraw(hmpuoff,"Pz", kOpenCircle, kRed+2);
    	leg->AddEntry(hmpuoff,Form("%s (not reweighted; %.2f #pm %.2f)",iov.c_str(), xmean_puoff, xmeanerr_puoff),"PLE");
		}
    
 

    c1->cd(2); //for ratio plot
    //gPad->SetLogx();
  } // for each iov

	//extratext->DrawLatex((x2-x1)*0.2, textpos, "Preliminary"); //reposition text




		//Now still adding the mu distribution for data (different file)
		//these were created with brilcalc
		if(name=="Mu"){
			c1->cd(1);
    	TFile *fdatapu(0);
			//open the manually (well, with my own script) created file containing all currently used pu distributions (i.a. MyDataPileupHisto...)
   		fdatapu = new TFile(Form("pileup/pileup_mc_data%s_%s.root",cdataera,id.c_str()), "READ"); //last without reweighting so far, w36
			assert(fdatapu && !fdatapu->IsZombie());
			TH1D *hdatapu = (TH1D*)fdatapu->Get(Form("pileup_%s_HLT_Photon50EB_TightID_TightIso", cdataera));  //only trigger i look at right now
			//TH1D *hdatapu = (TH1D*)fdatapu->Get(Form("pileup_%s_%s",ce,ct));  //with automatic reading of era (ce) and trigger (ct)
			//should I do fdatapu->Clone() instead?
			assert(hdatapu); 
			hdatapu->Scale(1./hdatapu->Integral()); //normalise to unity
			double xmean_data = hdatapu->GetMean(1); //get the average along x-axis
			double xmeanerr_data = hdatapu->GetMeanError(1); //get mean error
	


			tdrDraw(hdatapu,"Pz",kFullTriangleDown,kBlue+2,kSolid,-1,kNone);
			leg->AddEntry(hdatapu,Form("%s (from brilcalc); %.2f #pm %.2f",cdataera, xmean_data, xmeanerr_data),"PLE");
		}



  //finish everything.
  if (id!="")
		c1->SaveAs(Form("pdf/drawPU_%s_%s_%s.pdf",cdataera,name.c_str(),id.c_str()));
  else 
    c1->SaveAs(Form("pdf/drawPU_%s_%s.pdf",cdataera,name.c_str()));
} // void drawPU_2024


