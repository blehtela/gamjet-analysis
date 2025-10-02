// Purpose: draw photon+jet response vs run to check time stability
// This is the code for plotting 2025 only... could add such a switch to the general code.
#include "TFile.h"
#include "TProfile.h"
#include "TLatex.h"

#include "TArrow.h" //for some reason started to throw error with forward declaration

#include "../tdrstyle_mod22.C"

// Clean
void clean(TH1 *p, double maxerr = 0.005) {
  
  for (int i = 1; i != p->GetNbinsX()+1; ++i) {
    if (p->GetBinError(i)==0 || p->GetBinError(i)>maxerr) {
      p->SetBinContent(i, 0);
      p->SetBinError(i, 0);
    }
  }
} // clean

// Shift the x-axis position for higher run numbers, to avoid whitespace -- is this doing what it should?
void shiftX(TH1 *p, double maxrunnum = 378000){
    TAxis *a = p->GetXaxis();
    cout << "this is xmin: " << a->GetXmin();
    //cout << "getsize: " << a->GetXbins()->GetSize();
    //if(a->GetXmin()>maxrunnum){ //doesn't work, because they all start at 355000...

    if(a->GetXbins()->GetSize()){ //axis with variable bins
        cout << a << endl; 
        TArrayD X(*(a->GetXbins()));
        for(Int_t i=0; i<X.GetSize(); i++){
            X[i]-=3000; //shift the x-value 3000 to the left
        }

        //set new xbins
        a->Set((X.GetSize()-1), X.GetArray());
    }
    else{ //axis with fix bins
        a->Set(a->GetNbins(),(a->GetXmin())-3000,(a->GetXmax())-3000);
    }
    //}
}//shiftX

// Add offset to histogram, if bin content nonempty
void addOffset(TH1D *h, double off) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    if (h->GetBinContent(i)!=0 && h->GetBinError(i)!=0) {
      h->SetBinContent(i, h->GetBinContent(i) + off);
    }
  }
} // addOffset

TH1D *hadd(string name, TProfile *p1, TProfile *p2) {
  TH1D *h = p1->ProjectionX(name.c_str());
  for (int i = 1; i != p1->GetNbinsX()+1; ++i) {
    if (p2->GetBinContent(i)!=0) {
      h->SetBinContent(i, p2->GetBinContent(i));
      h->SetBinError(i, p2->GetBinError(i));
    }
  }
  return h;
} // hadd


void drawPFcompVsRun_2025only(string version, int rereco, bool zoom);
void drawPFcompEFBvsRun_2025only(string version, int rereco, bool zoom);



//void drawResponseVsRun_custom(string version = "w10") { //w10 for 2023 data, w11 for the new 2024 data
//void drawResponseVsRun_custom(string version = "w12", string year=2024) {//for plotting only one year
//void drawResponseVsRun_2025only(string version = "w39w40", int rereco = 0) { //switched to w15 and w16; w17 and w18; w19 and w20
//void drawResponseVsRun_2025only(string version = "w56", int rereco = 0) { //w48 added on 06.04.2025 (previous w44 was mid-Feb '25), rereco = 3 for CDE-rereco with FGHI-prompt
//void drawResponseVsRun_2025only(string version = "w58", int rereco = 0) { //w58 added on 12.08.2025, 25C and early 25D. Data until 6th of Aug approx.
void drawResponseVsRun_2025only(string version = "w62", int rereco = 0) { //w62 (30.09., processed a week earlier or so)

    //bool rereco = 1;

    //const char *cyear = year.c_str();
    const char *cv = version.c_str();

    setTDRStyle();
    TDirectory *curdir = gDirectory;

    // Open input files
    //TFile *f = new TFile(Form("rootfiles/GamHistosFill_data_Run3_%s.root",cv), "READ"); //is this the entire run 3 data with all runperiods?
    //TFile *f = new TFile(Form("rootfiles/GamHistosFill_data_Run3Summer23_%s.root",cv), "READ"); //for now: hadd on 2023 Cv123, Cv4, D
    
    TFile *f(0);
    
    if(rereco==0){
	cout << "Looking at 2025 prompt data." << endl;
	f = new TFile(Form("rootfiles/GamHistosFill_data_2025CDEF_%s.root",cv), "READ");
	//f = new TFile(Form("rootfiles/GamHistosFill_data_2025BC_%s.root",cv), "READ");
	//f = new TFile(Form("rootfiles/GamHistosFill_data_2025BC_noL2L3Res_%s_28May2025.root",cv), "READ");
    } 
    //if and else for future when we might have rereco cases etc.
/*
    else if(rereco==1){
	cout << "Looking at 2025 rereco (ECALRATIO) data." << endl;
    	f = new TFile(Form("rootfiles/GamHistosFill_data_2024-ECALRATIO_%s.root",cv), "READ"); //for now: hadd on 2024B and 2024C (need to redo this file with new daily json)
    }
    else if(rereco==2){
	cout << "Looking at 2024C 2nd rereco (ECALR-HCALDI) data." << endl;
    	f = new TFile(Form("rootfiles/GamHistosFill_data_2024C-ECALR-HCALDI_%s.root",cv), "READ"); //for now: hadd on 2024B and 2024C (need to redo this file with new daily json)
    }
    else if(rereco==3){
	cout << "Looking at 2024CDE-rereco and 2024FGHI-prompt data." << endl;
	f = new TFile(Form("rootfiles/GamHistosFill_data_2024CDE-rereco_2024FGHI-prompt_%s.root",cv), "READ"); //containing 2024CDE-rereco and 2024FGHI-prompt
    }
*/
    else{
	cout << "Check rereco int. Something's wrong." << endl;
	exit(0);
    }

    assert(f && !f->IsZombie());

    //f->cd("control"); // v29
    f->cd("runs"); // v30
    TDirectory *d = gDirectory;

    // Load input profiles + clean the errors (what should i set the maxerr to? tried to be generous now)
    double maxerr = 0.01; //apply same maxerr cleaning to all paths, so it is less confusing
    //double maxerr = 0.1; //TEST (12.02.2025)
    //TProfile *pr30b = (TProfile*)d->Get("pr30b"); clean(pr30b,maxerr); //b = balance
    //TProfile *pr30m = (TProfile*)d->Get("pr30m"); clean(pr30m,maxerr); //m = mpf
    TProfile *pr40b = (TProfile*)d->Get("pr40b"); clean(pr40b,maxerr);  //new in 25
    TProfile *pr40m = (TProfile*)d->Get("pr40m"); clean(pr40m,maxerr);  //new in 25
    TProfile *pr50b = (TProfile*)d->Get("pr50b"); clean(pr50b,maxerr); 
    TProfile *pr50m = (TProfile*)d->Get("pr50m"); clean(pr50m,maxerr); 
    //TProfile *pr110b = (TProfile*)d->Get("pr110b"); clean(pr110b,maxerr);
    TProfile *pr110m = (TProfile*)d->Get("pr110m"); clean(pr110m,maxerr);
    //TProfile *pr230b = (TProfile*)d->Get("pr230b"); //clean(pr230b,0.003);
    //TProfile *pr230m = (TProfile*)d->Get("pr230m"); //clean(pr230m,0.003);



 
    // Scale BAL and pT30 results
    double kbal = 1.13;
    //pr30b->Scale(kbal);
    //double kbal50 = 1.10; //first had this, Mikko suggests 1.107 to make it closer to the other ones
    double kbal50 = 1.10; //1.08; //1.15; //1.11 earlier before w29
    double kbal40 = 1.05; //1.03; //1.08;
    pr50b->Scale(kbal50); //now also scaling this
    pr40b->Scale(kbal40);

    double kmpf50 = 1.00; //1.06; 
    double kmpf110 = 1.00; //1.06; 
    pr50m->Scale(kmpf50);
    pr110m->Scale(kmpf110);


/*
    double kpt50 = 1.00; //0.99 earlier before w29
    pr50m->Scale(kpt50);
*/
    //pr110b->Scale(kbal);
    //double kpt30 = 0.98;
    //pr30b->Scale(kpt30);
    //pr30m->Scale(kpt30);
    //double kbal30 = 0.96;
    //pr30b->Scale(kbal30);
/*
    if (pr50b && pr50m) { //in 2024 always
        //pr50b->Scale(kbal);
        double kpt50 = 0.99;//0.98;
        //pr50b->Scale(kpt50);
        pr50m->Scale(kpt50);
        //double kbal50 = 0.91;//kbal30;
        //pr50b->Scale(kbal50);
    }
*/
    //NOTE: DO WE ALSO NEED TO SCALE 230GeV results? (in case we include those)

    // Setup canvas
    //TH1D *h = tdrHist("h","Response",0.92,1.08,"Run",366000,381300); //should start later when dropping 2022 data (here some before 2023C)
    ///TH1D *h = tdrHist("h","Response",0.92,1.06,"Run",378900,382000); //2024 only (standard was to display from 0.92 to 1.08, now zooming in for 2024
    //TH1D *h = tdrHist("h","Response",0.92,1.06,"Run",378900,387100); //2024 only (standard was to display from 0.92 to 1.08, now zooming in for 2024

    TH1D *h(0);
    if(rereco==3){ //starting from 24C
    	//h = tdrHist("h","Response",0.92,1.06,"Run",379300,387100); // after removing 24B (start from C-rereco), 24B start: 378900 approx.
    }
    else{ //starting from 25B, B start: 391531, but no data?
	//h = tdrHist("h","Response",0.92,1.06,"Run",391500,393000);//w53
	//h = tdrHist("h","Response",0.92,1.06,"Run",391500,393200);//w54
	//h = tdrHist("h","Response",0.92,1.06,"Run",391500,396000);//w58
	h = tdrHist("h","Response",0.92,1.06,"Run",391500,397300);//w62
    }
    //lumi_136TeV = Form("Photon+jet, Run 3, %s",cv);
    //lumi_136TeV = Form("Photon+jet, Run 3 2024, %s",cv); //for 2024-only version
    //extraText = "Private"; //would print it below CMS logo, but i want it to the right (save space)

    if(rereco==1){
			lumi_136TeV = Form("Photon+jet, 2024 ECALRATIO, %s",cv);
    }
    else if(rereco==2){
			lumi_136TeV = Form("Photon+jet, 2024 ECALR-HCALDI, %s",cv);
    }
    else if(rereco==3){
			//lumi_136TeV = Form("Photon+jet, 2024CDE-rereco 2024FGHI-prompt, %s",cv);
			lumi_136TeV = Form("#gamma+jet 2024: CDE-rereco FGHI-prompt, %s",cv);
    }
    else{
    	lumi_136TeV = Form("Photon+jet, 2025CDEF, %s",cv);
	//lumi_136TeV = Form("#gamma+jet, 2024CDEFGHI-prompt, %s",cv);
    }

    ////extraText = "Private"; //would print it below CMS logo, but i want it to the right (save space)
		//h->GetXaxis()->SetLabelSize(0.26);//
		h->GetXaxis()->SetLabelSize(0.06);
		h->GetYaxis()->SetLabelSize(0.036); //doesn't do it??


    //in case we look at high etas, write it into plot
    /*
    if(eta3==true){
        lumi_136TeV = Form("Photon+jet, Run3 2024, 3<#eta<4, %s",cv);
    }
    */

    //set the axis labels custom?
    //TAxis *xaxis = h->GetXaxis();
    //xaxis->SetAxisColor(kGreen);
    //gPad->Update();
 

    TCanvas *c1 = tdrCanvas("c1",h,8,11);
    TLine *l = new TLine();
    TLatex *t = new TLatex();
    //t->SetTextSize(0.045);// t->SetNDC();
    t->SetTextSize(0.036);


    //custom place for "Private" extraText (usually handled via tdrCanvas)
    TLatex *extratext = new TLatex();
    extratext->SetTextSize(0.039);
    extratext->SetTextFont(52);
    //extratext->DrawLatex(380200, 1.0465, "Preliminary");
    //extratext->DrawLatex(391730, 1.0460, "Preliminary"); //25B start 391531
    //extratext->DrawLatex(391590, 1.0400, "Preliminary"); //25B start 391531
    extratext->DrawLatex(391700, 1.0400, "Preliminary"); //25B start 391531, w58




    // Start drawing
    c1->cd();
    gPad->SetLogx();
    double x1(h->GetXaxis()->GetXmin()), x2(h->GetXaxis()->GetXmax());
    double y1(h->GetMinimum()), y2(h->GetMaximum());
    l->DrawLine(x1,1,x2,1);

    l->SetLineStyle(kDotted); 

    // 2024 (start and end of runs)
    // see: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis
    // this commented-out line took the numbers from json i think:
    //double run24b_start(378981), run24b_end(379356); //379355 end of 24B? 
 
    double run25b_start(391531), run25b_end(392158); //end of 24B? - but it based on start of C... Bstart:391531
    l->DrawLine(run25b_start,y1+0.045,run25b_start,y2-0.035);
    t->DrawLatex(run25b_start,0.960,"25B");
    l->DrawLine(run25b_end,y1+0.045,run25b_end,y2-0.035); //to separate visually 25B and 25C
 
    double run25cv1_start(392159), run25cv1_end(393108); //run25c_end(392524); //run25c_end(392441); //ends?
    l->DrawLine(run25cv1_start,y1+0.045,run25cv1_start,y2-0.035); //used to draw this -15 before the first data point of this era...
    l->DrawLine(run25cv1_end,y1+0.045,run25cv1_end,y2-0.035);
    t->DrawLatex(run25cv1_start,0.960,"25Cv1");

    double run25cv2_start(393111), run25cv2_end(393609); //run25c_end(392524); //run25c_end(392441); //ends?
    l->DrawLine(run25cv2_start,y1+0.045,run25cv2_start,y2-0.035); //used to draw this -15 before the first data point of this era...
    l->DrawLine(run25cv2_end,y1+0.045,run25cv2_end,y2-0.035);
    t->DrawLatex(run25cv2_start,0.960,"25Cv2");

    double run25d_start(394286), run25d_end(395967); //updated end of D (30.09.)
    l->DrawLine(run25d_start,y1+0.045,run25d_start,y2-0.035);
    t->DrawLatex(run25d_start,0.960,"25D");
    l->DrawLine(run25d_end,y1+0.045,run25d_end,y2-0.035); //to separate visually 25C and 25D
 
    double run25e_start(395968), run25e_end(396597); //
    l->DrawLine(run25e_start,y1+0.045,run25e_start,y2-0.035);
    t->DrawLatex(run25e_start,0.960,"25E");
    l->DrawLine(run25e_end,y1+0.045,run25e_end,y2-0.035);

    double run25f_start(396598), run25f_end(397097); // from my json used on 21st of Sep
    l->DrawLine(run25f_start,y1+0.045,run25f_start,y2-0.035);
    t->DrawLatex(run25f_start,0.960,"25F");
    l->DrawLine(run25f_end,y1+0.045,run25f_end,y2-0.035); //PRELIMINARY end of F
 
 
/*
    //removing B if rereco==3 (06.04.2025)
    double run24b_start(378971), run24b_end(379411); //379355 end of 24B? 
    if(rereco!=3){
      //l->DrawLine(run24b_start-15,y1+0.040,run24b_start-15,y2-0.045);
      l->DrawLine(run24b_start,y1+0.045,run24b_start,y2-0.035);
      //l->DrawLine(run24b_end,y1+0.035,run24b_end,y2-0.050);
      t->DrawLatex(run24b_start,0.960,"24B");
      l->DrawLine(run24b_end,y1+0.045,run24b_end,y2-0.035); //to separate visually 24B and 24C
    }

    double run24c_start(379412), run24c_end(380252); //actually starting at 379415 (first run#), ends?
    l->DrawLine(run24c_start,y1+0.045,run24c_start,y2-0.035); //used to draw this -15 before the first data point of this era...
    l->DrawLine(run24c_end,y1+0.045,run24c_end,y2-0.035);
    t->DrawLatex(run24c_start,0.960,"24C");

    double run24d_start(380253), run24d_end(380947); //(first run#), ends?
    l->DrawLine(run24d_start,y1+0.045,run24d_start,y2-0.035);
    l->DrawLine(run24d_end,y1+0.045,run24d_end,y2-0.035);
    t->DrawLatex(run24d_start,0.960,"24D");

    double run24ev1_start(380948), run24ev1_end(381383); //(first run#), ends?
    l->DrawLine(run24ev1_start,y1+0.045,run24ev1_start,y2-0.035);
    l->DrawLine(run24ev1_end,y1+0.035,run24ev1_end,y2-0.050);
    t->DrawLatex(run24ev1_start,0.960,"24Ev1");

    double run24ev2_start(381384), run24ev2_end(381943); //(first run#), ends?
    l->DrawLine(run24ev2_start,y1+0.025,run24ev2_start,y2-0.035);
    l->DrawLine(run24ev2_end,y1+0.035,run24ev2_end,y2-0.060);
    t->DrawLatex(run24ev2_start,0.948,"24Ev2"); //a bit lower than for the other eras, to make it look prettier
    //l->DrawLine(381594,y1+0.045,381594,y2-0.035); //"preliminary" end of Ev2, until we have real one

    double run24f_start(381944), run24f_end(383779); //(first run#), ends?
    l->DrawLine(run24f_start,y1+0.045,run24f_start,y2-0.035);
    l->DrawLine(run24f_end,y1+0.035,run24f_end,y2-0.050);
    t->DrawLatex(run24f_start,0.960,"24F");

    double run24g_start(383780), run24g_end(385813); //(first run#), ends?
    l->DrawLine(run24g_start,y1+0.045,run24g_start,y2-0.035);
    l->DrawLine(run24g_end,y1+0.035,run24g_end,y2-0.050);
    t->DrawLatex(run24g_start,0.960,"24G");
		//l->DrawLine(385619,y1+0.045,385619,y2-0.035); //"preliminary" end of G, until we have real one

    double run24h_start(385814), run24h_end(386402); //(first run#), ends?
    l->DrawLine(run24h_start,y1+0.045,run24h_start,y2-0.035);
    l->DrawLine(run24h_end,y1+0.035,run24h_end,y2-0.050);
    t->DrawLatex(run24h_start,0.960,"24H");

    double run24i_start(386403); //, run24i_end(); //(first run#), ends?
    l->DrawLine(run24i_start,y1+0.045,run24i_start,y2-0.035);
    //l->DrawLine(run24d_end,y1+0.035,run24d_end,y2-0.050);
    t->DrawLatex(run24i_start,0.960,"24I");
    //l->DrawLine(384491,y1+0.045,384491,y2-0.035); //"preliminary" end of I, until we have real one
		l->DrawLine(386951,y1+0.045,386951,y2-0.035); //"preliminary" end of I, until we have real one
*/



    //first run number for each era
    TLatex *runstart = new TLatex();
    runstart->SetTextSize(0.026);
    runstart->SetTextColor(kBlack);
    runstart->SetTextFont(52);
    runstart->DrawLatex(run25b_start, 0.955, Form("%.f",run25b_start));
    runstart->DrawLatex(run25cv1_start, 0.955, Form("%.f",run25cv1_start));
    runstart->DrawLatex(run25cv2_start, 0.955, Form("%.f",run25cv2_start));
    runstart->DrawLatex(run25d_start, 0.955, Form("%.f",run25d_start));
    runstart->DrawLatex(run25e_start, 0.955, Form("%.f",run25e_start));
    runstart->DrawLatex(run25f_start, 0.955, Form("%.f",run25f_start));


    //runstart->DrawLatex(395032, 0.960, "395032"); //preliminary end of 25D

/*
    if(rereco!=3){ runstart->DrawLatex(run24b_start, 0.955, Form("%.0f",run24b_start)); } //removed B on 06.04.2025 if rereco==3
    runstart->DrawLatex(run24c_start, 0.955, Form("%.f",run24c_start));
    runstart->DrawLatex(run24d_start, 0.955, Form("%.f",run24d_start));
    runstart->DrawLatex(run24ev1_start, 0.955, Form("%.f",run24ev1_start));
    runstart->DrawLatex(run24ev2_start, 0.943, Form("%.f",run24ev2_start));
    runstart->DrawLatex(run24f_start, 0.955, Form("%.f",run24f_start));
    runstart->DrawLatex(run24g_start, 0.955, Form("%.f",run24g_start));
    runstart->DrawLatex(run24h_start, 0.955, Form("%.f",run24h_start));
    runstart->DrawLatex(run24i_start, 0.955, Form("%.f",run24i_start));
*/



    //runstart->DrawLatex(380649, 0.960, "380649"); //preliminary end of 24D
    //runstart->DrawLatex(381594, 1.020, "381594"); //preliminary end of 24Ev2
		//runstart->DrawLatex(384491, 1.020, "384332"); //preliminary end of 24G
		//runstart->DrawLatex(385619, 1.020, "385619"); //preliminary end of 24G



    // ------------- STRANGE JUMPS ---------------------------//
    //strange things, where we try to find (and draw) the run-number to check back with hcal etc
    //i identified these visually and then put a line there and wrote the run number
    TLine *strangeline = new TLine();
    strangeline->SetLineStyle(kSolid);
    strangeline->SetLineColor(kRed-2);
    //double strangerun(379965); //379980
    //double strangerun(383200);
    double strangerun(384635);
    /////strangeline->DrawLine(strangerun,y1+0.030,strangerun,y2-0.030);     //strange jump
    //strangeline->DrawLine(384910,y1+0.040,384910,y2-0.030);     	//drop in response in G
    TLatex *strangeinfo = new TLatex();
    strangeinfo->SetTextSize(0.020); //0.026
    strangeinfo->SetTextColor(kRed-2);
    strangeinfo->SetTextFont(52);
    /////strangeinfo->DrawLatex(strangerun, 0.95, Form("%.f: change?",strangerun));

    strangeinfo->SetTextSize(0.020);
    //strangeinfo->DrawLatex(384730, 0.956, "384910"); //drop in response (G)
    //strangeinfo->DrawLatex(384730, 0.953, "change?"); //drop in response (G)


    //locating jump in F
    //strangeline->DrawLine(383219,y1+0.030,383219,y2-0.030); 	//HCALRespCorrs (run 383219) 
    //strangeline->DrawLine(383277,y1+0.030,383277,y2-0.030); //ECAL pedestals corrections (run 383186, w29) 


    //strangeinfo->DrawLatex(strangerun, 0.945, Form("%.f: What happens?",strangerun));


// UPDATE THE FOLLOWING FOR 2025
    //add important run numbers
    //HCAL Updates:     https://twiki.cern.ch/twiki/bin/view/CMS/HcalDPGRun3Updates
    //HCAL scans:       https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalPhaseScan
    //-------------------- HCAL -----------------------//
    // so far scans were reported for: Run 379349: 46.48/pb, and Run 379350: 38.75/pb 
    //could do this in a loop in case there will be more of these
    TLine *runinfoline = new TLine();
    runinfoline->SetLineStyle(kSolid);
    runinfoline->SetLineColor(kOrange+2);

    //Pedestal/LUT updates occurred twice in the last couple weeks: (got this info from Mikko, HCAL meeting)
    runinfoline->DrawLine(385221,y1+0.022,385221,y2-0.035);  	//2 September [elog], run 385221
    runinfoline->DrawLine(385550,y1+0.022,385550,y2-0.035);	//11 September [elog], run 385550


    TLatex *runinfo = new TLatex();
    runinfo->SetTextSize(0.020); //0.026
    runinfo->SetTextColor(kOrange+2);
    runinfo->SetTextFont(52);

    if(rereco!=3){
	runinfoline->DrawLine(379349,y1+0.018,379349,y2-0.035);     //hcal scan
	runinfoline->DrawLine(379350,y1+0.016,379350,y2-0.035);     //hcal scan

	runinfo->DrawLatex(379349, 0.932, "379349: HCAL scan");
	runinfo->DrawLatex(379350, 0.928, "379350: HCAL scan");
    }

    runinfo->SetTextSize(0.020);
    runinfo->DrawLatex(385221, 0.939, "385221 & 385550:");
    runinfo->DrawLatex(385221, 0.9353, "Pedestal/LUT updates");


    //HCALRespCorrs
    runinfoline->DrawLine(383219,y1+0.030,383219,y2-0.035); 	//HCALRespCorrs (run 383219) 
    runinfo->DrawLatex(383219, 0.946, "383219: HCALRespCorrs");


    //changes in middle of era G (various between 384635 and 384910, also in ECAL)
    runinfoline->DrawLine(384858,y1+0.035,384858,y2-0.016); 	//(three HCAL updates in same runnumber)
    //runinfo->SetTextSize(0.016); //temporarily go to smaller textsize, to fit a lot in small gap...
    //runinfo->DrawLatex(384858, 0.96, "384858:");
    runinfo->DrawLatex(384858, 1.0435, "384858:");
    runinfo->DrawLatex(384858, 1.040, "HcalL1TriggerObjects");	//(New Hcal L1 Trigger Objects AlCaDB cmsTalk) 
    runinfo->DrawLatex(384858, 1.0365, "HcalPedestalWidths");	//(New Hcal pedestals width conditions AlCaDB cmsTalk) 
    runinfo->DrawLatex(384858, 1.033, "HcalPedestals");		//(New Hcal pedestals conditions AlCaDB cmsTalk)




    //------------ ECAL --------------------//
    //(ECAL?) Updates: https://twiki.cern.ch/twiki/bin/viewauth/CMS/AlCaTSGConditionsUpdate
    TLine *ecalline = new TLine();
    ecalline->SetLineStyle(kSolid);
    ecalline->SetLineColor(kCyan+2);
    ecalline->DrawLine(379956,y1+0.022,379956,y2-0.035);     //

    //start end end of this ECalTimeCalibConstants update
    //Prompt | EcalIntercalibConstants | EcalIntercalibConstants_V1_prompt | 379956 | Update of EcalIntercalibConstants conditions from runs 378981-379616 ALCaDB cmsTalk.
    //ecalline->SetLineStyle(kDotted);
    //ecalline->DrawLine(378981,y1+0.030,378981,y2-0.030); 
    //ecalline->DrawLine(379616,y1+0.030,379616,y2-0.030); 

    TLatex *ecalinfo = new TLatex();
    ecalinfo->SetTextSize(0.020); //0.026
    ecalinfo->SetTextColor(kCyan+2);
    ecalinfo->SetTextFont(52);
    ecalinfo->DrawLatex(379956, 0.938, "379956: EcalIntercalibConstants");
 
    ecalline->DrawLine(383227,y1+0.0265,383227,y2-0.035); //(ECAL Intercalib constants: new payload with ICs up to run 382913 ALCaDB cmsTalk)
    ecalinfo->DrawLatex(383227, 0.9426, "383227: EcalIntercalibConstants");
	
    //at 383277: ECAL pedestals corrections (run 383186, w29), did not result in a major jump...


    //changes in middle of era G (various between 384635 and 384910, also in HCAL)
    /*
    ecalinfo->SetTextSize(0.016); //temporarily go to smaller textsize, to fit a lot in small gap...
    ecalline->DrawLine(384719, y1+0.0265, 384719, y2-0.030); //(Ecal Pedestals update (run 384719, w34) AlCaDB cmsTalk)
    ecalinfo->DrawLatex(384719, 0.9653, "384719: EcalPedestals");
    ecalline->DrawLine(384756, y1+0.0265, 384756, y2-0.030); //(EcalTimeCalibConstants update after front end delay change AlCaDB cmsTalk)
    ecalinfo->DrawLatex(384756, 0.963, "384756: EcalTimeCalibConstants");
    */ 
    //move this info to top of plot for space reasons...
    ecalline->DrawLine(384719, y1+0.035, 384719, y2-0.0065); //(Ecal Pedestals update (run 384719, w34) AlCaDB cmsTalk)
    ecalinfo->DrawLatex(384719, 1.053, "384719: EcalPedestals");
    ecalline->DrawLine(384756, y1+0.035, 384756, y2-0.010); //(EcalTimeCalibConstants update after front end delay change AlCaDB cmsTalk)
    ecalinfo->DrawLatex(384756, 1.050, "384756: EcalTimeCalibConstants");
 
	




    //------------ HCAL as arrow --------------------//
    //HCAL stuff taken care of before 2024B (mark with arrow to the left, remove later)
    //TArrow *hcalArrow = new TArrow(0.5, 0.3, 0.2, 0.1, 0.05, ">"); //1-4: position, 5:arrowsize, 6:option
    //hcalArrow->SetNDC(kTRUE); //set already here to put arrow more easily
    //runinfo->DrawLatex(0.4, 0.6, "hcal arrow");

    if(rereco!=3){
	TArrow *hcalArrow = new TArrow(379360, 1.033, 379000, 1.033, 0.02, ">"); //1-4: position, 5:arrowsize, 6:option
	hcalArrow->SetLineColor(kOrange+2);
	hcalArrow->Draw();

	runinfo->DrawLatex(379060, 1.036, "377804: phaseTuning");
    }


    //time (re-)alignments & other stuff, only put here when happening after start of 24B (378971 April 5)
    //379347 April 13 phaseTuning_HB_2024_v2, improves iphi-symmetry of timing in HB (all depths)  

    //HcalPedestals: https://twiki.cern.ch/twiki/bin/view/CMS/HcalPedestalsTagsRun3#Descriptions_of_the_HLT_Express
    //379349 (24 Apr 2024)
    ///runinfo->DrawLatex(379349, 0.930, "HcalPedestals"); //but this is same as scan?


    //info on the dotted lines
    TLatex *infotext = new TLatex();
    infotext->SetNDC(kTRUE);
    //infotext->SetTextSize(0.030); //0.036
		infotext->SetTextSize(0.020); //0.036
    infotext->SetTextColor(12);
    infotext->SetTextFont(52);
    infotext->DrawLatex(0.61,0.153, Form("Error cleaning applied (%.2f as maxerr).",maxerr));
    //infotext->DrawLatex(0.53,0.13, Form("The dotted line is: start-15 and end+15."));
    infotext->DrawLatex(0.61,0.13, Form("The dotted lines are: era start and end."));


    //hyperlinks for monitoring changes
    //(ECAL?) Updates: https://twiki.cern.ch/twiki/bin/viewauth/CMS/AlCaTSGConditionsUpdate
    //HCAL Updates:     https://twiki.cern.ch/twiki/bin/view/CMS/HcalDPGRun3Updates
    //HCAL scans:       https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalPhaseScan
    TLatex *linktext = new TLatex();
    linktext->SetNDC(kTRUE);
    linktext->SetTextSize(0.020); //0.036
    linktext->SetTextColor(kBlue);
    linktext->SetTextFont(12);
    //InfoLink: <https://twiki.cern.ch/>
    linktext->DrawLatex(0.13,0.06, "https://twiki.cern.ch/twiki/bin/viewauth/CMS/AlCaTSGConditionsUpdate");
    linktext->DrawLatex(0.13,0.04, "https://twiki.cern.ch/twiki/bin/view/CMS/HcalDPGRun3Updates");
    linktext->DrawLatex(0.13,0.02, "https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalPhaseScan");



    //information on JSON used
    linktext->SetTextColor(kRed+2);
    linktext->SetTextFont(82);
    //linktext->DrawLatex(0.13,0.09, "Created on 21.05.2024 using daily JSON: Collisions24_13p6TeV_378981_380883_DCSOnly_TkPx.json");
    //linktext->DrawLatex(0.13,0.09, "Created on 31.05.2024 using golden JSON: Cert_Collisions2024_378981_380649_Golden.json");
		//linktext->DrawLatex(0.13,0.09, "Created on 10.06.2024 using hybrid JSON = newest golden (5.6.) + additional daily (9.6.)");
		//linktext->DrawLatex(0.13,0.09, "Created on 09.10.2024 using hybrid JSON = newest golden (12.9.) + additional daily (12.9.)");
		//linktext->DrawLatex(0.13,0.09, "Created on 09.10.2024 using golden JSON (29.9.) and entire 24G data.");
		//linktext->DrawLatex(0.13,0.09, "Created on 28.10.2024 with w39 for BCDEFG and w40 for H and I.");
		//linktext->DrawLatex(0.13,0.09, "Created on 12.02.2025 with w44 for 2024BCDEFGHI Prompt data.");
		if(rereco==3){
			linktext->DrawLatex(0.13,0.09, "Created on 06.04.2025 with w48 for 2024CDE-rereco and 2024FGHI-prompt data (nib-based V8M).");
		}
		else{
			linktext->DrawLatex(0.13,0.09, "Created on 21.09.2025 with w62 for 2025CDEF data.");
		}













    //pr50m->GetXaxis()->SetAxisColor(kRed);


    //tdrDraw(pr30b,"Pz",kOpenSquare,kBlue,kSolid); pr30b->SetMarkerSize(0.7);
    //tdrDraw(pr30m,"Pz",kFullCircle,kBlue,kSolid); pr30m->SetMarkerSize(0.6);
    //tdrDraw(pr110b,"Pz",kOpenCircle,kRed,kSolid); pr110b->SetMarkerSize(0.7);
/*
    tdrDraw(pr50b,"Pz",kOpenSquare,kBlue,kSolid); pr50b->SetMarkerSize(0.5);            //balance --> instead of MPF30, showing DB50
    tdrDraw(pr110m,"Pz",kFullCircle,kRed,kSolid); pr110m->SetMarkerSize(0.5);           //print this one first, as 50EB captures also everything that goes through 110EB
    tdrDraw(pr50m,"Pz",kFullCircle,kGreen+2,kSolid); pr50m->SetMarkerSize(0.5);         //marker size 0.06 or 0.05?
*/

    //set the axis labels custom?
    TAxis *xaxis = h->GetXaxis();
    //need to create x-axis, as the h histogram is fixbin (i.e. only has xmin and xmax)
    cout << "This is x1: " << x1 << endl;
    cout << "This is x2: " << x2 << endl;
    int nbins = int(x2)-int(x1);
    cout << "This is nbins: " << nbins << endl;
    cout << "This is h->GetNbinsX(): " << h->GetNbinsX() << endl;
    cout << "Labels? " << xaxis->GetLabels() << endl;
    int xarr[nbins];
    int sizearr = sizeof(xarr);
    cout << "This is sizeof(xarr): " << sizearr << endl;
    

    //for(int i=0; i<sizeof(xarr); i++){
    for(int i=1; i<=nbins; i++){
        //xarr[i]=378900+i; //starting from mininmu x from h
        //xarr[i]=x1+i-1;
        //cout << Form("xarr[%d]=",i) << xarr[i] << ", ";
        string currlabel = Form("%.f",x1+i-1); //current bin's label
        /*
        cout << "label: " << currlabel << endl;
        */
        //xaxis->SetBinLabel(i,Form("%d",x1+i-1));
        xaxis->SetBinLabel(i,Form("%s",currlabel.c_str()));
    }

    tdrDraw(pr50b,"Pz",kOpenSquare,kBlue,kSolid); pr50b->SetMarkerSize(0.5);            //balance --> instead of MPF30, showing DB50
    tdrDraw(pr40b,"Pz",kOpenSquare,kCyan,kSolid); pr40b->SetMarkerSize(0.5);            //balance for Photon40EB_TightID_TightIso
    tdrDraw(pr110m,"Pz",kFullCircle,kRed,kSolid); pr110m->SetMarkerSize(0.5);           //print this one first, as 50EB captures also everything that goes through 110EB
    tdrDraw(pr50m,"Pz",kFullCircle,kGreen+2,kSolid); pr50m->SetMarkerSize(0.5);         //marker size 0.06 or 0.05?
    //tdrDraw(pr40m,"Pz",kFullCircle,kOrange+2,kSolid); pr40m->SetMarkerSize(0.5);        //new in 25
    tdrDraw(pr40m,"Pz",kFullCircle,kMagenta,kSolid); pr40m->SetMarkerSize(0.5);        //new in 25




    //manually set selected labels
    //xaxis->SetBinLabel(1,Form("%d",run24b_start));


    // NEW IDEA: WHAT ABOUT USING THE NUMERIC X-VALUES FROM THE TPROFILES AND SETTING THEM AS ALPHANUMERIC LABELS FOR THE BACKGROUND HISTOGRAM h (in bin middle)? 
    //cout << "this is pr50m's xaxis: " << pr50m->GetXaxis()->GetXbins();
    /*
    //for(int i=0; i<pr50m->GetNbinsX(); i++){
    //for(int i=378981; i<379000; i++){
    for(int i=0; i<10; i++){
        cout << pr50m->GetXaxis()->GetXbins()->At(i) << ","; //does not work here, fixbins
    }
    */
    //TAxis *xaxis = h->GetXaxis();
    //xaxis->SetAxisColor(kGreen); //--> works now
    //xaxis->SetBinLabel();
    //xaxis->LabelsOption("d");
    h->LabelsOption("d","X");
    c1->Update();
    gPad->RedrawAxis();
    gPad->Update();
    //xaxis->Draw();
    //xaxis->LabelsOption("d");
    //h->LabelsOption("d","X");
    //xaxis->Draw();
    //pr50m->SetXLabels();
    //h->LabelsOption("d","X");
    //pr50m->Draw(); //draw again to colour axis?

    //tdrDraw(pr230m,"Pz",kFullCircle,kGray,kSolid); pr230m->SetMarkerSize(0.6);


    // Add legend
    //pr50legentry = Form()
    //c1->cd(1);
    //TLegend *leg = tdrLeg(0.54,0.68,0.64,0.88);
    //TLegend *leg = tdrLeg(0.64,0.78,0.84,0.88);
    //TLegend *leg = tdrLeg(0.64,0.76,0.84,0.88);
    //
    //TLegend *leg = tdrLeg(0.62,0.78,0.82,0.90);
    //TLegend *leg = tdrLeg(0.42,0.80,0.62,0.88); //put it more to middle to make space for jumptexts
    TLegend *leg = tdrLeg(0.72,0.75,0.92,0.89); //for 2025
    leg->SetTextSize(0.030); //smaller size for squeezing more into the plot



    leg->SetTextFont(42);
    leg->AddEntry(pr110m,"MPF 110EB","PLE");
    ///leg->AddEntry(pr110m, Form("MPF 110EB x%.2f", kmpf110),"PLE");
    //leg->AddEntry(pr110b,"BAL 110EB","PLE");
    //leg->AddEntry(pr50m, Form("MPF 50EB x%.2f", kpt50), "PLE");
    //
    leg->AddEntry(pr50m, "MPF 50EB", "PLE");
    ///leg->AddEntry(pr50m, Form("MPF 50EB x%.2f", kmpf50), "PLE");
    leg->AddEntry(pr40m, "MPF 40EB", "PLE"); //new
    //leg->AddEntry(pr40m, Form("MPF 40EB x%.2f", kmpf40), "PLE"); //new
    leg->AddEntry(pr50b, Form("BAL 50EB x%.2f", kbal50), "PLE");
    leg->AddEntry(pr40b, Form("BAL 40EB x%.2f", kbal40), "PLE"); //new
    //leg->AddEntry(pr30m,"MPF 30EB","PLE");


    //leg->AddEntry(pr30b,"BAL 30EB","PLE");

    //h->Draw();
    c1->Update(); //updating helping?
    cout << "Labels? " << xaxis->GetLabels() << endl;


    if(rereco==0){
    	c1->SaveAs(Form("pdf/drawResponseVsRun_2025only_%s.pdf",cv));
    	//c1->SaveAs(Form("pdf/drawResponseVsRun_2025BC_noL2L3Res_%s.pdf",cv));
    }
    else if(rereco==1){
    	c1->SaveAs(Form("pdf/drawResponseVsRun_2024-ECALRATIO_%s.pdf",cv));
    }
    else if(rereco==2){
    	c1->SaveAs(Form("pdf/drawResponseVsRun_2024C-ECALR-HCALDI_%s.pdf",cv));
    }
    else if(rereco==3){
    	c1->SaveAs(Form("pdf/drawResponseVsRun_2024CDE-rereco_2024FGHI-prompt_%s.pdf",cv));
    }




    // Extra composition plots also modified.
		bool zoom = 0;
    drawPFcompVsRun_2025only(version, rereco, zoom);
    drawPFcompEFBvsRun_2025only(version, rereco, zoom);




} // drawResponseVsRun_2025only



//start of PF composition versus run, modified for 2024 only; for now focus on 50EB trigger
void drawPFcompVsRun_2025only(string version, int rereco, bool zoom) {

    const char *cv = version.c_str();

    setTDRStyle();
    TDirectory *curdir = gDirectory;

    // Open input files
    //TFile *f = new TFile(Form("rootfiles/GamHistosFill_data_2025only_%s.root",cv),"READ");
    TFile *f(0);

    if(rereco==0){
	cout << "Looking at 2025 prompt data." << endl;
	f = new TFile(Form("rootfiles/GamHistosFill_data_2025CDEF_%s.root",cv), "READ"); 
	//f = new TFile(Form("rootfiles/GamHistosFill_data_2025BC_noL2L3Res_%s_28May2025.root",cv), "READ"); //prompt, nib-based
    	//f = new TFile(Form("rootfiles/GamHistosFill_data_2025only_%s.root",cv), "READ"); //for now: hadd on 2024B and 2024C (need to redo this file with new daily json)
	//f = new TFile(Form("rootfiles/GamHistosFill_data_2025only_BCDEFG-w39_HI-w40.root"), "READ");//intermediate version where i mixed w39 and w40

    }
    else if(rereco==1){
	cout << "Looking at 2024 rereco (ECALRATIO) data." << endl;
    	f = new TFile(Form("rootfiles/GamHistosFill_data_2024-ECALRATIO_%s.root",cv), "READ"); //for now: hadd on 2024B and 2024C (need to redo this file with new daily json)
    }
    else if(rereco==2){
	cout << "Looking at 2024C 2nd rereco (ECALR-HCALDI) data." << endl;
    	f = new TFile(Form("rootfiles/GamHistosFill_data_2024C-ECALR-HCALDI_%s.root",cv), "READ"); //
    }
    else if(rereco==3){
	cout << "Looking at 2024CDE-rereco and 2024FGHI-prompt data." << endl;
	f = new TFile(Form("rootfiles/GamHistosFill_data_2024CDE-rereco_2024FGHI-prompt_%s.root",cv), "READ"); //containing 2024CDE-rereco and 2024FGHI-prompt
    }
    else {
	cout << "Check rereco boolian. Something's wrong." << endl;
	exit(0);
    }
    assert(f && !f->IsZombie());
    f->cd("runs");
    TDirectory *d = gDirectory;



    /*
    TFile *fr = new TFile("rootfiles/GamHistosFill_data_Run3_v29.root","READ");
    assert(fr && !fr->IsZombie());
    fr->cd("control");
    TDirectory *dr = gDirectory;
    */

    // Load input profiles; should adjust error cleaning and also add this to legend as info
    /*
    TProfile *pr30m = (TProfile*)dr->Get("pr30m"); clean(pr30m,0.006);
    TProfile *pr30chf = (TProfile*)d->Get("pr30chf"); clean(pr30chf,0.006);
    TProfile *pr30nhf = (TProfile*)d->Get("pr30nhf"); clean(pr30nhf,0.006);
    TProfile *pr30nef = (TProfile*)d->Get("pr30nef"); clean(pr30nef,0.006);
    TProfile *pr110m = (TProfile*)dr->Get("pr110m"); clean(pr110m,0.003);
    TProfile *pr110chf = (TProfile*)d->Get("pr110chf"); clean(pr110chf,0.003);
    TProfile *pr110nhf = (TProfile*)d->Get("pr110nhf"); clean(pr110nhf,0.003);
    TProfile *pr110nef = (TProfile*)d->Get("pr110nef"); clean(pr110nef,0.003);
    */
    TProfile *pr50m = (TProfile*)d->Get("pr50m"); clean(pr50m,0.006);
    TProfile *pr50chf = (TProfile*)d->Get("pr50chf"); clean(pr50chf,0.006);
    TProfile *pr50nhf = (TProfile*)d->Get("pr50nhf"); clean(pr50nhf,0.006);
    TProfile *pr50nef = (TProfile*)d->Get("pr50nef"); clean(pr50nef,0.006);

    // CAN NOW USE PHOTON40EB_TIGHTID_TIGHTISO
    TProfile *pr40m = (TProfile*)d->Get("pr40m"); clean(pr40m,0.006);
    TProfile *pr40chf = (TProfile*)d->Get("pr40chf"); clean(pr40chf,0.006);
    TProfile *pr40nhf = (TProfile*)d->Get("pr40nhf"); clean(pr40nhf,0.006);
    TProfile *pr40nef = (TProfile*)d->Get("pr40nef"); clean(pr40nef,0.006);
 
 

    curdir->cd();
    /*
    TH1D *hr30m = pr30m->ProjectionX("hr30m");
    TH1D *hr30chf = pr30chf->ProjectionX("hr30chf");
    TH1D *hr30nhf = pr30nhf->ProjectionX("hr30nhf");
    TH1D *hr30nef = pr30nef->ProjectionX("hr30nef");
    TH1D *hr110m = pr110m->ProjectionX("hr110m");
    TH1D *hr110chf = pr110chf->ProjectionX("hr110chf");
    TH1D *hr110nhf = pr110nhf->ProjectionX("hr110nhf");
    TH1D *hr110nef = pr110nef->ProjectionX("hr110nef");
    */
    TH1D *hr50m = pr50m->ProjectionX("hr50m");          //what happens here?!
    TH1D *hr50chf = pr50chf->ProjectionX("hr50chf");
    TH1D *hr50nhf = pr50nhf->ProjectionX("hr50nhf");
    TH1D *hr50nef = pr50nef->ProjectionX("hr50nef");

    //since 2025
    TH1D *hr40m = pr40m->ProjectionX("hr40m");          //what happens here?!
    TH1D *hr40chf = pr40chf->ProjectionX("hr40chf");
    TH1D *hr40nhf = pr40nhf->ProjectionX("hr40nhf");
    TH1D *hr40nef = pr40nef->ProjectionX("hr40nef");
 
 

    // Offset composition
    /*
    addOffset(hr30m,-1.020 +0.01);
    addOffset(hr30chf,-0.68);
    addOffset(hr30nef,-0.20);
    addOffset(hr30nhf,-0.11);
    addOffset(hr110m,-1.005 +0.01);
    addOffset(hr110chf,-0.61);
    addOffset(hr110nef,-0.26);
    addOffset(hr110nhf,-0.10);
    */
    ///*

/*
	  	double hr50moff = -1.00; //in V2M had this as -1.01
    	double hr50chfoff = -0.665;
    	double hr50nefoff = -0.23;
    	double hr50nhfoff = -0.075;
*/
	

		double hr50moff, hr50chfoff, hr50nefoff, hr50nhfoff;
		if(zoom){
			hr50moff = -1.00; //in V2M had this as -1.01
    	hr50chfoff = -0.665;
    	hr50nefoff = -0.23;
    	hr50nhfoff = -0.075;
		}
		else{
    	hr50moff = -1.00;
    	hr50chfoff = -0.68;
    	hr50nefoff = -0.20;
    	hr50nhfoff = -0.11;
		}

    addOffset(hr50m,hr50moff); //need to adjust this and the following still
    addOffset(hr50chf,hr50chfoff);
    addOffset(hr50nef,hr50nefoff);
    addOffset(hr50nhf,hr50nhfoff);

    //for now use same offsets in 25
    addOffset(hr40m,hr50moff); //need to adjust this and the following still
    addOffset(hr40chf,hr50chfoff);
    addOffset(hr40nef,hr50nefoff);
    addOffset(hr40nhf,hr50nhfoff);





    //*/
 

    // Setup canvas
    TH1D *h;
    TH1D *h40;
    //TH1D *h = tdrHist("h2","PF composition offset",-0.10,+0.10,"Run",378900,380800);
    /*
    if(zoom){
    	h = tdrHist("h2","PF composition offset",-0.025,+0.025,"Run",378900,387100);
    }
    else{
	h = tdrHist("h2","PF composition offset",-0.10,+0.10,"Run",378900,387100);
    }
    */
    if(rereco==3){ //starting from 24C
    	h = tdrHist("h","PF composition offset",-0.10,+0.10,"Run",379300,387100); // after removing 24B (start from C-rereco) 
    }
    else{ //starting from 24B
	//h = tdrHist("h","PF composition offset",-0.10,+0.10,"Run",378900,387100);
	//h = tdrHist("h","PF composition offset",-0.10,+0.10,"Run",391500,393200);//w54
	//h40 = tdrHist("h40","PF composition offset",-0.10,+0.10,"Run",391500,393200);//w54 for photon40
	h = tdrHist("h","PF composition offset",-0.10,+0.10,"Run",391500,397300);//w54
	h40 = tdrHist("h40","PF composition offset",-0.10,+0.10,"Run",391500,397300);//w54 for photon40
    }

    //TH1D *h = tdrHist("h2","PF composition offset",-0.80,+0.80,"Run",378900,380600);
    //h->GetXaxis()->SetLabelSize(0.15); //smaller labels
    h->GetXaxis()->SetLabelSize(0.26);
    h->GetYaxis()->SetLabelSize(0.036);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    gPad->Update();

    h40->GetXaxis()->SetLabelSize(0.26);
    h40->GetYaxis()->SetLabelSize(0.036);
    h40->GetXaxis()->SetTitleSize(0.05);
    h40->GetYaxis()->SetTitleSize(0.05);
    gPad->Update();


		

    //lumi_136TeV = Form("Photon+jet, Run 3, %s",cv);
    if(rereco==1){
	lumi_136TeV = Form("Photon+jet, 2024 ECALRATIO, %s",cv);
    }
    else if(rereco==2){
	lumi_136TeV = Form("Photon+jet, 2024 ECALR-HCALDI, %s",cv);
    }
    else if(rereco==3){
	lumi_136TeV = Form("#gamma+jet, 2024CDE-rereco 2024FGHI-prompt, %s",cv);
    }
    else{
    	//lumi_136TeV = Form("#gamma+jet, 2024BCDEFGHI-prompt, %s",cv);
    	lumi_136TeV = Form("#gamma+jet, 2025CDEF-prompt, %s",cv);
    }
    extraText = "Private";
    TCanvas *c1 = tdrCanvas("c2",h,8,11);
    TCanvas *c2 = tdrCanvas("c3",h40,8,11);


    TLine *l = new TLine();
    TLatex *t = new TLatex();
    //t->SetTextSize(0.045);
		t->SetTextSize(0.036);


    //added this to composition plot.
    TLatex *extratext = new TLatex();
    extratext->SetTextSize(0.039);
    extratext->SetTextFont(52);
    extratext->DrawLatex(391700, 0.0715, "Preliminary");

    h->GetXaxis()->SetLabelSize(0.15); //smaller labels
    gPad->Update();
	

    // Start drawing
    c1->cd();
    gPad->SetLogx();
    double x1(h->GetXaxis()->GetXmin()), x2(h->GetXaxis()->GetXmax());
    double y1(h->GetMinimum()), y2(h->GetMaximum());
    l->DrawLine(x1,0,x2,0);

    extratext->DrawLatex(391700, y2-0.0265, "Preliminary"); //TEST



    // 2024 Era definition
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#2024_Era_definition
    double run24b_start(378971), run24b_end(379411);
    double run24c_start(379412), run24c_end(380252);
    double run24d_start(380253), run24d_end(380947);
    double run24ev1_start(380948), run24ev1_end(381383);
    double run24ev2_start(381384), run24ev2_end(381943); //(first run#), ends?
    double run24f_start(381944), run24f_end(383779);
    double run24g_start(383780), run24g_end(385813); //(first run#), ends?
    double run24h_start(385814), run24h_end(386402); //(first run#), ends?
    double run24i_start(386403);

    //2025
    double run25b_start(391531), run25b_end(392158); //end of 24B? - but it based on start of C... Bstart:391531
    double run25cv1_start(392159), run25cv1_end(393108); //run25c_end(392524); //run25c_end(392441); //ends?
    double run25cv2_start(393111), run25cv2_end(393609);
    //double run25d_start(394286), run25d_end(395032);
    double run25d_start(394286), run25d_end(395967); //updated end of D (30.09.)
    double run25e_start(395968), run25e_end(396597); //
    double run25f_start(396598), run25f_end(397097); // from my json used on 21st of Sep
 

    //l->SetLineStyle(kDashed);
    l->SetLineColor(kGray);
    if(rereco!=3){
	l->DrawLine(run24b_start,y1,run24b_start,y2); 	//removed B on 06.04.2025
	l->DrawLine(run24b_end,y1,run24b_end,y2);	//removed B on 06.04.2025
	//t->DrawLatex(run24b_start+250,-0.095,"24B"); //+250 offset only when displaying in middle of era
    }

    //2025
    l->DrawLine(run25b_start,y1,run25b_start,y2);
    l->DrawLine(run25b_end,y1,run25b_end,y2); //to separate visually 25B and 25C
    l->DrawLine(run25cv1_start,y1,run25cv1_start,y2); //used to draw this -15 before the first data point of this era...
    l->DrawLine(run25cv1_end,y1,run25cv1_end,y2);
    l->DrawLine(run25cv2_start,y1,run25cv2_start,y2); //2025C-v2 start and end
    l->DrawLine(run25cv2_end,y1,run25cv2_end,y2);
    l->DrawLine(run25d_start,y1,run25d_start,y2); //2025D start and end
    l->DrawLine(run25d_end,y1,run25d_end,y2);
    l->DrawLine(run25e_start,y1,run25e_start,y2); //2025D start and end
    l->DrawLine(run25e_end,y1,run25e_end,y2);
    l->DrawLine(run25f_start,y1,run25f_start,y2); //2025D start and end
    l->DrawLine(run25f_end,y1,run25f_end,y2);



    //2024 (maybe remove, since this script is 2025only)
    l->DrawLine(run24c_start,y1,run24c_start,y2);
    l->DrawLine(run24c_end,y1,run24c_end,y2);
    l->DrawLine(run24d_start,y1,run24d_start,y2);
    l->DrawLine(run24d_end,y1,run24d_end,y2);
    l->DrawLine(run24ev1_start,y1,run24ev1_start,y2);
    l->DrawLine(run24ev2_start,y1,run24ev2_start,y2);
    l->DrawLine(run24f_start,y1,run24f_start,y2);
    l->DrawLine(run24g_start,y1,run24g_start,y2);
    l->DrawLine(run24h_start,y1,run24h_start,y2);
    l->DrawLine(run24i_start,y1,run24i_start,y2);




 
    //text for eras
    /*
    t->DrawLatex(run24b_start+10,-0.090,"24B");
    t->DrawLatex(run24c_start+10,-0.090,"24C");
    t->DrawLatex(run24d_start+10,-0.090,"24D");
    */

    double textposy;
    if(zoom){
	    textposy=-0.020;
    }
    else{
	    textposy=-0.090;
    }
    if(rereco!=3){ t->DrawLatex(run24b_start+10,textposy,"24B"); }

    //2025
    t->DrawLatex(run25b_start+10,textposy,"25B");
    t->DrawLatex(run25cv1_start+10,textposy,"25Cv1");
    t->DrawLatex(run25cv2_start+10,textposy,"25Cv2");
    t->DrawLatex(run25d_start+10,textposy,"25D");
    t->DrawLatex(run25e_start+10,textposy,"25E");
    t->DrawLatex(run25f_start+10,textposy,"25F");



    //2024
    t->DrawLatex(run24c_start+10,textposy,"24C");
    t->DrawLatex(run24d_start+10,textposy,"24D");
    //t->DrawLatex(run24ev1_start+10,textposy,"24Ev1,v2");
    t->DrawLatex(run24ev1_start+10,textposy,"24Ev1,");
    t->DrawLatex(run24ev1_start+450,textposy-0.008,"v2");
    //t->DrawLatex(run24ev2_start+10,textposy-0.01,"24Ev2");
    t->DrawLatex(run24f_start+10,textposy,"24F");
    t->DrawLatex(run24g_start+10,textposy,"24G");
    t->DrawLatex(run24h_start+10,textposy,"24H");
    t->DrawLatex(run24i_start+10,textposy,"24I");





    /*
    l->DrawLine(run23c1b,y1,run23c1b,y2);
    l->SetLineStyle(kDashed);
    l->DrawLine(run23c1e,y1+0.015,run23c1e,y2);
    l->DrawLine(run23c2b,y1+0.015,run23c2b,y2);
    l->DrawLine(run23c2e,y1+0.015,run23c2e,y2);
    l->DrawLine(run23c3b,y1+0.015,run23c3b,y2);
    l->SetLineStyle(kSolid);
    l->DrawLine(run23c3e,y1+0.015,run23c3e,y2);
    l->DrawLine(run23c4b,y1+0.015,run23c4b,y2);
    l->DrawLine(run23c4e,y1,run23c4e,y2);
    t->DrawLatex(run23c1b+250,-0.095,"23C");
    t->DrawLatex(run23c4b+250,-0.085,"v4");
    double run23d1b(369803), run23d1e(370602);
    double run23d2b(370603), run23d2e(372415);
    l->DrawLine(run23d1b,y1,run23d1b,y2);
    l->SetLineStyle(kDashed);
    l->DrawLine(run23d1e,y1+0.015,run23d1e,y2);
    l->DrawLine(run23d2b,y1+0.015,run23d2b,y2);
    l->SetLineStyle(kSolid);
    l->DrawLine(run23d2e,y1,run23d2e,y2);
    t->DrawLatex(run23d1b+250,-0.095,"23D");
    */

    /*
    l->SetLineWidth(2);
    l->SetLineColor(kMagenta+1);
    double run22ae(354640); // HB "ieta flattening" (2022A)
    //l->DrawLine(run22ae,y1+0.025,run22ae,y2-0.025);
    double run22f(362102); // HB realignment (2022F)
    l->DrawLine(run22f,y1+0.025,run22f,y2-0.025);
    l->SetLineColor(kOrange+1);
    double run22ab(352319); // (start of era A) no HB corrections (respCorr = 1.0)
    //l->DrawLine(run22ab,y1+0.025,run22ab,y2-0.025);
    double run22fb(360329); // 360329 (start of era F) big HB corrections (~ +17%)
    l->DrawLine(run22fb,y1+0.025,run22fb,y2-0.025);

    l->SetLineColor(kRed+1);
    double runHBG(367765); // response correction update using 2022G data, deployed to prompt only
    l->DrawLine(runHBG,y1+0.025,runHBG,y2-0.025);
    double runHB(368775); // HB realignment based on TDC timing
    //also: 368775 compensating response correction for the HB time alignment
    l->DrawLine(runHB,y1+0.025,runHB,y2-0.025);
    */

    //Drawing everything
    /*
    tdrDraw(hr30chf,"Pz",kOpenSquare,kRed,kSolid); hr30chf->SetMarkerSize(0.6);
    tdrDraw(hr30nef,"Pz",kOpenCircle,kBlue,kSolid); hr30nef->SetMarkerSize(0.6);
    tdrDraw(hr30nhf,"Pz",kOpenDiamond,kGreen+2,kSolid); hr30nhf->SetMarkerSize(0.7);
    tdrDraw(hr110chf,"Pz",kFullSquare,kRed,kSolid); hr110chf->SetMarkerSize(0.6);
    tdrDraw(hr110nef,"Pz",kFullCircle,kBlue,kSolid); hr110nef->SetMarkerSize(0.6);
    tdrDraw(hr110nhf,"Pz",kFullDiamond,kGreen+2,kSolid); hr110nhf->SetMarkerSize(0.7);

    tdrDraw(hr30m,"Pz",kOpenDiamond,kBlack,kSolid); hr30m->SetMarkerSize(0.7);
    tdrDraw(hr110m,"Pz",kFullDiamond,kBlack,kSolid); hr110m->SetMarkerSize(0.7);
    */
    tdrDraw(hr50chf,"Pz",kFullSquare,kRed,kSolid); hr50chf->SetMarkerSize(0.6);
    tdrDraw(hr50nef,"Pz",kFullCircle,kBlue,kSolid); hr50nef->SetMarkerSize(0.6);
    tdrDraw(hr50nhf,"Pz",kFullDiamond,kGreen+2,kSolid); hr50nhf->SetMarkerSize(0.7);
    
    tdrDraw(hr50m,"Pz",kOpenDiamond,kBlack,kSolid); hr50m->SetMarkerSize(0.7);

    hr50m->GetXaxis()->SetLabelSize(0.15);

    // Add legend
    c1->cd(1);
    //TLegend *leg = tdrLeg(0.54,0.88-8*0.04,0.64,0.88);
    TLegend *leg = tdrLeg(0.70,0.88-4*0.04,0.80,0.88); //4*0.4 for each leg entry 0.4
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    /*
    leg->AddEntry(hr110m,"MPF 110EB","PLE");
    leg->AddEntry(hr110chf,"CHF 110EB","PLE");
    leg->AddEntry(hr110nef,"NEF 110EB","PLE");
    leg->AddEntry(hr110nhf,"NHF 110EB","PLE");
    leg->AddEntry(hr30m,"MPF 30EB","PLE");
    leg->AddEntry(hr30chf,"CHF 30EB","PLE");
    leg->AddEntry(hr30nef,"NEF 30EB","PLE");
    leg->AddEntry(hr30nhf,"NHF 30EB","PLE");
    */

    leg->AddEntry(hr50m, Form("MPF 50EB %.3f",hr50moff), "PLE");
    leg->AddEntry(hr50chf, Form("CHF 50EB %.3f",hr50chfoff), "PLE");
    leg->AddEntry(hr50nef, Form("NEF 50EB %.3f",hr50nefoff), "PLE");
    leg->AddEntry(hr50nhf, Form("NHF 50EB %.3f",hr50nhfoff), "PLE");
 
    string ending;
    if(zoom){
	    ending="_zoomed";
    }
    else{
	    ending="";
    }

    if(rereco==1){
    //c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2025only_%s.pdf",cv));
    //c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2024-ECALRATIO_%s_zoomed.pdf",cv));
    	c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2024-ECALRATIO_%s%s.pdf",cv,ending.c_str()));
    //c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2025only_%s_new.pdf",cv));
    }
    if(rereco==2){
    	//c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2024C-ECALR-HCALDI_%s_zoomed.pdf",cv));
	c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2024C-ECALR-HCALDI_%s%s.pdf",cv,ending.c_str()));
    }
    if(rereco==3){
	c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2024CDE-rereco_2024FGHI-prompt_%s%s.pdf",cv,ending.c_str()));
    }
    else{
     	c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2025only_%s.pdf",cv));
    	//c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2025only_%s_zoomed.pdf",cv));
    	//c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2025only_%s%s.pdf",cv,ending.c_str()));
    	//c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2024BCDEFGHI-prompt_%s%s.pdf",cv,ending.c_str()));

    //c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2025only_%s_new.pdf",cv));
    }


} // drawPFcompVsRun_2025only


// TO DO: instead of using this extra function, should just call the PFcompVsRun and call it with a different dataset...
//start of PF composition ENERGY FRACTION BALANCE  versus run, modified for 2024 only; for now focus on 50EB trigger
void drawPFcompEFBvsRun_2025only(string version, int rereco, bool zoom) {

    const char *cv = version.c_str();

    setTDRStyle();
    TDirectory *curdir = gDirectory;

    // Open input files
    //TFile *f = new TFile(Form("rootfiles/GamHistosFill_data_2025only_%s.root",cv),"READ");
    TFile *f(0);

    if(rereco==0){
	cout << "Looking at 2025 prompt data." << endl;
	f = new TFile(Form("rootfiles/GamHistosFill_data_2025CDEF_%s.root",cv), "READ"); 
	    //f = new TFile(Form("rootfiles/GamHistosFill_data_2024BCDEFGHI-prompt_%s.root",cv), "READ"); //prompt, nib-based
    }
    else if(rereco==1){
	    cout << "Looking at 2024 rereco (ECALRATIO) data." << endl;
    	f = new TFile(Form("rootfiles/GamHistosFill_data_2024-ECALRATIO_%s.root",cv), "READ"); //for now: hadd on 2024B and 2024C (need to redo this file with new daily json)
    }
    else if(rereco==2){
	    cout << "Looking at 2024C 2nd rereco (ECALR-HCALDI) data." << endl;
    	f = new TFile(Form("rootfiles/GamHistosFill_data_2024C-ECALR-HCALDI_%s.root",cv), "READ"); //
    }
    else if(rereco==3){
	    cout << "Looking at 2024CDE-rereco and 2024FGHI-prompt data." << endl;
	    f = new TFile(Form("rootfiles/GamHistosFill_data_2024CDE-rereco_2024FGHI-prompt_%s.root",cv), "READ"); //containing 2024CDE-rereco and 2024FGHI-prompt
    }
    else {
	    cout << "Check rereco boolian. Something's wrong." << endl;
	    exit(0);
    }
    assert(f && !f->IsZombie());
    f->cd("runs");
    TDirectory *d = gDirectory;

    // Load input profiles; should adjust error cleaning and also add this to legend as info
    /*
    TProfile *pr50m = (TProfile*)d->Get("pr50m"); clean(pr50m,0.006);
    TProfile *pr50chf = (TProfile*)d->Get("pr50chf"); clean(pr50chf,0.006);
    TProfile *pr50nhf = (TProfile*)d->Get("pr50nhf"); clean(pr50nhf,0.006);
    TProfile *pr50nef = (TProfile*)d->Get("pr50nef"); clean(pr50nef,0.006);
    */

    TProfile *pr50m = (TProfile*)d->Get("pr50m"); clean(pr50m,0.006);
    TProfile *pr50efb_chf = (TProfile*)d->Get("pr50efb_chf"); clean(pr50efb_chf,0.006);
    TProfile *pr50efb_nhf = (TProfile*)d->Get("pr50efb_nhf"); clean(pr50efb_nhf,0.006);
    TProfile *pr50efb_nef = (TProfile*)d->Get("pr50efb_nef"); clean(pr50efb_nef,0.006);
 

    curdir->cd();
    TH1D *hr50m = pr50m->ProjectionX("hr50m");          //what happens here?!
    TH1D *hr50efb_chf = pr50efb_chf->ProjectionX("hr50efb_chf");
    TH1D *hr50efb_nhf = pr50efb_nhf->ProjectionX("hr50efb_nhf");
    TH1D *hr50efb_nef = pr50efb_nef->ProjectionX("hr50efb_nef");
 

    // Offset composition
		double hr50moff, hr50chfoff, hr50nefoff, hr50nhfoff;
		if(zoom){
			hr50moff = -1.00; //in V2M had this as -1.01
    	hr50chfoff = -0.665;
    	hr50nefoff = -0.23;
    	hr50nhfoff = -0.075;
		}
		else{
      /*
    	hr50moff = -1.00;
    	hr50chfoff = -0.68;
    	hr50nefoff = -0.20;
    	hr50nhfoff = -0.11;
      */
    	hr50moff = -1.00;
    	hr50chfoff = -0.7;
    	hr50nefoff = -0.25;
    	hr50nhfoff = -0.07;
	}

    addOffset(hr50m,hr50moff); //need to adjust this and the following still
    addOffset(hr50efb_chf,hr50chfoff);
    addOffset(hr50efb_nef,hr50nefoff);
    addOffset(hr50efb_nhf,hr50nhfoff);
 

    // Setup canvas
    TH1D *h;
    if(rereco==3){ //starting from 24C
      //should go from -0.10 to +0.10
    	//h = tdrHist("h","PF composition EFB offset",-0.10,+0.10,"Run",379300,387100); // after removing 24B (start from C-rereco)
    	//h = tdrHist("h","PF composition EFB offset",-0.10,+0.10,"Run",379300,387100); // after removing 24B (start from C-rereco)
    }
    else{ //starting from 24B
	//h = tdrHist("h","PF composition EFB offset",-0.10,+0.10,"Run",378900,387100);
	//h = tdrHist("h","PF composition EFB offset",-0.10,+0.10,"Run",391500,393200);//w54
	//h = tdrHist("h","PF composition EFB offset",-0.10,+0.10,"Run",391500,396000);//w58
	h = tdrHist("h","PF composition EFB offset",-0.10,+0.10,"Run",391500,397300);//w58
    }

    //TH1D *h = tdrHist("h2","PF composition offset",-0.80,+0.80,"Run",378900,380600);
    //h->GetXaxis()->SetLabelSize(0.15); //smaller labels
    h->GetXaxis()->SetLabelSize(0.26);
    h->GetYaxis()->SetLabelSize(0.036);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    gPad->Update();

		

    //lumi_136TeV = Form("Photon+jet, Run 3, %s",cv);
    if(rereco==1){
	    lumi_136TeV = Form("Photon+jet, 2024 ECALRATIO, %s",cv);
    }
    else if(rereco==2){
	    lumi_136TeV = Form("Photon+jet, 2024 ECALR-HCALDI, %s",cv);
    }
    else if(rereco==3){
	    lumi_136TeV = Form("#gamma+jet, 2024CDE-rereco 2024FGHI-prompt, %s",cv);
    }
    else{
    	//lumi_136TeV = Form("#gamma+jet, 2024BCDEFGHI-prompt, %s",cv);
    	lumi_136TeV = Form("#gamma+jet, 2025CDEF-prompt, %s",cv);
    }
    extraText = "Private";
    TCanvas *c1 = tdrCanvas("c2",h,8,11);
    TLine *l = new TLine();
    TLatex *t = new TLatex();
    //t->SetTextSize(0.045);
		t->SetTextSize(0.036);


    //added this to composition plot.
    TLatex *extratext = new TLatex();
    extratext->SetTextSize(0.039);
    extratext->SetTextFont(52);
    //extratext->DrawLatex(379575, 0.0715, "Preliminary");
    extratext->DrawLatex(391700, 0.0715, "Preliminary"); //25B start 391531, w58



    TLatex *efbinfo = new TLatex();
    efbinfo->SetTextSize(0.039);
    efbinfo->SetTextFont(52);
    efbinfo->SetTextColor(kMagenta);
    efbinfo->DrawLatex(381560, 0.0765, "Energy Fraction Balance");



    h->GetXaxis()->SetLabelSize(0.15); //smaller labels
    gPad->Update();
	

    // Start drawing
    c1->cd();
    gPad->SetLogx();
    double x1(h->GetXaxis()->GetXmin()), x2(h->GetXaxis()->GetXmax());
    double y1(h->GetMinimum()), y2(h->GetMaximum());
    l->DrawLine(x1,0,x2,0);

    // 2024 Era definition
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#2024_Era_definition
    double run24b_start(378971), run24b_end(379411);
    double run24c_start(379412), run24c_end(380252);
    double run24d_start(380253), run24d_end(380947);
    double run24ev1_start(380948), run24ev1_end(381383);
    double run24ev2_start(381384), run24ev2_end(381943); //(first run#), ends?
    double run24f_start(381944), run24f_end(383779);
    double run24g_start(383780), run24g_end(385813); //(first run#), ends?
    double run24h_start(385814), run24h_end(386402); //(first run#), ends?
    double run24i_start(386403);

    //2025
    double run25b_start(391531), run25b_end(392158); //end of 24B? - but it based on start of C... Bstart:391531
    double run25cv1_start(392159), run25cv1_end(393108); //run25c_end(392524); //run25c_end(392441); //ends?
    double run25cv2_start(393111), run25cv2_end(393609);
    double run25d_start(394286), run25d_end(395967); //updated end of D (30.09.)
    double run25e_start(395968), run25e_end(396597); //
    double run25f_start(396598), run25f_end(397097); // from my json used on 21st of Sep
 
 


    //l->SetLineStyle(kDashed);
    l->SetLineColor(kGray);
    if(rereco!=3){
	l->DrawLine(run24b_start,y1,run24b_start,y2); 	//removed B on 06.04.2025
	l->DrawLine(run24b_end,y1,run24b_end,y2);	//removed B on 06.04.2025
	//t->DrawLatex(run24b_start+250,-0.095,"24B"); //+250 offset only when displaying in middle of era
    }
    l->DrawLine(run24c_start,y1,run24c_start,y2);
    l->DrawLine(run24c_end,y1,run24c_end,y2);
    l->DrawLine(run24d_start,y1,run24d_start,y2);
    l->DrawLine(run24d_end,y1,run24d_end,y2);
    l->DrawLine(run24ev1_start,y1,run24ev1_start,y2);
    l->DrawLine(run24ev2_start,y1,run24ev2_start,y2);
    l->DrawLine(run24f_start,y1,run24f_start,y2);
    l->DrawLine(run24g_start,y1,run24g_start,y2);
    l->DrawLine(run24h_start,y1,run24h_start,y2);
    l->DrawLine(run24i_start,y1,run24i_start,y2);

    //2025
    l->DrawLine(run25b_start,y1,run25b_start,y2);
    l->DrawLine(run25b_end,y1,run25b_end,y2); //to separate visually 25B and 25C
    l->DrawLine(run25cv1_start,y1,run25cv1_start,y2); //used to draw this -15 before the first data point of this era...
    l->DrawLine(run25cv1_end,y1,run25cv1_end,y2);
    l->DrawLine(run25cv2_start,y1,run25cv2_start,y2); //used to draw this -15 before the first data point of this era...
    l->DrawLine(run25cv2_end,y1,run25cv2_end,y2);
    l->DrawLine(run25d_start,y1,run25d_start,y2); 
    l->DrawLine(run25d_end,y1,run25d_end,y2);
    l->DrawLine(run25e_start,y1,run25e_start,y2); 
    l->DrawLine(run25e_end,y1,run25e_end,y2);
    l->DrawLine(run25f_start,y1,run25f_start,y2); 
    l->DrawLine(run25f_end,y1,run25f_end,y2);







 
    //text for eras
    double textposy;
    if(zoom){
	    textposy=-0.020;
    }
    else{
	    textposy=-0.090;
    }
    if(rereco!=3){ t->DrawLatex(run24b_start+10,textposy,"24B"); }
    t->DrawLatex(run24c_start+10,textposy,"24C");
    t->DrawLatex(run24d_start+10,textposy,"24D");
    //t->DrawLatex(run24ev1_start+10,textposy,"24Ev1,v2");
    t->DrawLatex(run24ev1_start+10,textposy,"24Ev1,");
    t->DrawLatex(run24ev1_start+450,textposy-0.008,"v2");
    //t->DrawLatex(run24ev2_start+10,textposy-0.01,"24Ev2");
    t->DrawLatex(run24f_start+10,textposy,"24F");
    t->DrawLatex(run24g_start+10,textposy,"24G");
    t->DrawLatex(run24h_start+10,textposy,"24H");
    t->DrawLatex(run24i_start+10,textposy,"24I");

    //2025
    t->DrawLatex(run25b_start+10,textposy,"25B");
    t->DrawLatex(run25cv1_start+10,textposy,"25Cv1");
    t->DrawLatex(run25cv2_start+10,textposy,"25Cv2");
    t->DrawLatex(run25d_start+10,textposy,"25D");
    t->DrawLatex(run25e_start+10,textposy,"25E");
    t->DrawLatex(run25f_start+10,textposy,"25F");




    //Drawing everything
    tdrDraw(hr50efb_chf,"Pz",kFullSquare,kRed,kSolid); hr50efb_chf->SetMarkerSize(0.6);
    tdrDraw(hr50efb_nef,"Pz",kFullCircle,kBlue,kSolid); hr50efb_nef->SetMarkerSize(0.6);
    tdrDraw(hr50efb_nhf,"Pz",kFullDiamond,kGreen+2,kSolid); hr50efb_nhf->SetMarkerSize(0.7);
    
    tdrDraw(hr50m,"Pz",kOpenDiamond,kBlack,kSolid); hr50m->SetMarkerSize(0.7);

    hr50m->GetXaxis()->SetLabelSize(0.15);

    // Add legend
    c1->cd(1);
    //TLegend *leg = tdrLeg(0.54,0.88-8*0.04,0.64,0.88);
    TLegend *leg = tdrLeg(0.70,0.88-4*0.04,0.80,0.88); //4*0.4 for each leg entry 0.4
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);

    leg->AddEntry(hr50m, Form("MPF 50EB %.3f",hr50moff), "PLE");
    leg->AddEntry(hr50efb_chf, Form("CHF 50EB %.3f",hr50chfoff), "PLE");
    leg->AddEntry(hr50efb_nef, Form("NEF 50EB %.3f",hr50nefoff), "PLE");
    leg->AddEntry(hr50efb_nhf, Form("NHF 50EB %.3f",hr50nhfoff), "PLE");
 
    string ending;
    if(zoom){
	    ending="_zoomed";
    }
    else{
	    ending="";
    }

    if(rereco==1){
    	c1->SaveAs(Form("pdf/drawResponseVsRun_PFcompEFB_2024-ECALRATIO_%s%s.pdf",cv,ending.c_str()));
    }
    if(rereco==2){
	    c1->SaveAs(Form("pdf/drawResponseVsRun_PFcompEFB_2024C-ECALR-HCALDI_%s%s.pdf",cv,ending.c_str()));
    }
    if(rereco==3){
	    c1->SaveAs(Form("pdf/drawResponseVsRun_PFcompEFB_2024CDE-rereco_2024FGHI-prompt_%s%s.pdf",cv,ending.c_str()));
    }
    else{
    	//c1->SaveAs(Form("pdf/drawResponseVsRun_PFcompEFB_2024BCDEFGHI-prompt_%s%s.pdf",cv,ending.c_str()));
    	c1->SaveAs(Form("pdf/drawResponseVsRun_PFcompEFB_2025CDEF_%s%s.pdf",cv,ending.c_str()));
    }


} // drawPFcompVsRun_2025only

