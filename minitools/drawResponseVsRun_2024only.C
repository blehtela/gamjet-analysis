// Purpose: draw photon+jet response vs run to check time stability
// This is the code for plotting 2024 only... should add such a switch to the general code.
#include "TFile.h"
#include "TProfile.h"
#include "TLatex.h"

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


void drawPFcompVsRun_2024only(string version);


//void drawResponseVsRun_custom(string version = "w10") { //w10 for 2023 data, w11 for the new 2024 data
//void drawResponseVsRun_custom(string version = "w12", string year=2024) {//for plotting only one year
void drawResponseVsRun_2024only(string version = "w17") { //switched to w15 and w16; w17 and w18

    //const char *cyear = year.c_str();
    const char *cv = version.c_str();

    setTDRStyle();
    TDirectory *curdir = gDirectory;

    // Open input files
    //TFile *f = new TFile(Form("rootfiles/GamHistosFill_data_Run3_%s.root",cv), "READ"); //is this the entire run 3 data with all runperiods?
    //TFile *f = new TFile(Form("rootfiles/GamHistosFill_data_Run3Summer23_%s.root",cv), "READ"); //for now: hadd on 2023 Cv123, Cv4, D
    TFile *f = new TFile(Form("rootfiles/GamHistosFill_data_2024only_%s.root",cv), "READ"); //for now: hadd on 2024B and 2024C (need to redo this file with new daily json)

    assert(f && !f->IsZombie());

    //f->cd("control"); // v29
    f->cd("runs"); // v30
    TDirectory *d = gDirectory;

    // Load input profiles + clean the errors (what should i set the maxerr to? tried to be generous now)
    double maxerr = 0.01; //apply same maxerr cleaning to all paths, so it is less confusing
    //TProfile *pr30b = (TProfile*)d->Get("pr30b"); clean(pr30b,maxerr); //b = balance
    //TProfile *pr30m = (TProfile*)d->Get("pr30m"); clean(pr30m,maxerr); //m = mpf
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
    double kbal50 = 1.11;
    pr50b->Scale(kbal50); //now also scaling this
    double kpt50 = 0.99;
    pr50m->Scale(kpt50);
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
    TH1D *h = tdrHist("h","Response",0.92,1.06,"Run",378900,380800); //2024 only (standard was to display from 0.92 to 1.08, now zooming in for 2024
    //lumi_136TeV = Form("Photon+jet, Run 3, %s",cv);
    lumi_136TeV = Form("Photon+jet, Run 3 2024, %s",cv); //for 2024-only version
    extraText = "Private"; //would print it below CMS logo, but i want it to the right (save space)

    //set the axis labels custom?
    //TAxis *xaxis = h->GetXaxis();
    //xaxis->SetAxisColor(kGreen);
    //gPad->Update();
 

    TCanvas *c1 = tdrCanvas("c1",h,8,11);
    TLine *l = new TLine();
    TLatex *t = new TLatex();
    t->SetTextSize(0.045);// t->SetNDC();

    //custom place for "Private" extraText (usually handled via tdrCanvas)
    TLatex *extratext = new TLatex();
    extratext->SetTextSize(0.039);
    extratext->SetTextFont(52);
    extratext->DrawLatex(379226, 1.0465, "Preliminary");

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
    double run24b_start(378971), run24b_end(379411); //379355 end of 24B? 
    //l->DrawLine(run24b_start-15,y1+0.040,run24b_start-15,y2-0.045);
    l->DrawLine(run24b_start,y1+0.045,run24b_start,y2-0.035);
    //l->DrawLine(run24b_end,y1+0.035,run24b_end,y2-0.050);
    t->DrawLatex(run24b_start,0.960,"24B");
    l->DrawLine(run24b_end,y1+0.045,run24b_end,y2-0.035); //to separate visually 24B and 24C

    double run24c_start(379412), run24c_end(380252); //actually starting at 379415 (first run#), ends?
    l->DrawLine(run24c_start,y1+0.045,run24c_start,y2-0.035); //used to draw this -15 before the first data point of this era...
    l->DrawLine(run24c_end,y1+0.045,run24c_end,y2-0.035);
    t->DrawLatex(run24c_start,0.960,"24C");

    double run24d_start(380253); //, run24d_end(380600); //(first run#), ends?
    l->DrawLine(run24d_start,y1+0.045,run24d_start,y2-0.035);
    //l->DrawLine(run24d_end,y1+0.035,run24d_end,y2-0.050);
    t->DrawLatex(run24d_start,0.960,"24D");

    //first run number for each era
    TLatex *runstart = new TLatex();
    runstart->SetTextSize(0.026);
    runstart->SetTextColor(kBlack);
    runstart->SetTextFont(52);
    runstart->DrawLatex(run24b_start, 0.955, Form("%.0f",run24b_start));
    runstart->DrawLatex(run24c_start, 0.955, Form("%.f",run24c_start));
    runstart->DrawLatex(run24d_start, 0.955, Form("%.f",run24d_start));


    //add important run numbers
    //HCAL Updates:     https://twiki.cern.ch/twiki/bin/view/CMS/HcalDPGRun3Updates
    //HCAL scans:       https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalPhaseScan
    // so far scans were reported for: Run 379349: 46.48/pb, and Run 379350: 38.75/pb 
    //could do this in a loop in case there will be more of these
    TLine *runinfoline = new TLine();
    runinfoline->SetLineStyle(kSolid);
    runinfoline->SetLineColor(kOrange+2);
    runinfoline->DrawLine(379349,y1+0.030,379349,y2-0.030);     //hcal scan
    runinfoline->DrawLine(379350,y1+0.030,379350,y2-0.030);     //hcal scan
    TLatex *runinfo = new TLatex();
    runinfo->SetTextSize(0.026);
    runinfo->SetTextColor(kOrange+2);
    runinfo->SetTextFont(52);
    runinfo->DrawLatex(379349, 0.945, "379349: HCAL scan");
    runinfo->DrawLatex(379350, 0.940, "379350: HCAL scan");

    //strange things, where we try to find (and draw) the run-number to check back with hcal etc
    //i identified these visually and then put a line there and wrote the run number
    TLine *strangeline = new TLine();
    strangeline->SetLineStyle(kSolid);
    strangeline->SetLineColor(kRed-2);
    double strangerun(379965); //379980
    strangeline->DrawLine(strangerun,y1+0.030,strangerun,y2-0.030);     //hcal scan
    TLatex *strangeinfo = new TLatex();
    strangeinfo->SetTextSize(0.026);
    strangeinfo->SetTextColor(kRed-2);
    strangeinfo->SetTextFont(52);
    strangeinfo->DrawLatex(strangerun, 0.945, Form("%.f: What happens?",strangerun));

    //HCAL stuff taken care of before 2024B (mark with arrow to the left, remove later)
    //TArrow *hcalArrow = new TArrow(0.5, 0.3, 0.2, 0.1, 0.05, ">"); //1-4: position, 5:arrowsize, 6:option
    //hcalArrow->SetNDC(kTRUE); //set already here to put arrow more easily
    //runinfo->DrawLatex(0.4, 0.6, "hcal arrow");
    TArrow *hcalArrow = new TArrow(379360, 1.035, 379000, 1.035, 0.02, ">"); //1-4: position, 5:arrowsize, 6:option
    hcalArrow->SetLineColor(kOrange+2);
    hcalArrow->Draw();
    runinfo->DrawLatex(379060, 1.036, "377804: phaseTuning");


    //time (re-)alignments & other stuff, only put here when happening after start of 24B (378971 April 5)
    //379347 April 13 phaseTuning_HB_2024_v2, improves iphi-symmetry of timing in HB (all depths)  

    //HcalPedestals: https://twiki.cern.ch/twiki/bin/view/CMS/HcalPedestalsTagsRun3#Descriptions_of_the_HLT_Express
    //379349 (24 Apr 2024)
    ///runinfo->DrawLatex(379349, 0.930, "HcalPedestals"); //but this is same as scan?


    //info on the dotted lines
    TLatex *infotext = new TLatex();
    infotext->SetNDC(kTRUE);
    infotext->SetTextSize(0.030); //0.036
    infotext->SetTextColor(12);
    infotext->SetTextFont(52);
    infotext->DrawLatex(0.61,0.16, Form("Error cleaning applied (%.2f).",maxerr));
    //infotext->DrawLatex(0.53,0.13, Form("The dotted line is: start-15 and end+15."));
    infotext->DrawLatex(0.61,0.13, Form("The dotted line is: era start and end."));



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
    tdrDraw(pr110m,"Pz",kFullCircle,kRed,kSolid); pr110m->SetMarkerSize(0.5);           //print this one first, as 50EB captures also everything that goes through 110EB
    tdrDraw(pr50m,"Pz",kFullCircle,kGreen+2,kSolid); pr50m->SetMarkerSize(0.5);         //marker size 0.06 or 0.05?


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
    TLegend *leg = tdrLeg(0.64,0.76,0.84,0.88);

    leg->SetTextFont(42);
    leg->AddEntry(pr110m,"MPF 110EB","PLE");
    //leg->AddEntry(pr110b,"BAL 110EB","PLE");
    leg->AddEntry(pr50m, Form("MPF 50EB x%.2f", kpt50), "PLE");
    leg->AddEntry(pr50b, Form("BAL 50EB x%.2f", kbal50), "PLE");
    //leg->AddEntry(pr30m,"MPF 30EB","PLE");


    //leg->AddEntry(pr30b,"BAL 30EB","PLE");

    //h->Draw();
    c1->Update(); //updating helping?
    cout << "Labels? " << xaxis->GetLabels() << endl;


    c1->SaveAs(Form("pdf/drawResponseVsRun_2024only_%s.pdf",cv));


    // Extra composition plots also modified.
    drawPFcompVsRun_2024only(version);

} // drawResponseVsRun_2024only



//start of PF composition versus run, modified for 2024 only; for now focus on 50EB trigger
void drawPFcompVsRun_2024only(string version) {

    const char *cv = version.c_str();

    setTDRStyle();
    TDirectory *curdir = gDirectory;

    // Open input files
    TFile *f = new TFile(Form("rootfiles/GamHistosFill_data_2024only_%s.root",cv),"READ");
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


    double hr50moff = -1.01;
    double hr50chfoff = -0.665;
    double hr50nefoff = -0.23;
    double hr50nhfoff = -0.075;
/*
    double hr50moff = -1.01;
    double hr50chfoff = -0.68;
    double hr50nefoff = -0.20;
    double hr50nhfoff = -0.11;
*/ 
    addOffset(hr50m,hr50moff); //need to adjust this and the following still
    addOffset(hr50chf,hr50chfoff);
    addOffset(hr50nef,hr50nefoff);
    addOffset(hr50nhf,hr50nhfoff);
    //*/
 

    // Setup canvas
    //TH1D *h = tdrHist("h2","PF composition offset",-0.10,+0.10,"Run",378900,380800);
    TH1D *h = tdrHist("h2","PF composition offset",-0.025,+0.025,"Run",378900,380800);
    //TH1D *h = tdrHist("h2","PF composition offset",-0.80,+0.80,"Run",378900,380600);

    //lumi_136TeV = Form("Photon+jet, Run 3, %s",cv);
    lumi_136TeV = Form("Photon+jet, 2024only, %s",cv);
    extraText = "Private";
    TCanvas *c1 = tdrCanvas("c2",h,8,11);
    TLine *l = new TLine();
    TLatex *t = new TLatex();
    t->SetTextSize(0.045);

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
    double run24d_start(380253); //, run24d_end();

    //l->SetLineStyle(kDashed);
    l->SetLineColor(kGray);
    l->DrawLine(run24b_start,y1,run24b_start,y2);
    l->DrawLine(run24b_end,y1,run24b_end,y2);
    //t->DrawLatex(run24b_start+250,-0.095,"24B"); //+250 offset only when displaying in middle of era
    l->DrawLine(run24c_start,y1,run24c_start,y2);
    l->DrawLine(run24c_end,y1,run24c_end,y2);
    l->DrawLine(run24d_start,y1,run24d_start,y2);
 
    //text for eras
    /*
    t->DrawLatex(run24b_start+10,-0.090,"24B");
    t->DrawLatex(run24c_start+10,-0.090,"24C");
    t->DrawLatex(run24d_start+10,-0.090,"24D");
    */
    t->DrawLatex(run24b_start+10,-0.020,"24B");
    t->DrawLatex(run24c_start+10,-0.020,"24C");
    t->DrawLatex(run24d_start+10,-0.020,"24D");




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
 

    //c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2024only_%s.pdf",cv));
    c1->SaveAs(Form("pdf/drawResponseVsRun_PFcomp_2024only_%s_zoomed.pdf",cv));


} // drawPFcompVsRun_2024only

