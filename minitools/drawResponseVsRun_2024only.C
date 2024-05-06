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


//void drawResponseVsRun_custom(string version = "w10") { //w10 for 2023 data, w11 for the new 2024 data
//void drawResponseVsRun_custom(string version = "w12", string year=2024) {//for plotting only one year
void drawResponseVsRun_2024only(string version = "w13-1") {

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
    double maxerr = 0.3; //apply same maxerr cleaning to all paths, so it is less confusing
    //TProfile *pr30b = (TProfile*)d->Get("pr30b"); clean(pr30b,maxerr); //b = balance
    //TProfile *pr30m = (TProfile*)d->Get("pr30m"); clean(pr30m,maxerr); //m = mpf
    TProfile *pr50b = (TProfile*)d->Get("pr50b"); clean(pr50b,maxerr); 
    TProfile *pr50m = (TProfile*)d->Get("pr50m"); clean(pr50m,maxerr); 
    TProfile *pr110b = (TProfile*)d->Get("pr110b"); clean(pr110b,maxerr);
    TProfile *pr110m = (TProfile*)d->Get("pr110m"); clean(pr110m,maxerr);
    //TProfile *pr230b = (TProfile*)d->Get("pr230b"); //clean(pr230b,0.003);
    //TProfile *pr230m = (TProfile*)d->Get("pr230m"); //clean(pr230m,0.003);

 
    // Scale BAL and pT30 results
    double kbal = 1.13;
    //pr30b->Scale(kbal);
    pr110b->Scale(kbal);
    double kpt30 = 0.98;
    //pr30b->Scale(kpt30);
    //pr30m->Scale(kpt30);
    double kbal30 = 0.96;
    //pr30b->Scale(kbal30);
    if (pr50b && pr50m) { //in 2024 always
        pr50b->Scale(kbal);
        double kpt50 = 0.99;//0.98;
        pr50b->Scale(kpt50);
        pr50m->Scale(kpt50);
        double kbal50 = 0.91;//kbal30;
        pr50b->Scale(kbal50);
    }
    //NOTE: DO WE ALSO NEED TO SCALE 230GeV results? (in case we include those)

    // Setup canvas
    //TH1D *h = tdrHist("h","Response",0.92,1.08,"Run",366000,381300); //should start later when dropping 2022 data (here some before 2023C)
    TH1D *h = tdrHist("h","Response",0.86,1.06,"Run",378900,380200); //2024 only (standard was to display from 0.92 to 1.08, now zooming in for 2024
    //lumi_136TeV = Form("Photon+jet, Run 3, %s",cv);
    lumi_136TeV = Form("Photon+jet, Run 3 2024, %s",cv); //for 2024-only version
    extraText = "Private";
    TCanvas *c1 = tdrCanvas("c1",h,8,11);
    TLine *l = new TLine();
    TLatex *t = new TLatex();
    t->SetTextSize(0.045);// t->SetNDC();

    // Start drawing
    c1->cd();
    gPad->SetLogx();
    double x1(h->GetXaxis()->GetXmin()), x2(h->GetXaxis()->GetXmax());
    double y1(h->GetMinimum()), y2(h->GetMaximum());
    l->DrawLine(x1,1,x2,1);

    l->SetLineStyle(kDotted); 

    // 2024 (start and end of runs)
    double run24b_start(378981), run24b_end(379356); //379355 end of 24B? 
    l->DrawLine(run24b_start-15,y1+0.040,run24b_start-15,y2-0.045);
    //l->DrawLine(run24b_end,y1+0.035,run24b_end,y2-0.050);
    t->DrawLatex(run24b_start,0.960,"24B");
    l->DrawLine(run24b_end+15,y1+0.040,run24b_end+15,y2-0.045); //to separate visually 24B and 24C


    double run24c_start(379415), run24c_end(380100); //actually starting at 379415 (first run#), ends?
    l->DrawLine(run24c_start-15,y1+0.040,run24c_start-15,y2-0.045);
    l->DrawLine(run24c_end,y1+0.035,run24c_end,y2-0.050);
    t->DrawLatex(run24c_start,0.960,"24C");

    //info on the dotted lines
    TLatex *infotext = new TLatex();
    infotext->SetNDC(kTRUE);
    infotext->SetTextSize(0.036);
    infotext->SetTextColor(12);
    infotext->SetTextFont(52);
    infotext->DrawLatex(0.53,0.18, Form("Error cleaning applied (%.2f).",maxerr));
    infotext->DrawLatex(0.53,0.13, Form("The dotted line is: start-15 and end+15."));


    //pr50m->GetXaxis()->SetAxisColor(kRed);


    //tdrDraw(pr30b,"Pz",kOpenSquare,kBlue,kSolid); pr30b->SetMarkerSize(0.7);
    //tdrDraw(pr30m,"Pz",kFullCircle,kBlue,kSolid); pr30m->SetMarkerSize(0.6);
    //tdrDraw(pr110b,"Pz",kOpenCircle,kRed,kSolid); pr110b->SetMarkerSize(0.7);
    tdrDraw(pr50b,"Pz",kOpenSquare,kBlue,kSolid); pr50b->SetMarkerSize(0.5);            //balance --> instead of MPF30, showing DB50
    tdrDraw(pr110m,"Pz",kFullCircle,kRed,kSolid); pr110m->SetMarkerSize(0.5);           //print this one first, as 50EB captures also everything that goes through 110EB
    tdrDraw(pr50m,"Pz",kFullCircle,kGreen+2,kSolid); pr50m->SetMarkerSize(0.5);         //marker size 0.06 or 0.05?

    //set the axis labels custom?
    TAxis *xaxis = pr50b->GetXaxis();
    xaxis->SetAxisColor(kGreen);
    //xaxis->LabelsOption("d");
    //h->LabelsOption("d","X");
    //xaxis->Draw();
    //pr50m->SetXLabels();
    pr50m->LabelsOption("d","X");

    //tdrDraw(pr230m,"Pz",kFullCircle,kGray,kSolid); pr230m->SetMarkerSize(0.6);


    // Add legend
    //c1->cd(1);
    //TLegend *leg = tdrLeg(0.54,0.68,0.64,0.88);
    TLegend *leg = tdrLeg(0.74,0.78,0.84,0.88);
    leg->AddEntry(pr110m,"MPF 110EB","PLE");
    //leg->AddEntry(pr110b,"BAL 110EB","PLE");
    leg->AddEntry(pr50m,"MPF 50EB","PLE");
    leg->AddEntry(pr50b,"BAL 50EB","PLE");
    //leg->AddEntry(pr30m,"MPF 30EB","PLE");


    //leg->AddEntry(pr30b,"BAL 30EB","PLE");

    c1->SaveAs(Form("pdf/drawResponseVsRun_2024only_%s.pdf",cv));


    // Extra composition plots dropped --> see old script for those (drawPFcompVsRun(version);)
} // drawResponseVsRun_2024only
