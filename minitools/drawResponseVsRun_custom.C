// Purpose: draw photon+jet response vs run to check time stability
// _custom.C --> This is a copy of the original plotting script used to try out some stuff and removed unneeded things (can always look them up from the original code)
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

void drawPFcompVsRun(string version);

// v26->v29->(v30)->v31(after L2L3Res_V2)->v32(L2L3Res_V3) ==> these are Mikko's versions, mine are with "wX"
//void drawResponseVsRun_custom(string version = "w10") { //w10 for 2023 data, w11 for the new 2024 data
//void drawResponseVsRun_custom(string version = "w12", string year=2024) {//for plotting only one year
void drawResponseVsRun_custom(string version = "w12") {

    const char *cyear = year.c_str();
    const char *cv = version.c_str();

    setTDRStyle();
    TDirectory *curdir = gDirectory;

    // Open input files
    //TFile *f = new TFile(Form("rootfiles/GamHistosFill_data_Run3_%s.root",cv), "READ"); //is this the entire run 3 data with all runperiods?
    TFile *f = new TFile(Form("rootfiles/GamHistosFill_data_Run3Summer23_%s.root",cv), "READ"); //for now: hadd on 2023 Cv123, Cv4, D

    assert(f && !f->IsZombie());

    //f->cd("control"); // v29
    f->cd("runs"); // v30
    TDirectory *d = gDirectory;

    // Load input profiles + clean the errors
    TProfile *pr30b = (TProfile*)d->Get("pr30b"); clean(pr30b,0.006); //b = balance
    TProfile *pr30m = (TProfile*)d->Get("pr30m"); clean(pr30m,0.006); //m = mpf
    TProfile *pr110b = (TProfile*)d->Get("pr110b"); clean(pr110b,0.003);
    TProfile *pr110m = (TProfile*)d->Get("pr110m"); clean(pr110m,0.003);

    TProfile *pr50b, *pr50m;
    { // add 2024 --> this is so far added separately while the rest is read from one run3-file. could combine it.
        TFile *f24 = new TFile("rootfiles/GamHistosFill_data_2024B-PromptReco-v1_w11.root", "READ");
        assert(f24 && !f24->IsZombie());
        TDirectory *d24 = f24->GetDirectory("runs");

        TProfile *pr30b_24, *pr30m_24, *pr110b_24, *pr110m_24;
        pr30b_24 = (TProfile*)d24->Get("pr30b"); clean(pr30b_24,0.01); //shiftX(pr30m_24); //test
        pr30m_24 = (TProfile*)d24->Get("pr30m"); clean(pr30m_24,0.01);
        pr50b = (TProfile*)d24->Get("pr50b"); clean(pr50b,0.01); //shiftX(pr50b); //test
        pr50m = (TProfile*)d24->Get("pr50m"); clean(pr50m,0.01); //shiftX(pr50m); //test
        pr110b_24 = (TProfile*)d24->Get("pr110b"); clean(pr110b_24,0.01); //shiftX(pr110b_24); //test
        pr110m_24 = (TProfile*)d24->Get("pr110m"); clean(pr110m_24,0.01); //shiftX(pr110m_24); //test


        pr30b = (TProfile*)hadd("pr30b2",pr30b_24, pr30b);
        pr30m = (TProfile*)hadd("pr30m2",pr30m_24, pr30m);
        pr110b = (TProfile*)hadd("pr100b2",pr110b_24, pr110b);
        pr110m = (TProfile*)hadd("pr110m2",pr110m_24, pr110m);
    }
  
    // Scale BAL and pT30 results
    double kbal = 1.13;
    pr30b->Scale(kbal);
    pr110b->Scale(kbal);
    double kpt30 = 0.98;
    pr30b->Scale(kpt30);
    pr30m->Scale(kpt30);
    double kbal30 = 0.96;
    pr30b->Scale(kbal30);
    if (pr50b && pr50m) {
        pr50b->Scale(kbal);
        double kpt50 = 0.99;//0.98;
        pr50b->Scale(kpt50);
        pr50m->Scale(kpt50);
        double kbal50 = 0.91;//kbal30;
        pr50b->Scale(kbal50);
    }

    // Setup canvas
    //TH1D *h = tdrHist("h","Response",0.8,1.2,"Run",355300,371300);
    //TH1D *h = tdrHist("h","Response",0.92,1.08,"Run",355300,371300);
    //TH1D *h = tdrHist("h","Response",0.92,1.08,"Run",355300,381300); //wherefrom do i know these run numbers -> link
    TH1D *h = tdrHist("h","Response",0.92,1.08,"Run",366000,381300); //should start later when dropping 2022 data (here some before 2023C)
    //lumi_136TeV = Form("Photon+jet, Run 3, %s",cv);
    lumi_136TeV = Form("Photon+jet, Run 3 2023, %s",cv); //for 2023-only version
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

    // 2023 Era definition
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#2023_Era_definition
    //comment out 2022 FOR NOW --> TODO: run 2022 with w12 code and include all.
    /*
    double run22c1b(355794), run22c1e(357486);
    l->DrawLine(run22c1b,y1,run22c1b,y2);
    l->SetLineStyle(kDotted);
    l->DrawLine(run22c1e,y1,run22c1e,y2);
    t->DrawLatex(run22c1b+350,0.925,"22C");//"22C");
    double run22d1b(357487), run22d1e(357733);
    double run22d2b(357734), run22d2e(358219);
    double run22d3b(358220), run22d3e(359021);
    l->DrawLine(run22d1b,y1,run22d1b,y2);
    */

    //this was aleady commented out
    /*
    l->SetLineStyle(kDashed);
    l->DrawLine(run22d1e,y1+0.015,run22d1e,y2);
    l->DrawLine(run22d2b,y1+0.015,run22d2b,y2);
    l->DrawLine(run22d2e,y1+0.015,run22d2e,y2);
    l->DrawLine(run22d3b,y1+0.015,run22d3b,y2);
    l->SetLineStyle(kSolid);
    */

    //comment out 2022 FOR NOW --> TODO: run 2022 with w12 code and include all.
    /*
    l->DrawLine(run22d3e,y1,run22d3e,y2);
    t->DrawLatex(run22d1b+300,0.925,"D");//"22D");
    double run22e1b(359022), run22e1e(360331);
    l->DrawLine(run22e1b,y1,run22e1b,y2);
    l->DrawLine(run22e1e,y1,run22e1e,y2);
    t->DrawLatex(run22e1b+150,0.925,"E");//"22E");
    l->SetLineStyle(kSolid);
    double run22f1b(360332), run22f1e(362180);
    l->DrawLine(run22f1b,y1,run22f1b,y2);
    l->SetLineStyle(kDotted);
    l->DrawLine(run22f1e,y1,run22f1e,y2);
    t->DrawLatex(run22f1b+450,0.925,"22F");
    double run22g1b(362350), run22g1e(362760);
    l->DrawLine(run22g1b,y1,run22g1b,y2);
    l->SetLineStyle(kSolid);
    l->DrawLine(run22g1e,y1,run22g1e,y2);  
    t->DrawLatex(run22g1b+50,0.925,"G");
    */
    // 2023
    double run23c1b(367080), run23c1e(367515);
    double run23c2b(367516), run23c2e(367620);
    double run23c3b(367621), run23c3e(367763);
    double run23c4b(367765), run23c4e(369802);
    l->DrawLine(run23c1b,y1,run23c1b,y2);
    /*
    l->SetLineStyle(kDashed);
    l->DrawLine(run23c1e,y1+0.015,run23c1e,y2);
    l->DrawLine(run23c2b,y1+0.015,run23c2b,y2);
    l->DrawLine(run23c2e,y1+0.015,run23c2e,y2);
    l->DrawLine(run23c3b,y1+0.015,run23c3b,y2);
    l->SetLineStyle(kSolid);
    */
    l->SetLineStyle(kDotted);
    l->DrawLine(run23c3e,y1+0.015,run23c3e,y2);
    l->DrawLine(run23c4b,y1+0.015,run23c4b,y2);
    l->DrawLine(run23c4e,y1,run23c4e,y2);
    t->DrawLatex(run23c1b+250,0.925,"23C");//"23C");
    t->DrawLatex(run23c4b+250,0.935,"v4");
    double run23d1b(369803), run23d1e(370602);
    double run23d2b(370603), run23d2e(372415);
    l->DrawLine(run23d1b,y1,run23d1b,y2);
    /*
    l->SetLineStyle(kDashed);
    l->DrawLine(run23d1e,y1+0.015,run23d1e,y2);
    l->DrawLine(run23d2b,y1+0.015,run23d2b,y2);
    */
    l->SetLineStyle(kSolid);
    l->DrawLine(run23d2e,y1,run23d2e,y2);
    t->DrawLatex(run23d1b+250,0.925,"D");//"23D");

    double run24(378981), run24e(379154);
    l->DrawLine(run24,y1+0.035,run24,y2-0.050);
    l->DrawLine(run24e,y1+0.035,run24e,y2-0.050);
    t->DrawLatex(run24,0.945,"24B");
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
    //tdrDraw(pr30b,"Pz",kOpenSquare,kBlue,kSolid); pr30b->SetMarkerSize(0.7);
    tdrDraw(pr30m,"Pz",kFullSquare,kBlue,kSolid); pr30m->SetMarkerSize(0.6);
    //tdrDraw(pr110b,"Pz",kOpenCircle,kRed,kSolid); pr110b->SetMarkerSize(0.7);
    tdrDraw(pr110m,"Pz",kFullCircle,kRed,kSolid); pr110m->SetMarkerSize(0.6);

    // Add legend
    //c1->cd(1);
    //TLegend *leg = tdrLeg(0.54,0.68,0.64,0.88);
    TLegend *leg = tdrLeg(0.74,0.78,0.84,0.88);
    leg->AddEntry(pr110m,"MPF 110EB","PLE");
    //leg->AddEntry(pr110b,"BAL 110EB","PLE");
    leg->AddEntry(pr30m,"MPF 30EB","PLE");
    //leg->AddEntry(pr30b,"BAL 30EB","PLE");

    if (pr50b && pr50m) {
        //tdrDraw(pr50b,"Pz",kOpenDiamond,kGreen+2,kSolid);//pr50b->SetMarkerSize(0.7);
        tdrDraw(pr50m,"Pz",kFullDiamond,kGreen+2,kSolid);//pr50m->SetMarkerSize(0.6);
        tdrDraw(pr110m,"Pz",kFullCircle,kRed,kSolid); // Redraw
        //TLegend *leg2 = tdrLeg(0.74,0.13,0.84,0.23);
        //leg2->AddEntry(pr50m,"MPF 50EB","PLE");
        //leg2->AddEntry(pr50b,"BAL 50EB","PLE");
        leg->AddEntry(pr50m,"MPF 50EB","PLE");
        leg->SetY1(leg->GetY1()-0.05);
    }


  
    //c1->SaveAs(Form("pdf/drawResponseVsRun_custom_%s.pdf",cv));
    c1->SaveAs(Form("pdf/drawResponseVsRun_custom_%s_TEST.pdf",cv));


    // Extra composition plots
    drawPFcompVsRun(version);
} // drawResponseVsRun_custom


//for calling function and running this PFComp code with certain code version from command line:
// look up code in original version of this script, see function drawPFcompVsRun(string version), removed it from here
