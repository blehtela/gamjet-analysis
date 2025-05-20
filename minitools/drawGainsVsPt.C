// Purpose: draw Photon_seedGain responses vs pT to check their calibration
#include "TFile.h"
#include "TProfile.h"
#include "TLine.h"
#include "TF1.h"

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

//void drawGainsVsPt(string era = "2023D", string version = "w8") {
//void drawGainsVsPt(string era = "2023Cv4", string version = "w8") {
//void drawGainsVsPt(string era = "2023D", string version = "w9") {
//void drawGainsVsPt(string era = "2024BCD", string version = "w25") {
///void drawGainsVsPt(string era = "2024G", string version = "w42") {
///void drawGainsVsPt(string era = "2024F", string version = "w42") {
void drawGainsVsPt(string era = "2024I", string version = "w42") {


//void drawGainsVsPt(string era = "2023Cv4", string version = "w9") {
  //void drawGainsVsPt(string era = "Run2", string version = "v18") {

  const char *cera = era.c_str();
  const char *cv = version.c_str();
  //string eram = (era=="Run3"||era=="Run2" ? era : "2022");
  //string eram = "2023"; //hard-coded this for now...
  string eram = "2024"; //hard-coded this for now...
  const char *ceram = eram.c_str();

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Open input files
  TFile *f = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",
			    cera,cv),"READ");
  if (!f || f->IsZombie())
    f = new TFile(Form("files/GamHistosFill_data_%s_%s.root",
		       cera,cv),"READ");
  assert(f && !f->IsZombie());
  //TFile *fm = new TFile(Form("rootfiles/GamHistosFill_mc_%sP8_%s.root", //needed to hardcode for now
  TFile *fm = new TFile(Form("rootfiles/GamHistosFill_mc_winter%sP8_pu-%s_w41.root",
			     ceram,cera),"READ");

  if (!fm || fm->IsZombie())
    fm = new TFile(Form("files/GamHistosFill_mc_%sP8_%s.root",
			ceram,cv),"READ");
  assert(fm && !fm->IsZombie());
  
  // Get directories
  f->cd("control");
  TDirectory *d = gDirectory;
  fm->cd("control");
  TDirectory *dm = gDirectory;
  curdir->cd();

  // Load input profiles
  TProfile *pr1b = (TProfile*)d->Get("prbal1"); clean(pr1b);
  TProfile *pr1m = (TProfile*)d->Get("prmpf1"); clean(pr1m);
  TProfile *pr6b = (TProfile*)d->Get("prbal6"); clean(pr6b);
  TProfile *pr6m = (TProfile*)d->Get("prmpf6"); clean(pr6m);
  TProfile *pr12b = (TProfile*)d->Get("prbal12"); clean(pr12b);
  TProfile *pr12m = (TProfile*)d->Get("prmpf12"); clean(pr12m);
  //
  TProfile *pr1bm = (TProfile*)dm->Get("prbal1"); clean(pr1bm);
  TProfile *pr1mm = (TProfile*)dm->Get("prmpf1"); clean(pr1mm);
  TProfile *pr6bm = (TProfile*)dm->Get("prbal6"); clean(pr6bm);
  TProfile *pr6mm = (TProfile*)dm->Get("prmpf6"); clean(pr6mm);
  TProfile *pr12bm = (TProfile*)dm->Get("prbal12"); clean(pr12bm);
  TProfile *pr12mm = (TProfile*)dm->Get("prmpf12"); clean(pr12mm);

  // Scale BAL results
  double kbal = 1.13;
  pr1b->Scale(kbal);
  pr6b->Scale(kbal);
  pr12b->Scale(kbal);
  bool isRun2 = (era=="Run2");
  double kbalm = (isRun2 ? 1.13 : 1.10);
  pr1bm->Scale(kbalm);
  pr6bm->Scale(kbalm);
  pr12bm->Scale(kbalm);

  // Calculate gain ratios
  TH1D *hr6b = pr6b->ProjectionX("hr6b"); hr6b->Divide(pr12b);
  TH1D *hr1b = pr1b->ProjectionX("hr1b"); hr1b->Divide(pr6b);
  TH1D *hr6m = pr6m->ProjectionX("hr6m"); hr6m->Divide(pr12m);
  TH1D *hr1m = pr1m->ProjectionX("hr1m"); hr1m->Divide(pr6m);
  TH1D *hr6bm = pr6bm->ProjectionX("hr6bm"); hr6bm->Divide(pr12bm);
  TH1D *hr1bm = pr1bm->ProjectionX("hr1bm"); hr1bm->Divide(pr6bm);
  TH1D *hr6mm = pr6mm->ProjectionX("hr6mm"); hr6mm->Divide(pr12mm);
  TH1D *hr1mm = pr1mm->ProjectionX("hr1mm"); hr1mm->Divide(pr6mm);

  // Fit ratio
  TF1 *f1m = new TF1("f1m","[0]",230,700);
  hr1m->Fit(f1m,"QRN");
  TF1 *f1b = new TF1("f1b","[0]",230,700);
  hr1b->Fit(f1b,"QRN");
  TF1 *f6m = new TF1("f6m","[0]",110,500);
  hr6m->Fit(f6m,"QRN");
  TF1 *f6b = new TF1("f6b","[0]",110,500);
  hr6b->Fit(f6b,"QRN");

  // Setup canvas
  TH1D *h = tdrHist("h","Response",0.945,1.095);
  TH1D *hd = tdrHist("hd","Gain ratio",0.97,1.03);
  lumi_136TeV = Form("Photon+jet, %s, %s",cera,cv);
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas("c1",h,hd,8,11);
  TLine *l = new TLine();

  // Start drawing
  c1->cd(1);
  gPad->SetLogx();
  l->DrawLine(15,1,3500,1);

  tdrDraw(pr6mm,"HIST][",kNone,kMagenta+2,kSolid,-1,kNone,0);
  tdrDraw(pr12mm,"HIST][",kNone,kBlue,kSolid,-1,kNone,0);
  tdrDraw(pr1mm,"HIST][",kNone,kRed,kSolid,-1,kNone,0);

  tdrDraw(pr6bm,"HIST][",kNone,kMagenta+2,kDashed,-1,kNone,0);
  tdrDraw(pr12bm,"HIST][",kNone,kBlue,kDashed,-1,kNone,0);
  tdrDraw(pr1bm,"HIST][",kNone,kRed,kDashed,-1,kNone,0);

  tdrDraw(pr6m,"Pz",kFullSquare,kMagenta+2,kSolid,-1,kNone,0);
  tdrDraw(pr12m,"Pz",kFullDiamond,kBlue,kSolid,-1,kNone,0);
  tdrDraw(pr1m,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0);

  tdrDraw(pr6b,"Pz",kOpenSquare,kMagenta+2,kSolid,-1,kNone,0);
  tdrDraw(pr12b,"Pz",kOpenDiamond,kBlue,kSolid,-1,kNone,0);
  tdrDraw(pr1b,"Pz",kOpenCircle,kRed,kSolid,-1,kNone,0);

  c1->cd(2);
  gPad->SetLogx();
  l->DrawLine(15,1,3500,1);

  tdrDraw(hr6mm,"HIST][",kNone,kBlue,kSolid,-1,kNone,0);
  tdrDraw(hr1mm,"HIST][",kNone,kRed,kSolid,-1,kNone,0);
  tdrDraw(hr6bm,"HIST][",kNone,kBlue,kDashed,-1,kNone,0);
  tdrDraw(hr1bm,"HIST][",kNone,kRed,kDashed,-1,kNone,0);

  tdrDraw(hr6m,"Pz",kFullSquare,kBlue);
  tdrDraw(hr1m,"Pz",kFullCircle,kRed);
  tdrDraw(hr6b,"Pz",kOpenSquare,kBlue);
  tdrDraw(hr1b,"Pz",kOpenCircle,kRed);

  f6b->SetLineColor(kBlue);
  f6b->SetLineStyle(kDashed);
  f6b->Draw("SAME");
  f6m->SetLineColor(kBlue);
  f6m->Draw("SAME");

  f1b->SetLineColor(kRed);
  f1b->SetLineStyle(kDashed);
  f1b->Draw("SAME");
  f1m->SetLineColor(kRed);
  f1m->Draw("SAME");

  // Add legend
  c1->cd(1);
  TLatex *t = new TLatex();
  t->SetTextSize(0.045); t->SetNDC();
  t->DrawLatex(0.34,0.845,"|#eta_{jet}|<1.3, |#eta_{#gamma}|<1.3");
  t->DrawLatex(0.34,0.785,"TightID, r_{9}>0.94");
  t->DrawLatex(0.34,0.725,"AK4PFPUPPI");
  t->DrawLatex(0.18,0.08,Form("BAL #times %1.3f in MC",kbalm));
  t->DrawLatex(0.18,0.03,Form("BAL #times %1.3f in DATA",kbal));
  TLegend *leg1m = tdrLeg(0.65,0.63,0.85,0.88);
  leg1m->SetHeader("MPF");
  leg1m->AddEntry(pr12m," ","PLE");
  leg1m->AddEntry(pr6m," ","PLE");
  leg1m->AddEntry(pr1m," ","PLE");
  leg1m->AddEntry(pr1mm," ","L");
  TLegend *leg1b = tdrLeg(0.725,0.63,0.925,0.88);
  leg1b->SetHeader(" BAL");
  leg1b->AddEntry(pr12b,"  Gain12","PLE");
  leg1b->AddEntry(pr6b,"  Gain6","PLE");
  leg1b->AddEntry(pr1b,"  Gain1","PLE");
  leg1b->AddEntry(pr1bm,"  MC","L");

  c1->cd(2);
  TLegend *leg2m = tdrLeg(0.18,0.50,0.38,0.90);
  leg2m->SetTextSize(0.045*2.0);
  leg2m->SetHeader("MPF");
  leg2m->AddEntry(hr6m," ","PLE");
  leg2m->AddEntry(hr1m," ","PLE");
  leg2m->AddEntry(hr1mm," ","L");
  TLegend *leg2b = tdrLeg(0.255,0.50,0.455,0.90);
  leg2b->SetTextSize(0.045*2.0);
  leg2b->SetHeader(" BAL");
  leg2b->AddEntry(hr6b,"  Gain 6 / 12","PLE");
  leg2b->AddEntry(hr1b,"  Gain 1 / 6","PLE");
  leg2b->AddEntry(hr1bm,"  MC","L");

  t->SetTextSize(0.045*1.5);
  t->SetTextColor(kBlue);
  t->DrawLatex(0.65,0.80,Form("R_{6/12} = %1.2f#pm%1.2f%% (DB)",
			      100.*(f6b->GetParameter(0)-1),
			      100.*f6b->GetParError(0)));
  t->DrawLatex(0.65,0.72,Form("        = %1.2f#pm%1.2f%% (MPF)",
			      100.*(f6m->GetParameter(0)-1),
			      100.*f6m->GetParError(0)));
  t->SetTextColor(kRed);
  t->DrawLatex(0.725,0.48,Form("R_{1/6} = %1.2f#pm%1.2f%%",
			       100.*(f1b->GetParameter(0)-1),
			       100.*f1b->GetParError(0)));
  t->DrawLatex(0.725,0.40,Form("       = %1.2f#pm%1.2f%%",
			       100.*(f1m->GetParameter(0)-1),
			       100.*f1m->GetParError(0)));

  //c1->SaveAs(Form("pdf/drawGainVsPt_%s_%s.pdf",cera,cv));
  c1->SaveAs(Form("pdf/drawGainVsPt_%s_w41%s.pdf",cera,cv));


} // drawGainsVsPt
