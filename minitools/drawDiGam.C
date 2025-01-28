// Purpose: Draw DiGamHistos.Fill:/DiGam/hmgg* to find eta_t->gamma+gamma
// Date: 4 Oct 2024
// Author: miquork
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"

#include "../tdrstyle_mod22.C"

void drawDiGam() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/DiGamHistosFill_data_2024FG.root","READ");
  TFile *f = new TFile("rootfiles/DiGamHistosFill_data_2024F_vX.root","READ");
  //TFile *f = new TFile("rootfiles/DiGamHistosFill_data_2024G_vX.root","READ");
  assert(f && !f->IsZombie());

  //2024F:   29.558801123      | 27.757490090 
  //2024G:  40.399886502      | 37.772631997 
      //double lumi = (27.757+37.773)*1000; // 2024F+G, 1/fb->1/pb
  //double lumi = 37.773*1000; // 2024G, 1/fb->1/pb
  double lumi = 27.757*1000; // 2024G, 1/fb->1/pba
  //double lumi = (20+30)*1000; // 2024F+G, 1/fb->1/pb

  TH1D *hmgg = (TH1D*)f->Get("DiGam/hmgg1_13"); assert(hmgg);
  //TH1D *hmgg = (TH1D*)f->Get("DiGam/hmgg1"); assert(hmgg);
  TH1D *hmggn = (TH1D*)hmgg->Clone("hmggn");
  hmggn->Rebin(2.);
  hmggn->Scale(1./lumi,"width");

  double xmin = 140;//110;//120;
  double xmax = 1100;//2000;
  //TF1 *f1 = new TF1("f1","[0]*pow(x,[1]+[2]*log(x))",120,2000);
  //TF1 *f1 = new TF1("f1","[0]*pow(x,[1]+[2]*log(0.01*x))",120,2000);
  TF1 *f1 = new TF1("f1","[0]*pow(x,[1]+[2]*log(0.01*x))",xmin,xmax);
  //f1->SetParameters(1e5,-5,0);
  //f1->SetParameters(5.2294e1,-1.57839,-3.29205e-1);
  f1->SetParameters(5.22675e1,-1.57829,-3.29220e-1);
  //hmggn->Fit(f1,"QRN");
  hmggn->Fit(f1,"RN");

  double m_eta_t = 2*172.5; // GeV
  double gamma_total = 7.06;//2*1.35;//7.06; // GeV;
  double emres = 0.02;
  double resol = emres*m_eta_t;
  //double resol = 0.01*m_eta_t;
  //double resol = 0.03*m_eta_t; // vertex smearing?
  double br_gg = 1.33e-4;
  double width_eta_t = sqrt(pow(gamma_total,2) + pow(resol,2));
  double E_b = 1.53; // or 0.575? or 0.95 GeV?
  //TF1 *f2 = new TF1("f2","[0]*pow(x,[1]+[2]*log(x))+fabs([3])*TMath::Gaus(x,[4],[5],1)",120,2000);
  //TF1 *f2 = new TF1("f2","[0]*pow(x,[1]+[2]*log(0.01*x))+fabs([3])*TMath::Gaus(x,[4],[5],1)",120,2000);
  TF1 *f2 = new TF1("f2","[0]*pow(x,[1]+[2]*log(0.01*x))+fabs([3])*TMath::Gaus(x,[4]-[6],[5],1)",xmin,xmax);
  f2->SetParameters(f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2),0.,m_eta_t,width_eta_t,E_b);
  //f2->SetParLimits(3,0,1e10);
  f2->FixParameter(4,m_eta_t);
  f2->FixParameter(5,width_eta_t);
  f2->FixParameter(6,E_b);
  //hmggn->Fit(f2,"QRN");
  hmggn->Fit(f2,"RN");

  TH1D *hmggd = (TH1D*)hmggn->Clone("hmggd");
  TH1D *hmggd2 = (TH1D*)hmggn->Clone("hmggd2");
  TH1D *hmggr = (TH1D*)hmggn->Clone("hmggr");
  TH1D *hmggr2 = (TH1D*)hmggn->Clone("hmggr2");
  double chi2(0), chi22(0);
  for (int i = 1; i != hmggd->GetNbinsX()+1; ++i) {

    double pt = hmggd->GetBinCenter(i);
    double data = hmggd->GetBinContent(i); 
    double sigma = hmggd->GetBinError(i);
    double fit = f1->Eval(pt);
    double fit2 = f2->Eval(pt);

    if (sigma!=0) {

      hmggd->SetBinContent(i, (data-fit)/sigma);
      hmggd->SetBinError(i, 0.);
      hmggd2->SetBinContent(i, (data-fit2)/sigma);
      hmggd2->SetBinError(i, 0.);

      if (pt>200. && pt<1000.) {
	chi2 += (data-fit)/sigma;
	chi22 += (data-fit2)/sigma;
      }
      
      hmggr->SetBinContent(i, (data-fit)*20. * lumi);
      hmggr->SetBinError(i, 0.);
      hmggr2->SetBinContent(i, (data-fit2)*20. * lumi);
      hmggr2->SetBinError(i, 0.);
    }
  }
  
  //TH1D *h = tdrHist("h","N_{event} / GeV",0.02,1e6,"m_{#gamma#gamma} (GeV)",50,2400);
  TH1D *h = tdrHist("h","Cross section (pb / GeV)",2e-7,1e2,"m_{#gamma#gamma} (GeV)",50,2400);
  TH1D *hd = tdrHist("hd","(Data-fit)/#sigma",-5,5,"m_{#gamma#gamma} (GeV)",50,2400);
  lumi_136TeV = Form("2024FG, %1.1f fb^{-1}",lumi/1000.);
  //lumi_136TeV = Form("2024G, %1.1f fb^{-1}",lumi/1000.);
  //lumi_136TeV = Form("2024F, %1.1f fb^{-1}",lumi/1000.);
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas("c1",h,hd,8,11);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.035);

  
  c1->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();

  tdrDraw(hmggn,"HE",kNone,kBlack,kSolid,-1,1001,kBlue-9);

  f2->SetLineColor(kRed);
  f2->Draw("SAME");
  f1->SetLineColor(kBlue+2);
  f1->Draw("SAME");

  gPad->RedrawAxis();

  tex->DrawLatex(0.40,0.85,"Tight cut-based photon ID, R_{9}>0.94");
  //tex->DrawLatex(0.40,0.81,"|#eta_{#gamma,1}|<1.3 or |#eta_{#gamma,2}|<1.3 and p_{T,#gamma}>50 GeV");
  tex->DrawLatex(0.40,0.81,"|#eta_{#gamma,1,2}|<1.3 and p_{T,#gamma,1}>50 GeV");

  tex->DrawLatex(0.40,0.77,"f(m_{#gamma#gamma})=N#timespow(m_{#gamma#gamma},#alpha+#betalog(m_{#gamma#gamma}))+S#timesgaus");
  tex->DrawLatex(0.45,0.73,Form("Fit Range %1.0f<m_{#gamma#gamma}<%1.0f GeV",xmin,xmax));
  tex->DrawLatex(0.45,0.69,"#chi^{2} for 200<m_{#gamma#gamma}<1000 GeV");
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.45,0.64,Form("#chi^{2}(B) = %1.1f / %d",
				f1->GetChisquare(),f1->GetNDF()));
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.45,0.60,Form("#chi^{2}(S+B) = %1.1f / %d",
				f2->GetChisquare(),f2->GetNDF()));
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.45,0.56,Form("m(#eta_{t}) #equiv 2m_{t}-E_{b} = (%1.1f-%1.1f) GeV",m_eta_t,E_b));
  //tex->DrawLatex(0.50,0.52,Form("#sigma(#eta_{t}) #equiv 7.06 GeV #oplus 2%%#timesm(#eta_{t})"));
  //tex->DrawLatex(0.50,0.48,Form("          = %1.1f GeV",width_eta_t));
  tex->DrawLatex(0.45,0.52,Form("#sigma(#eta_{t}) #equiv %1.1f GeV#oplus%1.1f%%#timesm(#eta_{t}) = %1.1f GeV",gamma_total,emres*100,width_eta_t));
  //tex->DrawLatex(0.45,0.52,Form("#sigma(#eta_{t}) #equiv 7.06 GeV #oplus 1%%#timesm(#eta_{t}) = %1.1f GeV",width_eta_t));
  //tex->DrawLatex(0.45,0.52,Form("#sigma(#eta_{t}) #equiv 7.06 GeV #oplus 3%%#timesm(#eta_{t}) = %1.1f GeV",width_eta_t));
  tex->DrawLatex(0.45,0.48,Form("xsec(#eta_{t}) #equiv fidX / BR(#eta_{t}#rightarrow#gamma#gamma)"));
  //tex->DrawLatex(0.60,0.30,Form("= N/1.33#times10^{-4} = %1.1f pb",fabs(f2->GetParameter(3))/1.33e-4));
  tex->DrawLatex(0.56,0.44,Form("= S/1.33#times10^{-4} = %1.1f#pm%1.1f pb",fabs(f2->GetParameter(3))/1.33e-4,f2->GetParError(3)/1.33e-4));

  
  c1->cd(2);
  gPad->SetLogx();

  tdrDraw(hmggd,"HISTE",kNone,kBlack,kSolid,-1,1001,kBlue-9);
  tdrDraw(hmggd2,"HISTE",kNone,kRed,kSolid,-1,kNone,kBlue-9);

  gPad->RedrawAxis();

  
  c1->SaveAs("pdf/drawDiGam/drawDiGam.pdf");

  cout << "Observed " << f2->GetParameter(3) * lumi << " +/- "
       << f2->GetParError(3) * lumi << " signal events in "
       << lumi/1000. << "/fb" << endl;
  double B = f1->Integral(m_eta_t-width_eta_t,m_eta_t+width_eta_t)*lumi;
  cout << "Observed " << B << " background events" << endl;
  
  TH1D *h_2 = tdrHist("h_2","(Data - fit) / 20 GeV",-1200,2300,"m_{#gamma#gamma} (GeV)",200,1000);
  //TH1D *h_2d = tdrHist("h_2d","(Data - B+S fit) / 20 GeV",-1200,2300,"m_{#gamma#gamma} (GeV)",200,1000);
  lumi_136TeV = Form("2024FG, %1.1f fb^{-1}",lumi/1000.);
  //lumi_136TeV = Form("2024G, %1.1f fb^{-1}",lumi/1000.);
  //lumi_136TeV = Form("2024F, %1.1f fb^{-1}",lumi/1000.);
  extraText = "Private";
  //TCanvas *c2 = tdrDiCanvas("c2",h_2,h_2d,8,11);
  TCanvas *c2 = tdrCanvas("c2",h_2,8,11,kSquare);

  //c2->cd(1);

  //hmggr->Draw();
  tdrDraw(hmggr,"HISTE",kNone,kBlack,kSolid,-1,1001,kBlue-9);
  tdrDraw(hmggr2,"HISTE",kNone,kRed,kSolid,-1,kNone,kBlue-9);

  TLegend *leg = tdrLeg(0.5,0.85-2*0.05,0.85,0.85);
  leg->AddEntry(hmggr,"B-only fit","F");
  leg->AddEntry(hmggr2,"S+B fit","F");
  
  //c2->cd(2);

  //tdrDraw(hmggr2,"HISTE",kNone,kRed,kSolid,-1,kNone,kBlue-9);

  c2->SaveAs("pdf/drawDiGam/drawDiGam_diff.pdf");


  // Ravindra's interference in box diagram:
  // https://arxiv.org/abs/1605.07360
  // Original eta_t discovery paper:
  // https://cds.cern.ch/record/2911775
  // Toponium Deep Dive:
  // https://indico.cern.ch/event/1412949/
  TH1D *hmgg_13 = (TH1D*)f->Get("DiGam/hmgg_13"); assert(hmgg_13);
  TH1D *hmgg20 = (TH1D*)hmgg_13->Clone("hmgg20");
  hmgg20->Rebin(20);

  //TH1D *h_3 = tdrHist("h_3","Events / 20 GeV",20,1e5,"m_{#gamma#gamma} (GeV)",200,1000);
  TH1D *h_3 = tdrHist("h_3","Events / 20 GeV",5,2.5e4,"m_{#gamma#gamma} (GeV)",200,1000);
  TCanvas *c3 = tdrCanvas("c3",h_3,8,11,kSquare);
  gPad->SetLogy();
  tdrDraw(hmgg20,"PE",kFullCircle,kBlack);

  c3->SaveAs("pdf/drawDiGam/drawDiGam_AtlasEvents.pdf");


  TH1D *hmgg13 = (TH1D*)f->Get("DiGam/hmgg_13"); assert(hmgg13);
  TH1D *hmggn13 = (TH1D*)hmgg13->Clone("hmggn13");
  hmggn13->Rebin(10.);
  hmggn13->Scale(1./lumi,"width");
  double dpt = hmggn13->GetBinWidth(1);
  for (int i = 1; i != hmggn13->GetNbinsX()+1; ++i) {
    double dpt = hmggn13->GetBinWidth(i);
    double fit = f1->Integral(hmggn13->GetBinLowEdge(i),hmggn13->GetBinLowEdge(i+1)) / dpt;
    //double fit = f1->Eval(hmggn13->GetBinCenter(i));
    double nfit = fit * lumi;
    hmggn13->SetBinContent(i, (hmggn13->GetBinContent(i)-fit)*lumi*dpt);
    hmggn13->SetBinError(i, hmggn13->GetBinError(i)*lumi*dpt);
			 //sqrt(pow(hmggn13->GetBinError(i),2)+nfit/5.));
  } // for i

  TF1 *f3 = new TF1("f2","fabs([0])*TMath::Gaus(x,[1]-[3],[2],1)",xmin,xmax);
  f3->SetParameters(f2->GetParameter(3)*lumi*dpt,m_eta_t,width_eta_t,E_b);
  f3->FixParameter(0,f2->GetParameter(3)*lumi*dpt);
  f3->FixParameter(1,m_eta_t);
  f3->FixParameter(2,width_eta_t);
  f3->FixParameter(3,E_b);
  
  // Just focusing on the 300-400 GeV range
  TH1D *h_4 = tdrHist("h_4",Form("Bkg-substracted events / %1.0f GeV",dpt),
		      -40,60,"m_{#gamma#gamma} (GeV)",300-1e-3,400+1e-3);
  TCanvas *c4 = tdrCanvas("c4",h_4,8,11,kSquare);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(300,0,400,0);
  
  tdrDraw(hmggn13,"PE",kFullCircle,kBlack);
  f3->Draw("SAME");

  tex->SetTextSize(0.035);
  tex->DrawLatex(0.65,0.85,Form("N_{tot}(#eta_{t}) = %1.1f pb",
				fabs(f2->GetParameter(3))/1.33e-4));
  tex->DrawLatex(0.65,0.80,Form("BR(#eta_{t}#rightarrow#gamma#gamma) = %1.2f#times10^{%d}",
				1.33,-4));
  tex->DrawLatex(0.65,0.75,Form("m(#eta_{t}) = %1.1f GeV",
				m_eta_t-E_b));
  tex->DrawLatex(0.65,0.70,Form("#sigma(#eta_{t}) = %1.1f GeV",
				width_eta_t));
  
  c4->SaveAs("pdf/drawDiGam/drawDiGam_diff2.pdf");
  

}
