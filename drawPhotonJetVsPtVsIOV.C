// Purpose: draw stability of photon+jet results vs pT over time
//          compare data and MC for response, HOE, R9 etc.
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TLine.h"

#include "tdrstyle_mod22.C"

bool addMPFu2n = true;
bool addG1toMPF = false;//true;
bool addG12toMPF = false;
string id = "w15"; //change this always to current version (need to use w10 for comparison with the w10 files... change 2024B name manually to w10 although it is w11
//string id = "wX23"; //testing summer23 corrections with single files
//string id = "wX22full"; //testing summer22 corrections with all files
//string id = "wX22full-data_w5-mc"; //testing summer22 corrections with data files and summer23 corrections on mc files


bool drawFullIOVList = false;//true;

// Forward declaration of call
void drawPhotonJetVsPtVsIOVs(string so, string var, string name,
			     double y1, double y2, double z1, double z2);

//void drawPhotonJetVsPtVsIOV(string so = "resp_MPFchs_%s_a100_eta00_13") {
//void drawPhotonJetVsPtVsIOVs(string so = "control/phoevspt",
//			     string var = "H/E",
//			     double y1 = 0, double y2 = 0.01,
//			     double z1 = 0.7, double z2 = 1.4) {
			    
// Multiple calls to draw function
void drawPhotonJetVsPtVsIOV() {
  //drawPhotonJetVsPtVsIOVs("control/pphoj","Leakage","Leakage",0.,0.125,-0.005,0.03);//0.5,2.5); // before leakage cut
  //drawPhotonJetVsPtVsIOVs("control/pphoj","Leakage","Leakage",0.,0.125,-0.005,0.03);//0.5,2.5); // before leakage cut
  //drawPhotonJetVsPtVsIOVs("control/pphoj","Leakage","Leakage",0.005,0.020,-0.0024,0.004);// after leakage cut
  //drawPhotonJetVsPtVsIOVs("control/pphoj","Leakage","Leakage",0.00,0.025,-0.004,0.004);// before gain1
  drawPhotonJetVsPtVsIOVs("control/pphoj0","Leakage0","Leakage0",0.,0.025,-0.004,0.004);// after gain1, Run2 style
  //drawPhotonJetVsPtVsIOVs("control/pphoj0","Leakage0","Leakage0",0.,0.15,-0.07,0.07);// Run3
  //drawPhotonJetVsPtVsIOVs("control/pcorrvspt","Correction","Correction",0.99,1.01,-0.001,0.001);//0.5,2.5);
  //drawPhotonJetVsPtVsIOVs("control/phoevspt","H/E","HoE",0.,0.01,0.7,1.6);
  //drawPhotonJetVsPtVsIOVs("control/pr9vspt","R9","R9",0.83,1.03,0.9,1.14);

  //drawPhotonJetVsPtVsIOVs("control/pgain1vspt","Gain1","Gain1",0.,1.2,0.5,1.5);
  //drawPhotonJetVsPtVsIOVs("control/pgain6vspt","Gain6","Gain6",0.,1.2,0.5,1.5);
  //drawPhotonJetVsPtVsIOVs("control/pgain12vspt","Gain12","Gain12",0.,1.2,0.5,1.5);

  //drawPhotonJetVsPtVsIOVs("resp_MPFchs_%s_a100_eta00_13","MPF","MPF",0.88,1.23,0.88,1.08);
  drawPhotonJetVsPtVsIOVs("resp_MPFchs_%s_a100_eta00_13","MPF","MPF",0.90,1.15,0.89,1.09);//1.06); // Run3
  drawPhotonJetVsPtVsIOVs("resp_DBchs_%s_a100_eta00_13","DB","DB",0.75,1.20,0.89,1.09);//0.80,1.10); // Run3
  //drawPhotonJetVsPtVsIOVs("resp_DBchs_%s_a100_eta00_13","DB","DB",0.75,1.20,0.89,1.09);

  drawPhotonJetVsPtVsIOVs("resp_RES_%s_a100_eta00_13","RES","Residual JES",0.8,1.1,-0.1,0.1);
  
  //drawPhotonJetVsPtVsIOVs("resp_MPFchs_%s_a100_eta00_13","MPF","MPF",0.94,1.08,0.95,1.03);
  //drawPhotonJetVsPtVsIOVs("resp_MPFR1chs_%s_a100_eta00_13","MPF1","MPF1",0.8,1.10,0.95,1.05);
  //drawPhotonJetVsPtVsIOVs("resp_MPFRnchs_%s_a100_eta00_13","MPFn","MPFn",-0.1,0.2,-0.03,0.03);
  //drawPhotonJetVsPtVsIOVs("resp_MpfRuchs_%s_a100_eta00_13","MPFu","MPFu",-0.05,0.15,-0.03,0.02);
} // drawPhotonJetVsPtVsIOV

void drawPhotonJetVsPtVsIOVs(string so, string var, string name,
			     double y1, double y2, double z1, double z2) {
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //string iovs[] = {"2016BCDEF","2016FGH","2017BCDEF","2018ABCD","Run2"};
  //string mcs[] = {"2016APVP8","2016P8","2017P8","2018P8","Run2P8"};
  //string iovs[] = {"2016BCDEF","2016FGH","2017BCDEF","2018ABCD"};
  //string mcs[] = {"2016APVP8","2016P8","2017P8","2018P8"};
  string iovs_long[] = {
    "2022C","2022D","2022E","2022F","2022G",
    "2023Cv123","2023Cv4","2023D","2024B", "2024C"
  };
  string iovs_short[] = {
    //"2018ABCD","Run3",
    //"2022CD","2022E","2022FG", //hadd
    //"2023Cv123","2023Cv4D", //hadd
    "2023Cv123","2023Cv4","2023D","2024B", "2024C" //hadd <-- change back to this after testing single files
    //"2023Cv123X","2023Cv4X","2023DX" //testing single files with different jec
  };

  string mcs_long[] = {
    "2022P8","2022P8","2022EEP8","2022EEP8","2022EEP8", //orig
    "2022P8","2022P8","2022P8"
  };
  string mcs_short[] = {
    //"2018P8","Run3P8",
    //"2022P8","2022EEP8","2022EEP8", //hadd
    //"2022P8","2022P8"
    //"2022P8","2022P8","2022P8"
    //"2023P8QCD","2023P8QCD","2023-BPixP8QCD" //in principle use QCD
    "2023P8","2023P8","2023P8-BPix", "2023P8-BPix", "2023P8-BPix" //<-- change back to this after testing single files; use 2023P8-BPix also for 2024B and C
    //"2023P8X","2023P8X","2023P8-BPixX"//testing single files with different jec
  };
  const int niov_long = sizeof(iovs_long)/sizeof(iovs_long[0]);
  const int nmc_long = sizeof(mcs_long)/sizeof(mcs_long[0]);
  const int niov_short = sizeof(iovs_short)/sizeof(iovs_short[0]);
  const int nmc_short = sizeof(mcs_short)/sizeof(mcs_short[0]);

  string *iovs = (drawFullIOVList ? &iovs_long[0] : &iovs_short[0]);
  string *mcs = (drawFullIOVList ? &mcs_long[0] : &mcs_short[0]);
  const int niov = (drawFullIOVList ? niov_long : niov_short);
  const int nmc = (drawFullIOVList ? nmc_long : nmc_short);

  assert(niov==nmc);
  
  map<string,int> mcolor;
  mcolor["2016BCDEF"] = kBlue;
  mcolor["2016FGH"] = kCyan+2;
  mcolor["2017BCDEF"] = kGreen+2;
  mcolor["2018ABCD"] = kGray+2;//kRed;
  mcolor["Run2"] = kBlack;
  //
  mcolor["2022C"] = kGreen+2;//kBlue;
  mcolor["2022D"] = kGreen+3;//kCyan+1;
  mcolor["2022E"] = kCyan+2;
  mcolor["2022CD"] = kGreen+2;//kBlue;
  mcolor["2022CDE"] = kGreen+2;//kBlue;
  mcolor["2022F"] = kRed;
  mcolor["2022G"] = kRed+1;//kOrange+2;
  mcolor["2022FG"] = kRed;
  mcolor["2023Cv123"] = kOrange+1;//kYellow+2;
  mcolor["2023Cv4"] = kBlue;//kGreen+2;
  mcolor["2023D"] = kMagenta;//+2;
  mcolor["2023Cv4D"] = kBlue;//kGreen+2;
  mcolor["2024B"] = kRed;
  mcolor["2024C"] = kGreen;
  mcolor["Run3"] = kBlack;

  //for my investigations on the mpf issue
  mcolor["2023Cv123X"] = kOrange+1;//kYellow+2;
  mcolor["2023Cv4X"] = kBlue;//kGreen+2;
  mcolor["2023DX"] = kGreen+2;//+2;
  mcolor["2023D"] = kGreen+2;//+2; --> overwriting old setting (keep only for wX22full version)
  mcolor["2024B"] = kRed;
  mcolor["2024C"] = kGreen;

 

  map<string,int> mmarker;
  mmarker["2022C"] = kFullSquare;
  mmarker["2022D"] = kOpenSquare;
  mmarker["2022E"] = kOpenSquare;
  mmarker["2022CD"] = kFullSquare;
  mmarker["2022CDE"] = kFullSquare;
  mmarker["2022F"] = kFullTriangleUp;
  mmarker["2022G"] = kOpenTriangleUp;
  mmarker["2022FG"] = kFullTriangleUp;
  mmarker["2023Cv123"] = kFullDiamond;//kFullCircle;
  mmarker["2023BCv123"] = kFullDiamond;//kFullCircle;
  mmarker["2023Cv4"] = kFullTriangleDown;//kFullDiamond;
  mmarker["2023D"] = kOpenTriangleDown;//kOpenDiamond;
  mmarker["2023Cv4D"] = kFullTriangleDown;//kFullDiamond;
  mmarker["2024B-PromptReco-v1"] = kFullCircle;
  mmarker["Run3"] = kFullSquare;

  //for my investigation on the mpf issue
  mmarker["2023Cv123X"] = kFullDiamond;//kFullCircle;
  mmarker["2023Cv4X"] = kFullTriangleDown;//kFullDiamond;
  mmarker["2023DX"] = kOpenTriangleDown;//kOpenDiamond;
 
  
  const char *cvar = var.c_str();
  const char *cname = name.c_str();

  TH1D *h = tdrHist("h",cvar,y1,y2);
  TH1D *h2 = tdrHist("h2","Data/MC",z1,z2);
  //lumi_13TeV = "Run2 v16";
  //lumi_13TeV = "Run3 v21";
  lumi_13TeV = "Run2"; // 4=13 TeV
  if (id!="") lumi_136TeV = Form("Run3 %s",id.c_str()); // 8=13.6 TeV
  TCanvas *c1 = tdrDiCanvas(Form("c1_%s",cname),h,h2,8,11);


  c1->cd(1);
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(h->GetXaxis()->GetXmin(),1,h->GetXaxis()->GetXmax(),1);
  l->DrawLine(h->GetXaxis()->GetXmin(),0,h->GetXaxis()->GetXmax(),0);
  gPad->SetLogx();

  //TLegend *leg = tdrLeg(0.60,0.90-niov*0.06,0.80,0.90);
  TLegend *leg = tdrLeg(0.60,0.90-niov*0.045,0.80,0.90);

  c1->cd(2);
  gPad->SetLogx();
  l->DrawLine(h->GetXaxis()->GetXmin(),1,h->GetXaxis()->GetXmax(),1);
  l->DrawLine(h->GetXaxis()->GetXmin(),0,h->GetXaxis()->GetXmax(),0);

  for (int i = 0; i != niov; ++i) {

    string iov = iovs[i];
    const char *ciov = iov.c_str();
    const char *cmc = mcs[i].c_str();
    const char *cid = id.c_str();

    TFile *fd(0), *fm(0);
    /*
    if (iovs[i]=="2018ABCD") {
    //fd = new TFile(Form("files/GamHistosFill_data_%s_v16.root",ciov));
      fd = new TFile(Form("files/GamHistosFill_data_%s_v20.root",ciov));
    }
    else {
      //fd = new TFile(Form("files/GamHistosFill_data_%s_v21.root",ciov));
      fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",ciov,cid));
    }
    assert(fd && !fd->IsZombie());
    if (iovs[i]=="2018ABCD") {
      //fm = new TFile(Form("files/GamHistosFill_mc_%s_v16.root",cmc));      
      fm = new TFile(Form("files/GamHistosFill_mc_%s_v20.root",cmc));
    }
    else {
      //fm = new TFile(Form("files/GamHistosFill_mc_%s_v21.root",cmc));
      fm = new TFile(Form("rootfiles/GamHistosFill_mc_%s_%s.root",cmc,cid));
    }
    assert(fm && !fm->IsZombie());
    */

    //add this just for short time to check w2-w3 samples (w3 needed for BPix, otherwise same)
    /*
    if (iovs[i]=="2023D") { //bpix sample, needs w3+
      fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",ciov,cid));
      fm = new TFile(Form("rootfiles/GamHistosFill_mc_%s_%s.root",cmc,cid));
    }
    else { //can use w2 samples, no bpix issue here
      fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_w2.root",ciov));
      fm = new TFile(Form("rootfiles/GamHistosFill_mc_%s_w2.root",cmc));
    }
    */

    //fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_wX22full.root",ciov)); //only to test combination of 2022-jec for data and 2023-jec for MC
    //fm = new TFile(Form("rootfiles/GamHistosFill_mc_%s_w5.root",cmc)); //only to test combination of 2022-jec for data and 2023-jec for MC
    fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",ciov,cid));
    fm = new TFile(Form("rootfiles/GamHistosFill_mc_%s_%s.root",cmc,cid));
 
    assert(fd && !fd->IsZombie());
    assert(fm && !fm->IsZombie());


    curdir->cd();
    
    TObject *od = fd->Get(Form(so.c_str(),"DATA")); assert(od);
    TObject *om = fm->Get(Form(so.c_str(),"MC")); assert(om);
    
    TH1D *hd(0), *hm(0);
    if (od->InheritsFrom("TProfile")) {
      hd = ((TProfile*)od)->ProjectionX(Form("hd_%s_%s",cvar,ciov));
      hm = ((TProfile*)om)->ProjectionX(Form("hm_%s_%s",cvar,ciov));

      if (name=="MPFn" && addMPFu2n) {
	const char *co = "resp_MpfRuchs_%s_a100_eta00_13";
	TProfile *pd = (TProfile*)fd->Get(Form(co,"DATA")); assert(pd);
	TProfile *pm = (TProfile*)fm->Get(Form(co,"MC")); assert(pm);
	hd->Add(pd,1.2);
	hm->Add(pm,1.2);
      }
      if ((name=="MPF" || name=="MPF1") && addG1toMPF) {
	const char *co = "control/pgain1vspt";
	TProfile *pd = (TProfile*)fd->Get(Form(co,"DATA")); assert(pd);
	TProfile *pm = (TProfile*)fm->Get(Form(co,"MC")); assert(pm);
	hd->Add(pd,0.02);
      }
      if ((name=="MPF" || name=="MPF1") && addG12toMPF) {
	const char *co = "control/pgain12vspt";
	TProfile *pd = (TProfile*)fd->Get(Form(co,"DATA")); assert(pd);
	TProfile *pm = (TProfile*)fm->Get(Form(co,"MC")); assert(pm);
	hd->Add(pd,0.0);
      }
    }
    else {
      hd = (TH1D*)od;
      hm = (TH1D*)om;
    }
    assert(hd);
    assert(hm);
    
    TH1D *hr = (TH1D*)hd->Clone(Form("hr_%s_%s",cvar,ciov));
    if (name=="MPFn" || name=="MPFu" || name=="Leakage" || name=="Leakage0") {
      hr->Add(hm,-1);
      h2->SetYTitle("Data-MC");
    }
    else
      hr->Divide(hm);
    
    c1->cd(1);
    gPad->SetLogx();
    tdrDraw(hm,"H",kNone,(mcolor[iov] ? mcolor[iov] : kBlack),kSolid,-1,kNone);
    tdrDraw(hd,"Pz",(mmarker[iov] ? mmarker[iov] : kFullCircle),
	    (mcolor[iov] ? mcolor[iov] : kBlack));
    
    leg->AddEntry(hd,ciov,"PLE");

    c1->cd(2);
    gPad->SetLogx();
    tdrDraw(hr,"Pz",(mmarker[iov] ? mmarker[iov] : kFullCircle),
	    (mcolor[iov] ? mcolor[iov] : kBlack));
  } // for iov

  if (id!="")
      c1->SaveAs(Form("pdf/drawPhotonJetVsPtVsIOVs_%s_%s.pdf",
		      name.c_str(),id.c_str()));
  else
    c1->SaveAs(Form("pdf/drawPhotonJetVsPtVsIOVs_%s.pdf",name.c_str()));
} // void drawPhotonJetVsPtVsIOVs
