// Purpose: draw stability of photon+jet results vs pT over time
//          compare data and MC for response, HOE, R9 etc.
// Note: This is essentially a copy of the original drawPhotonJetVsPtVsIOV.C, adjusted to our needs while deriving new correctiongs for 2023. I removed everything not needed at this point.
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TLine.h"

#include "tdrstyle_mod22.C"

bool addMPFu2n = true;
bool addG1toMPF = false;//true;
bool addG12toMPF = false;
//string id = "w5"; //change this back to current version!!
//string id = "wX23"; //testing summer23 corrections with single files
//string id = "wX22full"; //testing summer22 corrections with all files
//string id = "wX22full-data_w5-mc"; //testing summer22 corrections with data files and summer23 corrections on mc files
//string id = "wX22full-data_w5-mc_plus-extra"; //displaying even more data and mc results.
//string id = "w7-29feb2024"; //comparing: Cv123, Cv4, D, P8, P8-BPix, P8QCD and P8QCD-BPix for 2023
//string id = "w8-09mar2024"; //comparing: Cv123, Cv4, D, P8, P8-BPix, P8QCD and P8QCD-BPix for 2023 (after bug-fix)
//string id = "w9-13mar2024"; //comparing: Cv123, Cv4, D, P8, P8-BPix, P8QCD and P8QCD-BPix for 2023 (using jetvetomaps for photons)
string id = "w17"; //comparing 2024B and 2024C
//string id = "various comparions"; //displaying even more data and mc results.



bool drawFullIOVList = false;//true;

// Forward declaration of call
void drawPhotonJetVsPtVsIOVscomparisons(string so, string var, string name,
			     double y1, double y2, double z1, double z2);

//void drawPhotonJetVsPtVsIOV(string so = "resp_MPFchs_%s_a100_eta00_13") {
//void drawPhotonJetVsPtVsIOVs(string so = "control/phoevspt",
//			     string var = "H/E",
//			     double y1 = 0, double y2 = 0.01,
//			     double z1 = 0.7, double z2 = 1.4) {
			    
// Multiple calls to draw function
void drawPhotonJetVsPtVsIOVcomparisons() {
  drawPhotonJetVsPtVsIOVscomparisons("control/pphoj0","Leakage0","Leakage0",0.,0.025,-0.004,0.004);// after gain1, Run2 style
  //drawPhotonJetVsPtVsIOVs("control/pphoj0","Leakage0","Leakage0",0.,0.15,-0.07,0.07);// Run3
  //drawPhotonJetVsPtVsIOVs("control/pcorrvspt","Correction","Correction",0.99,1.01,-0.001,0.001);//0.5,2.5);
  drawPhotonJetVsPtVsIOVscomparisons("resp_MPFchs_%s_a100_eta00_13","MPF","MPF",0.90,1.15,0.89,1.09);//1.06); // Run3
  drawPhotonJetVsPtVsIOVscomparisons("resp_DBchs_%s_a100_eta00_13","DB","DB",0.75,1.20,0.89,1.09);//0.80,1.10); // Run3
  drawPhotonJetVsPtVsIOVscomparisons("resp_RES_%s_a100_eta00_13","RES","Residual JES",0.8,1.1,-0.1,0.1);
} // drawPhotonJetVsPtVsIOVcomparisons

void drawPhotonJetVsPtVsIOVscomparisons(string so, string var, string name,
			     double y1, double y2, double z1, double z2) {
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //string iovs[] = {"2016BCDEF","2016FGH","2017BCDEF","2018ABCD","Run2"};
  //string mcs[] = {"2016APVP8","2016P8","2017P8","2018P8","Run2P8"};
  //string iovs[] = {"2016BCDEF","2016FGH","2017BCDEF","2018ABCD"};
  //string mcs[] = {"2016APVP8","2016P8","2017P8","2018P8"};
  string iovs_long[] = {
    //"2022C","2022D","2022E","2022F","2022G", "2023Cv123","2023Cv4","2023D"
    "2024B", "2024C", "2024D"
  };
  string iovs_short[] = {
    //"2018ABCD","Run3",
    //"2023Cv123","2023Cv4","2023D"//, //hadd <-- change back to this after testing single files
    "2024B", "2024C", "2024D"
    //"2018ABCD" //added UL2018, use v20 for this
  };

  string mcs_long[] = {
    "2022P8","2022P8","2022EEP8","2022EEP8","2022EEP8", //orig
    "2022P8","2022P8","2022P8"
  };
  string mcs_short[] = {
    //"2023P8QCD","2023P8QCD","2023-BPixP8QCD" //in principle use QCD
    //"2023P8","2023P8","2023P8-BPix"//, //<-- change back to this after testing single files
    "2023P8-BPix", "2023P8-BPix", "2023P8-BPix" //for 2024
    //"2018P8" //use v20 for this
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
  mcolor["2018ABCD"] = kGray+2;//kRed;
  mcolor["Run2"] = kBlack;
  //
  mcolor["2022E"] = kCyan+2;
  mcolor["2022CD"] = kGreen+2;//kBlue;
  mcolor["2022CDE"] = kGreen+2;//kBlue;
  mcolor["2022F"] = kRed;
  mcolor["2022G"] = kRed+1;//kOrange+2;
  mcolor["2023D"] = kMagenta;//+2;
  mcolor["Run3"] = kBlack;

  //for my investigations on the mpf issue
  mcolor["2023Cv123"] = kOrange+1;//kYellow+2;
  mcolor["2023Cv4"] = kBlue;//kGreen+2;
  mcolor["2023D"] = kGreen+2;//+2; --> overwriting old setting (keep only for wX22full version)
  mcolor["2018ABCD"] = kGray+2;//kRed;
  mcolor["2022P8"] = kRed;


  //NOTE: TO DO - SHOULD START HANDLING THIS AS A MAP OR DICT
 

  map<string,int> mmarker;
  mmarker["2022C"] = kFullSquare;
  mmarker["2022F"] = kFullTriangleUp;
  mmarker["2022G"] = kOpenTriangleUp;
  mmarker["2023Cv4D"] = kFullTriangleDown;//kFullDiamond;
  mmarker["Run3"] = kFullSquare;

  //for my investigation on the mpf issue
  mmarker["2023Cv123"] = kFullDiamond;//kFullCircle;
  mmarker["2023Cv4"] = kFullTriangleDown;//kFullDiamond;
  mmarker["2023D"] = kOpenTriangleDown;//kOpenDiamond;
  mmarker["2018ABCD"] = kFullSquare;

  //2024
  mmarker["2024B"] = kFullDiamond;//kFullCircle;
  mmarker["2024C"] = kFullTriangleDown;//kFullDiamond;
  mmarker["2024D"] = kOpenCircle;
  mcolor["2024B"] = kBlue;//kGreen+2;
  mcolor["2024C"] = kGreen+2;
  mcolor["2024D"] = kOrange+2;
 
 
  
  const char *cvar = var.c_str();
  const char *cname = name.c_str();

  TH1D *h = tdrHist("h",cvar,y1,y2);
  //TH1D *h2 = tdrHist("h2","Data/MC(noQCD)",z1,z2);
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
  //TLegend *leg = tdrLeg(0.60,0.90-niov*0.045,0.80,0.90);
  TLegend *leg = tdrLeg(0.60,0.90-niov*0.08,0.80,0.90);
  //TLegend *leg = tdrLeg(0.60,0.90-niov*0.06,0.80,0.90);

  leg->SetTextFont(42);




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

    //getting the correct input files (as they also differ in code version)
    if (iovs[i]=="2023Cv123" || iovs[i]=="2023Cv4" || iovs[i]=="2023D" || 
        iovs[i]=="2024B" || iovs[i]=="2024C" || iovs[i]=="2024D") {
      fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_w17.root",ciov)); //newest ones
      fm = new TFile(Form("rootfiles/GamHistosFill_mc_%s_w13.root",cmc));    //newest ones
      //fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_wX22full.root",ciov)); //newest data but with 2022 correctiongs
      //fm = new TFile(Form("rootfiles/GamHistosFill_mc_%s_wX22full.root",cmc));    //newest  data but with 2022 correctiongs
    }
/*
    if (iovs[i]=="2018ABCD") {
      fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_v20.root",ciov));
      fm = new TFile(Form("rootfiles/GamHistosFill_mc_%s_v20.root",cmc));
    }
*/

    //fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",ciov,cid));
    //fm = new TFile(Form("rootfiles/GamHistosFill_mc_%s_%s.root",cmc,cid));
 
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
    
    //to add one info on mc?
    if(iovs[i]=="2024D"){
        leg->AddEntry(hm,"2023 bpix mc","PLE");
    }

    c1->cd(2);
    gPad->SetLogx();
    tdrDraw(hr,"Pz",(mmarker[iov] ? mmarker[iov] : kFullCircle),
	    (mcolor[iov] ? mcolor[iov] : kBlack));
  } // for iov

    //after all the iov specific curves have been drawn, still add the 2022 Monte Carlo (or whatever else) --> TO DO: generalise this to a loop
    //for example, add "-" to the data list and when there is "-" the corresponding MC will only be added to upper plot and not to ratio. <-- TO DO
    c1->cd(1);
    TFile *fmc(0);
    TFile *fmcbpix(0);
    //fmc = new TFile(Form("rootfiles/GamHistosFill_mc_2022P8_v32.root"));
    fmc = new TFile(Form("rootfiles/GamHistosMix_mc_2023P8QCD_w13.root")); //add P8QCD (mix)
    fmcbpix = new TFile(Form("rootfiles/GamHistosMix_mc_2023-BPixP8QCD_w13.root")); //add P8QCD-BPix (mix)
    TObject *omc = fmc->Get(Form(so.c_str(),"MC")); assert(omc);
    TObject *omcbpix = fmcbpix->Get(Form(so.c_str(),"MC")); assert(omcbpix);

    TH1D *hmc(0);
    hmc = (TH1D*)omc;
    assert(hmc);
    //tdrDraw(hmc,"H",kNone,(mcolor["2022P8"] ? mcolor["2022P8"] : kBlack),kSolid,-1,kNone);
    tdrDraw(hmc,"H",kNone,kMagenta+2,kSolid,-1,kNone);


    TH1D *hmcbpix(0);
    hmcbpix = (TH1D*)omcbpix;
    assert(hmcbpix);
    tdrDraw(hmcbpix,"H",kNone,kRed+2,kSolid,-1,kNone);


    //leg->AddEntry(hmc,"2022P8","PLE");
    leg->AddEntry(hmc,"2023P8QCD","PLE");
    leg->AddEntry(hmcbpix,"2023P8QCD-BPix","PLE");




  //finish everything.
  if (id!="")
      c1->SaveAs(Form("pdf/drawPhotonJetVsPtVsIOVs-comparisons_%s_%s.pdf",
                     name.c_str(),id.c_str()));
  else 
    c1->SaveAs(Form("pdf/drawPhotonJetVsPtVsIOVs-comparisons_%s.pdf",name.c_str()));
} // void drawPhotonJetVsPtVsIOVs
