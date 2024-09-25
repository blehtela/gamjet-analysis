// Purpose: draw stability of photon+jet results vs pT over time
//          compare data and MC for response, HOE, R9 etc.
// Note: This is essentially a copy of the original drawPhotonJetVsPtVsIOV.C, adjusted to our needs while deriving new correctiongs for 2023. I removed everything not needed at this point.
// Note: starting point for this was my version drawPhotonJetVsPtVsIOV_2024only.C
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TLine.h"

#include "tdrstyle_mod22.C"

bool addMPFu2n = true;
bool addG1toMPF = false;//true;
bool addG12toMPF = false;
string id = "w37"; //comparing 2024 stuff
//string id = "various comparions"; //displaying even more data and mc results.


bool drawFullIOVList = false;//true; //need to use the short list, haven't updated the long lists properly

// Forward declaration of call
void drawRhoNpvMu_2024only(string so, string var, string name, double y1, double y2);
			     //double y1, double y2, double z1, double z2); //y1, y2 upper plot, z1, z2 lower plot

//void drawPhotonJetVsPtVsIOV(string so = "resp_MPFchs_%s_a100_eta00_13") {
//void drawPhotonJetVsPtVsIOVs(string so = "control/phoevspt",
//			     string var = "H/E",
//			     double y1 = 0, double y2 = 0.01,
//			     double z1 = 0.7, double z2 = 1.4) {
			    
// Multiple calls to draw function
void drawRhoNpvMu_2024only() {
  drawRhoNpvMu_2024only("control/prhovspt","Rho","Rho",20.,60.); //0.,50.
  drawRhoNpvMu_2024only("control/pnpvgoodvspt","NPVgood","NPVgood",20.0,60.0);
  drawRhoNpvMu_2024only("control/pnpvallvspt","NPVall","NPVall",20.0,60.0);
	drawRhoNpvMu_2024only("control/pmuvspt","Mu","Mu",0.0,100.0);
} // drawPhotonJetVsPtVsIOVcomparisons

void drawRhoNpvMu_2024only(string so, string var, string name, double y1, double y2){
			     //double y1, double y2, double z1, double z2) {
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  string iovs[] = {"2024C", "2024D", "2024E", "2024F", "2024G", "winter2024P8"}; 
//, "2024C-ECALCC-HCALDI", "winter2024P8", "2023P8-BPix"};
	//string reweighted_mc = "w37";

  const int niov = sizeof(iovs)/sizeof(iovs[0]);
  //string *iovs = &iovs[0];

  map<string,int> mcolor;
  map<string,int> mmarker;
  //2024
  mmarker["2024B"] = kFullDiamond;//kFullCircle;
  mmarker["2024C"] = kFullTriangleDown;//kFullDiamond;
  mmarker["2024D"] = kFullCircle;//kOpenCircle; 
	mmarker["2024Ev1"] = kFullSquare;//kOpenSquare;
	mmarker["2024Ev2"] = kOpenTriangleUp;
	mmarker["2024E"] = kFullSquare;
	mmarker["2024F"] = kFullTriangleUp;
	mmarker["2024G"] = kFullDiamond;

  mmarker["2024B-ECALRATIO"] = kFullDiamond;
  mmarker["2024C-ECALRATIO"] = kFullTriangleDown;
  mmarker["2024C-ECALR-HCALDI"] = kFullTriangleDown;
  mmarker["2024C-ECALCC-HCALDI"] = kFullTriangleDown;

  mmarker["winter2024P8"] = kOpenCircle;
  mmarker["2023P8-BPix"] = kOpenSquare;


  mcolor["2024B"] = kBlue;//kGreen+2;
  mcolor["2024C"] = kGreen+2;
  mcolor["2024D"] = kOrange+2;
 	mcolor["2024Ev1"] = kCyan+2;
  mcolor["2024Ev2"] = kMagenta+2;
 	mcolor["2024E"] = kCyan+2;
  mcolor["2024F"] = kMagenta+2;
  mcolor["2024G"] = kBlue+2;

 	mcolor["2024B-ECALRATIO"] = kBlue+3;
  mcolor["2024C-ECALRATIO"] = kGreen+3;
  mcolor["2024C-ECALR-HCALDI"] = kGreen+4;
	mcolor["2024C-ECALCC-HCALDI"] = kGreen+4;

  mcolor["winter2024P8"] = kGreen+2; //reweighted //kBlack;
  mcolor["2023P8-BPix"] = kGray+2;


 
  
  const char *cvar = var.c_str();
  const char *cname = name.c_str();

  TH1D *h = tdrHist("h",cvar,y1,y2);
  //TH1D *h2 = tdrHist("h2","Data/MC(noQCD)",z1,z2);
  //TH1D *h2 = tdrHist("h2","Data/MC",z1,z2);

  //lumi_13TeV = "Run2 v16";
  //lumi_13TeV = "Run3 v21";
  lumi_13TeV = "Run3"; // 4=13 TeV
  if (id!="") lumi_136TeV = Form("Run3 %s",id.c_str()); // 8=13.6 TeV
  //TCanvas *c1 = tdrDiCanvas(Form("c1_%s",cname),h,h2,8,11);
  TCanvas *c1 = tdrCanvas(Form("c1_%s",cname),h,8,11);

	//custom place for "Private" extraText (usually handled via tdrCanvas)
	TLatex *extratext = new TLatex();
	extratext->SetTextSize(0.039);
	extratext->SetTextFont(52);
	extratext->DrawLatex(40, 56, "Preliminary");

  c1->cd(1);
/*
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(h->GetXaxis()->GetXmin(),1,h->GetXaxis()->GetXmax(),1);
  l->DrawLine(h->GetXaxis()->GetXmin(),0,h->GetXaxis()->GetXmax(),0);
*/
  gPad->SetLogx();

  //TLegend *leg = tdrLeg(0.60,0.90-niov*0.06,0.80,0.90);
  //TLegend *leg = tdrLeg(0.60,0.90-niov*0.045,0.80,0.90);

  //TLegend *leg = tdrLeg(0.60,0.90-niov*0.08,0.80,0.90); //this is good for fewer graphs
	//TLegend *leg = tdrLeg(0.40,0.90-niov*0.03,0.60,0.90); //moved a bit to the left and smaller font
	TLegend *leg = tdrLeg(0.50,0.90-niov*0.03,0.60,0.90); //moved a bit to the left and smaller font
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);


  //TLegend *leg = tdrLeg(0.60,0.90-niov*0.06,0.80,0.90);

  leg->SetTextFont(42);




/*
  c1->cd(2);
  gPad->SetLogx();
  l->DrawLine(h->GetXaxis()->GetXmin(),1,h->GetXaxis()->GetXmax(),1);
  l->DrawLine(h->GetXaxis()->GetXmin(),0,h->GetXaxis()->GetXmax(),0);
*/

  for (int i = 0; i != niov; ++i) {

    string iov = iovs[i];
    const char *ciov = iov.c_str();
    //const char *cmc = mcs[i].c_str();
    const char *cid = id.c_str();

    TFile *fd(0), *fm(0);

    //getting the correct input files (as they also differ in code version)
/*
    if (iovs[i]=="2023Cv123" || iovs[i]=="2023Cv4" || iovs[i]=="2023D"){ //|| 
        //iovs[i]=="2024B" || iovs[i]=="2024C" || iovs[i]=="2024D") {
      fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_w17.root",ciov)); //newest ones
      fm = new TFile(Form("rootfiles/GamHistosFill_mc_%s_w13.root",cmc));    //newest ones
      //fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_wX22full.root",ciov)); //newest data but with 2022 correctiongs
      //fm = new TFile(Form("rootfiles/GamHistosFill_mc_%s_wX22full.root",cmc));    //newest  data but with 2022 correctiongs
    }
		else{
      fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_w29.root",ciov)); //newest ones
      fm = new TFile(Form("rootfiles/GamHistosFill_mc_%s_w13.root",cmc));    //newest ones
		}
*/

/*
    if (iovs[i]=="2018ABCD") {
      fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_v20.root",ciov));
      fm = new TFile(Form("rootfiles/GamHistosFill_mc_%s_v20.root",cmc));
    }
*/

		//if(TString(ciov.c_str()).Contains("P8")){// Monte Carlo
		if(iov=="winter2024P8" || iov=="2023P8-BPix"){// Monte Carlo
    	fd = new TFile(Form("rootfiles/GamHistosFill_mc_%s_%s.root",ciov,cid));
		}
		else{ //data
    	fd = new TFile(Form("rootfiles/GamHistosFill_data_%s_%s.root",ciov,cid));
		}
 
    assert(fd && !fd->IsZombie());
    //assert(fm && !fm->IsZombie());


    curdir->cd();
    
    TObject *od = fd->Get(Form(so.c_str())); assert(od);
    //TObject *om = fm->Get(Form(so.c_str(),"MC")); assert(om);
    
    TH1D *hd(0), *hm(0);
/*
    if (od->InheritsFrom("TProfile")) {
      hd = ((TProfile*)od)->ProjectionX(Form("hd_%s_%s",cvar,ciov));
      //hm = ((TProfile*)om)->ProjectionX(Form("hm_%s_%s",cvar,ciov));

      if (name=="MPFn" && addMPFu2n) {
	const char *co = "resp_MpfRuchs_%s_a100_eta00_13";
	TProfile *pd = (TProfile*)fd->Get(Form(co,"DATA")); assert(pd);
	hd->Add(pd,1.2);
      }
*/
    
    hd = (TH1D*)od;
		cout << "Got histo from file " << fd << endl << flush;
    assert(hd);
   

    c1->cd(1);
    gPad->SetLogx();
    //tdrDraw(hm,"H",kNone,(mcolor[iov] ? mcolor[iov] : kBlack),kSolid,-1,kNone); //draw it extra later in red, as only bpix used
    tdrDraw(hd,"Pz",(mmarker[iov] ? mmarker[iov] : kFullCircle), (mcolor[iov] ? mcolor[iov] : kBlack));
    
    leg->AddEntry(hd,ciov,"PLE");
    
 

    c1->cd(2);
    gPad->SetLogx();
  } // for iov

/*
    //after all the iov specific curves have been drawn, still add the 2022 Monte Carlo (or whatever else) --> TO DO: generalise this to a loop
    //for example, add "-" to the data list and when there is "-" the corresponding MC will only be added to upper plot and not to ratio. <-- TO DO
// TO DO (note from 20.08.2024): could add the winter2024P8-QCD mc here.
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

		//NOTE: I DO NOT DRAW THIS FOR 2024 BECAUSE WE USE BPIX FOR ALL OF THEM
    //tdrDraw(hmc,"H",kNone,kMagenta+2,kSolid,-1,kNone);


    TH1D *hmcbpix(0);
    hmcbpix = (TH1D*)omcbpix;
    assert(hmcbpix);
    tdrDraw(hmcbpix,"H",kNone,kRed+2,kSolid,-1,kNone);


    //leg->AddEntry(hmc,"2022P8","PLE");
    //leg->AddEntry(hmc,"2023P8QCD","PLE");
    leg->AddEntry(hmcbpix,"2023P8QCD-BPix","PLE");

*/

		//Now still adding the MC without reweighting (w36)
		c1->cd(1);
    TFile *fmc_noreweight(0);
   	fmc_noreweight = new TFile(Form("rootfiles/GamHistosFill_mc_winter2024P8-v14_w36.root")); //last without reweighting so far, w36
    //TObject *omc = fd->Get("control/prhovspt"); assert(omc);
		TH1D *hmc = (TH1D*)fmc_noreweight->Get("control/prhovspt"); assert(hmc);

		tdrDraw(hmc,"Pz",kOpenCircle,kRed+2,kSolid,-1,kNone);
		leg->AddEntry(hmc,"winter2024P8 without PU reweighting (w36)","PLE");


  //finish everything.
  if (id!="")
      c1->SaveAs(Form("pdf/drawRhoNpvMu-2024only_%s_%s.pdf",name.c_str(),id.c_str()));
  else 
    c1->SaveAs(Form("pdf/drawRhoNpvMu-2024only_%s.pdf",name.c_str()));
} // void drawRhoNpvMu_2024only
