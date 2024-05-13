HOW TO RUN on Hefaistos:
------------------------
(try mosh to not drop connection)
- from local: 'rsync -rutP GamHistosFill.C GamHistosFill.h Hefaistos:/media/storage/gamjet/'
[- source /work/data/root/bin/thisroot.sh]
- source /work/data/rootbinaries/root/bin/thisroot.sh
- (rm *.d *.so *.pcm)
- root -l -b -q mk_CondFormats.C
- #define GPU in mk_GamHistosFill.C
=> edit (version, IOV_list) and execute 'python runAllIOVs.py'
[- nohup root -l -b -q mk_GamHistosFill.C > log.txt & ]
+ runtime about 12h, up from 9h (v19: 09:24 -> 21:18 18D2)
=> edit (version, IOV_list) and execute 'python renameAllIOVs.py'

+ tail -f log.txt
+ starting up takes quite a bit (~100 sec) due to GetEntries() call
   => code TStopWatch to ignore this startup time, or skip GetEntries
      (1M entries is just 0.1% ?! => 1.39B entries, about 3.4h runtime)

Local analysis:
- rsync -rutP files from Hefaistos
[  - rename GamHistosFill_mc_2016*APV to mc_2016APV* (automatize later)]
- python addAllIOVs.py
- root -l -b GamHistosMix.C+g
- root -l -b GamHistosRatio.C+g

Extra plotting
- root -l -b -q drawPhotonJetVsPtVsIOV.C
- root -l -b -q minitools/drawResponseVsRun.C
- root -l -b -q minitools/drawL2Res.C+g\(\"batchRun\"\)


TIMING estimates, one at a time
-------------------------------
2023C 5762->4851s* (9650s)
2023D 1466->1323s (3816s)
2022E 2670->2355s (6691s)
2023F 8104->7116s*** (23000s, or 6.3h)
2022G 1746->1215s (4304s)
2023B 843->687s (2482s)
2023Cv123 3556->2920s (10462s)
2023Cv4 6172->5109s** (15983s, or 4.4h)
2023D empty
2022P8 (14866s, or 4.1h)
2022EEP8 (scaled 52943s or 14.7h)

CHANGELOG (latest first):
----------------------------
To-do later: 
- update h2r9vspt r9 binning. Check if r9 cut makes sense at pT>230 GeV for Photon200 trigger, and if efficiency drops too much

w18: Golden JSON (03.05.), but with new version V2 of L2L3Residuals and jetvetomap. Otherwise like w17.
w17: New DAILY json (12.5.), more data; new L2L3Residual corrections (V2); new jetvetomap (V2). (data + json from 12.05.2024)
w16: Code unchanged (same as w15), but using new GOLDEN json. (data + json from 05./06.05.2024).
w15: New DAILY json and more data; new L2L3Residual corrections (Prompt24); new jetveto map (2024BC). (data + json from 05./06.05.2024).
w14: For runnng the code with Golden JSON (which was released in the morning of 30.04.2024).
w13: Using 2024Winter corrections for L2Relative jec, using daily JSON from 29.04.2024 (is the same as 28.04.2024)
w12: changed eta-binning in Gamjet2 folder (using now dijet narrow eta binning)
w11: first time running with 2024 data
w10: upper pT bins in some histograms were combined for better statistics. (might want to change this back?)
w9: applying jetvetomap also to photons (for all data periods: Cv123, Cv4, D)
w8: fixed bug --> using correct L2L3Res for the RunD data
w7: first run with new L2L3Res (29.02.2024)
w6: for running with L2Rel + L2Res (27.02.2024) on new QCD monte carlo
w5: used L2Rel + L2L3Res --> turned out to be not good (issue with L3-residuals)
w4: for running with L2Rel + L2Res (16.02.2024)
w3: used for running the bpix samples. Added line for setting jet phi, so that FactorizedJetCorrector does not complain. 
w2: changed JEC to Summer23Run3_V1_MC_L2Relative_AK4PUPPI.txt and for BPix samples Summer23BPixRun3_V1_MC_L2Relative_AK4PUPPI.txt respectively.
w1: first version including my updates - used to run on 2023 DATA (Cv123, Cv4, D) and 2023 MC (GJ-4Jets_HT* Summer23, and the same with BPix fix) (blehtela)


v32b: rerun QCD adter filtering more bad files.
v32: Update to 2022Sep2023 JEC V3. Add debugFiles to catch bad files.

v31b - Add debugFiles printout to catch bad files better. Replace photon idx assert with error printout and continue. Add printout for photon dR also. 
v31 - Increase Jet vector size from 100 to 200 just in case. Fix Int_t to Short_t for Photon_jetIdx and GenJet_partonFlavour in NanoV12.

v30b - Fix mistyped jec switch for 2022EEQCD.
v30 - Add new L2L3Res for 22Sep2023. Add doGamjet with vs pTave, pTjet. Fix if(pass) to doGamjet2. Add pr30chf etc. for PFcomp vs run. Add pass_filt to more trigger plots.

v29b - Remove 2022G from is22v10 patch.
v29 - Remove 2022[EE]P8 and 2022[EE]QCD from is22v10 patch.

v28 - Add prgamorigqcd[2]. Move METfilters earlier befor trigger plots (try to clean 700 GeV bump out in 2022FG). Update all but 22G (23B) to 22Sep2023-NanoV12.

v27c => should clean also 2022C and 2022F more, same for 2022EEQCD.
v27b - 2022F failed due to bad files, clean and try again.
v27 - Scale Gain1 path by 1.7% based on new minitools/drawGainVsPt.C. Add control/pr230b, prg1b etc.

v26 - Fix Run3 trigger for QCD bkg EM-jets. Add gam.DeltaR(jeti) for 2022CDE re-reco safety. Fix deltaPhi>2.7 to |deltaPhi|>2.7. Update MET filters list. Update doGamJet2 L2Res binning. Add control/pcutflow. Add control/pr110b etc. Update 2022F to 22Sep2023.v12. (Summer22(EE) MC still downloading, same with 2023 DATA.)

v25 - Flag_METFilters on (pass_filt) for Run3. Update 2023D file list. Add QCD to runAllIOVs.py.

v24(patch3) - Fix MC weighing for MadGraph QCD and add it
v24(patch1,2) - Add is22v10 to apply phoj patch to only 2022 prompt data and MC in v10/v11, not to new 2022CDE re-reco samples (NanoV12).
v24 - Fix fixedGridRho->fixedGridRhoFastjetAll, add safety for NaN. Add versioning to files. Sync JEC to Mikel (L2Relative only). Update 2022D, 2022E file lists to 22Sep2023 partial re-reco

v23 - Fix rawmet to use RawPuppiMET. Patch rawgam in PuppiJet in 2022 data and MC for phoj0 (MET ok). Return pass_leak. Add re-JEC + pdjes (DeltaJES). Add Gamjet2 folder (doGamJet2). Extend h2gametaphi to EC, although may not help with just 110EB trigger active. Remove h2gametaphi3,4 as redundant.

v22 - Switch HLT_Photon20 off for Run3. Switch 1% scaling to gain1 off for 22-23. Switch footprint correction off for 22-23. Switch off L1RC (corrl1rc was for AK4PFchs). Todo: override 75m,90m with 50m for 2023C. Override 50m with 30t for others. Next round: determine footprint from control/pphoj0.

v21e - Bybass data trigger to check trigger menu
v21 - Start implementing Run3. Add new JSON files. Add new jetvetomaps. Add new L1RC (and switch it off). Add new triggers. Disable _pu and PU reweighing for 22-23 (temporarily). Disable pass_leak and pass_filt for 22-23 (temporarily). Add code for MadGraph HT bins from dijet and add G-4Jet xsec.
=> HLT_Photon20 often missing, lots of errors on it

v20a - Rename 2016QCDAPV and 2016P8APV to 2016APVQCD and 2016APVP8
v20 - add pileup and PS reweighing for MC and lumi normalization for data
19f - comment out branches not available for 2017QCD => fixed in v20
19f - UL2017QCD files re-updated Dec 7
19f - add 'rho' in flavor, reset genjet and genjet2 for every event, require iGenJet!=-1 for 'rjet' and 'gjet', safety vs gen match to both recjet1 and recjet2
19f - add 'Rho_CHS' into BasicHistos
19f - add 'hgamtrg' with 'hgam' bins: wider bins and higher pt than 'hgamtrig'
v19e - add separate "ud" and "s" in 'flavor' folder in addition to "q"
v19d - change QCD photon match 'hasorig' from dR<0.4 to dR<0.2
v19c - add windowing for prgamqcd2, pfakeqcd2 to [0.8,1.1]
v19b - Add pfakeqcd, prgamqcd, h2rgamqcd to better estimate QCD bkg. Expect x5.5 drop in fake rate and fake photon scale close to unity
v19 - Add PF composition variables (+6). Add folder 'pf' with 1D composition and 2D response+composition. Fix missing |eta,gamma|<1.3 cut. Add DB(=MPF1).
v19x - Scale UE correction by jet aread (~0.5). Remove UE from phoj. Fix nFilesMax from 8 to 100 for MC in mk_GamHistosFill (had up to 24).
v18 - Remove UE from QCD gam and fox. Rerun everything to get data mu's fixed
v17 - Add QCD bkg, back to medium b-jet ID, gain1 to 1.01 (from 1.02)
v16 - Change b-jet ID to medium to match ttbar analysis. Extra 2D with yet finer phi granularity
=> reduce gain1 by 1% on-line (now 2%)
v15 - Add flavor breakup for data and MC. Runtime now 9.1h, compared to  6.7h in v12
v14 - Add missing pileup file, extra photon 2D with finer phi granularity. Photon 2D map to use any active triggers
v13 - Extra plots on gains, non-zero suppressed footprint plots.
v12 - Cut leak<0.06*ptgam, add 2D plots of H/E and R9 vs pTgam
v11 - add R9>0.94 cuts to avoid steps between triggers and data/MC difference; add eeBadScFilter for MC; remove ecalBadCalibFilter from UL16; implement jet veto also with 2016, 2017 maps
    => R9vsPt huge improvement, H/E similar (still 2018 step), Leakage improved at pT<60 GeV
    => MPF down by ~0.5% overall, MPF1 also at 200-300 GeV, MPFn slightly down (MC up), MPFu quite stable
v8 - extending to 2016 and 2018 data =>
read privileges missing on some files
missing Photon_eCorr on some files?
missing HLT_Photon20, 33, 200, 100EBTight, 110EBTight, 120EBtight, 
v7 - add ecalBadCalibFilter (+reweight MC pileup), back to old pT binning
v6 - reoptimized trigger selection and pT binnig
v5 - rerun of v4 since update failed
v4 - add hotjets veto
v3 - add JSON filter
v2 - add MET filters
v1 - initial version, keeping all events from EGamma for controls

LUMINOSITIES from Laura
------------------------
https://cernbox.cern.ch/index.php/s/MDX4TdtJZTqQCem
tail -n3 *_out.csv | grep '==>' -a1 | awk {'print $2" "$4}'

2016
HLT_Photon22_R9Id90_HE10_IsoM_out.csv recorded 0.019308130
HLT_Photon30_R9Id90_HE10_IsoM_out.csv recorded 0.098493381
HLT_Photon36_R9Id90_HE10_IsoM_out.csv recorded 0.220584079
HLT_Photon50_R9Id90_HE10_IsoM_out.csv recorded 0.506332692
HLT_Photon75_R9Id90_HE10_IsoM_out.csv recorded 2.599433386
HLT_Photon90_R9Id90_HE10_IsoM_out.csv recorded 5.230487368
HLT_Photon120_R9Id90_HE10_IsoM_out.csv recorded 14.673422774
HLT_Photon165_HE10_out.csv recorded 36.329675572
HLT_Photon165_R9Id90_HE10_IsoM_out.csv recorded 36.329675572
HLT_Photon175_out.csv recorded 36.329675572
2017
HLT_Photon20_HoverELoose_out.csv recorded 0.000447966
HLT_Photon20_out.csv recorded 0.000000601
HLT_Photon30_HoverELoose_out.csv recorded 0.002159100
HLT_Photon33_out.csv recorded 0.014177882
HLT_Photon50_R9Id90_HE10_IsoM_out.csv recorded 0.307208422
HLT_Photon75_R9Id90_HE10_IsoM_out.csv recorded 1.326146145
HLT_Photon90_R9Id90_HE10_IsoM_out.csv recorded 2.640888480
HLT_Photon120_R9Id90_HE10_IsoM_out.csv recorded 7.797325137
HLT_Photon165_R9Id90_HE10_IsoM_out.csv recorded 28.315865389
HLT_Photon200_out.csv recorded 41.477735690

TRIGGER MENUs
------------
[2016] https://docs.google.com/spreadsheets/d/1bII_92pCrgk20A9FMIIHsOsYYj3f-lLjsoRkP_ZNQW4/edit?usp=sharing
[2017] https://docs.google.com/spreadsheets/d/1SqSKATiHcZdSB-t7rMb0SRSUL4Ot1jr2TExjnMrDUpA/edit?usp=sharing
[2018] https://docs.google.com/spreadsheets/d/1D_og1_J6Hp4uALUg-R4Hkq3CF0VN6IK5komHPHv5R-o/edit?usp=sharing

TO-DO LIST
-----------
- Complete implementation of PU reweighing. Add PU variations (systematic)
- Implement JER SF (once pT-dependent?). Add JER SF variations (systematic)
- Add parton shower weights (systematic)
- Add flavor corrections and its variations (systematic)


-- Correctgamscale seems critical for photon+jet. At low pT we probably don't have pT,e ~ pT,Z/2 anymore, and low pT photons may require more correction like high pT ones (quadratic shape with minimum at pT,Z~40 GeV). Should also calibrate with superclusters instead of electrons (track+EM more like track at low pT?)
** add veto on bad photon regions?

TBD:
- Investigate 1.5% leakage at high pT, photon vs PFphoton inconsistent? (h2phoj)
  - Could explain MPF lower than MPF1. Leakage also higher in data
- Add cut on maximum leakage at low pT to reduce QCD bkg? (h2phoj in data,MC)
  - do data/MC ratio of normalized leakage in pT bins 
- Investigate non-continuous H/E vs pT
- Veto photons near ECAL holes? And near SC boundaries? (h2gametaphi1,2)
  => at least ECAL holes, some missing from MC

tbd
- fix prhovsmu for data_Run2 => missing pileup_ASCII file from GPU machine?

tbd
+ reweight MC pileup (brilcalc for target profiles)
+ add JER SF

Brilcalc luminosity
cd /afs/cern.ch/user/v/voutila/scratch0/CMSSW_10_2_25
cmsenv
setenv PATH $HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i [your json] --hltpath [your !HLT path] -o output.csv

Luminosity recommendations:
https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#SummaryTable
https://twiki.cern.ch/twiki/bin/viewauth/CMS/BrilcalcQuickStart

Plan:
- fix running over 2018A so get complete sample => done
  - increased nPhotonMax to 200, crashed with one event having 191
- implement JSON, rerun v3
  - /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt
  - /eos/home-c/cmsdqm/www/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json
  - /eos/home-c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json
- implement hotzones, rerun v4,v5
- implement JER SF, rerun v6 (MC only?)

To-do list:
- Check that runs on data => done
- Implement merging of 2018A,B,C,D => hadd (but in the future, give data luminosity weight for correct averaging?)
- Implement data/MC ratio => ok
- Run on full 2018 data and MC (flat?) => v1 ok, v2 on-going
- Implement photon corrections and uncertainties
   => 2018 MPFdata bit low, bal less different => photon footprint?
   - Look into photon jet minus photon minus L1RC. Separate CHF, NEF, NHF
- Jet and MET smearing? => careful, check shift on jet scales
- Rho vs pT mu etc.
- Select exactly one photon? => make plot of these vs pTgamma
- MET cleaning => done; add control plot
- Study photon seedGain vs pT, energyErr vs pT, hoe vs pT, r9 vs pT
- Check photon eta-phi distribution to locate hot/cold zones => draft, update
- Compare to Higgs to gamma gamma workflow

 Events->Draw("Photon_energyErr[0]/Photon_pt[0]:Photon_pt[0]","fabs(Photon_eta[0])<1.3 && Photon_pt[0]>30 && Photon_cutBased[0]==3","")

Lucas Torterotot, latest talk on gamma+jet I could find (29 July 2019, EOY18):
- https://indico.cern.ch/event/837707/#10-gamma-jet-l3res-for-2018-v1
Code through Stephanie:
- https://github.com/lucastorterotot/DijetRootTreeMaker/blob/master/instructions/GammaJetTree_Instruction.md

Latest runs v2, v3 (trig+JSON; all in parallel, except 2018a in v2 and MC in v3)
- 2018A: runtime 3453 sec. => 6387 s (alone) / 9307 s (~14/fb)
- 2018B: runtime 5169 sec. => 6167 s / 6262 s (~7.1/fb)
- 2018C: runtime 5110 sec. => 6102 s / 6217 s (~6.9/fb)
- 2018D: runtime 14631 sec. => 14834 s / 16220 s (~31/fb)
- 2018MC: runtime 1075 sec. => 910 s / 452 s (alone)
=> would be justified to split 2018D in two to go from 4.1h to 2.1h

Sami's code:
- https://gitlab.cern.ch/slehti/ZbAnalysis/-/blob/master/src/Analysis.h
  - return false; for event selection
- https://gitlab.cern.ch/slehti/ZbAnalysis/-/blob/master/analyse.py
  - inputList = ROOT.TList() kommunikoi C++:aan

Observations:
- bit trick to compile and run on Hefaistos, differet ROOT version?
- HLT_Photon[100,110,120]EB_TightID_TightIso correspond to a lot of statistics in midrange, but are apparently not active in all of 2018 data => anything replacing them?
- can remove HLT_Photon100EBHE10, HLT_Photon100EEHE10, HLT_Photon100EE_TightID_TightIso also, not always active and not used at the moment
- data has 1.4B events, takes a lot of time to run (16.6h). Photon triggers are only a fraction of this. As a quick fix, activate trigger pre-selection for data
- first run took 51326 sec (14.5h), crashed before starting MC, but writing file succeeded (or not?)
- HLT_PhotonXX pass through some noise at 200 GeV that 110T does not?


// NanoAOD info: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD
// Variables: https://cms-nanoaod-integration.web.cern.ch/integration/cms-swCMSSW_10_6_19/mc102X_doc.html

	// Tight cut-based photon ID:
	// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Offline_selection_criteria
	// "Do not implement your ID by hand, instead use VID" (94X-V2 ID)
	//
	// NB: if using generic photon ID for high pT region (>200 GeV),
	// some detailed studies of efficiency and SF are needed:
	// https://twiki.cern.ch/twiki/bin/view/CMS/EGMPhotonIDHighPtPhotons
	//
	// NanoAOD information:
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Photons
	// https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/data102X_doc.html#Photon
	// Photon_cutBased includes isolation. Checked w/ cutBased==3 && pT>100 GeV
	// However, could only check charged, as neutral and photon not separated
	// Also checked sieie. The HoverE is not applied!? 2.4% leakage

	// Select photon matching with LHE photon (always #2)
	// Only possible in MadGraph samples
	// PythiaFlat has GenIsolatedPhotons (MG has <50% of these, P8 >90%)


	/*
	for (unsigned int j = 0; j != nLHEPart; ++j) {
	  if (LHEPart_status[j]==1 && LHEPart_pdgId[j]==22) {
	    assert(j==2);
	    assert(lhe.Pt()==0);
	    lhe.SetPtEtaPhiM(LHEPart_pt[j],LHEPart_eta[j],
			     LHEPart_phi[j],LHEPart_mass[j]);
	    gam.SetPtEtaPhiM(Photon_pt[i],Photon_eta[i],
			     Photon_phi[i],Photon_mass[i]);
	    if (lhe.DeltaR(gam)<0.4 && iGamGen==-1) {
	      iGamGen = i;
	    }
	  } // LHE photon
	}
	*/
