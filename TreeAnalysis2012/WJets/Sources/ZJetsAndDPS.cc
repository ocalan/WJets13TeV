
#define PI 3.14159265359
#define BARREDEPROGRESSION 0
#define DEBUG 0
#define PRINTEVENT 0

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <RooUnfoldResponse.h>
#include <TDatime.h>
#include <TMath.h>
#include <TRandom3.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "LHAPDF/LHAPDF.h"

#include "functions.h"
#include "getFilesAndHistograms.h"
#include "standalone_LumiReWeighting.h"
#include "funcReweightResp.h"
#include "HistoSet.h"
#include "ZJetsAndDPS.h"

using namespace std;

TRandom3* RandGen = new TRandom3();
TRandom3* Rand_MER_Gen = new TRandom3();

ClassImp(ZJetsAndDPS);

void ZJetsAndDPS::Loop(bool hasRecoInfo, bool hasGenInfo, int doQCD, bool doSSign, bool doInvMassCut, 
        int doBJets, int doPUStudy, bool doFlat, bool useRoch, bool doVarWidth,  bool hasPartonInfo, string pdfSet, int pdfMember)
{

    //--- Initialize PDF from LHAPDF if needed ---
    if (pdfSet != "") {
        LHAPDF::initPDFSet(1, pdfSet.c_str(), pdfMember);
        LHAPDF::initPDFSet(2, "CT10.LHgrid");
        const int numberPDFS(LHAPDF::numberPDF() + 1);
        if (pdfMember > numberPDFS) {
            std::cout << "Warning pdfMember to high" << std::endl;
            return;
        }
    }
    //--------------------------------------------

    //--- Check weither it is 7 or 8 TeV ---
    string energy = getEnergy();
    //--------------------------------------
 if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    //--- Counters to check the yields ---
    unsigned int nEvents(0), nEventsIncl0Jets(0), nEventsUNFOLDIncl0Jets(0);
    unsigned int nEventsWithTwoGoodLeptonsNoChargeNoMass(0), nEventsWithTwoGoodLeptonsNoMass(0), nEventsWithTwoGoodLeptons(0);
    unsigned int nEventsExcl0Jets(0), nEventsExcl1Jets(0), nEventsExcl2Jets(0), nEventsExcl3Jets(0),nEventsIncBJets(0);
    unsigned int GENnEventsIncl0Jets(0), GENnEventsIncl1Jets(0), GENnEventsIncl2Jets(0), GENnEventsIncl3Jets(0);
    double TotalGenWeight(0.), TotalGenWeightPassGEN(0.), TotalGenWeightPassGENPU(0.), TotalGenWeightPassRECO(0.), TotalRecoWeightPassRECO(0.);
    //------------------------------------
    if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    //==========================================================================================================//
    double MTCut(50.);
    double ZMCutLow(71), ZMCutHigh(111);
    //------------------------------------
    bool doZ(true), doW(false), doTT(false), doDR(false), doTTreweighting(false);
    if (leptonFlavor == "SingleElectron" || leptonFlavor == "SingleMuon"){
        doW = true; 
        doDR = true;
    }
    if (leptonFlavor == "TTMuE") doTT = true; 
    if (doW || doTT) doZ = false;
    if (fileName.find("_dR_") != string::npos) doDR = true;
    if (fileName.find("TopReweighting") != string::npos) { hasGenInfo = true ; doTTreweighting = true;} // we don't want to use gen plots for ttbar, we just need to load the lepton branch to read the t qurak pt 
    if ( doZ ) METcut = 0; // no need for MET cut on Z+jets analysis 
    //if ( doW ) METcut = 30; // remove this if no MET cut
    
    // additional muons variables
    double leptonMass(0.00051);  // 
    int LeptonID(11);
    if (leptonFlavor == "Muons" || leptonFlavor == "SingleMuon"){
        leptonMass = 0.105658;
        LeptonID = 13;
    }
    //------------------------------------
    cout << " begin: "<< hasRecoInfo <<"  " << hasGenInfo <<"  " << doQCD<<"  " << doSSign<<"  " << doInvMassCut << "  " << METcut << "  " <<doBJets <<"  " <<doPUStudy << endl;
    
    
    //==========================================================================================================//
    //         Output file name           //
    //===================================//
    string command = "mkdir -p " + outputDirectory;
    system(command.c_str());
    string outputFileName = CreateOutputFileName(useRoch, doFlat, doPUStudy, doVarWidth, doBJets, doQCD, doSSign , doInvMassCut, pdfSet, pdfMember);
    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    //TFile *outputFile = new TFile("TEST.root", "RECREATE");
    //==========================================================================================================//

    
    if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    //==========================================================================================================//
    //       Load efficiency tables        //
    //====================================//
    table LeptIso, LeptID, LeptTrig, Ele_Rec;
    table TableJESunc("EfficiencyTables/JESUnce_FT_53_V21_AN4_Uncertainty_AK5PFchs.txt");
    if (energy == "8TeV"){
        if (leptonFlavor == "Electrons" || leptonFlavor == "SingleElectron"){
            /// electron SF
            table Ele_Rec_8TeV("EfficiencyTables/Ele_SF_Reconstruction_2012.txt");
            table SC_Ele_2012EA("EfficiencyTables/Ele_SF_EA2012.txt");
            
            Ele_Rec = Ele_Rec_8TeV ;
            LeptID = SC_Ele_2012EA;
        }
        if (leptonFlavor == "SingleMuon")  {
            // Single Muon SF for 2012 data (Rereco on Jan2013)
            table SF_Muon_IDTight_Rereco("EfficiencyTables2012/Muon_IDTight_Efficiencies_ReReco2012_Eta_Pt.txt");
            table SF_Muon_ISOTight_ReReco("EfficiencyTables2012/Muon_ISOTight_forTight_Efficiencies_ReReco2012_Eta_Pt.txt");
            table SF_TrigIsoMu24eta2p1_ReReco("EfficiencyTables2012/Efficiency_SF_ReReco2012_IsoMu24_eta2p1.txt");
            
            LeptID = SF_Muon_IDTight_Rereco;
            LeptIso = SF_Muon_ISOTight_ReReco;
            LeptTrig = SF_TrigIsoMu24eta2p1_ReReco;
        }
    }
    //==========================================================================================================//
    cout << "Phase space cuts -- jet pt:" << jetPtCutMin <<"  " << jetPtCutMax<<"  -- jet eta : " << jetEtaCutMin<< "  " << jetEtaCutMax<< "  " << "  -- Z eta: " << ZEtaCutMin<<"   " << ZEtaCutMax<< "  -- MET cut: " << METcut << "    "   << endl;
    cout << " other selections:  " <<endl;
    cout << " doQCD: " << doQCD <<"  do SS: " << doSSign <<" inv. mass cut: " << doInvMassCut <<"  use MET cut: " << METcut<<"  use B jets: " << doBJets <<" do PU study: " << doPUStudy << endl;


    if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    //==========================================================================================================//
    //     Systematics: jec, pu, xsec     //
    //====================================//
    cout << "Lepton Flavor: " << leptonFlavor << endl;
    int puYear(2011); 
    if (energy == "8TeV") puYear = 2013;
    cout << "Pile Up Distribution: " << puYear << endl;
    standalone_LumiReWeighting puWeight(leptonFlavor, puYear), puUp(leptonFlavor, puYear, 1), puDown(leptonFlavor, puYear, -1);
    cout << "systematics: " << systematics << "  direction: " << direction << endl;
    if (systematics == 1 && direction ==  1) puWeight = puUp;
    if (systematics == 1 && direction == -1) puWeight = puDown;

    int scale(0);//0,+1,-1; (keep 0 for noJEC shift study)
    if (systematics == 2 && direction ==  1) scale =  1;
    if (systematics == 2 && direction == -1) scale = -1;

    double xsec(1.);
    if (systematics == 3 && direction ==  1) xsec = 1. + xsecfactor;
    if (systematics == 3 && direction == -1) xsec = 1. - xsecfactor;

    int smearJet(0);
    if (systematics == 4 && direction ==  1) smearJet =  1;
    if (systematics == 4 && direction == -1) smearJet = -1;
    
    int sysLepSF(0);
    if (systematics == 5 && direction ==  1) sysLepSF =  1;
    if (systematics == 5 && direction == -1) sysLepSF = -1;
    
    int sysBtagSF(0);
    if (systematics == 6 && direction ==  1) sysBtagSF =  1;
    if (systematics == 6 && direction == -1) sysBtagSF = -1;
    
    double muScale(1.0);
    if (systematics == 7 && direction ==  1) muScale = 1.002;
    if (systematics == 7 && direction == -1) muScale = 0.998;
    
    bool doMer(false); // the number used for MER : 0.006
    double merUncer(0); Rand_MER_Gen->SetSeed(0); // set random seed; 0.006
    if (systematics == 8 && direction ==  1) doMer = true;
    
    // Wb study
    bool doWbsyst(false);
    double WbSystSF(1.3);
    if (systematics == 9 && direction ==  1) doWbsyst = true;
    
    bool doRespSyst(false);
    if (systematics == 10 && direction ==  1) doRespSyst = true;
    
    //int smearLepSF(0);
    //if ((systematics == 5 || systematics == 6) && direction ==  1) smearLepSF = 1;
    //if ((systematics == 5 || systematics == 6) && direction == -1) smearLepSF = -1;

    //==========================================================================================================//


    // initialize rochester corrrection
    //rochcor2012 *rmcor = new rochcor2012(); // make the pointer of rochcor class
    //REMARK : Need to call "rochcor(seed)" to assign the systematic error
    //rochcor2012 *rmcor = new rochcor2012(seed); //where "seed" is the random seed number
 if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    //---  Retreive the NVtx comparison histogram to have the exact weight to re-weight for pile-up
    // and have a flat NVtx distribution
    if (doFlat){
        string nVtxFileName = "Includes/";
        if (leptonFlavor == "Muons") nVtxFileName += "DMu_NVtx.root";
        else if (leptonFlavor == "Electrons") nVtxFileName += "DE_NVtx.root";
        TFile *DataMCRawComparison = new TFile(nVtxFileName.c_str());
        TCanvas *can = (TCanvas*) DataMCRawComparison->Get("NVtx");
        TPad *pad = (TPad*) can->FindObject("pad2");
        FlatNVtxWeight = (TH1D*) pad->FindObject("NVtx");
    }

    /// for Z counting
    int NZtotal = 0 ;

    double sumSherpaW = 0. ;
    cout << " Sherpa initial weight :  " << sumSherpaW <<endl;

    double sumTrig[4] = {0};
    int countTrigSum[4] = {0};

    // setting weight when running on MIX of exclusive DY/WJets files to match number of parton events
    double mixingWeightsDY[4] = {0.1926862,  0.07180968,  0.04943502,  0.03603373 }; // here we match all partons, and combine electron and muon side
    //double mixingWeightsDY[4] = {  0.1927289,  0.07199641,  0.04966403,  0.03631544 }; // here we match  partons that pas gen cuts, and combine electron and muon side
    //double mixingWeightsDY_DMu[4] = {0.1925615791, 0.07927772, 0.04974768769, 0.03640484898};// OLDDDD  here we match only those partons that pass the gen cuts
    //double mixingWeightsDY_DMu[4] = { 0.192623 , 0.0719199 , 0.0495369 , 0.03617}; //  here we match only those partons that pass the gen    cuts
    //double mixingWeightsDY_DE[4] = {0.192749004, 0.07880291002, 0.04933310117, 0.03589736768};// 
    //double mixingWeightsDY_DE[4] = {0.1929798803, 0.07896770167, 0.04960555955, 0.03619125322}; // here we match only those partons that pass the gen cuts
    double mixingWeightsWJ_SMu[4] ={0.366713,  0.1119323,  0.07641136,  0.03803325};
    double mixingWeightsWJ_SE[4] ={0.3667127984048746, 0.111932213229137, 0.076411344088767, 0.0380331330318}; // this need to be updated
    // Bugras weights: The weights: {2.4630522,0.473822,1.091304,0.121762,0.092715}
    //                  Events:      9566350,  2820770 ,969224   ,307529  ,133303

    //==========================================================================================================//
    // Start looping over all the events //
    //===================================//
    cout << endl;
    printf("\nProcessing : %s    -->   %s \n", fileName.c_str(), outputFileName.c_str()); 
    double weightSum(0.), weightSumNoTopRew(0.);

    //--- Initialize the tree branches ---
    Init(hasRecoInfo, hasGenInfo, hasPartonInfo);
    if (fChain == 0) return;
    Long64_t nbytes(0), nb(0);
    Long64_t nentries = fChain->GetEntriesFast();
    if (nEvents_10000) {
        nentries = 10000;
        std::cout << "We plane to run on 100000 events" << std::endl;
    }
    std::cout << "We will run on " << nentries << " events" << std::endl;
    //------------------------------------
    
    //--- get weight and bin edges for systematic for reweighted response
    vector<double> vecFwHT1; vector<double> binEdgeHT1;
    vector<double> vecFwHT2; vector<double> binEdgeHT2;
    vector<double> vecFwRap12; vector<double> binEdgeRap12;
    vector<double> vecFwRapFB; vector<double> binEdgeRapFB;
    if (doRespSyst){
        vecFwHT1     = getFW("MeanNJetsHT_Zinc1jet");
        binEdgeHT1   = getXBin("MeanNJetsHT_Zinc1jet");
        vecFwHT2     = getFW("MeanNJetsHT_Zinc2jet");
        binEdgeHT2   = getXBin("MeanNJetsHT_Zinc2jet");
        vecFwRap12   = getFW("MeanNJetsdRapidity_Zinc2jet");
        binEdgeRap12 = getXBin("MeanNJetsdRapidity_Zinc2jet");
        vecFwRapFB   = getFW("MeanNJetsdRapidityFB_Zinc2jet");
        binEdgeRapFB = getXBin("MeanNJetsdRapidityFB_Zinc2jet");
    }
    //------------------------------------
    
    cout << " run on " << nentries << "events" << endl;
    //--- Begin Loop All Entries --
    //KOfor (Long64_t jentry(0); jentry < nentries; jentry++){
    for (Long64_t jentry(0); jentry < 100000; jentry++){
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;

        if (jentry % 100000 == 0) std::cout << jentry << std::endl;
        nb = fChain->GetEntry(jentry);  
        nbytes += nb;
        nEvents++;

        if (BARREDEPROGRESSION && !DEBUG) barre_de_progression((int) 100*(jentry+1)/nentries);
 if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        //=======================================================================================================//
        //         Continue Statements        //
        //====================================//
        //if (jentry % 2 == 0) continue;
        //if (EvtInfo_NumVtx <= 14) continue;
        //if (PRINTEVENT && EvtInfo_RunNum != 190733) continue ;

        //=======================================================================================================//


        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        //=======================================================================================================//
        //         Computing weight           //
        //====================================//
        //--- weight variable ---
        double weight(1.);
        double genWeight(1.);
        //-----------------------
        double weightTop = 1.; 
        // for data PU_npT == -2
        // line below is to see distributions as provided with default MC PU distribution
        double reweighting(1);
        
/* APICHART
        if (hasRecoInfo && !isData){
            weight *= (double)puWeight.weight(int(PU_npT));
            //-- reweight again to IMPOSE FLAT #VTX DATA/MC RATIO
            if (doFlat){
                reweighting = FlatNVtxWeight->GetBinContent(EvtVtxCnt + 1);
                //-- for safety check the value of the weight...
                if (reweighting <= 0 || reweighting > 1000) reweighting = 1;
                weight *= reweighting;
            }
        }
        if (weight > 10000 || weight < 0) weight = 1;
//APICHART */
        
        weight = weight * lumiScale * xsec;
        
        if (fileName.find("MIX") != string::npos && nup_ > 5) {
            if (fileName.find("DYJets") != string::npos) weight *= mixingWeightsDY[nup_ - 6]; 
            else if (fileName.find("WJets") != string::npos && fileName.find("SMu_") != string::npos) weight *= mixingWeightsWJ_SMu[nup_ - 6]; 
            else if (fileName.find("WJets") != string::npos && fileName.find("SE_") != string::npos) weight *= mixingWeightsWJ_SE[nup_ - 6]; 
            else weight *= mixingWeightsDY[nup_ - 6];
        }
        // To compare generator distributions for comparison between Sherpa and MadGraph.
        // if below is mainly for Sherpa since events can have enhancement factors of different jet multiplicities
        if (fileName.find("mcEveWeight") != string::npos || fileName.find("MiNLO") != string::npos) {
            weight *= mcEveWeight_;
        }
        if (fileName.find("HepMC") != string::npos) {
            weight *= mcEveWeight_;
            //sumSherpaW += mcSherpaSumWeight3_ ;
        }
        if (fileName.find("Sherpa2") != string::npos){
            weight *= mcSherpaWeights_->at(0);
            sumSherpaW += mcSherpaWeights_->at(4);
        }

        //==========================================================================================================//
        // Compute the weight for PDF syst   //
        //===================================//
        //-- get the pdgId of the two colliding partons 
        double wPdf(1);
        if (pdfSet != "") {
            int id1 = pdfInfo_->at(2);
            int id2 = pdfInfo_->at(3);
            if (id1 == 21) id1 = 0;
            if (id2 == 21) id2 = 0;

            LHAPDF::usePDFMember(2, 0);
            double pdf1 = LHAPDF::xfx(1, pdfInfo_->at(2), pdfInfo_->at(4), id1);
            double pdf2 = LHAPDF::xfx(1, pdfInfo_->at(3), pdfInfo_->at(4), id2);
            double pdf01 = LHAPDF::xfx(2, pdfInfo_->at(2), pdfInfo_->at(4), id1);
            double pdf02 = LHAPDF::xfx(2, pdfInfo_->at(3), pdfInfo_->at(4), id2);

            if (pdfInfo_->at(2) * pdfInfo_->at(3) > 0) {
                wPdf = pdf1 * pdf2;
                if (pdf01*pdf02 <= 0 || pdf1*pdf2 <= 0) {
                    wPdf = 1;
                }
                else {
                    wPdf /= (pdf01 * pdf02);
                }
            }
        }
        //==========================================================================================================//

        //--- There is no pile-up so no need to reweight for that ---
        genWeight = weight * wPdf;
        double genWeightBackup(genWeight);
        TotalGenWeight += genWeightBackup;
        //=======================================================================================================//


        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        //=======================================================================================================//
        //         Retrieving leptons           //
        //====================================//
        bool doMuons(leptonFlavor == "Muons" || doW || doTT);
        bool doElectrons(leptonFlavor == "Electrons" || doW || doTT);
        bool passesLeptonCut(0);
        unsigned short nTotLeptons(0), nLeptons(0), nMuons(0), nElectrons(0);
        vector<leptonStruct> leptons, muons, electrons, mets, allLooseLeptons, selLeptons;
        TLorentzVector lep1, lep2, Z;
        leptonStruct lepton1 = {0, 0, 0, 0, 0, 0, 0};
        leptonStruct lepton2 = {0, 0, 0, 0, 0, 0, 0};
        double METphi(0), METpt(0), MT(0);
        //double METpt(0), MT(0);
        int whichMet(2); //  0 - pfMETPFlow, 1 - pfMet, 2 - pfType1CorrectedMet, 3 - pfType1p2CorrectedMet

        int sumLepCharge(1000);
        //APICHART if (hasRecoInfo && doTT && METPt->at(whichMet) < METcut) continue;
        if (hasRecoInfo) {

            bool eventTrigger = false;
            //--- DO MUONS ---
            if (doMuons){
                nTotLeptons = MuEta->size();
                // test if any of the leptons is atached to trigger
                // if we don't really care to match both leptons to trigger
                // we also have event trigger variables --> we should at least match one of the leptons to trigger
              /*  for (unsigned short i(0); i < nTotLeptons; i++) {
                    int whichTrigger(MuTrig->at(i));
                    if (energy == "7TeV" && whichTrigger > 0) eventTrigger = true;
                    //if (energy == "8TeV" && (whichTrigger % 2) == 1 && doW) eventTrigger = true;
                    if (energy == "8TeV" && (whichTrigger & 0x1) && doW) eventTrigger = true;
                    if (energy == "8TeV" && doTT && whichTrigger >= 16) eventTrigger = true; // for TT background
                }*/ //KO

                for (unsigned short i(0); i < nTotLeptons; i++) {
                    if (doMer) merUncer = Rand_MER_Gen->Gaus(0, (MuPt->at(i) * 0.006));
                    leptonStruct mu = {(MuPt->at(i) * muScale) + merUncer, MuEta->at(i), MuPhi->at(i), MuE->at(i), MuCh->at(i), MuPfIso->at(i), 0};
                    //KOint whichTrigger(MuTrig->at(i));
                    bool muPassesPtCut(( (doZ || doTT) && mu.pt >= 20.) || (doW && mu.pt >= 25.));
                    //bool muPassesPtCut(( (doZ || doTT) && mu.pt >= 20.) || (doW && mu.pt >= 30.));
                    
                    bool muPassesEtaLooseCut(fabs(mu.eta) <= 2.4);
                    bool muPassesEtaCut( ((doZ || doTT) && muPassesEtaLooseCut) || (doW && fabs(mu.eta) <= 2.1) );
                    
                    // We use Tight muon, only for tight the patMuonCombId_ is odd
                    //kobool muPassesIdCut(int(MuIdTight->at(i)) % 2 == 1); // this is for tight ID --> odd number
                    bool muPassesIdCut((MuIdTight->at(i) & 1)); // this is for tight ID --> odd number
                    //bool muPassesIdCut(int(patMuonCombId_->at(i)) >= 1 ); // this is for Loose ID

                    //kobool muPassesDxyCut(MuDxy->at(i) < 0.2);
                    bool muPassesIsoCut((!doW && MuPfIso->at(i) < 0.2) || (doW && MuPfIso->at(i) < 0.12));  
                    bool muPassesQCDIsoCut(doW && MuPfIso->at(i) >= 0.2); // use 0.12 if you want to cover full Iso space
                    //KObool muPassesEMuAndWJetsTrig( whichTrigger == 1 || whichTrigger == 16 || whichTrigger == 17 || whichTrigger == 32 || whichTrigger == 33 || whichTrigger == 48 || whichTrigger ==  49  ) ;
                    //KObool muPassesAnyTrig((doZ && ((energy == "7TeV" && whichTrigger > 0) || (energy == "8TeV" && whichTrigger > 7 && !muPassesEMuAndWJetsTrig))) ||
                           //KO (doW && whichTrigger % 2 == 1) || (doTT && whichTrigger >= 16)); // 8TeV comment: Mu17Mu8Tk = 4; Mu17Mu8 = 8 
                    /// for files obtained form bugra
                   //KO if (fileName.find("DYJets_Sherpa_UNFOLDING_dR_5311") != string::npos && whichTrigger > 0) muPassesAnyTrig = 1; 

                    // select the good muons only
                    //-- no Isolation Cut
                    //koif (!doTT && muPassesEtaLooseCut && MuIdTight->at(i) > 0 && mu.pt >= 15) muons.push_back(mu);
                    if (!doTT && muPassesEtaLooseCut && mu.pt >= 15) muons.push_back(mu);
                    
                    //koif (doW && fabs(mu.eta) > 2.1) muPassesEtaCut = false;
                    //koif (doTT && fabs(mu.eta) > 2.4) muPassesEtaCut = false;
                    
                    //KOif (muPassesPtCut && muPassesEtaCut && muPassesIdCut && muPassesDxyCut && (!useTriggerCorrection || muPassesAnyTrig || eventTrigger)){
                    if (muPassesPtCut && muPassesEtaCut && muPassesIdCut){
                        // fill isolation histograms for control    
                        //koMuDetIsoRhoCorr->Fill(MuPfIso->at(i), weight);
                        //koMuPFIsoDBetaCorr->Fill(MuPfIso->at(i), weight);
                        //-- isolation Cut
                        if (doQCD > 1 && muPassesQCDIsoCut && leptonFlavor != "SingleElectron") leptons.push_back(mu);
                        if (muPassesIsoCut){  
                            if (doQCD < 2 && leptonFlavor != "SingleElectron") leptons.push_back(mu); 
                            if (doTT && fabs(mu.eta) < 2.4) muons.push_back(mu); 
                        } 
                    }
                }//End of loop over all the muons
            }
            nMuons = muons.size();
 if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
            //------ DO ELECTRONS -------
            if (doElectrons) {
                nTotLeptons = 0;
                nTotLeptons = ElEta->size();
                // if we don't really care to match both leptons to trigger
                if (doW) eventTrigger = false;
               /*//KO for (unsigned short i(0); i < nTotLeptons; i++){
                //KO    int whichTrigger(ElTrig->at(i));
                  //KO  if (energy == "7TeV" && whichTrigger > 0) eventTrigger = true; ///
                   //KO if (energy == "8TeV" && (whichTrigger & 0x1) && doW) eventTrigger = true; ///
                    if (energy == "8TeV" && doTT && whichTrigger >= 16) eventTrigger = true; // for TT background ///
                }*///KO
                
                for (unsigned short i(0); i < nTotLeptons; i++){
                    leptonStruct ele = {ElPt->at(i), ElEta->at(i), ElPhi->at(i), ElE->at(i),  ElCh->at(i), 0., ElEtaSc->at(i)};
                   //KO int whichTrigger(ElTrig->at(i));
                    bool elePassesPtCut( ( !doW && ele.pt >= 20.)  || ( doW && ele.pt >= 30.));
                    
                    bool elePassesEtaLooseCut(fabs(ElEtaSc->at(i)) <= 1.4442 || (fabs(ElEtaSc->at(i)) >= 1.566 && fabs(ElEtaSc->at(i)) <= 2.4));
                    bool elePassesEtaCut( ((doZ || doTT) && elePassesEtaLooseCut) || (doW && elePassesEtaLooseCut && fabs(ElEtaSc->at(i)) <= 2.1) );
                    
                    // We use medium electron id
                    //bool elePassesIdCut(int(ElId->at(i)) >= 4); /// >=4 is medium ID; >=2 is Loose ID
                    bool elePassesIdCut((ElId->at(i) & 1<<2)); ///// APICHART: medium ID
                    
                    bool elePassesIsoCut(ElPfIsoRho->at(i) < 0.15 );

                    //KObool elePassesEMuAndWJetsTrig(whichTrigger == 1 || whichTrigger == 16 || whichTrigger == 17 || whichTrigger == 32 || whichTrigger == 33 || whichTrigger == 48 || whichTrigger == 49 ) ;
                   //KO bool elePassesAnyTrig(  (doZ && (whichTrigger >= 2 && !elePassesEMuAndWJetsTrig )) || ( doTT && whichTrigger >= 16 ) || ( doW && whichTrigger % 2 == 1) );
                    //KOif ( DEBUG ) cout << EvtNum << "  lepton loop: "<<elePassesAnyTrig <<"   " << ele.pt <<"   " << ele.eta <<"  " <<"  " << ElE->at(i) <<"  " <<elePassesIdCut<<"  SIZE  " << nTotLeptons <<  endl;
                    // elePassesAnyTrig = true ;

                    // select the good electrons only
                    // APICHART: Loose ID
                    if (!doTT && elePassesEtaLooseCut && (ElId->at(i) & 1<<1) && ele.pt >= 15. && ElPfIsoRho->at(i) < 0.2 )  electrons.push_back(ele); /// DO I WANT THIS !!!!!!
                    
                    if (doW && fabs(ElEtaSc->at(i)) > 2.1) elePassesEtaCut = false ;
                    
                    //KOif (elePassesPtCut && elePassesEtaCut && elePassesIdCut && (!useTriggerCorrection || elePassesAnyTrig || eventTrigger)){
                    if (elePassesPtCut && elePassesEtaCut && elePassesIdCut){
                        //-- isolation Cut
                        if (doQCD > 1  && !elePassesIsoCut && leptonFlavor != "SingleMuon") leptons.push_back(ele);
                        if ( elePassesIsoCut ) {
                                if (doQCD < 2 && leptonFlavor != "SingleMuon") leptons.push_back(ele);
                                if (doTT) electrons.push_back(ele);
                            }
                    }
                }//End of loop over all the electrons
            }

            nElectrons = electrons.size();
            int diLepCharge(0);
            nLeptons = leptons.size();
            vector<leptonStruct> tempVec;
            for ( int iLep = 0 ; iLep < nLeptons ; iLep++){
                tempVec.push_back(leptons[iLep]);
            }
            selLeptons = tempVec ;
            if (nLeptons == 2) diLepCharge = abs(leptons[0].charge) + abs(leptons[1].charge);

 if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
            //=======================================================================================================//
            //          Retrieving MET             //
            //====================================//
            //APICHART
            //cout << "METPt->size() " << METPt->size() << endl;
            //APICHART
            //cout << "METPx->size() " << METPx->size() << endl;
            
            if (doW && !(METPt->size() > 0)) continue;
            
                //KOMETphi = METPhi->at(whichMet);
            METpt = METPt->at(whichMet);

            
            //cout << " jentry: " << jentry << " nMuons: " << nMuons << " nLeptons: " << nLeptons << " nElectrons: " << nElectrons <<  endl;
            if (doW && ((leptonFlavor == "SingleMuon" && nMuons == 1 && nLeptons == 1 && nElectrons == 0) || (leptonFlavor == "SingleElectron" && nMuons == 0 && nLeptons == 1 && nElectrons == 1))) {
                
                lepton1 = leptons[0];
                
                
                TLorentzVector tmpVecMet;
                tmpVecMet.SetPxPyPzE(METPx->at(whichMet), METPy->at(whichMet), METPz->at(whichMet), METE->at(whichMet));
                METphi = tmpVecMet.Phi();
                
                
                
                // build the TLorentzVectors, the Z candidate and the kinematic
                lep1.SetPtEtaPhiM(lepton1.pt, lepton1.eta, lepton1.phi, leptonMass);
                lep2.SetPtEtaPhiM(METpt, 0, METphi, 0);
                
                Z = lep1 + lep2;

                leptonStruct tempMet = {METpt, 0., METphi, METpt, 0, 0, 0};
                lepton2 = tempMet;

                MT = sqrt(2 * METpt * lepton1.pt * (1 - cos(METphi - lepton1.phi)));
                //MT = sqrt(2 * METpt * lepton1.pt * (1 - cos(0 - lepton1.phi)));

                /// 2D histograms for ABCD method to extract the QCD abckround ?
                /* APICHART
                fullMET->Fill(METpt, weight);
                fullMET_pfMETPFlow->Fill(METPt->at(0), weight);
                fullMET_pfMet->Fill(METPt->at(1), weight);
                fullMET_pfType1CorrectedMet->Fill(METPt->at(2), weight);
                fullMET_pfType1p2CorrectedMet->Fill(METPt->at(3), weight);
                fullMT->Fill(MT, weight);
                METvslepIso->Fill(METpt, lepton1.iso, weight);
                MTvslepIso->Fill(MT, lepton1.iso, weight);
                 APICHART */
                // APICHART double effWeight = 1.;
                
                if (METpt >= METcut && (((doQCD % 2) == 0 && MT >= MTCut) || ((doQCD % 2) == 1 && MT < MTCut))) {
                    passesLeptonCut = true;
                    if (fabs(scale) == 0) nEventsWithTwoGoodLeptons++;
                    // correct for identification and isolation efficiencies if required by useEfficiencyCorrection
                    // apply scale factors only on MC
                    
                    /* APICHART
                    if (fabs(scale) == 0 && useEfficiencyCorrection) {
                        if (leptonFlavor == "SingleMuon") {
                            effWeight = LeptID.getEfficiency(lepton1.pt, fabs(lepton1.eta), sysLepSF);
                            effWeight *= LeptIso.getEfficiency(lepton1.pt, fabs(lepton1.eta), sysLepSF);
                            if (useTriggerCorrection) effWeight *= LeptTrig.getEfficiency(fabs(lepton1.pt), fabs(lepton1.eta), sysLepSF);
                        }
                        else if (leptonFlavor == "SingleElectron") {
                            effWeight *= LeptID.getEfficiency(lepton1.pt, fabs(lepton1.eta), sysLepSF);
                            effWeight *= Ele_Rec.getEfficiency(lepton1.pt, fabs(lepton1.scEta), sysLepSF);
                            
                        }
                    }
                    // APICHART */
                }
                // APICHART if (isData) weight /= effWeight;
                // APICHART else weight *= effWeight;
                
            }  // END IF RECO FOR MET AND LEPTONS
            //=======================================================================================================//
        }// end has reco info

 if (DEBUG) cout << "Stop after line " << __LINE__ << endl;

        //=======================================================================================================//
        //       Retrieving gen leptons        //
        //====================================//
        bool passesGenLeptonCut(0);
        unsigned short nTotGenLeptons(0), nGenLeptons(0), nTotGenPhotons(0);
        vector<leptonStruct> genLeptons;
        vector<int> usedGenPho;
        TLorentzVector genLep1, genLep2, genZ;
        leptonStruct genLepton1, genLepton2;
        int countTauS3 = 0;

        // to use the TOP PAG TTBAR reweighting recommendation
        // line below is to check contribution from mainy tau decay : use passesLeptonCut = 0  only if you want to have RECO events that originate from tau ; countTauS3 is used in passesGenLeptonCut
        // I also put top quarks in this collection
        // ...DELETED
        /// end top reweighting

        countTauS3 = 0; 
        if (hasGenInfo) {
            if (hasRecoInfo) countTauS3 = 2;
            if (hasRecoInfo && doW) countTauS3 = 1;
            nTotGenPhotons = genPhoEta_->size();
            nTotGenLeptons = genLepEta_->size();
            //-- retriveing generated leptons with status 1
            for (unsigned short i(0); i < nTotGenLeptons; i++) {
                // line below is to check contribution from mainy tau decay : use passesLeptonCut = 0  only if you want to have RECO events that originate from tau ; countTauS3 is used in passesGenLeptonCut
                bool lepSelector( 
                        (doZ && abs(genLepId_->at(i)) == LeptonID) || 
                        (doW && (abs(genLepId_->at(i)) == LeptonID || abs(genLepId_->at(i)) == 12 || abs(genLepId_->at(i)) == 14)));
                
                
                // following two lines should give the same result
                if (genLepSt_->at(i) == 3 && abs(genLepId_->at(i)) != LeptonID && (abs(genLepId_->at(i)) == 15 || abs(genLepId_->at(i)) == 13 || abs(genLepId_->at(i)) == 11)) countTauS3++;
                if (genLepSt_->at(i) == 3 && abs(genLepId_->at(i)) == LeptonID ) countTauS3--;
                
                
                if (!lepSelector) continue ;
                double charge(genLepQ_->at(i)); 
                if (abs(genLepId_->at(i)) == 12 || abs(genLepId_->at(i)) == 14 || abs(genLepId_->at(i)) == 16) charge = 0.;
                leptonStruct genLep = {genLepPt_->at(i), genLepEta_->at(i), genLepPhi_->at(i), genLepE_->at(i), charge, 0., 0.};
                leptonStruct genLepNoFSR = {genLepPt_->at(i), genLepEta_->at(i), genLepPhi_->at(i), genLepE_->at(i), charge, 0., 0. };
                
                //-- dress the leptons with photon (cone size = 0.1). Only for status 1 leptons (after FSR)
                if ( ( genLepSt_->at(i) == 1 && lepSelector && abs(genLepId_->at(i)) == LeptonID) || ( doW && charge == 0 ) ){
                    // only charged lepton(s) will be dressed
                    if( fabs(genLep.charge) > 0 ){
                        TLorentzVector tmpGenLep;
                        tmpGenLep.SetPtEtaPhiM(genLep.pt, genLep.eta, genLep.phi, leptonMass);
                        // loop over all photons
                        for (unsigned short j(0); j < nTotGenPhotons; j++){
                            TLorentzVector tmpGenPho;
                            tmpGenPho.SetPtEtaPhiM(genPhoPt_->at(j), genPhoEta_->at(j), genPhoPhi_->at(j), 0.);
                            int used(0);
                            for (unsigned short k(0); k < usedGenPho.size(); k++){
                                if (j == usedGenPho[k]) used = 1;
                            }
                            if (deltaR(tmpGenPho.Phi(), tmpGenPho.Eta(), genLepNoFSR.phi, genLepNoFSR.eta) <= 0.1 && !used){
                                tmpGenLep += tmpGenPho;
                                usedGenPho.push_back(j);
                            }
                        }
                        genLep.pt = tmpGenLep.Pt();
                        genLep.eta = tmpGenLep.Eta();
                        genLep.phi = tmpGenLep.Phi();
                        genLep.energy = tmpGenLep.E();
                    }
                    
                    //-- store lepton in the collection
                    if (doZ && genLep.pt >= 20 && fabs(genLep.eta) <= 2.4 && fabs(genLep.charge) > 0){
                        genLeptons.push_back(genLep);
                    }
                    
                    //if (doW && ((fabs(genLep.charge) > 0 && genLep.pt >= 25 && fabs(genLep.eta) <= 2.1)||(fabs(genLep.charge) == 0))){
                    if (doW && ((fabs(genLep.charge) > 0 && genLep.pt >= 25 && fabs(genLep.eta) <= 2.1) || (fabs(genLep.charge) == 0 && genLep.pt >= METcut))){
                        genLeptons.push_back(genLep); 
                    }
                }
            }
            nGenLeptons = genLeptons.size();

            if (countTauS3 == 0 && fileName.find("UNFOLDING") != string::npos){
                partonsN->Fill(nup_-5);
                partonsNWeighted->Fill(nup_-5, genWeight);
            }
            
            /// --- if there are taus, but we do not run on the Tau file, thus we run on the WJets file,
            //    then we don't count the event at reco.
            // APICHART if (useEfficiencyCorrection && countTauS3 > 0 && fileName.find("Tau") == string::npos  ) passesLeptonCut = 0 ;
            
            
            if ( ( fileName.find("Tau") == string::npos &&  countTauS3 > 0  ) || ( fileName.find("Tau") != string::npos &&  countTauS3 == 0) ){
                ZMassAllPassLep->Fill(Z.M(),weight);
                AllPassLepID->Fill(sumLepCharge,weight);
                if (Z.M() > 50 )  AllPassWithMassCutLepID->Fill(sumLepCharge,weight);
            }
            if (passesLeptonCut)      AllPassWithMassCutLepIDCharge->Fill(sumLepCharge,weight);
            
if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
            //-- determine if the event passes the leptons requirements
            if (nGenLeptons >= 2){
                // sort leptons by descending pt
                sort(genLeptons.begin(), genLeptons.end(), LepDescendingOrder);
                genLepton1 = genLeptons[0];
                genLepton2 = genLeptons[1];
                
                if (doW){
                    double genMT;
                    if (abs(genLeptons[0].charge) > 0 && genLeptons[1].charge == 0){
                        genLepton1 = genLeptons[0];
                        genLepton2 = genLeptons[1];
                        genMT = sqrt(2 *genLepton2.pt * genLepton1.pt * (1 - cos(genLepton2.phi - genLepton1.phi)));
                    }
                    else if (abs(genLeptons[1].charge) > 0 && genLeptons[0].charge == 0){
                        genLepton1 = genLeptons[1];
                        genLepton2 = genLeptons[0];
                        genMT = sqrt(2 *genLepton2.pt * genLepton1.pt * (1 - cos(genLepton2.phi - genLepton1.phi)));
                    }
                    else genMT = -99.0;
                    
                    if (genMT >= MTCut && genLepton2.pt >= METcut) passesGenLeptonCut = 1;
                }

                //----- For Z+jets -------
                // select the first two leptons with opposite charge
                if (doZ && genLepton1.charge*genLepton2.charge > 0 && nGenLeptons > 2) {
                    genLepton2 = genLeptons[2];
                }
                // build the TLorentzVectors, the Z candidate and the kinematic
                genLep1.SetPtEtaPhiE(genLepton1.pt, genLepton1.eta, genLepton1.phi, genLepton1.energy);
                genLep2.SetPtEtaPhiE(genLepton2.pt, genLepton2.eta, genLepton2.phi, genLepton2.energy);
                genZ = genLep1 + genLep2;
                // apply charge, mass and eta cut
                if (doZ && genLepton1.charge*genLepton2.charge < 0 
                        && genZ.M() > ZMCutLow && genZ.M() < ZMCutHigh 
                        && genZ.Eta()*100 > ZEtaCutMin && genZ.Eta()*100 < ZEtaCutMax && genZ.Pt()>= ZPtCutMin) {
                    
                    passesGenLeptonCut = 1;
                }
                //----- End For Z+jets -------
                
                //--- if there are taus we don't want the gen level
                // APICHART if (countTauS3 > 0) passesGenLeptonCut = 0;
            }
        }
        if (passesGenLeptonCut) {
            TotalGenWeightPassGEN += genWeightBackup; 
            TotalGenWeightPassGENPU += weight;
            partonsNAfterGenCut->Fill(nup_ - 5);
            partonsNAfterGenCutWeighted->Fill(nup_ - 5, genWeight);
        }
        //=======================================================================================================//



        //=======================================================================================================//
        //          Retrieving jets           //
        //====================================//
        bool passesJetCut(1), passesEWKJetPt(0), passesEWKJetFwdEta(0);
        unsigned short nGoodJets(0), nTotJets(0), nJetsAdd(0);
        double jetsHT(0);
        //KOdouble METscale(0.);
        //KOdouble XMETscale(0.), YMETscale(0.);            // for calculating METscale
        //KOdouble TempMETpt(0.), TempMETphi(0.), XMETpt(0.), YMETpt(0.) ; // for calculating METscale

        vector<jetStruct> jets, jetsAdditional;
        TLorentzVector leadJ, secondJ, jet1Plus2, jet1Minus2;
        
        //*************************************** begin edit *************************************************************//
        TLorentzVector newLeadJ, newSecondJ, newThirdJ, newFourthJ;
        double ForwardJetRapidity(0), BackwardJetRapidity(0);
        vector<TLorentzVector> vJetYOrdered;
        //**************************************** end edit **************************************************************//
        
        int countBJets = 0;
        int countWbBjets = 0; // Wb study

        if (hasRecoInfo) {
            int countNJetsVSBeta[10] = {0};
            nTotJets = JetAk04Eta->size();
            
if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
            
            //--- loop over all the jets ----------
            for (unsigned short i(0); i < nTotJets; i++) {
                double jetPtTemp(0.); // for calculating METscale
                bool passBJets(0);
                // APICHART if (JetAk04BTagCsv->at(i) >= 0.679) passBJets = true;
                if (JetAk04BTagCsv->at(i) >= 0.890) passBJets = true;
                
                

                //************************* B-tag Veto Correcction *******************************//
       /*         float this_rand = RandGen->Rndm(); // Get a random number.
                float pt= JetAk04Pt->at(i);
                float eta= JetAk04Eta->at(i);
                float x = 0.679;     ///discrim_cut;
                // --------- MC-only
                if (isData == false){
                    
                    bool passBJets_SFB_sys_up = passBJets;     // Initialize the systematic_up as the central value
                    bool passBJets_SFB_sys_down = passBJets; // Initialize the systematic_down as the central value
                    
                    int jetflavour= JetAk04PartFlav->at(i);
                if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
                    if (abs(jetflavour)==5){
                        float effb = -1.73338329789*x*x*x*x +  1.26161794785*x*x*x +  0.784721653518*x*x +  -1.03328577451*x +  1.04305075822;
                        
                        float SFb = (0.938887+(0.00017124*pt))+(-2.76366e-07*(pt*pt));
                        if (pt < 20.) SFb = (0.938887+(0.00017124*20.))+(-2.76366e-07*(20.*20.));
                        if (pt > 800.) SFb = (0.938887+(0.00017124*800.))+(-2.76366e-07*(800.*800.));
                        
                        float SFb_error = 0.0;
                        if (pt < 20.)                SFb_error = 0.0415707*2.;
                        if (pt >= 20. && pt < 30.)   SFb_error = 0.0415707;
                        if (pt >= 30. && pt < 40.)   SFb_error = 0.0204209;
                        if (pt >= 40. && pt < 50.)   SFb_error = 0.0223227;
                        if (pt >= 50. && pt < 60.)   SFb_error = 0.0206655;
                        if (pt >= 60. && pt < 70.)   SFb_error = 0.0199325;
                        if (pt >= 70. && pt < 80.)   SFb_error = 0.0174121;
                        if (pt >= 80. && pt < 100.)  SFb_error = 0.0202332;
                        if (pt >= 100. && pt < 120.) SFb_error = 0.0182446;
                        if (pt >= 120. && pt < 160.) SFb_error = 0.0159777;
                        if (pt >= 160. && pt < 210.) SFb_error = 0.0218531;
                        if (pt >= 210. && pt < 260.) SFb_error = 0.0204688;
                        if (pt >= 260. && pt < 320.) SFb_error = 0.0265191;
                        if (pt >= 320. && pt < 400.) SFb_error = 0.0313175;
                        if (pt >= 400. && pt < 500.) SFb_error = 0.0415417;
                        if (pt >= 500. && pt < 600.) SFb_error = 0.0740446;
                        if (pt >= 600. && pt < 800.) SFb_error = 0.0596716;
                        if (pt >= 800.)              SFb_error = 0.0596716*2.;
                        
                        float SFb_up = SFb + SFb_error;
                        float SFb_down = SFb - SFb_error;
                        
                        // F values for rand comparison
                        float f = 0.0;
                        float f_up = 0.0;
                        float f_down = 0.0;
                        
                        if (SFb <1.0) f = (1.0 - SFb);
                        if (SFb_up <1.0) f_up = (1.0 - SFb_up);
                        if (SFb_down <1.0) f_down = (1.0 - SFb_down);
                        
                        if (SFb > 1.0) f = (1.0 - SFb)/(1.0 - 1.0/effb);
                        if (SFb_up > 1.0) f_up = (1.0 - SFb_up)/(1.0 - 1.0/effb);
                        if (SFb_down > 1.0) f_down = (1.0 - SFb_down)/(1.0 - 1.0/effb);
                        
                        passBJets_SFB_sys_up = passBJets;     // Initialize the systematic_up as the central value
                        passBJets_SFB_sys_down = passBJets; // Initialize the systematic_down as the central value
                        
                        // Untag a tagged jet
                        if ((passBJets==true) && (SFb<1.0) && (this_rand < f)) passBJets = false; // for central value
                        if ((passBJets_SFB_sys_up==true)   && (SFb_up<1.0) && (this_rand < f_up))   passBJets_SFB_sys_up = false; // for systematic_up
                        if ((passBJets_SFB_sys_down==true) && (SFb_down<1.0) && (this_rand < f_down)) passBJets_SFB_sys_down = false; // for sytematic_down
                        
                        
                        // Tag an untagged jet
                        if ((passBJets==false) && (SFb>1.0) && (this_rand < f)) passBJets = true; // for central value
                        if ((passBJets_SFB_sys_up==false)   && (SFb_up>1.0) && (this_rand < f_up))   passBJets_SFB_sys_up = true; // for systematic_up
                        if ((passBJets_SFB_sys_down==false) && (SFb_down>1.0) && (this_rand < f_down)) passBJets_SFB_sys_down = true; // for sytematic_down
                        
                    }


 if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
                    // ---------------- For Real C-jets--------------- //
                    if (abs(jetflavour)==4){
                        float effc = -1.5734604211*x*x*x*x +  1.52798999269*x*x*x +  0.866697059943*x*x +  -1.66657942274*x +  0.780639301724;
                        
                        float SFc = (0.938887+(0.00017124*pt))+(-2.76366e-07*(pt*pt));
                        if (pt < 20.) SFc = (0.938887+(0.00017124*20.))+(-2.76366e-07*(20.*20.));
                        if (pt > 800.) SFc = (0.938887+(0.00017124*800.))+(-2.76366e-07*(800.*800.));
                        
                        float SFc_error = 0.0;
                        if (pt < 20.)                SFc_error = 0.0415707*4.;
                        if (pt >= 20. && pt < 30.)   SFc_error = 0.0415707*2.;
                        if (pt >= 30. && pt < 40.)   SFc_error = 0.0204209*2.;
                        if (pt >= 40. && pt < 50.)   SFc_error = 0.0223227*2.;
                        if (pt >= 50. && pt < 60.)   SFc_error = 0.0206655*2.;
                        if (pt >= 60. && pt < 70.)   SFc_error = 0.0199325*2.;
                        if (pt >= 70. && pt < 80.)   SFc_error = 0.0174121*2.;
                        if (pt >= 80. && pt < 100.)  SFc_error = 0.0202332*2.;
                        if (pt >= 100. && pt < 120.) SFc_error = 0.0182446*2.;
                        if (pt >= 120. && pt < 160.) SFc_error = 0.0159777*2.;
                        if (pt >= 160. && pt < 210.) SFc_error = 0.0218531*2.;
                        if (pt >= 210. && pt < 260.) SFc_error = 0.0204688*2.;
                        if (pt >= 260. && pt < 320.) SFc_error = 0.0265191*2.;
                        if (pt >= 320. && pt < 400.) SFc_error = 0.0313175*2.;
                        if (pt >= 400. && pt < 500.) SFc_error = 0.0415417*2.;
                        if (pt >= 500. && pt < 600.) SFc_error = 0.0740446*2.;
                        if (pt >= 600. && pt < 800.) SFc_error = 0.0596716*2.;
                        if (pt >= 800.)              SFc_error = 0.0596716*4.;
                        
                        float SFc_up = SFc + SFc_error;
                        float SFc_down = SFc - SFc_error;
                        
                        // F values for rand comparison
                        float f = 0.0;
                        float f_up = 0.0;
                        float f_down = 0.0;
                        
                        if (SFc <1.0) f = (1.0 - SFc);
                        if (SFc_up <1.0) f_up = (1.0 - SFc_up);
                        if (SFc_down <1.0) f_down = (1.0 - SFc_down);
                        
                        if (SFc > 1.0) f = (1.0 - SFc)/(1.0 - 1.0/effc);
                        if (SFc_up > 1.0) f_up = (1.0 - SFc_up)/(1.0 - 1.0/effc);
                        if (SFc_down > 1.0) f_down = (1.0 - SFc_down)/(1.0 - 1.0/effc);
                        
                        passBJets_SFB_sys_up = passBJets;     // Initialize the systematic_up as the central value
                        passBJets_SFB_sys_down = passBJets; // Initialize the systematic_down as the central value
                        
                        // Untag a tagged jet
                        if ((passBJets==true) && (SFc<1.0) && (this_rand < f)) passBJets = false; // for central value
                        if ((passBJets_SFB_sys_up==true)   && (SFc_up<1.0) && (this_rand < f_up))   passBJets_SFB_sys_up = false; // for systematic_up
                        if ((passBJets_SFB_sys_down==true) && (SFc_down<1.0) && (this_rand < f_down)) passBJets_SFB_sys_down = false; // for sytematic_down
                        
                        // Tag an untagged jet
                        if ((passBJets==false) && (SFc>1.0) && (this_rand < f)) passBJets = true; // for central value
                        if ((passBJets_SFB_sys_up==false)   && (SFc_up>1.0) && (this_rand < f_up))   passBJets_SFB_sys_up = true; // for systematic_up
                        if ((passBJets_SFB_sys_down==false) && (SFc_down>1.0) && (this_rand < f_down)) passBJets_SFB_sys_down = true; // for sytematic_down
                        
                    }
 if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
                    // ---------------- For REAL Light-jets --------------- //
                    if (abs(jetflavour)<4){
                        float SFlight=1.0;
                        float SFlight_up=1.0;
                        float SFlight_down=1.0;
                        float eff_l = 0.0;
                        float pt_temp(0.0), xpt(0.0);
                        
                        xpt = pt;
                        if (pt < 20. ) xpt = 20.;
                        if (pt > 670.) xpt = 670.;
                        
                        if ((fabs(eta)>=0.0) && (fabs(eta)<=0.8)){
                            eff_l = ((0.00967751+(2.54564e-05*xpt))+(-6.92256e-10*(xpt*xpt)));
                            
                            pt_temp = pt;
                            if (pt < 20.) pt = 20.;
                            if (pt > 1000.) pt = 1000.;
                            
                            SFlight = (((1.07541+(0.00231827*pt))+(-4.74249e-06*(pt*pt)))+(2.70862e-09*(pt*(pt*pt))));
                            SFlight_up = (((1.18638+(0.00314148*pt))+(-6.68993e-06*(pt*pt)))+(3.89288e-09*(pt*(pt*pt))));
                            SFlight_down = (((0.964527+(0.00149055*pt))+(-2.78338e-06*(pt*pt)))+(1.51771e-09*(pt*(pt*pt))));
                            
                            if (pt_temp < 20.|| pt_temp > 1000.) {
                                SFlight_up =  SFlight + 2.0*fabs(SFlight_up - SFlight);
                                SFlight_down = SFlight - 2.0*fabs(SFlight - SFlight_down);
                            }
                            pt = pt_temp;
                        }
                        if ((fabs(eta)>0.8) && (fabs(eta)<=1.6)){
                            eff_l = ((0.00974141+(5.09503e-05*xpt))+(2.0641e-08*(xpt*xpt)));
                            
                            pt_temp = pt;
                            if (pt < 20.) pt = 20.;
                            if (pt > 1000.) pt = 1000.;
                            
                            SFlight = (((1.05613+(0.00114031*pt))+(-2.56066e-06*(pt*pt)))+(1.67792e-09*(pt*(pt*pt))));
                            SFlight_up = (((1.16624+(0.00151884*pt))+(-3.59041e-06*(pt*pt)))+(2.38681e-09*(pt*(pt*pt))));
                            SFlight_down = (((0.946051+(0.000759584*pt))+(-1.52491e-06*(pt*pt)))+(9.65822e-10*(pt*(pt*pt))));
                            
                            if (pt_temp < 20.|| pt_temp > 1000.) {
                                SFlight_up =  SFlight + 2.0*fabs(SFlight_up - SFlight);
                                SFlight_down = SFlight - 2.0*fabs(SFlight - SFlight_down);
                            }
                            pt = pt_temp;
                        }
                        if ((fabs(eta)>1.6) && (fabs(eta)<=2.4)){
                            eff_l = ((0.013595+(0.000104538*xpt))+(-1.36087e-08*(xpt*xpt)));
                            
                            pt_temp = pt;
                            if (pt < 20.) pt = 20.;
                            if (pt > 850.) pt = 850.;
                            
                            SFlight = (((1.05625+(0.000487231*pt))+(-2.22792e-06*(pt*pt)))+(1.70262e-09*(pt*(pt*pt))));
                            SFlight_up = (((1.15575+(0.000693344*pt))+(-3.02661e-06*(pt*pt)))+(2.39752e-09*(pt*(pt*pt))));
                            SFlight_down = (((0.956736+(0.000280197*pt))+(-1.42739e-06*(pt*pt)))+(1.0085e-09*(pt*(pt*pt))));
                            
                            if (pt_temp < 20.|| pt_temp > 850.) {
                                SFlight_up =  SFlight + 2.0*fabs(SFlight_up - SFlight);
                                SFlight_down = SFlight - 2.0*fabs(SFlight - SFlight_down);
                            }
                            pt = pt_temp;
                        }
                        if (fabs(eta)>2.4){
                            eff_l = ((0.013595+(0.000104538*xpt))+(-1.36087e-08*(xpt*xpt)));
                        }
                        
                        // F values for rand comparison
                        float f = 0.0;
                        float f_up = 0.0;
                        float f_down = 0.0;
                        
                        if (SFlight <1.0) f = (1.0 - SFlight);
                        if (SFlight_up <1.0) f_up = (1.0 - SFlight_up);
                        if (SFlight_down <1.0) f_down = (1.0 - SFlight_down);
                        
                        if (SFlight > 1.0) f = (1.0 - SFlight)/(1.0 - 1.0/eff_l);
                        if (SFlight_up > 1.0) f_up = (1.0 - SFlight_up)/(1.0 - 1.0/eff_l);
                        if (SFlight_down > 1.0) f_down = (1.0 - SFlight_down)/(1.0 - 1.0/eff_l);
                        
                        passBJets_SFB_sys_up = passBJets;     // Initialize the systematic_up as the central value
                        passBJets_SFB_sys_down = passBJets; // Initialize the systematic_down as the central value
                        
                        // Untag a tagged jet
                        if ((passBJets==true) && (SFlight<1.0) && (this_rand < f)) passBJets = false; // for central value
                        if ((passBJets_SFB_sys_up==true)   && (SFlight_up<1.0) && (this_rand < f_up))   passBJets_SFB_sys_up = false; // for systematic_up
                        if ((passBJets_SFB_sys_down==true) && (SFlight_down<1.0) && (this_rand < f_down)) passBJets_SFB_sys_down = false; // for sytematic_down
                        
                        // Tag an untagged jet
                        if ((passBJets==false) && (SFlight>1.0) && (this_rand < f)) passBJets = true; // for central value
                        if ((passBJets_SFB_sys_up==false)   && (SFlight_up>1.0) && (this_rand < f_up))   passBJets_SFB_sys_up = true; // for systematic_up
                        if ((passBJets_SFB_sys_down==false) && (SFlight_down>1.0) && (this_rand < f_down)) passBJets_SFB_sys_down = true; // for sytematic_down
                    }   ////////flavour lop
                    
                    if (sysBtagSF ==  1) passBJets = passBJets_SFB_sys_up;
                    if (sysBtagSF == -1) passBJets = passBJets_SFB_sys_down;
                    
                    // Wb study
                    if (abs(jetflavour)==5) countWbBjets++ ;
                } *///KO
                // --------- End MC-only
                //************************* End B-tag Veto Correcction *******************************//                
                int jetflavour = int(JetAk04PartFlav->at(i));
               
                jetStruct jet = {JetAk04Pt->at(i), JetAk04Eta->at(i), JetAk04Phi->at(i), JetAk04E->at(i), i, passBJets};

                //-- apply jet energy scale uncertainty (need to change the scale when initiating the object)
                double jetEnergyCorr = 0.; 
                bool jetPassesPtCut(jet.pt >= 10); // for MET uncertainty should the cut be before or aftes adding unc.?????
                jetEnergyCorr = TableJESunc.getEfficiency(jet.pt, jet.eta);

                jetPtTemp = jet.pt; // for calculating METscale
                jet.pt *= (1 + scale * jetEnergyCorr);
                jet.energy *= (1 + scale * jetEnergyCorr);

                bool jetPassesEtaCut((jet.eta >= jetEtaCutMin / 10.) && (jet.eta <= jetEtaCutMax / 10.)); 
                bool jetPassesIdCut(JetAk04Id->at(i) > 0);
                //KObool jetPassesBetaCut(JetAk04JetBeta->at(i) > 0.1 * doPUStudy);
                //KObool jetPassesBetaStarCut(JetAk04JetBetaStar->at(i) < 1);
                double tempMVA = JetAk04PuMva->at(i);
                bool jetPassesMVACut(0);
                if (energy == "7TeV") {
                    jetPassesMVACut = ((tempMVA > -0.9 && jet.pt <= 20) || (jet.pt > 20 && (
                                    (tempMVA > -0.4  && fabs(jet.eta) <= 2.5) ||
                                    (tempMVA > -0.85 && fabs(jet.eta) > 2.5  && fabs(jet.eta) <= 2.75) ||
                                    (tempMVA > -0.7  && fabs(jet.eta) > 2.75 && fabs(jet.eta) <= 3.) ||
                                    (tempMVA > -0.6  && fabs(jet.eta) > 3.   && fabs(jet.eta) <= 5.))));  
                }
                if (energy == "8TeV") {
                    jetPassesMVACut = (
                            (tempMVA > -0.89 && fabs(jet.eta) <= 2.5) || 
                            (tempMVA > -0.77 && fabs(jet.eta) > 2.5  && fabs(jet.eta) <= 2.75) ||
                            (tempMVA > -0.69 && fabs(jet.eta) > 2.75 && fabs(jet.eta) <= 3.) ||
                            (tempMVA > -0.75 && fabs(jet.eta) > 3.   && fabs(jet.eta) <= 5.)) ;  
                    //  new training does not work for 22Jan rereco, we use simple loose PU ID : -1  - does not pass, 1 passes
                    if (tempMVA > 0) jetPassesMVACut = true ;
                    else jetPassesMVACut = false ;
                }
                //bool jetPassesMVACut(patJetPfAk05jetpuMVA_->at(i) >= - 0.4); // -0.4 set for the cut was for 44x training. for 53x chs loose jet id is set to -0.89? 
                //bool jetPassesMVACut(patJetPfAk05jetpuMVA_->at(i) >= - 0.0 && patJetPfAk05jetBZ_->at(i) > 0.05 ); // 0 set for the cut was for 44x training. for 53x chs loose jet id is set to -0.89? 

                bool jetPassesdRCut(1);
                double dPtJetMuon(-100.);
                unsigned short nRemovedLep = min(int(nLeptons), doW ? 1:2);
                
                //if (jentry % 100000 == 0) cout << "--- nLeptons: " << nLeptons << "  nRemovedLep: " << nRemovedLep << " ---" << endl;

if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
                for (unsigned short j(0); j < nRemovedLep; j++) {
                    // determine if passes dRCut
                    if (deltaR(jet.phi, jet.eta, selLeptons[j].phi, selLeptons[j].eta) < 0.5) {
                        if (doDR) jetPassesdRCut = 0;
                       dPtJetMuon = jet.pt-selLeptons[j].pt;
                        if (jet.pt >= jetPtCutMin && jetPassesMVACut && passesLeptonCut && jetPassesEtaCut && jetPassesIdCut) { 
                            ZMass_lowDeltaR->Fill(Z.M(), weight);
                          //  deltaPtjetMu->Fill(dPtJetMuon,weight);
                        }
                    }
                    // for MET scale
                    //KOif (fabs(scale) > 0. && jetPassesPtCut && jetPassesMVACut && jetPassesIdCut){
                        
                        //METscale -= scale * jetEnergyCorr * jetPtTemp ;
                        //-------- my calculation --------
                 //KO       XMETscale += ( scale * jetEnergyCorr * jetPtTemp * cos(jet.phi) ) ;
                  //KO      YMETscale += ( scale * jetEnergyCorr * jetPtTemp * sin(jet.phi) ) ;
                        //-------- end my calculation --------
                        
                   //KO }
                    if (jet.pt >= jetPtCutMin && jetPassesMVACut && passesLeptonCut && jetPassesEtaCut && jetPassesIdCut){
                       //KO deltaRjetMu->Fill(deltaR(jet.phi, jet.eta, selLeptons[j].phi, selLeptons[j].eta), weight);
                    }
                }
                //cout << " jet passes  :" <<jet.pt <<"     " << jetPassesIdCut <<"  " << jetPassesEtaCut <<"   " << jetPassesdRCut <<"   " << jetPassesMVACut << "     " << tempMVA << endl;
                
if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
                if ( jetPassesEtaCut && jetPassesIdCut && jetPassesdRCut) {
                //if ( jetPassesEtaCut && jetPassesIdCut) {
                    if (jetPassesPtCut){
                        if (jet.pt >= jetPtCutMin && passesLeptonCut){
                            /* APICHART
                            Beta->Fill(JetAk04JetBeta->at(i), weight);
                            BetaStar->Fill(JetAk04JetBetaStar->at(i), weight);
                            puMVA->Fill(JetAk04PuMva->at(i), weight);
                            puMVAvsBeta->Fill(JetAk04PuMva->at(i),JetAk04JetBetaStar->at(i), weight);
                             APICHART*/
                        }

                        //if ( fabs(doBJets) > 0 && patJetPfAk05OCSV_->at(i) >=  0.679 )  countBJets++ ;// count BJets, used for BVeto
                        if ( fabs(doBJets) > 0 && passBJets == true) countBJets++ ;// count BJets, used for BVeto
                        
                        // We apply only the pu MVA variable for the time being.
                        // This is the recommended one.
                        // if (jetPassesBetaCut && jetPassesBetaStarCut && jetPassesMVACut) 
                        if ( jet.pt >= 50 ) passesEWKJetPt = true ;
                        if ( fabs(jet.eta) > 2.4 ) passesEWKJetFwdEta = true ;
                        if (doZ){	
                        if (energy == "8TeV" && doPUStudy < 0  && jetPassesMVACut) jets.push_back(jet);
                        //KOif (energy == "8TeV" && doPUStudy >= 0 && jetPassesBetaCut && jetPassesBetaStarCut) jets.push_back(jet);
                        //KOif (energy == "7TeV" && jetPassesBetaCut && jetPassesBetaStarCut) jets.push_back(jet);
                        for ( int k = 0 ; k < 10 ; k++){
                            if ( JetAk04JetBeta->at(i) >= 0.1 * k )  countNJetsVSBeta[k]++;
                        }

                        }
                        else jets.push_back(jet);
                    }
                    if (jet.pt >=  15.){
                        jetsAdditional.push_back(jet);	
                    }
                }
if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
            }
            //--- End of loop over all the jets ---

            for ( int k = 0 ; k < 10 ; k++){
                ZNGoodJetsBeta_Zexc->Fill(countNJetsVSBeta[k], k  , weight);
            }

            nGoodJets = jets.size();
            nJetsAdd = jetsAdditional.size();
            
            // if do JES systematic variations (Data only), we recalculate MET
/* APICHART
            if (doW && fabs(scale) > 0. && ((leptonFlavor == "SingleMuon" && nMuons == 1 && nLeptons == 1 && nElectrons == 0) || (leptonFlavor == "SingleElectron" && nMuons == 0 && nLeptons == 1 && nElectrons == 1))){
            
                passesLeptonCut = false ;
                
                // recalculate METpt and METphi
                    //METphi = patMetPhi_->at(whichMet);
                    //METpt = patMetPt_->at(whichMet) + METscale;
                
                //-------- my calculation --------
               // TempMETphi = METPhi->at(whichMet);
              //KO  TempMETpt = METPt->at(whichMet);
               //KO XMETpt = (TempMETpt * cos(TempMETphi)) - XMETscale;
               //KO YMETpt = (TempMETpt * sin(TempMETphi)) - YMETscale;
                
              //KO  TVector2 METvec;
              //KO  METvec.Set(XMETpt, YMETpt);
              //KO  METpt  = METvec.Mod();
               //KO METphi = METvec.Phi_mpi_pi(METvec.Phi());
                //-------- end my calculation -------- 
                
                //kolep2.SetPtEtaPhiM(METpt, 0., METphi, 0.);
                lep2.SetPtEtaPhiM(METpt, 0., 0., 0.);
                Z = lep1 + lep2;
                //koleptonStruct tempMet = { METpt , 0., METphi, METpt, 0 , 0, 0.};
                leptonStruct tempMet = { METpt , 0., 0., METpt, 0 , 0, 0.};
                lepton2 = tempMet;
                //koMT = sqrt( 2 * METpt * lepton1.pt * (1 - cos(METphi - lepton1.phi)));
                MT = sqrt( 2 * METpt * lepton1.pt * (1 - cos(0 - lepton1.phi)));
                double effWeight = 1.;
                if (METpt >= METcut && (((doQCD % 2) == 0 && MT >= MTCut) || ((doQCD % 2) == 1 && MT < MTCut))) {
                    passesLeptonCut = true; // This is ok since we do JES for Data only. So, we do not overide passesLeptonCut determined by countTauS3 (GEN info).
                    nEventsWithTwoGoodLeptons++;
                    // correct for identification and isolation efficiencies are not needed since this is Data
                    // So, lines below are not needed
                    
                }
            }
APICHART */
            
            // line below to test reco events that originate from TAU
            
/* APICHART
            if (fileName.find("Tau") != string::npos && countTauS3 == 0 && hasGenInfo ){
                passesLeptonCut = 0;
            }
APICHART */
            
            if (doBJets < 0 && countBJets >= fabs(doBJets) )  { nEventsIncBJets++; 
                passesLeptonCut = 0  ; }

            if (doBJets > 0 && doBJets < 99 && countBJets < fabs(doBJets) ) passesLeptonCut = 0  ;
            if (doBJets == 101 && countBJets != 1) passesLeptonCut = 0  ;
            
        }  // END IF HAS RECO
        //=======================================================================================================//



        //=======================================================================================================//
        //        Print Event Information     //
        //====================================//
        if (PRINTEVENT && passesLeptonCut){
            vector<jetStruct> tmpJets;
            for (unsigned short i(0); i < nGoodJets; i++){
                if (jets[i].pt >= jetPtCutMin) tmpJets.push_back(jets[i]);
            }
            unsigned short tempnGoodJets(tmpJets.size());
            NZtotal++;
            cout << "event info: " << EvtRunNum << "  " << EvtNum << endl;
            cout << "Z event #" << NZtotal << "  Zmass : " << Z.M() << "  Zpt : " << Z.Pt() << " NJets : " << tempnGoodJets <<"    " <<weight << endl;
            if (nGoodJets > 0) cout << "JETS:"<< endl;
            for (unsigned short i(0); i < tempnGoodJets; i++) 
                cout << " jet #" << i + 1 << "  pt: " << tmpJets[i].pt << "  eta:"<<tmpJets[i].eta << "   " << endl;
            cout << "-----------------------------------------------------------------------------------------"<< endl;
        }
        //=======================================================================================================//



        //=======================================================================================================//
        //        Retrieving gen jets         //
        //====================================//
        bool passesGenJetCut(1), passesGenEWKJetPt(0), passesGenEWKJetFwdEta(0);
        unsigned short nGoodGenJets(0), nGenJetsAdd(0), nTotGenJets(0);
        double genJetsHT(0);
        vector<jetStruct> genJets, genJetsAdditional;
        TLorentzVector genLeadJ, genSecondJ, genJet1Plus2, genJet1Minus2;
        
        //*************************************** begin edit *************************************************************//
        TLorentzVector genNewLeadJ, genNewSecondJ, genNewThirdJ, genNewFourthJ;
        double genForwardJetRapidity(0), genBackwardJetRapidity(0);
        vector<TLorentzVector> genvJetYOrdered;
        //**************************************** end edit **************************************************************//

        if (hasGenInfo){
            nTotGenJets = genJetEta_->size();
            //-- retrieving generated jets
            for (unsigned short i(0); i < nTotGenJets; i++){
                jetStruct genJet = {genJetPt_->at(i), genJetEta_->at(i), genJetPhi_->at(i), genJetE_->at(i), i, 0};
                bool genJetPassesdRCut(1);
                double dRmin = 9999.;
                for (unsigned short j(0); j < nGenLeptons; j++){ 
                    if (genJet.pt >=  jetPtCutMin){
                        gendeltaRjetMu->Fill(deltaR(genJet.phi, genJet.eta, genLeptons[j].phi, genLeptons[j].eta), genWeight);
                    }
                    if ( deltaR(genJet.phi, genJet.eta, genLeptons[j].phi, genLeptons[j].eta) < dRmin ) dRmin = deltaR(genJet.phi, genJet.eta, genLeptons[j].phi, genLeptons[j].eta);
                    // I need this line because for to me unknown reason I CAN NO REMOVE ELECTRONS FROM Z IN SHERPA !!!!
                    if ((genLeptons[j].charge != 0)
                        && (doDR || (leptonFlavor == "Electrons" && fileName.find("HepMC") != string::npos))
                        && (deltaR(genJet.phi, genJet.eta, genLeptons[j].phi, genLeptons[j].eta) < 0.5)){
                        genJetPassesdRCut = 0;
                    }
                }
                //if (genJet.pt >= 10 && genJet.pt < 1000. && fabs(genJet.eta) <= 4.7 && genJetPassesdRCut){
                if (genJetPassesdRCut && genJet.pt >= 10 && fabs(genJet.eta) <= 4.7){
                    passesGenEWKJetPt = (genJet.pt >= 50);
                    passesGenEWKJetFwdEta = (fabs(genJet.eta) > 2.4);
                    genJets.push_back(genJet);                  
                    if (genJet.pt >=  15.) genJetsAdditional.push_back(genJet);
                }
            }
            nGoodGenJets = genJets.size();
            nGenJetsAdd = genJetsAdditional.size();
        }
        //=======================================================================================================//

if (DEBUG) cout << "Stop after line " << __LINE__ << "   " << hasGenInfo <<"    gen Wgh = " << genWeight << "  pass gen cuts = " << passesGenLeptonCut <<"  nGenJets = " << nGoodGenJets <<  endl;

        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        //=======================================================================================================//
        //     Matching gen and reco jets     //
        //====================================//
        vector<int> genJetsIndex(nGoodGenJets, 0);
        vector<vector<int> > matchingTable(nGoodJets, genJetsIndex);
        
        /* APICHART
        if (hasRecoInfo && hasGenInfo){
            for (unsigned short i(0); i < nGoodJets; i++){
                double mindR(0.5);
                int index(-1);
                double dR(9999);
                for (unsigned short j(0); j < nGoodGenJets; j++){
                    dR = deltaR(genJets[j].phi, genJets[j].eta, jets[i].phi, jets[i].eta);
                    if (dR < mindR){
                        mindR = dR;
                        index = j;
                    }
                }
                if (index > -1 ){
                    matchingTable[i][index] = 1;
                    //cout << " " << jentry << " " << i <<  " " << index << "  " << jets[i].pt << "  " << genJets[index].pt;
                    //cout << " " << jentry << " " << i <<  " " << index << "  " << jets[i].pt << "  " << jets[i].energy ;
                    //jets[i].pt = SmearJetPt(jets[i].pt, genJets[index].pt, jets[i].eta, smearJet);
                    double oldJetPt = jets[i].pt;
                    double newJetPt = SmearJetPt(oldJetPt, genJets[index].pt, jets[i].eta, smearJet);
                    jets[i].pt = newJetPt;
                    jets[i].energy = jets[i].energy * (newJetPt / oldJetPt);
                    //cout << "  new : " << jets[i].pt << "  " << jets[i].energy << endl;
                    puMVA_JetsMatchGenJets->Fill(JetAk04PuMva->at(jets[i].patIndex), weight);
                    puBeta_JetsMatchGenJets->Fill(JetAk04JetBeta->at(jets[i].patIndex), weight);
                    puBetaStar_JetsMatchGenJets->Fill(JetAk04JetBetaStar->at(jets[i].patIndex), weight);
                    jetsEta_JetsMatchGenJets->Fill(JetAk04Eta->at(jets[i].patIndex), weight);
                }
                else {
                    puMVA_JetsNoMatchGenJets->Fill(JetAk04PuMva->at(jets[i].patIndex), weight);
                    puBeta_JetsNoMatchGenJets->Fill(JetAk04JetBeta->at(jets[i].patIndex), weight);
                    puBetaStar_JetsNoMatchGenJets->Fill(JetAk04JetBetaStar->at(jets[i].patIndex), weight);
                    jetsEta_JetsNoMatchGenJets->Fill(JetAk04Eta->at(jets[i].patIndex), weight);
                }
            }

            //-- print the mathcing table
            //cout << "\n mathcing Table: \n" << endl; 
            //for (int i = 0; i < int(matchingTable.size()); i++){
            //  for (int j = 0; j < int(matchingTable[i].size()); j++){
            //    cout << matchingTable[i][j] << "  ";
            //  }
            //  cout << endl;
            //}

        }
         APICHART */
        
        
        // Re-analyze the jets collections and Cut on the Pt
        // we can do it only now since we needed to smear 
        // the jet pt distribution for the MC

        if (hasRecoInfo){     
            vector<jetStruct> tmpJets;
            for (unsigned short i(0); i < nGoodJets; i++){
                if (jets[i].pt >= jetPtCutMin) tmpJets.push_back(jets[i]);
            }
            jets.clear(); 
            jets = tmpJets; 
            tmpJets.clear(); 
            nGoodJets = jets.size();
            if (nGoodJets >= 1){
                sort(jets.begin(), jets.end(), JetDescendingOrder);
                sort(jetsAdditional.begin(), jetsAdditional.end(), JetDescendingOrder);
                leadJ.SetPtEtaPhiE(jets[0].pt, jets[0].eta, jets[0].phi, jets[0].energy);               
                //*************************************** begin edit *************************************************************//
                newLeadJ.SetPtEtaPhiE(jets[0].pt, jets[0].eta, jets[0].phi, jets[0].energy);
                vector<double> vJetRapidity;
                for (unsigned short i(0); i < nGoodJets; i++) {
                    TLorentzVector LVJet;
                    LVJet.SetPtEtaPhiE(jets[i].pt, jets[i].eta, jets[i].phi, jets[i].energy);
                    vJetRapidity.push_back(LVJet.Rapidity());
                    vJetYOrdered.push_back(LVJet);
                }
                ForwardJetRapidity = *max_element(vJetRapidity.begin(), vJetRapidity.end());
                BackwardJetRapidity = *min_element(vJetRapidity.begin(), vJetRapidity.end());
                sort(vJetYOrdered.begin(), vJetYOrdered.end(), JetYDescendingOrder);
                //**************************************** end edit ************************************************************//               
                if (jetPtCutMax > jetPtCutMin){
                    passesJetCut = jets[0].pt < jetPtCutMax;
                }
            }
            if (nGoodJets >= 2){
                secondJ.SetPtEtaPhiE(jets[1].pt, jets[1].eta, jets[1].phi, jets[1].energy);               
                //*************************************** begin edit *************************************************************//
                newSecondJ.SetPtEtaPhiE(jets[1].pt, jets[1].eta, jets[1].phi, jets[1].energy);
                //**************************************** end edit ************************************************************// 
		        jet1Plus2 = leadJ + secondJ;
                jet1Minus2 = leadJ - secondJ;
            }
            //*************************************** begin edit *************************************************************//
            if (nGoodJets >= 3){
                newThirdJ.SetPtEtaPhiE(jets[2].pt, jets[2].eta, jets[2].phi, jets[2].energy);
            }
            if (nGoodJets >= 4){
                newFourthJ.SetPtEtaPhiE(jets[3].pt, jets[3].eta, jets[3].phi, jets[3].energy);
            }
            //**************************************** end edit ************************************************************//
	        jetsHT = 0;
            for (unsigned short i(0); i < nGoodJets; i++){
                jetsHT += jets[i].pt;  
            }
        }

        if (hasGenInfo){
            vector< jetStruct> tmpJets;
            for (unsigned short i(0); i < nGoodGenJets; i++){
                if (genJets[i].pt >= jetPtCutMin && genJets[i].eta >= double (jetEtaCutMin/10.) && genJets[i].eta <= (double )(jetEtaCutMax/10.) ){
                    tmpJets.push_back(genJets[i]);
                }
            }
            genJets.clear();
            genJets = tmpJets; 
            tmpJets.clear(); 
            nGoodGenJets = genJets.size();
            if (nGoodGenJets >= 1){
                sort(genJets.begin(), genJets.end(), JetDescendingOrder);
                genLeadJ.SetPtEtaPhiE(genJets[0].pt, genJets[0].eta, genJets[0].phi, genJets[0].energy);
                //*************************************** begin edit ***********************************************************//
                genNewLeadJ.SetPtEtaPhiE(genJets[0].pt, genJets[0].eta, genJets[0].phi, genJets[0].energy);
                vector<double> vJetRapidity;
                for (unsigned short i(0); i < nGoodGenJets; i++) {
                    TLorentzVector LVJet;
                    LVJet.SetPtEtaPhiE(genJets[i].pt, genJets[i].eta, genJets[i].phi, genJets[i].energy);
                    vJetRapidity.push_back(LVJet.Rapidity());
                    genvJetYOrdered.push_back(LVJet);
                }
                genForwardJetRapidity = *max_element(vJetRapidity.begin(), vJetRapidity.end());
                genBackwardJetRapidity = *min_element(vJetRapidity.begin(), vJetRapidity.end());
                sort(genvJetYOrdered.begin(), genvJetYOrdered.end(), JetYDescendingOrder);
                //**************************************** end edit ************************************************************//
                if (jetPtCutMax > jetPtCutMin){
                    passesGenJetCut = genJets[0].pt < jetPtCutMax;
                }
            }
            if (nGoodGenJets >= 2){
                genSecondJ.SetPtEtaPhiE(genJets[1].pt, genJets[1].eta, genJets[1].phi, genJets[1].energy);
                //*************************************** begin edit *************************************************************//
                genNewSecondJ.SetPtEtaPhiE(genJets[1].pt, genJets[1].eta, genJets[1].phi, genJets[1].energy);
                //**************************************** end edit ************************************************************//
                genJet1Plus2 = genLeadJ + genSecondJ;
                genJet1Minus2 = genLeadJ - genSecondJ;
            }
            //*************************************** begin edit *************************************************************//
            if (nGoodGenJets >= 3){
                genNewThirdJ.SetPtEtaPhiE(genJets[2].pt, genJets[2].eta, genJets[2].phi, genJets[2].energy);
            }
            if (nGoodGenJets >= 4){
                genNewFourthJ.SetPtEtaPhiE(genJets[3].pt, genJets[3].eta, genJets[3].phi, genJets[3].energy);
            }
            //**************************************** end edit ************************************************************//
            genJetsHT = 0.;
            for (unsigned short i(0); i < nGoodGenJets; i++){
                genJetsHT += genJets[i].pt;  
            }

        }
        //=======================================================================================================//


        //=======================================================================================================//
        // Select the best pair of jets for DPS  //
        //=======================================//
        pair<TLorentzVector, TLorentzVector> bestTwoJets;
        TLorentzVector bestJet1Plus2, bestJet1Minus2;
        if (hasRecoInfo){
            bestTwoJetsCandidatesPt(jets, bestTwoJets);
            //bestTwoJetsCandidatesPhi(jets, bestTwoJets);
            if (nGoodJets >= 2){
                bestJet1Plus2 = bestTwoJets.first + bestTwoJets.second;
                bestJet1Minus2 = bestTwoJets.first - bestTwoJets.second;
            }
        }

        pair<TLorentzVector, TLorentzVector> genBestTwoJets;
        TLorentzVector genBestJet1Plus2, genBestJet1Minus2;
        if (hasGenInfo){
            bestTwoJetsCandidatesPt(genJets, genBestTwoJets);
            //bestTwoJetsCandidatesPhi(genJets, genBestTwoJets);
            if (nGoodGenJets >= 2){
                genBestJet1Plus2 = genBestTwoJets.first + genBestTwoJets.second; 
                genBestJet1Minus2 = genBestTwoJets.first - genBestTwoJets.second; 
            }
        }
        //=======================================================================================================//


        //=======================================================================================================//
        //   Retrieving DPS partons and jets  //
        //====================================//
        bool passesDPSPartonCut(1); 
        unsigned short nTotDPSPartons(0);
        vector<pair<TLorentzVector, TLorentzVector> > genDPSAndItsJets;
        if (hasPartonInfo){
            nTotDPSPartons = dpsParton_Eta->size();
            short jetInd(0);
            for (unsigned short i(0); i < nTotDPSPartons; i++){
                TLorentzVector genDPSParton, genDPSJet;
                genDPSParton.SetPtEtaPhiE(dpsParton_Pt->at(i), dpsParton_Eta->at(i), dpsParton_Phi->at(i), dpsParton_E->at(i));
                jetInd = genMatchDPSpar->at(i); 
                if (jetInd >= 0){
                    genDPSJet.SetPtEtaPhiE(genJetPt_->at(jetInd), genJetEta_->at(jetInd), genJetPhi_->at(jetInd), genJetE_->at(jetInd));
                    pair<TLorentzVector, TLorentzVector> genDPSAndItsJet(genDPSParton, genDPSJet);
                    genDPSAndItsJets.push_back(genDPSAndItsJet);
                }
            }
            // remove DPS from inclusive sample
            if (fileName.find("PartonInfo") != string::npos && fileName.find("enriched") != string::npos &&
                    (nTotDPSPartons < 2 || fabs(dpsParton_Eta->at(0)) > 2.4 || fabs(dpsParton_Eta->at(1)) > 2.4)){
                passesDPSPartonCut = 0;
            }
        }
        //=======================================================================================================//

        // Wb study
        if ( doWbsyst && countWbBjets == 1 && hasGenInfo) {
            weight = weight * WbSystSF;
            genWeight = genWeight * WbSystSF;
        }

        if (DEBUG) cout << "Stop after line " << __LINE__ << "   " << hasGenInfo <<"    gen Wgh = " << genWeight << "  pass gen cuts = " << passesGenLeptonCut <<"  nGenJets = " << nGoodGenJets <<  endl;
        //=======================================================================================================//
        //   Filling gen and parton histos    //
        //====================================//
        if (hasGenInfo){
            if (passesGenLeptonCut && passesDPSPartonCut && passesGenJetCut){
                if (fileName.find("HepMC") == string::npos && pdfInfo_->size()>3){
                    partonX2D->Fill(pdfInfo_->at(2),pdfInfo_->at(3),genWeight);
                }
                GENnEventsIncl0Jets++;
                genZNGoodJets_Zexc->Fill(nGoodGenJets, genWeight);
                genZNGoodJetsFull_Zexc->Fill(nGoodGenJets, genWeight);
                genZNGoodJets_Zinc->Fill(0., genWeight);
                genZNGoodJetsFull_Zinc->Fill(0., genWeight);
                genZMass_Zinc0jet->Fill(genZ.M(), genWeight);
                genZPt_Zinc0jet->Fill(genZ.Pt(), genWeight);
                genZRapidity_Zinc0jet->Fill(genZ.Rapidity(), genWeight);
                genZEta_Zinc0jet->Fill(genZ.Eta(), genWeight);
                genlepPt_Zinc0jet->Fill(genLep1.Pt(), genWeight);
                genlepEta_Zinc0jet->Fill(genLep1.Eta(), genWeight);
                if (doZ || doTT){
                    genlepPt_Zinc0jet->Fill(genLep2.Pt(), genWeight);
                    genlepEta_Zinc0jet->Fill(genLep2.Eta(), genWeight);
                }
                if (nGoodGenJets >= 1){
                    GENnEventsIncl1Jets++;
                    genZNGoodJets_Zinc->Fill(1., genWeight);
                    genZNGoodJetsFull_Zinc->Fill(1., genWeight);
                    genZMass_Zinc1jet->Fill(genZ.M(), genWeight);
                    genlepPt_Zinc1jet->Fill(genLep1.Pt(), genWeight);
                    genlepEta_Zinc1jet->Fill(genLep1.Eta(), genWeight);
                    if (doZ || doTT){
                        genlepPt_Zinc1jet->Fill(genLep2.Pt(), genWeight);
                        genlepEta_Zinc1jet->Fill(genLep2.Eta(), genWeight);
                    }
                    genZPt_Zinc1jet->Fill(genZ.Pt(), genWeight);
                    genZRapidity_Zinc1jet->Fill(genZ.Rapidity(), genWeight);
                    genZEta_Zinc1jet->Fill(genZ.Eta(), genWeight);
                    genFirstJetPt_Zinc1jet->Fill(genLeadJ.Pt(), genWeight);
                    genFirstJetPt_1_Zinc1jet->Fill(genLeadJ.Pt(), genWeight);
                    genFirstJetPt_2_Zinc1jet->Fill(genLeadJ.Pt(), genWeight);
                    
                    genFirstJetEta_Zinc1jet->Fill(fabs(genLeadJ.Eta()), genWeight);
                    genFirstJetEta_2_Zinc1jet->Fill(fabs(genLeadJ.Eta()), genWeight);
                    //*************************************** begin edit *************************************************************//
                    genFirstJetRapidity_Zinc1jet->Fill(fabs(genNewLeadJ.Rapidity()), genWeight);
                    genFirstJetRapidityFull_Zinc1jet->Fill(genNewLeadJ.Rapidity(), genWeight);
                    genMeanNJetsHT_1D_Zinc1jet->Fill(genJetsHT, genWeight*nGoodGenJets);
                    genMeanNJetsHT_Zinc1jet->Fill(genJetsHT, nGoodGenJets, genWeight);
                    //*************************************** end edit ***************************************************************//
                    genFirstJetPtEta_Zinc1jet->Fill(genLeadJ.Pt(), fabs(genLeadJ.Eta()), genWeight);
                    genFirstHighestJetPt_Zinc1jet->Fill(genLeadJ.Pt(), genWeight);
                    genJetsHT_Zinc1jet->Fill(genJetsHT, genWeight);
                    genJetsHT_1_Zinc1jet->Fill(genJetsHT, genWeight);
                    genJetsHT_2_Zinc1jet->Fill(genJetsHT, genWeight);
                    
                    for ( int i =0 ; i < NbinsEta2D - 1 ; i++){
                        if ( fabs(genLeadJ.Eta()) >= j_Y_range[i] &&  fabs(genLeadJ.Eta()) < j_Y_range[i+1] )                                genFirstJetPt_Zinc1jet_Eta[i]->Fill(fabs(genLeadJ.Pt()), genWeight);
                    }
                    if ( doW ) gendEtaBosonJet_Zinc1jet->Fill(fabs(genLeadJ.Eta() - genLep1.Eta()), genWeight);
                    else gendEtaBosonJet_Zinc1jet->Fill(fabs(genLeadJ.Eta()-genZ.Eta()), genWeight);
                    if (nGoodGenJets == 1){
                        genFirstJetPt_Zexc1jet->Fill(genLeadJ.Pt(), genWeight);
                        //    gendEtaBosonJet_Zexc1jet->Fill(fabs(genLeadJ.Eta()), genWeight);
                        if ( doW ) gendEtaBosonJet_Zexc1jet->Fill(fabs(genLeadJ.Eta() - genLep1.Eta()), genWeight);
                        else gendEtaBosonJet_Zexc1jet->Fill(fabs(genLeadJ.Eta()-genZ.Eta()), genWeight);

                    }
                }
                if (nGoodGenJets >= 2){
                    TLorentzVector genJet1Plus2PlusZ = genJet1Plus2 + genZ;
                    GENnEventsIncl2Jets++;
                    genZNGoodJets_Zinc->Fill(2., genWeight);
                    genZNGoodJetsFull_Zinc->Fill(2., genWeight);
                    genZMass_Zinc2jet->Fill(genZ.M(), genWeight);
                    genTwoJetsPtDiff_Zinc2jet->Fill(genJet1Minus2.Pt(), genWeight);
                    genBestTwoJetsPtDiff_Zinc2jet->Fill(genBestJet1Minus2.Pt(), genWeight);
                    genJetsMass_Zinc2jet->Fill(genJet1Plus2.M(), genWeight);
                    genlepPt_Zinc2jet->Fill(genLep1.Pt(), genWeight);
                    genlepEta_Zinc2jet->Fill(genLep1.Eta(), genWeight);
                    if (doZ || doTT){
                        genlepPt_Zinc2jet->Fill(genLep2.Pt(), genWeight);
                        genlepEta_Zinc2jet->Fill(genLep2.Eta(), genWeight);
                    }
                    genZPt_Zinc2jet->Fill(genZ.Pt(), genWeight);
                    genZRapidity_Zinc2jet->Fill(genZ.Rapidity(), genWeight);
                    genZEta_Zinc2jet->Fill(genZ.Eta(), genWeight);
                    genFirstHighestJetPt_Zinc2jet->Fill(genLeadJ.Pt(), genWeight);
                    genSecondHighestJetPt_Zinc2jet->Fill(genSecondJ.Pt(), genWeight);
                    genSecondJetPt_Zinc2jet->Fill(genSecondJ.Pt(), genWeight);
                    genSecondJetPt_1_Zinc2jet->Fill(genSecondJ.Pt(), genWeight);
                    genSecondJetPt_2_Zinc2jet->Fill(genSecondJ.Pt(), genWeight);
                    
                    genSecondJetEta_Zinc2jet->Fill(fabs(genSecondJ.Eta()), genWeight);
                    genSecondJetEta_2_Zinc2jet->Fill(fabs(genSecondJ.Eta()), genWeight);
                    //*************************************** begin edit *******************************************************//
                    genSecondJetRapidity_Zinc2jet->Fill(fabs(genNewSecondJ.Rapidity()), genWeight);
                    genSecondJetRapidityFull_Zinc2jet->Fill(genNewSecondJ.Rapidity(), genWeight);
                    gendRapidityJets_Zinc2jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight);
                    gendRapidityJets_2_Zinc2jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight);
                    
                    gendRapidityJetsFB_Zinc2jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight);
                    gendRapidityJetsFB_2_Zinc2jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight);
                    
                    gendiJetMass_Zinc2jet->Fill(genJet1Plus2.M(), genWeight);
                    gendiJetMass_2_Zinc2jet->Fill(genJet1Plus2.M(), genWeight);
                    
                    gendPhiJetsFB_Zinc2jet->Fill(deltaPhi(genvJetYOrdered[0], genvJetYOrdered[genvJetYOrdered.size() - 1] ), genWeight);
                    gendPhiJetsFB_2_Zinc2jet->Fill(deltaPhi(genvJetYOrdered[0], genvJetYOrdered[genvJetYOrdered.size() - 1] ), genWeight);
                    
                    gendRJets_Zinc2jet->Fill(deltaRYPhi(genNewLeadJ, genNewSecondJ), genWeight);
                    gendRJets_2_Zinc2jet->Fill(deltaRYPhi(genNewLeadJ, genNewSecondJ), genWeight);
                    
                    gendiJetPt_Zinc2jet->Fill(genJet1Plus2.Pt(), genWeight);
                    gendiJetPt_2_Zinc2jet->Fill(genJet1Plus2.Pt(), genWeight);
                    
                    genMeanNJetsHT_1D_Zinc2jet->Fill(genJetsHT, genWeight*nGoodGenJets);
                    genMeanNJetsdRapidity_1D_Zinc2jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight*nGoodGenJets);
                    genMeanNJetsdRapidityFB_1D_Zinc2jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight*nGoodGenJets);
                    
                    genMeanNJetsHT_Zinc2jet->Fill(genJetsHT, nGoodGenJets, genWeight);
                    genMeanNJetsdRapidity_Zinc2jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), nGoodGenJets, genWeight);
                    genMeanNJetsdRapidityFB_Zinc2jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, nGoodGenJets, genWeight);
                    //*************************************** end edit *********************************************************//
                    genSecondJetPtEta_Zinc2jet->Fill(genSecondJ.Pt(), fabs(genSecondJ.Eta()), genWeight);
                    genRatioJetPt21_Zinc2jet->Fill(genJets[1].pt/genJets[0].pt, genWeight);
                    genJetsHT_Zinc2jet->Fill(genJetsHT, genWeight);
                    genJetsHT_1_Zinc2jet->Fill(genJetsHT, genWeight);
                    genJetsHT_2_Zinc2jet->Fill(genJetsHT, genWeight);
                    
                    genptBal_Zinc2jet->Fill(genJet1Plus2PlusZ.Pt(), genWeight);
                    
                    gendPhiJets_Zinc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                    gendPhiJets_2_Zinc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                    
                    genBestdPhiJets_Zinc2jet->Fill(deltaPhi(genBestTwoJets.first, genBestTwoJets.second), genWeight);
                    gendEtaJets_Zinc2jet->Fill(genLeadJ.Eta() - genSecondJ.Eta(), genWeight);
                    gendEtaFirstJetZ_Zinc2jet->Fill(genLeadJ.Eta() - genZ.Eta(), genWeight);
                    gendEtaSecondJetZ_Zinc2jet->Fill(genSecondJ.Eta() - genZ.Eta(), genWeight);
                    gendEtaJet1Plus2Z_Zinc2jet->Fill(genJet1Plus2.Eta() - genZ.Eta(), genWeight);
                    genPHI_Zinc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    genBestPHI_Zinc2jet->Fill(PHI(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                    genPHI_T_Zinc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    genBestPHI_T_Zinc2jet->Fill(PHI_T(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                    genSpT_Zinc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    genBestSpT_Zinc2jet->Fill(SpT(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                    genSpTJets_Zinc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                    genBestSpTJets_Zinc2jet->Fill(SpTsub(genBestTwoJets.first, genBestTwoJets.second), genWeight);
                    genSpTLeptons_Zinc2jet->Fill(SpTsub(genLep1, genLep2), genWeight);
                    genSPhi_Zinc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    genBestSPhi_Zinc2jet->Fill(SPhi(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                    for ( int i =0 ; i < NbinsEta2D - 1 ; i++){
                        if ( fabs(genSecondJ.Eta()) >= j_Y_range[i] &&  fabs(genSecondJ.Eta()) < j_Y_range[i+1] )                                                genSecondJetPt_Zinc2jet_Eta[i]->Fill(fabs(genSecondJ.Pt()), genWeight);
                    }


                    //--- START EWK ---
                    if (passesGenEWKJetPt){
                        // met histos
                        //
                        //METEWK_Zinc2jet->Fill(METpt, weight);
                        //METphiEWK_Zinc2jet->Fill(METphi,  weight);
                        //MTEWK_Zinc2jet->Fill(MT,  weight);
                        //          
                        // jet histos
                        genZNGoodJetsEWK_Zexc->Fill(nGoodGenJets, genWeight);
                        genFirstJetPtEWK_Zinc2jet->Fill(genJets[0].pt, genWeight);
                        genFirstJetEtaEWK_Zinc2jet->Fill(genJets[0].eta, genWeight);
                        genFirstJetPhiEWK_Zinc2jet->Fill(genJets[0].phi, genWeight);


                        genSecondJetPtEWK_Zinc2jet->Fill(genJets[1].pt, genWeight);
                        genSecondJetEtaEWK_Zinc2jet->Fill(genJets[1].eta, genWeight);
                        genSecondJetPhiEWK_Zinc2jet->Fill(genJets[1].phi, genWeight);
                        int temIND(0), temIND1(1);
                        if (fabs(genJets[1].eta) > fabs(genJets[0].eta)){ 
                            temIND = 1; 
                            temIND1 = 0;
                        }
                        genForwardJetPtEWK_Zinc2jet->Fill(genJets[temIND].pt, genWeight);
                        genForwardJetEtaEWK_Zinc2jet->Fill(genJets[temIND].eta, genWeight);
                        genForwardJetPhiEWK_Zinc2jet->Fill(genJets[temIND].phi, genWeight);


                        genCentralJetPtEWK_Zinc2jet->Fill(genJets[temIND1].pt, genWeight);
                        genCentralJetEtaEWK_Zinc2jet->Fill(genJets[temIND1].eta, genWeight);
                        genCentralJetPhiEWK_Zinc2jet->Fill(genJets[temIND1].phi, genWeight);


                        genJetsHTEWK_Zinc2jet->Fill(genJetsHT, genWeight);
                        genJetsMassEWK_Zinc2jet->Fill(genJet1Plus2.M(), genWeight);

                        // multi jet variables
                        genSumEtaJetsEWK_Zinc2jet->Fill(fabs(genLeadJ.Eta() + genSecondJ.Eta()),  genWeight);
                        genAbsSumEtaJetsEWK_Zinc2jet->Fill(fabs(genLeadJ.Eta()) + fabs(genSecondJ.Eta()),  genWeight);
                        genTwoJetsPtDiffEWK_Zinc2jet->Fill(genJet1Minus2.Pt(), genWeight);
                        genptBalEWK_Zinc2jet->Fill(genJet1Plus2PlusZ.Pt(), genWeight);
                        gendPhiJetsEWK_Zinc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                        gendEtaJetsEWK_Zinc2jet->Fill(fabs(genLeadJ.Eta() - genSecondJ.Eta()), genWeight);
                        genSpTJetsEWK_Zinc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                        gendPhiJetsEWK_Zinc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);

                        // find jet properties of the third jet that is  between the two leading jets
                        int nGoodJetsBtw(0.);
                        double jetsHTBtw(0.);
                        if (nGenJetsAdd > 2){
                            genThirdJetPtEWKadd_Zinc2jet->Fill(genJetsAdditional[2].pt, genWeight);
                            genThirdJetEtaEWKadd_Zinc2jet->Fill(genJetsAdditional[2].eta, genWeight);
                            genThirdJetPhiEWKadd_Zinc2jet->Fill(genJetsAdditional[2].phi, genWeight);
                            for (unsigned short i(2); i < nGenJetsAdd; i++) {
                                if (genJetsAdditional[i].eta < max(genJets[0].eta,genJets[1].eta) -0.5 && genJetsAdditional[i].eta > min(genJets[0].eta,genJets[1].eta) + 0.5){
                                    jetsHTBtw += genJetsAdditional[i].pt;
                                    nGoodJetsBtw++;
                                    genAllJetPtEWKbtw_Zinc2jet->Fill(genJetsAdditional[i].pt, genWeight);
                                    genAllJetEtaEWKbtw_Zinc2jet->Fill(genJetsAdditional[i].eta, genWeight);
                                    genAllJetPhiEWKbtw_Zinc2jet->Fill(genJetsAdditional[i].phi, genWeight);
                                    if (nGoodJetsBtw == 1){
                                        genThirdJetPtEWKbtw_Zinc2jet->Fill(genJetsAdditional[i].pt, genWeight);
                                        genThirdJetEtaEWKbtw_Zinc2jet->Fill(genJetsAdditional[i].eta, genWeight);
                                        genThirdJetPhiEWKbtw_Zinc2jet->Fill(genJetsAdditional[i].phi, genWeight);
                                    }
                                }
                            }
                            genJetsHTEWKbtw_Zinc2jet->Fill(jetsHTBtw, genWeight);
                            genZNGoodJetsEWKbtw_Zexc->Fill(nGoodJetsBtw, genWeight);
                        }

                    }
                    //////////////////////// STOP EWK
                    if (genZ.Pt() < 25){
                        genptBal_LowPt_Zinc2jet->Fill(genJet1Plus2PlusZ.Pt(), genWeight);
                        gendPhiJets_LowPt_Zinc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                        genBestdPhiJets_LowPt_Zinc2jet->Fill(deltaPhi(genBestTwoJets.first, genBestTwoJets.second), genWeight);
                        gendPhiLeptons_LowPt_Zinc2jet->Fill(deltaPhi(genLep1, genLep2), genWeight);
                        genPHI_T_LowPt_Zinc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genBestPHI_T_LowPt_Zinc2jet->Fill(PHI_T(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                        genPHI_LowPt_Zinc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genBestPHI_LowPt_Zinc2jet->Fill(PHI(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                        genSpTJets_LowPt_Zinc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                        genBestSpTJets_LowPt_Zinc2jet->Fill(SpTsub(genBestTwoJets.first, genBestTwoJets.second), genWeight);
                        genSpTLeptons_LowPt_Zinc2jet->Fill(SpTsub(genLep1, genLep2), genWeight);
                        genSpT_LowPt_Zinc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genBestSpT_LowPt_Zinc2jet->Fill(SpT(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                        genSPhi_LowPt_Zinc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genBestSPhi_LowPt_Zinc2jet->Fill(SPhi(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                        if (SpT(genLep1, genLep2, genLeadJ, genSecondJ) < 0.5){ 
                            genPHI_LowSpT_LowPt_Zinc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            genSPhi_LowSpT_LowPt_Zinc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        }
                        else {
                            genPHI_HighSpT_LowPt_Zinc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            genSPhi_HighSpT_LowPt_Zinc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        }
                        if (SPhi(genLep1, genLep2, genLeadJ, genSecondJ) < 0.5){
                            genSpT_LowSPhi_LowPt_Zinc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        }
                        else {
                            genSpT_HighSPhi_LowPt_Zinc2jet ->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        }
                    }
                    else {
                        genptBal_HighPt_Zinc2jet->Fill(genJet1Plus2PlusZ.Pt(), genWeight);
                        gendPhiJets_HighPt_Zinc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                        gendPhiLeptons_HighPt_Zinc2jet->Fill(deltaPhi(genLep1, genLep2), genWeight);
                        genPHI_HighPt_Zinc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genPHI_T_HighPt_Zinc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genSpTJets_HighPt_Zinc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                        genSpTLeptons_HighPt_Zinc2jet->Fill(SpTsub(genLep1, genLep2), genWeight);
                        genSpT_HighPt_Zinc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genSPhi_HighPt_Zinc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    }
                    if (nGoodGenJets == 2){
                        genTwoJetsPtDiff_Zexc2jet->Fill(genJet1Minus2.Pt(), genWeight);
                        genJetsMass_Zexc2jet->Fill(genJet1Plus2.M(), genWeight);
                        genSecondJetPt_Zexc2jet->Fill(genSecondJ.Pt(), genWeight);
                        gendPhiJets_Zexc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                        genPHI_Zexc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genPHI_T_Zexc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        gendEtaJets_Zexc2jet->Fill(genLeadJ.Eta() - genSecondJ.Eta(), genWeight);
                        gendEtaFirstJetZ_Zexc2jet->Fill(genLeadJ.Eta() - genZ.Eta(), genWeight);
                        gendEtaSecondJetZ_Zexc2jet->Fill(genSecondJ.Eta() - genZ.Eta(), genWeight);
                        gendEtaJet1Plus2Z_Zexc2jet->Fill(genJet1Plus2.Eta() - genZ.Eta(), genWeight);
                        genSpT_Zexc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genSpTJets_Zexc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                        genSpTLeptons_Zexc2jet->Fill(SpTsub(genLep1, genLep2), genWeight);
                        genSPhi_Zexc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genptBal_Zexc2jet->Fill(genJet1Plus2PlusZ.Pt(), genWeight);

                        if (genZ.Pt() < 25){
                            genptBal_LowPt_Zexc2jet->Fill(genJet1Plus2PlusZ.Pt(), genWeight);
                            gendPhiJets_LowPt_Zexc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                            gendPhiLeptons_LowPt_Zexc2jet->Fill(deltaPhi(genLep1, genLep2), genWeight);
                            genPHI_T_LowPt_Zexc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            genPHI_LowPt_Zexc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            genSpTJets_LowPt_Zexc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                            genSpTLeptons_LowPt_Zexc2jet->Fill(SpTsub(genLep1, genLep2), genWeight);
                            genSpT_LowPt_Zexc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            genSPhi_LowPt_Zexc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            if (SpT(genLep1, genLep2, genLeadJ, genSecondJ) < 0.5) { 
                                genPHI_LowSpT_LowPt_Zexc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                                genSPhi_LowSpT_LowPt_Zexc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            }
                            else {
                                genPHI_HighSpT_LowPt_Zexc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                                genSPhi_HighSpT_LowPt_Zexc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            }
                            if (SPhi(genLep1, genLep2, genLeadJ, genSecondJ) < 0.5) {
                                genSpT_LowSPhi_LowPt_Zexc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            }
                            else {
                                genSpT_HighSPhi_LowPt_Zexc2jet ->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            }
                        }
                        else {
                            genptBal_HighPt_Zexc2jet->Fill(genJet1Plus2PlusZ.Pt(), genWeight);
                            gendPhiJets_HighPt_Zexc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                            gendPhiLeptons_HighPt_Zexc2jet->Fill(deltaPhi(genLep1, genLep2), genWeight);
                            genPHI_HighPt_Zexc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            genPHI_T_HighPt_Zexc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            genSpTJets_HighPt_Zexc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                            genSpTLeptons_HighPt_Zexc2jet->Fill(SpTsub(genLep1, genLep2), genWeight);
                            genSpT_HighPt_Zexc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            genSPhi_HighPt_Zexc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        }
                    }

                    if (hasPartonInfo){
                        if (genDPSAndItsJets.size() == 2) {
                            if ((deltaR(genDPSAndItsJets[0].first, genLeadJ) < 0.5 || deltaR(genDPSAndItsJets[0].first, genSecondJ) < 0.5) 
                                    && (deltaR(genDPSAndItsJets[1].first, genLeadJ) < 0.5 || deltaR(genDPSAndItsJets[1].first, genSecondJ) < 0.5)) {
                                genSpTJetsDeltaR_Zexc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                                gendPhiJetsDeltaR_Zexc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                                genPHI_TDeltaR_Zexc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                                genSpTDeltaR_Zexc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                                if (hasRecoInfo) {
                                    resSpTJetsDeltaR_Zexc2jet->Fill(SpTsub(genLeadJ, genSecondJ) - SpTsub(leadJ, secondJ) , genWeight);
                                    resdPhiJetsDeltaR_Zexc2jet->Fill(deltaPhi(genLeadJ, genSecondJ) - deltaPhi(leadJ, secondJ), genWeight);
                                    resPHI_TDeltaR_Zexc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ) - PHI_T(lep1, lep2, leadJ, secondJ), genWeight);
                                    resSpTDeltaR_Zexc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ) - SpT(lep1, lep2, leadJ, secondJ), genWeight);
                                }
                            }
                            gendPhiJetsDPS_Zexc2jet->Fill(deltaPhi(genDPSAndItsJets[0].second, genDPSAndItsJets[1].second), genWeight);
                            genPHI_TDPS_Zexc2jet->Fill(PHI_T(genLep1, genLep2, genDPSAndItsJets[0].second, genDPSAndItsJets[1].second), genWeight);
                            genSpTDPSPartons_Zexc2jet->Fill(SpTsub(genDPSAndItsJets[0].first, genDPSAndItsJets[1].first), genWeight);
                            genSpTJetsDPS_Zexc2jet->Fill(SpTsub(genDPSAndItsJets[0].second, genDPSAndItsJets[1].second), genWeight);
                            genSpTDPS_Zexc2jet->Fill(SpT(genLep1, genLep2, genDPSAndItsJets[0].second, genDPSAndItsJets[1].second), genWeight);
                            if (deltaR(genDPSAndItsJets[0].first, genDPSAndItsJets[0].second) < 0.5 && deltaR(genDPSAndItsJets[1].first, genDPSAndItsJets[1].second)) {
                                gendPhiJetsDPSDeltaR_Zexc2jet->Fill(deltaPhi(genDPSAndItsJets[0].second, genDPSAndItsJets[1].second), genWeight);
                                genSpTJetsDPSDeltaR_Zexc2jet->Fill(SpTsub(genDPSAndItsJets[0].second, genDPSAndItsJets[1].second), genWeight);
                                genSpTDPSDeltaR_Zexc2jet->Fill(SpT(genLep1, genLep2, genDPSAndItsJets[0].second, genDPSAndItsJets[1].second), genWeight);
                                genPHI_TDPSDeltaR_Zexc2jet->Fill(PHI_T(genLep1, genLep2, genDPSAndItsJets[0].second, genDPSAndItsJets[1].second), genWeight);
                                gendPhiJetsDPSDeltaR_ZpT_Zexc2jet->Fill(deltaPhi(genDPSAndItsJets[0].second, genDPSAndItsJets[1].second), genZ.Pt(), genWeight);
                            }
                        }
                    }
                }
                if (nGoodGenJets >= 3){
                    GENnEventsIncl3Jets++;
                    genZNGoodJets_Zinc->Fill(3., genWeight);
                    genZNGoodJetsFull_Zinc->Fill(3., genWeight);
                    genZMass_Zinc3jet->Fill(genZ.M(), genWeight);
                    genZPt_Zinc3jet->Fill(genZ.Pt(), genWeight);
                    genZRapidity_Zinc3jet->Fill(genZ.Rapidity(), genWeight);
                    genZEta_Zinc3jet->Fill(genZ.Eta(), genWeight);
                    genThirdJetPt_Zinc3jet->Fill(genJets[2].pt, genWeight);
                    genThirdJetPt_1_Zinc3jet->Fill(genJets[2].pt, genWeight);
                    genThirdJetPt_2_Zinc3jet->Fill(genJets[2].pt, genWeight);
                    
                    genThirdJetEta_Zinc3jet->Fill(fabs(genJets[2].eta), genWeight);
                    genThirdJetEta_2_Zinc3jet->Fill(fabs(genJets[2].eta), genWeight);
                    //*************************************** begin edit *************************************************************//
                    genThirdJetRapidity_Zinc3jet->Fill(fabs(genNewThirdJ.Rapidity()), genWeight);
                    genThirdJetRapidityFull_Zinc3jet->Fill(genNewThirdJ.Rapidity(), genWeight);
                    gendRapidityJets_Zinc3jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight);
                    gendRapidityJets_2_Zinc3jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight);
                    
                    gendRapidityJets_First_Third_Zinc3jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewThirdJ.Rapidity()), genWeight);
                    gendRapidityJets_2_First_Third_Zinc3jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewThirdJ.Rapidity()), genWeight);
                    
                    gendRapidityJets_Second_Third_Zinc3jet->Fill(fabs(genNewSecondJ.Rapidity() - genNewThirdJ.Rapidity()), genWeight);
                    gendRapidityJets_2_Second_Third_Zinc3jet->Fill(fabs(genNewSecondJ.Rapidity() - genNewThirdJ.Rapidity()), genWeight);
                    
                    gendRapidityJetsFB_Zinc3jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight);
                    gendRapidityJetsFB_2_Zinc3jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight);
                    
                    gendiJetMass_Zinc3jet->Fill(genJet1Plus2.M(), genWeight);
                    gendiJetMass_2_Zinc3jet->Fill(genJet1Plus2.M(), genWeight);
                    
                    gendiJetPt_Zinc3jet->Fill(genJet1Plus2.Pt(), genWeight);
                    gendiJetPt_2_Zinc3jet->Fill(genJet1Plus2.Pt(), genWeight);
                    //*************************************** end edit ***************************************************************//
                    genThirdJetPtEta_Zinc3jet->Fill(genJets[2].pt, fabs(genJets[2].eta), genWeight);
                    genRatioJetPt32_Zinc3jet->Fill(genJets[2].pt/genJets[1].pt, genWeight);
                    genFirstHighestJetPt_Zinc3jet->Fill(genJets[0].pt, genWeight);
                    genSecondHighestJetPt_Zinc3jet->Fill(genJets[1].pt, genWeight);
                    genThirdHighestJetPt_Zinc3jet->Fill(genJets[2].pt, genWeight);
                    genJetsHT_Zinc3jet->Fill(genJetsHT, genWeight);
                    genJetsHT_1_Zinc3jet->Fill(genJetsHT, genWeight);
                    genJetsHT_2_Zinc3jet->Fill(genJetsHT, genWeight);
                    
                }
                if (nGoodGenJets >= 4){
                    genZNGoodJets_Zinc->Fill(4., genWeight);
                    genZNGoodJetsFull_Zinc->Fill(4., genWeight);
                    genZMass_Zinc4jet->Fill(genZ.M(), genWeight);
                    genZPt_Zinc4jet->Fill(genZ.Pt(), genWeight);
                    genZRapidity_Zinc4jet->Fill(genZ.Rapidity(), genWeight);
                    genZEta_Zinc4jet->Fill(genZ.Eta(), genWeight);
                    genFourthJetPt_Zinc4jet->Fill(genJets[3].pt, genWeight);
                    genFourthJetPt_1_Zinc4jet->Fill(genJets[3].pt, genWeight);
                    genFourthJetPt_2_Zinc4jet->Fill(genJets[3].pt, genWeight);
                    genFourthJetEta_Zinc4jet->Fill(fabs(genJets[3].eta), genWeight);
                    genFourthJetEta_2_Zinc4jet->Fill(fabs(genJets[3].eta), genWeight);
                    //*************************************** begin edit *************************************************************//
                    genFourthJetRapidity_Zinc4jet->Fill(fabs(genNewFourthJ.Rapidity()), genWeight);
                    genFourthJetRapidityFull_Zinc4jet->Fill(genNewFourthJ.Rapidity(), genWeight);
                    gendRapidityJets_Zinc4jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight);
                    gendRapidityJets_2_Zinc4jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight);
                    
                    gendRapidityJetsFB_Zinc4jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight);
                    gendRapidityJetsFB_2_Zinc4jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight);
                    
                    gendiJetMass_Zinc4jet->Fill(genJet1Plus2.M(), genWeight);
                    gendiJetMass_2_Zinc4jet->Fill(genJet1Plus2.M(), genWeight);
                    
                    gendiJetPt_Zinc4jet->Fill(genJet1Plus2.Pt(), genWeight);
                    gendiJetPt_2_Zinc4jet->Fill(genJet1Plus2.Pt(), genWeight);
                    //*************************************** end edit ***************************************************************//
                    genFourthJetPtEta_Zinc4jet->Fill(genJets[3].pt, fabs(genJets[3].eta), genWeight);
                    genFirstHighestJetPt_Zinc4jet->Fill(genJets[0].pt, genWeight);
                    genSecondHighestJetPt_Zinc4jet->Fill(genJets[1].pt, genWeight);
                    genThirdHighestJetPt_Zinc4jet->Fill(genJets[2].pt, genWeight);
                    genJetsHT_Zinc4jet->Fill(genJetsHT, genWeight);
                    genJetsHT_1_Zinc4jet->Fill(genJetsHT, genWeight);
                    genJetsHT_2_Zinc4jet->Fill(genJetsHT, genWeight);
                }
                if (nGoodGenJets >= 5){
                    genZNGoodJets_Zinc->Fill(5., genWeight);
                    genZNGoodJetsFull_Zinc->Fill(5., genWeight);
                    genZMass_Zinc5jet->Fill(genZ.M(), genWeight);
                    genZPt_Zinc5jet->Fill(genZ.Pt(), genWeight);
                    genZRapidity_Zinc5jet->Fill(genZ.Rapidity(), genWeight);
                    genZEta_Zinc5jet->Fill(genZ.Eta(), genWeight);
                    genFifthJetPt_Zinc5jet->Fill(genJets[4].pt, genWeight);
                    genFifthJetPt_1_Zinc5jet->Fill(genJets[4].pt, genWeight);
                    genFifthJetPt_2_Zinc5jet->Fill(genJets[4].pt, genWeight);
                    genFifthJetEta_Zinc5jet->Fill(fabs(genJets[4].eta), genWeight);
                    genFifthJetEta_2_Zinc5jet->Fill(fabs(genJets[4].eta), genWeight);
                    genFifthJetPtEta_Zinc5jet->Fill(genJets[4].pt, fabs(genJets[4].eta), genWeight);
                    genFirstHighestJetPt_Zinc5jet->Fill(genJets[0].pt, genWeight);
                    genSecondHighestJetPt_Zinc5jet->Fill(genJets[1].pt, genWeight);
                    genThirdHighestJetPt_Zinc5jet->Fill(genJets[2].pt, genWeight);
                    genJetsHT_Zinc5jet->Fill(genJetsHT, genWeight);
                    genJetsHT_1_Zinc5jet->Fill(genJetsHT, genWeight);
                    genJetsHT_2_Zinc5jet->Fill(genJetsHT, genWeight);
                }
                if (nGoodGenJets >= 6){
                    genZNGoodJets_Zinc->Fill(6., genWeight);
                    genZNGoodJetsFull_Zinc->Fill(6., genWeight);
                    genZMass_Zinc6jet->Fill(genZ.M(), genWeight);
                    genZPt_Zinc6jet->Fill(genZ.Pt(), genWeight);
                    genZRapidity_Zinc6jet->Fill(genZ.Rapidity(), genWeight);
                    genZEta_Zinc6jet->Fill(genZ.Eta(), genWeight);
                    genSixthJetPt_Zinc6jet->Fill(genJets[5].pt, genWeight);
                    genSixthJetPt_1_Zinc6jet->Fill(genJets[5].pt, genWeight);
                    genSixthJetEta_Zinc6jet->Fill(fabs(genJets[5].eta), genWeight);
                    genSixthJetPtEta_Zinc6jet->Fill(genJets[5].pt, fabs(genJets[5].eta), genWeight);
                    genFirstHighestJetPt_Zinc6jet->Fill(genJets[0].pt, genWeight);
                    genSecondHighestJetPt_Zinc6jet->Fill(genJets[1].pt, genWeight);
                    genThirdHighestJetPt_Zinc6jet->Fill(genJets[2].pt, genWeight);
                    genJetsHT_Zinc6jet->Fill(genJetsHT, genWeight);
                    genJetsHT_1_Zinc6jet->Fill(genJetsHT, genWeight);
                }
                if (nGoodGenJets >= 7){
                    genZNGoodJets_Zinc->Fill(7., genWeight);
                    genZNGoodJetsFull_Zinc->Fill(7., genWeight);
                }
                
                if (nGoodGenJets >= 8) genZNGoodJets_Zinc->Fill(8., genWeight);
                if (doW && nGoodGenJets >= 9) genZNGoodJets_Zinc->Fill(9., genWeight);
                if (doW && nGoodGenJets >= 10) genZNGoodJets_Zinc->Fill(10., genWeight);
            }
        }
        
        //=======================================================================================================//

        //=======================================================================================================//
        //      Selection for Reco Histos      //
        //====================================//
        if (hasRecoInfo && passesLeptonCut && passesJetCut && passesDPSPartonCut) {
            //=======================================================================================================//
            
            if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
            //=======================================================================================================//
            //      Start filling histograms      //
            //====================================//
            
            
            TotalRecoWeightPassRECO+=weight;
            TotalGenWeightPassRECO+=genWeightBackup;
            NVtx->Fill(EvtVtxCnt, weight);

            //NVtx->Fill(EvtInfo_NumVtx + 1000 , weight);
            if (fileName.find("Sherpa") != string::npos && fileName.find("UNFOL") == string::npos ) PUWeight->Fill(puWeight.weight(int(PU_npT)) * reweighting * mcEveWeight_, 1);
            else PUWeight->Fill(puWeight.weight(int(PU_npT)) * reweighting, 1);
            if (nGoodJets == 0){
                PUWeight0->Fill(puWeight.weight(int(PU_npT)) * reweighting, 1);
            }
            else {
                PUWeight1->Fill(puWeight.weight(int(PU_npT)) * reweighting, 1);
            }
            
            if (lepton1.charge > 0){
                MuPlusPt->Fill(lepton1.pt, weight);
                MuPlusEta->Fill(lepton1.eta, weight);
                if (doZ || doTT){
                    MuMinusPt->Fill(lepton2.pt, weight);
                    MuMinusEta->Fill(lepton2.eta, weight);
                }
            }
            else {
                MuMinusPt->Fill(lepton1.pt, weight);
                MuMinusEta->Fill(lepton1.eta, weight);
                if (doZ || doTT){
                    MuPlusPt->Fill(lepton2.pt, weight);
                    MuPlusEta->Fill(lepton2.eta, weight);
                }
            }
            nEventsIncl0Jets++;
            ZNGoodJetsNVtx_Zexc->Fill(nGoodJets, EvtVtxCnt  , weight);
            ZNGoodJets_Zinc->Fill(0., weight);
            ZNGoodJetsFull_Zinc->Fill(0., weight);
            ZNGoodJets_Zexc->Fill(nGoodJets, weight);
            ZNGoodJetsFull_Zexc->Fill(nGoodJets, weight);
            ZNGoodJets_Zinc_NoWeight->Fill(0.);
            ZMass_Zinc0jet->Fill(Z.M(), weight);
            MET_Zinc0jet->Fill(METpt, weight);
            //koMETphi_Zinc0jet->Fill(METphi, weight);
            MT_Zinc0jet->Fill(MT, weight);
            
            MuPFIso_Zinc0jet->Fill(lepton1.iso, weight);
            MuPFIso_2ndZinc0jet->Fill(lepton1.iso, weight);
            MuPFIso_3rdZinc0jet->Fill(lepton1.iso, weight);
            
            ZPt_Zinc0jet->Fill(Z.Pt(), weight);
            ZRapidity_Zinc0jet->Fill(Z.Rapidity(), weight);
            ZEta_Zinc0jet->Fill(Z.Eta(), weight);
            lepPt_Zinc0jet->Fill(lepton1.pt, weight);
            lepEta_Zinc0jet->Fill(lepton1.eta, weight);
            lepPhi_Zinc0jet->Fill(lepton1.phi, weight);
            if (doZ || doTT){
                lepPt_Zinc0jet->Fill(lepton2.pt, weight);
                lepEta_Zinc0jet->Fill(lepton2.eta, weight);
                lepPhi_Zinc0jet->Fill(lepton2.phi, weight);
            }
            lepEtaEta_Zinc0jet->Fill(lepton1.eta, lepton2.eta, weight);
            dPhiLeptons_Zinc0jet->Fill(deltaPhi(lep1, lep2), weight);
            dEtaLeptons_Zinc0jet->Fill(lepton1.eta - lepton2.eta, weight);
            dRLeptons_Zinc0jet->Fill(deltaR(lepton1.phi, lepton1.eta, lepton2.phi, lepton2.eta), weight);
            SpTLeptons_Zinc0jet->Fill(SpTsub(lep1, lep2), weight);
            if (nGoodJets == 0){
                nEventsExcl0Jets++;
                ZNGoodJets_Zexc_NoWeight->Fill(0.);
                ZMass_Zexc0jet->Fill(Z.M(), weight);
                ZPt_Zexc0jet->Fill(Z.Pt(), weight);
                ZRapidity_Zexc0jet->Fill(Z.Rapidity(), weight);
                ZEta_Zexc0jet->Fill(Z.Eta(), weight);
                lepPt_Zexc0jet->Fill(lepton1.pt, weight);
                lepEta_Zexc0jet->Fill(lepton1.eta, weight);
                if ( doZ || doTT){
                    lepPt_Zexc0jet->Fill(lepton2.pt, weight);
                    lepEta_Zexc0jet->Fill(lepton2.eta, weight);
                }
                dPhiLeptons_Zexc0jet->Fill(deltaPhi(lep1, lep2), weight);
                dEtaLeptons_Zexc0jet->Fill(lepton1.eta - lepton2.eta, weight);
                SpTLeptons_Zexc0jet->Fill(SpTsub(lep1, lep2), weight);
            }
            
            if (nGoodJets >= 1){
                ZNGoodJets_Zinc->Fill(1., weight);
                ZNGoodJetsFull_Zinc->Fill(1., weight);
                ZNGoodJets_Zinc_NoWeight->Fill(1.);
                ZMass_Zinc1jet->Fill(Z.M(), weight);
                MET_Zinc1jet->Fill(METpt, weight);
                //koMETphi_Zinc1jet->Fill(METphi, weight);
                MT_Zinc1jet->Fill(MT, weight);
                ZPt_Zinc1jet->Fill(Z.Pt(), weight);
                ZRapidity_Zinc1jet->Fill(Z.Rapidity(), weight);
                ZEta_Zinc1jet->Fill(Z.Eta(), weight);
                lepPt_Zinc1jet->Fill(lepton1.pt, weight);
                lepEta_Zinc1jet->Fill(lepton1.eta, weight);
                lepPhi_Zinc1jet->Fill(lepton1.phi, weight);
                if (doZ || doTT){
                    lepPt_Zinc1jet->Fill(lepton2.pt, weight);
                    lepEta_Zinc1jet->Fill(lepton2.eta, weight);
                    lepPhi_Zinc1jet->Fill(lepton2.phi, weight);
                }
                dPhiLeptons_Zinc1jet->Fill(deltaPhi(lep1, lep2), weight);
                dEtaLeptons_Zinc1jet->Fill(lepton1.eta - lepton2.eta, weight);
                dRLeptons_Zinc1jet->Fill(deltaR(lepton1.phi, lepton1.eta, lepton2.phi, lepton2.eta), weight);
                SpTLeptons_Zinc1jet->Fill(SpTsub(lep1, lep2), weight);
                FirstJetPt_Zinc1jet->Fill(jets[0].pt, weight);
                FirstJetPt_1_Zinc1jet->Fill(jets[0].pt, weight);
                FirstJetPt_2_Zinc1jet->Fill(jets[0].pt, weight);
                
                FirstHighestJetPt_Zinc1jet->Fill(jets[0].pt, weight);
                FirstJetEta_Zinc1jet->Fill(fabs(jets[0].eta), weight);
                FirstJetEta_2_Zinc1jet->Fill(fabs(jets[0].eta), weight);
                FirstJetEtaFull_Zinc1jet->Fill(jets[0].eta, weight);
                //*************************************** begin edit *************************************************************//
                FirstJetRapidity_Zinc1jet->Fill(fabs(newLeadJ.Rapidity()), weight);
                FirstJetRapidityFull_Zinc1jet->Fill(newLeadJ.Rapidity(), weight);
                FirstJetmass_Zinc1jet->Fill(newLeadJ.M(), weight);
                FirstJetmass_1_Zinc1jet->Fill(newLeadJ.M(), weight);
                MeanNJetsHT_1D_Zinc1jet->Fill(jetsHT, weight*nGoodJets);
                MeanNJetsHT_Zinc1jet->Fill(jetsHT, nGoodJets, weight);
                //*************************************** end edit ***************************************************************//
                FirstJetPtEta_Zinc1jet->Fill(jets[0].pt, fabs(jets[0].eta), weight);
                FirstJetPhi_Zinc1jet->Fill(jets[0].phi, weight);
                JetsHT_Zinc1jet->Fill(jetsHT, weight);
                JetsHT_1_Zinc1jet->Fill(jetsHT, weight);
                JetsHT_2_Zinc1jet->Fill(jetsHT, weight);
                
                for (unsigned short j(0); j < nGoodJets; j++){
                    AllJetPt_Zinc1jet->Fill(jets[j].pt, weight);
                    AllJetEta_Zinc1jet->Fill(jets[j].eta, weight);
                    AllJetPhi_Zinc1jet->Fill(jets[j].phi, weight);
                }
                if ( doW ) dEtaBosonJet_Zinc1jet->Fill(fabs(jets[0].eta - lepton1.eta), weight);
                else dEtaBosonJet_Zinc1jet->Fill(fabs(jets[0].eta-Z.Eta()), weight);
                for ( int i =0 ; i < NbinsEta2D - 1 ; i++){
                    if ( fabs(jets[0].eta) >= j_Y_range[i] &&  fabs(jets[0].eta) < j_Y_range[i+1] )                                                FirstJetPt_Zinc1jet_Eta[i]->Fill(fabs(jets[0].pt), weight);
                }
                
                if (nGoodJets == 1){
                    // compute Delta pt between Z and jets
                    if (Z.Pt() > 0.8 * jetPtCutMin  && jets[0].pt/Z.Pt() < 1.2 && jets[0].pt/Z.Pt() > 0.8 && deltaPhi(leadJ, Z) > 2.7){
                        hPtEtaBackJet_Zexc1jet->Fill(leadJ.Pt(), leadJ.Eta(), weight);
                        if (JetAk04PuMva->at(jets[0].patIndex) > 0 ) hPtEtaBackJetMVA_Zexc1jet->Fill(leadJ.Pt(), leadJ.Eta(), weight);
                    }
           
                    nEventsExcl1Jets++;
                    ZNGoodJets_Zexc_NoWeight->Fill(1.);
                    ZMass_Zexc1jet->Fill(Z.M(), weight);
                    ZPt_Zexc1jet->Fill(Z.Pt(), weight);
                    ZRapidity_Zexc1jet->Fill(Z.Rapidity(), weight);
                    ZEta_Zexc1jet->Fill(Z.Eta(), weight);
                    lepPt_Zexc1jet->Fill(lepton1.pt, weight);
                    lepEta_Zexc1jet->Fill(lepton1.eta, weight);
                    if (doZ || doTT){
                        lepPt_Zexc1jet->Fill(lepton2.pt, weight);
                        lepEta_Zexc1jet->Fill(lepton2.eta, weight);
                    }
                    dPhiLeptons_Zexc1jet->Fill(deltaPhi(lep1, lep2), weight);
                    dEtaLeptons_Zexc1jet->Fill(lepton1.eta - lepton2.eta, weight);
                    SpTLeptons_Zexc1jet->Fill(SpTsub(lep1, lep2), weight);
                    FirstJetPt_Zexc1jet->Fill(jets[0].pt, weight);
                    FirstJetEta_Zexc1jet->Fill(jets[0].eta, weight);
                    FirstJetPhi_Zexc1jet->Fill(jets[0].phi, weight);
                    if ( doW ) dEtaBosonJet_Zexc1jet->Fill(fabs(jets[0].eta - lepton1.eta), weight);
                    else dEtaBosonJet_Zexc1jet->Fill(fabs(jets[0].eta-Z.Eta()), weight);
                    
                }
            }
            if (nGoodJets >= 2){
                TLorentzVector jet1Plus2PlusZ = jet1Plus2 + Z;
                ZNGoodJets_Zinc->Fill(2., weight);
                ZNGoodJetsFull_Zinc->Fill(2., weight);
                ZNGoodJets_Zinc_NoWeight->Fill(2.);
                ZMass_Zinc2jet->Fill(Z.M(), weight);
                MET_Zinc2jet->Fill(METpt, weight);
                //koMETphi_Zinc2jet->Fill(METphi, weight);
                MT_Zinc2jet->Fill(MT, weight);
                TwoJetsPtDiff_Zinc2jet->Fill(jet1Minus2.Pt(), weight);
                BestTwoJetsPtDiff_Zinc2jet->Fill(bestJet1Minus2.Pt(), weight);
                JetsMass_Zinc2jet->Fill(jet1Plus2.M(), weight);
                ZPt_Zinc2jet->Fill(Z.Pt(), weight);
                ZRapidity_Zinc2jet->Fill(Z.Rapidity(), weight);
                ZEta_Zinc2jet->Fill(Z.Eta(), weight);
                lepPt_Zinc2jet->Fill(lepton1.pt, weight);
                lepEta_Zinc2jet->Fill(lepton1.eta, weight);
                if (doZ || doTT){
                    lepPt_Zinc2jet->Fill(lepton2.pt, weight);
                    lepEta_Zinc2jet->Fill(lepton2.eta, weight);
                }
                dPhiLeptons_Zinc2jet->Fill(deltaPhi(lep1, lep2), weight);
                dEtaLeptons_Zinc2jet->Fill(lepton1.eta - lepton2.eta, weight);
                dRLeptons_Zinc2jet->Fill(deltaR(lepton1.phi, lepton1.eta, lepton2.phi, lepton2.eta), weight);
                SpTLeptons_Zinc2jet->Fill(SpTsub(lep1, lep2), weight);
                FirstHighestJetPt_Zinc2jet->Fill(jets[0].pt, weight);
                SecondHighestJetPt_Zinc2jet->Fill(jets[1].pt, weight);
                SecondJetPt_Zinc2jet->Fill(jets[1].pt, weight);
                SecondJetPt_1_Zinc2jet->Fill(jets[1].pt, weight);
                SecondJetPt_2_Zinc2jet->Fill(jets[1].pt, weight);
                SecondJetEta_Zinc2jet->Fill(fabs(jets[1].eta), weight);
                SecondJetEta_2_Zinc2jet->Fill(fabs(jets[1].eta), weight);
                
                SecondJetEtaFull_Zinc2jet->Fill(jets[1].eta, weight);
                //*************************************** begin edit *************************************************************//
                SecondJetRapidity_Zinc2jet->Fill(fabs(newSecondJ.Rapidity()), weight);
                SecondJetRapidityFull_Zinc2jet->Fill(newSecondJ.Rapidity(), weight);
                SecondJetmass_Zinc2jet->Fill(newSecondJ.M(), weight);
                SecondJetmass_1_Zinc2jet->Fill(newSecondJ.M(), weight);
                dRapidityJets_Zinc2jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight);
                dRapidityJets_2_Zinc2jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight);
                
                dRapidityJetsFB_Zinc2jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight);
                dRapidityJetsFB_2_Zinc2jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight);
                
                diJetMass_Zinc2jet->Fill(jet1Plus2.M(), weight);
                diJetMass_2_Zinc2jet->Fill(jet1Plus2.M(), weight);
                
                dPhiJetsFB_Zinc2jet->Fill(deltaPhi(vJetYOrdered[0], vJetYOrdered[vJetYOrdered.size() - 1] ), weight);
                dPhiJetsFB_2_Zinc2jet->Fill(deltaPhi(vJetYOrdered[0], vJetYOrdered[vJetYOrdered.size() - 1] ), weight);
                
                dRJets_Zinc2jet->Fill(deltaRYPhi(newLeadJ, newSecondJ), weight);
                dRJets_2_Zinc2jet->Fill(deltaRYPhi(newLeadJ, newSecondJ), weight);
                
                diJetPt_Zinc2jet->Fill(jet1Plus2.Pt(), weight);
                diJetPt_2_Zinc2jet->Fill(jet1Plus2.Pt(), weight);
                
                MeanNJetsHT_1D_Zinc2jet->Fill(jetsHT, weight*nGoodJets);
                MeanNJetsdRapidity_1D_Zinc2jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight*nGoodJets);
                MeanNJetsdRapidityFB_1D_Zinc2jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight*nGoodJets);
                
                MeanNJetsHT_Zinc2jet->Fill(jetsHT, nGoodJets, weight);
                MeanNJetsdRapidity_Zinc2jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), nGoodJets, weight);
                MeanNJetsdRapidityFB_Zinc2jet->Fill(ForwardJetRapidity - BackwardJetRapidity, nGoodJets, weight);
                //*************************************** end edit ***************************************************************//
                SecondJetPtEta_Zinc2jet->Fill(jets[1].pt, fabs(jets[1].eta), weight);
                RatioJetPt21_Zinc2jet->Fill(jets[1].pt/jets[0].pt, weight);
                SecondJetPhi_Zinc2jet->Fill(jets[1].phi, weight);
                JetsHT_Zinc2jet->Fill(jetsHT, weight);
                JetsHT_1_Zinc2jet->Fill(jetsHT, weight);
                JetsHT_2_Zinc2jet->Fill(jetsHT, weight);
                
                ptBal_Zinc2jet->Fill(jet1Plus2PlusZ.Pt(), weight);
                
                dPhiJets_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                dPhiJets_2_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                
                BestdPhiJets_Zinc2jet->Fill(deltaPhi(bestTwoJets.first, bestTwoJets.second), weight);
                dEtaJets_Zinc2jet->Fill(leadJ.Eta() - secondJ.Eta(), weight);
                dEtaFirstJetZ_Zinc2jet->Fill(leadJ.Eta() - Z.Eta(), weight);
                dEtaSecondJetZ_Zinc2jet->Fill(secondJ.Eta() - Z.Eta(), weight);
                dEtaJet1Plus2Z_Zinc2jet->Fill(jet1Plus2.Eta() - Z.Eta(), weight);
                PHI_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                BestPHI_Zinc2jet->Fill(PHI(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                PHI_T_Zinc2jet->Fill(PHI_T(lep1, lep2, leadJ, secondJ), weight);
                BestPHI_T_Zinc2jet->Fill(PHI_T(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                SpT_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                BestSpT_Zinc2jet->Fill(SpT(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                SpTJets_Zinc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                BestSpTJets_Zinc2jet->Fill(SpTsub(bestTwoJets.first, bestTwoJets.second), weight);
                SPhi_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                BestSPhi_Zinc2jet->Fill(SPhi(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                
                
                for ( int i =0 ; i < NbinsEta2D - 1 ; i++){
                    if ( fabs(jets[1].eta) >= j_Y_range[i] &&  fabs(jets[1].eta) < j_Y_range[i+1]                                   )                                                SecondJetPt_Zinc2jet_Eta[i]->Fill(fabs(jets[0].pt), weight);
                }
                //--- V + 2 jets EWK histograms
                if (passesEWKJetPt){
                    // met histos
                    METEWK_Zinc2jet->Fill(METpt, weight);
                    //koMETphiEWK_Zinc2jet->Fill(METphi,  weight);
                    MTEWK_Zinc2jet->Fill(MT,  weight);
                    
                    // jet histos
                    ZNGoodJetsEWK_Zexc->Fill(nGoodJets, weight);
                    FirstJetPtEWK_Zinc2jet->Fill(jets[0].pt, weight);
                    FirstJetEtaEWK_Zinc2jet->Fill(jets[0].eta, weight);
                    FirstJetPhiEWK_Zinc2jet->Fill(jets[0].phi, weight);
                    
                    
                    SecondJetPtEWK_Zinc2jet->Fill(jets[1].pt, weight);
                    SecondJetEtaEWK_Zinc2jet->Fill(jets[1].eta, weight);
                    SecondJetPhiEWK_Zinc2jet->Fill(jets[1].phi, weight);
                    
                    int temIND (0), temIND1(1);
                    if (fabs(jets[1].eta) > fabs(jets[0].eta)){
                        temIND = 1;
                        temIND1 = 0;
                    }
                    ForwardJetPtEWK_Zinc2jet->Fill(jets[temIND].pt, weight);
                    ForwardJetEtaEWK_Zinc2jet->Fill(jets[temIND].eta, weight);
                    ForwardJetPhiEWK_Zinc2jet->Fill(jets[temIND].phi, weight);
                    
                    CentralJetPtEWK_Zinc2jet->Fill(jets[temIND1].pt, weight);
                    CentralJetEtaEWK_Zinc2jet->Fill(jets[temIND1].eta, weight);
                    CentralJetPhiEWK_Zinc2jet->Fill(jets[temIND1].phi, weight);
                    
                    JetsHTEWK_Zinc2jet->Fill(jetsHT, weight);
                    JetsMassEWK_Zinc2jet->Fill(jet1Plus2.M(), weight);
                    
                    // multi jet variables
                    SumEtaJetsEWK_Zinc2jet->Fill(fabs(leadJ.Eta() + secondJ.Eta()),  weight);
                    AbsSumEtaJetsEWK_Zinc2jet->Fill(fabs(leadJ.Eta()) + fabs(secondJ.Eta()),  weight);
                    TwoJetsPtDiffEWK_Zinc2jet->Fill(jet1Minus2.Pt(), weight);
                    ptBalEWK_Zinc2jet->Fill(jet1Plus2PlusZ.Pt(), weight);
                    dPhiJetsEWK_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                    dEtaJetsEWK_Zinc2jet->Fill(fabs(leadJ.Eta() - secondJ.Eta()), weight);
                    SpTJetsEWK_Zinc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                    dPhiJetsEWK_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                    
                    // find jet properties of the third jet that is between the two leading jets
                    int nGoodJetsBtw(0.);
                    double jetsHTBtw(0.);
                    if (nJetsAdd > 2){
                        ThirdJetPtEWKadd_Zinc2jet->Fill(jetsAdditional[2].pt, weight);
                        ThirdJetEtaEWKadd_Zinc2jet->Fill(jetsAdditional[2].eta, weight);
                        ThirdJetPhiEWKadd_Zinc2jet->Fill(jetsAdditional[2].phi, weight);
                        
                        for (unsigned short i(2); i < nJetsAdd; i++){
                            int coutBtw = 0 ;
                            if (jetsAdditional[i].eta < max(jets[0].eta,jets[1].eta) - 0.5
                                && jetsAdditional[i].eta > min(jets[0].eta,jets[1].eta) + 0.5){
                                jetsHTBtw += jetsAdditional[i].pt ;
                                nGoodJetsBtw++;
                                AllJetPtEWKbtw_Zinc2jet->Fill(jetsAdditional[i].pt, weight);
                                AllJetEtaEWKbtw_Zinc2jet->Fill(jetsAdditional[i].eta, weight);
                                AllJetPhiEWKbtw_Zinc2jet->Fill(jetsAdditional[i].phi, weight);
                                if (coutBtw == 0){
                                    ThirdJetPtEWKbtw_Zinc2jet->Fill(jetsAdditional[i].pt, weight);
                                    ThirdJetEtaEWKbtw_Zinc2jet->Fill(jetsAdditional[i].eta, weight);
                                    ThirdJetPhiEWKbtw_Zinc2jet->Fill(jetsAdditional[i].phi, weight);
                                    coutBtw = 1;
                                }
                            }
                        }
                    }
                    ZNGoodJetsEWKbtw_Zexc->Fill(nGoodJetsBtw, weight);
                    JetsHTEWKbtw_Zinc2jet->Fill(jetsHTBtw, weight);
                    
                    // at least one forward jet
                    if (passesEWKJetFwdEta){
                        METEWKfwd_Zinc2jet->Fill(METpt,  weight);
                        //KoMETphiEWKfwd_Zinc2jet->Fill(METphi,  weight);
                        MTEWKfwd_Zinc2jet->Fill(MT,  weight);
                        
                        // jet histos
                        ZNGoodJetsEWKfwd_Zexc->Fill(nGoodJets, weight);
                        FirstJetPtEWKfwd_Zinc2jet->Fill(jets[0].pt, weight);
                        FirstJetEtaEWKfwd_Zinc2jet->Fill(jets[0].eta, weight);
                        FirstJetPhiEWKfwd_Zinc2jet->Fill(jets[0].phi, weight);
                        
                        SecondJetPtEWKfwd_Zinc2jet->Fill(jets[1].pt, weight);
                        SecondJetEtaEWKfwd_Zinc2jet->Fill(jets[1].eta, weight);
                        SecondJetPhiEWKfwd_Zinc2jet->Fill(jets[1].phi, weight);
                        
                        JetsHTEWKfwd_Zinc2jet->Fill(jetsHT, weight);
                        JetsMassEWKfwd_Zinc2jet->Fill(jet1Plus2.M(), weight);
                        
                        // multi jet variables
                        SumEtaJetsEWKfwd_Zinc2jet->Fill(leadJ.Eta() + secondJ.Eta(), weight);
                        TwoJetsPtDiffEWKfwd_Zinc2jet->Fill(jet1Minus2.Pt(), weight);
                        ptBalEWKfwd_Zinc2jet->Fill(jet1Plus2PlusZ.Pt(), weight);
                        dPhiJetsEWKfwd_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                        dEtaJetsEWKfwd_Zinc2jet->Fill(fabs(leadJ.Eta() - secondJ.Eta()), weight);
                        SpTJetsEWKfwd_Zinc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                        dPhiJetsEWKfwd_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                        
                        
                    }//--- end at least one forward jet
                    
                    if (jet1Plus2.M() > 1000.){
                        METEWKmjj_Zinc2jet->Fill(METpt,  weight);
                        //KoMETphiEWKmjj_Zinc2jet->Fill(METphi,  weight);
                        MTEWKmjj_Zinc2jet->Fill(MT,  weight);
                        short nGoodJetsAdd(nJetsAdd-2);
                        double jetsHTAdd(0);
                        // jet histos
                        ZNGoodJetsEWKmjj_Zexc->Fill(nGoodJetsAdd, weight);
                        FirstJetPtEWKmjj_Zinc2jet->Fill(jets[0].pt, weight);
                        FirstJetEtaEWKmjj_Zinc2jet->Fill(jets[0].eta, weight);
                        FirstJetPhiEWKmjj_Zinc2jet->Fill(jets[0].phi, weight);
                        
                        
                        SecondJetPtEWKmjj_Zinc2jet->Fill(jets[1].pt, weight);
                        SecondJetEtaEWKmjj_Zinc2jet->Fill(jets[1].eta, weight);
                        SecondJetPhiEWKmjj_Zinc2jet->Fill(jets[1].phi, weight);
                        
                        JetsHTEWKmjj_Zinc2jet->Fill(jetsHT, weight);
                        JetsMassEWKmjj_Zinc2jet->Fill(jet1Plus2.M(), weight);
                        
                        // multi jet variables
                        SumEtaJetsEWKmjj_Zinc2jet->Fill(leadJ.Eta() + secondJ.Eta(),  weight);
                        TwoJetsPtDiffEWKmjj_Zinc2jet->Fill(jet1Minus2.Pt(), weight);
                        ptBalEWKmjj_Zinc2jet->Fill(jet1Plus2PlusZ.Pt(), weight);
                        dPhiJetsEWKmjj_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                        dEtaJetsEWKmjj_Zinc2jet->Fill(leadJ.Eta() - secondJ.Eta(), weight);
                        SpTJetsEWKmjj_Zinc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                        dPhiJetsEWKmjj_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                        
                        // find jet properties of the third jet that is  between the two leading jets
                        int nGoodJetsmjjBtw(0.);
                        if (nJetsAdd > 2){
                            double jetsHTmjjBtw(0.);
                            for (unsigned short i(2); i < nJetsAdd; i++){
                                if (jetsAdditional[i].eta < max(jets[0].eta,jets[1].eta) - 0.5
                                    && jetsAdditional[i].eta > min(jets[0].eta,jets[1].eta) + 0.5 ){
                                    jetsHTmjjBtw += jetsAdditional[i].pt ;
                                    nGoodJetsmjjBtw++;
                                }
                                jetsHTAdd += jetsAdditional[i].pt ;
                            }
                            
                            ThirdJetPtEWKmjj_Zinc3jet->Fill(jetsAdditional[2].pt, weight);
                            JetsHTEWKmjjAdd_Zinc2jet->Fill(jetsHTAdd, weight);
                            TLorentzVector thirdJAdd;
                            thirdJAdd.SetPtEtaPhiE(jetsAdditional[2].pt, jetsAdditional[2].eta, jetsAdditional[2].phi, jetsAdditional[2].energy);
                            double tempRapidiy3Jet = thirdJAdd.Rapidity() - 0.5 * (leadJ.Rapidity() + secondJ.Rapidity());
                            ThirdJetEtaEWKmjj_Zinc3jet->Fill(tempRapidiy3Jet, weight);
                        }
                    }//--- end dijet mass > 1000
                    
                    //--- higher jet properties
                    if (nGoodJets > 2){
                        METEWK_Zinc3jet->Fill(METpt, weight);
                        //koMETphiEWK_Zinc3jet->Fill(METphi, weight);
                        MTEWK_Zinc3jet->Fill(MT, weight);
                        
                        TLorentzVector thirdJ;
                        thirdJ.SetPtEtaPhiE(jets[2].pt, jets[2].eta, jets[2].phi, jets[2].energy);
                        double tempRapidiy3Jet = thirdJ.Rapidity() - 0.5 * (leadJ.Rapidity() + secondJ.Rapidity());
                        EtaThirdJetsEWK_Zinc3jet->Fill(tempRapidiy3Jet, weight);
                    }
                }
                
                for (unsigned short j(0); j < nGoodJets; j++){
                    AllJetPt_Zinc2jet->Fill(jets[j].pt, weight);
                    AllJetEta_Zinc2jet->Fill(jets[j].eta, weight);
                    AllJetPhi_Zinc2jet->Fill(jets[j].phi, weight);
                }
                if (Z.Pt() < 25){
                    ptBal_LowPt_Zinc2jet->Fill(jet1Plus2PlusZ.Pt(), weight);
                    dPhiJets_LowPt_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                    BestdPhiJets_LowPt_Zinc2jet->Fill(deltaPhi(bestTwoJets.first, bestTwoJets.second), weight);
                    dPhiLeptons_LowPt_Zinc2jet->Fill(deltaPhi(lep1, lep2), weight);
                    PHI_LowPt_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                    BestPHI_LowPt_Zinc2jet->Fill(PHI(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                    PHI_T_LowPt_Zinc2jet->Fill(PHI_T(lep1, lep2, leadJ, secondJ), weight);
                    BestPHI_T_LowPt_Zinc2jet->Fill(PHI_T(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                    SpT_LowPt_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    BestSpT_LowPt_Zinc2jet->Fill(SpT(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                    SpTJets_LowPt_Zinc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                    BestSpTJets_LowPt_Zinc2jet->Fill(SpTsub(bestTwoJets.first, bestTwoJets.second), weight);
                    SpTLeptons_LowPt_Zinc2jet->Fill(SpTsub(lep1, lep2), weight);
                    SPhi_LowPt_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    BestSPhi_LowPt_Zinc2jet->Fill(SPhi(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                    if (SpT(lep1, lep2, leadJ, secondJ) < 0.5){
                        PHI_LowSpT_LowPt_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        SPhi_LowSpT_LowPt_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    }
                    else {
                        PHI_HighSpT_LowPt_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        SPhi_HighSpT_LowPt_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    }
                    if (SPhi(lep1, lep2, leadJ, secondJ) < 0.5){
                        SpT_LowSPhi_LowPt_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    }
                    else {
                        SpT_HighSPhi_LowPt_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    }
                }
                else {
                    ptBal_HighPt_Zinc2jet->Fill(jet1Plus2PlusZ.Pt(),weight);
                    dPhiJets_HighPt_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                    dPhiLeptons_HighPt_Zinc2jet->Fill(deltaPhi(lep1, lep2), weight);
                    PHI_HighPt_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                    PHI_T_HighPt_Zinc2jet->Fill(PHI_T(lep1, lep2, leadJ, secondJ), weight);
                    SpT_HighPt_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    SpTJets_HighPt_Zinc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                    SpTLeptons_HighPt_Zinc2jet->Fill(SpTsub(lep1, lep2), weight);
                    SPhi_HighPt_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    if (SpT(lep1, lep2, leadJ, secondJ) < 0.5){
                        PHI_LowSpT_HighPt_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        SPhi_LowSpT_HighPt_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    }
                    else {
                        PHI_HighSpT_HighPt_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        SPhi_HighSpT_HighPt_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    }
                    if (SPhi(lep1, lep2, leadJ, secondJ) < 0.5){
                        SpT_LowSPhi_HighPt_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    }
                    else {
                        SpT_HighSPhi_HighPt_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    }
                }
                if (SPhi(lep1, lep2, leadJ, secondJ) < 0.5){
                    SpT_LowSPhi_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                }
                else {
                    SpT_HighSPhi_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                }
                if (SpT(lep1, lep2, leadJ, secondJ) < 0.5){
                    PHI_LowSpT_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                    SPhi_LowSpT_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                }
                else {
                    PHI_HighSpT_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                    SPhi_HighSpT_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                }
                if (nGoodJets == 2){
                    nEventsExcl2Jets++;
                    ZNGoodJets_Zexc_NoWeight->Fill(2.);
                    ZMass_Zexc2jet->Fill(Z.M(), weight);
                    ZPt_Zexc2jet->Fill(Z.Pt(), weight);
                    ZRapidity_Zexc2jet->Fill(Z.Rapidity(), weight);
                    ZEta_Zexc2jet->Fill(Z.Eta(), weight);
                    lepPt_Zexc2jet->Fill(lepton1.pt, weight);
                    lepEta_Zexc2jet->Fill(lepton1.eta, weight);
                    if (doZ || doTT){
                        lepPt_Zexc2jet->Fill(lepton2.pt, weight);
                        lepEta_Zexc2jet->Fill(lepton2.eta, weight);
                    }
                    dPhiLeptons_Zexc2jet->Fill(deltaPhi(lep1, lep2), weight);
                    dEtaLeptons_Zexc2jet->Fill(lepton1.eta - lepton2.eta, weight);
                    SpTLeptons_Zexc2jet->Fill(SpTsub(lep1, lep2), weight);
                    SecondJetPt_Zexc2jet->Fill(jets[1].pt, weight);
                    SecondJetEta_Zexc2jet->Fill(jets[1].eta, weight);
                    SecondJetPhi_Zexc2jet->Fill(jets[1].phi, weight);
                    
                    //-- DPS Histograms
                    TwoJetsPtDiff_Zexc2jet->Fill(jet1Minus2.Pt(), weight);
                    JetsMass_Zexc2jet->Fill(jet1Plus2.M(), weight);
                    ptBal_Zexc2jet->Fill(jet1Plus2PlusZ.Pt(), weight);
                    dPhiJets_Zexc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                    dEtaJets_Zexc2jet->Fill(leadJ.Eta() - secondJ.Eta(), weight);
                    dEtaFirstJetZ_Zexc2jet->Fill(leadJ.Eta() - Z.Eta(), weight);
                    dEtaSecondJetZ_Zexc2jet->Fill(secondJ.Eta() - Z.Eta(), weight);
                    dEtaJet1Plus2Z_Zexc2jet->Fill(jet1Plus2.Eta() - Z.Eta(), weight);
                    PHI_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                    PHI_T_Zexc2jet->Fill(PHI_T(lep1, lep2, leadJ, secondJ), weight);
                    SpT_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    SpTJets_Zexc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                    SPhi_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    if (Z.Pt() < 25){
                        ptBal_LowPt_Zexc2jet->Fill(jet1Plus2PlusZ.Pt(), weight);
                        dPhiJets_LowPt_Zexc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                        dPhiLeptons_LowPt_Zexc2jet->Fill(deltaPhi(lep1, lep2), weight);
                        PHI_LowPt_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        PHI_T_LowPt_Zexc2jet->Fill(PHI_T(lep1, lep2, leadJ, secondJ), weight);
                        SpT_LowPt_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                        SpTJets_LowPt_Zexc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                        SpTLeptons_LowPt_Zexc2jet->Fill(SpTsub(lep1, lep2), weight);
                        SPhi_LowPt_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                        if (SpT(lep1, lep2, leadJ, secondJ) < 0.5){
                            PHI_LowSpT_LowPt_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                            SPhi_LowSpT_LowPt_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                        }
                        else {
                            PHI_HighSpT_LowPt_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                            SPhi_HighSpT_LowPt_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                        }
                        if (SPhi(lep1, lep2, leadJ, secondJ) < 0.5){
                            SpT_LowSPhi_LowPt_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                        }
                        else {
                            SpT_HighSPhi_LowPt_Zexc2jet ->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                        }
                    }
                    else {
                        ptBal_HighPt_Zexc2jet->Fill(jet1Plus2PlusZ.Pt(),weight);
                        dPhiJets_HighPt_Zexc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                        dPhiLeptons_HighPt_Zexc2jet->Fill(deltaPhi(lep1, lep2), weight);
                        PHI_HighPt_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        PHI_T_HighPt_Zexc2jet->Fill(PHI_T(lep1, lep2, leadJ, secondJ), weight);
                        SpT_HighPt_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                        SpTJets_HighPt_Zexc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                        SpTLeptons_HighPt_Zexc2jet->Fill(SpTsub(lep1, lep2), weight);
                        SPhi_HighPt_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                        if (SpT(lep1, lep2, leadJ, secondJ) < 0.5){
                            PHI_LowSpT_HighPt_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                            SPhi_LowSpT_HighPt_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                        }
                        else {
                            PHI_HighSpT_HighPt_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                            SPhi_HighSpT_HighPt_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                        }
                        if (SPhi(lep1, lep2, leadJ, secondJ) < 0.5){
                            SpT_LowSPhi_HighPt_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                        }
                        else {
                            SpT_HighSPhi_HighPt_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                        }
                    }
                    if (SPhi(lep1, lep2, leadJ, secondJ) < 0.5){
                        SpT_LowSPhi_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    }
                    else {
                        SpT_HighSPhi_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    }
                    if (SpT(lep1, lep2, leadJ, secondJ) < 0.5){
                        PHI_LowSpT_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        SPhi_LowSpT_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    }
                    else {
                        PHI_HighSpT_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        SPhi_HighSpT_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    }
                }
            }
            if (nGoodJets >= 3) {
                ZNGoodJets_Zinc->Fill(3., weight);
                ZNGoodJetsFull_Zinc->Fill(3., weight);
                ZNGoodJets_Zinc_NoWeight->Fill(3.);
                ZMass_Zinc3jet->Fill(Z.M(), weight);
                MET_Zinc3jet->Fill(METpt, weight);
                //koMETphi_Zinc3jet->Fill(METphi, weight);
                MT_Zinc3jet->Fill(MT, weight);
                ZPt_Zinc3jet->Fill(Z.Pt(), weight);
                ZRapidity_Zinc3jet->Fill(Z.Rapidity(), weight);
                ZEta_Zinc3jet->Fill(Z.Eta(), weight);
                lepPt_Zinc3jet->Fill(lepton1.pt, weight);
                lepEta_Zinc3jet->Fill(lepton1.eta, weight);
                if ( doZ || doTT){
                    lepPt_Zinc3jet->Fill(lepton2.pt, weight);
                    lepEta_Zinc3jet->Fill(lepton2.eta, weight);
                }
                dPhiLeptons_Zinc3jet->Fill(deltaPhi(lep1, lep2), weight);
                dEtaLeptons_Zinc3jet->Fill(lepton1.eta - lepton2.eta, weight);
                dRLeptons_Zinc3jet->Fill(deltaR(lepton1.phi, lepton1.eta, lepton2.phi, lepton2.eta), weight);
                SpTLeptons_Zinc3jet->Fill(SpTsub(lep1, lep2), weight);
                ThirdJetPt_Zinc3jet->Fill(jets[2].pt, weight);
                ThirdJetPt_1_Zinc3jet->Fill(jets[2].pt, weight);
                ThirdJetPt_2_Zinc3jet->Fill(jets[2].pt, weight);
                
                FirstHighestJetPt_Zinc3jet->Fill(jets[0].pt, weight);
                SecondHighestJetPt_Zinc3jet->Fill(jets[1].pt, weight);
                ThirdHighestJetPt_Zinc3jet->Fill(jets[2].pt, weight);
                ThirdJetEta_Zinc3jet->Fill(fabs(jets[2].eta), weight);
                ThirdJetEta_2_Zinc3jet->Fill(fabs(jets[2].eta), weight);
                ThirdJetEtaFull_Zinc3jet->Fill(jets[2].eta, weight);
                //*************************************** begin edit *************************************************************//
                ThirdJetRapidity_Zinc3jet->Fill(fabs(newThirdJ.Rapidity()), weight);
                ThirdJetRapidityFull_Zinc3jet->Fill(newThirdJ.Rapidity(), weight);
                ThirdJetmass_Zinc3jet->Fill(newThirdJ.M(), weight);
                ThirdJetmass_1_Zinc3jet->Fill(newThirdJ.M(), weight);
                dRapidityJets_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight);
                dRapidityJets_2_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight);
                
                dRapidityJets_First_Third_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newThirdJ.Rapidity()), weight);
                dRapidityJets_2_First_Third_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newThirdJ.Rapidity()), weight);
                
                dRapidityJets_Second_Third_Zinc3jet->Fill(fabs(newSecondJ.Rapidity() - newThirdJ.Rapidity()), weight);
                dRapidityJets_2_Second_Third_Zinc3jet->Fill(fabs(newSecondJ.Rapidity() - newThirdJ.Rapidity()), weight);
                
                dRapidityJetsFB_Zinc3jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight);
                dRapidityJetsFB_2_Zinc3jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight);
                
                diJetMass_Zinc3jet->Fill(jet1Plus2.M(), weight);
                diJetMass_2_Zinc3jet->Fill(jet1Plus2.M(), weight);
                
                diJetPt_Zinc3jet->Fill(jet1Plus2.Pt(), weight);
                diJetPt_2_Zinc3jet->Fill(jet1Plus2.Pt(), weight);
                //*************************************** end edit ***************************************************************//
                ThirdJetPtEta_Zinc3jet->Fill(jets[2].pt, fabs(jets[2].eta), weight);
                RatioJetPt32_Zinc3jet->Fill(jets[2].pt/jets[1].pt, weight);
                ThirdJetPhi_Zinc3jet->Fill(jets[2].phi, weight);
                JetsHT_Zinc3jet->Fill(jetsHT, weight);
                JetsHT_1_Zinc3jet->Fill(jetsHT, weight);
                JetsHT_2_Zinc3jet->Fill(jetsHT, weight);
                
                for (unsigned short j(0); j < nGoodJets; j++){
                    AllJetPt_Zinc3jet->Fill(jets[j].pt, weight);
                    AllJetEta_Zinc3jet->Fill(jets[j].eta, weight);
                    AllJetPhi_Zinc3jet->Fill(jets[j].phi, weight);
                }
                if (nGoodJets == 3){
                    nEventsExcl3Jets++;
                    ZNGoodJets_Zexc_NoWeight->Fill(3.);
                    ZMass_Zexc3jet->Fill(Z.M(), weight);
                    ZPt_Zexc3jet->Fill(Z.Pt(), weight);
                    ZRapidity_Zexc3jet->Fill(Z.Rapidity(), weight);
                    ZEta_Zexc3jet->Fill(Z.Eta(), weight);
                    lepPt_Zexc3jet->Fill(lepton1.pt, weight);
                    lepEta_Zexc3jet->Fill(lepton1.eta, weight);
                    if (doZ || doTT){
                        lepPt_Zexc3jet->Fill(lepton2.pt, weight);
                        lepEta_Zexc3jet->Fill(lepton2.eta, weight);
                    }
                    dPhiLeptons_Zexc3jet->Fill(deltaPhi(lep1, lep2), weight);
                    dEtaLeptons_Zexc3jet->Fill(lepton1.eta - lepton2.eta, weight);
                    SpTLeptons_Zexc3jet->Fill(SpTsub(lep1, lep2), weight);
                }
            }
            if (nGoodJets >= 4){
                ZNGoodJets_Zinc->Fill(4., weight);
                ZNGoodJetsFull_Zinc->Fill(4., weight);
                ZNGoodJets_Zinc_NoWeight->Fill(4.);
                ZMass_Zinc4jet->Fill(Z.M(), weight);
                ZPt_Zinc4jet->Fill(Z.Pt(), weight);
                ZRapidity_Zinc4jet->Fill(Z.Rapidity(), weight);
                ZEta_Zinc4jet->Fill(Z.Eta(), weight);
                lepPt_Zinc4jet->Fill(lepton1.pt, weight);
                lepEta_Zinc4jet->Fill(lepton1.eta, weight);
                if (doZ || doTT){
                    lepEta_Zinc4jet->Fill(lepton2.eta, weight);
                    lepPt_Zinc4jet->Fill(lepton2.pt, weight);
                }
                dPhiLeptons_Zinc4jet->Fill(deltaPhi(lep1, lep2), weight);
                dEtaLeptons_Zinc4jet->Fill(lepton1.eta - lepton2.eta, weight);
                dRLeptons_Zinc4jet->Fill(deltaR(lepton1.phi, lepton1.eta, lepton2.phi, lepton2.eta), weight);
                SpTLeptons_Zinc4jet->Fill(SpTsub(lep1, lep2), weight);
                FourthJetPt_Zinc4jet->Fill(jets[3].pt, weight);
                FourthJetPt_1_Zinc4jet->Fill(jets[3].pt, weight);
                FourthJetPt_2_Zinc4jet->Fill(jets[3].pt, weight);
                FirstHighestJetPt_Zinc4jet->Fill(jets[0].pt, weight);
                SecondHighestJetPt_Zinc4jet->Fill(jets[1].pt, weight);
                ThirdHighestJetPt_Zinc4jet->Fill(jets[2].pt, weight);
                FourthJetEta_Zinc4jet->Fill(fabs(jets[3].eta), weight);
                FourthJetEta_2_Zinc4jet->Fill(fabs(jets[3].eta), weight);
                FourthJetEtaFull_Zinc4jet->Fill(jets[3].eta, weight);
                //*************************************** begin edit *************************************************************//
                FourthJetRapidity_Zinc4jet->Fill(fabs(newFourthJ.Rapidity()), weight);
                FourthJetRapidityFull_Zinc4jet->Fill(newFourthJ.Rapidity(), weight);
                FourthJetmass_Zinc4jet->Fill(newFourthJ.M(), weight);
                FourthJetmass_1_Zinc4jet->Fill(newFourthJ.M(), weight);
                dRapidityJets_Zinc4jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight);
                dRapidityJets_2_Zinc4jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight);
                
                dRapidityJetsFB_Zinc4jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight);
                dRapidityJetsFB_2_Zinc4jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight);
                
                diJetMass_Zinc4jet->Fill(jet1Plus2.M(), weight);
                diJetMass_2_Zinc4jet->Fill(jet1Plus2.M(), weight);
                
                diJetPt_Zinc4jet->Fill(jet1Plus2.Pt(), weight);
                diJetPt_2_Zinc4jet->Fill(jet1Plus2.Pt(), weight);
                //*************************************** end edit ***************************************************************//
                FourthJetPtEta_Zinc4jet->Fill(jets[3].pt, fabs(jets[3].eta), weight);
                FourthJetPhi_Zinc4jet->Fill(jets[3].phi, weight);
                JetsHT_Zinc4jet->Fill(jetsHT, weight);
                JetsHT_1_Zinc4jet->Fill(jetsHT, weight);
                JetsHT_2_Zinc4jet->Fill(jetsHT, weight);
                for (unsigned short j(0); j < nGoodJets; j++){
                    AllJetPt_Zinc4jet->Fill(jets[j].pt, weight);
                    AllJetEta_Zinc4jet->Fill(jets[j].eta, weight);
                    AllJetPhi_Zinc4jet->Fill(jets[j].phi, weight);
                }
                if (nGoodJets == 4){
                    ZNGoodJets_Zexc_NoWeight->Fill(4.);
                    ZMass_Zexc4jet->Fill(Z.M(), weight);
                    ZPt_Zexc4jet->Fill(Z.Pt(), weight);
                    ZRapidity_Zexc4jet->Fill(Z.Rapidity(), weight);
                    ZEta_Zexc4jet->Fill(Z.Eta(), weight);
                    lepPt_Zexc4jet->Fill(lepton1.pt, weight);
                    lepEta_Zexc4jet->Fill(lepton1.eta, weight);
                    if (doZ || doTT){
                        lepPt_Zexc4jet->Fill(lepton2.pt, weight);
                        lepEta_Zexc4jet->Fill(lepton2.eta, weight);
                    }
                    dPhiLeptons_Zexc4jet->Fill(deltaPhi(lep1, lep2), weight);
                    dEtaLeptons_Zexc4jet->Fill(lepton1.eta - lepton2.eta, weight);
                    SpTLeptons_Zexc4jet->Fill(SpTsub(lep1, lep2), weight);
                }
            }
            if (nGoodJets >= 5){
                ZNGoodJets_Zinc->Fill(5., weight);
                ZNGoodJetsFull_Zinc->Fill(5., weight);
                ZNGoodJets_Zinc_NoWeight->Fill(5.);
                ZMass_Zinc5jet->Fill(Z.M(), weight);
                ZPt_Zinc5jet->Fill(Z.Pt(), weight);
                ZRapidity_Zinc5jet->Fill(Z.Rapidity(), weight);
                ZEta_Zinc5jet->Fill(Z.Eta(), weight);
                lepPt_Zinc5jet->Fill(lepton1.pt, weight);
                lepEta_Zinc5jet->Fill(lepton1.eta, weight);
                if (doZ || doTT){
                    lepPt_Zinc5jet->Fill(lepton2.pt, weight);
                    lepEta_Zinc5jet->Fill(lepton2.eta, weight);
                }
                dPhiLeptons_Zinc5jet->Fill(deltaPhi(lep1, lep2), weight);
                dEtaLeptons_Zinc5jet->Fill(lepton1.eta - lepton2.eta, weight);
                dRLeptons_Zinc5jet->Fill(deltaR(lepton1.phi, lepton1.eta, lepton2.phi, lepton2.eta), weight);
                SpTLeptons_Zinc5jet->Fill(SpTsub(lep1, lep2), weight);
                FifthJetPt_Zinc5jet->Fill(jets[4].pt, weight);
                FifthJetPt_1_Zinc5jet->Fill(jets[4].pt, weight);
                FifthJetPt_2_Zinc5jet->Fill(jets[4].pt, weight);
                FirstHighestJetPt_Zinc5jet->Fill(jets[0].pt, weight);
                SecondHighestJetPt_Zinc5jet->Fill(jets[1].pt, weight);
                ThirdHighestJetPt_Zinc5jet->Fill(jets[2].pt, weight);
                FifthJetEta_Zinc5jet->Fill(fabs(jets[4].eta), weight);
                FifthJetEta_2_Zinc5jet->Fill(fabs(jets[4].eta), weight);
                FifthJetEtaFull_Zinc5jet->Fill(jets[4].eta, weight);
                FifthJetPtEta_Zinc5jet->Fill(jets[4].pt, fabs(jets[4].eta), weight);
                FifthJetPhi_Zinc5jet->Fill(jets[4].phi, weight);
                JetsHT_Zinc5jet->Fill(jetsHT, weight);
                JetsHT_1_Zinc5jet->Fill(jetsHT, weight);
                JetsHT_2_Zinc5jet->Fill(jetsHT, weight);
                if (nGoodJets == 5){
                    ZNGoodJets_Zexc_NoWeight->Fill(5.);
                    ZMass_Zexc5jet->Fill(Z.M(), weight);
                    ZPt_Zexc5jet->Fill(Z.Pt(), weight);
                    ZRapidity_Zexc5jet->Fill(Z.Rapidity(), weight);
                    ZEta_Zexc5jet->Fill(Z.Eta(), weight);
                    lepPt_Zexc5jet->Fill(lepton1.pt, weight);
                    lepEta_Zexc5jet->Fill(lepton1.eta, weight);
                    if (doZ || doTT){
                        lepPt_Zexc5jet->Fill(lepton2.pt, weight);
                        lepEta_Zexc5jet->Fill(lepton2.eta, weight);
                    }
                    dPhiLeptons_Zexc5jet->Fill(deltaPhi(lep1, lep2), weight);
                    dEtaLeptons_Zexc5jet->Fill(lepton1.eta - lepton2.eta, weight);
                    SpTLeptons_Zexc5jet->Fill(SpTsub(lep1, lep2), weight);
                }
            }
            if (nGoodJets >= 6){
                ZNGoodJets_Zinc->Fill(6., weight);
                ZNGoodJetsFull_Zinc->Fill(6., weight);
                ZNGoodJets_Zinc_NoWeight->Fill(6.);
                ZMass_Zinc6jet->Fill(Z.M(), weight);
                ZPt_Zinc6jet->Fill(Z.Pt(), weight);
                ZRapidity_Zinc6jet->Fill(Z.Rapidity(), weight);
                ZEta_Zinc6jet->Fill(Z.Eta(), weight);
                SixthJetPt_Zinc6jet->Fill(jets[5].pt, weight);
                SixthJetPt_1_Zinc6jet->Fill(jets[5].pt, weight);
                FirstHighestJetPt_Zinc6jet->Fill(jets[0].pt, weight);
                SecondHighestJetPt_Zinc6jet->Fill(jets[1].pt, weight);
                ThirdHighestJetPt_Zinc6jet->Fill(jets[2].pt, weight);
                SixthJetEta_Zinc6jet->Fill(fabs(jets[5].eta), weight);
                SixthJetEtaFull_Zinc6jet->Fill(jets[5].eta, weight);
                SixthJetPtEta_Zinc6jet->Fill(jets[5].pt, fabs(jets[5].eta), weight);
                SixthJetPhi_Zinc6jet->Fill(jets[5].phi, weight);
                JetsHT_Zinc6jet->Fill(jetsHT, weight);
                JetsHT_1_Zinc6jet->Fill(jetsHT, weight);
                if (nGoodJets == 6){
                    ZNGoodJets_Zexc_NoWeight->Fill(6.);
                    ZMass_Zexc6jet->Fill(Z.M(), weight);
                    ZPt_Zexc6jet->Fill(Z.Pt(), weight);
                    ZRapidity_Zexc6jet->Fill(Z.Rapidity(), weight);
                    ZEta_Zexc6jet->Fill(Z.Eta(), weight);
                }
            }
            if (nGoodJets >= 7){
                ZNGoodJets_Zinc->Fill(7., weight);
                ZNGoodJetsFull_Zinc->Fill(7., weight);
            }
            
            //more bins in Single Lepton dataset than Double --> xsec bigger
            if (nGoodJets >= 8) ZNGoodJets_Zinc->Fill(8., weight);
            if (doW && nGoodJets >= 9)  ZNGoodJets_Zinc->Fill(9., weight);
            if (doW && nGoodJets >= 10) ZNGoodJets_Zinc->Fill(10., weight);
        }
        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        
        //=======================================================================================================//
        //             Unfolding              //
        //====================================//
        if (hasRecoInfo && hasGenInfo && passesGenLeptonCut && passesLeptonCut){
            //-- jet multiplicity exc
            hresponseZNGoodJets_Zexc->Fill(nGoodJets, nGoodGenJets, weight);
            hresponseZNGoodJetsFull_Zexc->Fill(nGoodJets, nGoodGenJets, weight);
            hresponseZNGoodJets_Zinc->Fill(0., 0., weight);
            hresponseZNGoodJetsFull_Zinc->Fill(0., 0., weight);
            //-- First Jet Pt
            if (nGoodGenJets >= 1 && nGoodJets >= 1){
                hresponseZNGoodJets_Zinc->Fill(1., 1., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(1., 1., weight);
                hresponseFirstJetPt_Zinc1jet->Fill(jets[0].pt, genJets[0].pt, weight);
                hresponseFirstJetEta_Zinc1jet->Fill(fabs(jets[0].eta), fabs(genJets[0].eta), weight);
                hresponseJetsHT_Zinc1jet->Fill(jetsHT, genJetsHT, weight);
                
                hresponseFirstJetPt_1_Zinc1jet->Fill(jets[0].pt, genJets[0].pt, weight);
                hresponseFirstJetPt_2_Zinc1jet->Fill(jets[0].pt, genJets[0].pt, weight);
                hresponseFirstJetEta_2_Zinc1jet->Fill(fabs(jets[0].eta), fabs(genJets[0].eta), weight);
                hresponseJetsHT_1_Zinc1jet->Fill(jetsHT, genJetsHT, weight);
                hresponseJetsHT_2_Zinc1jet->Fill(jetsHT, genJetsHT, weight);
            }
            //-- Second Jet Pt
            if (nGoodGenJets >= 2 && nGoodJets >= 2){
                hresponseZNGoodJets_Zinc->Fill(2., 2., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(2., 2., weight);
                hresponseSecondJetPt_Zinc2jet->Fill(jets[1].pt, genJets[1].pt, weight);
                hresponseSecondJetEta_Zinc2jet->Fill(fabs(jets[1].eta), fabs(genJets[1].eta), weight);
                hresponseJetsHT_Zinc2jet->Fill(jetsHT, genJetsHT, weight);
                
                hresponseSecondJetPt_1_Zinc2jet->Fill(jets[1].pt, genJets[1].pt, weight);
                hresponseSecondJetPt_2_Zinc2jet->Fill(jets[1].pt, genJets[1].pt, weight);
                hresponseSecondJetEta_2_Zinc2jet->Fill(fabs(jets[1].eta), fabs(genJets[1].eta), weight);
                hresponseJetsHT_1_Zinc2jet->Fill(jetsHT, genJetsHT, weight);
                hresponseJetsHT_2_Zinc2jet->Fill(jetsHT, genJetsHT, weight);
                
                hresponsedRapidityJets_Zinc2jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), weight);
                hresponsedRapidityJets_2_Zinc2jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), weight);
                
                hresponsedRapidityJetsFB_Zinc2jet->Fill(ForwardJetRapidity - BackwardJetRapidity, genForwardJetRapidity - genBackwardJetRapidity, weight);
                hresponsedRapidityJetsFB_2_Zinc2jet->Fill(ForwardJetRapidity - BackwardJetRapidity, genForwardJetRapidity - genBackwardJetRapidity, weight);
                
                hresponsedPhiJets_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), deltaPhi(genLeadJ, genSecondJ), weight);
                hresponsedPhiJets_2_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), deltaPhi(genLeadJ, genSecondJ), weight);
                
                hresponsedPhiJetsFB_Zinc2jet->Fill(deltaPhi(vJetYOrdered[0], vJetYOrdered[vJetYOrdered.size() - 1] ), deltaPhi(genvJetYOrdered[0], genvJetYOrdered[genvJetYOrdered.size() - 1] ), weight);
                hresponsedPhiJetsFB_2_Zinc2jet->Fill(deltaPhi(vJetYOrdered[0], vJetYOrdered[vJetYOrdered.size() - 1] ), deltaPhi(genvJetYOrdered[0], genvJetYOrdered[genvJetYOrdered.size() - 1] ), weight);
                
                hresponsedRJets_Zinc2jet->Fill(deltaRYPhi(newLeadJ, newSecondJ), deltaRYPhi(genNewLeadJ, genNewSecondJ), weight);
                hresponsedRJets_2_Zinc2jet->Fill(deltaRYPhi(newLeadJ, newSecondJ), deltaRYPhi(genNewLeadJ, genNewSecondJ), weight);
                
                hresponsediJetMass_Zinc2jet->Fill(jet1Plus2.M(), genJet1Plus2.M(), weight);
                hresponsediJetMass_2_Zinc2jet->Fill(jet1Plus2.M(), genJet1Plus2.M(), weight);
                
                hresponsediJetPt_Zinc2jet->Fill(jet1Plus2.Pt(), genJet1Plus2.Pt(), weight);
                hresponsediJetPt_2_Zinc2jet->Fill(jet1Plus2.Pt(), genJet1Plus2.Pt(), weight);
            }
            //-- Third Jet Pt
            if (nGoodGenJets >= 3 && nGoodJets >= 3){
                hresponseZNGoodJets_Zinc->Fill(3., 3., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(3., 3., weight);
                hresponseThirdJetPt_Zinc3jet->Fill(jets[2].pt, genJets[2].pt, weight);
                hresponseThirdJetEta_Zinc3jet->Fill(fabs(jets[2].eta), fabs(genJets[2].eta), weight);
                hresponseJetsHT_Zinc3jet->Fill(jetsHT, genJetsHT, weight);
                
                hresponseThirdJetPt_1_Zinc3jet->Fill(jets[2].pt, genJets[2].pt, weight);
                hresponseThirdJetPt_2_Zinc3jet->Fill(jets[2].pt, genJets[2].pt, weight);
                hresponseThirdJetEta_2_Zinc3jet->Fill(fabs(jets[2].eta), fabs(genJets[2].eta), weight);
                hresponseJetsHT_1_Zinc3jet->Fill(jetsHT, genJetsHT, weight);
                hresponseJetsHT_2_Zinc3jet->Fill(jetsHT, genJetsHT, weight);
                
                hresponsedRapidityJets_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), weight);
                hresponsedRapidityJets_2_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), weight);
                
                hresponsedRapidityJetsFB_Zinc3jet->Fill(ForwardJetRapidity - BackwardJetRapidity, genForwardJetRapidity - genBackwardJetRapidity, weight);
                hresponsedRapidityJetsFB_2_Zinc3jet->Fill(ForwardJetRapidity - BackwardJetRapidity, genForwardJetRapidity - genBackwardJetRapidity, weight);
                
                hresponsedRapidityJets_First_Third_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newThirdJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewThirdJ.Rapidity()), weight);
                hresponsedRapidityJets_2_First_Third_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newThirdJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewThirdJ.Rapidity()), weight);
                
                hresponsedRapidityJets_Second_Third_Zinc3jet->Fill(fabs(newSecondJ.Rapidity() - newThirdJ.Rapidity()), fabs(genNewSecondJ.Rapidity() - genNewThirdJ.Rapidity()), weight);
                hresponsedRapidityJets_2_Second_Third_Zinc3jet->Fill(fabs(newSecondJ.Rapidity() - newThirdJ.Rapidity()), fabs(genNewSecondJ.Rapidity() - genNewThirdJ.Rapidity()), weight);
                
                hresponsediJetMass_Zinc3jet->Fill(jet1Plus2.M(),genJet1Plus2.M(), weight);
                hresponsediJetMass_2_Zinc3jet->Fill(jet1Plus2.M(),genJet1Plus2.M(), weight);
                
                hresponsediJetPt_Zinc3jet->Fill(jet1Plus2.Pt(), genJet1Plus2.Pt(), weight);
                hresponsediJetPt_2_Zinc3jet->Fill(jet1Plus2.Pt(), genJet1Plus2.Pt(), weight);
            }
            //-- Fourth Jet Pt
            if (nGoodGenJets >= 4 && nGoodJets >= 4){
                hresponseZNGoodJets_Zinc->Fill(4., 4., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(4., 4., weight);
                hresponseFourthJetPt_Zinc4jet->Fill(jets[3].pt, genJets[3].pt, weight);
                hresponseFourthJetEta_Zinc4jet->Fill(fabs(jets[3].eta), fabs(genJets[3].eta), weight);
                hresponseJetsHT_Zinc4jet->Fill(jetsHT, genJetsHT, weight);
                
                hresponseFourthJetPt_1_Zinc4jet->Fill(jets[3].pt, genJets[3].pt, weight);
                hresponseFourthJetPt_2_Zinc4jet->Fill(jets[3].pt, genJets[3].pt, weight);
                hresponseFourthJetEta_2_Zinc4jet->Fill(fabs(jets[3].eta), fabs(genJets[3].eta), weight);
                hresponseJetsHT_1_Zinc4jet->Fill(jetsHT, genJetsHT, weight);
                hresponseJetsHT_2_Zinc4jet->Fill(jetsHT, genJetsHT, weight);
                
                
                hresponsedRapidityJets_Zinc4jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), weight);
                hresponsedRapidityJets_2_Zinc4jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), weight);
                
                hresponsedRapidityJetsFB_Zinc4jet->Fill(ForwardJetRapidity - BackwardJetRapidity, genForwardJetRapidity - genBackwardJetRapidity, weight);
                hresponsedRapidityJetsFB_2_Zinc4jet->Fill(ForwardJetRapidity - BackwardJetRapidity, genForwardJetRapidity - genBackwardJetRapidity, weight);
                
                hresponsediJetMass_Zinc4jet->Fill(jet1Plus2.M(), genJet1Plus2.M(), weight);
                hresponsediJetMass_2_Zinc4jet->Fill(jet1Plus2.M(), genJet1Plus2.M(), weight);
                
                hresponsediJetPt_Zinc4jet->Fill(jet1Plus2.Pt(), genJet1Plus2.Pt(), weight);
                hresponsediJetPt_2_Zinc4jet->Fill(jet1Plus2.Pt(), genJet1Plus2.Pt(), weight);
            }
            //-- Fifth Jet Pt
            if (nGoodGenJets >= 5 && nGoodJets >= 5){
                hresponseZNGoodJets_Zinc->Fill(5., 5., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(5., 5., weight);
                hresponseFifthJetPt_Zinc5jet->Fill(jets[4].pt, genJets[4].pt, weight);
                hresponseFifthJetEta_Zinc5jet->Fill(fabs(jets[4].eta), fabs(genJets[4].eta), weight);
                hresponseJetsHT_Zinc5jet->Fill(jetsHT, genJetsHT, weight);
                
                hresponseFifthJetPt_1_Zinc5jet->Fill(jets[4].pt, genJets[4].pt, weight);
                hresponseFifthJetPt_2_Zinc5jet->Fill(jets[4].pt, genJets[4].pt, weight);
                hresponseFifthJetEta_2_Zinc5jet->Fill(fabs(jets[4].eta), fabs(genJets[4].eta), weight);
                hresponseJetsHT_1_Zinc5jet->Fill(jetsHT, genJetsHT, weight);
                hresponseJetsHT_2_Zinc5jet->Fill(jetsHT, genJetsHT, weight);
            }
            //-- Sixth Jet Pt
            if (nGoodGenJets >= 6 && nGoodJets >= 6){
                hresponseZNGoodJets_Zinc->Fill(6., 6., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(6., 6., weight);
            }
            //-- inc 7 jets case
            if (nGoodGenJets >= 7 && nGoodJets >= 7){
                hresponseZNGoodJets_Zinc->Fill(7., 7., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(7., 7., weight);
            }
            //-- inc 8 jets case
            if (nGoodGenJets >= 8 && nGoodJets >= 8){
                hresponseZNGoodJets_Zinc->Fill(8., 8., weight);
            }
            //-- inc 9 jets case
            if (nGoodGenJets >= 9 && nGoodJets >= 9){
                hresponseZNGoodJets_Zinc->Fill(9., 9., weight);
            }
            //-- inc 9 jets case
            if (nGoodGenJets >= 10 && nGoodJets >= 10){
                hresponseZNGoodJets_Zinc->Fill(10., 10., weight);
            }
            
        }
        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        
        // MeanNJ use old method of filling response.
        if (hasRecoInfo && hasGenInfo){
            
            //--- Set weight for systematic for reweighted response
            double WforMeanNHT1       = weight;
            double genWforMeanNHT1    = genWeight;
            double WforMeanNHT2       = weight;
            double genWforMeanNHT2    = genWeight;
            double WforMeanNdRap12    = weight;
            double genWforMeanNdRap12 = genWeight;
            double WforMeanNdRapFB    = weight;
            double genWforMeanNdRapFB = genWeight;
            if (doRespSyst){
                WforMeanNHT1 = ReweightForResp(vecFwHT1, binEdgeHT1, weight, jetsHT);
                genWforMeanNHT1 = ReweightForResp(vecFwHT1, binEdgeHT1, genWeight, jetsHT);
                WforMeanNHT2 = ReweightForResp(vecFwHT2, binEdgeHT2, weight, jetsHT);
                genWforMeanNHT2 = ReweightForResp(vecFwHT2, binEdgeHT2, genWeight, jetsHT);
                WforMeanNdRap12 = ReweightForResp(vecFwRap12, binEdgeRap12, weight, fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()));
                genWforMeanNdRap12 = ReweightForResp(vecFwRap12, binEdgeRap12, genWeight, fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()));
                WforMeanNdRapFB = ReweightForResp(vecFwRapFB, binEdgeRapFB, weight, fabs(ForwardJetRapidity - BackwardJetRapidity));
                genWforMeanNdRapFB = ReweightForResp(vecFwRapFB, binEdgeRapFB, genWeight, fabs(ForwardJetRapidity - BackwardJetRapidity));
            }
            //------------------------------

            //-- one jet inclusive
            if (nGoodGenJets >= 1 && passesGenLeptonCut){
                if (nGoodJets >= 1 && passesLeptonCut){
                    responseMeanNJetsHT_Zinc1jet->Fill(jetsHT, nGoodJets, genJetsHT, nGoodGenJets, WforMeanNHT1);
                }
                else{
                    responseMeanNJetsHT_Zinc1jet->Miss(genJetsHT, nGoodGenJets, genWforMeanNHT1);
                }
            }
            if (nGoodJets >= 1 && passesLeptonCut){
                if (!(nGoodGenJets >= 1 && passesGenLeptonCut)){
                    responseMeanNJetsHT_Zinc1jet->Fake(jetsHT, nGoodJets, WforMeanNHT1);
                }
            }
            
            //-- two jets inclusive
            if (nGoodGenJets >= 2 && passesGenLeptonCut){
                if (nGoodJets >= 2 && passesLeptonCut){
                    responseMeanNJetsHT_Zinc2jet->Fill(jetsHT, nGoodJets, genJetsHT, nGoodGenJets, WforMeanNHT2);
                    responseMeanNJetsdRapidity_Zinc2jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), nGoodJets, fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), nGoodGenJets, WforMeanNdRap12);
                    responseMeanNJetsdRapidityFB_Zinc2jet->Fill(ForwardJetRapidity - BackwardJetRapidity, nGoodJets, genForwardJetRapidity - genBackwardJetRapidity, nGoodGenJets, WforMeanNdRapFB);
                }
                else {
                    responseMeanNJetsHT_Zinc2jet->Miss(genJetsHT, nGoodGenJets, genWforMeanNHT2);
                    responseMeanNJetsdRapidity_Zinc2jet->Miss(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), nGoodGenJets, genWforMeanNdRap12);
                    responseMeanNJetsdRapidityFB_Zinc2jet->Miss(genForwardJetRapidity - genBackwardJetRapidity, nGoodGenJets, genWforMeanNdRapFB);
                }
            }
            if (nGoodJets >= 2 && passesLeptonCut){
                if (!(nGoodGenJets >= 2 && passesGenLeptonCut)){
                    responseMeanNJetsHT_Zinc2jet->Fake(jetsHT, nGoodJets, WforMeanNHT2);
                    responseMeanNJetsdRapidity_Zinc2jet->Fake(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), nGoodJets, WforMeanNdRap12);
                    responseMeanNJetsdRapidityFB_Zinc2jet->Fake(ForwardJetRapidity - BackwardJetRapidity, nGoodJets, WforMeanNdRapFB);
                }
            }
        }
        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        
        //=======================================================================================================//

    } //End of loop over all the events
    
    //==========================================================================================================//

    if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    //==========================================================================================================//
    //         Writing file             //
    //==================================//

    outputFile->cd();

    //--- Save all the histograms ---
    unsigned short numbOfHistograms = listOfHistograms.size();
    for (unsigned short i(0); i < numbOfHistograms; i++){
        string hName = listOfHistograms[i]->GetName();
        if ( (!hasGenInfo && hName.find("gen") != string::npos ) || (!hasRecoInfo && hName.find("gen") == string::npos )) continue; 
        if (fileName.find("HepMC") != string::npos){
            //if (sumSherpaW > 0) listOfHistograms[i]->Scale(1/sumSherpaW);
            if (doTTreweighting) listOfHistograms[i]->Scale(weightSumNoTopRew / weightSum );
        }
        if (fileName.find("Sherpa2") != string::npos){
            if (sumSherpaW > 0) listOfHistograms[i]->Scale(1/sumSherpaW);
        }
        listOfHistograms[i]->Write();
    }
    //--- Save all the RooUnfoldResponses ---
    if ( hasGenInfo && hasRecoInfo ){
        unsigned short numbOfResponses = listOfResponses.size();
        for (unsigned short i(0); i < numbOfResponses; i++){
            string currentName = listOfResponses[i]->GetName();
            currentName = currentName.substr(0, currentName.find("_gen"));
            string savingName = "response" + currentName;
            outputFile->WriteTObject(listOfResponses[i], savingName.c_str());
        }
    }
    
    // let's delete all histograms, just to be safe
    for (unsigned short i(0); i < numbOfHistograms; i++){
        delete listOfHistograms[i];
    }
    
    outputFile->Write();
    outputFile->Close();
    
    //==========================================================================================================//


    cout << "Number of events                               : " << nEvents << endl;
    cout << "Total GEN weight of all events                 : " << TotalGenWeight << endl;
    cout << "Number with two good leptons no charge no mass : " << nEventsWithTwoGoodLeptonsNoChargeNoMass << endl;
    cout << "Number with two good leptons no mass           : " << nEventsWithTwoGoodLeptonsNoMass << endl;
    cout << "Total GEN pass: RECO weight of all events      : " << TotalGenWeightPassGENPU << endl;
    cout << "Total GEN pass: GEN weight of all events       : " << TotalGenWeightPassGEN << endl;
    cout << "Total RECO pass: RECO weight of all events     : " << TotalRecoWeightPassRECO << endl;
    cout << "Total RECO pass: GEN weight of all events      : " << TotalGenWeightPassRECO << endl;
    cout << "Number with two good leptons                   : " << nEventsWithTwoGoodLeptons << endl;
    cout << "How many times do we visit unfolding 0 jets    : " << nEventsUNFOLDIncl0Jets << endl;
    cout << "Number Inclusif 0 jets                         : " << nEventsIncl0Jets << endl;
    cout << "Number Exclusif 0 jets                         : " << nEventsExcl0Jets << endl;
    cout << "Number Exclusif 1 jets                         : " << nEventsExcl1Jets << endl;
    cout << "Number Exclusif 2 jets                         : " << nEventsExcl2Jets << endl;
    cout << "Number Exclusif 3 jets                         : " << nEventsExcl3Jets << endl;
    cout << "Number Inclusive 1 B-jet                       : " << nEventsIncBJets << endl;
    cout << "Number GEN Inclusif 0 jets                     : " << GENnEventsIncl0Jets << endl;
    cout << "Number GEN Inclusif 1 jets                     : " << GENnEventsIncl1Jets << endl;
    cout << "Number GEN Inclusif 2 jets                     : " << GENnEventsIncl2Jets << endl;
    cout << "Number GEN Inclusif 3 jets                     : " << GENnEventsIncl3Jets << endl;
    cout << "Sherpa weight                                  : " << sumSherpaW << endl;
    if (doTTreweighting)       cout << "We run to TTbar with reweighting :   " << weightSum << "  and the original weight is :" << weightSumNoTopRew << endl;
    cout << " Trigger summary"<< endl;
    for (unsigned short k(0); k < 4; k++) {
        if (countTrigSum[k] > 0) cout << sumTrig[k] << "    " << countTrigSum[k]  << "    " << sumTrig[k]/countTrigSum[k] << endl;
    }
}











ZJetsAndDPS::ZJetsAndDPS(string fileName_, float lumiScale_, float puScale_, bool useTriggerCorrection_, bool useEfficiencyCorrection_, 
        int systematics_, int direction_, float xsecfactor_, int jetPtCutMin_, int jetPtCutMax_, int ZPtCutMin_, int ZEtaCutMin_, int ZEtaCutMax_, int METcut_, bool nEvents_10000_, int jetEtaCutMin_, int jetEtaCutMax_): 
    HistoSet(fileName_.substr(0, fileName_.find("_"))), nEvents_10000(nEvents_10000_), outputDirectory("HistoFiles/"),
    fileName(fileName_), lumiScale(lumiScale_), puScale(puScale_), useTriggerCorrection(useTriggerCorrection_), useEfficiencyCorrection(useEfficiencyCorrection_), 
    systematics(systematics_), direction(direction_), xsecfactor(xsecfactor_), jetPtCutMin(jetPtCutMin_), jetPtCutMax(jetPtCutMax_), jetEtaCutMin(jetEtaCutMin_), jetEtaCutMax(jetEtaCutMax_), ZPtCutMin(ZPtCutMin_), ZEtaCutMin(ZEtaCutMin_), ZEtaCutMax(ZEtaCutMax_), METcut(METcut_)
{

    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.

    TChain *chain = new TChain("", "");
    isData = false;
    string fullFileName =  "../Data_Z_5311_New/" + fileName;


    if (fileName.find("DMu_") == 0) leptonFlavor = "Muons";
    else if (fileName.find("DE_") == 0)  leptonFlavor = "Electrons"; 
    else if (fileName.find("SMu_") == 0) leptonFlavor = "SingleMuon";
    else if (fileName.find("SE_") == 0)  leptonFlavor = "SingleElectron";
    else if (fileName.find("SMuE_") == 0){
        leptonFlavor = "TTMuE";
        fullFileName =  "../DataTTbarEMu/" + fileName;
    }
    if (fileName.find("Data") != string::npos ) isData = true;
    if ( fileName.find("SMu_") == 0 || fileName.find("SE_") == 0 ) fullFileName =  "../DataW/" + fileName;
    if (fileName.find("Sherpa2") != string::npos) fullFileName =  "../DataSherpa2/" + fileName;
    if (fileName.find("List") == string::npos){
        if (fileName.find("Sherpa2") != string::npos){
            fullFileName += ".root";
            string treePath = fullFileName + "/tree";
            cout << "Loading file: " << fullFileName << endl;
            chain->Add(treePath.c_str());
        }
        else{
            fullFileName += ".root";
            string treePath = fullFileName + "/tupel/EventTree";
            cout << "Loading file: " << fullFileName << endl;
            chain->Add(treePath.c_str());
        }
    }
    else {
        fullFileName += ".txt";
        ifstream infile(fullFileName.c_str());
        string line; 
        int countFiles(0);
        while (getline(infile, line)){
            countFiles++;
            string treePath =  line + "/tree/tree";
            chain->Add(treePath.c_str());       
        }
    }
    fChain = chain;
}

ZJetsAndDPS::~ZJetsAndDPS(){
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

string ZJetsAndDPS::CreateOutputFileName(bool useRoch, bool doFlat, int doPUStudy, bool doVarWidth, int doBJets, int doQCD, bool doSSign, bool doInvMassCut, string pdfSet, int pdfMember)
{
    ostringstream result;
    result << outputDirectory << fileName;
    result << "_EffiCorr_" << useEfficiencyCorrection;
    result << "_TrigCorr_" << useTriggerCorrection;
    result << "_Syst_" << systematics;
    if (direction == 1) result << "_Up";
    else if (direction == -1) result << "_Down";
    result << "_JetPtMin_" << jetPtCutMin;
    if (jetPtCutMax > jetPtCutMin) result << "_JetPtMax_" << jetPtCutMax;
    if (ZPtCutMin > 0) result << "_ZPtMin" << abs(ZPtCutMin);
    if (ZEtaCutMin > -999999 && ZEtaCutMin <  0) result << "_ZEtaMin_m" << abs(ZEtaCutMin);
    if (ZEtaCutMin > -999999 && ZEtaCutMin >= 0) result << "_ZEtaMin_"  << abs(ZEtaCutMin);
    if (ZEtaCutMax <  999999 && ZEtaCutMax >= 0) result << "_ZEtaMax_"  << abs(ZEtaCutMax);
    if (ZEtaCutMax <  999999 && ZEtaCutMax <  0) result << "_ZEtaMax_m" << abs(ZEtaCutMax);

    if (useRoch) result << "_rochester";
    if (!isData && doFlat) result << "_Flat";
    if (doPUStudy >= 0) result << "_Beta" << doPUStudy;
    if (doVarWidth) result << "_VarWidth";
    if (doInvMassCut) result << "_InvMass";
    if (doSSign) result << "_SS";
    if (doBJets > 0) result << "_BJets";
    if (doBJets < 0) result << "_BVeto";
    if (doQCD>0) result << "_QCD" << doQCD;
    if (METcut > 0) result << "_MET" << METcut;
    if (pdfSet != "") result << "_PDF_" << pdfSet << "_" << pdfMember;

    //--- Add your test names here ---
    //result << "_NoPUCut";
    //result << "_LooseID";
    //result << "_SRANJE";

    result << ".root";
    return result.str();
}

Int_t ZJetsAndDPS::GetEntry(Long64_t entry){
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t ZJetsAndDPS::LoadTree(Long64_t entry){
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();Notify();  }
    return centry;
}

void ZJetsAndDPS::Init(bool hasRecoInfo, bool hasGenInfo, bool hasPartonInfo){
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    mcSherpaWeights_ = 0 ;
    pdfInfo_ = 0;
    genLepPt_ = 0;
    genLepEta_ = 0;
    genLepPhi_ = 0;
    genLepE_ = 0;
   // genLepQ_ = 0;
    genLepId_ = 0;
    genLepSt_ = 0;
    genPhoPt_ = 0;
    genPhoEta_ = 0;
    genPhoPhi_ = 0;
    genJetPt_ = 0;
    genJetEta_ = 0;
    genJetPhi_ = 0;
    genJetE_ = 0;
    dpsParton_Pt = 0;
    dpsParton_Eta = 0;
    dpsParton_Phi = 0;
    dpsParton_E = 0;
    genMatchDPSpar = 0;
    dpsParton_dR = 0;

    //gsfElecPt_ = 0;
    //gsfElecEta_ = 0;
    //gsfElecPhi_ = 0;
    //gsfElecEnergy_ = 0;

    ElPt = 0;
    ElEta = 0;
    ElPhi = 0;
    ElE = 0;
    ElCh = 0;
    ElId = 0;
    //ElTrig = 0;
    ElPfIsoRho = 0;
    ElEtaSc = 0;

    MuPt = 0;
    MuEta = 0;
    MuPhi = 0;
    MuVtxZ = 0;
    MuE = 0;
    MuCh = 0;
    MuDxy = 0;
    MuIdTight = 0;
    //MuTrig = 0;
    MuPfIso = 0;

    JetAk04Pt = 0;   
    JetAk04Eta = 0;   
    JetAk04Phi = 0;   
    JetAk04E = 0;   
    JetAk04Id = 0;   
    JetAk04PuMva = 0;   
    JetAk04BTagCsv = 0;   
    JetAk04PartFlav = 0;   
    JetAk04JetBeta = 0;   
    JetAk04JetBetaStar = 0;   
  
    METPt = 0 ;
    METPx = 0 ;
    METPy = 0 ;
    METPz = 0 ;
    METE = 0 ;
    //koMETPhi = 0 ;
    //METsig = 0 ;
 if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    // Set branch addresses and branch pointers
    fCurrent = -1;
    fChain->SetMakeClass(1);
    if (fileName.find("Data") == string::npos) fChain->SetBranchAddress("PU_npT", &PU_npT, &b_PU_npT);
    if (fileName.find("UNFOLDING") != string::npos) fChain->SetBranchAddress("nup_", &nup_, &b_nup_);
    if (hasRecoInfo){
        fChain->SetBranchAddress("EvtVtxCnt", &EvtVtxCnt, &b_EvtVtxCnt);
        fChain->SetBranchAddress("EvtRunNum", &EvtRunNum, &b_EvtRunNum); // not used
        fChain->SetBranchAddress("EvtNum", &EvtNum, &b_EvtNum); // not used

        fChain->SetBranchAddress("JetAk04E", &JetAk04E, &b_JetAk04E);
        fChain->SetBranchAddress("JetAk04Pt", &JetAk04Pt, &b_JetAk04Pt);
        fChain->SetBranchAddress("JetAk04Eta", &JetAk04Eta, &b_JetAk04Eta);
        fChain->SetBranchAddress("JetAk04Phi", &JetAk04Phi, &b_JetAk04Phi);
        fChain->SetBranchAddress("JetAk04Id", &JetAk04Id, &b_JetAk04Id);
        fChain->SetBranchAddress("JetAk04PuMva", &JetAk04PuMva, &b_JetAk04PuMva);
        fChain->SetBranchAddress("JetAk04BTagCsv", &JetAk04BTagCsv, &b_JetAk04BTagCsv);
        fChain->SetBranchAddress("JetAk04PartFlav", &JetAk04PartFlav, &b_JetAk04PartFlav);
        fChain->SetBranchAddress("JetAk04JetBeta", &JetAk04JetBeta, &b_JetAk04JetBeta);
        fChain->SetBranchAddress("JetAk04JetBetaStar", &JetAk04JetBetaStar, &b_JetAk04JetBetaStar);
        fChain->SetBranchAddress("METPt", &METPt, &b_METPt);
        fChain->SetBranchAddress("METPx", &METPx, &b_METPx);
        fChain->SetBranchAddress("METPy", &METPy, &b_METPy);
        fChain->SetBranchAddress("METPz", &METPz, &b_METPz);
        fChain->SetBranchAddress("METE", &METE, &b_METE);
        //KOfChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
        //fChain->SetBranchAddress("METsig", &METsig, &b_METsig); // not used
  
        if (leptonFlavor != "Muons"){
            //fChain->SetBranchAddress("gsfElecPt_", &gsfElecPt_, &b_gsfElecPt_); // not used
            //fChain->SetBranchAddress("gsfElecEta_", &gsfElecEta_, &b_gsfElecEta_); // not used
            //fChain->SetBranchAddress("gsfElecPhi_", &gsfElecPhi_, &b_gsfElecPhi_); // not used
            //fChain->SetBranchAddress("gsfElecEnergy_", &gsfElecEnergy_, &b_gsfElecEnergy_); // not used
            fChain->SetBranchAddress("ElPt", &ElPt, &b_ElPt);
            fChain->SetBranchAddress("ElEta", &ElEta, &b_ElEta);
            fChain->SetBranchAddress("ElPhi", &ElPhi, &b_ElPhi);
            fChain->SetBranchAddress("ElE", &ElE, &b_ElE);
            fChain->SetBranchAddress("ElCh", &ElCh, &b_ElCh);
            fChain->SetBranchAddress("ElId", &ElId, &b_ElId);
            //fChain->SetBranchAddress("ElTrig", &ElTrig, &b_ElTrig);
            //fChain->SetBranchAddress("patElecDetIso_", &patElecDetIso_, &b_patElecDetIso_); // not used
            fChain->SetBranchAddress("ElPfIsoRho", &ElPfIsoRho, &b_ElPfIsoRho); 
            fChain->SetBranchAddress("ElEtaSc", &ElEtaSc, &b_ElEtaSc);
            //fChain->SetBranchAddress("patElecIsPF_", &patElecIsPF_, &b_patElecIsPF_); // not used

        }
        if (leptonFlavor != "Electrons"){
            fChain->SetBranchAddress("MuPt", &MuPt, &b_MuPt);
            fChain->SetBranchAddress("MuEta", &MuEta, &b_MuEta);
            fChain->SetBranchAddress("MuPhi", &MuPhi, &b_MuPhi);
            //fChain->SetBranchAddress("patMuonVtxZ_", &patMuonVtxZ_, &b_patMuonVtxZ_); // not used
            fChain->SetBranchAddress("MuE", &MuE, &b_MuE);
            fChain->SetBranchAddress("MuCh", &MuCh, &b_MuCh);
            fChain->SetBranchAddress("MuDxy", &MuDxy, &b_MuDxy);
            fChain->SetBranchAddress("MuIdTight", &MuIdTight, &b_MuIdTight);
          //  fChain->SetBranchAddress("MuTrig", &MuTrig, &b_MuTrig);
            fChain->SetBranchAddress("MuPfIso", &MuPfIso, &b_MuPfIso);
        }
    }
    if (hasGenInfo){
        fChain->SetBranchAddress("genLepPt_", &genLepPt_, &b_genLepPt_);
        fChain->SetBranchAddress("genLepEta_", &genLepEta_, &b_genLepEta_);
        fChain->SetBranchAddress("genLepPhi_", &genLepPhi_, &b_genLepPhi_);
        fChain->SetBranchAddress("genLepE_", &genLepE_, &b_genLepE_);
        //fChain->SetBranchAddress("genLepQ_", &genLepQ_, &b_genLepQ_);
        fChain->SetBranchAddress("genJetPt_", &genJetPt_, &b_genJetPt_);
        fChain->SetBranchAddress("genJetEta_", &genJetEta_, &b_genJetEta_);
        fChain->SetBranchAddress("genJetPhi_", &genJetPhi_, &b_genJetPhi_);
        fChain->SetBranchAddress("genJetE_", &genJetE_, &b_genJetE_);

        if (fileName.find("Sherpa") != string::npos || 
                (fileName.find("UNFOLDING") != string::npos && hasGenInfo ) || 
                fileName.find("Powheg") != string::npos || 
                fileName.find("P8") != string::npos || 
                fileName.find("TopReweighting") != string::npos ||
                fileName.find("Z2") != string::npos){

            fChain->SetBranchAddress("pdfInfo_", &pdfInfo_, &b_pdfInfo_);
            fChain->SetBranchAddress("genLepId_", &genLepId_, &b_genLepId_);
            fChain->SetBranchAddress("genLepSt_", &genLepSt_, &b_genLepSt_);
            fChain->SetBranchAddress("genPhoPt_", &genPhoPt_, &b_genPhoPt_);
            fChain->SetBranchAddress("genPhoEta_", &genPhoEta_, &b_genPhoEta_);
            fChain->SetBranchAddress("genPhoPhi_", &genPhoPhi_, &b_genPhoPhi_);

            if (fileName.find("MiNLO") != string::npos || 
                    fileName.find("mcEveWeight") != string::npos || 
                    fileName.find("HepMC") != string::npos){
                fChain->SetBranchAddress("mcEveWeight_", &mcEveWeight_, &b_mcEveWeight_);
            }

            if (fileName.find("HepMC") != string::npos){
                fChain->SetBranchAddress("mcSherpaSumWeight3_", &mcSherpaSumWeight3_, &b_mcSherpaSumWeight3_);
                fChain->SetBranchAddress("mcSherpaWeights_", &mcSherpaWeights_, &b_mcSherpaWeights_);
            }
            
            if (fileName.find("Sherpa2") != string::npos){
                fChain->SetBranchAddress("mcSherpaWeights_", &mcSherpaWeights_, &b_mcSherpaWeights_);
            }

        }
    }
    if (hasPartonInfo){
        fChain->SetBranchAddress("dpsParton_Pt", &dpsParton_Pt, &b_dpsParton_Pt);
        fChain->SetBranchAddress("dpsParton_Eta", &dpsParton_Eta, &b_dpsParton_Eta);
        fChain->SetBranchAddress("dpsParton_Phi", &dpsParton_Phi, &b_dpsParton_Phi);
        fChain->SetBranchAddress("dpsParton_E", &dpsParton_E, &b_dpsParton_E);
        fChain->SetBranchAddress("genMatchDPSpar", &genMatchDPSpar, &b_genMatchDPSpar);
        fChain->SetBranchAddress("dpsParton_dR", &dpsParton_dR, &b_dpsParton_dR);
    }
    Notify();
    cout << "Branches are properly initialized." << endl;
}

Bool_t ZJetsAndDPS::Notify(){
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void ZJetsAndDPS::Show(Long64_t entry){
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}

Int_t ZJetsAndDPS::Cut(Long64_t entry){
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    printf("entry %lld", entry);
    return 1;
}
