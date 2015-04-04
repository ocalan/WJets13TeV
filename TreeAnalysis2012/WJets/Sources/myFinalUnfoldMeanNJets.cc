#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <vector>
#include <string.h>

#include <TROOT.h>
#include <TFile.h>
#include <TDirectory.h>

#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TList.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>

#include <RooUnfoldResponse.h>
#include <RooUnfold.h>
#include <RooUnfoldBayes.h>
#include <RooUnfoldSvd.h>
#include <RooUnfoldBinByBin.h>
#include <RooUnfoldInvert.h>
#include <TSVDUnfold.h>

#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TPaveStats.h>
#include <TText.h>
#include <THStack.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TDatime.h>

#include "getFilesAndHistograms.h"
#include "unfoldingFunctions.h"
#include "writeFastPlotsTex.h"

#include "myUnfoldingSyst.h"
#include "myFinalUnfoldMeanNJets.h"

using namespace std;

//--  Setting global variables --------------------------------------------------------------
#include <TRandom3.h>
#include <TMath.h>
#include "fileNames.h"
#include "variablesOfInterestVarWidth.h"
TRandom3* RandGen = new TRandom3();
const int NFILESWJETS_M(15);
const int FilesWJets_M[NFILESWJETS_M] = {0, 2, 3, 4, 23, 5, 6, 7, 8, 9, 10, 11, 12, 13, 24};
//-------------------------------------------------------------------------------------------


void myFinalUnfoldMeanNJets(string var, string leptonFlavor, string unfAlg)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    //--- load the proper VAROFINTERST structure for W+jets study---
    TList *varList = new TList();
    int nVarStruct(0);
    variableStruct *varStruct = NULL;
    
    nVarStruct = NVAROFINTERESTWJETS;
    varStruct = new variableStruct[nVarStruct];
    varStruct = VAROFINTERESTWJETS;
    
    makeLISTOFVAROFINTERESTWJETS();
    varList = LISTOFVAROFINTERESTWJETS;
    //---------------------------------------------
    
    //--- check whether muon or electron is used ---
    bool isMuon;
    if (leptonFlavor.find("Mu") != string::npos) isMuon = 1;
    else isMuon = 0;
    cout << " Are you doing Muon? : " << isMuon << endl;
    //---------------------------------------------
    
    TObjString *ObjNameVar = new TObjString(var.c_str());
    unsigned short j(varList->IndexOf(ObjNameVar));
    cout << " Processing variable: " << varStruct[j].name << " with BayesKterm : " << varStruct[j].MuBayeskterm << endl;
    
    
    if (leptonFlavor == "SMu" || leptonFlavor == "SE"){
        if (isMuon) FuncUnfold(leptonFlavor, isMuon, unfAlg, varStruct[j].name, varStruct[j].MuBayeskterm, varStruct[j].MuSVDkterm);
        else FuncUnfold(leptonFlavor, isMuon, unfAlg, varStruct[j].name, varStruct[j].EBayeskterm, varStruct[j].ESVDkterm);
    }

    /*
    else if (leptonFlavor == "DMu" || leptonFlavor == "DE"){
        if (isMuon) FuncUnfold(leptonFlavor, isMuon, unfAlg, VAROFINTERESTZJETS[i].name, VAROFINTERESTZJETS[i].MuBayeskterm, VAROFINTERESTZJETS[i].MuSVDkterm);
        else FuncUnfold(leptonFlavor, isMuon, unfAlg, VAROFINTERESTZJETS[i].name, VAROFINTERESTZJETS[i].EBayeskterm, VAROFINTERESTZJETS[i].ESVDkterm);
    }
    */
    
    
}


void FuncUnfold(string leptonFlavor, bool isMuon, string unfAlg, string variable, int UsedKtermBayes, int UsedKtermSVD, bool doFlat, bool doVarWidth)
{
    int JetPtMin(30);
    int JetPtMax(0);
    
    string energy = getEnergy();
    string outputDirectory = "PNGFiles/FinalUnfoldMeanNJets/";
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    string command = "mkdir -p " + outputDirectory;
    system(command.c_str());
    //if (variable.find("JetPt")!= string::npos && variable.find("Highest") == string::npos && JetPtMin == 20) JetPtMin = 15;
    //if (variable.find("JetPt")!= string::npos && variable.find("Highest") == string::npos && JetPtMin == 30) JetPtMin = 20;
    
    
    //--- set number of toys used in getErrorRespToy2D()
    int NumberOfToys(2000), oppNumberOfToys(4);
    int UsedKterm = UsedKtermBayes, oppUsedKterm = UsedKtermSVD;
    string oppUnfAlg = "SVD";
    if (unfAlg == "SVD") {
        NumberOfToys = 4;
        oppNumberOfToys = 2000;
        oppUnfAlg = "Bayes";
        UsedKterm = UsedKtermSVD;
        oppUsedKterm = UsedKtermBayes;
    }
    
    
    //====== set efficiency errors ===============
    double muonIDIsoHLTError = 0.025;// 2012 : id,iso = 0.005, 0.002 --> 1.4% for dimuon + HLT: 2% ??????
    double electronIDIsoHLTError = 0.005; // SF factors for 2012: 0.002
    double EffError = muonIDIsoHLTError;
    if (!isMuon) EffError = electronIDIsoHLTError;
    //--------------------------------------------
    
    cout << " We now unfold:  " << variable << " with jet pt cut of:" << JetPtMin << " - " << JetPtMax << endl;
    
    //-- fetch the data files and histograms --------------
    cout << " Fetch Data files " << endl;
    TFile *fData[5] = {NULL};  // 0 = central, 1 = JES Up, 2 = JES Down , 3 - SF Up , 4 - SF down
    TH2D *hData[5] = {NULL};
    getFilesNew(FILESDIRECTORY, fData, leptonFlavor, energy, ProcessInfo[DATAFILENAME].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, 0, 0, 0, 0, -1, false);
    get2DHistos(hData, fData, variable, false);
    
    
    
    //-- fetch the WJets files and histograms ----------------
    cout << " Fetch WJets files " << endl;
    TFile *fDYMadGraph[8] = {NULL}; // WJets File : 0 = central, 1 = PU Up,  2 = PU Down,  3 = JER Up, 4 = BtagEff Up, 5 = BtagEff Down, 6 = Wb, 7 = Resp
    TH2D *hDY[8] = {NULL}; // Reco Wjets Hist => actually these are not used
    TH2D *hDYGenMadGraph = NULL; // GenMadGraph Hist
    
    getFilesNew(FILESDIRECTORY, fDYMadGraph, leptonFlavor, energy, ProcessInfo[24].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, 0, 0, 0, 0, -1, false);
    get2DHistos(hDY, fDYMadGraph, variable, false); // Reco Wjets Hist => actually these are not used
    string genName = "gen" + variable;
    hDYGenMadGraph = get2DHisto(fDYMadGraph[0], genName);
    cout << "line th " << __LINE__ << endl;
    //-----------------------------------------------------
    
    
    //--- Get MadGraph Resp Matrix ------------------------
    RooUnfoldResponse *resDY[8] = {NULL};
    getRespsNew(resDY, fDYMadGraph, "response" + variable);
    cout << "line th " << __LINE__ << endl;
    //-----------------------------------------------------
    
    
    //--- fetch the BG files and histograms ----------------
    cout << " Fetch background files " << endl;
    
    const int NBkgSyst(7);
    int nFilesBkg(0);
    
    TFile *fBG[15][NBkgSyst] = {{NULL}};  // 0 = central, 1 = PU Up,  2 = PU Down,  3 = XS Up,  4 = XS Down,  5 = BtagEff Up,  6 = BtagEff Down
    TH2D *hBG[15][NBkgSyst] = {{NULL}}, *hSumBG[NBkgSyst] = {NULL} ;
    
    int nFiles = NFILESDYJETS; //
    bool isDoubleLep(1);
    if (leptonFlavor == "SMuE" || leptonFlavor == "SMu" || leptonFlavor == "Muon") {
        isDoubleLep = 0;
        //nFiles = NFILESTTBAR; //from "Includes/fileNames.h" NFILESTTBAR(9)=>(used the merged files for Top and DY)
        nFiles = NFILESWJETS_M;
    }
    
    int countFiles = 0;
    
    for (int i(0); i < nFiles; i++){
        int fileSelect = FilesDYJets[i];
        //if (!isDoubleLep) fileSelect = FilesTTbarWJets[i]; //FilesTTbarWJets[NFILESTTBAR] = {0, 4, 3, 2, 23, 20, 22, 11, 24};
        if (!isDoubleLep) fileSelect = FilesWJets_M[i];
        string fileNameTemp = ProcessInfo[fileSelect].filename;
        
        // Tau is not considered as background. It is taken into account by the unfolding
        if (fileNameTemp.find("Tau") != string::npos) continue;
        
        // we want to fetch background here, so ignore data files
        if ((fileNameTemp.find("Data") != string::npos) && (fileNameTemp.find("QCD") == string::npos)) continue;
        
        // DYJets is not a background, it is the signal (for double lepton)
        if (isDoubleLep && fileNameTemp.find("DYJets") != string::npos) continue;
        
        // WJets is not a background, it is the signal (for single lepton)
        if (!isDoubleLep && fileNameTemp.find("WJets") != string::npos && fileNameTemp.find("UNF") != string::npos) continue;
        
        // change the name of QCD files got from ProcessInfo[fileSelect].filename
        if ((fileNameTemp.find("QCD") != string::npos) && (variable.find("MeanNJets") != string::npos)) fileNameTemp = "DataQCD_ForMeanNJets_dR_5311";
        
        nFilesBkg++;
        
        //--- fetch the BG file and histos ---
        getFilesNew(FILESDIRECTORY, fBG[countFiles], leptonFlavor, energy, fileNameTemp, JetPtMin, JetPtMax, doFlat, doVarWidth, 0, 0, 0, 0, -1, false);
        get2DHistos(hBG[countFiles], fBG[countFiles], variable, false);
        cout << "line th " << __LINE__ << endl;
        //------------------------------------
        
        for (int j(0); j < NBkgSyst; j++){
            if (countFiles == 0) hSumBG[j] = (TH2D*) hBG[countFiles][j]->Clone();
            else hSumBG[j]->Add(hBG[countFiles][j]);
            cout << " Adding 2D Background of " << fileNameTemp << " : Integral : " << hSumBG[j]->Integral() << endl;
        }
        countFiles++;
    }
    //--- End fetch the BG files and histograms ------------------
    
    cout << " loaded all files: 1 Data : 1 Wjets : " <<  nFilesBkg << " Background files" << endl;
    
    // set the output
    string outputRootFileName = outputDirectory + leptonFlavor + "_" + energy + "_unfolded_" + variable + "_histograms_" + unfAlg;
    if (doVarWidth) outputRootFileName += "_VarWidth";
    outputRootFileName += ".root";
    TFile* outputRootFile = new TFile(outputRootFileName.c_str(), "RECREATE");
    
    /*
    //-- save data reco
    outputRootFile->cd(); hData[0]->Write("Data"); hData[1]->Write("DataJESup"); hData[2]->Write("DataJESdown");
    //-- save SumBG reco
    outputRootFile->cd(); hSumBG[0]->Write("RecoBg"); hSumBG[1]->Write("RecoBgPUup"); hSumBG[2]->Write("RecoBgPUdown");
                          hSumBG[3]->Write("RecoBgXSECup"); hSumBG[4]->Write("RecoBgXSECdown");
    //-- save WJETs reco
    outputRootFile->cd(); hDY[0]->Write("RecoWJETs"); hDY[1]->Write("RecoWJETsPUup"); hDY[2]->Write("RecoWJETsPUdown");
                          hDY[3]->Write("RecoWJETsJERup");
    cout << " saved data reco and madgraph files" << endl;
     */
    
    //--- produce projected 1D hist for gen ----
    TH1D* genMadMeanJ = NULL;
    genMadMeanJ = computeProject1DMeanNJets(hDYGenMadGraph, genName);
    outputRootFile->cd();
    genMadMeanJ->Write("genMad"); //-- save madgraph gen into the file
    //--- End produce projected 1D hist for gen ----
    
    
    
    //--- unfolding 12 different cases
    cout << "Start unfolding" << endl;
    string hNames[12] = {"Central", "JESup", "JESdown", "PUup", "PUdown", "XSECup", "XSECdown", "JERup", "BtagEFFup", "BtagEFFdown", "WBup", "RESPup"};
    int SelData[12] = {0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int SelDY[12]   = {0, 0, 0, 1, 2, 0, 0, 3, 4, 5, 6, 7};
    int SelBG[12]   = {0, 0, 0, 1, 2, 3, 4, 0, 5, 6, 0, 0};
    for ( int n(0); n < 12; n++){
        //TH1D* hUnfolded = Unfold(unfAlg, resDY[SelDY[i]], hData[SelData[i]], hSumBG[SelBG[i]],  UsedKterm, hNames[i]);
        cout << "Start unfolding histograms on : " << hNames[n] <<endl;
        TH2D *hMeas = (TH2D*) hData[SelData[n]]->Clone();
        hMeas->Add(hSumBG[SelBG[n]], -1);
        // set negative bin to zero
        setNegBinZero(hMeas);

        TH2D *hUnfoldedC = NULL;
        RooUnfoldBayes unfoldC(resDY[SelDY[n]], hMeas, UsedKterm);
        unfoldC.SetVerbose(0);
        hUnfoldedC = (TH2D*) unfoldC.Hreco();
        
        cout << "Start computing MeanNjets histograms on : " << hNames[n] <<endl;
        TH1D* hUnfoldedCMeanJ = NULL;
        hUnfoldedCMeanJ = computeProject1DMeanNJets(hUnfoldedC, variable);
        hUnfoldedCMeanJ->SetName(hNames[n].c_str());
        outputRootFile->cd();
        hUnfoldedCMeanJ->Write();
    }
    
    
    //----- get uncertainty contribution from Toy experiments on Response matrix (MC)
    TH2D *hMeasForMC = (TH2D*) hData[0]->Clone();
    hMeasForMC->Add(hSumBG[0], -1);
    // set negative bin to zero
    setNegBinZero(hMeasForMC);
    TH1D* hErrorFromMC = getErrorRespToy2D(NumberOfToys, hMeasForMC, resDY[0], UsedKterm, variable);
    hErrorFromMC->SetName("MC");
    outputRootFile->cd();
    hErrorFromMC->Write();
    //----- End get uncertainty (MC) ----------------------------
    
    
    
    //----- get uncertainty contribution from TTbar modeling
    // Get Scale Factor
    TH2D* hisSF = SFAnalyzeTTbarMeanNJets("SMu", variable, JetPtMin, JetPtMax, doFlat, doVarWidth, 0, false, false, 0, 2);
    
    cout << "Start unfolding histograms on : " << "TTBARup" << endl;
    TH2D *hMeasForTT = (TH2D*) hData[0]->Clone();
    TH2D *hSumBGForTT = NULL;
    
    for (int i(0); i < nFilesBkg; i++){
        
        if (i == 0) hSumBGForTT = (TH2D*) hBG[0][0]->Clone();
        else if (i == 10){
            TH2D *TempScaledBG = (TH2D*) hBG[i][0]->Clone();
            TempScaledBG->Multiply(hisSF);
            hSumBGForTT->Add(TempScaledBG);
        }
        else hSumBGForTT->Add(hBG[i][0]);
    }
    
    hMeasForTT->Add(hSumBGForTT, -1);
    // set negative bin to zero
    setNegBinZero(hMeasForTT);
    
    TH2D *hUnfoldedCtt = NULL;
    RooUnfoldBayes unfoldCtt(resDY[0], hMeasForTT, UsedKterm);
    unfoldCtt.SetVerbose(0);
    hUnfoldedCtt = (TH2D*) unfoldCtt.Hreco();
    
    cout << "Start computing MeanNjets histograms on : " << "TTBARup" << endl;
    TH1D* hUnfoldedCMeanJtt = NULL;
    hUnfoldedCMeanJtt = computeProject1DMeanNJets(hUnfoldedCtt, variable);
    hUnfoldedCMeanJtt->SetName("TTBARup");
    outputRootFile->cd();
    hUnfoldedCMeanJtt->Write();
    //----- End get uncertainty (ttbar) ----------------------------
    
    
    //-- Close all the files ------------------------------
    outputRootFile->Close();
    closeFilesNew(fData);
    closeFilesNew(fDYMadGraph);
    for (int i(0); i < nFilesBkg; i++) closeFilesNew(fBG[i]);
    //-----------------------------------------------------
}




TH1D* getErrorRespToy2D(int myNToys, TH2D* hMeas2D, RooUnfoldResponse* response, int kterm, string variable)
{
    cout << endl;
    cout << "--- Begin Toy experiment on Response matrix ( for MeanNJets Dist ) ------" << endl;
    
    //--- test distribution of unfolded hist
    int whichBin = 8; // which bin to be test
    RooUnfoldBayes unfoldCTest(response, hMeas2D, kterm);
    unfoldCTest.SetVerbose(0);
    TH2D *hUnfoldedCTest2D = (TH2D*) unfoldCTest.Hreco();
    TH1D* hUnfoldedCTest = NULL;
    hUnfoldedCTest = computeProject1DMeanNJets(hUnfoldedCTest2D, "test");
    
    double meanC = hUnfoldedCTest->GetBinContent(whichBin);
    double sigmaC = hUnfoldedCTest->GetBinError(whichBin);
    for (int i = 1; i <= hUnfoldedCTest->GetNbinsX(); i++){
        cout << "bin: " << i << " " << hUnfoldedCTest->GetBinContent(i) << " " << hUnfoldedCTest->GetBinError(i) << endl;
    }
    
    TH1D *distriUnfolded = new TH1D ("ErrTest for Smeared Resp", "ErrTest for Smeared Resp", 20, -3.0, 3.0);
    //-----------
    
    
    //--- output setup
    double uncerFromToyResp[50];
    TH1D * hUncerFromToyResp = (TH1D*) hUnfoldedCTest->Clone("Error from Resp Toy (MC)");
    
    double xiSum[50];
    double xiSquareSum[50];
    for(int i = 0; i < 50; i++){
        xiSum[i] = 0;
        xiSquareSum[i] = 0;
    }
    
    int numBinXFinal = hUnfoldedCTest->GetNbinsX();
    
    for (int nthToy = 0; nthToy < myNToys; nthToy ++){
        
        if (nthToy % 100 == 0) cout << "looping nthToy: " << nthToy << endl;
        
        //--- generate new RooUnfoldBayes object with smeared Response matrix
        RooUnfoldResponse *newRespToy =  (RooUnfoldResponse*) response->RunToy();
        RooUnfoldBayes unfoldToyResp(newRespToy, hMeas2D, kterm);
        unfoldToyResp.SetVerbose(0);
        TH2D *hUnfoldedToyResp2D = NULL;
        hUnfoldedToyResp2D = (TH2D*) unfoldToyResp.Hreco();
        
        TH1D* hUnfoldedToyResp = NULL;
        hUnfoldedToyResp = computeProject1DMeanNJets(hUnfoldedToyResp2D, "tempt");
        
        //--- next calculate the uncertainty
        for (int i = 1; i <= numBinXFinal; i++){
            double xi = hUnfoldedToyResp->GetBinContent(i);
            //cout << " k: " << nthToy << " binth : " << i << " xi: " << xi << endl;
            if (xi < 0){
                cout << " unfolded value less than Zero " << nthToy << " " << i << endl;
                xi = 0;
            }
            xiSum[i] += xi;
            xiSquareSum[i] += xi * xi;
            
            //--- test distribution of unfolded hist ( whichBin = 8 only)
            if (i == whichBin) {
                double distriVal = (xi - meanC)/sigmaC;
                distriUnfolded->Fill(distriVal);
            }
        }
        
    }
    
    //--- test distribution of unfolded hist
    TCanvas* canvasErrTest = new TCanvas("canvasErrTest for Smeared Resp", "canvasErrTest for Smeared Resp");
    canvasErrTest->cd();
    distriUnfolded->Draw();
    
    string MCoutputDir = "PNGFiles/MCsystematics/";
    string command2 = "mkdir -p " + MCoutputDir;
    system(command2.c_str());
    string testNamePDF = "PNGFiles/MCsystematics/TestErrorToyResp_" + variable + ".pdf";
    canvasErrTest->Print(testNamePDF.c_str());
    //-------------
    
    //--- output
    for (int i = 1; i <= numBinXFinal; i++){
        //uncerFromToyResp[i] = sqrt( fabs( (xiSquareSum[i]/myNToys) - (pow( xiSum[i]/myNToys, 2)) ) );
        uncerFromToyResp[i] = sqrt( (xiSquareSum[i]/myNToys) - (pow( xiSum[i]/myNToys, 2)) );
        hUncerFromToyResp->SetBinContent(i, uncerFromToyResp[i]);
        hUncerFromToyResp->SetBinError(i, 0.0);
        //cout << " sumSquare: " << xiSquareSum[i]/myNToys << " squareOfSum: " << pow( xiSum[i]/myNToys, 2) <<  " diff: " << (xiSquareSum[i]/myNToys) - (pow( xiSum[i]/myNToys, 2)) << endl;
    }
    
    cout << " stat uncer of central (8th binth): " << sigmaC << endl;
    cout << " syts uncer from MC (8th binth)   : " << uncerFromToyResp[whichBin] << endl;
    cout << "....." << endl;
    for (int i = 1; i <= hUncerFromToyResp->GetNbinsX(); i++){
        cout << "absolute error for binth " << i << " " << hUncerFromToyResp->GetBinContent(i) << endl;
    }
    
    cout << "--- End Toy experiment on Response matrix ---" << endl;
    cout << endl;
    
    return hUncerFromToyResp;
}

TH1D* computeProject1DMeanNJets( TH2D *hUnfoldedC, string variable )
{
    
    int nBinsX(hUnfoldedC->GetNbinsX());
    int nBinsY(hUnfoldedC->GetNbinsY());
    double xmin = hUnfoldedC->GetXaxis()->GetXmin();
    double xmax = hUnfoldedC->GetXaxis()->GetXmax();
    //cout << " nBinsX: " << nBinsX << endl;
    //cout << " nBinsY: " << nBinsY << endl;
    
    const double *xBins = new double[nBinsX+1];
    xBins = hUnfoldedC->GetXaxis()->GetXbins()->GetArray();
    
    TH1D* hUnfoldedCMeanJ;
    
    if (xBins == 0){
        hUnfoldedCMeanJ = new TH1D(variable.c_str(), variable.c_str(), nBinsX, xmin, xmax);
    }
    else{
        hUnfoldedCMeanJ = new TH1D(variable.c_str(), variable.c_str(), nBinsX, xBins);
    }
    
    //--- produce projected 1D hist for hUnfoldedC ---
    for (int i=0; i<nBinsX; i++){
        // calculate mean # of jets
        double meanNJets(0);
        double NEvents(0);
        double NJetsTot(0);
        for(int j=0; j<nBinsY; j++){
            double NEventsTemp(0);
            double NJetsTemp(0);
            NEventsTemp = hUnfoldedC->GetBinContent(i+1,j+1);
            NEvents = NEvents + NEventsTemp;
            
            NJetsTemp = (hUnfoldedC->GetBinContent(i+1,j+1)) * (j+1);
            NJetsTot = NJetsTot + NJetsTemp;
        }
        if (NEvents > 0) meanNJets = NJetsTot/NEvents;
        
        hUnfoldedCMeanJ->SetBinContent(i+1, meanNJets);
        
        // calculate error
        double errOfMean(0);
        if (NEvents > 0) {
            for(int k=0; k<nBinsY; k++){
                double errTemp(0);
                errTemp = (((NEvents * (k+1)) - NJetsTot) / pow(NEvents, 2)) * (hUnfoldedC->GetBinError(i+1,k+1));
                errTemp = pow(errTemp, 2);
                errOfMean = errOfMean + errTemp;
            }
            errOfMean = sqrt(errOfMean);
        }
        hUnfoldedCMeanJ->SetBinError(i+1, errOfMean);
        //cout << " loop " << i << endl;
    }
    
    return hUnfoldedCMeanJ;
    
}

TH2D* SFAnalyzeTTbarMeanNJets(string leptonFlavor, string variable, int JetPtMin, int JetPtMax, bool doFlat, bool doVarWidth, int doQCD, bool doSSign, bool doInvMassCut, int MET, int doBJets)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("1.2f");
    
    cout << "...>>>" << endl;
    cout << "...>>> Processing Data and MC from TTbar Systematic Control Region ; Variable is: " << variable << endl;
    
    string energy = getEnergy();
    string HistoFilesDirTT = "HistoFilesTT/";
    
    int nsel = 0;
    const int *sel = NULL;
    if (leptonFlavor == "DMu" || leptonFlavor == "DE") {
        nsel = NFILESDYJETS;
        sel = new int[nsel];
        sel = FilesDYJets;
    }
    else if (leptonFlavor == "SMu" || leptonFlavor == "SE") {
        //nsel = NFILESTTBAR;
        nsel = NFILESWJETS_M;
        sel = new int[nsel];
        //sel = FilesTTbarWJets;
        sel = FilesWJets_M;
    }
    
    //--- fetch the files and histograms of the selected variable --------------
    // note that now we need 2D histo for MeanNJets
    TH2D **histoArr = new TH2D*[nsel];
    for (int i(0); i < nsel; i++) {
        TFile *fCurrent;
        string fileNameTemp = ProcessInfo[sel[i]].filename;
        
        // change the name of QCD files got from ProcessInfo[sel[i]].filename
        if ((fileNameTemp.find("QCD") != string::npos) && (variable.find("MeanNJets") != string::npos) ) fileNameTemp = "DataQCD_ForMeanNJets_dR_5311";
        
        fCurrent = getFile(HistoFilesDirTT, leptonFlavor, energy, fileNameTemp, JetPtMin, JetPtMax, doFlat, doVarWidth, doQCD , doSSign,  doInvMassCut, MET, doBJets,  "", "0");
        
        std::cout << "opened :  " << i << "   " << sel[i] << "   " << FilesWJets_M[i] << "  " << fileNameTemp << "   " << std::endl;
        
        histoArr[i] = get2DHisto(fCurrent, variable);
        // close all input root files
        fCurrent->Close();
    }
    std::cout << "Closed all files for ttbar" << std::endl;
    
    //--- calculate hMeasTTbar ---------
    TH2D *hMeasTTbar = (TH2D*) histoArr[0]->Clone();
    int nBinsForX(hMeasTTbar->GetNbinsX());
    int nBinsForY(hMeasTTbar->GetNbinsY());
    
    int nselTTbar = 11; // index 11 in sel[] is ttbar
    for (int i = 1; i < nsel; i++) {
        if (i == nselTTbar) continue; // we want data minus WjetsMC and all bgs but ttbar
        hMeasTTbar->Add(histoArr[i], -1);
        cout << " subtracting : " << ProcessInfo[sel[i]].filename << endl;
        //cout << hMeasTTbar->GetBinContent(9,4) << " " << hMeasTTbar->GetBinError(9,4) << endl;
    }
    setNegBinZero(hMeasTTbar);
    
    TH2D *histScaleF = (TH2D*) hMeasTTbar->Clone();
    histScaleF->Divide(histoArr[nselTTbar]);
    
    // The (0-jet and) 1-jet bins are irrelevant. So set them to 1
    for (int i = 1; i <= nBinsForX; i++) {
        histScaleF->SetBinContent(i, 1, 1.0);
        histScaleF->SetBinError(i, 1, 0.0);
    }
    
    // Print out SF for inspection
    /*
     for (int j = 1; j <= nBinsForY; j++) {
     for (int i = 1; i <= nBinsForX; i++) {
     //cout << hMeasTTbar->GetBinContent(j)  << "\t" << hMeasTTbar->GetBinError(j)  << "\t";
     //cout << histoArr[7]->GetBinContent(j) << "\t" << histoArr[7]->GetBinError(j) << "\t";
     //cout << histScaleF->GetBinContent(j)  << "\t" << histScaleF->GetBinError(j)  << "\t";
     //cout << hMeasTTbar->GetBinContent(j)/histoArr[7]->GetBinContent(j) ;
     cout << histScaleF->GetBinContent(i,j)  << "\t";
     }
     cout << endl;
     }
     */
    
    cout << "...>>> Scale Factors are obtained ..." << endl;
    cout << "...>>>" << endl;
    
    TCanvas *canSF = new TCanvas("can1", "can1", 900, 600);
    canSF->cd();
    histScaleF->GetZaxis()->SetRangeUser(0.51, 1.49);
    histScaleF->GetZaxis()->SetTitle("Scale Factor");
    histScaleF->SetTitle(("TTbarSyst for " + variable).c_str());
    histScaleF->SetFillColor(kYellow);
    histScaleF->Draw("LEGO1 0");
    
    string TTbaroutputDir = "PNGFiles/TTbarSystematics/";
    string command1 = "mkdir -p " + TTbaroutputDir;
    system(command1.c_str());
    
    string SFplotName = "PNGFiles/TTbarSystematics/SFplots_" + variable + ".pdf";
    canSF->Print(SFplotName.c_str());
    
    return histScaleF;
    
}

void setNegBinZero(TH2D *histograms2D)
{
    int nBinsX(histograms2D->GetNbinsX());
    int nBinsY(histograms2D->GetNbinsY());
    
    for (int i = 1; i <= nBinsX; i++){
        for(int j = 1; j <= nBinsY; j++){
            if ( histograms2D->GetBinContent(i,j) < 0 ){
                histograms2D->SetBinContent(i, j, 0.);
                histograms2D->SetBinError(i, j, 0.);
            }
        }
    }
    
}



void get2DHistos(TH2D *histograms2D[], TFile *Files[], string variable, bool isDoubleLep)
{
    string fileName = Files[0]->GetName();
    int nFiles;
    if ((fileName.find("Data") != string::npos) && (fileName.find("QCD") == string::npos)){
        nFiles = 3;
        if (fileName.find("DE") != string::npos) nFiles = 5;
    }
    else if (((isDoubleLep && fileName.find("DYJets") != string::npos) || (!isDoubleLep && fileName.find("WJets") != string::npos)) && fileName.find("UNFOLDING") != string::npos) nFiles = 8;
    // for BG
    else nFiles = 7;
    
    for (int i(0); i < nFiles; i++){
        histograms2D[i] = (TH2D*) Files[i]->Get(variable.c_str());
        cout << "Integral() " << histograms2D[i]->Integral() << endl;
        //histograms2D[i]->SetDirectory(0);
    }
}

void getRespsNew(RooUnfoldResponse *responses[], TFile *Files[], string variable)
{
    string fileName = Files[0]->GetName();
    int nFiles;
    //if (fileName.find("Data") != string::npos) nFiles = 3;
    if ((fileName.find("Data") != string::npos) && (fileName.find("QCD") == string::npos)) nFiles = 3;
    else if (fileName.find("UNFOLDING") != string::npos) nFiles = 8;
    // for BG
    else nFiles = 7;
    
    for (int i(0); i < nFiles; i++){
        responses[i] = (RooUnfoldResponse*) Files[i]->Get(variable.c_str());
        cout << "Getting Response with Dimension (Mea x Truth): " << responses[i]->GetDimensionMeasured() << " x " << responses[i]->GetDimensionTruth() << endl;
    }
}

void getFilesNew(string histoFilesDirectory, TFile *Files[], string leptonFlavor, string energy, string Name, int JetPtMin, int JetPtMax, bool doFlat, bool doVarWidth, int doQCD, bool doSSign, bool doInvMassCut, int MET, int doBJets, bool useUnfoldingFiles)
{
    
    //--- make sure leptonFlavor is short version ---
    if (leptonFlavor == "Muons" || leptonFlavor == "DMu_") leptonFlavor = "DMu";
    else if (leptonFlavor == "Electrons" || leptonFlavor == "DE_") leptonFlavor = "DE";
    else if (leptonFlavor == "Electron" || leptonFlavor == "SE_") leptonFlavor = "SE";
    else if (leptonFlavor == "Muon" || leptonFlavor == "SMu_") leptonFlavor = "SMu";
    else if (leptonFlavor == "MuonElectron" || leptonFlavor == "SMuE_") leptonFlavor = "SMuE";
    //-----------------------------------------------
    
    //--- set the double lepton flag ---
    bool isDoubleLep(0);
    if (leptonFlavor == "DMu" || leptonFlavor == "DE") isDoubleLep = 1;
    //----------------------------------
    
    vector<string> Syst;
    if ((Name.find("Data") != string::npos) && (Name.find("QCD") == string::npos)) { // for data we have:
        Syst.push_back("0");                 //   0: central
        Syst.push_back("2_Up");              //   2 up: JES up
        Syst.push_back("2_Down");            //   2 down: JES down
        if (leptonFlavor == "DE") {          // additionaly for electron:
            Syst.push_back("5_Up");          //   5 up: scale factor up
            Syst.push_back("5_Down");        //   5 down: scale factor down
        }
    }
    else if (Name.find("UNFOLDING") != string::npos && ((isDoubleLep && Name.find("DYJets") != string::npos) || (!isDoubleLep && Name.find("WJets") != string::npos))) {
        // for DYJets in case of Z+Jets or for WJets in case of W+Jets analysis we have:
        Syst.push_back("0");         // 0:      central
        Syst.push_back("1_Up");      // 1 up:   PU up
        Syst.push_back("1_Down");    // 1 down: PU down
        Syst.push_back("4_Up");      // 4 up:   JER up
        Syst.push_back("6_Up");      // 6 up:   BtagEff up
        Syst.push_back("6_Down");    // 6 down: BtagEff down
        Syst.push_back("7_Up");      // 7 up:   Wb
        Syst.push_back("8_Up");      // 8 up:   Resp
    }
    else { // for background we have
        Syst.push_back("0");         // 0:      central
        Syst.push_back("1_Up");      // 1 up:   PU up
        Syst.push_back("1_Down");    // 1 down: PU down
        Syst.push_back("3_Up");      // 3 up:   XSec up
        Syst.push_back("3_Down");    // 3 down: Xsec down
        Syst.push_back("6_Up");      // 6 up:   BtagEff up
        Syst.push_back("6_Down");    // 6 down: BtagEff down
    };
    
    //--- determnie how many files we have and open them all ---
    int nSyst(Syst.size());
    for (int i(0); i < nSyst; i++) {
        Files[i] = getFile(histoFilesDirectory, leptonFlavor, energy, Name, JetPtMin, JetPtMax, doFlat, doVarWidth, doQCD, doSSign, doInvMassCut, MET, doBJets, "", Syst[i], false, useUnfoldingFiles);
    }
    //----------------------------------------------------------
}

void closeFilesNew(TFile *Files[])
{
    if (Files[0]) {
        string fileName = Files[0]->GetName();
        int nFiles;
        //if (fileName.find("Data") != string::npos) {
        if ((fileName.find("Data") != string::npos) && (fileName.find("QCD") == string::npos)){
            nFiles = 3;
            if (fileName.find("DE") != string::npos) nFiles = 5;
        }
        else if (fileName.find("UNFOLDING") != string::npos) nFiles = 8;
        // for BG
        else nFiles = 7;
        
        for (int i(0); i < nFiles; i++){
            closeFile(Files[i]);
        }
    }
}

TH2D* get2DHisto(TFile *File, string variable)
{
    TH2D *histo2D = (TH2D*) File->Get(variable.c_str());
    histo2D->SetDirectory(0);
    return histo2D;
}





