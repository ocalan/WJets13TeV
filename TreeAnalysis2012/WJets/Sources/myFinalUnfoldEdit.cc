#include <vector>
#include <TDirectory.h>
#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TPaveStats.h>
#include <TText.h>
#include <TLegend.h>
#include <iostream>
#include <sstream>
#include <TLorentzVector.h>
#include <RooUnfoldResponse.h>
#include <RooUnfold.h>
#include <RooUnfoldBayes.h>
#include <RooUnfoldSvd.h>
#include <RooUnfoldBinByBin.h>
#include <RooUnfoldInvert.h>
#include <TSVDUnfold.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TDatime.h>
#include <TRandom3.h>
#include <TMatrixD.h>
#include "myUnfoldingSyst.h"
#include "getFilesAndHistograms.h"
#include "myFinalUnfoldEdit.h"


using namespace std;

//--  Setting global variables --------------------------------------------------------------
#include <TRandom3.h>
#include <TMath.h>
#include "fileNames.h"
#include "variablesOfInterestVarWidth.h"
TRandom3* RandGen = new TRandom3();
//-------------------------------------------------------------------------------------------

TH1D* getErrorRespToy(int myNToys = 1000, TH1D* hMeas = NULL, RooUnfoldResponse* response = NULL, int kterm = 4);
TH1D* getErrorMeasToy(int myNToys = 1000, TH1D* hMeas = NULL, RooUnfoldResponse* response = NULL, int kterm = 4);

void myFinalUnfoldEdit(int start, int end, string leptonFlavor, string unfAlg)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    //--- check whether muon or electron is used ---
    bool isMuon;
    if (leptonFlavor.find("Mu") != string::npos) isMuon = 1;
    else isMuon = 0;
    //---------------------------------------------
    
    if (end == -1) end = start + 1;
    for (int i(start); i < end /*NVAROFINTERESTZJETS*/; i++) {
        for (int j(0); j < 1; j++){
            if (leptonFlavor == "DMu" || leptonFlavor == "DE"){
                if (isMuon) FuncUnfold(leptonFlavor, isMuon, unfAlg, VAROFINTERESTZJETS[i].name, VAROFINTERESTZJETS[i].MuBayeskterm, VAROFINTERESTZJETS[i].MuSVDkterm);
                else FuncUnfold(leptonFlavor, isMuon, unfAlg, VAROFINTERESTZJETS[i].name, VAROFINTERESTZJETS[i].EBayeskterm, VAROFINTERESTZJETS[i].ESVDkterm);
            }
            else if (leptonFlavor == "SMu" || leptonFlavor == "SE"){
                if (isMuon) FuncUnfold(leptonFlavor, isMuon, unfAlg, VAROFINTERESTWJETS[i].name, VAROFINTERESTWJETS[i].MuBayeskterm, VAROFINTERESTWJETS[i].MuSVDkterm);
                else FuncUnfold(leptonFlavor, isMuon, unfAlg, VAROFINTERESTWJETS[i].name, VAROFINTERESTWJETS[i].EBayeskterm, VAROFINTERESTWJETS[i].ESVDkterm);
            }
        }
    }
    
}


void FuncUnfold(string leptonFlavor, bool isMuon, string unfAlg, string variable, int UsedKtermBayes, int UsedKtermSVD, bool doFlat, bool doVarWidth)
{
    int JetPtMin(30);
    int JetPtMax(0);
    
    string energy = getEnergy();
    string outputDirectory = "PNGFiles/FinalUnfold/";
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    if (JetPtMin == 30) outputDirectory = "PNGFiles/FinalUnfold_inc_30_1000_Toys/";
    else if (JetPtMin == 20) outputDirectory = "PNGFiles/FinalUnfold_inc_20_1000_Toys/";
    string command = "mkdir -p " + outputDirectory;
    system(command.c_str());
    //if (variable.find("JetPt")!= string::npos && variable.find("Highest") == string::npos && JetPtMin == 20) JetPtMin = 15;
    //if (variable.find("JetPt")!= string::npos && variable.find("Highest") == string::npos && JetPtMin == 30) JetPtMin = 20;
    
    
    int NumberOfToys(2000), oppNumberOfToys(4);
    int UsedKterm = UsedKtermBayes, oppUsedKterm = UsedKtermSVD;
    string oppUnfAlg = "SVD";
    if (unfAlg == "SVD") {
        UsedKterm = UsedKtermSVD;
        NumberOfToys = 3;
        oppNumberOfToys = 400;
        oppUnfAlg = "Bayes";
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
    TH1D *hData[5] = {NULL};
    getFiles(FILESDIRECTORY, fData, leptonFlavor, energy, ProcessInfo[DATAFILENAME].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, 0, 0, 0, 0, -1, false);
    getHistos(hData, fData, variable, false);
    
    
    //-- fetch the WJets files and histograms ----------------
    cout << " Fetch WJets files " << endl;
    TFile *fDYMadGraph[4] = {NULL}; // WJets File : 0 = central, 1 = PU Up,  2 = PU Down,  3 = JER Up
    TH1D *hDY[4] = {NULL}; // Reco Hist
    TH1D *hDYGenMadGraph = NULL; // GenMadGraph Hist
    
    getFiles(FILESDIRECTORY, fDYMadGraph, leptonFlavor, energy, ProcessInfo[24].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, 0, 0, 0, 0, -1, false);
    getHistos(hDY, fDYMadGraph, variable, false);
    hDYGenMadGraph = getHisto(fDYMadGraph[0], "gen" + variable);
    cout << "line th " << __LINE__ << endl;
    //-----------------------------------------------------
    
    
    //--- Get MadGraph Resp Matrix ------------------------
    RooUnfoldResponse *resDY[4] = {NULL};
    getResps(resDY, fDYMadGraph, "response" + variable);
    cout << "line th " << __LINE__ << endl;
    //-----------------------------------------------------
    
    
    //-- fetch the BG files and histograms ----------------
    cout << " Fetch background files " << endl;
    const int NBkgGroups(6);
    const int NBkgSyst(5);
    int nFilesBkg(0);
    
    TFile *fBG[15][NBkgSyst] = {{NULL}};  // 0 = central, 1 = PU Up,  2 = PU Down,  3 = XS Up,  4 = XS Down
    TH1D *hBG[15][NBkgSyst] = {{NULL}}, *hSumBG[NBkgSyst] = {NULL}, *hSumBGgroup[15][10] = {{NULL}};
    
    int group = -1;
    double sumEve = 0.;
    
    int nFiles = NFILESDYJETS; // from "Includes/fileNames.h"
    bool isDoubleLep(1);
    if (leptonFlavor == "SMuE" || leptonFlavor == "SMu" || leptonFlavor == "Muon") {
        isDoubleLep = 0;
        nFiles = NFILESTTBAR; //NFILESTTBAR(9);
    }
    
    int countFiles = 0;
    int countBGp[2] = {-1,-1};
    
    for (int i(0); i < nFiles; i++){
        int fileSelect = FilesDYJets[i];
        if (!isDoubleLep) fileSelect = FilesTTbarWJets[i]; //FilesTTbarWJets[NFILESTTBAR] = {0, 4, 3, 2, 23, 20, 22, 11, 24};
        string fileNameTemp = ProcessInfo[fileSelect].filename;
        
        // Tau is not considered as background. It is taken into account by the unfolding
        if (fileNameTemp.find("Tau") != string::npos) continue;
        
        if ((fileNameTemp.find("Data") != string::npos) && (fileNameTemp.find("QCD") == string::npos)) continue;
        //if (fileNameTemp.find("Data") != string::npos) continue;
        
        // DYJets is not a background, it is the signal (for double lepton)
        if (isDoubleLep && fileNameTemp.find("DYJets") != string::npos) continue;
        
        // WJets is not a background, it is the signal (for single lepton)
        if (!isDoubleLep && fileNameTemp.find("WJets") != string::npos && fileNameTemp.find("UNF") != string::npos) continue;
        
        nFilesBkg++;
        
        //--- fetch the BG file and histos ---
        getFiles(FILESDIRECTORY, fBG[countFiles], leptonFlavor, energy, fileNameTemp, JetPtMin, JetPtMax, doFlat, doVarWidth, 0, 0, 0, 0, -1, false);
        getHistos(hBG[countFiles], fBG[countFiles], variable, false);
        cout << "line th " << __LINE__ << endl;
        //------------------------------------
        
        if (isDoubleLep){
            if (fileNameTemp.find("ZZ") != string::npos){
                if (countBGp[0] == -1) group++;
            }
            else if (fileNameTemp.find("WZ") != string::npos){
                if (countBGp[1] == -1) group++;
            }
            else group++;
        }
        else group++;
        
        for (int j(0); j < NBkgSyst; j++){
            if (countFiles == 0) hSumBG[j] = (TH1D*) hBG[countFiles][j]->Clone();
            else hSumBG[j]->Add(hBG[countFiles][j]);
            
            //assign hSumBGgroup[j][group]
            if (fileNameTemp.find("ZZ") != string::npos ) {
                if (countBGp[0] == -1) {
                    if (j == 4) countBGp[0] = group;
                    hSumBGgroup[j][group] = (TH1D*) hBG[countFiles][j]->Clone();
                }
                else {
                    hSumBGgroup[j][countBGp[0]]->Add(hBG[countFiles][j]);
                }
            }
            else if (fileNameTemp.find("WZ") != string::npos) {
                if (countBGp[1] == -1) {
                    if (j == 4) countBGp[1] = group;
                    hSumBGgroup[j][group] = (TH1D*) hBG[countFiles][j]->Clone();
                }
                else {
                    hSumBGgroup[j][countBGp[1]]->Add(hBG[countFiles][j]);
                }
            }
            else {
                hSumBGgroup[j][group] = (TH1D*) hBG[countFiles][j]->Clone();
            }
            
            if ( j == 0 ) {
                sumEve += hBG[countFiles][j]->GetBinContent(1);
            }
        }
        countFiles++;
    }
    //-----------------------------------------------------
    
    cout << " loaded all files"<<endl;
    
    // set the output
    string outputRootFileName = outputDirectory + leptonFlavor + "_" + energy + "_unfolded_" + variable + "_histograms_" + unfAlg;
    if (doVarWidth) outputRootFileName += "_VarWidth";
    outputRootFileName += ".root";
    TFile* outputRootFile = new TFile(outputRootFileName.c_str(), "RECREATE");
    
    //-- save data reco
    outputRootFile->cd(); hData[0]->Write("Data"); hData[1]->Write("DataJESup"); hData[2]->Write("DataJESdown");
    //-- save SumBG reco
    outputRootFile->cd(); hSumBG[0]->Write("RecoBg"); hSumBG[1]->Write("RecoBgPUup"); hSumBG[2]->Write("RecoBgPUdown");
                          hSumBG[3]->Write("RecoBgXSECup"); hSumBG[4]->Write("RecoBgXSECdown");
    //-- save WJETs reco
    outputRootFile->cd(); hDY[0]->Write("RecoWJETs"); hDY[1]->Write("RecoWJETsPUup"); hDY[2]->Write("RecoWJETsPUdown");
                          hDY[3]->Write("RecoWJETsJERup");
    
    //-- save madgraph gen also in the file
    outputRootFile->cd(); hDYGenMadGraph->Write("genMad");
    
    cout << " saved data reco and madgraph files" << endl;
    
    // unfolding 8 different cases
    cout << "Start unfolding"<<endl;
    string hNames[8] = {"Central", "JESup", "JESdown", "PUup", "PUdown", "XSECup", "XSECdown", "JERup"};
    int SelData[8] = {0, 1, 2, 0, 0, 0, 0, 0};
    int SelDY[8]   = {0, 0, 0, 1, 2, 0, 0, 3};
    int SelBG[8]   = {0, 0, 0, 1, 2, 3, 4, 0};
    for ( int i(0); i < 8; i++){
        TH1D* hUnfolded = Unfold(unfAlg, resDY[SelDY[i]], hData[SelData[i]], hSumBG[SelBG[i]],  UsedKterm, hNames[i]);
        outputRootFile->cd();  hUnfolded->Write();
        cout << "Start unfolding offseted histograms on RECO " << hNames[i] <<endl;
    }
    
    // get uncertainty contribution from Toy experiments on Response matrix (MC)
    TH1D *hMeas = (TH1D*) hData[0]->Clone();
    hMeas->Add(hSumBG[0], -1);
    TH1D* hErrorFromToyResponse = getErrorRespToy(NumberOfToys, hMeas, resDY[0], UsedKterm);
    outputRootFile->cd();
    hErrorFromToyResponse->Write("MC");
    
    //-- Close all the files ------------------------------
    outputRootFile->Close();
    closeFiles(fData);
    closeFiles(fDYMadGraph);
    for (int i(0); i < nFilesBkg; i++) closeFiles(fBG[i]);
    //-----------------------------------------------------
}


TH1D* getErrorRespToy(int myNToys, TH1D* hMeas, RooUnfoldResponse* response, int kterm)
{
    cout << "--- Begin Toy experiment on Response matrix ------" << endl;
    double uncerFromToyResp[50]; // output
    TH1D * hUncerFromToyResp = (TH1D*) hMeas->Clone("Error from Resp Toy");
    
    //TRandom3* RandGen = new TRandom3();
    RandGen->SetSeed(0); // set random seed
    double xiSum[50];
    double xiSquareSum[50];
    for(int i = 0; i<50; i++){
        xiSum[i] = 0;
        xiSquareSum[i] = 0;
    }
    int numBinXFinal = response->GetNbinsMeasured();
    
    //--- test distribution of unfolded hist
    RooUnfoldBayes unfoldCTest(response, hMeas, kterm);
    unfoldCTest.SetVerbose(0);
    TH1D *hUnfoldedCTest = (TH1D*) unfoldCTest.Hreco();
    double meanC = hUnfoldedCTest->GetBinContent(2);
    double sigmaC = hUnfoldedCTest->GetBinError(2);
    TH1D *distriUnfolded = new TH1D ("ErrTest for Smeared Resp", "ErrTest for Smeared Resp", 20, -0.2, 0.2);
    //-----------
    
    int myCount = 0;
    
    for (int nthToy = 0; nthToy < myNToys; nthToy ++){
        
        // generate smeared Response matrix and resulting unfolded histo. Inputs are "response", hMeas, kterm
        
        // creat new response matrix
        TH2D *hRespToy = (TH2D*) response->Hresponse();
        TH2D *hNewRespToy = (TH2D*) hRespToy->Clone();
        
        //int numBinUsed(0);
        int binthStart = 1;
        int binthEnd = numBinXFinal;
        //numBinUsed = numBinXFinal;
        if (response->UseOverflowStatus()){
            binthStart = 0;
            binthEnd = numBinXFinal +1;
            //numBinUsed = numBinXFinal +2 ;
        }
        
        for (int i = binthStart; i <= binthEnd; i++) {
            for(int j = binthStart; j <= binthEnd; j++){
                double errRespTemp (0);
                double contentRespTemp (0);
                errRespTemp = hRespToy->GetBinError(i,j);
                contentRespTemp = hRespToy->GetBinContent(i,j);
                if (errRespTemp > 0.0 && contentRespTemp > 0){
                    
                    //double tempora = contentRespTemp;
                    
                    contentRespTemp += RandGen->Gaus(0,errRespTemp);
                    if (contentRespTemp < 0){
                        //cout << " New response metrix less than zero. Error: " << tempora << " " << errRespTemp << " " << contentRespTemp  << endl;
                        contentRespTemp = 0;
                        myCount += 1;
                    }
                }
                hNewRespToy->SetBinContent(i,j,contentRespTemp);
            }
        }
        
        // creat new RooUnfoldResponse object
        RooUnfoldResponse *newRespToy = new RooUnfoldResponse(response->Hmeasured(), response->Htruth(), hNewRespToy);
        
        // creat new RooUnfold object and Unfolded histo
        RooUnfoldBayes unfoldToyResp(newRespToy, hMeas, kterm);
        unfoldToyResp.SetVerbose(0);
        TH1D *hUnfoldedToyResp = NULL;
        hUnfoldedToyResp = (TH1D*) unfoldToyResp.Hreco();
        
        // next calculate the uncertainty
        for (int i = 1; i <= numBinXFinal; i++){
            double xi = hUnfoldedToyResp->GetBinContent(i);
            //cout << " k: " << nthToy << " binth : " << i << " xi: " << xi << endl;
            if (xi < 0){
                cout << " unfolded value less than Zero " << nthToy << " " << i << endl;
                xi = 0;
            }
            xiSum[i] += xi;
            xiSquareSum[i] += xi * xi;
            
            //--- test distribution of unfolded hist ( second bin only [exclusive 1 jet] )
            if (i == 2) {
                double distriVal = (xi - meanC)/sigmaC;
                distriUnfolded->Fill(distriVal);
            }
            //---------------
        }
        
    }
    
    cout << " Counts for bins of New Resp < zero : " << myCount << endl;
    
    //--- test distribution of unfolded hist
    //TCanvas* canvasErrTest = new TCanvas("canvasErrTest for Smeared Resp", "canvasErrTest for Smeared Resp");
    //canvasErrTest->cd();
    //distriUnfolded->Draw();
    //canvasErrTest->Print("PNGFiles/TestToyErrorResp.pdf");
    //-------------
    
    // set up output
    for (int i = 1; i <= numBinXFinal; i++){
        uncerFromToyResp[i] = sqrt( (xiSquareSum[i]/myNToys) - (pow( xiSum[i]/myNToys, 2)) );
        //uncerFromToyResp[i] = sqrt( fabs( (xiSquareSum[i]/myNToys) - (pow( xiSum[i]/myNToys, 2)) ) );
        hUncerFromToyResp->SetBinContent(i, uncerFromToyResp[i]);
        hUncerFromToyResp->SetBinError(i, 0.0);
        
        cout << " sumSquare: " << xiSquareSum[i]/myNToys << " squareOfSum: " << pow( xiSum[i]/myNToys, 2) <<  " diff: " << (xiSquareSum[i]/myNToys) - (pow( xiSum[i]/myNToys, 2)) << endl;
        
    }
    
    cout << " stat uncer of central (2nd binth): " << sigmaC << endl;
    cout << " syts uncer from MC (2nd binth)   : " << uncerFromToyResp[2] << endl;
    
    cout << "--- End Toy experiment on Response matrix ---" << endl;
    
    return hUncerFromToyResp;
}

TH1D* getErrorMeasToy(int myNToys, TH1D* hMeas, RooUnfoldResponse* response, int kterm)
{
    //--- Toy experiment on hMeas ------
    double uncerFromToyMeas[50]; // output
    TH1D * hUncerFromToyMeas = (TH1D*) hMeas->Clone("Error from Meas Toy");
    
    //TRandom3* RandGen = new TRandom3();
    RandGen->SetSeed(0); // set random seed
    double xiSum[50];
    double xiSquareSum[50];
    for(int i = 0; i<50; i++){
        xiSum[i] = 0;
        xiSquareSum[i] = 0;
    }
    int numBinXFinal = response->GetNbinsMeasured();
    
    //--- test distribution of unfolded hist
    RooUnfoldBayes unfoldCTest(response, hMeas, kterm);
    unfoldCTest.SetVerbose(0);
    TH1D *hUnfoldedCTest = (TH1D*) unfoldCTest.Hreco();
    double meanC = hUnfoldedCTest->GetBinContent(2);
    double sigmaC = hUnfoldedCTest->GetBinError(2);
    TH1D *distriUnfoldedMeas = new TH1D ("ErrTest for Smeared Meas", "ErrTest for Smeared Meas", 20, -3.0, 3.0);
    //-----------
    
    
    for (int nthToy = 0; nthToy < myNToys; nthToy ++){
        
        // creat new Meas
        TH1D *hMeasToy = (TH1D*) hMeas->Clone();
        TH1D *hNewMeasToy = (TH1D*) hMeas->Clone();
        
        //int numBinUsed(0);
        //int numBinX = response->GetNbinsMeasured();
        int binthStart = 1;
        int binthEnd = numBinXFinal;
        //numBinUsed = numBinXFinal;
        if (response->UseOverflowStatus()){
            binthStart = 0;
            binthEnd = numBinXFinal +1;
            //numBinUsed = numBinXFinal +2 ;
        }
        
        for (int i = binthStart; i <= binthEnd; i++) {
            
            double errMeasTemp (0);
            double contentMeasTemp (0);
            errMeasTemp = hMeasToy->GetBinError(i);
            contentMeasTemp = hMeasToy->GetBinContent(i);
            if (errMeasTemp > 0.0){
                contentMeasTemp += RandGen->Gaus(0,errMeasTemp);
                if (contentMeasTemp < 0){
                    contentMeasTemp = 0;
                }
            }
            hNewMeasToy->SetBinContent(i, contentMeasTemp);
            
        }
        
        
        // creat new RooUnfold object and Unfolded histo
        RooUnfoldBayes unfoldToyMeas(response, hNewMeasToy, kterm);
        unfoldToyMeas.SetVerbose(0);
        TH1D *hUnfoldedToyMeas = NULL;
        hUnfoldedToyMeas = (TH1D*) unfoldToyMeas.Hreco();
        
        // next calculate the uncertainty
        //cout << endl;
        for (int i = 1; i <= numBinXFinal; i++){
            double xi = hUnfoldedToyMeas->GetBinContent(i);
            //cout << " k: " << nthToy << " binth : " << i << " xi: " << xi << endl;
            if (xi < 0){
                cout << " unfolded value less than Zero " << nthToy << " " << i << endl;
                xi = 0;
            }
            xiSum[i] += xi;
            xiSquareSum[i] += xi * xi;
            
            //--- test distribution of unfolded hist ( second bin only [exclusive 1 jet] )
            if (i == 2) {
                double distriVal = (xi - meanC)/sigmaC;
                distriUnfoldedMeas->Fill(distriVal);
            }
            //---------------
            
        }
        
        //cout << " numBinUsedInReap: " << numBinUsed << endl;
    }
    
    
    //--- test distribution of unfolded hist
    //TCanvas* canvasErrTestMeas = new TCanvas("canvasErrTest for Smeared Meas", "canvasErrTest for Smeared Meas");
    //canvasErrTestMeas->cd();
    //distriUnfoldedMeas->Draw();
    //canvasErrTestMeas->Print("PNGFiles/TestToyErrorMeas.pdf");
    //-------------
    
    // set up output
    for (int i = 1; i <= numBinXFinal; i++){
        uncerFromToyMeas[i] = sqrt( (xiSquareSum[i]/myNToys) - (pow( xiSum[i]/myNToys, 2)) );
        hUncerFromToyMeas->SetBinContent(i, uncerFromToyMeas[i]);
        hUncerFromToyMeas->SetBinError(i, 0.0);
        
        cout << " sumSquare: " << xiSquareSum[i]/myNToys << " squareOfSum: " << pow( xiSum[i]/myNToys, 2) <<  " diff: " << (xiSquareSum[i]/myNToys) - (pow( xiSum[i]/myNToys, 2)) << endl;
    }
    
    cout << " sigmaC:   " << sigmaC << endl;
    cout << " sigmaToy: " << uncerFromToyMeas[2] << endl;
    
    return hUncerFromToyMeas;
    
    //--- End Toy experiment on hMeas ------
    
}



