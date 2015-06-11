// History
//---- 2015_05_24
// Fix the function GetSysErrorTable()
//
//---- 2015_05_16
// BH is done in myFinalUnfoldMeanNJets.cc so that the function of PlottingEditMeanNJets is only to do "plotting". However the part of code that do BH is included here but is comment out. This is useful when we want to test only BH procedure since running myFinalUnfoldMeanNJets.cc is slow.
// For pre approval I do not want to include BH so I comment out lines related to "genBhatALL" and  "gen2" in the plotting section. If I want to include then I would uncomment these lines.
//---- 2015_05_10
// MES and MER is now included


#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include "getFilesAndHistograms.h"

//-- Setting global variables -------------
#include "fileNames.h"
#include "variablesOfInterestVarWidth.h"

using namespace std;
//--------------------------------


//--- Declaring functions -------------
void PlottingEditMeanNJets(int start = 0, int end = -1);
void FuncPlot(string variable = "ZNGoodJets_Zexc", bool log = 1, bool decrease = 1);
void GetSysErrorTable (string variable, string outputTableName, TH1D* dataCenBackUp = NULL, TH1D* hDiffer[] = {NULL}, TH1D* hStatError = NULL, TH1D* hTotalError = NULL, int nGroup = 10);
void plotSystematicBreakdown (string outputDirectory, string variable, TH1D* dataCentral = NULL, TH1D* hSyst[] = NULL, TH1D* hErrorFromToyResponse = NULL);
//--------------------------------


//--- More global variables -------------

double Luminosity(19.244 * 1000); // for SMu
bool isMuon = 1;
string leptonFlavor = "SMu";
string unfAlg = "Bayes";
bool doVarWidth = true ;
//--------------------------------
const int ZJetsFillStyle = 1001;



void PlottingEditMeanNJets(int start, int end)
{
    //NVAROFINTEREST = 41;
    //for (int i(0); i < 1 /*NVAROFINTEREST*/; i++) {
    if (end == -1) end = start + 1;
    for (int i(start); i < end /*NVAROFINTERESTZJETS*/; i++) {
        FuncPlot(VAROFINTERESTWJETS[i].name, VAROFINTERESTWJETS[i].log, VAROFINTERESTWJETS[i].decrease);
    }
    
}

void FuncPlot(string variable, bool logZ, bool decrease)
{
    cout << "processing variable: " << variable << endl;
    string energy = getEnergy();
    if ( !isMuon ) leptonFlavor = "DE";
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    
    // set file name for openning
    string fileName = "PNGFiles/FinalUnfoldMeanNJets/" + leptonFlavor + "_" + energy +   "_unfolded_" + variable + "_histograms_" + unfAlg ;
    if (doVarWidth) fileName += "_VarWidth";
    fileName += ".root";
    cout << " opening : " << fileName <<endl;
    TFile *f  = new TFile(fileName.c_str());
    
    // get cental and gen
    
    TH1D *dataCentral       = (TH1D*) f->Get("Central"); // unfolded central (data -BG)
    //TH1D *data              = (TH1D*) f->Get("Data"); // reco data
    TH1D *data              = (TH1D*) dataCentral->Clone();
    TH1D *genMad            = (TH1D*) f->Get("genMad");
    TH1D *dataCentralBackUp = (TH1D*) dataCentral->Clone();
    
    // get unfolded histogram for systematics
    int nSyst(15);
    TH1D *hSyst[15];
    hSyst[0] = (TH1D*) f->Get("JESup");
    hSyst[1] = (TH1D*) f->Get("JESdown");
    hSyst[2] = (TH1D*) f->Get("PUup");
    hSyst[3] = (TH1D*) f->Get("PUdown");
    hSyst[4] = (TH1D*) f->Get("XSECup");
    hSyst[5] = (TH1D*) f->Get("XSECdown");
    hSyst[6] = (TH1D*) f->Get("JERup");
    hSyst[7] = (TH1D*) f->Get("BtagEFFup");
    hSyst[8] = (TH1D*) f->Get("BtagEFFdown");
    hSyst[9] = (TH1D*) f->Get("WBup");
    hSyst[10] = (TH1D*) f->Get("TTBARup");
    hSyst[11] = (TH1D*) f->Get("RESPup");
    hSyst[12] = (TH1D*) f->Get("MESup");
    hSyst[13] = (TH1D*) f->Get("MESdown");
    hSyst[14] = (TH1D*) f->Get("MER");
    
    //--- get uncer from Toy Resp (MC)
    TH1D *hErrorFromToyResponse = (TH1D*) f->Get("MC");
    
    //--- Get BlackHat ------
    TH1D *genBhatALL = (TH1D*) f->Get("genBhat");
    
    /*
     // if you want to do BH here then uncomment this part
    TFile *fBhat[3];
    fBhat[0] = new TFile("PNGFiles/BlackHat/W1j_all.root", "READ");
    fBhat[1] = new TFile("PNGFiles/BlackHat/W2j_all.root", "READ");
    fBhat[2] = new TFile("PNGFiles/BlackHat/W3j_all.root", "READ");
    cout << " Opening: " << "W1j_all.root" << "   --->   Opened ? " << fBhat[0]->IsOpen() << endl;
    cout << " Opening: " << "W2j_all.root" << "   --->   Opened ? " << fBhat[1]->IsOpen() << endl;
    cout << " Opening: " << "W3j_all.root" << "   --->   Opened ? " << fBhat[2]->IsOpen() << endl;
    
    //-- matching variable name
    string varBlackHatName, varBlackHatName2;
    int useBHFile(0);
    if (variable == "MeanNJetsHT_Zinc1jet") {
        varBlackHatName = "totNJhtjet1";
        varBlackHatName2 = "htjet1";
        useBHFile = 0;
    }
    else if (variable == "MeanNJetsHT_Zinc2jet") {
        varBlackHatName = "totNJhtjet2";
        varBlackHatName2 = "htjet2";
        useBHFile = 1;
    }
    else if (variable == "MeanNJetsdRapidity_Zinc2jet") {
        varBlackHatName = "totNJdyj1j2";
        varBlackHatName2 = "dyj1j2jet2";
        useBHFile = 1;
    }
    else if (variable == "MeanNJetsdRapidityFB_Zinc2jet") {
        varBlackHatName = "totNJdyjFjB";
        varBlackHatName2 = "dyjFjBjet2";
        useBHFile = 1;
    }
    
    TH1D *genBhat[3], *genBhatDeno[3];
    genBhat[0] = (TH1D*) fBhat[0]->Get(varBlackHatName.c_str()); genBhat[0]->Scale(0.001);
    genBhat[1] = (TH1D*) fBhat[1]->Get(varBlackHatName.c_str()); genBhat[1]->Scale(0.001);
    genBhat[2] = (TH1D*) fBhat[2]->Get(varBlackHatName.c_str()); genBhat[2]->Scale(0.001);
    
    genBhatDeno[0] = (TH1D*) fBhat[0]->Get(varBlackHatName2.c_str()); genBhatDeno[0]->Scale(0.001);
    genBhatDeno[1] = (TH1D*) fBhat[1]->Get(varBlackHatName2.c_str()); genBhatDeno[1]->Scale(0.001);
    genBhatDeno[2] = (TH1D*) fBhat[2]->Get(varBlackHatName2.c_str()); genBhatDeno[2]->Scale(0.001);
    
    TH1D *genBhatMeanJ = NULL;
    genBhatMeanJ = (TH1D*) genBhat[useBHFile]->Clone();
    genBhatMeanJ->Divide(genBhatDeno[useBHFile]);
    
    TH1D *genBhatALL = (TH1D*) genBhatMeanJ->Clone();
    */
    //--- End Get BlackHat ------
    
    
    //--- Get Sherpa ------
    TH1D *genShe = (TH1D*) f->Get("genShe");
    
    // print out integral value for testing
    const int nBins(dataCentral->GetNbinsX());
    cout << " central integral:  " << variable << "   " << dataCentral->Integral(1, nBins) << " with under/over:  " << dataCentral->Integral(0, nBins+1) << endl;
    for (int syst(0); syst < nSyst; syst++) {
        double tempScale = hSyst[syst]->Integral(1, nBins) ;
        cout << " original integral for " << syst << " : " << tempScale << endl;
    }
    
    // no need to do cross section or normalized
    bool doXSec(1);
    bool doNormalize(0);
    if (doXSec) doNormalize = false;
    if (doNormalize) doXSec = false;
    
    //--- systematic error computation --------------------------
    cout << "go to uncertainty computation " << endl;
    double xCoor[nBins], yCoor[nBins], xErr[nBins], ySystDown[nBins], ySystUp[nBins];
    for (int bin(1); bin <= nBins; bin++) {
        
        double centralValue (dataCentral->GetBinContent(bin));
        double totalStatistics (dataCentral->GetBinError(bin));
        double totalSystematicsUp (0.);
        double totalSystematicsDown (0.);
        
        // set effSF = 0. So it has no contribution.
        double effSF = 0.00;
        totalSystematicsUp += pow( pow((1.+effSF), 2) - 1 , 2);
        totalSystematicsDown += pow( pow((1.+effSF), 2) - 1 , 2);
        
        for(int i = 0; i< 10; i++) {
            // contribution from JES PU and Xsec
            if (i < 3) {
                double diffFromUp(hSyst[2*i]->GetBinContent(bin) - centralValue);
                double diffFromDown(hSyst[(2*i)+1]->GetBinContent(bin) - centralValue);
                double SysDiffUPTemp(0);
                double SysDiffDownTemp(0);
                
                if (diffFromUp >= diffFromDown) {
                    SysDiffUPTemp = diffFromUp;
                    SysDiffDownTemp = diffFromDown;
                }
                else {
                    SysDiffUPTemp = diffFromDown;
                    SysDiffDownTemp = diffFromUp;
                }
                
                if (SysDiffUPTemp >= 0) {
                    totalSystematicsUp += SysDiffUPTemp * SysDiffUPTemp;
                }
                if (SysDiffDownTemp <= 0) {
                    totalSystematicsDown += SysDiffDownTemp * SysDiffDownTemp;
                }
            }
            // contribution from JER
            if (i == 3) {
                double difference(hSyst[2*i]->GetBinContent(bin) - centralValue);
                totalSystematicsUp += difference * difference;
                totalSystematicsDown += difference * difference;
            }
            // contribution from BtagEff
            if (i == 4){
                double diffFromUp(hSyst[7]->GetBinContent(bin) - centralValue);
                double diffFromDown(hSyst[8]->GetBinContent(bin) - centralValue);
                double SysDiffUPTemp(0);
                double SysDiffDownTemp(0);
                
                if (diffFromUp >= diffFromDown) {
                    SysDiffUPTemp = diffFromUp;
                    SysDiffDownTemp = diffFromDown;
                }
                else {
                    SysDiffUPTemp = diffFromDown;
                    SysDiffDownTemp = diffFromUp;
                }
                
                if (SysDiffUPTemp >= 0) {
                    totalSystematicsUp += SysDiffUPTemp * SysDiffUPTemp;
                }
                if (SysDiffDownTemp <= 0) {
                    totalSystematicsDown += SysDiffDownTemp * SysDiffDownTemp;
                }
            }
            // contribution from Wb
            if (i == 5) {
                double difference(hSyst[9]->GetBinContent(bin) - centralValue);
                totalSystematicsUp += difference * difference;
                totalSystematicsDown += difference * difference;
            }
            // contribution from ttbar
            if (i == 6) {
                double difference(hSyst[10]->GetBinContent(bin) - centralValue);
                totalSystematicsUp += difference * difference;
                totalSystematicsDown += difference * difference;
            }
            // contribution from Resp
            if (i == 7) {
                double difference(hSyst[11]->GetBinContent(bin) - centralValue);
                totalSystematicsUp += difference * difference;
                totalSystematicsDown += difference * difference;
            }
            // contribution from MES
            if (i == 8) {
                double diffFromUp(hSyst[12]->GetBinContent(bin) - centralValue);
                double diffFromDown(hSyst[13]->GetBinContent(bin) - centralValue);
                double SysDiffUPTemp(0);
                double SysDiffDownTemp(0);
                
                if (diffFromUp >= diffFromDown) {
                    SysDiffUPTemp = diffFromUp;
                    SysDiffDownTemp = diffFromDown;
                }
                else {
                    SysDiffUPTemp = diffFromDown;
                    SysDiffDownTemp = diffFromUp;
                }
                
                if (SysDiffUPTemp >= 0) {
                    totalSystematicsUp += SysDiffUPTemp * SysDiffUPTemp;
                }
                if (SysDiffDownTemp <= 0) {
                    totalSystematicsDown += SysDiffDownTemp * SysDiffDownTemp;
                }
            }
            // contribution from MER
            if (i == 9) {
                double difference(hSyst[14]->GetBinContent(bin) - centralValue);
                totalSystematicsUp += difference * difference;
                totalSystematicsDown += difference * difference;
            }
        }
        
        // contribution from integrated lumi = 2.6%
        totalSystematicsUp += pow( ((2.6/100.0) * centralValue), 2 );
        totalSystematicsDown += pow( ((2.6/100.0) * centralValue), 2 );
        
        
        // contribution from Toy Resp (MC)
            // (what written in the file is the absolute error value in "counts", so we need the dataCentralBackUp for normalization to xsec) => this is not true for meanNjets >> I do not need this step anymore >> FIXME
        hErrorFromToyResponse->SetBinContent( bin, ((hErrorFromToyResponse->GetBinContent(bin))/(dataCentralBackUp->GetBinContent(bin))) * centralValue );
        totalSystematicsUp += pow( hErrorFromToyResponse->GetBinContent(bin)  ,  2);
        totalSystematicsDown += pow( hErrorFromToyResponse->GetBinContent(bin)  ,  2);
        
        
        
        //--- Set up for systematic plot on pad 1
            // set x-coordinate and error in x
        xCoor[bin-1]    = dataCentral->GetBinCenter(bin);
        xErr[bin-1]     = 0.5 * dataCentral->GetBinWidth(bin);
            // set y-coordinate
        yCoor[bin-1]    = centralValue;
            // set systematics errors in y => combine both statistical and systematics
        ySystUp[bin-1]     = sqrt(totalStatistics * totalStatistics + totalSystematicsUp );
        ySystDown[bin-1]   = sqrt(totalStatistics * totalStatistics + totalSystematicsDown );
    }
    cout << " finished uncertainty computation " << endl;
    
    //--- Set up for systematic plot on pad 1
    TGraphAsymmErrors *grCentralSyst = new TGraphAsymmErrors(nBins, xCoor, yCoor, xErr, xErr, ySystDown, ySystUp);
    
    
    //--- plotting ---
    TH1D *hisUnfolded = (TH1D*) dataCentral->Clone();
    
    TCanvas *can = new TCanvas(variable.c_str(), variable.c_str(), 600, 800);
    can->cd();
    //--- Set pad1.
    TPad *pad1 = new TPad("pad1","pad1", 0.01, 0.55, 0.99, 0.99);
    pad1->SetTopMargin(0.11);
    pad1->SetBottomMargin(0);
    pad1->SetRightMargin(0.04);
    pad1->SetLeftMargin(0.11);
    pad1->SetTicks();
    pad1->Draw();
    pad1->cd();
    //pad1->SetLogy();
    
    // set y axis title
    string temp = "<N_{jets}>";
    cout << temp << endl;
    //-------
    //Set x axis title
    string tempXTitle;
    double RyMin, RyMax;
    if (variable == "MeanNJetsHT_Zinc1jet"){
        tempXTitle = "H_{T}(jets) [GeV]";
        RyMin = 0.9;
        RyMax = 4.6;
    }
    else if (variable == "MeanNJetsHT_Zinc2jet"){
        tempXTitle = "H_{T}(jets) [GeV]";
        RyMin = 1.9;
        RyMax = 4.2;
    }
    else if (variable == "MeanNJetsdRapidity_Zinc2jet"){
        tempXTitle = "#Delta y(j_{1}j_{2})";
        RyMin = 1.9;
        RyMax = 2.7;
    }
    else if (variable == "MeanNJetsdRapidityFB_Zinc2jet"){
        tempXTitle = "#Delta y(j_{F}j_{B})";
        RyMin = 1.9;
        RyMax = 3.2;
    }
    
    
    //double MineYMax = genMad->GetMaximum();
    //double MineYMin = genMad->GetMinimum();
    //if (hisUnfolded->GetMaximum() > MineYMax) MineYMax = hisUnfolded->GetMaximum();
    
    hisUnfolded->SetTitle("");
    if (variable.find("ZNGood") != string::npos) {
        hisUnfolded->GetXaxis()->SetRange(2, 8);
    }
    hisUnfolded->GetYaxis()->SetTitle(temp.c_str());
    hisUnfolded->GetYaxis()->SetTitleSize(0.04);
    hisUnfolded->GetYaxis()->SetTitleOffset(1.45);
    hisUnfolded->GetYaxis()->SetTitleFont(42);
    hisUnfolded->GetYaxis()->SetRangeUser(RyMin, RyMax);
    //hisUnfolded->SetMaximum(MineYMax*1.2);
    //hisUnfolded->SetMinimum(MineYMin*0.8);
    cout << "Stop after line " << __LINE__ << endl;
    
    hisUnfolded->SetFillStyle(3354);
    hisUnfolded->SetFillColor(12);
    hisUnfolded->SetLineColor(kBlack);
    hisUnfolded->SetLineWidth(2);
    hisUnfolded->SetMarkerColor(kBlack);
    hisUnfolded->SetMarkerStyle(20);
    hisUnfolded->SetMarkerSize(0.8);
    hisUnfolded->Draw("E1");
    cout << "Stop after line " << __LINE__ << endl;
    
    grCentralSyst->SetFillStyle(3354);
    grCentralSyst->SetFillColor(12);
    grCentralSyst->Draw("2");
    
    genMad->SetFillStyle(ZJetsFillStyle);
    genMad->SetFillColor(kBlue-10);
    genMad->SetLineColor(kBlue);
    genMad->SetLineWidth(2);
    genMad->SetMarkerColor(kBlue);
    genMad->SetMarkerStyle(24);
    genMad->Draw("E2 same");
    genMad->Draw("E same");
    cout << "Stop after line " << __LINE__ << endl;
    
    genBhatALL->SetFillStyle(ZJetsFillStyle);
    genBhatALL->SetFillColor(kOrange-2);
    genBhatALL->SetLineColor(kOrange+10);
    genBhatALL->SetLineWidth(2);
    genBhatALL->SetMarkerColor(kOrange+10);
    genBhatALL->SetMarkerStyle(25);
    //genBhatALL->Draw("E2 same");
    //genBhatALL->Draw("E same");
    cout << "Stop after line " << __LINE__ << endl;
    
    genShe->SetFillStyle(ZJetsFillStyle);
    genShe->SetFillColor(kGreen-8);
    genShe->SetLineColor(kGreen+3);
    genShe->SetLineWidth(2);
    genShe->SetMarkerColor(kGreen+3);
    genShe->SetMarkerStyle(26);
    genShe->Draw("E2 same");
    genShe->Draw("E same");
    cout << "Stop after line " << __LINE__ << endl;
    
    pad1->Draw();
    //--- TLegend ---
    TLegend *legend = new TLegend(0.47, 0.74, 0.99, 0.98);
    legend->SetX1(0.44);
    legend->SetY1(0.77);
    legend->SetX2(0.95);
    legend->SetY2(0.98);
    legend->SetTextSize(.034);
    legend->SetFillColor(0);
    legend->SetFillStyle(1001);
    legend->SetBorderSize(1);
    legend->AddEntry(hisUnfolded, "Data", "PLEF");
    legend->AddEntry(genMad, "MG+Py6 at NNLO Stat. Error", "plef");
    legend->AddEntry(genBhatALL, "BlackHat at NLO Stat. Error", "plef");
    legend->AddEntry(genShe, "Sherpa2 Stat. Error", "plef");
    legend->Draw();
    //------------------
    cout << "Stop after line " << __LINE__ << endl;
    
    //--- calculate ratio >> pad2, pad3
    
    // for data statistical uncer
    TH1D  *data1 = (TH1D*) hisUnfolded->Clone();
    for(int i = 1; i <= hisUnfolded->GetNbinsX(); i++) {
        data1->SetBinContent(i, 1.0);
        if (hisUnfolded->GetBinContent(i) > 0.) {
            data1->SetBinError(i, hisUnfolded->GetBinError(i)/hisUnfolded->GetBinContent(i));
        }
        else data1->SetBinError(i, 0.);
        
    }
    TH1D  *data2 = (TH1D*) data1->Clone();
    TH1D  *data3 = (TH1D*) data1->Clone();
    
    // for gen statistical uncer
    TH1D  *gen1 = (TH1D*) genMad->Clone();
    for(int i = 1; i <= hisUnfolded->GetNbinsX(); i++) {
        if (hisUnfolded->GetBinContent(i) > 0.) {
            gen1->SetBinContent(i, genMad->GetBinContent(i)/hisUnfolded->GetBinContent(i));
            gen1->SetBinError(i, genMad->GetBinError(i)/hisUnfolded->GetBinContent(i));
        }
        else {
            gen1->SetBinContent(i, 0.);
            gen1->SetBinError(i, 0.);
        }
    }
    TH1D  *gen2 = (TH1D*) hisUnfolded->Clone();
    for(int i = 1; i <= hisUnfolded->GetNbinsX(); i++) {
        if (hisUnfolded->GetBinContent(i) > 0.) {
            gen2->SetBinContent(i, genBhatALL->GetBinContent(i)/hisUnfolded->GetBinContent(i));
            gen2->SetBinError(i, genBhatALL->GetBinError(i)/hisUnfolded->GetBinContent(i));
        }
        else {
            gen2->SetBinContent(i, 0.);
            gen2->SetBinError(i, 0.);
        }
    }
    TH1D  *gen3 = (TH1D*) hisUnfolded->Clone();
    for(int i = 1; i <= hisUnfolded->GetNbinsX(); i++) {
        if (hisUnfolded->GetBinContent(i) > 0.) {
            gen3->SetBinContent(i, genShe->GetBinContent(i)/hisUnfolded->GetBinContent(i));
            gen3->SetBinError(i, genShe->GetBinError(i)/hisUnfolded->GetBinContent(i));
        }
        else {
            gen3->SetBinContent(i, 0.);
            gen3->SetBinError(i, 0.);
        }
        //cout << " sherpa " << gen3->GetBinContent(i) << endl;
    }
    
    //--- Set up for systematic plot on pad 2, 3
    double xCoorCenRatio[nBins], yCoorCenRatio[nBins], xErrCenRatio[nBins], ySystDownCenRatio[nBins], ySystUpCenRatio[nBins];
    for (int bin(1); bin <= nBins; bin++) {
        // set x-coordinate and error in x
        xCoorCenRatio[bin-1] = hisUnfolded->GetBinCenter(bin);
        xErrCenRatio[bin-1]  = 0.5 * hisUnfolded->GetBinWidth(bin);
        // set y-coordinate
        yCoorCenRatio[bin-1] = 1.;
        // set systematics errors in y => combine both statistical and systematics
        if (hisUnfolded->GetBinContent(bin) > 0) {
            ySystUpCenRatio[bin-1]     = ySystUp[bin-1]/hisUnfolded->GetBinContent(bin);
            ySystDownCenRatio[bin-1]   = ySystDown[bin-1]/hisUnfolded->GetBinContent(bin);
        }
        else {
            ySystUpCenRatio[bin-1] = 0.;
            ySystDownCenRatio[bin-1] = 0.;
        }
        
    }
    TGraphAsymmErrors *grCentralSystRatio = new TGraphAsymmErrors(nBins, xCoorCenRatio, yCoorCenRatio, xErrCenRatio, xErrCenRatio, ySystDownCenRatio, ySystUpCenRatio);
    

    //--- NEW TPAD => pad2 ---//
    can->cd();
    TPad *pad2 = new TPad("pad2","pad2", 0.01, 0.39, 0.99, 0.55);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0);
    pad2->SetRightMargin(0.04);
    pad2->SetLeftMargin(0.11);
    pad2->SetGridy();
    pad2->SetTicks();
    pad2->Draw();
    pad2->cd();
    
    data1->SetTitle("");
    data1->GetYaxis()->SetRangeUser(0.41,1.59);
    data1->GetYaxis()->SetTitle("MC / Data");
    data1->GetYaxis()->SetTitleSize(0.135);
    data1->GetYaxis()->SetTitleOffset(0.44);
    data1->GetYaxis()->CenterTitle();
    data1->GetYaxis()->SetLabelSize(0.12);
    data1->GetYaxis()->SetLabelOffset(0.01);
    
    data1->GetXaxis()->SetTickLength(0.06);
    
    data1->SetLineColor(kBlack);
    data1->SetLineWidth(1);
    data1->SetMarkerColor();
    data1->SetMarkerStyle();
    data1->Draw("E1");
    
    grCentralSystRatio->SetFillStyle(3354);
    grCentralSystRatio->SetFillColor(12);
    grCentralSystRatio->Draw("2");
    
    gen1->SetFillStyle(ZJetsFillStyle);
    gen1->SetFillColor(kBlue-10);
    gen1->SetLineColor(kBlue);
    gen1->SetLineWidth(2);
    gen1->SetMarkerColor(kBlue);
    gen1->SetMarkerStyle(24);
    gen1->Draw("E2 same");
    gen1->Draw("E1 same");
    
    TLatex *latexLabel = new TLatex();
    latexLabel->SetTextSize(0.15);
    latexLabel->SetTextFont(42);
    latexLabel->SetLineWidth(2);
    latexLabel->SetTextColor(kAzure-5);
    latexLabel->SetNDC();
    latexLabel->DrawLatex(0.15,0.09,"Madgraph");
    
    pad2->Draw();
    
    
    //--- NEW TPAD => pad3 ---//
    can->cd();
    TPad *pad3 = new TPad("pad3","pad3", 0.01, 0.23, 0.99, 0.39);
    pad3->SetTopMargin(0);
    pad3->SetBottomMargin(0);
    pad3->SetRightMargin(0.04);
    pad3->SetLeftMargin(0.11);
    pad3->SetGridy();
    pad3->SetTicks();
    pad3->Draw();
    pad3->cd();
    
    data2->SetTitle("");
    data2->GetYaxis()->SetRangeUser(0.41,1.59);
    data2->GetYaxis()->SetTitle("MC / Data");
    data2->GetYaxis()->SetTitleSize(0.135);
    data2->GetYaxis()->SetTitleOffset(0.44);
    data2->GetYaxis()->CenterTitle();
    data2->GetYaxis()->SetLabelSize(0.12);
    data2->GetYaxis()->SetLabelOffset(0.01);
    
    data2->GetXaxis()->SetTickLength(0.06);
    
    data2->SetLineColor(kBlack);
    data2->SetLineWidth(1);
    data2->SetMarkerColor();
    data2->SetMarkerStyle();
    data2->Draw("E1");
    
    grCentralSystRatio->SetFillStyle(3354);
    grCentralSystRatio->SetFillColor(12);
    grCentralSystRatio->Draw("2");
    
    gen2->SetFillStyle(ZJetsFillStyle);
    gen2->SetFillColor(kOrange-2);
    gen2->SetLineColor(kOrange+10);
    gen2->SetLineWidth(2);
    gen2->SetMarkerColor(kOrange+10);
    gen2->SetMarkerStyle(25);
    //gen2->Draw("E2 same");
    //gen2->Draw("E1 same");
    
    // Text BlackHat
    TLatex *latexLabel1 = new TLatex();
    latexLabel1->SetTextSize(0.15);
    latexLabel1->SetTextFont(42);
    latexLabel1->SetLineWidth(2);
    latexLabel1->SetTextColor(kOrange+2);
    latexLabel1->SetNDC();
    latexLabel1->DrawLatex(0.15,0.09,"BlackHat");
    
    pad3->Draw();
    
    //--- NEW TPAD => pad4 ---//
    can->cd();
    TPad *pad4 = new TPad("pad3","pad3", 0.01, 0.01, 0.99, 0.23);
    pad4->SetTopMargin(0);
    pad4->SetBottomMargin(0.25);
    pad4->SetRightMargin(0.04);
    pad4->SetLeftMargin(0.11);
    pad4->SetGridy();
    pad4->SetTicks();
    pad4->Draw();
    pad4->cd();
    
    data3->SetTitle("");
    data3->GetYaxis()->SetRangeUser(0.41,1.59);
    data3->GetYaxis()->SetTitle("MC / Data");
    data3->GetYaxis()->SetTitleSize(0.1);
    data3->GetYaxis()->SetTitleOffset(0.6);
    data3->GetYaxis()->CenterTitle();
    data3->GetYaxis()->SetLabelSize(0.085);
    data3->GetYaxis()->SetLabelOffset(0.012);
    
    data3->GetXaxis()->SetTitle(tempXTitle.c_str());
    data3->GetXaxis()->SetTitleSize(0.1);
    data3->GetXaxis()->SetTitleOffset(1.05);
    data3->GetXaxis()->SetLabelSize(0.09);
    data3->GetXaxis()->SetLabelOffset(0.014);
    data3->GetXaxis()->SetTickLength(0.06);
    
    data3->SetLineColor(kBlack);
    data3->SetLineWidth(1);
    data3->SetMarkerColor();
    data3->SetMarkerStyle();
    data3->Draw("E1");
    
    grCentralSystRatio->SetFillStyle(3354);
    grCentralSystRatio->SetFillColor(12);
    grCentralSystRatio->Draw("2");
    
    gen3->SetFillStyle(ZJetsFillStyle);
    gen3->SetFillColor(kGreen-8);
    gen3->SetLineColor(kGreen+3);
    gen3->SetLineWidth(2);
    gen3->SetMarkerColor(kGreen+3);
    gen3->SetMarkerStyle(26);
    gen3->Draw("E2 same");
    gen3->Draw("E1 same");
    
    // Text Sherpa
    TLatex *latexLabel2 = new TLatex();
    latexLabel2->SetTextSize(0.11);
    latexLabel2->SetTextFont(42);
    latexLabel2->SetLineWidth(2);
    latexLabel2->SetTextColor(kGreen+3);
    latexLabel2->SetNDC();
    latexLabel2->DrawLatex(0.15,0.35,"Sherpa");
    
    pad4->Draw();
    pad4->cd();
    can->cd();
    
    //// Set output .pdf file in the direcctory "PNGFiles/FinalUnfold/"
    string outputDirectory = "PNGFiles/FinalUnfoldMeanNJets/";
    string outputFileNamePNG = outputDirectory + leptonFlavor + "_" + energy + "_RunPlotting_";
    
    if (doXSec) outputFileNamePNG += "XSec_";
    else if (doNormalize) outputFileNamePNG += "Normalized_";
    
    outputFileNamePNG += variable + "_" + unfAlg ;
    if (doVarWidth)  outputFileNamePNG += "_VarWidth";
    outputFileNamePNG += ".pdf";
    can->Print(outputFileNamePNG.c_str());
    
    
    //--- Generate new plot to see Contribution : all input histo are already scaled to the final unit ---
    plotSystematicBreakdown(outputDirectory, variable, dataCentral, hSyst, hErrorFromToyResponse);

}


void plotSystematicBreakdown (string outputDirectory, string variable, TH1D* dataCentral, TH1D* hSyst[], TH1D* hErrorFromToyResponse)
{
    TH1D *dataCenBackUp = (TH1D*) dataCentral->Clone();
    
    const int nGroup = 12;
    const int nBins(dataCentral->GetNbinsX());
    
    int group = nGroup -1;
    if (variable != "ZNGoodJets_Zexc") {
        group = nGroup;
    }
    
    TH1D *hUpDiffer[nGroup] = {NULL}, *hDownDiffer[nGroup] = {NULL}, *hStatError = {NULL}, *hTotalErrorUp = {NULL}, *hTotalErrorDown = {NULL} ; //nSyst(12);
    
    // set up histogram
    hStatError =  (TH1D*) dataCentral->Clone();
    hStatError->SetTitle("");
    hStatError->GetYaxis()->SetTitle("");
    
    hTotalErrorUp =  (TH1D*) dataCentral->Clone();
    hTotalErrorUp->SetTitle("");
    hTotalErrorUp->GetYaxis()->SetTitle("");
    
    hTotalErrorDown =  (TH1D*) dataCentral->Clone();
    hTotalErrorDown->SetTitle("");
    hTotalErrorDown->GetYaxis()->SetTitle("");
    
    // Initialize hUpDiffer hDownDiffer
    for (int syst(0); syst < group; syst++){
        hUpDiffer[syst] = (TH1D*) dataCentral->Clone();
        hDownDiffer[syst] = (TH1D*) dataCentral->Clone();
        
        hUpDiffer[syst]->SetTitle("");
        hUpDiffer[syst]->GetYaxis()->SetTitle("");
        hDownDiffer[syst]->SetTitle("");
        hDownDiffer[syst]->GetYaxis()->SetTitle("");
    }
    
    // Assign Value
    for (int bin(1); bin <= nBins; bin++){
        for(int i = 0; i < group; i++){
            
            hUpDiffer[i]->SetBinContent(bin, 0.);
            hDownDiffer[i]->SetBinContent(bin, 0.);
            
            if( dataCentral->GetBinContent(bin) > 0 ){
                // contribution from JES PU and Xsec
                if (i < 3){
                    double SysDiffUPTemp(0);
                    double SysDiffDownTemp(0);
                    
                    double diffFromUp( hSyst[2*i]->GetBinContent(bin) - (dataCentral->GetBinContent(bin)) );
                    double diffFromDown( hSyst[(2*i)+1]->GetBinContent(bin) - (dataCentral->GetBinContent(bin)) );
                    
                    
                    if (diffFromUp >= diffFromDown) {
                        SysDiffUPTemp = diffFromUp;
                        SysDiffDownTemp = diffFromDown;
                    }
                    else{
                        SysDiffUPTemp = diffFromDown;
                        SysDiffDownTemp = diffFromUp;
                    }
                    
                    if (SysDiffUPTemp >= 0){
                        hUpDiffer[i]->SetBinContent( bin, SysDiffUPTemp * (100 / dataCentral->GetBinContent(bin)) ) ;
                    }
                    if (SysDiffDownTemp <= 0){
                        hDownDiffer[i]->SetBinContent( bin, (fabs(SysDiffDownTemp)) * (100 / dataCentral->GetBinContent(bin)) ) ;
                    }
                }
                // contribution from JER
                if (i == 3){
                    double difference(hSyst[2*i]->GetBinContent(bin) - (dataCentral->GetBinContent(bin)));
                    hUpDiffer[i]  ->SetBinContent( bin, (fabs(difference)) * (100 / dataCentral->GetBinContent(bin)) ) ;
                    hDownDiffer[i]->SetBinContent( bin, (fabs(difference)) * (100 / dataCentral->GetBinContent(bin)) ) ;
                }
                // contribution from MC
                if (i == 4){
                    double uncer(hErrorFromToyResponse->GetBinContent(bin));
                    hUpDiffer[i]  ->SetBinContent( bin, uncer * (100 / dataCentral->GetBinContent(bin)) ) ;
                    hDownDiffer[i]->SetBinContent( bin, uncer * (100 / dataCentral->GetBinContent(bin)) ) ;
                }
                // contribution from ttbar
                if (i == 5){
                    double difference(hSyst[10]->GetBinContent(bin) - (dataCentral->GetBinContent(bin)));
                    hUpDiffer[i]  ->SetBinContent( bin, (fabs(difference)) * (100 / dataCentral->GetBinContent(bin)) ) ;
                    hDownDiffer[i]->SetBinContent( bin, (fabs(difference)) * (100 / dataCentral->GetBinContent(bin)) ) ;
                }
                // contribution from Wb
                if (i == 6){
                    double difference(hSyst[9]->GetBinContent(bin) - (dataCentral->GetBinContent(bin)));
                    hUpDiffer[i]  ->SetBinContent( bin, (fabs(difference)) * (100 / dataCentral->GetBinContent(bin)) ) ;
                    hDownDiffer[i]->SetBinContent( bin, (fabs(difference)) * (100 / dataCentral->GetBinContent(bin)) ) ;
                }
                // contribution from integrated lumi
                if (i == 7){
                    hUpDiffer[i]  ->SetBinContent( bin, 2.6 ) ;
                    hDownDiffer[i]->SetBinContent( bin, 2.6 ) ;
                }
                // contribution from Btag CorrFactor
                if (i == 8){
                    double BdiffFromUp(hSyst[7]->GetBinContent(bin) - (dataCentral->GetBinContent(bin)));
                    double BdiffFromDown(hSyst[8]->GetBinContent(bin) - (dataCentral->GetBinContent(bin)));
                    double BSysDiffUPTemp(0);
                    double BSysDiffDownTemp(0);
                    
                    if (BdiffFromUp >= BdiffFromDown) {
                        BSysDiffUPTemp = BdiffFromUp;
                        BSysDiffDownTemp = BdiffFromDown;
                    }
                    else {
                        BSysDiffUPTemp = BdiffFromDown;
                        BSysDiffDownTemp = BdiffFromUp;
                    }
                    
                    if (BSysDiffUPTemp >= 0) {
                        hUpDiffer[i]->SetBinContent( bin, BSysDiffUPTemp * (100. / dataCentral->GetBinContent(bin)) ) ;
                    }
                    if (BSysDiffDownTemp <= 0) {
                        hDownDiffer[i]->SetBinContent( bin, (fabs(BSysDiffDownTemp)) * (100. / dataCentral->GetBinContent(bin)) ) ;
                    }
                }
                // contribution from MES
                if (i == 9){
                    double BdiffFromUp(hSyst[12]->GetBinContent(bin) - (dataCentral->GetBinContent(bin)));
                    double BdiffFromDown(hSyst[13]->GetBinContent(bin) - (dataCentral->GetBinContent(bin)));
                    double BSysDiffUPTemp(0);
                    double BSysDiffDownTemp(0);
                    
                    if (BdiffFromUp >= BdiffFromDown) {
                        BSysDiffUPTemp = BdiffFromUp;
                        BSysDiffDownTemp = BdiffFromDown;
                    }
                    else {
                        BSysDiffUPTemp = BdiffFromDown;
                        BSysDiffDownTemp = BdiffFromUp;
                    }
                    
                    if (BSysDiffUPTemp >= 0) {
                        hUpDiffer[i]->SetBinContent( bin, BSysDiffUPTemp * (100. / dataCentral->GetBinContent(bin)) ) ;
                    }
                    if (BSysDiffDownTemp <= 0) {
                        hDownDiffer[i]->SetBinContent( bin, (fabs(BSysDiffDownTemp)) * (100. / dataCentral->GetBinContent(bin)) ) ;
                    }
                }
                // contribution from MER
                if (i == 10){
                    double difference(hSyst[14]->GetBinContent(bin) - (dataCentral->GetBinContent(bin)));
                    hUpDiffer[i]  ->SetBinContent( bin, (fabs(difference)) * (100 / dataCentral->GetBinContent(bin)) ) ;
                    hDownDiffer[i]->SetBinContent( bin, (fabs(difference)) * (100 / dataCentral->GetBinContent(bin)) ) ;
                }
                // contribution from Resp
                if (i == nGroup-1){
                    double difference(hSyst[11]->GetBinContent(bin) - (dataCentral->GetBinContent(bin)));
                    hUpDiffer[i]  ->SetBinContent( bin, (fabs(difference)) * (100 / dataCentral->GetBinContent(bin)) ) ;
                    hDownDiffer[i]->SetBinContent( bin, (fabs(difference)) * (100 / dataCentral->GetBinContent(bin)) ) ;
                }
            }
        }
        
        hStatError->SetBinContent(bin, 0.);
        if (dataCentral->GetBinContent(bin) > 0){
            hStatError->SetBinContent( bin, (dataCentral->GetBinError(bin)) * (100 / dataCentral->GetBinContent(bin)) ) ;
            //cout << " Binth: " << bin << " Stat  difference = " << dataCentral->GetBinError(bin)  << "  % difference = " << (dataCentral->GetBinError(bin)) * (100 / dataCentral->GetBinContent(bin))  << endl;
        }
        //cout << endl;
    }
    
    // sum all
    for (int bin(1); bin <= nBins; bin++){
        hTotalErrorUp->SetBinContent(bin, 0.);
        hTotalErrorDown->SetBinContent(bin, 0.);
        
        double sumErrorUp(0.);
        double sumErrorDown(0.);
        for (int syst(0); syst < group; syst++){
            sumErrorUp += pow(hUpDiffer[syst]->GetBinContent(bin), 2);
            sumErrorDown += pow(hDownDiffer[syst]->GetBinContent(bin), 2);
        }
        sumErrorUp += pow(hStatError->GetBinContent(bin), 2);
        sumErrorDown += pow(hStatError->GetBinContent(bin), 2);
        sumErrorUp = sqrt(sumErrorUp);
        sumErrorDown = sqrt(sumErrorDown);
        
        hTotalErrorUp->SetBinContent(bin, sumErrorUp) ;
        hTotalErrorDown->SetBinContent(bin, sumErrorDown) ;
    }
    
    // canvas1
    TCanvas* canvas1 = new TCanvas("c1", "canvas1", 900, 900);
    canvas1->cd();
    canvas1->SetLogy();
    canvas1->SetTicks();
    
    string tempXTitle;
    if (variable == "MeanNJetsHT_Zinc1jet"){
        tempXTitle = "H_{T}(jets) [GeV]";
    }
    else if (variable == "MeanNJetsHT_Zinc2jet"){
        tempXTitle = "H_{T}(jets) [GeV]";
    }
    else if (variable == "MeanNJetsdRapidity_Zinc2jet"){
        tempXTitle = "#Delta y(j_{1}j_{2})";
    }
    else if (variable == "MeanNJetsdRapidityFB_Zinc2jet"){
        tempXTitle = "#Delta y(j_{F}j_{B})";
    }
    
    if (variable.find("ZNGood") != string::npos){
        hUpDiffer[0]->GetXaxis()->SetRange(2, 8);
    }
    hUpDiffer[0]->SetMaximum(1000);
    hUpDiffer[0]->SetMinimum(0.1);
    hUpDiffer[0]->SetTitle("Systematic Up");
    hUpDiffer[0]->GetYaxis()->SetTitle("Uncertainty (%)");
    hUpDiffer[0]->GetXaxis()->SetTitle(tempXTitle.c_str());
    hUpDiffer[0]->SetLineColor(kRed);
    hUpDiffer[0]->SetLineWidth(2);
    hUpDiffer[0]->SetLineStyle(2);
    hUpDiffer[0]->DrawCopy("HIST");
    
    hUpDiffer[1]->SetLineColor(kBlue);
    hUpDiffer[1]->SetLineWidth(2);
    hUpDiffer[1]->SetLineStyle(5);
    hUpDiffer[1]->DrawCopy("HIST same");
    
    hUpDiffer[2]->SetLineColor(kGreen);
    hUpDiffer[2]->SetLineWidth(2);
    hUpDiffer[2]->SetLineStyle(5);
    hUpDiffer[2]->DrawCopy("HIST same");
    
    hUpDiffer[3]->SetLineColor(kGray+2);
    hUpDiffer[3]->SetLineWidth(3);
    hUpDiffer[3]->SetLineStyle(3);
    hUpDiffer[3]->DrawCopy("HIST same");
    
    hUpDiffer[4]->SetLineColor(kCyan);
    hUpDiffer[4]->SetLineWidth(3);
    hUpDiffer[4]->SetLineStyle(3);
    hUpDiffer[4]->DrawCopy("HIST same");
    
    hUpDiffer[5]->SetLineColor(kGreen+3);
    hUpDiffer[5]->SetLineWidth(2);
    hUpDiffer[5]->SetLineStyle(5);
    hUpDiffer[5]->DrawCopy("HIST same");
    
    hUpDiffer[6]->SetLineColor(kOrange);
    hUpDiffer[6]->SetLineWidth(3);
    hUpDiffer[6]->SetLineStyle(3);
    hUpDiffer[6]->DrawCopy("HIST same");
    
    hUpDiffer[7]->SetLineColor(kViolet+5);
    hUpDiffer[7]->SetLineWidth(2);
    hUpDiffer[7]->SetLineStyle(5);
    hUpDiffer[7]->DrawCopy("HIST same");
    
    hUpDiffer[8]->SetLineColor(kPink+2);
    hUpDiffer[8]->SetLineWidth(2);
    hUpDiffer[8]->SetLineStyle(5);
    hUpDiffer[8]->DrawCopy("HIST same");
    
    hUpDiffer[9]->SetLineColor(kSpring+9);
    hUpDiffer[9]->SetLineWidth(2);
    hUpDiffer[9]->SetLineStyle(3);
    hUpDiffer[9]->DrawCopy("HIST same");
    
    hUpDiffer[10]->SetLineColor(kAzure+4);
    hUpDiffer[10]->SetLineWidth(2);
    hUpDiffer[10]->SetLineStyle(2);
    hUpDiffer[10]->DrawCopy("HIST same");
    
    hUpDiffer[nGroup-1]->SetLineColor(kYellow);
    hUpDiffer[nGroup-1]->SetLineWidth(3);
    hUpDiffer[nGroup-1]->SetLineStyle(3);
    hUpDiffer[nGroup-1]->DrawCopy("HIST same");
    
    hStatError->SetLineColor(kMagenta);
    hStatError->SetLineWidth(2);
    hStatError->SetLineStyle(2);
    hStatError->DrawCopy("HIST same");
    
    hTotalErrorUp->SetLineColor(kBlack);
    hTotalErrorUp->SetLineWidth(2);
    hTotalErrorUp->SetLineStyle(1);
    hTotalErrorUp->DrawCopy("HIST same");
    
    TLegend *leg = new TLegend(0.15, 0.63, 0.34, 0.89);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hUpDiffer[0]   , "JES" , "l");
    leg->AddEntry(hUpDiffer[1]   , "PU"  , "l");
    leg->AddEntry(hUpDiffer[2]   , "XSEC", "l");
    leg->AddEntry(hUpDiffer[3]   , "JER" , "l");
    leg->AddEntry(hUpDiffer[4]   , "MC"  , "l");
    leg->AddEntry(hUpDiffer[5]   , "ttbar"    , "l");
    leg->AddEntry(hUpDiffer[6]   , "Wb"    , "l");
    leg->AddEntry(hUpDiffer[7]   , "Int.Lumi"    , "l");
    leg->AddEntry(hUpDiffer[8]   , "Btag Corr"  , "l");
    leg->AddEntry(hUpDiffer[9]  , "MES"    , "l");
    leg->AddEntry(hUpDiffer[10] , "MER"    , "l");
    leg->AddEntry(hUpDiffer[nGroup-1]  , "Resp", "l");
    
    leg->AddEntry(hStatError     , "Statistical", "l");
    leg->AddEntry(hTotalErrorUp  , "Total"      , "l");
    leg->Draw();
    
    string outputFileNameErrorUpPNG = outputDirectory + variable + "ErrorUp.pdf";
    canvas1->Print(outputFileNameErrorUpPNG.c_str());
    
    // canvas2
    TCanvas* canvas2 = new TCanvas("c2", "canvas2", 900, 900);
    canvas2->cd();
    canvas2->SetLogy();
    canvas2->SetTicks();
    
    if (variable.find("ZNGood") != string::npos){
        hDownDiffer[0]->GetXaxis()->SetRange(2, 8);
    }
    hDownDiffer[0]->SetMaximum(1000);
    hDownDiffer[0]->SetMinimum(0.1);
    hDownDiffer[0]->SetTitle("Systematic Down");
    hDownDiffer[0]->GetYaxis()->SetTitle("Uncertainty (%)");
    hDownDiffer[0]->GetXaxis()->SetTitle(tempXTitle.c_str());
    hDownDiffer[0]->SetLineColor(kRed);
    hDownDiffer[0]->SetLineWidth(2);
    hDownDiffer[0]->SetLineStyle(2);
    hDownDiffer[0]->DrawCopy("HIST");
    
    hDownDiffer[1]->SetLineColor(kBlue);
    hDownDiffer[1]->SetLineWidth(2);
    hDownDiffer[1]->SetLineStyle(5);
    hDownDiffer[1]->DrawCopy("HIST same");
    
    hDownDiffer[2]->SetLineColor(kGreen);
    hDownDiffer[2]->SetLineWidth(2);
    hDownDiffer[2]->SetLineStyle(5);
    hDownDiffer[2]->DrawCopy("HIST same");
    
    hDownDiffer[3]->SetLineColor(kGray+2);
    hDownDiffer[3]->SetLineWidth(3);
    hDownDiffer[3]->SetLineStyle(3);
    hDownDiffer[3]->DrawCopy("HIST same");
    
    hDownDiffer[4]->SetLineColor(kCyan);
    hDownDiffer[4]->SetLineWidth(3);
    hDownDiffer[4]->SetLineStyle(3);
    hDownDiffer[4]->DrawCopy("HIST same");
    
    hDownDiffer[5]->SetLineColor(kGreen+3);
    hDownDiffer[5]->SetLineWidth(2);
    hDownDiffer[5]->SetLineStyle(5);
    hDownDiffer[5]->DrawCopy("HIST same");
    
    hDownDiffer[6]->SetLineColor(kOrange);
    hDownDiffer[6]->SetLineWidth(3);
    hDownDiffer[6]->SetLineStyle(3);
    hDownDiffer[6]->DrawCopy("HIST same");
    
    hDownDiffer[7]->SetLineColor(kViolet+5);
    hDownDiffer[7]->SetLineWidth(2);
    hDownDiffer[7]->SetLineStyle(5);
    hDownDiffer[7]->DrawCopy("HIST same");
    
    hDownDiffer[8]->SetLineColor(kPink+2);
    hDownDiffer[8]->SetLineWidth(2);
    hDownDiffer[8]->SetLineStyle(5);
    hDownDiffer[8]->DrawCopy("HIST same");
    
    hDownDiffer[9]->SetLineColor(kSpring+9);
    hDownDiffer[9]->SetLineWidth(2);
    hDownDiffer[9]->SetLineStyle(3);
    hDownDiffer[9]->DrawCopy("HIST same");
    
    hDownDiffer[10]->SetLineColor(kAzure+4);
    hDownDiffer[10]->SetLineWidth(2);
    hDownDiffer[10]->SetLineStyle(2);
    hDownDiffer[10]->DrawCopy("HIST same");
    
    hDownDiffer[nGroup-1]->SetLineColor(kYellow);
    hDownDiffer[nGroup-1]->SetLineWidth(3);
    hDownDiffer[nGroup-1]->SetLineStyle(3);
    hDownDiffer[nGroup-1]->DrawCopy("HIST same");
    
    hStatError->SetLineColor(kMagenta);
    hStatError->SetLineWidth(2);
    hStatError->SetLineStyle(2);
    hStatError->DrawCopy("HIST same");
    
    hTotalErrorDown->SetLineColor(kBlack);
    hTotalErrorDown->SetLineWidth(2);
    hTotalErrorDown->SetLineStyle(1);
    hTotalErrorDown->DrawCopy("HIST same");
    
    TLegend *leg2 = new TLegend(0.15, 0.63, 0.34, 0.89);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->AddEntry(hDownDiffer[0]   , "JES" , "l");
    leg2->AddEntry(hDownDiffer[1]   , "PU"  , "l");
    leg2->AddEntry(hDownDiffer[2]   , "XSEC", "l");
    leg2->AddEntry(hDownDiffer[3]   , "JER" , "l");
    leg2->AddEntry(hDownDiffer[4]   , "MC"  , "l");
    leg2->AddEntry(hDownDiffer[5]   , "ttbar"    , "l");
    leg2->AddEntry(hDownDiffer[6]   , "Wb"    , "l");
    leg2->AddEntry(hDownDiffer[7]   , "Int.Lumi"    , "l");
    leg2->AddEntry(hDownDiffer[8]   , "Btag Corr"  , "l");
    leg2->AddEntry(hDownDiffer[9]  , "MES"    , "l");
    leg2->AddEntry(hDownDiffer[10] , "MER"    , "l");
    leg2->AddEntry(hDownDiffer[nGroup-1]   , "Resp", "l");
    
    leg2->AddEntry(hStatError       , "Statistical", "l");
    leg2->AddEntry(hTotalErrorDown  , "Total"      , "l");
    leg2->Draw();
    
    string outputFileNameErrorDownPNG = outputDirectory + variable + "ErrorDown.pdf";
    canvas2->Print(outputFileNameErrorDownPNG.c_str());
    
    string outputTableErrorDown = outputDirectory + "Table/" + variable + "TableErrorDown.tex";
    string outputTableErrorUp   = outputDirectory + "Table/" + variable + "TableErrorUp.tex";
    
    GetSysErrorTable(variable, outputTableErrorDown, dataCenBackUp, hDownDiffer, hStatError, hTotalErrorDown, nGroup);
    GetSysErrorTable(variable, outputTableErrorUp,   dataCenBackUp, hUpDiffer,   hStatError, hTotalErrorUp, nGroup);
    
}


void GetSysErrorTable (string variable, string outputTableName, TH1D* dataCenBackUp, TH1D* hDiffer[], TH1D* hStatError, TH1D* hTotalError, const int nGroup)
{
    cout << "Now producing table of " << outputTableName << endl;
    string command = "mkdir -p PNGFiles/FinalUnfoldMeanNJets/Table/" ;
    system(command.c_str());
    
    string tableTitle;
    if (outputTableName.find("TableErrorDown") != string::npos){
        tableTitle = variable + "TableErrorDown";
    }
    if (outputTableName.find("TableErrorUp") != string::npos){
        tableTitle = variable + "TableErrorUp";
    }
    
    int nXbin = hStatError->GetNbinsX();
    int nbinStart = 1;
    int group = nGroup;
    if (variable.find("ZNGood") != string::npos){
        nbinStart = 2;
        group = group -1;
        nXbin = 8;
    }
    
    //------------ set column headers
    string xtitle;
    string temp = "<N_{jets}>";
    
    if (variable == "MeanNJetsHT_Zinc1jet"){
        xtitle = "H_{T}(jets) [GeV]";
    }
    else if (variable == "MeanNJetsHT_Zinc2jet"){
        xtitle = "H_{T}(jets) [GeV]";
    }
    else if (variable == "MeanNJetsdRapidity_Zinc2jet"){
        xtitle = "#Delta y(j_{1}j_{2})";
    }
    else if (variable == "MeanNJetsdRapidityFB_Zinc2jet"){
        xtitle = "#Delta y(j_{F}j_{B})";
    }
    
    
    // formating text to match latex style
    size_t pos_m;
    if (xtitle.find("#") != string::npos) {
        pos_m = xtitle.find("#");
        xtitle.replace(pos_m, 1, "\\");
    }
    if (temp.find("#") != string::npos) {
        pos_m = temp.find("#");
        temp.replace(pos_m, 1, "\\");
    }
    if (tableTitle.find("_") != string::npos) {
        pos_m = tableTitle.find("_");
        tableTitle.replace(pos_m, 1, "\\_");
    }
    //-------------------
    
    
    
    FILE *outFile = fopen(outputTableName.c_str(),"w");
    fprintf(outFile, "\\input{header.tex} \n ");
    fprintf(outFile, "\\usepackage[cm]{fullpage} \n ");
    fprintf(outFile, "\\begin{document} \n ");
    fprintf(outFile, "%s \\\\ \n ", tableTitle.c_str());
    fprintf(outFile, "\\\\ \n ");
    fprintf(outFile, "\\resizebox{\\textwidth}{!}     {\n\\begin{tabular}{l|ccccccccccccccc} \n ");
    fprintf(outFile, "$%s$ & ", xtitle.c_str());
    fprintf(outFile, "$%s$ & ", temp.c_str());
    fprintf(outFile, "$Tot. Err.(\\%%)$ & $Stat(\\%%)$ & $JES(\\%%)$ & $PU(\\%%)$ & $XSEC(\\%%)$ & $JER(\\%%)$ & $MC(\\%%)$ & $ttbar(\\%%)$ & $Wb(\\%%)$ & $Int.Lumi(\\%%)$ & $BtagCorr(\\%%)$ & $MES(\\%%)$ & $MER(\\%%)$ & $Resp(\\%%)$ \\\\ \\hline \n ");
    
    // get array of bin edges
    double arrX[50];
    const double *arr = new double[nXbin+1];
    arr = dataCenBackUp->GetXaxis()->GetXbins()->GetArray();
    if (arr == 0){ // non-variable binning is defined
        if (variable.find("ZNGood") != string::npos) {
        }
        else{
            for (int j = 0; j < nXbin ; j++){
                arrX[j] = dataCenBackUp->GetXaxis()->GetBinLowEdge(j+1);
            }
            arrX[nXbin] = dataCenBackUp->GetXaxis()->GetBinUpEdge(nXbin);
        }
    }
    else{ // variable binning is defined
        for (int j = 0; j <= nXbin ; j++){
            arrX[j] = dataCenBackUp->GetXaxis()->GetXbins()->At(j);
        }
    }
    
    // fill out the table
    for (int i = nbinStart; i<= nXbin ; i++){
        if (variable.find("ZNGood") != string::npos) {
            fprintf(outFile, " $N_{\\text{jets}} =$ %d  & ", i-1);
        }
        else{
            fprintf(outFile, " %g - %g  & ", arrX[i-1], arrX[i]);
        }
        fprintf(outFile, "%.3g  &", double(dataCenBackUp->GetBinContent(i)));
        fprintf(outFile, "%.3g  &", double(hTotalError->GetBinContent(i)));
        fprintf(outFile, "%.3g  &", double(hStatError->GetBinContent(i)));
        
        for (int j = 0 ; j < group  ; j++ ) {
            if (j != group -1){
                fprintf(outFile, "%.3g & ", double(hDiffer[j]->GetBinContent(i)));
            }
            else {
                fprintf(outFile, "%.3g ", double(hDiffer[j]->GetBinContent(i)));
            }
        }
        fprintf(outFile, "\\\\ \n");
    }
    fprintf( outFile, "\\end{tabular}} \\\\ \\\\  \n");
    fprintf( outFile, "\\end{document}");
    
    fclose(outFile);
}
//cout << "Stop after line " << __LINE__ << endl;
////////////////////////////////////////////////////////////////////////////////////


