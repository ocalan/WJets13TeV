#include <TFile.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <RooUnfoldResponse.h>
#include <RooUnfoldSvd.h>
#include <RooUnfoldBayes.h>
#include <RooUnfoldBinByBin.h>
#include <RooUnfoldInvert.h>
#include <TSVDUnfold.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TList.h>
#include "getFilesAndHistograms.h"
#include "unfoldingFunctions.h"
#include "writeFastPlotsTex.h"
#include "myFastPlots.h"

//-- Setting global variables --------------------------------------------------------------
#include "fileNames.h"
#include "variablesOfInterestVarWidth.h"

void myFastPlots(string leptonFlavor, string var, bool closureTest)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("1.2f");

    //--- load the proper VAROFINTERST structure: Z+jets or W+jets ---
    int nsel = 0;
    const int *sel = NULL;
    variableStruct *varStruct = NULL;
    TList *varList = new TList();
    int nVarStruct;
    if (leptonFlavor == "DMu" || leptonFlavor == "DE"){
        nVarStruct = NVAROFINTERESTZJETS;
        varStruct = new variableStruct[nVarStruct];
        varStruct = VAROFINTERESTZJETS;
        makeLISTOFVAROFINTERESTZJETS();
        varList = LISTOFVAROFINTERESTZJETS;
        nsel = NFILESDYJETS;
        sel = new int[nsel];
        sel = FilesDYJets;
    }
    else if (leptonFlavor == "SMu" || leptonFlavor == "SE"){
        nVarStruct = NVAROFINTERESTWJETS;
        varStruct = new variableStruct[nVarStruct];
        varStruct = VAROFINTERESTWJETS;
        makeLISTOFVAROFINTERESTWJETS();
        varList = LISTOFVAROFINTERESTWJETS;
        nsel = NFILESTTBAR;
        sel = new int[nsel];
        sel = FilesTTbarWJets;
    }

    //---------------------------------------------

    //--- check wether muon or electron is used ---
    bool isMuon;
    if (leptonFlavor.find("Mu") != string::npos) isMuon = 1;
    else isMuon = 0;
    //---------------------------------------------
    
    if (var == ""){
        //nVarStruct = 1;
        for (int i(40); i < 50/*nVarStruct*/; i++){
            cout << setw(3) << i << ") Processing variable: " << varStruct[i].name << endl; 
            if (!isMuon) FastPlotsRun(sel, nsel, leptonFlavor, varStruct[i].name, varStruct[i].log, varStruct[i].decrease, varStruct[i].ESVDkterm, varStruct[i].EBayeskterm);
            if (isMuon)  FastPlotsRun(sel, nsel, leptonFlavor, varStruct[i].name, varStruct[i].log, varStruct[i].decrease, varStruct[i].MuSVDkterm, varStruct[i].MuBayeskterm);
            cout << "Type \"n\" for next" << endl;
            int a = 0;
            cin >> a;
            cin.clear();
            cin.ignore();
        }
    }
    else{
        if(varList->Contains(var.c_str())){
            TObjString *name = new TObjString(var.c_str());
            unsigned short j(varList->IndexOf(name));
            cout << setw(3) << " " << ") Processing variable: " << varStruct[j].name << endl; 
            if (!isMuon) FastPlotsRun(sel, nsel, leptonFlavor, varStruct[j].name, varStruct[j].log, varStruct[j].decrease, varStruct[j].ESVDkterm, varStruct[j].EBayeskterm, closureTest);
            if (isMuon)  FastPlotsRun(sel, nsel, leptonFlavor, varStruct[j].name, varStruct[j].log, varStruct[j].decrease, varStruct[j].MuSVDkterm, varStruct[j].MuBayeskterm, closureTest);
        }
        else{
            if (!isMuon) FastPlotsRun(sel, nsel, leptonFlavor, var.c_str(),1,1,5,5);
            if (isMuon)  FastPlotsRun(sel, nsel, leptonFlavor, var.c_str(),1,1,5,5);
        }
    }

}


void FastPlotsRun(const int *sel, int nsel, string leptonFlavor, string variable, bool logZ, bool decrease, int SVDkterm, int Bayeskterm, bool closureTest, int JetPtMin, int JetPtMax, bool doFlat, bool doVarWidth)
{
    //--- make sure we have trailing "_" at the end of leptonFlavor
    if (leptonFlavor == "Muons" || leptonFlavor == "DMu") leptonFlavor = "DMu_";
    else if (leptonFlavor == "Electrons" || leptonFlavor == "DE") leptonFlavor = "DE_";
    else if (leptonFlavor == "Muon" || leptonFlavor == "SE") leptonFlavor = "SE_";
    else if (leptonFlavor == "Electron" || leptonFlavor == "SMu") leptonFlavor = "SMu_";

    //--- retreive the year from current location and set the energy accordingly
    string energy = getEnergy();

    //if (variable.find("JetPt")!= string::npos && variable.find("Highest") == string::npos && JetPtMin == 30)  JetPtMin = 20 ;
    ostringstream jetPtMin; jetPtMin << JetPtMin;
    string outputDirectory = "PNGFiles/FastPlots_" + jetPtMin.str() ;
    if (doVarWidth)  outputDirectory+= "_VarWidth";
    string command = "mkdir -p " + outputDirectory;
    system(command.c_str());
    outputDirectory +=  "/";

    string outputRootFileName = outputDirectory + leptonFlavor + energy + "_" + variable + "_fastplots";
    if (closureTest) outputRootFileName += "_ClosureTest";
    if (doVarWidth) outputRootFileName += "_VarWidth";
    string outputTexFileName = outputRootFileName;
    outputRootFileName += ".root";
    outputTexFileName += ".tex";

    TFile *outputRootFile = new TFile(outputRootFileName.c_str(), "RECREATE");
    ofstream outputTexFile(outputTexFileName.c_str());
    writeHeaderPage(outputTexFile, leptonFlavor, variable);

    //-- get data histogram (or MC Even for closure test) ----------------------------
    TFile *DataFile = NULL;
    if (!closureTest) DataFile = getFile(FILESDIRECTORY, leptonFlavor, energy, ProcessInfo[sel[0]].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, 0, 0, 0, 0, -1, "", "0", false, 0);
    else DataFile = getFile(FILESDIRECTORY, leptonFlavor, energy, ProcessInfo[sel[nsel-1]].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, 0, 0, 0, 0, -1, "", "0", false, 0);
    
    TH1D *meas = getHisto(DataFile, variable);
    

    //-- get MC madgraph histograms (or MC Odd for closure test) ---------------------
    TFile *DY = NULL;
    if (!closureTest) DY = getFile(FILESDIRECTORY, leptonFlavor, energy, ProcessInfo[sel[nsel-1]].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, 0, 0, 0, 0, -1, "", "0", false, 0);
    else DY = getFile(FILESDIRECTORY, leptonFlavor, energy, ProcessInfo[sel[nsel-1]].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, 0, 0, 0, 0, -1, "", "0", false, 0);
    
    TH1D *recoMad = (TH1D*) getHisto(DY, variable);
    outputRootFile->cd();
    recoMad->Write("recoMad");

    //-- get Powheg histograms -------------------------------------------------------
    //TFile *Powheg = getFile(FILESDIRECTORY, leptonFlavor, energy, DYPOWHEGFILENAME, JetPtMin, JetPtMax, doFlat);

    //-- get Powheg histograms -------------------------------------------------------
    //TFile *Sherpa = getFile(FILESDIRECTORY, leptonFlavor, energy, DYSHERPAFILENAME, JetPtMin, JetPtMax, doFlat);

    //-- get BG histograms -----------------------------------------------------------
    int nBG(nsel-2);
    if (closureTest) nBG = 0;
    TFile **BG = new TFile*[nBG];
    TH1D **hBG = new TH1D*[nBG];
    for (int i(0); i < nBG; i++){
        BG[i] = getFile(FILESDIRECTORY, leptonFlavor, energy, ProcessInfo[sel[i+1]].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, 0, 0, 0, 0, -1, "", "0", false, 0);
        hBG[i] = getHisto(BG[i], variable);
                
    }
    //TH1D *measSubBG = (TH1D*) meas->Clone();//== some functions will need this as an input == need to add code to build this)//
    

    //-- build the repsonse matrix, MCtruth, MCreco ------------------------
    string respName = "response" + variable;
    RooUnfoldResponse *response = (RooUnfoldResponse*) DY->Get(respName.c_str());
    response->UseOverflow();
    string genName = "gen" + variable;
    TH1D *genMad = NULL;
    if (!closureTest) genMad = (TH1D*) DY->Get(genName.c_str());
    else genMad = (TH1D*) DY->Get(genName.c_str());

    TH2D *hResponse = (TH2D*) response->Hresponse();
    TH2D *hNormResp = (TH2D*) hResponse->Clone(); 
    normalizeTH2D(hNormResp);
    //TH1D *MCtruth   = (TH1D*) response->Htruth();
    //TH1D *MCreco    = (TH1D*) response->Hmeasured();
    

    //-- get TH2D from Powheg -------------------------------------------------------
    //RooUnfoldResponse *respPow = (RooUnfoldResponse*) Powheg->Get(respName.c_str());
    //TH2D *hRespPow = (TH2D*) respPow->Hresponse();

    //-- get TH2D from Sherpa -------------------------------------------------------
    //RooUnfoldResponse *respShe = (RooUnfoldResponse*) Sherpa->Get(respName.c_str());
    //TH2D *hRespShe = (TH2D*) respShe->Hresponse();

    
    //-- draw normalized response matrix --------------------------------------------
    plotHNormResp(hNormResp, leptonFlavor, variable, energy, outputDirectory, outputRootFile, closureTest);
    
    
    //-- draw Bayes selected kTerm and +/- 1 ----------------------------------------
    myplotSelectedMethod2("Bayes", response, genMad, Bayeskterm, meas, hBG, nBG, leptonFlavor, variable, logZ, decrease, outputDirectory, outputRootFile, closureTest, true);
    
    //-- draw chi2 of change --------------------------------------------------------
    myplotChi2OfChange(response, Bayeskterm, meas, hBG, nBG, leptonFlavor, variable, energy, outputDirectory, outputRootFile, closureTest);
    
    
    //-- Bin by bin method ----------------------------------------------------------
    //myplotSelectedMethod2("BinByBin", response, genMad, 0, meas, hBG, nBG, leptonFlavor, variable, logZ, decrease, outputDirectory, outputRootFile, closureTest, true);
    
    
    //-- draw SVD selected kTerm and +/- 1 ------------------------------------------
    //myplotSelectedMethod2("SVD", response, genMad, SVDkterm, meas, hBG, nBG, leptonFlavor, variable, logZ, decrease, outputDirectory, outputRootFile, closureTest, true);

    //-- close every files ----------------------------------------------------------
    outputRootFile->Close();
    closeFile(DataFile);
    closeFile(DY);
    //closeFile(Powheg);
    //closeFile(Sherpa);
    for (int i(0); i < nBG; i++) closeFile(BG[i]);
    outputTexFile.close();
}


void myplotSelectedMethod2(string method, RooUnfoldResponse *response, TH1D *genMad, int kterm, TH1D *hData, TH1D *hBG[], int nBG, string leptonFlavor, string variable, bool logZ, bool decrease, string outputDirectory, TFile *outputRootFile, bool closureTest, bool save)
{
    cout << " Start processing Unfolding using method: " << method << endl;
    
    if (leptonFlavor == "Muons" || leptonFlavor == "DMu") leptonFlavor = "DMu_";
    else if (leptonFlavor == "Electrons" || leptonFlavor == "DE") leptonFlavor = "DE_";
    else if (leptonFlavor == "Muon" || leptonFlavor == "SE") leptonFlavor = "SE_";
    else if (leptonFlavor == "Electron" || leptonFlavor == "SMu") leptonFlavor = "SMu_";
    
    string energy = getEnergy();
    if (kterm < 1) kterm = 2;
    ostringstream ktermStr; ktermStr << kterm;
    
    TH2D *hResponse = (TH2D*) response->Hresponse();
    TH1D *MCreco    = (TH1D*) response->Hmeasured();
    
    TH1D *gen       = (TH1D*) genMad->Clone();
    
    string yTitle = MCreco->GetXaxis()->GetTitle();
    string respTitle = hResponse->GetTitle();
    string niceTitle = respTitle.substr(0, respTitle.find("#right"));
    niceTitle = niceTitle.substr(niceTitle.find("Response")+9);
    if (niceTitle.find("gen") != string::npos) niceTitle = niceTitle.substr(niceTitle.find("gen")+4);
    
    string title = method + " Unfolded Distributions for " + niceTitle + " (kTerm = " + ktermStr.str() + " #pm 1)";
    if (closureTest) title += " C.T.";
    string name = leptonFlavor + energy + "_" + variable + "_Unfolded_Distribution_" + method;
    if (closureTest) name += "_CT";
    
    TH1D *hMeas = (TH1D*) hData->Clone();
    int nBins(hMeas->GetNbinsX());
    
    for (int i(0); i < nBG; i++) {
        hMeas->Add(hBG[i], -1);
    }
    
    
    TH1D *hUnfoldedC = NULL, *hUnfoldedD = NULL, *hUnfoldedU = NULL;
    if (method == "SVD") {
        cout << " Unfolding using Center Kterm " << endl;
        RooUnfoldSvd unfoldC(response, hMeas, kterm);
        unfoldC.SetVerbose(0);
        hUnfoldedC = (TH1D*) unfoldC.Hreco();
        if (kterm > 1){
            cout << " Unfolding using Down Kterm " << endl;
            RooUnfoldSvd unfoldD(response, hMeas, kterm - 1);
            unfoldD.SetVerbose(0);
            hUnfoldedD = (TH1D*) unfoldD.Hreco();
        }
        else cout << "       Impossible to show kterm - 1" << endl;
        
        if (kterm < nBins){
            cout << " Unfolding using Up Kterm" << endl;
            RooUnfoldSvd unfoldU(response, hMeas, kterm + 1);
            unfoldU.SetVerbose(0);
            hUnfoldedU = (TH1D*) unfoldU.Hreco();
        }
        else cout << "       Impossible to show kterm + 1" << endl;
        
    }
    else if (method == "Bayes") {
        cout << " Unfolding using Center Kterm " << endl;
        RooUnfoldBayes unfoldC(response, hMeas, kterm);
        unfoldC.SetVerbose(0);
        hUnfoldedC = (TH1D*) unfoldC.Hreco();
        if (kterm > 1){
            cout << " Unfolding using Down Kterm " << endl;
            RooUnfoldBayes unfoldD(response, hMeas, kterm - 1);
            unfoldD.SetVerbose(0);
            hUnfoldedD = (TH1D*) unfoldD.Hreco();
        }
        else cout << "       Impossible to show kterm - 1" << endl;
        if (kterm < nBins){
            cout << " Unfolding using Up Kterm" << endl;
            RooUnfoldBayes unfoldU(response, hMeas, kterm + 1);
            unfoldU.SetVerbose(0);
            hUnfoldedU = (TH1D*) unfoldU.Hreco();
        }
        else cout << "       Impossible to show kterm + 1" << endl;
        
    }
    else if (method == "BinByBin"){
        RooUnfoldBinByBin unfoldC(response, hMeas);
        unfoldC.SetVerbose(0);
        hUnfoldedC = (TH1D*) unfoldC.Hreco();
    }
    
    // save unfolded histo to root file
    outputRootFile->cd();
    hMeas->Write("DataSubBG");
    hUnfoldedC->Write("Central");
    genMad->Write("genMad");
    
    TLatex *cmsColl = new TLatex();
    cmsColl->SetTextSize(0.025);
    cmsColl->SetTextFont(42);
    cmsColl->SetLineWidth(2);
    cmsColl->SetTextColor(kBlack);
    cmsColl->SetNDC();
    cmsColl->SetTextAlign(11);
    
    TLatex *sqrtXTeV = new TLatex();
    sqrtXTeV->SetTextSize(0.025);
    sqrtXTeV->SetTextFont(42);
    sqrtXTeV->SetLineWidth(2);
    sqrtXTeV->SetTextColor(kBlack);
    sqrtXTeV->SetNDC();
    sqrtXTeV->SetTextAlign(21);
    
    TLatex *intLumi = new TLatex();
    intLumi->SetTextSize(0.022);
    intLumi->SetTextFont(42);
    intLumi->SetLineWidth(2);
    intLumi->SetTextColor(kBlack);
    intLumi->SetNDC();
    intLumi->SetTextAlign(31);
    
    
    TCanvas *canUnfoldedDistribution = new TCanvas(name.c_str(), name.c_str(), 700, 900);
    canUnfoldedDistribution->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
    pad1->SetTopMargin(0.11);
    pad1->SetBottomMargin(0);
    
    if (variable == "FirstJetPt_2_Zinc1jet") pad1->SetRightMargin(0.04);
    else pad1->SetRightMargin(0.02);
    
    pad1->Draw();
    pad1->cd();
    if (logZ) pad1->SetLogy();
    gen->SetTitle(title.c_str());
    gen->GetYaxis()->SetTitleOffset(1.4);
    gen->SetLineColor(kBlack);
    gen->SetLineWidth(2);
    gen->SetFillStyle(3003);
    gen->SetFillColor(kGray+2);
    gen->DrawCopy("HIST");
    
    
    if (kterm > 1 && method != "BinByBin") {
            hUnfoldedD->SetTitle("");
            hUnfoldedD->SetMarkerStyle(32);
            hUnfoldedD->SetMarkerColor(kRed+1);
            hUnfoldedD->SetLineColor(kRed+1);
            hUnfoldedD->SetLineWidth(2);
            hUnfoldedD->DrawCopy("E1same");
    }
    if (kterm < nBins && method != "BinByBin") {
            hUnfoldedU->SetTitle("");
            hUnfoldedU->SetMarkerStyle(26);
            hUnfoldedU->SetMarkerColor(kGreen+2);
            hUnfoldedU->SetLineColor(kGreen+2);
            hUnfoldedU->SetLineWidth(2);
            hUnfoldedU->DrawCopy("E1same");
    }
    
    
    hUnfoldedC->SetMarkerStyle(24);
    hUnfoldedC->SetMarkerColor(kBlue+1);
    hUnfoldedC->SetLineColor(kBlue+1);
    hUnfoldedC->SetLineWidth(2);
    hUnfoldedC->DrawCopy("E1same");
    
    double xLowLeg(0.76), xHighLeg(0.95);
    if (!decrease){
        xLowLeg = 0.15;
        xHighLeg = 0.34;
    }
    
    //=== legend ====//
    TLegend *leg = new TLegend(xLowLeg, 0.63, xHighLeg, 0.89);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    if (!closureTest) leg->AddEntry(gen, "MC Gen.", "fl");
    else leg->AddEntry(gen, "MC Gen.", "fl");
    
    string legEntry = "kTerm = " + ktermStr.str();
    
    if (method != "BinByBin") {leg->AddEntry(hUnfoldedC, legEntry.c_str(), "pl");}
    else if (method == "BinByBin"){leg->AddEntry(hUnfoldedC, "BinByBin", "pl");}
    
    if (kterm < nBins && method != "BinByBin") leg->AddEntry(hUnfoldedU, "kTerm + 1", "pl");
    if (kterm > 1 && method != "BinByBin") leg->AddEntry(hUnfoldedD, "kTerm - 1", "pl");
    
    leg->Draw();
    //=== end legend ====//
    
    if (!closureTest){
        cmsColl->DrawLatex(0.11,0.905, "CMS Preliminary");
        if (energy == "7TeV") sqrtXTeV->DrawLatex(0.54,0.905, "#sqrt{s} = 7 TeV");
        if (energy == "8TeV") sqrtXTeV->DrawLatex(0.54,0.905, "#sqrt{s} = 8 TeV");
        if (energy == "7TeV") intLumi->DrawLatex(0.98,0.91, "#int L dt = 5.05 fb^{-1}");
        if (energy == "8TeV") intLumi->DrawLatex(0.98,0.91, "#int L dt = 19.6 fb^{-1}");
    }
    pad1->Draw();
    canUnfoldedDistribution->cd();
    
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    
    if (variable == "FirstJetPt_2_Zinc1jet") pad2->SetRightMargin(0.04);
    else pad2->SetRightMargin(0.02);
    
    pad2->SetTicks();
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
    hUnfoldedC->Divide(genMad);
    if (kterm > 1 && method != "BinByBin") hUnfoldedD->Divide(genMad);
    if (kterm < nBins && method != "BinByBin") hUnfoldedU->Divide(genMad);
    
    hUnfoldedC->SetTitle("");
    hUnfoldedC->GetXaxis()->SetTickLength(0.03);
    hUnfoldedC->GetXaxis()->SetTitleSize(0.1);
    hUnfoldedC->GetXaxis()->SetLabelSize(0.1);
    hUnfoldedC->GetXaxis()->SetLabelOffset(0.018);
    hUnfoldedC->GetXaxis()->SetTitleOffset(1.2);
    
    hUnfoldedC->GetYaxis()->SetRangeUser(0.51, 1.49);
    hUnfoldedC->GetYaxis()->SetNdivisions(5,5,0);
    if (!closureTest) {
        hUnfoldedC->GetYaxis()->SetTitle("Gen / Unfolded");
        hUnfoldedC->GetYaxis()->SetTitleSize(0.1);
        hUnfoldedC->GetYaxis()->SetTitleOffset(0.5);
    }
    else {
        hUnfoldedC->GetYaxis()->SetTitle("Gen. / Unf. C.T.");
        hUnfoldedC->GetYaxis()->SetTitleSize(0.07);
        hUnfoldedC->GetYaxis()->SetTitleOffset(0.7);
    }
    hUnfoldedC->GetYaxis()->CenterTitle();
    hUnfoldedC->GetYaxis()->SetLabelSize(0.08);
    
    for (int j(1); j <= nBins; j++){
        double content(hUnfoldedC->GetBinContent(j));
        double contentErr(hUnfoldedC->GetBinError(j));
        if (content > 0){
            hUnfoldedC->SetBinContent(j, 1./content);
            hUnfoldedC->SetBinError(j, contentErr/(content * content));
        }
        
        if(method != "BinByBin"){
            if(kterm > 1){
                content = (hUnfoldedD->GetBinContent(j));
                contentErr = (hUnfoldedD->GetBinError(j));
                if (content > 0){
                    hUnfoldedD->SetBinContent(j, 1./content);
                    hUnfoldedD->SetBinError(j, contentErr/(content * content));                    
                }
            }
            if(kterm < nBins){
                content = (hUnfoldedU->GetBinContent(j));
                contentErr = (hUnfoldedU->GetBinError(j));
                if (content > 0){
                    hUnfoldedU->SetBinContent(j, 1./content);
                    hUnfoldedU->SetBinError(j, contentErr/(content * content));
                }
            }
        }
    }
    
    hUnfoldedC->DrawCopy("E1");
    if (kterm > 1 && method != "BinByBin") hUnfoldedD->DrawCopy("E1same");
    if (kterm < nBins && method != "BinByBin") hUnfoldedU->DrawCopy("E1same");
    pad2->Draw();
    
    string outputFileNamePNG = outputDirectory + name.c_str() + ".pdf";
    if (save){
        canUnfoldedDistribution->Print(outputFileNamePNG.c_str());
        outputRootFile->cd();
        canUnfoldedDistribution->Write();
    }
}

void myplotChi2OfChange(RooUnfoldResponse *response, int kterm, TH1D *hData, TH1D *hBG[], int nBG, string leptonFlavor, string variable, string energy, string outputDirectory, TFile *outputRootFile, bool closureTest, bool save)
{
    if (leptonFlavor == "Muons" || leptonFlavor == "DMu") leptonFlavor = "DMu_";
    else if (leptonFlavor == "Electrons" || leptonFlavor == "DE") leptonFlavor = "DE_";
    else if (leptonFlavor == "Muon" || leptonFlavor == "SE") leptonFlavor = "SE_";
    else if (leptonFlavor == "Electron" || leptonFlavor == "SMu") leptonFlavor = "SMu_";
    
    TH2D *hResponse = (TH2D*) response->Hresponse();
    TH1D *hMeas = (TH1D*) hData->Clone();
    int nBins(hMeas->GetNbinsX()  );
    
    // not allow this to go futher than kterm = 20
    int iterNum = nBins;
    if (iterNum > 20) {iterNum = 20;}
    
    for (int i(0); i < nBG; i++){
        hMeas->Add(hBG[i], -1);
    }
    
    cout << " Start processing Chi2OfChange Plot " << endl;
    
    RooUnfoldBayes unfoldBayesChi2(response, hMeas, iterNum);
    unfoldBayesChi2.SetVerbose(0);
    TH1D *tmp = (TH1D*) unfoldBayesChi2.Hreco();
            //nBins = tmp->GetNbinsX();
    TH1D *Chi2OfChange = unfoldBayesChi2.Chi2OfChange;
    
    double hmax(Chi2OfChange->GetMaximum());
    double kCont(Chi2OfChange->GetBinContent(kterm));
    kCont *= 2;
    
    TArrow *arrow = new TArrow(kterm+0.5, hmax, kterm+0.5, kCont, 0.03, "|>");
    arrow->SetAngle(40);
    arrow->SetLineWidth(2);
    arrow->SetLineColor(kRed+2);
    arrow->SetFillColor(kRed+2);
    
    
    string respTitle = hResponse->GetTitle();
    string niceTitle = respTitle.substr(0, respTitle.find("#right"));
    niceTitle = niceTitle.substr(niceTitle.find("Response")+9);
    if (niceTitle.find("gen") != string::npos) niceTitle = niceTitle.substr(niceTitle.find("gen")+4);
    string chi2Title = "#chi^{2} of change for " + niceTitle;
    if (closureTest) chi2Title += " C.T.";
    string canNameChi2 = leptonFlavor + energy + "_" + variable + "_Chi2";
    if (closureTest) canNameChi2 += "_CT";
    
    Chi2OfChange->SetTitle(chi2Title.c_str());
    Chi2OfChange->SetLineWidth(3);
    Chi2OfChange->SetLineColor(kBlack);
    Chi2OfChange->GetXaxis()->SetTitle("kterm");
    Chi2OfChange->GetXaxis()->SetTitleSize(0.04);
    Chi2OfChange->GetXaxis()->SetTitleOffset(1.0);
    Chi2OfChange->GetYaxis()->SetTitleSize(0.04);
    Chi2OfChange->GetYaxis()->SetTitleOffset(1.15);
    Chi2OfChange->GetXaxis()->CenterLabels();
    Chi2OfChange->GetXaxis()->SetLabelSize(0.1/sqrt(iterNum));
    Chi2OfChange->GetXaxis()->SetNdivisions(iterNum, 0, 0);
    Chi2OfChange->GetYaxis()->SetTitle("#chi^{2}");
    TCanvas *canChi2 = new TCanvas(canNameChi2.c_str(), canNameChi2.c_str(), 700, 700);
    canChi2->cd();
    canChi2->SetLeftMargin(0.1);
    canChi2->SetRightMargin(0.05);
    canChi2->SetLogy();
    canChi2->SetGrid();
    canChi2->SetTicks();
    Chi2OfChange->DrawCopy("TEXT60");
    arrow->Draw();
    string outputFileNamePNG = outputDirectory + canNameChi2 + ".pdf";
    if (save){
        canChi2->Print(outputFileNamePNG.c_str());
        outputRootFile->cd();
        canChi2->Write();
    }
}



