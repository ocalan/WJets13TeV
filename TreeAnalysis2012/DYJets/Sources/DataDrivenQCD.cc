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
#include "getFilesAndHistograms.h"
#include "DataDrivenQCD.h"
using namespace std;

//--  Setting global variables --------------------------------------------------------------
#include "fileNames.h"
#include "variablesOfInterestVarWidth.h"
string outputDirectory="DataDrivenQCD/";

const int NQCD = 4 ;
TFile *fData[NQCD],*fSignal[NQCD];             // 0 = central, 1 = JES Up, 2 = JES Down 
TH1D *hData[NQCD],*hSignal[NQCD],*hBack[NQCD];
TFile *fMC[NQCD][NFILESWJETS];

bool isMuon(0);
string unfAlg = "Bayes";
string energy = getEnergy();
int JetPtMin(30);
int JetPtMax(0);
TFile *fOut;
//TFile *fOut = new TFile("DataDrivenQCD/out.root","recreate");
// TFile* outputRootFile = new TFile(outputRootFileName.c_str(), "RECREATE");

//-------------------------------------------------------------------------------------------
void DataDrivenQCD(  string leptonFlavor, int METcut , int doBJets ){
    FuncOpenAllFiles(leptonFlavor, METcut, false, true, doBJets);
    vector<string> histoNameRun = getVectorOfHistoNames();
  //  for (int i(200); i < 323/*NVAROFINTEREST*/; i++){
    for (int i(0); i < histoNameRun.size()/*NVAROFINTEREST*/; i++){
//        cout << i<<"   "<<histoNameRun[i] << endl;
        FuncDataDrivenQCD(histoNameRun[i],  leptonFlavor, METcut, doBJets);
    }
//    fData[0]->ls();    
    fOut ->Close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void FuncOpenAllFiles(string leptonFlavor,int METcut, bool doFlat , bool doVarWidth, int doBJets){
    for ( int i=0 ; i < NQCD ; i++){
        //getFiles(FILESDIRECTORY, fData, leptonFlavor, energy, DATAFILENAME, JetPtMin, JetPtMax, doFlat, doVarWidth);
                fData[i] = getFile(FILESDIRECTORY,  leptonFlavor, energy, ProcessInfo[DATAFILENAME].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, i , 0, 0, METcut,doBJets, "","0");
        if ( i == 0 ) {
              string nameQCDout = fData[i]->GetName();
              nameQCDout.insert(nameQCDout.find("Data") + 4,"QCD");
            fOut = new TFile(nameQCDout.c_str(),"recreate");
        }
    }
    /// get MC sampples
    for ( int i=0 ; i < NQCD ; i++){
        for ( int j = 1 ; j < NFILESWJETS ; j++){
            int sel = j ;
            if ( j == 1 ) sel = 24 ;
            cout << endl;
            fMC[i][j] = getFile(FILESDIRECTORY,  leptonFlavor, energy, ProcessInfo[sel].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, i , 0, 0, METcut, doBJets, "","0");
            TH1D *hTemp = getHisto(fMC[i][j], "ZNGoodJets_Zexc");
             cout << " going to fetch " << ProcessInfo[sel].filename << "   " << hTemp ->Integral()<<endl;
        }
    }
cout << " opened all data and MC files "<< endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FuncDataDrivenQCD(string variable, string leptonFlavor, int METcut, bool doFlat , bool doVarWidth, int doBjets){
    cout << " test" << endl;
    //-- fetch the data files and histograms --------------
    /* const int NQCD = 4 ;
       TFile *fData[NQCD],*fSignal[NQCD];             // 0 = central, 1 = JES Up, 2 = JES Down 
       TH1D *hData[NQCD],*hSignal[NQCD],*hBack[NQCD];
       */
    for ( int i=0 ; i < NQCD ; i++){
        //getFiles(FILESDIRECTORY, fData, leptonFlavor, energy, DATAFILENAME, JetPtMin, JetPtMax, doFlat, doVarWidth);
//        fData[i] = getFile(FILESDIRECTORY,  leptonFlavor, energy, ProcessInfo[DATAFILENAME].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, i , 0, 0, "","0");
        cout << " got data " << variable << endl;
        TH1D *hTemp = getHisto(fData[i], variable);
        hData[i] = (TH1D *) hTemp->Clone();
        cout << " got data histos " << ProcessInfo[i].filename << "    " << hData[i]->Integral()<<endl;
        //	  fData[i]->Close();
    }
    //-----------------------------------------------------
    cout << " now open all the MC of interest " << endl;
    // ----------------- get all MC histos ---
    for ( int i=0 ; i < NQCD ; i++){
        for ( int j = 1 ; j < NFILESWJETS ; j++){
            int sel = j ;
            if ( j == 1 ) sel = 24 ;
            cout << endl;
  //          TFile *fTemp = getFile(FILESDIRECTORY,  leptonFlavor, energy, ProcessInfo[sel].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, i , 0, 0,"","0");
            TH1D *hTemp = getHisto(fMC[i][j], variable);
            cout << " going to fetch " << ProcessInfo[sel].filename << "   " << hTemp ->Integral()<<endl;
            //fOut ->cd(); hTemp->Write(); 
            if ( ProcessInfo[sel].filename.find("WJetsALL") != -1 ) {
                cout << ProcessInfo[sel].filename << endl;
                hSignal[i] = (TH1D *) hTemp->Clone();
            }
            else{
                if ( j == 2 ) {hBack[i] = (TH1D *) hTemp->Clone();}
                else { cout << " adding " << endl; hBack[i]->Add(hTemp);  }
            }
            //fTemp->Close();
            cout << " got histo" << i << "    " << j << "   " <<  hTemp->Integral() <<  endl; 
        }

    }
    // QCD estiamtion procedure
    cout << " to QCD estimation " << endl;
    // followed as descirbed in AN2012_331 pages 10,11
    TH1D* scaledMC[NQCD], *hQCD[NQCD];
    // step 0 : initial normalization of Wjets and data  
    cout << " normalization of signal "  <<hData[0]->Integral() <<"   " << hSignal[0]->Integral() <<"  " <<  hBack[0]->Integral() << endl;
    //  TH1D *hSigInData = (TH1D*) hData[0]->Clone() ;
    //  hSigInData->Add(hBack[0], -1 );
    double NormFactor = (hData[0]->Integral() - hBack[0]->Integral() ) / hSignal[0]->Integral() ;
    cout << " normalization of signal " << NormFactor <<hData[0]->Integral() <<"   " << hSignal[0]->Integral() <<"  " <<  hBack[0]->Integral() << endl;
    for ( int j = 0  ; j < 12 ; j++){
    if ( j > 0 )  NormFactor = (hData[0]->Integral() - hBack[0]->Integral() - hQCD[0] ->Integral() ) / hSignal[0]->Integral() ;
    // step 1: 
    for ( int i = 0  ; i < NQCD ; i++){
        scaledMC[i] = (TH1D*) hSignal[i]->Clone();
        scaledMC[i]->Scale(NormFactor) ;
        hQCD[i] = (TH1D*) hData[i]->Clone();
        //hQCD[i] = (TH1D*) hData[i]->Clone("QCD");
        hQCD[i]->Add(scaledMC[i],-1);
        hQCD[i]->Add(hBack[i],-1);
        // set the negative bins to 0
        for ( int j = 0 ; j < hData[i]->GetNbinsX()+1 ; j++ ){
        if (  hQCD[i]->GetBinContent(j) < 0 ) hQCD[i]->SetBinContent(j, 0. ) ;
        }
        double NormTemp =  (hData[i]->Integral() - hBack[i]->Integral() ) / hSignal[i]->Integral() ;
        cout << " data vs MC : " << i <<"  " << NormTemp << "   " << hData[i]->Integral() << "   " << scaledMC[i] ->Integral() + hBack[i]->Integral() <<  endl;
    }

    cout << " QCD intgral  A to D " << hQCD[0]->Integral() << "   " << hQCD[1]->Integral() <<"   " << hQCD[2]->Integral() <<"  " << hQCD[3]->Integral()<<  endl;
    // step 2: isoaltion fake rate
    double NormFactorISO =  hQCD[1]->Integral(); //  hQCD[3]->Integral();
    if ( hQCD[3]->Integral() > 0 ) NormFactorISO /= hQCD[3]->Integral();
    else NormFactorISO = 0;
    cout << " ratio of regions B to C " << NormFactorISO << "   " << hQCD[1]->Integral() <<"   " << hQCD[3]->Integral() <<  endl;
    //   hQCD[2] ->Scale(NormFactorISO);
    // step 3 : isolation fake-rate from step 2 is aplied to QCD[0]
    hQCD[0] = (TH1D*) hQCD[2]->Clone();
    hQCD[0] ->Scale(NormFactorISO) ;
    }
    // ---- save to file
    cout << " now save results to file of " << variable << endl; 
    fOut ->cd();
    hQCD[0] ->Write();

    /*
    cout << " now save results to file " << endl; 
    hData[0]->Write();
    cout << " now save results to file " << endl; 
    hSignal[0]->Write();
    cout << " now save results to file " << endl; 
    hBack[0]->Write();
    cout << " now save results to file " << endl; 

    for ( int i = 0  ; i < NQCD ; i++){
        hQCD[i] -> Write();
    }
    */
}

vector<string> getVectorOfHistoNames(){
    //TFile* file =  new TFile("HistoFiles/SMu_8TeV_Data_dR_5311_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth.root", "READ");
    unsigned short nHist = fData[0]->GetListOfKeys()->GetEntries();


    //vector<string> histoName(nHist);
    vector<string> histoName;

    int countHist = 0 ;
    for (unsigned short i(0); i < nHist; i++) {
        string histoNameTemp = fData[0]->GetListOfKeys()->At(i)->GetName();
        TH1D* histTemp = (TH1D*) fData[0]->Get(histoNameTemp.c_str());
        //if ( histoNameTemp.find("ZMass_low") == -1 ) continue;
        if (histTemp->GetEntries() < 1) continue;
        //if (!histTemp->InheritsFrom(TH1D::Class())) cout << histoNameTemp << endl; 
        //if (!histTemp->InheritsFrom(TH1D::Class())) continue; 
        if ( histTemp->InheritsFrom(TH1D::Class())) {
                //histoName[countHist] = file->GetListOfKeys()->At(i)->GetName();
                histoName.push_back(fData[0]->GetListOfKeys()->At(i)->GetName());
            countHist++;
        }
    }
    cout << " we will produce "<<countHist << " QCD histograms  "<<  histoName.size() << endl;
//    file->Close();
    return histoName;

}
