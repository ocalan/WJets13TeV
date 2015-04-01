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

const int NQCD = 4 ;
TFile *fData[NQCD], *fSignal[NQCD];
TH1D  *hData[NQCD], *hSignal[NQCD], *hBack[NQCD];
TFile *fMC[NQCD][NFILESWJETS];

bool isMuon(0);
string unfAlg = "Bayes";
string energy = getEnergy();
int JetPtMin(30);
int JetPtMax(0);
TFile *fOut;

//-------------------------------------------------------------------------------------------
void DataDrivenQCD(  string leptonFlavor, int METcut , int doBJets ){
    
    TH1::SetDefaultSumw2();
    
    FuncOpenAllFiles(leptonFlavor, METcut, false, true, doBJets);
    vector<string> histoNameRun = getVectorOfHistoNames();
    
    for (int i(0); i < histoNameRun.size() ; i++){
        cout << endl; cout << endl;
        cout << " --- processing histogram: " << i << " : " << histoNameRun[i] << endl;
        FuncDataDrivenQCD(histoNameRun[i]);
    }
    
    //-- Close all the files ------------------------------
    cout << "I m closing the files" << endl;
    for (int i(0); i < NQCD; i++) {
        closeFile(fData[i]);
        closeFile(fSignal[i]);
        for (int j(0); j < NFILESWJETS; j++){
            closeFile(fMC[i][j]);
        }
    }
    fOut ->Close();
    cout << "QCD calculation is done. All files are closed." << endl;
    //-----------------------------------------------------
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void FuncOpenAllFiles(string leptonFlavor,int METcut, bool doFlat , bool doVarWidth, int doBJets){
    // Get data files
    for ( int i = 0 ; i < NQCD ; i++){
        fData[i] = getFile(FILESDIRECTORY,  leptonFlavor, energy, ProcessInfo[DATAFILENAME].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, i, 0, 0, METcut, doBJets, "", "0");
        if ( i == 0 ) {
            string nameQCDout = fData[i]->GetName();
            nameQCDout.insert(nameQCDout.find("Data") + 4,"QCD");
            fOut = new TFile(nameQCDout.c_str(),"recreate");
        }
    }
    /// get MC files
    for ( int i=0 ; i < NQCD ; i++){
        cout << endl;
        for ( int j = 1 ; j < NFILESWJETS ; j++){
            int sel = j ;
            if ( j == 1 ) sel = 24 ;
            cout << endl;
            fMC[i][j] = getFile(FILESDIRECTORY,  leptonFlavor, energy, ProcessInfo[sel].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, i , 0, 0, METcut, doBJets, "", "0");
            TH1D *hTemp = getHisto(fMC[i][j], "ZNGoodJets_Zexc");
            cout << " going to fetch " << ProcessInfo[sel].filename << "   " << hTemp ->Integral()<<endl;
        }
    }
    cout << endl;
    cout << "--- opened all data and MC files ---"<< endl;
}

vector<string> getVectorOfHistoNames(){
    
    unsigned short nHist = fData[0]->GetListOfKeys()->GetEntries();
    vector<string> histoName;
    
    int countHist = 0 ;
    for (unsigned short i(0); i < nHist; i++) {
        string histoNameTemp = fData[0]->GetListOfKeys()->At(i)->GetName();
        TH1D* histTemp = (TH1D*) fData[0]->Get(histoNameTemp.c_str());
        //if ( histoNameTemp.find("ZMass_low") == -1 ) continue;
        if (histTemp->GetEntries() < 1) continue;
        if ( histTemp->InheritsFrom(TH1D::Class())) {
            histoName.push_back(fData[0]->GetListOfKeys()->At(i)->GetName());
            countHist++;
        }
    }
    cout << endl;
    cout << "--- checking the histogram names to be studied ---" << endl;
    cout << " we will produce "<< countHist << " QCD histograms  " <<  histoName.size() << endl;
    cout << "--------------------------------------------------" << endl;
    return histoName;
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FuncDataDrivenQCD(string variable){
    cout << " test" << endl;
    cout << " --- now opening histograms from Data files ---" << endl;
    //--- fetch the data histograms -------------
    for ( int i=0 ; i < NQCD ; i++){
        
        TH1D *hTemp = getHisto(fData[i], variable);
        hData[i] = (TH1D *) hTemp->Clone();
        cout << " got data of histo name: " << variable << " for QCD: " << i << " " << ProcessInfo[DATAFILENAME].filename << endl;
        cout << "   integral: " << hData[i]->Integral() << endl;
    }
    cout << " --- got Data histos ---" << endl;
    
    //--- fetch the MC histograms --------------
    cout << endl;
    cout << " --- now opening histograms from all the MC files ---" << endl;
    for ( int i=0 ; i < NQCD ; i++){
        for ( int j = 1 ; j < NFILESWJETS ; j++){
            int sel = j ;
            if ( j == 1 ) sel = 24 ;
            
            cout << endl;
            TH1D *hTemp = getHisto(fMC[i][j], variable);
            cout << " going to fetch " << ProcessInfo[sel].filename << "   integral: " << hTemp ->Integral()<<endl;
            if ( ProcessInfo[sel].filename.find("WJetsALL") != -1 ) {
                cout << " This is signal: " << ProcessInfo[sel].filename << endl;
                hSignal[i] = (TH1D *) hTemp->Clone();
            }
            else{
                if ( j == 2 ) {
                    cout << " This is background ..." << endl;
                    hBack[i] = (TH1D *) hTemp->Clone();
                }
                else {
                    cout << " adding... " << endl;
                    hBack[i]->Add(hTemp);
                }
            }
            cout << " got MC histo " << "  QCD: " << i << "   MC type: " << j << "   integral: " << hTemp->Integral() << endl;
        }
        
    }
    cout << " --- got MC histos ---" << endl;
    
    //--- QCD estiamtion procedure -------------
    cout << endl;
    cout << " --- now do QCD estimation --- " << endl;
    
    // followed as descirbed in AN2012_331 pages 10,11
    TH1D* scaledMC[NQCD], *hQCD[NQCD];
    if (variable.find("ZNGoodJets") == string::npos){
        // step 0 : initial normalization of Wjets and data
        double NormFactor = (hData[0]->Integral() - hBack[0]->Integral() ) / hSignal[0]->Integral() ;
        cout << " normalization of signal " << NormFactor << "   " << hData[0]->Integral() << "   " << hSignal[0]->Integral() << "   " <<  hBack[0]->Integral() << endl;
        for ( int j = 0  ; j < 12 ; j++){
            if ( j > 0 )  NormFactor = (hData[0]->Integral() - hBack[0]->Integral() - hQCD[0] ->Integral() ) / hSignal[0]->Integral() ;
            
            // step 1:
            for ( int i = 0  ; i < NQCD ; i++){
                scaledMC[i] = (TH1D*) hSignal[i]->Clone();
                scaledMC[i]->Scale(NormFactor) ;
                hQCD[i] = (TH1D*) hData[i]->Clone();
                hQCD[i]->Add(scaledMC[i],-1);
                hQCD[i]->Add(hBack[i],-1);
                // set the negative bins to 0
                for ( int k = 0 ; k < hData[i]->GetNbinsX()+1 ; k++ ){
                    if ( hQCD[i]->GetBinContent(k) <= 0 ) {
                        hQCD[i]->SetBinContent(k, 0. ) ;
                        hQCD[i]->SetBinError(k, 0. ) ;
                    }
                }
                double NormTemp =  hData[i]->Integral() / (scaledMC[i] ->Integral() + hBack[i]->Integral()) ;
                cout << " data vs MC : " << i << "  " << NormTemp << "   " << hData[i]->Integral() << "   " << scaledMC[i] ->Integral() + hBack[i]->Integral() <<  endl;
            }
            cout << " QCD intgral A B C and D   " << hQCD[0]->Integral() << "   " << hQCD[1]->Integral() << "   " << hQCD[2]->Integral() << "  " << hQCD[3]->Integral() <<  endl;
            
            // step 2: isoaltion fake rate
            double NormFactorISO =  hQCD[1]->Integral(); //  hQCD[3]->Integral();
            if ( hQCD[3]->Integral() > 0 ) NormFactorISO /= hQCD[3]->Integral();
            else NormFactorISO = 0;
            cout << " ratio of regions B to D   " << NormFactorISO << "   " << hQCD[1]->Integral() << "   " << hQCD[3]->Integral() <<  endl;
            
            // step 3 : isolation fake-rate from step 2 is aplied to QCD[2] to get QCD[0]
            hQCD[0] = (TH1D*) hQCD[2]->Clone();
            hQCD[0] ->Scale(NormFactorISO);
        }
    }
    else{
        // for ZNGoodJets dist.
        for ( int i = 0  ; i < NQCD ; i++){
            scaledMC[i] = (TH1D*) hSignal[i]->Clone();
            hQCD[i] = (TH1D*) hData[i]->Clone();
        }
        
        for (int m = 1; m <= hData[0]->GetNbinsX(); m++){
            cout << " -------- processing binth : " << m << endl;
            // step 0 : initial normalization of Wjets and data
            double NormFactor(1);
            if(hSignal[0]->GetBinContent(m) > 0){NormFactor = (hData[0]->GetBinContent(m) - hBack[0]->GetBinContent(m)) / hSignal[0]->GetBinContent(m) ;}
            else {NormFactor = 1.;}
            
            for ( int j = 0  ; j < 12 ; j++){
                if ( j > 0 )  {
                    if(hSignal[0]->GetBinContent(m) > 0 ){NormFactor = (hData[0]->GetBinContent(m) - hBack[0]->GetBinContent(m) - hQCD[0]->GetBinContent(m) ) / hSignal[0]->GetBinContent(m) ;}
                    else {NormFactor = 1. ;}
                }
                cout << " NormFactor: " << NormFactor << endl;
                // step 1:
                for ( int i = 0  ; i < NQCD ; i++){
                    scaledMC[i]->SetBinContent(m, NormFactor * hSignal[i]->GetBinContent(m));
                    scaledMC[i]->SetBinError(m, NormFactor * hSignal[i]->GetBinError(m));
                    
                    hQCD[i]->SetBinContent(m, hData[i]->GetBinContent(m) - (scaledMC[i]->GetBinContent(m) + hBack[i]->GetBinContent(m)) );
                    hQCD[i]->SetBinError(m,  sqrt(pow(hData[i]->GetBinError(m), 2) + pow(scaledMC[i]->GetBinError(m), 2) + pow(hBack[i]->GetBinError(m), 2)) );
                    
                    // set the negative bins to 0
                    if ( hQCD[i]->GetBinContent(m) <= 0 ) {
                        hQCD[i]->SetBinContent(m, 0. ) ;
                        hQCD[i]->SetBinError(m, 0. ) ;
                    }
                }
                cout << " QCD intgral A B C and D   " << hQCD[0]->GetBinContent(m) << "   " << hQCD[1]->GetBinContent(m) << "   " << hQCD[2]->GetBinContent(m) << "  " << hQCD[3]->GetBinContent(m) <<  endl;
                
                // step 2: isoaltion fake rate
                double NormFactorISO(0);
                if ( hQCD[3]->GetBinContent(m) > 0 ) NormFactorISO = hQCD[1]->GetBinContent(m) / hQCD[3]->GetBinContent(m);
                else NormFactorISO = 0;
                cout << " ratio of regions B to D   " << NormFactorISO << "   " << hQCD[1]->GetBinContent(m) << "   " << hQCD[3]->GetBinContent(m) <<  endl;
                
                // step 3 : isolation fake-rate from step 2 is aplied to QCD[2] to get QCD[0]
                hQCD[0]->SetBinContent(m, NormFactorISO * hQCD[2]->GetBinContent(m));
                hQCD[0]->SetBinError(m, NormFactorISO * hQCD[2]->GetBinError(m));
                
                // Note: the only quantity that keeps updating is NormFactor and we treat it as constant without uncertainty. So, the uncertainty in QCD[i] is determined by the operation in the last loop, although we do calculate it in every loop. For QCD[0], the uncertainty is further scaled by NormFactorISO.
            }
            
        }
    }
    
    // ---- save to file
    cout << " --- now save results to file of: " << variable << endl;
    fOut ->cd();
    hQCD[0] ->Write();
    
}


