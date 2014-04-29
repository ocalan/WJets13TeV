#ifndef _DataDrivenQCD_h_
#define _DataDrivenQCD_h_

#include <iostream>

using namespace std;

void DataDrivenQCD( string leptonFlavor = "SMu", int METcut = 0,  int doBJets = 0);
vector<string> getVectorOfHistoNames();
void FuncDataDrivenQCD(string variable="ZNGoodJets_Zexc" , string leptonFlavor = "SMu", int METcut = 0, bool doFlat = false, bool doVarWidth = true, int doBJets = 0 );
void FuncOpenAllFiles(string leptonFlavor = "SMu",  int METcut = 0, bool doFlat = false, bool doVarWidth = true , int doBJets = 0 );
#endif

