#ifndef __variablesOfInterest__
#define __variablesOfInterest__

#include <TObject.h>
#include <TObjString.h>
#include <TList.h>

struct variableStruct{
   string name;
   bool log, decrease;  
   int MuSVDkterm, MuBayeskterm, ESVDkterm, EBayeskterm, BayesktermCombine;
};

//--- Z Jets ------------------------------------

const int NVAROFINTERESTZJETS(62);
variableStruct VAROFINTERESTZJETS[] = {
   //--  Name  ------------------ log - decrease - Mu SVD kterm - Mu Bayes kterm - E SVD kterm - E Bayes kterm  - E Bayes kterm COMBINE
   {"ZNGoodJets_Zexc",             1,      1,           6,              4,             3,            3 ,               3    },

   //-- jet Pt distribution
   {"FirstJetPt_Zinc1jet",         1,      1,           9,              5,             9,            9 ,               6    }, 
   {"FirstJetPt_2_Zinc1jet",       1,      1,          14,              6,            14,           10 ,               7    }, 
   {"SecondJetPt_Zinc2jet",        1,      1,          11,              8,            11,            7 ,               7    }, 
   {"SecondJetPt_2_Zinc2jet",      1,      1,          15,              8,            14,            8 ,               8    }, 
   {"ThirdJetPt_Zinc3jet",         1,      1,           8,              9,             8,            7 ,               7    }, 
   {"FourthJetPt_Zinc4jet",        1,      1,           6,              4,             6,            5 ,               4    }, 
   {"FifthJetPt_Zinc5jet",         1,      1,           4,              2,             4,            3 ,               2    },

   {"FirstHighestJetPt_Zinc1jet",  1,      1,           7,              4,             8,            3 ,               7    },
   {"FirstHighestJetPt_Zinc2jet",  1,      1,           9,              4,             8,            4 ,               7    },
   {"FirstHighestJetPt_Zinc3jet",  1,      1,           8,              4,             8,            3 ,               7    },
   {"FirstHighestJetPt_Zinc4jet",  1,      0,           5,              2,             7,            5 ,               7    },
   {"FirstHighestJetPt_Zinc5jet",  1,      0,           9,              3,             5,            5 ,               7    },

   {"SecondHighestJetPt_Zinc2jet", 1,      1,           9,              3,            11,            3 ,               7    },
   {"SecondHighestJetPt_Zinc3jet", 1,      1,          10,              3,            10,            4 ,               7    },
   {"SecondHighestJetPt_Zinc4jet", 1,      1,           9,              3,             9,            3 ,               7    },
   {"SecondHighestJetPt_Zinc5jet", 1,      0,           4,              2,             6,            4 ,               7    },

   {"ThirdHighestJetPt_Zinc3jet",  1,      1,           6,              2,             6,            2 ,               7    },
   {"ThirdHighestJetPt_Zinc4jet",  1,      1,           6,              3,             5,            2 ,               7    },
   {"ThirdHighestJetPt_Zinc5jet",  1,      1,           5,              3,             3,            3 ,               7    },

   //-- jets HT distribution
   {"JetsHT_Zinc1jet",             1,      1,           7,              4,             7,            3 ,               3    }, 
   {"JetsHT_2_Zinc1jet",           1,      1,           9,              4,             8,            4 ,               4    }, 
   {"JetsHT_Zinc2jet",             1,      1,           7,              3,             9,            4 ,               2    },
   {"JetsHT_2_Zinc2jet",           1,      1,           7,              3,             8,            4 ,               2    },
   {"JetsHT_Zinc3jet",             1,      1,           6,              3,             6,            6 ,               3    },
   {"JetsHT_Zinc4jet",             1,      1,           5,              5,             4,            4 ,               3    },
   {"JetsHT_Zinc5jet",             1,      1,           3,              2,             4,            3 ,               1    },/// ?????

   //-- jet Eta distribution
   {"FirstJetEta_Zinc1jet",        0,      1,          28,              4,            18,            4 ,               5    },
   {"SecondJetEta_Zinc2jet",       0,      1,          11,              5,            16,            5 ,               5    }, 
   {"ThirdJetEta_Zinc3jet",        0,      1,          10,              4,             8,            4 ,               5    }, 
   {"FourthJetEta_Zinc4jet",       0,      1,           4,              4,             4,            5 ,               4    }, 
   {"FifthJetEta_Zinc5jet",        0,      1,           3,              2,             2,            2 ,               1    }, // ???

   //-- Z Mass distribution
   //{"ZMass_Zinc0jet",              0,      1,          19,             17,             7,          6 ,               7    },

   //-- Jets Mass distribution
   {"JetsMass_Zexc2jet",           0,      1,           7,              5,            13,            6 ,               4    },
   {"JetsMass_Zinc2jet",           0,      1,           9,              5,            13,            7 ,               5    },

   //-- DPS distribution
   {"SpTJets_Zexc2jet",            0,      0,           7,              4,             6,            4 ,               3    },
   {"SpTJets_Zinc2jet",            0,      0,           7,              4,             8,            4 ,               3    },
   {"SpT_Zexc2jet",                0,      0,          14,              3,             8,            3 ,               2    },
   {"SpT_Zinc2jet",                0,      0,          16,              3,             5,            3 ,               2    },
   {"dPhiJets_Zexc2jet",           0,      0,           4,              4,             4,            4 ,               3    },
   {"dPhiJets_Zinc2jet",           0,      0,           6,              4,             8,            4 ,               3    },
   {"PHI_Zexc2jet",                0,      0,           5,              5,             6,            4 ,               3    },
   {"PHI_Zinc2jet",                0,      0,          14,              4,            10,            4 ,               3    },
   {"PHI_T_Zexc2jet",              0,      0,           6,              6,             4,            6 ,               5    },
   {"PHI_T_Zinc2jet",              0,      0,           6,              6,             7,            5 ,               5    },
   {"SPhi_Zexc2jet",               0,      0,          18,              4,             9,            4 ,               3    },
   {"SPhi_Zinc2jet",               0,      0,           7,              4,             9,            4 ,               3    },

   {"SpTJets_LowPt_Zexc2jet",      0,      1,           9,             12,             7,            9 ,               5    },
   {"SpTJets_LowPt_Zinc2jet",      0,      1,           8,             12,             9,            8 ,               6    },
   {"SpT_LowPt_Zexc2jet",          0,      1,           6,             13,             4,            9 ,               7    },
   {"SpT_LowPt_Zinc2jet",          0,      1,           6,             12,            14,            9 ,               7    },
   {"dPhiJets_LowPt_Zexc2jet",     0,      0,           6,             13,             4,            9 ,               5    },
   {"dPhiJets_LowPt_Zinc2jet",     0,      0,           4,             12,             8,            9 ,               6    },
   {"PHI_LowPt_Zexc2jet",          0,      1,           7,              5,             4,            8 ,               3    },
   {"PHI_LowPt_Zinc2jet",          0,      1,           5,              5,             2,           15 ,               3    },
   {"PHI_T_LowPt_Zexc2jet",        0,      0,           3,              5,             7,            6 ,               4    },
   {"PHI_T_LowPt_Zinc2jet",        0,      0,           3,              4,             3,            6 ,               4    },
   {"SPhi_LowPt_Zexc2jet",         0,      0,          11,              7,             3,            6 ,               7    },
   {"SPhi_LowPt_Zinc2jet",         0,      0,           8,              7,             6,            6 ,               7    },
   // additional correlations  =----- NEED TP CHECK THE UNFOLDING PARAMTERS !!!!
   {"dEtaJets_Zinc2jet",           0,      0,           8,             10,             6,            6 ,               4    },
   {"dEtaJets_Zexc2jet",           0,      0,           8,             10,             6,            6 ,               3    },
   {"dEtaBosonJet_Zinc1",          0,      0,           8,             10,             6,            6 ,               8    },
   {"dEtaBosonJet_Zexc1",          0,      0,           8,             10,             6,            6 ,               8    },
//   {"dEtaFirstJetZ_Zinc2jet",      0,      0,           8,             10,             6,            6 ,               7    },
//   {"dEtaSecondJetZ_Zinc2jet",     0,      0,           8,             10,             6,            6 ,               7    }

}; 

TList *LISTOFVAROFINTERESTZJETS = new TList();

void makeLISTOFVAROFINTERESTZJETS(){
    for (int i(0); i < NVAROFINTERESTZJETS; i++){
        TObjString *name = new TObjString(VAROFINTERESTZJETS[i].name.c_str());
        LISTOFVAROFINTERESTZJETS->Add(name);
    }
}

//--- W Jets ------------------------------------

const int NVAROFINTERESTWJETS(66);
variableStruct VAROFINTERESTWJETS[] = {
    //--  Name  ------------------ log - decrease - Mu SVD kterm - Mu Bayes kterm - E SVD kterm - E Bayes kterm
    {"ZNGoodJets_Zexc",             1,      1,           4,              4,             4,            3,           7},

    //-- jet Pt distribution // I add these later: without carefully study the optimum kterm
    {"FirstJetPt_Zinc1jet",         1,      1,          14,              6,            14,           10,           7},
    {"FirstJetPt_1_Zinc1jet",       1,      1,          14,              6,            14,           10 ,          7},
    {"FirstJetPt_2_Zinc1jet",       1,      1,          14,              6,            14,           10 ,          7},
    {"SecondJetPt_Zinc2jet",        1,      1,          15,              8,            14,            8,           8},
    {"SecondJetPt_1_Zinc2jet",      1,      1,          15,              8,            14,            8 ,          8},
    {"ThirdJetPt_Zinc3jet",         1,      1,           5,              4,             5,            5,           7},
    {"ThirdJetPt_1_Zinc3jet",       1,      1,           5,              4,             5,            5,           7},
    {"FourthJetPt_Zinc4jet",        1,      1,           5,              3,             4,            4,           7},
    {"FourthJetPt_1_Zinc4jet",      1,      1,           5,              3,             4,            4,           7},
    {"FifthJetPt_Zinc5jet",         1,      1,           4,              3,             3,            2,           7},
    {"FifthJetPt_1_Zinc5jet",       1,      1,           4,              3,             3,            2,           7},

    {"FirstHighestJetPt_Zinc1jet",  1,      1,           9,              4,             9,            3,           7},
    {"FirstHighestJetPt_Zinc2jet",  1,      1,           9,              7,             9,            4,           7},
    {"FirstHighestJetPt_Zinc3jet",  1,      1,           3,              4,             5,            4,           7},
    {"FirstHighestJetPt_Zinc4jet",  1,      0,           3,              5,             3,            5,           7},
    {"FirstHighestJetPt_Zinc5jet",  1,      0,           3,              2,             2,            2,           7},

    {"SecondHighestJetPt_Zinc2jet", 1,      1,           8,              3,             8,            3,           7},
    {"SecondHighestJetPt_Zinc3jet", 1,      1,           5,              4,             6,            4,           7},
    {"SecondHighestJetPt_Zinc4jet", 1,      1,           3,              4,             3,            5,           7},
    {"SecondHighestJetPt_Zinc5jet", 1,      0,           3,              2,             2,            4,           7},

    {"ThirdHighestJetPt_Zinc3jet",  1,      1,           5,              2,             5,            4,           7},
    {"ThirdHighestJetPt_Zinc4jet",  1,      1,           4,              3,             5,            2,           7},
    {"ThirdHighestJetPt_Zinc5jet",  1,      1,           3,              3,             3,            3,           7},

    //-- jets HT distribution
    {"JetsHT_Zinc1jet",             1,      1,           9,              4,             8,            4,           4},
    {"JetsHT_1_Zinc1jet",           1,      1,           9,              4,             8,            4 ,          4},
    
    {"JetsHT_Zinc2jet",             1,      1,           7,              4,             8,            4,           2},
    {"JetsHT_1_Zinc2jet",           1,      1,           7,              4,             8,            4 ,          2},
    
    {"JetsHT_Zinc3jet",             1,      1,           7,              3,             7,            5,           7},
    {"JetsHT_1_Zinc3jet",           1,      1,           7,              3,             7,            5,           7},
    
    {"JetsHT_Zinc4jet",             1,      1,           6,              3,             5,            9,           7},
    {"JetsHT_1_Zinc4jet",           1,      1,           6,              3,             5,            9,           7},
    
    {"JetsHT_Zinc5jet",             1,      1,           3,              2,             4,            2,           7},
    {"JetsHT_1_Zinc5jet",           1,      1,           3,              2,             4,            2,           7},

    //-- jet Eta distribution
    {"FirstJetEta_Zinc1jet",        0,      1,          13,              3,            15,            3,           7},
    {"SecondJetEta_Zinc2jet",       0,      1,          11,              5,            11,            4,           7}, 
    {"ThirdJetEta_Zinc3jet",        0,      1,          10,              4,             8,            4,           7}, 
    {"FourthJetEta_Zinc4jet",       0,      1,           5,              4,             3,            4,           7}, 
    {"FifthJetEta_Zinc5jet",        0,      1,           3,              2,             2,            2,           7}, 

    //-- Z Mass distribution
    {"ZMass_Zinc0jet",              0,      1,          19,              4,             7,           36,           7},

    //-- Jets Mass distribution
    {"JetsMass_Zexc2jet",           0,      1,           9,              6,            10,            7,           7},
    {"JetsMass_Zinc2jet",           0,      1,           9,              7,             9,            6,           7},

    //-- DPS distribution
    {"SpTJets_Zexc2jet",            0,      0,           7,              4,             5,            4,           7},
    {"SpTJets_Zinc2jet",            0,      0,           7,              5,             5,            4,           7},
    {"SpT_Zexc2jet",                0,      0,          10,              4,             6,            4,           7},
    {"SpT_Zinc2jet",                0,      0,          11,              5,             6,            3,           7},
    {"dPhiJets_Zexc2jet",           0,      0,           7,              4,             5,            3,           3},
    {"dPhiJets_Zinc2jet",           0,      0,           7,              4,             8,            3,           3},
    {"PHI_Zexc2jet",                0,      0,           5,              5,             6,            6,           7},
    {"PHI_Zinc2jet",                0,      0,          14,              6,            10,            5,           7},
    {"PHI_T_Zexc2jet",              0,      0,           5,              7,             5,            6,           7},
    {"PHI_T_Zinc2jet",              0,      0,           5,              7,             4,            6,           7},
    {"SPhi_Zexc2jet",               0,      0,          18,              5,             9,            6,           7},
    {"SPhi_Zinc2jet",               0,      0,           7,              5,             9,            5,           7},

    {"SpTJets_LowPt_Zexc2jet",      0,      1,           4,              7,             6,           11,           7},
    {"SpTJets_LowPt_Zinc2jet",      0,      1,           4,              7,             9,            8,           7},
    {"SpT_LowPt_Zexc2jet",          0,      1,           4,              7,             3,            5,           7},
    {"SpT_LowPt_Zinc2jet",          0,      1,           6,              8,             8,            5,           7},
    {"dPhiJets_LowPt_Zexc2jet",     0,      0,           4,              6,             6,            7,           7},
    {"dPhiJets_LowPt_Zinc2jet",     0,      0,           4,              7,             2,            6,           7},
    {"PHI_LowPt_Zexc2jet",          0,      1,           7,              6,             4,           14,           7},
    {"PHI_LowPt_Zinc2jet",          0,      1,           5,              5,             2,           15,           7},
    {"PHI_T_LowPt_Zexc2jet",        0,      0,           3,              5,             7,            6,           7},
    {"PHI_T_LowPt_Zinc2jet",        0,      0,           3,              4,             3,            6,           7},
    {"SPhi_LowPt_Zexc2jet",         0,      0,          11,              7,             3,            6,           7},
    {"SPhi_LowPt_Zinc2jet",         0,      0,           8,              7,             6,            6,           7},
}; 

TList *LISTOFVAROFINTERESTWJETS = new TList();

void makeLISTOFVAROFINTERESTWJETS(){
    for (int i(0); i < NVAROFINTERESTWJETS; i++){
        TObjString *name = new TObjString(VAROFINTERESTWJETS[i].name.c_str());
        LISTOFVAROFINTERESTWJETS->Add(name);
    }
}
#endif




