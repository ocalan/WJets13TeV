
#ifndef __xsecVariableNames__
#define __xsecVariableNames__

struct xsecVariableStruct{
    string ZJetsAndDPSName;
    string RivetName;
    string BlackHatName;
};

const xsecVariableStruct xsecVarNames[] = {
    // ZJetsAndDPSName             RivetName                  BlackHatName

    {"ZNGoodJets_Zexc",             "njetWJet_excl",            "njet_WMuNu"}, //0
    
    {"FirstJetPt_Zinc1jet",         "addjet_Pt1jetcase",        "ptjet1"},     //1
    {"SecondJetPt_Zinc2jet",        "addjet_Pt2jetcase",        "ptjet2"},
    {"ThirdJetPt_Zinc3jet",         "addjet_Pt3jetcase",        "ptjet3"},
    {"FourthJetPt_Zinc4jet",        "addjet_Pt4jetcase",        "ptjet4"},
    {"FifthJetPt_Zinc5jet",         "addjet_Pt5jetcase",        "NA"},         //5
    
    
    {"FirstJetPt_1_Zinc1jet",       "jet_Pt1jetcase",           "NA"},
    {"SecondJetPt_1_Zinc2jet",      "jet_Pt2jetcase",           "NA"},
    {"ThirdJetPt_1_Zinc3jet",       "jet_Pt3jetcase",           "NA"},
    {"FourthJetPt_1_Zinc4jet",      "jet_Pt4jetcase",           "NA"},
    {"FifthJetPt_1_Zinc5jet",       "jet_Pt5jetcase",           "NA"},  //10
    
    //-- jets HT distribution
    {"JetsHT_Zinc1jet",             "addJetsHT_inc1jet",        "htjet1"},
    {"JetsHT_Zinc2jet",             "addJetsHT_inc2jet",        "htjet2"},
    {"JetsHT_Zinc3jet",             "addJetsHT_inc3jet",        "htjet3"},
    {"JetsHT_Zinc4jet",             "addJetsHT_inc4jet",        "htjet4"},
    {"JetsHT_Zinc5jet",             "addJetsHT_inc5jet",        "NA"},  //15
    
    
    {"JetsHT_1_Zinc1jet",           "JetsHT_inc1jet",           "NA"},
    {"JetsHT_1_Zinc2jet",           "JetsHT_inc2jet",           "NA"},
    {"JetsHT_1_Zinc3jet",           "JetsHT_inc3jet",           "NA"},
    {"JetsHT_1_Zinc4jet",           "JetsHT_inc4jet",           "NA"},
    {"JetsHT_1_Zinc5jet",           "JetsHT_inc5jet",           "NA"},  //20
    
    //-- jet Eta distribution
    {"FirstJetEta_Zinc1jet",        "jet_eta1jetcase",      "etajet1"},
    {"SecondJetEta_Zinc2jet",       "jet_eta2jetcase",      "etajet2"},
    {"ThirdJetEta_Zinc3jet",        "jet_eta3jetcase",      "etajet3"},
    {"FourthJetEta_Zinc4jet",       "jet_eta4jetcase",      "etajet4"},
    {"FifthJetEta_Zinc5jet",        "jet_eta5jetcase",      "NA"},      //25
    
    
    // ZJetsAndDPSName                          RivetName        BlackHatName
    
    //-- angular variables
    {"dRapidityJets_Zinc2jet",                  "NA",           "NA"},
    {"dRapidityJets_Zinc3jet",                  "NA",           "NA"},
    {"dRapidityJets_Zinc4jet",                  "NA",           "NA"}, //28
    
    {"dRapidityJetsFB_Zinc2jet",                "NA",           "NA"},
    {"dRapidityJetsFB_Zinc3jet",                "NA",           "NA"},
    {"dRapidityJetsFB_Zinc4jet",                "NA",           "NA"}, //31

    {"dRapidityJets_First_Third_Zinc3jet",      "NA",           "NA"},
    {"dRapidityJets_Second_Third_Zinc3jet",     "NA",           "NA"},
    
    {"dPhiJets_Zinc2jet",                       "NA",           "NA"},
    {"dPhiJetsFB_Zinc2jet",                     "NA",           "NA"},
    {"dRJets_Zinc2jet",                         "NA",           "NA"}, //36
    
    //-- dijet properties
    {"diJetMass_Zinc2jet",                      "NA",           "NA"},
    {"diJetMass_Zinc3jet",                      "NA",           "NA"},
    {"diJetMass_Zinc4jet",                      "NA",           "NA"},
    
    {"diJetPt_Zinc2jet",                        "NA",           "NA"},
    {"diJetPt_Zinc3jet",                        "NA",           "NA"},
    {"diJetPt_Zinc4jet",                        "NA",           "NA"}, //42
    
    //-- Mean number of Jets
    {"MeanNJetsHT_Zinc1jet",                    "NA",           "NA"},
    {"MeanNJetsHT_Zinc2jet",                    "NA",           "NA"},
    {"MeanNJetsdRapidity_Zinc2jet",             "NA",           "NA"},
    {"MeanNJetsdRapidityFB_Zinc2jet",           "NA",           "NA"}, //46

};

#endif
