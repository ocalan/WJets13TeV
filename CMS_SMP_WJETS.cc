#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, WFinder...
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/AnalysisLoader.hh"
#include <iostream>

//using namespace std;

namespace Rivet {


  class CMS_SMP_WJETS: public Analysis {
  public:

    CMS_SMP_WJETS()
      : Analysis("CMS_SMP_WJETS")
    {
      setBeams(PROTON, PROTON);
      setNeedsCrossSection(true);
    }

  public:


    // Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      WFinder wfinder_mu(fs, -2.4, 2.4, 0*GeV, MUON, 50*GeV, 100*GeV, 0*GeV, 0.1, true, false);
      addProjection(wfinder_mu, "WFinder_mu");


     // Define veto FS 
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(wfinder_mu);
      vfs.addVetoPairId(MUON);      
      vfs.vetoNeutrinos();

      FastJets fastjets(vfs, FastJets::ANTIKT, 0.5);
      addProjection(fastjets, "Jets");

      vector<double> addjetPt_Winc1jet ;

      vector<double> jetPt_Winc1jet ;
      vector<double> jetPt_Winc2jet ;
      vector<double> jetPt_Winc3jet ;
      vector<double> jetPt_Winc4jet ;
      vector<double> jetPt_Winc5jet ;

      vector<double> jetHT_Winc1jet ;
      vector<double> jetHT_Winc2jet ;
      vector<double> jetHT_Winc3jet ;
      vector<double> jetHT_Winc4jet ;
      vector<double> jetHT_Winc5jet ;

      jetPt_Winc1jet += 20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 700, 900, 1400;
      jetPt_Winc2jet += 20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 700, 1000;
      jetPt_Winc3jet += 20, 24, 30, 39, 49, 62, 78, 105, 142, 185, 235, 300, 500;
      jetPt_Winc4jet += 20, 24, 30, 39, 49, 62, 78, 96, 120, 160, 250;
      jetPt_Winc5jet += 20, 24, 30, 39, 49, 62, 90, 140;

  
      jetHT_Winc1jet += 30, 39, 49, 62, 78, 96, 118, 150, 190, 240, 300, 370, 450, 540, 650, 800, 1000, 1300, 2000;
      jetHT_Winc2jet += 60, 78, 96, 118, 150, 190, 240, 300, 370, 450, 540, 650, 800, 1000, 1300, 2000;
      jetHT_Winc3jet += 90, 118, 150, 190, 240, 300, 370, 450, 540, 650, 800, 1100, 1800;
      jetHT_Winc4jet += 120, 140, 165, 200, 250, 310, 380, 480, 600, 800, 1100, 1800;
      jetHT_Winc5jet += 150, 180, 222, 282, 365, 485, 650, 880, 1300;

      addjetPt_Winc1jet += 20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 700, 1000;

     _hist_addJetPt1j =bookHistogram1D("addjet_Pt1jetcase", addjetPt_Winc1jet);


     _hist_JetPt1j =bookHistogram1D("jet_Pt1jetcase", jetPt_Winc1jet);
     _hist_JetPt2j =bookHistogram1D("jet_Pt2jetcase", jetPt_Winc2jet);
     _hist_JetPt3j =bookHistogram1D("jet_Pt3jetcase", jetPt_Winc3jet);
     _hist_JetPt4j =bookHistogram1D("jet_Pt4jetcase", jetPt_Winc4jet);
     _hist_JetPt5j =bookHistogram1D("jet_Pt5jetcase", jetPt_Winc5jet);

     _hist_Jeteta1j =bookHistogram1D("jet_eta1jetcase", 32, 0, 2.4);
     _hist_Jeteta2j =bookHistogram1D("jet_eta2jetcase", 32, 0, 2.4);
     _hist_Jeteta3j =bookHistogram1D("jet_eta3jetcase", 24, 0, 2.4);
     _hist_Jeteta4j =bookHistogram1D("jet_eta4jetcase", 12, 0, 2.4);
     _hist_Jeteta5j =bookHistogram1D("jet_eta5jetcase", 6,  0, 2.4);

     _hist_Ht_1j =bookHistogram1D("JetsHT_inc1jet", jetHT_Winc1jet);
     _hist_Ht_2j =bookHistogram1D("JetsHT_inc2jet", jetHT_Winc2jet);
     _hist_Ht_3j =bookHistogram1D("JetsHT_inc3jet", jetHT_Winc3jet);
     _hist_Ht_4j =bookHistogram1D("JetsHT_inc4jet", jetHT_Winc4jet);
     _hist_Ht_5j =bookHistogram1D("JetsHT_inc5jet", jetHT_Winc5jet);

     _hist_inc_WJetMult = bookHistogram1D("njetWJet_incl", 8, -0.5, 7.5);
     _hist_excl_WJetMult= bookHistogram1D("njetWJet_excl", 8, -0.5, 7.5);

  }

 void Fill(AIDA::IHistogram1D*& _histJetMult, const double& weight, std::vector<FourMomentum>& finaljet_list){
       _histJetMult->fill(0, weight);
       for (size_t i=0 ; i<finaljet_list.size() ; ++i) {
         if (i==6) break;
         _histJetMult->fill(i+1, weight);  // inclusive multiplicity
   
        }
   }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
   
      const double weight = event.weight();
      const WFinder& wfinder_mu = applyProjection<WFinder>(event, "WFinder_mu");

      if (wfinder_mu.bosons().size() != 1) {
          vetoEvent;
      }

      if (wfinder_mu.bosons().size()==1) {

       const FourMomentum& boson = wfinder_mu.bosons().front().momentum();
       const FourMomentum& neutrino = wfinder_mu.constituentNeutrinos()[0].momentum();
       const FourMomentum& lepton0 = wfinder_mu.constituentLeptons()[0].momentum();
       const FourMomentum& lepton1 = wfinder_mu.constituentLeptons()[1].momentum();

       //double mT = sqrt(2.0 * lepton.pT() * p_miss.Et() * (1.0 - cos( lepton.phi()-p_miss.phi()) ) );

       double pt0 = lepton0.pT()/GeV;
       double pt1 = lepton1.pT()/GeV;
       double eta0 = lepton0.eta();
       double eta1 = lepton1.eta();

       bool inAcceptance = fabs(eta0)<2.1 && fabs(eta1)<2.1 && pt0>25 && pt1>15;

       if (!inAcceptance)return;

       //Obtain the jets::::::::::::::

       vector<FourMomentum> finaljet_list;
       double HT = 0.0;
     
      foreach (const Jet& j, applyProjection<FastJets>(event, "Jets").jetsByPt(30.0*GeV)) {
        const double jeta = j.momentum().eta();
        const double jpt = j.momentum().pT();

        if (fabs(jeta) < 2.4)
        if(jpt>30){
                finaljet_list.push_back(j.momentum());
                HT += j.momentum().pT();     
      }
 }

  //Multiplicity plots.     

     _hist_excl_WJetMult->fill(finaljet_list.size(), weight);
      Fill(_hist_inc_WJetMult, weight, finaljet_list);
      
      if(finaljet_list.size()>=1) {
            _hist_addJetPt1j->fill(finaljet_list[0].pT(), weight);

            _hist_JetPt1j->fill(finaljet_list[0].pT(), weight);
            _hist_Jeteta1j->fill(finaljet_list[0].eta(), weight);
            _hist_Ht_1j->fill(HT, weight);
      }

      if(finaljet_list.size()>=2) {
            _hist_JetPt2j->fill(finaljet_list[1].pT(), weight);
            _hist_Jeteta2j->fill(finaljet_list[1].eta(), weight);
            _hist_Ht_2j->fill(HT, weight);
      }

      if(finaljet_list.size()>=3) {
           _hist_JetPt3j->fill(finaljet_list[2].pT(), weight);
           _hist_Jeteta3j->fill(finaljet_list[2].eta(), weight);
           _hist_Ht_3j->fill(HT, weight);
      }

      if(finaljet_list.size()>=4) {
           _hist_JetPt4j->fill(finaljet_list[3].pT(), weight);
           _hist_Jeteta4j->fill(finaljet_list[3].eta(), weight);
           _hist_Ht_4j->fill(HT, weight);
     }


     if(finaljet_list.size()>=5) {
           _hist_JetPt5j->fill(finaljet_list[4].pT(), weight);
           _hist_Jeteta5j->fill(finaljet_list[4].eta(), weight);
           _hist_Ht_5j->fill(HT, weight);
     }

  }// close the muon loop

} //void loop


/// Normalise histograms etc., after the run
    void finalize() {

        double crossSec=36702;

        scale(_hist_inc_WJetMult, crossSec/sumOfWeights());
        scale(_hist_excl_WJetMult, crossSec/sumOfWeights());
 
        scale(_hist_addJetPt1j, crossSec/sumOfWeights());

        scale(_hist_JetPt1j, crossSec/sumOfWeights());
        scale(_hist_JetPt2j, crossSec/sumOfWeights());
        scale(_hist_JetPt3j, crossSec/sumOfWeights());
        scale(_hist_JetPt4j, crossSec/sumOfWeights());
        scale(_hist_JetPt5j, crossSec/sumOfWeights());

        scale(_hist_Jeteta1j, crossSec/sumOfWeights());
        scale(_hist_Jeteta2j, crossSec/sumOfWeights());
        scale(_hist_Jeteta3j, crossSec/sumOfWeights());
        scale(_hist_Jeteta4j, crossSec/sumOfWeights());
        scale(_hist_Jeteta5j, crossSec/sumOfWeights());

        scale(_hist_Ht_1j, crossSec/sumOfWeights());
        scale(_hist_Ht_2j, crossSec/sumOfWeights());
        scale(_hist_Ht_3j, crossSec/sumOfWeights());
        scale(_hist_Ht_4j, crossSec/sumOfWeights());
        scale(_hist_Ht_5j, crossSec/sumOfWeights());


}

  private:

    // Data members like post-cuts event weight counters go here

  private:

        AIDA::IHistogram1D* _hist_inc_WJetMult;
        AIDA::IHistogram1D* _hist_excl_WJetMult;

        AIDA::IHistogram1D* _hist_addJetPt1j;

        AIDA::IHistogram1D* _hist_JetPt1j;
        AIDA::IHistogram1D* _hist_JetPt2j;
        AIDA::IHistogram1D* _hist_JetPt3j;
        AIDA::IHistogram1D* _hist_JetPt4j;
        AIDA::IHistogram1D* _hist_JetPt5j;

        AIDA::IHistogram1D* _hist_Jeteta1j;
        AIDA::IHistogram1D* _hist_Jeteta2j;
        AIDA::IHistogram1D* _hist_Jeteta3j;
        AIDA::IHistogram1D* _hist_Jeteta4j;
        AIDA::IHistogram1D* _hist_Jeteta5j;
        AIDA::IHistogram1D* _hist_Jeteta6j;

        AIDA::IHistogram1D* _hist_Ht_1j;
        AIDA::IHistogram1D* _hist_Ht_2j;
        AIDA::IHistogram1D* _hist_Ht_3j;
        AIDA::IHistogram1D* _hist_Ht_4j;
        AIDA::IHistogram1D* _hist_Ht_5j;


  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CMS_SMP_WJETS> plugin_CMS_SMP_WJETS;     

}



