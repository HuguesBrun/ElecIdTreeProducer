// -*- C++ -*-
//
// Package:    hugues/ElecIdTreeProducer
// Class:      ElecIdTreeProducer
//
/**\class ElecIdTreeProducer ElecIdTreeProducer.cc hugues/ElecIdTreeProducer/plugins/ElecIdTreeProducer.cc
 
 Description: [one line class summary]
 
 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Hugues Louis Brun
//         Created:  Fri, 30 May 2014 13:29:42 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"


// root stuff !
#include "TH1D.h"
#include <map>
#include "TFile.h"
#include <math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TString.h"
#include "TTree.h"



//
// class declaration
//

class ElecIdTreeProducer : public edm::EDAnalyzer {
public:
    explicit ElecIdTreeProducer(const edm::ParameterSet&);
    ~ElecIdTreeProducer();
    
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    
private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    virtual void beginEvent();
    virtual void endEvent();
    
    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    
    
    // ----------member data ---------------------------
    bool isMC_;
    edm::InputTag electronsCollection_;
    std::string outputFile_; // output file
    
    
    // ---------- output ROOT file
    TFile*  rootFile_;
    
    // ---------- tree declaration
    TTree *mytree_;
    
    // -----------tree variables
    //Events
    int T_Event_RunNumber;
    int T_Event_EventNumber;
    int T_Event_LuminosityBlock;
    // MC info if need
    int T_Event_processID;
    int T_Event_ptHat;
    // all PU info
    int T_Event_nPU;
    float T_Event_nTruePU;
    int T_Event_nPUm;
    int T_Event_nPUp;
    float T_Event_AveNTruePU;
    float T_Event_Rho;

    
    // electron variables
    // kinematics
    std::vector<float> *T_Elec_Eta;
    std::vector<float> *T_Elec_Phi;
    std::vector<float> *T_Elec_Px;
    std::vector<float> *T_Elec_Py;
    std::vector<float> *T_Elec_Pz;
    std::vector<float> *T_Elec_Pt;
    std::vector<float> *T_Elec_Energy;
    std::vector<int> *T_Elec_Charge;

    //position in ECAL
    std::vector<int> *T_Elec_isEB;
    std::vector<int> *T_Elec_isEE;

    
    // tracks
    std::vector<float> *T_Elec_vz;
    std::vector<float> *T_Elec_vy;
    std::vector<float> *T_Elec_vx;
    std::vector<int>   *T_Elec_nLost;
    std::vector<int>   *T_Elec_nHits;
    std::vector<float> *T_Elec_kfchi2;
    std::vector<int>   *T_Elec_kfhits;
    std::vector<int>   *T_Elec_gsfhits;
    std::vector<float> *T_Elec_gsfchi2;
    std::vector<float> *T_Elec_fbrem;
    std::vector<int>   *T_Elec_nbrems;
    std::vector<int>   *T_Elec_missingHits;
    std::vector<float> *T_Elec_Dist;
    std::vector<float> *T_Elec_Dcot;
    std::vector<float> *T_Elec_D0;
    std::vector<float> *T_Elec_Dz;
    std::vector<float> *T_Elec_ip3d;
    std::vector<float> *T_Elec_ip3ds;

    
    // tracks calo matching
    std::vector<float> *T_Elec_detacalo;
    std::vector<float> *T_Elec_eledeta;
    std::vector<float> *T_Elec_dphicalo;
    std::vector<float> *T_Elec_deltaPhiIn;
    std::vector<float> *T_Elec_deltaEtaIn;
    std::vector<float> *T_Elec_EoP;
    std::vector<float> *T_Elec_EoPin;
    std::vector<float> *T_Elec_ESeedoP;
    std::vector<float> *T_Elec_ESeedoPout;
    std::vector<float> *T_Elec_EEleoPout;
    std::vector<float> *T_Elec_IoEmIoP;
    std::vector<float> *T_Elec_eleEoPout;

    
    
    //SC
    std::vector<float> *T_Elec_SC_Et;
    std::vector<float> *T_Elec_SC_Eta;
    std::vector<float> *T_Elec_SC_Phi;
    std::vector<float> *T_Elec_SC_RawEnergy;
    std::vector<float> *T_Elec_EcalEnergy;
    std::vector<float> *T_Elec_EsEnergy;
    std::vector<float> *T_Elec_PreShowerOverRaw;
    std::vector<int>   *T_Elec_NClusters;

    
    
    //seed Basic clusters
    std::vector<float> *T_Elec_EtaSeed;
    std::vector<float> *T_Elec_PhiSeed;
    std::vector<float> *T_Elec_ESeed;
    std::vector<int> *T_Elec_IEta;
    std::vector<int> *T_Elec_IPhi;

    
    //seed Xtal
    std::vector<float> *T_Elec_EtaSeedXtal;
    std::vector<float> *T_Elec_PhiSeedXtal;
    std::vector<float> *T_Elec_IEtaSeedXtal;
    std::vector<float> *T_Elec_IPhiSeedXtal;
    
    // SC shape
    std::vector<float> *T_Elec_HtoE;
    std::vector<float> *T_Elec_EmaxSeed;
    std::vector<float> *T_Elec_EtopSeed;
    std::vector<float> *T_Elec_EbottomSeed;
    std::vector<float> *T_Elec_EleftSeed;
    std::vector<float> *T_Elec_ErightSeed;
    std::vector<float> *T_Elec_E2ndSeed;
    
    std::vector<float> *T_Elec_E2x5RightSeed;
    std::vector<float> *T_Elec_E2x5LeftSeed;
    std::vector<float> *T_Elec_E2x5TopSeed;
    std::vector<float> *T_Elec_E2x5BottomSeed;
    std::vector<float> *T_Elec_E2x5MaxSeed;
    
    std::vector<float> *T_Elec_E1x5Seed;
    std::vector<float> *T_Elec_E2x2Seed;
    std::vector<float> *T_Elec_E3x3Seed;
    std::vector<float> *T_Elec_E5x5Seed;
    
    
    std::vector<float> *T_Elec_see;
    std::vector<float> *T_Elec_spp;
    std::vector<float> *T_Elec_sep;
    std::vector<float> *T_Elec_etawidth;
    std::vector<float> *T_Elec_phiwidth;
    std::vector<float> *T_Elec_e1x5e5x5;
    std::vector<float> *T_Elec_s9e25;
    std::vector<float> *T_Elec_R9;
    
    
    
    
    
    //IP infos
    
    //conversion rejection
    std::vector<int> *T_Elec_MatchConv;
    std::vector<int> *T_Elec_EcalDriven;

    
    
    
    
};

