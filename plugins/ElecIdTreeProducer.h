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

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"


#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "EgammaAnalysis/ElectronTools/interface/SuperClusterHelper.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"


#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

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
    virtual bool hasWZasMother(const reco::GenParticle);
    virtual bool isMatchedWithGen(reco::GenParticle, const reco::GsfElectron &);

    
    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    
    
    // ----------member data ---------------------------
    bool isGeomInitialized_;
    const CaloTopology * ecalTopology_;
    const CaloGeometry * caloGeometry_;

    HLTConfigProvider hltConfig; 
    
    bool isMC_;
    edm::InputTag electronsCollection_;
    edm::InputTag primaryVertexInputTag_;
    edm::InputTag EBRecHitsLabel_;
    edm::InputTag EERecHitsLabel_;
    edm::InputTag conversionsInputTag_;
    edm::InputTag beamSpotInputTag_;
    edm::InputTag metTag_;
    edm::InputTag jetCollectionTag_;
    edm::InputTag triggerResultsTag_;
    edm::InputTag triggerSummaryLabel_;
    std::vector<edm::InputTag> rhoInputTags_;
    std::vector<std::string> pathsToSave_;
    std::vector<std::string> filterToMatch_;
    std::string HLTprocess_;
    std::string outputFile_; // output file
    
    
    std::vector<int> triggerBits_;
    
    
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
    
    
    
    std::vector<float> * T_Event_Rho;
    
    std::vector<int> *T_Event_pathsFired;
    
    
    
    
    // electron variables
    
    // gen informations
    // gen info on the electron
    std::vector<float> *T_Gen_Elec_Px;
    std::vector<float> *T_Gen_Elec_Py;
    std::vector<float> *T_Gen_Elec_Pz;
    std::vector<float> *T_Gen_Elec_Energy;
    std::vector<int> *T_Gen_Elec_PDGid;
    std::vector<int> *T_Gen_Elec_MotherID;
    std::vector<int> *T_Gen_Elec_GndMotherID;
    
    //trigger leg
    std::vector<int> *T_Elec_TriggerLeg;

    
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
    std::vector<float> *T_Elec_dphicalo;
    std::vector<float> *T_Elec_eledeta;
    std::vector<float> *T_Elec_eledphi;
    std::vector<float> *T_Elec_deltaPhiIn;
    std::vector<float> *T_Elec_deltaEtaIn;
    std::vector<float> *T_Elec_EoP;
    std::vector<float> *T_Elec_ESeedoP;
    std::vector<float> *T_Elec_ESeedoPout;
    std::vector<float> *T_Elec_EEleoPout;
    std::vector<float> *T_Elec_IoEmIoP;
    
    
    
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
    
    
    std::vector<float> *T_Elec_noZSsee;
    std::vector<float> *T_Elec_noZSspp;
    std::vector<float> *T_Elec_noZSsep;
    std::vector<float> *T_Elec_noZSr9;
    std::vector<float> *T_Elec_noZSe1x5;
    std::vector<float> *T_Elec_noZSe2x5MaxSeed;
    std::vector<float> *T_Elec_noZSe5x5;
    
    
    
    // isolation stuff
    std::vector<float> *T_Elec_puChargedIso;
    std::vector<float> *T_Elec_allChargedHadronIso;
    std::vector<float> *T_Elec_chargedHadronIso;
    std::vector<float> *T_Elec_neutralHadronIso;
    std::vector<float> *T_Elec_photonIso;
    std::vector<float> *T_Elec_puChargedIso04;
    std::vector<float> *T_Elec_allChargedHadronIso04;
    std::vector<float> *T_Elec_chargedHadronIso04;
    std::vector<float> *T_Elec_neutralHadronIso04;
    std::vector<float> *T_Elec_photonIso04;
    
    // det based isolation stuff
    std::vector<float> *T_Elec_ECALiso;
    std::vector<float> *T_Elec_HCALiso;
    std::vector<float> *T_Elec_TKiso;
    
    
    // infos on the BCs
    std::vector<int> *T_Elec_nbBC;
    std::vector<float> *T_Elec_BC1_eta;
    std::vector<float> *T_Elec_BC1_phi;
    std::vector<float> *T_Elec_BC1_energy;
    std::vector<float> *T_Elec_BC2_eta;
    std::vector<float> *T_Elec_BC2_phi;
    std::vector<float> *T_Elec_BC2_energy;
    std::vector<float> *T_Elec_BC3_eta;
    std::vector<float> *T_Elec_BC3_phi;
    std::vector<float> *T_Elec_BC3_energy;

    
    //conversion rejection
    std::vector<int> *T_Elec_MatchConv;
    std::vector<int> *T_Elec_EcalDriven;
    
    //jets and met
    std::vector<float> *T_Jet_Px;
    std::vector<float> *T_Jet_Py;
    std::vector<float> *T_Jet_Pz;
    std::vector<float> *T_Jet_Et;
    std::vector<float> *T_Jet_Eta;
    std::vector<float> *T_Jet_Energy;
    std::vector<float> *T_Jet_Phi;
    

    
    //met of the event
    float T_METPF_ET;
    float T_METPF_px;
    float T_METPF_py;
    float T_METPF_Phi;
    float T_METPF_Sig;
    float T_METPFTypeI_ET;
    float T_METPFTypeI_Phi;
    
};

typedef std::vector< edm::Handle< double > >   rhoHandles;


