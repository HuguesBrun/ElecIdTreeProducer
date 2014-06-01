#include "ElecIdTreeProducer.h"


ElecIdTreeProducer::ElecIdTreeProducer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   isMC_                   = iConfig.getParameter<bool>("isMC");
    
   electronsCollection_      = iConfig.getParameter<edm::InputTag>("electronsCollection");
    
   outputFile_   = iConfig.getParameter<std::string>("outputFile");
   rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");
}


ElecIdTreeProducer::~ElecIdTreeProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElecIdTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    beginEvent(); //create the vectors
    
    ///load the collections:
    edm::Handle<reco::GsfElectronCollection> electronsCollection;
    iEvent.getByLabel(electronsCollection_ , electronsCollection);
    
    

    T_Event_RunNumber = iEvent.id().run();
    T_Event_EventNumber = iEvent.id().event();
    T_Event_LuminosityBlock = iEvent.id().luminosityBlock();
    
    
    /// now start the loop on the electrons:
    for(reco::GsfElectronCollection::const_iterator eleIt = electronsCollection->begin(); eleIt != electronsCollection->end(); eleIt++){
	cout << "pt=" << eleIt->pt() << endl;
    }

    mytree_->Fill();

    endEvent();    
}


// ------------ method called once each job just before starting event loop  ------------
void 
ElecIdTreeProducer::beginJob()
{
    mytree_ = new TTree("eventsTree","");
    mytree_->Branch("T_Event_RunNumber", &T_Event_RunNumber, "T_Event_RunNumber/I");
    mytree_->Branch("T_Event_EventNumber", &T_Event_EventNumber, "T_Event_EventNumber/I");
    mytree_->Branch("T_Event_LuminosityBlock", &T_Event_LuminosityBlock, "T_Event_LuminosityBlock/I");
   

mytree_->Branch("T_Elec_Eta", "std::vector<float>", &T_Elec_Eta);
mytree_->Branch("T_Elec_Px", "std::vector<float>", &T_Elec_Px);
mytree_->Branch("T_Elec_Py", "std::vector<float>", &T_Elec_Py);
mytree_->Branch("T_Elec_Pz", "std::vector<float>", &T_Elec_Pz);
mytree_->Branch("T_Elec_Pt", "std::vector<float>", &T_Elec_Pt);
mytree_->Branch("T_Elec_Energy", "std::vector<float>", &T_Elec_Energy);
mytree_->Branch("T_Elec_Charge", "std::vector<int>", &T_Elec_Charge);
mytree_->Branch("T_Elec_isEB", "std::vector<int>", &T_Elec_isEB);
mytree_->Branch("T_Elec_isEE", "std::vector<int>", &T_Elec_isEE);
mytree_->Branch("T_Elec_vz", "std::vector<float>", &T_Elec_vz);
mytree_->Branch("T_Elec_vy", "std::vector<float>", &T_Elec_vy);
mytree_->Branch("T_Elec_vx", "std::vector<float>", &T_Elec_vx);
mytree_->Branch("T_Elec_nLost", "std::vector<int>", &T_Elec_nLost);
mytree_->Branch("T_Elec_nHits", "std::vector<int>", &T_Elec_nHits);
mytree_->Branch("T_Elec_kfchi2", "std::vector<float>", &T_Elec_kfchi2);
mytree_->Branch("T_Elec_kfhits", "std::vector<int>", &T_Elec_kfhits);
mytree_->Branch("T_Elec_gsfhits", "std::vector<int>", &T_Elec_gsfhits);
mytree_->Branch("T_Elec_gsfchi2", "std::vector<float>", &T_Elec_gsfchi2);
mytree_->Branch("T_Elec_fbrem", "std::vector<float>", &T_Elec_fbrem);
mytree_->Branch("T_Elec_nbrems", "std::vector<int>", &T_Elec_nbrems);
mytree_->Branch("T_Elec_missingHits", "std::vector<int>", &T_Elec_missingHits);
mytree_->Branch("T_Elec_Dist", "std::vector<float>", &T_Elec_Dist);
mytree_->Branch("T_Elec_Dcot", "std::vector<float>", &T_Elec_Dcot);
mytree_->Branch("T_Elec_D0", "std::vector<float>", &T_Elec_D0);
mytree_->Branch("T_Elec_Dz", "std::vector<float>", &T_Elec_Dz);
mytree_->Branch("T_Elec_ip3d", "std::vector<float>", &T_Elec_ip3d);
mytree_->Branch("T_Elec_ip3ds", "std::vector<float>", &T_Elec_ip3ds);
mytree_->Branch("T_Elec_detacalo", "std::vector<float>", &T_Elec_detacalo);
mytree_->Branch("T_Elec_eledeta", "std::vector<float>", &T_Elec_eledeta);
mytree_->Branch("T_Elec_dphicalo", "std::vector<float>", &T_Elec_dphicalo);
mytree_->Branch("T_Elec_deltaPhiIn", "std::vector<float>", &T_Elec_deltaPhiIn);
mytree_->Branch("T_Elec_deltaEtaIn", "std::vector<float>", &T_Elec_deltaEtaIn);
mytree_->Branch("T_Elec_EoP", "std::vector<float>", &T_Elec_EoP);
mytree_->Branch("T_Elec_EoPin", "std::vector<float>", &T_Elec_EoPin);
mytree_->Branch("T_Elec_ESeedoP", "std::vector<float>", &T_Elec_ESeedoP);
mytree_->Branch("T_Elec_ESeedoPout", "std::vector<float>", &T_Elec_ESeedoPout);
mytree_->Branch("T_Elec_EEleoPout", "std::vector<float>", &T_Elec_EEleoPout);
mytree_->Branch("T_Elec_IoEmIoP", "std::vector<float>", &T_Elec_IoEmIoP);
mytree_->Branch("T_Elec_eleEoPout", "std::vector<float>", &T_Elec_eleEoPout);
mytree_->Branch("T_Elec_SC_Et", "std::vector<float>", &T_Elec_SC_Et);
mytree_->Branch("T_Elec_SC_Eta", "std::vector<float>", &T_Elec_SC_Eta);
mytree_->Branch("T_Elec_SC_Phi", "std::vector<float>", &T_Elec_SC_Phi);
mytree_->Branch("T_Elec_SC_RawEnergy", "std::vector<float>", &T_Elec_SC_RawEnergy);
mytree_->Branch("T_Elec_EcalEnergy", "std::vector<float>", &T_Elec_EcalEnergy);
mytree_->Branch("T_Elec_EsEnergy", "std::vector<float>", &T_Elec_EsEnergy);
mytree_->Branch("T_Elec_PreShowerOverRaw", "std::vector<float>", &T_Elec_PreShowerOverRaw);
mytree_->Branch("T_Elec_NClusters", "std::vector<int>", &T_Elec_NClusters);
mytree_->Branch("T_Elec_EtaSeed", "std::vector<float>", &T_Elec_EtaSeed);
mytree_->Branch("T_Elec_PhiSeed", "std::vector<float>", &T_Elec_PhiSeed);
mytree_->Branch("T_Elec_ESeed", "std::vector<float>", &T_Elec_ESeed);
mytree_->Branch("T_Elec_IEta", "std::vector<int>", &T_Elec_IEta);
mytree_->Branch("T_Elec_IPhi", "std::vector<int>", &T_Elec_IPhi);
mytree_->Branch("T_Elec_EtaSeedXtal", "std::vector<float>", &T_Elec_EtaSeedXtal);
mytree_->Branch("T_Elec_PhiSeedXtal", "std::vector<float>", &T_Elec_PhiSeedXtal);
mytree_->Branch("T_Elec_IEtaSeedXtal", "std::vector<float>", &T_Elec_IEtaSeedXtal);
mytree_->Branch("T_Elec_IPhiSeedXtal", "std::vector<float>", &T_Elec_IPhiSeedXtal);
mytree_->Branch("T_Elec_HtoE", "std::vector<float>", &T_Elec_HtoE);
mytree_->Branch("T_Elec_EmaxSeed", "std::vector<float>", &T_Elec_EmaxSeed);
mytree_->Branch("T_Elec_EtopSeed", "std::vector<float>", &T_Elec_EtopSeed);
mytree_->Branch("T_Elec_EbottomSeed", "std::vector<float>", &T_Elec_EbottomSeed);
mytree_->Branch("T_Elec_EleftSeed", "std::vector<float>", &T_Elec_EleftSeed);
mytree_->Branch("T_Elec_ErightSeed", "std::vector<float>", &T_Elec_ErightSeed);
mytree_->Branch("T_Elec_E2ndSeed", "std::vector<float>", &T_Elec_E2ndSeed);
mytree_->Branch("T_Elec_E2x5RightSeed", "std::vector<float>", &T_Elec_E2x5RightSeed);
mytree_->Branch("T_Elec_E2x5LeftSeed", "std::vector<float>", &T_Elec_E2x5LeftSeed);
mytree_->Branch("T_Elec_E2x5TopSeed", "std::vector<float>", &T_Elec_E2x5TopSeed);
mytree_->Branch("T_Elec_E2x5BottomSeed", "std::vector<float>", &T_Elec_E2x5BottomSeed);
mytree_->Branch("T_Elec_E2x5MaxSeed", "std::vector<float>", &T_Elec_E2x5MaxSeed);
mytree_->Branch("T_Elec_E1x5Seed", "std::vector<float>", &T_Elec_E1x5Seed);
mytree_->Branch("T_Elec_E2x2Seed", "std::vector<float>", &T_Elec_E2x2Seed);
mytree_->Branch("T_Elec_E3x3Seed", "std::vector<float>", &T_Elec_E3x3Seed);
mytree_->Branch("T_Elec_E5x5Seed", "std::vector<float>", &T_Elec_E5x5Seed);
mytree_->Branch("T_Elec_see", "std::vector<float>", &T_Elec_see);
mytree_->Branch("T_Elec_spp", "std::vector<float>", &T_Elec_spp);
mytree_->Branch("T_Elec_sep", "std::vector<float>", &T_Elec_sep);
mytree_->Branch("T_Elec_etawidth", "std::vector<float>", &T_Elec_etawidth);
mytree_->Branch("T_Elec_phiwidth", "std::vector<float>", &T_Elec_phiwidth);
mytree_->Branch("T_Elec_e1x5e5x5", "std::vector<float>", &T_Elec_e1x5e5x5);
mytree_->Branch("T_Elec_s9e25", "std::vector<float>", &T_Elec_s9e25);
mytree_->Branch("T_Elec_R9", "std::vector<float>", &T_Elec_R9);
mytree_->Branch("T_Elec_MatchConv", "std::vector<int>", &T_Elec_MatchConv);
mytree_->Branch("T_Elec_EcalDriven", "std::vector<int>", &T_Elec_EcalDriven); 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElecIdTreeProducer::endJob() 
{
    rootFile_->Write();
    rootFile_->Close();
}

void 
ElecIdTreeProducer::beginEvent()
{
T_Elec_Eta = new std::vector<float>;
T_Elec_Px = new std::vector<float>;
T_Elec_Py = new std::vector<float>;
T_Elec_Pz = new std::vector<float>;
T_Elec_Pt = new std::vector<float>;
T_Elec_Energy = new std::vector<float>;
T_Elec_Charge = new std::vector<int>;
T_Elec_isEB = new std::vector<int>;
T_Elec_isEE = new std::vector<int>;
T_Elec_vz = new std::vector<float>;
T_Elec_vy = new std::vector<float>;
T_Elec_vx = new std::vector<float>;
T_Elec_nLost = new std::vector<int>;
T_Elec_nHits = new std::vector<int>;
T_Elec_kfchi2 = new std::vector<float>;
T_Elec_kfhits = new std::vector<int>;
T_Elec_gsfhits = new std::vector<int>;
T_Elec_gsfchi2 = new std::vector<float>;
T_Elec_fbrem = new std::vector<float>;
T_Elec_nbrems = new std::vector<int>;
T_Elec_missingHits = new std::vector<int>;
T_Elec_Dist = new std::vector<float>;
T_Elec_Dcot = new std::vector<float>;
T_Elec_D0 = new std::vector<float>;
T_Elec_Dz = new std::vector<float>;
T_Elec_ip3d = new std::vector<float>;
T_Elec_ip3ds = new std::vector<float>;
T_Elec_detacalo = new std::vector<float>;
T_Elec_eledeta = new std::vector<float>;
T_Elec_dphicalo = new std::vector<float>;
T_Elec_deltaPhiIn = new std::vector<float>;
T_Elec_deltaEtaIn = new std::vector<float>;
T_Elec_EoP = new std::vector<float>;
T_Elec_EoPin = new std::vector<float>;
T_Elec_ESeedoP = new std::vector<float>;
T_Elec_ESeedoPout = new std::vector<float>;
T_Elec_EEleoPout = new std::vector<float>;
T_Elec_IoEmIoP = new std::vector<float>;
T_Elec_eleEoPout = new std::vector<float>;
T_Elec_SC_Et = new std::vector<float>;
T_Elec_SC_Eta = new std::vector<float>;
T_Elec_SC_Phi = new std::vector<float>;
T_Elec_SC_RawEnergy = new std::vector<float>;
T_Elec_EcalEnergy = new std::vector<float>;
T_Elec_EsEnergy = new std::vector<float>;
T_Elec_PreShowerOverRaw = new std::vector<float>;
T_Elec_NClusters = new std::vector<int>;
T_Elec_EtaSeed = new std::vector<float>;
T_Elec_PhiSeed = new std::vector<float>;
T_Elec_ESeed = new std::vector<float>;
T_Elec_IEta = new std::vector<int>;
T_Elec_IPhi = new std::vector<int>;
T_Elec_EtaSeedXtal = new std::vector<float>;
T_Elec_PhiSeedXtal = new std::vector<float>;
T_Elec_IEtaSeedXtal = new std::vector<float>;
T_Elec_IPhiSeedXtal = new std::vector<float>;
T_Elec_HtoE = new std::vector<float>;
T_Elec_EmaxSeed = new std::vector<float>;
T_Elec_EtopSeed = new std::vector<float>;
T_Elec_EbottomSeed = new std::vector<float>;
T_Elec_EleftSeed = new std::vector<float>;
T_Elec_ErightSeed = new std::vector<float>;
T_Elec_E2ndSeed = new std::vector<float>;
T_Elec_E2x5RightSeed = new std::vector<float>;
T_Elec_E2x5LeftSeed = new std::vector<float>;
T_Elec_E2x5TopSeed = new std::vector<float>;
T_Elec_E2x5BottomSeed = new std::vector<float>;
T_Elec_E2x5MaxSeed = new std::vector<float>;
T_Elec_E1x5Seed = new std::vector<float>;
T_Elec_E2x2Seed = new std::vector<float>;
T_Elec_E3x3Seed = new std::vector<float>;
T_Elec_E5x5Seed = new std::vector<float>;
T_Elec_see = new std::vector<float>;
T_Elec_spp = new std::vector<float>;
T_Elec_sep = new std::vector<float>;
T_Elec_etawidth = new std::vector<float>;
T_Elec_phiwidth = new std::vector<float>;
T_Elec_e1x5e5x5 = new std::vector<float>;
T_Elec_s9e25 = new std::vector<float>;
T_Elec_R9 = new std::vector<float>;
T_Elec_MatchConv = new std::vector<int>;
T_Elec_EcalDriven = new std::vector<int>;
}
void
ElecIdTreeProducer::endEvent()
{
delete T_Elec_Eta;
delete T_Elec_Px;
delete T_Elec_Py;
delete T_Elec_Pz;
delete T_Elec_Pt;
delete T_Elec_Energy;
delete T_Elec_Charge;
delete T_Elec_isEB;
delete T_Elec_isEE;
delete T_Elec_vz;
delete T_Elec_vy;
delete T_Elec_vx;
delete T_Elec_nLost;
delete T_Elec_nHits;
delete T_Elec_kfchi2;
delete T_Elec_kfhits;
delete T_Elec_gsfhits;
delete T_Elec_gsfchi2;
delete T_Elec_fbrem;
delete T_Elec_nbrems;
delete T_Elec_missingHits;
delete T_Elec_Dist;
delete T_Elec_Dcot;
delete T_Elec_D0;
delete T_Elec_Dz;
delete T_Elec_ip3d;
delete T_Elec_ip3ds;
delete T_Elec_detacalo;
delete T_Elec_eledeta;
delete T_Elec_dphicalo;
delete T_Elec_deltaPhiIn;
delete T_Elec_deltaEtaIn;
delete T_Elec_EoP;
delete T_Elec_EoPin;
delete T_Elec_ESeedoP;
delete T_Elec_ESeedoPout;
delete T_Elec_EEleoPout;
delete T_Elec_IoEmIoP;
delete T_Elec_eleEoPout;
delete T_Elec_SC_Et;
delete T_Elec_SC_Eta;
delete T_Elec_SC_Phi;
delete T_Elec_SC_RawEnergy;
delete T_Elec_EcalEnergy;
delete T_Elec_EsEnergy;
delete T_Elec_PreShowerOverRaw;
delete T_Elec_NClusters;
delete T_Elec_EtaSeed;
delete T_Elec_PhiSeed;
delete T_Elec_ESeed;
delete T_Elec_IEta;
delete T_Elec_IPhi;
delete T_Elec_EtaSeedXtal;
delete T_Elec_PhiSeedXtal;
delete T_Elec_IEtaSeedXtal;
delete T_Elec_IPhiSeedXtal;
delete T_Elec_HtoE;
delete T_Elec_EmaxSeed;
delete T_Elec_EtopSeed;
delete T_Elec_EbottomSeed;
delete T_Elec_EleftSeed;
delete T_Elec_ErightSeed;
delete T_Elec_E2ndSeed;
delete T_Elec_E2x5RightSeed;
delete T_Elec_E2x5LeftSeed;
delete T_Elec_E2x5TopSeed;
delete T_Elec_E2x5BottomSeed;
delete T_Elec_E2x5MaxSeed;
delete T_Elec_E1x5Seed;
delete T_Elec_E2x2Seed;
delete T_Elec_E3x3Seed;
delete T_Elec_E5x5Seed;
delete T_Elec_see;
delete T_Elec_spp;
delete T_Elec_sep;
delete T_Elec_etawidth;
delete T_Elec_phiwidth;
delete T_Elec_e1x5e5x5;
delete T_Elec_s9e25;
delete T_Elec_R9;
delete T_Elec_MatchConv;
delete T_Elec_EcalDriven;
}


// ------------ method called when starting to processes a run  ------------
/*
void 
ElecIdTreeProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ElecIdTreeProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ElecIdTreeProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ElecIdTreeProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElecIdTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters

}

//define this as a plug-in
DEFINE_FWK_MODULE(ElecIdTreeProducer);
