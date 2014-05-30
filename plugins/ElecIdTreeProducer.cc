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
    
}


// ------------ method called once each job just before starting event loop  ------------
void 
ElecIdTreeProducer::beginJob()
{
    mytree_ = new TTree("eventsTree","");
    mytree_->Branch("T_Event_RunNumber", &T_Event_RunNumber, "T_Event_RunNumber/I");
    mytree_->Branch("T_Event_EventNumber", &T_Event_EventNumber, "T_Event_EventNumber/I");
    mytree_->Branch("T_Event_LuminosityBlock", &T_Event_LuminosityBlock, "T_Event_LuminosityBlock/I");
    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElecIdTreeProducer::endJob() 
{
    rootFile_->Write();
    rootFile_->Close();
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
