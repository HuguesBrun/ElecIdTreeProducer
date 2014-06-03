#include "ElecIdTreeProducer.h"


ElecIdTreeProducer::ElecIdTreeProducer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   isMC_                   = iConfig.getParameter<bool>("isMC");
    
   electronsCollection_      = iConfig.getParameter<edm::InputTag>("electronsCollection");
   primaryVertexInputTag_    = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");

    EBRecHitsLabel_           = iConfig.getParameter<edm::InputTag>("rechitCollectionEB");
    EERecHitsLabel_           = iConfig.getParameter<edm::InputTag>("rechitCollectionEE");
    
    conversionsInputTag_ = iConfig.getParameter<edm::InputTag>("conversionsCollection");
    beamSpotInputTag_ = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
    rhoInputTags_ = iConfig.getParameter<std::vector<edm::InputTag> >("rhoTags");
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
    
    // load the vertices collection
    edm::Handle<reco::VertexCollection> vtx_h;
    iEvent.getByLabel(primaryVertexInputTag_, vtx_h);
    
    //load the conversion collection
    edm::Handle<reco::ConversionCollection> conversions_h;
    iEvent.getByLabel(conversionsInputTag_, conversions_h);
    
    // get the beam spot
    edm::Handle<reco::BeamSpot> beamspot_h;
    iEvent.getByLabel(beamSpotInputTag_, beamspot_h);
    const reco::BeamSpot &beamSpot = *(beamspot_h.product());
    
    
    //load the RecHits
    edm::Handle< EcalRecHitCollection > EBRecHits;
    iEvent.getByLabel(EBRecHitsLabel_ , EBRecHits);
    edm::Handle< EcalRecHitCollection > EERecHits;
    iEvent.getByLabel(EERecHitsLabel_ , EERecHits);
    
    
    // handles to gen infos
    edm::Handle<GenEventInfoProduct> genEvent;
    edm::Handle <reco::GenParticleCollection> genParticles;

    
    
    //load the transiant tracks builder
    edm::ESHandle<TransientTrackBuilder> builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
    TransientTrackBuilder thebuilder = *(builder.product());
    
    
    

    //initialize the ECAL geom, if not yet done
    if (!(isGeomInitialized_)){
             edm::ESHandle<CaloTopology> theCaloTopology;
             iSetup.get<CaloTopologyRecord>().get(theCaloTopology);
             ecalTopology_ = & (*theCaloTopology);
    
             edm::ESHandle<CaloGeometry> theCaloGeometry;
             iSetup.get<CaloGeometryRecord>().get(theCaloGeometry);
             caloGeometry_ = & (*theCaloGeometry);
             isGeomInitialized_ = true;
    }
    
    

    T_Event_RunNumber = iEvent.id().run();
    T_Event_EventNumber = iEvent.id().event();
    T_Event_LuminosityBlock = iEvent.id().luminosityBlock();
    
    //fill rho
    for (int iteRho = 0 ; iteRho)
    
    float truePu=0.;
    if (isMC_){
        iEvent.getByLabel("generator", genEvent);
        iEvent.getByLabel( "genParticles", genParticles );
        
        T_Event_processID = genEvent->signalProcessID();
        if ( genEvent->binningValues().size()>0 ) T_Event_ptHat = genEvent->binningValues()[0];
        
        Handle<std::vector< PileupSummaryInfo > > puInfo;
        try {
            iEvent.getByLabel("addPileupInfo",puInfo);
            std::vector<PileupSummaryInfo>::const_iterator PVI;
            //The in-time crossing is getBunchCrossing = 0; negative ones are early, positive ones are late.
            for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) {
                
                //    std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
                if(PVI->getBunchCrossing()==0){
                    T_Event_nPU =PVI->getPU_NumInteractions();
                    T_Event_nTruePU=PVI->getTrueNumInteractions();
                    
                }
                
                else if(PVI->getBunchCrossing()==-1){
                    T_Event_nPUm=PVI->getPU_NumInteractions();
                }
                else if(PVI->getBunchCrossing()==1){
                    T_Event_nPUp=PVI->getPU_NumInteractions();
                }
                truePu += PVI->getTrueNumInteractions();
            }
        } catch (...) {}
        T_Event_AveNTruePU=truePu;
    }
    
    
    
    /// now start the loop on the electrons:
    for(reco::GsfElectronCollection::const_iterator eleIt = electronsCollection->begin(); eleIt != electronsCollection->end(); eleIt++){
	     cout << "pt=" << eleIt->pt() << endl;
        T_Elec_Eta->push_back(eleIt->eta());
        T_Elec_Phi->push_back(eleIt->phi());
        T_Elec_Px->push_back(eleIt->px());
        T_Elec_Py->push_back(eleIt->py());
        T_Elec_Pz->push_back(eleIt->pz());
        T_Elec_Pt->push_back(eleIt->pt());
        T_Elec_Energy->push_back(eleIt->energy());
        T_Elec_Charge->push_back(eleIt->charge());
        
        T_Elec_vx->push_back(eleIt->vx());
        T_Elec_vy->push_back(eleIt->vy());
        T_Elec_vz->push_back(eleIt->vz());

        T_Elec_nLost->push_back(eleIt->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits());
        T_Elec_nHits->push_back(eleIt->gsfTrack()->trackerExpectedHitsInner().numberOfHits());
        T_Elec_gsfhits->push_back(eleIt->gsfTrack()->numberOfValidHits());
        T_Elec_gsfchi2->push_back(eleIt->gsfTrack()->normalizedChi2());
        
        T_Elec_fbrem->push_back(eleIt->fbrem());
        T_Elec_nbrems->push_back(eleIt->numberOfBrems());
        
        
        T_Elec_Dist->push_back(eleIt->convDist());
        T_Elec_Dcot->push_back(eleIt->convDcot());
        
        T_Elec_D0->push_back(eleIt->gsfTrack()->d0());
        T_Elec_Dz->push_back(eleIt->gsfTrack()->dz());
        
        // ip stuff
        float pv_ip3d = -999.0;
        float pv_ip3dSig = 0.0;
        if (eleIt->gsfTrack().isNonnull()) {
            const double gsfsign   = ( (-eleIt->gsfTrack()->dxy(vtx_h->at(0).position()))   >=0 ) ? 1. : -1.;
            
            const reco::TransientTrack &tt = thebuilder.build(eleIt->gsfTrack());
            const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,vtx_h->at(0));
            if (ip3dpv.first) {
                double ip3d = gsfsign*ip3dpv.second.value();
                double ip3derr = ip3dpv.second.error();
                pv_ip3d = ip3d;
                pv_ip3dSig = ip3d/ip3derr;
            }
        }
        
        
        T_Elec_ip3d->push_back(pv_ip3d);
        T_Elec_ip3ds->push_back(pv_ip3dSig);
        
        
        
        //all kfStuff
        bool validKF= false;
        reco::TrackRef myTrackRef = eleIt->closestTrack();
        validKF = (myTrackRef.isNonnull());
        T_Elec_kfchi2->push_back((validKF) ? myTrackRef->normalizedChi2() : 0 );
        T_Elec_kfhits->push_back((validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : 0 );
        
        
        //track cluster matching variables
        T_Elec_detacalo->push_back(eleIt->deltaEtaSeedClusterTrackAtCalo());
        T_Elec_dphicalo->push_back(eleIt->deltaPhiSeedClusterTrackAtCalo());
        
        T_Elec_deltaPhiIn->push_back(eleIt->deltaPhiSuperClusterTrackAtVtx());
        T_Elec_deltaEtaIn->push_back(eleIt->deltaEtaSuperClusterTrackAtVtx());
        
        T_Elec_eledeta->push_back(eleIt->deltaEtaEleClusterTrackAtCalo());
        T_Elec_eledphi->push_back(eleIt->deltaPhiEleClusterTrackAtCalo());
        
        T_Elec_EoP->push_back(eleIt->eSuperClusterOverP());
        T_Elec_ESeedoP->push_back(eleIt->eSeedClusterOverP());
        T_Elec_ESeedoPout->push_back(eleIt->eSeedClusterOverPout());
        T_Elec_EEleoPout->push_back(eleIt->eEleClusterOverPout());
        
        
        float IoEmIo = -1;
        if (!((eleIt->ecalEnergy()==0) || (eleIt->p()==0))) IoEmIo = 1.0/eleIt->ecalEnergy() - 1.0/eleIt->p();
        T_Elec_IoEmIoP->push_back(IoEmIo);
        
        T_Elec_SC_Et->push_back((TMath::CosH(eleIt->superCluster()->eta())!=0)? eleIt->superCluster()->energy()/TMath::CosH(eleIt->superCluster()->eta()): -1);
        T_Elec_SC_Eta->push_back(eleIt->superCluster()->eta());
        T_Elec_SC_Phi->push_back(eleIt->superCluster()->phi());
        T_Elec_SC_RawEnergy->push_back(eleIt->superCluster()->rawEnergy());
        T_Elec_EcalEnergy->push_back(eleIt->superCluster()->energy());
        T_Elec_EsEnergy->push_back(eleIt->superCluster()->preshowerEnergy());
        T_Elec_PreShowerOverRaw->push_back((eleIt->superCluster()->rawEnergy()>0) ? eleIt->superCluster()->preshowerEnergy() / eleIt->superCluster()->rawEnergy() : -1);
        T_Elec_NClusters->push_back(eleIt->superCluster()->clustersSize());
        
        T_Elec_EtaSeed->push_back(eleIt->superCluster()->seed()->eta());
        T_Elec_PhiSeed->push_back(eleIt->superCluster()->seed()->phi());
        T_Elec_ESeed->push_back(eleIt->superCluster()->seed()->energy());

        
        T_Elec_HtoE->push_back(eleIt->hadronicOverEm());
        
        const EcalRecHitCollection * recHits=0;
        if( eleIt->isEB() ) recHits = EBRecHits.product();
        else recHits = EERecHits.product();
     
        SuperClusterHelper mySCHelper( &(*eleIt), recHits, ecalTopology_, caloGeometry_ );
        
        T_Elec_EmaxSeed->push_back((eleIt->superCluster()->energy()>0) ? mySCHelper.eMax() / eleIt->superCluster()->energy() : -1);
        T_Elec_EtopSeed->push_back((eleIt->superCluster()->energy()>0) ? mySCHelper.eTop() / eleIt->superCluster()->energy() : -1);
        T_Elec_EbottomSeed->push_back((eleIt->superCluster()->energy()>0) ? mySCHelper.eBottom() / eleIt->superCluster()->energy() : -1);
        T_Elec_EleftSeed->push_back((eleIt->superCluster()->energy()>0) ? mySCHelper.eLeft() / eleIt->superCluster()->energy() : -1);
        T_Elec_ErightSeed->push_back((eleIt->superCluster()->energy()>0) ? mySCHelper.eRight() / eleIt->superCluster()->energy() : -1);
        T_Elec_E2ndSeed->push_back((eleIt->superCluster()->energy()>0) ? mySCHelper.e2nd() / eleIt->superCluster()->energy() : -1);
        
        
        
        T_Elec_E2x5RightSeed->push_back((eleIt->superCluster()->energy()>0) ? mySCHelper.e2x5Right() / eleIt->superCluster()->energy() : -1);
        T_Elec_E2x5LeftSeed->push_back((eleIt->superCluster()->energy()>0) ? mySCHelper.e2x5Left() / eleIt->superCluster()->energy() : -1);
        T_Elec_E2x5TopSeed->push_back((eleIt->superCluster()->energy()>0) ? mySCHelper.e2x5Top() / eleIt->superCluster()->energy() : -1);
        T_Elec_E2x5BottomSeed->push_back((eleIt->superCluster()->energy()>0) ? mySCHelper.e2x5Bottom() / eleIt->superCluster()->energy() : -1);
        T_Elec_E2x5MaxSeed->push_back((eleIt->superCluster()->energy()>0) ? mySCHelper.e2x5Max() / eleIt->superCluster()->energy() : -1);
        
        // remove 2x2 1x5
        T_Elec_E3x3Seed->push_back((eleIt->superCluster()->energy()>0) ? mySCHelper.e3x3() / eleIt->superCluster()->energy() : -1);
        T_Elec_E5x5Seed->push_back((eleIt->superCluster()->energy()>0) ? mySCHelper.e5x5() / eleIt->superCluster()->energy() : -1);
        
        T_Elec_see->push_back(mySCHelper.sigmaIetaIeta());
        T_Elec_sep->push_back(mySCHelper.sep());
        T_Elec_spp->push_back(mySCHelper.spp());
        
        T_Elec_etawidth->push_back(mySCHelper.etaWidth());
        T_Elec_phiwidth->push_back(mySCHelper.phiWidth());
        
        T_Elec_s9e25->push_back((mySCHelper.e5x5()>0) ? mySCHelper.e3x3()/mySCHelper.e5x5() : -1);
        
        T_Elec_e1x5e5x5->push_back((eleIt->e5x5()) !=0. ? 1.-(eleIt->e1x5()/eleIt->e5x5()) : -1.);
        
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(*eleIt, conversions_h, beamSpot.position());
        T_Elec_MatchConv->push_back(vtxFitConversion);
        T_Elec_EcalDriven->push_back(eleIt->ecalDrivenSeed());

        
    }

    mytree_->Fill();

    endEvent();    
}


// ------------ method called once each job just before starting event loop  ------------
void 
ElecIdTreeProducer::beginJob()
{
    
    isGeomInitialized_ = false;
    
    mytree_ = new TTree("eventsTree","");
    mytree_->Branch("T_Event_RunNumber", &T_Event_RunNumber, "T_Event_RunNumber/I");
    mytree_->Branch("T_Event_EventNumber", &T_Event_EventNumber, "T_Event_EventNumber/I");
    mytree_->Branch("T_Event_LuminosityBlock", &T_Event_LuminosityBlock, "T_Event_LuminosityBlock/I");
   

    mytree_->Branch("T_Event_Rho", "std::vector<float>", &T_Event_Rho);
    mytree_->Branch("T_Elec_Eta", "std::vector<float>", &T_Elec_Eta);
    mytree_->Branch("T_Elec_Phi", "std::vector<float>", &T_Elec_Phi);
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
    mytree_->Branch("T_Elec_eledphi", "std::vector<float>", &T_Elec_eledphi);
mytree_->Branch("T_Elec_dphicalo", "std::vector<float>", &T_Elec_dphicalo);
mytree_->Branch("T_Elec_deltaPhiIn", "std::vector<float>", &T_Elec_deltaPhiIn);
mytree_->Branch("T_Elec_deltaEtaIn", "std::vector<float>", &T_Elec_deltaEtaIn);
mytree_->Branch("T_Elec_EoP", "std::vector<float>", &T_Elec_EoP);
mytree_->Branch("T_Elec_ESeedoP", "std::vector<float>", &T_Elec_ESeedoP);
mytree_->Branch("T_Elec_ESeedoPout", "std::vector<float>", &T_Elec_ESeedoPout);
mytree_->Branch("T_Elec_EEleoPout", "std::vector<float>", &T_Elec_EEleoPout);
mytree_->Branch("T_Elec_IoEmIoP", "std::vector<float>", &T_Elec_IoEmIoP);
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
    T_Event_Rho = new std::vector<float>;
    T_Elec_Eta = new std::vector<float>;
    T_Elec_Phi = new std::vector<float>;
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
    T_Elec_eledphi = new std::vector<float>;
T_Elec_dphicalo = new std::vector<float>;
T_Elec_deltaPhiIn = new std::vector<float>;
T_Elec_deltaEtaIn = new std::vector<float>;
T_Elec_EoP = new std::vector<float>;
T_Elec_ESeedoP = new std::vector<float>;
T_Elec_ESeedoPout = new std::vector<float>;
T_Elec_EEleoPout = new std::vector<float>;
T_Elec_IoEmIoP = new std::vector<float>;
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
    delete T_Event_Rho;
    delete T_Elec_Eta;
    delete T_Elec_Phi;
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
    delete T_Elec_eledphi;
delete T_Elec_dphicalo;
delete T_Elec_deltaPhiIn;
delete T_Elec_deltaEtaIn;
delete T_Elec_EoP;
delete T_Elec_ESeedoP;
delete T_Elec_ESeedoPout;
delete T_Elec_EEleoPout;
delete T_Elec_IoEmIoP;
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
