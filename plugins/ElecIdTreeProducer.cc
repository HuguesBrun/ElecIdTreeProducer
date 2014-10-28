#include "ElecIdTreeProducer.h"


ElecIdTreeProducer::ElecIdTreeProducer(const edm::ParameterSet& iConfig)

{
    //now do what ever initialization is needed
    isMC_                   = iConfig.getParameter<bool>("isMC");
    
    electronsCollection_      = iConfig.getParameter<edm::InputTag>("electronsCollection");
    primaryVertexInputTag_    = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
    
    EBRecHitsLabel_           = iConfig.getParameter<edm::InputTag>("rechitCollectionEB");
    EERecHitsLabel_           = iConfig.getParameter<edm::InputTag>("rechitCollectionEE");

    ecalRechitEBToken_ = consumes<EcalRecHitCollection>(EBRecHitsLabel_);
    ecalRechitEEToken_ = consumes<EcalRecHitCollection>(EERecHitsLabel_);
    
    conversionsInputTag_ = iConfig.getParameter<edm::InputTag>("conversionsCollection");
    beamSpotInputTag_ = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
    rhoInputTags_ = iConfig.getParameter<std::vector<edm::InputTag> >("rhoTags");
    metTag_ = iConfig.getParameter<edm::InputTag>("metTag");
    jetCollectionTag_= iConfig.getParameter<edm::InputTag>("jetCollectionTag");
    muonProducers_			= iConfig.getParameter<vtag>("muonProducer");
    triggerResultsTag_ = iConfig.getParameter<edm::InputTag>("triggerResultTag");
    triggerSummaryLabel_= iConfig.getParameter<edm::InputTag>("triggerSummaryTag");
    pathsToSave_ = iConfig.getParameter<std::vector<std::string> >("pathsToSave");
    filterToMatch_ = iConfig.getParameter<std::vector<std::string> >("filterToMatch");
    HLTprocess_   = iConfig.getParameter<std::string>("HLTprocess");
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
  /*  edm::Handle< EcalRecHitCollection > EBRecHits;
    iEvent.getByLabel(EBRecHitsLabel_ , EBRecHits);
    edm::Handle< EcalRecHitCollection > EERecHits;
    iEvent.getByLabel(EERecHitsLabel_ , EERecHits);*/
    
    
    // handles to gen infos
    edm::Handle<GenEventInfoProduct> genEvent;
    edm::Handle <reco::GenParticleCollection> genParticles;
    
    
    
    //load the transiant tracks builder
    edm::ESHandle<TransientTrackBuilder> builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
    TransientTrackBuilder thebuilder = *(builder.product());
    
    
    // Jet
    edm::Handle < std::vector <reco::PFJet> > recoPFJets;
    iEvent.getByLabel(jetCollectionTag_, recoPFJets);
    int nJets = recoPFJets->size();
    
    //mEt collection
    edm::Handle<reco::PFMETCollection> metPF;
    iEvent.getByLabel(metTag_,metPF);
    const reco::PFMET * metsPF= &((metPF.product())->front());
    
    //muon collection :
		edm::Handle < std::vector <reco::Muon> > recoMuons;
    edm::InputTag muonProducer = muonProducers_.at(0);
	iEvent.getByLabel(muonProducer, recoMuons);
    
    
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
     rhoHandles rhos(rhoInputTags_.size());
     for (unsigned int iteRho = 0 ; iteRho < rhoInputTags_.size() ; iteRho++){
         iEvent.getByLabel(rhoInputTags_[iteRho], rhos[iteRho]);
         T_Event_Rho->push_back(*rhos[iteRho]);
     }
    

    noZS::EcalClusterLazyTools lazyToolsNoZS(iEvent, iSetup, ecalRechitEBToken_, ecalRechitEBToken_);
    EcalClusterLazyTools lazyTools(iEvent, iSetup, ecalRechitEBToken_, ecalRechitEEToken_);

    
    float truePu=0.;
    int theNbOfGenParticles = 0;
    if (isMC_){
        iEvent.getByLabel("generator", genEvent);
        iEvent.getByLabel( "genParticles", genParticles );
        theNbOfGenParticles = genParticles->size();
        
        
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
    
    
    reco::Vertex dummy;
    const reco::Vertex *pv = &dummy;
    if (vtx_h->size() != 0) {
        pv = &*vtx_h->begin();
    } else { // create a dummy PV
        reco::Vertex::Error e;
        e(0, 0) = 0.0015 * 0.0015;
        e(1, 1) = 0.0015 * 0.0015;
        e(2, 2) = 15. * 15.;
        reco::Vertex::Point p(0, 0, 0);
        dummy = reco::Vertex(p, e, 0, 0, 0);
    }
    
    
    // now look at the trigger info !

    
    bool changedConfig = false;
    if (!hltConfig.init(iEvent.getRun(), iSetup, HLTprocess_.c_str(), changedConfig)) {
        cout << "Initialization of HLTConfigProvider failed!!" << endl;
        return;
    }
    if (changedConfig){
        unsigned int nbPaths = pathsToSave_.size();

        std::cout << "the curent menu is " << hltConfig.tableName() << std::endl;
        for (unsigned int itePath=0 ; itePath<nbPaths ; itePath++){
            for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
                if (TString(hltConfig.triggerNames()[j]).Contains(pathsToSave_.at(itePath))){
                    triggerBits_.push_back(j);
                    cout << "found the path " << pathsToSave_.at(itePath) << endl;
                }

            }
        }
        if (triggerBits_.size() <nbPaths) cout << "an HLT paths is not found ! ! " << endl;
        
    }
    
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerResultsTag_, triggerResults);
    


    
    for (unsigned int itePath = 0 ; itePath < triggerBits_.size() ; itePath++){
        if (triggerResults->accept(triggerBits_.at(itePath))) {
		T_Event_pathsFired->push_back(1);
	}
        else T_Event_pathsFired->push_back(0);
    }
    
    
    edm::Handle<trigger::TriggerEvent> triggerSummary;
    iEvent.getByLabel(triggerSummaryLabel_, triggerSummary);
    trigger::TriggerObjectCollection allTriggerObjects = triggerSummary->getObjects();
    trigger::TriggerObjectCollection legObjects;
    std::vector<unsigned int> legRefs;
    // find the ref of the legs
    for (size_t iteFilter=0; iteFilter<filterToMatch_.size(); iteFilter++) {
	edm::InputTag filterTag = edm::InputTag(filterToMatch_.at(iteFilter), "", HLTprocess_);
        size_t filterIndex = (*triggerSummary).filterIndex(filterTag);
        if (filterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
            //save the trigger objects corresponding to muon leg
            const trigger::Keys &keys = (*triggerSummary).filterKeys(filterIndex);
            for (size_t j = 0; j < keys.size(); j++) {
                trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
                legObjects.push_back(foundObject);
                legRefs.push_back(iteFilter);
            }
        }
    }
    
    // cout << "nbGen=" << theNbOfGenParticles << endl;
    
    /// now start the loop on the electrons:
    for(reco::GsfElectronCollection::const_iterator eleIt = electronsCollection->begin(); eleIt != electronsCollection->end(); eleIt++){
          //cout << "pt=" << eleIt->pt() << endl;
        
        // if in MC then do the matching with MC particles
        if (isMC_){
            int theGenPartRef = -1;
            float minDr = 1000;
            int iteMinDr=-1;
            for (int iteGen = 0 ; iteGen < theNbOfGenParticles ; iteGen++){
                const reco::GenParticle & genElectron = (*genParticles)[iteGen];
                if (fabs(genElectron.pdgId())!=11) continue;
                float deltaR = sqrt(pow(eleIt->eta()-genElectron.eta(),2)+ pow(acos(cos(eleIt->phi()-genElectron.phi())),2)) ;
                if (deltaR<minDr){
                    minDr = deltaR;
                    iteMinDr = iteGen;
                }
            }
            if (minDr<0.1) theGenPartRef = iteMinDr;
            
 
            if (theGenPartRef>=0){
                const reco::GenParticle & genElectron = (*genParticles)[theGenPartRef];
                T_Gen_Elec_softElectron->push_back(isNonPromptElectron(genElectron));
                T_Gen_Elec_Px->push_back(genElectron.px());
                T_Gen_Elec_Py->push_back(genElectron.py());
                T_Gen_Elec_Pz->push_back(genElectron.pz());
                T_Gen_Elec_Energy->push_back(genElectron.energy());
                T_Gen_Elec_PDGid->push_back(genElectron.pdgId());
                T_Gen_Elec_status->push_back(genElectron.status());
                int  isFromTau = hasTauasMother(genElectron);
                T_Gen_Elec_fromTAU->push_back(isFromTau);
                if (genElectron.numberOfMothers()>0) {
                    const reco::Candidate  *part = (genElectron.mother());
                    const reco::Candidate  *MomPart =(genElectron.mother()); //dummy initialisation :)
                    // loop on the  particles to check if is has a W has mother
                    while ((part->numberOfMothers()>0)) {
                        MomPart =part->mother();
                        if ((fabs(MomPart->pdgId())>=22)&&(fabs(MomPart->pdgId())<=24)){
                            break;
                        }
                        part = MomPart;
                    }
                    T_Gen_Elec_MotherID->push_back(MomPart->pdgId());
                    if (MomPart->numberOfMothers()>0) {
                        const reco::Candidate  *grandMa = MomPart->mother();
                        T_Gen_Elec_GndMotherID->push_back(grandMa->pdgId());
                    }
                    else T_Gen_Elec_GndMotherID->push_back(-1);
                }
                else {
                    T_Gen_Elec_MotherID->push_back(-1);
                    T_Gen_Elec_GndMotherID->push_back(-1);
                }
            }
            else{
                T_Gen_Elec_Px->push_back(-1);
                T_Gen_Elec_Py->push_back(-1);
                T_Gen_Elec_Pz->push_back(-1);
                T_Gen_Elec_Energy->push_back(-1);
                T_Gen_Elec_PDGid->push_back(-1);
                T_Gen_Elec_MotherID->push_back(-1);
                T_Gen_Elec_GndMotherID->push_back(-1);
                T_Gen_Elec_status->push_back(-1);
                T_Gen_Elec_fromTAU->push_back(-1);
                T_Gen_Elec_softElectron->push_back(-1);
            }
            
            
        }
        
        //now look if the electron is matched with trigger
        //FIXME remove the double couting
        int theLegInfo = 0;
        for (unsigned int iteTrigObj = 0 ; iteTrigObj < filterToMatch_.size() ; iteTrigObj++){
            bool foundTheLeg = false;
            for (unsigned int i = 0 ; i < legObjects.size() ; i++){
                if (legRefs.at(i)==iteTrigObj) continue;
                float deltaR = sqrt(pow(legObjects[i].eta()-eleIt->eta(),2)+ pow(acos(cos(legObjects[i].phi()-eleIt->phi())),2)) ;
                if (deltaR<0.1) {
                    foundTheLeg = true;
                }
            }
            if (foundTheLeg){
            theLegInfo += std::pow(2,iteTrigObj);
            }
        }
        T_Elec_TriggerLeg->push_back(theLegInfo);
        
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
        
        T_Elec_nLost->push_back(eleIt->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
        //T_Elec_nHits->push_back(eleIt->gsfTrack()d$d$->trackerExpectedHitsInner().numberOfHits()); // the same as the previous one
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
        
        /*const EcalRecHitCollection * recHits=0;
        if( eleIt->isEB() ) recHits = EBRecHits.product();
        else recHits = EERecHits.product();*/
        
     /*   SuperClusterHelper mySCHelper( &(*eleIt), recHits, ecalTopology_, caloGeometry_ );*/
        const auto & seedCluster = eleIt->superCluster()->seed();

      
      
        T_Elec_EmaxSeed->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.eMax(*seedCluster) / eleIt->superCluster()->energy() : -1);
        T_Elec_EtopSeed->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.eTop(*seedCluster) / eleIt->superCluster()->energy() : -1);
        T_Elec_EbottomSeed->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.eBottom(*seedCluster) / eleIt->superCluster()->energy() : -1);
        T_Elec_EleftSeed->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.eLeft(*seedCluster) / eleIt->superCluster()->energy() : -1);
        T_Elec_ErightSeed->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.eRight(*seedCluster) / eleIt->superCluster()->energy() : -1);
        T_Elec_E2ndSeed->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.e2nd(*seedCluster) / eleIt->superCluster()->energy() : -1);
        
        
        
        T_Elec_E2x5RightSeed->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.e2x5Right(*seedCluster) / eleIt->superCluster()->energy() : -1);
        T_Elec_E2x5LeftSeed->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.e2x5Left(*seedCluster) / eleIt->superCluster()->energy() : -1);
        T_Elec_E2x5TopSeed->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.e2x5Top(*seedCluster) / eleIt->superCluster()->energy() : -1);
        T_Elec_E2x5BottomSeed->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.e2x5Bottom(*seedCluster) / eleIt->superCluster()->energy() : -1);
        T_Elec_E2x5MaxSeed->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.e2x5Max(*seedCluster) / eleIt->superCluster()->energy() : -1);
        
        // remove 2x2 1x5
        T_Elec_E3x3Seed->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.e3x3(*seedCluster) / eleIt->superCluster()->energy() : -1);
        T_Elec_E5x5Seed->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.e5x5(*seedCluster) / eleIt->superCluster()->energy() : -1);
        
        std::vector<float> vCovN = lazyTools.localCovariances(*seedCluster);
        T_Elec_see->push_back(sqrt(vCovN[0]));
        T_Elec_sep->push_back(sqrt(vCovN[2]));
        T_Elec_spp->push_back(sqrt(vCovN[1]));
        
        T_Elec_etawidth->push_back(eleIt->superCluster()->etaWidth());
        T_Elec_phiwidth->push_back(eleIt->superCluster()->phiWidth());
        
        T_Elec_s9e25->push_back((lazyTools.e5x5(*seedCluster)>0) ? lazyTools.e3x3(*seedCluster)/lazyTools.e5x5(*seedCluster) : -1);
        
        T_Elec_e1x5e5x5->push_back((lazyTools.e5x5(*seedCluster)) !=0. ? 1.-(eleIt->e1x5()/eleIt->e5x5()) : -1.);
        
        
        T_Elec_R9->push_back((eleIt->superCluster()->energy()>0) ? lazyTools.e3x3(*seedCluster) / eleIt->superCluster()->energy() : -1);
   
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(*eleIt, conversions_h, beamSpot.position());
        T_Elec_MatchConv->push_back(vtxFitConversion);
        T_Elec_EcalDriven->push_back(eleIt->ecalDrivenSeed());
        
        std::vector<float> vCov = lazyToolsNoZS.localCovariances(*seedCluster);
        
        T_Elec_noZSsee->push_back(  sqrt(vCov[0]));
        T_Elec_noZSspp->push_back(  sqrt(vCov[2]));
        T_Elec_noZSsep->push_back(  sqrt(vCov[1]));
        T_Elec_noZSr9->push_back((eleIt->superCluster()->energy()>0) ? lazyToolsNoZS.e3x3(*seedCluster)/eleIt->superCluster()->energy(): -1);
        T_Elec_noZSe1x5->push_back(lazyToolsNoZS.e1x5(*seedCluster));
        T_Elec_noZSe2x5MaxSeed->push_back(lazyToolsNoZS.e2x5Max(*seedCluster));
        T_Elec_noZSe5x5->push_back(lazyToolsNoZS.e5x5(*seedCluster));
      
       
        
        T_Elec_puChargedIso->push_back((*eleIt).pfIsolationVariables().sumChargedParticlePt);
        T_Elec_allChargedHadronIso->push_back((*eleIt).pfIsolationVariables().sumPUPt);
        T_Elec_chargedHadronIso->push_back((*eleIt).pfIsolationVariables().sumChargedHadronPt);
        T_Elec_neutralHadronIso->push_back((*eleIt).pfIsolationVariables().sumNeutralHadronEt);
        T_Elec_photonIso->push_back((*eleIt).pfIsolationVariables().sumPhotonEt);
        
        
        T_Elec_ECALiso->push_back(eleIt->dr03EcalRecHitSumEt());
        T_Elec_HCALiso->push_back(eleIt->dr03HcalTowerSumEt());
        T_Elec_TKiso->push_back(eleIt->dr03TkSumPt());
        
        
     //   fill the infos on the 3 first BC
        int nbBC = 1;
        for (reco::CaloCluster_iterator subBC = eleIt->superCluster()->clustersBegin(); subBC != eleIt->superCluster()->clustersEnd(); ++subBC){
            switch (nbBC) {
                case 1:
                    T_Elec_BC1_eta->push_back((*subBC)->eta());
                    T_Elec_BC1_phi->push_back((*subBC)->phi());
                    T_Elec_BC1_energy->push_back((*subBC)->energy());
                    break;
                    
                case 2:
                    T_Elec_BC2_eta->push_back((*subBC)->eta());
                    T_Elec_BC2_phi->push_back((*subBC)->phi());
                    T_Elec_BC2_energy->push_back((*subBC)->energy());
                    break;
                    
                case 3:
                    T_Elec_BC3_eta->push_back((*subBC)->eta());
                    T_Elec_BC3_phi->push_back((*subBC)->phi());
                    T_Elec_BC3_energy->push_back((*subBC)->energy());
                    break;
                    
                default:
                    break;
            }
            nbBC++;
        }
        
        T_Elec_nbBC->push_back(nbBC);
        
        if (nbBC==1) {
            T_Elec_BC2_eta->push_back(-1);
            T_Elec_BC2_phi->push_back(-1);
            T_Elec_BC2_energy->push_back(-1);
            
            T_Elec_BC3_eta->push_back(-1);
            T_Elec_BC3_phi->push_back(-1);
            T_Elec_BC3_energy->push_back(-1);
        }
        
        if (nbBC==2) {
            T_Elec_BC3_eta->push_back(-1);
            T_Elec_BC3_phi->push_back(-1);
            T_Elec_BC3_energy->push_back(-1);
        }
        
        
    }
    
    //save the jets
    for (int k = 0 ; k < nJets ; k++){
        const reco::Jet* jet = (const reco::Jet*) ( & ((*recoPFJets)[k]) );
        T_Jet_Px->push_back(jet->px());
        T_Jet_Py->push_back(jet->py());
        T_Jet_Pz->push_back(jet->pz());
        T_Jet_Et->push_back(jet->et());
        T_Jet_Eta->push_back(jet->eta());
        T_Jet_Energy->push_back(jet->energy());
        T_Jet_Phi->push_back(jet->phi());
    }
    
    
    //save the met
    T_METPF_ET = metsPF[0].pt();
    T_METPF_px = metsPF[0].px();
    T_METPF_py = metsPF[0].py();
    T_METPF_Phi = metsPF[0].phi();
    T_METPF_Sig = metsPF[0].significance();
   // T_METPFTypeI_ET =
 //   T_METPFTypeI_Phi =
    
    
    int nbMuons = recoMuons->size();
    //cout << "il y a " << nbMuons << " muons " << endl;
    //loop on the muons in the event
    for (int k = 0 ; k < nbMuons ; k++){
        
        const reco::Muon* muon = &((*recoMuons)[k]);
        //  cout << "le muon : eta=" << muon->eta() << " phi=" << muon->phi() << endl;
        
        if (isMC_){
            int theGenPartRef = -1;
            float minDr = 1000;
            int iteMinDr=-1;
            for (int iteGen = 0 ; iteGen < theNbOfGenParticles ; iteGen++){
                const reco::GenParticle & genMuon = (*genParticles)[iteGen];
                if (fabs(genMuon.pdgId())!=13) continue;
                float deltaR = sqrt(pow(muon->eta()-genMuon.eta(),2)+ pow(acos(cos(muon->phi()-genMuon.phi())),2)) ;
                if (deltaR<minDr){
                    minDr = deltaR;
                    iteMinDr = iteGen;
                }
            }
            if (minDr<0.1) theGenPartRef = iteMinDr;
            
            
            if (theGenPartRef>=0){
                const reco::GenParticle & genMuon = (*genParticles)[theGenPartRef];
                T_Gen_Muon_softMuon->push_back(isNonPromptElectron(genMuon));
                T_Gen_Muon_Px->push_back(genMuon.px());
                T_Gen_Muon_Py->push_back(genMuon.py());
                T_Gen_Muon_Pz->push_back(genMuon.pz());
                T_Gen_Muon_Energy->push_back(genMuon.energy());
                T_Gen_Muon_PDGid->push_back(genMuon.pdgId());
                T_Gen_Muon_status->push_back(genMuon.status());
                int  isFromTau = hasTauasMother(genMuon);
                T_Gen_Muon_fromTAU->push_back(isFromTau);
                if (genMuon.numberOfMothers()>0) {
                    const reco::Candidate  *part = (genMuon.mother());
                    const reco::Candidate  *MomPart =(genMuon.mother()); //dummy initialisation :)
                    // loop on the  particles to check if is has a W has mother
                    while ((part->numberOfMothers()>0)) {
                        MomPart =part->mother();
                        if ((fabs(MomPart->pdgId())>=22)&&(fabs(MomPart->pdgId())<=24)){
                            break;
                        }
                        part = MomPart;
                    }
                    T_Gen_Muon_MotherID->push_back(MomPart->pdgId());
                    if (MomPart->numberOfMothers()>0) {
                        const reco::Candidate  *grandMa = MomPart->mother();
                        T_Gen_Muon_GndMotherID->push_back(grandMa->pdgId());
                    }
                    else T_Gen_Muon_GndMotherID->push_back(-1);
                }
                else {
                    T_Gen_Muon_MotherID->push_back(-1);
                    T_Gen_Muon_GndMotherID->push_back(-1);
                }
            }
            else{
                T_Gen_Muon_Px->push_back(-1);
                T_Gen_Muon_Py->push_back(-1);
                T_Gen_Muon_Pz->push_back(-1);
                T_Gen_Muon_Energy->push_back(-1);
                T_Gen_Muon_PDGid->push_back(-1);
                T_Gen_Muon_MotherID->push_back(-1);
                T_Gen_Muon_GndMotherID->push_back(-1);
                T_Gen_Muon_status->push_back(-1);
                T_Gen_Muon_fromTAU->push_back(-1);
                T_Gen_Muon_softMuon->push_back(-1);
            }
            
            
        }
        //now look if the electron is matched with trigger
        //FIXME remove the double couting
        int theLegInfo = 0;
        for (unsigned int iteTrigObj = 0 ; iteTrigObj < filterToMatch_.size() ; iteTrigObj++){
            bool foundTheLeg = false;
            for (unsigned int i = 0 ; i < legObjects.size() ; i++){
                if (legRefs.at(i)==iteTrigObj) continue;
                float deltaR = sqrt(pow(legObjects[i].eta()-muon->eta(),2)+ pow(acos(cos(legObjects[i].phi()-muon->phi())),2)) ;
                if (deltaR<0.1) {
                    foundTheLeg = true;
                }
            }
            if (foundTheLeg){
                theLegInfo += std::pow(2,iteTrigObj);
            }
        }
        T_Muon_TriggerLeg->push_back(theLegInfo);
        
        
        T_Muon_Eta->push_back(muon->eta());
        T_Muon_Phi->push_back(muon->phi());
        T_Muon_IsGlobalMuon->push_back(muon->isGlobalMuon());
        T_Muon_IsPFMuon->push_back(muon->isPFMuon());
        T_Muon_IsTrackerMuon->push_back(muon->isTrackerMuon());
        T_Muon_IsCaloMuon->push_back(muon->isCaloMuon());
        T_Muon_IsStandAloneMuon->push_back(muon->isStandAloneMuon());
        T_Muon_IsMuon->push_back(muon->isMuon());
        T_Muon_Energy->push_back(muon->energy());
        T_Muon_Et->push_back(muon->et());
        T_Muon_Pt->push_back(muon->pt());
        T_Muon_Px->push_back(muon->px());
        T_Muon_Py->push_back(muon->py());
        T_Muon_Pz->push_back(muon->pz());
        T_Muon_Mass->push_back(muon->mass());
        T_Muon_charge->push_back(muon->charge());
        
        T_Muon_numberOfChambers->push_back(muon->numberOfChambers());
        T_Muon_numberOfChambersRPC->push_back(muon->numberOfChambersNoRPC());
        T_Muon_numberOfMatches->push_back(muon->numberOfMatches());
        T_Muon_numberOfMatchedStations->push_back(muon->numberOfMatchedStations());
        bool isMatchTheStation = muon::isGoodMuon(*muon, muon::TMOneStationTight);
        bool isGlobalMuonPT = muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight);
        bool isGlobalMuonArbitrated = muon::isGoodMuon(*muon, muon::TrackerMuonArbitrated);
        T_Muon_TMLastStationTight->push_back(isMatchTheStation);
        T_Muon_IsGlobalMuon_PromptTight->push_back(isGlobalMuonPT);
        T_Muon_IsTrackerMuonArbitrated->push_back(isGlobalMuonArbitrated);
        
        if (muon->globalTrack().isNull()) T_Muon_globalTrackChi2->push_back(-1); else T_Muon_globalTrackChi2->push_back(muon->globalTrack()->normalizedChi2());
        if (muon->globalTrack().isNull()) T_Muon_validMuonHits->push_back(-1); else T_Muon_validMuonHits->push_back(muon->globalTrack()->hitPattern().numberOfValidMuonHits());
        T_Muon_trkKink->push_back(muon->combinedQuality().trkKink);
        if (muon->muonBestTrack().isNull()) {
            T_Muon_trkNbOfTrackerLayers->push_back(-1);
            T_Muon_trkError->push_back(-1);
            T_Muon_dB->push_back(-1);
            T_Muon_dBstop->push_back(-1);
            T_Muon_dzPV->push_back(-1);
            T_Muon_dzstop->push_back(-1);
            T_Muon_trkValidPixelHits->push_back(-1);
            T_Muon_trkNbOfValidTrackeHits->push_back(-1);
        }
        else {
            T_Muon_trkNbOfTrackerLayers->push_back(muon->muonBestTrack()->hitPattern().trackerLayersWithMeasurement());
            T_Muon_trkError->push_back(muon->muonBestTrack()->ptError());
            T_Muon_trkValidPixelHits->push_back(muon->muonBestTrack()->hitPattern().numberOfValidPixelHits());
            T_Muon_dB->push_back(fabs(muon->muonBestTrack()->dxy(pv->position())));
            T_Muon_dBstop->push_back(fabs(muon->muonBestTrack()->dxy(pv->position())));
            T_Muon_dzPV->push_back(fabs(muon->muonBestTrack()->dz(pv->position())));
            T_Muon_dzstop->push_back(fabs(muon->muonBestTrack()->dz(pv->position())));
            T_Muon_trkNbOfValidTrackeHits->push_back(muon->muonBestTrack()->hitPattern().numberOfValidTrackerHits());
        }
        T_Muon_isoR03_emEt->push_back(muon->isolationR03().emEt);
        T_Muon_isoR03_hadEt->push_back(muon->isolationR03().hadEt);
        T_Muon_isoR03_hoEt->push_back(muon->isolationR03().hoEt);
        T_Muon_isoR03_sumPt->push_back(muon->isolationR03().sumPt);
        T_Muon_isoR03_nTracks->push_back(muon->isolationR03().nTracks);
        T_Muon_isoR03_nJets->push_back(muon->isolationR03().nJets);
        T_Muon_chargedHadronIsoR04->push_back(muon->pfIsolationR04().sumChargedHadronPt);
        T_Muon_neutralHadronIsoR04->push_back(muon->pfIsolationR04().sumNeutralHadronEt);
        T_Muon_photonIsoR04->push_back(muon->pfIsolationR04().sumPhotonEt);
        T_Muon_chargedHadronIsoPUR04->push_back(muon->pfIsolationR04().sumPUPt);
        T_Muon_chargedHadronIsoR03->push_back(muon->pfIsolationR03().sumChargedHadronPt);
        T_Muon_neutralHadronIsoR03->push_back(muon->pfIsolationR03().sumNeutralHadronEt);
        T_Muon_photonIsoR03->push_back(muon->pfIsolationR03().sumPhotonEt);
        T_Muon_chargedHadronIsoPUR03->push_back(muon->pfIsolationR03().sumPUPt);

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
    
    mytree_->Branch("T_Event_nPU", &T_Event_nPU, "T_Event_nPU/I");
    mytree_->Branch("T_Event_nTruePU", &T_Event_nTruePU, "T_Event_nTruePU/F");
    mytree_->Branch("T_Event_nPUm", &T_Event_nPUm, "T_Event_nPUm/I");
    mytree_->Branch("T_Event_nPUp", &T_Event_nPUp, "T_Event_nPUp/I");
    mytree_->Branch("T_Event_AveNTruePU", &T_Event_AveNTruePU, "T_Event_AveNTruePU/F");

    
    
    mytree_->Branch("T_Event_Rho", "std::vector<float>", &T_Event_Rho);
    
    mytree_->Branch("T_Event_pathsFired", "std::vector<int>", &T_Event_pathsFired);
    
    
    mytree_->Branch("T_Gen_Elec_Px", "std::vector<float>", &T_Gen_Elec_Px);
    mytree_->Branch("T_Gen_Elec_Py", "std::vector<float>", &T_Gen_Elec_Py);
    mytree_->Branch("T_Gen_Elec_Pz", "std::vector<float>", &T_Gen_Elec_Pz);
    mytree_->Branch("T_Gen_Elec_Energy", "std::vector<float>", &T_Gen_Elec_Energy);
    mytree_->Branch("T_Gen_Elec_status", "std::vector<int>", &T_Gen_Elec_status);
    mytree_->Branch("T_Gen_Elec_PDGid", "std::vector<int>", &T_Gen_Elec_PDGid);
    mytree_->Branch("T_Gen_Elec_MotherID", "std::vector<int>", &T_Gen_Elec_MotherID);
    mytree_->Branch("T_Gen_Elec_GndMotherID", "std::vector<int>", &T_Gen_Elec_GndMotherID);
    mytree_->Branch("T_Gen_Elec_fromTAU", "std::vector<int>", &T_Gen_Elec_fromTAU);
    mytree_->Branch("T_Gen_Elec_softElectron", "std::vector<int>", &T_Gen_Elec_softElectron);
    
    
    mytree_->Branch("T_Elec_TriggerLeg", "std::vector<int>", &T_Elec_TriggerLeg);
    
    mytree_->Branch("T_Elec_puChargedIso", "std::vector<float>", &T_Elec_puChargedIso);
    mytree_->Branch("T_Elec_allChargedHadronIso", "std::vector<float>", &T_Elec_allChargedHadronIso);
    mytree_->Branch("T_Elec_chargedHadronIso", "std::vector<float>", &T_Elec_chargedHadronIso);
    mytree_->Branch("T_Elec_neutralHadronIso", "std::vector<float>", &T_Elec_neutralHadronIso);
    mytree_->Branch("T_Elec_photonIso", "std::vector<float>", &T_Elec_photonIso);
    mytree_->Branch("T_Elec_puChargedIso04", "std::vector<float>", &T_Elec_puChargedIso04);
    mytree_->Branch("T_Elec_allChargedHadronIso04", "std::vector<float>", &T_Elec_allChargedHadronIso04);
    mytree_->Branch("T_Elec_chargedHadronIso04", "std::vector<float>", &T_Elec_chargedHadronIso04);
    mytree_->Branch("T_Elec_neutralHadronIso04", "std::vector<float>", &T_Elec_neutralHadronIso04);
    mytree_->Branch("T_Elec_photonIso04", "std::vector<float>", &T_Elec_photonIso04);
    
    
    
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
    
    
    mytree_->Branch("T_Elec_noZSsee", "std::vector<float>", &T_Elec_noZSsee);
    mytree_->Branch("T_Elec_noZSspp", "std::vector<float>", &T_Elec_noZSspp);
    mytree_->Branch("T_Elec_noZSsep", "std::vector<float>", &T_Elec_noZSsep);
    mytree_->Branch("T_Elec_noZSr9", "std::vector<float>", &T_Elec_noZSr9);
    mytree_->Branch("T_Elec_noZSe1x5", "std::vector<float>", &T_Elec_noZSe1x5);
    mytree_->Branch("T_Elec_noZSe2x5MaxSeed", "std::vector<float>", &T_Elec_noZSe2x5MaxSeed);
    mytree_->Branch("T_Elec_noZSe5x5", "std::vector<float>", &T_Elec_noZSe5x5);

    
    
    mytree_->Branch("T_Elec_ECALiso", "std::vector<float>", &T_Elec_ECALiso);
    mytree_->Branch("T_Elec_HCALiso", "std::vector<float>", &T_Elec_HCALiso);
    mytree_->Branch("T_Elec_TKiso", "std::vector<float>", &T_Elec_TKiso);
    
    
    mytree_->Branch("T_Elec_nbBC", "std::vector<int>", &T_Elec_nbBC);
    mytree_->Branch("T_Elec_BC1_eta", "std::vector<float>", &T_Elec_BC1_eta);
    mytree_->Branch("T_Elec_BC1_phi", "std::vector<float>", &T_Elec_BC1_phi);
    mytree_->Branch("T_Elec_BC1_energy", "std::vector<float>", &T_Elec_BC1_energy);

    mytree_->Branch("T_Elec_BC2_eta", "std::vector<float>", &T_Elec_BC1_eta);
    mytree_->Branch("T_Elec_BC2_phi", "std::vector<float>", &T_Elec_BC1_phi);
    mytree_->Branch("T_Elec_BC2_energy", "std::vector<float>", &T_Elec_BC1_energy);
    
    mytree_->Branch("T_Elec_BC3_eta", "std::vector<float>", &T_Elec_BC3_eta);
    mytree_->Branch("T_Elec_BC3_phi", "std::vector<float>", &T_Elec_BC3_phi);
    mytree_->Branch("T_Elec_BC3_energy", "std::vector<float>", &T_Elec_BC3_energy);
    
    //muons
    
    mytree_->Branch("T_Muon_TriggerLeg", "std::vector<int>", &T_Muon_TriggerLeg);

    mytree_->Branch("T_Muon_Eta", "std::vector<float>", &T_Muon_Eta);
    mytree_->Branch("T_Muon_Phi", "std::vector<float>", &T_Muon_Phi);
    mytree_->Branch("T_Muon_Energy", "std::vector<float>", &T_Muon_Energy);
    mytree_->Branch("T_Muon_Et", "std::vector<float>", &T_Muon_Et);
    mytree_->Branch("T_Muon_Pt", "std::vector<float>", &T_Muon_Pt);
    mytree_->Branch("T_Muon_Px", "std::vector<float>", &T_Muon_Px);
    mytree_->Branch("T_Muon_Py", "std::vector<float>", &T_Muon_Py);
    mytree_->Branch("T_Muon_Pz", "std::vector<float>", &T_Muon_Pz);
    mytree_->Branch("T_Muon_Mass", "std::vector<float>", &T_Muon_Mass);
    mytree_->Branch("T_Muon_IsGlobalMuon", "std::vector<bool>", &T_Muon_IsGlobalMuon);
    mytree_->Branch("T_Muon_IsTrackerMuon", "std::vector<bool>", &T_Muon_IsTrackerMuon);
    mytree_->Branch("T_Muon_IsPFMuon", "std::vector<bool>", &T_Muon_IsPFMuon);
    mytree_->Branch("T_Muon_IsCaloMuon", "std::vector<bool>", &T_Muon_IsCaloMuon);
    mytree_->Branch("T_Muon_IsStandAloneMuon", "std::vector<bool>", &T_Muon_IsStandAloneMuon);
    mytree_->Branch("T_Muon_IsMuon", "std::vector<bool>", &T_Muon_IsMuon);
    mytree_->Branch("T_Muon_IsGlobalMuon_PromptTight", "std::vector<bool>", &T_Muon_IsGlobalMuon_PromptTight);
    mytree_->Branch("T_Muon_IsTrackerMuonArbitrated", "std::vector<bool>", &T_Muon_IsTrackerMuonArbitrated);
    mytree_->Branch("T_Muon_numberOfChambers", "std::vector<int>", &T_Muon_numberOfChambers);
    mytree_->Branch("T_Muon_numberOfChambersRPC", "std::vector<int>", &T_Muon_numberOfChambersRPC);
    mytree_->Branch("T_Muon_numberOfMatches", "std::vector<int>", &T_Muon_numberOfMatches);
    mytree_->Branch("T_Muon_numberOfMatchedStations", "std::vector<int>", &T_Muon_numberOfMatchedStations);
    mytree_->Branch("T_Muon_charge", "std::vector<int>", &T_Muon_charge);
    mytree_->Branch("T_Muon_TMLastStationTight", "std::vector<bool>", &T_Muon_TMLastStationTight);
    mytree_->Branch("T_Muon_globalTrackChi2", "std::vector<float>", &T_Muon_globalTrackChi2);
    mytree_->Branch("T_Muon_validMuonHits", "std::vector<int>", &T_Muon_validMuonHits);
    mytree_->Branch("T_Muon_trkKink", "std::vector<float>", &T_Muon_trkKink);
    mytree_->Branch("T_Muon_trkNbOfTrackerLayers", "std::vector<int>", &T_Muon_trkNbOfTrackerLayers);
    mytree_->Branch("T_Muon_trkNbOfValidTrackeHits", "std::vector<int>", &T_Muon_trkNbOfValidTrackeHits);
    mytree_->Branch("T_Muon_trkValidPixelHits", "std::vector<int>", &T_Muon_trkValidPixelHits);
    mytree_->Branch("T_Muon_trkError", "std::vector<float>", &T_Muon_trkError);
    mytree_->Branch("T_Muon_dB", "std::vector<float>", &T_Muon_dB);
    mytree_->Branch("T_Muon_dzPV", "std::vector<float>", &T_Muon_dzPV);
    mytree_->Branch("T_Muon_dBstop", "std::vector<float>", &T_Muon_dBstop);
    mytree_->Branch("T_Muon_dzstop", "std::vector<float>", &T_Muon_dzstop);
    mytree_->Branch("T_Muon_chargedHadronIsoR04", "std::vector<float>", &T_Muon_chargedHadronIsoR04);
    mytree_->Branch("T_Muon_neutralHadronIsoR04", "std::vector<float>", &T_Muon_neutralHadronIsoR04);
    mytree_->Branch("T_Muon_photonIsoR04", "std::vector<float>", &T_Muon_photonIsoR04);
    mytree_->Branch("T_Muon_chargedHadronIsoPUR04", "std::vector<float>", &T_Muon_chargedHadronIsoPUR04);
    mytree_->Branch("T_Muon_chargedHadronIsoR03", "std::vector<float>", &T_Muon_chargedHadronIsoR03);
    mytree_->Branch("T_Muon_neutralHadronIsoR03", "std::vector<float>", &T_Muon_neutralHadronIsoR03);
    mytree_->Branch("T_Muon_photonIsoR03", "std::vector<float>", &T_Muon_photonIsoR03);
    mytree_->Branch("T_Muon_chargedHadronIsoPUR03", "std::vector<float>", &T_Muon_chargedHadronIsoPUR03);
    mytree_->Branch("T_Muon_isoR03_emEt", "std::vector<float>", &T_Muon_isoR03_emEt);
    mytree_->Branch("T_Muon_isoR03_hadEt", "std::vector<float>", &T_Muon_isoR03_hadEt);
    mytree_->Branch("T_Muon_isoR03_hoEt", "std::vector<float>", &T_Muon_isoR03_hoEt);
    mytree_->Branch("T_Muon_isoR03_sumPt", "std::vector<float>", &T_Muon_isoR03_sumPt);
    mytree_->Branch("T_Muon_isoR03_nTracks", "std::vector<int>", &T_Muon_isoR03_nTracks);
    mytree_->Branch("T_Muon_isoR03_nJets", "std::vector<int>", &T_Muon_isoR03_nJets);
    mytree_->Branch("T_Muon_isoRingsMVA", "std::vector<float>", &T_Muon_isoRingsMVA);
    
    
    mytree_->Branch("T_Gen_Muon_Px", "std::vector<float>", &T_Gen_Muon_Px);
    mytree_->Branch("T_Gen_Muon_Py", "std::vector<float>", &T_Gen_Muon_Py);
    mytree_->Branch("T_Gen_Muon_Pz", "std::vector<float>", &T_Gen_Muon_Pz);
    mytree_->Branch("T_Gen_Muon_Energy", "std::vector<float>", &T_Gen_Muon_Energy);
    mytree_->Branch("T_Gen_Muon_status", "std::vector<int>", &T_Gen_Muon_status);
    mytree_->Branch("T_Gen_Muon_PDGid", "std::vector<int>", &T_Gen_Muon_PDGid);
    mytree_->Branch("T_Gen_Muon_MotherID", "std::vector<int>", &T_Gen_Muon_MotherID);
    mytree_->Branch("T_Gen_Muon_GndMotherID", "std::vector<int>", &T_Gen_Muon_GndMotherID);
    mytree_->Branch("T_Gen_Muon_fromTAU", "std::vector<int>", &T_Gen_Muon_fromTAU);
    mytree_->Branch("T_Gen_Muon_softMuon", "std::vector<int>", &T_Gen_Muon_softMuon);
    
    //jets
    mytree_->Branch("T_Jet_Px", "std::vector<float>", &T_Jet_Px);
    mytree_->Branch("T_Jet_Py", "std::vector<float>", &T_Jet_Py);
    mytree_->Branch("T_Jet_Pz", "std::vector<float>", &T_Jet_Pz);
    mytree_->Branch("T_Jet_Et", "std::vector<float>", &T_Jet_Et);
    mytree_->Branch("T_Jet_Eta", "std::vector<float>", &T_Jet_Eta);
    mytree_->Branch("T_Jet_Energy", "std::vector<float>", &T_Jet_Energy);
    mytree_->Branch("T_Jet_Phi", "std::vector<float>", &T_Jet_Phi);
    
    
    //MET
    mytree_->Branch("T_METPF_ET", &T_METPF_ET, "T_METPF_ET/F");
    mytree_->Branch("T_METPF_px", &T_METPF_px, "T_METPF_px/F");
    mytree_->Branch("T_METPF_py", &T_METPF_py, "T_METPF_py/F");
    mytree_->Branch("T_METPF_Phi", &T_METPF_Phi, "T_METPF_Phi/F");
    mytree_->Branch("T_METPF_Sig", &T_METPF_Sig, "T_METPF_Sig/F");
    mytree_->Branch("T_METPFTypeI_ET", &T_METPFTypeI_ET, "T_METPFTypeI_ET/F");
    mytree_->Branch("T_METPFTypeI_Phi", &T_METPFTypeI_Phi, "T_METPFTypeI_Phi/F");


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
    
    T_Event_pathsFired = new std::vector<int>;
    
    
    T_Elec_puChargedIso = new std::vector<float>;
    T_Elec_allChargedHadronIso = new std::vector<float>;
    T_Elec_chargedHadronIso = new std::vector<float>;
    T_Elec_neutralHadronIso = new std::vector<float>;
    T_Elec_photonIso = new std::vector<float>;
    
    T_Elec_puChargedIso04 = new std::vector<float>;
    T_Elec_allChargedHadronIso04 = new std::vector<float>;
    T_Elec_chargedHadronIso04 = new std::vector<float>;
    T_Elec_neutralHadronIso04 = new std::vector<float>;
    T_Elec_photonIso04 = new std::vector<float>;
    
    
    T_Gen_Elec_Px = new std::vector<float>;
    T_Gen_Elec_Py = new std::vector<float>;
    T_Gen_Elec_Pz = new std::vector<float>;
    T_Gen_Elec_Energy = new std::vector<float>;
    T_Gen_Elec_PDGid = new std::vector<int>;
    T_Gen_Elec_MotherID = new std::vector<int>;
    T_Gen_Elec_GndMotherID = new std::vector<int>;
    T_Gen_Elec_status = new std::vector<int>;
    T_Gen_Elec_fromTAU = new std::vector<int>;
    T_Gen_Elec_softElectron = new std::vector<int>;
    
    T_Elec_TriggerLeg = new std::vector<int>;
    
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
    
    
    T_Elec_noZSsee = new std::vector<float>;
    T_Elec_noZSspp = new std::vector<float>;
    T_Elec_noZSsep = new std::vector<float>;
    T_Elec_noZSr9 = new std::vector<float>;
    T_Elec_noZSe1x5 = new std::vector<float>;
    T_Elec_noZSe2x5MaxSeed = new std::vector<float>;
    T_Elec_noZSe5x5 = new std::vector<float>;

    
    
    T_Elec_ECALiso = new std::vector<float>;
    T_Elec_HCALiso = new std::vector<float>;
    T_Elec_TKiso = new std::vector<float>;
    
    T_Elec_nbBC = new std::vector<int>;
    T_Elec_BC1_eta = new std::vector<float>;
    T_Elec_BC1_phi = new std::vector<float>;
    T_Elec_BC1_energy = new std::vector<float>;

    T_Elec_BC2_eta = new std::vector<float>;
    T_Elec_BC2_phi = new std::vector<float>;
    T_Elec_BC2_energy = new std::vector<float>;
    
    T_Elec_BC3_eta = new std::vector<float>;
    T_Elec_BC3_phi = new std::vector<float>;
    T_Elec_BC3_energy = new std::vector<float>;
    
    
    T_Muon_TriggerLeg = new std::vector<int>;
    
    T_Muon_Eta = new std::vector<float>;
    T_Muon_Phi = new std::vector<float>;
    T_Muon_Energy = new std::vector<float>;
    T_Muon_Et = new std::vector<float>;
    T_Muon_Pt = new std::vector<float>;
    T_Muon_Px = new std::vector<float>;
    T_Muon_Py = new std::vector<float>;
    T_Muon_Pz = new std::vector<float>;
    T_Muon_Mass = new std::vector<float>;
    T_Muon_IsGlobalMuon = new std::vector<bool>;
    T_Muon_IsTrackerMuon = new std::vector<bool>;
    T_Muon_IsPFMuon = new std::vector<bool>;
    T_Muon_IsCaloMuon = new std::vector<bool>;
    T_Muon_IsStandAloneMuon = new std::vector<bool>;
    T_Muon_IsMuon = new std::vector<bool>;
    T_Muon_IsGlobalMuon_PromptTight = new std::vector<bool>;
    T_Muon_IsTrackerMuonArbitrated = new std::vector<bool>;
    T_Muon_numberOfChambers = new std::vector<int>;
    T_Muon_numberOfChambersRPC = new std::vector<int>;
    T_Muon_numberOfMatches = new std::vector<int>;
    T_Muon_numberOfMatchedStations = new std::vector<int>;
    T_Muon_charge = new std::vector<int>;
    T_Muon_TMLastStationTight = new std::vector<bool>;
    T_Muon_globalTrackChi2 = new std::vector<float>;
    T_Muon_validMuonHits = new std::vector<int>;
    T_Muon_trkKink = new std::vector<float>;
    T_Muon_trkNbOfTrackerLayers = new std::vector<int>;
    T_Muon_trkNbOfValidTrackeHits = new std::vector<int>;
    T_Muon_trkValidPixelHits = new std::vector<int>;
    T_Muon_trkError = new std::vector<float>;
    T_Muon_dB = new std::vector<float>;
    T_Muon_dzPV = new std::vector<float>;
    T_Muon_dBstop = new std::vector<float>;
    T_Muon_dzstop = new std::vector<float>;
    T_Muon_chargedHadronIsoR04 = new std::vector<float>;
    T_Muon_neutralHadronIsoR04 = new std::vector<float>;
    T_Muon_photonIsoR04 = new std::vector<float>;
    T_Muon_chargedHadronIsoPUR04 = new std::vector<float>;
    T_Muon_chargedHadronIsoR03 = new std::vector<float>;
    T_Muon_neutralHadronIsoR03 = new std::vector<float>;
    T_Muon_photonIsoR03 = new std::vector<float>;
    T_Muon_chargedHadronIsoPUR03 = new std::vector<float>;
    T_Muon_isoR03_emEt = new std::vector<float>;
    T_Muon_isoR03_hadEt = new std::vector<float>;
    T_Muon_isoR03_hoEt = new std::vector<float>;
    T_Muon_isoR03_sumPt = new std::vector<float>;
    T_Muon_isoR03_nTracks = new std::vector<int>;
    T_Muon_isoR03_nJets = new std::vector<int>;
    T_Muon_isoRingsMVA = new std::vector<float>;
    
    T_Gen_Muon_Px = new std::vector<float>;
    T_Gen_Muon_Py = new std::vector<float>;
    T_Gen_Muon_Pz = new std::vector<float>;
    T_Gen_Muon_Energy = new std::vector<float>;
    T_Gen_Muon_status = new std::vector<int>;
    T_Gen_Muon_PDGid = new std::vector<int>;
    T_Gen_Muon_MotherID = new std::vector<int>;
    T_Gen_Muon_GndMotherID = new std::vector<int>;
    T_Gen_Muon_fromTAU = new std::vector<int>;
    T_Gen_Muon_softMuon = new std::vector<int>;
    
    
    T_Jet_Px = new std::vector<float>;
    T_Jet_Py = new std::vector<float>;
    T_Jet_Pz = new std::vector<float>;
    T_Jet_Et = new std::vector<float>;
    T_Jet_Eta = new std::vector<float>;
    T_Jet_Energy = new std::vector<float>;
    T_Jet_Phi = new std::vector<float>;

    
}
void
ElecIdTreeProducer::endEvent()
{
    delete T_Event_Rho;
    
    delete T_Event_pathsFired;
    
    delete T_Elec_puChargedIso;
    delete T_Elec_allChargedHadronIso;
    delete T_Elec_chargedHadronIso;
    delete T_Elec_neutralHadronIso;
    delete T_Elec_photonIso;
    delete T_Elec_puChargedIso04;
    delete T_Elec_allChargedHadronIso04;
    delete T_Elec_chargedHadronIso04;
    delete T_Elec_neutralHadronIso04;
    delete T_Elec_photonIso04;
    
    
    delete T_Gen_Elec_Px;
    delete T_Gen_Elec_Py;
    delete T_Gen_Elec_Pz;
    delete T_Gen_Elec_Energy;
    delete T_Gen_Elec_PDGid;
    delete T_Gen_Elec_MotherID;
    delete T_Gen_Elec_GndMotherID;
    delete T_Gen_Elec_status;
    delete T_Gen_Elec_fromTAU;
    delete T_Gen_Elec_softElectron;
    
    delete T_Elec_TriggerLeg;
    
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
    
    delete T_Elec_noZSsee;
    delete T_Elec_noZSspp;
    delete T_Elec_noZSsep;
    delete T_Elec_noZSr9;
    delete T_Elec_noZSe1x5;
    delete T_Elec_noZSe2x5MaxSeed;
    delete T_Elec_noZSe5x5;

    
    delete T_Elec_ECALiso;
    delete T_Elec_HCALiso;
    delete T_Elec_TKiso;
    
    delete T_Elec_nbBC;
    delete T_Elec_BC1_eta;
    delete T_Elec_BC1_phi;
    delete T_Elec_BC1_energy;
    delete T_Elec_BC2_eta;
    delete T_Elec_BC2_phi;
    delete T_Elec_BC2_energy;
    delete T_Elec_BC3_eta;
    delete T_Elec_BC3_phi;
    delete T_Elec_BC3_energy;
    
    delete T_Muon_TriggerLeg;
    
    delete T_Muon_Eta;
    delete T_Muon_Phi;
    delete T_Muon_Energy;
    delete T_Muon_Et;
    delete T_Muon_Pt;
    delete T_Muon_Px;
    delete T_Muon_Py;
    delete T_Muon_Pz;
    delete T_Muon_Mass;
    delete T_Muon_IsGlobalMuon;
    delete T_Muon_IsTrackerMuon;
    delete T_Muon_IsPFMuon;
    delete T_Muon_IsCaloMuon;
    delete T_Muon_IsStandAloneMuon;
    delete T_Muon_IsMuon;
    delete T_Muon_IsGlobalMuon_PromptTight;
    delete T_Muon_IsTrackerMuonArbitrated;
    delete T_Muon_numberOfChambers;
    delete T_Muon_numberOfChambersRPC;
    delete T_Muon_numberOfMatches;
    delete T_Muon_numberOfMatchedStations;
    delete T_Muon_charge;
    delete T_Muon_TMLastStationTight;
    delete T_Muon_globalTrackChi2;
    delete T_Muon_validMuonHits;
    delete T_Muon_trkKink;
    delete T_Muon_trkNbOfTrackerLayers;
    delete T_Muon_trkNbOfValidTrackeHits;
    delete T_Muon_trkValidPixelHits;
    delete T_Muon_trkError;
    delete T_Muon_dB;
    delete T_Muon_dzPV;
    delete T_Muon_dBstop;
    delete T_Muon_dzstop;
    delete T_Muon_chargedHadronIsoR04;
    delete T_Muon_neutralHadronIsoR04;
    delete T_Muon_photonIsoR04;
    delete T_Muon_chargedHadronIsoPUR04;
    delete T_Muon_chargedHadronIsoR03;
    delete T_Muon_neutralHadronIsoR03;
    delete T_Muon_photonIsoR03;
    delete T_Muon_chargedHadronIsoPUR03;
    delete T_Muon_isoR03_emEt;
    delete T_Muon_isoR03_hadEt;
    delete T_Muon_isoR03_hoEt;
    delete T_Muon_isoR03_sumPt;
    delete T_Muon_isoR03_nTracks;
    delete T_Muon_isoR03_nJets;
    delete T_Muon_isoRingsMVA;
    
    delete T_Gen_Muon_Px;
    delete T_Gen_Muon_Py;
    delete T_Gen_Muon_Pz;
    delete T_Gen_Muon_Energy;
    delete T_Gen_Muon_status;
    delete T_Gen_Muon_PDGid;
    delete T_Gen_Muon_MotherID;
    delete T_Gen_Muon_GndMotherID;
    delete T_Gen_Muon_fromTAU;
    delete T_Gen_Muon_softMuon;
    
    delete T_Jet_Px;
    delete T_Jet_Py;
    delete T_Jet_Pz;
    delete T_Jet_Et;
    delete T_Jet_Eta;
    delete T_Jet_Energy;
    delete T_Jet_Phi;

    
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

bool
ElecIdTreeProducer::hasWZasMother(const reco::GenParticle  p)
{
    bool foundW = false;
    if (p.numberOfMothers()==0) return foundW;
    const reco::Candidate  *part = (p.mother());
    // loop on the mother particles to check if is has a W has mother
    while ((part->numberOfMothers()>0)) {
        const reco::Candidate  *MomPart =part->mother();
        if ((fabs(MomPart->pdgId())>=22)&&(fabs(MomPart->pdgId())<=24)){
            foundW = true;
            break;
        }
        part = MomPart;
    }
    return foundW;
}

bool
ElecIdTreeProducer::hasTauasMother(const reco::GenParticle  p)
{
    bool foundW = false;
    if (p.numberOfMothers()==0) return foundW;
    const reco::Candidate  *part = (p.mother());
    // loop on the mother particles to check if is has a W has mother
    while ((part->numberOfMothers()>0)) {
        const reco::Candidate  *MomPart =part->mother();
        if ((fabs(MomPart->pdgId())==15)){
            foundW = true;
            break;
        }
        part = MomPart;
    }
    return foundW;
}
bool
ElecIdTreeProducer::isNonPromptElectron(const reco::GenParticle  p)
{
    bool foundW = false;
    if (p.numberOfMothers()==0) return foundW;
    const reco::Candidate  *part = (p.mother());
    // loop on the mother particles to check if is has a W has mother
    while ((part->numberOfMothers()>0)&&(part->status()==1)) {
        const reco::Candidate  *MomPart =part->mother();
        part = MomPart;
    }
   // std::cout << "status=" << part->status()  << " pdgID=" << part->pdgId() << std::endl;
    if (( part->status()==2)&&(fabs(part->pdgId())>50)) foundW=1;
    return foundW;
}


bool
ElecIdTreeProducer::isMatchedWithGen(reco::GenParticle  p, const reco::GsfElectron &recoElec)
{
    float deltaR = sqrt(pow(recoElec.eta()-p.eta(),2)+ pow(acos(cos(recoElec.phi()-p.phi())),2)) ;
    if (deltaR<0.1) return true;
    return false;
}



//define this as a plug-in
DEFINE_FWK_MODULE(ElecIdTreeProducer);
