import FWCore.ParameterSet.Config as cms


isMC = True

process = cms.Process("runAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration/StandardSequences/MagneticField_38T_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Services_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.MessageLogger.cerr.FwkReport.reportEvery = 10


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
                                      'file:/tmp/hbrun/theDYJetFile.root'
    )
)


process.GlobalTag.globaltag = 'FT_R_70_V1::All'
if (isMC):
    process.GlobalTag.globaltag = 'PHYS14_25_V1::All'


typeProcess = "reRECO"
if (isMC):
    typeProcess = "RECO"

process.ElecIdTreeProducer = cms.EDAnalyzer('ElecIdTreeProducer',
    isMC                        = cms.bool(False),
    doMuon                      = cms.bool(False),
    electronsCollection       	= cms.InputTag("gedGsfElectrons","",typeProcess),
    elecIdName                  = cms.VInputTag(cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-CSA14-PU20bx25-V0-standalone-veto"),
                                                cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-CSA14-PU20bx25-V0-standalone-loose"),
                                                cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-CSA14-PU20bx25-V0-standalone-medium"),
                                                cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight")),
    muonProducer 	         	= cms.VInputTag(cms.InputTag("muons")),
    primaryVertexInputTag   	= cms.InputTag("offlinePrimaryVertices","",typeProcess),
    rechitCollectionEB   	= cms.InputTag("reducedEcalRecHitsEB","",typeProcess),
    rechitCollectionEE   	= cms.InputTag("reducedEcalRecHitsEE","",typeProcess),
    conversionsCollection   = cms.InputTag("allConversions","",typeProcess),
    beamSpotInputTag   = cms.InputTag("offlineBeamSpot","",typeProcess),
    rhoTags =               cms.VInputTag(cms.InputTag("kt6PFJetsForIsolation","rho","runAnalyzer")),#,cms.InputTag("fixedGridRhoAll","","RECO"), cms.InputTag("fixedGridRhoFastjetAll","","RECO"), cms.InputTag("fixedGridRhoFastjetAllCalo","","RECO"), cms.InputTag("fixedGridRhoFastjetCentralCalo","","RECO"), cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp","","RECO"), cms.InputTag("fixedGridRhoFastjetCentralNeutral","","RECO")),
    metTag     = cms.InputTag("pfMet", "", typeProcess),
    jetCollectionTag     = cms.InputTag("ak5PFJets", "", typeProcess),
    triggerResultTag     = cms.InputTag("TriggerResults", "", "HLT"),
    triggerSummaryTag    = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                                            pathsToSave           =cms.vstring("HLT_Ele27_eta2p1_WP85_Gsf_v1",
                                                                               "HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1",
                                                                               "HLT_IsoMu20_eta2p1_IterTrk02_v1",
                                                                               "HLT_IsoMu24_eta2p1_IterTrk02_v1",
                                                                               "HLT_IsoTkMu20_eta2p1_IterTrk02_v1",
                                                                               "HLT_IsoTkMu24_eta2p1_IterTrk02_v1",
                                                                               "HLT_Mu17_Mu8_v1",
                                                                               "HLT_Mu17_TkMu8_v1",
                                                                               "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1",
                                                                               "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1",
                                                                               "HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1",
                                                                               "HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1"),
                                            filterToMatch           =cms.vstring("hltEle27WP85GsfTrackIsoFilter",
                                                                                 "hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter",
                                                                                 "hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter",
                                                                                 "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f20QL3crIsoRhoFiltered0p15IterTrk02",
                                                                                 "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15IterTrk02",
                                                                                 "hltL3fL1sMu16L1Eta2p1f0TkFiltered20QL3crIsoRhoFiltered0p15IterTrk02",
                                                                                 "hltL3fL1sMu16L1Eta2p1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02",
                                                                                 "hltDiMuonGlb17Glb8DzFiltered0p2",
                                                                                 "hltDiMuonGlb17Trk8DzFiltered0p2",
                                                                                 "hltDiMuonGlb17Glb8DzFiltered0p2RelTrkIsoFiltered0p4",
                                                                                 "hltDiMuonGlb17Trk8DzFiltered0p2RelTrkIsoFiltered0p4",
                                                                                 "hltL1sL1Mu5EG20ORL1Mu5IsoEG18L3IsoFiltered8",
                                                                                 "hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter",
                                                                                 "hltL1Mu12EG7L3IsoMuFiltered23",
                                                                                 "hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter"),
    HLTprocess            = cms.string("HLT"),
    outputFile		        = cms.string("ElecIDtree.root")
)


# to compute FastJet rho to correct isolation (note: EtaMax restricted to 2.5)
# all the analyses
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                           )


process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )


#keep only events passing the single Electron path
process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")
process.triggerResultsFilter.triggerConditions = cms.vstring("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
                                                             "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
                                                             "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
                                                             "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*",
                                                             "HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*")
process.triggerResultsFilter.l1tResults = ''
process.triggerResultsFilter.throw = False
process.triggerResultsFilter.hltResults = cms.InputTag( "TriggerResults", "", "HLT" )



if (isMC):
    process.ElecIdTreeProducer.isMC = cms.bool(True)



process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cff")


process.electronIDValueMapProducer.ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB")
process.electronIDValueMapProducer.eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE")
process.electronIDValueMapProducer.esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES")



#process.p = cms.Path(process.triggerResultsFilter * process.primaryVertexFilter * process.noscraping * process.kt6PFJetsForIsolation*process.ElecIdTreeProducer)
#process.p = cms.Path(process.primaryVertexFilter * process.noscraping * process.kt6PFJetsForIsolation*process.ElecIdTreeProducer)
process.p = cms.Path(process.kt6PFJetsForIsolation*process.egmGsfElectronIDSequence*process.ElecIdTreeProducer)
#process.p = cms.Path(process.kt6PFJetsForIsolation*process.ElecIdTreeProducer)





