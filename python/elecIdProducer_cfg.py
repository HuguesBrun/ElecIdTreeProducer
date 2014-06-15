import FWCore.ParameterSet.Config as cms


isMC = True

process = cms.Process("runAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration/StandardSequences/MagneticField_38T_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Services_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.MessageLogger.cerr.FwkReport.reportEvery = 10


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
                                      #'file:/tmp/hbrun/theLocalReco.root'
                                      'file:/tmp/hbrun/theDY_70_file.root'
                                      #'file:/tmp/hbrun/theDYfile_new.root'
    )
)


process.GlobalTag.globaltag = 'FT_R_70_V1::All'
if (isMC):
    process.GlobalTag.globaltag = 'POSTLS170_V5::All'


typeProcess = "reRECO"
if (isMC):
    typeProcess = "RECO"

process.ElecIdTreeProducer = cms.EDAnalyzer('ElecIdTreeProducer',
    isMC                        = cms.bool(False),
    electronsCollection       	= cms.InputTag("gedGsfElectrons","",typeProcess),
    primaryVertexInputTag   	= cms.InputTag("offlinePrimaryVertices","",typeProcess),
    rechitCollectionEB   	= cms.InputTag("reducedEcalRecHitsEB","",typeProcess),
    rechitCollectionEE   	= cms.InputTag("reducedEcalRecHitsEE","",typeProcess),
    conversionsCollection   = cms.InputTag("allConversions","",typeProcess),
    beamSpotInputTag   = cms.InputTag("offlineBeamSpot","",typeProcess),
    rhoTags =               cms.VInputTag(cms.InputTag("kt6PFJetsForIsolation","rho","runAnalyzer")),
    metTag     = cms.InputTag("pfMet", "", typeProcess),
    jetCollectionTag     = cms.InputTag("ak5PFJets", "", typeProcess),
    triggerResultTag     = cms.InputTag("TriggerResults", "", "HLT"),
    triggerSummaryTag    = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    pathsToSave           =cms.vstring("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
                                       "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
                                       "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
                                       "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v",
                                       "HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v"),
    filterToMatch           =cms.vstring("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",
                                         "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter",
                                         "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",
                                         "hltEle8TightIdLooseIsoTrackIsoFilter",
                                         "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter",
                                         "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter",
                                         "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter",
                                         "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter"),
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


#process.p = cms.Path(process.triggerResultsFilter * process.primaryVertexFilter * process.noscraping * process.kt6PFJetsForIsolation*process.ElecIdTreeProducer)
process.p = cms.Path(process.primaryVertexFilter * process.noscraping * process.kt6PFJetsForIsolation*process.ElecIdTreeProducer)
