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
                                      #'file:/tmp/hbrun/theDY_70_file.root'
					'/store/relval/CMSSW_7_1_0/RelValZMM_13/GEN-SIM-RECO/POSTLS171_V15-v1/00000/6650F961-99FB-E311-BA90-0025905A48BC.root'                                      

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
    muonProducer 	         	= cms.VInputTag(cms.InputTag("muons")),
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
    pathsToSave           =cms.vstring("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v20",
                                       "HLT_Mu17_Mu8_v23",
                                       "HLT_Mu17_TkMu8_v15",
                                       "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",
                                       "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10"),
    filterToMatch           =cms.vstring("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",
                                         "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter",
                                         "hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8",
                                         "hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17",
                                         "hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17",
                                         "hltDiMuonGlbFiltered17TrkFiltered8",
                                         "hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8",
                                         "hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",
                                         "hltL1Mu12EG7L3MuFiltered17",
                                         "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter"),
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
