import FWCore.ParameterSet.Config as cms


isMC = True

process = cms.Process("runAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration/StandardSequences/MagneticField_38T_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Services_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
                                      #'file:/tmp/hbrun/theLocalReco.root'
                                      #'file:/tmp/hbrun/theDY_70_file.root'
                                      'file:/tmp/hbrun/theDYfile_new.root'
    )
)


process.GlobalTag.globaltag = 'FT_R_70_V1::All'

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
    rhoTags =               cms.VInputTag(cms.InputTag("ak5CaloJets","rho",""),
                                          cms.InputTag("ak5PFJets","rho","")),
    triggerResultTag     = cms.InputTag("TriggerResults", "", "HLT"),
    triggerSummaryTag    = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    pathsToSave           =cms.vstring("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
                                       "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v"),
    HLTprocess            = cms.string("HLT"),
    outputFile		        = cms.string("ElecIDtree.root")
)








if (isMC):
    process.ElecIdTreeProducer.isMC = cms.bool(True)


process.p = cms.Path(process.ElecIdTreeProducer)
