import FWCore.ParameterSet.Config as cms

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
                                      'file:/tmp/hbrun/theLocalReco.root'
                                      #'file:/tmp/hbrun/theDY_70_file.root'
    )
)


process.GlobalTag.globaltag = 'FT_R_70_V1::All'

process.ElecIdTreeProducer = cms.EDAnalyzer('ElecIdTreeProducer',
    isMC                        = cms.bool(False),
    electronsCollection       	= cms.InputTag("gedGsfElectrons","","reRECO"),
    primaryVertexInputTag   	= cms.InputTag("offlinePrimaryVertices","","reRECO"),
    rechitCollectionEB   	= cms.InputTag("reducedEcalRecHitsEB","","reRECO"),
    rechitCollectionEE   	= cms.InputTag("reducedEcalRecHitsEE","","reRECO"),
    conversionsCollection   = cms.InputTag("allConversions","","reRECO"),
    beamSpotInputTag   = cms.InputTag("offlineBeamSpot","","reRECO"),
    outputFile		        = cms.string("ElecIDtree.root")
)


process.p = cms.Path(process.ElecIdTreeProducer)
