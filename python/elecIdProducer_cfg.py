import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/tmp/hbrun/theLocalReco.root'
    )
)

process.ElecIdTreeProducer = cms.EDAnalyzer('ElecIdTreeProducer',
    isMC                        = cms.bool(False),
    electronsCollection       	= cms.InputTag("gedGsfElectrons","","reRECO"),
    primaryVertexInputTag   	= cms.InputTag("offlinePrimaryVertices","","reRECO"),
    outputFile		        = cms.string("ElecIDtree.root")
)


process.p = cms.Path(process.ElecIdTreeProducer)
