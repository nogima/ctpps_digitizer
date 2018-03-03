import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'file:/afs/cern.ch/work/p/polme/public/FullSimPPS/test19/CMSSW_10_1_0_pre1/src/step3_RAW2DIGI_RECO.root'
    )
)

process.validation = cms.EDAnalyzer('CTPPSValidation',
	MCEvent = cms.untracked.InputTag("LHCTransport"),
        tagPixelTrack = cms.InputTag("ctppsPixelLocalTracks"),
        tagPixelRecHit = cms.InputTag("ctppsPixelRecHits"),
        tagPixelDigiBefRaw = cms.InputTag("ctppsPixelDigis"),
        tagPixelDigiAftRaw = cms.InputTag("RPixDetDigitizer"),
        psimHitTag = cms.untracked.InputTag('g4SimHits',"CTPPSPixelHits")
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.TFileService = cms.Service("TFileService",
                fileName = cms.string('outValidation10_LHCT21955.root')
)

process.p = cms.Path(process.validation)
