# Auto generated configuration file
# using: 
# Revision: 1.19 
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process('RECO',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("CondFormats.CTPPSReadoutObjects.CTPPSPixelDAQMappingESSourceXML_cfi")
process.load('Configuration.StandardSequences.CTPPSRawToDigi_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("EmptyIOVSource",
    timetype = cms.string('runnumber'),
    firstValue = cms.uint64(1),
    lastValue = cms.uint64(1),
    interval = cms.uint64(1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step2_DIGI_DIGI2RAW.root'),
    secondaryFileNames = cms.untracked.vstring()
)

# Track memory leaks
process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Output definition

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('step3_RAW2DIGI_RECO.root'),
    outputCommands = cms.untracked.vstring("drop *","keep PSimHits*_*_*_*","keep CTPPS*_*_*_*","keep *_RP*_*_*",'keep *_LHCTransport_*_*')
)


# Additional output definition
# Other statements

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')

process.load("CondCore.CondDB.CondDB_cfi")

process.CondDB.connect = 'sqlite_file:CTPPSPixel_Mapping_Mask_Gain_Jan23rd2018.db'

process.PoolDBESSource = cms.ESSource("PoolDBESSource",
        process.CondDB,
        DumpStat = cms.untracked.bool(True),
        toGet = cms.VPSet(
          cms.PSet(
            record = cms.string('CTPPSPixelGainCalibrationsRcd'),
            tag = cms.string("CTPPSPixelGainCalibNew_v4")
          ),
          cms.PSet(
            record = cms.string('CTPPSPixelAnalysisMaskRcd'),
            tag = cms.string("PixelAnalysisMask")
          ),
          cms.PSet(
            record = cms.string('CTPPSPixelDAQMappingRcd'),
            tag = cms.string("PixelDAQMapping")
          )
        )
)

from RecoCTPPS.PixelLocal.ctppsPixelLocalReconstruction_cff import *
process.load("RecoCTPPS.PixelLocal.ctppsPixelLocalReconstruction_cff")

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reco_step = cms.Path(process.ctppsPixelLocalReconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.output_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reco_step,process.endjob_step,process.output_step)

# filter all path with the production filter sequence
for path in process.paths:
    #  getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq
    getattr(process,path)._seq = getattr(process,path)._seq

