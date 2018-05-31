# Auto generated configuration file
# using: 
# Revision: 1.19 
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process('RECO',eras.Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Geometry.VeryForwardGeometry.geometryCTPPSFromDD_2016_cfi')


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
    fileNames = cms.untracked.vstring('file:step2_DIGI_DIGI2RAW2016.root'),
    secondaryFileNames = cms.untracked.vstring()
)

# Track memory leaks
process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Output definition

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('step3_RAW2DIGI_RECO2016.root'),
    outputCommands = cms.untracked.vstring("drop *","keep PSimHits*_*_*_*","keep CTPPS*_*_*_*","keep *_*RP*_*_*",'keep *_LHCTransport_*_*')
)


# Additional output definition
# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


process.load('CondFormats.CTPPSReadoutObjects.totemDAQMappingESSourceXML_cfi')
process.totemDAQMappingESSourceXML.configuration = cms.VPSet(
    cms.PSet(
      validityRange = cms.EventRange("1:min - 999999999:max"),
      mappingFileNames = cms.vstring("CondFormats/CTPPSReadoutObjects/xml/mapping_tracking_strip_2016_from_fill_5330.xml"),
      maskFileNames = cms.vstring()
    )
)


# modify CTPPS 2018 raw-to-digi modules  
process.load("EventFilter.CTPPSRawToDigi.pgctppsRawToDigi_cfi")     
from EventFilter.CTPPSRawToDigi.pgctppsRawToDigi_cfi import *
from Configuration.Eras.Modifier_ctpps_2016_cff import ctpps_2016
ctpps_2016.toReplaceWith(ctppsRawToDigi, ctppsRawToDigi.copyAndExclude([ctppsDiamondRawToDigi]))
ctpps_2016.toReplaceWith(ctppsRawToDigi, ctppsRawToDigi.copyAndExclude([ctppsPixelDigis]))

from RecoCTPPS.Configuration.recoCTPPS_DD_2016_cff import *
process.load("RecoCTPPS.Configuration.recoCTPPS_DD_2016_cff")
               
# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.ctppsRawToDigi) 
process.reco_step = cms.Path(process.recoCTPPS)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.output_step = cms.EndPath(process.output)


# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reco_step,process.endjob_step,process.output_step)

# filter all path with the production filter sequence
for path in process.paths:
    #  getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq
    getattr(process,path)._seq = getattr(process,path)._seq

