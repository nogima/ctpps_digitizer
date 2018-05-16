import FWCore.ParameterSet.Config as cms

# ---------- trigger data ----------
from EventFilter.CTPPSRawToDigi.totemTriggerRawToDigi_cfi import totemTriggerRawToDigi
totemTriggerRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")

# ---------- Si strips ----------
totemDAQMappingESSourceXML_TrackingStrip = cms.ESSource("TotemDAQMappingESSourceXML",
  verbosity = cms.untracked.uint32(0),
  subSystem = cms.untracked.string("TrackingStrip"),
  configuration = cms.VPSet(
    # 2017
    cms.PSet(
      validityRange = cms.EventRange("1:min - 999999999:max"),
      mappingFileNames = cms.vstring("CondFormats/CTPPSReadoutObjects/xml/mapping_tracking_strip_2017.xml"),
      maskFileNames = cms.vstring()
    )
  )
)

from EventFilter.CTPPSRawToDigi.totemRPRawToDigi_cfi import totemRPRawToDigi
totemRPRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")

# various error/warning/info output may be enabled with these flags
#  totemRPRawToDigi.RawUnpacking.verbosity = 1
#  totemRPRawToDigi.RawToDigi.verbosity = 1 # or higher number for more output
#  totemRPRawToDigi.RawToDigi.printErrorSummary = 1
#  totemRPRawToDigi.RawToDigi.printUnknownFrameSummary = 1

# ---------- diamonds ----------
totemDAQMappingESSourceXML_TimingDiamond = cms.ESSource("TotemDAQMappingESSourceXML",
  verbosity = cms.untracked.uint32(0),
  subSystem = cms.untracked.string("TimingDiamond"),
  configuration = cms.VPSet(
    # 2017
    cms.PSet(
      validityRange = cms.EventRange("1:min - 999999999:max"),
      mappingFileNames = cms.vstring("CondFormats/CTPPSReadoutObjects/xml/mapping_timing_diamond_2017.xml"),
      maskFileNames = cms.vstring()
    )
  )
)

#from EventFilter.CTPPSRawToDigi.ctppsDiamondRawToDigi_cfi import ctppsDiamondRawToDigi
#ctppsDiamondRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")

# ---------- pixels ----------

from EventFilter.CTPPSRawToDigi.ctppsPixelRawToDigi_cfi import ctppsPixelDigis
ctppsPixelDigis.InputLabel = cms.InputTag("ctppsPixelRawData")

# raw-to-digi sequence
ctppsRawToDigi = cms.Sequence(
  totemTriggerRawToDigi *
  totemRPRawToDigi *
  #ctppsDiamondRawToDigi*
  ctppsPixelDigis
)
