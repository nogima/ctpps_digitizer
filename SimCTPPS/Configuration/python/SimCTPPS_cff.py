import FWCore.ParameterSet.Config as cms

# CTPPS Digitization
from SimCTPPS.CTPPSPixelDigiProducer.RPixDetConf_cfi import *
from SimCTPPS.RPDigiProducer.RPSiDetConf_cfi import *

ctppsDigi = cms.Sequence(RPixDetDigitizer+RPSiDetDigitizer)

# add CTPPS 2016 digi modules
from Configuration.Eras.Modifier_ctpps_2016_cff import ctpps_2016

_ctpps_2016_Digi = ctppsDigi.copy() 
_ctpps_2016_Digi =cms.Sequence(RPSiDetDigitizer)
ctpps_2016.toReplaceWith(ctppsDigi,_ctpps_2016_Digi)

# add CTPPS 2017 digi modules
from Configuration.Eras.Modifier_ctpps_2017_cff import ctpps_2017

_ctpps_2017_Digi = ctppsDigi.copy() 
_ctpps_2017_Digi = cms.Sequence(RPixDetDigitizer+RPSiDetDigitizer)
ctpps_2017.toReplaceWith(ctppsDigi,_ctpps_2017_Digi)

# add CTPPS 2018 digi modules
from Configuration.Eras.Modifier_ctpps_2018_cff import ctpps_2018

_ctpps_2018_Digi = ctppsDigi.copy() 
_ctpps_2018_Digi = cms.Sequence(RPixDetDigitizer)
ctpps_2018.toReplaceWith(ctppsDigi,_ctpps_2018_Digi)


