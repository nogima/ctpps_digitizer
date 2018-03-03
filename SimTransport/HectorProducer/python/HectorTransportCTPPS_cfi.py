import FWCore.ParameterSet.Config as cms

from SimG4Core.Application.hectorParameter_cfi import *

from SimTransport.HectorProducer.CTPPSOpticsParameters_cfi import *

LHCTransport = cms.EDProducer("CTPPSHectorProducer",
    HepMCProductLabel = cms.string('generatorSmeared'),  ## HepMC source to be processed
    CTPPSTransport = cms.bool(True), 
    Verbosity = cms.bool(False),
    CTPPSHector = cms.PSet(
        HectorEtaCut,
        Nominal_TDR, #Validated_PreTS2_2016, #Nominal_2017_beta40cm,
        BeamLineLengthCTPPS = cms.double(250.0),
        CTPPSf = cms.double(212.55),    ##in meters
        CTPPSb = cms.double(212.55),    ##in meters
        smearEnergy = cms.bool(True),       ## if False: no Energy smearing(i.e. sigmaEnergy =0.0)
        smearAng = cms.bool(True),       ## if False: no Angle smearing(i.e. sigmaSTX(Y) =0.0)
        CrossAngleCorr = cms.bool(True),
        MomentumMin = cms.double(3.000),	
    )
)
