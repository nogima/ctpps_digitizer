import FWCore.ParameterSet.Config as cms

from Geometry.VeryForwardGeometry.geometryCTPPSFromDD_2018_cfi import *

from RecoCTPPS.Configuration.recoCTPPS_sequences_cff import *

recoCTPPS = cms.Sequence(recoCTPPSdets)

