#ifndef SimPPS_CTPPSPixelDigiProducer_RPix_Charge_Share_H
#define SimPPS_CTPPSPixelDigiProducer_RPix_Charge_Share_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <vector>
#include "Geometry/VeryForwardGeometry/interface/CTPPSPixelSimTopology.h"
#include "SimCTPPS/CTPPSPixelDigiProducer/interface/RPixSignalPoint.h"
//class TFile;
//class TH2D;
class RPixChargeShare
{
public:
  RPixChargeShare(const edm::ParameterSet &params, uint32_t det_id);
//    SimRP::strip_charge_map 
  std::map<unsigned short, double, std::less<unsigned short> > Share(const std::vector<RPixSignalPoint> &charge_map);
private:
  uint32_t det_id_;
  std::vector<double> signalCoupling_;
  std::map<unsigned short, double, std::less<unsigned short> > thePixelChargeMap;
  CTPPSPixelSimTopology theRPixDetTopology_;
  double sqrt_2;
  double chargeMap2E[4][60][40];
  int no_of_pixels_;
  int verbosity_;
  std::string ChargeMapFile2E_[4];
};

#endif
