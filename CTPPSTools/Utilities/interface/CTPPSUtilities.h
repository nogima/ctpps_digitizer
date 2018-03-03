#ifndef CTPPSTOOLS_UTILITIES
#define CTPPSTOOLS_UTILITIES
#include <cmath>
#include <string>
#include "H_BeamParticle.h"
#include <CLHEP/Vector/LorentzVector.h>
#include "CLHEP/Units/PhysicalConstants.h"

namespace CTPPSTools {

typedef CLHEP::HepLorentzVector LorentzVector;

double fCrossingAngleBeam1;
double fCrossingAngleBeam2;
double fBeamMomentum;
double fBeamEnergy;
const double ProtonMass = 0.93827;
const double ProtonMassSQ = pow(ProtonMass,2);
const double   urad     = 1./1000000.; 

CLHEP::HepLorentzVector HectorParticle2LorentzVector(H_BeamParticle hp,int =1);

H_BeamParticle LorentzVector2HectorParticle(CLHEP::HepLorentzVector p);

void LorentzBoost(H_BeamParticle& h_p, const string& frame);

void LorentzBoost(CLHEP::HepLorentzVector& p_out, const std::string& frame);

void Get_t_and_xi(const CLHEP::HepLorentzVector* proton,double& t,double& xi) ;

};
#endif
