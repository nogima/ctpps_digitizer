#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "SimTransport/HectorProducer/interface/CTPPSHector.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "HepMC/SimpleVector.h"

#include "CLHEP/Random/RandGauss.h"
#include "CTPPSTools/Utilities/interface/CTPPSUtilities.h"

#include "TRandom3.h"
#include <TMath.h>

#include "H_Parameters.h"
#include "H_RecRPObject.h"

#include <math.h>
#include <iomanip>
#include <cstdlib>


CTPPSHector::CTPPSHector(const edm::ParameterSet & param, bool verbosity,bool CTPPSTransport) : 
    m_smearAng(false),m_sig_e(0.),m_smearE(false),m_sigmaSTX(0.),m_sigmaSTY(0.),m_sigmaSX(0.),m_sigmaSY(0.),
    fCrossAngleCorr(false),fCrossingAngleBeam1(0.),fCrossingAngleBeam2(0.),fBeamMomentum(0),fBeamEnergy(0),
    fVtxMeanX(0.),fVtxMeanY(0.),fVtxMeanZ(0.),fMomentumMin(0.),
    m_verbosity(verbosity), 
    m_CTPPSTransport(CTPPSTransport),NEvent(0)

{
    // Create LHC beam line
    edm::ParameterSet hector_par = param.getParameter<edm::ParameterSet>("CTPPSHector");

    // User definitons
    lengthctpps     = hector_par.getParameter<double>("BeamLineLengthCTPPS" );
    m_f_ctpps_f     = (float) hector_par.getParameter<double>("CTPPSf");
    m_b_ctpps_b     = (float) hector_par.getParameter<double>("CTPPSb");

    beam1filename   = hector_par.getParameter<string>("Beam1");
    beam2filename   = hector_par.getParameter<string>("Beam2");  
    m_smearAng      = hector_par.getParameter<bool>("smearAng");
    m_sigmaSTX      = hector_par.getParameter<double>("sigmaSTX" );
    m_sigmaSTY      = hector_par.getParameter<double>("sigmaSTY" );
    m_sigmaSX       = hector_par.getParameter<double>("sigmaSX");
    m_sigmaSY       = hector_par.getParameter<double>("sigmaSY");
    m_smearE        = hector_par.getParameter<bool>("smearEnergy");
    m_sig_e         = hector_par.getParameter<double>("sigmaEnergy");
    etacut          = hector_par.getParameter<double>("EtaCutForHector" );
    //CTPPS
    fCrossAngleCorr = hector_par.getParameter<bool>("CrossAngleCorr");
    fCrossingAngleBeam1  = hector_par.getParameter<double>("CrossingAngleBeam1");
    fCrossingAngleBeam2  = hector_par.getParameter<double>("CrossingAngleBeam2");
    fBeamEnergy     = hector_par.getParameter<double>("BeamEnergy"); // beam energy in GeV
    fVtxMeanX       = hector_par.getParameter<double>("VtxMeanX");
    fVtxMeanY       = hector_par.getParameter<double>("VtxMeanY");
    fVtxMeanZ       = hector_par.getParameter<double>("VtxMeanZ");
    fMomentumMin    = hector_par.getParameter<double>("MomentumMin"); 
    fBeamXatIP      = hector_par.getParameter<double>("BeamXatIP"); // position in mm
    fBeamYatIP      = hector_par.getParameter<double>("BeamYatIP"); // position in mm

    fBeamMomentum = sqrt(fBeamEnergy*fBeamEnergy - ProtonMassSQ);

    theCorrespondenceMap.clear();
    
    CTPPSTools::fBeamMomentum=fBeamMomentum;
    CTPPSTools::fBeamEnergy=fBeamEnergy;
    CTPPSTools::fCrossingAngleBeam1=-fCrossingAngleBeam1;
    CTPPSTools::fCrossingAngleBeam2=-fCrossingAngleBeam2;

    // construct beam line for CTPPS (forward 1 backward 2):                                                                                           
    if(!SetBeamLine()) {
        if ( m_verbosity ) LogDebug("CTPPSHectorSetup") << "CTPPSHector: WARNING: lengthctpps=  " << lengthctpps;
    } 
}

CTPPSHector::~CTPPSHector(){

    for (std::map<unsigned int,H_BeamParticle*>::iterator it = m_beamPart.begin(); it != m_beamPart.end(); ++it ) {
        delete (*it).second;
    }

}

void CTPPSHector::clearApertureFlags(){
    m_isStoppedctpps.clear();
}

void CTPPSHector::clear(){
    for ( std::map<unsigned int,H_BeamParticle*>::iterator it = m_beamPart.begin(); it != m_beamPart.end(); ++it ) {
        delete (*it).second;
    };
    m_beamPart.clear();
    m_direct.clear();
    m_eta.clear();
    m_pdg.clear();
    m_pz.clear();
    m_isCharged.clear();  
}

void CTPPSHector::add( const HepMC::GenEvent * evt ,const edm::EventSetup & iSetup, CLHEP::HepRandomEngine * engine) {

    H_BeamParticle* h_p  = NULL;
    unsigned int line;

    for (HepMC::GenEvent::particle_const_iterator eventParticle =evt->particles_begin();
            eventParticle != evt->particles_end();
            ++eventParticle ) {
        if ( (*eventParticle)->status() == 1 && (*eventParticle)->pdg_id()==2212 ){
            if ( abs( (*eventParticle)->momentum().eta())>etacut && abs( (*eventParticle)->momentum().pz())>fMomentumMin){
                line = (*eventParticle)->barcode();
                if ( m_beamPart.find(line) == m_beamPart.end() ) {
                    double charge=1.;
                    m_isCharged[line] = false;// neutrals
                    HepMC::GenParticle * g = (*eventParticle);	
                    iSetup.getData( pdt );
                    const ParticleData * part = pdt->particle( g->pdg_id() );
                    if (part){
                        charge = part->charge();
                    }
                    if(charge !=0) m_isCharged[line] = true;//charged
                    double mass = (*eventParticle)->generatedMass();

                    h_p = new H_BeamParticle(mass,charge);

                    double px,py,pz,e;
                    double TXforPosition=0.0, TYforPosition=0.0;//urad

                    px = (*eventParticle)->momentum().px();	  
                    py = (*eventParticle)->momentum().py();	  
                    pz = (*eventParticle)->momentum().pz();	  

                    e = sqrt(pow(mass,2)+pow(px,2)+pow(py,2)+pow(pz,2));

                    TXforPosition=0.;
                    TYforPosition=0.;

                    // Apply Beam and Crossing Angle Corrections
                    LorentzVector p_out(px,py,pz,e);
                    ApplyBeamCorrection(p_out, engine);
                    CTPPSTools::LorentzBoost(p_out,"LAB");

                    // from mm to cm        
                    double XforPosition = (*eventParticle)->production_vertex()->position().x()/cm;//cm
                    double YforPosition = (*eventParticle)->production_vertex()->position().y()/cm;//cm
                    double ZforPosition = (*eventParticle)->production_vertex()->position().z()/cm;//cm

                    if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") <<
                                    " fVtxMeanX: " << fVtxMeanX << " fVtxMeanY: " << fVtxMeanY << " fVtxMeanZ: "  << fVtxMeanZ ;

                    h_p->set4Momentum(-p_out.px(), p_out.py(), abs(p_out.pz()), p_out.e());
                    h_p->setPosition(-((XforPosition-fVtxMeanX)*cm_to_um+fBeamXatIP*mm_to_um),(YforPosition-fVtxMeanY)*cm_to_um+fBeamYatIP*mm_to_um,
                                        h_p->getTX()+TXforPosition,h_p->getTY()+TYforPosition,-(ZforPosition-fVtxMeanZ)*cm_to_m);
                    
                    m_beamPart[line] = h_p;
                    m_direct[line] = 0;
                    m_direct[line] = ( pz > 0 ) ? 1 : -1;
                    m_eta[line] = p_out.eta();
                    m_pdg[line] = (*eventParticle)->pdg_id();
                    m_pz[line]  = p_out.pz();
                    if(m_verbosity) { 
                        LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector:add: barcode = " << line 
                            << " status = " << g->status() 
                            << " PDG Id = " << g->pdg_id() 
                            << " mass = " << mass 
                            << " pz = " << pz 
                            << " charge = " << charge 
                            << " m_isCharged[line] = " << m_isCharged[line];
                    } 
                }// if find line
            }// if eta > 8.2
        }// if status
    }// for loop

}

void CTPPSHector::filterCTPPS(TRandom3* rootEngine){

    unsigned int line;
    H_BeamParticle * part = NULL;
 
    std::map< unsigned int, H_BeamParticle* >::iterator it;

    bool is_stop;
    int direction;

    float x1_ctpps;
    float y1_ctpps;

    if ( m_beamPart.size() && lengthctpps>0. ) {

        for (it = m_beamPart.begin(); it != m_beamPart.end(); ++it ) {
            line = (*it).first;
            part = (*it).second;

            if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector:filterCTPPS: barcode = " << line;
            if ( (*m_isCharged.find( line )).second ) {
                direction = (*m_direct.find( line )).second;
                if ( direction == 1 && m_beamlineCTPPS1 != 0 ) {
                    
 		    part->computePath(&*m_beamlineCTPPS1);

                    is_stop = part->stopped(&* m_beamlineCTPPS1 );
                    if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << 
                                    "CTPPSHector:filterCTPPS: barcode = " << line << " positive is_stop=  "<< is_stop;
                }
                else if ( direction == -1 && m_beamlineCTPPS2 != 0 ){

                    part->computePath(&*m_beamlineCTPPS2);

                    is_stop = part->stopped(&*m_beamlineCTPPS2 );
                    if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << 
                                    "CTPPSHector:filterCTPPS: barcode = " << line << " negative is_stop=  "<< is_stop;
                }
                else {
                    is_stop = true;
                    if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector:filterCTPPS: barcode = " << line << " 0      is_stop=  "<< is_stop;
                }

                //propagating
                m_isStoppedctpps[line] = is_stop;
                if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << 
                                "CTPPSHector:filterCTPPS: barcode = " << line << " isStopped=" << (*m_isStoppedctpps.find(line)).second;

                if (!is_stop) {
                    if ( direction == 1 ) part->propagate( m_f_ctpps_f ); 
                    if ( direction == -1 ) part->propagate( m_b_ctpps_b );  
                    x1_ctpps = -part->getX()/millimeter;
                    y1_ctpps = part->getY()/millimeter;
                    if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << 
                                   "CTPPSHector:filterCTPPS: barcode = " << line << " x=  "<< x1_ctpps <<" y= " << y1_ctpps;
                    m_xAtTrPoint[line]  = x1_ctpps;
                    m_yAtTrPoint[line]  = y1_ctpps;
                    m_TxAtTrPoint[line] = -part->getTX(); // needs to be reflected due to the way phi is calculated here
                    m_TyAtTrPoint[line] = part->getTY();
                    m_eAtTrPoint[line]  = part->getE();

                }
            }// if isCharged
            else {
                m_isStoppedctpps[line] = true;// imply that neutral particles stopped to reach 420m
                if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << 
                                "CTPPSHector:filterCTPPS: barcode = " << line << " isStopped=" << (*m_isStoppedctpps.find(line)).second;
            }

        } // for (it = m_beamPart.begin(); it != m_beamPart.end(); it++ ) 
    } // if ( m_beamPart.size() )

}//

int  CTPPSHector::getDirect( unsigned int part_n ) const {
    std::map<unsigned int, int>::const_iterator it = m_direct.find( part_n );
    if ( it != m_direct.end() ){
        return (*it).second;
    }
    return 0;
}

void CTPPSHector::print() const {
    for (std::map<unsigned int,H_BeamParticle*>::const_iterator it = m_beamPart.begin(); it != m_beamPart.end(); ++it ) {
        (*it).second->printProperties();
    };
}

void CTPPSHector::ApplyBeamCorrection(LorentzVector& p_out, CLHEP::HepRandomEngine* engine)
{

    double theta  = p_out.theta();
    double thetax = atan(p_out.px()/abs(p_out.pz()));
    double thetay = atan(p_out.py()/abs(p_out.pz()));
    double energy = p_out.e();

    int direction = (p_out.pz()>0)?1:-1;

    if (p_out.pz()<0) theta=CLHEP::pi-theta;

    double dtheta_x = (double)(m_smearAng)?CLHEP::RandGauss::shoot(engine,0.,m_sigmaSTX):0;
    double dtheta_y = (double)(m_smearAng)?CLHEP::RandGauss::shoot(engine,0.,m_sigmaSTY):0;
    double denergy  = (double)(m_smearE)?CLHEP::RandGauss::shoot(engine,0.,m_sig_e):0.;

    double s_theta = sqrt(pow(thetax+dtheta_x*urad,2)+pow(thetay+dtheta_y*urad,2)); 
    double s_phi = atan2(thetay+dtheta_y*urad,thetax+dtheta_x*urad);
    energy+=denergy;
    double p = sqrt(pow(energy,2)-ProtonMassSQ);

    p_out.setPx((double)p*sin(s_theta)*cos(s_phi));
    p_out.setPy((double)p*sin(s_theta)*sin(s_phi));
    p_out.setPz((double)p*(cos(s_theta))*direction);
    p_out.setE(energy);
}

HepMC::GenEvent * CTPPSHector::addPartToHepMC( HepMC::GenEvent * evt ){
    NEvent++;
    theCorrespondenceMap.clear();

    unsigned int line;

    HepMC::GenParticle * gpart;
    double tx,ty,theta,fi,energy,time = 0;
    std::map< unsigned int, H_BeamParticle* >::iterator it;

    for (it = m_beamPart.begin(); it != m_beamPart.end(); ++it ) {
        line = (*it).first;
        if(!m_CTPPSTransport) m_isStoppedctpps[line] = true;
        if(m_verbosity) {
            LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector:addPartToHepMC: barcode = " << line << "\n"
                << "CTPPSHector:addPartToHepMC: isStoppedctpps=" << (*m_isStoppedctpps.find(line)).second;
        }
        if (!((*m_isStoppedctpps.find(line)).second)){

            gpart = evt->barcode_to_particle( line );
            if ( gpart ) {
                tx     = (*m_TxAtTrPoint.find(line)).second / 1000000.;
                ty     = (*m_TyAtTrPoint.find(line)).second / 1000000.;
                theta  = sqrt((tx*tx) + (ty*ty));
                double ddd = 0.;
                double fi_  = 0.; 
                if ( !((*m_isStoppedctpps.find(line)).second)) {
                    if( (*m_direct.find( line )).second >0 ) {
                        ddd = m_f_ctpps_f;
                        fi_    = std::atan2(tx,ty); // tx, ty never == 0?
                    }
                    else if((*m_direct.find( line )).second <0 ) {
                        ddd = m_b_ctpps_b;
                        theta= CLHEP::pi-theta;
                        fi_    = std::atan2(tx,ty); // tx, ty never == 0?
                    }
                } 
                fi = fi_; 
                energy = (*m_eAtTrPoint.find(line)).second;
                time = ( ddd*meter - gpart->production_vertex()->position().z()*mm ); // mm

                if(ddd != 0.) {
                    if(m_verbosity) {
                        LogDebug("CTPPSHectorEventProcessing") <<"CTPPSHector:: x= "<< (*(m_xAtTrPoint.find(line))).second*0.001<< "\n"
                            <<"CTPPSHector:: y= "<< (*(m_yAtTrPoint.find(line))).second*0.001<< "\n"
                            <<"CTPPSHector:: z= "<< ddd * (*(m_direct.find( line ))).second*1000. << "\n"
                            <<"CTPPSHector:: t= "<< time;
                    }

                    HepMC::GenVertex * vert = new HepMC::GenVertex( HepMC::FourVector( (*(m_xAtTrPoint.find(line))).second*0.001,
                                (*(m_yAtTrPoint.find(line))).second*0.001,
                                ddd * (*(m_direct.find( line ))).second*1000.,
                                time + .001*time ) );

                    gpart->set_status( 2 );
                    vert->add_particle_in( gpart );
                    double Pmom = sqrt(energy*energy - ProtonMassSQ);
                    vert->add_particle_out( new HepMC::GenParticle( HepMC::FourVector(Pmom*std::sin(theta)*std::sin(fi),
                                    Pmom*std::sin(theta)*std::cos(fi),
                                    Pmom*std::cos(theta),
                                    energy ),gpart->pdg_id(), 1, gpart->flow() ) );
                    evt->add_vertex( vert );

                    int ingoing = (*vert->particles_in_const_begin())->barcode();
                    int outgoing = (*vert->particles_out_const_begin())->barcode();
                    LHCTransportLink theLink(ingoing,outgoing);
                    if (m_verbosity) LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector:addPartToHepMC: LHCTransportLink " << theLink;
                    theCorrespondenceMap.push_back(theLink);

                    if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector::TRANSPORTED pz= " << gpart->momentum().pz()  
                        << " eta= "<< gpart->momentum().eta() << " status= "<< gpart->status();
                }// ddd
            }// if gpart
        }// if !isStopped

        else {
            gpart = evt->barcode_to_particle( line );
            if ( gpart ) {
                gpart->set_status( 2 );
                if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector::NON-transp. pz= " << gpart->momentum().pz()  
                    << " eta= "<< gpart->momentum().eta()  << " status= "<< gpart->status();
            }
        }
    }//for 

    return evt;
} 
bool CTPPSHector::SetBeamLine()
{
    edm::FileInPath b1(beam1filename.c_str());
    edm::FileInPath b2(beam2filename.c_str());
    if(m_verbosity) {
        edm::LogInfo("CTPPSHectorSetup") << "===================================================================\n"  
            << " * * * * * * * * * * * * * * * * * * * * * * * * * * * *           \n"  
            << " *                                                         *       \n"  
            << " *                   --<--<--  A fast simulator --<--<--     *     \n"  
            << " *                 | --<--<--     of particle   --<--<--     *     \n"  
            << " *  ----HECTOR----<                                          *     \n"  
            << " *                 | -->-->-- transport through-->-->--      *     \n"   
            << " *                   -->-->-- generic beamlines -->-->--     *     \n"  
            << " *                                                           *     \n"   
            << " * JINST 2:P09005 (2007)                                     *     \n"  
            << " *      X Rouby, J de Favereau, K Piotrzkowski (CP3)         *     \n"  
            << " *       http://www.fynu.ucl.ac.be/hector.html               *     \n"  
            << " *                                                           *     \n"  
            << " * Center for Cosmology, Particle Physics and Phenomenology  *     \n"  
            << " *              Universite catholique de Louvain             *     \n"  
            << " *                 Louvain-la-Neuve, Belgium                 *     \n"  
            << " *                                                         *       \n"  
            << " * * * * * * * * * * * * * * * * * * * * * * * * * * * *           \n"   
            << " CTPPSHector configuration: \n" 
            << " m_CTPPSTransport   = " << m_CTPPSTransport << "\n"
            << " lengthctpps      = " << lengthctpps << "\n"
            << " m_f_ctpps_f      =  " << m_f_ctpps_f << "\n"
            << " m_b_ctpps_b      =  " << m_b_ctpps_b << "\n"
            << "===================================================================\n";
    }  
    m_beamlineCTPPS1=NULL;
    m_beamlineCTPPS2=NULL;

    // construct beam line for CTPPS (forward 1 backward 2):                                                                                           
    if(m_CTPPSTransport && lengthctpps>0. ) {
        m_beamlineCTPPS1 = std::unique_ptr<H_BeamLine>(new H_BeamLine(-1, lengthctpps + 0.1 )); // (direction, length)
        m_beamlineCTPPS1->fill( b2.fullPath(), 1, "IP5");
        m_beamlineCTPPS2 = std::unique_ptr<H_BeamLine>(new H_BeamLine( 1, lengthctpps + 0.1 )); //
        m_beamlineCTPPS2->fill( b1.fullPath(), 1, "IP5");
        //m_beamlineCTPPS1->offsetElements( 120, 0.097 );
        //m_beamlineCTPPS2->offsetElements( 120,-0.097 );
    }
    else {
        if ( m_verbosity ) LogDebug("CTPPSHectorSetup") << "CTPPSHector: WARNING: lengthctpps=  " << lengthctpps;        
        return false;
    }
    //if (m_verbosity) {
          std::cout  << "====================================================================\n"
                     << "                  Forward beam line elements \n";
          m_beamlineCTPPS1->showElements();
          std::cout << "====================================================================\n"
                    << "                 Backward beam line elements \n";
          m_beamlineCTPPS2->showElements();
    //}

    return true;
}
