// -*- C++ -*-
//
// Package:    Validation/CTPPSValidation
// Class:      CTPPSValidation
// 
/**\class CTPPSValidation CTPPSValidation.cc Validation/CTPPSValidation/plugins/CTPPSValidation.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Maria Elena Pol
//         Created:  Mon, 20 Nov 2017 09:35:19 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelLocalTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "RecoCTPPS/PixelLocal/interface/CTPPSPixelRecHitProducer.h" 
#include "DataFormats/CTPPSDetId/interface/CTPPSPixelDetId.h"
#include "DataFormats/CTPPSDigi/interface/CTPPSPixelDigi.h"
#include "EventFilter/CTPPSRawToDigi/interface/CTPPSPixelDataFormatter.h"

#include <TROOT.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TMath.h>
#include <TStyle.h>
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include <cmath>
#include <vector>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class CTPPSValidation : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit CTPPSValidation(const edm::ParameterSet&);
		~CTPPSValidation();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

		typedef std::vector<CTPPSPixelDigi> DetDigis;
		typedef std::map<cms_uint32_t,DetDigis> Digis;

	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
		edm::EDGetTokenT< edm::HepMCProduct > mcEventToken; // label of MC event
		edm::Handle< edm::HepMCProduct > EvtHandle ;

		edm::EDGetTokenT< edm::DetSetVector<CTPPSPixelLocalTrack> > pixelTrackToken_;

		edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelRecHit>> tokenCTPPSPixelRecHit_;

		edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelDigi>> tokenCTPPSPixelDigi_; 
		edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelDigi>> tokenCTPPSRPix_; 

		edm::EDGetTokenT<edm::PSimHitContainer> psimHitToken;
		edm::Handle<edm::PSimHitContainer> PSimHit;

		const double   ProtonMass = 0.93827;
		const double ProtonMassSQ = pow(ProtonMass,2);

		int EventKind, evNumber= 0;
		//int sh45size = 0 , lt45size = 0; 
		//int sh56size = 0 , lt56size = 0; 
		//int  lt45size1 = 0; 
		//int  lt56size1 = 0; 
		//TFileService
		unsigned int tracks45=0;
		unsigned int tracks56=0;
		unsigned int simhits45=0;  
		unsigned int simhits56=0;  
		unsigned int digiBefRaw45=0;  
		unsigned int digiBefRaw56=0;  
		unsigned int digiAftRaw45=0;  
		unsigned int digiAftRaw56=0;  
		unsigned int rechits45=0;  
		unsigned int rechits56=0;  
		edm::Service<TFileService> fs;
		TTree* AnalysisTree;
		// Proton Histograms
		TH2F *h_LHCT_xVsy_ARMFs1;
		TH2F *h_LHCT_xVsy_ARMBs1;
		TH1F *h_LHCT_x_ARMFs1;
		TH1F *h_LHCT_y_ARMFs1;
		TH1F *h_LHCT_x_ARMBs1;
		TH1F *h_LHCT_y_ARMBs1;
		TH1F *h_x_ARMF, *h_x_ARMB, *h_y_ARMF, *h_y_ARMB;
		TH2F *h_xx_ARMF, *h_xx_ARMB, *h_yy_ARMF, *h_yy_ARMB;
		TH1F *h_Dx_ARMF, *h_Dx_ARMB, *h_Dy_ARMF, *h_Dy_ARMB;
		TH1F *hdiff_x_recsim_ARMF, *hdiff_y_recsim_ARMF;
		TH1F *hdiff_x_recsim_ARMB, *hdiff_y_recsim_ARMB;
		TH1F *h_LHCT_eta_ARMF, *h_LHCT_eta_ARMB, *h_LHCT_phimom_ARMF, *h_LHCT_phimom_ARMB;
		TH1F *h_LHCT_px_ARMF, *h_LHCT_px_ARMB;
		TH1F *h_LHCT_py_ARMF, *h_LHCT_py_ARMB;
		TH1F *h_LHCT_pz_ARMF, *h_LHCT_pz_ARMB;
		TH1F *h_LHCT_pt_ARMF, *h_LHCT_pt_ARMB, *h_LHCT_p_ARMF, *h_LHCT_p_ARMB;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CTPPSValidation::CTPPSValidation(const edm::ParameterSet& iConfig)

{
	//now do what ever initialization is needed
	usesResource("TFileService");
	mcEventToken    = mayConsume<edm::HepMCProduct>(iConfig.getUntrackedParameter<edm::InputTag>("MCEvent",std::string("")));

	auto tagPixelTrack_ = iConfig.getParameter<edm::InputTag>("tagPixelTrack"); 
	auto tagPixelRecHit_ = iConfig.getParameter<edm::InputTag>("tagPixelRecHit");
	auto tagPixelDigiBefRaw_ = iConfig.getParameter<edm::InputTag>("tagPixelDigiBefRaw");
	auto tagPixelDigiAftRaw_ = iConfig.getParameter<edm::InputTag>("tagPixelDigiAftRaw");
	auto tagPixelSimHit_ = iConfig.getUntrackedParameter<edm::InputTag>("psimHitTag"); 
	psimHitToken = mayConsume<edm::PSimHitContainer>(tagPixelSimHit_);
	tokenCTPPSPixelDigi_ = consumes<edm::DetSetVector<CTPPSPixelDigi> >(tagPixelDigiBefRaw_); 	
	tokenCTPPSRPix_ = consumes<edm::DetSetVector<CTPPSPixelDigi> >(tagPixelDigiAftRaw_); 	
	tokenCTPPSPixelRecHit_ = consumes<edm::DetSetVector<CTPPSPixelRecHit> >(edm::InputTag(tagPixelRecHit_));
	if (not tagPixelTrack_.label().empty()) pixelTrackToken_   = consumes< edm::DetSetVector<CTPPSPixelLocalTrack> >  (tagPixelTrack_);
}


CTPPSValidation::~CTPPSValidation()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
CTPPSValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	iEvent.getByToken( mcEventToken, EvtHandle );
	const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
	EventKind = Evt->signal_process_id();

	edm::Handle< edm::DetSetVector<CTPPSPixelLocalTrack> > inputPixelTracks;
	iEvent.getByToken( psimHitToken, PSimHit);

	double x45 = 999., y45 = 999.;
	double x56 = 999., y56 = 999.;
	//double z45 = 0.  , z56 = 0. ; 
/*	for(HepMC::GenEvent::particle_const_iterator i=Evt->particles_begin(); i != Evt->particles_end();++i) {
		int myId = (*i)->pdg_id();
		if(myId!=2212) continue;

		HepMC::GenVertex* pv = (*i)->production_vertex();
		HepMC::FourVector vertex = pv->position();
		//status() == 1 -> vertex at detector entrance
		if(vertex.eta() >= 8. && (*i)->status() == 1) { //Arm F
			h_LHCT_xVsy_ARMFs1->Fill(vertex.x(),vertex.y());
			h_LHCT_x_ARMFs1->Fill(vertex.x());
			h_LHCT_y_ARMFs1->Fill(vertex.y());
			x45 = vertex.x(); 
			y45 = vertex.y(); 
			//z45 = 219550.;

		}
		if(vertex.eta() <= -8. && (*i)->status() == 1) { //Arm B
			h_LHCT_xVsy_ARMBs1->Fill(vertex.x(),vertex.y());
			h_LHCT_x_ARMBs1->Fill(vertex.x());
			h_LHCT_y_ARMBs1->Fill(vertex.y());
			x56 = vertex.x(); 
			y56 = vertex.y(); 
			//z56 = -219550.;
		}
	}
*/	if(PSimHit.isValid()) {
		edm::PSimHitContainer::const_iterator itPSimHit;
		unsigned int simsize56 =0, simsize45 =0;  
		for (itPSimHit = PSimHit->begin(); itPSimHit != PSimHit->end(); ++itPSimHit) {
			const uint32_t simHitId = itPSimHit->detUnitId(); 
			if(simHitId & 0x1000000) simsize56++;
			else simsize45++;
		}
		if(simsize56 != 0 ) simhits56++;
		if(simsize45 != 0 ) simhits45++;
	} 

	edm::Handle<edm::DetSetVector<CTPPSPixelRecHit> > recHits;
	iEvent.getByToken(tokenCTPPSPixelRecHit_, recHits);
	edm::DetSetVector<CTPPSPixelRecHit>::const_iterator rhDSViter = recHits->begin();
	//std::cout<< " (*rhDSViter).size() " << (*rhDSViter).size()  << std::endl;
	//        if((*rhDSViter).size()!=0){
	unsigned int recsize56 =0, recsize45 =0;  
	for (; rhDSViter != recHits->end(); rhDSViter++) {
		edm::PSimHitContainer::const_iterator itPSimHit;
		if(PSimHit.isValid()) {
			for (itPSimHit = PSimHit->begin(); itPSimHit != PSimHit->end(); ++itPSimHit) {
				const uint32_t simHitId = itPSimHit->detUnitId(); 
				// std::cout << " simHitId =  " << simHitId << " SimHit x,y " << itPSimHit->entryPoint().x()<< ", " << itPSimHit->entryPoint().y() <<  std::endl; 
				edm::DetSet<CTPPSPixelRecHit>::const_iterator begin = (*rhDSViter).begin();
				edm::DetSet<CTPPSPixelRecHit>::const_iterator end = (*rhDSViter).end();
				for (edm::DetSet<CTPPSPixelRecHit>::const_iterator rh = begin; rh != end; rh++) {
					const uint32_t recHitId = rhDSViter->detId(); 
					//std::cout << " recHitId =  " <<  rhDSViter->detId() << " RecHit x,y " << (*rh).getPoint().x()<< ", "<<  (*rh).getPoint().y() << std::endl;
					if (simHitId == recHitId ) {
						if(recHitId < 2030e+6) { 
							hdiff_x_recsim_ARMF->Fill((*rh).getPoint().x()-itPSimHit->entryPoint().x());
							hdiff_y_recsim_ARMF->Fill((*rh).getPoint().y()-itPSimHit->entryPoint().y());
							recsize45++;
						} 
						else {
							hdiff_x_recsim_ARMB->Fill((*rh).getPoint().x()-itPSimHit->entryPoint().x());
							hdiff_y_recsim_ARMB->Fill((*rh).getPoint().y()-itPSimHit->entryPoint().y());
							recsize56++;
						}
					} 
				}
			}
		}
	}
	if(recsize56 != 0 ) rechits56++;
	if(recsize45 != 0 ) rechits45++;

	if (not pixelTrackToken_.isUninitialized()){
		iEvent.getByToken( pixelTrackToken_, inputPixelTracks );
		// process tracks from pixels  
		for ( const auto& rpv : *inputPixelTracks ) {
			const uint32_t rpId = rpv.detId();
			// std::cout << rpv.size() <<  " ------------- rpId = "  << rpId << std::endl; 
			if(rpv.size()==1) for ( const auto& trk : rpv ) {
				//std::cout << " X Y Z "  << trk.getX0() <<  " " << trk.getY0()  <<  " " << trk.getZ0()  << std::endl; 
				//std::cout << " X Y Z 45 "  << x45 <<  " " << y45  <<  " " << std::endl; 
				//std::cout << " X Y Z 56 "  << x56 <<  " " << y56  <<  " " << std::endl; 
				if ( !trk.isValid() ) continue;
				for(HepMC::GenEvent::particle_const_iterator i=Evt->particles_begin(); i != Evt->particles_end();++i) {
					int myId = (*i)->pdg_id();
					if(myId!=2212) continue;

					HepMC::GenVertex* pv = (*i)->production_vertex();
					HepMC::FourVector vertex = pv->position();
					HepMC::FourVector momentum=(*i)->momentum();
					const HepMC::FourVector p((*i)->momentum());
					double px = momentum.x();
					double py = momentum.y();
					double pz = momentum.z();
					double e = sqrt(px*px+py*py+pz*pz+ProtonMassSQ);
					TLorentzVector* proton = new TLorentzVector(px,py,pz,e);

					//status() == 1 -> vertex at detector entrance
					if(rpId & 0x1000000){
						if(vertex.eta() <= -8. && (*i)->status() == 1) { //Arm B
							h_LHCT_xVsy_ARMBs1->Fill(vertex.x(),vertex.y());
							h_LHCT_x_ARMBs1->Fill(vertex.x());
							h_LHCT_y_ARMBs1->Fill(vertex.y());
							x56 = vertex.x(); 
							y56 = vertex.y(); 
							h_x_ARMB->Fill(trk.getX0());
							h_y_ARMB->Fill(trk.getY0());
							h_xx_ARMB->Fill(trk.getX0(),x56);
							h_yy_ARMB->Fill(trk.getY0(),y56);
							h_Dx_ARMB->Fill(trk.getX0()-x56);
							h_Dy_ARMB->Fill(trk.getY0()-y56);
							tracks56++;
							//z45 = 219550.;
							h_LHCT_pt_ARMB->Fill(proton->Pt());
							h_LHCT_p_ARMB->Fill(proton->P());
							h_LHCT_px_ARMB->Fill(proton->Px());
							h_LHCT_py_ARMB->Fill(proton->Py());
							h_LHCT_pz_ARMB->Fill(proton->Pz());
							h_LHCT_eta_ARMB->Fill(proton->Eta());
							h_LHCT_phimom_ARMB->Fill(proton->Phi());
						}
					}
					else { 
						if(vertex.eta() >= 8. && (*i)->status() == 1) { //Arm F
							h_LHCT_xVsy_ARMFs1->Fill(vertex.x(),vertex.y());
							h_LHCT_x_ARMFs1->Fill(vertex.x());
							h_LHCT_y_ARMFs1->Fill(vertex.y());
							x45 = vertex.x(); 
							y45 = vertex.y(); 
							h_x_ARMF->Fill(trk.getX0());
							h_y_ARMF->Fill(trk.getY0());
							h_xx_ARMF->Fill(trk.getX0(),x45);
							h_yy_ARMF->Fill(trk.getY0(),y45);
							h_Dx_ARMF->Fill(trk.getX0()-x45);
							h_Dy_ARMF->Fill(trk.getY0()-y45);
							h_LHCT_pt_ARMF->Fill(proton->Pt());
							h_LHCT_p_ARMF->Fill(proton->P());
							h_LHCT_px_ARMF->Fill(proton->Px());
							h_LHCT_py_ARMF->Fill(proton->Py());
							h_LHCT_pz_ARMF->Fill(proton->Pz());
							h_LHCT_eta_ARMF->Fill(proton->Eta());
							h_LHCT_phimom_ARMF->Fill(proton->Phi());
							//std::cout << "---  "  << (trk.getX0()-x45) << " " << trk.getX0() << " "  << x45 << std::endl;  
							//std::cout << "---  "  << (trk.getY0()-y45) << " " << trk.getY0() << " "  << y45 << std::endl;  
							tracks45++;
							//z56 = -219550.;
						}
					}
				}				
			}
			//std::cout << " localTracks45 "  << tracks45 << std::endl; 
			//std::cout << " localTracks56 "  << tracks56 << std::endl;
		} 
	}

	edm::Handle< edm::DetSetVector<CTPPSPixelDigi> > digiCollection;
	iEvent.getByToken( tokenCTPPSPixelDigi_, digiCollection);
	typedef std::vector< edm::DetSet<CTPPSPixelDigi> >::const_iterator DI;
	CTPPSPixelDataFormatter::Digis digis;
	int digiCounter = 0; 
	for (DI di=digiCollection->begin(); di != digiCollection->end(); di++) {
		digiCounter += (di->data).size(); 
		digis[ di->id] = di->data;
	}
	unsigned int pixsize45 =0, pixsize56 =0; 
	for (Digis::const_iterator im = digis.begin(); im != digis.end(); im++) {
		cms_uint32_t pixelId = im->first;
		//std::cout << " DigiBefRaw rawId = " << pixelId << std::endl;
		if(pixelId < 2030e+6) pixsize45++; 
		else  pixsize56++; 
	}
	if(pixsize56 != 0 ) digiBefRaw56++;
	if(pixsize45 != 0 ) digiBefRaw45++;

	//RpixDetDigitizer
	edm::Handle< edm::DetSetVector<CTPPSPixelDigi> > digiCollection2;
	iEvent.getByToken( tokenCTPPSRPix_, digiCollection2);
	typedef std::vector< edm::DetSet<CTPPSPixelDigi> >::const_iterator DI2;
	CTPPSPixelDataFormatter::Digis digis2;
	for (DI2 di=digiCollection2->begin(); di != digiCollection2->end(); di++) digis2[ di->id] = di->data;
	unsigned int pixAftsize45 =0, pixAftsize56 =0; 
	for (Digis::const_iterator im = digis2.begin(); im != digis2.end(); im++) {
		cms_uint32_t pixelId = im->first;
		if(pixelId < 2030e+6) pixAftsize45++; 
		else  pixAftsize56++; 
	}
	if(pixAftsize56 != 0 ) digiAftRaw56++;
	if(pixAftsize45 != 0 ) digiAftRaw45++;


	//if(simhits45 > 0) sh45size++;
	//if(simhits56 > 0) sh56size++;
	//if(tracks45 > 0) lt45size++;
	//if(tracks56 > 0) lt56size++;

	//if(tracks45  == 1 ) lt45size1++;
	//if(tracks56  == 1 ) lt56size1++;
	//std::cout << " simhits56 e tracks56 =  "  << simhits56 << " " << tracks56 << std::endl; 
	//std::cout << " simhits45 e tracks45 =  "  << simhits45 << " " << tracks45 << std::endl; 

}

// ------------ method called once each job just before starting event loop  ------------
	void 
CTPPSValidation::beginJob()
{
	h_LHCT_xVsy_ARMFs1 = fs->make<TH2F>( "LHCT_xVsyARMFs1" , "LHCT_xVsyARMFs1; x [mm] ; y [mm]" , 100,-30,0.,100,-15,15);
	h_LHCT_xVsy_ARMBs1 = fs->make<TH2F>( "LHCT_xVsyARMBs1" , "LHCT_xVsyARMBs1; x [mm] ; y [mm]" , 100,-30,0.,100,-15,15);
	h_LHCT_x_ARMFs1 = fs->make<TH1F>( "LHCT_xARMFs1" , "LHCT_xARMFs1; x [mm]" , 100,-30,0.);
	h_LHCT_y_ARMFs1 = fs->make<TH1F>( "LHCT_yARMFs1" , "LHCT_yARMFs1; y [mm]" , 100,-15,15.);
	h_LHCT_x_ARMBs1 = fs->make<TH1F>( "LHCT_xARMBs1" , "LHCT_xARMBs1; x [mm]" , 100,-30,0.);
	h_LHCT_y_ARMBs1 = fs->make<TH1F>( "LHCT_yARMBs1" , "LHCT_yARMBs1; y [mm]" , 100,-15,15.);
	h_x_ARMF = fs->make<TH1F>( "xARMF" , "global coordinates ARMF ; x_track [mm]" , 100,-30,0.);
	h_x_ARMB = fs->make<TH1F>( "xARMB" , "global coordinates ARMB ; x_track [mm]" , 100,-30,0.);
	h_y_ARMF = fs->make<TH1F>( "yARMF" , "global coordinates ARMF ; y_track [mm]" , 100,-15,15);
	h_y_ARMB = fs->make<TH1F>( "yARMB" , "global coordinates ARMB ; y_track [mm]" , 100,-15,15);
	h_xx_ARMF = fs->make<TH2F>( "xxARMF" , "global coordinates ARMF ; x_track [mm] ; x_lhct [mm]" , 100,-20,0.,100,-20,0);
	h_xx_ARMB = fs->make<TH2F>( "xxARMB" , "global coordinates ARMB ; x_track [mm] ; x_lhct [mm]" , 100,-20,0.,100,-20,0);
	h_yy_ARMF = fs->make<TH2F>( "yyARMF" , "global coordinates ARMF ; y_track [mm] ; y_lhct [mm]" , 100,-10,10.,100,-10,10);
	h_yy_ARMB = fs->make<TH2F>( "yyARMB" , "global coordinates ARMB ; y_track [mm] ; y_lhct [mm]" , 100,-10,10.,100,-10,10);
	h_Dx_ARMF = fs->make<TH1F>( "DxARMF" , "global coordinates ARMF ; x_track - x_lhct [mm]" , 100,-.3,.3);
	h_Dx_ARMB = fs->make<TH1F>( "DxARMB" , "global coordinates ARMB ; x_track - x_lhct [mm]" , 100,-.3,.3);
	h_Dy_ARMF = fs->make<TH1F>( "DyARMF" , "global coordinates ARMF ; y_track - y_lhct [mm]" , 100,-.3,.3);
	h_Dy_ARMB = fs->make<TH1F>( "DyARMB" , "global coordinates ARMB ; y_track - y_lhct [mm]" , 100,-.3,.3);
	hdiff_x_recsim_ARMF = fs->make<TH1F>( "hdiff_x_recsimARMF" , "local coordinates ARMF ; x_rec - x_sim [mm]" , 100,-1.,1.);
	hdiff_y_recsim_ARMF = fs->make<TH1F>( "hdiff_y_recsimARMF" , "local coordinates ARMF ; y_rec - y_sim [mm]" , 100,-1.,1.);
	hdiff_x_recsim_ARMB = fs->make<TH1F>( "hdiff_x_recsimARMB" , "local coordinates ARMB ; x_rec - x_sim [mm]" , 100,-1.,1.);
	hdiff_y_recsim_ARMB = fs->make<TH1F>( "hdiff_y_recsimARMB" , "local coordinates ARMB ; y_rec - y_sim [mm]" , 100,-1.,1.);
	h_LHCT_pt_ARMF = fs->make<TH1F>( "LHCT_pt_armF" , "LHCT_ptARMF; pt [GeV] ; nEvents", 100, .0, 2.0 );
	h_LHCT_px_ARMF = fs->make<TH1F>( "LHCT_px_armF" , "LHCT_pxARMF; px [GeV] ; nEvents", 100, -3.0, 3.0 );
	h_LHCT_py_ARMF = fs->make<TH1F>( "LHCT_py_armF" , "LHCT_pyARMF; py [GeV] ; nEvents", 100, -3.0, 3.0 );
	h_LHCT_pz_ARMF = fs->make<TH1F>( "LHCT_pz_armF" , "LHCT_pzARMF; pz [GeV] ; nEvents", 300, 4500.0, 6500.0 );
	h_LHCT_p_ARMF = fs->make<TH1F>( "LHCT_p_armF" , "LHCT_pARMF; p [GeV] ; nEvents", 300, 4500.0, 6500.0 );
	h_LHCT_eta_ARMF = fs->make<TH1F>( "LHCT_etaARMF" , "LHCT_etaARMF; #eta; nEvents" , 100, 8., 15.);
	h_LHCT_phimom_ARMF = fs->make<TH1F>( "LHCT_phimomARMF" , "LHCT_phimomARMF; #phi_P; nEvents" , 100, -3.2, 3.2 );
	h_LHCT_pt_ARMB = fs->make<TH1F>( "LHCT_pt_armB" , "LHCT_ptARMB; pt [GeV] ; nEvents", 100, .0, 2.0 );
	h_LHCT_px_ARMB = fs->make<TH1F>( "LHCT_px_armB" , "LHCT_pxARMB; px [GeV] ; nEvents", 100, -3.0, 3.0 );
	h_LHCT_py_ARMB = fs->make<TH1F>( "LHCT_py_armB" , "LHCT_pyARMB; py [GeV] ; nEvents", 100, -3.0, 3.0 );
	h_LHCT_p_ARMB = fs->make<TH1F>( "LHCT_p_armB" , "LHCT_pARMB; p [GeV] ; nEvents", 300, 4500.0, 6500.0 );
	h_LHCT_pz_ARMB = fs->make<TH1F>( "LHCT_pz_armB" , "LHCT_pzARMB; p [GeV] ; nEvents", 300, -6500.0, -4500.0 );
	h_LHCT_eta_ARMB  = fs->make<TH1F>( "LHCT_etaARMB" , "LHCT_etaARMB; #eta; nEvents" , 100, -15.,-8.);
	h_LHCT_phimom_ARMB = fs->make<TH1F>( "LHCT_phimomARMB" , "LHCT_phimomARMB; #phi_P; nEvents" , 100, -3.2, 3.2 );
	gStyle->SetPalette(1);
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
CTPPSValidation::endJob() 
{
	std::cout << " -----------------------------------------------------------" << std::endl;
	std::cout << " track45 , track56 =  "  << tracks45 << " " << tracks56 << std::endl; 
	std::cout << " rechits45 , rechits56 =  "  << rechits45 << " " << rechits56 << std::endl; 
	std::cout << " digiAftRaw45 , digiAftRaw56 =  "  << digiAftRaw45 << " " << digiAftRaw56 << std::endl; 
	std::cout << " digiBefRaw45 , digiBefRaw56 =  "  << digiBefRaw45 << " " << digiBefRaw56 << std::endl; 
	std::cout << " simhits45 , simhits56 =  "  << simhits45 << " " << simhits56 << std::endl; 
	//std::cout << " sh45size , sh56size =  "  << sh45size << " " << sh56size << std::endl; 
	//std::cout << " lt45size , lt56size =  "  << lt45size << " " << lt56size << std::endl; 
	//std::cout << " lt45size1 , lt56size1 =  "  << lt45size1 << " " << lt56size1 << std::endl; 
	std::cout << " -----------------------------------------------------------" << std::endl;
	TCanvas histo ("histo");
	h_LHCT_xVsy_ARMFs1->Draw();
	histo.Print("histo.pdf(");
	histo.Print("h_LHCT_xVsy_ARMFs1.png");
	h_LHCT_xVsy_ARMBs1->Draw();
	histo.Print("histo.pdf");   
	histo.Print("h_LHCT_xVsy_ARMBs1.png");

	h_LHCT_x_ARMFs1->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_LHCT_x_ARMFs1.png");
	h_LHCT_y_ARMFs1->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_LHCT_y_ARMFs1.png");
	h_LHCT_x_ARMBs1->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_LHCT_x_ARMBs1.png");
	h_LHCT_y_ARMBs1->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_LHCT_y_ARMBs1.png");
	h_x_ARMF->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_x_ARMF.png");
	h_x_ARMB->Draw();
	histo.Print("h_x_ARMB.png");
	histo.Print("histo.pdf");
	h_y_ARMF->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_y_ARMF.png");
	h_y_ARMB->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_y_ARMB.png");
	h_xx_ARMF->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_xx_ARMF.png");
	h_xx_ARMB->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_xx_ARMB.png");
	h_yy_ARMF->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_yy_ARMF.png");
	h_yy_ARMB->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_yy_ARMB.png");
	h_Dx_ARMF->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_Dx_ARMF.png");
	h_Dx_ARMB->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_Dx_ARMB.png");
	h_Dy_ARMF->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_Dy_ARMF.png");
	h_Dy_ARMB->Draw();
	histo.Print("histo.pdf");
	histo.Print("h_Dy_ARMB.png");
	hdiff_x_recsim_ARMF->Draw();
	histo.Print("histo.pdf"); 
	histo.Print("hdiff_x_recsim_ARMF.png");
	hdiff_y_recsim_ARMF->Draw();
	histo.Print("histo.pdf");
	histo.Print("hdiff_y_recsim_ARMF.png");
	hdiff_x_recsim_ARMB->Draw();
	histo.Print("histo.pdf"); 
	histo.Print("hdiff_x_recsim_ARMB.png");
	hdiff_y_recsim_ARMB->Draw();
	histo.Print("histo.pdf)");
	histo.Print("hdiff_y_recsim_ARMB.png");

	h_LHCT_pt_ARMF->Draw(); 
	histo.Print("h_LHCT_pt_ARMF.png"); 
	h_LHCT_px_ARMF->Draw();
	histo.Print("h_LHCT_px_ARMF.png"); 
	h_LHCT_py_ARMF->Draw();
	histo.Print("h_LHCT_py_ARMF.png"); 
	h_LHCT_pz_ARMF->Draw();
	histo.Print("h_LHCT_pz_ARMF.png"); 
	h_LHCT_p_ARMF->Draw();
	histo.Print("h_LHCT_p_ARMF.png"); 
	h_LHCT_eta_ARMF->Draw();
	histo.Print("h_LHCT_eta_ARMF.png"); 
	h_LHCT_phimom_ARMF->Draw();
	histo.Print("h_LHCT_phimom_ARMF.png"); 
	h_LHCT_pt_ARMB->Draw();
	histo.Print("h_LHCT_pt_ARMB.png"); 
	h_LHCT_px_ARMB->Draw();
	histo.Print("h_LHCT_px_ARMB.png"); 
	h_LHCT_py_ARMB->Draw();
	histo.Print("h_LHCT_py_ARMB.png"); 
	h_LHCT_p_ARMB->Draw();
	histo.Print("h_LHCT_p_ARMB.png"); 
	h_LHCT_pz_ARMB->Draw();
	histo.Print("h_LHCT_pz_ARMB.png"); 
	h_LHCT_eta_ARMB->Draw();
	histo.Print("h_LHCT_eta_ARMB.png"); 
	h_LHCT_phimom_ARMB->Draw();
	histo.Print("h_LHCT_phimom_ARMB.png");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CTPPSValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CTPPSValidation);
