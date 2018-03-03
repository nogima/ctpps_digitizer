#include "CTPPSTotemDigiToRaw.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSDigi/interface/TotemRPDigi.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"

#include "CondFormats/CTPPSReadoutObjects/interface/TotemDAQMapping.h"
#include "CondFormats/DataRecord/interface/TotemReadoutRcd.h"

#include "EventFilter/CTPPSRawToDigi/interface/CTPPSTotemDataFormatter.h"

//raw test
#include "DataFormats/DetId/interface/DetIdCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"

using namespace std;

CTPPSTotemDigiToRaw::CTPPSTotemDigiToRaw( const edm::ParameterSet& pset ) :
  config_(pset)
{

  tTotemRPDigi = consumes<edm::DetSetVector<TotemRPDigi> >(config_.getParameter<edm::InputTag>("InputLabel")); 

  // Define EDProduct type
  produces<FEDRawDataCollection>();

  // start the counters
  eventCounter = 0;
  allDigiCounter = 0;
  allWordCounter = 0;

}

// -----------------------------------------------------------------------------
CTPPSTotemDigiToRaw::~CTPPSTotemDigiToRaw() {
  edm::LogInfo("CTPPSTotemDigiToRaw")  << " CTPPSTotemDigiToRaw destructor!";

}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void CTPPSTotemDigiToRaw::produce( edm::Event& ev,
    const edm::EventSetup& es)
{
  eventCounter++;
  edm::LogInfo("CTPPSTotemDigiToRaw") << "[CTPPSTotemDigiToRaw::produce] "
    << "event number: " << eventCounter;

  edm::Handle< edm::DetSetVector<TotemRPDigi> > digiCollection;
  label_ = config_.getParameter<edm::InputTag>("InputLabel");
  ev.getByToken( tTotemRPDigi, digiCollection);

  CTPPSTotemDataFormatter::RawData rawdata;
  CTPPSTotemDataFormatter::Digis digis;

  //typedef edm::DetSet<TotemRPDigi>::const_iterator DI;
  int digiCounter = 0; 
  typedef vector< edm::DetSet<TotemRPDigi> >::const_iterator DI;
  for (DI di=digiCollection->begin(); di != digiCollection->end(); di++) {
    digiCounter += (di->data).size();
    //std::cout << "print digi size =  "  <<  (di->data).size() <<  std::endl;
    //std::cout << "D2R - detId strip =  "  <<  di->detId()  << " TotemRPDetId = " << TotemRPDetId(di->detId()) << std::endl;
    digis[ di->detId()] = di->data;
    for (edm::DetSet<TotemRPDigi>::const_iterator dit = di->begin(); dit != di->end(); ++dit) {
      //std::cout << "\t strip number  =  " << dit->getStripNumber() << std::endl;
    }
  }
  allDigiCounter += digiCounter;
  //std::cout << "print all digi size =  "  <<  allDigiCounter << std::endl; 
  edm::ESHandle<TotemDAQMapping> mapping;
  // label of the CTPPS sub-system
  std::string subSystemName = "TrackingStrip";
  map<TotemFramePosition, Record> records;
  if (recordWatcher.check( es )) {
    es.get<TotemReadoutRcd>().get(mapping);
    //std::cout << " recordWatcher TotemReadoutRcd "  << std::endl; 
    for (const auto &p : mapping->VFATMapping) {
      TotemVFATStatus st;
      st.setMissing(true);
      records[p.first] = { &p.second, nullptr, st };
      //std::cout << " first === " << p.first << " second === " << p.second   << std::endl; 
      const uint32_t pID = (p.second.symbolicID).symbolicID;  
      unsigned int phwID = p.second.hwID; 
      std::map<const uint32_t, unsigned int> mapSymb; 
      mapSymb.insert(std::pair<const uint32_t, unsigned int>(pID,phwID)); 

      short unsigned int pFediD = p.first.getFEDId();
      short unsigned int pIdxInFiber = p.first.getIdxInFiber();	
      std::map<short unsigned int, short unsigned int> mapFedFiber; 
      mapFedFiber.insert(std::pair<short unsigned int, short unsigned int>(pFediD,pIdxInFiber)); 
      iDdet2fed_.insert(std::pair<std::map<const uint32_t, unsigned int>, std::map<short unsigned int, short unsigned int>>(mapSymb,mapFedFiber));
      //std::cout <<    "---------------- " << " pid = " << pID << " hwID =  " << phwID <<  std::endl;
    }
  }

  for (const auto &p : records)
  {
    //     if (p.second.status.isMissing())

    Record record = p.second;
    TotemRPDetId chipId(record.info->symbolicID.symbolicID);
    uint8_t chipPosition = chipId.chip();
    //std::cout <<    "---------------- " << record.info->symbolicID.symbolicID << " frame for VFAT " << p.first << " chipPosition = " << unsigned(chipPosition) <<  std::endl;
  }

  CTPPSTotemDataFormatter formatter(mapping->VFATMapping);

  //formatter.formatRawData( ev.id().event(), rawdata, digis, iDdet2fed_, records_);
  formatter.formatRawData( ev.id().event(), rawdata, digis, iDdet2fed_);

//88888888888888888888888888888888888888888888888888888888888888888888
  /*
  // pack raw data into collection
  for (auto it = fedIds_.begin(); it != fedIds_.end(); it++) { 
  FEDRawData& fedRawData = buffers->FEDData( *it );
  CTPPSTotemDataFormatter::RawData::iterator fedbuffer = rawdata.find( *it );
  if( fedbuffer != rawdata.end() ) fedRawData = fedbuffer->second;

  int nWords = fedRawData.size()/sizeof(Word64);
  if (nWords!=0) data_exist = true; 
  }
  if(data_exist) { 
  allWordCounter += formatter.nWords();

  if (debug) LogDebug("CTPPSTotemDigiToRaw") 
  << "Words/Digis this ev: "<<digiCounter<<"(fm:"<<formatter.nDigis()<<")/"
  <<formatter.nWords()
  <<"  all: "<< allDigiCounter <<"/"<<allWordCounter;

  ev.put(std::move(buffers));
  }
  */
}
