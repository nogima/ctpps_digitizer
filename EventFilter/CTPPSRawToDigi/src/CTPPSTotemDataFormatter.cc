/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 ****************************************************************************/

#include "EventFilter/CTPPSRawToDigi/interface/CTPPSTotemDataFormatter.h"

#include "EventFilter/CTPPSRawToDigi/interface/CounterChecker.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"

#include "EventFilter/CTPPSRawToDigi/interface/DiamondVFATFrame.h"
#include "DataFormats/CTPPSDigi/interface/TotemRPDigi.h"
#include "CondFormats/DataRecord/interface/TotemReadoutRcd.h"

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;
using namespace ctppstotemobjects;

//----------------------------------------------------------------------------------------------------
CTPPSTotemDataFormatter::CTPPSTotemDataFormatter(std::map<TotemFramePosition, TotemVFATInfo> const &mapping)  :  theWordCounter(0), theDigiCounter(0), mapping_(mapping)
{
  int s32 = sizeof(Word32);
  int s64 = sizeof(Word64);
  int s8  = sizeof(char);
  if ( s8 != 1 || s32 != 4*s8 || s64 != 2*s32) {
    LogError("UnexpectedSizes")
      <<" unexpected sizes: "
      <<"  size of char is: " << s8
      <<", size of Word32 is: " << s32
      <<", size of Word64 is: " << s64
      <<", send exception" ;
  }

  allDetDigis = 0;
  hasDetDigis = 0;
}
//----------------------------------------------------------------------------------------------------
void CTPPSTotemDataFormatter::formatRawData(unsigned int lvl1_ID, RawData & fedRawData, const 
    Digis & digis, std::map<std::map<const uint32_t, unsigned int>, std::map<short unsigned int, short unsigned int>> iDdet2fed)
//void CTPPSTotemDataFormatter::formatRawData(unsigned int lvl1_ID, RawData & fedRawData, const 
//    Digis & digis, std::map<std::map<const uint32_t, unsigned int>, std::map<short unsigned int, short unsigned int>> iDdet2fed,  std::map<TotemFramePosition, CTPPSTotemDigiToRaw::Record> &record)
{
  std::map<int, vector<Word32> > words;
  for (Digis::const_iterator im = digis.begin(); im != digis.end(); im++) {
    allDetDigis++;
    cms_uint32_t rawId = im->first;
    edm::LogInfo("--- RPSi") << " \t\t digi rawId = " << rawId;
    //std::cout << " \t digi rawId = " << rawId <<  " - theDigiCounter = " << theDigiCounter << std::endl; //" ---- TotemRPDetID = " << TotemRPDetId(rawId) << std::endl; 
    hasDetDigis++;
    const DetDigis & detDigis = im->second;
    for (DetDigis::const_iterator it = detDigis.begin(); it != detDigis.end(); it++) {
      theDigiCounter++;
      const TotemRPDigi & digi = (*it);
      int chipPosition = -9;  
      const int nCH = 128; 	
      int nStrip = digi.getStripNumber();  
      //std::cout << " \t\t  StripNumber  = " << nStrip << "  "  << theDigiCounter << std::endl;
      if (nStrip<= nCH ) chipPosition = 0; //std::cout << " \t\t  StripNumber  = " << nStrip << "  "  << theDigiCounter << " chipPosition = 0 - ch = " <<  nStrip<<  std::endl; 
      if (nStrip> nCH && nStrip<= 2*nCH ) chipPosition = 1; // std::cout << " \t\t  StripNumber  = " << nStrip << "  "  << theDigiCounter << " chipPosition = 1 - ch = " <<  (nStrip-nCH)  << std::endl; 
      if (nStrip> 2*nCH && nStrip<= 3*nCH ) chipPosition = 2; //std::cout << " \t\t  StripNumber  = " << nStrip << "  "  << theDigiCounter << " chipPosition = 2 - ch = " <<  (nStrip-2*nCH) << std::endl; 
      if (nStrip> 3*nCH ) chipPosition = 3; //std::cout << " \t\t  StripNumber  = " << nStrip << "  "  << theDigiCounter << " chipPosition = 3 - ch = " <<  (nStrip-3*nCH)  << std::endl; 
      //channel form DIGI		
      int channel = nStrip - chipPosition*nCH; 
	
     //std::cout << "----------------  channel = " << channel << std::endl;  		
      for (auto &p : iDdet2fed) {
	for (auto &pf : p.first){
	  cms_uint32_t pfrawId = pf.first;
	  TotemRPDetId chipId(pfrawId);
          int hwid = pf.second; 
	    //uint8_t chipPosition = chipId.chip();
            //unsigned short offset = chipPosition * 128;
	    //std::cout <<  " \t\t\t\t---  pfrawId = " << pf.first << " - chipPosition = " << unsigned(chipPosition) << " - offset = " << offset  <<  " - chipId =  " << chipId << std::endl;
	    cms_uint32_t newrawId = rawId+8192*chipPosition; 

	    if (pfrawId == newrawId){
	    //std::cout <<  " \n\t\t\t---  pfrawId = " << pf.first << " - chipPosition = " << chipPosition <<  " newRawID =  " << newrawId << " RAWiDoRIGINAL = " << rawId  << " detId" << pf.first << " hwid = "  << pf.second << endl ; 
		for (auto &ps : p.second){
	            int idxInFiber = ps.second; 		
		    //std::cout <<  " \t\t\t\t ---  channel =  " << channel <<  " hwid = " << pf.second << " newRawID =  " << newrawId << " FedId = " << ps.first << " - IdxInFiber = " << ps.second << std::endl;
	 	     Word64 word = 
  			(idxInFiber<< 0 ) 
			| (channel << 9 ) 
			| (rawId   << 10)  		
			| ( 0      << 11)  			
			| ( 0      << 12);  			
	//std::cout << " word  = " <<  print(word) << " ++++++++++++++++++++++++ " << std::endl;
		}
	   }//rawId
	} 
      }

    }//detDigis
    //std::cout << " theDigiCounter = " << theDigiCounter << std::endl; 
  }//digis
}

std::string CTPPSTotemDataFormatter::print(const  Word64 & word) const
{
  ostringstream str;
  str  <<"word64:  " << reinterpret_cast<const bitset<64>&> (word);
  return str.str();
}
//---------------------------------------------------------------------------------------------------
//iDdet2fed (detId,hwID;FedId,IdxInFiber)
/*
 * Raw data frame as sent by electronics.
 * The container is organized as follows (reversed Figure 8 at page 23 of VFAT2 manual):
 * \verbatim
 * buffer index   content       size
 * ---------------------------------------------------------------
 *   0            CRC           16 bits
 *   1->8         Channel data  128 bits, channel 0 first
 *   9            ChipID        4 constant bits (1110) + 12 bits
 *   10           EC, Flags     4 constant bits (1100) + 8, 4 bits
 *   11           BC            4 constant bits (1010) + 12 bits
 * \endverbatim
 */
