/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors: 
*  Jan Ka??par (jan.kaspar@gmail.com) 
*    
****************************************************************************/

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/CTPPSDigi/interface/TotemRPDigi.h"
//#include "DataFormats/CTPPSDigi/interface/TotemRPDigiCollection.h"
#include "DataFormats/CTPPSDigi/interface/TotemTriggerCounters.h"
#include "DataFormats/CTPPSDigi/interface/TotemVFATStatus.h"
#include "DataFormats/CTPPSDigi/interface/TotemFEDInfo.h"
#include "DataFormats/CTPPSDigi/interface/CTPPSDiamondDigi.h"

#include "DataFormats/CTPPSDigi/interface/CTPPSPixelDigi.h"
#include "DataFormats/CTPPSDigi/interface/CTPPSPixelDigiCollection.h"

#include "DataFormats/CTPPSDigi/interface/RPStripDigi.h"
#include "DataFormats/CTPPSDigi/interface/RPDetTrigger.h"
#include "DataFormats/CTPPSDigi/interface/RPDigCluster.h"

#include <vector>

namespace DataFormats_CTPPSDigi {
  struct dictionary {
    TotemRPDigi rp_str_dig;
    edm::DetSet<TotemRPDigi> ds_rp_str_dig;
    std::vector<TotemRPDigi> vec_rp_str_dig;
    edm::DetSetVector<TotemRPDigi> dsv_rp_str_dig;
    std::vector<edm::DetSet<TotemRPDigi> > vec_ds_rp_str_dig;
    edm::Wrapper<edm::DetSet<TotemRPDigi> > wds_rp_str_dig;
    edm::Wrapper<edm::DetSetVector<TotemRPDigi> > wdsv_rp_str_dig;

    RPStripDigi rp_str_digT;
    edm::DetSet<RPStripDigi> ds_rp_str_digT;
    std::vector<RPStripDigi> vec_rp_str_digT;
    edm::DetSetVector<RPStripDigi> dsv_rp_str_digT;
    std::vector<edm::DetSet<RPStripDigi> > vec_ds_rp_str_digT;
    edm::Wrapper<edm::DetSet<RPStripDigi> > wds_rp_str_digT;
    edm::Wrapper<edm::DetSetVector<RPStripDigi> > wdsv_rp_str_digT;

    RPDetTrigger rp_str_tri;
    edm::DetSet<RPDetTrigger> ds_rp_str_tri;
    std::vector<RPDetTrigger> vec_rp_str_tri;
    std::vector<edm::DetSet<RPDetTrigger> > vec_ds_rp_str_tri;
    edm::DetSetVector<RPDetTrigger> dsv_rp_str_tri;
    edm::Wrapper<edm::DetSet<RPDetTrigger> > wds_rp_str_tri;
    edm::Wrapper<edm::DetSetVector<RPDetTrigger> > wdsv_rp_str_tri;

    RPDigCluster dc;
    edm::DetSet<RPDigCluster> dsdc;
    std::vector<RPDigCluster> svdc;
    std::vector<edm::DetSet<RPDigCluster> > svdsdc;
    edm::DetSetVector<RPDigCluster> dsvdc;
    edm::Wrapper<edm::DetSetVector<RPDigCluster> > wdsvdc;

    TotemTriggerCounters dummy10;
    edm::Wrapper<TotemTriggerCounters> dummy11;

    std::map<unsigned int, uint64_t> dummy27;

    TotemVFATStatus dummy30;
    edm::Wrapper< TotemVFATStatus > dummy31;
    edm::DetSetVector<TotemVFATStatus> dummy32;
    edm::Wrapper< edm::DetSetVector<TotemVFATStatus> > dummy33;

    std::bitset<8> dummy50;
    edm::Wrapper< std::bitset<8> > dummy51;

    TotemFEDInfo fi;
    std::vector<TotemFEDInfo> v_fi;
    edm::Wrapper<std::vector<TotemFEDInfo>> w_v_fi;

    CTPPSDiamondDigi rm_diamo_dig;
    edm::DetSet<CTPPSDiamondDigi> ds_rp_diamo_dig;
    std::vector<CTPPSDiamondDigi> vec_rp_diamo_dig;
    edm::DetSetVector<CTPPSDiamondDigi> dsv_rp_diamo_dig;
    std::vector<edm::DetSet<CTPPSDiamondDigi> > vec_ds_rp_diamo_dig;
    edm::Wrapper<edm::DetSet<CTPPSDiamondDigi> > wds_rp_diamo_dig;
    edm::Wrapper<edm::DetSetVector<CTPPSDiamondDigi> > wdsv_rp_diamo_dig;

    HPTDCErrorFlags rm_hptdcerr;
    CTPPSPixelDigi ff0;
    CTPPSPixelDigiCollection ffc0;
    std::vector<CTPPSPixelDigi>  ff1;
    edm::DetSet<CTPPSPixelDigi>  ff2;
    std::vector<edm::DetSet<CTPPSPixelDigi> >  ff3;
    edm::DetSetVector<CTPPSPixelDigi> ff4;


    edm::Wrapper<CTPPSPixelDigi> wff0;
    edm::Wrapper<CTPPSPixelDigiCollection> wffc0;
    edm::Wrapper< std::vector<CTPPSPixelDigi>  > wff1;
    edm::Wrapper< edm::DetSet<CTPPSPixelDigi> > wff2;
    edm::Wrapper< std::vector<edm::DetSet<CTPPSPixelDigi> > > wff3;
    edm::Wrapper< edm::DetSetVector<CTPPSPixelDigi> > wff4;

  };
}
