#include "SimCTPPS/CTPPSPixelDigiProducer/interface/RPixChargeShare.h"
#include <iostream>
#include <fstream>

RPixChargeShare::RPixChargeShare(const edm::ParameterSet &params, uint32_t det_id)
  : det_id_(det_id), theRPixDetTopology_(), sqrt_2(sqrt(2.0))
{
  verbosity_ = params.getParameter<int>("RPixVerbosity");
  signalCoupling_.clear();
  ChargeMapFile2E_[0] = params.getParameter<std::string>("ChargeMapFile2E");
  ChargeMapFile2E_[1] = params.getParameter<std::string>("ChargeMapFile2E_2X");
  ChargeMapFile2E_[2] = params.getParameter<std::string>("ChargeMapFile2E_2Y");
  ChargeMapFile2E_[3] = params.getParameter<std::string>("ChargeMapFile2E_2X2Y");
  
  double coupling_constant_ = params.getParameter<double>("RPixCoupling");

  signalCoupling_.push_back(coupling_constant_);
  signalCoupling_.push_back( (1.0-coupling_constant_)/2 );
  
  no_of_pixels_ = theRPixDetTopology_.detPixelNo();

  double xMap, yMap;
  double chargeprobcollect;
  int xUpper[]={75,150,75,150};
  int yUpper[]={50,50,100,100};
  int ix,iy;
  for(int i=0;i<4;i++){
   std::ifstream fChargeMap (ChargeMapFile2E_[i]);
   if (fChargeMap.is_open())
   {
     while (fChargeMap>>xMap>>yMap>>chargeprobcollect )
     {
      ix = int((xMap+xUpper[i])/5);
      iy = int((yMap+yUpper[i])/5);
      chargeMap2E[i][ix][iy]=chargeprobcollect;     
     }
     fChargeMap.close();
   }
   else std::cout << "Unable to open charge Map file"; 
  }
}

std::map<unsigned short, double, std::less<unsigned short> >  RPixChargeShare::Share(
     const std::vector<RPixSignalPoint> &charge_map)
{
  thePixelChargeMap.clear();
  if(verbosity_>1)
    std::cout<<"RPixChargeShare "<<det_id_<<" : Clouds to be induced= "<<charge_map.size()<<std::endl;

  double CH =0;

  for(std::vector<RPixSignalPoint>::const_iterator i=charge_map.begin(); 
      i!=charge_map.end(); ++i)
    {
      double hit_pos_x,hit_pos_y;
// Used to avoid the abort due to hits out of detector 
      if (((*i).X()+16.6/2)<0||((*i).X()+16.6/2)>16.6) {
	std::cout << "**** Attention ((*i).X()+simX_width_/2.)<0||((*i).X()+simX_width_/2.)>simX_width  " << std::endl;
	std::cout << "(*i).X() = " << (*i).X() << std::endl;
	continue;
      }
      if (((*i).Y()+24.4/2.)<0||((*i).Y()+24.4/2.)>24.4) {
	std::cout << "**** Attention ((*i).Y()+simY_width_/2.)<0||((*i).Y()+simY_width_/2.)>simY_width  " << std::endl;
	std::cout << "(*i).Y() = " << (*i).Y() << std::endl;
	continue;
      }

      std::vector<pixel_info> relevant_pixels = theRPixDetTopology_.getPixelsInvolved(
										      (*i).X(), (*i).Y(), (*i).Sigma(), hit_pos_x, hit_pos_y);
      if(verbosity_>1)
	{
	  std::cout<<"RPixChargeShare "<<det_id_<<" : relevant_pixels = "
		   <<relevant_pixels.size()<<std::endl;
	}
      for(std::vector<pixel_info>::const_iterator j = relevant_pixels.begin(); 
	  j!=relevant_pixels.end(); ++j)
	{

	  double effic = (*j).effFactor();
 
	  unsigned short pixel_no = (*j).pixelIndex();
      
	  double charge_in_pixel =  (*i).Charge()*effic;
      
	  CH += charge_in_pixel;

	  if(verbosity_>1)
	    std::cout<<"Efficiency in detector "<<det_id_<< " and pixel no " << pixel_no << "  : " <<effic<< "  ch: "<< charge_in_pixel << "   CHtot: "<< CH << std::endl;

//        QUI SI POTREBBE INTRODURRE IL CHARGE SHARING TRA I PIXELS ..................................       

	  if (signalCoupling_[0]==0.){
	    thePixelChargeMap[pixel_no] += charge_in_pixel;
	  } else {
	    int pixel_row = (*j).pixelRowNo();
	    int pixel_col = (*j).pixelColNo();
	    double pixel_lower_x=0;
	    double pixel_lower_y=0;
	    double pixel_upper_x=0;
	    double pixel_upper_y=0;
	    int    psize=0;
	    theRPixDetTopology_.pixelRange(pixel_row,pixel_col,pixel_lower_x,pixel_upper_x,pixel_lower_y,pixel_upper_y);
	    double pixel_width_x = pixel_upper_x-pixel_lower_x;
	    double pixel_width_y = pixel_upper_y-pixel_lower_y;
            if(pixel_row==0 || pixel_row==pxlRowSize-1) pixel_width_x = 0.1; // Correct edge pixel width
            if(pixel_col==0 || pixel_col==pxlColSize-1) pixel_width_y = 0.15; //
	    double pixel_center_x = pixel_lower_x + (pixel_width_x)/2.;
	    double pixel_center_y = pixel_lower_y + (pixel_width_y)/2.;
// xbin and ybin are coordinates (um) ínside the pixel as in the test beam, swapped wrt plane coordinates.
	    int xbin=int(round((((*i).Y()-pixel_center_y)+pixel_width_y/2)*1.e3/5));
            int ybin=int(round((((*i).X()-pixel_center_x)+pixel_width_x/2)*1.e3/5));
	    if(pixel_width_x<0.11&&pixel_width_y<0.151) {  // pixel 100x150 um^2
             psize = 0; 				    
	     if(xbin>xBinMax[psize]||ybin>yBinMax[psize]) continue; 
	    }
	    if(pixel_width_x>0.11&&pixel_width_y<0.151) {  // pixel 200x150 um^2
             psize = 2;
             if(xbin>xBinMax[psize]||ybin>yBinMax[psize]) continue; 
            } 
	    if(pixel_width_x<0.11&&pixel_width_y>0.151) {  // pixel 100x300 um^2
             psize = 1; 
             if(xbin>xBinMax[psize]||ybin>yBinMax[psize]) continue; 
            }
	    if(pixel_width_x>0.11&&pixel_width_y>0.151) { // pixel 200x300 um^2
             psize = 3;
             if(xbin>xBinMax[psize]||ybin>yBinMax[psize]) continue; 
            }

	    double hit2neighbour[8];
            double collect_prob = chargeMap2E[psize][xbin][ybin];
	    thePixelChargeMap[pixel_no] += charge_in_pixel*collect_prob;
	    unsigned short neighbour_no[8];
	    unsigned short m=0;
	    double closer_neighbour=0;
	    unsigned short closer_no=0;
//      Considering the 8 neighbours to share charge
	    for (int k = pixel_row - 1; k <= pixel_row + 1; k++){
	      for (int l = pixel_col - 1; l <= pixel_col + 1; l++){
		if ((k<0) || k>pxlRowSize-1 || l<0 || l>pxlColSize-1) continue;
		if ((k==pixel_row) && (l==pixel_col)) continue;
		double neighbour_pixel_lower_x=0;
		double neighbour_pixel_lower_y=0;
		double neighbour_pixel_upper_x=0;
		double neighbour_pixel_upper_y=0;
		double neighbour_pixel_center_x=0;
		double neighbour_pixel_center_y=0;
//         Check the hit approach to the neighbours
		theRPixDetTopology_.pixelRange(k,l,neighbour_pixel_lower_x,neighbour_pixel_upper_x,neighbour_pixel_lower_y,neighbour_pixel_upper_y);
		neighbour_pixel_center_x = neighbour_pixel_lower_x + 
		  (neighbour_pixel_upper_x - neighbour_pixel_lower_x)/2.;
		neighbour_pixel_center_y = neighbour_pixel_lower_y + 
		  (neighbour_pixel_upper_y - neighbour_pixel_lower_y)/2.;
		hit2neighbour[m] = 1/sqrt(pow((*i).X()-neighbour_pixel_center_x,2.)+pow((*i).Y()-neighbour_pixel_center_y,2.));
		neighbour_no[m]=l*pxlRowSize+k;
		if (hit2neighbour[m]>closer_neighbour){
		  closer_neighbour=hit2neighbour[m];
		  closer_no = neighbour_no[m];
		}
		m++;
	      }
	    }
	    double chargetransfereff = (1-collect_prob)*signalCoupling_[0];
//	    std::cout << "chargetransfereff : " << chargetransfereff << " " << 1-collect_prob << " "<< signalCoupling_[0]<< std::endl;
	    thePixelChargeMap[closer_no] += charge_in_pixel*chargetransfereff;
//	  thePixelChargeMap[pixel_no] += charge_in_pixel*(1-chargetransfereff);
	  }
	}
    }


//  double CHS=0;
//  for(auto &IO: thePixelChargeMap){
//    std::cout << "carica nella mappa post sharing  "<<IO.second << std::endl;
//    CHS += IO.second;
//  }
//  std::cout << "carica TOTALE post sharing  "<<CHS << std::endl;
  return thePixelChargeMap;
}
