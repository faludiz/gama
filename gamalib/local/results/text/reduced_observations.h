/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2002, 2003  Jan Pytel  <pytel@gama.fsv.cvut.cz>
                        2003  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU GaMa / GaMaLib C++ Library.
    
    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
 *  $Id: reduced_observations.h,v 1.1 2003/01/20 18:02:23 cepek Exp $
 */

#ifndef GaMaLib_local_results_text_reduced_observations_h
#define GaMaLib_local_results_text_reduced_observations_h

#include <gamalib/local/results/text/underline.h>
#include <gamalib/local/network.h>
#include <gamalib/local/acord.h>
#include <cctype>
#include <iomanip>

namespace GaMaLib {

template <class OutStream>
void ReducedObservationsText(GaMaLib::LocalNetwork* IS,
    	                GaMaLib::ReducedObservations* reduced, OutStream& out)    
{
   using namespace std;
   using namespace GaMaLib;

   if ( !reduced->size() )
       return;
   
   out << T_GaMa_reduced_Review_of_reduced_observations << "\n"
       << underline(T_GaMa_reduced_Review_of_reduced_observations, '*') 
       << "\n\n";

   if ( reduced->size_nonexist() )
       ; // size_nonexist
   

   int minval_obs = 12;
   int maxval_obs = minval_obs;

   int minval_dh  = 8;
   int maxval_dh  = minval_dh;
   
   {
       for (ReducedObservations::ListReducedObs_iter it = reduced->begin();
	                                            it != reduced->end(); it++)
       {
	   Double value = it->orig_value();
	   Double dh    = it->ptr_obs->to_dh() - it->ptr_obs->from_dh();
	   
	   int oz = 0;
	   int dz = 0;
	   
	   if (value < 0)
           {
             oz = 1;
             value = -value;
           }
	   
	   if (value < 1e5) continue;
	   
	   oz += 6;   // ... decimal point plus 5 digits
	   
	   do {
	       oz++; 
	       value /= 10;
	   } while (value >= 1);
	   
	   if (oz > maxval_obs) maxval_obs = oz;

	   if (dh < 0)
	   {
	       dz = 1;
	       dh = -dh;
	   }

	   if (dh < 1e2) continue;

	   dz += 6; // ... decimal point plus 5 digits

	   do {
	       dz++;
	       dh /= 10;
	   } while (dh >=1);
       }
   }

   // out.width(IS->maxw_obs());
   //   out << "i" << " ";

   out.width(IS->maxw_id());
   out << T_GaMa_standpoint << " ";
   out.width(IS->maxw_id());
   out << T_GaMa_target << "       ";
   out.width(maxval_obs);
   out << T_GaMa_reduced_observed << " ";
   out.width(maxval_obs); 
   out << T_GaMa_reduced_reduced << " ";
   out.width(maxval_dh);
   out << T_GaMa_reduced_header1;

   {  
     int kk = 12 + 2*IS->maxw_id() + maxval_obs - minval_obs;
     for (int i=0; i < kk; i++) out << "=";
   }

   out << T_GaMa_reduced_value;

   {
     int kk = maxval_obs - minval_obs + 1;
     for (int i=0; i < kk; i++) out << "=";
   }  

   out << "==== [m/g] ";

   {
       int kk = maxval_dh - minval_dh;
       for (int i=0; i < kk; i++) out << "=";
   }

   out << "=== [mm] =\n\n";

   out.flush();
   
   
   PointID predcs = "";   // provious standpoint ID

   for (ReducedObservations::ListReducedObs_iter it = reduced->begin();
  	                                        it != reduced->end(); it++)
   {
      Observation* pm = it->ptr_obs;
      PointID cs = pm->from();
      out.width(IS->maxw_id());
      if (cs != predcs)
         out << cs.c_str();
      else
         out << " ";
      out << " ";
      PointID cc = pm->to();
      out.width(IS->maxw_id());
      out << cc.c_str();
      out.setf(ios::fixed, ios::floatfield);

      string str_nonexist("");
      for (int i=0; i < maxval_obs - 2; i++)
	  str_nonexist+="-";
      
      {   // ***************************************************
	 if (S_Distance* sd = dynamic_cast<S_Distance*>(pm))
          {
            out << T_GaMa_s_distance; 
            out.precision(5);
            out.width(maxval_obs);
	    out << it->orig_value() << " ";
            out.width(maxval_obs);
	    if ( it->reduced() )
		out << pm->value();
	    else
		out << str_nonexist.c_str();
            out << " ";
          }

        else if (Z_Angle* za = dynamic_cast<Z_Angle*>(pm))
          {
	    out << T_GaMa_z_angle;
            out.precision(6); 
            out.width(maxval_obs); 
            out << R2G*(it->orig_value()) << " ";
            out.width(maxval_obs);
	    if ( it->reduced() )
		out << R2G*pm->value();
	    else
		out << str_nonexist.c_str();
	    out << " ";
          }
        else if (H_Diff* h = dynamic_cast<H_Diff*>(pm))
          {
            out << T_GaMa_levell;
            out.precision(5);
            out.width(maxval_obs);
            out << it->orig_value() << " ";
            out.width(maxval_obs);
	    if ( it->reduced() )
		out << pm->value();
	    else
		out << str_nonexist.c_str();
	    out << " ";
          }
        else  
          {
            throw GaMaLib::Exception("review/reduced_observations.h - "
                                     "unknown observation type");
          }
      }   // ***************************************************
      out << " ";
      out.precision(4);
      out.width(maxval_dh);
      out << pm->to_dh() - pm->from_dh() << "\n";
      out.flush();

      predcs = cs;  // previous standpoint ID
      }

   out << "\n\n";
   out.flush();
}

}

#endif



