////////////////////////////////////////////////////////////////////////
// Implementation of the class Earth
////////////////////////////////////////////////////////////////////////
//#include "CLHEP/config/CLHEP.h"
//#include "CLHEP/config/iostream.h"
//#include "CLHEP/config/fstream.h"
#include "Earth.hh"
#include "math.h"
#include "stdio.h"

////////////////////////////////////////////////////////////////////////
//  Earth::Earth(){ }
////////////////////////////////////////////////////////////////////////
//  Distribution of the matter density in the Earth measured in g/cm^3
//  as the function of distance from the center R measured in meters
//  Earth model was referred in Gandhi et al.  astro-ph/9512364 p.23
////////////////////////////////////////////////////////////////////////
double Earth::GetDensity(double R_cm){
     double dens = 0;
     double R_meters = R_cm/100.;
     double x = R_meters/EARTHRADIUS;
     double h = R_meters/EARTHRADIUS-1.;
     double depth_new_layer_meters = depth_new_layer*1.e3;  // meters
       /*if (R_meters < 0) {
        cerr << "  Attention, R < 0 !!! " << endl;
        exit(EXIT_FAILURE); 
       }*/


     //double dens_new_layer  = water_layer_density_g_cm3;  // density of new layer (g/cm^3) 
     //double depth_new_layer = water_layer_thickness_km;   // depth of new layer (km)
     //double dens_new_layer=1.027;  // density of new layer (g/cm^3) 
     //double depth_new_layer=3.688;   // depth of new layer (km)
                                   // maximum value = 22 km 
     //printf("depth_new_layer %f\n", depth_new_layer);

       if (R_meters < 1221500) dens = 13.0885 - 8.8381*x*x;
       else if (R_meters < 3480000) dens = 12.5815 - x*(1.2638 + x*(3.6426 + x*5.5281));
       else if (R_meters < 5701000) dens = 7.9565 - x*(6.4761 - x*(5.5283 - x*3.0807));
       else if (R_meters < 5771000) dens = 5.3197 - 1.4836*x;
       else if (R_meters < 5971000) dens = 11.2494 - 8.0298*x;
       else if (R_meters < 6151000) dens = 7.1089 - 3.8045*x;
       else if (R_meters < 6346600) dens = 2.691 + 0.6924*x;
       else if (R_meters < 6356000) dens = 2.9;
//old       else if (R_meters < (EARTHRADIUS - 3000) ) dens = 2.6;
       else if (R_meters < (EARTHRADIUS - depth_new_layer_meters) ) dens = 2.6;
       else if (R_meters <= EARTHRADIUS) dens = dens_new_layer;
//old       else if (R_meters <= EARTHRADIUS) dens = 1.02;
       else if (R_meters <= ATMORADIUS)  dens = 1.29e-3*exp(-h/8.5e3);

       /*else if (R_meters > ATMORADIUS){
       cerr << " Attention, R = " << R_meters 
            << " is higher than ATMOSRADIUS =" << ATMORADIUS 
            << " m !!! " << endl;
       cerr << " R - ATMORADIUS = " << R_meters - ATMORADIUS 
            << " m" << endl;
       exit(EXIT_FAILURE);
       }*/

     return dens;
 }


////////////////////////////////////////////////////////////////////////




