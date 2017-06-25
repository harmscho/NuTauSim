/////////////////////
// Earth mode      //
/////////////////////

#ifndef EARTH
#define EARTH
#define EARTHRADIUS 6.378e6       // Earth radius in meter
#define ATMORADIUS  6.478e6       // Atmos + radius in meter

using namespace std;

/////////////////////////////////
  class Earth {
//private:

  public:
    double dens_new_layer;  // density of new layer (g/cm^3) 
    double depth_new_layer;   // depth of new layer (km)
    Earth(double water_layer_thickness_km, double water_layer_density_g_cm3){
        dens_new_layer  = water_layer_density_g_cm3;  // density of new layer (g/cm^3) 
        depth_new_layer = water_layer_thickness_km*1.e3;   // depth of new layer (km)
    }

    ~Earth(){}

    double GetDensity(double R);
    double earthdens(double distance_along_chord_cm, double chord_length_cm); // this is the distance along the chord and the chord length
  };
                                   
#endif 
