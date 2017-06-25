///////////////////////////////////////////////////////////////////////////
// here you find various specialized table classes for parametrization of 
// crossections, fluxes (both table) and final states (TableFinal)
///////////////////////////////////////////////////////////////////////////

#ifndef TABLE
#define TABLE
#include <fstream>
#include <cstdlib>
#include <iostream>

using namespace std; 

//////////////////////////////////////////////////////////////////////
//  table class - loads a table from a file and provides            //
//  reading and interpolation utilities                             //
//////////////////////////////////////////////////////////////////////
 class Table {
   public:
      Table(){}
     ~Table();
      void   InitTable (char *filename);  
      double Evaluate  (double x, double y);  // function to evaluate 2dim table
      double Evaluate  (double x);            // function to evaluate 1dim table
      int    IsInit()  { return Init;}  
   private:
      int      Init;       // Initialized 
      int      Imax;       // size of table a[i][j]
      int      Jmax;
      double   Xmin;       // value of i=0
      double   Ymin;       // value of j=0
      double   Xmax;       // value of i=imax
      double   Ymax;       // value of j=jmax
      double **Tablearray;
      char    *Inputfile;
 };

/////////////////////////////////////////////////////////////////////
// FinalTable assumes a table of final states, of the structure     /
// [energy][entriesperenergy][2], where [2] could be theta and y.   /
// ThrowFinal(energy) draws one of the final states                 /
/////////////////////////////////////////////////////////////////////
 class FinalTable {
   public:
      FinalTable();
     ~FinalTable();
      int  IsInit()  {return Init;}  
      int  InitTable (char *filename);  
      void ThrowFinal(double energy, double final[]);  
      void ThrowFinal(double final[]);  
   private:
      int       Init;      // Initialized?
      int       Imax;      // number of steps in energy
      int       Nfinal;    // number of final states per energy
      int       Ndim;
      double    Emin;
      double    Emax;
      double ***Tablearray;
 };

//////////////////////////////////////////////////////////////////////
#endif 









