//////////////////////////////////////////////////////////////////////
//  table class - loads a table from a file and provides various    //
//  reading utilities                                               //
//////////////////////////////////////////////////////////////////////
 
#include "Table.hh"

//////////////////////////////////////////////////////////////////////
 void Table::InitTable(char *file_name){ 
    ifstream tablefile(file_name);   // opening input file
    if (!tablefile){
      cerr << " Failed to open file "<<file_name<<" for reading table" << endl;
      exit(EXIT_FAILURE);
    }
    tablefile >> Imax >> Jmax;
    tablefile >> Xmin >> Xmax;

    if (Jmax > 1) {   // if 2 dimensional table
      tablefile >> Ymin >> Ymax;
    } else  {         // 1 dimensional table
      Ymin = 0;       // filled with dummy values
      Ymax = 1;       // filled with dummy values
    }

    Tablearray = new double*[Imax];   // allocate memory for table
    for (int i=0; i<Imax; i++)
      Tablearray[i] = new double[Jmax]; 

    for (int j=0; j<Jmax; j++)
      for (int i=0; i<Imax; i++)
        tablefile >> Tablearray[i][j];
      
    Init = 1;   // set status
    if (tablefile.fail()){
      cerr << " something has gone wrong with reading file " << endl;
      exit(EXIT_FAILURE);
    }
    tablefile.close();
}

//////////////////////////////////////////////////////////////////////
 Table::~Table(){
    if (Init){
	for (int i=0; i<Imax; i++)
	    delete[] Tablearray[i];
	delete[] Tablearray;
    }
 }    

//////////////////////////////////////////////////////////////////////
// reading and interpolating a table
 double Table::Evaluate(double x, double y){
    int    ilow, ihigh, jlow, jhigh;
    double djlow, dilow;
    double eval;

    dilow = (x-Xmin)/(Xmax-Xmin)*(Imax - 1);  // calculating lower index 
    djlow = (y-Ymin)/(Ymax-Ymin)*(Jmax - 1);  
    ilow  = int (dilow);	       
    jlow  = int (djlow);	       
    ihigh = ilow + 1;                          
    jhigh = jlow + 1;
    if (dilow == Imax-1)  // if x exactly on upper boundary 
         ihigh = ilow;
    if (djlow == Jmax-1)  // if y exactly on upper boundary 
         jhigh = jlow;
    if (dilow < 0 || djlow < 0 || dilow > Imax-1 || djlow > Jmax-1) {
      eval = 0;           //checking if index outside array boundaries 
    }  else {             // interpolating table
// distance to interpolation 
        double di = dilow - ilow;  
        double dj = djlow - jlow;
// interpolating the values
        eval  =  Tablearray[ilow][jlow];
        eval += (Tablearray[ihigh][jlow]-Tablearray[ilow][jlow])*di;
        eval += (Tablearray[ilow][jhigh]-Tablearray[ilow][jlow])*dj;
   }
   return eval;
 }

//////////////////////////////////////////////////////////////////////
 double Table::Evaluate(double x) {  // for one dimensional tables
     return Evaluate(x,Ymin);
 }

//////////////////////////////////////////////////////////////////////
  FinalTable::FinalTable() {
        Init = 0;
 }
//////////////////////////////////////////////////////////////////////
 int  FinalTable::InitTable(char *infile){
     ifstream tablefile(infile);      // opening input file
     if (!tablefile){
         cerr << " Failed to open file for reading final table" << endl;
	 exit(EXIT_FAILURE);
     }

    tablefile >> Nfinal;              // number of final states per energy
    tablefile >> Ndim;                // number of dimensions of final state
    tablefile >> Imax;                // number of energy steps and
    tablefile >> Emin >> Emax;

    Tablearray = new double**[Imax];  // allocate memory for table
    for (int i=0; i < Imax; i++){
	Tablearray[i] = new double*[Nfinal];
	for (int j=0; j<Nfinal; j++)
	    Tablearray[i][j] = new double[Ndim];
    }
    for (int i=0;i<Imax;i++){
	for (int j = 0; j < Nfinal; j++){
	    for (int k=0; k < Ndim; k++)
		tablefile >> Tablearray[i][j][k];
	}
    }
    Init = 1;   // set status
    if (tablefile.fail()) {
	cerr << " something has gone wrong with reading the file" << endl;
	exit(EXIT_FAILURE);
    }
    tablefile.close();
    return Init;
 }    

//////////////////////////////////////////////////////////////////////
 void  FinalTable::ThrowFinal(double energy, double final[]) {
    int    ilow, energyindex;
    double dilow;
    int    randomentry;
    dilow = (energy-Emin)/(Emax-Emin)*(Imax - 1);  // calculating lower index  
    ilow = int (dilow);
    if (ilow<0){
	energyindex = 0;
    } else if (ilow >= Imax-1){
	energyindex = ilow;
    } else {
	energyindex = ilow+1;
	if (drand48() > dilow - ilow)
	    energyindex = ilow;
    }
    if (energy <= Emax && energy >= Emin){
	randomentry = int (drand48()*Nfinal);
	for (int i=0; i < Ndim; i++) {
	    final[i]=Tablearray[energyindex][randomentry][i];
	}
    } else {
	for (int i=0; i < Ndim; i++) 
	    final[i] = 0;    // set finalstate to 0
    }
 }

//////////////////////////////////////////////////////////////////////
 void  FinalTable::ThrowFinal(double final[]){
   int randomentry;
   int energyindex=0;
   randomentry=int (drand48()*Nfinal);
   for (int i=0; i < Ndim; i++) 
     final[i] = Tablearray[energyindex][randomentry][i];
 } 

//////////////////////////////////////////////////////////////////////
 FinalTable::~FinalTable(){
     if (Init){
	 for (int i=0; i < Imax; i++){
	    for (int j=0; j < Nfinal; j++)
		delete[] Tablearray[i][j];
	    delete[] Tablearray[i];
	 }
         delete[] Tablearray;
     }
 }    
//////////////////////////////////////////////////////////////////////












