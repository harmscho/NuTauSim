// Déclarations de toutes les constantes utilisées par le programme

const double PI = 3.1415926535898;	// PI :-)

const double alpha=1/137.03599976;	// Constante de structure fine
const double M=0.938271998;		// Masse du proton en GeV
const double m=1.77699;		// Masse du tau en GeV
const double le=3.8616*1e-11;		// Longueur d'onde de Compton de l'électron / 2PI
const double Navo=6.02214199*1e+23;	// Nombre d'Avogadro
const double me=0.510998902*1e-3;	// Masse de l'électron en GeV
const double e=2.718;
const double mpi=0.1349766;		// Masse du PI0 en GeV
const double taudl=86.93*1e-4; 	// longueur de désintégration du tau en cm
const double R0=6.378e8; 		// Rayon de la terre en cm

// Constantes secondaires

const double M2=M*M;
const double m2=m*m;
const double a4=alpha*alpha*alpha*alpha;
const double a2=alpha*alpha;
const double a3=alpha*alpha*alpha;
const double me2=me*me;
const double R02=R0*R0;

// Paramètres Standard rock

const double A=22.;		//
const double Z=11.;		//
const double X0=0.049;	//
const double X1=3.055;	//
const double aa=0.083;	//   Paramètres pour le calcul de l'ionisation
const double mm=3.412;	//
const double I=136.4*1e-9;	//
const double I2=I*I;		//
const double CC=-3.774;	//

const double densrock=2.65; 	// densité moyenne du milieu en g/cm^3

// Constantes utilisées dans le calcul de la perte par ionisation

const double Cbb1=a2*2*PI*Navo*le*le*Z*me/A;
const double Cbb2=2*me;

// Paramètres de la paramétrisation pour le beta total (effet photonucléaire + brem + production de paire): 
// Polynome de degré 9 entre 100 GeV et 10e14 GeV  (pour la roche standard)

const double p0=3.45062e-07;   
const double p1=-3.72936e-07;  
const double p2=2.05637e-07; 
const double p3=-4.84147e-08;
const double p4=5.07199e-09;   
const double p5=1.94132e-11;   
const double p6=-6.13689e-11;   
const double p7=6.59919e-12;
const double p8=-3.09033e-13;   
const double p9=5.86727e-15;


// Paramètres pour le calcul de la section efficace totale de l'interaction CC: 
// Polynome de degré 7 entre 10e5 et 10e12 GeV

const double pCC0=2.20473e-30;   
const double pCC1=-2.25731e-30;  
const double pCC2=9.73692e-31; 
const double pCC3=-2.29314e-31;
const double pCC4=3.18342e-32;   
const double pCC5=-2.60396e-33;   
const double pCC6=1.16138e-34;   
const double pCC7=-2.17449e-36;   


// Paramètres pour le calcul de la section efficace totale de l'interaction NC: 
// Polynome de degré 9 entre 10e5 et 10e12 GeV

const double pNC0=-1.96378e-30;   
const double pNC1=1.79257e-30;  
const double pNC2=-6.17758e-31; 
const double pNC3=7.85698e-32;
const double pNC4=7.90663e-33;   
const double pNC5=-4.20927e-33;   
const double pNC6=6.44764e-34;   
const double pNC7=-5.11796e-35;   
const double pNC8=2.12964e-36;   
const double pNC9=-3.68115e-38;   

// Parameterization values for tau eloss pair production and bremsstrahlung
// There are 3 layer compositions (iron, rock, water) and 4 parameters for the fit
const double pBrem[3][4] = {{1.56995122e-08, 4.28188407e+00, 2.94450001e+00, 6.98412618e-01}, 
			    {7.88982596e-09, 4.48475227e+00, 3.08880391e+00, 6.09300006e-01}, 
                            {6.21177609e-09, 4.64573538e+00, 3.21763680e+00, 4.94713279e-01}
                            };
const double pPair[3][4] = {{2.81202672e-07, 5.42597316e+00, 5.79540809e+00, 3.57248989e-02}, 
                            {1.34819038e-07, 1.60062940e+00, 1.52179335e+00, 8.98529667e+00}, 
                            {9.26738275e-08, 1.70211443e+00, 1.57079809e+00, 7.81591550e+00}};

