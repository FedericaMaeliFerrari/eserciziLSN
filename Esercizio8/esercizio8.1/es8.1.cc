#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Quantum.h"
#include "random.h"

using namespace std;

int main(){
	
//SALVA LE VARIABILI DA FILE PRIMES E SEED.IN PER GENERARE NUMERI CASUALI TRA 0 E 1
	Random *rnd=new Random();
	int p1, p2;
	int seed[4];
   	ifstream Primes("Primes");
   	Primes >> p1 >> p2 ;
   	Primes.close();

   	ifstream input("seed.in");
   	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   	rnd->SetRandom(seed,p1,p2);
   	input.close();

//Inizializzazione parametri
	double mu=1.;
	double sigma=0.5;
	double delta=1.3;   //Passo per Metropolis
	double x0=0.;
	int Nblk=101;	   //numero blocchi--> il primo blocco servirà solo per stabilizzare il sistema, non verrà utilizzato nel codice python
	int Nstep=1000;	   //numero step in ogni blocco

	Quantum Q(rnd,x0,mu,sigma,delta);
	for(int iblk=1; iblk <= Nblk; ++iblk){ //Simulation 
    		Q.Reset(iblk);   //Reset block averages
    		for(int istep=1; istep <= Nstep; ++istep){
      			Q.Metropolis();
      			Q.Measure();
		}
    		Q.Media(iblk);   //Print results for current block
  	}
 

return 0;
}

