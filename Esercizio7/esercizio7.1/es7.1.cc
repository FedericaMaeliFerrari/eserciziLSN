#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "NVT.h"
#include "random.h"

using namespace std;

int main(int argc, char **argv){
	if(argc!=4){
	cerr<<"per l'esecuzione del programma" <<argv[0]<< " bisogna inserire: 1)il nome del programma, 2)i files da cui leggere i valori di input: input.dat,input.liquid,input.solid" <<endl;
	}
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

//SIMULAZIONE INSIEME CANONICO GAS
	cout<<"SIMULAZIONE GAS"<<endl<<endl;	
	NVT gas(argv[1],rnd,1); //Inizialization
  	//int nconf = 1;
	int nblk=gas.GetNblk();
	int nstep=gas.GetNstep();
  	for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation 
    		gas.Reset(iblk);   //Reset block averages
    		for(int istep=1; istep <= nstep; ++istep){
      			gas.Move();
      			gas.Measure(1);
      			gas.Accumulate(); //Update block averages
      			//if(istep%10 == 0){
        			//Ngas.ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        			//nconf += 1;
      			//}
		}
    		gas.Averages(iblk);   //Print results for current block
  	}
  	gas.ConfFinal(); //Write final configuration

//SIMULAZIONE INSIEME CANONICO LIQUIDO
	cout<<"°°°°°°°°°°°°°°°°°°°°°°°°°°°°°"<<endl<<endl<<endl;
	cout<<"SIMULAZIONE LIQUIDO"<<endl<<endl;
	NVT liquid(argv[2],rnd,2);
	nblk=liquid.GetNblk();
	nstep=liquid.GetNstep();
  	for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation 
    		liquid.Reset(iblk);   //Reset block averages
    		for(int istep=1; istep <= nstep; ++istep){
      			liquid.Move();
      			liquid.Measure(2);
      			liquid.Accumulate(); //Update block averages
      			
		}
    		liquid.Averages(iblk);   //Print results for current block
  	}
  	liquid.ConfFinal(); //Write final configuration


//SIMULAZIONE INSIEME CANONICO SOLIDO
	cout<<"°°°°°°°°°°°°°°°°°°°°°°°°°°°°°"<<endl<<endl<<endl;
	cout<<"SIMULAZIONE SOLIDO"<<endl<<endl;
	NVT solid(argv[3],rnd,3);
	nblk=solid.GetNblk();
	nstep=solid.GetNstep();
  	for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation 
    		solid.Reset(iblk);   //Reset block averages
    		for(int istep=1; istep <= nstep; ++istep){
      			solid.Move();
      			solid.Measure(3);
      			solid.Accumulate(); //Update block averages
      			
		}
    		solid.Averages(iblk);   //Print results for current block
  	}
  	solid.ConfFinal(); //Write final configuration



return 0;
}
