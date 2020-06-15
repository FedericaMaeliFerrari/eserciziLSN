#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Ising.h"
#include "random.h"

using namespace std;

int main(){

//SALVA LE VARIABILI DA FILE PRIMES E SEED.IN PER GENERARE NUMERI CASUALI TRA 0 E 1
	Random *rnd=new Random();
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
   	if (Primes.is_open()){
      	Primes >> p1 >> p2 ;
   	} else cerr << "PROBLEM: Unable to open Primes" << endl;
   	Primes.close();

   	ifstream input("seed.in");
   	string property;
   	if (input.is_open()){
      	while ( !input.eof() ){
	 	input >> property;
	 	if( property == "RANDOMSEED" ){
	    	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	    	rnd->SetRandom(seed,p1,p2);
	 	}
      	}
      	input.close();
   	} else cerr << "PROBLEM: Unable to open seed.in" << endl;


//PROGRAMMA ISING METROPOLIS
	Ising I(rnd);	 //Inizialization
	I.SetSpin(0);    	//Genero il vettore di spin iniziale in modo randomico
	int nblk=I.GetNblk();
	int nstep=I.GetNstep();
	int metro=I.GetMetro();
	double temp=I.GetTemp();
	double Temp=I.GetTemp();   //mi servirà per riinizializzare la classe Ising con il valore di temperature iniziale per la simulazione con Gibbs

	ofstream OutT;		//File in cui salverò i valori della temperatura con cui ho effettuato le misure
	OutT.open("temperatura.dat");
	cout<<"ALGORITMO DI METROPOLIS"<<endl;
	while(temp < 2.05){
	  	for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
	    		I.Reset(iblk);   //Reset block averages
	    		for(int istep=1; istep <= nstep; ++istep){
	      			I.Move(metro);
	      			I.Measure();
	      			I.Accumulate(); //Update block averages
	    		}
	    		I.Averages(iblk,metro);   //Print results for current block
	  	}
	  	I.ConfFinal(); //Write final configuration
		I.SetSpin(1);  //Uso gli spin della configurazione precedente (non li genero più randomicamente)
		OutT<<temp<<endl;
		cout<<" LA TEMPERATURA E': "<<temp<<endl;
		temp=temp+0.05;
		I.SetTemp(temp);	
		
	}
	OutT.close();


//PROGRAMMA ISING GIBBS
	I.SetSpin(0);	 //Inizialization
	metro=0;
	I.SetMetro(0);
	I.SetTemp(Temp);  // riinizializzo la temperature con il valore iniziale del file input
	temp=I.GetTemp();
	cout<<"°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°"<<endl;
	cout<<"ALGORITMO DI GIBBS"<<endl;
	while(temp < 2.05){
	  	for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
	    		I.Reset(iblk);   //Reset block averages
	    		for(int istep=1; istep <= nstep; ++istep){
	      			I.Move(metro);
	      			I.Measure();
	      			I.Accumulate(); //Update block averages
	    		}
	    		I.Averages(iblk,metro);   //Print results for current block
	  	}
	  	I.ConfFinal(); //Write final configuration
		I.SetSpin(1);
		cout<<" LA TEMPERATURA E': "<<temp<<endl;
		temp=temp+0.05;
		I.SetTemp(temp);
	}


//SELEZIONO I VALORI CHE ANDRANNO PLOTTATI, DUNQUE SEPARO GIBBS DA METROPOLIS E SELEZIONO LE MEDIE DOPO 20 BLOCCHI
	double e[1][5];
	double c[1][5];
	double x[1][5];
	double m[1][5];
	
	ifstream E,C,X,M;
	E.open("output.ene.0");
	C.open("output.heat.0");
	X.open("output.chi.0");
	M.open("output.mag.0");
	ofstream Out,OutG;
	Out.open("datiMetr.dat");
	OutG.open("datiGibbs.dat");
	cout<<"°°°°°°°°°°°°°°°°°°°°°°°°°°°"<<endl;
	cout<<"SCRITTURA SU FILE"<<endl;
	for(int i=0;i<1240;i++){
		E>>e[0][0]>>e[0][1]>>e[0][2]>>e[0][3]>>e[0][4];
		C>>c[0][0]>>c[0][1]>>c[0][2]>>c[0][3]>>c[0][4];
		X>>x[0][0]>>x[0][1]>>x[0][2]>>x[0][3]>>x[0][4];
		M>>m[0][0]>>m[0][1]>>m[0][2]>>m[0][3]>>m[0][4];
		
		if(e[0][0]==20){
			if(e[0][4]==1){
				Out<<e[0][2]<<" "<<e[0][3]<<" "<<c[0][2]<<" "<<c[0][3]<<" "<<x[0][2]<<" "<<x[0][3]<<" "<<m[0][2]<<" "<<m[0][3]<<endl;
			}
			else{
				OutG<<e[0][2]<<" "<<e[0][3]<<" "<<c[0][2]<<" "<<c[0][3]<<" "<<x[0][2]<<" "<<x[0][3]<<" "<<m[0][2]<<" "<<m[0][3]<<endl;
			}
		}
	}
	E.close();
	C.close();
	X.close();
	M.close();
	Out.close();
	OutG.close();
	

  return 0;
}

