#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "MD.h"

using namespace std;

int main(int argc, char **argv){
	if(argc!=3){
	cerr<<"per l'esecuzione del programma" <<argv[0]<< " bisogna inserire: 1)il nome del programma, 2)il file da cui leggere i valori (rho, num.particelle ecc.) e 3)il file da cui caricare le posizioni iniziali delle particelle"<<endl;
	}
	int opz;
	cout<<"Se si vuole fare la simulazione con l'opzione Restart digitare 1, altrimenti digitare 0"<<endl;
	cin>>opz;

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



//INIZIALIZZO I VALORI PER LA DINAMICA MOLECOLARE;
	MD molecole(argv[1],rnd);
	molecole.inputPos(argv[2]);
	molecole.inputV();


//MOTO DELLE MOLECOLE SECONDO ALGORITMO DI VELVET
	int nstep=molecole.getnstep();
	int iprint=molecole.getiprint();
	//int nconf = 1;

	for(int i=1;i<=nstep;i++){
		molecole.Move();
		if(i%iprint==0){cout<<"Numero di Time-Step: "<<i<<endl;}
		if(i%10==0){
			molecole.Misura(0);
			//nconf += 1;
		}
	}

	molecole.OldConf();
	molecole.ConfFinal();
	
//MOTO MOLECOLE CON OPZIONE RESTART	
	if(opz==1){
		for(int k=1;k<=10;k++){
			cout<<"Ciclo Simulazione: "<<k<<endl;
			molecole.inputPosRestart();
			for(int i=1;i<=1000;i++){
				molecole.Move();
				molecole.Misura(opz);
					
			}
			molecole.OldConf();
		}
	}
	




return 0;
}
