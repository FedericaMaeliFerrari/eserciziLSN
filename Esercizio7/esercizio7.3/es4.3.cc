#include <iostream>
#include <fstream>
#include<stdlib.h>
#include <cmath>
#include <string>
#include "MD.h"
#include "random.h"


using namespace std;

int main(int argc, char **argv){
	if(argc!=5){
	cerr<<"per l'esecuzione del programma" <<argv[0]<< " bisogna inserire: 1)il nome del programma, 2,3,4)i files da cui leggere i valori (input.gas,input.liquid,input.solid) e 3)il file da cui caricare le posizioni iniziali delle particelle (config.0)"<<endl;
	}
	int nblk=100;
	cout<<"QUESTA SIMULAZIONE VERRÀ CONDOTTA UTILIZZANDO LE UNITA' DI MISURA DEL SISTEMA INTERNAZIONALE"<<endl;
	//double epskb=120.;
	
	//double eps=epskb*1.380*pow(10.,-23.);
	
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

//INIZIALIZZO I VETTORI PER LA MEDIA A BLOCCHI
	for(int num=1;num<=3;num++){
		cout<<"SIMULAZIONE DEL FILE "<<num<<endl;
	//INIZIALIZZO I VALORI PER LA DINAMICA MOLECOLARE;
		MD molecole(argv[num],rnd);
		molecole.inputPos(argv[4]);
		molecole.inputV();

	//MOTO DELLE MOLECOLE SECONDO ALGORITMO DI VELVET
		int nstep=molecole.getnstep();
		int iprint=molecole.getiprint();
		//int nconf = 1;

		for(int i=1;i<=nstep;i++){
			molecole.Move();
			if(i%iprint==0){cout<<"Numero di Time-Step: "<<i<<endl;}
			if(i%10==0){
				//molecole.Misura(0);
				//nconf += 1;
			}
		}

		molecole.OldConf();
		molecole.ConfFinal();

	//MISURA CON OPZIONE RESTART
		nstep=1000;
		molecole.inputPosRestart();
		for(int iblk=1;iblk<=nblk;iblk++){
			molecole.Reset(iblk);
			if(iblk%10==0){cout<<"Blocco Simulazione: "<<iblk<<endl;}
			for(int j=0;j<nstep;j++){
				molecole.Move();
				molecole.Misura_blocchi();
				molecole.Accumulate();
				
			}
			molecole.Averages(iblk);
				
		}
		molecole.OldConf();
		molecole.Gavefinal();
	}
		


	return 0;
}
	
