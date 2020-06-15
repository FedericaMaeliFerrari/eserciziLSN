#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "funzioneBase.h"
#include "funzione.h"
#include "integrale.h"

using namespace std;

int main(int argc, char *argv[]){
	int M;
	int N;
	double a;
	double b;
	cout<<"Inserire quanti numeri generare il numero di intervalli per calcolare la media progressiva"<<endl;
	cin>>M>>N;
	cout<<"Inserire estremi integrale"<<endl;
	cin>>a>>b;
	int L=M/N;

//SALVA LE VARIABILI DA FILE PRIMES E SEED.IN PER GENERARE NUMERI CASUALI TRA 0 E 1
	Random *rnd;
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

// CREO I PUNTATORI ALLE CLASSI E UTILIZZO IL METODO INTEGRALEBYMEDIA PER CALCOLARE IL VALORE DELL'INTEGRALE E DELL'ERRORE ASSOCIATO SIGMA
	coseno *funz;
	coseno *funzerr;
	funz->cos_costruttore(M_PI/2.,M_PI/2.,1.);          //funzione da utilizzare per il calcolo dell'integrale
	funzerr->cos_costruttore(M_PI/2.,M_PI/2.,2.);       //funzione da utilizzare per il calcolo dell'errore (g(x)^2)
	integrale I(a,b,funz,rnd);
	integrale Ierr(a,b,funzerr,rnd);
	double vett[N];
	double vett2[N];

	for (int i=0;i<N;i++){
		vett[i]=I.integraleByMedia(a,b,L);
		vett2[i]=Ierr.integraleByMedia(a,b,L);	
	}


//CALCOLO LE MEDIE PROGRESSIVE E GLI ERRORI E LI STAMPO SU UN FILE	
	ofstream outData;
	outData.open ("Integrale_unif.dat");

	double I_prog[N];
	double Ierr_prog[N];
	double err[N];
	for(int i=0;i<N;i++){
		for(int j=0;j<i+1;j++){
			I_prog[i]+=vett[j];
			Ierr_prog[i]+=vett2[j];
		}
		I_prog[i]=I_prog[i]/(i+1);
		Ierr_prog[i]=Ierr_prog[i]/(i+1);
		if (i==0){err[i]=0;}
		else{err[i]=sqrt((Ierr_prog[i]-pow(I_prog[i],2.))/i);}
		outData<<I_prog[i]<<" "<<err[i]<<endl;
	}
	
	



	return 0;
}
