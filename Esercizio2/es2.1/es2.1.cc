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
	cout<<"Inserire quanti numeri generare (M) e il numero di intervalli (N) per calcolare la media progressiva"<<endl;
	cin>>M>>N;
	cout<<"Inserire estremi integrale"<<endl;
	cin>>a>>b;
	int L=M/N;

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

// CREO I PUNTATORI ALLE CLASSI E UTILIZZO IL METODO INTEGRALEBYMEDIA PER CALCOLARE IL VALORE DELL'INTEGRALE E DELL'ERRORE ASSOCIATO SIGMA
	funzioneBase *funz=new coseno(M_PI/2.,0.5,1.);        //funzione da utilizzare per il calcolo dell'integrale
	funzioneBase *funzerr=new coseno(M_PI/2.,0.5,2.);     //funzione da utilizzare per il calcolo dell'errore (g(x)^2)  
	funzioneBase *dx_trasf=new retta(-1.,1.,0.5,1.,-1.);
	funzioneBase *f=new cosretta(-2.,2.,1.,0.,1.,M_PI/2.,0.5,1.);
	funzioneBase *ferr_dx=new cosretta(-2.,2.,2.,0.,1.,M_PI/2.,0.5,2.);           
	integrale I(a,b,funz,rnd);
	integrale Idx(a,b,f,rnd);
	integrale Ierr(a,b,funzerr,rnd);
	integrale Idxerr(a,b,ferr_dx,rnd);
	double vett[N];
	double vett2[N];
	double vettdx[N];
	double vett2dx[N];

	for (int i=0;i<N;i++){
		vett[i]=I.integraleByMedia(a,b,L);
		vettdx[i]=Idx.integraleByMedia_distrib(a,b,L,dx_trasf);
		vett2[i]=Ierr.integraleByMedia(a,b,L);
		vett2dx[i]=Idxerr.integraleByMedia_distrib(a,b,L,dx_trasf);
		cout<<vett[i]<<" "<<vettdx[i]<<" "<<vett2[i]<<" "<<vett2dx[i]<<endl;
	}
	

//CALCOLO LE MEDIE PROGRESSIVE E GLI ERRORI E LI STAMPO SU UN FILE	
	ofstream outData;
	outData.open ("Integrale.dat");

	double I_prog[N]={0};
	double Idx_prog[N]={0};
	double Ierr_prog[N]={0};
	double Idxerr_prog[N]={0};
	double err[N]={0};
	double errdx[N]={0};
	for(int i=0;i<N;i++){
		for(int j=0;j<i+1;j++){
			I_prog[i]+=vett[j];
			Idx_prog[i]+=vettdx[j];
			Ierr_prog[i]+=vett2[j];
			Idxerr_prog[i]+=vett2dx[j];
		}
		I_prog[i]=I_prog[i]/(i+1);
		Ierr_prog[i]=Ierr_prog[i]/(i+1);
		Idx_prog[i]=Idx_prog[i]/(i+1);
		Idxerr_prog[i]=Idxerr_prog[i]/(i+1);
		if (i==0){
			err[i]=0;
			errdx[i]=0;
		}
		else{
			err[i]=sqrt((Ierr_prog[i]-pow(I_prog[i],2.))/i);
			errdx[i]=sqrt((Idxerr_prog[i]-pow(Idx_prog[i],2.))/i);

		}
		outData<<I_prog[i]<<" "<<err[i]<<" "<<Idx_prog[i]<<" "<<errdx[i]<<endl;
	}
	outData.close();



	
	

	return 0;
}
