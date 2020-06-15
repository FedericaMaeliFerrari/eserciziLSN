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


/*
	double delta=1.3;
	double x0=0.;
	double mu[20];
	double sigma[20];
	int Nblk=100;
	int Nstep=1000;
//creo dei vettori di mu e sigma da 0 a 2 (vedendo l'andamento del potenziale e della funzione d'onda sicuramente i parametri che minimizzano l'energia non si trovano all'infuori di questo range)	
	for(int i=0;i<20;i++){
		mu[i]=0.1*i;
		sigma[i]=0.1*i;
	}

//creo una matrice con i valori di <H> trovati al variare di mu e sigma 
	double stimaH[20][20];
	Quantum Q(rnd,x0,mu[0],sigma[0],delta);
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			Q.SetParametres(mu[i],sigma[j]);
			for(int iblk=1; iblk <= Nblk; ++iblk){ //Simulation 
		    		Q.Reset(iblk);   //Reset block averages
		    		for(int istep=1; istep <= Nstep; ++istep){
		      			Q.Metropolis(0);
		      			Q.Measure();
				}
		    		Q.Media(iblk,0); 
		  	}
			stimaH[i][j]=Q.GetGlobH()/Nblk;
		}	
	}


//trovo il valore minimo della matrice
	double Hmin=stimaH[0][0];
	int posmumin=0;					//posizione in cui si trova la mu che minimizza <H>
	int possigmamin=0;				//posizione in cui si trova il sigma che minimizza <H>
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			if(stimaH[i][j]<Hmin){
				Hmin=stimaH[i][j];
				posmumin=i;
				possigmamin=j;
			}
		}
	}
	cout<<"mu _min= "<<posmumin*0.1<<endl;
	cout<<"sigma _min= "<<possigmamin*0.1<<endl;
	double mu_min=posmumin*0.1;
	double sigma_min=possigmamin*0.1;

*/

//Faccio partire la simulazione con i parametri che minimizzano <H>

	double mu_min=0.8;
	double sigma_min=0.8;
	double delta=2.95;   //Passo per Metropolis
	double x0=0.;
	int Nblk=101;	   //numero blocchi--> il primo blocco servirà solo per stabilizzare il sistema, non verrà utilizzato nel codice python
	int Nstep=1000;	   //numero step in ogni blocco

	Quantum Qmin(rnd,x0,mu_min,sigma_min,delta);
	for(int iblk=1; iblk <= Nblk; ++iblk){ //Simulation 
    		Qmin.Reset(iblk);   //Reset block averages
    		for(int istep=1; istep <= Nstep; ++istep){
      			Qmin.Metropolis(1);
      			Qmin.Measure();
		}
    		Qmin.Media(iblk,1);   //Print results for current block
  	} 

return 0;
}

