#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "random.h"

using namespace std;


/************************************************
Programma che genera un set di numeri casualitra 0 e 1
e stima la loro media e incertezza e la media degli errori con la loro relativa incertezza
a seconda del numero di blocchi considerati per il calcolo N.
In ingresso viene chiesto quanti numeri casuali generare(M) e in quanti blocchi suddividerli
per calcolare le medie(N). Stampa infine in un file di testo "Dati.dat" quanto calcolato.
*************************************************/


int main (int argc, char *argv[]){

   int M;
   int N;
   cout<<"Inserire quanti numeri casuali si intende generare (numero intero)"<<endl;
   cin>>M;
   cout<<"Inserire in quanti blocchi si vuole suddividere i valori per il calcolo della media (numero intero)"<<endl;
   cin>>N;
   double vettR[M];

//SALVA LE VARIABILI DA FILE PRIMES E SEED.IN PER GENERARE NUMERI CASUALI TRA 0 E 1
   Random rnd;
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
	    rnd.SetRandom(seed,p1,p2);
	 }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


// GENERA M NUMERI TRA 0 E 1
   for(int i=0; i<M; i++){
      vettR[i]=rnd.Rannyu();
   }

   rnd.SaveSeed();

// CALCOLA LA MEDIA IN BLOCCHI DI 100 LANCI, LA SOMMA DELLE MEDIE PROGRESSIVE E L'ERRORE TRAMUTE IL METODO A BLOCCHI. STESSA COSA PER MEDIE DEGLI ERRORI E RELATIVA INCERTEZZA
   int L=M/N;
   double vettMedie[N];		// <A_i>
   double vettMedie_err[N];
   double vettMedie2[N];	// <A_i>^2
   double vettMedie2_err[N];
   double sum_prog[N];		//vett con somme progressive delle medie
   double sum_prog_err[N];
   double sum_prog2[N];		// vett con somme progressive delle medie al quadrato
   double sum_prog2_err[N];
   double err[N];		//vett con errori delle medie
   double err_err[N];

   for(int i=0;i<N;i++){
	double sum=0;
	double sum_err=0;
	for(int j=0;j<L;j++){
		int k=j+i*L;
		sum=sum+vettR[k];
		sum_err=sum_err+pow((vettR[k]-0.5),2.);
	}
	vettMedie[i]=sum/L;
	vettMedie_err[i]=sum_err/L;
	vettMedie2[i]=pow(vettMedie[i],2.);
	vettMedie2_err[i]=pow(vettMedie_err[i],2.);
   }

   ofstream outData;		
   outData.open ("Dati.dat");	// apro un file di testo in cui andrò a salvare la somma delle medie progressive e l'errore associato, la media degli errori e l'icertezza associata

   for(int i=0;i<N;i++){
	for(int j=0;j<i+1;j++){
		sum_prog[i]+=vettMedie[j];
		sum_prog_err[i]+=vettMedie_err[j];
		sum_prog2[i]+=vettMedie2[j];
		sum_prog2_err[i]+=vettMedie2_err[j];
	}
	sum_prog[i]=sum_prog[i]/(i+1);
	sum_prog_err[i]=sum_prog_err[i]/(i+1);
	sum_prog2[i]=sum_prog2[i]/(i+1);
	sum_prog2_err[i]=sum_prog2_err[i]/(i+1);
	if(i==0){
		err[i]=0;
		err_err[i]=0;
	}
	else{
		err[i]=sqrt((sum_prog2[i]-pow(sum_prog[i],2.))/i);
		err_err[i]=sqrt((sum_prog2_err[i]-pow(sum_prog_err[i],2.))/i);
	}
	outData<<setprecision(6);
	outData<<sum_prog[i]<<setw(12)<<err[i]<<setw(12)<<sum_prog_err[i]<<setw(12)<<err_err[i]<<endl;
	
   }
   outData.close();










   return 0;
}

