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
   int Run;
   cout<<"Inserire in quanti blocchi si vuole suddividere l'intervallo [0,1]"<<endl;
   cin>>M;
   cout<<"Inserire quanti numeri casuali si intende generare (numero intero)"<<endl;
   cin>>N;
   cout<<"Inserire quante volte si vuole ripetere l'operazione"<<endl;
   cin>>Run;

   double vettNumCas[N];
   double vettChi[Run];

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


// VETTORE CONTENENTE GLI INTERVALLI CON CUI HO DIVISO L'INTERVALLO [0,1]
   double m=1./M;
   double num=0.;
   double vettPos[M];
   for(int i=0;i<M;i++){
	vettPos[i]=num+m;
	num=num+m;		
   }


   ofstream outData;		
   outData.open ("DatiChi.dat");           //APRO UN FILE CHE CONTERRÀ I VALORI DI CHI QUADRO PER OGNI RIPETIZIONE

   for(int j=0; j<Run; j++){
	double sumChi=0.;
	for(int i=0; i<N; i++){
		vettNumCas[i]=rnd.Rannyu();       // GENERA M NUMERI TRA 0 E 1
	}
	for(int k=0; k<M; k++){
		int cont=0;
		for(int w=0; w<N; w++){
			if(k==0){
				if(vettNumCas[w]<=vettPos[k]){        //CONTA QUANTI NUMERI SONO CADUTI ALL'INTERNO DI OGNI INTERVALLO
					cont=cont+1;
				}
				
			}
			else{
				if(vettNumCas[w]>vettPos[k-1] && vettNumCas[w]<=vettPos[k]){
					cont=cont+1;
				}
			}
		
		}	
		sumChi=sumChi+(pow((cont-(N/M)),2)/(N/M));		//CALCOLA IL CHI QUADRO
	}
	vettChi[j]=sumChi;
	outData<<vettChi[j]<<endl;
   }

   outData.close();
   rnd.SaveSeed();


   return 0;
}

