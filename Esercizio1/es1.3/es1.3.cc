#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"


using namespace std;

int main (int argc, char *argv[]){

	int l;          //lunghezza bastoncino
	int d;          //distanza tra due linee
	int M;          // valori da generare random
	int N;          // numero di blochhi
	cout<<"Inserire lunghezza del bastoncino (l) e distanza tra le righe (d)."<<endl<<" Ricorda che d non deve essere troppo più grande di l e che devono essere numeri interi."<<endl;
	cin>>l>>d;
	cout<<"Inserire quanti valori generare (M) e il nuemro di blocchi da considerare per il calcolo della media(N)."<<endl;
	cin>>M>>N;


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

	double vettPos[M];    //vettore che tiene conto in che coordinata y del piano è caduto il baricentro del legnetto
	double vettL[M];      //vettore che tiene conto con che inclinazione è caduto il legnetto sul piano.
	double Nhit[M];
	

// GENERA M NUMERI TRA 0 E 1, CALCOLA LA POSIZIONE DEL BARICENTRO DEL LEGNETTO E LA POSIZIONE DI UNO DEGLI ESTREMI RISPETTO AD ESSO (PROBLEMA SIMMETRICO).
// CALCOLA DUNQUE SE, CON LE COORDINATE OTTENUTE, IL LEGNETTO HA COLPITO LA RIGA (Nhit=1) O SE NON L'HA COLPITA (Nhit=0).
	for(int i=0; i<M; i++){
		double ybar=rnd.Rannyu(0.,d);
		vettPos[i]=ybar;     //posizione in cui è caduto il baricentro del legnetto
		int cont=0;
		while(cont==0){	    			   //metodo accetpt-reject per calcolare l'inclinazione del legnetto: genero numeri all'interno del quadrato di lato l/2
			double x=rnd.Rannyu(0.,l/2.);      //in cui il baricentro si trova in uno degli angoli. Se il punto generato si trova all'interno della circonferenza di raggio l/2
			double y=rnd.Rannyu(0.,l/2.);      // allora accetto il dato. 
			if(y<sqrt(l/2.-pow(x,2.))){
				cont=1;
				vettL[i]=l/2.*(y/sqrt(pow(x,2.)+pow(y,2.)));   //distanza sull'asse delle y tra il baricentro e il vertice del legnetto.
			}
		}
		if((vettPos[i]+vettL[i])>d || (vettPos[i]-vettL[i])<0.){
			Nhit[i]=1;

		}
		else{
			Nhit[i]=0;
		}
	}
	rnd.SaveSeed();

// CALCOLA PI GRECO IN BLOCCHI DA N	
	int L=M/N;
   	double vettPi[N];			
   	double vettPi2[N];
	double Pi_prog[N];
	double Pi2_prog[N];
	double err[N];
	
	for(int i=0;i<N;i++){
		double sum_hit=0;
		for(int j=0;j<L;j++){
			int k=j+i*L;
			if(Nhit[k]==1){
				sum_hit=sum_hit+1;
			}
		}
		vettPi[i]=(2.*l*L)/(sum_hit*d);
		vettPi2[i]=pow(vettPi[i],2.);
	}

	ofstream outData;		
   	outData.open ("DatiPi.dat");

	for(int i=0;i<N;i++){
		for(int j=0;j<i+1;j++){
			Pi_prog[i]+=vettPi[j];
			Pi2_prog[i]+=vettPi2[j];
		}
		Pi_prog[i]=Pi_prog[i]/(i+1);
		Pi2_prog[i]=Pi2_prog[i]/(i+1);
		if (i==0){
			err[i]=0;
		}
		else{
			err[i]=sqrt((Pi2_prog[i]-pow(Pi_prog[i],2.))/i);
		}
		outData<<Pi_prog[i]<<" "<<err[i]<<endl;
	}

	outData.close();
		
		


return 0;
}
