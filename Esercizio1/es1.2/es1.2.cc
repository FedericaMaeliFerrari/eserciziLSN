#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "randomDistrib.h"


using namespace std;

int main (int argc, char *argv[]){

	int n;
	int numMedie;
	cout<<"Il programma svolge gli istogrammi dei valori medi dei nuemeri generati in modo casuale"<<endl<<" secondo tre distribuzioni: normale, esponenziale, lorentziana."<<endl<<"per tutte e tre le tipologie gli istogrammi verranno calcolati al variare del numero di dati"<<endl<<"utilizzati per la media."<<endl;
	cout<<"Inserire quante medie calcolare"<<endl;
	cin>>numMedie;
	cout<<"Inserire quanti istogrammi calcolare per ogni distribuzione"<<endl;
	cin>>n;
	double vettMedie[numMedie][n];
	double N[n];
	cout<<"Inserire il numero di valori N con cui fare la media andando ogni volta a capo"<<endl;
	for (int i=0;i<n;i++){
		cin>>N[i];
	}
	/*N[0]=1;
	N[1]=2;
	N[2]=10;
	N[3]=100;*/

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


//GENERAZIONE NUMERI CASUALI TRA [0,1]

	for(int i=0;i<n;i++){
		for(int j=0;j<numMedie;j++){
			double sum=0.;
			for(int k=0;k<N[i];k++){
				double xi=rnd.Rannyu();
				sum=sum+xi;
			}
			vettMedie[j][i]=sum/N[i];
		}
	}


	ofstream outData;
	outData.open("Dati[0,1].dat");
	outData<<" ";
	for(int i=0;i<numMedie;i++){
		for(int j=0;j<n;j++){
			outData<<setprecision(6);
			outData<<vettMedie[i][j]<<" ";
			
			}
		outData<<endl;
		}
		
	outData.close();


//GENERAZIONE NUMERI SECONDO DISTRIBUZIONE ESPONENZIALE

	RandomDistrib rde;
	double l;
	cout<<"Inserire il parametro lambda per la distribuzione esponenziale exp(-lambda*x)"<<endl;
	cin>>l;
	for(int i=0;i<n;i++){
		for(int j=0;j<numMedie;j++){
			double sum=0.;
			for(int k=0;k<N[i];k++){
				double xi=rnd.Rannyu();
				double x=rde.exp_trasformata(l,xi);
				sum=sum+x;
			}
			vettMedie[j][i]=sum/N[i];
		}
	}

	ofstream outDataExp;
	outDataExp.open("DatiExp.dat");
	outDataExp<<" ";
	for(int i=0;i<numMedie;i++){
		for(int j=0;j<n;j++){
			outDataExp<<setprecision(6);
			outDataExp<<vettMedie[i][j]<<" ";
		
			}
		outDataExp<<endl;
		}
	
	outDataExp.close();

	

//GENERAZIONE NUMERI SECONDO DISTRIBUZIONE LORENTZIANA

	RandomDistrib rdl;
	double x0;
	double g;
	cout<<"Inserire i paramentri x0 e gamma per la distribuzione lorentziana "<<endl;
	cin>>x0>>g;

	for(int i=0;i<n;i++){
		for(int j=0;j<numMedie;j++){
			double sum=0.;
			for(int k=0;k<N[i];k++){
				double xi=rnd.Rannyu();
				double x=rdl.lorentz_trasformata(x0,g,xi);
				sum=sum+x;
			}
			vettMedie[j][i]=sum/N[i];
		}
	}
	
	ofstream outDataL;
	outDataL.open("DatiLorentz.dat");
	outDataL<<" ";
	for(int i=0;i<numMedie;i++){
		for(int j=0;j<n;j++){
			outDataL<<setprecision(6);
			outDataL<<vettMedie[i][j]<<" ";
		
			}
		outDataL<<endl;
		}
	
	outDataL.close();




return 0;

}
