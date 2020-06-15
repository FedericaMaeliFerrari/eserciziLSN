#include <iostream>
#include <fstream>
#include <cmath>
#include "integrale.h"
#include "funzioneBase.h"
#include "funzione.h"
#include "random.h"
#include "finanza.h"
#include "randomDistrib.h"

using namespace std;

int main(int argc, char *argv[]){

//VALORI AL TEMPO T=0
	double t=0.;
	double S0=100.;
	double T=1.;
	double K=100.;
	double r=0.1;
	double sigma=0.25;

	int M;
	int N;
	cout<<"Quante simulazioni condurre (M)? E in quanti passi (N,usati per media a blocchi)?"<<endl;
	cin>>M>>N;
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


//CALCOLO IN MOTO ANALITICO CALL E PUT
	funzioneBase *funz=new esponenziale(1.,-1.,2.,0.);  //funzione esponenziale che andrà integrata per trovare N(d1) e N(d2)
	integrale *I=new integrale(0.,1.,funz,rnd);     
	Finanza f(t,S0,T,K,r,sigma,I);                  //calcola C e P a seconda dei parametri iniziali e tramite l'integrale (metodo by media).				
						
	double c=f.getc();
	double p=f.getp();

	ofstream out;
	out.open("CP_attesi.dat");
	out<<c<<" "<<p<<endl;
	out.close();

//CALCOLO CALL E PUT IN MODO DIRETTO
	double CP_direct[N][2];
	double CP2_direct[N][2];
	RandomDistrib *g=new RandomDistrib();
	for(int i=0;i<N;i++){
		CP_direct[i][0]=f.direct_sampleC(L,0.,g,rnd);
		CP_direct[i][1]=f.direct_sampleP(L,0.,g,rnd);
		CP2_direct[i][0]=pow(CP_direct[i][0],2.);
		CP2_direct[i][1]=pow(CP_direct[i][1],2.);
	}

//CALCOLO CALL E PUT IN MODO DISCRETO
	double CP_discret[N][2];
	double CP2_discret[N][2];
	for(int i=0;i<N;i++){
		CP_discret[i][0]=f.discret_sampleC(100.,L,g,rnd);
		CP_discret[i][1]=f.discret_sampleP(100.,L,g,rnd);
		CP2_discret[i][0]=pow(CP_discret[i][0],2.);
		CP2_discret[i][1]=pow(CP_discret[i][1],2.);
	}


//MEDIA A BLOCCHI E SALVATAGGIO SU UN FILE
	ofstream outData;
	ofstream outData2;
	outData2.open("Call-Put_discreto.dat");
	outData.open ("Call-Put_direct.dat");
	double CP_d[N][2]={0};
	double CPerr_d[N][2]={0};
	double CP_disc[N][2]={0};
	double CPerr_disc[N][2]={0};
	double err_dir[N][2]={0};
	double err_disc[N][2]={0};
	for(int i=0;i<N;i++){
		for(int j=0;j<i+1;j++){
			CP_d[i][0]+=CP_direct[j][0];
			CP_d[i][1]+=CP_direct[j][1];
			CPerr_d[i][0]+=CP2_direct[j][0];
			CPerr_d[i][1]+=CP2_direct[j][1];
			CP_disc[i][0]+=CP_discret[j][0];
			CP_disc[i][1]+=CP_discret[j][1];
			CPerr_disc[i][0]+=CP2_discret[j][0];
			CPerr_disc[i][1]+=CP2_discret[j][1];
		}
		CP_d[i][0]=CP_d[i][0]/(i+1);		
		CP_d[i][1]=CP_d[i][1]/(i+1);
		CPerr_d[i][0]=CPerr_d[i][0]/(i+1);
		CPerr_d[i][1]=CPerr_d[i][1]/(i+1);
		CP_disc[i][0]=CP_disc[i][0]/(i+1);
		CP_disc[i][1]=CP_disc[i][1]/(i+1);
		CPerr_disc[i][0]=CPerr_disc[i][0]/(i+1);
		CPerr_disc[i][1]=CPerr_disc[i][1]/(i+1);
		if (i==0){
			err_dir[i][0]=0;
			err_dir[i][1]=0;
			err_disc[i][0]=0;
			err_disc[i][1]=0;
		}
		else{
			err_dir[i][0]=sqrt((CPerr_d[i][0]-pow(CP_d[i][0],2.))/i);
			err_dir[i][1]=sqrt((CPerr_d[i][1]-pow(CP_d[i][1],2.))/i);
			err_disc[i][0]=sqrt((CPerr_disc[i][0]-pow(CP_disc[i][0],2.))/i);
			err_disc[i][1]=sqrt((CPerr_disc[i][1]-pow(CP_disc[i][1],2.))/i);

		}
		outData<<CP_d[i][0]<<" "<<err_dir[i][0]<<" "<<CP_d[i][1]<<" "<<err_dir[i][1]<<endl;
		outData2<<CP_disc[i][0]<<" "<<err_disc[i][0]<<" "<<CP_disc[i][1]<<" "<<err_disc[i][1]<<endl;
	}
	outData.close();
	outData2.close();

return 0;
}
