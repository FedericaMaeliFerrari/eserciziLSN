#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "randomDistrib.h"
#include "funzione.h"
#include "vettore.h"
#include "funzioneBase.h"
#include "metropolis.h"

using namespace std;

int main(int argc, char *argv[]){
	cout<<"Calcolo della posizione nello STATO FONDAMENTALE dell'idrogeno utilizzando una matrice di Transizione uniforme (coordinate sferiche)"<<endl;
	int n=100;
	int l=1000;
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


//INIZIALIZZO VETTORI E FUNZIONI CHE MI SERVIRANNO PER IL CALCOLO
	Vettore vettr(n);
	Vettore vettrgauss(n);
	double xcoord[n*l];      //vettori coordinate punti generati con matrice transizione uniforme
	double ycoord[n*l];
	double zcoord[n*l];
	double xcoord_g[n*l];	 //vettori coordinate punti generati con matrice transizione gaussiana
	double ycoord_g[n*l];
	double zcoord_g[n*l];
	funzioneBase *prob_fond= new esponenziale(1./M_PI,-2.,1.,0.);         //|Psi(r,t)|^2 dello stato fondamentale
	funzioneBase *prob_ecc1= new esponenziale(1./(32.*M_PI),-1.,1.,0.);   //|Psi(r,t)|^2 dello stato eccitato: ho spezzato la funzione in due funzioni, una dipendente da r (prob_ecc1)
	funzioneBase *prob_ecc2= new parabola(1.,0.,0.);                           // e l'altra dipendente da theta (prob_ecc2)				
	Metropolis M(rnd,prob_fond);
	Metropolis M2(rnd,prob_ecc1,prob_ecc2);


//ALGORITMO METROPOLIS, MATRICE TRANSIZIONE UNIFORME E GAUSSIANA, STATO FONDAMENTALE
	double deltaunif=1.25;
	double deltagauss=0.75;
	double xunif,yunif,zunif;
	double xgauss,ygauss,zgauss;
	double sumrunif;
	double sumrgauss;
	for(int i=0;i<n;i++){
		xunif=0.;       //Inizializzazione coordinate per il calcolo
		yunif=0.;
		zunif=0.;
		sumrunif=0.;
		xgauss=1.;
		ygauss=1.;
		zgauss=1.;
		sumrgauss=0.;
		M.setcont();
		for(int j=0;j<l;j++){
			//stato fondamentale,uniforme
			M.MetrUnif(xunif,yunif,zunif,deltaunif);
			xunif=M.getX();
			xcoord[i*l+j]=xunif;       //Riempio il vettore dell coordinate (che mi servirà per lo scatterplot)
			yunif=M.getY();
			ycoord[i*l+j]=yunif;
			zunif=M.getZ();
			zcoord[i*l+j]=zunif;
			sumrunif+=sqrt(pow(xunif,2.)+pow(yunif,2.)+pow(zunif,2.));

			//stato fondamentale,gauss
			M.MetrGauss(xgauss,ygauss,zgauss,deltagauss);
			xgauss=M.getX();
			xcoord_g[i*l+j]=xgauss;
			ygauss=M.getY();
			ycoord_g[i*l+j]=ygauss;
			zgauss=M.getZ();
			zcoord_g[i*l+j]=zgauss;
			sumrgauss+=sqrt(pow(xgauss,2.)+pow(ygauss,2.)+pow(zgauss,2.));
		}
		cout<<M.getcont()/l<<endl;

		vettr.SetComponent(i,sumrunif/l);

		vettrgauss.SetComponent(i,sumrgauss/l);

	}
	
//MEDIA A BLOCCHI E SCRIVO SU FILE
	vettr.mediaByBlocchi();
	vettrgauss.mediaByBlocchi();
	ofstream Output;
	ofstream Out;
	Out.open("Fond_Unif_coord.dat");
	Output.open("Fond_Unif.dat");
	for(int i=0;i<n*l;i++){
		Out<<xcoord[i]<<" "<<ycoord[i]<<" "<<zcoord[i]<<endl;
	}
	Out.close();
	for(int i=0;i<n;i++){
		Output<<vettr.GetComp_media(i)<<" "<<vettr.GetComp_err(i)<<endl;
	}
	Output.close();


	ofstream OutputGauss;
	ofstream OutGauss;
	OutGauss.open("Fond_Gauss_coord.dat");
	OutputGauss.open("Fond_Gauss.dat");
	for(int i=0;i<n*l;i++){
		OutGauss<<xcoord_g[i]<<" "<<ycoord_g[i]<<" "<<zcoord_g[i]<<endl;
	}
	OutGauss.close();
	for(int i=0;i<n;i++){
		OutputGauss<<vettrgauss.GetComp_media(i)<<" "<<vettrgauss.GetComp_err(i)<<endl;
	}
	OutputGauss.close();
	


//ALGORITMO METROPOLIS, MATRICE TRANSIZIONE UNIFORME E GAUSSIANA, STATO ECCITATO
	deltaunif=3.0;
	deltagauss=1.85;
	for(int i=0;i<n;i++){
		xunif=0.;
		yunif=0.;
		zunif=0.;
		sumrunif=0.;
		xgauss=1.;
		ygauss=1.;
		zgauss=1.;
		sumrgauss=0.;
		M2.setcont();
		for(int j=0;j<l;j++){
			//stato eccitato,uniforme
			M2.MetrUnif_ecc(xunif,yunif,zunif,deltaunif);
			xunif=M2.getX();
			xcoord[i*l+j]=xunif;
			yunif=M2.getY();
			ycoord[i*l+j]=yunif;
			zunif=M2.getZ();
			zcoord[i*l+j]=zunif;
			sumrunif+=sqrt(pow(xunif,2.)+pow(yunif,2.)+pow(zunif,2.));

			//stato eccitato,gauss
			M2.MetrGauss_ecc(xgauss,ygauss,zgauss,deltagauss);
			xgauss=M2.getX();
			xcoord_g[i*l+j]=xgauss;
			ygauss=M2.getY();
			ycoord_g[i*l+j]=ygauss;
			zgauss=M2.getZ();
			zcoord_g[i*l+j]=zgauss;
			sumrgauss+=sqrt(pow(xgauss,2.)+pow(ygauss,2.)+pow(zgauss,2.));
		}

		vettr.SetComponent(i,sumrunif/l);

		vettrgauss.SetComponent(i,sumrgauss/l);

	}


//MEDIA A BLOCCHI E SCRIVO SU FILE
	vettr.mediaByBlocchi();
	vettrgauss.mediaByBlocchi();
	ofstream OutputEcc;
	ofstream OutEcc;
	OutEcc.open("Ecc_Unif_coord.dat");
	OutputEcc.open("Ecc_Unif.dat");
	for(int i=0;i<n*l;i++){
		OutEcc<<xcoord[i]<<" "<<ycoord[i]<<" "<<zcoord[i]<<endl;
	}
	OutEcc.close();
	for(int i=0;i<n;i++){
		OutputEcc<<vettr.GetComp_media(i)<<" "<<vettr.GetComp_err(i)<<endl;
	}
	OutputEcc.close();


	ofstream OutputGaussEcc;
	ofstream OutGaussEcc;
	OutGaussEcc.open("Ecc_Gauss_coord.dat");
	OutputGaussEcc.open("Ecc_Gauss.dat");
	for(int i=0;i<n*l;i++){
		OutGaussEcc<<xcoord_g[i]<<" "<<ycoord_g[i]<<" "<<zcoord_g[i]<<endl;
	}
	OutGaussEcc.close();
	for(int i=0;i<n;i++){
		OutputGaussEcc<<vettrgauss.GetComp_media(i)<<" "<<vettrgauss.GetComp_err(i)<<endl;
	}
	OutputGaussEcc.close();
		




return 0;
}
