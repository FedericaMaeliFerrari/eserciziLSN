#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "MD.h"
//#include "vettore.h"

using namespace std;

int main(int argc, char **argv){
	if(argc!=3){
	cerr<<"per l'esecuzione del programma" <<argv[0]<< " bisogna inserire: 1)il nome del programma, 2)il file da cui leggere i valori (rho, num.particelle ecc.) e 3)il file da cui caricare le posizioni iniziali delle particelle"<<endl;
	}
	//int n=100;
	//int l=100;

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
/*
//INIZIALIZZO I VETTORI PER LA MEDIA A BLOCCHI
	Vettore T(n);
	Vettore Ekin(n);
	Vettore Epot(n);
	Vettore Etot(n);
*/
	
//INIZIALIZZO I VALORI PER LA DINAMICA MOLECOLARE;
	MD molecole(argv[1],rnd);
	molecole.inputPos(argv[2]);
	molecole.inputV();

//MOTO DELLE MOLECOLE SECONDO ALGORITMO DI VELVET
	int nstep=molecole.getnstep();
	int iprint=molecole.getiprint();
	//int nconf = 1;

	for(int i=1;i<=nstep;i++){
		molecole.Move();
		if(i%iprint==0){cout<<"Numero di Time-Step: "<<i<<endl;}
		if(i%10==0){
			molecole.Misura(0);
			//nconf += 1;
		}
	}

	molecole.OldConf();
	molecole.ConfFinal();


	
//MISURA CON OPZIONE RESTART
	int nblk=100;
	nstep=1000;
	molecole.inputPosRestart();
	for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
    		molecole.Reset(iblk);   //Reset block averages
    		for(int istep=1; istep <= nstep; ++istep){
      			molecole.Move();
      			molecole.Misura_blocchi();
      			molecole.Accumulate(); //Update block averages
    		}
    		molecole.Averages(iblk);   //Print results for current block
  	}
  	molecole.ConfFinal(); //Write final configuration









/*
	double t,ekin,epot,etot;
	int nconf=0.;
	for(int k=0;k<n;k++){
		t=0.;
		ekin=0.;
		epot=0.;
		etot=0.;
		cout<<"Blocco Simulazione: "<<k<<endl;
		molecole.inputPosRestart();
		for(int j=0;j<l;j++){
			molecole.Move();
			if(j%10==0){
				molecole.Misura_blocchi();
				molecole.ConfXYZ(nconf);
				nconf++;
				t+=molecole.getT();
				ekin+=molecole.getEkin();
				epot+=molecole.getEpot();
				etot+=molecole.getEtot();
				
			}
		}
		cout<<t<<" "<<ekin<<" "<<epot<<" "<<etot<<endl;
		t=t/l;
		ekin=ekin/l;
		epot=epot/l;
		etot=etot/l;
		T.SetComponent(k,t);
		Ekin.SetComponent(k,ekin);
		Epot.SetComponent(k,epot);
		Etot.SetComponent(k,etot);
		molecole.OldConf();
			
	}

	cout<<"HO FINITO LE SIMULAZIONI SUI 100 BLOCCHI"<<endl;
	cout<<T.GetComponent(0)<<endl;
	cout<<T.GetComponent(99)<<endl;
//MEDIA A BLOCCHI E RELATIVO ERRORE
	T.mediaByBlocchi();
	Ekin.mediaByBlocchi();
	Epot.mediaByBlocchi();
	Etot.mediaByBlocchi();
	cout<<"CALCOLATA MEDIA A BLOCCHI"<<endl;

//STAMPO FILE CON MEDIE E ERRORI
	ofstream Output;
	Output.open("Simulazione.dat");
	for(int i=0;i<n;i++){
		Output<<T.GetComp_media(i)<<" "<<T.GetComp_err(i)<<" "<<Ekin.GetComp_media(i)<<" "<<Ekin.GetComp_err(i)<<" "<<Epot.GetComp_media(i)<<" "<<Epot.GetComp_err(i)<<" "<<Etot.GetComp_media(i)<<" "<<Etot.GetComp_err(i)<<endl;
	}
	Output.close();


*/

	return 0;
}		
	


