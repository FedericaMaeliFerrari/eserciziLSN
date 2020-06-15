#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "randomDistrib.h"

using namespace std;

int main(int argc, char *argv[]){
	int M;
	int N;
	cout<<"Inserire quante simulazioni fare e quanti passi fare"<<endl;
	cin>>M>>N;
	double rWalk[N];
	double err[N];
	double rWalk_s[N];
	double err_s[N];
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

//RANDOM WALK DISCRETO SU UN CUBO + RANDOM WALK CONTINUO SU UNA SFERA
	RandomDistrib rds;
	ofstream outData;
	outData.open("RandomWalk.dat");
	for(int i=0;i<N;i++){
		double percorso=0;
		double percorso2=0;
		double p_sfera=0;
		double p2_sfera=0;
		for(int j=0;j<M;j++){
			double x=0;
			double y=0;
			double z=0;
			double x_sfera=0;
			double y_sfera=0;
			double z_sfera=0;
			for(int k=0;k<i+1;k++){
				//random walk discreto
				int a=0;
				double passo=rnd.Rannyu(-1.,1.);
				if (passo>=0.){a =1;}
				if (passo<0.){a =-1;}
				int dir=rnd.Rannyu(1.,4.);
				if(dir==1){x=x+a;}
				if(dir==2){y=y+a;}
				if(dir==3){z=z+a;}
				//random walk continuo
				double y_t=rnd.Rannyu(0.,2.);
				double y_f=rnd.Rannyu(0.,2.);
				double theta=rds.sen_trasformata(y_t);
				double phi=rds.sen_trasformata(y_f);
				if (a==-1){phi=phi+M_PI;}
				x_sfera=x_sfera+cos(phi)*sin(theta);
				y_sfera=y_sfera+sin(phi)*sin(theta);
				z_sfera=cos(theta);
			}	
			percorso=percorso+pow(x,2.)+pow(y,2.)+pow(z,2.);
			percorso2=percorso2+pow(pow(x,2.)+pow(y,2.)+pow(z,2.),2.);
			p_sfera=p_sfera+pow(x_sfera,2.)+pow(y_sfera,2.)+pow(z_sfera,2.);
			p2_sfera=p2_sfera+pow(pow(x_sfera,2.)+pow(y_sfera,2.)+pow(z_sfera,2.),2.);
		}
		rWalk[i]=pow(percorso/M,0.5);
		rWalk_s[i]=pow(p_sfera/M,0.5);
		err[i]=(sqrt((percorso2/M-pow(percorso/M,2.))/M));
		err_s[i]=(sqrt((p2_sfera/M-pow(p_sfera/M,2.))/M));
		outData<<rWalk[i]<<" "<<err[i]<<" "<<rWalk_s[i]<<" "<<err_s[i]<<endl;
	}
	

			


return 0;
}
