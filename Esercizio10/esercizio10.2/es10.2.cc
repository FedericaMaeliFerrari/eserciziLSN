#include <iostream>
#include <cmath>
#include <fstream>
#include "mpi.h"
#include "random.h"
#include "GA.h"


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

//CIRCONFERENZA
	cout<<"SIMULAZIONE CITTÀ SU CIRCONFERENZA DI RAGGIO 1"<<endl;
	int nstep=1000;
	int nblk=100;
	GA Circ(rnd,0);
	Circ.Coord();
	for(int iblk=1; iblk <= nblk; ++iblk){
		Circ.Reset(iblk);
		for(int istep=1; istep <= nstep; ++istep){
			Circ.Mutazione();
			Circ.NewGeneration();
			Circ.Accumulate();
		}
		Circ.Averages(iblk);
	}
	Circ.BestPath();
	
//QUADRATO
	cout<<"SIMULAZIONE CITTÀ SU QUADRATO DI LATO 2"<<endl;
	GA Quad(rnd,1);
	Quad.Coord();
	for(int iblk=1; iblk <= nblk; ++iblk){
		Quad.Reset(iblk);
		for(int istep=1; istep <= nstep; ++istep){
			Quad.Mutazione();
			Quad.NewGeneration();
			Quad.Accumulate();
		}
		Quad.Averages(iblk);
	}
	Quad.BestPath();
	

return 0;
}
