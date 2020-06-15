#ifndef _Quantum_h_
#define _Quantum_h
#include "random.h"

class Quantum {

public:	
	Quantum(Random* rand,double x0,double mu,double sigma,double delta);
	double funz_prob(double x);
	void Metropolis(int k);
	double Vpot(double x);
	void Measure();
	void Reset(int iblk);
	void Media(int iblk,int k);
	double Error(double sum, double sum2, int iblk);

	double GetGlobH(){return _globH;}
	void SetParametres(double mu, double sigma);


protected:
	Random* _rand;
	double _mu,_sigma,_delta;
	double _x;
	int _blknorm;
	double _globH,_globH2,_blkH,accepted,attempted;



};
#endif
