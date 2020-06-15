#ifndef _finanza_h
#define _finanza_h
#include "integrale.h"
#include "randomDistrib.h"

class Finanza {

public:
	Finanza(double t,double S0,double T,double K,double r,double sigma,integrale* i);
	double getd1(){return _d1;}
	double getd2(){return _d2;}
	double getc(){return _c;}
	double getp(){return _p;}
	double direct_sampleC(int L,double mu, RandomDistrib* g,Random* r);
	double direct_sampleP(int L,double mu, RandomDistrib* g,Random* r);
	double discret_sampleC(double passi,int L, RandomDistrib* g,Random* r);
	double discret_sampleP(double passi,int L, RandomDistrib* g,Random* r);


protected:
	double _t,_S0,_T,_K,_r,_sigma;
	double _d1,_d2;
	double _c,_p;
	integrale* _i;

};
#endif
