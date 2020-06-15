#ifndef _INTEGRALE_H
#define _INTEGRALE_H
#include "random.h"
#include "funzioneBase.h"
//#include "vettore.h"


class integrale {

public:
	//costruttore
	integrale (double a, double b, funzioneBase *funzione,Random *r);
	//distruttore
	~integrale();
	
	//metodi
	double midpoint (int nrStep);
	double simpson (int nrStep);
	double trapezoidi(int nrStep);
	double integraleByMedia (double a, double b, int n);
	double integraleHitorMiss(double a,double b,  double fmax, int n);
	//double integraleMultiDim( double *vettorea, double *vettoreb, int n);
	double integraleByMedia_distrib(double a, double b, int n, funzioneBase *f);



private:
	double _a;
	double _b;
	int _segno;
	double _h;
	double _somma;
	double _integrale;
	//int _dim;
	funzioneBase *_integranda; 
	Random* _rand;
};

#endif 
