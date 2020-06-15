#ifndef _Metropolis_h_
#define _Metropolis_h_
#include "random.h"
#include "randomDistrib.h"
#include "funzioneBase.h"


class Metropolis {

public:
	Metropolis(Random *rand, funzioneBase *funz);
	Metropolis(Random *rand, funzioneBase *funz, funzioneBase *funz2);
	void MetrUnif(double x,double y,double z,double delta);
	void MetrGauss(double x,double y,double z,double delta);
	void MetrUnif_ecc(double x,double y,double z,double delta);
	void MetrGauss_ecc(double x,double y,double z,double delta);
	double getX(){return _x;}
	double getY(){return _y;}
	double getZ(){return _z;}
	double getcont(){return _cont;}
	void setcont(){_cont=0;}

protected:
	Random* _rand;
	funzioneBase* _funz;
	funzioneBase* _funz2;
	int _cont;
	
	double _x, _y, _z;






};

#endif
