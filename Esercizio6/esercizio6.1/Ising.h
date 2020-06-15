#ifndef _Ising_h_
#define _Ising_h_
#include "random.h"

class Ising {

public:
	Ising(Random *r);
	void SetSpin(int opz);
	void Reset(int);
	void Accumulate();
	void Averages(int,int);
	void Move(int);
	void ConfFinal();
	void Measure();
	double Boltzmann(int, int);
	double BoltzmannM(int, int);
	int Pbc(int);
	double Error(double,double,int);

	int GetNblk(){return _nblk;}
	int GetNstep(){return _nstep;}
	int GetMetro(){return _metro;}
	double GetTemp(){return _temp;}
	void SetTemp(double temp);
	void SetMetro(int metro){_metro=metro;}

protected:
	Random* _rnd;

	int _nspin,_nstep,_nblk,_metro;
	double _beta,_temp,_J,_h;
	int _iu,_im,_ix,_ic;
	int n_props;
	double _kb;

	int* _s;
	int* _sm;
	double* _walker;

	double* glob_av;
	double* glob_av2;
	double* blk_av;
	int blk_norm,attempted,accepted;

};

#endif
