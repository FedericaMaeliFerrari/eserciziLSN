#ifndef _NVT_h_
#define _NVT_h_
#include "random.h"

class NVT{

public:
	NVT(const char*filename,Random* r,int k);	
	void Move();
	void Measure(int k);

	void Reset(int);
	void Accumulate();
	void Averages(int);

	double Boltzmann(double, double, double, int);
	double Pbc(double);
	double Error(double,double,int);

	void ConfFinal();
	void ConfXYZ(int);

	int GetNblk(){return _nblk;}
	int GetNstep(){return _nstep;}

protected:

	Random* _rnd;

	int _npart,_nstep,_nblk;
	double _temp,_rcut,_rho,_beta,_vol,_box,_delta;
	double _vtail,_ptail;
	int _iv,_iw,n_props,blk_norm;
	int _igofr,_nbins;
	double bin_size;
	double attempted,accepted;

	double* _x;
	double* _y;
	double* _z;
	double* _walker;

	double _xold,_yold,_zold;
	double _xnew,_ynew,_znew;


	double* glob_av;
	double* glob_av2;
	double* blk_av;
	


	



};
#endif

