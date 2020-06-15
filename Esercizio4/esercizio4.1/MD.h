#ifndef _MD_h_
#define _MD_h_
#include "random.h"

class MD {

public:
	//Costruttore e inizializzazione valori
	MD(const char*filename,Random *r);
	void inputPos(const char*filename);
	void inputV();
	int getnstep(){return _nstep;}
	int getiprint(){return _iprint;}

	//Condizioni al contorno periodiche
	double pbc(double r);
	//Potenziale di Lennard-Jones 
	double LennJon(int i,int asse);

	//Moto delle molecole
	void Move();

	//Misura dei paramteri
	void Misura(int k);
	
	//Scrittura su files
	void ConfFinal();
	void ConfXYZ(int nconf);
	void OldConf();

	//Funzioni per l'opzione re-start
	void inputPosRestart();
	

protected:
	double _t,_rho,_rcut,_delta;
	int _npart,_nstep,_iprint;
	double _vol,_box;

	double* _posx;
	double* _posy;
	double* _posz;

	double* _vx;
	double* _vy;
	double* _vz;

	double* _xold;
	double* _yold;
	double* _zold;

	Random* _rand;


};
#endif
