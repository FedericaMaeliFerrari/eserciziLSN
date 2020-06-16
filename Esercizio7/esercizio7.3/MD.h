#ifndef _MD_h_
#define _MD_h_
#include "random.h"
//#include "vettore.h"

class MD {

public:
	//Costruttore e inizializzazione valori
	MD(const char*filename,Random *r);
	void inputPos(const char*filename);
	void inputV();
	void Accumulate();
	void Reset(int);
	void Averages(int);

	//funzioni che restituscono variabili protected
	int getnstep(){return _nstep;}
	int getiprint(){return _nblk;}
	//double getT(){return _temp;}
	//double getEkin(){return _ekin;}
	//double getEpot(){return _epot;}
	//double getEtot(){return _etot;}
	
	//Condizioni al contorno periodiche
	double pbc(double r);
	//Potenziale di Lennard-Jones 
	double LennJon(int i,int asse);

	//Moto delle molecole
	void Move();

	//Misura dei paramteri
	void Misura();
	void Misura_blocchi();
	double Error(double,double,int);
	
	//Scrittura su files
	void ConfFinal();
	void ConfXYZ(int nconf);
	void OldConf();
	void Gavefinal();

	//Funzioni per l'opzione re-start
	void inputPosRestart();
	

protected:
	double _t,_rho,_rcut,_delta;
	int _npart,_nstep,_nblk;
	double _vol,_box;
	double _vtail,_ptail;
	int _igofr,_nbins,_iv,_iw;
	double bin_size;


	double* _posx;
	double* _posy;
	double* _posz;

	double* _vx;
	double* _vy;
	double* _vz;

	double* _xold;
	double* _yold;
	double* _zold;

	//double _temp;
	//double _ekin;
	//double _epot;
	//double _etot;
	int _it,_iek,_iep,_iet,n_props;

	double* _walker;

	double* glob_av;
	double* glob_av2;
	double* blk_av;
	int blk_norm;


	Random* _rand;


};
#endif
