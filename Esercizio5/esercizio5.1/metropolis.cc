#include "metropolis.h"
#include "funzioneBase.h"
#include "random.h"
#include "randomDistrib.h"
#include "iostream"
#include <cmath>

using namespace std;

Metropolis::Metropolis(Random *rand, funzioneBase *funz){
	_rand=rand;
	_funz=funz;
	_cont=0;
}




Metropolis::Metropolis(Random *rand, funzioneBase *funz, funzioneBase *funz2){
	_rand=rand;
	_funz=funz;
	_funz2=funz2;
	_cont=0;
}




void Metropolis::MetrUnif(double x,double y,double z,double delta){
	double x0=x;
	double y0=y;
	double z0=z;
	double r0=sqrt(pow(x0,2.)+pow(y0,2.)+pow(z0,2.));
	
	double xnew=x0+_rand->Rannyu(-delta,delta);
	double ynew=y0+_rand->Rannyu(-delta,delta);
	double znew=z0+_rand->Rannyu(-delta,delta);
	
	double rnew=sqrt(pow(xnew,2.)+pow(ynew,2.)+pow(znew,2.));

	//Valuto la probabilità di accettazione
	double px=_funz->eval(r0);
	double px2=_funz->eval(rnew);
	double a=min(1.,px2/px);
		
	
	//Accettare o Rifiutare?
	double val=_rand->Rannyu();
	if(val<=a){
		_x=xnew;
		_y=ynew;
		_z=znew;
		_cont+=1;
	}
	else{
		_x=x0;
		_y=y0;
		_z=z0;
	}


	return;
}




void Metropolis::MetrGauss(double x,double y,double z,double delta){
	double x0=x;
	double y0=y;
	double z0=z;
	double r0=sqrt(pow(x0,2.)+pow(y0,2.)+pow(z0,2.));
	
	double xnew=_rand->Gauss(x0,pow(delta,2.));
	double ynew=_rand->Gauss(y0,pow(delta,2.));
	double znew=_rand->Gauss(z0,pow(delta,2.));
	double rnew=sqrt(pow(xnew,2.)+pow(ynew,2.)+pow(znew,2.));

	//Valuto la probabilità di accettazione

	double px=_funz->eval(r0);
	double px2=_funz->eval(rnew);
	double a=min(1.,px2/px);
			
	//Accettare o Rifiutare?
	double val=_rand->Rannyu();
	if(val<=a){
		_x=xnew;
		_y=ynew;
		_z=znew;
	}
	else{
		_x=x0;
		_y=y0;
		_z=z0;
	}


	return;
}





void Metropolis::MetrUnif_ecc(double x,double y,double z,double delta){
	double x0=x;
	double y0=y;
	double z0=z;
	double r0=sqrt(pow(x0,2.)+pow(y0,2.)+pow(z0,2.));
	
	double xnew=x0+_rand->Rannyu(-delta,delta);
	double ynew=y0+_rand->Rannyu(-delta,delta);
	double znew=z0+_rand->Rannyu(-delta,delta);
	double rnew=sqrt(pow(xnew,2.)+pow(ynew,2.)+pow(znew,2.));

	//Valuto la probabilità di accettazione
	double px1=_funz->eval(rnew);
	double px2=_funz2->eval(znew);
	double pxtot=px1*px2;
	double px1old=_funz->eval(r0);
	double px2old=_funz2->eval(z0);
	double pxold=px1old*px2old;
	double a=min(1.,pxtot/pxold);
	
	//Accettare o Rifiutare?
	double val=_rand->Rannyu();
	if(val<=a){
		_x=xnew;
		_y=ynew;
		_z=znew;
		_cont+=1;
	}
	else{
		_x=x0;
		_y=y0;
		_z=z0;
	}


	return;
}




void Metropolis::MetrGauss_ecc(double x,double y,double z,double delta){
	double x0=x;
	double y0=y;
	double z0=z;
	double r0=sqrt(pow(x0,2.)+pow(y0,2.)+pow(z0,2.));
	
	double xnew=_rand->Gauss(x0,pow(delta,2.));
	double ynew=_rand->Gauss(y0,pow(delta,2.));
	double znew=_rand->Gauss(z0,pow(delta,2.));
	double rnew=sqrt(pow(xnew,2.)+pow(ynew,2.)+pow(znew,2.));

	//Valuto la probabilità di accettazione
	double px1=_funz->eval(rnew);
	double px2=_funz2->eval(znew);
	double pxtot=px1*px2;
	double px1old=_funz->eval(r0);
	double px2old=_funz2->eval(z0);
	double pxold=px1old*px2old;
	double a=min(1.,pxtot/pxold);
			
	//Accettare o Rifiutare?
	double val=_rand->Rannyu();
	if(val<=a){
		_x=xnew;
		_y=ynew;
		_z=znew;
	}
	else{
		_x=x0;
		_y=y0;
		_z=z0;
	}


	return;
}

