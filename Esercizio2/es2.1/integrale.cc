#include <iostream>
#include "random.h"
#include "integrale.h"
#include <algorithm>
#include <cmath>

using namespace std;

//costruttore
integrale::integrale (double a, double b,funzioneBase *funzione,Random *r){

	_integranda=funzione;
	_a=min(a,b);
	_b=max(a,b);
	//_dim=funzione->getdim();

   // _rand= new Random(3);					//in questo caso la funzione random viene passata dal main, se no si deve creare al suo interno come nel punto nerettato
	_rand=r;
	if (a>b){
		_segno=-1;
		}
	else{
		_segno=1;
		}

}

//distruttore
integrale::~integrale(){}


//metodi

double integrale::midpoint(int nrStep){

	_somma=0.;
	_h=(_b-_a)/nrStep;

	for (int i=0; i<nrStep; ++i){
		double x= _a+(i+0.5)*_h;
		_somma=_somma+_integranda->eval(x);
	}

	_integrale=_segno*_somma*_h;
	return _integrale;

}


double integrale::simpson(int nrStep){

	_h= (_b-_a)/nrStep/2.;
	_somma= 0.;
	double f=0.;
	
	for( int i=0; i<=2*nrStep; i++){
		double x= _a+i*_h;
		if(i==0 || i==2*nrStep){
			f=1./3.*_integranda->eval(x);
		}
		else{
			if(i%2==0){
				f=2./3.*(_integranda->eval(x));
			}
			if(i%2==1){
				f=4./3.*(_integranda->eval(x));
			}
		}
		_somma=_somma+f;
	}
	_integrale=_segno*(_h)*_somma;
	return _integrale;

}


double integrale:: trapezoidi(int nrStep){
	_h=(_b-_a)/nrStep;
	_somma= 0.;
	double f=0.;

	for(int i=0;i<=nrStep;i++){
		double x= _a+i*_h;
		if(i==0 || i==nrStep){
			f=_integranda->eval(x)/2.;
		}
		else {
			f=_integranda->eval(x);
		}
		_somma=_somma+f;
	}

	_integrale=_segno*_h*_somma;
	return _integrale;
}



double integrale:: integraleByMedia(double a, double b, int n){
		double sumfunz=0;
		for(int j=0;j<n;j++){
			double s=_rand->Rannyu();
			double x= a+(b-a)*s;
			double y=_integranda->eval(x);
			sumfunz=sumfunz+y;
		}
		double integrale=0;
		integrale=sumfunz*(b-a)/n;
		return integrale;
	
}

double integrale:: integraleByMedia_distrib(double a, double b, int n, funzioneBase *f){
	double sumfunz=0;
	for(int j=0;j<n;j++){
		double s=_rand->Rannyu();
		double s2=f->eval(s);
		double x=a+(b-a)*s2;
		double y=_integranda->eval(x);
		sumfunz=sumfunz+y;
	}
	double integrale=0;
	integrale=sumfunz*(b-a)/n;
	return integrale;
}


double integrale:: integraleHitorMiss(double a,double b, double fmax, int n){
		double Nint=0;
		double integranda;
		for(int i=0;i<n;i++){
			double s=_rand->Rannyu();
			double t=_rand->Rannyu();
			double x=a+(b-a)*s;
			double y=fmax*t;
			double fx=_integranda->eval(x);
			if(fx>y){Nint++;}
			else continue;
		}
		integranda=(b-a)*fmax*(Nint/n);
		return integranda;
	
}



/*double integrale::integraleMultiDim( double *vettorea, double *vettoreb, int n){
	double sum=0;
	double *array=new double[n];
	int dim=_dim;
	double *vettorex=new double [dim];
	double volume=1;
	for (int i=0; i<n; i++) {
		//volume=1;
		for (int k=0; k<dim; k++) {
			double s=_rand->Rannyu();
			vettorex[k]=vettorea[k]+(vettoreb[k]-vettorea[k])*s;
			volume=volume*(vettoreb[k]-vettorea[k]);
		}
		double y=_integranda->Eval(vettorex);
		sum=sum+y;
		array[i]=y;
		}
	double integrale=sum*volume/n;
	delete [] array;
	
	return integrale;

}*/





