#include "Quantum.h"
#include "random.h"
#include "iostream"
#include <cmath>
#include <fstream>

using namespace std;

Quantum::Quantum(Random* rand,double x0,double mu,double sigma,double delta){
	_rand=rand;
	_mu=mu;
	_sigma=sigma;
	_delta=delta;
	_x=x0;
}


double Quantum::funz_prob(double x){
	return exp(-pow((x-_mu),2.)/pow(_sigma,2.))+2*exp(-(pow(x,2.)+pow(_mu,2.))/pow(_sigma,2.))+exp(-pow((x+_mu),2.)/pow(_sigma,2.));
}


double Quantum::Vpot(double x){
	return pow(x,4.)-2.5*pow(x,2.);
}



void Quantum::Metropolis(int k){

	double x0=_x;
	double xnew=x0+_rand->Rannyu(-_delta,_delta);

	double pxold=funz_prob(x0);
	double pxnew=funz_prob(xnew);
	double a=min(1.,pxnew/pxold);
	double val=_rand->Rannyu();
	if(val<=a){
		_x=xnew;
		accepted=accepted+1;
	}
	else{_x=x0;}
	attempted=attempted+1;
	if(k==1){
		ofstream Pos;
		Pos.open("posX.dat",ios::app);
		Pos<<_x<<endl;
		Pos.close();
	}

	return;
}


void Quantum::Measure(){
	double x=_x;
	_blkH+=Vpot(x);
	_blknorm+=1;
	return;
}


void Quantum::Reset(int iblk){
	if(iblk==1){
		_globH=0;
		_globH2=0;
	}
	_blkH=0;
	_blknorm=0;
	attempted = 0;
   	accepted = 0;
	return;
}




void Quantum::Media(int iblk,int k){
	if(k==1){	
		cout << "Block number " << iblk << endl;
   		cout << "Acceptance rate " <<accepted/attempted<< endl << endl;
	}
	
	
	double stima_H=_blkH/_blknorm;
	_globH+=stima_H;
	_globH2+=stima_H*stima_H;
	double err_H=Error(_globH,_globH2,iblk);
	if(k==1){
		ofstream Out;
		Out.open("stima_H.dat",ios::app);
		Out<<iblk<<" "<< stima_H<<" "<< _globH/(double)iblk <<" "<< err_H << endl;
		Out.close();
	}
	return;
}


double Quantum::Error(double sum, double sum2, int iblk){
	 if( iblk == 1 ) return 0.0;
	 else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}




void Quantum::SetParametres(double mu, double sigma){
	_mu=mu;
	_sigma=sigma;
	return;
}
	
	
	

	 
	

