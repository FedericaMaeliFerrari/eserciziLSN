#include "finanza.h"
#include <iostream>
#include <cmath>

using namespace std;

Finanza::Finanza(double t,double S0,double T,double K,double r,double sigma,integrale* i){
	_t=t;
	_S0=S0;
	_T=T;
	_K=K;
	_r=r;
	_sigma=sigma;
	_i=i;
	_d1=1./(_sigma*sqrt(_T-_t))*(log(_S0/_K)+_r+pow(_sigma,2.)/2.*(_T-_t));
	_d2=_d1-_sigma*sqrt(_T-_t);
	double integraled1=(_i->integraleByMedia(0., _d1/sqrt(2.),1000))*2./sqrt(M_PI);
	double integraled2=(_i->integraleByMedia(0., _d2/sqrt(2.),1000))*2./sqrt(M_PI);
	double nd1=1./2.*(1.+integraled1);
	double nd2=1./2.*(1.+integraled2);
	double c=_S0*nd1-_K*exp(-_r*(_T-_t))*nd2;
	double p=_S0*(nd1-1.)-_K*exp(-_r*(_T-_t))*(nd2-1.);
	_c=c;
	_p=p;
}


double Finanza::direct_sampleC(int L,double mu,RandomDistrib* g,Random* r){
	double sumCi=0.;
	for(int i=0;i<L;i++){
		double xi=r->Rannyu();
		double xi2=r->Rannyu();
		double x=g->gauss_trasformata(mu,_T,xi,xi2);
		double St=_S0*exp((_r-1./2.*pow(_sigma,2.))*_T+_sigma*x);
		if(St-_K>0.){sumCi=sumCi+exp(-_r*_T)*(St-_K);}
		else{sumCi=sumCi+0.;}	
	}
	double C_media=sumCi/L;
	return C_media;	

}


double Finanza::direct_sampleP(int L,double mu,RandomDistrib* g,Random* r){
	double sumPi=0.;
	for(int i=0;i<L;i++){
		double xi=r->Rannyu();
		double xi2=r->Rannyu();
		double x=g->gauss_trasformata(mu,_T,xi,xi2);
		double St=_S0*exp((_r-1./2.*pow(_sigma,2.))*_T+_sigma*x);
		if(_K-St>0.){sumPi=sumPi+exp(-_r*_T)*(_K-St);}
		else{sumPi=sumPi+0.;}	
	}
	double P_media=sumPi/L;
	return P_media;	

}


double Finanza::discret_sampleC(double passi,int L,RandomDistrib* g,Random* r){
	double dt=_T/passi;
	double sumCi=0.;
	for(int i=0;i<L;i++){
		double St_bef=_S0;
		for(int j=0;j<passi;j++){
			double xi=r->Rannyu();
			double xi2=r->Rannyu();
			double x=g->gauss_trasformata(0.,1.,xi,xi2);
			double St=St_bef*exp((_r-1./2.*pow(_sigma,2.))*dt+_sigma*x*sqrt(dt));
			St_bef=St;
		}
		if(St_bef-_K>0.){sumCi=sumCi+exp(-_r*_T)*(St_bef-_K);}
		else{sumCi=sumCi+0.;}	
	}
	double C_media=sumCi/L;
	return C_media;	
}

double Finanza::discret_sampleP(double passi,int L,RandomDistrib* g,Random* r){
	double dt=_T/passi;
	double sumPi=0.;
	for(int i=0;i<L;i++){
		double St_bef=_S0;
		for(int j=0;j<passi;j++){
			double xi=r->Rannyu();
			double xi2=r->Rannyu();
			double x=g->gauss_trasformata(0.,1.,xi,xi2);
			double St=St_bef*exp((_r-1./2.*pow(_sigma,2.))*dt+_sigma*x*sqrt(dt));
			St_bef=St;
		}
		if(_K-St_bef>0.){sumPi=sumPi+exp(-_r*_T)*(_K-St_bef);}
		else{sumPi=sumPi+0.;}	
	}
	double P_media=sumPi/L;
	return P_media;	
}
