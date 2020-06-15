
#include "funzioneBase.h"
#include "funzione.h"
#include <cmath>

using namespace std;

double coseno::eval(double x)const{
	double f=pow(_a*cos(_b*x*M_PI),_c);
	return f;
	
}
coseno::coseno(double a,double b,double c){
	_a=a;
	_b=b;
	_c=c;

}

coseno::coseno(){
	_a=1.;
	_b=1.;
	_c=1.;

}
//double funzione::Eval (double *vettore)const{}      //metodi che servono per costruire volumi o comunque spazi in più dimensioni (vedi sfera)
//double funzione::getdim() {}
//double funzione::getr() {}



double seno::eval(double x)const{
	double f=pow(_a*sin(_b*x*M_PI),_c);
	return f;
}

seno::seno(double a,double b,double c){
	_a=a;
	_b=b;
	_c=c;
}

seno::seno(){
	_a=1.;
	_b=1.;
	_c=1.;
}


double parabola::eval(double x)const{
	double f=_a*pow(x,2.)+_b*x+_c;
	return f;
}

parabola::parabola(double a,double b,double c){
	_a=a;
	_b=b;
	_c=c;
}


double retta::eval(double x)const{
	double f=_d*pow((_a*x+_b),_pot)+_c;
	return f;
}

retta::retta(double a, double b, double pot, double c, double d){
	_a=a;
	_b=b;
	_pot=pot;
	_c=c;
	_d=d;
}


retta::retta(double a, double b){
	_a=a;
	_b=b;
	_pot=1.;
	_c=0.;
	_d=1.;
}

esponenziale::esponenziale(double a,double b,double pot,double c){
	_a=a;
	_b=b;
	_c=c;
	_pot=pot;
	
}

double esponenziale::eval(double x)const{
	double f=_a*exp(_b*pow(x,_pot)+_c);
	return f;
}


double cosretta::eval(double x)const{
	double f=pow(_acos*cos(_bcos*x*M_PI),_ccos)/(_d*(pow((_a*x+_b),_pot)+_c));
	return f;
}

cosretta::cosretta(double a,double b,double pot,double c,double d,double acos,double bcos,double ccos){
	_a=a;
	_b=b;
	_pot=pot;
	_c=c;
	_d=d;
	_acos=acos;
	_bcos=bcos;
	_ccos=ccos;
}

cosretta::cosretta(double a, double b){
	_a=a;
	_b=b;
	_pot=1.;
	_c=0.;
	_d=1.;
	_acos=1.;
	_bcos=1.;
	_ccos=1.;
}




