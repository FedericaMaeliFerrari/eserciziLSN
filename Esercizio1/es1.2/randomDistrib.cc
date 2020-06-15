#include "randomDistrib.h"
#include <cmath>

double RandomDistrib::exp_trasformata(double lambda, double xi){
	_lambda=lambda;
	double x= -1./_lambda*log(1.-xi);
	return x;
}

double RandomDistrib::lorentz_trasformata(double mu, double gamma, double xi){
	_mu=mu;
	_gamma=gamma;
	double x=_gamma*tan(M_PI*(xi-0.5))+_mu;
	return x;
}
