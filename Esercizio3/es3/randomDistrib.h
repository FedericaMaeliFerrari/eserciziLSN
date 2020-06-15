#ifndef _randomDistrib_h_
#define _randomDistrib_h_


class RandomDistrib {

public:
	  RandomDistrib(){};
	  double exp_trasformata(double lambda,double xi);
	  double lorentz_trasformata(double mu, double gamma, double xi);
	  double gauss_trasformata(double mu,double sigma,double xi,double xi2);

private:
	double _xi;
	double _mu;
	double _gamma;
	double _lambda;	
	double _sigma;
};

#endif
