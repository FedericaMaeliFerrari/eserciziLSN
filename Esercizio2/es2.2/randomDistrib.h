#ifndef _randomDistrib_h_
#define _randomDistrib_h_


class RandomDistrib {

public:
	  double exp_trasformata(double lambda,double xi);
	  double lorentz_trasformata(double mu, double gamma, double xi);
	  double sen_trasformata(double xi);

private:
	double _xi;
	double _mu;
	double _gamma;
	double _lambda;
	
};

#endif
