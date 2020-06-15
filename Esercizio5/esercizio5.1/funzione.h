#ifndef _funzione_h_
#define _funzione_h_
#include "funzioneBase.h"

class coseno:public funzioneBase{

public:
	double eval (double x)const;
	coseno(double a,double b,double c);
	coseno();

	

    	//virtual double Eval (double *vettore)const;       //metodi che servono per costruire volumi o comunque spazi in più dimensioni (vedi sfera)
	//virtual double getdim() ;
	//virtual double getr() ;
protected:
	double _a,_b,_c;
	
};

#endif


#ifndef _seno_h
#define _seno_h

class seno:public funzioneBase{

public:

	double eval (double x)const;
	seno(double a, double b,double c);
	seno();
protected:
	double _a,_b,_c;
};
#endif



#ifndef _parabola_h
#define _parabola_h

class parabola:public funzioneBase{

public:
	double eval (double x)const;
	parabola(double a,double b,double c);
protected:
	double _a,_b,_c;
};
#endif



#ifndef _retta_h
#define _retta_h

class retta:public funzioneBase{

public:
	double eval (double x)const;
	retta(double a,double b, double pot, double c,double d);
	retta(double a,double b);
protected:
	double _a,_b,_pot,_c,_d;
};
#endif


#ifndef _esponenziale_h
#define _esponenziale_h

class esponenziale:public funzioneBase{

public:
	double eval (double x)const;
	esponenziale(double a,double b,double pot, double c);
protected:
	double _a,_b,_c,_pot;
	
};
#endif


#ifndef _cosretta_h
#define _cosretta_h

class cosretta:public funzioneBase{

public:
	double eval (double x)const;
	cosretta(double a,double b,double pot,double c,double d,double acos,double bcos, double ccos);
	cosretta(double a,double b);
protected:
	double _a,_b,_pot,_c,_d,_acos,_bcos,_ccos;
};
#endif 
	
