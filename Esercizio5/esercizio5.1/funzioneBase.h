#ifndef _funzioneBase_h_
#define _funzioneBase_h_

//classe astratta
class funzioneBase {

public:
	virtual double eval(double x)const=0;

    	//virtual double Eval (double *vettore)const=0;       //metodi che servono per costruire volumi o comunque spazi in più dimensioni (vedi sfera)
	//virtual double getdim() =0;
	//virtual double getr() =0;

	
//il metodo virtual posto a zero implica che non posso creare oggetti di questa classe...le classi figlie implementano tale metodo

};
#endif 
