#ifndef _GA_h_
#define _GA_h_
#include "random.h"

class GA {

public:
	GA(Random* rand,int metro);
	int check();
	void NewGeneration();
	void NewGeneration_Metropolis(double temp);
	void Mutazione();

	void crossover(int idx,int idx2, int cont, int cont2);
	void permPari(int ind);		//permutazioni pari
	void permTrasl(int ind);	//permutazioni contigue con traslazione
	void permCont(int ind);		//permutazioni contigue
	void inversione(int ind);
	int Select();
	
	double Lx(int ind);     //Calcola la funzione di costo per il percorso numero ind;
	void Ordina();		//Ordina il vettore delle distanze in ordine cfrescente--->utilizzerò ciò per ordinare i percorsi di conseguenza

	void Coord();		//Scrive in un file di testo le coordinate delle città e i relativi indici
	void BestPath();	//Scrive gli indici delle città che minimizzano la distanza
	void Reset(int);
	void Accumulate();
	void Averages(int);
	double Error(double,double,int);

protected:
	Random* _rnd;
	int _num;
	int _metro;

	double *_x;
	double *_y;
	int *_indici;

	double *_distanza;
	double *_distanzaOld;
	int **_percorsi;
	int **_percorsiOld;

	double* _walker;
	double* glob_av;
	double* glob_av2;
	double* blk_av;
	int blk_norm;




};

#endif
