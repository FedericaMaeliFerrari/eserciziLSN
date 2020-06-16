#include "GA.h"
#include "random.h"
#include "iostream"
#include<fstream>
#include <cmath>

using namespace std;


GA::GA(Random* rand,int metro){
	_rnd=rand;
	_num=32;
	_metro=metro;
	
	_x=new double [32];
	_y=new double [32];
	_indici= new int [32];
	_distanza=new double [100];
	_percorsi= new int* [100];
	for(int i=0;i<100;i++){
		_percorsi[i]=new int [32];
	}
	_percorsiOld= new int* [100];
	for(int i=0;i<100;i++){
		_percorsiOld[i]=new int [32];
	}

	_walker=new double [2];
	glob_av=new double [2];
	glob_av2=new double [2];
	blk_av=new double [2];
//CREO CITTA' VINCOLATE ALLA CIRCONFERENZA DI RAGGIO 1 CENTRATA IN 0 E NE SALVO LE COORDINATE
	if(_metro==0){
		for(int i=0;i<_num;i++){
			_indici[i]=i;
			double r=_rnd->Rannyu(0.,2.);
			double theta=acos(1.-r);
			int opz=(int)_rnd->Rannyu(0.,2.);
			if(opz==0){
				_x[i]=cos(theta);
				_y[i]=sin(theta);
			}
			if(opz==1){
				_x[i]=cos(theta+M_PI);
				_y[i]=sin(theta+M_PI);
			}
		}
	}
//CREO CITTÀ ISCRITTE IN UN QUADRATO DI LATO 2
	else{
		for(int i=0;i<_num;i++){
			_indici[i]=i;
			_x[i]=_rnd->Rannyu(0.,2.);
			_y[i]=_rnd->Rannyu(0.,2.);
		}
	}

//CONTROLLO CHE TUTTE LE CITTA' SIANO IN COORDINATE DIVERSE
	for(int i=0;i<_num-1;i++){
		for(int j=i+1;j<_num;j++){
			if(_x[i]==_x[j]){
				if(_y[i]==_y[j]){cout<<"ERRORE! DUE CITTÀ HANNO LE STESSE COORDINATE!!"<<endl;}
			}
		}
	}
		

//CREO 100 PERCORSI SFRUTTANDO L'OPERATORE PERMUTAZIONI PARI--> POPOLAZIONE INIZIALE
	for(int k=0;k<100;k++){
		for(int i=0;i<_num;i++){
			_percorsi[k][i]=_indici[i];     //crea 100 percorsi tutti uguali
		}
	}

	for(int j=1;j<100;j++){
		permPari(j);				//modifica i percorsi creati
	}

//CALCOLO LA DISTANZA PER OGNI PERCORSO
	for(int i=0;i<100;i++){
		_distanza[i]=Lx(i);
	}
	
//ORDINO I PERCORSI DA QUELLI CON DISTANZA (L(x)) MINORE A L(x) MAGGIORE
	Ordina();		
	
}


void GA::Mutazione(){
	double pari,cont,trasl,inv;
	for(int i=0;i<100;i++){
		pari=_rnd->Rannyu();
		cont=_rnd->Rannyu();
		trasl=_rnd->Rannyu();
		inv=_rnd->Rannyu();
		if(pari<0.1){ permPari(i);}	//faccio la mutazione Pari sul percorso i-esimo
		if(cont<0.1){ permCont(i);}
		if(trasl<0.1){ permTrasl(i);}
		if(inv<0.1){ inversione(i);}
	}
	int controllo=check();
	if(controllo !=1){cout<<"MUTAZIONE HA AVUTO UN ERRORE!"<<endl;}
	//if(controllo==1){cout<<"MUTAZIONE ESEGUITA PERFETTAMENTE"<<endl;}
	if(controllo==-1){cout<<"ERRORE! IL PRIMO ELEMENTO DEL PERCORSO NON E' ZERO(CITTÀ INIZIALE)"<<endl;}
	if(controllo==2){cout<<"ERRORE! CI SONO CITTÀ UGUALI NELLO STESSO STESSO PERCORSO"<<endl;}	
	if(controllo==-2){cout<<"ERRORE! HAI SBAGLIATO TUTTO"<<endl;}
	
	for(int i=0;i<100;i++){
		_distanza[i]=Lx(i);
	}
	
	Ordina();	

}


void GA::NewGeneration(){
	for(int i=0;i<100;i++){					//matrice di appoggio in cui salvo i valori nuovi della popolazione man mano;
		for(int j=0;j<_num;j++){
			_percorsiOld[i][j]=_percorsi[i][j];
		}
	}

	int idx,idx2;
	double prob; 
	for(int i=0;i<50;i++){
		idx=Select();
		idx2=Select();
		prob=_rnd->Rannyu();
		if(prob<0.7){
			crossover(idx,idx2,(int)i*2,(int)(i*2)+1);
		}
		else{
			for(int j=0;j<_num;j++){
				_percorsi[(int)i*2][j]=_percorsiOld[idx][j];
				_percorsi[(int)(i*2)+1][j]=_percorsiOld[idx2][j];
			}
		}
	}

	
	int controllo=check();
	if(controllo !=1){cout<<"CROSSOVER HA AVUTO UN ERRORE!"<<endl;}
	//if(controllo==1){cout<<"CROSSOVER ESEGUITO PERFETTAMENTE"<<endl;}
	if(controllo==-1){cout<<"ERRORE! IL PRIMO ELEMENTO DEL PERCORSO NON E' ZERO(CITTÀ INIZIALE)"<<endl;}
	if(controllo==2){cout<<"ERRORE! CI SONO CITTÀ UGUALI NELLO STESSO STESSO PERCORSO"<<endl;}	
	if(controllo==-2){cout<<"ERRORE! HAI SBAGLIATO TUTTO"<<endl;}
	
	
	for(int i=0;i<100;i++){
		_distanza[i]=Lx(i);
	}
	
	Ordina();

	_walker[0]=_distanza[0];		//salvo il valore della distanza del percorso più breve
	double sum=0;
	for(int i=0;i<100;i++){
		sum+=_distanza[i];
	}
	_walker[1]=sum;				//salvo la somma di tutti i percorsi;
	return;
}	
	

void GA::permPari(int ind){
	int r=(int)_rnd->Rannyu(1.,_num-1.);	//genero casualmente il punto del vettore in cui devo scambiare il valore selezionato è quello successivo.
	int a=_percorsi[ind][r];
	int b=_percorsi[ind][r+1];
	_percorsi[ind][r]=b;
	_percorsi[ind][r+1]=a;

	return;
}



double GA::Lx(int ind){
	double sum=0.;
	double dx,dy;
	for(int i=0;i<_num;i++){
		if(i==_num-1){
			dx=_x[_percorsi[ind][i]]-_x[_percorsi[ind][0]];	
			dy=_y[_percorsi[ind][i]]-_y[_percorsi[ind][0]];
		}
		else{
			dx=_x[_percorsi[ind][i]]-_x[_percorsi[ind][i+1]];	//calcolo la distanza sull'asse x dei due indici consecutivi del percorso numero "ind"
			dy=_y[_percorsi[ind][i]]-_y[_percorsi[ind][i+1]];	//idem per l'asse y					
		}
		sum+=pow((pow(dx,2.)+pow(dy,2.)),0.5);		//calcolo la distanza tra le due città e le sommo alle sistanze precedenti-->funzione di costo
	}
	return sum;
}



void GA::Ordina(){
	int indOrdinamento[100];	//creo un vettore di indici per tener conto di come si sposta il vettore _distanza per ordinarsi in modo crescente
	double appoggio;
	int appoggioInd;

	for(int i=0;i<100;i++){
		indOrdinamento[i]=i;   //vettore degli indici [0,1,2,.....99]
	}

			
	for(int i=0;i<99;i++){
		for(int j=i+1;j<100;j++){
			if(_distanza[j]<_distanza[i]){
				appoggio=_distanza[j];
				appoggioInd=indOrdinamento[j];
				_distanza[j]=_distanza[i];
				_distanza[i]=appoggio;
				indOrdinamento[j]=indOrdinamento[i];
				indOrdinamento[i]=appoggioInd;
			}
		}
	}

//CONTROLLO SE EFFETTIVAMENTE IL VETTORE DELLE DISTANZE E' ORDINATO
	for(int i=0;i<99;i++){
		if(_distanza[i]>_distanza[i+1]){cout<<"ERRORE IN POSIZIONE "<<i<<endl;}
		//cout<<_distanza[i]<<endl;
	}
	

//A questo punto dovrei avere il vettore delle distanze L(x) ordianto  e un vettore di indici che tiene conto di come è stato ordinato--->seguendo l'ordine degli indici dovrei poter ordinare 
	int appPercorsi[100][_num];
	//creo una matrice di appoggio per i percorsi uguale a _percorsi
	for(int i=0;i<100;i++){
		for(int j=0;j<_num;j++){
			appPercorsi[i][j]=_percorsi[i][j];
		}
	}
	
	//ordino i percorsi secondo gli indici che tengono conto di in che modo era stato ordinato il vettore delle distanze
	for(int i=0;i<100;i++){
		appoggioInd=indOrdinamento[i];
		for(int j=0;j<_num;j++){
			_percorsi[i][j]=appPercorsi[appoggioInd][j];
		}
	}

	return;
		
}
			

void GA::permTrasl(int ind){
	int n=(int)_rnd->Rannyu(1.,10.);		//quanti step traslare il vettore
	int m=(int)_rnd->Rannyu(1.,_num-1.);	//quante città traslare, partendo dalla città che occupa la pos 1 nel vettore percorsi (non quella iniziale) a quella che occupa la posizione m
	if(n+m<_num-1){
		int vettApp[n];
		for(int i=0;i<n;i++){
			vettApp[i]=_percorsi[ind][m+i+1];   //salvo i numeri non appartenenti alle m città che però finiranno sovrascritti durante la traslazione
		}
		for(int j=m;j>=1;j--){
			_percorsi[ind][j+n]=_percorsi[ind][j];	//traslo le m città di n step
		}
		for(int i=0;i<n;i++){
			_percorsi[ind][i+1]=vettApp[i];        //inserisco in testa le città sovrascritte durante la traslazione
		}
	}
	/*
	else{
		int posEccesso=n+m-(_num-1);
		int posEccesso2=m-posEccesso;
		int vettApp2[posEccesso2];
		if(posEccesso2 != 0){
			for(int i=0;i<posEccesso2;i++){
				vettApp2[i]=_percorsi[ind][_num-1-i]; //salvo le città che verrebbero sovrascritte in quanto occupano posti in cui andranno alcune delle m città (quelle che non si riavvolgono)
			}
		}
		int vettApp[posEccesso];
		for(int i=0;i<posEccesso;i++){
			vettApp[i]=_percorsi[ind][i+1];      //salvo le città che verrebbero sovrascritte in quanto ho un riavvolgimento del vettore
		}
		for(int j=m;j>posEccesso;j--){
			_percorsi[ind][j+n-(_num-1)]=_percorsi[ind][j];		//sovrascrivo le città che si riavvolgono nel vettore
		}
		for(int i=0;i<posEccesso;i++){
			_percorsi[ind][i+1+n]=vettApp[i];			//sovrascrivo le città che si traslano senza riavvolgersi (occupando i posti finali del vettore)
		}
		if(posEccesso2 !=0){
			for(int i=0;i<posEccesso2;i++){
				_percorsi[ind][m-posEccesso2+1+i]=vettApp2[i];	//sovrascrivo le città sovrascritte che non appartenevano alle m città
			}
		}
	}
	*/
	return;				

}



void GA:: permCont(int ind){
	int m =(int)_rnd->Rannyu(1.,_num/2.);
	int app;
	for(int j=1;j<=m;j++){
		app=_percorsi[ind][m+j];
		_percorsi[ind][m+j]=_percorsi[ind][j];
		_percorsi[ind][j]=app;
	}
	return;
}




void GA::inversione(int ind){
	int m=(int)_rnd->Rannyu(1.,_num);
	int app;
	if(m%2==0){
		int cont=1;	
		for(int j=m;j>m/2;j--){
			app=_percorsi[ind][j];
			_percorsi[ind][j]=_percorsi[ind][cont];
			_percorsi[ind][cont]=app;
			cont=cont+1;
		}
	}
	else{
		int cont2=1;
		for(int j=m;j>(m+1)/2;j--){
			app=_percorsi[ind][j];
			_percorsi[ind][j]=_percorsi[ind][cont2];
			_percorsi[ind][cont2]=app;
			cont2=cont2+1;
		}
	}
	return;	
}



int GA::check(){
	int test=1;
	int test2=1;
	for(int i=0;i<100;i++){
		if(_percorsi[i][0] != 0){ test=-1;}
		for(int j=0;j<_num-1;j++){
			for(int k=j+1;k<_num;k++){
				if(_percorsi[i][j]==_percorsi[i][k]){test2=2;}
			}
		}
	}
	return test*test2;
}
	

int GA::Select(){
	double r=_rnd->Rannyu();
	int j=(int)(100*pow(r,2));
	return j;
}


void GA:: crossover(int idx, int idx2,int cont, int cont2){
	int cut=(int)_rnd->Rannyu(1,_num-2);     //indice in cui andranno tagliati i vettori
	int temp[(int)_num-cut];      		//vettore che terrà conto dell'ordine in cui appaiono nell'altro vettore i valori da crossoverare
	int temp2[(int)_num-cut];
	for(int i=cut;i<_num;i++){
		for(int j=0;j<_num;j++){
			if(_percorsiOld[idx][i]==_percorsiOld[idx2][j]){
				temp[(int)i-cut]=j;			//Adesso in temp e in temp2 ho salvato l'rodine  (ovvero l'indice) in cui appaiono i valori che voglio crossoverare dei due vettori
			}
			if(_percorsiOld[idx2][i]==_percorsiOld[idx][j]){
				temp2[(int)i-cut]=j;
			}
		}
	}
	
	//ordino i vettori temp e temp2 in ordine crescente, in modo da sapere in che ordine dovranno essere reinseriti i valori tagliati dai valori iniziali
	int appoggio,appoggio2;
	for(int i=0;i<(int)_num-cut-1;i++){
		for(int j=i+1;j<(int)_num-cut;j++){
			if(temp[j]<temp[i]){
				appoggio=temp[j];
				temp[j]=temp[i];
				temp[i]=appoggio;
			}
			if(temp2[j]<temp2[i]){
				appoggio2=temp2[j];
				temp2[j]=temp2[i];
				temp2[i]=appoggio2;
			}
		}
	}

	//faccio il crossover
	for(int i=0;i<_num;i++){
		_percorsi[cont][i]=_percorsiOld[idx][i];
		_percorsi[cont2][i]=_percorsiOld[idx2][i];
	}
	for(int i=cut;i<_num;i++){
		_percorsi[cont][i]=_percorsiOld[idx2][temp[(int)i-cut]];
		_percorsi[cont2][i]=_percorsiOld[idx][temp2[(int)i-cut]];
	}
	
	return;
}




void GA::Coord(){
	ofstream Out;
	if(_metro==0){
		Out.open("Città_circonferenza.dat");
	}
	else{
		Out.open("Città_quadrato.dat");
	}
	for(int i=0;i<_num;i++){
		Out<<_indici[i]<<" "<<_x[i]<<" "<<_y[i]<<endl;
	}
	Out.close();
	return;
}


void GA::BestPath(){
	ofstream BP;
	if(_metro==0){
		BP.open("bestpath_circ.dat");
	}
	else{
		BP.open("bestpath_quadr.dat");
	}
	for(int i=0;i<_num;i++){
		BP<<_percorsi[0][i]<<endl;
	}
	BP.close();
	return;
}
	
	
		
void GA::Reset(int iblk){
	if(iblk == 1){
		for(int i=0; i<2; ++i){
           		glob_av[i] = 0;
           		glob_av2[i] = 0;
       		}
	}

	for(int i=0; i<2; ++i){
		blk_av[i] = 0;
	}
	blk_norm = 0;

	return;
}


void GA::Accumulate(){
	for(int i=0; i<2; ++i){
		blk_av[i] = blk_av[i] + _walker[i];
	}
	blk_norm = blk_norm + 1.0;
	return;
}



void GA::Averages(int iblk){
	ofstream Lx,Lxm;
	if(_metro==0){
		Lx.open("Lx_circ.dat",ios::app);
	}
	else{
		Lx.open("Lx_quadr.dat",ios::app);
	}
		
	double stima_l=blk_av[0]/blk_norm;
	glob_av[0]+=stima_l;
	glob_av2[0]+=stima_l*stima_l;
	double err_l=Error(glob_av[0],glob_av2[0],iblk);
	Lx <<  iblk <<"  "<< stima_l <<" "<< glob_av[0]/(double)iblk <<" "<< err_l<<endl;
	Lx.close();

	if(_metro==0){
		Lxm.open("Lx_media_circ.dat",ios::app);
	}
	else{
		Lxm.open("Lx_media_quadr.dat",ios::app);
	}
	double stima_lm=blk_av[1]/blk_norm/100;
	glob_av[1]+=stima_lm;
	glob_av2[1]+=stima_lm*stima_lm;
	double err_lm=Error(glob_av[1],glob_av2[1],iblk);
	Lxm <<  iblk <<"  "<< stima_lm <<" "<< glob_av[1]/(double)iblk <<" "<< err_lm<<endl;
	Lxm.close();

	return;
}

	

double GA::Error(double sum, double sum2, int iblk){
	return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}
			
		
	
			
			
	
	



		
