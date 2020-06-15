#include<fstream>
#include<iostream>
#include<cmath>
#include "random.h"
#include "Ising.h"

using namespace std;

Ising::Ising(Random *r){
	ifstream ReadInput;

  	cout << "Classic 1D Ising model             " << endl;
  	cout << "Monte Carlo simulation             " << endl << endl;
  	cout << "Nearest neighbour interaction      " << endl << endl;
  	cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  	cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read input informations
	_kb=1.380649*pow(10.,-23.);
	_rnd=r;
  	ReadInput.open("input.dat");

  	ReadInput >> _temp;
  	_beta = 1.0/_temp;
  	cout << "Temperature = " << _temp << endl;

  	ReadInput >> _nspin;
  	cout << "Number of spins = " << _nspin << endl;

  	ReadInput >> _J;
  	cout << "Exchange interaction = " << _J << endl;

  	ReadInput >> _h;
  	cout << "External field = " << _h << endl << endl;
    
  	ReadInput >> _metro; // if=1 Metropolis else Gibbs

  	ReadInput >> _nblk;

  	ReadInput >> _nstep;

  	if(_metro==1) cout << "The program perform Metropolis moves" << endl;
  	else cout << "The program perform Gibbs moves" << endl;
  	cout << "Number of blocks = " << _nblk << endl;
  	cout << "Number of steps in one block = " << _nstep << endl << endl;
  	ReadInput.close();

//initial configuration(T=infinito)
	_s = new int [_nspin];
	_sm= new int [_nspin];
						
	

//Prepare arrays for measurements
  	_iu = 0; //Energy
  	_ic = 1; //Heat capacity
  	_im = 2; //Magnetization
  	_ix = 3; //Magnetic susceptibility
	n_props = 4; //Number of observables

	_walker=new double [n_props];
	glob_av=new double [n_props];
	glob_av2=new double [n_props];
	blk_av=new double [n_props];

	return;


}	




void Ising::SetSpin(int opz){
	if (opz==0){
		for (int i=0; i<_nspin; ++i){
    		if(_rnd->Rannyu() >= 0.5) _s[i] = 1;
    		else _s[i] = -1;
		_sm[i]=_s[i];

 	 	}
	}
	if(opz==1){
		ifstream Input;
		Input.open("config.final");
		for(int i=0;i<_nspin;i++){
			Input>>_s[i]>>_sm[i];
		}
		Input.close();
	}
	
	return;
}




void Ising::Move(int _metro){
	int o;
  	double oene,alpha,r;   //variabili per simulazione con h=0
	double pUp,pDown;
	double oenem,alpham;   //variabili per simulazione con h=0.02 (che serviranno per calcolare la magnetizzazione)
	double pUpm,pDownm;

  	for(int i=0; i<_nspin; ++i){
//Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    		o = (int)(_rnd->Rannyu()*_nspin);

   		if(_metro==1){          //Metropolis
			int sm;
			int random=int(_rnd->Rannyu(0.,2.));  //per decidere come flippare lo spin
			if(random==0){sm=-1;}
			else{sm=1;}
			oene=Boltzmann(sm,o)-Boltzmann(_s[o],o);     //differenza di energia
			oenem=(Boltzmann(sm,o)-0.02*sm) - (BoltzmannM(_sm[o],o)-0.02*_sm[o]);
			alpha=min(1.,exp(-_beta*oene));
			alpham=min(1.,exp(-_beta*oenem));
			r=_rnd->Rannyu();
			if(r<=alpha){_s[o]=sm;}
			if(r<=alpham){_sm[o]=sm;}
	
    		}
    		else {              //Gibbs sampling
		
    			pUp=1./(1.+exp(_beta*(Boltzmann(-1,o)-Boltzmann(1,o))));
			pDown=1./(1.+exp(_beta*(Boltzmann(1,o)-Boltzmann(-1,o))));
			pUpm=1./(1.+exp(-_beta*(Boltzmann(-1,o)-Boltzmann(1,o)+0.04)));
			pDownm=1./(1.+exp(_beta*(Boltzmann(1,o)-Boltzmann(-1,o)-0.04)));
			//if(pUp+pDown!=1.){cout<<"ERRORE! La somma di pUp e pDown in Gibbs è diversa da 1!!"<<endl;}
			r=_rnd->Rannyu();
			if(pDown<=pUp){
				if(r<pDown){ _s[o]=1;}
				if(r>=pDown){ _s[o]=-1;}
			}
			if(pDownm<=pUpm){
				if(r<pDownm){ _sm[o]=1;}
				if(r>=pDownm){ _sm[o]=-1;}
			}
			if(pDown>pUp){
				if(r<pUp){ _s[o]=-1;}
				if(r>=pUp){ _s[o]=1;}	
			}
			if(pDownm>pUpm){
				if(r<pUpm){_sm[o]=-1;}
				if(r>=pUpm){_sm[o]=1;}
		 
			
   			}
  		}
	}
	return;
}



double Ising::Boltzmann(int sm, int ip){
	double ene = -_J * sm * ( _s[Pbc(ip-1)] + _s[Pbc(ip+1)] ) - _h * sm;
  	return ene;
}



double Ising::BoltzmannM(int sm, int ip){
	double ene = -_J * sm * ( _sm[Pbc(ip-1)] + _sm[Pbc(ip+1)] ) - _h * sm;
  	return ene;
}




int Ising::Pbc(int i){
	if(i >= _nspin) i = i - _nspin;
    	else if(i < 0) i = i + _nspin;
    	return i;
}




void Ising::Measure(){
	//int bin;
  	double u = 0.0, m2 = 0.0;
	double u2=0.0,m=0.0;
	double en;
	

	for (int i=0; i<_nspin; ++i){
		en=-_J * _s[i] * _s[Pbc(i+1)] - 0.5 * _h * (_s[i] + _s[Pbc(i+1)]);
     		u +=en;
		m2+=_s[i];
		m+=_sm[i];
  	}
	u2=pow(u,2.);
  	_walker[_iu] = u;
	_walker[_ic] = u2;
	_walker[_ix] = pow(m2,2.);
	_walker[_im] = m; 
	// INCLUDE YOUR CODE HERE

	return;
}




void Ising::Reset(int iblk){
	if(iblk == 1){
		for(int i=0; i<n_props; ++i){
           		glob_av[i] = 0;
           		glob_av2[i] = 0;
       		}
	}

	for(int i=0; i<n_props; ++i){
		blk_av[i] = 0;
	}
	blk_norm = 0;
   	attempted = 0;
   	accepted = 0;

	return;
}





void Ising:: Accumulate(){
	for(int i=0; i<n_props; ++i){
		blk_av[i] = blk_av[i] + _walker[i];
	}
	blk_norm = blk_norm + 1.0;
}




void Ising::Averages(int iblk,int metro){
	ofstream Ene, Heat, Mag, Chi;
	//cout << "Block number " << iblk << endl;
    	//cout << "Acceptance rate " << accepted/attempted << endl << endl;

	Ene.open("output.ene.0",ios::app);
    	double stima_u = blk_av[_iu]/blk_norm/(double)_nspin; //Energy
    	glob_av[_iu]  += stima_u;
    	glob_av2[_iu] += stima_u*stima_u;
    	double err_u=Error(glob_av[_iu],glob_av2[_iu],iblk);
    	Ene <<  iblk <<"  "<< stima_u <<" "<< glob_av[_iu]/(double)iblk <<" "<< err_u <<" "<<metro<<endl;
    	Ene.close();
	

	Heat.open("output.heat.0",ios::app);
	double stima_u2=blk_av[_ic]/blk_norm/(double)_nspin; //Energy^2
	double stima_heat=pow(_beta,2.)*(stima_u2-pow(stima_u*(double)_nspin,2.)/(double)_nspin);   //Heat Capacity
	glob_av[_ic]  += stima_heat;
    	glob_av2[_ic] += stima_heat*stima_heat;
	double err_heat=Error(glob_av[_ic],glob_av2[_ic],iblk);
	Heat<<  iblk <<"  "<< stima_heat <<" "<< glob_av[_ic]/(double)iblk <<" "<< err_heat <<" "<<metro<<endl;
	Heat.close();


	Chi.open("output.chi.0",ios::app);
	double stima_chi=_beta*(blk_av[_ix]/blk_norm/(double)_nspin);  //Suscettibilità magnetica
	glob_av[_ix]  += stima_chi;
    	glob_av2[_ix] += stima_chi*stima_chi;
	double err_chi=Error(glob_av[_ix],glob_av2[_ix],iblk);
	Chi<<  iblk <<"  "<< stima_chi <<" "<< glob_av[_ix]/(double)iblk <<" "<< err_chi <<" "<<metro<< endl;
	Chi.close();

	Mag.open("output.mag.0",ios::app);
	double stima_mag=blk_av[_im]/blk_norm/(double)_nspin;
	glob_av[_im]  += stima_mag;
    	glob_av2[_im] += stima_mag*stima_mag;
	double err_mag=Error(glob_av[_im],glob_av2[_im],iblk);
	Mag<<  iblk <<"  "<< stima_mag <<" "<< glob_av[_im]/(double)iblk <<" "<< err_mag <<" "<<metro<< endl;
	Mag.close();
	

		
}




void Ising::ConfFinal(){
	ofstream WriteConf;

  	//cout << "Print final configuration to file config.final " << endl << endl;
  	WriteConf.open("config.final");
  	for (int i=0; i<_nspin; ++i){
		WriteConf << _s[i] <<" "<< _sm[i] << endl;
	}

	WriteConf.close();
}





double Ising::Error(double sum, double sum2, int iblk){
	return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}
	

void Ising::SetTemp(double temp){
	_temp=temp;
	_beta=1./temp;
	return;
}
	
