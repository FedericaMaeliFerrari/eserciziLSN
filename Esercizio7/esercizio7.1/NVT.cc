#include "NVT.h"
#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

NVT::NVT(const char*filename,Random* r,int k){
	  ifstream ReadInput,ReadConf;

	  cout << "Classic Lennard-Jones fluid        " << endl;
	  cout << "Monte Carlo simulation             " << endl << endl;
	  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
	  cout << "The program uses Lennard-Jones units " << endl;

	  _rnd=r;
//Read input informations
	  ReadInput.open(filename);

	  ReadInput >> _temp;
	  _beta = 1.0/_temp;
	  cout << "Temperature = " << _temp << endl;

	  ReadInput >> _npart;
	  cout << "Number of particles = " << _npart << endl;

	  ReadInput >> _rho;
	  cout << "Density of particles = " << _rho << endl;
	  _vol = (double)_npart/_rho;
	  _box = pow(_vol,1.0/3.0);
	  cout << "Volume of the simulation box = " << _vol << endl;
	  cout << "Edge of the simulation box = " << _box << endl;

	  ReadInput >> _rcut;
	  cout << "Cutoff of the interatomic potential = " << _rcut << endl << endl;
	  ReadInput >> _delta;
	  ReadInput >> _nblk;
	  ReadInput >> _nstep;

//Tail corrections for potential energy and pressure
	  _vtail = (8.0*M_PI*_rho)/(9.0*pow(_rcut,9)) - (8.0*M_PI*_rho)/(3.0*pow(_rcut,3));
	  _ptail = (32.0*M_PI*_rho)/(9.0*pow(_rcut,9)) - (16.0*M_PI*_rho)/(3.0*pow(_rcut,3));
	  cout << "Tail correction for the potential energy = " << _vtail << endl;
	  cout << "Tail correction for the virial           = " << _ptail << endl; 
	  cout << "The program perform Metropolis moves with uniform translations" << endl;
	  cout << "Moves parameter = " << _delta << endl;
	  cout << "Number of blocks = " << _nblk << endl;
	  cout << "Number of steps in one block = " << _nstep << endl << endl;

	  ReadInput.close();

	  _x = new double [_npart];
	  _y = new double [_npart];
	  _z = new double [_npart];
	  
//Prepare arrays for measurements
	  _iv = 0; //Potential energy
	  _iw = 1; //Virial
	 
	  n_props = 2; //Number of observables
	
//measurement of g(r)
	  _igofr = 2;
	  _nbins = 100;
	  n_props = n_props + _nbins;
	  bin_size = (_box/2.0)/(double)_nbins;

	  _walker = new double [n_props];
	  glob_av=new double [n_props];
	  glob_av2=new double [n_props];
	  blk_av=new double [n_props];
//Read initial configuration
	  cout << "Read initial configuration from file config.0 " << endl << endl;
	  ReadConf.open("config.0");
	  for (int i=0; i<_npart; ++i){
		    ReadConf >> _x[i] >> _y[i] >> _z[i];
		    _x[i] = Pbc( _x[i] * _box );
		    _y[i] = Pbc( _y[i] * _box );
		    _z[i] = Pbc( _z[i] * _box );
	  }
	  ReadConf.close();
  
//Evaluate potential energy and virial of the initial configuration
 	  Measure(k);

//Print initial values for the potential energy and virial
	  cout << "Initial potential energy (with tail corrections) = " << _walker[_iv]/(double)_npart + _vtail << endl;
	  cout << "Virial                   (with tail corrections) = " << _walker[_iw]/(double)_npart + _ptail << endl;
	  cout << "Pressure                 (with tail corrections) = " << _rho * _temp + (_walker[_iw] + (double)_npart * _ptail) / _vol << endl << endl;

	return;
}




void NVT::Move(){
	  int o;
  	  double p, energy_old, energy_new;

	  for(int i=0; i<_npart; ++i){
  //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
		    o = (int)(_rnd->Rannyu()*_npart);
	  //Old
		    _xold = _x[o];
		    _yold = _y[o];
		    _zold = _z[o];
		    energy_old = Boltzmann(_xold,_yold,_zold,o);
	  //New
		    _xnew = Pbc( _x[o] + _delta*(_rnd->Rannyu() - 0.5) );
		    _ynew = Pbc( _y[o] + _delta*(_rnd->Rannyu() - 0.5) );
		    _znew = Pbc( _z[o] + _delta*(_rnd->Rannyu() - 0.5) );
		    energy_new = Boltzmann(_xnew,_ynew,_znew,o);
	  //Metropolis test
		    p = exp(_beta*(energy_old-energy_new));
		    if(p >= _rnd->Rannyu()){
			       _x[o] = _xnew;
			       _y[o] = _ynew;
			       _z[o] = _znew;
			       accepted = accepted + 1.0;
	    	    }
	   	    attempted = attempted + 1.0;
	    }
	return;
}




double NVT::Boltzmann(double xx,double yy,double zz, int ip){
	double ene=0.0;
  	double dx, dy, dz, dr;
	for(int i=0;i<_npart;i++){
		if(i != ip){		// distance ip-i in pbc
			dx = Pbc(xx - _x[i]);
      			dy = Pbc(yy - _y[i]);
      			dz = Pbc(zz - _z[i]);
			dr = dx*dx + dy*dy + dz*dz;
      			dr = sqrt(dr);
			if(dr < _rcut){
				ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
			}
		}
	}
	return 4.0*ene;
}




double NVT::Pbc(double r){
	return r-_box*rint(r/_box);
}




void NVT::Measure(int k){
	ofstream EP;
	if(k==1){
		EP.open("Gas.dat",ios::app);  //stampa energia potenziale e pressione istantanea
		
	}
	if(k==2){
		EP.open("Liquido.dat",ios::app);  //stampa energia potenziale e pressione istantanea
	}
	if(k==3){
		EP.open("Solido.dat",ios::app);  //stampa energia potenziale e pressioneistantanea
	}
	//int bin;
  	double v = 0.0, w = 0.0;
  	double vij, wij;
  	double dx, dy, dz, dr;

	//reset the hystogram of g(r)
  	for (int k=_igofr; k<_igofr+_nbins; ++k){_walker[k]=0.0;}
	//cycle over pairs of particles
 	for (int i=0; i<_npart-1; ++i){
    		for (int j=i+1; j<_npart; ++j){          // distance i-j in pbc
     			dx = Pbc(_x[i] - _x[j]);
			dy = Pbc(_y[i] - _y[j]);
			dz = Pbc(_z[i] - _z[j]);
			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);
	//CODICE PER RIEMPIRE L'ISTOGRAMMA DI G(R)
			if(dr < _rcut){
				vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
      				wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
				v += vij;
       				w += wij;

			}
		}
	}
	_walker[_iv] = 4.0 * v;
  	_walker[_iw] = 48.0 * w / 3.0;

//Salvo valori istantanei di energia e pressione con le loro relative incertezze
	double emedia=_walker[_iv]/_npart+_vtail;
	double pmedia=_rho*_temp+(_walker[_iw]+(double)_npart*_ptail)/_vol;
	EP<<emedia<<" "<<pmedia<<endl;
	EP.close();

	return;
}




void NVT::Reset(int iblk){
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




void NVT::Accumulate(){
	for(int i=0; i<n_props; ++i){
		blk_av[i] = blk_av[i] + _walker[i];
	}
	 blk_norm = blk_norm + 1.0;
	return;
}




void NVT::Averages(int iblk){
	//double r, gdir;
   	ofstream Gofr, Gave, Epot, Pres;
	 
	cout << "Block number " << iblk << endl;
   	cout << "Acceptance rate " << accepted/attempted << endl << endl;

	Epot.open("output.epot.0",ios::app);
    	Pres.open("output.pres.0",ios::app);
    	Gofr.open("output.gofr.0",ios::app);
    	Gave.open("output.gave.0",ios::app);
    	
	double stima_pot = blk_av[_iv]/blk_norm/(double)_npart + _vtail; //Potential energy
    	glob_av[_iv] += stima_pot;
    	glob_av2[_iv] += stima_pot*stima_pot;
   	double err_pot=Error(glob_av[_iv],glob_av2[_iv],iblk);

	double stima_pres = _rho * _temp + (blk_av[_iw]/blk_norm + _ptail * (double)_npart) / _vol; //Pressure
        glob_av[_iw] += stima_pres;
    	glob_av2[_iw] += stima_pres*stima_pres;
   	double err_press=Error(glob_av[_iw],glob_av2[_iw],iblk);

//Potential energy per particle
    	Epot <<iblk<<" "<< stima_pot<<" "<< glob_av[_iv]/(double)iblk <<" "<< err_pot << endl;
//Pressure
    	Pres << iblk <<" "<< stima_pres <<" "<< glob_av[_iw]/(double)iblk <<" "<<err_press << endl;
 //CALCOLO G(R)

	 cout << "----------------------------" << endl << endl;

    	Epot.close();
    	Pres.close();
    	Gofr.close();
	return;
}




double NVT::Error(double sum, double sum2, int iblk){
	 if( iblk == 1 ) return 0.0;
	 else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}




void NVT::ConfFinal(){
	ofstream WriteConf;

  	cout << "Print final configuration to file config.final " << endl << endl;
  	WriteConf.open("config.final");
  	for (int i=0; i<_npart; ++i){
    		WriteConf << _x[i]/_box << "   " <<  _y[i]/_box << "   " << _z[i]/_box << endl;
 	 }
  	WriteConf.close();
}

void NVT::ConfXYZ(int nconf){
	ofstream WriteXYZ;
  	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  	WriteXYZ << _npart << endl;
  	WriteXYZ << "This is only a comment!" << endl;
  	for (int i=0; i<_npart; ++i){
    		WriteXYZ << "LJ  " << Pbc(_x[i]) << "   " <<  Pbc(_y[i]) << "   " << Pbc(_z[i]) << endl;
  	}
  	WriteXYZ.close();
}


			    














