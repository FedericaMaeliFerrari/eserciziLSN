#include "MD.h"
#include<fstream>
#include<iostream>
#include<cmath>
#include "random.h"
//#include "vettore.h"

using namespace std;

MD::MD(const char*filename,Random* r){
	ifstream input;
	input.open(filename);

	input>>_t;
	input>>_npart;
	input>>_rho;
	input>>_rcut;
	input>>_delta;
	input>>_nstep;
	input>>_nblk;
	_vol=_npart/_rho;
	_box=pow(_vol,1./3.);

	input.close();
	_rand=r;
	//_temp=0.;
	//_ekin=0.;
	//_etot=0.;
	//_epot=0.;
	 _iv = 0; //Potential energy
	 _iw = 1; //Virial
	 
	  //n_props = 2; //Number of observables
	
//measurement of g(r)
	  _igofr = 2;
	  _nbins = 100;
	  n_props = _igofr + _nbins;
	  bin_size = (_box/2.0)/(double)_nbins;


	_walker=new double [n_props];
	glob_av=new double [n_props];
	glob_av2=new double [n_props];
	blk_av=new double [n_props];

	 _vtail = (8.0*M_PI*_rho)/(9.0*pow(_rcut,9)) - (8.0*M_PI*_rho)/(3.0*pow(_rcut,3));
	 _ptail = (32.0*M_PI*_rho)/(9.0*pow(_rcut,9)) - (16.0*M_PI*_rho)/(3.0*pow(_rcut,3));
	 cout << "Tail correction for the potential energy = " << _vtail << endl;
	 cout << "Tail correction for the virial           = " << _ptail << endl; 
	 cout << "The program perform Metropolis moves with uniform translations" << endl;
	 cout << "Moves parameter = " << _delta << endl;
	 cout << "Number of blocks = " << _nblk << endl;
	 cout << "Number of steps in one block = " << _nstep << endl << endl;

}




void MD::inputPos(const char*filename){
	ifstream input;
	input.open(filename);
	_posx = new double [_npart];
	_posy = new double [_npart];
	_posz = new double [_npart];
	for(int i=0;i<_npart;i++){
		input>>_posx[i]>>_posy[i]>>_posz[i];
		_posx[i]=_posx[i]*_box;
		_posy[i]=_posy[i]*_box;
		_posz[i]=_posz[i]*_box;
	}
	input.close();
	return;
}




void MD::inputV(){
	_vx = new double [_npart];
	_vy = new double [_npart];
	_vz = new double [_npart];
	_xold = new double [_npart];
	_yold = new double [_npart];
	_zold = new double [_npart];
	double sumv[3]={0.};
	double fs;
	for(int i=0;i<_npart;i++){
		_vx[i]=_rand->Rannyu()/double(RAND_MAX)-0.5;
		_vy[i]=_rand->Rannyu()/double(RAND_MAX)-0.5;
		_vz[i]=_rand->Rannyu()/double(RAND_MAX)-0.5;
		sumv[0]+=_vx[i];
		sumv[1]+=_vy[i];
		sumv[2]+=_vz[i];
	}
	for(int i=0;i<3;i++){sumv[i]/=double(_npart);}
	double sumv2=0.;
	for(int i=0;i<_npart;i++){
		_vx[i]=_vx[i]-sumv[0];
		_vy[i]=_vy[i]-sumv[1];
		_vz[i]=_vz[i]-sumv[2];
		sumv2+=pow(_vx[i],2.)+pow(_vy[i],2.)+pow(_vz[i],2.);
	}

	sumv2/=double(_npart);
	fs=sqrt(3.*_t/sumv2);
	for(int i=0;i<_npart;i++){
		_vx[i]*=fs;
		_vy[i]*=fs;
		_vz[i]*=fs;

		_xold[i]=pbc(_posx[i]-_vx[i]*_delta);
		_yold[i]=pbc(_posy[i]-_vy[i]*_delta);
		_zold[i]=pbc(_posz[i]-_vz[i]*_delta);
	}
	return;
}




double MD::pbc(double r){
	return r-_box*rint(r/_box);
}




double MD::LennJon(int i,int asse){
	double f=0.;
	double vett[3];
	double dr;
	for(int j=0;j<_npart;j++){
		if(j != i){
			vett[0]=pbc(_posx[i]-_posx[j]);
			vett[1]=pbc(_posy[i]-_posy[j]);
			vett[2]=pbc(_posz[i]-_posz[j]);

			dr=sqrt(pow(vett[0],2.)+pow(vett[1],2.)+pow(vett[2],2.));

			if(dr<_rcut){
				f+=vett[asse]*(48.0/pow(dr,14.)-24.0/pow(dr,8.));
			}
		}
	}
	return f;
}




void MD::Move(){
	double fx[_npart];
	double fy[_npart];
	double fz[_npart];
	double x,y,z;
	for(int i=0;i<_npart;i++){
		fx[i]=LennJon(i,0);
		fy[i]=LennJon(i,1);
		fz[i]=LennJon(i,2);
	}

	for(int i=0;i<_npart;i++){
		x=pbc(2.*_posx[i]-_xold[i]+fx[i]*pow(_delta,2.));
		y=pbc(2.*_posy[i]-_yold[i]+fy[i]*pow(_delta,2.));
		z=pbc(2.*_posz[i]-_zold[i]+fz[i]*pow(_delta,2.));
		
		_vx[i]=pbc(x-_xold[i])/(2.*_delta);
		_vy[i]=pbc(y-_yold[i])/(2.*_delta);
		_vz[i]=pbc(z-_zold[i])/(2.*_delta);
	
		_xold[i]=_posx[i];
		_yold[i]=_posy[i];
		_zold[i]=_posz[i];

		_posx[i]=x;
		_posy[i]=y;
		_posz[i]=z;
	}
	return;
}




void MD::Misura_blocchi(){
	//int bin;
	int bin;
  	double v = 0.0, w = 0.0;
  	double vij, wij;
  	double dx, dy, dz, dr;

	//reset the hystogram of g(r)
  	for (int k=_igofr; k<_igofr+_nbins; ++k){_walker[k]=0.0;}
		
	//cycle over pairs of particles
 	for (int i=0; i<_npart-1; ++i){
    		for (int j=i+1; j<_npart; ++j){          // distance i-j in pbc
     			dx = pbc(_posx[i] - _posx[j]);
			dy = pbc(_posy[i] - _posy[j]);
			dz = pbc(_posz[i] - _posz[j]);
			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);
	//CODICE PER RIEMPIRE L'ISTOGRAMMA DI G(R)
			if(dr<=_box*0.5){
				bin=0;
				for(float l=0;l+bin_size<dr;l+=bin_size){
					bin++;
				}
			//bin=int(dr/bin_size);
			_walker[_igofr+bin]+=2;
			}

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

	return;
}



void MD::ConfFinal(){
	ofstream WriteConf;
	cout << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("config.final");
	for (int i=0; i<_npart; ++i){
		WriteConf << _posx[i]/_box << "   " <<  _posy[i]/_box << "   " << _posz[i]/_box << endl;
	}
	WriteConf.close();
	return;
}



void MD::ConfXYZ(int nconf){
	ofstream WriteXYZ;
  	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  	// << _npart << endl;
  	//WriteXYZ << "This is only a comment!" << endl;
  	for (int i=0; i<_npart; ++i){
    		WriteXYZ << pbc(_posx[i]) << "   " <<  pbc(_posy[i]) << "   " << pbc(_posz[i]) << endl;
  	}
  	WriteXYZ.close();
	return;
}




void MD::inputPosRestart(){
	ifstream input;
	ifstream input2;
	input.open("old.0");
	input2.open("old.final");
	for(int i=0;i<_npart;i++){
		input>>_posx[i]>>_posy[i]>>_posz[i];
		input2>>_xold[i]>>_yold[i]>>_zold[i];	
		_posx[i]=_posx[i]*_box;
		_posy[i]=_posy[i]*_box;
		_posz[i]=_posz[i]*_box;
		_xold[i]=_xold[i]*_box;
		_yold[i]=_yold[i]*_box;
		_zold[i]=_zold[i]*_box;
	}
	input.close();
	input2.close();
	
	Move();
	
	double sumv2=0.;
	for(int i=0;i<_npart;i++){
		_vx[i]=pbc(_posx[i]-_xold[i])/(_delta);
		_vy[i]=pbc(_posy[i]-_yold[i])/(_delta);
		_vz[i]=pbc(_posz[i]-_zold[i])/(_delta);
		sumv2+=pow(_vx[i],2.)+pow(_vy[i],2.)+pow(_vz[i],2.);
	}
	
	sumv2/=double(_npart);
	double Tmis=sumv2/3.;
	double fs=sqrt(_t/Tmis);
	for(int i=0;i<_npart;i++){
		_vx[i]*=fs;
		_vy[i]*=fs;
		_vz[i]*=fs;
		_xold[i]=pbc(_posx[i]-_delta*_vx[i]);
		_yold[i]=pbc(_posy[i]-_delta*_vy[i]);
		_zold[i]=pbc(_posz[i]-_delta*_vz[i]);
	}
	return;
}




void MD::OldConf(){
	ofstream OldConf;
	OldConf.open("old.0");
	for(int i=0;i<_npart;i++){
		OldConf<<_posx[i]/_box<<" "<<_posy[i]/_box<<" "<<_posz[i]/_box<<endl;
	}
	OldConf.close();
	
	OldConf.open("old.final");
	for(int i=0;i<_npart;i++){
		OldConf<<_xold[i]/_box<<" "<<_yold[i]/_box<<" "<<_zold[i]/_box<<endl;
	}
	OldConf.close();
	return;
}




void MD:: Accumulate(){
	for(int i=0; i<n_props; ++i){
		blk_av[i] = blk_av[i] + _walker[i];
	}
	blk_norm = blk_norm + 1.0;
}





void MD::Reset(int iblk){
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


	return;
}





void MD::Averages(int iblk){
double r, gdir;
   	ofstream Gofr, Epot, Pres;
	 
	//cout << "Block number " << iblk << endl;
   	//cout << "Acceptance rate " << accepted/attempted << endl << endl;


	Epot.open("output.epot.4",ios::app);
    	Pres.open("output.pres.4",ios::app);
    	Gofr.open("output.gofr.4",ios::app);
    	
	double stima_pot = blk_av[_iv]/blk_norm/(double)_npart + _vtail; //Potential energy
    	glob_av[_iv] += stima_pot;
    	glob_av2[_iv] += stima_pot*stima_pot;
	
   	double err_pot=Error(glob_av[_iv],glob_av2[_iv],iblk);

	double stima_pres = _rho * _t + (blk_av[_iw]/blk_norm + _ptail * (double)_npart) / _vol; //Pressure
        glob_av[_iw] += stima_pres;
    	glob_av2[_iw] += stima_pres*stima_pres;
   	double err_press=Error(glob_av[_iw],glob_av2[_iw],iblk);

//Potential energy per particle
    	Epot <<iblk<<" "<< stima_pot<<" "<< glob_av[_iv]/(double)iblk <<" "<< err_pot << endl;
//Pressure
    	Pres << iblk <<" "<< stima_pres <<" "<< glob_av[_iw]/(double)iblk <<" "<<err_press << endl;
 //CALCOLO G(R)
	Gofr<<iblk<<" ";
	for(int i=0;i<_nbins;i++){
		double norm=0.;
		r=i*bin_size;
		norm=_rho*_npart*(4*M_PI*(pow(r+bin_size,3.)-pow(r,3.)))/3.;
		gdir=blk_av[_igofr+i]/blk_norm/norm;
		glob_av[_igofr+i]+=gdir;
		glob_av2[_igofr+i]+=gdir*gdir;
		Gofr<<gdir<<" ";	
	}
	Gofr<<endl;
	

	 //cout << "----------------------------" << endl << endl;

    	Epot.close();
    	Pres.close();
    	Gofr.close();
	return;
		
}




double MD::Error(double sum, double sum2, int iblk){
	return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}


void MD::Gavefinal(){
	ofstream Gave;
	Gave.open("output.gave.4",ios::app);
	double err=0,r=0;
	for(int i=0;i<_nbins;i++){
		r=i*bin_size;
		err=Error(glob_av[_igofr+i],glob_av2[_igofr+i],_nblk);
		Gave<<r<<" "<<glob_av[_igofr+i]/(double)_nblk<<" "<<err<<endl;
	}
	Gave.close();
}














		




















		
