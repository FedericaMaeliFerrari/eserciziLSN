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
	input>>_iprint;
	_vol=_npart/_rho;
	_box=pow(_vol,1./3.);

	input.close();
	_rand=r;
	//_temp=0.;
	//_ekin=0.;
	//_etot=0.;
	//_epot=0.;
	_it = 0; //Energy
  	_iek= 1; //Heat capacity
  	_iep = 2; //Magnetization
  	_iet= 3; //Magnetic susceptibility
	n_props = 4; //Number of observables

	_walker=new double [n_props];
	glob_av=new double [n_props];
	glob_av2=new double [n_props];
	blk_av=new double [n_props];

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




void MD::Misura(int k){
	//int bin;
	double v=0.;
	double e=0.;
	double vij;
	double dx,dy,dz,dr;
	ofstream Epot, Ekin, Etot, T;
	
	if(k==0){
		Epot.open("Potenziale.dat",fstream::app);
		Ekin.open("EnCinetica.dat",fstream::app);
		Etot.open("EnergiaTot.dat",fstream::app);
		T.open("Temperatura.dat",fstream::app);
	}
	else {
		Epot.open("Potenziale_Restart.dat",fstream::app);
		Ekin.open("EnCinetica_Restart.dat",fstream::app);
		Etot.open("EnergiaTot_Restart.dat",fstream::app);
		T.open("Temperatura_Restart.dat",fstream::app);	
	}


	for(int i=0;i<_npart-1;i++){
		for(int j=i+1;j<_npart;j++){
			dx=pbc(_posx[i]-_posx[j]);
			dy=pbc(_posy[i]-_posy[j]);
			dz=pbc(_posz[i]-_posz[j]);
			dr=sqrt(pow(dx,2.)+pow(dy,2.)+pow(dz,2.));
			if(dr<_rcut){
				vij=4./pow(dr,12.)-4./pow(dr,6.);
				v+=vij;     //potenziale
			}
		}
	}
	
	for(int i=0;i<_npart;i++){e+=0.5*(pow(_vx[i],2.)+pow(_vy[i],2.)+pow(_vz[i],2.));}     //energia cinetica
	
	double pot=v/(double)_npart;
	double ekin=e/(double)_npart;
	double t=2./3.*ekin;
	double etot=pot+ekin;
	

	Epot<<pot<<endl;
	Ekin<<ekin<<endl;
	T<<t<<endl;
	Etot<<etot<<endl;

	Epot.close();
	Ekin.close();
	T.close();
	Etot.close();

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



void MD::Misura_blocchi(){

	double v=0.;
	double e=0.;
	double vij;
	double dx,dy,dz,dr;
	
	for(int i=0;i<_npart-1;i++){
		for(int j=i+1;j<_npart;j++){
			dx=pbc(_posx[i]-_posx[j]);
			dy=pbc(_posy[i]-_posy[j]);
			dz=pbc(_posz[i]-_posz[j]);
			dr=sqrt(pow(dx,2.)+pow(dy,2.)+pow(dz,2.));
			if(dr<_rcut){
				vij=4./pow(dr,12.)-4./pow(dr,6.);
				v+=vij;     //potenziale
			}
		}
	}

	for(int i=0;i<_npart;i++){e+=0.5*(pow(_vx[i],2.)+pow(_vy[i],2.)+pow(_vz[i],2.));}     //energia cinetica

	//_epot=v/(double)_npart;
	//_ekin=e/(double)_npart;
	//_temp=2./3.*_ekin;
	//_etot=_epot+_ekin;
	_walker[_it]=2./3.*e/(double)_npart;
	_walker[_iek]=e/(double)_npart;
	_walker[_iep]=v/(double)_npart;
	_walker[_iet]=e/(double)_npart+v/(double)_npart;

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
	ofstream Temp, Ek, Ep, Etot;
	//cout << "Block number " << iblk << endl;
    	//cout << "Acceptance rate " << accepted/attempted << endl << endl;

	Temp.open("output.Temp.0",ios::app);
    	double stima_t = blk_av[_it]/blk_norm; //Energy
    	glob_av[_it]  += stima_t;
    	glob_av2[_it] += stima_t*stima_t;
    	double err_t=Error(glob_av[_it],glob_av2[_it],iblk);
    	Temp<<  iblk <<"  "<< stima_t <<" "<< glob_av[_it]/(double)iblk <<" "<< err_t <<endl;
    	Temp.close();
	

	Ek.open("output.Ekin.0",ios::app);
	double stima_ek=blk_av[_iek]/blk_norm; //Energy^2
	glob_av[_iek]  += stima_ek;
    	glob_av2[_iek] += stima_ek*stima_ek;
	double err_ek=Error(glob_av[_iek],glob_av2[_iek],iblk);
	Ek<<  iblk <<"  "<< stima_ek <<" "<< glob_av[_iek]/(double)iblk <<" "<< err_ek <<endl;
	Ek.close();


	Ep.open("output.Epot.0",ios::app);
	double stima_ep=blk_av[_iep]/blk_norm;  //Suscettibilità magnetica
	glob_av[_iep]  += stima_ep;
    	glob_av2[_iep] += stima_ep*stima_ep;
	double err_ep=Error(glob_av[_iep],glob_av2[_iep],iblk);
	Ep<<  iblk <<"  "<< stima_ep <<" "<< glob_av[_iep]/(double)iblk <<" "<< err_ep <<endl;
	Ep.close();

	Etot.open("output.Etot.0",ios::app);
	double stima_et=blk_av[_iet]/blk_norm;
	glob_av[_iet]  += stima_et;
    	glob_av2[_iet] += stima_et*stima_et;
	double err_et=Error(glob_av[_iet],glob_av2[_iet],iblk);
	Etot<<  iblk <<"  "<< stima_et <<" "<< glob_av[_iet]/(double)iblk <<" "<< err_et <<endl;
	Etot.close();
	

		
}




double MD::Error(double sum, double sum2, int iblk){
	return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}















		




















		
