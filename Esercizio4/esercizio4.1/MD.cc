#include "MD.h"
#include<fstream>
#include<iostream>
#include<cmath>
#include "random.h"

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
	cout<<"sumv2: "<<sumv2<<endl;
	fs=sqrt(3.*_t/sumv2);
	cout<<"fs: "<<fs<<endl;
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
	cout<<"sumv2: "<<sumv2<<endl;
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
	cout<<"fs: "<<fs<<endl;
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
	
	

	










		




















		
