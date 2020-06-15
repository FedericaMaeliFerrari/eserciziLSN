#include "vettore.h"
#include <iostream>
#include <cmath>
using namespace std;
Vettore::Vettore(){
	m_N=0;
	//*m_v=NULL;
	//*m_media=NULL;
}

Vettore::Vettore(int N){
	m_N=N;
	m_v = new double [N];
	m_media=new double [N];
	m_err=new double [N];
	//*m_v=NULL;
}

Vettore::~Vettore(){
	delete[] m_v;
}

void Vettore::SetComponent(int n, double a){
	if(n>m_N) cout << "non esiste questa componente" <<endl;
	else m_v[n]=a;
}

double Vettore::GetComponent(int n) const { 
	if(n>=m_N) cout<<"non esiste questa componente"<<endl;
	return m_v[n];
}

double Vettore::GetComp_media(int n){
	if(n>=m_N) cout<<"non esiste questa componente"<<endl;
	return m_media[n];
}

double Vettore::GetComp_err(int n){
	if(n>=m_N) cout<<"non esiste questa componente"<<endl;
	return m_err[n];
}

Vettore::Vettore(const Vettore &v){ //copy constructor
	m_N = v.m_N;
	m_v = new double [m_N];
	for(int i = 0;i<m_N; i++)	m_v[i]=v.m_v[i];
}

Vettore& Vettore::operator=(const Vettore &v){
	m_N = v.m_N;
	if(m_v) delete[] m_v;
	m_v = new double[m_N];
	for(int i = 0; i<m_N; i++) m_v[i]=v.m_v[i];
	return *this;
}

double Vettore::media(){
	double sum=0;
	for(int i=0;i<m_N;i++){
		sum=sum+m_v[i];
	}
	double m=sum/m_N;
	return m;
}

double Vettore::sigma(){
	double p=0;
	double m=media();
	for(int i=0;i<m_N;i++){
		p=p+pow((m_v[i]-m),2);
	}
	double s=pow((p/(m_N-1)),0.5);
	return s;
}

void Vettore::mediaByBlocchi(){
	double vett2[m_N]={0};
	double media2[m_N]={0};
	double media[m_N]={0};
	for(int i=0;i<m_N;i++){
		vett2[i]=pow(m_v[i],2.);
	}
	
	for(int i=0;i<m_N;i++){
		for(int j=0;j<i+1;j++){
			media[i]+=m_v[j];
			media2[i]+=vett2[j];
		}
		media[i]=media[i]/(i+1);
		m_media[i]=media[i];
		media2[i]=media2[i]/(i+1);
		if (i==0){
			m_err[i]=0;
		}
		else{
			m_err[i]=sqrt((media2[i]-pow(m_media[i],2.))/i);
		}
	}
	return;
}


	
			
	




