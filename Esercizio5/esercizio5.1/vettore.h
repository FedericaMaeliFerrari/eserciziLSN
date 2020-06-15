#ifndef __Vettore_h__
#define __Vettore_h__
class Vettore{
	public:
		Vettore();
		Vettore(int N);
		Vettore(const Vettore&);
		Vettore& operator=(const Vettore&);
		~Vettore();
		unsigned int GetN() const { return m_N;}
		void SetComponent(int, double);
		double GetComponent(int) const;
		double GetComp_media(int);
		double GetComp_err(int);
		double media();
		double sigma();
		void mediaByBlocchi();
		
	protected:
		int m_N;
		double* m_v;
		double* m_media;
		double* m_err;
};
	
#endif
