///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//   This program saves spikes for one neuron based in                               //
//                                                                                   //
//   Christof Koch. Biophysics of Computation - Information Processing               //
//   in Single Neurons. Oxford University Press, 2004                                // 
//                                                                                   //
//   NEURONAL MODEL: Hodgkin Huxley with Runge-Kutta's Method                        //
//                                                                                   //
//   COMPILE: g++ spike.cpp -o main -lm                                             // 
//   RUN: ./main                                                                     //  
//                                                                                   //
//   Written by Luana J. N. Ferreira, 2020.                                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

// step size
#define step 0.01

class Neuron
{
	private:
		// 2011 PRE -- condutance: mS -- potential: mV
		double Gna = 120.0; 
		double Gk = 36.0;
		double Gm = 0.3;
		double Vrest = 10.6;
		double Ena = 115.0;
		double Ek = -12;
		double P = step; // RK4 integration step
		double v1,c1,k1,r1;
		double v2,c2,k2,r2;
		double v3,c3,k3,r3;
		double v4,c4,k4,r4;
		double v0,m0,n0,h0;
		double Vmais,mmais,nmais,hmais;		
	public:
		double Euler(double, double,double,double);
		int Spike(int);
		void initializate();
		void RK4(void);
		double dot_v(double, double,double,double);
		double dot_m(double, double,double,double);
	  	double dot_n(double,double,double,double);
	  	double dot_h(double,double,double,double);
	  	double I;
		double Isynapse;
  		double V,N,M,H;
		double Vpre, Vpos, Vatual;
		double spike;               
};

double Neuron::dot_v(double v,double m,double n,double h)
{return Gna*(pow(m,3))*h*(Ena-v) + Gk*(pow(n,4))*(Ek-v) + Gm*(Vrest-v) + I/(9*(M_PI));}

double Neuron::dot_m(double v,double m,double n,double h)
{return ((25-v)/(10*(exp((25-v)/10)-1)))*(1-m) - 4*(exp(-v/18))*m;}	

double Neuron::dot_n(double v,double m,double n,double h)
{return ((10-v)/(100*(exp((10-v)/10)-1)))*(1-n) - 0.125*(exp(-v/80))*n;}

double Neuron::dot_h(double v,double m,double n,double h)
{return 0.07*(exp(-v/20))*(1-h) - h/(exp((30-v)/10)+1);}

double Neuron::Euler(double v,double m,double n,double h)
{        
    v0 = v;
    m0 = m;
    n0 = n;
    h0 = h;
	v0 = v0 + step*dot_v(v,m,n,h);
	m0 = m0 + step*dot_m(v,m,n,h);
	n0 = n0 + step*dot_n(v,m,n,h);
	h0 = h0 + step*dot_h(v,m,n,h);
	V = v0;
    M = m0;
    N = n0;
    H = h0;
}

void Neuron::RK4()
{
	v0 = V;
    m0 = M;
    n0 = N;
    h0 = H;	

  v1 = dot_v(V,M,N,H);  
  k1 = dot_m(V,M,N,H);
  c1 = dot_n(V,M,N,H);
  r1 = dot_h(V,M,N,H);
  Vmais = V+P*0.5*v1;  
  mmais = M+P*0.5*k1;
  nmais = N+P*0.5*c1;
  hmais = H+P*0.5*r1;

  v2 = dot_v(Vmais,mmais,nmais,hmais);  
  k2 = dot_m(Vmais,mmais,nmais,hmais);
  c2 = dot_n(Vmais,mmais,nmais,hmais);
  r2 = dot_h(Vmais,mmais,nmais,hmais);
  Vmais = V+P*0.5*v2;  
  mmais = M+P*0.5*k2;
  nmais = N+P*0.5*c2;
  hmais = H+P*0.5*r2;
  
  v3 = dot_v(Vmais,mmais,nmais,hmais);  
  k3 = dot_m(Vmais,mmais,nmais,hmais);
  c3 = dot_n(Vmais,mmais,nmais,hmais);
  r3 = dot_h(Vmais,mmais,nmais,hmais);
  Vmais = V+P*v3;
  nmais = N+P*c3;
  mmais = M+P*k3;
  hmais = H+P*r3;
  
  v4 = dot_v(Vmais,mmais,nmais,hmais);
  k4 = dot_m(Vmais,mmais,nmais,hmais);
  c4 = dot_n(Vmais,mmais,nmais,hmais);
  r4 = dot_h(Vmais,mmais,nmais,hmais);
  V = v0+(P/6.0)*(v1+2.0*v2+2.0*v3+v4);
  M = m0+(P/6.0)*(k1+2.0*k2+2.0*k3+k4);
  N = n0+(P/6.0)*(c1+2.0*c2+2.0*c3+c4);  
  H = h0+(P/6.0)*(r1+2.0*r2+2.0*r3+r4);
}

int Neuron::Spike(int i)
{
	Vpre = Vatual;
	Vatual = Vpos;
	Vpos = V;

	if(Vpre<Vatual && Vatual>Vpos && i>1)
		{return 1;}

	return 0;

}

void Neuron::initializate()
{
	//Como inicializar v, m, n e h?
	V = 0.0;
	M = 0.0;
	N = 0.0;
	H = 0.0;
	
	Vpre = 0.0;
	Vpos = 0.0;
	Vatual = 0.0;
}

/****************************MAIN STARTS HERE*******************************/

int main()
{
	
	Neuron neuron;
	cout << "This is running Hodgkin-Huxley Model"<< endl;
	cout << "Enter the current I(pA): 280 or 300"<< endl;
	cin >> neuron.I;	

	ofstream file;
	if(neuron.I == 300)
	   file.open("1neuronI300.csv");
	else if(neuron.I == 280)
		file.open("1neuronI280.csv");
    else
        file.open("1neuron.csv");
    

	int iterations = 40000;
	double time = 0.0;
	neuron.initializate();
	

	for(int i = 0; i <= iterations; i++)
	{
		time = step*i;		
		neuron.RK4();
		neuron.spike = neuron.Spike(i);

		//neuron.Euler(neuron.V,neuron.M,neuron.N,neuron.H);
		file<<time<<","<<neuron.V<<","<<neuron.M<<","<<neuron.N<<","<<neuron.H<<","<<neuron.spike<<endl;
	}
	cout << "End of execution!"<<endl;
	file.close();
	return 0;
}
