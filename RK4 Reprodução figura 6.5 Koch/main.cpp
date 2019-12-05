///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//   This program reproduces figures 6.5 (a), (b) and (c) in                         //
//                                                                                   //
//   Christof Koch. Biophysics of Computation - Information Processing               //
//   in Single Neurons. Oxford University Press, 2004                                // 
//                                                                                   //
//                                                                                   //
//   NEURONAL MODEL: Hodgkin Huxley with Runge-Kutta's Method                        //
//                                                                                   //
//   COMPILE: g++ main.cpp -o main -lm                                               // 
//   RUN: ./main                                                                     //  
//                                                                                   //
//   Written by Luana J. N. Ferreira, 2019.                                          //
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
	public:
		double Euler(double, double,double,double);
		void RK4(void);
		double dot_v(double, double,double,double);
		double dot_m(double, double,double,double);
	  	double dot_n(double,double,double,double);
	  	double dot_h(double,double,double,double);
	  	double I;
  		double V,N,M,H;
                
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
    double Vo = v;
    double Mo = m;
    double No = n;
    double Ho = h;
	Vo = Vo + step*dot_v(v,m,n,h);
	Mo = Mo + step*dot_m(v,m,n,h);
	No = No + step*dot_n(v,m,n,h);
	Ho = Ho + step*dot_h(v,m,n,h);
	V = Vo;
    M = Mo;
    N = No;
    H = Ho;
}

void Neuron::RK4()
{
	//implementar
	double Vo = V;
    double Mo = M;
    double No = N;
    double Ho = H;

	double v1,c1,k1,r1;
	double v2,c2,k2,r2;
	double v3,c3,k3,r3;
	double v4,c4,k4,r4;
	double v0,m0,n0,h0;

  v1 = dot_v(V,N,M,H);  
  k1 = dot_m(V,N,M,H);
  c1 = dot_n(V,N,M,H);
  r1 = dot_h(V,N,M,H);
  v0 = V+P*0.5*v1;  
  m0 = M+P*0.5*k1;
  n0 = N+P*0.5*c1;
  h0 = H+P*0.5*r1;
  
  v2 = dot_v(v0,m0,n0,h0);  
  k2 = dot_m(v0,m0,n0,h0);
  c2 = dot_n(v0,m0,n0,h0);
  r2 = dot_h(v0,m0,n0,h0);
  v0 = V+P*0.5*v2;  
  m0 = M+P*0.5*k2;
  n0 = N+P*0.5*c2;
  h0 = H+P*0.5*r2;
  
  v3 = dot_v(v0,m0,n0,h0);  
  k3 = dot_m(v0,m0,n0,h0);
  c3 = dot_n(v0,m0,n0,h0);
  r3 = dot_h(v0,m0,n0,h0);
  v0 = V+P*v3;
  n0 = N+P*c3;
  m0 = M+P*k3;
  h0 = H+P*r3;
  
  v4 = dot_v(v0,m0,n0,h0);
  k4 = dot_m(v0,m0,n0,h0);
  c4 = dot_n(v0,m0,n0,h0);
  r4 = dot_h(v0,m0,n0,h0);
  V = Vo+(P/6.0)*(v1+2.0*v2+2.0*v3+v4);
  M = Mo+(P/6.0)*(k1+2.0*k2+2.0*k3+k4);
  N = No+(P/6.0)*(c1+2.0*c2+2.0*c3+c4);  
  H = Ho+(P/6.0)*(r1+2.0*r2+2.0*r3+r4);
}

int main()
{
	
	Neuron neuron;
	cout << "This is running Hodgkin-Huxley Model"<< endl;
	cout << "The current I(pA): 280"<< endl;
	// cin >> neuron.I;
	neuron.I = 280;

	ofstream file;
	file.open("1neuronI280.csv");

	int iterations = 40000;
	double time = 0.0;

	//Como inicializar v, m, n e h?
	neuron.V = 0.0;
	neuron.M = 0.0;
	neuron.N = 0.0;
	neuron.H = 0.0;

	for(int i = 0; i <= iterations; i++)
	{
		time = step*i;
		neuron.RK4();
		//neuron.Euler(neuron.V,neuron.M,neuron.N,neuron.H);
		file<<time<<","<<neuron.V<<","<<neuron.M<<","<<neuron.N<<","<<neuron.H<<endl;
	}

	file.close();
	return 0;
}
