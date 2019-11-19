///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//   This program reproduces figures 6.5 (a), (b) and (c) in                         //
//                                                                                   //
//   Christof Koch. Biophysics of Computation - Information Processing               //
//   in Single Neurons. Oxford University Press, 2004                                // 
//                                                                                   //
//                                                                                   //
//   NEURONAL MODEL: Hodgkin Huxley with Euler's Method                              //
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
#define step 0.001

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
	public:
		double Euler(double, double,double,double);
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

        double Vo, Mo,No,Ho; 
        Vo = v;
        Mo = m;
        No = n;
        Ho = h;
	Vo = Vo + step*dot_v(v,m,n,h);
	Mo = Mo + step*dot_m(v,m,n,h);
	No = No + step*dot_n(v,m,n,h);
	Ho = Ho + step*dot_h(v,m,n,h);
        V = Vo;
        M = Mo;
        N = No;
        H = Ho;



}

int main()
{
	
	Neuron neuron;

	cout << "Insert the current I(pA):";
	cin >> neuron.I;

	ofstream file;
	file.open("figure6.5.csv");

	int iterations = 200000;
	double time = 0.0;

	//Como inicializar v, m, n e h?
	neuron.V = 0.0;
	neuron.M = 0.0;
	neuron.N = 0.0;
	neuron.H = 0.0;

	for(int i = 0; i <= iterations; i++)
	{
		time = step*i;
		neuron.Euler(neuron.V,neuron.M,neuron.N,neuron.H);
		file<<time<<","<<neuron.V<<","<<neuron.M<<","<<neuron.N<<","<<neuron.H<<endl;
	}

	file.close();
	return 0;
}
