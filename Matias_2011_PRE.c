///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//   This program simulates the 3-neuron motif described in                           //
//                                                                                   //
//   Matias, F. S., Carelli, P. V., Mirasso, C. R., & Copelli, M. (2011).            //
//   Anticipated synchronization in a biologically plausible                         //
//   model of neuronal motifs. Physical Review E, 84(2), 021922.                     // 
//                                                                                   //
//   This reference should be cited whenever this code is used.                      //
//                                                                                   //
//   NEURONAL MODEL: Hodgkin Huxley                                                  //
//   SYNAPTIC MODEL: Excitatory (AMPA) and Inhibitory (GABA_A)                       //
//                                                                                   //
//   COMPILE: g++ -O3 Matias_2011_PRE.c -o Matias_2011_PRE -lm                       // 
//   RUN: ./Matias_2011_PRE gms gis gsi Im alpha_exc beta_exc alpha_inh beta_inh     //  
//   EX:  ./Matias_2011_PRE 10 10 10 280 1.1 0.19 5 0.3                              //
//                                                                                   //
//   Written by Fernanda S. Matias, 2011.                                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 5000000      // number of steps to integrate
#define P (double)0.01 // RK4 integration step
#define OUTPUT_STEP 1  // saves data every OUTPUT_STEP steps
#define nneuronios 4   // number of neurons

double Vm[nneuronios]; 
double time;
unsigned long int i;

//Hodgkin Huxley neuron model
double Ena= 115.0, Ek=-12.0; //mV
double an,bn,am,bm,ah,bh;


/////////////////// RAN2 Random number generator/////////////////////////////
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS) 

float ran2(long *idum)
{  
  int j;
  long k;
  static long idum2 = 123456789, iy = 0, iv[NTAB];
  float temp;
  
  if(*idum <= 0)
    {
      if(-(*idum) < 1)   *idum = 1;
      else   *idum = -(*idum);
      idum2 = (*idum);
      for(j = NTAB+7;j >= 0;j--)
	{
	  k = (*idum)/IQ1;
	  *idum = IA1*(*idum-k*IQ1)-k*IR1;
	  if(*idum < 0)   *idum += IM1;
	  if(j < NTAB)   iv[j] = *idum;
	}  
      iy = iv[0];
    }  
  k = (*idum)/IQ1;
  *idum = IA1*(*idum-k*IQ1)-k*IR1;
  if(*idum < 0)   *idum += IM1; 
  
  k = idum2/IQ2;
  idum2 = IA2*(idum2-k*IQ2)-k*IR2;
  if(idum2 < 0)   idum2 += IM2;
  
  j = iy/NDIV;
  iy = iv[j]-idum2;
  iv[j] = *idum;
  if(iy < 1)   iy += IMM1;
  
  if((temp = AM*iy) > RNMX)   return RNMX;
  else   return temp;
}  
//////////////////////////// RAN2 (ENDS) ////////////////////////////////////


//Hodgkin-Huxley neuronal model
class HH_neuron{
  //private
  double dot_V(double,double,double,double);
  double dot_n(double,double,double,double);
  double dot_m(double,double,double,double);
  double dot_h(double,double,double,double);
  double Vo,mo,no,ho,Vmais,mmais,nmais,hmais;
  double c1,c2,c3,c4,k1,k2,k3,k4,r1,r2,r3,r4,v1,v2,v3,v4;
  
 public:
  double I,I_syn;
  double V,n,m,h;
  void RK4_step(void);
  
};
void HH_neuron::RK4_step(){
  
  Vo=V;
  no=n;
  mo=m;
  ho=h;
  
  v1=dot_V(V,n,m,h);
  c1=dot_n(V,n,m,h);
  k1=dot_m(V,n,m,h);
  r1=dot_h(V,n,m,h);
  Vmais=V+P*0.5*v1;
  nmais=n+P*0.5*c1;
  mmais=m+P*0.5*k1;
  hmais=h+P*0.5*r1;
  
  v2=dot_V(Vmais,nmais,mmais,hmais);
  c2=dot_n(Vmais,nmais,mmais,hmais);
  k2=dot_m(Vmais,nmais,mmais,hmais);
  r2=dot_h(Vmais,nmais,mmais,hmais);
  Vmais=V+P*0.5*v2;
  nmais=n+P*0.5*c2;
  mmais=m+P*0.5*k2;
  hmais=h+P*0.5*r2;
  
  v3=dot_V(Vmais,nmais,mmais,hmais);
  c3=dot_n(Vmais,nmais,mmais,hmais);
  k3=dot_m(Vmais,nmais,mmais,hmais);
  r3=dot_h(Vmais,nmais,mmais,hmais);
  Vmais=V+P*v3;
  nmais=n+P*c3;
  mmais=m+P*k3;
  hmais=h+P*r3;
  
  v4=dot_V(Vmais,nmais,mmais,hmais);
  c4=dot_n(Vmais,nmais,mmais,hmais);
  k4=dot_m(Vmais,nmais,mmais,hmais);
  r4=dot_h(Vmais,nmais,mmais,hmais);
  V=Vo+(P/6.0)*(v1+2.0*v2+2.0*v3+v4);
  n=no+(P/6.0)*(c1+2.0*c2+2.0*c3+c4);
  m=mo+(P/6.0)*(k1+2.0*k2+2.0*k3+k4);
  h=ho+(P/6.0)*(r1+2.0*r2+2.0*r3+r4);
}

double HH_neuron::dot_V(double ve, double ene, double eme, double aga){
  return( 120*(eme*eme*eme)*aga*(Ena-ve)    
	  +36*(ene*ene*ene*ene)*(Ek-ve)
	  +0.3*(10.613-ve)
	  +(I + I_syn)/(9*3.1415)); 
}

double HH_neuron::dot_n(double ve, double ene, double eme, double aga){
  return( ((10-ve)/(100*(exp((10-ve)/10)-1)))*(1-ene)
	  -(0.125* (exp(-ve/80)))*ene );
}


double HH_neuron::dot_m(double ve, double ene, double eme, double aga){
  return( ((25-ve)/(10*(exp ((25-ve)/10)-1)))*(1-eme)
	  -(4* (exp(-ve/18)))*eme);
}

double HH_neuron::dot_h(double ve, double ene, double eme, double aga){
  return( 
	 (0.07*(exp (-ve/20)))*(1-aga)
	 -(1/(exp ((30-ve)/10)+1))*aga);
}

//Excitatory synapse (AMPA)
class synapseAMPA{
  //private
  double dot_s(double);
  double so,smais;
  double c1,c2,c3,c4;
  
 public:
  unsigned pre_cell,post_cell;
  double G_syn, V_syn;
  double s,alpha,beta;
  void parametros(void);
  void RK4_step(void);
  void update_syn(void);
  double get_I_syn(void);
  double s_infinity,T;
};

void synapseAMPA::parametros(){
  V_syn = 60.0; //Vrest= 0.0 if Vrest = -60.0
}

void synapseAMPA::RK4_step(){
  
  so=s;
  c1=dot_s(s);
  smais=s+P*0.5*c1;
  c2=dot_s(smais);
  smais=s+P*0.5*c2;
  c3=dot_s(smais);
  smais=s+P*c3;
  c4=dot_s(smais);
  s=so+(P/6.0)*(c1+2.0*c2+2.0*c3+c4);
}

void synapseAMPA::update_syn(){	
  T=1/(1+exp((62-Vm[pre_cell])/5));
  //  T=1/(1+exp((2-Vm[pre_cell])/2.5)); if Vrest = -60.0
  
}

double synapseAMPA::dot_s(double esse){
  return(alpha*T*(1-esse)-beta*esse );  
}

double synapseAMPA::get_I_syn(){
  return(G_syn*s*(V_syn-Vm[post_cell]));
}

//Inhibitory synapse
class synapseGABAa{
  //private
  double dot_s(double);
  double so,smais;
  double c1,c2,c3,c4;
  
  
 public:
  unsigned pre_cell,post_cell;
  double G_syn, V_syn;
  double s,alpha,beta;
  void parametros(void);	
  void RK4_step(void);
  void update_syn(void);
  double get_I_syn(void);
  double s_infinity,T;
};

void synapseGABAa::parametros(){
  V_syn = -20.0; //V_syn =-80; // if Vrest = -60.0 
}

void synapseGABAa::RK4_step(){
  
  so=s;
  c1=dot_s(s);
  smais=s+P*0.5*c1;
  c2=dot_s(smais);
  smais=s+P*0.5*c2;
  c3=dot_s(smais);
  smais=s+P*c3;
  c4=dot_s(smais);
  s=so+(P/6.0)*(c1+2.0*c2+2.0*c3+c4);
}

void synapseGABAa::update_syn(){	
   T=1/(1+exp((62-Vm[pre_cell])/5));
  // T=1/(1+exp((2-Vm[pre_cell])/2.5));	// if Vrest = -60.0 
  
}
double synapseGABAa::dot_s(double esse){
  return(alpha*T*(1-esse)-beta*esse );  
}

double synapseGABAa::get_I_syn(){
  return(G_syn*s*(V_syn-Vm[post_cell]));
}


//objects
HH_neuron neuron[nneuronios];
synapseAMPA  syn_0_1,syn_1_2, syn_3_0,syn_3_1,syn_3_2;
synapseGABAa syn_2_1;



////////////////////////////////////////////////////////////////
//                     MAIN STARTS HERE                       //
////////////////////////////////////////////////////////////////

int main(int argc,char *argv[]){

  double gms, gis, gsi, Im, alpha_exc, beta_exc, alpha_inh, beta_inh;
  double gMm=0.0, gMs=0.0, gMi=0.0;
  double IM = 0.0;
  double transiente =49900.0;

  //inputs
  gms=atof(argv[1]);
  gis=atof(argv[2]);
  gsi=atof(argv[3]);
  Im=atof(argv[4]);
  alpha_exc=atof(argv[5]);
  beta_exc=atof(argv[6]);
  alpha_inh=atof(argv[7]);
  beta_inh=atof(argv[8]);

  syn_0_1.G_syn = gms;
  syn_2_1.G_syn = gis;
  syn_1_2.G_syn = gsi;

  syn_0_1.alpha = alpha_exc;
  syn_0_1.beta  = beta_exc;
  syn_1_2.alpha = alpha_exc;
  syn_1_2.beta  = beta_exc;
  syn_2_1.alpha = alpha_inh;
  syn_2_1.beta  = beta_inh;
 
  syn_3_0.alpha = alpha_exc;
  syn_3_0.beta  = beta_exc;
  syn_3_1.alpha = alpha_exc;
  syn_3_1.beta  = beta_exc;
  syn_3_2.alpha = alpha_exc;
  syn_3_2.beta  = beta_exc;
 
  FILE *Vmemb, *I;
  char strVmemb[100], strI[100];

  //open output files
  sprintf(strVmemb,"ZVmemb_g%.1f_%.1f_%.1fI%.1f_Im%.1f_gM%.1f_%.1f_%.1fabab%.2f_%.2f_%.2f_%.2f.dat",gms,gis,gsi,IM,Im,gMm,gMs,gMi,syn_0_1.alpha,syn_0_1.beta,syn_2_1.alpha,syn_2_1.beta);
  remove(strVmemb);  
  Vmemb = fopen(strVmemb, "w");

  sprintf(strI,"ZIsyn_g%.1f_%.1f_%.1fI%.1f_Im%.1f_gM%.1f_%.1f_%.1fabab%.2f_%.2f_%.2f_%.2f.dat",gms,gis,gsi,IM,Im,gMm,gMs,gMi,syn_0_1.alpha,syn_0_1.beta,syn_2_1.alpha,syn_2_1.beta);
  remove(strI);  
  I = fopen(strI, "w");
  
  time=0.0;
  div_t q;
  long semente = 98458813;
  // Random number generator Warmup:

  for(i=0;i<1000;i++) ran2(&semente);
  
  //neurons parameters
  neuron[0].I = Im;
  neuron[1].I = Im; //Is;
  neuron[2].I = Im; //Ii;
  neuron[3].I = IM;
  
  //synaptic parameters
  syn_0_1.pre_cell=0;
  syn_0_1.post_cell=1;
  //syn_0_1.G_syn = gms;
  //syn_0_1.alpha = 1.1;
  //syn_0_1.beta = 0.19;  
   
  syn_1_2.pre_cell=1;
  syn_1_2.post_cell=2;
  //syn_1_2.G_syn = gsi;
  //syn_1_2.alpha = 1.1;
  //syn_1_2.beta = 0.19;
 
  syn_2_1.pre_cell=2;
  syn_2_1.post_cell=1;
  // syn_2_1.G_syn = gis;
  // syn_2_1.alpha = alpha_inh;
  // syn_2_1.beta = 0.18;

  syn_3_0.pre_cell=3;
  syn_3_0.post_cell=0;
  syn_3_0.G_syn = gMm;

  syn_3_1.pre_cell=3;
  syn_3_1.post_cell=1;
  syn_3_1.G_syn = gMs;
 
  syn_3_2.pre_cell=3;
  syn_3_2.post_cell=2;
  syn_3_2.G_syn = gMi;
  
  //initial conditions:
  //neurons
  for(i=0;i<nneuronios;i++){
    neuron[i].V=10.0*(ran2(&semente)-0.5);  
    
    an=((10-neuron[i].V)/(100*(exp((10-neuron[i].V)/10)-1)));  
    bn=(0.125* (exp(-neuron[i].V/80)));                  
    am=((25-neuron[i].V)/(10*(exp ((25-neuron[i].V)/10)-1)));
    bm=(4* (exp(-neuron[i].V/18))); 
    ah=(0.07*(exp (-neuron[i].V/20)));
    bh=(1/(exp ((30-neuron[i].V)/10)+1));
    
    neuron[i].n=(an/(an+bn))*ran2(&semente) ;   
    neuron[i].m=(am/(am+bm))*ran2(&semente) ;  
    neuron[i].h=(ah/(ah+bh))*ran2(&semente) ;
  }  

  neuron[0].I_syn=0.0;
  neuron[1].I_syn=0.0;
  neuron[2].I_syn=0.0;
  neuron[3].I_syn=0.0;

  //synapses
  syn_0_1.s=0.0;
  syn_1_2.s=0.0;
  syn_2_1.s=0.0;
  syn_3_0.s=0.0;
  syn_3_1.s=0.0;
  syn_3_2.s=0.0;
  
  syn_0_1.parametros();	
  syn_1_2.parametros();
  syn_2_1.parametros();
  syn_3_0.parametros();
  syn_3_1.parametros();
  syn_3_2.parametros();


  ///////////////////////////// Time integration ////////////////////////////
  for(i=0;i<N;i++){
    
    Vm[0]=neuron[0].V;
    Vm[1]=neuron[1].V;
    Vm[2]=neuron[2].V;
    Vm[3]=neuron[3].V;
    
    syn_0_1.update_syn();
    syn_0_1.RK4_step();		
    syn_1_2.update_syn();
    syn_1_2.RK4_step();
    syn_2_1.update_syn();
    syn_2_1.RK4_step();
    syn_3_0.update_syn();
    syn_3_0.RK4_step();
    syn_3_1.update_syn();
    syn_3_1.RK4_step();
    syn_3_2.update_syn();
    syn_3_2.RK4_step();
    
    neuron[0].I_syn=syn_3_0.get_I_syn();//syn_2_0.get_I_syn()+syn_1_0.get_I_syn();
    neuron[1].I_syn=syn_0_1.get_I_syn()+syn_2_1.get_I_syn()+syn_3_1.get_I_syn();
    neuron[2].I_syn=syn_1_2.get_I_syn()+syn_3_2.get_I_syn();
    neuron[3].I_syn=0.0;
    
    neuron[0].RK4_step();
    neuron[1].RK4_step();
    neuron[2].RK4_step();
    neuron[3].RK4_step();
	  
    time=(double)i*P;

    //outputs
    q=div(i,OUTPUT_STEP);
    if((q.rem==0)&&(time>transiente)) {
      fprintf(Vmemb,"%lf %lf %lf %lf %lf\n",
      time, Vm[0], Vm[1], Vm[2], Vm[3]);
   
      // fprintf(I,"%lf %lf %lf %lf %lf %lf %lf %lf %lf \n",time, syn_0_1.get_I_syn(),syn_1_2.get_I_syn(), syn_2_1.get_I_syn(),syn_0_1.s,syn_1_2.s,syn_2_1.s, syn_0_1.get_I_syn() + syn_2_1.get_I_syn() ,syn_3_1.get_I_syn(), syn_0_1.get_I_syn() + syn_2_1.get_I_syn() + syn_3_1.get_I_syn() );
   
    }
  }	       
  //////////////////// End of time integration /////////////////////	

  //close output files
  fclose(Vmemb);
  fclose(I);
  
  return(0);
} //END OF MAIN
