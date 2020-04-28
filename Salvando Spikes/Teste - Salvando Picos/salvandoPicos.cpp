#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
#define step 0.1

double spike(double i);
double pre = 0.0;
double atual = 0.0;
double pos = 0.0;

double spike(double value){
    pre = atual;
    atual = pos;
    pos = value;

    if(pre < atual && atual > pos)
        return 1;
    
    return 0;
}

int main(){
    ofstream file;
    double x, y;
    int i;
    file.open("seno.csv");    

    for(i=0; i<200; i++)
    {
        x = step*i;
        y = sin(x);

        //lembrar que o spike é calculado utilizando a posição posterior
        file<<x<<","<<y<<","<<spike(y)<<endl; 
        
    }
    cout << "Succesfully executed! :)"<<endl;
    file.close();    
    return 0;
}

