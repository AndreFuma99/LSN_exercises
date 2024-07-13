#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

Random :: Random(){}
// Default constructor, does not perform any action

Random :: ~Random(){}
// Default destructor, does not perform any action

void Random :: SaveSeed(){
   // This function saves the current state of the random number generator to a file "seed.out"
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << "RANDOMSEED	" << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
   WriteSeed.close();
   return;
}

double Random :: Gauss(double mean, double sigma) {
   // This function generates a random number from a Gaussian distribution with given mean and sigma
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   // This function generates a random number in the range [min, max)
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  // This function generates a random number in the range [0,1)
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  // This function sets the seed and parameters of the random number generator
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

double error(vector <double> AV, vector <double> AV2, int n){
  if (n==0){
      return 0;
   } else return sqrt((AV2[n]-(AV[n]*AV[n]))/n);
}

double psi_trial(double x, double u, double s){
  return exp((-1.)*pow(x-u,2)/(2*pow(s,2)))+exp((-1.)*pow(x+u,2)/(2*pow(s,2)));
}

double potential(double x){
  return pow(x,4)-(5./2.)*pow(x,2);
}

double kinetic_energy(double x, double u, double s){
  return -0.5*((((x-u)*(x-u)-s*s)/pow(s,4))*exp((-1.)*(x-u)*(x-u)/(2.*s*s)) + (((x+u)*(x+u)-s*s)/pow(s,4))*exp((-1.)*(x+u)*(x+u)/(2.*s*s)));
}

pair<double, double> calculate_H(Random &rnd, double mu, double sigma, int M, int N, double step) {
    int L = M / N;
    vector<double> ave(N, 0.);
    vector<double> av2(N, 0.);
    vector<double> sum_prog(N, 0.);
    vector<double> su2_prog(N, 0.);
    vector<double> err_prog(N, 0.);
    
    double x_n = 1.0;
    
    // Equilibrazione
    for(int k=0; k<2000; k++){
      
         double x_p = rnd.Rannyu(x_n-step,x_n+step);
      	 
      	 double A = min(1., pow(psi_trial(x_p,mu,sigma),2)/pow(psi_trial(x_n,mu,sigma),2));
      	 
      	 double s = rnd.Rannyu();
      	 if (s<=A) x_n = x_p;
      }
    
    for(int i = 0; i < N; i++) {
        double sum = 0.0;
        for(int j = 0; j < L; j++) {
            double x_p = rnd.Rannyu(x_n - step, x_n + step);
            double A = min(1.0, pow(psi_trial(x_p, mu, sigma), 2) / pow(psi_trial(x_n, mu, sigma), 2));
            double s = rnd.Rannyu();
            if (s <= A) x_n = x_p;
            sum += kinetic_energy(x_n, mu, sigma) + potential(x_n);
        }
        ave[i] = sum / L;
        av2[i] = ave[i] * ave[i];
    }
    
    for(int i = 0; i < N; i++) {
        for(int j = 0; j <= i; j++) {
            sum_prog[i] += ave[j];
            su2_prog[i] += av2[j];
        }
        sum_prog[i] /= (i + 1);
        su2_prog[i] /= (i + 1);
        err_prog[i] = error(sum_prog, su2_prog, i);
    }
    
    return make_pair(sum_prog[N-1], err_prog[N-1]);
}


