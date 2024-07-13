#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "random.h"
#include <vector>

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   
// L'obiettivo è calcolare il valore di aspettazione dell'Hamiltoniana che deve essere usato per trovare i parametri che minimizzano questa quantità

   int M = 10000;   // M rappresenta il numero di "lanci"
   int N = 100;   // N rappresenta il numero di blocchi in cui divido M
   int L = M/N;   // L rappresenta il numero di lanci in ogni blocco
   vector<double> ave(N, 0.);
   vector<double> av2(N, 0.);
   vector<double> sum_prog(N, 0.);
   vector<double> su2_prog(N, 0.);
   vector<double> err_prog(N, 0.);
   
   
   double x_n = 1.;   // punto di partenza
   double step = 0.8;   // step della probabilità di transizione uniforme
   double mu = 0.75;
   double sigma = 0.7;
   // mu e sigma sono i due parametri variazionali da ottimizzare
   
   // Prima di tutto devo trovare uno step della probabilità di transizione uniforme T: la probabilità di accettazione A è circa = 0.5;
   /*
   double SUM = 0;
   
   for(int k=0; k<M; k++){
      
      double p = rnd.Rannyu(-step,step);
      double x_p = x_n + p;
      	 
      double A = min(1., (psi_trial(x_p, mu, sigma)*psi_trial(x_p, mu, sigma))/(psi_trial(x_n, mu, sigma)*psi_trial(x_n, mu, sigma)));
      	 
      double s = rnd.Rannyu();
      if (s<=A) x_n = x_p;
      
      SUM += A;
      
   }
   cout << SUM/M << endl;
   */
   cout << "Utilizzando uno step della probabilità di transizione uniforme pari a " << step << ", la probabilità di accettazione è circa del 50%" << endl;


   for(int i=0; i<N; i++){
   
      double sum = 0;
      
      for(int j=0; j<L; j++){
      
         double x_p = rnd.Rannyu(x_n-step,x_n+step);
      	 
      	 double A = min(1., pow(psi_trial(x_p,mu,sigma),2)/pow(psi_trial(x_n,mu,sigma),2));
      	 
      	 double s = rnd.Rannyu();
      	 if (s<=A) x_n = x_p;
      	 
      	 sum += kinetic_energy(x_n, mu, sigma)/psi_trial(x_n, mu, sigma) + potential(x_n);
      }
      
      ave[i] = sum/L;
      av2[i] = ave[i]*ave[i];
      
   }
 
   
   for(int i=0; i<N; i++){
      for(int j=0; j<(i+1); j++){
         sum_prog[i] += ave[j];
         su2_prog[i] += av2[j];
      }
      sum_prog[i] /= (i+1);
      su2_prog[i] /= (i+1);
      err_prog[i] = error(sum_prog,su2_prog,i);
   }
   
   ofstream output;
   output.open("Es08_01.csv");
   if (output.is_open()){
      output << "CumulativeAverage" << "," << "StatisticalUncertainty" << endl;
      for(int i=0; i<N; i++){
         output << sum_prog[i] << "," << err_prog[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es08_01.csv" << endl;
   output.close();
   
   rnd.SaveSeed();
   return 0;
}


