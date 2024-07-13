#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
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
   
// Ora voglio calcolare al tempo t=0 attraverso Monte Carlo i prezzi europei Call-option e Put-option campionando direttamente il prezzo finale dell'asset S(T) per un GBM con r,s >0

   int M = 10000;   // M rappresenta il numero di "lanci"
   int N = 100;   // N rappresenta il numero di blocchi in cui divido M
   int L = M/N;   // L rappresenta il numero di lanci in ogni blocco
   std::vector<double> aveC(N);   // vettore che contiene i valori medi per ogni blocco per Call-option
   std::vector<double> aveP(N);   // vettore che contiene i valori medi per ogni blocco per Put-option
   std::vector<double> av2C(N);   // vettore che contiene i valori quadratici medi per ogni blocco per Call-option
   std::vector<double> av2P(N);   // vettore che contiene i valori quadratici medi per ogni blocco per Put-option
   std::vector<double> sum_prog(N);   // vettore che contiene la media cumulativa per ogni blocco
   std::vector<double> su2_prog(N);   // vettore che contiene la media quadratica cumulativa per ogni blocco
   std::vector<double> err_prog(N);   // vettore che contiene l'incertezza statistica per ogni blocco
   double S_0 = 100.;   // prezzo dell'asset a t=0
   double T = 1.;   // tempo di consegna
   double K = 100.;   // prezzo d'esercizio
   double r = 0.1;   // tasso d'interesse a rischio libero
   double s = 0.25;   // volatilità  
   
   for(int i=0; i<N; i++){   // per ogni blocco
      double sumC = 0;
      double sumP = 0;
      for(int j=0; j<L; j++){   // per ogni lancio (all'interno di un blocco)
      	 double x = rnd.Gauss(0., 1.);
      	 double S_T = asset(S_0, r, s, T, x);
      	 double C = exp(-r*T)*max(0., S_T-K);
      	 double P = exp(-r*T)*max(0., K-S_T);
      	 
         sumC += C;
         sumP += P;
      }
      aveC[i] = sumC/L;
      av2C[i] = aveC[i]*aveC[i];
      
      aveP[i] = sumP/L;
      av2P[i] = aveP[i]*aveP[i];
   }
   
   for(int i=0; i<N; i++){   // per ogni blocco
      for(int j=0; j<(i+1); j++){   // per ogni blocco fino all'i-esimo
         sum_prog[i] += aveC[j];
         su2_prog[i] += av2C[j];
      }
      sum_prog[i] /= (i+1);
      su2_prog[i] /= (i+1);
      err_prog[i] = error(sum_prog,su2_prog,i);
   }
   
   ofstream output_1;
   output_1.open("Es03_1C.csv");
   if (output_1.is_open()){
      output_1 << "CumulativeAverage" << "," << "StatisticalUncertainty" << endl;
      for(int i=0; i<N; i++){
         output_1 << sum_prog[i] << "," << err_prog[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es03_1C.out" << endl;
   output_1.close();
   
   
   
   for(int i=0; i<N; i++){   // per ogni blocco
      sum_prog[i] = 0;
      su2_prog[i] = 0;
      err_prog[i] = 0;
      
      for(int j=0; j<(i+1); j++){   // per ogni blocco fino all'i-esimo
         sum_prog[i] += aveP[j];
         su2_prog[i] += av2P[j];
      }
      sum_prog[i] /= (i+1);
      su2_prog[i] /= (i+1);
      err_prog[i] = error(sum_prog,su2_prog,i);
   }
   
   ofstream output_2;
   output_2.open("Es03_1P.csv");
   if (output_2.is_open()){
      output_2 << "CumulativeAverage" << "," << "StatisticalUncertainty" << endl;
      for(int i=0; i<N; i++){
         output_2 << sum_prog[i] << "," << err_prog[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es03_1P.out" << endl;
   output_2.close();
   
   rnd.SaveSeed();
   return 0;
}


