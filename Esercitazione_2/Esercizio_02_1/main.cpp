#include <iostream>
#include <fstream>
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
   
// Ora voglio calcolare un integrale 1D via Monte Carlo campionando una distribuzione uniforme nell'intervallo [0,1].

   int M = 10000;   // M rappresenta il numero di "lanci"
   int N = 100;   // N rappresenta il numero di blocchi in cui divido M
   int L = M/N;   // L rappresenta il numero di lanci in ogni blocco
   std::vector<double> ave(N);   // ave è un vettore che contiene i valori medi di ogni blocco
   std::vector<double> av2(N);   // av2 è un vettore che contiene i valori quadratici medi di ogni blocco
   std::vector<double> sum_prog(N);   // sum_prog è un vettore che contiene la media cumulativa per ogni blocco
   std::vector<double> su2_prog(N);   // su2_prog è un vettore che contiene la media quadratica cumulativa per ogni blocco
   std::vector<double> err_prog(N);   // err_prog è un vettore che contiene l'incertezza statistica per ogni blocco
   
   for(int i=0; i<N; i++){   // per ogni blocco
      double sum = 0;
      for(int j=0; j<L; j++){   // per ogni lancio (all'interno di un blocco)
      	 double k = rnd.Rannyu();
      	 my_function(k);
         sum += my_function(k);
      }
      ave[i] = sum/L;
      av2[i] = ave[i]*ave[i];
   }   // A questo punto ho calcolato per ogni blocco il valore medio e il valore quadratico medio di f(x).
   
   for(int i=0; i<N; i++){   // per ogni blocco
      for(int j=0; j<(i+1); j++){   // per ogni blocco fino all'i-esimo
         sum_prog[i] += ave[j];
         su2_prog[i] += av2[j];
      }
      sum_prog[i] /= (i+1);
      su2_prog[i] /= (i+1);
      err_prog[i] = error(sum_prog,su2_prog,i);
   }   // A questo punto ho calcolato per ogni blocco la media cumulativa, la media quadratica cumulativa e l'incertezza statistica
   
   // Ora dovrò stampare in un file di output la media cumulativa e l'incertezza statistica per ogni blocco
   
   ofstream output_1;
   output_1.open("Es02_1_01.csv");
   if (output_1.is_open()){
      output_1 << "Cumulative_Average" << "," << "Statistical_Uncertainty" << endl;
      for(int i=0; i<N; i++){
         output_1 << sum_prog[i] << "," << err_prog[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es02_1_01.csv" << endl;
   output_1.close();
   
// Ora voglio calcolare lo stesso integrale 1D via Monte Carlo utilizzando importance sampling.

   for(int i=0; i<N; i++){   // per ogni blocco
      ave[i] = 0;
      av2[i] = 0;
      double sum = 0;
      for(int j=0; j<L; j++){   // per ogni lancio (all'interno di un blocco)
      	 double k = rnd.Rannyu();
      	 double q = rnd.Rannyu();
      	 while (q>=(1-k*k)){   // per generare dei numeri casuali da p(x) utilizzo il rejection method
      	    k = rnd.Rannyu();
      	    q = rnd.Rannyu();
      	 }
      	 sum += MY_g(k);
      }
      ave[i] = sum/L;
      av2[i] = ave[i]*ave[i];
   }   // A questo punto ho calcolato per ogni blocco il valore medio e il valore quadratico medio di g(x)
   
   for(int i=0; i<N; i++){   // per ogni blocco
      sum_prog[i] = 0;
      su2_prog[i] = 0;
      err_prog[i] = 0;
      for(int j=0; j<(i+1); j++){   // per ogni blocco fino all'i-esimo
         sum_prog[i] += ave[j];
         su2_prog[i] += av2[j];
      }
      sum_prog[i] /= (i+1);
      su2_prog[i] /= (i+1);
      err_prog[i] = error(sum_prog,su2_prog,i);
   }   // A questo punto ho calcolato per ogni blocco la media cumulativa, la media quadrata cumulativa e l'incertezza statistica
   
   // Ora dovrò stampare in un file di output la media cumulativa e l'incertezza statistica per ogni blocco
   
   ofstream output_2;
   output_2.open("Es02_1_02.csv");
   if (output_2.is_open()){
      output_2 << "Cumulative_Average" << "," << "Statistical_Uncertainty" << endl;
      for(int i=0; i<N; i++){
         output_2 << sum_prog[i] << "," << err_prog[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es02_1_02.csv" << endl;
   output_2.close();
     
   rnd.SaveSeed();
   return 0;
}


