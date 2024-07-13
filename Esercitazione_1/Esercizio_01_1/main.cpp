/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>

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
   
// Ora voglio testare questo generatore di numeri pseudo-casuali stimando il valore medio di r, ossia la media della distribuzione uniforme e la sua incertezza

   int M = 100000;   // M rappresenta il numero di "lanci"
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
         sum += k;
      }
      ave[i] = sum/L;
      av2[i] = ave[i]*ave[i];
   }   // A questo punto ho calcolato per ogni blocco il valore medio e il valore quadratico medio
   
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
   
   ofstream output;
   output.open("Es01_1_01.csv");
   if (output.is_open()){
      output << "Cumulative_Average" << "," << "Statistical_Uncertainty" << endl;
      for(int i=0; i<N; i++){
         output << sum_prog[i] << "," << err_prog[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es01_1_01.out" << endl;
   output.close();
   
   
   
   // Ora voglio testare questo generatore di numeri pseudo-casuali stimando la varianza sigma^2, ossia la varianza della distribuzione uniforme e la sua incertezza
   
   for(int i=0; i<N; i++){   // per ogni blocco
      ave[i] = 0;
      av2[i] = 0;
      
      double sum = 0;   
      for(int j=0; j<L; j++){   // per ogni lancio (all'interno di un blocco)
      	 double k = rnd.Rannyu();
         sum += (k - 0.5)*(k - 0.5);
      }
      ave[i] = sum/L;
      av2[i] = ave[i]*ave[i];
   }   // A questo punto ho calcolato per ogni blocco il valore medio e il valore quadratico medio della varianza.
   
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
   }   // A questo punto ho calcolato per ogni blocco la media cumulativa, la media quadratica cumulativa e l'incertezza statistica relative alla varianza
   
   // Ora dovrò stampare in un file di output la media cumulativa e l'incertezza statistica, legate alla varianza, per ogni blocco
   
   ofstream output_2;
   output_2.open("Es01_1_02.csv");
   if (output_2.is_open()){
      output_2 << "Cumulative_Average" << "," << "Statistical_Uncertainty" << endl;
      for(int i=0; i<N; i++){
         output_2 << sum_prog[i] << "," << err_prog[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es01_1_02.out" << endl;
   output_2.close();
   
   
   
   // In questo ultimo test io voglio stimare la probabilità che i numeri non sono tratti da una distribuzione uniforme
   
   // N=100 ora è il numero di sotto-intervalli identici in cui divido [0,1]
   int n = 10000;   // numero di "lanci"
   int E = n/N;   // numero di eventi attesi osservati per ogni sotto-intervallo dopo n lanci
   std::vector<int> counts(N);   // vettore in cui verranno inseriti i conteggi per ogni blocco
   std::vector<double> chi_square(N);   // vettore in cui verranno inseriti i valori di chi^quadrato
   
   for(int i=0; i<100; i++){
      double chi = 0;
      for(int j=0; j<n; j++){
         double k = rnd.Rannyu();
         double K = 100*k;
         counts[round(K)]++;
      }
      for(int p=0; p<N; p++){
         chi += ((counts[p] - E)*(counts[p] - E));
         counts[p] = 0;
      }
      chi/=E;
      chi_square[i] = chi;
   }
   
   // Ora dovrò stampare in un file di output i 100 valori di chi^quadrato contenuti nel vettore chi_square
   
   ofstream output_3;
   output_3.open("Es01_1_03.csv");
   if (output_3.is_open()){
      output_3 << "Chi_Square" << endl;
      for(int i=0; i<N; i++){
         output_3 << chi_square[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es01_1_03.out" << endl;
   output_3.close();
   
   rnd.SaveSeed();
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
