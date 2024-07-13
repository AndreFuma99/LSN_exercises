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
   
// L'obiettivo è quello di simulare l'esperimento di Buffon: stimare pi-greco dai lanci di un ago

   double l = 1.0;   // lunghezza dell'ago
   double d = 2.0;   // distanza tra le linee del piano orizzontale

   int M = 100000;   // M rappresenta il numero di trials ma anche il numero di lanci dell'ago
   int N = 100;   // N rappresenta il numero di blocchi in cui divido M
   int L = M/N;   // L rappresenta il numero di lanci in ogni blocco
   std::vector<double> ave(N);   // ave è un vettore che contiene i valori medi di ogni blocco
   std::vector<double> av2(N);   // av2 è un vettore che contiene i valori quadratici medi di ogni blocco
   std::vector<double> sum_prog(N);   // sum_prog è un vettore che contiene la media cumulativa per ogni blocco
   std::vector<double> su2_prog(N);   // su2_prog è un vettore che contiene la media quadratica cumulativa per ogni blocco
   std::vector<double> err_prog(N);   // err_prog è un vettore che contiene l'incertezza statistica per ogni blocco
   
   for(int i=0; i<N; i++){   // per ogni blocco
      double sum = 0;   
      for(int j=0; j<L; j++){   // per ogni trial (all'interno di un blocco)
         double counts = 0.;
      	 for(int k=0; k<M; k++){   // per ogni lancio dell'ago
      	    double c = rnd.Rannyu(0., d/2);   // genero una posizione casuale del centro dell'ago rispetto alla distanza tra le linee
      	    double theta = rnd.Rannyu(0., 3.14159/2);   // genero casualmente l'angolo tra l'ago e le linee parallele
      	    if(c<=(l*sin(theta))/2)   counts++;
      	 }
      	 double pi_greco = (2*l*M)/(counts*d);
         sum += pi_greco;
      }
      ave[i] = sum/L;
      av2[i] = ave[i]*ave[i];
   }
   
   for(int i=0; i<N; i++){   // per ogni blocco
      for(int j=0; j<(i+1); j++){   // per ogni blocco fino all'i-esimo
         sum_prog[i] += ave[j];
         su2_prog[i] += av2[j];
      }
      sum_prog[i] /= (i+1);
      su2_prog[i] /= (i+1);
      err_prog[i] = error(sum_prog,su2_prog,i);
   }
   
   ofstream output;
   output.open("Es01_3.csv");
   if (output.is_open()){
      output << "Cumulative_Average" << "," << "Statistical_Uncertainty" << endl;
      for(int i=0; i<N; i++){
         output << sum_prog[i] << "," << err_prog[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es01_3.csv" << endl;
   output.close();
   
   rnd.SaveSeed();
   return 0;
}


