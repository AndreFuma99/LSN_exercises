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
   
// L'obiettivo è utilizzare l'algoritmo di Simulated Annealing per trovare i parametri che minimizzano il valore di aspettazione dell'Hamiltoniana

   int M = 10000;   // numero di "lanci"
   int N = 100;   // numero di blocchi
   double step = 0.8;   // step della probabilità di transizione uniforme
   double mu = 0.75;
   double sigma = 0.7;
   // mu e sigma sono i due parametri variazionali da ottimizzare
    
   double T = 2.0;   // temperatura iniziale
   double alpha = 0.99;   // fattore di raffreddamento
   int N_SA = 700;   // numero massimo di iterazioni dell'algoritmo
   ofstream output;
   output.open("SA_optimization.csv");   // file di output in cui salvare i risultati di ogni passo dell'algoritmo
   if (!output.is_open()) {
      cerr << "PROBLEM: Unable to open SA_optimization.csv" << endl;
      return -1;
   }
   output << "Step" << "," << "Mu" << "," << "Sigma" << "," << "Hamiltonian" << "," << "Error" << endl;
   
   double mu_min = mu;
   double sigma_min = sigma;   // queste sono tre variabili in cui salverò il set di parametri ottimali
   double H_min = 2.0;
   
   for (int i = 0; i < N_SA; i++) {
      auto [H, error] = calculate_H(rnd, mu, sigma, M, N, step);
      output << i << "," << mu << "," << sigma << "," << H << "," << error << endl;
        
      double new_mu = mu + rnd.Rannyu(-0.1, 0.1);   // propongo un nuovo valore di mu
      double new_sigma = sigma + rnd.Rannyu(-0.1, 0.1);   // propongo un nuovo valore di sigma
      auto [new_H, new_error] = calculate_H(rnd, new_mu, new_sigma, M, N, step);   // H ed errore per i nuovi parametri
        
      double delta_H = new_H - H;
      // l'algoritmo accetta il nuovo set di parametri se il nuovo valore di H è minore del precedente altrimenti accetta con una probabilità che dipende dalla temperatura
      if (delta_H < 0 || rnd.Rannyu() < exp(-delta_H / T)) {
         mu = new_mu;
         sigma = new_sigma;
         H = new_H;
         error = new_error;
      }
      T *= alpha;   // riduco la temperatura
      
      if (H < H_min){   // aggiorno i parametri ottimali solo per valori di energia più piccoli
         mu_min = mu;
         sigma_min = sigma;
         H_min = H;
      }
   }
   
   output.close();
   
   cout << "Valore ottimale di mu = " << mu_min << endl;
   cout << "Valore ottimale di sigma = " << sigma_min << endl;
   
   // Faccio andare un'ultima volta la simulazione per fare data blocking utilizzando il set di parametri ottimali
   
   int L = M/N;   // L rappresenta il numero di lanci in ogni blocco
   vector<double> ave(N, 0.);
   vector<double> av2(N, 0.);
   vector<double> sum_prog(N, 0.);
   vector<double> su2_prog(N, 0.);
   vector<double> err_prog(N, 0.);
   vector<double> mod_quad(M, 0.);   // vettore che riempirò con il modulo quadro della funzione d'onda di trial campionata
   
   double x_n = 1.;   // punto di partenza
   
   // Equilibrazione  
   for(int k=0; k<2000; k++){
      
      double x_p = rnd.Rannyu(x_n-step,x_n+step);
      	 
      double A = min(1., pow(psi_trial(x_p,mu_min,sigma_min),2)/pow(psi_trial(x_n,mu_min,sigma_min),2));
      	 
      double s = rnd.Rannyu();
      if (s<=A) x_n = x_p;
   }
   
   for(int i=0; i<N; i++){
   
      double sum = 0;
      
      for(int j=0; j<L; j++){
      
         double x_p = rnd.Rannyu(x_n-step,x_n+step);
      	 
      	 double A = min(1., pow(psi_trial(x_p,mu_min,sigma_min),2)/pow(psi_trial(x_n,mu_min,sigma_min),2));
      	 
      	 double s = rnd.Rannyu();
      	 if (s<=A) x_n = x_p;
      	 
      	 sum += kinetic_energy(x_n, mu_min, sigma_min)/psi_trial(x_n, mu_min, sigma_min) + potential(x_n);
      	 mod_quad[i*L+j] = x_n;
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
   
   ofstream Output;
   Output.open("Es08_02.csv");
   if (Output.is_open()){
      Output << "CumulativeAverage" << "," << "StatisticalUncertainty" << endl;
      for(int i=0; i<N; i++){
         Output << sum_prog[i] << "," << err_prog[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es08_02.csv" << endl;
   Output.close();
   
   ofstream OutPut;
   OutPut.open("|psi|2.csv");
   if (OutPut.is_open()){
      for(int i=0; i<M; i++){
         OutPut << mod_quad[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open |psi|2.csv" << endl;
   OutPut.close();
   
   rnd.SaveSeed();
   return 0;
}


