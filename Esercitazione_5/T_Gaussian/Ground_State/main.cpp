#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"

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
   
// Ora utilizzo l'algoritmo di Metropolis per campionare il modulo quadro della funzione d'onda relativa allo stato fondamentale dell'atomo di H utilizzando una probabilità di transizione gaussiana   

   int M = 1000000;   // M rappresenta il numero di "lanci"
   int N = 100;   // N rappresenta il numero di blocchi in cui divido M
   int L = M/N;   // L rappresenta il numero di lanci in ogni blocco
   int dim = 100000;   // numero di iterazioni che utilizzo per l'equilibrazione
   std::vector<double> ave(N);   // ave è un vettore che contiene i valori medi di ogni blocco
   std::vector<double> av2(N);   // av2 è un vettore che contiene i valori quadratici medi di ogni blocco
   std::vector<double> sum_prog(N);   // sum_prog è un vettore che contiene la media cumulativa per ogni blocco
   std::vector<double> su2_prog(N);   // su2_prog è un vettore che contiene la media quadratica cumulativa per ogni blocco
   std::vector<double> err_prog(N);   // err_prog è un vettore che contiene l'incertezza statistica per ogni blocco
   std::vector<double> R(dim);   // R è un vettore di distanze che utilizzo durante la fase di equilibrazione
   std::vector<double> X(M);   // X è un vettore di ascisse che utilizzo per visualizzare il campionamento
   std::vector<double> Y(M);   // Y è un vettore di ordinate che utilizzo per visualizzare il campionamento
   std::vector<double> Z(M);   // Z è un vettore di altezze che utilizzo per visualizzare il campionamento
   
   double x_n = sqrt(3.0);   // ascissa del punto di partenza
   double y_n = sqrt(3.0);   // ordinata del punto di partenza
   double z_n = sqrt(3.0);   // altezza del punto di partenza
   double step = 0.755;   // step della probabilità di transizione gaussiana T(x|y)
   
   
   // Prima di tutto devo trovare uno step della probabilità di transizione gaussiana T(x|y) : la probabilità di accettazione A è circa = 0.5;
   /*
   double SUM = 0;
   
   for(int k=0; k<dim; k++){
      
      double p = rnd.Gauss(0.,step);
      double x_p = x_n + p;
      double q = rnd.Gauss(0.,step);
      double y_p = y_n + q;
      double r = rnd.Gauss(0.,step);
      double z_p = z_n + r;
      	 
      double A = min(1., psi_100(x_p, y_p, z_p)/psi_100(x_n, y_n, z_n));
      	 
      double s = rnd.Rannyu();
      if (s<=A){
      	  x_n = x_p;
      	  y_n = y_p;
      	  z_n = z_p;
      }
      
      SUM += A;
      
   }
   */
   cout << "Utilizzando uno step della probabilità di transizione gaussiana pari a " << step << ", la probabilità di accettazione è circa del 50%" << endl;

   for(int k=0; k<dim; k++){   // Questa è la fase di equilibrazione
      
      double p = rnd.Gauss(0.,step);
      double x_p = x_n + p;
      double q = rnd.Gauss(0.,step);
      double y_p = y_n + q;
      double r = rnd.Gauss(0.,step);
      double z_p = z_n + r;
      	 
      double A = min(1., psi_100(x_p, y_p, z_p)/psi_100(x_n, y_n, z_n));
      	 
      double s = rnd.Rannyu();
      if (s<=A){
      	  x_n = x_p;
      	  y_n = y_p;
      	  z_n = z_p;
      }
      
      R[k] = sqrt(x_n*x_n + y_n*y_n + z_n*z_n);
      
   }

   ofstream output_1;
   output_1.open("Es05_E0.csv");
   if (output_1.is_open()){
      output_1 << "Radius" << endl;
      for(int i=0; i<dim; i++){
         output_1 << R[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es05_E0.csv" << endl;
   output_1.close();

   for(int i=0; i<N; i++){
   
      double sum = 0;
      
      for(int j=0; j<L; j++){
      
         double p = rnd.Gauss(0.,step);
      	 double x_p = x_n + p;
      	 double q = rnd.Gauss(0.,step);
      	 double y_p = y_n + q;
      	 double r = rnd.Gauss(0.,step);
      	 double z_p = z_n + r;
      	 
      	 double A = min(1., psi_100(x_p, y_p, z_p)/psi_100(x_n, y_n, z_n));
      	 
      	 double s = rnd.Rannyu();
      	 if (s<=A){
      	    x_n = x_p;
      	    y_n = y_p;
      	    z_n = z_p;
      	 }
      	 
         sum += sqrt(x_n*x_n + y_n*y_n + z_n*z_n);
         
         X[i*L+j] = x_n;
         Y[i*L+j] = y_n;
         Z[i*L+j] = z_n;
      }
      
      ave[i] = sum/L;
      av2[i] = ave[i]*ave[i];
      
   }
   
   ofstream coordinates_1;
   coordinates_1.open("Es05_N0_xyz.csv");
   if (coordinates_1.is_open()){
      coordinates_1 << "X" << "," << "Y" << "," << "Z" << endl;
      for(int i=0; i<M; i++){
         coordinates_1 << X[i] << "," << Y[i] << "," << Z[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es05_N0_xyz.csv" << endl;
   coordinates_1.close();
   
   for(int i=0; i<N; i++){
      for(int j=0; j<(i+1); j++){
         sum_prog[i] += ave[j];
         su2_prog[i] += av2[j];
      }
      sum_prog[i] /= (i+1);
      su2_prog[i] /= (i+1);
      err_prog[i] = error(sum_prog,su2_prog,i);
   }
   
   ofstream output_2;
   output_2.open("Es05_N0.csv");
   if (output_2.is_open()){
      output_2 << "CumulativeAverage" << "," << "StatisticalUncertainty" << endl;
      for(int i=0; i<N; i++){
         output_2 << sum_prog[i] << "," << err_prog[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es05_N0.csv" << endl;
   output_2.close();

   return 0;
}


